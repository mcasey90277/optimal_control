function out = campaign_heartbeat(action, hbDir, varargin)
%% Purpose:
%
%   HEARTBEATS as the only liveness signal for campaign workers, and the
%   five states a monitor must be able to tell apart.
%
%   Why this exists (2026-09-13). Two monitor bugs in one campaign, both of
%   which reported health that had not been verified:
%     * liveness was a `pgrep` on the worker's tag, but the tag lived in the
%       ENVIRONMENT and never appeared in the command line, so four working
%       processes were reported dead;
%     * a MISSING log file read as "alive, age 0 s", because the age
%       expression fell back to the current time -- so three workers that
%       never launched at all looked exactly like workers that had just
%       checked in.
%   Both failures share a cause: the monitor inferred a state it had no
%   evidence for. A heartbeat is evidence. No heartbeat is NOT evidence of
%   health, and this reports that as its own state.
%
%   THE FIVE STATES, which a caller must never collapse into two:
%     never   -- no heartbeat file: the worker did not start. THE STATE MY
%                MONITORS KEPT LOSING.
%     running -- beat within staleSec
%     stalled -- beat older than staleSec, no done marker
%     done    -- the worker wrote its own completion
%     failed  -- the worker wrote a failure
%
%% Inputs:
%
%  action                   char                    'beat' | 'done' |
%                                                   'fail' | 'read'
%  hbDir                    char                    heartbeat directory
%  varargin                 beat(hbDir, tag, msg)   msg char, free text
%                           done/fail(hbDir, tag, msg)
%                           read(hbDir, tags, opts) opts .staleSec [1800]
%
%% Outputs:
%
%  out                      struct array            for 'read', one per tag:
%                                                   .tag .state .age .msg
%                                                   .when
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if ~isfolder(hbDir), mkdir(hbDir); end
switch lower(action)
    case {'beat', 'done', 'fail', 'ready'}
        tag = varargin{1};  msg = '';
        if numel(varargin) > 1, msg = varargin{2}; end
        % READY is its own PERSISTENT file. It used to be a beat record,
        % which the worker overwrote milliseconds later with its first
        % claim, so a launcher polling every 2 s never saw it and reported
        % workers that had finished a whole small campaign as "exited
        % before ready" (test_campaign_processes, 2026-09-13).
        if strcmpi(action, 'ready'), f = fullfile(hbDir, [tag '.ready']);
        else,                              f = fullfile(hbDir, [tag '.hb']); end
        % ATOMIC: write beside, then move. Truncating the live file in place
        % let a reader see an empty record between the truncate and the
        % write -- and an empty record used to read as RUNNING.
        tmp = sprintf('%s.%s.tmp', f, char(java.util.UUID.randomUUID()));
        fid = fopen(tmp, 'w');
        assert(fid >= 0, 'campaign_heartbeat:write', 'cannot write heartbeat %s', f);
        fprintf(fid, '%s|%s|%d|%s\n', lower(action), char(datetime('now')), matlabProcessID, msg);
        fclose(fid);
        publish_atomic(tmp, f);
        out = f;

    case 'read'
        tags = varargin{1};
        opts = struct();  if numel(varargin) > 1, opts = varargin{2}; end
        staleSec = fieldd(opts, 'staleSec', 1800);
        if ischar(tags), tags = {tags}; end
        out = repmat(struct('tag', '', 'state', '', 'age', NaN, 'msg', '', 'when', '', 'pid', NaN, 'ready', false), 1, numel(tags));
        for k = 1:numel(tags)
            out(k).tag = tags{k};
            f = fullfile(hbDir, [tags{k} '.hb']);
            if ~isfile(f)
                % NO FILE IS NOT HEALTH. The whole point of this function.
                out(k).state = 'never';  out(k).msg = 'no check-in observed (never wrote a heartbeat)';
                continue
            end
            d = dir(f);
            out(k).age = seconds(datetime('now') - datetime(d.datenum, 'ConvertFrom', 'datenum'));
            try
                txt = strtrim(fileread(f));
            catch
                txt = '';
            end
            parts = strsplit(txt, '|');
            % ONLY a complete, well-formed record can mean anything:
            % kind|time|pid|msg with a known kind, a parseable time and a
            % numeric pid. A truncated 'beat' alone, an empty file, or a
            % record from another format is UNKNOWN and alarms.
            kind = parts{1};  wellFormed = numel(parts) >= 4 && any(strcmp(kind, {'beat', 'done', 'fail'}));
            out(k).ready = isfile(fullfile(hbDir, [tags{k} '.ready']));
            if wellFormed
                out(k).when = parts{2};  out(k).pid = str2double(parts{3});
                out(k).msg = strjoin(parts(4:end), '|');
                wellFormed = ~isnan(out(k).pid) && ~isempty(regexp(parts{2}, '^\d{2}-\w{3}-\d{4} \d{2}:\d{2}:\d{2}$', 'once'));
            end
            if ~wellFormed
                out(k).state = 'unknown';
                out(k).msg = sprintf('unreadable heartbeat record "%s"', txt);
                continue
            end
            switch kind
                case 'done', out(k).state = 'done';
                case 'fail', out(k).state = 'failed';
                case 'beat'
                    if out(k).age > staleSec, out(k).state = 'stalled';
                    else,                      out(k).state = 'running'; end
            end
        end

    otherwise
        error('campaign_heartbeat:action', 'unknown action "%s"', action);
end
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
