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
    case {'beat', 'done', 'fail'}
        tag = varargin{1};  msg = '';
        if numel(varargin) > 1, msg = varargin{2}; end
        f = fullfile(hbDir, [tag '.hb']);
        % ATOMIC: write beside, then move. Truncating the live file in place
        % let a reader see an empty record between the truncate and the
        % write -- and an empty record used to read as RUNNING.
        tmp = sprintf('%s.%d.tmp', f, matlabProcessID);
        fid = fopen(tmp, 'w');
        assert(fid >= 0, 'campaign_heartbeat:write', 'cannot write heartbeat %s', f);
        fprintf(fid, '%s|%s|%s\n', lower(action), char(datetime('now')), msg);
        fclose(fid);
        movefile(tmp, f);
        out = f;

    case 'read'
        tags = varargin{1};
        opts = struct();  if numel(varargin) > 1, opts = varargin{2}; end
        staleSec = fieldd(opts, 'staleSec', 1800);
        if ischar(tags), tags = {tags}; end
        out = repmat(struct('tag', '', 'state', '', 'age', NaN, 'msg', '', 'when', ''), 1, numel(tags));
        for k = 1:numel(tags)
            out(k).tag = tags{k};
            f = fullfile(hbDir, [tags{k} '.hb']);
            if ~isfile(f)
                % NO FILE IS NOT HEALTH. The whole point of this function.
                out(k).state = 'never';  out(k).msg = 'no heartbeat: the worker never started';
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
            kind = parts{1};
            if numel(parts) > 1, out(k).when = parts{2}; end
            if numel(parts) > 2, out(k).msg = strjoin(parts(3:end), '|'); end
            % ONLY a well-formed 'beat' record can mean running. Anything
            % else -- empty, truncated, a token this reader does not know --
            % is UNKNOWN and alarms; it must never pass as health.
            switch kind
                case 'done', out(k).state = 'done';
                case 'fail', out(k).state = 'failed';
                case 'beat'
                    if out(k).age > staleSec, out(k).state = 'stalled';
                    else,                      out(k).state = 'running'; end
                otherwise
                    out(k).state = 'unknown';
                    out(k).msg = sprintf('unreadable heartbeat record "%s"', txt);
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
