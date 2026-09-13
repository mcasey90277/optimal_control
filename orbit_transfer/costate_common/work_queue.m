function out = work_queue(action, qDir, varargin)
%% Purpose:
%
%   A DISK WORK QUEUE for campaign units, so that N workers pull the next
%   unclaimed unit instead of being handed static ranges.
%
%   Why this exists (2026-09-13, doc/CAMPAIGN_DISCIPLINE.md). The 24x24
%   library's departure axis was split by giving each worker a fixed list of
%   columns. Columns run from 90 minutes to over 5 hours, so workers with a
%   slow column queued behind another sat idle while their queue waited --
%   about five hours lost in one day. Worse, a worker killed mid-unit took
%   its whole remaining list with it, because nothing else could pick the
%   work up.
%
%   The unit here is whatever the caller says it is (a rib column, a catalog
%   entry, a thrust rung). The queue owns three facts per unit: is it done,
%   is someone working on it, and when did that someone last say so.
%
%  ATOMIC CLAIMING. A claim is a DIRECTORY created with mkdir, which fails
%  when it already exists -- there is no window in which two workers both
%  believe they own a unit. (MATLAB's fopen has no exclusive-create mode,
%  so the obvious 'wx' is not available; mkdir is the portable primitive.)
%  Inside it a 'beat' file carries the worker tag and its mtime is the
%  progress signal, because changing a file does not touch its directory.
%
%  STALE CLAIMS ARE RECLAIMED. A claim whose file has not been touched for
%  staleSec is treated as abandoned and may be taken by another worker: a
%  worker that is killed does not take its unit out of circulation. This is
%  the property that made a watchdog kill cost nine hours of column 13.
%
%  DONE IS AN ARTIFACT, NOT A FLAG. A unit is done when its OUTPUT FILE
%  exists. The queue never records completion separately, so it cannot
%  disagree with the results on disk -- the failure mode where a sidecar
%  says measured and the catalog says empty.
%
%  ATTEMPTS ARE COUNTED, AND WRITTEN BEFORE THE WORK. A unit that always
%  fails is released and immediately re-claimed, forever: the first
%  end-to-end test of this queue livelocked on exactly that. After maxAtt
%  attempts a unit is RETIRED and handed out no more. The count is written
%  when the claim is taken, not when the work returns, so a unit that hangs
%  hard enough to kill the process still spends an attempt -- otherwise the
%  next run retries it and hangs identically.
%
%% Inputs:
%
%  action                   char                    'init' | 'claim' |
%                                                   'beat' | 'release' |
%                                                   'status'
%  qDir                     char                    queue directory
%  varargin                 per action:
%    init(qDir, units, outFcn)   units [1 x n] ids; outFcn(id) -> the output
%                                path whose existence means DONE. Saved as a
%                                HANDLE: func2str would drop the variables it
%                                closes over.
%    claim(qDir, tag, opts)      -> out.id (NaN if nothing left), out.claim
%                                opts .staleSec [1800] .maxAtt [3]
%    beat(qDir, claimFile)       touch a claim so it is not reclaimed
%    release(qDir, claimFile)    give a unit back (failed, not done)
%    status(qDir)                -> out .done .running .todo .stale .units
%
%% Outputs:
%
%  out                      struct                  per action, see above
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

switch lower(action)
    case 'init'
        units = varargin{1};  outFcn = varargin{2};
        if ~isfolder(qDir), mkdir(qDir); end
        % the handle is SAVED, not func2str'd: func2str throws away an
        % anonymous function's captured workspace, so a rebuilt handle could
        % not see the queue directory it was closed over (caught by the test)
        S = struct('units', units(:).', 'outFcn', outFcn, 'created', char(datetime('now')));
        save(fullfile(qDir, 'queue.mat'), '-struct', 'S');
        % clear any claim left by a previous run: claims are per-run state,
        % the OUTPUTS are what persists
        delete(fullfile(qDir, '*.att'));                  % attempts are per-run
        old = dir(fullfile(qDir, '*.claim'));
        for k = 1:numel(old)
            p_ = fullfile(qDir, old(k).name);
            if isfolder(p_), rmdir(p_, 's'); else, delete(p_); end
        end
        out = S;

    case 'claim'
        tag = varargin{1};
        opts = struct();  if numel(varargin) > 1, opts = varargin{2}; end
        staleSec = fieldd(opts, 'staleSec', 1800);
        maxAtt = fieldd(opts, 'maxAtt', 3);
        Q = load(fullfile(qDir, 'queue.mat'));
        outFcn = Q.outFcn;
        out = struct('id', NaN, 'claim', '');
        for id = Q.units
            if isfile(outFcn(id)), continue, end              % DONE: artifact exists
            if attempts(qDir, id) >= maxAtt, continue, end    % RETIRED: never again
            cf = fullfile(qDir, sprintf('%d.claim', id));
            if isfolder(cf)
                if claimAge(cf) < staleSec, continue, end     % someone live owns it
                rmdir(cf, 's');                               % abandoned: reclaim
            end
            [okMk, ~, idMk] = mkdir(cf);                      % ATOMIC: create-or-exists
            if ~okMk || strcmp(idMk, 'MATLAB:MKDIR:DirectoryExists')
                continue                                      % lost the race, try the next
            end
            writeBeat(cf, tag);
            bumpAttempt(qDir, id);        % BEFORE the work, not after
            out.id = id;  out.claim = cf;
            return
        end

    case 'beat'
        cf = varargin{1};
        if ~isempty(cf) && isfolder(cf), writeBeat(cf, 'beat'); end
        out = [];

    case 'release'
        cf = varargin{1};
        if ~isempty(cf) && isfolder(cf), rmdir(cf, 's'); end
        out = [];

    case 'status'
        Q = load(fullfile(qDir, 'queue.mat'));
        outFcn = Q.outFcn;
        n = numel(Q.units);
        maxAtt = 3;  if numel(varargin) > 0 && isstruct(varargin{1}), maxAtt = fieldd(varargin{1}, 'maxAtt', 3); end
        isDone = false(1, n);  isRun = false(1, n);  isStale = false(1, n);  isRet = false(1, n);
        att = zeros(1, n);
        for k = 1:n
            id = Q.units(k);
            isDone(k) = isfile(outFcn(id));
            att(k) = attempts(qDir, id);
            isRet(k) = ~isDone(k) && att(k) >= maxAtt;
            cf = fullfile(qDir, sprintf('%d.claim', id));
            if ~isDone(k) && isfolder(cf)
                age = claimAge(cf);
                isRun(k) = age < 1800;  isStale(k) = ~isRun(k);
            end
        end
        todoMask = ~isDone & ~isRun & ~isStale & ~isRet;
        out = struct('units', Q.units, 'done', Q.units(isDone), 'running', Q.units(isRun), ...
                     'stale', Q.units(isStale & ~isRet), 'retired', Q.units(isRet), ...
                     'todo', Q.units(todoMask), 'attempts', att, ...
                     'nDone', nnz(isDone), 'nRunning', nnz(isRun), 'nStale', nnz(isStale & ~isRet), ...
                     'nRetired', nnz(isRet), 'nTodo', nnz(todoMask), 'n', n);

    otherwise
        error('work_queue:action', 'unknown action "%s"', action);
end
end

% ------------------------------------------------------------------------
function n = attempts(qDir, id)
% ATTEMPTS  How many times unit id has been claimed.  INPUTS: qDir; id.
% OUTPUTS: n.
f = fullfile(qDir, sprintf('%d.att', id));
if ~isfile(f), n = 0;  return, end
v = str2double(strtrim(fileread(f)));
if isnan(v), n = 0; else, n = v; end
end

% ------------------------------------------------------------------------
function bumpAttempt(qDir, id)
% BUMPATTEMPT  Record a claim BEFORE the work, so a unit that hangs hard
% enough to kill the process still spends an attempt.
% INPUTS: qDir; id.  OUTPUTS: none.
f = fullfile(qDir, sprintf('%d.att', id));
n = attempts(qDir, id) + 1;
fid = fopen(f, 'w');  if fid >= 0, fprintf(fid, '%d', n);  fclose(fid); end
end

% ------------------------------------------------------------------------
function writeBeat(cf, tag)
% WRITEBEAT  (Re)write the beat file inside a claim: its mtime is the
% progress signal.  INPUTS: cf; tag.  OUTPUTS: none.
fid = fopen(fullfile(cf, 'beat'), 'w');
if fid >= 0, fprintf(fid, '%s|%s\n', tag, char(datetime('now')));  fclose(fid); end
end

% ------------------------------------------------------------------------
function a = claimAge(cf)
% CLAIMAGE  Seconds since the claim last reported progress.  Reads the beat
% FILE, because modifying a file does not update its directory's mtime.
% INPUTS: cf.  OUTPUTS: a.
d = dir(fullfile(cf, 'beat'));
if isempty(d), d = dir(cf); end
if isempty(d), a = inf;  return, end
a = seconds(datetime('now') - datetime(d(1).datenum, 'ConvertFrom', 'datenum'));
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
