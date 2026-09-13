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
%  Inside it a 'beat' file carries the OWNER TOKEN and its mtime is the
%  progress signal, because changing a file does not touch its directory.
%
%  OWNERSHIP IS A TOKEN, NOT A PATH. Every claim is issued with a random
%  token; beat and release act only if the claim still carries that token.
%  Without it a worker that was reclaimed while blocked in a solver would,
%  on returning, refresh the NEW owner's claim and then delete it.
%  (Astra chain review 2026-09-13.)
%
%  STALE CLAIMS ARE RECLAIMED, BY ONE RECLAIMER. A claim whose beat has not
%  been touched for staleSec is treated as abandoned. Takeover is a RENAME
%  of the stale directory -- atomic, so of two workers that both saw it
%  stale exactly one wins the rename and the other moves on. A plain
%  stat-then-rmdir let both "reclaim" and both compute the unit.
%
%  DONE IS AN ARTIFACT, NOT A FLAG. A unit is done when its OUTPUT FILE
%  exists. The queue never records completion separately, so it cannot
%  disagree with the results on disk. Callers publish that file ATOMICALLY
%  (write to a temp name, then move), so a partial save is never mistaken
%  for a finished unit. After a claim is won the done and attempt checks
%  are made AGAIN, because they were made before ownership and the world
%  may have moved.
%
%  ATTEMPTS ARE COUNTED, WRITTEN BEFORE THE WORK, AND PERSIST. A unit that
%  always fails is released and immediately re-claimed, forever: the first
%  end-to-end test of this queue livelocked on exactly that. After maxAtt
%  attempts a unit is RETIRED. The count is written when the claim is
%  taken, so a unit that hangs hard enough to kill the process still spends
%  an attempt. It is NOT reset by opening the queue again -- the first
%  version wiped attempts on every 'init', which handed the livelock back
%  to the next run. Only an explicit 'reset' clears them. An attempt record
%  that cannot be read or written FAILS CLOSED: unreadable means blocked,
%  unwritable means no claim.
%
%% Inputs:
%
%  action                   char                    'init' | 'open' |
%                                                   'create' | 'reset' |
%                                                   'claim' | 'beat' |
%                                                   'release' | 'status'
%  qDir                     char                    queue directory
%  varargin                 per action:
%    init(qDir, units, outFcn)   SAFE: 'open' if the queue exists, else
%                                'create'. Never destroys claims or attempts.
%    create(qDir, units, outFcn) a NEW queue; refuses if one exists
%    open(qDir, units, outFcn)   re-open: keep claims and attempts, add any
%                                new units, update outFcn (a HANDLE:
%                                func2str drops what it closes over)
%    reset(qDir, opts)           EXPLICIT wipe of claims and attempts;
%                                refuses while any claim is fresh unless
%                                opts.force
%    claim(qDir, tag, opts)      -> out.id (NaN if nothing available now),
%                                out.claim, out.token, out.reason
%                                opts .staleSec [1800] .maxAtt [3]
%    beat(qDir, claim, token)    -> true if still the owner
%    release(qDir, claim, token) -> true if released (was the owner)
%    status(qDir, opts)          -> one state per unit: done | running |
%                                stale | retired | todo; opts .staleSec
%                                .maxAtt as for claim
%
%% Outputs:
%
%  out                      per action (above)
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

qMat = fullfile(qDir, 'queue.mat');
switch lower(action)
    case 'init'
        if isfile(qMat), out = work_queue('open', qDir, varargin{:});
        else,            out = work_queue('create', qDir, varargin{:}); end

    case 'create'
        units = varargin{1};  outFcn = varargin{2};
        assert(~isfile(qMat), 'work_queue:exists', ...
               ['a queue already exists in %s: use ''open'' to resume it, or ' ...
                '''reset'' to start over deliberately'], qDir);
        if ~isfolder(qDir), mkdir(qDir); end
        S = struct('units', units(:).', 'outFcn', outFcn, 'created', char(datetime('now')));
        saveAtomic(qMat, S);
        out = S;

    case 'open'
        units = varargin{1};  outFcn = varargin{2};
        assert(isfile(qMat), 'work_queue:missing', 'no queue in %s to open', qDir);
        S = load(qMat);
        added = setdiff(units(:).', S.units);
        S.units = [S.units, added];  S.outFcn = outFcn;
        S.opened = char(datetime('now'));
        saveAtomic(qMat, S);
        out = S;  out.added = added;

    case 'reset'
        opts = struct();  if ~isempty(varargin), opts = varargin{1}; end
        force = fieldd(opts, 'force', false);  staleSec = fieldd(opts, 'staleSec', 1800);
        cl = dir(fullfile(qDir, '*.claim'));
        fresh = cl(arrayfun(@(c) c.isdir && claimAge(fullfile(qDir, c.name)) < staleSec, cl));
        assert(force || isempty(fresh), 'work_queue:live', ...
               ['reset refused: %d claim(s) are FRESH (a worker is on them). Stop the ' ...
                'workers first, or pass .force to discard their ownership.'], numel(fresh));
        delete(fullfile(qDir, '*.att'));
        for k = 1:numel(cl)
            p_ = fullfile(qDir, cl(k).name);
            if isfolder(p_), rmdir(p_, 's'); else, delete(p_); end
        end
        out = numel(cl);

    case 'claim'
        tag = varargin{1};
        opts = struct();  if numel(varargin) > 1, opts = varargin{2}; end
        staleSec = fieldd(opts, 'staleSec', 1800);
        maxAtt = fieldd(opts, 'maxAtt', 3);
        Q = load(qMat);
        outFcn = Q.outFcn;
        out = struct('id', NaN, 'claim', '', 'token', '', 'reason', 'nothing available now');
        for id = Q.units
            if isfile(outFcn(id)), continue, end              % DONE: artifact exists
            if attempts(qDir, id) >= maxAtt, continue, end    % RETIRED (or unreadable)
            cf = fullfile(qDir, sprintf('%d.claim', id));
            if isfolder(cf)
                if claimAge(cf) < staleSec, continue, end     % someone live owns it
                if ~takeover(cf), continue, end               % another reclaimer won
            end
            [okMk, ~, idMk] = mkdir(cf);                      % ATOMIC: create-or-exists
            if ~okMk || strcmp(idMk, 'MATLAB:MKDIR:DirectoryExists')
                continue                                      % lost the race, try the next
            end
            % RE-CHECK UNDER OWNERSHIP: both tests above were made before
            % the claim existed, and the previous owner may have finished
            % or spent the last attempt in between
            if isfile(outFcn(id)) || attempts(qDir, id) >= maxAtt
                rmdir(cf, 's');  continue
            end
            token = sprintf('%s:%s', tag, randomHex(8));
            try
                bumpAttempt(qDir, id);        % BEFORE the work; fails closed
                writeBeat(cf, token, tag);
            catch ME
                rmdir(cf, 's');
                rethrow(ME);
            end
            out.id = id;  out.claim = cf;  out.token = token;  out.reason = 'claimed';
            return
        end

    case 'beat'
        cf = varargin{1};  token = varargin{2};
        out = false;
        if ~isempty(cf) && isfolder(cf) && ownedBy(cf, token)
            writeBeat(cf, token, tagOf(token));  out = true;
        end

    case 'release'
        cf = varargin{1};  token = varargin{2};
        out = false;
        if ~isempty(cf) && isfolder(cf) && ownedBy(cf, token)
            rmdir(cf, 's');  out = true;
        end

    case 'status'
        opts = struct();  if ~isempty(varargin) && isstruct(varargin{1}), opts = varargin{1}; end
        staleSec = fieldd(opts, 'staleSec', 1800);
        maxAtt = fieldd(opts, 'maxAtt', 3);
        Q = load(qMat);
        outFcn = Q.outFcn;
        n = numel(Q.units);
        state = repmat({'todo'}, 1, n);  att = zeros(1, n);  owner = repmat({''}, 1, n);
        for k = 1:n
            id = Q.units(k);
            att(k) = attempts(qDir, id);
            cf = fullfile(qDir, sprintf('%d.claim', id));
            % ONE state per unit, in priority order. A live final attempt
            % is RUNNING, not retired; a unit is retired only when nobody
            % is on it and its attempts are spent.
            if isfile(outFcn(id))
                state{k} = 'done';
            elseif isfolder(cf)
                owner{k} = tagOf(readToken(cf));
                if claimAge(cf) < staleSec, state{k} = 'running';
                else,                         state{k} = 'stale'; end
            elseif att(k) >= maxAtt
                state{k} = 'retired';
            end
        end
        is = @(s) strcmp(state, s);
        out = struct('units', Q.units, 'state', {state}, 'owner', {owner}, 'attempts', att, ...
                     'done', Q.units(is('done')), 'running', Q.units(is('running')), ...
                     'stale', Q.units(is('stale')), 'retired', Q.units(is('retired')), ...
                     'todo', Q.units(is('todo')), ...
                     'nDone', nnz(is('done')), 'nRunning', nnz(is('running')), ...
                     'nStale', nnz(is('stale')), 'nRetired', nnz(is('retired')), ...
                     'nTodo', nnz(is('todo')), 'n', n);
        % what the caller usually wants to know, stated instead of inferred
        out.nOpen = out.nTodo + out.nStale;                 % dispatchable now
        out.complete = out.nDone == n;                       % every artifact exists
        out.finished = out.nDone + out.nRetired == n;        % nothing more will change

    otherwise
        error('work_queue:action', 'unknown action "%s"', action);
end
end

% ------------------------------------------------------------------------
function ok = takeover(cf)
% TAKEOVER  Remove a stale claim so that exactly ONE reclaimer proceeds:
% rename it to a unique name first (atomic), then delete the renamed copy.
% A second reclaimer's rename fails because the source is gone.
% INPUTS: cf.  OUTPUTS: ok logical.
tomb = sprintf('%s.stale.%s', cf, randomHex(6));
[ok, ~] = movefile(cf, tomb);
if ok, try rmdir(tomb, 's'); catch, end, end
end

% ------------------------------------------------------------------------
function n = attempts(qDir, id)
% ATTEMPTS  How many times unit id has been claimed. FAILS CLOSED: a record
% that exists but cannot be parsed is reported as Inf (blocked), never as
% zero, so corruption cannot re-open an exhausted unit.
% INPUTS: qDir; id.  OUTPUTS: n.
f = fullfile(qDir, sprintf('%d.att', id));
if ~isfile(f), n = 0;  return, end
try
    v = str2double(strtrim(fileread(f)));
catch
    v = NaN;
end
if isnan(v) || v < 0 || v ~= round(v), n = inf; else, n = v; end
end

% ------------------------------------------------------------------------
function bumpAttempt(qDir, id)
% BUMPATTEMPT  Record a claim BEFORE the work, so a unit that hangs hard
% enough to kill the process still spends an attempt. Written atomically
% (temp + move); an I/O failure is an ERROR, not a silent zero.
% INPUTS: qDir; id.  OUTPUTS: none.
f = fullfile(qDir, sprintf('%d.att', id));
n = attempts(qDir, id);
assert(isfinite(n), 'work_queue:attempts', 'attempt record %s is unreadable; refusing to claim', f);
tmp = sprintf('%s.%s.tmp', f, randomHex(4));
fid = fopen(tmp, 'w');
assert(fid >= 0, 'work_queue:attempts', 'cannot write the attempt record %s', f);
fprintf(fid, '%d', n + 1);  fclose(fid);
[ok, msg] = movefile(tmp, f);
assert(ok, 'work_queue:attempts', 'cannot publish the attempt record %s: %s', f, msg);
end

% ------------------------------------------------------------------------
function writeBeat(cf, token, tag)
% WRITEBEAT  (Re)write the beat file inside a claim: its mtime is the
% progress signal and its first field is the owner token. Written atomically
% so a reader never sees a truncated record.
% INPUTS: cf; token; tag.  OUTPUTS: none.
f = fullfile(cf, 'beat');
tmp = sprintf('%s.%s.tmp', f, randomHex(4));
fid = fopen(tmp, 'w');
assert(fid >= 0, 'work_queue:beat', 'cannot write the beat file in %s', cf);
fprintf(fid, '%s|%s|%s\n', token, tag, char(datetime('now')));  fclose(fid);
movefile(tmp, f);
end

% ------------------------------------------------------------------------
function t = readToken(cf)
% READTOKEN  The owner token recorded in a claim ('' if none yet).
% INPUTS: cf.  OUTPUTS: t char.
t = '';
f = fullfile(cf, 'beat');
if ~isfile(f), return, end
try
    parts = strsplit(strtrim(fileread(f)), '|');
    t = parts{1};
catch
end
end

% ------------------------------------------------------------------------
function ok = ownedBy(cf, token)
% OWNEDBY  True if the claim still carries this token.  INPUTS: cf; token.
% OUTPUTS: ok.
ok = ~isempty(token) && strcmp(readToken(cf), token);
end

% ------------------------------------------------------------------------
function tag = tagOf(token)
% TAGOF  The worker tag part of a token ('tag:hex').  INPUTS: token.
% OUTPUTS: tag.
c = find(token == ':', 1);
if isempty(c), tag = token; else, tag = token(1:c-1); end
end

% ------------------------------------------------------------------------
function a = claimAge(cf)
% CLAIMAGE  Seconds since the claim last reported progress.  Reads the beat
% FILE, because modifying a file does not update its directory's mtime. A
% claim with no beat yet (a worker between mkdir and its first write) is
% aged from the directory itself.
% INPUTS: cf.  OUTPUTS: a.
d = dir(fullfile(cf, 'beat'));
if isempty(d), d = dir(cf); end                          % '.' carries the dir mtime
if isempty(d), a = inf;  return, end
a = seconds(datetime('now') - datetime(min([d.datenum]), 'ConvertFrom', 'datenum'));
end

% ------------------------------------------------------------------------
function saveAtomic(f, S)
% SAVEATOMIC  save() to a temp name, then move: a worker loading queue.mat
% never sees a half-written file.  INPUTS: f; S struct.  OUTPUTS: none.
tmp = sprintf('%s.%s.tmp', f, randomHex(4));
save(tmp, '-struct', 'S');
movefile(tmp, f);
end

% ------------------------------------------------------------------------
function h = randomHex(n)
% RANDOMHEX  n random hex characters (independent of the global RNG state).
% INPUTS: n.  OUTPUTS: h char.
seed = mod(round(posixtime(datetime('now'))*1e6) + matlabProcessID, 2^31);
s = RandStream('mt19937ar', 'Seed', seed);
h = lower(dec2hex(randi(s, [0 15], 1, n), 1)).';
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
