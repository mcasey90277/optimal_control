function out = work_queue(action, qDir, varargin)
%% Purpose:
%
%   A DISK WORK QUEUE for campaign units, so that N worker processes on one
%   host pull the next unclaimed unit instead of being handed static
%   ranges. The unit is whatever the caller says it is (a rib column, a
%   catalog entry, a thrust rung).
%
%   Why this exists (2026-09-13, doc/CAMPAIGN_DISCIPLINE.md): static ranges
%   idled workers behind their own slow column, a watchdog inside a unit
%   destroyed nine hours of work twice, and a unit that always failed
%   livelocked the first version of this queue.
%
%  OWNERSHIP IS A PROCESS-HELD LOCK (unit_lock: a kernel file lock the
%  owning MATLAB process holds until it releases it or dies). It is never
%  transferred on heartbeat age. The first version reclaimed a claim whose
%  beat was 30 minutes old, which let a slow-but-alive owner and its
%  reclaimer both compute and both publish one unit (Astra chain review
%  pass 2, A/B). Now a live owner cannot be stolen from; a hung owner is
%  killed by the launcher's inactivity watchdog, the kernel frees the lock,
%  and only then can another worker take the unit. Heartbeat age is an
%  ALARM, not an authority.
%
%  DONE IS A PUBLISHED ARTIFACT. A unit is done when its output file
%  exists. Outputs reach their names only through 'publish', which the
%  owner calls WHILE HOLDING THE LOCK, from an attempt-specific temporary
%  file, after an optional validation, in one rename (publish_atomic). No
%  separate completion flag can disagree with the file on disk, and a
%  process that lost its lock cannot publish.
%
%  ATTEMPTS ARE COUNTED UNDER THE LOCK, BEFORE THE WORK, AND PERSIST. After
%  maxAtt attempts a unit is RETIRED: terminal, reported, never handed out.
%  A unit whose owner died on its last attempt is retired the moment the
%  lock is free -- not left half-open. Only an explicit 'reset' clears
%  attempts, and it cannot touch a unit whose lock is held. An attempt
%  record that cannot be read is BLOCKED (never zero); one that cannot be
%  written means no claim.
%
%  THE UNIT SET IS FIXED AT CREATION. 'open' validates that the caller's
%  units and outputs are the ones on disk and refuses otherwise: a resumed
%  campaign must mean the same thing as the one it resumes.
%
%% Inputs:
%
%  action                   char                    see below
%  qDir                     char                    queue directory
%  varargin                 per action:
%    create(qDir, units, outputs)  NEW queue. units [1 x n] ids; outputs
%                                  cellstr of artifact paths (one per unit)
%                                  or a handle outputs(id) evaluated NOW.
%                                  Refuses if a queue exists.
%    open(qDir, units, outputs)    resume: READ-ONLY, validates units and
%                                  outputs against the queue on disk
%    init(qDir, units, outputs)    open if a queue exists, else create
%    reset(qDir, opts)             clear attempts (and dead-owner records)
%                                  of units whose lock is NOT held;
%                                  opts.units [all] restricts it
%    claim(qDir, tag, opts)        -> c: .id (NaN if nothing claimable now)
%                                  .token .output .tmpOut (attempt-specific
%                                  temporary name in the output's folder)
%                                  .lock (KEEP IT: releasing it is
%                                  releasing the unit) .reason
%                                  opts .maxAtt [3]
%    beat(qDir, c)                 -> true if still the owner
%    publish(qDir, c, src, validateFcn) -> [ok, msg]: validate src (fcn
%                                  returns [ok, msg]; optional), then move
%                                  it onto c.output. Requires the lock.
%    release(qDir, c)              give the unit back (done or not)
%    status(qDir, opts)            -> one state per unit: done | running |
%                                  abandoned | retired | todo, plus .held
%                                  (lock probe) .beatAge .owner .stalled
%                                  (age > opts.staleSec [1800], ALARM only)
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
        [units, outputs] = unitsAndOutputs(varargin{:});
        assert(~isfile(qMat), 'work_queue:exists', ...
               'a queue already exists in %s: ''open'' resumes it; a different campaign needs a different directory', qDir);
        if ~isfolder(qDir), mkdir(qDir); end
        S = struct('units', units, 'outputs', {outputs}, 'created', char(datetime('now')));
        saveAtomic(qMat, S);
        out = S;

    case 'open'
        [units, outputs] = unitsAndOutputs(varargin{:});
        assert(isfile(qMat), 'work_queue:missing', 'no queue in %s to open', qDir);
        S = load(qMat);
        assert(isequal(sort(S.units), sort(units)), 'work_queue:units', ...
               ['the queue in %s holds units %s, not %s: a campaign''s unit set is fixed when it is ' ...
                'created. Use a different directory for a different grid.'], qDir, mat2str(S.units), mat2str(units));
        [~, ia] = ismember(units, S.units);
        assert(isequal(S.outputs(ia), outputs), 'work_queue:outputs', ...
               'the queue in %s was created with different output paths; refusing to reinterpret it', qDir);
        out = S;

    case 'reset'
        opts = struct();  if ~isempty(varargin), opts = varargin{1}; end
        S = load(qMat);
        which = fieldd(opts, 'units', S.units);
        out = [];
        for id = which(:).'
            if unit_lock('probe', lockFile(qDir, id)).held
                warning('work_queue:held', 'unit %d is held by a live process; not reset', id);
                continue
            end
            f = attFile(qDir, id);   if isfile(f), delete(f); end
            f = ownerFile(qDir, id); if isfile(f), delete(f); end
            out(end+1) = id; %#ok<AGROW>
        end

    case 'claim'
        tag = varargin{1};
        opts = struct();  if numel(varargin) > 1, opts = varargin{2}; end
        maxAtt = fieldd(opts, 'maxAtt', 3);
        Q = load(qMat);
        out = struct('id', NaN, 'token', '', 'output', '', 'tmpOut', '', 'lock', [], ...
                     'owner', '', 'reason', 'nothing claimable now');
        for k = 1:numel(Q.units)
            id = Q.units(k);  art = Q.outputs{k};
            if isfile(art), continue, end                     % DONE
            if attempts(qDir, id) >= maxAtt, continue, end    % RETIRED (or unreadable)
            L = unit_lock('try', lockFile(qDir, id));
            if ~L.held, continue, end                         % a live process owns it
            % RE-CHECK UNDER THE LOCK: the world may have moved since the
            % two tests above
            if isfile(art) || attempts(qDir, id) >= maxAtt
                unit_lock('release', '', L);  continue
            end
            token = sprintf('%s:%s', tag, char(java.util.UUID.randomUUID()));
            try
                bumpAttempt(qDir, id);                        % under the lock; fails closed
                writeOwner(qDir, id, token, tag);
            catch ME
                unit_lock('release', '', L);
                rethrow(ME);
            end
            [ad, an, ae] = fileparts(art);
            out.id = id;  out.token = token;  out.output = art;  out.lock = L;
            out.tmpOut = fullfile(ad, sprintf('%s%s.%s.part', an, ae, token(numel(tag)+2:end)));
            out.owner = ownerFile(qDir, id);  out.reason = 'claimed';
            return
        end

    case 'beat'
        c = varargin{1};
        out = holds(c);
        if out, writeOwner(qDir, c.id, c.token, tagOf(c.token)); end

    case 'publish'
        c = varargin{1};  src = varargin{2};
        validateFcn = [];  if numel(varargin) > 2, validateFcn = varargin{3}; end
        out = false;  msg = '';
        if ~holds(c), msg = 'this process no longer holds the unit''s lock';
        elseif ~isfile(src), msg = sprintf('nothing to publish: %s does not exist', src);
        else
            okV = true;
            if ~isempty(validateFcn)
                try
                    [okV, msg] = validateFcn(src);
                catch ME
                    okV = false;  msg = ['validation threw: ' ME.message];
                end
            end
            if okV
                publish_atomic(src, c.output);
                out = true;
            end
        end
        out = struct('ok', out, 'msg', msg);

    case 'release'
        c = varargin{1};
        if holds(c)
            f = ownerFile(qDir, c.id);  if isfile(f), delete(f); end
            if isfield(c, 'tmpOut') && isfile(c.tmpOut), delete(c.tmpOut); end
        end
        if isstruct(c) && isfield(c, 'lock'), unit_lock('release', '', c.lock); end
        out = [];

    case 'status'
        opts = struct();  if ~isempty(varargin) && isstruct(varargin{1}), opts = varargin{1}; end
        staleSec = fieldd(opts, 'staleSec', 1800);
        maxAtt = fieldd(opts, 'maxAtt', 3);
        Q = load(qMat);
        n = numel(Q.units);
        state = repmat({'todo'}, 1, n);  att = zeros(1, n);  owner = repmat({''}, 1, n);
        held = false(1, n);  age = nan(1, n);  stalled = false(1, n);
        for k = 1:n
            id = Q.units(k);
            att(k) = attempts(qDir, id);
            held(k) = unit_lock('probe', lockFile(qDir, id)).held;
            [owner{k}, age(k)] = readOwner(qDir, id);
            % ONE state per unit. The lock decides who is live; the owner
            % record without a lock is a process that died holding it.
            if isfile(Q.outputs{k})
                state{k} = 'done';
            elseif held(k)
                state{k} = 'running';  stalled(k) = age(k) > staleSec;
            elseif att(k) >= maxAtt
                state{k} = 'retired';
            elseif ~isempty(owner{k})
                state{k} = 'abandoned';
            end
        end
        is = @(s) strcmp(state, s);
        out = struct('units', Q.units, 'outputs', {Q.outputs}, 'state', {state}, 'owner', {owner}, ...
                     'attempts', att, 'held', held, 'beatAge', age, 'stalled', stalled, ...
                     'done', Q.units(is('done')), 'running', Q.units(is('running')), ...
                     'abandoned', Q.units(is('abandoned')), 'retired', Q.units(is('retired')), ...
                     'todo', Q.units(is('todo')), ...
                     'nDone', nnz(is('done')), 'nRunning', nnz(is('running')), ...
                     'nAbandoned', nnz(is('abandoned')), 'nRetired', nnz(is('retired')), ...
                     'nTodo', nnz(is('todo')), 'nHeld', nnz(held), 'nStalled', nnz(stalled), 'n', n);
        out.nOpen = out.nTodo + out.nAbandoned;              % claimable now
        out.complete = out.nDone == n;                        % every artifact exists
        out.finished = (out.nDone + out.nRetired == n) && out.nHeld == 0;   % nothing will change

    otherwise
        error('work_queue:action', 'unknown action "%s"', action);
end
end

% ------------------------------------------------------------------------
function [units, outputs] = unitsAndOutputs(units, outputs)
% UNITSANDOUTPUTS  Normalise the (units, outputs) pair: a handle is
% evaluated now so the queue stores explicit paths, not a closure.
% INPUTS: units; outputs handle or cellstr.  OUTPUTS: units row; outputs cell row.
units = units(:).';
if isa(outputs, 'function_handle')
    outputs = arrayfun(outputs, units, 'UniformOutput', false);
end
outputs = outputs(:).';
assert(iscellstr(outputs) && numel(outputs) == numel(units), 'work_queue:outputs', ...
       'outputs must be one path per unit');
end

% ------------------------------------------------------------------------
function ok = holds(c)
% HOLDS  True if c is a live claim this process still holds.  INPUTS: c.
% OUTPUTS: ok.
ok = isstruct(c) && isfield(c, 'lock') && isstruct(c.lock) && c.lock.held;
end

% ------------------------------------------------------------------------
function f = lockFile(qDir, id),  f = fullfile(qDir, sprintf('%d.lock', id));  end
function f = attFile(qDir, id),   f = fullfile(qDir, sprintf('%d.att', id));   end
function f = ownerFile(qDir, id), f = fullfile(qDir, sprintf('%d.owner', id)); end

% ------------------------------------------------------------------------
function n = attempts(qDir, id)
% ATTEMPTS  How many times unit id has been claimed. FAILS CLOSED: a record
% that exists but cannot be parsed is Inf (blocked), never zero.
% INPUTS: qDir; id.  OUTPUTS: n.
f = attFile(qDir, id);
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
% BUMPATTEMPT  Record a claim BEFORE the work, under the unit's lock, so a
% unit that hangs hard enough to kill the process still spends an attempt.
% INPUTS: qDir; id.  OUTPUTS: none (throws if it cannot be recorded).
f = attFile(qDir, id);
n = attempts(qDir, id);
assert(isfinite(n), 'work_queue:attempts', 'attempt record %s is unreadable; refusing to claim', f);
writeText(f, sprintf('%d', n + 1));
end

% ------------------------------------------------------------------------
function writeOwner(qDir, id, token, tag)
% WRITEOWNER  (Re)write the owner record; its mtime is the beat.
% INPUTS: qDir; id; token; tag.  OUTPUTS: none.
writeText(ownerFile(qDir, id), sprintf('%s|%s|%d|%s', token, tag, matlabProcessID, char(datetime('now'))));
end

% ------------------------------------------------------------------------
function [tag, age] = readOwner(qDir, id)
% READOWNER  Owner tag and seconds since the last beat ('' / NaN if none).
% INPUTS: qDir; id.  OUTPUTS: tag; age.
tag = '';  age = NaN;
f = ownerFile(qDir, id);
d = dir(f);
if isempty(d), return, end
age = seconds(datetime('now') - datetime(d.datenum, 'ConvertFrom', 'datenum'));
try
    parts = strsplit(strtrim(fileread(f)), '|');
    if numel(parts) >= 2, tag = parts{2}; else, tag = tagOf(parts{1}); end
catch
    tag = '?';
end
end

% ------------------------------------------------------------------------
function tag = tagOf(token)
% TAGOF  The worker tag part of a token ('tag:uuid').  INPUTS: token.
% OUTPUTS: tag.
c = find(token == ':', 1);
if isempty(c), tag = token; else, tag = token(1:c-1); end
end

% ------------------------------------------------------------------------
function writeText(f, txt)
% WRITETEXT  Write a small record atomically: exclusive temp, then publish.
% INPUTS: f; txt.  OUTPUTS: none (throws on any I/O failure).
tmp = sprintf('%s.%s.tmp', f, char(java.util.UUID.randomUUID()));
fid = fopen(tmp, 'w');
assert(fid >= 0, 'work_queue:io', 'cannot write %s', tmp);
nw = fprintf(fid, '%s', txt);
st = fclose(fid);
assert(nw == numel(txt) && st == 0, 'work_queue:io', 'short write to %s', tmp);
publish_atomic(tmp, f);
end

% ------------------------------------------------------------------------
function saveAtomic(f, S)
% SAVEATOMIC  save() to an exclusive temp name, then publish: a worker
% loading queue.mat never sees a half-written file.  INPUTS: f; S.
% OUTPUTS: none.
tmp = sprintf('%s.%s.tmp', f, char(java.util.UUID.randomUUID()));
save(tmp, '-struct', 'S');
publish_atomic(tmp, f);
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
