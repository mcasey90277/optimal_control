function out = campaign_worker(qDir, hbDir, tag, unitFcn, opts)
%% Purpose:
%
%   One campaign worker: claim a unit, run it, PUBLISH its validated result
%   under the unit's lock, release, repeat, until the campaign is finished
%   or the budget is spent. Launched as its own MATLAB process by
%   run_campaign_workers.sh, N at a time, against one work_queue.
%
%   THE UNIT CONTRACT. unitFcn(id, beat, tmpOut) must write its result to
%   tmpOut (an attempt-specific temporary name in the output's folder) and
%   call beat() at its natural inner cadence. It must NOT write the final
%   output itself: the worker validates tmpOut (opts.validateFcn) and moves
%   it onto the unit's output in one rename, while it still holds the lock.
%   A function that returns without a valid tmpOut is a FAILED attempt, not
%   a done unit -- "returned normally" used to count as success (Astra pass
%   2, P1-17).
%
%   READY. The worker writes its first heartbeat only after it has opened
%   the queue, so the launcher's verification means "attached to the
%   queue", not "a file appeared".
%
%   CHECKPOINTS. A unit function may keep a resumable checkpoint at
%   <output>.ckpt (walk_checkpoint); it belongs to the UNIT, so the next
%   attempt after a kill resumes from it. The worker deletes it after the
%   unit is published.
%
%   THE TAIL. When nothing is claimable but units are still held by other
%   workers, this one WAITS (beating 'idle') until the campaign is finished
%   -- so if the last long column's owner dies, someone is still there to
%   pick it up. opts.idleSec bounds the wait.
%
%% Inputs:
%
%  qDir, hbDir              char                    queue and heartbeat dirs
%  tag                      char                    worker identity from the
%                                                   launcher (required)
%  unitFcn                  function_handle         unitFcn(id, beat, tmpOut)
%  opts                     struct (optional)
%   .validateFcn []  [ok, msg] = validateFcn(tmpOut, id), run before publishing
%   .budgetSec [inf] stop cleanly between units after this
%   .idleSec [inf]   longest wait for held units before leaving
%   .maxAtt [3]      attempts after which a unit is retired (queue policy)
%   .staleSec [1800] heartbeat age reported as STALLED (alarm only)
%   .maxUnits [inf] .logFile ''
%
%% Outputs:
%
%  out                      struct                  .nDone .nFailed .units
%                                                   .failed .wallSec .stopped
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 5, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
logFile = d('logFile', '');
lg = @(varargin) safeLog(logFile, varargin{:});
assert(ischar(tag) && ~isempty(tag), 'campaign_worker:tag', 'a worker needs a real tag from its launcher');

% THE WHOLE LIFECYCLE IS GUARDED: a failure anywhere leaves a 'fail'
% heartbeat, never a stale 'running' one
try
    out = lifecycle(qDir, hbDir, tag, unitFcn, d, lg);
catch ME
    safeCall(@() campaign_heartbeat('fail', hbDir, tag, ME.message));
    lg('worker %s: FATAL %s', tag, ME.message);
    rethrow(ME);
end
end

% ------------------------------------------------------------------------
function out = lifecycle(qDir, hbDir, tag, unitFcn, d, lg)
% LIFECYCLE  The worker loop proper.  INPUTS: as campaign_worker; d option
% getter; lg logger.  OUTPUTS: out.
budgetSec = d('budgetSec', inf);  idleSec = d('idleSec', inf);
maxAtt = d('maxAtt', 3);  staleSec = d('staleSec', 1800);
maxUnits = d('maxUnits', inf);  validateFcn = d('validateFcn', []);
policy = struct('staleSec', staleSec, 'maxAtt', maxAtt);

t0 = tic;  nDone = 0;  nFail = 0;  units = [];  failed = [];  stopped = 'finished';
st = work_queue('status', qDir, policy);            % READY means: the queue opened
campaign_heartbeat('ready', hbDir, tag, sprintf('%d units, %d done', st.n, st.nDone));   % persistent marker
campaign_heartbeat('beat', hbDir, tag, sprintf('ready %d units', st.n));
lg('worker %s: ready, queue %s (%d units, %d done)', tag, qDir, st.n, st.nDone);
tIdle = [];
while true
    if toc(t0) > budgetSec, stopped = 'budget';  break, end
    if nDone + nFail >= maxUnits, stopped = 'maxUnits';  break, end
    c = work_queue('claim', qDir, tag, struct('maxAtt', maxAtt));
    if isnan(c.id)
        st = work_queue('status', qDir, policy);
        if st.finished, stopped = 'finished';  break, end
        % held by others: WAIT, so the tail is never unattended
        if isempty(tIdle), tIdle = tic; end
        if toc(tIdle) > idleSec, stopped = 'idle timeout';  break, end
        campaign_heartbeat('beat', hbDir, tag, sprintf('idle: %d unit(s) held elsewhere', st.nHeld));
        pause(15);
        continue
    end
    tIdle = [];

    % THE UNIT IS RELEASED WHATEVER HAPPENS BELOW: an error in a heartbeat,
    % a log line or the cleanup itself must not leave the lock held by a
    % worker that has stopped working (release is idempotent, so the
    % ordinary release further down makes this a no-op)
    guard = onCleanup(@() work_queue('release', qDir, c));
    campaign_heartbeat('beat', hbDir, tag, sprintf('unit %d', c.id));
    lg('worker %s: unit %d claimed (attempt record written)', tag, c.id);
    tU = tic;
    beat = @() beatBoth(qDir, c, hbDir, tag);
    % ACCOUNTING HAPPENS ONCE, at the commit boundary
    try
        unitFcn(c.id, beat, c.tmpOut);
        if isempty(validateFcn), vf = [];
        else, vf = @(f) validateFcn(f, c.id); end        % the validator knows WHICH unit
        P = work_queue('publish', qDir, c, c.tmpOut, vf);
        okUnit = P.ok;  err = P.msg;
    catch ME
        okUnit = false;  err = sprintf('%s (%s)', ME.message, ME.identifier);
    end
    if okUnit
        nDone = nDone + 1;  units(end+1) = c.id; %#ok<AGROW>
        lg('worker %s: unit %d DONE in %.0f s -> %s', tag, c.id, toc(tU), c.output);
        % a unit may keep a resumable checkpoint at <output>.ckpt while it
        % runs (walk_checkpoint); once published it is no longer needed
        if isfile([c.output '.ckpt']), delete([c.output '.ckpt']); end
    else
        nFail = nFail + 1;  failed(end+1) = c.id; %#ok<AGROW>
        lg('worker %s: unit %d FAILED after %.0f s -- %s', tag, c.id, toc(tU), err);
        if isfile(c.tmpOut)           % keep the evidence, out of the publisher's way
            [okMv, msgMv] = movefile(c.tmpOut, [c.tmpOut '.failed']);
            if ~okMv, lg('worker %s: could not preserve %s: %s', tag, c.tmpOut, msgMv); end
        end
    end
    work_queue('release', qDir, c);   % the artifact, not the claim, says done
    clear guard
end

out = struct('nDone', nDone, 'nFailed', nFail, 'units', units, 'failed', failed, ...
             'wallSec', toc(t0), 'stopped', stopped);
msg = sprintf('%d done, %d failed, stopped: %s', nDone, nFail, stopped);
% 'done' means THIS WORKER exited cleanly; the queue's artifacts say
% whether the campaign is complete
campaign_heartbeat('done', hbDir, tag, msg);
lg('worker %s: %s', tag, msg);
end

% ------------------------------------------------------------------------
function beatBoth(qDir, c, hbDir, tag)
% BEATBOTH  Refresh the owner record AND the heartbeat. Best-effort: a beat
% that cannot be written must not abort the solver it reports on.
% INPUTS: qDir; c; hbDir; tag.  OUTPUTS: none.
safeCall(@() work_queue('beat', qDir, c));
safeCall(@() campaign_heartbeat('beat', hbDir, tag, sprintf('unit %d', c.id)));
end

% ------------------------------------------------------------------------
function safeCall(f)
% SAFECALL  Run f, swallowing any error.  INPUTS: f handle.  OUTPUTS: none.
try f(); catch, end
end

% ------------------------------------------------------------------------
function safeLog(logFile, varargin)
% SAFELOG  Append one line to the log (or stdout) without ever throwing.
% INPUTS: logFile; fmt, args.  OUTPUTS: none.
try
    s = sprintf(varargin{:});
    if isempty(logFile), fprintf('%s\n', s);
    else, fid = fopen(logFile, 'a'); fprintf(fid, '%s %s\n', char(datetime('now', 'Format', 'HH:mm:ss')), s); fclose(fid); end
catch
end
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
