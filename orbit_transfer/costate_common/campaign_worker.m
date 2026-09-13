function out = campaign_worker(qDir, hbDir, tag, unitFcn, opts)
%% Purpose:
%
%   ONE campaign worker: pull the next unclaimed unit, do it, repeat until
%   the queue is empty or the budget runs out. Several run at once, as
%   separate processes, needing no coordination beyond the queue.
%
%   This replaces hand-assigned unit ranges, the most expensive mistake of
%   the 24x24 campaign: units vary from 90 minutes to over five hours, so a
%   worker with a slow unit queued behind another sat idle while its list
%   waited, and a worker killed mid-list took the whole remainder with it.
%   A worker here owns exactly one unit at a time.
%
%  THE CONTRACT WITH unitFcn:
%
%   unitFcn(id, beat) does the work for unit `id` and WRITES ITS ARTIFACT --
%   the file whose existence the queue reads as done. It should call beat()
%   as it progresses (cheap; it touches the claim) so a long unit is not
%   mistaken for a dead worker. If it throws, the unit is RELEASED, not
%   marked done: the next worker, or the next run, picks it up.
%
%  THE BUDGET STOPS BETWEEN UNITS, never inside one. A watchdog that fires
%  mid-unit destroys that unit's work, which is what happened to columns 13
%  and 16 after nine hours each. An external watchdog must therefore be
%  longer than one unit.
%
%% Inputs:
%
%  qDir                     char                    work_queue directory
%  hbDir                    char                    heartbeat directory
%  tag                      char                    this worker's name
%  unitFcn                  fhandle                 unitFcn(id, beat)
%  opts                     struct (optional)
%   .budgetSec [inf] stop cleanly between units after this
%   .staleSec [1800] when another worker's claim counts as abandoned
%   .maxAtt [3]      attempts after which a unit is retired (forwarded to
%                    the queue, so worker and monitor apply ONE policy)
%   .idleSec [0]     when nothing is claimable but units are still held by
%                    others, wait up to this long for one to be released
%                    or go stale instead of exiting -- so the LAST live
%                    worker is not gone when a peer dies at the tail
%   .maxUnits [inf] .logFile ''
%  The tag must be a real launcher-supplied identity: 'w0' is refused.
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
budgetSec = d('budgetSec', inf);  staleSec = d('staleSec', 1800);
maxAtt = d('maxAtt', 3);  idleSec = d('idleSec', 0);
maxUnits = d('maxUnits', inf);  logFile = d('logFile', '');
lg = @(varargin) safeLog(logFile, varargin{:});
assert(~isempty(tag) && ~strcmp(tag, 'w0'), 'campaign_worker:tag', ...
       'a worker needs a real tag from its launcher (got "%s")', tag);

% THE WHOLE LIFECYCLE IS GUARDED. A failure in claiming, logging or release
% used to escape the per-unit try and leave the heartbeat saying 'running'
% (or nothing at all) for a process that had died.
try
    out = lifecycle(qDir, hbDir, tag, unitFcn, budgetSec, staleSec, maxAtt, idleSec, maxUnits, lg);
catch ME
    safeCall(@() campaign_heartbeat('fail', hbDir, tag, ME.message));
    lg('worker %s: FATAL %s', tag, ME.message);
    rethrow(ME);
end
end

% ------------------------------------------------------------------------
function out = lifecycle(qDir, hbDir, tag, unitFcn, budgetSec, staleSec, maxAtt, idleSec, maxUnits, lg)
% LIFECYCLE  The worker loop proper.  INPUTS: as campaign_worker.
% OUTPUTS: out.
t0 = tic;  nDone = 0;  nFail = 0;  units = [];  failed = [];  stopped = 'queue empty';
campaign_heartbeat('beat', hbDir, tag, 'starting');
lg('worker %s: start', tag);
tIdle = [];
while true
    if toc(t0) > budgetSec, stopped = 'budget';  break, end
    if nDone + nFail >= maxUnits, stopped = 'maxUnits';  break, end
    c = work_queue('claim', qDir, tag, struct('staleSec', staleSec, 'maxAtt', maxAtt));
    if isnan(c.id)
        % NOTHING AVAILABLE NOW is not QUEUE EMPTY: the other units may be
        % held by live workers that could still die. With idleSec > 0 the
        % worker waits for that instead of leaving the tail unattended.
        st = work_queue('status', qDir, struct('staleSec', staleSec, 'maxAtt', maxAtt));
        if st.finished || idleSec <= 0, stopped = 'queue empty';  break, end
        if isempty(tIdle), tIdle = tic; end
        if toc(tIdle) > idleSec, stopped = 'idle timeout';  break, end
        campaign_heartbeat('beat', hbDir, tag, 'idle: waiting for held units');
        pause(min(60, staleSec/4));
        continue
    end
    tIdle = [];

    campaign_heartbeat('beat', hbDir, tag, sprintf('unit %d', c.id));
    lg('worker %s: unit %d claimed', tag, c.id);
    tU = tic;
    beat = @() beatBoth(qDir, c.claim, c.token, hbDir, tag, c.id);
    % ACCOUNTING HAPPENS ONCE, at the boundary of the work; the log lines
    % are best-effort and cannot turn a saved unit into a failed one.
    try
        unitFcn(c.id, beat);
        okUnit = true;  err = '';
    catch ME
        okUnit = false;  err = sprintf('%s (%s)', ME.message, ME.identifier);
    end
    if okUnit
        nDone = nDone + 1;  units(end+1) = c.id; %#ok<AGROW>
        lg('worker %s: unit %d DONE in %.0f s', tag, c.id, toc(tU));
    else
        nFail = nFail + 1;  failed(end+1) = c.id; %#ok<AGROW>
        lg('worker %s: unit %d FAILED after %.0f s -- %s', tag, c.id, toc(tU), err);
    end
    % RELEASED either way: the ARTIFACT, not the claim, says done. A unit
    % that threw must return to the queue rather than look finished. The
    % release is conditional on the token: if this unit was reclaimed
    % while we were blocked, the new owner's claim is left alone.
    if ~work_queue('release', qDir, c.claim, c.token)
        lg('worker %s: unit %d was RECLAIMED by another worker while we held it', tag, c.id);
    end
end

out = struct('nDone', nDone, 'nFailed', nFail, 'units', units, 'failed', failed, ...
             'wallSec', toc(t0), 'stopped', stopped);
msg = sprintf('%d done, %d failed, stopped: %s', nDone, nFail, stopped);
% 'done' means THIS WORKER exited cleanly. It says nothing about the
% campaign; the queue's artifacts do.
campaign_heartbeat('done', hbDir, tag, msg);
lg('worker %s: %s', tag, msg);
end

% ------------------------------------------------------------------------
function beatBoth(qDir, claim, token, hbDir, tag, id)
% BEATBOTH  Report progress to the claim (so it is not reclaimed) AND to the
% heartbeat (so the monitor sees it). Best-effort: a beat that cannot be
% written must not abort the solver it is reporting on.
% INPUTS: qDir; claim; token; hbDir; tag; id.  OUTPUTS: none.
safeCall(@() work_queue('beat', qDir, claim, token));
safeCall(@() campaign_heartbeat('beat', hbDir, tag, sprintf('unit %d', id)));
end

% ------------------------------------------------------------------------
function safeCall(f)
% SAFECALL  Run f, swallowing any error.  INPUTS: f handle.  OUTPUTS: none.
try f(); catch, end
end

% ------------------------------------------------------------------------
function safeLog(logFile, varargin)
% SAFELOG  logmsg that cannot throw.  INPUTS: logFile; fmt, args.
% OUTPUTS: none.
try logmsg(logFile, sprintf(varargin{:})); catch, end
end

% ------------------------------------------------------------------------
function logmsg(f, s)
% LOGMSG  Append one line to a log file, or to stdout.  INPUTS: f; s.
% OUTPUTS: none.
line = sprintf('%s  %s', char(datetime('now', 'Format', 'HH:mm:ss')), s);
if isempty(f), fprintf('%s\n', line);
else, fid = fopen(f, 'a');  fprintf(fid, '%s\n', line);  fclose(fid); end
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
