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
budgetSec = d('budgetSec', inf);  staleSec = d('staleSec', 1800);
maxUnits = d('maxUnits', inf);  logFile = d('logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

t0 = tic;  nDone = 0;  nFail = 0;  units = [];  failed = [];  stopped = 'queue empty';
campaign_heartbeat('beat', hbDir, tag, 'starting');
lg('worker %s: start', tag);

while true
    if toc(t0) > budgetSec, stopped = 'budget';  break, end
    if nDone + nFail >= maxUnits, stopped = 'maxUnits';  break, end
    c = work_queue('claim', qDir, tag, struct('staleSec', staleSec));
    if isnan(c.id), stopped = 'queue empty';  break, end

    campaign_heartbeat('beat', hbDir, tag, sprintf('unit %d', c.id));
    lg('worker %s: unit %d claimed', tag, c.id);
    tU = tic;
    beat = @() beatBoth(qDir, c.claim, hbDir, tag, c.id);
    try
        unitFcn(c.id, beat);
        nDone = nDone + 1;  units(end+1) = c.id; %#ok<AGROW>
        lg('worker %s: unit %d DONE in %.0f s', tag, c.id, toc(tU));
    catch ME
        nFail = nFail + 1;  failed(end+1) = c.id; %#ok<AGROW>
        lg('worker %s: unit %d FAILED after %.0f s -- %s', tag, c.id, toc(tU), ME.message);
    end
    % RELEASED either way: the ARTIFACT, not the claim, says done. A unit
    % that threw must return to the queue rather than look finished.
    work_queue('release', qDir, c.claim);
end

out = struct('nDone', nDone, 'nFailed', nFail, 'units', units, 'failed', failed, ...
             'wallSec', toc(t0), 'stopped', stopped);
msg = sprintf('%d done, %d failed, stopped: %s', nDone, nFail, stopped);
campaign_heartbeat('done', hbDir, tag, msg);
lg('worker %s: %s', tag, msg);
end

% ------------------------------------------------------------------------
function beatBoth(qDir, claim, hbDir, tag, id)
% BEATBOTH  Report progress to the claim (so it is not reclaimed) AND to the
% heartbeat (so the monitor sees it).  INPUTS: qDir; claim; hbDir; tag; id.
% OUTPUTS: none.
work_queue('beat', qDir, claim);
campaign_heartbeat('beat', hbDir, tag, sprintf('unit %d', id));
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
