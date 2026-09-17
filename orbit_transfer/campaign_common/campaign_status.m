function S = campaign_status(qDir, hbDir, tags, opts)
%% Purpose:
%
%   What a campaign is actually doing, read from ARTIFACTS ONLY: the queue's
%   outputs and the workers' heartbeats. No process matching, because that
%   is what reported four live workers dead (their tag was in the
%   environment, not the command line).
%
%   Prints, and returns, the two things worth acting on: units not being
%   worked on, and workers that are not working.
%
%% Inputs:
%
%  qDir, hbDir              char                    campaign directories
%  tags                     cellstr                 worker names expected
%  opts                     struct (optional)
%   .staleSec [1800] .quiet [false]
%
%% Outputs:
%
%  S                        struct                  .queue (work_queue
%                                                   status) .workers
%                                                   (heartbeat states)
%                                                   .alarm (cellstr: what
%                                                   needs attention)
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 4, opts = struct(); end
staleSec = fieldd(opts, 'staleSec', 1800);  maxAtt = fieldd(opts, 'maxAtt', 3);
% ONE policy for the queue and the heartbeats: status used to hand the
% queue a hardcoded 1800 s while forwarding staleSec only to the workers
S.queue = work_queue('status', qDir, struct('staleSec', staleSec, 'maxAtt', maxAtt));
S.workers = campaign_heartbeat('read', hbDir, tags, struct('staleSec', staleSec));
S.alarm = {};

never = {S.workers(strcmp({S.workers.state}, 'never')).tag};
stall = {S.workers(strcmp({S.workers.state}, 'stalled')).tag};
fail  = {S.workers(strcmp({S.workers.state}, 'failed')).tag};
unk   = {S.workers(strcmp({S.workers.state}, 'unknown')).tag};
live  = nnz(strcmp({S.workers.state}, 'running'));
if ~isempty(unk),   S.alarm{end+1} = sprintf('UNREADABLE HEARTBEAT: %s', strjoin(unk, ' ')); end
if ~isempty(never), S.alarm{end+1} = sprintf('NEVER STARTED: %s', strjoin(never, ' ')); end
if ~isempty(stall), S.alarm{end+1} = sprintf('STALLED: %s', strjoin(stall, ' ')); end
if ~isempty(fail),  S.alarm{end+1} = sprintf('FAILED: %s', strjoin(fail, ' ')); end
if S.queue.nRetired > 0
    S.alarm{end+1} = sprintf('%d unit(s) RETIRED after repeated failure: %s', ...
                             S.queue.nRetired, mat2str(S.queue.retired));
end
if S.queue.nStalled > 0
    % held by a LIVE process (the lock says so) that has not beaten for
    % staleSec: a hang, or a solve longer than the policy expects. The
    % launcher's inactivity watchdog is what acts on it; this only says so.
    S.alarm{end+1} = sprintf('%d running unit(s) with no beat for > %d s (owner alive, possibly hung): %s', ...
                             S.queue.nStalled, staleSec, mat2str(S.queue.units(S.queue.stalled)));
end
if S.queue.nAbandoned > 0
    S.alarm{end+1} = sprintf('%d unit(s) whose owner died holding them (claimable): %s', ...
                             S.queue.nAbandoned, mat2str(S.queue.abandoned));
end
% JOIN: every running unit's owner must be a worker we can see beating
runIdx = find(strcmp(S.queue.state, 'running'));
wtags = {S.workers.tag};
for k = runIdx
    ow = S.queue.owner{k};
    w = find(strcmp(wtags, ow), 1);
    if isempty(w)
        S.alarm{end+1} = sprintf('unit %d is held by "%s", which is not a worker of this launch', S.queue.units(k), ow);
    elseif ~any(strcmp(S.workers(w).state, {'running', 'stalled'}))
        S.alarm{end+1} = sprintf('unit %d is held by %s whose heartbeat says %s', S.queue.units(k), ow, S.workers(w).state);
    end
end
if live == 0 && S.queue.nOpen > 0
    S.alarm{end+1} = sprintf('NO LIVE WORKER but %d unit(s) still claimable', S.queue.nOpen);
end

if ~fieldd(opts, 'quiet', false)
    fprintf('campaign: %d/%d done, %d running, %d todo, %d abandoned, %d retired%s\n', ...
            S.queue.nDone, S.queue.n, S.queue.nRunning, S.queue.nTodo, S.queue.nAbandoned, ...
            S.queue.nRetired, tern(S.queue.finished, ' -- FINISHED', ''));
    for k = 1:numel(S.workers)
        w = S.workers(k);
        fprintf('  %-6s %-8s %s%s\n', w.tag, w.state, w.msg, ...
                tern(isnan(w.age), '', sprintf(' (%.0f s ago)', w.age)));
    end
    if isempty(S.alarm), fprintf('  nothing needs attention\n');
    else, fprintf('  ** %s\n', S.alarm{:}); end
end
end

% ------------------------------------------------------------------------
function v = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: v.
if c, v = a; else, v = b; end
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
