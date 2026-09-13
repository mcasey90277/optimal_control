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
staleSec = fieldd(opts, 'staleSec', 1800);
S.queue = work_queue('status', qDir);
S.workers = campaign_heartbeat('read', hbDir, tags, struct('staleSec', staleSec));
S.alarm = {};

never = {S.workers(strcmp({S.workers.state}, 'never')).tag};
stall = {S.workers(strcmp({S.workers.state}, 'stalled')).tag};
fail  = {S.workers(strcmp({S.workers.state}, 'failed')).tag};
live  = nnz(strcmp({S.workers.state}, 'running'));
if ~isempty(never), S.alarm{end+1} = sprintf('NEVER STARTED: %s', strjoin(never, ' ')); end
if ~isempty(stall), S.alarm{end+1} = sprintf('STALLED: %s', strjoin(stall, ' ')); end
if ~isempty(fail),  S.alarm{end+1} = sprintf('FAILED: %s', strjoin(fail, ' ')); end
if isfield(S.queue, 'nRetired') && S.queue.nRetired > 0
    S.alarm{end+1} = sprintf('%d unit(s) RETIRED after repeated failure: %s', ...
                             S.queue.nRetired, mat2str(S.queue.retired));
end
if S.queue.nStale > 0
    S.alarm{end+1} = sprintf('%d unit(s) held by a dead claim (reclaimable): %s', ...
                             S.queue.nStale, mat2str(S.queue.stale));
end
if live == 0 && S.queue.nTodo > 0
    S.alarm{end+1} = sprintf('NO LIVE WORKER but %d unit(s) still to do', S.queue.nTodo);
end

if ~fieldd(opts, 'quiet', false)
    fprintf('campaign: %d/%d done, %d running, %d todo, %d stale, %d retired\n', ...
            S.queue.nDone, S.queue.n, S.queue.nRunning, S.queue.nTodo, S.queue.nStale, ...
            fieldd(S.queue, 'nRetired', 0));
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
