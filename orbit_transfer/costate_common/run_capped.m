function [ok, varargout] = run_capped(pool, fcn, nout, capSec, varargin)
%% Purpose:
%
%   Runs fcn(args) on a parfeval worker under a HARD wall-clock cap; on
%   timeout or worker error the future is CANCELLED (the worker is killed
%   and restarted), so no single stuck computation can stall the caller.
%
%   This is the only fence that bounds a CRAWLING integration: in-process
%   wall checks (fsolve OutputFcn, loop-top toc guards) fire only between
%   iterations, and a CR3BP segment propagation whose iterate parks a
%   junction near a primary can grind for hours inside ONE function
%   evaluation (measured twice on the 0.5 N extension campaign,
%   2026-08-31, >4 h each). Promoted from extend_thrust_ladder on its
%   second consumer (probe_deep_rungs), per the consolidation queue.
%
%% Inputs:
%
%  pool                     parallel.Pool           From gcp (workers
%                                                   inherit the client
%                                                   path at pool creation)
%
%  fcn                      function_handle         The function to run
%
%  nout                     Integer                 Number of outputs to
%                                                   request from fcn
%
%  capSec                   double                  Hard wall-clock cap [s]
%
%  varargin                 any                     Arguments passed to fcn
%
%% Outputs:
%
%  ok                       logical                 true if fcn finished
%                                                   in time without error
%
%  varargout                any                     fcn's nout outputs
%                                                   ([] each when ~ok)
%
%  cause (opt-in)           char                    ONLY when the caller
%                                                   requests nout + 2
%                                                   outputs: '' when ok,
%                                                   'capped' on timeout,
%                                                   else the worker error
%                                                   'identifier: message'
%                                                   (fut.Error, otherwise
%                                                   discarded). Callers
%                                                   asking for nout + 1
%                                                   outputs are unchanged
%
%% Revision History:
%  M. Casey                                                   (c) 09/01/2026
%  M. Casey  opt-in cause output (worker error text)           (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

fut  = parfeval(pool, fcn, nout, varargin{:});
done = wait(fut, 'finished', capSec);
if ~done
    % A CANCEL IS A REQUEST. A worker stuck in native code may not honour
    % it, and then every later call on this one-worker pool would time out
    % behind it while the client kept beating -- activity without work
    % (Astra pass 3, 4.3). So the cancellation is VERIFIED: if the future
    % has not finished within 30 s, the pool is deleted (which kills its
    % worker process) and the caller gets an infrastructure error, not a
    % numerical failure.
    cancel(fut);
    stopped = wait(fut, 'finished', 30);
    if ~stopped || strcmp(fut.State, 'running')
        try delete(pool); catch, end
        error('run_capped:poolUnusable', ...
              ['a capped call did not stop within 30 s of cancellation (cap %g s); the pool was ' ...
               'deleted so that a fresh attempt starts with a fresh worker'], capSec);
    end
    ok = false;  varargout = repmat({[]}, 1, nout);  cause = 'capped';
elseif ~isempty(fut.Error)
    ok = false;  varargout = repmat({[]}, 1, nout);
    cause = error_text(fut.Error);
else
    ok = true;   [varargout{1:nout}] = fetchOutputs(fut);  cause = '';
end
if nargout > nout + 1, varargout{nout + 1} = cause; end
end

% ------------------------------------------------------------------------
function t = error_text(E)
% ERROR_TEXT  'identifier: message' of a future's error, innermost remote
% cause first when parfeval wrapped it.  INPUTS: E MException.  OUTPUTS: t char.
try
    if ~isempty(E.remotecause), E = E.remotecause{1}; end
catch
end
try
    t = sprintf('%s: %s', E.identifier, E.message);
catch
    t = 'worker error (no MException available)';
end
end
