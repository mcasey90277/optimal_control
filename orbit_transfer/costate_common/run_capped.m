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
%% Revision History:
%  M. Casey                                                   (c) 09/01/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

fut  = parfeval(pool, fcn, nout, varargin{:});
done = wait(fut, 'finished', capSec);
if ~done || ~isempty(fut.Error)
    cancel(fut);
    ok = false;  varargout = repmat({[]}, 1, nout);
else
    ok = true;   [varargout{1:nout}] = fetchOutputs(fut);
end
end
