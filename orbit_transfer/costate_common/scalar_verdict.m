function [ok, v] = scalar_verdict(x)
%% Purpose:
%
%   Read an external result as a gate verdict, or refuse to read it.
%
%   A gate that interprets whatever a solver hands back can be passed by
%   malformed data rather than by satisfied physics. In MATLAB every one of
%   these slips through a plain `if` (Astra chain review, 2026-09-10):
%
%     []      `if v ~= 1` is FALSE on empty  -> an empty verdict acts like a pass
%     NaN     every comparison is false      -> a NaN residual trips no bound
%     Inf     `Inf > 0` is true              -> a degenerate gate reads satisfied
%     [1 1]   `if` on a vector means ALL     -> not the test anyone intended
%
%   So: a verdict is a REAL FINITE SCALAR or it is not a verdict. Callers
%   test `ok` first and fail the candidate with a named reason when it is
%   false -- never fall through to interpreting `v`.
%
%% Inputs:
%
%  x                        any                     the external result
%
%% Outputs:
%
%  ok                       logical                 x is a real finite
%                                                   numeric/logical scalar
%  v                        double                  double(x) when ok, NaN
%                                                   otherwise (so an
%                                                   unchecked use compares
%                                                   false rather than true)
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = (isnumeric(x) || islogical(x)) && isscalar(x) && isreal(x) && isfinite(x);
if ok, v = double(x); else, v = NaN; end
end
