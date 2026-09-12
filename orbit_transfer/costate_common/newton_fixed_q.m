function [p, converged, normR, nCalls] = newton_fixed_q(resFactory, q, p, Dx, tol, nMax)
%% Purpose:
%
%   Scaled Newton at a FIXED continuation parameter, used to land exactly on
%   a requested level.
%
%   Extracted from arclength_ms on its second consumer (crossings_from_arc,
%   2026-09-12), per the library rule. The two must not diverge: a crossing
%   re-scanned out of a saved arc has to be corrected by the SAME step the
%   walk would have applied at that level, or it is not the same crossing.
%
%  ASSUMPTIONS / NOTES:
%
% • resFactory(q) returns a handle h(p) -> [R, J]: the residual at that
%   fixed parameter and its Jacobian in p.
% • Dx scales the unknowns; the step is taken in scaled coordinates and
%   mapped back, so the tolerance is on the residual, not on p.
% • A non-finite residual or step stops the iteration and reports
%   converged = false: a Newton that has left the domain has not converged,
%   and saying so is the caller's cue to discard the candidate.
%
%% Inputs:
%
%  resFactory               fhandle                 q -> h(p) -> [R, J]
%  q                        double                  the fixed parameter
%  p                        [n x 1]                 starting guess
%  Dx                       [n x 1]                 unknown scaling
%  tol                      double                  |R|_inf to accept
%  nMax                     double                  iteration cap
%
%% Outputs:
%
%  p                        [n x 1]                 corrected point
%  converged                logical                 |R|_inf < tol, finite
%  normR                    double                  final |R|_inf
%  nCalls                   double                  residual evaluations
%
%% Revision History:
%  M. Casey                                                   (c) 09/12/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

converged = false;  nCalls = 0;  normR = inf;
h = resFactory(q);
for it = 1:nMax
    [R, J] = h(p);  nCalls = nCalls + 1;
    normR = norm(R, inf);
    if ~all(isfinite(R)), return, end
    if normR < tol, converged = true; return, end
    dx = -((J .* Dx(:)') \ R);
    if ~all(isfinite(dx)), return, end
    p = p + dx .* Dx;
end
R = h(p);  nCalls = nCalls + 1;  normR = norm(R, inf);
converged = all(isfinite(R)) && normR < tol;
end
