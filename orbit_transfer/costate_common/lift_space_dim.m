function [dimS, tol, gap] = lift_space_dim(sv, nullResid, rankTol)
%% Purpose:
%
%   Numerical dimension of the lift space S from the singular values of the
%   stacked lift-constraint matrix C (mintime_hypothesis_gates):
%
%       dim S = #{ k : sv(k) < tol },
%       tol   = min( max(rankTol*sv(1), 10*nullResid), 1e-3*sv(1) ).
%
%   Why the residual enters: the accepted normal lift lam0 is KNOWN to be in
%   S, and it satisfies C*lam0 = 0 only to nullResid = |C lam0|/|lam0| (the
%   accuracy of the flight + adjoint integration). A null space cannot be
%   resolved finer than its known member's residual, so a fixed 1e-8 rank
%   tolerance under-counts whenever nullResid > 1e-8 -- the 2026-09-07
%   catalog pass flagged 512 entries "dim S = 0" for exactly that reason
%   (sv(7)/sv(6) ~ 1e-7..2.5e-6, i.e. a clean one-dimensional null space).
%   The 1e-3*sv(1) cap keeps a noisy lift from inflating the count: with
%   the cap, a poor lift reports dim S = 0 (unresolved), never > 1.
%
%% Inputs:
%
%  sv                       [n x 1]                 singular values of C,
%                                                   descending (svd order)
%
%  nullResid                double                  |C lam0| / |lam0|
%
%  rankTol                  double                  relative rank tolerance
%                                                   (mintime_hypothesis_gates
%                                                   default 1e-8)
%
%% Outputs:
%
%  dimS                     double                  numerical dim S
%
%  tol                      double                  the threshold used
%
%  gap                      double                  sv(end)/sv(end-1): the
%                                                   spectral gap that makes
%                                                   dim S = 1 decisive
%
%% Revision History:
%  M. Casey                                                   (c) 09/07/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

sv  = sv(:);
tol = min(max(rankTol*sv(1), 10*nullResid), 1e-3*sv(1));
dimS = nnz(sv < tol);
gap  = sv(end) / max(sv(end-1), realmin);
end
