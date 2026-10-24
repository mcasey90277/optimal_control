function pb = antipodal_pair_bound(x0, xT, T, r)
%% Purpose:
%
%   Rigorous lower bound on the hard-clearance cost: relax the problem to
%   its antipodal SWAP pairs and drop every other constraint. Agents i and
%   j form a pair only when the swap is mutual (x0(:,j) == xT(:,i) AND
%   x0(:,i) == xT(:,j)), and each agent sits in at most one pair, so the
%   pairs' variables are disjoint and their costs add. Each pair is then the
%   two-agent head-on problem, whose cost is at least ell^2/(4T) with
%       ell = 2 sqrt(D^2 - R^2) + R (pi - 2 acos(R/D)),
%   D the pair distance and R = 2r (tests/test_direct_oracle.m checks the
%   closed form against a solve). Agents in no pair contribute zero, so the
%   bound is 0 on a layout without swaps. Used by drone_table1 and
%   run_drones; tests/test_table1_fixture.m keeps its own copy.
%
%% Inputs:
%
%  x0                       [d x M]                 Start positions
%
%  xT                       [d x M]                 Target positions
%
%  T                        double                  Horizon (s)
%
%  r                        double                  Agent radius (R = 2r)
%
%% Outputs:
%
%  pb                       double                  Sum over the pairs of
%                                                   ell^2/(4T); errors with
%                                                   antipodal_pair_bound:
%                                                   geometry if a pair has
%                                                   D <= R
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the paper's M = 4 circle (two head-on pairs at D = 1, r = 0.02):
     here = fileparts(mfilename('fullpath'));
     addpath(fullfile(fileparts(here), 'drone_dynamics'));
     [a, b, g] = drone_layouts('paper', 4);
     pb = antipodal_pair_bound(a, b, 10, g.r);
     fprintf('antipodal_pair_bound demo (paper M = 4, T = 10): %.7f (expected 0.2003202)\n', pb);
     return
end

tol  = 1e-12;
R    = 2*r;
M    = size(x0, 2);
used = false(1, M);
pb   = 0;
for ia = 1:M
    if used(ia), continue, end
    for ib = ia+1:M
        if used(ib), continue, end
        swapAB = max(abs(x0(:, ib) - xT(:, ia))) < tol;
        swapBA = max(abs(x0(:, ia) - xT(:, ib))) < tol;
        if swapAB && swapBA
            D = sqrt(sum((x0(:, ia) - x0(:, ib)).^2));
            if D <= R
                error('antipodal_pair_bound:geometry', ...
                      'agents %d and %d: pair distance %.3g <= R = %.3g', ia, ib, D, R);
            end
            ell = 2*sqrt(D^2 - R^2) + R*(pi - 2*acos(R/D));
            pb  = pb + ell^2/(4*T);
            used([ia ib]) = true;
            break
        end
    end
end
end
