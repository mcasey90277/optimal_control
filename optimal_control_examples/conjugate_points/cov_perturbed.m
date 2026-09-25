function P = cov_perturbed(prob, ext, delta, tEnd)
%% Purpose:
%
%   Fly the NEIGHBOURING extremal that leaves the same left point with slope
%   p + delta, and find where it crosses the original extremal again.
%
%   For small delta the neighbour is y0 + delta h + O(delta^2), so its
%   crossings approach the zeros of the Jacobi field h -- the conjugate
%   points. For a finite delta on a nonlinear problem the crossing is only
%   NEAR the conjugate point (off by O(delta)); shrinking delta shows it
%   converge. On a linear problem (quadratic F) the neighbour is exactly
%   y0 + delta h and the crossing is exactly the conjugate point.
%
%% Inputs:
%
%  prob                     struct                  cov_problem output
%  ext                      struct                  one cov_extremals entry
%  delta                    scalar                  slope change, ~= 0
%  tEnd                     scalar                  how far to fly, > a
%
%% Outputs:
%
%  P                        struct                  .S (cov_shoot of the
%                                                   neighbour) .delta .tCross
%                                                   (crossings of y0 in
%                                                   (a, tMax], sorted) .tMax
%                                                   (where both flights
%                                                   exist)
%
%% Revision History:
%  M. Casey                                                   (c) 09/24/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

assert(delta ~= 0, 'cov_perturbed:delta', 'delta must be nonzero');
a = prob.a;
S0 = cov_shoot(prob, ext.p, tEnd);
S1 = cov_shoot(prob, ext.p + delta, tEnd);
P = struct('S', S1, 'delta', delta, 'tCross', [], 'tMax', a);
if isempty(S0.sol) || isempty(S1.sol), return, end
tMax = min(S0.tStop, S1.tStop);
P.tMax = tMax;
if tMax <= a, return, end

% sample the difference, skipping the shared start point, then refine
t = linspace(a, tMax, 4001);
t = t(2:end);
d = diffAt(t);
tc = [];
for k = 1:numel(t)-1
    if d(k) == 0
        tc(end+1) = t(k);                                            %#ok<AGROW>
    elseif sign(d(k)) ~= sign(d(k+1))
        tc(end+1) = fzero(@diffAt, t([k k+1]), optimset('TolX', 1e-13)); %#ok<AGROW>
    end
end
P.tCross = tc;

    function d = diffAt(tt)
        % neighbour minus original, at the times tt
        z1 = deval(S1.sol, tt);  z0 = deval(S0.sol, tt);
        d = z1(1,:) - z0(1,:);
    end
end
