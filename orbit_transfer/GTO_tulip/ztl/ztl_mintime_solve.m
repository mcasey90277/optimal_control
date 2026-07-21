function [z, resNorm, out] = ztl_mintime_solve(rv0, rvf, zGuess, Tmax, c, muStar, maxEvals)
% ZTL_MINTIME_SOLVE  Min-time CR3BP TPBVP: pumpkyn's analytic-STM residual,
% real solver budgets.
%
% Rebuilds exactly the residual/Jacobian of pumpkyn.cr3bp.tfMin's private
% shootingResidual -- rendezvous (6), lambda_m(tf)=0, H(tf)=0, with dR/dx
% assembled from the propagated 14x14 STM (pumpkyn.cr3bp.tfMinProp, ode45
% RelTol 1e-10 with throttle-switch event restarts) -- and solves it with
% lsqnonlin trust-region-reflective at REAL budgets. pumpkyn's own wrapper
% hardcodes MaxFunctionEvaluations = MaxIterations = 100, which stalls every
% rung of a thrust march (ZTL P0d finding).
%
% CAVEAT: the STM is integrated ACROSS throttle-switch events without a
% saltation correction, so J is exact only while the min-time arc is all-burn
% (S < 0 throughout); with interior coast arcs it is approximate. out.nSwitch
% reports the converged arc's switch count so the caller can judge.
%
% INPUTS:
%   rv0      - initial position/velocity, ND rotating frame [1x6]
%   rvf      - target position/velocity, ND [1x6]
%   zGuess   - [8x1]: [lambda_r(3); lambda_v(3); lambda_m; tf]
%   Tmax     - max thrust acceleration at m = 1 (ND) [scalar]
%   c        - exhaust velocity (ND) [scalar]
%   muStar   - Earth-Moon mass ratio [scalar]
%   maxEvals - (optional) MaxFunctionEvaluations [default 1500]
%
% OUTPUTS:
%   z       - solution [8x1] (costates + tf_min)
%   resNorm - ||R|| at z [scalar]
%   out     - struct: .flag .nSwitch .R
%
% REFERENCES:
%   [1] pumpkyn.cr3bp.tfMin (Koblick 2025) -- source of residual/J assembly.
%   [2] Zhang et al., JGCD 38(8), 2015 (STM-based indirect shooting).

if nargin < 7 || isempty(maxEvals), maxEvals = 1500; end
m0 = 1;

% tf lower bound: impulsive minimum dV at infinite thrust (as in pumpkyn)
dV    = pumpkyn.cr3bp.minDeltaV(rv0(:)', rvf(:)', muStar);
tfLow = dV/Tmax;
LB = [-1e3*ones(6,1); -100; tfLow];
UB = [+1e3*ones(6,1); +100; 100*tfLow];

opts = optimoptions('lsqnonlin', ...
    'Display', 'off', ...
    'Algorithm', 'trust-region-reflective', ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance', 1e-24, ...
    'StepTolerance', 1e-14, ...
    'MaxIterations', 1e3, ...
    'MaxFunctionEvaluations', maxEvals);

[z, res2, ~, flag] = lsqnonlin(@resfun, zGuess(:), LB, UB, [], [], [], [], [], opts);
resNorm = sqrt(res2);

% converged-arc switch count (S sign changes along the trajectory)
[~, Y] = pumpkyn.cr3bp.tfMinProp(z(8), [rv0(:); m0; z(1:7)], Tmax, c, muStar);
Ssgn = sign(-sqrt(sum(Y(:,11:13).^2, 2))*c./Y(:,7) - Y(:,14));
out  = struct('flag', flag, 'nSwitch', sum(abs(diff(Ssgn > 0))), 'R', resfun(z));

% ---------------------------------------------------------------------------
    function [R, J] = resfun(x)
        lam = x(1:7);  tf = x(8);
        if nargout > 1
            PHI0 = eye(14);
            y0 = [rv0(:); m0; lam; PHI0(:)];
        else
            y0 = [rv0(:); m0; lam];
        end
        [~, Yaug] = pumpkyn.cr3bp.tfMinProp(tf, y0, Tmax, c, muStar);
        yf = Yaug(end, 1:14).';
        [Ff, Hf, dHdy] = pumpkyn.cr3bp.tfMinEoM(tf, yf, Tmax, c, muStar);
        R = [yf(1:6) - rvf(:); yf(14); Hf];
        if nargout > 1
            PHIf = reshape(Yaug(end, 15:210), 14, 14);
            J = [PHIf(1:6, 8:14),      Ff(1:6);
                 PHIf(14,  8:14),      Ff(14);
                 dHdy*PHIf(:, 8:14),   dHdy*Ff];
        end
    end
end
