function [zSol, resNorm, exitflag] = solve_tfmin_indirect(rv0, rvf, zGuess, Tmax, c, muStar)
% SOLVE_TFMIN_INDIRECT  Solve the min-time TPBVP by single shooting.
%
% Wraps fsolve (trust-region-dogleg, complex-step Jacobian supplied by
% SHOOT_RESIDUAL_TF) to find the 7 initial costates and the transfer
% time that satisfy the minimum-time necessary conditions. Finite
% differencing is NOT viable here -- see the sensitivity note in
% SHOOT_RESIDUAL_TF.
%
% INPUTS:
%   rv0    - initial position/velocity (ND, rotating frame) [1x6]
%   rvf    - target position/velocity (ND) [1x6]
%   zGuess - initial guess [8x1]: [lambda_r0(3); lambda_v0(3); lambda_m0; tf]
%   Tmax   - max thrust acceleration at m = 1 (ND) [scalar]
%   c      - exhaust velocity (ND) [scalar]
%   muStar - Earth-Moon mass ratio [scalar]
%
% OUTPUTS:
%   zSol     - converged decision vector [8x1]
%   resNorm  - 2-norm of the residual at zSol [scalar]
%   exitflag - fsolve exit flag [scalar]
%
% REFERENCES:
%   [1] pumpkyn.cr3bp.tfMin (analytic-Jacobian counterpart).

opts = optimoptions('fsolve', ...
    'Display', 'iter', ...
    'Algorithm', 'trust-region-dogleg', ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance', 1e-24, ...
    'StepTolerance', 1e-12, ...
    'MaxIterations', 100, ...
    'MaxFunctionEvaluations', 400);

resFun = @(z) shoot_residual_tf(z, rv0, rvf, Tmax, c, muStar);
[zSol, Rf, exitflag] = fsolve(resFun, zGuess(:), opts);
resNorm = sqrt(sum(Rf.^2));
end
