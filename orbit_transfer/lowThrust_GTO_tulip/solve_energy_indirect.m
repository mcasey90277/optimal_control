function [lamSol, resNorm, flag] = solve_energy_indirect(rv0, m0, rvf, tf, lamGuess, Tmax, c, muStar)
% SOLVE_ENERGY_INDIRECT  Min-energy TPBVP by single shooting (no homotopy).
%
% Solves the fixed-tf minimum-energy problem with LEVENBERG-MARQUARDT on
% the 7-condition shooting residual (complex-step Jacobian). Unlike
% min-fuel, the min-energy control is already continuous, so NO smoothing
% continuation is required -- a single solve from a decent seed suffices.
% The remaining difficulty on the full multi-revolution spiral is the
% integrator sensitivity (~40 perigee passes), the same wall the min-time
% shooter faced; a covector-mapped seed from the direct solution is the
% robust source (see GTO_tulip driver).
%
% INPUTS:
%   rv0      - initial position/velocity (ND) [1x6]
%   m0       - initial mass fraction [scalar]
%   rvf      - target position/velocity (ND) [1x6]
%   tf       - fixed transfer time (ND) [scalar]
%   lamGuess - initial costate guess [7x1]
%   Tmax     - max thrust acceleration at m = 1 (ND) [scalar]
%   c        - exhaust velocity (ND) [scalar]
%   muStar   - Earth-Moon mass ratio [scalar]
%
% OUTPUTS:
%   lamSol  - converged initial costates [7x1]
%   resNorm - residual 2-norm at lamSol [scalar]
%   flag    - lsqnonlin exit flag [scalar]
%
% REFERENCES:
%   [1] Caillau, Gergaud, Noailles, JOTA 2003.

opts = optimoptions('lsqnonlin', ...
    'Display', 'off', ...
    'Algorithm', 'levenberg-marquardt', ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance', 1e-24, ...
    'StepTolerance', 1e-14, ...
    'MaxIterations', 300, ...
    'MaxFunctionEvaluations', 2000);

resFun = @(lam0) shoot_residual_energy(lam0, tf, rv0, m0, rvf, Tmax, c, muStar);
[lamSol, res2, ~, flag] = lsqnonlin(resFun, lamGuess(:), [], [], opts);
resNorm = sqrt(res2);
fprintf('  energy shoot: ||R|| = %.3g, flag %d\n', resNorm, flag);
end
