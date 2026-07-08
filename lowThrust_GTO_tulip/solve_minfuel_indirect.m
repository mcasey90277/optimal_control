function [lamSol, resNorm, epsHist] = solve_minfuel_indirect(rv0, m0, rvf, tf, lamGuess, Tmax, c, muStar, epsSchedule)
% SOLVE_MINFUEL_INDIRECT  Min-fuel TPBVP by shooting + smoothing continuation.
%
% Solves the fixed-tf minimum-fuel problem by walking the Bertrand-Epenoy
% smoothing parameter down a schedule, re-solving the 7-costate shooting
% problem at each step from the previous converged solution. Large eps
% makes the throttle gentle and the convergence basin wide; small eps
% sharpens toward the true bang-bang solution.
%
% INPUTS:
%   rv0         - initial position/velocity (ND) [1x6]
%   m0          - initial mass fraction [scalar]
%   rvf         - target position/velocity (ND) [1x6]
%   tf          - fixed transfer time (ND) [scalar]
%   lamGuess    - initial costate guess [7x1] (e.g. the min-time solution)
%   Tmax        - max thrust acceleration at m = 1 (ND) [scalar]
%   c           - exhaust velocity (ND) [scalar]
%   muStar      - Earth-Moon mass ratio [scalar]
%   epsSchedule - (optional) decreasing smoothing values
%                 [default 1 -> 0.3 -> 0.1 -> 0.03 -> 0.01 -> 3e-3 -> 1e-3]
%
% OUTPUTS:
%   lamSol  - converged initial costates at the final (smallest) eps [7x1]
%   resNorm - residual 2-norm at lamSol [scalar]
%   epsHist - struct array: per-eps converged costates, residual, flag
%
% REFERENCES:
%   [1] Bertrand & Epenoy, OCAM 23(4), 2002.

if nargin < 9 || isempty(epsSchedule)
    epsSchedule = [1, 0.3, 0.1, 0.03, 0.01, 3e-3, 1e-3];
end

% Levenberg-Marquardt handles this landscape far better than fsolve's
% trust-region-dogleg (measured: 4x deeper residual decrease per budget).
opts = optimoptions('lsqnonlin', ...
    'Display', 'off', ...
    'Algorithm', 'levenberg-marquardt', ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance', 1e-24, ...
    'StepTolerance', 1e-14, ...
    'MaxIterations', 200, ...
    'MaxFunctionEvaluations', 1000);

lam     = lamGuess(:);
epsHist = struct('epsSmooth', {}, 'lambda0', {}, 'resNorm', {}, 'flag', {});

for epsSmooth = epsSchedule
    resFun = @(lam0) shoot_residual_minfuel(lam0, tf, rv0, m0, rvf, ...
                                            Tmax, c, muStar, epsSmooth);
    [lam, res2, ~, flag] = lsqnonlin(resFun, lam, [], [], opts);
    resNorm = sqrt(res2);
    epsHist(end+1) = struct('epsSmooth', epsSmooth, 'lambda0', lam, ...
                            'resNorm', resNorm, 'flag', flag); %#ok<AGROW>
    fprintf('  eps = %-8.3g  ||R|| = %-10.3g  flag %d\n', ...
            epsSmooth, resNorm, flag);
    if resNorm > 1e-2
        warning('solve_minfuel_indirect:stall', ...
                'continuation stalled at eps = %g', epsSmooth);
        break;
    end
end
lamSol = lam;
end
