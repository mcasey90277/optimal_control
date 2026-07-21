function [Z, out] = solve_minfuel_nlp(Z0, sigma, tf, rv0, m0, rvf, Tmax, c, muStar)
% SOLVE_MINFUEL_NLP  Direct min-fuel solve via fmincon (interior-point).
%
% Maximizes final mass (objective -m_{N+1}, linear) subject to trapezoidal
% defects, the throttle cone w'*w = s^2 with s in [0,1], endpoint pinning,
% and fixed transfer time. Unlike the min-time solve, throttle bounds WILL
% go active at the solution (bang-bang structure) -- that is fine for an
% interior-point method as long as the WARM START is strictly interior
% (clip the guess throttle away from 0/1; see BUILD_GUESS_MINFUEL).
%
% INPUTS:
%   Z0     - initial decision vector [11*(N+1) x 1]
%   sigma  - normalized node times [(N+1)x1]
%   tf     - fixed transfer time (ND) [scalar]
%   rv0    - initial position/velocity (ND) [1x6]
%   m0     - initial mass fraction (< 1 for a mid-transfer leg) [scalar]
%   rvf    - target position/velocity (ND) [1x6]
%   Tmax   - max thrust acceleration at m = 1 (ND) [scalar]
%   c      - exhaust velocity (ND) [scalar]
%   muStar - Earth-Moon mass ratio [scalar]
%
% OUTPUTS:
%   Z      - converged decision vector
%   out    - struct: .X [7x(N+1)], .U [4x(N+1)], .mf (final mass fraction),
%            .exitflag, .maxDefect, .fmincon, .eqnonlin (KKT multipliers of
%            [defects; cone] -- the discrete costates, see covector mapping)
%
% REFERENCES:
%   [1] Betts, SIAM 2010, Ch. 4.

sigma  = sigma(:);
N      = numel(sigma) - 1;
nNodes = N + 1;
nZ     = 11*nNodes;

assert(numel(Z0) == nZ, 'solve_minfuel_nlp:z0', ...
       'Z0 has %d elements; expected 11*(N+1) = %d', numel(Z0), nZ);
assert(sigma(1) == 0 && sigma(end) == 1 && all(diff(sigma) > 0), ...
       'solve_minfuel_nlp:sigma', 'sigma must increase strictly 0 -> 1');

% --- bounds ---------------------------------------------------------------
lb = -inf(nZ, 1);  ub = inf(nZ, 1);
for k = 1:nNodes
    xIdx = (k-1)*7 + (1:7);
    lb(xIdx) = [-3; -3; -3; -12; -12; -12; 0.3];
    ub(xIdx) = [ 3;  3;  3;  12;  12;  12; 1.0];
    uIdx = 7*nNodes + (k-1)*4 + (1:4);
    lb(uIdx) = [-1.1; -1.1; -1.1; 0];
    ub(uIdx) = [ 1.1;  1.1;  1.1; 1];
end
lb(1:7)   = [rv0(:); m0];  ub(1:7)   = [rv0(:); m0];
xfIdx     = (nNodes-1)*7 + (1:6);
lb(xfIdx) = rvf(:);       ub(xfIdx) = rvf(:);

% --- objective: minimize -m(tf) (maximize final mass), linear --------------
idxMf = 7*nNodes;                       % m at the last node
gradJ = sparse(nZ, 1);  gradJ(idxMf) = -1;
objFun = @(Z) deal(-Z(idxMf), gradJ);

conFun = @(Z) nlp_constraints_minfuel(Z, sigma, tf, Tmax, c, muStar);

opts = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true, ...
    'InitBarrierParam', 1e-4, ...
    'HessianApproximation', 'lbfgs', ...
    'MaxIterations', 3000, ...
    'MaxFunctionEvaluations', 1e5, ...
    'ConstraintTolerance', 1e-10, ...
    'OptimalityTolerance', 1e-7, ...
    'StepTolerance', 1e-12, ...
    'Display', 'iter');

[Z, ~, exitflag, output, lambdaS] = fmincon(objFun, Z0, [], [], [], [], ...
                                             lb, ub, conFun, opts);

X = reshape(Z(1:7*nNodes), 7, nNodes);
U = reshape(Z(7*nNodes + (1:4*nNodes)), 4, nNodes);
[~, ceq] = nlp_constraints_minfuel(Z, sigma, tf, Tmax, c, muStar);

out = struct('X', X, 'U', U, 'mf', X(7, end), 'exitflag', exitflag, ...
             'maxDefect', max(abs(ceq(1:7*N))), 'fmincon', output, ...
             'eqnonlin', lambdaS.eqnonlin);
end
