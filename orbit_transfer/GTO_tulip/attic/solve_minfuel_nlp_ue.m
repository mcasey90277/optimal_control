function [Z, out] = solve_minfuel_nlp_ue(Z0, sigma, tf, rv0, m0, rvf, Tmax, c, muStar, maxIter)
% SOLVE_MINFUEL_NLP_UE  Direct min-fuel solve, cone-ELIMINATED (unit-dir).
%
% Same fixed-tf max-final-mass problem as SOLVE_MINFUEL_NLP, but with the
% direction+throttle control [alpha; s] and the unit-sphere constraint
% ||alpha|| = 1 (NLP_CONSTRAINTS_MINFUEL_UE) instead of the (w, s) cone.
% Decoupling direction from throttle removes the w->0 degeneracy at coasts,
% so the interior-point method can actually drive the throttle onto its
% bounds (sharp bang-bang) rather than stalling on the cone.
%
% INPUTS:  as SOLVE_MINFUEL_NLP, plus
%   maxIter - (optional) fmincon MaxIterations [default 5000]
%
% OUTPUTS:
%   Z, out - out has .X, .U [alpha;s], .mf, .exitflag, .maxDefect (ALL
%            constraints), .maxDynDefect, .maxUnit, .eqnonlin
%
% REFERENCES:
%   [1] Betts, SIAM 2010, Ch. 4.

if nargin < 10 || isempty(maxIter), maxIter = 5000; end
sigma  = sigma(:);
N      = numel(sigma) - 1;
nNodes = N + 1;
nZ     = 11*nNodes;

assert(numel(Z0) == nZ, 'solve_minfuel_nlp_ue:z0', ...
       'Z0 has %d elements; expected 11*(N+1) = %d', numel(Z0), nZ);
assert(sigma(1) == 0 && sigma(end) == 1 && all(diff(sigma) > 0), ...
       'solve_minfuel_nlp_ue:sigma', 'sigma must increase strictly 0 -> 1');

% --- bounds: alpha in [-1.1,1.1] (unit dir), s in [0,1] --------------------
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

idxMf = 7*nNodes;
gradJ = sparse(nZ, 1);  gradJ(idxMf) = -1;
objFun = @(Z) deal(-Z(idxMf), gradJ);
conFun = @(Z) nlp_constraints_minfuel_ue(Z, sigma, tf, Tmax, c, muStar);

opts = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true, ...
    'InitBarrierParam', 1e-3, ...
    'HessianApproximation', 'lbfgs', ...
    'MaxIterations', maxIter, ...
    'MaxFunctionEvaluations', 5e5, ...
    'ConstraintTolerance', 1e-10, ...
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance', 1e-14, ...
    'Display', 'iter');

[Z, ~, exitflag, output, lambdaS] = fmincon(objFun, Z0, [], [], [], [], ...
                                            lb, ub, conFun, opts);

X = reshape(Z(1:7*nNodes), 7, nNodes);
U = reshape(Z(7*nNodes + (1:4*nNodes)), 4, nNodes);
[~, ceq] = nlp_constraints_minfuel_ue(Z, sigma, tf, Tmax, c, muStar);

out = struct('X', X, 'U', U, 'mf', X(7, end), 'exitflag', exitflag, ...
             'maxDefect', max(abs(ceq)), 'maxDynDefect', max(abs(ceq(1:7*N))), ...
             'maxUnit', max(abs(ceq(7*N+1:end))), 'fmincon', output, ...
             'eqnonlin', lambdaS.eqnonlin);
end
