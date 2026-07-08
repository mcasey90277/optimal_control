function [Z, out] = solve_minfuel_nlp_hs(Z0, sigma, tf, rv0, m0, rvf, Tmax, c, muStar, maxIter)
% SOLVE_MINFUEL_NLP_HS  Direct min-fuel via separated Hermite-Simpson.
%
% 4th-order collocation (NLP_CONSTRAINTS_MINFUEL_HS) of the cone-eliminated
% min-fuel problem. The higher order accommodates sharp throttle transitions
% with far smaller defects than trapezoidal, so the optimizer is no longer
% forced to smear the throttle -- the path to a crisp bang-bang solution.
%
% INPUTS:  as SOLVE_MINFUEL_NLP_UE. Z0 is the HS decision vector
%          [X(:); U(:); Xc(:); Uc(:)] (see NLP_CONSTRAINTS_MINFUEL_HS).
%   maxIter - (optional) fmincon MaxIterations [default 5000]
%
% OUTPUTS:
%   Z, out - out has .X, .U (endpoints), .Xc, .Uc (midpoints), .mf,
%            .exitflag, .maxDefect (all constraints), .maxDynDefect,
%            .maxUnit, .eqnonlin
%
% REFERENCES:
%   [1] Betts, SIAM 2010, Ch. 4.  [2] Kelly, SIAM Review 59(4), 2017.

if nargin < 10 || isempty(maxIter), maxIter = 5000; end
sigma = sigma(:);
N  = numel(sigma) - 1;
nN = N + 1;
oX = 0;  oU = 7*nN;  oXc = 7*nN + 4*nN;  oUc = oXc + 7*N;
nZ = oUc + 4*N;

assert(numel(Z0) == nZ, 'solve_minfuel_nlp_hs:z0', ...
       'Z0 has %d elements; expected 11*nN + 11*N = %d', numel(Z0), nZ);
assert(sigma(1) == 0 && sigma(end) == 1 && all(diff(sigma) > 0), ...
       'solve_minfuel_nlp_hs:sigma', 'sigma must increase strictly 0 -> 1');

xLo = [-3; -3; -3; -12; -12; -12; 0.3];  xHi = [3; 3; 3; 12; 12; 12; 1.0];
uLo = [-1.1; -1.1; -1.1; 0];             uHi = [1.1; 1.1; 1.1; 1];

lb = -inf(nZ, 1);  ub = inf(nZ, 1);
for k = 1:nN
    xi = oX + (k-1)*7 + (1:7);  lb(xi) = xLo;  ub(xi) = xHi;
    ui = oU + (k-1)*4 + (1:4);  lb(ui) = uLo;  ub(ui) = uHi;
end
for k = 1:N
    xi = oXc + (k-1)*7 + (1:7);  lb(xi) = xLo;  ub(xi) = xHi;
    ui = oUc + (k-1)*4 + (1:4);  lb(ui) = uLo;  ub(ui) = uHi;
end
lb(oX+(1:7))          = [rv0(:); m0];  ub(oX+(1:7))          = [rv0(:); m0];
xfIdx = oX + (nN-1)*7 + (1:6);
lb(xfIdx) = rvf(:);   ub(xfIdx) = rvf(:);

idxMf = 7*nN;                        % X(7, nN) = final mass
gradJ = sparse(nZ, 1);  gradJ(idxMf) = -1;
objFun = @(Z) deal(-Z(idxMf), gradJ);
conFun = @(Z) nlp_constraints_minfuel_hs(Z, sigma, tf, Tmax, c, muStar);

opts = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true, ...
    'InitBarrierParam', 1e-3, ...
    'HessianApproximation', 'lbfgs', ...
    'MaxIterations', maxIter, ...
    'MaxFunctionEvaluations', 5e5, ...
    'ConstraintTolerance', 1e-11, ...
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance', 1e-14, ...
    'Display', 'iter');

[Z, ~, exitflag, output, lambdaS] = fmincon(objFun, Z0, [], [], [], [], ...
                                            lb, ub, conFun, opts);

Xr  = reshape(Z(oX  + (1:7*nN)), 7, nN);
Ur  = reshape(Z(oU  + (1:4*nN)), 4, nN);
Xcr = reshape(Z(oXc + (1:7*N )), 7, N );
Ucr = reshape(Z(oUc + (1:4*N )), 4, N );
[~, ceq] = nlp_constraints_minfuel_hs(Z, sigma, tf, Tmax, c, muStar);

out = struct('X', Xr, 'U', Ur, 'Xc', Xcr, 'Uc', Ucr, 'mf', Xr(7,end), ...
             'exitflag', exitflag, 'maxDefect', max(abs(ceq)), ...
             'maxDynDefect', max(abs(ceq(1:14*N))), ...
             'maxUnit', max(abs(ceq(14*N+1:end))), 'fmincon', output, ...
             'eqnonlin', lambdaS.eqnonlin);
end
