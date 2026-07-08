function [Z, out] = solve_tfmin_nlp(Z0, sigma, rv0, rvf, Tmax, c, muStar)
% SOLVE_TFMIN_NLP  Direct min-time solve via fmincon (interior-point).
%
% Minimizes tf subject to the trapezoidal defect equalities, the
% unit-sphere control constraint w'*w = 1, loose safety bounds, and
% endpoint pinning (via lb = ub). All bounds are strictly inactive at the
% warm start (deliberate: an interior-point method forced off an active
% bound at the start will wander).
% Objective and constraint gradients are analytic and sparse; the Hessian
% is limited-memory BFGS (a dense quasi-Newton Hessian is infeasible at
% ~10^4 variables).
%
% INPUTS:
%   Z0     - initial decision vector [10*(N+1)+1 x 1] (see BUILD_GUESS)
%   sigma  - normalized node times [(N+1)x1] (see BUILD_GUESS)
%   rv0    - initial position/velocity (ND) [1x6]
%   rvf    - target position/velocity (ND) [1x6]
%   Tmax   - max thrust acceleration at m = 1 (ND) [scalar]
%   c      - exhaust velocity (ND) [scalar]
%   muStar - Earth-Moon mass ratio [scalar]
%
% OUTPUTS:
%   Z      - converged decision vector
%   out    - struct: .tf, .X [7x(N+1)], .W [3x(N+1)], .exitflag,
%            .maxDefect (inf-norm of defects at Z), .fmincon (solver output)
%
% REFERENCES:
%   [1] Betts, "Practical Methods for Optimal Control and Estimation Using
%       Nonlinear Programming," 2nd ed., SIAM, 2010.

N      = numel(sigma) - 1;
nNodes = N + 1;
nZ     = 10*nNodes + 1;

assert(numel(Z0) == nZ, 'solve_tfmin_nlp:z0', ...
       'Z0 has %d elements; expected 10*(N+1)+1 = %d', numel(Z0), nZ);
assert(sigma(1) == 0 && sigma(end) == 1 && all(diff(sigma(:)) > 0), ...
       'solve_tfmin_nlp:sigma', ...
       'sigma must increase strictly from 0 to 1');

% --- bounds ---------------------------------------------------------------
lb = -inf(nZ, 1);  ub = inf(nZ, 1);
for k = 1:nNodes
    xIdx = (k-1)*7 + (1:7);
    lb(xIdx) = [-3; -3; -3; -12; -12; -12; 0.3];
    ub(xIdx) = [ 3;  3;  3;  12;  12;  12; 1.0];
    uIdx = 7*nNodes + (k-1)*3 + (1:3);
    lb(uIdx) = [-1.5; -1.5; -1.5];
    ub(uIdx) = [ 1.5;  1.5;  1.5];
end
% pin the endpoints (cleaner than extra equality constraints)
lb(1:7)   = [rv0(:); 1];  ub(1:7)   = [rv0(:); 1];
xfIdx     = (nNodes-1)*7 + (1:6);
lb(xfIdx) = rvf(:);       ub(xfIdx) = rvf(:);
lb(end)   = 2;            ub(end)   = 15;          % tf window (ND)

% --- objective: J = tf, with sparse gradient ------------------------------
gradJ = sparse(nZ, 1);  gradJ(end) = 1;
objFun = @(Z) deal(Z(end), gradJ);

conFun = @(Z) nlp_constraints(Z, sigma, Tmax, c, muStar);

% InitBarrierParam: we start from a near-feasible warm start; the default
% (0.1) lets the barrier push iterates far off it (observed: tf plunges
% while feasibility decays). A small initial barrier keeps the IPM local.
opts = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true, ...
    'InitBarrierParam', 1e-6, ...
    'HessianApproximation', 'lbfgs', ...
    'MaxIterations', 3000, ...
    'MaxFunctionEvaluations', 1e5, ...
    'ConstraintTolerance', 1e-10, ...
    'OptimalityTolerance', 1e-7, ...
    'StepTolerance', 1e-12, ...
    'Display', 'iter');

[Z, ~, exitflag, output] = fmincon(objFun, Z0, [], [], [], [], lb, ub, ...
                                   conFun, opts);

[X, W, tf] = unpack_z(Z, N);
[~, ceq]   = nlp_constraints(Z, sigma, Tmax, c, muStar);

out = struct('tf', tf, 'X', X, 'W', W, 'exitflag', exitflag, ...
             'maxDefect', max(abs(ceq(1:7*N))), 'fmincon', output);
end
