function [Z, out] = solve_energy_nlp(Z0, sigma, tf, rv0, m0, rvf, Tmax, c, muStar)
% SOLVE_ENERGY_NLP  Direct min-ENERGY solve via fmincon (interior-point).
%
% Minimizes the quadratic control energy J = int (1/2) s^2 dt (trapezoidal
% quadrature) subject to the SAME transcription used by min-fuel:
% trapezoidal defects, the throttle cone w'*w = s^2 with s in [0,1],
% endpoint pinning, and fixed transfer time. Only the OBJECTIVE differs
% from SOLVE_MINFUEL_NLP (energy is smooth and strictly convex in s, so the
% solution has NO bound-active bang-bang structure -- the throttle rides a
% saturated ramp and the NLP is far better conditioned).
%
% INPUTS:
%   Z0     - initial decision vector [11*(N+1) x 1]
%   sigma  - normalized node times [(N+1)x1]
%   tf     - fixed transfer time (ND) [scalar]
%   rv0    - initial position/velocity (ND) [1x6]
%   m0     - initial mass fraction [scalar]
%   rvf    - target position/velocity (ND) [1x6]
%   Tmax   - max thrust acceleration at m = 1 (ND) [scalar]
%   c      - exhaust velocity (ND) [scalar]
%   muStar - Earth-Moon mass ratio [scalar]
%
% OUTPUTS:
%   Z    - converged decision vector
%   out  - struct: .X [7x(N+1)], .U [4x(N+1)], .mf, .energy (J value),
%          .exitflag, .maxDefect, .fmincon, .eqnonlin (KKT multipliers)
%
% REFERENCES:
%   [1] Betts, SIAM 2010, Ch. 4.

sigma  = sigma(:);
N      = numel(sigma) - 1;
nNodes = N + 1;
nZ     = 11*nNodes;

assert(numel(Z0) == nZ, 'solve_energy_nlp:z0', ...
       'Z0 has %d elements; expected 11*(N+1) = %d', numel(Z0), nZ);
assert(sigma(1) == 0 && sigma(end) == 1 && all(diff(sigma) > 0), ...
       'solve_energy_nlp:sigma', 'sigma must increase strictly 0 -> 1');

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

% --- objective: J = sum (h_k/2)(0.5 s_k^2 + 0.5 s_{k+1}^2) -----------------
% trapezoidal node weights omega_j (endpoints half); gradient = omega.*s in
% the s-slots only.
dSig  = diff(sigma).';
h     = tf.*dSig;                       % 1 x N segment widths
omega = zeros(nNodes, 1);
omega(1:N)     = omega(1:N)     + h(:)/2;
omega(2:N+1)   = omega(2:N+1)   + h(:)/2;   % node weights (trapezoid)
sIdx  = 7*nNodes + (1:nNodes)*4;            % the s-slot of each control node
    function [J, g] = energyObj(Z)
        s = Z(sIdx);
        J = 0.5*sum(omega.*s.^2);
        if nargout > 1
            g = sparse(nZ, 1);
            g(sIdx) = omega.*s;
        end
    end

conFun = @(Z) nlp_constraints_minfuel(Z, sigma, tf, Tmax, c, muStar);

% Standard interior-point returns a good-enough solution to SEED the
% indirect (primer directions + rough scale) far faster than feasibility
% mode + CG, which grinds pathologically on this 40-rev problem. Tight
% feasibility is the indirect's job, not the NLP's.
opts = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true, ...
    'InitBarrierParam', 1e-4, ...
    'HessianApproximation', 'lbfgs', ...
    'MaxIterations', 1200, ...
    'MaxFunctionEvaluations', 1e5, ...
    'ConstraintTolerance', 1e-8, ...
    'OptimalityTolerance', 1e-7, ...
    'StepTolerance', 1e-12, ...
    'Display', 'iter');

[Z, ~, exitflag, output, lambdaS] = fmincon(@energyObj, Z0, [], [], [], [], ...
                                            lb, ub, conFun, opts);

X = reshape(Z(1:7*nNodes), 7, nNodes);
U = reshape(Z(7*nNodes + (1:4*nNodes)), 4, nNodes);
[~, ceq] = nlp_constraints_minfuel(Z, sigma, tf, Tmax, c, muStar);
Jval = energyObj(Z);

out = struct('X', X, 'U', U, 'mf', X(7, end), 'energy', Jval, ...
             'exitflag', exitflag, 'maxDefect', max(abs(ceq(1:7*N))), ...
             'fmincon', output, 'eqnonlin', lambdaS.eqnonlin);
end
