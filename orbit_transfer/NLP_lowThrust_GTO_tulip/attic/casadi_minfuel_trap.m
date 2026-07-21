function out = casadi_minfuel_trap(sigma, tf, rv0, m0, rvf, Tmax, c, muStar, X0, U0, maxIter)
% CASADI_MINFUEL_TRAP  Cone-eliminated trapezoidal min-fuel via CasADi + IPOPT.
%
% Re-poses the min-fuel NLP symbolically in CasADi and solves with IPOPT
% (exact sparse Jacobian AND Hessian via automatic differentiation) -- the
% exact-Hessian solver that fmincon interior-point + lbfgs could not match
% on the many-switch bang-bang problem. Control is [alpha(3); s] with the
% unit-direction constraint ||alpha|| = 1 (cone-eliminated) and s in [0,1].
%
% INPUTS:
%   sigma   - normalized node times [(N+1)x1], 0 -> 1
%   tf      - fixed transfer time (ND) [scalar]
%   rv0     - initial position/velocity (ND) [1x6]
%   m0      - initial mass fraction [scalar]
%   rvf     - target position/velocity (ND) [1x6]
%   Tmax,c,muStar - dynamics constants
%   X0, U0  - warm-start node states/controls [7x(N+1)], [4x(N+1)] (U=[alpha;s])
%   maxIter - (optional) IPOPT max iterations [default 3000]
%
% OUTPUTS:
%   out - struct: .X, .U, .mf, .maxDefect, .maxUnit, .switches, .edge
%         (fraction of nodes at a throttle bound = bang-bang sharpness),
%         .success, .ipoptStatus
%
% REFERENCES:
%   [1] Andersson et al., "CasADi," Math. Prog. Comp. 11(1), 2019.
%   [2] Waechter & Biegler, "IPOPT," Math. Prog. 106(1), 2006.

if nargin < 11 || isempty(maxIter), maxIter = 3000; end
addpath('/Users/msc/casadi-3.7.0');
import casadi.*

sigma = sigma(:);  N = numel(sigma) - 1;  nN = N + 1;
h = (tf*diff(sigma)).';                       % 1 x N

% --- symbolic dynamics (cone-eliminated: thrust = s*Tmax*alpha/m) ----------
x = MX.sym('x', 7);  u = MX.sym('u', 4);
r = x(1:3);  v = x(4:6);  m = x(7);  al = u(1:3);  s = u(4);
dd = [r(1)+muStar; r(2); r(3)];
rr = [r(1)-1+muStar; r(2); r(3)];
d3 = (dd.'*dd)^1.5;  r3 = (rr.'*rr)^1.5;
gr = [r(1); r(2); 0] - (1-muStar)*dd/d3 - muStar*rr/r3;
hv = [2*v(2); -2*v(1); 0];
xdot = [v; gr + hv + (s*Tmax/m)*al; -(Tmax/c)*s];
fdyn = Function('f', {x, u}, {xdot});
Fmap = fdyn.map(nN);

% --- NLP via Opti ----------------------------------------------------------
opti = Opti();
X = opti.variable(7, nN);
U = opti.variable(4, nN);
F = Fmap(X, U);                               % 7 x nN

% trapezoidal defects
D = X(:,2:end) - X(:,1:end-1) - (repmat(h,7,1)/2).*(F(:,1:end-1) + F(:,2:end));
opti.subject_to(D(:) == 0);

% unit-direction constraints
unit = sum(U(1:3,:).^2, 1) - 1;
opti.subject_to(unit(:) == 0);

% bounds -- explicit two-sided (MATLAB does NOT chain a<=x<=b correctly;
% write each side so the box constraints are unambiguous)
lbX = repmat([-3;-3;-3;-12;-12;-12;0.3], 1, nN);
ubX = repmat([ 3; 3; 3; 12; 12; 12;1.0], 1, nN);
opti.subject_to(X(:) >= lbX(:));
opti.subject_to(X(:) <= ubX(:));
lbU = repmat([-1.1;-1.1;-1.1;0], 1, nN);
ubU = repmat([ 1.1; 1.1; 1.1;1], 1, nN);
opti.subject_to(U(:) >= lbU(:));
opti.subject_to(U(:) <= ubU(:));

% boundary conditions
opti.subject_to(X(:,1)  == [rv0(:); m0]);
opti.subject_to(X(1:6,nN) == rvf(:));

opti.minimize(-X(7, nN));                     % maximize final mass
opti.set_initial(X, X0);
opti.set_initial(U, U0);

p = struct;
p.print_time      = true;
p.ipopt.max_iter  = maxIter;
p.ipopt.tol       = 1e-7;
p.ipopt.constr_viol_tol = 1e-8;
p.ipopt.acceptable_tol  = 1e-5;
p.ipopt.acceptable_iter = 15;
p.ipopt.mu_strategy = 'monotone';             % stable barrier schedule
p.ipopt.nlp_scaling_method = 'gradient-based';
p.ipopt.linear_solver = 'mumps';
p.ipopt.print_level = 5;
opti.solver('ipopt', p);

success = true;  status = '';
try
    sol = opti.solve();
    Xs = sol.value(X);  Us = sol.value(U);
    status = char(opti.return_status());
catch solveErr
    Xs = opti.debug.value(X);  Us = opti.debug.value(U);   % last iterate
    success = false;  status = solveErr.message;
end

% --- metrics ---------------------------------------------------------------
Fs = full(Fmap(Xs, Us));
Dd = Xs(:,2:end) - Xs(:,1:end-1) - (repmat(h,7,1)/2).*(Fs(:,1:end-1) + Fs(:,2:end));
ss = Us(4,:);
out = struct('X', Xs, 'U', Us, 'mf', Xs(7,end), ...
             'maxDefect', max(abs(Dd(:))), ...
             'maxUnit', max(abs(sum(Us(1:3,:).^2,1) - 1)), ...
             'switches', sum(abs(diff(ss > 0.5))), ...
             'edge', mean(ss > 0.95 | ss < 0.05), ...
             'success', success, 'ipoptStatus', status);
end
