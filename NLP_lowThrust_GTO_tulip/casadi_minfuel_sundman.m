function out = casadi_minfuel_sundman(sigma, tf, rv0, rvf, Tmax, c, muStar, X0, U0, tauf0, pSund, maxIter, epsilon)
% CASADI_MINFUEL_SUNDMAN  Sundman-regularized min-fuel collocation (CasADi+IPOPT).
%
% Path A: fix the near-perigee ill-conditioning that stalled the plain
% collocation by changing the independent variable from time t to a
% regularizing variable tau, with the Sundman relation
%     dt/dtau = kappa(r) = r1^pSund,     r1 = ||r - r_Earth||.
% Every state ODE is multiplied by kappa, so the near-perigee gravity Hessian
% terms (~1/r^3) that blew up IPOPT's exact Hessian become r1^(pSund-3)
% (bounded for pSund ~ 3, mild for pSund ~ 1.5). A UNIFORM mesh in tau also
% concentrates nodes near perigee in TIME automatically. Time t is carried as
% an 8th state; the fixed transfer time is the terminal constraint
% t(tau_f) = tf. The total regularized length tau_f is held FIXED (= tauf0, the
% warm-start value), NOT a decision variable: a free scalar tau_f multiplies
% every collocation defect, producing one dense KKT column -> catastrophic
% MUMPS fill-in / OOM at large N. Fixing tau_f and enforcing t(tau_f)=tf on the
% carried time state instead keeps the Jacobian sparse; the trajectory adjusts
% so that int(kappa dtau) = tf.
%
% The objective is the Bertrand-Epenoy energy->fuel homotopy in epsilon:
%   J(eps) = Int[s]dt - eps*Int[s(1-s)]dt   (physical-time measure dt=kappa dtau)
%   eps=1 -> Int[s^2]dt (energy, strictly convex, smooth ramp)
%   eps=0 -> Int[s]dt   (fuel, linear -> bang-bang; equals propellant up to a
%                        positive constant, since m(tf)=1-(Tmax/c)Int[s]dt).
% Sweep eps 1->0, warm-starting each solve from the last (see run_sundman_*).
%
% State  x = [r(3); v(3); m; t]  (8).   Control u = [alpha(3); s]  (4).
% Cone-eliminated: thrust = s*Tmax*alpha/m, ||alpha|| = 1, s in [0,1].
%
% INPUTS:
%   sigma   - normalized independent-variable nodes [(N+1)x1], 0 -> 1
%   tf      - fixed transfer TIME (ND) [scalar]
%   rv0,rvf - initial / target position-velocity (ND) [1x6]
%   Tmax,c,muStar - dynamics constants [scalars]
%   X0      - warm-start states [8x(N+1)] ([r;v;m;t])
%   U0      - warm-start controls [4x(N+1)] ([alpha;s])
%   tauf0   - fixed total regularized length [scalar]
%   pSund   - Sundman power [scalar, default 1.5]
%   maxIter - IPOPT max iterations [scalar, default 3000]
%   epsilon - homotopy parameter in [0,1]: 0=fuel, 1=energy [scalar, default 0]
%
% OUTPUTS:
%   out - struct: .X [8x(N+1)] .U [4x(N+1)] .tauf .mf .maxDefect .maxUnit
%         .switches .edge (bang-bang node fraction) .success .ipoptStatus
%
% REFERENCES:
%   [1] Bertrand & Epenoy, "New smoothing techniques for solving bang-bang
%       optimal control problems," Optim. Control Appl. Methods 23 (2002).
%   [2] Sundman regularization of the two/three-body problem; e.g. dt = r dtau.
%   [3] Andersson et al., "CasADi," Math. Prog. Comp. 11 (2019); Wachter &
%       Biegler (IPOPT), Math. Prog. 106 (2006).

if nargin < 11 || isempty(pSund),  pSund  = 1.5;  end
if nargin < 12 || isempty(maxIter), maxIter = 3000; end
if nargin < 13 || isempty(epsilon), epsilon = 0;   end   % 0=fuel, 1=energy
addpath('/Users/msc/casadi-3.7.0');
import casadi.*

sigma = sigma(:);  N = numel(sigma) - 1;  nN = N + 1;
dsig  = diff(sigma).';                       % 1 x N

% --- symbolic Sundman-regularized dynamics dX/dtau = kappa * f -------------
x = MX.sym('x', 8);  u = MX.sym('u', 4);
r = x(1:3);  v = x(4:6);  m = x(7);  al = u(1:3);  s = u(4);
dd = [r(1)+muStar; r(2); r(3)];              % vector from Earth
rr = [r(1)-1+muStar; r(2); r(3)];            % vector from Moon
r1 = sqrt(dd.'*dd + 1e-12);                  % Earth distance (guarded)
d3 = (dd.'*dd + 1e-12)^1.5;  r3 = (rr.'*rr + 1e-12)^1.5;   % guarded denoms
gr = [r(1); r(2); 0] - (1-muStar)*dd/d3 - muStar*rr/r3;
hv = [2*v(2); -2*v(1); 0];
accel = gr + hv + (s*Tmax/m)*al;
mdot  = -(Tmax/c)*s;
kappa = r1^pSund;
fdyn  = Function('f', {x,u}, {kappa*[v; accel; mdot; 1]});   % dX/dtau (8x1)
Fmap  = fdyn.map(nN);
% energy->fuel homotopy integrands (Bertrand-Epenoy), on the physical-time
% measure dt = kappa dtau:  q_fuel = s*kappa,  q_smooth = s(1-s)*kappa.
gint  = Function('g', {x,u}, {[s*kappa; s*(1-s)*kappa]});
Gmap  = gint.map(nN);

% --- NLP ------------------------------------------------------------------
opti = Opti();
X    = opti.variable(8, nN);
U    = opti.variable(4, nN);
% tau-length is FIXED (from the warm start), NOT a decision variable: a free
% scalar tau_f couples to every defect -> a dense KKT column -> catastrophic
% MUMPS fill-in / OOM at large N. Fixed transfer time is still enforced
% exactly by the t-state terminal condition t(tau_end) = tf below; the
% trajectory adjusts so that int(kappa dtau) = tf.
tauf = tauf0;
F    = Fmap(X, U);                           % 8 x nN, = dX/dtau

% trapezoidal defects in sigma: dX/dsigma = tauf * dX/dtau
D = X(:,2:end) - X(:,1:end-1) - tauf*(repmat(dsig,8,1)/2).*(F(:,1:end-1) + F(:,2:end));
opti.subject_to(D(:) == 0);

% unit-direction
opti.subject_to((sum(U(1:3,:).^2, 1) - 1).' == 0);

% bounds (explicit two-sided)
lbX = repmat([-3;-3;-3;-12;-12;-12;0.3;0], 1, nN);
ubX = repmat([ 3; 3; 3; 12; 12; 12;1.0; 2*tf], 1, nN);
opti.subject_to(X(:) >= lbX(:));   opti.subject_to(X(:) <= ubX(:));
lbU = repmat([-1.1;-1.1;-1.1;0], 1, nN);
ubU = repmat([ 1.1; 1.1; 1.1;1], 1, nN);
opti.subject_to(U(:) >= lbU(:));   opti.subject_to(U(:) <= ubU(:));

% boundary conditions (fixed transfer TIME via the t-state)
opti.subject_to(X(1:6,1) == rv0(:));   opti.subject_to(X(7,1) == 1);   opti.subject_to(X(8,1) == 0);
opti.subject_to(X(1:6,nN) == rvf(:));  opti.subject_to(X(8,nN) == tf);

% objective: J(eps) = Int[s]dt - eps*Int[s(1-s)]dt  (trapezoid in tau)
%   eps=0 -> Int s dt   = fuel (linear in s -> bang-bang)
%   eps=1 -> Int s^2 dt = energy (strictly convex -> smooth ramp, no restoration)
G    = Gmap(X, U);                            % 2 x nN  [q_fuel; q_smooth]
IntF = tauf * sum((dsig/2).*(G(1,1:end-1) + G(1,2:end)));
IntS = tauf * sum((dsig/2).*(G(2,1:end-1) + G(2,2:end)));
opti.minimize(IntF - epsilon*IntS);
opti.set_initial(X, X0);
opti.set_initial(U, U0);

p = struct;
p.print_time      = true;
p.ipopt.max_iter  = maxIter;
p.ipopt.tol       = 1e-7;
p.ipopt.constr_viol_tol = 1e-7;
p.ipopt.acceptable_tol  = 1e-5;
p.ipopt.acceptable_iter = 15;
p.ipopt.mu_strategy = 'monotone';
p.ipopt.nlp_scaling_method = 'gradient-based';
p.ipopt.linear_solver = 'mumps';
p.ipopt.print_level = 5;
% Warm-start from a bang-bang seed: the throttle is PINNED at its bounds
% (s = 0 on coasts, s = 1 on burns). IPOPT's default bound_push (0.01) shoves
% every s off its bound at startup -- a large disruption that immediately
% triggers the restoration phase and a false "locally infeasible" exit. Keep
% the initial point hugging its bounds and start from a small barrier so the
% seed is honored.
p.ipopt.warm_start_init_point       = 'yes';
p.ipopt.mu_init                     = 1e-4;
p.ipopt.warm_start_bound_push       = 1e-9;
p.ipopt.warm_start_bound_frac       = 1e-9;
p.ipopt.warm_start_slack_bound_push = 1e-9;
p.ipopt.warm_start_slack_bound_frac = 1e-9;
p.ipopt.warm_start_mult_bound_push  = 1e-9;
opti.solver('ipopt', p);

success = true;  status = 'solved';
try
    sol = opti.solve();
    Xs = sol.value(X);  Us = sol.value(U);
    status = char(opti.return_status());
catch solveErr
    Xs = opti.debug.value(X);  Us = opti.debug.value(U);
    success = false;  status = solveErr.message;
end

% metrics
Fs = full(Fmap(Xs, Us));
Dd = Xs(:,2:end) - Xs(:,1:end-1) - tauf*(repmat(dsig,8,1)/2).*(Fs(:,1:end-1) + Fs(:,2:end));
ss = Us(4,:);
out = struct('X', Xs, 'U', Us, 'tauf', tauf, 'mf', Xs(7,end), ...
             'maxDefect', max(abs(Dd(:))), ...
             'maxUnit', max(abs(sum(Us(1:3,:).^2,1) - 1)), ...
             'switches', sum(abs(diff(ss > 0.5))), ...
             'edge', mean(ss > 0.95 | ss < 0.05), ...
             'success', success, 'ipoptStatus', status);
end
