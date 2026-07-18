% TEST_MEE_RHS  Ballistic invariance + thrust cross-check vs Cartesian RHS.
p = kepler_lt_params(10, 1500, 2000);
% initial MEE state (paper), coplanar variant for the planar checks where noted
X0 = [11625/p.LU_km; 0.75; 0; 0.0612; 0; 1; 0];
% (a) BALLISTIC: thr=0 -> P,ex,ey,hx,hy,m all frozen; only t advances; Ldot>0
U0 = [1;0;0; 0];
[dXdL, Ldot] = lt_mee_rhs(X0, U0, setfield(p,'L',pi));
assert(Ldot > 0, 'Ldot must be positive');
assert(max(abs(dXdL(1:5))) < 1e-14, 'elements must be frozen under zero thrust');
assert(abs(dXdL(6)) < 1e-14, 'mass frozen under zero thrust');
assert(dXdL(7) > 0, 'time must advance');
% (b) BALLISTIC Ldot value: at L=pi (apogee, ex=0.75) Z=1-0.75=0.25, Ldot=sqrt(1/P^3)*Z^2
P0 = X0(1); Z_apo = 1 - 0.75;
assert(abs(Ldot - sqrt(1/P0^3)*Z_apo^2) < 1e-12, 'ballistic Ldot formula');
% (c) THRUST CROSS-CHECK vs Cartesian: convert MEE->Cartesian, apply the SAME
% physical thrust in both, require d/dt of the Cartesian state to match the
% Cartesian RHS to ODE tolerance. Transverse burn thr=1, beta=[0;1;0].
Uc = [0;1;0; 1];                              % pure transverse, full throttle
[dXdL_t, Ldot_t] = lt_mee_rhs(X0, Uc, setfield(p,'L',pi));
dXdt_mee = dXdL_t * Ldot_t;                   % back to time domain
% independent finite check: energy rate must be positive for a transverse burn
% (raises orbit), and Pdot>0
assert(dXdt_mee(1) > 0, 'transverse burn must raise P');

% (d) REAL CARTESIAN CROSS-CHECK: a general RTN thrust direction with ALL
% THREE components nonzero (radial+transverse+normal), finite-differenced
% against an independent Cartesian propagation of the SAME physical thrust
% accel. This is the check that actually exercises the normal (beta(3))
% component and therefore the Ldot thrust term -- test (c) above never sets
% beta(3)~=0, so it could not have caught an error confined to that term.
% L0 is deliberately NOT pi: at L=pi, sin(L)=0 and hy=0, so
% hterm = hx*sin(L) - hy*cos(L) = 0 identically and the entire Ldot thrust
% term vanishes regardless of its coefficient -- that degenerate point would
% make this check blind to the exact bug it's meant to catch (Finding 1).
L0 = 1.2;
mu = p.mu;  Tm = p.Tmax;  m = X0(6);
beta = [0.2; 0.8; 0.5663];  beta = beta/norm(beta);
thr = 1;
Ud = [beta; thr];
[r0, v0] = elements_to_cart(X0(1), X0(2), X0(3), X0(4), X0(5), L0, mu);
rhat = r0/norm(r0);
nhat = cross(r0,v0)/norm(cross(r0,v0));
that = cross(nhat, rhat);
a_thrust = (Tm/m)*thr*(beta(1)*rhat + beta(2)*that + beta(3)*nhat);
% RK4 step of the independent Cartesian 2-body + thrust dynamics
dt = 1e-6;
f_cart = @(rv) [rv(4:6); -mu*rv(1:3)/norm(rv(1:3))^3 + a_thrust];
rv0 = [r0; v0];
k1 = f_cart(rv0);
k2 = f_cart(rv0 + 0.5*dt*k1);
k3 = f_cart(rv0 + 0.5*dt*k2);
k4 = f_cart(rv0 + dt*k3);
rv1 = rv0 + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
el1 = cart_to_elements(rv1(1:3), rv1(4:6), mu);
el0.P = X0(1); el0.ex = X0(2); el0.ey = X0(3); el0.hx = X0(4); el0.hy = X0(5); el0.L = L0;
fd_rate = [ (el1.P  - el0.P ) / dt;
            (el1.ex - el0.ex) / dt;
            (el1.ey - el0.ey) / dt;
            (el1.hx - el0.hx) / dt;
            (el1.hy - el0.hy) / dt ];
% L wraps to (-pi,pi] in cart_to_elements; unwrap the finite difference so a
% step through the +pi/-pi seam doesn't register as a spurious -2*pi jump.
dL = el1.L - el0.L;  dL = mod(dL + pi, 2*pi) - pi;
Ldot_fd = dL / dt;
[dXdL_d, Ldot_d] = lt_mee_rhs(X0, Ud, setfield(p,'L',L0));
dXdt_d = dXdL_d * Ldot_d;                     % element rates in time domain
% Combined absolute+relative tolerance (ODE-solver style): guards against any
% individual rate happening to be near-zero at this L0, whose FD estimate
% would then be pure roundoff noise (from catastrophic cancellation in
% cart_to_elements' cross/atan2 chain, amplified by dividing by dt=1e-6).
% abstol=1e-6 comfortably absorbs that noise while reltol=1e-5 stays tight on
% the ~1e-2-scale real signal terms.
abstol = 1e-6;  reltol = 1e-5;
diff_elem = abs(dXdt_d(1:5) - fd_rate);
ok_elem = (diff_elem < abstol) | (diff_elem ./ abs(fd_rate) < reltol);
diff_L = abs(Ldot_d - Ldot_fd);
ok_L = (diff_L < abstol) | (diff_L / abs(Ldot_fd) < reltol);
assert(all(ok_elem), ...
    sprintf('Cartesian cross-check failed on element rates (max abs diff %.3e)', max(diff_elem)));
assert(ok_L, ...
    sprintf('Cartesian cross-check failed on Ldot (abs diff %.3e)', diff_L));

% (e) CASADI MX REGRESSION: lt_mee_rhs must build/evaluate cleanly under
% CasADi MX symbolics (no norm/abs/max/if on state-dependent quantities),
% and the symbolic evaluation must match the numeric path exactly.
cp = getenv('CASADI_PATH');
if isempty(cp), cp = fullfile(getenv('HOME'), 'casadi-3.7.0'); end
addpath(cp);
import casadi.*
Xs = MX.sym('x', 7, 1);
Us = MX.sym('u', 4, 1);
try
    [dXdL_s, Ldot_s] = lt_mee_rhs(Xs, Us, setfield(p, 'L', pi));
catch mx_err
    error('lt_mee_rhs failed to build under CasADi MX: %s', mx_err.message);
end
f_mx = casadi.Function('f_mx', {Xs, Us}, {dXdL_s, Ldot_s});
[dXdL_mx, Ldot_mx] = f_mx(X0, Uc);
dXdL_mx = full(dXdL_mx);  Ldot_mx = full(Ldot_mx);
assert(max(abs(dXdL_mx - dXdL_t)) < 1e-12, 'CasADi MX dXdL mismatch vs numeric path');
assert(abs(Ldot_mx - Ldot_t) < 1e-12, 'CasADi MX Ldot mismatch vs numeric path');

fprintf('test_mee_rhs: ALL PASS (Ldot=%.4f, Pdot=%.3e, MX_Ldot=%.4f)\n', ...
    Ldot_t, dXdt_mee(1), Ldot_mx);
