function ok = test_stm_variational()
%% Purpose:
%
%   INDEPENDENT check of the variational generator behind S4. X2 validates
%   the 14-vector FIELD row by row; it says nothing about the 196 STM
%   entries the conjugate instruments chain, which need the derivative of
%   the OPTIMISED Hamiltonian flow, including d alpha*/d lam_v =
%   -(I - alpha alpha')/|lam_v| (Astra review #2, 2026-09-11). Two checks
%   on one shooting segment of the certified anchor:
%     1. DIRECTIONAL DERIVATIVES: every column of PHI against a central
%        finite difference of two flights, to the FD's own accuracy;
%     2. SYMPLECTICITY: the 14-state PMP flow is canonical in ([r v m],
%        [lam_r lam_v lam_m]), so PHI' J PHI = J with J = [0 I; -I 0], to
%        integration accuracy. Measured 2026-09-11: this does NOT
%        discriminate the control-law derivative -- a FROZEN-control
%        generator is the linearisation of H(x, lam, alpha) at fixed alpha,
%        itself Hamiltonian, and is symplectic to the same 1e-11. It is a
%        check of the integration and of the canonical structure only.
%     3. DISCRIMINATION: the frozen-control STM is integrated beside the
%        true one and the two must DIFFER in the lam_v columns by O(1)
%        relative, while the FD comparison (1) sides with the true one.
%        That is what shows the generator carries d alpha*/d lam_v.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));
lib = dro_tulip_library();
E = lib(strcmp({lib.src}, 'anchor'));  E = E(1);
[B, ~] = arclength_arrival('setup');
rv0 = B.stateD(E.sD);  z = E.z(:);
y0 = [rv0(1:6); 1; z(1:7)];
dt = z(8)/24;                                  % one shooting segment

[yh, PHI] = mintime_prop_seg(dt, y0, true, B.Tnd, B.cnd, B.mu);
ok = chk(ok, isequal(size(PHI), [14 14]) && all(isfinite(PHI(:))), 'a 14x14 finite STM is returned');

% 1. directional derivatives, column by column
h = 1e-6*max(1, abs(y0));
errCol = zeros(1, 14);
for jj = 1:14
    e = zeros(14, 1);  e(jj) = h(jj);
    yp = mintime_prop_seg(dt, y0 + e, false, B.Tnd, B.cnd, B.mu);
    ym = mintime_prop_seg(dt, y0 - e, false, B.Tnd, B.cnd, B.mu);
    fd = (yp - ym)/(2*h(jj));
    errCol(jj) = norm(fd - PHI(:, jj))/max(norm(PHI(:, jj)), 1);
end
ok = chk(ok, max(errCol) < 1e-6, ...
         sprintf('all 14 columns match central differences: worst relative %.1e (FD accuracy ~1e-8..1e-6)', max(errCol)));
ok = chk(ok, max(errCol(11:13)) < 1e-6, ...
         sprintf('including the lam_v columns that carry d alpha*/d lam_v: %.1e', max(errCol(11:13))));

% 2. symplecticity of the canonical 14-state flow
J = [zeros(7), eye(7); -eye(7), zeros(7)];
sympErr = norm(PHI.'*J*PHI - J)/norm(J);
ok = chk(ok, sympErr < 1e-8, sprintf('PHI'' J PHI = J to %.1e', sympErr));

% 3. a generator WITHOUT the control-law derivative: integrate the
%    variational equation with the FROZEN-control Jacobian (the one the lift
%    space uses). It is symplectic too (a fixed-alpha Hamiltonian), so that
%    identity cannot tell the two apart; the lam_v columns can.
Af = frozenJacobian();
tspan = [0 dt];
rhs = @(t, w) [mintime_rhs_point(w(1:14), B.Tnd, B.cnd, B.mu); ...
               reshape(full(Af(w(1:7), w(8:14), -w(11:13)/norm(w(11:13)), B.Tnd, B.cnd, B.mu)) * reshape(w(15:end), 14, 14), [], 1)];
[~, W] = ode113(rhs, tspan, [y0; reshape(eye(14), [], 1)], odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
PHIf = reshape(W(end, 15:end), 14, 14);
sympErrF = norm(PHIf.'*J*PHIf - J)/norm(J);
ok = chk(ok, sympErrF < 1e-8, ...
         sprintf('a frozen-control generator is ALSO symplectic (%.1e): symplecticity does not discriminate', sympErrF));
dLam = norm(PHIf(:, 11:13) - PHI(:, 11:13))/norm(PHI(:, 11:13));
dPos = norm(PHIf(:, 1:3) - PHI(:, 1:3))/norm(PHI(:, 1:3));
ok = chk(ok, dLam > 1e-3, ...
         sprintf('but its lam_v columns differ from the true STM by %.2e relative (position columns %.1e): the true generator carries d alpha*/d lam_v, and FD sided with it', dLam, dPos));

if ok, fprintf('TEST_STM_VARIATIONAL: ALL PASS\n'); else, fprintf('TEST_STM_VARIATIONAL: FAIL\n'); end
end

function Af = frozenJacobian()
% FROZENJACOBIAN  CasADi Jacobian of the full 14-state PMP field with the
% control FROZEN at an input alpha: d[f(x,alpha); -H_x(x,lam,alpha)] /
% d[x; lam]. This is the generator a lift-space calculation uses; it omits
% d alpha*/d lam_v and so is NOT the derivative of the optimised flow.
% INPUTS: none.  OUTPUTS: Af (CasADi Function of x, lam, alpha, T, c, mu).
if isempty(which('casadi.SX')), addpath(fullfile(getenv('HOME'), 'casadi-3.7.0')); end
xs = casadi.SX.sym('x', 7);  ls = casadi.SX.sym('l', 7);  as = casadi.SX.sym('a', 3);
Ts = casadi.SX.sym('T');  cs = casadi.SX.sym('c');  mus = casadi.SX.sym('mu');
r = xs(1:3);  v = xs(4:6);  m = xs(7);
dd = sqrt((r(1)+mus)^2 + r(2)^2 + r(3)^2);
rr = sqrt((r(1)-1+mus)^2 + r(2)^2 + r(3)^2);
gr = [r(1) - (1-mus)*(r(1)+mus)/dd^3 - mus*(r(1)-1+mus)/rr^3;
      r(2) - (1-mus)*r(2)/dd^3       - mus*r(2)/rr^3;
           - (1-mus)*r(3)/dd^3       - mus*r(3)/rr^3];
hv = [2*v(2); -2*v(1); 0];
f7 = [v; gr + hv + (Ts/m)*as; -Ts/cs];
H = 1 + ls.'*f7;
F14 = [f7; -jacobian(H, xs).'];
Af = casadi.Function('Afrozen', {xs, ls, as, Ts, cs, mus}, {jacobian(F14, [xs; ls])}).expand();
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
