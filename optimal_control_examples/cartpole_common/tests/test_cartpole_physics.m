function ok = test_cartpole_physics()
%% Purpose:
%
%   THE INDEPENDENT PHYSICAL ORACLE for the cart-pole dynamics. Every other
%   test in this folder checks the field against itself or against the ex2
%   helpers, so an error shared by all of them is invisible -- and one was:
%   until 2026-09-17 the pendulum row carried the wrong sign in ex2 AND here
%   (found by GPT-6 Astra's review; the direct and indirect solves agreed to
%   0.075% because they solved the same wrong problem).
%
%   This test derives nothing from the acceleration formulas. It uses only
%   the GEOMETRY: the bob sits at (q1 + L sin q2, -L cos q2), so
%       T = 1/2 m1 q1dot^2 + 1/2 m2 |v_bob|^2,    V = -m2 g L cos q2,
%   and mechanics then demands, at EVERY state and control,
%       dE/dt = grad E . xdot = u * q1dot        (power in = force x velocity)
%   -- a pointwise identity, no integration needed. A sign error in any
%   acceleration row breaks it.
%
%   Checks:
%     1. POWER BALANCE for cartpole_field at a scatter of states, controls.
%     2. POWER BALANCE for the ex2 helpers (the tutorial is pinned by
%        physics too, not only by agreement with this folder).
%     3. HANGING IS STABLE: linearised about q2 = 0 the pendulum mode is
%        oscillatory, omega^2 = (m1+m2) g / (L m1).
%     4. UPRIGHT IS UNSTABLE: about q2 = pi the same mode is a real pair
%        +/- sqrt((m1+m2) g / (L m1)).
%     5. ENERGY IS CONSERVED along an unforced flight (u = 0), to integrator
%        tolerance.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
addpath(fullfile(fileparts(here), 'ex2_cart_pole_swing_up', 'try2'));
p = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);

rngWas = rng(17);  restore = onCleanup(@() rng(rngWas));
worstField = 0;  worstEx2 = 0;  scale = 0;
for k = 1:30
    x = [2*randn; pi*randn; 2*randn; 2*randn];
    u = 40*randn;
    gradE = energy_gradient(x, p);

    [F, G] = cartpole_field(x, p);
    xdot = F + G*u;
    assert(all(isfinite(xdot)), 'test_cartpole_physics: non-finite field');
    worstField = max(worstField, abs(gradE.'*xdot - u*x(3)));

    xdot2 = [x(3); x(4); ...
             cart_accel(x(2), x(4), u, p.L, p.m1, p.m2, p.g); ...
             pendulum_accel(x(2), x(4), u, p.L, p.m1, p.m2, p.g)];
    worstEx2 = max(worstEx2, abs(gradE.'*xdot2 - u*x(3)));
    scale = max(scale, abs(u*x(3)));
end
ok = chk(ok, worstField < 1e-10*scale, ...
         sprintf('power balance, cartpole_field: worst |dE/dt - u q1dot| = %.2e (scale %.1f W)', worstField, scale));
ok = chk(ok, worstEx2 < 1e-10*scale, ...
         sprintf('power balance, ex2 helpers:    worst |dE/dt - u q1dot| = %.2e', worstEx2));

%% 3-4. Equilibria: hanging stable, upright unstable
w2 = (p.m1 + p.m2)*p.g/(p.L*p.m1);
A0 = cartpole_state_jac([0; 0; 0; 0], 0, p);
ev0 = eig(A0);
osc = ev0(abs(imag(ev0)) > 1e-8);
w2meas = NaN;  if ~isempty(osc), w2meas = abs(imag(osc(1)))^2; end
ok = chk(ok, numel(osc) == 2 && max(abs(real(ev0))) < 1e-10 && abs(w2meas - w2) < 1e-10, ...
         sprintf('hanging (q2 = 0) is STABLE: oscillatory mode omega^2 = %.6f (theory %.6f), max Re = %.1e', ...
                 w2meas, w2, max(real(ev0))));
Api = cartpole_state_jac([0; pi; 0; 0], 0, p);
evp = eig(Api);
grow = evp(real(evp) > 1e-8);
ok = chk(ok, numel(grow) == 1 && abs(grow - sqrt(w2)) < 1e-10 && max(abs(imag(evp))) < 1e-10, ...
         sprintf('upright (q2 = pi) is UNSTABLE: growth rate %.6f (theory %.6f)', max(real(evp)), sqrt(w2)));

%% 5. Unforced flight conserves energy
oo = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[~, X] = ode113(@(t, x) cartpole_field(x, p), [0 3], [0; 0.7; 0.3; -0.2], oo);
Ev = zeros(size(X, 1), 1);
for k = 1:size(X, 1), Ev(k) = energy(X(k,:).', p); end
ok = chk(ok, max(abs(Ev - Ev(1))) < 1e-9, ...
         sprintf('unforced flight conserves energy: drift %.2e J over 3 s', max(abs(Ev - Ev(1)))));

if ok, fprintf('TEST_CARTPOLE_PHYSICS: ALL PASS\n');
else,  fprintf('TEST_CARTPOLE_PHYSICS: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function E = energy(x, p)
%% Purpose:
%
%   Total mechanical energy FROM THE GEOMETRY ALONE: cart on a rail, bob at
%   (q1 + L sin q2, -L cos q2). Complex-step safe.
%
vbx = x(3) + p.L*cos(x(2))*x(4);
vby = p.L*sin(x(2))*x(4);
E = 0.5*p.m1*x(3)^2 + 0.5*p.m2*(vbx^2 + vby^2) - p.m2*p.g*p.L*cos(x(2));
end

function gE = energy_gradient(x, p)
%% Purpose:
%
%   grad E by complex step through the geometric energy.
%
gE = zeros(4, 1);
h = 1e-20;
for col = 1:4
    xp = x;  xp(col) = xp(col) + 1i*h;
    gE(col) = imag(energy(xp, p))/h;
end
end

function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
