function ok = test_cartpole_field()
%% Purpose:
%
%   The indirect demo must solve the SAME problem as the direct example, so
%   this pins cartpole_field against that example's own helpers:
%     1. F + G*u reproduces cart_accel / pendulum_accel at a scatter of
%        states and controls;
%     2. G is the exact d(xdot)/du, verified against independent helpers;
%     3. the helpers are truly affine in u (second differences vanish);
%     4. the first two rows are the kinematics;
%     5. it is complex-step safe (an imaginary perturbation propagates).
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
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
addpath(fullfile(fileparts(here), 'ex2_cart_pole_swing_up', 'try2'));
p = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);

rngWas = rng(11);  restore = onCleanup(@() rng(rngWas));
worst = 0;  worstGRef = 0;  worstAffine2 = 0;
for k = 1:25
    x = [2*randn; pi*randn; randn; 2*randn];
    u = 50*randn;
    [F, G] = cartpole_field(x, p);
    dx = F + G*u;
    a1 = cart_accel(x(2), x(4), u, p.L, p.m1, p.m2, p.g);
    a2 = pendulum_accel(x(2), x(4), u, p.L, p.m1, p.m2, p.g);
    worst = max(worst, max(abs(dx(3:4) - [a1; a2])));

    % G is the true control derivative (independent reference)
    a1p = cart_accel(x(2), x(4), u+1, p.L, p.m1, p.m2, p.g);
    a1m = cart_accel(x(2), x(4), u,   p.L, p.m1, p.m2, p.g);
    a2p = pendulum_accel(x(2), x(4), u+1, p.L, p.m1, p.m2, p.g);
    a2m = pendulum_accel(x(2), x(4), u,   p.L, p.m1, p.m2, p.g);
    Gref = [0; 0; a1p - a1m; a2p - a2m];
    worstGRef = max(worstGRef, max(abs(G - Gref)));

    % Helpers are truly affine: second difference should vanish
    a1pp = cart_accel(x(2), x(4), u+1, p.L, p.m1, p.m2, p.g);
    a1zm = cart_accel(x(2), x(4), u,   p.L, p.m1, p.m2, p.g);
    a1mm = cart_accel(x(2), x(4), u-1, p.L, p.m1, p.m2, p.g);
    a2pp = pendulum_accel(x(2), x(4), u+1, p.L, p.m1, p.m2, p.g);
    a2zm = pendulum_accel(x(2), x(4), u,   p.L, p.m1, p.m2, p.g);
    a2mm = pendulum_accel(x(2), x(4), u-1, p.L, p.m1, p.m2, p.g);
    worstAffine2 = max(worstAffine2, max(abs([a1pp - 2*a1zm + a1mm; ...
                                              a2pp - 2*a2zm + a2mm])));
end
ok = chk(ok, worst < 1e-12, sprintf('F + G*u equals the example helpers (worst %.1e)', worst));
ok = chk(ok, worstGRef < 1e-12, sprintf('G is the exact d(xdot)/du from independent helpers (worst %.1e)', worstGRef));
ok = chk(ok, worstAffine2 < 1e-12, sprintf('the helpers are truly affine in u: second difference (worst %.1e)', worstAffine2));

x = [0.3; 1.1; -0.7; 0.4];
[F, G] = cartpole_field(x, p);
ok = chk(ok, isequal(F(1:2), x(3:4)) && isequal(G(1:2), [0;0]), ...
         'rows 1-2 are the kinematics, with no control');

h = 1e-20;
[Fc, Gc] = cartpole_field(x + [0; 1i*h; 0; 0], p);
ok = chk(ok, any(imag(Fc) ~= 0) && any(imag(Gc) ~= 0) && all(isfinite([Fc; Gc])), ...
         'complex-step safe: an imaginary perturbation propagates through both outputs');

if ok, fprintf('TEST_CARTPOLE_FIELD: ALL PASS\n');
else,  fprintf('TEST_CARTPOLE_FIELD: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
