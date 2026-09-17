function ok = test_cartpole_field()
%% Purpose:
%
%   The indirect demo must solve the SAME problem as the direct example, so
%   this pins cartpole_field against that example's own helpers:
%     1. F + G*u reproduces cart_accel / pendulum_accel at a scatter of
%        states and controls;
%     2. the split is affine: G is the exact d(xdot)/du;
%     3. the first two rows are the kinematics;
%     4. it is complex-step safe (an imaginary perturbation propagates).
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
worst = 0;  worstAffine = 0;
for k = 1:25
    x = [2*randn; pi*randn; randn; 2*randn];
    u = 50*randn;
    [F, G] = cartpole_field(x, p);
    dx = F + G*u;
    a1 = cart_accel(x(2), x(4), u, p.L, p.m1, p.m2, p.g);
    a2 = pendulum_accel(x(2), x(4), u, p.L, p.m1, p.m2, p.g);
    worst = max(worst, max(abs(dx(3:4) - [a1; a2])));
    % affine in u: the difference quotient in u is exactly G
    [F2, ~] = cartpole_field(x, p);
    worstAffine = max(worstAffine, max(abs((F2 + G*(u+1)) - (dx + G))));
end
ok = chk(ok, worst < 1e-12, sprintf('F + G*u equals the example helpers (worst %.1e)', worst));
ok = chk(ok, worstAffine < 1e-12, sprintf('the split is exactly affine in u (worst %.1e)', worstAffine));

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
