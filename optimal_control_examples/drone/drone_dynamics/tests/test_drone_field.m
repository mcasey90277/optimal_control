function ok = test_drone_field()
%% Purpose:
%
%   Shape and structure of drone_field for each model, and complex-step
%   safety: a complex perturbation of one state must produce a finite
%   complex result whose imaginary part is a derivative (compared against
%   central differences of F here; the full Jacobian test is
%   test_drone_state_jac). These are CONSISTENCY checks on the field's
%   form; the physics is pinned by test_drone_physics.
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
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
rngWas = rng(7);  restore = onCleanup(@() rng(rngWas));

models = {'kinematic', 'pointmass', 'rigid6dof'};
for mdl = 1:3
    p = drone_params(models{mdl});
    x = randn(p.nx, 1);  if p.nx == 12, x(7:8) = 0.3*x(7:8); end
    [F, B] = drone_field(x, p);
    ok = chk(ok, isequal(size(F), [p.nx 1]) && isequal(size(B), [p.nx p.nu]), ...
             sprintf('%s: F is %dx1, B is %dx%d', p.model, p.nx, p.nx, p.nu));
    ok = chk(ok, all(isfinite(F(:))) && all(isfinite(B(:))), sprintf('%s: finite', p.model));
    % complex-step safety: imag part of F(x + i h e_k)/h vs central FD
    worst = 0;  hc = 1e-20;  hd = 1e-6;
    for k = 1:p.nx
        xc = x;  xc(k) = xc(k) + 1i*hc;
        Fc = drone_field(xc, p);
        xp = x;  xp(k) = xp(k) + hd;  xm = x;  xm(k) = xm(k) - hd;
        Fd = (drone_field(xp, p) - drone_field(xm, p))/(2*hd);
        worst = max(worst, max(abs(imag(Fc)/hc - Fd)));
        ok = ok && all(isfinite(Fc(:)));
    end
    ok = chk(ok, worst < 1e-7, sprintf('%s: complex step of F matches central differences (%.1e)', p.model, worst));
end

p = drone_params('rigid6dof');
x = randn(12,1);  x(7:8) = 0.3*x(7:8);
[~, B] = drone_field(x, p);
ok = chk(ok, max(max(abs(B(4:6,:) - B(4:6,1)*ones(1,4)))) < 1e-15, ...
         'rigid6dof: every rotor pushes along the same body z axis (B(4:6,:) columns equal)');
ok = chk(ok, isequal(B(1:3,:), zeros(3,4)) && isequal(B(7:9,:), zeros(3,4)), ...
         'rigid6dof: thrusts do not enter the position or Euler-rate rows');
Mix = [0 p.l 0 -p.l; -p.l 0 p.l 0; p.kappa -p.kappa p.kappa -p.kappa];
ok = chk(ok, max(max(abs(B(10:12,:) - diag(1./p.I)*Mix))) < 1e-12, 'rigid6dof: B(10:12,:) = I^-1 Mix (the "+" mixer)');
xT = x;  xT(8) = p.thetaMax;
F = drone_field(xT, p);
ok = chk(ok, all(isfinite(F)), 'rigid6dof: finite at theta = thetaMax');
xS = x;  xS(8) = pi/2;
F = drone_field(xS, p);
ok = chk(ok, ~all(isfinite(F)) || max(abs(F)) > 1e10, 'rigid6dof: blows up at theta = pi/2 (the singularity is not hidden)');

if ok, fprintf('TEST_DRONE_FIELD: ALL PASS\n'); else, fprintf('TEST_DRONE_FIELD: FAILURE\n'); end
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
