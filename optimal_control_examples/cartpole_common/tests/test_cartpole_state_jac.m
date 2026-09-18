function ok = test_cartpole_state_jac()
%% Purpose:
%
%   The generated state Jacobian must be the derivative it claims to be, and
%   must be usable INSIDE a complex step (that is the whole reason it is
%   generated rather than complex-stepped):
%     1. it matches a complex-step derivative of cartpole_field to 1e-12;
%     2. it matches central finite differences to 1e-7 (an independent
%        method, in case the complex step and the symbolic derivation share
%        a mistake in the field);
%     3. it is real and finite for real inputs, and complex-step clean: a
%        complex-perturbed x gives a finite A with no NaN.
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
p = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);

rngWas = rng(3);  restore = onCleanup(@() rng(rngWas));
worstCS = 0;  worstFD = 0;
for k = 1:15
    x = [randn; 2*randn; randn; randn];
    u = 30*randn;
    A = cartpole_state_jac(x, u, p);

    Acs = zeros(4);
    hcs = 1e-20;
    for col = 1:4
        xp = x;  xp(col) = xp(col) + 1i*hcs;
        [Fp, Gp] = cartpole_field(xp, p);
        Acs(:,col) = imag(Fp + Gp*u)/hcs;
    end
    worstCS = max(worstCS, max(abs(A(:) - Acs(:))));

    Afd = zeros(4);
    hfd = 1e-6;
    for col = 1:4
        xp = x;  xp(col) = xp(col) + hfd;
        xm = x;  xm(col) = xm(col) - hfd;
        [Fp, Gp] = cartpole_field(xp, p);
        [Fm, Gm] = cartpole_field(xm, p);
        Afd(:,col) = ((Fp + Gp*u) - (Fm + Gm*u))/(2*hfd);
    end
    worstFD = max(worstFD, max(abs(A(:) - Afd(:))));
end
ok = chk(ok, worstCS < 1e-12, sprintf('matches a complex step of the field (worst %.1e)', worstCS));
ok = chk(ok, worstFD < 1e-7,  sprintf('matches central differences (worst %.1e)', worstFD));

x = [0.2; 1.0; -0.3; 0.6];
A = cartpole_state_jac(x, 12, p);
ok = chk(ok, isreal(A) && all(isfinite(A(:))), 'real and finite for real inputs');
Ac = cartpole_state_jac(x + [0; 1i*1e-20; 0; 0], 12, p);
ok = chk(ok, all(isfinite(Ac(:))), 'survives a complex-perturbed state (no nested-step damage)');

if ok, fprintf('TEST_CARTPOLE_STATE_JAC: ALL PASS\n');
else,  fprintf('TEST_CARTPOLE_STATE_JAC: FAILURE (see lines above)\n');
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
