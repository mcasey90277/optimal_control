function ok = test_cartpole_pmp_rhs()
%% Purpose:
%
%   The PMP field must BE the Pontryagin conditions, not merely resemble
%   them:
%     1. the control is the stationary point: 2u + lam'G = 0 exactly;
%     2. it is the MINIMISER, not just a stationary point: H at u* is below
%        H at u* +/- d for a scatter of d (Legendre, d2H/du2 = 2 > 0);
%     3. the state rows are F + G u*;
%     4. the costate rows are -dH/dx, checked against a complex step of H
%        taken with u HELD FIXED (which is legitimate: dH/du = 0 at u*);
%     5. zero costates give zero costate rates and a drifting state (the
%        uncontrolled flow), so lam = 0 is the uncontrolled solution.
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

rngWas = rng(5);  restore = onCleanup(@() rng(rngWas));
worstStat = 0;  worstState = 0;  worstCostate = 0;  minGap = inf;
for k = 1:15
    y = [randn; 2*randn; randn; randn; 3*randn; 3*randn; 3*randn; 3*randn];
    x = y(1:4);  lam = y(5:8);
    [dy, u] = cartpole_pmp_rhs(y, p);
    [F, G] = cartpole_field(x, p);

    worstStat  = max(worstStat, abs(2*u + lam.'*G));
    worstState = max(worstState, max(abs(dy(1:4) - (F + G*u))));

    Hstar = u^2 + lam.'*(F + G*u);
    for d = [-1 -0.1 0.1 1]
        Hd = (u+d)^2 + lam.'*(F + G*(u+d));
        minGap = min(minGap, Hd - Hstar);
    end

    % -dH/dx by complex step, u held fixed
    dHdx = zeros(4,1);
    hcs = 1e-20;
    for col = 1:4
        xp = x;  xp(col) = xp(col) + 1i*hcs;
        [Fp, Gp] = cartpole_field(xp, p);
        Hp = u^2 + lam.'*(Fp + Gp*u);
        dHdx(col) = imag(Hp)/hcs;
    end
    worstCostate = max(worstCostate, max(abs(dy(5:8) + dHdx)));
end
ok = chk(ok, worstStat < 1e-12, sprintf('u is stationary: max |2u + lam''G| = %.1e', worstStat));
ok = chk(ok, minGap > 0, sprintf('and a MINIMISER: worst H(u*+d) - H(u*) = %.3e > 0', minGap));
ok = chk(ok, worstState < 1e-12, sprintf('state rows are F + G u* (worst %.1e)', worstState));
ok = chk(ok, worstCostate < 1e-9, sprintf('costate rows are -dH/dx (worst %.1e)', worstCostate));

y0 = [0.1; 0.2; 0; 0; 0; 0; 0; 0];
[dy0, u0] = cartpole_pmp_rhs(y0, p);
ok = chk(ok, u0 == 0 && all(dy0(5:8) == 0) && any(dy0(1:4) ~= 0), ...
         'lam = 0 gives u = 0, still costates, and the uncontrolled drift');

if ok, fprintf('TEST_CARTPOLE_PMP_RHS: ALL PASS\n');
else,  fprintf('TEST_CARTPOLE_PMP_RHS: FAILURE (see lines above)\n');
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
