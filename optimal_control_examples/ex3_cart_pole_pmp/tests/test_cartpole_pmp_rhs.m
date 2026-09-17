function ok = test_cartpole_pmp_rhs()
%% Purpose:
%
%   The PMP field must BE the Pontryagin conditions, not merely resemble
%   them:
%     1. the control is the stationary point: 2u + lam'G = 0 exactly;
%     2. POINTWISE Hamiltonian minimisation: H at u* is below H at u* +/- d.
%        Honest scope: given check 1 this gap is algebraically d^2, so it adds
%        no independent evidence; it documents the Legendre condition
%        (d2H/du2 = 2 > 0). It says NOTHING about the trajectory being a
%        minimum -- H_uu > 0 gives the unique pointwise minimiser, and the BVP
%        can still have several normal extremals, saddles included;
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
addpath(here, fullfile(fileparts(here), 'cartpole_common'));
p = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);

rngWas = rng(5);  restore = onCleanup(@() rng(rngWas));
worstStat = 0;  worstState = 0;  worstCostate = 0;  minGap = inf;  lamMax = 0;
for k = 1:15
    % half the scatter at O(1) costates, half at the SOLUTION scale (the
    % converged arc has |lam| ~ 1e3): a term that only matters at scale
    % would be invisible to unit-sized costates
    lamScale = 3;  if k > 8, lamScale = 1500; end
    y = [randn; 2*randn; randn; randn; lamScale*randn(4,1)];
    assert(all(isfinite(y)), 'test_cartpole_pmp_rhs: non-finite sample');
    x = y(1:4);  lam = y(5:8);
    [dy, u] = cartpole_pmp_rhs(y, p);
    [F, G] = cartpole_field(x, p);

    assert(all(isfinite(dy)) && isfinite(u), 'test_cartpole_pmp_rhs: non-finite field');
    lamMax     = max(lamMax, max(abs(lam)));
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
ok = chk(ok, worstStat < 1e-12*max(1, lamMax), sprintf('u is stationary: max |2u + lam''G| = %.1e', worstStat));
ok = chk(ok, minGap > 0, sprintf('pointwise H-minimiser (Legendre): worst H(u*+d) - H(u*) = %.3e > 0', minGap));
ok = chk(ok, worstState < 1e-12*max(1, lamMax), sprintf('state rows are F + G u* (worst %.1e)', worstState));
ok = chk(ok, worstCostate < 1e-9*max(1, lamMax), ...
         sprintf('costate rows are -dH/dx (worst %.1e at |lam| up to %.0f)', worstCostate, lamMax));

y0 = [0.1; 0.2; 0; 0; 0; 0; 0; 0];
[dy0, u0] = cartpole_pmp_rhs(y0, p);
F0 = cartpole_field(y0(1:4), p);
ok = chk(ok, u0 == 0 && all(dy0(5:8) == 0) && isequal(dy0(1:4), F0), ...
         'lam = 0 gives u = 0, still costates, and EXACTLY the uncontrolled drift F');

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
