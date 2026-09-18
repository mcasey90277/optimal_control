%% CARTPOLE_MINENERGY_STUDY  Minimum-effort cart-pole swing-up, step by step.
%
%  Edit the parameter blocks, press Run. The STEPS are done here rather than
%  behind the front door (run_cartpole_pmp), so the machinery is visible.
%  Two computations stay in the library because reimplementing them here
%  would create a second, unverified copy: the physics ORACLE (V1,
%  test_cartpole_physics -- power balance from geometry, not from the
%  equations under test) and the Jacobi determinant sweep (S2,
%  oc.ms_conjugate_test). Section 8 asserts the inline numbers against that
%  instrument for the same reason.
%
%  Prerequisite: optimal_control_examples/cartpole_common and oclib on the
%  path (this script adds both), and data/cartpole_direct_ref.mat present --
%  the committed direct-collocation solve that seeds the shooting. Nothing
%  here touches orbit_transfer.
%
%  Sections:
%    0  tolerances          one struct, feeding every printed line AND verdict
%    1  the plant           constants, and the physics ORACLE run here
%    2  the problem         cost, boundary conditions, the control law
%    3  the seed            direct multipliers -> costates, sign resolved
%    4  solve               oc.ms_bvp at fixed t_f, STMs kept
%    5  independent check   a perturbed seed must find the same root
%    6  NECESSARY  (N1-N6)  Pontryagin, one condition at a time
%    7  SUFFICIENCY (S1-S2) strict Legendre, conjugate point
%    8  assert vs library   the inline numbers against oc.ms_conjugate_test
%    9  plot                states, control, switching-free sigma, H(t)
%
%  Diagnostic IDs:  N necessary | S sufficiency | V validity | X cross-check.
%  A failed N or X gate means the root or the implementation is broken. A
%  failed or unresolved S gate is a FINDING about this trajectory: reported,
%  and thrown only when selfCheck.strict is set.
%
%  ONE FLIGHT. Section 6 flies the converged costates once, into (t, X, Lam,
%  U, H); every number in sections 6-9 and the figure come from that flight,
%  so the reader never has to ask which trajectory owns a reported number.
%
%  M. Casey                                                   (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------------------------------------------------------

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'cartpole_common'), ...
        fullfile(fileparts(here), 'cartpole_common', 'tests'), ...
        fullfile(fileparts(fileparts(here)), 'oclib'));

%% 0. Tolerances -- every gate below reads this struct and nothing else:
tol = struct( ...
    'R',        1e-10, ...   % ms_bvp residual. ABSOLUTE, and the costates are
    ...                      % O(1e3), so this is ~1e-13 relative; the
    ...                      % propagator's RelTol floors the achievable
    ...                      % residual near 1e-11 regardless of iterations.
    ...                      % MEASURED 2026-09-17 on this root: 9.607e-12,
    ...                      % identical to run_cartpole_pmp's own solve.
    'H',        1e-9,  ...   % spread of H along the arc, relative. MEASURED
    ...                      % 3.381e-12 on the 2001-point re-flight at RelTol
    ...                      % 2.5e-14 -- this is integration quality, not
    ...                      % optimality.
    'miss',     1e-9,  ...   % terminal miss of the re-flown arc. MEASURED
    ...                      % 1.475e-10 at the RelTol 2.5e-14 this script
    ...                      % actually flies. (The 3.67e-11 in older notes
    ...                      % came from a probe at RelTol 1e-14 and was never
    ...                      % re-taken at the shipped setting.)
    'flown',    1e-6,  ...   % miss when flying the CLOSED-FORM control.
    ...                      % MEASURED 1.972e-07.
    'gap',      -1e-10, ...  % minimum-principle gap floor: H(u*+d) - H(u*)
    ...                      % must exceed this (i.e. be non-negative to
    ...                      % round-off) for every probe d. MEASURED worst
    ...                      % gap +1.000e-02, which is exactly d^2 at the
    ...                      % smallest probe (d = 0.1): H is quadratic in u
    ...                      % with H_uu = 2, so the gap IS d^2 and the floor
    ...                      % is never approached on this problem.
    'adj',      1e-9,  ...   % adjoint equations, relative, vs complex-step
    ...                      % differentiation of H. MEASURED 6.091e-16.
    'agree',    1e-6,  ...   % inline numbers vs the library instrument.
    'xJ',       2e-3,  ...   % direct vs indirect cost, relative. MEASURED
    ...                      % 6.214e-04 -- this is the DISCRETISATION gap
    ...                      % between the two methods, not solver error, so
    ...                      % it does not shrink with tolerances.
    'xU',       2e-2,  ...   % direct vs indirect control, RMS difference over
    ...                      % RMS level. MEASURED 3.752e-03, same origin as
    ...                      % xJ (the 201-node trapezoid rule), so the gate
    ...                      % is set an order above the measurement rather
    ...                      % than at solver precision.
    'sigRatio', 1e-10);      % conjugate test: below this the determinant's
                             % sign is not trusted and the sample is reported
                             % unresolved rather than counted.

% THE SCRIPT'S SELF-CHECK. A study script that prints FAIL and exits 0 is a
% report, not a check. `strict` decides what an unsuccessful study DOES:
%   false (default)  a failed NECESSARY, VALIDITY or CROSS-CHECK gate throws
%                    -- those mean the root or the implementation is broken.
%                    A failed or unresolved SUFFICIENCY gate is a legitimate
%                    finding about this trajectory and is reported, not thrown.
%   true             any gate short of a full claim throws. Use it in a
%                    regression harness, where "unresolved" is a regression.
selfCheck = struct('strict', false);

%% 1. The plant, and the ORACLE that says the physics is right (V1):
p = cartpole_params();
fprintf('\n=== 1. plant: m1 %g kg, m2 %g kg, L %g m, g %g m/s^2\n', p.m1, p.m2, p.L, p.g);
V1 = test_cartpole_physics();     % power balance from GEOMETRIC energy, plus
                                  % equilibrium character and unforced
                                  % conservation -- an oracle OUTSIDE the
                                  % equations it checks. A sign error lived in
                                  % this plant until 2026-09-17 because every
                                  % other test compared the code with a copy
                                  % of itself.
fprintf('  V1 physics oracle                                       %s\n', pf(V1));
assert(V1, 'study:physics', 'the plant fails its own oracle: nothing below means anything');

%% 2. The problem:
R  = load(fullfile(here, 'data', 'cartpole_direct_ref.mat'));
tf = R.tf;
x0 = [0; 0; 0; 0];
xf = [0; pi; 0; 0];
fprintf('\n=== 2. minimise J = int_0^%g u^2 dt,  x(0) = [0 0 0 0],  x(t_f) = [0 pi 0 0]\n', tf);
fprintf('  control UNBOUNDED, so H = u^2 + lam''(F + G u) is minimised at its\n');
fprintf('  stationary point: u* = -lam''G/2, with H_uu = 2 > 0 exactly.\n');
fprintf('  t_f FIXED, so H is constant along the arc but NOT zero, and there\n');
fprintf('  is no transversality condition to check (N4 does not apply).\n');

%% 3. The seed: the direct solve's multipliers become costates
%  duals_to_costates resolves the GLOBAL SIGN by a primer vote over 3-vector
%  thrust directions. This problem's control is a scalar force, so there is
%  no such vote: the mapping is asked for the costates unsigned, and the sign
%  is fixed here by the same principle in scalar form -- the implied control
%  u = -lam'G/2 must agree in sign with the control the direct solve used.
[lamS, tS] = oc.duals_to_costates(struct('scheme', 'trapezoid', 'mu', R.muDefect, ...
                                         'tNodes', R.tN, 'velRows', 3:4));
tS       = tS(:).';
uImplied = zeros(1, numel(tS));
Xs       = interp1(R.tN, R.X.', tS, 'pchip').';
for k = 1:numel(tS)
    [~, Gk] = cartpole_field(Xs(:,k), p);
    uImplied(k) = -(lamS(:,k).'*Gk)/2;
end
uDirect  = reshape(interp1(R.tN, R.U, tS, 'pchip'), 1, []);
signCorr = sum(uImplied.*uDirect) / sqrt(sum(uImplied.^2)*sum(uDirect.^2));
if signCorr < 0, lamS = -lamS;  uImplied = -uImplied; end
ampRatio = sqrt(sum(uImplied.^2)/sum(uDirect.^2));
fprintf('\n=== 3. seed from the direct multipliers\n');
fprintf('  sign correlation %+.4f (flipped: %d), amplitude ratio %.4f\n', ...
        signCorr, signCorr < 0, ampRatio);
fprintf('  the correlation votes the SIGN only; the amplitude ratio is printed\n');
fprintf('  beside it because a wrong costate SCALE is invisible to that vote.\n');

K     = 8;
tGrid = linspace(0, tf, K+1);
Xg    = interp1(R.tN, R.X.', tGrid, 'pchip').';
Lg    = interp1(tS, lamS.', tGrid, 'pchip', 'extrap').';
seed  = struct('tf', tf, 'tGrid', tGrid, 'Y', [Xg; Lg]);
seed.Y(1:4,1) = x0;

%% 4. Solve: four unknowns lam(0) against four terminal conditions
prob = struct('ny', 8, 'freeIdx0', 5:8, ...
    'prop', @(dt, y0, needSTM) cartpole_pmp_prop(dt, y0, needSTM, p), ...
    'rhs',  @(y) cartpole_pmp_rhs(y, p), ...
    'terminal', @(y, needJ) studyTerminal(y, xf, needJ));
% keepSTMs: section 7 needs info.PHI. Without it the conjugate test has
% nothing to read.
[~, it] = oc.ms_bvp(prob, seed, struct('fixedTf', true, 'tolR', tol.R, ...
                                       'maxIter', 60, 'polishMax', 5, ...
                                       'keepSTMs', true));
lam0 = it.Y(5:8, 1);
fprintf('\n=== 4. multiple shooting: %d segments, converged %d, |R| = %.3e\n', ...
        K, it.converged, it.normR);
fprintf('  lam(0) = [%.6f %.6f %.6f %.6f]\n', lam0);

%% 5. Independent checks
%  X2 first, because it is about the solver and costs one more solve: a
%  different seed must find the SAME root. This guards the SOLVER (one basin,
%  not two); it cannot guard the physics, which is V1's job. X1 -- direct
%  against indirect -- needs the flight, so it is printed in section 6.
seed2 = seed;  seed2.Y(5:8,:) = 1.05*seed.Y(5:8,:);
[~, it2] = oc.ms_bvp(prob, seed2, struct('fixedTf', true, 'tolR', tol.R, ...
                                         'maxIter', 60, 'polishMax', 5));
dLam0 = max(abs(it2.Y(5:8,1) - lam0));
X2    = it2.converged && dLam0 < 1e-6;
fprintf('\n=== 5. perturbed seed (costates x 1.05)\n');
fprintf('  X2 same root: max|dlam(0)| = %.3e / 1.0e-06                %s\n', dLam0, pf(X2));

%% 6. NECESSARY conditions
%  ONE flight, read by sections 6 through 9. RelTol 2.5e-14, not the
%  propagator's 1e-12: this single shot re-flies the whole horizon from
%  lam0, so its own integration error lands directly in the terminal miss
%  below. Measured on this root: miss 6.88e-09 at RelTol 1e-12, 8.20e-10 at
%  1e-13, and 1.475e-10 at the 2.5e-14 used here, while the engine's own
%  last-arc residual is 7.268e-14 -- the loose figure was the MEASUREMENT,
%  not the solution. (The first two figures are inherited from the front
%  door's own probe; 1.475e-10 and 7.268e-14 were re-measured here.)
nPts   = 2001;
[t, Y] = ode113(@(tt, y) cartpole_pmp_rhs(y, p), linspace(0, tf, nPts), ...
                [x0; lam0], odeset('RelTol', 2.5e-14, 'AbsTol', 1e-16));
X   = Y(:,1:4).';
Lam = Y(:,5:8).';
U   = zeros(1, nPts);
H   = zeros(1, nPts);
sig = zeros(1, nPts);          % lam'G -- no switching structure here, but it
for k = 1:nPts                 % is the same object the bounded problems switch
    [Fk, Gk] = cartpole_field(X(:,k), p);
    sig(k)   = Lam(:,k).'*Gk;
    U(k)     = -sig(k)/2;
    H(k)     = U(k)^2 + Lam(:,k).'*(Fk + Gk*U(k));
end
J = trapz(t, U.^2);

fprintf('\n=== 6. NECESSARY\n');
N1 = it.converged && it.normR < tol.R;
fprintf('  N1 BVP residual         %.3e / %.1e                    %s\n', it.normR, tol.R, pf(N1));

Hrel = (max(H) - min(H))/max(abs(H));
N2   = Hrel < tol.H;
fprintf('  N2 H constant (rel)     %.3e / %.1e                    %s\n', Hrel, tol.H, pf(N2));
fprintf('     H = %.6f, constant but NOT zero: t_f is fixed.\n', mean(H));

missTerminal = max(abs(X(:,end) - xf));
N3 = missTerminal < tol.miss;
fprintf('  N3 terminal miss        %.3e / %.1e                    %s\n', missTerminal, tol.miss, pf(N3));

fprintf('  N4 transversality       n/a -- t_f and x(t_f) both fixed\n');

%  N5: the adjoint equation, lamdot = -dH/dx at fixed u, checked against
%  COMPLEX-STEP differentiation of H -- code that does not share a line with
%  the analytic Jacobian inside cartpole_pmp_rhs.
adjWorst = 0;
for k = round(linspace(1, nPts, 41))
    dy   = cartpole_pmp_rhs([X(:,k); Lam(:,k)], p);
    gH   = zeros(4,1);
    for m = 1:4
        xc     = complex(X(:,k));
        xc(m)  = xc(m) + 1e-30i;
        [Fc, Gc] = cartpole_field(xc, p);
        gH(m)  = imag(Lam(:,k).'*(Fc + Gc*U(k)))/1e-30;
    end
    adjWorst = max(adjWorst, max(abs(dy(5:8) + gH))/max(1, max(abs(gH))));
end
N5 = adjWorst < tol.adj;
fprintf('  N5 adjoint equations    %.3e / %.1e (rel)              %s\n', adjWorst, tol.adj, pf(N5));

%  N6: the minimum principle itself -- H(u* + d) - H(u*) >= 0 for every
%  probe d, at a scatter of times. Not a derivative check: a probe at finite
%  d would catch a stationary point that is a MAXIMUM, which dH/du = 0 would
%  not.
gapWorst = inf;
for k = round(linspace(1, nPts, 41))
    [Fk, Gk] = cartpole_field(X(:,k), p);
    for d = [-10 -1 -0.1 0.1 1 10]
        uP  = U(k) + d;
        HP  = uP^2 + Lam(:,k).'*(Fk + Gk*uP);
        gapWorst = min(gapWorst, HP - H(k));
    end
end
N6 = gapWorst > tol.gap;
fprintf('  N6 min-principle gap    %+.3e / %+.1e (worst)         %s\n', gapWorst, tol.gap, pf(N6));

%  The flown-control check: the CLOSED-FORM control on the flown state, not
%  an interpolant of U (interpolating u put the resampling error straight
%  into the miss: 8.9e-07 against a 1e-06 gate).
ppL  = pchip(t, Lam);
zEnd = oc.fly_control(x0, [0 tf], ...
    @(tt, x) studyFlyRhs(x, studyU(ppL, min(max(tt, t(1)), t(end)), x, p), p), ...
    struct('mode', 'span', 'solver', @ode113, 'RelTol', 1e-13, 'AbsTol', 1e-15));
missFlown = max(abs(zEnd - xf));
V2a = missFlown < tol.flown;
fprintf('  V2 flown-control miss   %.3e / %.1e                    %s\n', missFlown, tol.flown, pf(V2a));

%  X1: the DIRECT solve and this INDIRECT one are two different methods on
%  one problem. They must agree on the cost and on the control history. They
%  will NOT agree to solver tolerance -- the direct solve discretises then
%  optimises, so its cost carries the trapezoid rule's own error; the gate is
%  set above the measured agreement (J to 0.062%, control RMS to 0.38% of
%  the RMS level), not at machine precision.
dJrel = abs(J - R.J)/R.J;
uDrms = sqrt(trapz(t, (U - reshape(interp1(R.tN, R.U, t, 'pchip'), 1, [])).^2)/tf);
uRms  = sqrt(trapz(t, U.^2)/tf);
X1    = (dJrel < tol.xJ) && (uDrms/uRms < tol.xU);
fprintf('  X1 direct vs indirect   dJ %.3e / %.1e,  du/u %.3e / %.1e   %s\n', ...
        dJrel, tol.xJ, uDrms/uRms, tol.xU, pf(X1));
fprintf('     J = %.6f (indirect) vs %.6f (direct, 201-node trapezoid)\n', J, R.J);

%% 7. SUFFICIENCY
fprintf('\n=== 7. SUFFICIENCY\n');
S1 = true;      % H_uu = d2/du2 (u^2 + lam'(F + Gu)) = 2, EXACTLY, everywhere
fprintf('  S1 strict Legendre      H_uu = 2 > 0 (exact, not measured)    %s\n', pf(S1));

%  S2: the Jacobi / conjugate-point test. quotientDir [] and freeTime false:
%  t_f is FIXED, and the running cost u^2 breaks the costate-scaling
%  invariance the orbit min-time problems quotient out -- the same
%  convention costate_common/tests/test_conj_fixedtf.m uses. stateRows 1:4
%  is ALL state rows, as the fixed-t_f rule in the instrument's header
%  requires (the interior Jacobi field must vanish in the full state).
conjOut = oc.ms_conjugate_test(it, struct('stateRows', 1:4, 'costateCols', 5:8, ...
                                       'quotientDir', [], 'freeTime', false, ...
                                       'resolvedTol', tol.sigRatio));
S2 = (conjOut.nInterior == 0) && conjOut.covered && (conjOut.nUnresolved == 0);
fprintf('  S2 conjugate point      interior %d, touch %d, unresolved %d, covered %d  %s\n', ...
        conjOut.nInterior, conjOut.nTouch, conjOut.nUnresolved, conjOut.covered, pf(S2));
fprintf('     instrument verdict %s: %s\n', conjOut.verdict, conjOut.reason);
if ~S2
    fprintf('     A conjugate point, or an unresolved sample, is a FINDING about\n');
    fprintf('     THIS trajectory -- not a tolerance to loosen. Report it.\n');
end

%% 8. ASSERT the inline numbers against the library instrument
fprintf('\n=== 8. inline vs library\n');
%  What the instrument computed, re-derived here where it is cheap: the
%  determinant it reports at the FINAL sample must match a determinant built
%  from the same STM by this script's own two lines. The ORDER is the
%  instrument's own -- it accumulates PhiCum = PHI{k}*PhiCum from k = 1, so
%  Phi(t_f, 0) = PHI{K}*...*PHI{1}, left-multiplying as k advances.
PHIend = it.PHI{end};
for k = numel(it.PHI)-1:-1:1, PHIend = PHIend*it.PHI{k}; end
Mend    = PHIend(1:4, 5:8);
detMine = det(Mend);
dMine   = sign(detMine);
dLib    = sign(conjOut.detScaled(end));
%  MAGNITUDE TOO, not the sign alone. The instrument reports
%  sign(det)*|det|^(1/m), so |det(Mend)|^(1/4) is directly comparable -- and
%  the comparison is what gives this gate teeth. MEASURED: a deliberately
%  REVERSED product (PHI{1}*...*PHI{K} instead of PHI{K}*...*PHI{1}) still
%  came out sign +1 on this trajectory, so a sign-only check would have
%  passed the very error this gate exists to catch. Its |det|^(1/4) is
%  2.246689e-01 against the instrument's 1.091147e-01 (106% off), so the
%  magnitude comparison fires on it. Agreement on the correct ordering:
%  2.72e-13 relative.
magMine = abs(detMine)^(1/4);
magLib  = abs(conjOut.detScaled(end));
dMag    = abs(magMine - magLib)/max(magLib, realmin);
V2b     = isequal(dMine, dLib) && dMine ~= 0 && dMag < tol.agree;
fprintf('  V2 det at t_f           sign %+d vs %+d, |det|^(1/4) %.6e vs %.6e (rel %.2e / %.1e)  %s\n', ...
        dMine, dLib, magMine, magLib, dMag, tol.agree, pf(V2b));

%  Where re-deriving would mean a second unverified copy of delicate
%  machinery (the sweep, the bracketing, the equilibration), the instrument's
%  own bookkeeping is checked instead -- shapes, the coverage claim against
%  the time it says it sampled through, and nUnresolved re-counted from the
%  sigma ratios it reports. Note .tested is a LOGICAL ("something was
%  testable"), not a sample count: the live samples are counted here.
nSamp = numel(conjOut.t);
live  = conjOut.firstFullRank:nSamp;
nLive = numel(live);
V2c   = islogical(conjOut.tested) && conjOut.tested && islogical(conjOut.covered) && ...
        nSamp == numel(conjOut.detScaled) && nSamp == numel(conjOut.sigRatio) && ...
        isequal(conjOut.covered, abs(conjOut.sampledThrough - tf) <= tol.agree) && ...
        isequal(conjOut.pass, strcmp(conjOut.verdict, 'PASS')) && ...
        conjOut.nUnresolved == nnz(conjOut.sigRatio(live) <= tol.sigRatio | ...
                                   conjOut.detScaled(live) == 0);
fprintf('  V2 instrument self-consistent (%d samples, %d live from t = %.4f s)  %s\n', ...
        nSamp, nLive, conjOut.tFirstFullRank, pf(V2c));
V2 = V2a && V2b && V2c;

%% 9. Plot
figure('color', [1 1 1]);
subplot(4,1,1); plot(t, X(1,:), 'b', t, X(2,:), 'k'); grid on
ylabel('states'); legend({'q_1 (m)', 'q_2 (rad)'}, 'Location', 'best');
title(sprintf('cart-pole minimum effort: J = %.6f, t_f = %g s', J, tf));
subplot(4,1,2); plot(t, U, 'k'); grid on; ylabel('u (N)');
subplot(4,1,3); plot(t, sig, 'k'); grid on; ylabel('\sigma = \lambda''G');
subplot(4,1,4); plot(t, H - mean(H), 'k'); grid on
xlabel('t (s)'); ylabel('H - mean(H)');

%% Verdict
verdict = struct('necessary', all([N1 N2 N3 N5 N6 X1 X2 V1 V2]), ...
                 'sufficiency', S1 && S2, 'claim', false, 'why', '');
verdict.claim = verdict.necessary && verdict.sufficiency;
if ~verdict.necessary
    verdict.why = 'a necessary or validity gate failed: the root or the implementation is wrong';
elseif ~verdict.sufficiency
    verdict.why = 'necessary conditions hold; the second-order test did not resolve or found a conjugate point';
else
    verdict.why = ['consistent with a strict strong local minimum for FIXED endpoints ' ...
                   'and FIXED t_f; numerical and sampled, not a certificate'];
end
fprintf('\n=== VERDICT: necessary %d, sufficiency %d, claim %d\n', ...
        verdict.necessary, verdict.sufficiency, verdict.claim);
fprintf('  %s\n', verdict.why);
if selfCheck.strict
    assert(verdict.necessary, 'study:necessary', '%s', verdict.why);
    assert(verdict.sufficiency, 'study:sufficiency', '%s', verdict.why);
else
    assert(verdict.necessary, 'study:necessary', '%s', verdict.why);
end

% ------------------------------------------------------------------------
function s = pf(c)
%% Purpose:
%
%   'PASS' if c, else 'FAIL'.
%
if c, s = 'PASS'; else, s = 'FAIL'; end
end

function [g, dgdy] = studyTerminal(y, xf, needJ)
%% Purpose:
%
%   The four terminal conditions x(t_f) = xf, and their Jacobian.
%
g = y(1:4) - xf;
dgdy = [];
if needJ, dgdy = [eye(4), zeros(4)]; end
end

function u = studyU(ppL, tt, x, p)
%% Purpose:
%
%   The PMP control at time tt: the costate from the converged arc, G from
%   the FLOWN state, u = -lam'G/2. No interpolation of u itself.
%
[~, G] = cartpole_field(x, p);
u = -(ppval(ppL, tt).'*G)/2;
end

function dx = studyFlyRhs(x, u, p)
%% Purpose:
%
%   The TRUE dynamics with a supplied control -- what the flown check flies.
%
[F, G] = cartpole_field(x, p);
dx = F + G*u;
end
