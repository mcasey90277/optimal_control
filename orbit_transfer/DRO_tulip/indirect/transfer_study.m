%% TRANSFER_STUDY  One minimum-time low-thrust transfer: DRO-to-Tulip
%
%   Edit the parameter blocks, press Run. The STEPS are done here rather than
%   behind a front door, so the machinery is visible. Three computations stay
%   in the library because reimplementing them here would create a second,
%   unverified copy: the lift-space rank (S3), the Jacobi determinant (S4) and
%   the H6 margin (V1). Section 7 asserts the inline numbers against that
%   library for the same reason.
%
%   PREREQUISITE. Section 4 needs a converged SEED. The seed library
%   (dro_tulip_library, with the certified 70 mN catalog included) covers ONE
%   operating point -- the tau = 1 DRO to the 7-petal tulip (pm = -1) at
%   70 mN, Isp 900 s, 150 kg -- at its 115 certified phase pairs: the grid
%   sD = k/12, sA = 0.0754 + j/12, less the cells the catalog does not hold.
%   For any other case, walk to it first with run_dro_tulip and study its
%   result here.
%
%     0  tolerances                       (one source, set before anything runs)
%     1  generate the DEPARTURE orbit      (family + parameters -> a real orbit)
%     2  generate the TARGET orbit         (same)
%     3  engine, phases, endpoints         (nondimensionalisation spelled out)
%     4  solve                             (seed -> multiple shooting -> costates)
%     5  independent verification          (a second solver must not move them)
%     6  NECESSARY conditions              (Pontryagin, one at a time)
%     7  SUFFICIENCY hypotheses            (Bonnard-Caillau-Trelat, one at a time)
%     8  interactive 3D plot               (drag to rotate)
%
%   Diagnostic IDs are stable and grouped. N1-N6 are Pontryagin's necessary
%   conditions. S1-S4 are the hypotheses of the BCT sufficiency theorem. V1-V2
%   are the validity of the instruments that evaluated them. X1 is a
%   cross-check of OUR implementation, not a condition of the theory.
%
%   Sections 6 and 7 are separate because the theory separates them: the
%   first-order conditions are HYPOTHESES of the sufficiency theorem, not
%   consequences of it, and the conjugate test is computed ALONG the extremal
%   -- off one, its determinant means nothing. On this problem 12 of 14
%   candidates satisfied every first-order condition and were then refuted.
%
%   Section 5 is not an optimality condition at all. It guards against OUR
%   solver: a bug in our shooting could give a self-consistent answer to the
%   wrong problem, and only a second implementation catches that.
%
%   ONE FLIGHT. Section 4 flies the converged costates once, into `flight`;
%   every number in sections 6-8 and the figure come from that object, so the
%   reader never has to ask which trajectory owns a reported number.
%
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.

%% paths
clear; clc
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));

%% CR3BP constants
muStar = 0.012150585609624;                 % Earth-Moon mass ratio
lStar  = 389703.264829278;                  % km
tStar  = 382981.289129055;                  % s
day    = @(tnd) tnd*tStar/86400;

%% ========================================================================
%  0. TOLERANCES -- one source, used by every printed line AND the verdict.
%     Written twice they drift, and a report then disagrees with itself.
%% ========================================================================
tol = struct( ...
    'closure',   1e-7,  ...   % periodicity of each generated orbit, and of the
    ...                       %   interpolant's VALUE across s = 0 (it cannot be
    ...                       %   more periodic than its data)
    'seamDeriv', 1e-6,  ...   % interpolant DERIVATIVE mismatch across s = 0
    'R',         1e-8,  ...   % N1 BVP residual, inf-norm (solve must reach it)
    'H',         1e-6,  ...   % N2 |H| along the arc (exactly 0 in theory)
    'km',        100,   ...   % N3 flown arrival, position
    'ms',        10,    ...   % N3 flown arrival, velocity
    'lamm',      1e-6,  ...   % N4 transversality lambda_m(t_f) = 0
    'adj',       1e-7,  ...   % N5 adjoint equations, relative
    'gap',       1e-12, ...   % N6 minimum-principle gap of the applied control
    'dz',        1e-6,  ...   % X1 second solver must agree to this
    'h6Margin',  1.0,   ...   % V1 required (c/T)/lambda_m(0)
    'agree',     1e-6);       % V2 inline numbers vs the library instruments

%% ========================================================================
%  1. DEPARTURE ORBIT -- a distant retrograde orbit (DRO)
%     Parameter: the PERIOD tau (ND). Larger tau = larger orbit.
%       tau 0.5 -> 2.22 d      tau 1.0 -> 4.43 d      tau 2.0 -> 8.87 d
%% ========================================================================
dep.family = 'dro';
dep.tau    = 1.0;                            % <-- change the departure orbit here

% Generate DRO
[tD, rvD, infoD] = get_family_orbit(dep.family, struct('tau', dep.tau, 'muStar', muStar));
departure = struct('params', dep, 'timeND', tD, 'stateND', rvD, ...
                   'periodND', infoD.periodND, 'closure', norm(rvD(end,:) - rvD(1,:)));

fprintf('1. DEPARTURE %s: period %.4f ND (%.3f d), %d samples\n', ...
        upper(dep.family), infoD.periodND, day(infoD.periodND), numel(tD));
fprintf('   max radius from the barycentre %.4f ND (%.0f km)\n', ...
        max(vecnorm(rvD(:,1:3),2,2)), max(vecnorm(rvD(:,1:3),2,2))*lStar);
fprintf('   periodicity closure |x(T) - x(0)| = %.2e / %-7.0e %s\n', departure.closure, ...
        tol.closure, pass(departure.closure < tol.closure));
assert(departure.closure < tol.closure, 'the departure orbit is not periodic to tolerance');

%% ========================================================================
%  2. TARGET ORBIT -- a tulip
%     Parameters: the PETAL COUNT Np and the branch pm (+1 is the z-mirror
%     of -1). The period is NOT free: it is locked to 2*pi*(Np-2)/(Np-1).
%       Np 3 -> 13.93 d    Np 5 -> 20.89 d    Np 7 -> 23.21 d    Np 9 -> 24.37 d
%% ========================================================================
arr.family = 'tulip';
arr.Np     = 7;                              % <-- change the target orbit here
arr.pm     = -1;                             %     branch: -1 or +1

% Generate tulip
[tT, rvT, infoT] = get_family_orbit(arr.family, struct('Np', arr.Np, 'pm', arr.pm, 'muStar', muStar));
arrival = struct('params', arr, 'timeND', tT, 'stateND', rvT, ...
                 'periodND', infoT.periodND, 'closure', norm(rvT(end,:) - rvT(1,:)));

fprintf('\n2. TARGET %s: %d petals, branch %+d, period %.4f ND (%.3f d) LOCKED by Np\n', ...
        upper(arr.family), arr.Np, arr.pm, infoT.periodND, day(infoT.periodND));
fprintf('   out-of-plane extent max|z| = %.4f ND (%.0f km)\n', ...
        max(abs(rvT(:,3))), max(abs(rvT(:,3)))*lStar);
fprintf('   periodicity closure |x(T) - x(0)| = %.2e / %-7.0e %s\n', arrival.closure, ...
        tol.closure, pass(arrival.closure < tol.closure));
assert(arrival.closure < tol.closure, 'the target orbit is not periodic to tolerance');

%% ========================================================================
%  3. ENGINE, PHASES, ENDPOINTS
%     A phase is a fraction of its orbit's period. The endpoint state is read
%     from a PERIODIC SPLINE through the propagated orbit -- that interpolant,
%     not the raw table and not the dynamics, is what the boundary condition
%     matches, so its seam is checked in value AND derivative.
%% ========================================================================
thrustN = 0.070;      ispS = 900;      m0kg = 150;
sD      = 0.0;                                % departure phase, fraction of tau_D
sA      = 0.0754;                             % arrival phase,  fraction of tau_A
                                              % certified pairs: sD = k/12,
                                              %   sA = 0.0754 + j/12 (section 4)

g0  = 9.80665*tStar^2/(1000*lStar);           % ND gravity at sea level
cnd = (ispS/tStar)*g0;                        % ND exhaust speed
Tnd = (thrustN/m0kg)*tStar^2/(lStar*1000);    % ND thrust acceleration at m = 1

% PERIODIC interpolant, from the shared endpoint rule
% (costate_common/phase_state): an ordinary not-a-knot spline is not C1
% across the seam at s = 0 -- exactly where a phase near zero is evaluated
% -- so a periodic cubic is REQUIRED, not preferred, and phase_state
% refuses to build anything else unless asked. seam.value is the ORBIT's
% own closure (no interpolant can improve it); seam.deriv is the
% interpolant's derivative mismatch across the seam.
[stateD, seamD] = phase_state(tD, rvD);
[stateA, seamA] = phase_state(tT, rvT);
departure.stateAtPhase = stateD;   arrival.stateAtPhase = stateA;
rv0 = stateD(sD);   rvf = stateA(sA);
departure.endpointND = rv0;        arrival.endpointND = rvf;

fprintf('\n3. ENGINE %.0f mN, Isp %g s, m0 %g kg  ->  T_nd = %.6e, c_nd = %.4f\n', ...
        thrustN*1000, ispS, m0kg, Tnd, cnd);
fprintf('   departure phase %.4f -> r = [%+.5f %+.5f %+.5f]\n', sD, rv0(1:3));
fprintf('   arrival   phase %.4f -> r = [%+.5f %+.5f %+.5f]\n', sA, rvf(1:3));
fprintf('   periodic cubic seam: value %.1e / %.1e, derivative %.1e / %.1e   %s\n', ...
        seamD.value, seamA.value, seamD.deriv, seamA.deriv, ...
        pass(max(seamD.value, seamA.value) < tol.closure && max(seamD.deriv, seamA.deriv) < tol.seamDeriv));
assert(max(seamD.value, seamA.value) < tol.closure && max(seamD.deriv, seamA.deriv) < tol.seamDeriv, ...
       'the endpoint interpolant is not periodic across its seam');

%% ========================================================================
%  4. SOLVE -- seed, then multiple shooting on the PMP boundary-value problem
%     Unknowns: the seven initial costates and the final time.
%% ========================================================================
% The seed library was built at ONE operating point. Check the WHOLE of it --
% period, petal count, branch, thrust, Isp and mass -- not just two of the six,
% or a seed from a different engine would silently start the solve.
% The library includes the certified 70 mN CATALOG -- 115 phase pairs, not
% only the ten solutions kept in the anchor and sweep files.
lib = dro_tulip_library([], struct('includeCatalog', true));
match = find(abs([lib.sD] - mod(sD,1)) < 1e-6 & abs([lib.sA] - mod(sA,1)) < 1e-6, 1);
libOp = struct('tau', 1.0, 'Np', 7, 'pm', -1, 'thrustN', 0.070, 'ispS', 900, 'm0kg', 150);
sameOp = abs(dep.tau - libOp.tau) < 1e-12 && arr.Np == libOp.Np && arr.pm == libOp.pm && ...
         abs(thrustN - libOp.thrustN) < 1e-12 && abs(ispS - libOp.ispS) < 1e-9 && ...
         abs(m0kg - libOp.m0kg) < 1e-9;
assert(~isempty(match) && sameOp, 'transfer_study:noSeed', '%s', sprintf( ...
    ['no seed for this operating point and phase pair. The library covers the\n' ...
     'tau = 1 DRO to the 7-petal tulip (pm = -1) at 70 mN, Isp 900 s, 150 kg, at\n' ...
     'its 115 certified phase pairs (sD = k/12, sA = 0.0754 + j/12). For any other\n' ...
     'case, walk to it with   T = run_dro_tulip(sD, sA, opts)   and study T here.']));

if isempty(lib(match).Y)
    % a CATALOG entry stores z8 but no junction states: fly z8 once and cut
    % the flight into junctions (seed_from_z8), which puts the seed AT the root
    seed = seed_from_z8(lib(match).z, rv0(1:6), 24, Tnd, cnd, muStar);
else
    K = size(lib(match).Y, 2);
    seed = struct('tf', lib(match).z(8), 'tGrid', linspace(0, lib(match).z(8), K+1), ...
                  'Y', [lib(match).Y, lib(match).Y(:,end)]);
    seed.Y(1:7,1) = [rv0(1:6); 1];   seed.Y(8:14,1) = lib(match).z(1:7);
end
fprintf('\n4. SEED from the %s entry at phase pair (%.4f, %.4f), %d segments\n', ...
        lib(match).src, lib(match).sD, lib(match).sA, size(seed.Y, 2) - 1);

% Multiple-Shooting min-time solve
[z8, it] = ms_tfmin(rv0(1:6), rvf(1:6), seed, Tnd, cnd, muStar, ...
                    struct('tolR', 3e-11, 'wallSec', 600, 'conjTest', true));

% ENFORCE convergence before anything downstream. A best-iterate from an
% unsuccessful solve would otherwise be verified, plotted and described in
% the language of optimality.
assert(it.converged && isfinite(it.normR) && it.normR < tol.R, ...
       'the shooting solve did not converge (|R| = %.2e): nothing below is meaningful', it.normR);

% ONE FLIGHT, then its admissibility through the SAME validator the
% certifier uses: reached t_f, finite, all-burn mass law, clear of both
% primaries. A returned array is not a completed flight.
[tu, Y] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(1:6); 1; z8(1:7)], Tnd, cnd, muStar);
flight = struct('t', tu, 'Y', Y, 'z8', z8, 'rv0', rv0(1:6), 'rvf', rvf(1:6), ...
                'Tnd', Tnd, 'cnd', cnd, 'muStar', muStar, 'tStar', tStar, ...
                'lStar', lStar, 'm0kg', m0kg, 'nSamples', numel(tu));
VF = validate_flight(flight.t, flight.Y, z8(8), Tnd, cnd, muStar, lStar);
assert(VF.ok, 'the flight is inadmissible: %s', VF.reason);
flight.admissibility = VF;
mf = flight.Y(end,7);
dV = cnd*log(1/mf)*lStar/tStar;
fprintf('\n4. SOLVED: %d Newton iterations, |R| = %.2e, converged = %d\n', ...
        it.iters, it.normR, it.converged);
fprintf('   t_f = %.6f ND (%.4f d)   Delta-V = %.4f km/s   propellant %.2f kg\n', ...
        z8(8), day(z8(8)), dV, m0kg*(1-mf));
fprintf('   lambda_0 = [%s]\n', strjoin(compose('%+.6g', z8(1:7)'), ' '));
fprintf('   flight: %d samples, %s; closest approach %.0f km (Moon) / %.0f km (Earth)\n', ...
        flight.nSamples, VF.reason, VF.moonKm, VF.earthKm);

%% ========================================================================
%  5. INDEPENDENT VERIFICATION -- a second solver must not move the costates,
%     and its own answer must fly to the target too.
%% ========================================================================
B = struct('problem', struct('lStar', lStar, 'tStar', tStar, 'muStar', muStar, ...
                             'thrustN', thrustN, 'ispS', ispS, 'm0kg', m0kg, ...
                             'tauDRO', dep.tau, 'NpTulip', arr.Np, 'pmTulip', arr.pm), ...
           'Tnd', Tnd, 'cnd', cnd, 'mu', muStar, 'stateD', stateD, 'stateA', stateA);
V = verify_with_pumpkyn(struct('z', z8, 'sD', sD, 'sA', sA), B, ...
                        struct('tolDz', tol.dz, 'gateKm', tol.km, 'gateVms', tol.ms));

%% ========================================================================
%  6. NECESSARY CONDITIONS (Pontryagin, first order) -- one at a time
%% ========================================================================
fprintf('\n6. NECESSARY CONDITIONS      (value / threshold)\n');
tf = flight.t(end);
Yf = flight.Y;   tfl = flight.t;
lamV = Yf(:, 11:13);   lamM = Yf(:, 14);   mass = Yf(:, 7);
lamVmag = vecnorm(lamV, 2, 2);                % |lambda_v|, the Legendre quantity

% N1  the boundary-value residual: costate equations, terminal matching,
%     transversality. This IS the statement that the first variation vanishes.
[~, chk] = ms_tfmin(rv0(1:6), rvf(1:6), seed, Tnd, cnd, muStar, struct('assembleOnly', true));
R1 = chk.residual([z8(1:7); reshape(it.Y(:,2:end), [], 1); z8(8)]);
resid = norm(R1, inf);
fprintf('   N1 BVP residual   |R|_inf   %9.2e / %-9.0e  %s\n', resid, tol.R, pass(resid < tol.R));

% N2  the Hamiltonian. Autonomous problem, free final time  =>  H == 0.
%     H = 1 + lam_r.v + lam_v.(g + (T/m) alpha) - lam_m T/c,  alpha = -lam_v/|lam_v|
Hval = zeros(size(tfl));
for k = 1:numel(tfl)
    F = mintime_rhs_point(Yf(k,:).', Tnd, cnd, muStar);     % [xdot; lamdot]
    Hval(k) = 1 + Yf(k,8:14)*F(1:7);
end
Hmax = max(abs(Hval));
fprintf('   N2 Hamiltonian    max|H|    %9.2e / %-9.0e  %s\n', Hmax, tol.H, pass(Hmax < tol.H));

% N3  the flight actually reaches the target, in position AND velocity
missKm  = norm(Yf(end,1:3) - rvf(1:3)')*lStar;
missVms = norm(Yf(end,4:6) - rvf(4:6)')*lStar/tStar*1000;
fprintf('   N3 arrival        %.4f km / %-4.0f km, %.4f m/s / %-3.0f m/s  %s\n', ...
        missKm, tol.km, missVms, tol.ms, pass(missKm < tol.km && missVms < tol.ms));

% N4  transversality on the free mass: lam_m(t_f) = 0
lamMf = abs(Yf(end,14));
fprintf('   N4 transversality |lam_m(t_f)| %6.2e / %-9.0e  %s\n', lamMf, tol.lamm, pass(lamMf < tol.lamm));

% N5  the ADJOINT equations themselves: lambda-dot = -dH/dx. A small shooting
%     residual says the pieces MATCH each other; it does not say the costate
%     equations are the right ones. Differentiate H in the STATE at fixed
%     costate, at two step sizes scaled to each coordinate, combine them
%     (Richardson) and compare with the costate rate the field returns.
%     Differencing the propagator's OUTPUT instead measures its sample
%     spacing: on this arc that read 5e-4, all truncation at the lunar pass.
kk = unique(round(linspace(2, numel(tfl)-1, 120)));
hRel = 1e-6;  adjErr = 0;  fdAgree = 0;
for k = kk
    yk = Yf(k,:).';   lam = yk(8:14);
    F  = mintime_rhs_point(yk, Tnd, cnd, muStar);
    h1 = hRel*max(1, abs(yk(1:7)));
    g1 = zeros(7,1);  g2 = zeros(7,1);
    for jj = 1:7
        for step = 1:2
            h = h1(jj)/step;
            yp = yk;  yp(jj) = yp(jj) + h;
            ym = yk;  ym(jj) = ym(jj) - h;
            Fp = mintime_rhs_point(yp, Tnd, cnd, muStar);
            Fm = mintime_rhs_point(ym, Tnd, cnd, muStar);
            gj = (lam.'*Fp(1:7) - lam.'*Fm(1:7))/(2*h);
            if step == 1, g1(jj) = gj; else, g2(jj) = gj; end
        end
    end
    gR = (4*g2 - g1)/3;
    adjErr  = max(adjErr,  norm(F(8:14) + gR)/max(norm(gR), 1));
    fdAgree = max(fdAgree, norm(g1 - g2)/max(norm(gR), 1));
end
fprintf('   N5 adjoint eqns   rel err   %9.2e / %-9.0e  %s   (FD steps agree to %.1e)\n', ...
        adjErr, tol.adj, pass(adjErr < tol.adj), fdAgree);

% N6  the MINIMUM principle, EXACTLY, for the control the propagator APPLIED.
%     H is affine in the direction alpha with coefficient (T/m) lam_v, so its
%     minimiser over the unit sphere is -lam_v/|lam_v| ANALYTICALLY: sampling
%     the sphere against that formula tests nothing. What CAN be wrong is the
%     propagator's control law. Recover the thrust acceleration it applied as
%     (powered field - coasting field) and evaluate the gap
%        H(applied) - min_alpha H = (T/m)(lam_v . alpha_applied + |lam_v|) >= 0
%     which is zero exactly when the applied direction is the minimiser. The
%     throttle enters H linearly with slope -T Q, Q = |lam_v|/m + lam_m/c, so
%     u = 1 minimises H iff Q >= 0 (weak, necessary); S2 asks for Q > 0.
dirGap = 0;  throttleErr = 0;
for k = kk
    yk = Yf(k,:).';   m = yk(7);   lv = yk(11:13);
    F  = mintime_rhs_point(yk, Tnd, cnd, muStar);
    F0 = mintime_rhs_point(yk, 0,   cnd, muStar);
    alphaApplied = (F(4:6) - F0(4:6))*m/Tnd;               % u * alpha as flown
    throttleErr  = max(throttleErr, abs(norm(alphaApplied) - 1));
    dirGap = max(dirGap, (Tnd/m)*(lv.'*alphaApplied + norm(lv)));
end
Qmt = lamVmag./mass + lamM/cnd;
fprintf('   N6 min principle  gap       %9.2e / %-9.0e  %s   (throttle 1 to %.1e; min Q = %.3f >= 0 %s)\n', ...
        dirGap, tol.gap, pass(dirGap <= tol.gap && throttleErr < 1e-10), throttleErr, min(Qmt), pass(min(Qmt) >= 0));

necessary = resid < tol.R && Hmax < tol.H && missKm < tol.km && missVms < tol.ms && ...
            lamMf < tol.lamm && adjErr < tol.adj && dirGap <= tol.gap && ...
            throttleErr < 1e-10 && min(Qmt) >= 0;

% X1  the independent solve is NOT one of the above: it is a cross-check on
%     our implementation, reported apart and excluded from `necessary`, but
%     a failure still blocks the claim -- an implementation in doubt cannot
%     certify anything.
fprintf('   X1 cross-check    second solver |dz| %6.2e / %-7.0e, its flight %.4f km / %.4f m/s  %s\n', ...
        V.dz, tol.dz, V.flyKm, V.flyVms, pass(V.ok));
crossCheck = V.ok;

%% ========================================================================
%  7. SUFFICIENCY HYPOTHESES (Bonnard-Caillau-Trelat) -- one at a time
%     With section 6, these give a strict STRONG local minimizer among
%     trajectories with the same endpoints.
%% ========================================================================
fprintf('\n7. SUFFICIENCY HYPOTHESES\n');

% S1  strengthened Legendre. For a direction on the unit sphere the second
%     derivative of H in the control, restricted to that sphere, is
%     (T/m)|lam_v| I -- positive definite exactly when |lam_v| > 0.
[minLamV, iMin] = min(lamVmag);
fprintf('   S1 Legendre      min|lam_v| = %.4e at t/t_f = %.3f   %s\n', ...
        minLamV, tfl(iMin)/tf, pass(minLamV > 0));

% S2  STRICT bang: the switching function stays strictly positive, so the
%     throttle is determined (no singular arc) -- the strict form of N6's
%     weak condition.
fprintf('   S2 strict bang   min Q = %.4e                        %s\n', ...
        min(Qmt), pass(min(Qmt) > 0));

% S3  normality: no abnormal lift of the SAME trajectory (dim S = 1)
gates = mintime_hypothesis_gates(z8, rv0(1:6), Tnd, cnd, muStar, struct());
fprintf('   S3 normality     dim S = %d (1 = no abnormal lift)        %s\n', ...
        gates.dimS, pass(gates.dimS == 1));
% the dim S number is only as good as the lift it was measured around, so
% show the two self-consistency residuals it rests on rather than hiding them
fprintf('      lift residual %.1e, |lambda.f + 1| %.1e, sv gap %.1e\n', ...
        gates.nullResid, gates.Hresid, gates.svRatio);

% S4  no conjugate time in (0, t_f], by the free-time quotiented Jacobi test.
%     PASS/FAIL is not the whole vocabulary: the instrument can also return
%     ENDPOINT, which is UNRESOLVED rather than either. A verdict is only
%     interpretable beside its COVERAGE, reported here as t/t_f.
cj = it.conj;
s4Status = 'FAIL';
if cj.pass == 1
    s4Status = 'PASS';
elseif strcmpi(cj.verdict, 'ENDPOINT')
    s4Status = 'UNRESOLVED';
end
fprintf('   S4 conjugate     %s, %d crossing(s), min|det| %.2e       %s\n', ...
        cj.verdict, gv(cj, 'nCrossings'), min(abs(cj.detScaled)), s4Status);
fprintf(['      coverage: %d junctions, first full-rank at t/t_f = %.3f, sampled through\n' ...
         '      t/t_f = %.3f (uncovered final interval %.3f); a sign test cannot see an\n' ...
         '      even-order zero or two zeros inside one segment\n'], ...
        numel(cj.detScaled), cj.t(cj.firstFullRank)/tf, cj.sampledThrough/tf, 1 - cj.sampledThrough/tf);

% V1  H6: the reduced conjugate instrument's determinant can vanish
%     spuriously when lambda_m(0) >= c/T (FINDINGS 40). Not a hypothesis of
%     the theorem -- the VALIDITY condition of S4's instrument, with margin.
if isfield(gates, 'h6Margin') && isfinite(gates.h6Margin)
    v1Status = pass(gates.h6Margin >= tol.h6Margin);
    fprintf('   V1 H6 validity   lambda_m(0) %.3f vs c/T %.3f, margin %5.1fx / %.1fx  %s\n', ...
            gates.h6LamM0, gates.h6Threshold, gates.h6Margin, tol.h6Margin, v1Status);
else
    v1Status = 'NOT CHECKED';
    fprintf('   V1 H6 validity   NOT CHECKED (the gates returned no margin): S4 is not interpretable\n');
end

% V2  CONSISTENCY of the inline numbers with the library instruments. Exposing
%     the scaffolding risks growing a second, unverified implementation of the
%     tests; this binds the two -- and runs BEFORE any verdict is printed.
PW = pmp_pointwise_checks(flight.t, flight.Y, Tnd, cnd, muStar);
agreeErr = max([abs(minLamV - gates.minLamV)/max(gates.minLamV, 1), ...
                abs(min(Qmt) - gates.minQmt)/max(gates.minQmt, 1), ...
                abs(Hmax - PW.Hmax)/max(PW.Hmax, 1e-9), ...
                abs(adjErr - PW.adjErr)/max(PW.adjErr, 1e-9), ...
                abs(dirGap - PW.dirGap)]);
assert(agreeErr < tol.agree, ...
       'the inline checks disagree with the library instruments (worst %.1e)', agreeErr);
fprintf('   V2 consistency   inline S1/S2/N2/N5/N6 vs library instruments, worst %.1e / %.0e  %s\n', ...
        agreeErr, tol.agree, pass(agreeErr < tol.agree));

sufficient = minLamV > 0 && min(Qmt) > 0 && gates.dimS == 1 && strcmp(s4Status, 'PASS') && ...
             strcmp(v1Status, 'PASS');
unresolved = strcmp(s4Status, 'UNRESOLVED') || strcmp(v1Status, 'NOT CHECKED');

%% ------------------------------------------------------------------------
%  VERDICT. The same wording as report_optimality, and the same four-way
%  outcome: claim / not an extremal / unresolved / refuted. Strong local
%  minimality is a property of the whole trajectory, so the sampling
%  qualifies the EVIDENCE, not the claim.
%% ------------------------------------------------------------------------
fprintf('\n   VERDICT: ');
if necessary && sufficient && crossCheck
    fprintf(['All required numerical checks passed. If the stated hypotheses hold\n' ...
             '            exactly, this arc is a strict strong local minimizer among trajectories\n' ...
             '            with the same endpoints and phases.\n']);
    fprintf(['            EVIDENCE IS NUMERICAL AND SAMPLED: positivity is tested at the\n' ...
             '            sampled times, dim S is a numerical rank, and the conjugate test\n' ...
             '            is a sign test at %d junctions -- an even-order zero, or two zeros\n' ...
             '            inside one segment, would not be seen. This is a strong numerical\n' ...
             '            audit, not a proof. See doc/mintime_second_order_audit.tex.\n'], ...
             size(it.Y, 2));
elseif ~necessary
    fprintf(['NOT an extremal to tolerance. The second-order test is meaningless\n' ...
             '            off an extremal, so no minimality is claimed.\n']);
elseif ~crossCheck
    fprintf(['every PMP and sufficiency line passed, but the CROSS-CHECK failed: the\n' ...
             '            implementation is in doubt, so no claim is made.\n']);
elseif unresolved
    fprintf(['an extremal; sufficiency is UNRESOLVED (%s / %s): no minimality is\n' ...
             '            claimed and none is refuted.\n'], s4Status, v1Status);
else
    fprintf('an extremal, but a sufficiency hypothesis fails: no minimality claimed.\n');
end

%% ========================================================================
%  8. INTERACTIVE 3D PLOT  (drag to rotate, scroll to zoom) -- drawn from the
%     SAME flight as every number above.
%% ========================================================================
T = struct('z', z8, 'sD', sD, 'sA', sA, 'tfDays', day(z8(8)), 'dvKms', dV, ...
           'propellantKg', m0kg*(1-mf), 'finalMassKg', m0kg*mf);
P = plot_transfer_3d(T, B, struct('flight', flight));
fprintf('\n8. Figure %d is rotatable (flight supplied: %d).\n', P.fig.Number, P.flightSupplied);

%% ------------------------------------------------------------------------
function s = pass(c)
% PASS  Verdict text.  INPUTS: c.  OUTPUTS: s.
if c, s = 'PASS'; else, s = 'FAIL'; end
end

function v = gv(s, f)
% GV  Field or NaN.  INPUTS: s; f.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = NaN; end
end
