%% TRANSFER_STUDY  One minimum-time low-thrust transfer: DRO-to-Tulip
%
%   Edit the parameter blocks, press Run. The STEPS are done here rather than
%   behind a front door, so the machinery is visible. Five computations stay
%   in the library because reimplementing them here would create a second,
%   unverified copy: the lift-space rank and its Eckart-Young margin (S3),
%   the Jacobi determinant and the dense spectrum scan (S4) and the H6
%   margin (V1). Section 8 asserts the inline numbers against that library
%   for the same reason.
%
%   PREREQUISITE. The solve needs a converged SEED, and section 4 finds it
%   (dro_tulip_seed). The certified library covers ONE operating point --
%   the tau = 1 DRO to the 7-petal tulip (pm = -1) at 70 mN, Isp 900 s,
%   150 kg -- at its 115 certified phase pairs: the grid sD = k/12,
%   sA = 0.0754 + j/12, less the cells the catalog does not hold. Anything
%   else is refused by name, :noSeed for the phase pair and :operatingPoint
%   for the engine. For any other case, walk to it first with run_dro_tulip
%   and study its result here.
%
%     0  tolerances                       (one source, set before anything runs)
%     1  generate the DEPARTURE orbit      (family + parameters -> a real orbit)
%     2  generate the TARGET orbit         (same)
%     3  engine, phases, endpoints         (nondimensionalisation spelled out)
%     4  get the seed                      (which certified solution, and is
%                                          it the right engine?)
%     5  solve                             (multiple shooting -> costates)
%     6  independent verification          (a second solver must not move them)
%     7  NECESSARY conditions              (Pontryagin, one at a time)
%     8  SUFFICIENCY hypotheses            (Bonnard-Caillau-Trelat, one at a time)
%     9  interactive 3D plot               (drag to rotate)
%
%   Diagnostic IDs are stable and grouped, and the letters mean:
%     N  NECESSARY   -- N1-N6, Pontryagin's necessary conditions. Fail one
%                      and this is not an extremal, so nothing below it
%                      means anything.
%     S  SUFFICIENCY -- S1-S4, the HYPOTHESES of the Bonnard-Caillau-Trelat
%                      sufficiency theorem (its hypotheses, not its
%                      conclusions: satisfying them, with N1-N6, is what
%                      buys the strict strong local minimum).
%     V  VALIDITY    -- V1-V2, whether a verdict above can be believed at
%                      all: V1 is the conjugate instrument's own
%                      precondition (H6), V2 is a wiring check that the
%                      instruments were handed the same flight and physics.
%     X  CROSS-CHECK -- X1, a second SHOOTING implementation against ours
%                      (both propagate pumpkyn's field, so it checks the
%                      solver, not the physics); X2, pumpkyn's field against
%                      an independently written CR3BP + thrust field, row by
%                      row (that IS the physics check). Neither is a
%                      condition of the theory; a failure still blocks the
%                      claim.
%
%   Sections 7 and 8 are separate because the theory separates them: the
%   first-order conditions are HYPOTHESES of the sufficiency theorem, not
%   consequences of it, and the conjugate test is computed ALONG the extremal
%   -- off one, its determinant means nothing. On this problem 12 of 14
%   candidates satisfied every first-order condition and were then refuted.
%
%   Section 6 is not an optimality condition at all. It guards against OUR
%   solver: a bug in our shooting could give a self-consistent answer to the
%   wrong problem, and only a second implementation catches that. It shares
%   pumpkyn's equations with us, so a wrong FIELD would pass both; X2 in
%   section 7 is the check against an independently written field.
%
%   ONE FLIGHT. Section 5 flies the converged costates once, into `flight`;
%   every number in sections 7-9 and the figure come from that object, so the
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
    'interp',    1e-6,  ...   % endpoint interpolant vs a propagation to the
    ...                       %   same phase (ND; 1e-6 = 0.39 km) -- this is
    ...                       %   the physical error of the boundary condition
    'km',        100,   ...   % N3 flown arrival, position (a SCREEN on the
    'ms',        10,    ...   %   single-shot flight; the fixed-endpoint
    ...                       %   statement is N1's residual)
    'lamm',      1e-6,  ...   % N4 transversality lambda_m(t_f) = 0
    'adj',       1e-7,  ...   % N5 adjoint equations, relative
    'gap',       1e-12, ...   % N6 |full minimum-principle gap| of the applied control
    'throttle',  1e-10, ...   % N6 applied throttle, acceleration AND mass rows
    'dz',        1e-6,  ...   % X1 second solver must agree to this
    'field',     1e-10, ...   % X2 pumpkyn field vs independent field, relative
    'lift',      1e-4,  ...   % S3 accepted lift's own residual |C lam0|/|lam0|;
    ...                       %   = lift_margin's liftTol (one value, handed to it).
    ...                       %   The floor is the pchip interpolation of the flown
    ...                       %   arc that the frozen-control adjoint is integrated
    ...                       %   over: 2.2e-6 on the anchor at BOTH relTol 1e-12
    ...                       %   and 1e-10, so it is not the integration tolerance
    'liftMargin', 10,   ...   % S3 Eckart-Young margin sigma_6 / |dC| (lift_margin)
    'h6Margin',  1.0,   ...   % V1 required (c/T)/lambda_m(0), STRICT
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

% The seam check is a property of the interpolant's CONSTRUCTION; it says
% nothing about how far the interpolated endpoint sits from the orbit
% (Astra review 2026-09-11). So propagate each orbit from its first sample
% to the chosen phase and compare: that difference is the physical error of
% the boundary condition the shooting then matches to 1e-11.
interpErr = [phase_interp_error(stateD, rvD, sD, infoD.periodND, muStar), ...
             phase_interp_error(stateA, rvT, sA, infoT.periodND, muStar)];
fprintf('   interpolant vs propagation to the phase: %.1e / %.1e ND (%.4f / %.4f km) / %-6.0e %s\n', ...
        interpErr, interpErr*lStar, tol.interp, pass(max(interpErr) < tol.interp));
assert(max(interpErr) < tol.interp, 'the endpoint interpolant is off the orbit by more than tol.interp');

%% ========================================================================
%  4. GET THE SEED -- which certified solution starts the solve, and is it
%     the right engine? Both questions belong together, so they are one
%     call: dro_tulip_seed matches the phase pair EXACTLY (never nearest --
%     a neighbouring phase can seed a slower branch that then passes every
%     first-order check), checks the operating point in all six fields, and
%     builds the seed through costate_common/seed_from_entry.
%% ========================================================================
op = struct('tau', dep.tau, 'Np', arr.Np, 'pm', arr.pm, ...
            'thrustN', thrustN, 'ispS', ispS, 'm0kg', m0kg, ...
            'rv0', rv0, 'Tnd', Tnd, 'cnd', cnd, 'muStar', muStar, ...
            'lStar', lStar, 'tStar', tStar, 'K', 24);
[seed, seedInfo] = dro_tulip_seed(sD, sA, op);

fprintf('\n4. SEED from the %s entry at (%.4f, %.4f): %d segments, its t_f %.4f d\n', ...
        seedInfo.src, seedInfo.sD, seedInfo.sA, size(seed.Y, 2) - 1, seedInfo.tfDays);
fprintf('   operating point checked in all six fields, %d certified pairs on file\n', ...
        seedInfo.nLibrary);

%% ========================================================================
%  5. SOLVE -- multiple shooting on the PMP boundary-value problem.
%     Unknowns: the seven initial costates and the final time.
%% ========================================================================
[z8, it] = ms_tfmin(rv0(1:6), rvf(1:6), seed, Tnd, cnd, muStar, ...
                    struct('tolR', 3e-11, 'wallSec', 600, 'conjTest', true));


% TWO POLICY GATES, kept here because they are the SCRIPT's policy, not the
% library's: nothing below is meaningful from a best-iterate of an
% unsuccessful solve, or from a flight that is not admissible.
assert(it.converged && isfinite(it.normR) && it.normR < tol.R, ...
       'the shooting solve did not converge (|R| = %.2e): nothing below is meaningful', it.normR);

% ONE FLIGHT (costate_common/fly_transfer): flown once, validated through
% the SAME validator the certifier uses -- reached t_f, finite, all-burn
% mass law, clear of both primaries -- with the mass and Delta-V that follow
% from it. Every number in sections 7-9 and the figure come from this object.
flight = fly_transfer(z8, rv0(1:6), rvf(1:6), op);
assert(flight.admissibility.ok, 'the flight is inadmissible: %s', flight.admissibility.reason);

fprintf('\n5. SOLVED: %d Newton iterations, |R| = %.2e, converged = %d\n', ...
        it.iters, it.normR, it.converged);
print_transfer_summary(flight);
fprintf('   lambda_0 = [%s]\n', strjoin(compose('%+.6g', z8(1:7)'), ' '));

%% ========================================================================
%  6. INDEPENDENT VERIFICATION -- a second solver must not move the costates,
%     and its own answer must fly to the target too.
%% ========================================================================
B = struct('problem', struct('lStar', lStar, 'tStar', tStar, 'muStar', muStar, ...
                             'thrustN', thrustN, 'ispS', ispS, 'm0kg', m0kg, ...
                             'tauDRO', dep.tau, 'NpTulip', arr.Np, 'pmTulip', arr.pm), ...
           'Tnd', Tnd, 'cnd', cnd, 'mu', muStar, 'stateD', stateD, 'stateA', stateA);
V = verify_with_pumpkyn(struct('z', z8, 'sD', sD, 'sA', sA), B, ...
                        struct('tolDz', tol.dz, 'gateKm', tol.km, 'gateVms', tol.ms));

%% ========================================================================
%  7. NECESSARY CONDITIONS (Pontryagin, first order) -- one at a time
%% ========================================================================
fprintf('\n7. NECESSARY CONDITIONS      (value / threshold)\n');
% THE production instrument -- the same call certify_root makes on every
% certified entry: H along the arc, transversality, the adjoint equations
% and the exact minimum-principle gap of the control the propagator applied.
% The script used to recompute all four and compare (V2); that comparison
% measured 0.0e+00, because a copy of the same arithmetic on the same
% samples is not an independent implementation. Running the instrument the
% catalogs are certified with is the stronger statement.
PW = pmp_pointwise_checks(flight.t, flight.Y, Tnd, cnd, muStar);

tf = flight.t(end);                           % S4 reports its coverage against it

% N1  the boundary-value residual: junction matching, terminal state,
%     transversality and H(t_f) = 0 -- the PMP SHOOTING EQUATIONS, satisfied
%     to tolerance. Not "the first variation vanishes": on a boundary
%     control the first-order condition in the control is the minimum
%     principle, an inequality, which is N6 (Astra review 2026-09-11).
%     It is the solve's OWN residual, at the point it returned -- section 5
%     already computed it, and re-solving to get it back would report a
%     different number for a reason worth knowing: the residual depends on
%     HOW the arc is propagated. Measured on this anchor, the same point
%     gives 7.84e-12 when the Jacobian is requested (the integrator carries
%     210 states and takes finer steps, which is the mode the solver ran in
%     and the mode every campaign's tolR gate is evaluated in) and 1.41e-09
%     under plain 14-state propagation. Quoting the second as "the BVP
%     residual" would understate the solution by 180x in the wrong
%     direction.
resid = it.normR;
fprintf('   N1 BVP residual   |R|_inf   %9.2e / %-9.0e  %s\n', resid, tol.R, pass(resid < tol.R));

% N2  the Hamiltonian. Autonomous problem, free final time  =>  H == 0.
%     H = 1 + lam_r.v + lam_v.(g + (T/m) alpha) - lam_m T/c,  alpha = -lam_v/|lam_v|
%     On an exact solution of N1 this follows from H(t_f) = 0 plus canonical
%     propagation, and N4 is one of N1's equations: they are RE-EVALUATIONS
%     on an independent flight, which exposes propagation defects, not
%     additional restrictions on the extremal.
fprintf('   N2 Hamiltonian    max|H|    %9.2e / %-9.0e  %s\n', PW.Hmax, tol.H, pass(PW.Hmax < tol.H));

% N3  the flight actually reaches the target, in position AND velocity.
%     fly_transfer measured this in section 5 from the same flight; there is
%     no second copy of the formula here. The tolerances are a SCREEN on the
%     single-shot flight (the propagated miss grows with the arc); the
%     fixed-endpoint statement itself is N1's residual.
fprintf('   N3 arrival        %.4f km / %-4.0f km, %.4f m/s / %-3.0f m/s  %s\n', ...
        flight.flyKm, tol.km, flight.flyVms, tol.ms, ...
        pass(flight.flyKm < tol.km && flight.flyVms < tol.ms));

% N4  transversality on the free mass: lam_m(t_f) = 0
fprintf('   N4 transversality |lam_m(t_f)| %6.2e / %-9.0e  %s\n', PW.lamMf, tol.lamm, pass(PW.lamMf < tol.lamm));

% N5  the ADJOINT equations themselves: lambda-dot = -dH/dx. A small shooting
%     residual says the pieces MATCH each other; it does not say the costate
%     equations are the right ones. Differentiate H in the STATE at fixed
%     costate, at two step sizes scaled to each coordinate, combine them
%     (Richardson) and compare with the costate rate the field returns.
%     Differencing the propagator's OUTPUT instead measures its sample
%     spacing: on this arc that read 5e-4, all truncation at the lunar pass.
fprintf('   N5 adjoint eqns   rel err   %9.2e / %-9.0e  %s   (%d samples; FD steps agree to %.1e)\n', ...
        PW.adjErr, tol.adj, pass(PW.adjErr < tol.adj), PW.nSample, PW.fdAgree);

% N6  the MINIMUM principle, EXACTLY, for the control the propagator APPLIED.
%     H is affine in the direction alpha with coefficient (T/m) lam_v, so its
%     minimiser over the unit sphere is -lam_v/|lam_v| ANALYTICALLY: sampling
%     the sphere against that formula tests nothing. What CAN be wrong is the
%     propagator's control law. Recover the control it applied from the field
%     itself, on BOTH rows it enters -- b = u alpha from (powered - coasting)
%     acceleration, and u again from the mass row, -c F_m / T -- and evaluate
%     the gap against the FULL control minimum (direction AND throttle):
%        H(u, alpha) - min H = (T/m)(lam_v . b + |b||lam_v|) + T(max(Q,0) - |b| Q)
%     Q = |lam_v|/m + lam_m/c. Zero exactly at the minimiser; its signed
%     MINIMUM is kept too, so an over-unit thrust along the minimiser (a
%     negative gap) is not clipped away. u = 1 minimises H iff Q >= 0 (weak,
%     necessary); S2 asks for Q > 0.
n6 = PW.fullGap <= tol.gap && PW.throttleErr < tol.throttle && PW.minQmt >= 0;
fprintf('   N6 min principle  |gap|     %9.2e / %-9.0e  %s   (signed %+.1e..%+.1e; throttle 1 to %.1e acc / %.1e mass; min Q = %.3f >= 0 %s)\n', ...
        PW.fullGap, tol.gap, pass(n6), PW.gapMin, PW.gapMax, ...
        PW.throttleAccErr, PW.throttleMassErr, PW.minQmt, pass(PW.minQmt >= 0));

necessary = resid < tol.R && PW.Hmax < tol.H && flight.flyKm < tol.km && ...
            flight.flyVms < tol.ms && PW.lamMf < tol.lamm && PW.adjErr < tol.adj && n6;

% ONE call to the hypothesis gates -- the instrument that judges S1, S2, S3
% and V1 (section 8), the same one gates_catalog_pass ran over all 18,360
% catalog entries -- made HERE because it also carries X2. It flies its own
% dense trajectory from z8, so its numbers are not a re-reading of the
% script's flight. keepC at the tight tolerance, for lift_margin below.
gates = mintime_hypothesis_gates(z8, rv0(1:6), Tnd, cnd, muStar, ...
                                 struct('keepC', true, 'relTol', 1e-12));

% X1  the independent solve is NOT one of the above: it is a cross-check on
%     our implementation, reported apart and excluded from `necessary`, but
%     a failure still blocks the claim -- an implementation in doubt cannot
%     certify anything. It is a second SHOOTING implementation on pumpkyn's
%     field; agreement to |dz| pins its answer to ours (including the two
%     terminal multiplier conditions, which arrival alone would not).
fprintf('   X1 cross-check    second solver |dz| %6.2e / %-7.0e, its flight %.4f km / %.4f m/s  %s\n', ...
        V.dz, tol.dz, V.flyKm, V.flyVms, pass(V.ok));
% X2  the PHYSICS: pumpkyn's field (which every solver and every check above
%     propagates) against an independently written CR3BP + thrust field and
%     its CasADi adjoint, ROW BY ROW along the arc. N2 and N5 are consistency
%     checks of one field with itself; this is the seven-row comparison.
x2 = gates.fieldErr < tol.field && gates.adjErrRef < tol.field;
fprintf('   X2 physics        pumpkyn vs independent field: state rows %.1e, adjoint rows %.1e / %-6.0e %s\n', ...
        gates.fieldErr, gates.adjErrRef, tol.field, pass(x2));
crossCheck = V.ok && x2;

%% ========================================================================
%  8. SUFFICIENCY HYPOTHESES (Bonnard-Caillau-Trelat) -- one at a time
%     With section 6, these give a strict STRONG local minimizer among
%     trajectories with the same endpoints. Every one of them is a
%     CONTINUOUS statement tested on SAMPLES (the propagator's own dense
%     output, or the junctions): no between-sample bound is claimed here.
%% ========================================================================
fprintf('\n8. SUFFICIENCY HYPOTHESES\n');

% S1  strengthened Legendre. For a direction on the unit sphere the second
%     derivative of H in the control, restricted to that sphere, is
%     (T/m)|lam_v| I -- positive definite exactly when |lam_v| > 0.
fprintf('   S1 Legendre      min|lam_v| = %.4e at t/t_f = %.3f   %s\n', ...
        gates.minLamV, gates.tMinLamV/z8(8), pass(gates.minLamV > 0));

% S2  STRICT bang: the switching function stays strictly positive, so the
%     throttle is determined (no singular arc) -- the strict form of N6's
%     weak condition. On an EXACT lift this is not independent of the rest:
%     lam_m(t) = integral_t^tf T|lam_v|/m^2 >= 0 from N4 + N5, so S1 gives
%     Q > 0. It is kept as the margin diagnostic; a failure here means one
%     of those premises is not established.
fprintf('   S2 strict bang   min Q = %.4e at t/t_f = %.3f        %s\n', ...
        gates.minQmt, gates.tMinQmt/z8(8), pass(gates.minQmt > 0));

% S3  normality: no abnormal lift of the SAME trajectory. The gate builds
%     the space of STATIONARY lifts (frozen-control adjoint, lam_v || alpha,
%     lam_m(t_f) = 0) and reports its numerical dimension; dim S = 1 is the
%     corank-one condition that EXCLUDES an abnormal lift (the sufficient
%     direction; it is stronger than "no abnormal minimising lift"). The
%     rank is only as good as (a) the accepted lift's own residual, (b) the
%     Hamiltonian residual of the independent fixed field, and (c) the
%     Eckart-Young separation sigma_6 against the measured error in C --
%     the threshold rule alone forces at least one small singular value
%     (sigma_min <= nullResid by construction; Astra review 2026-09-11).
%     All three are REQUIRED here, not merely printed.
gatesLoose = mintime_hypothesis_gates(z8, rv0(1:6), Tnd, cnd, muStar, ...
                                      struct('keepC', true, 'relTol', 1e-9));
LM = lift_margin(gates.C, gatesLoose.C, z8(1:7), struct('marginMin', tol.liftMargin, 'liftTol', tol.lift));
s3 = gates.dimS == 1 && gates.nullResid < tol.lift && gates.Hresid < tol.H && LM.certified;
fprintf('   S3 normality     dim S = %d (1 = no abnormal lift)        %s\n', gates.dimS, pass(s3));
fprintf('      lift residual %.1e / %.0e, |lambda.f + 1| %.1e / %.0e, sv gap %.1e\n', ...
        gates.nullResid, tol.lift, gates.Hresid, tol.H, gates.svRatio);
fprintf('      Eckart-Young: sigma_6 %.2e vs |dC| %.2e between two integrations, margin %.0fx / %gx  %s\n', ...
        LM.sigma6, LM.errEst, LM.margin, tol.liftMargin, pass(LM.certified));

% S4  no conjugate time in (0, t_f], by the free-time quotiented Jacobi test
%     at the junctions, AND the dense singular-spectrum scan that closes the
%     sign test's blind spots (two zeros inside one segment, an even-order
%     touch). The junction test can return ENDPOINT (UNRESOLVED) and must
%     have COVERED the arc through t_f; the interval before its first
%     full-rank junction is uncovered by it and is sampled by the scan.
%     Two exact identities of the quotiented form are measured beside it.
cj = it.conj;
s4Status = 'FAIL';                                       % three outcomes, not two
if cj.pass == 1 && cj.covered,            s4Status = 'PASS';
elseif strcmpi(cj.verdict, 'ENDPOINT') || strcmpi(cj.verdict, 'UNDETERMINED') || ~cj.covered
                                          s4Status = 'UNRESOLVED';
end
fprintf('   S4 conjugate     %s (%s), %d crossing(s), min|det| %.2e       %s\n', ...
        cj.verdict, cj.reason, gv(cj, 'nCrossings'), min(abs(cj.detScaled)), s4Status);
fprintf(['      coverage: %d junctions, first full-rank at t/t_f = %.3f (before it: uncovered),\n' ...
         '      sampled through t/t_f = %.3f%s; identities J p(0) = 0 to %.1e, p(t)''J = 0 to %.1e\n'], ...
        numel(cj.detScaled), cj.tFirstFullRank/tf, cj.sampledThrough/tf, ...
        tern(cj.covered, '', ' (FINAL SEGMENT NOT COVERED)'), cj.kernelRight, cj.kernelLeft);
% the dense scan: 8 samples per segment, interior candidates located and
% refined twice (a zero keeps falling, a near-miss plateaus)
CS = conj_spectrum(z8, rv0(1:6), Tnd, cnd, muStar, struct('K', 24, 'nSub', 8));
s4Dense = CS.nInterior == 0 && CS.nZero == 0 && CS.multiplicity == 0;
fprintf(['      dense scan (%d samples): %d interior sign change(s), %d interior candidate(s):\n' ...
         '      %d zero, %d near-miss, %d multiplicity                                %s\n'], ...
        numel(CS.t), CS.nInterior, CS.nInteriorCand, CS.nZero, CS.nNearMiss, CS.multiplicity, pass(s4Dense));
if strcmp(s4Status, 'PASS') && ~s4Dense, s4Status = 'FAIL'; end

% V1  H6: the reduced conjugate instrument's determinant can vanish
%     spuriously when lambda_m(0) >= c/T (FINDINGS 40). Not a hypothesis of
%     the theorem -- the VALIDITY condition of S4's instrument. The gate is
%     the helper's OWN verdict (strict margin AND clearance above the
%     Hamiltonian residual), then the script's margin policy on top; the
%     margin ratio alone let equality through and ignored the clearance
%     (Astra review 2026-09-11). Valid for the normalisation p_0 = 1 this
%     solver uses.
v1Status = 'NOT CHECKED';                                % a missing margin BLOCKS S4
if isfield(gates, 'h6Ok') && isfield(gates, 'h6Margin') && isfinite(gates.h6Margin)
    v1Status = pass(gates.h6Ok && gates.h6Margin > tol.h6Margin);
    fprintf('   V1 H6 validity   lambda_m(0) %.3f vs c/T %.3f, margin %5.1fx > %.1fx, clearance %.3f > |H| %.1e  %s\n', ...
            gates.h6LamM0, gates.h6Threshold, gates.h6Margin, tol.h6Margin, ...
            gates.h6Clearance, gates.Hresid, v1Status);
else
    fprintf('   V1 H6 validity   NOT CHECKED (the gates returned no verdict): S4 is not interpretable\n');
end

% V2  the two INSTRUMENTS against each other on min Q: pmp_pointwise_checks
%     reads the flight this script flew, mintime_hypothesis_gates flies its
%     own from z8. A WIRING check, and worth calling it that -- both
%     propagate the same z8 deterministically, so it agrees to 0.0e+00 and
%     cannot detect an implementation divergence. What it CAN catch is the
%     script handing one of them the wrong flight or exhaust speed (min Q
%     depends on the flight and on c only); a wrong thrust or mass ratio
%     handed to the pointwise checks shows up in N2 instead, not here.
agreeErr = abs(PW.minQmt - gates.minQmt)/max(gates.minQmt, 1);
assert(agreeErr < tol.agree, ...
       'the pointwise checks and the gates disagree on min Q (%.1e)', agreeErr);
fprintf('   V2 wiring        min Q: pointwise vs gates, %.1e / %.0e       %s\n', ...
        agreeErr, tol.agree, pass(agreeErr < tol.agree));

sufficient = gates.minLamV > 0 && gates.minQmt > 0 && s3 && strcmp(s4Status, 'PASS') && ...
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
    fprintf(['The numerical diagnostics are consistent with a regular normal extremal\n' ...
             '            and with strict strong local minimality for fixed departure and\n' ...
             '            arrival position and velocity, initial mass fraction one, and free\n' ...
             '            terminal mass and time (phases are FIXED inputs, not optimised).\n' ...
             '            No determinant zero was detected on the %d-junction sign test or the\n' ...
             '            %d-sample spectrum scan. If the stated hypotheses hold exactly, the\n' ...
             '            Bonnard-Caillau-Trelat theorem gives the strict strong local minimum.\n'], ...
             size(it.Y, 2), numel(CS.t));
    fprintf(['            THIS IS NUMERICAL EVIDENCE, NOT A CERTIFICATE: positivity is tested\n' ...
             '            at sampled times with no between-sample bound, dim S is a numerical\n' ...
             '            rank with a measured (not proven) error, and the conjugate scan\n' ...
             '            locates and refines candidates rather than enclosing roots. See\n' ...
             '            doc/mintime_second_order_audit.tex.\n']);
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
%  9. INTERACTIVE 3D PLOT  (drag to rotate, scroll to zoom) -- drawn from the
%     SAME flight as every number above.
%% ========================================================================
T = struct('z', z8, 'sD', sD, 'sA', sA, 'tfDays', flight.tfDays, 'dvKms', flight.dvKms, ...
           'propellantKg', flight.propellantKg, 'finalMassKg', flight.finalMassKg);
P = plot_transfer_3d(T, B, struct('flight', flight));
fprintf('\n9. Figure %d is rotatable (flight supplied: %d).\n', P.fig.Number, P.flightSupplied);

%% ------------------------------------------------------------------------
function s = pass(c)
% PASS  Verdict text.  INPUTS: c.  OUTPUTS: s.
if c, s = 'PASS'; else, s = 'FAIL'; end
end

function v = gv(s, f)
% GV  Field or NaN.  INPUTS: s; f.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = NaN; end
end

function s = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end

function e = phase_interp_error(stateFun, rvTable, s, periodND, muStar)
% PHASE_INTERP_ERROR  Distance between the periodic interpolant at phase s
% and a CR3BP propagation of the orbit's first sample to the same phase.
% INPUTS: stateFun (phase -> [6x1]); rvTable [N x 6]; s [scalar];
% periodND; muStar.  OUTPUTS: e [scalar, ND, 6-vector norm].
if s <= 0
    rvProp = rvTable(1, 1:6).';
else
    [~, rvp] = pumpkyn.cr3bp.prop(s*periodND, rvTable(1, 1:6).', muStar);
    rvProp = rvp(end, 1:6).';
end
x = stateFun(s);
e = norm(x(1:6) - rvProp);
end
