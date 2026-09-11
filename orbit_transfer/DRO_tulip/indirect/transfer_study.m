%% TRANSFER_STUDY  One minimum-time low-thrust transfer: DRO-to-Tulip
%
%   Edit the parameter blocks, press Run. The STEPS are done here rather than
%   behind a front door, so the machinery is visible. Two computations stay in
%   the library because reimplementing them here would create a second,
%   unverified copy: the lift-space rank (S3) and the Jacobi determinant (S4).
%   Section 7 asserts its inline numbers against that library for the same
%   reason.
%
%     1  generate the DEPARTURE orbit      (family + parameters -> a real orbit)
%     2  generate the TARGET orbit         (same)
%     3  engine, phases, endpoints         (nondimensionalisation spelled out)
%     4  solve                             (seed -> multiple shooting -> costates)
%     5  independent verification          (a second solver must not move them)
%     6  NECESSARY conditions              (Pontryagin, one at a time)
%     7  SUFFICIENCY hypotheses            (Bonnard-Caillau-Trelat, one at a time)
%     8  interactive 3D plot               (drag to rotate)
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
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.

clear; clc
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
muStar = 0.012150585609624;                 % Earth-Moon mass ratio
lStar  = 389703.264829278;                  % km
tStar  = 382981.289129055;                  % s
day    = @(tnd) tnd*tStar/86400;

%% ========================================================================
%  1. DEPARTURE ORBIT -- a distant retrograde orbit (DRO)
%     Parameter: the PERIOD tau (ND). Larger tau = larger orbit.
%       tau 0.5 -> 2.22 d      tau 1.0 -> 4.43 d      tau 2.0 -> 8.87 d
%% ========================================================================
dep.family = 'dro';
dep.tau    = 1.0;                            % <-- change the departure orbit here

[tD, rvD, infoD] = get_family_orbit(dep.family, struct('tau', dep.tau, 'muStar', muStar));
departure = struct('params', dep, 'timeND', tD, 'stateND', rvD, ...
                   'periodND', infoD.periodND, 'closure', norm(rvD(end,:) - rvD(1,:)));
closureD = departure.closure;
fprintf('1. DEPARTURE %s: period %.4f ND (%.3f d), %d samples\n', ...
        upper(dep.family), infoD.periodND, day(infoD.periodND), numel(tD));
fprintf('   max radius from the barycentre %.4f ND (%.0f km)\n', ...
        max(vecnorm(rvD(:,1:3),2,2)), max(vecnorm(rvD(:,1:3),2,2))*lStar);
fprintf('   periodicity closure |x(T) - x(0)| = %.2e   %s\n', closureD, ...
        pass(closureD < 1e-8));

%% ========================================================================
%  2. TARGET ORBIT -- a tulip
%     Parameters: the PETAL COUNT Np and the branch pm (+1 is the z-mirror
%     of -1). The period is NOT free: it is locked to 2*pi*(Np-2)/(Np-1).
%       Np 3 -> 13.93 d    Np 5 -> 20.89 d    Np 7 -> 23.21 d    Np 9 -> 24.37 d
%% ========================================================================
arr.family = 'tulip';
arr.Np     = 7;                              % <-- change the target orbit here
arr.pm     = -1;                             %     branch: -1 or +1

[tT, rvT, infoT] = get_family_orbit(arr.family, struct('Np', arr.Np, 'pm', arr.pm, 'muStar', muStar));
arrival = struct('params', arr, 'timeND', tT, 'stateND', rvT, ...
                 'periodND', infoT.periodND, 'closure', norm(rvT(end,:) - rvT(1,:)));
closureT = arrival.closure;
fprintf('\n2. TARGET %s: %d petals, branch %+d, period %.4f ND (%.3f d) LOCKED by Np\n', ...
        upper(arr.family), arr.Np, arr.pm, infoT.periodND, day(infoT.periodND));
fprintf('   out-of-plane extent max|z| = %.4f ND (%.0f km)\n', ...
        max(abs(rvT(:,3))), max(abs(rvT(:,3)))*lStar);
fprintf('   periodicity closure |x(T) - x(0)| = %.2e   %s\n', closureT, ...
        pass(closureT < 1e-7));

%% ========================================================================
%  3. ENGINE, PHASES, ENDPOINTS
%     A phase is a fraction of its orbit's period. The target state is read
%     from a PERIODIC SPLINE through the propagated orbit -- that interpolant,
%     not the raw table and not the dynamics, is what the boundary condition
%     matches.
%% ========================================================================
thrustN = 0.070;      ispS = 900;      m0kg = 150;
sD      = 0.0;                                % departure phase, fraction of tau_D
sA      = 0.0754;                             % arrival phase,  fraction of tau_A

g0  = 9.80665*tStar^2/(1000*lStar);           % ND gravity at sea level
cnd = (ispS/tStar)*g0;                        % ND exhaust speed
Tnd = (thrustN/m0kg)*tStar^2/(lStar*1000);    % ND thrust acceleration at m = 1

% PERIODIC interpolant. `spline` builds a not-a-knot cubic; wrapping the
% argument in mod() does not make it periodic, so its derivative jumps across
% the seam at s = 0 -- exactly where a phase near zero is evaluated.
[ppD, kindD] = periodicPP(tD, rvD);
[ppA, kindA] = periodicPP(tT, rvT);
stateD = @(s) ppval(ppD, mod(s,1)*tD(end));
stateA = @(s) ppval(ppA, mod(s,1)*tT(end));
departure.stateAtPhase = stateD;   arrival.stateAtPhase = stateA;
rv0 = stateD(sD);   rvf = stateA(sA);
departure.endpointND = rv0;        arrival.endpointND = rvf;
seamD = norm(ppval(ppD, 0) - ppval(ppD, tD(end)));
seamA = norm(ppval(ppA, 0) - ppval(ppA, tT(end)));

fprintf('\n3. ENGINE %.0f mN, Isp %g s, m0 %g kg  ->  T_nd = %.6e, c_nd = %.4f\n', ...
        thrustN*1000, ispS, m0kg, Tnd, cnd);
fprintf('   departure phase %.4f -> r = [%+.5f %+.5f %+.5f]\n', sD, rv0(1:3));
fprintf('   arrival   phase %.4f -> r = [%+.5f %+.5f %+.5f]\n', sA, rvf(1:3));
fprintf('   interpolant: %s (departure), %s (target); seam mismatch %.1e / %.1e\n', ...
        kindD, kindA, seamD, seamA);

%% ========================================================================
%  4. SOLVE -- seed, then multiple shooting on the PMP boundary-value problem
%     Unknowns: the seven initial costates and the final time. The seed comes
%     from the certified library when this operating point is in it; a new
%     orbit pair needs a continuation from a solved one (run_dro_tulip).
%% ========================================================================
% The seed library was built at ONE operating point. Check the WHOLE of it --
% period, petal count, branch, thrust, Isp and mass -- not just two of the six,
% or a seed from a different engine would silently start the solve.
lib = dro_tulip_library();
match = find(abs([lib.sD] - mod(sD,1)) < 1e-6 & abs([lib.sA] - mod(sA,1)) < 1e-6, 1);
libOp = struct('tau', 1.0, 'Np', 7, 'pm', -1, 'thrustN', 0.070, 'ispS', 900, 'm0kg', 150);
sameOp = abs(dep.tau - libOp.tau) < 1e-12 && arr.Np == libOp.Np && arr.pm == libOp.pm && ...
         abs(thrustN - libOp.thrustN) < 1e-12 && abs(ispS - libOp.ispS) < 1e-9 && ...
         abs(m0kg - libOp.m0kg) < 1e-9;
assert(~isempty(match) && sameOp, ...
    ['no seed for this operating point. The library covers the tau = 1 DRO to the\n' ...
     '7-petal tulip (pm = -1) at 70 mN, Isp 900 s, 150 kg. For any other case,\n' ...
     'walk to it with   T = run_dro_tulip(sD, sA, opts)   and study T here.']);

K = size(lib(match).Y, 2);
seed = struct('tf', lib(match).z(8), 'tGrid', linspace(0, lib(match).z(8), K+1), ...
              'Y', [lib(match).Y, lib(match).Y(:,end)]);
seed.Y(1:7,1) = [rv0(1:6); 1];   seed.Y(8:14,1) = lib(match).z(1:7);

[z8, it] = ms_tfmin(rv0(1:6), rvf(1:6), seed, Tnd, cnd, muStar, ...
                    struct('tolR', 3e-11, 'wallSec', 600, 'conjTest', true));

% ENFORCE convergence before anything downstream. A best-iterate from an
% unsuccessful solve would otherwise be verified, plotted and described in
% the language of optimality.
assert(it.converged && isfinite(it.normR) && it.normR < 1e-8, ...
       'the shooting solve did not converge (|R| = %.2e): nothing below is meaningful', it.normR);

% ONE FLIGHT. Every number in sections 6-8 comes from this object. Flying
% again elsewhere with the same routine is repetition, not independence, and
% it leaves the reader asking which trajectory owns the reported numbers.
% (Astra script review 2026-09-10.)
[tu, Y] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(1:6); 1; z8(1:7)], Tnd, cnd, muStar);
flight = struct('t', tu, 'Y', Y, 'z8', z8, 'rv0', rv0(1:6), 'rvf', rvf(1:6), ...
                'Tnd', Tnd, 'cnd', cnd, 'muStar', muStar, 'tStar', tStar, ...
                'lStar', lStar, 'm0kg', m0kg, 'nSamples', numel(tu));
% ADMISSIBILITY, before any quantity that assumes it: positive finite mass
% (log(1/mf) and every division by m), and no approach to a primary.
assert(all(isfinite(Y(:))), 'the flight returned non-finite states');
assert(min(Y(:,7)) > 0, 'mass reached %.3g: the arc leaves the model''s domain', min(Y(:,7)));
altMoon = min(vecnorm(Y(:,1:3) - [1-muStar 0 0], 2, 2))*lStar - 1737;
altEarth = min(vecnorm(Y(:,1:3) - [-muStar 0 0], 2, 2))*lStar - 6378;
mf   = Y(end,7);
dV   = cnd*log(1/mf)*lStar/tStar;
fprintf('\n4. SOLVED: %d Newton iterations, |R| = %.2e, converged = %d\n', ...
        it.iters, it.normR, it.converged);
fprintf('   t_f = %.6f ND (%.4f d)   Delta-V = %.4f km/s   propellant %.2f kg\n', ...
        z8(8), day(z8(8)), dV, m0kg*(1-mf));
fprintf('   lambda_0 = [%s]\n', strjoin(compose('%+.6g', z8(1:7)'), ' '));
fprintf('   admissible: min mass %.4f, closest approach %.0f km (Moon) / %.0f km (Earth)\n', ...
        min(Y(:,7)), altMoon, altEarth);

%% ========================================================================
%  5. INDEPENDENT VERIFICATION -- a second solver must not move the costates
%% ========================================================================
B = struct('problem', struct('lStar', lStar, 'tStar', tStar, 'muStar', muStar, ...
                             'thrustN', thrustN, 'ispS', ispS, 'm0kg', m0kg, ...
                             'tauDRO', dep.tau, 'NpTulip', arr.Np, 'pmTulip', arr.pm), ...
           'Tnd', Tnd, 'cnd', cnd, 'mu', muStar, 'stateD', stateD, 'stateA', stateA);
V = verify_with_pumpkyn(struct('z', z8, 'sD', sD, 'sA', sA), B);

%% ========================================================================
%  6. NECESSARY CONDITIONS (Pontryagin, first order) -- one at a time
%% ========================================================================
% THRESHOLDS ONCE, used by both the printed lines and the verdict. Written
% twice they drift, and the report then disagrees with its own conclusion.
tol = struct('R', 1e-8, ...       % BVP residual, inf-norm
             'H', 1e-6, ...       % |H| along the arc (exactly 0 in theory)
             'km', 100, ...       % flown arrival, position
             'ms', 10, ...        % flown arrival, velocity
             'lamm', 1e-6, ...    % transversality lam_m(t_f) = 0
             'adj', 1e-7, ...     % adjoint equations, relative
             'min', 1e-12, ...    % H minimality slack over the sphere
             'dz', 1e-6);         % cross-check, second solver
fprintf('\n6. NECESSARY CONDITIONS      (value / threshold)\n');

% N1  the boundary-value residual: costate equations, terminal matching,
%     transversality. This IS the statement that the first variation vanishes.
[~, chk] = ms_tfmin(rv0(1:6), rvf(1:6), seed, Tnd, cnd, muStar, struct('assembleOnly', true));
R1 = chk.residual([z8(1:7); reshape(it.Y(:,2:end), [], 1); z8(8)]);
resid = norm(R1, inf);
fprintf('   N1 BVP residual   |R|_inf   %9.2e / %-9.0e  %s\n', resid, tol.R, pass(resid < tol.R));

% N2  the Hamiltonian. Autonomous problem, free final time  =>  H == 0.
%     H = 1 + lam_r.v + lam_v.(g + (T/m) alpha) - lam_m T/c,  alpha = -lam_v/|lam_v|
lamR = Y(:, 8:10);   lamV = Y(:, 11:13);   lamM = Y(:, 14);   mass = Y(:, 7);
lamVmag = vecnorm(lamV, 2, 2);      % |lambda_v|, the Legendre quantity
Hval = zeros(size(tu));
for k = 1:numel(tu)
    F = mintime_rhs_point(Y(k,:).', Tnd, cnd, muStar);     % [xdot; lamdot]
    Hval(k) = 1 + Y(k,8:14)*F(1:7);
end
Hmax = max(abs(Hval));
fprintf('   N2 Hamiltonian    max|H|    %9.2e / %-9.0e  %s\n', Hmax, tol.H, pass(Hmax < tol.H));

% N3  the flight actually reaches the target, in position AND velocity
missKm  = norm(Y(end,1:3) - rvf(1:3)')*lStar;
missVms = norm(Y(end,4:6) - rvf(4:6)')*lStar/tStar*1000;
fprintf('   N3 arrival        %.4f km / %-4.0f km, %.4f m/s / %-3.0f m/s  %s\n', ...
        missKm, tol.km, missVms, tol.ms, pass(missKm < tol.km && missVms < tol.ms));

% N4  transversality on the free mass: lam_m(t_f) = 0
fprintf('   N4 transversality |lam_m(t_f)| %6.2e / %-9.0e  %s\n', ...
        abs(Y(end,14)), tol.lamm, pass(abs(Y(end,14)) < tol.lamm));

% N5  the ADJOINT equations themselves: lambda-dot = -dH/dx. A small shooting
%     residual says the pieces MATCH each other; it does not say the costate
%     equations are the right ones. This tests them directly.
%     Done by differentiating H in the STATE at fixed costate and comparing
%     with the costate rate the propagator uses. Differencing the propagator's
%     OUTPUT instead would measure its sample spacing: on this arc that gave
%     5e-4, all of it truncation near the 4,674 km lunar passage.
kk = round(linspace(2, numel(tu)-1, 120));
hFD = 1e-6;  adjErr = 0;
for k = kk
    yk = Y(k,:).';
    F  = mintime_rhs_point(yk, Tnd, cnd, muStar);
    dHdx = zeros(7,1);
    for jj = 1:7
        yp = yk;  yp(jj) = yp(jj) + hFD;
        ym = yk;  ym(jj) = ym(jj) - hFD;
        Fp = mintime_rhs_point(yp, Tnd, cnd, muStar);
        Fm = mintime_rhs_point(ym, Tnd, cnd, muStar);
        dHdx(jj) = (yk(8:14).'*Fp(1:7) - yk(8:14).'*Fm(1:7))/(2*hFD);
    end
    adjErr = max(adjErr, norm(F(8:14) + dHdx)/max(norm(dHdx), 1));
end
fprintf('   N5 adjoint eqns   rel err   %9.2e / %-9.0e  %s\n', adjErr, tol.adj, pass(adjErr < tol.adj));

% N6  the MINIMUM principle itself: H must be MINIMISED over the admissible
%     control set, not merely stationary. H is affine in the direction, so
%     compare the flown alpha against a dense sample of the unit sphere.
alpha = -lamV ./ max(lamVmag, realmin);
nS = 400;  rng(0);
Asamp = randn(nS, 3);  Asamp = Asamp ./ vecnorm(Asamp, 2, 2);
worst = 0;
for k = kk
    hStar  = (Tnd/mass(k)) * (lamV(k,:) * alpha(k,:).');        % the flown term
    hOther = (Tnd/mass(k)) * (Asamp * lamV(k,:).');             % every sampled one
    worst  = max(worst, hStar - min(hOther));                  % must be <= 0
end
fprintf('   N6 H minimised    slack     %9.2e / %-9.0e  %s\n', worst, tol.min, pass(worst <= tol.min));

necessary = resid < tol.R && Hmax < tol.H && missKm < tol.km && missVms < tol.ms && ...
            abs(Y(end,14)) < tol.lamm && adjErr < tol.adj && worst <= tol.min;

% The independent solve is NOT one of the above. It is a cross-check on our
% implementation, not a condition of the theory, so it is reported apart and
% excluded from `necessary`. (Astra script review 2026-09-10: the header
% already said this and the code counted it anyway.)
fprintf('   -- cross-check (NOT a PMP condition): second solver |dz| %6.2e / %-9.0e  %s\n', ...
        V.dz, tol.dz, pass(~V.moved));

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
        minLamV, tu(iMin)/tu(end), pass(minLamV > 0));

% S2  all-burn really is the PMP control: the switching function stays positive
Qmt = lamVmag./mass + lamM/cnd;
fprintf('   S2 switching     min Q = %.4e                        %s\n', ...
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
%     ENDPOINT, which is UNRESOLVED rather than either. And a verdict is only
%     interpretable beside its COVERAGE -- how many samples, over what span.
cj = it.conj;
s4Status = 'FAIL';
if cj.pass == 1
    s4Status = 'PASS';
elseif strcmpi(cj.verdict, 'ENDPOINT')
    s4Status = 'UNRESOLVED';
end
fprintf('   S4 conjugate     %s, %d crossing(s), min|det| %.2e       %s\n', ...
        cj.verdict, gv(cj, 'nCrossings'), min(abs(cj.detScaled)), s4Status);
fprintf('      coverage: %d samples on (0, t_f], first full-rank at %d; a sign\n', ...
        numel(cj.detScaled), gv(cj, 'firstFullRank'));
fprintf('      test cannot see an even-order zero or two zeros in one segment\n');

% S5  H6: exclude the REDUCED problem's spurious-zero mechanism. Not a
%     hypothesis of the theorem -- a hazard of OUR reduction (FINDINGS 40).
if isfield(gates, 'h6Margin')
    fprintf('   S5 H6 spurious   lambda_m(0) %.3f vs c/T %.3f, margin %5.1fx  %s\n', ...
            gates.h6LamM0, gates.h6Threshold, gates.h6Margin, pass(gates.h6Ok));
end

h6ok = ~isfield(gates, 'h6Ok') || gates.h6Ok;
sufficient = minLamV > 0 && min(Qmt) > 0 && gates.dimS == 1 && ...
             strcmp(s4Status, 'PASS') && h6ok;

% CONSISTENCY FIRST, verdict second. Exposing the scaffolding risks growing a
% second, unverified implementation of the tests; this binds the inline
% numbers to the shared instrument -- and it has to run BEFORE any verdict is
% printed, or a mismatch is announced too late to matter.
assert(abs(minLamV - gates.minLamV) < 1e-6*max(gates.minLamV,1) && ...
       abs(min(Qmt) - gates.minQmt) < 1e-6*max(gates.minQmt,1), ...
       'the inline gates disagree with mintime_hypothesis_gates');
fprintf('   (inline gates agree with mintime_hypothesis_gates)\n');

fprintf('\n   VERDICT: ');
if necessary && sufficient
    % Deliberately NOT "certified at the sampled times": strong local
    % minimality is a property of the whole trajectory, so attaching the
    % sampling to the CLAIM is ill-posed. The sampling qualifies the
    % EVIDENCE, which is what the second line says. (Astra 2026-09-10.)
    fprintf(['every hypothesis of the BCT sufficiency theorem was checked and holds.\n' ...
             '            If they hold exactly, this arc is a strict strong local minimizer\n' ...
             '            among trajectories with the same endpoints and phases.\n']);
    fprintf(['            EVIDENCE IS NUMERICAL AND SAMPLED: positivity is tested at the\n' ...
             '            sampled times, dim S is a numerical rank, and the conjugate test\n' ...
             '            is a sign test at %d junctions -- an even-order zero, or two zeros\n' ...
             '            inside one segment, would not be seen. This is a strong numerical\n' ...
             '            audit, not a proof. See doc/mintime_second_order_audit.tex.\n'], ...
             size(it.Y, 2));
elseif ~necessary
    fprintf(['NOT an extremal to tolerance. The second-order test is meaningless\n' ...
             '            off an extremal, so no minimality is claimed.\n']);
else
    fprintf('an extremal, but a sufficiency hypothesis fails: no minimality claimed.\n');
end


%% ========================================================================
%  8. INTERACTIVE 3D PLOT  (drag to rotate, scroll to zoom)
%% ========================================================================
T = struct('z', z8, 'sD', sD, 'sA', sA, 'tfDays', day(z8(8)), 'dvKms', dV, ...
           'propellantKg', m0kg*(1-mf), 'finalMassKg', m0kg*mf);
P = plot_transfer_3d(T, B);
fprintf('\n8. Figure %d is rotatable.\n', P.fig.Number);

%% ------------------------------------------------------------------------
function [pp, kind] = periodicPP(tt, yy)
% PERIODICPP  Piecewise-polynomial interpolant of one period of an orbit,
% PERIODIC where the Curve Fitting Toolbox allows. An ordinary not-a-knot
% spline is not C1 across the seam at s = 0, and a phase near zero is
% evaluated exactly there.  INPUTS: tt [N x 1]; yy [N x 6].
% OUTPUTS: pp; kind (which one was built).
tt = tt(:).';  Y = yy.';
if exist('csape', 'file') == 2
    try
        pp = csape(tt, Y, 'periodic');  kind = 'periodic cubic';  return
    catch
    end
end
pp = spline(tt, Y);  kind = 'not-a-knot cubic (NOT periodic)';
end

function s = pass(c)
% PASS  Verdict text.  INPUTS: c.  OUTPUTS: s.
if c, s = 'PASS'; else, s = 'FAIL'; end
end

function v = gv(s, f)
% GV  Field or NaN.  INPUTS: s; f.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = NaN; end
end
