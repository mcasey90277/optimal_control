function G = leg1_gates(L, P, ctx, opts, varargin)
%% Purpose:
%
%   THE GATES OF ONE LEG-1 TRANSFER (GTO -> DRO, 5 N, minimum time):
%   given one polished transfer P (polish_gto_dro) and the leg-1 context
%   ctx, decide whether it is acceptable, report its second-order
%   verdicts, and compute its mass and power bookkeeping. Everything is
%   measured on an INDEPENDENT re-flight of P.z8 from P.rv0
%   (fly_transfer); P.flight is not used.
%
%   arrival  flown miss against P.rvf: missKm < flyKmMax [100] and
%            missMps < flyVmsMax [10] (the leg-2 inner rule).
%   oneGto   the departure state recomputed from ctx and P.tDep ALONE:
%            orientation ctx.orient0 + L.sgn*(180/pi)*P.tDep deg, phase
%            mod(ctx.gtoPhase0 + P.tDep/L.tauG, 1), state L.gtoState of
%            those; max |difference| with P.rv0 < oneGtoTolNd [1e-9 ND].
%            P.orientDeg and P.sD are NOT used.
%   flight   F.admissibility of the re-flight (fly_transfer's
%            validate_flight: reached t_f, mass law, Moon/Earth centre
%            distances); G.flight.ok is required by G.pass.
%   handoff  (Review Focus 5) the DRO phase of the FLOWN arrival state
%            (the last state of the re-flight, not P.rvf): nearest point
%            of the DRO in position (dense search over nDro [4000] phases
%            of L.slotState, then fminbnd on the bracketing cell), its
%            distance distKm < handoffKmMax [100]; advanced from the
%            FLOWN arrival epoch tArrFlown = P.tDep + F.t(end) by
%            (ctx.T0 - tArrFlown)/L.tauD and compared with (ctx.kSat-1)/6
%            by a wrapped phase difference errPhase in [-1/2, 1/2):
%            |errPhase| < handoffTolPhase [1e-5 of a period, about 4.4 s
%            or 1.9 km along track on the DRO], and tArrFlown must agree
%            with P.tArr to 1e-9 ND. Independent of the schedule's
%            bookkeeping. TOLERANCE (controller ruling 2026-11-11): the
%            gate exists to catch clock and numbering errors, whose scale
%            is 1/6 of a period; the arrival gate already bounds the miss.
%            Measured on the test transfer: errPhase 7.26e-7, which is the
%            along-track part (0.139 km at 0.50 km/s) of its 0.144 km
%            arrival miss, so 1e-6 would have rejected good transfers
%            with misses above about 0.2 km.
%   earth,   (Review Focus 2) the LOCATED minimum altitude over the Earth
%   moon     (radius 6378.137 km, floor earthFloorKm [200]) and the Moon
%            (radius 1737.4 km, floor moonFloorKm [500]); see LOCATED_MIN.
%            The sample minimum (earth_clearance, the GTO_DRO instrument)
%            is reported beside it.
%   second   (reported, NOT a blocker) the leg-2 instruments: a re-solve
%            ms_tfmin(P.rv0, P.rvf, seed_from_z8(P.z8, P.rv0, K [24]),
%            tolR 3e-11, wallSec 600, conjTest) counted as reproduced when
%            converged and max |z - P.z8| < reproTol [1e-6]; conjPass and
%            conjCovered from its info.conj; mintime_hypothesis_gates
%            (relTol 1e-12): minLamV, minQmt, dimS. pass = reproduced &&
%            conjPass == 1 && conjCovered == 1 && minLamV > 0 && minQmt > 0
%            && dimS == 1. second.status is 'pass', 'refuted' (reproduced
%            and a conjugate or hypothesis gate fails), 'notReproduced',
%            'error' (a throw inside, caught) or 'skipped' (doSecond
%            false).
%   mass     mfFrac = final-to-initial mass of the flight (its 7th state,
%            m0 = 1 ND), f = 1 - mfFrac; transfer propellant per m0
%            propKgPer150 = m0*f; stage propellant mp = m0*f/(1 - f -
%            dryFrac*f) [dryFrac 0.15], stage dry d = dryFrac*mp, stack wet
%            m0 + d + mp, thrust L.thrustN*wet/m0 (same acceleration as
%            the 5 N, 150 kg regime), jet power thrust*Isp*g0/2 [kW],
%            dvKms = -Isp*g0*log(mfFrac)/1000, with g0 = 9.80665 m/s^2.
%            GUARD: mfFrac not finite or not in (0, 1], or 1 - f -
%            dryFrac*f <= 0 (the stage cannot exist), or an inadmissible
%            flight: every derived field NaN, feasible = false, reason.
%
%   G.pass = flight.ok && arrival && oneGto && handoff && earth && moon
%            && mass.feasible.
%   G.refuted = second.status is 'refuted' (ONLY then).
%   G.secondUnknown = second.status is 'notReproduced' or 'error'.
%
%   LOCATED_MIN (minimum distance to a primary BETWEEN samples). The
%   samples (t, X) are refined densely by a cubic Hermite interpolant of
%   position with the sampled velocities, nSub [16] points per sample
%   interval; every place where the range rate (r - r_b).v changes sign
%   from negative to non-negative brackets a local minimum, and its root
%   is found by fzero on the range rate of the state PROPAGATED from the
%   left sample of the interval that holds the query time (prop handle;
%   a sample interval of a tfMinProp flight holds no switch), or on the
%   Hermite interpolant when no propagator is given. If the propagated
%   range rate does not change sign on the dense bracket, the bracket
%   widens to [t(k-1), t(k+2)] (clamped to the ends); if still no sign
%   change, fminbnd on the distance there. hMinKm is the smallest
%   altitude over all located minima and the two end points. A located
%   minimum that cannot be computed (NaN state: a re-integration that
%   stopped before its end time, or a throw) FAILS the gate with a
%   reason; it is never dropped. Callable on its own:
%
%     C = leg1_gates('located_min', t, X, prop, p)
%       t [M x 1] ND; X [M x n >= 6] ND rotating-frame states (rows);
%       prop handle (tA, xA [n x 1], tB) -> [n x 1] or [] for Hermite
%       only; p .center [1 x 3] ND .radiusKm .lStar .muStar .floorKm
%       (.nSub [16]). C: .hMinKm .tAtMin .pass .reason .sampleMinKm
%       .sampleTAtMin .samplePass .minima [nMin x 2] [t hKm] .nBrackets
%       .nFailed .method.
%
%   TEST SEAMS (opts.inject, empty by default): .truncateFlight q keeps
%   the flight samples with t <= q*t(end) and re-validates them with
%   validate_flight; .propStopFrac q makes every fine re-integration stop
%   at tA + q*(tB - tA) (its end-time check must catch it);
%   .secondError true throws inside the second-order block.
%
%     Ms = leg1_gates('mass', mfFrac, L, opts)   (opts.dryFrac [0.15])
%
%   Self-demo (nargin == 0): the test_polish_gto_dro case (2) transfer
%   (cell (1,1), +10 deg, slot +0.03, tDep 0.25) with a consistent ctx,
%   verbose.
%
%% Inputs:
%
%  L                        struct                  leg1_inputs(struct())
%
%  P                        struct                  polish_gto_dro output
%                                                   with P.ok (z8 [8 x 1],
%                                                   rv0, rvf [6 x 1], tDep,
%                                                   tArr [ND])
%
%  ctx                      struct                  .orient0 GTO orientation
%                                                   at leg-1 time 0 [deg]
%                                                   .gtoPhase0 GTO phase at
%                                                   leg-1 time 0 [fraction]
%                                                   .T0 hand-off epoch,
%                                                   leg-1 clock [ND]
%                                                   .kSat leg-2 satellite
%                                                   index 1..6
%
%  opts                     struct (optional)       .flyKmMax [100]
%                                                   .flyVmsMax [10]
%                                                   .oneGtoTolNd [1e-9]
%                                                   .handoffTolPhase [1e-5]
%                                                   .handoffKmMax [100]
%                                                   .earthFloorKm [200]
%                                                   .moonFloorKm [500]
%                                                   .dryFrac [0.15]
%                                                   .doSecond [true]
%                                                   .K [24] .reproTol [1e-6]
%                                                   .nSub [16] .nDro [4000]
%                                                   .inject [struct()] test
%                                                   seams (see above)
%                                                   .verbose [false]
%                                                   unknown names throw
%                                                   leg1_gates:opts
%
%% Outputs:
%
%  G                        struct                  .flight (validate_flight
%                                                   output: .ok .reason ...)
%                                                   .arrival (.missKm
%                                                   .missMps .pass)
%                                                   .oneGto (.dStateNd
%                                                   .orientDeg .phase .pass)
%                                                   .handoff (.phaseArr
%                                                   .distKm .dVelMps
%                                                   .tfFlown .tArrFlown
%                                                   .dtArr .phaseAtT0
%                                                   .target .errPhase
%                                                   .reason .pass)
%                                                   .earth, .moon (.hMinKm
%                                                   .tAtMin .sampleMinKm
%                                                   .sampleTAtMin .floorKm
%                                                   .minima .nBrackets
%                                                   .nFailed .reason .pass)
%                                                   .second (.status
%                                                   .conjPass .conjCovered
%                                                   .conjVerdict .minLamV
%                                                   .minQmt .dimS
%                                                   .reproduced .dz .resid
%                                                   .reason .skipped .wallS
%                                                   .pass)
%                                                   .mass (.mfFrac
%                                                   .propKgPer150
%                                                   .stagePropKg .stageDryKg
%                                                   .wetKg .thrustN .jetKW
%                                                   .dvKms .dryFrac
%                                                   .feasible .reason)
%                                                   .pass .refuted
%                                                   .secondUnknown .opts
%                                                   .wallS
%
%% Revision History:
%  M. Casey                                                   (c) 11/11/2026
%  M. Casey  fix round 1: flight admissibility, mass guard,   (c) 11/11/2026
%            located-min failures, flown arrival epoch,
%            second.status, handoff tolerance 1e-5
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    run(fullfile(here, 'setup_lunar_pnt_paths.m'));
    L = leg1_inputs(struct());
    cand = struct('iD', 1, 'iA', 1, 'iS', 1, 'tDep', 0.25, ...
                  'orientDeg', L.sheetOrient(1) + 10, 'slotPhaseAtDep', ...
                  mod(L.sA(1) + 0.03 - L.tf(1, 1, 1)/L.tauD, 1));
    P = polish_gto_dro(L, cand, struct());
    ctx.orient0   = P.orientDeg - L.sgn*(180/pi)*P.tDep;
    ctx.gtoPhase0 = mod(P.sD - P.tDep/L.tauG, 1);
    ctx.kSat      = floor(6*P.arrPhase) + 2;
    ctx.T0        = P.tArr + mod((ctx.kSat - 1)/6 - P.arrPhase, 1)*L.tauD;
    G = leg1_gates(L, P, ctx, struct('verbose', true));
    return
end
if ischar(L) || isstring(L)
    switch char(L)
        case 'located_min'
            G = located_min(P, ctx, opts, varargin{1});
        case 'mass'
            if nargin < 4, opts = struct(); end
            G = mass_block(P, ctx, fdef(opts, 'dryFrac', 0.15));
        otherwise
            error('leg1_gates:mode', 'unknown mode "%s"', char(L));
    end
    return
end
if nargin < 4, opts = struct(); end
o = knobs(opts);
tStart = tic;
if ~(isfield(P, 'ok') && P.ok && numel(P.z8) == 8)
    error('leg1_gates:notPolished', 'P is not a polished transfer (ok false or z8 empty)');
end

%% Re-flight of z8 from rv0 (and its admissibility):
F = fly_transfer(P.z8, P.rv0, P.rvf, L.phys);
if isfield(o.inject, 'truncateFlight')
    keep = F.t <= o.inject.truncateFlight*F.t(end);
    F.t = F.t(keep);  F.Y = F.Y(keep, :);
    F.admissibility = validate_flight(F.t, F.Y, P.z8(8), L.phys.Tnd, L.phys.cnd, ...
                                      L.muStar, L.lStar);
    F.flyKm  = L.lStar*sqrt(sum((F.Y(end, 1:3).' - P.rvf(1:3)).^2));
    F.flyVms = 1000*L.lStar/L.tStar*sqrt(sum((F.Y(end, 4:6).' - P.rvf(4:6)).^2));
end
G.flight = F.admissibility;
rvArr = F.Y(end, 1:6).';

%% Arrival:
G.arrival = struct('missKm', F.flyKm, 'missMps', F.flyVms, ...
                   'pass', F.flyKm < o.flyKmMax && F.flyVms < o.flyVmsMax);

%% One GTO (from ctx and tDep only):
orientDeg = ctx.orient0 + L.sgn*(180/pi)*P.tDep;
gPhase    = mod(ctx.gtoPhase0 + P.tDep/L.tauG, 1);
dState    = max(abs(L.gtoState(orientDeg, gPhase) - P.rv0(:)));
G.oneGto  = struct('dStateNd', dState, 'orientDeg', orientDeg, 'phase', gPhase, ...
                   'pass', dState < o.oneGtoTolNd);

%% Hand-off (DRO phase of the flown arrival state, from the flown epoch):
[sArr, dKm, dVms] = dro_phase_of(rvArr, L, o.nDro);
tfFlown   = F.t(end);
tArrFlown = P.tDep + tfFlown;
dtArr     = tArrFlown - P.tArr;
phaseAtT0 = mod(sArr + (ctx.T0 - tArrFlown)/L.tauD, 1);
target    = (ctx.kSat - 1)/6;
errPhase  = mod(phaseAtT0 - target + 0.5, 1) - 0.5;
why = {};
if ~(abs(errPhase) < o.handoffTolPhase)
    why{end+1} = sprintf('phase error %+.3g (tol %.3g)', errPhase, o.handoffTolPhase);
end
if ~(dKm < o.handoffKmMax)
    why{end+1} = sprintf('arrival %.3g km off the DRO (max %.3g)', dKm, o.handoffKmMax);
end
if ~(abs(dtArr) < 1e-9)
    why{end+1} = sprintf('flown arrival epoch differs from P.tArr by %.3g ND', dtArr);
end
G.handoff = struct('phaseArr', sArr, 'distKm', dKm, 'dVelMps', dVms, ...
                   'tfFlown', tfFlown, 'tArrFlown', tArrFlown, 'dtArr', dtArr, ...
                   'phaseAtT0', phaseAtT0, 'target', target, 'errPhase', errPhase, ...
                   'reason', strjoin(why, '; '), 'pass', isempty(why));

%% Clearance on the located minimum:
stopFrac = fdef(o.inject, 'propStopFrac', 1);
prop = @(tA, xA, tB) fly14(tA, xA, tB, L.phys, stopFrac);
pE = struct('center', [-L.muStar 0 0], 'radiusKm', 6378.137, 'lStar', L.lStar, ...
            'muStar', L.muStar, 'floorKm', o.earthFloorKm, 'nSub', o.nSub);
pM = pE;  pM.center = [1 - L.muStar 0 0];  pM.radiusKm = 1737.4;  pM.floorKm = o.moonFloorKm;
G.earth = clearance_out(located_min(F.t, F.Y, prop, pE));
G.moon  = clearance_out(located_min(F.t, F.Y, prop, pM));

%% Second order (reported):
G.second = second_order(L, P, o);

%% Mass and power (NaN / infeasible on an inadmissible flight):
if G.flight.ok
    G.mass = mass_block(F.Y(end, 7), L, o.dryFrac);
else
    G.mass = mass_block(NaN, L, o.dryFrac);
    G.mass.reason = ['flight not admissible: ' char(string(G.flight.reason))];
end

G.pass    = G.flight.ok && G.arrival.pass && G.oneGto.pass && G.handoff.pass && ...
            G.earth.pass && G.moon.pass && G.mass.feasible;
G.refuted = strcmp(G.second.status, 'refuted');
G.secondUnknown = any(strcmp(G.second.status, {'notReproduced', 'error'}));
G.opts    = o;
G.wallS   = toc(tStart);
if o.verbose
    print_gates(G);
end
end

% ------------------------------------------------------------------------
function C = clearance_out(C)
% CLEARANCE_OUT  The per-primary gate fields.  INPUTS: C located_min
%   output.  OUTPUTS: C with the field order of G.earth / G.moon.
C = struct('hMinKm', C.hMinKm, 'tAtMin', C.tAtMin, 'sampleMinKm', C.sampleMinKm, ...
           'sampleTAtMin', C.sampleTAtMin, 'floorKm', C.floorKm, ...
           'minima', C.minima, 'nBrackets', C.nBrackets, 'nFailed', C.nFailed, ...
           'reason', C.reason, 'pass', C.pass);
end

% ------------------------------------------------------------------------
function C = located_min(t, X, prop, p)
% LOCATED_MIN  Minimum altitude over a primary between samples (see the
%   header). INPUTS: t [M x 1]; X [M x n >= 6]; prop handle or [];
%   p struct. OUTPUTS: C struct.
nSub = fdef(p, 'nSub', 16);
t = t(:);
M = numel(t);
pc = struct('muStar', p.muStar, 'lStar', p.lStar, 'center', p.center, ...
            'radiusKm', p.radiusKm);
[hS, tS, sPass, sOut] = earth_clearance(t, X, p.floorKm, pc);

% dense Hermite refinement of the range rate
s   = (0:nSub - 1)/nSub;
h00 = 2*s.^3 - 3*s.^2 + 1;   h10 = s.^3 - 2*s.^2 + s;
h01 = -2*s.^3 + 3*s.^2;      h11 = s.^3 - s.^2;
d00 = 6*s.^2 - 6*s;          d10 = 3*s.^2 - 4*s + 1;
d01 = -6*s.^2 + 6*s;         d11 = 3*s.^2 - 2*s;
nD  = (M - 1)*nSub + 1;
tt  = zeros(nD, 1);  rr = zeros(nD, 1);  kInt = zeros(nD, 1);
for k = 1:M - 1
    dt  = t(k+1) - t(k);
    r0  = X(k, 1:3) - p.center;    v0 = X(k, 4:6);
    r1  = X(k+1, 1:3) - p.center;  v1 = X(k+1, 4:6);
    pos = h00.'*r0 + (dt*h10.')*v0 + h01.'*r1 + (dt*h11.')*v1;
    vel = (d00.'*r0 + (dt*d10.')*v0 + d01.'*r1 + (dt*d11.')*v1)/dt;
    idx = (k - 1)*nSub + (1:nSub);
    tt(idx)   = t(k) + s.'*dt;
    rr(idx)   = sum(pos.*vel, 2);
    kInt(idx) = k;
end
tt(nD)   = t(M);
rr(nD)   = (X(M, 1:3) - p.center)*X(M, 4:6).';
kInt(nD) = M - 1;

br = find(rr(1:end-1) < 0 & rr(2:end) >= 0);
mins = nan(numel(br), 2);
why  = cell(1, numel(br));
for q = 1:numel(br)
    k = kInt(br(q));
    a = tt(br(q));  b = tt(br(q) + 1);
    try
        if isempty(prop)
            g = @(tq) hermite_rr(tq, t(k), t(k+1), X(k, :), X(k+1, :), p.center);
            tR = fzero(g, [a b], optimset('TolX', 1e-15));
            rv = hermite_state(tR, t(k), t(k+1), X(k, :), X(k+1, :));
        else
            g  = @(tq) rr_of(piece_state(tq, t, X, prop), p.center);
            ga = g(a);  gb = g(b);
            if ~(ga < 0 && gb >= 0) && all(isfinite([ga gb]))
                % Hermite bracket off: widen to [t(k-1), t(k+2)], clamped
                a = t(max(k - 1, 1));  b = t(min(k + 2, M));  ga = g(a);  gb = g(b);
            end
            if ~all(isfinite([ga gb]))
                error('leg1_gates:prop', 're-integration stopped early or non-finite near t = %.6g', a);
            end
            if sign(ga) ~= sign(gb)
                tR = fzero(g, [a b], optimset('TolX', 1e-15));
            else                                 % no sign change in the exact flight
                tR = fminbnd(@(tq) dist_of(piece_state(tq, t, X, prop), p.center), a, b, ...
                             optimset('TolX', 1e-15));
            end
            rv = piece_state(tR, t, X, prop);
        end
        hq = p.lStar*sqrt(sum((rv(1:3).' - p.center).^2)) - p.radiusKm;
        if ~isfinite(hq)
            error('leg1_gates:prop', 're-integration stopped early or non-finite at t = %.6g', tR);
        end
        mins(q, :) = [tR, hq];
    catch err
        mins(q, :) = [a, NaN];
        why{q} = sprintf('minimum near t = %.6g not located: %s', a, err.message);
    end
end
cand = [mins; t(1), sOut.hKm(1); t(M), sOut.hKm(M)];
[hMin, kMin] = min(cand(:, 2));
nFailed = nnz(isnan(mins(:, 2)));
why = why(~cellfun(@isempty, why));
pass = hMin >= p.floorKm && nFailed == 0;
if hMin < p.floorKm
    why = [{sprintf('located minimum %.4g km below the %.4g km floor', hMin, p.floorKm)}, why];
end
C = struct('hMinKm', hMin, 'tAtMin', cand(kMin, 1), 'pass', pass, ...
           'reason', strjoin(why, '; '), ...
           'sampleMinKm', hS, 'sampleTAtMin', tS, 'samplePass', sPass, ...
           'floorKm', p.floorKm, 'minima', mins, 'nBrackets', numel(br), ...
           'nFailed', nFailed, 'method', ternary(isempty(prop), 'hermite', 'propagated'));
end

% ------------------------------------------------------------------------
function x = piece_state(tq, t, X, prop)
% PIECE_STATE  State at tq propagated from the left sample of the sample
%   interval that holds tq (the last interval for tq = t(end)).
%   INPUTS: tq ND; t [M x 1]; X [M x n]; prop handle.  OUTPUTS: x [n x 1].
k = find(t <= tq, 1, 'last');
if isempty(k), k = 1; end
k = min(k, numel(t) - 1);
x = prop(t(k), X(k, :).', tq);
end

% ------------------------------------------------------------------------
function g = rr_of(x, center)
% RR_OF  Range rate (r - r_b).v of a state.  INPUTS: x [n >= 6 x 1];
%   center [1 x 3].  OUTPUTS: g scalar.
g = (x(1:3).' - center)*x(4:6);
end

% ------------------------------------------------------------------------
function d = dist_of(x, center)
% DIST_OF  Distance to the primary.  INPUTS: x [n >= 6 x 1]; center
%   [1 x 3].  OUTPUTS: d scalar ND.
d = sqrt(sum((x(1:3).' - center).^2));
end

% ------------------------------------------------------------------------
function x = hermite_state(tq, ta, tb, xa, xb)
% HERMITE_STATE  Cubic Hermite position and its derivative at tq from the
%   end states of one interval.  INPUTS: tq, ta, tb ND; xa, xb [1 x n >= 6].
%   OUTPUTS: x [6 x 1] (position, velocity of the interpolant).
dt = tb - ta;  s = (tq - ta)/dt;
pos = (2*s^3 - 3*s^2 + 1)*xa(1:3) + dt*(s^3 - 2*s^2 + s)*xa(4:6) + ...
      (-2*s^3 + 3*s^2)*xb(1:3) + dt*(s^3 - s^2)*xb(4:6);
vel = ((6*s^2 - 6*s)*xa(1:3) + dt*(3*s^2 - 4*s + 1)*xa(4:6) + ...
       (-6*s^2 + 6*s)*xb(1:3) + dt*(3*s^2 - 2*s)*xb(4:6))/dt;
x = [pos, vel].';
end

% ------------------------------------------------------------------------
function g = hermite_rr(tq, ta, tb, xa, xb, center)
% HERMITE_RR  Range rate of the Hermite interpolant.  INPUTS: as
%   HERMITE_STATE plus center [1 x 3].  OUTPUTS: g scalar.
g = rr_of(hermite_state(tq, ta, tb, xa, xb), center);
end

% ------------------------------------------------------------------------
function x = fly14(tA, xA, tB, phys, stopFrac)
% FLY14  Transfer state (state + costates, 14) at tB from (tA, xA) on
%   pumpkyn's minimum-time field, no event handling (one sample interval
%   of a tfMinProp flight holds no switch). The integration must REACH
%   tB (relative 1e-12 of the span); otherwise the state is NaN, which
%   fails the clearance gate. stopFrac < 1 is the test seam that stops
%   it early.  INPUTS: tA, tB ND; xA [14 x 1]; phys (.Tnd .cnd .muStar);
%   stopFrac [1].  OUTPUTS: x [14 x 1] (NaN if tB not reached).
if tB == tA
    x = xA(:);
    return
end
tE = tA + stopFrac*(tB - tA);
w  = warning('off', 'MATLAB:ode113:IntegrationTolNotMet');
cw = onCleanup(@() warning(w));
[tO, Y] = ode113(@pumpkyn.cr3bp.tfMinEoM, [tA, 0.5*(tA + tE), tE], xA(1:14), ...
                 odeset('RelTol', 1e-12, 'AbsTol', 1e-13), phys.Tnd, phys.cnd, phys.muStar);
if abs(tO(end) - tB) > 1e-12*max(abs(tB - tA), eps) || ~all(isfinite(Y(end, :)))
    x = NaN(14, 1);
else
    x = Y(end, :).';
end
end

% ------------------------------------------------------------------------
function [sArr, dKm, dVms] = dro_phase_of(rv, L, nDro)
% DRO_PHASE_OF  Phase of the nearest point (in position) of the DRO to rv:
%   dense search over nDro phases, then fminbnd on the bracketing cell.
%   INPUTS: rv [6 x 1]; L leg1_inputs; nDro.  OUTPUTS: sArr [fraction];
%   dKm distance to the orbit [km]; dVms velocity difference there [m/s].
sg = (0:nDro - 1)/nDro;
Xd = zeros(6, nDro);
for kk = 1:nDro
    Xd(:, kk) = L.slotState(sg(kk));
end
[~, kMin] = min(sum((Xd(1:3, :) - rv(1:3)).^2, 1));
dfun = @(s) pos_dist2(L.slotState(mod(s, 1)), rv);
sStar = fminbnd(dfun, sg(kMin) - 1/nDro, sg(kMin) + 1/nDro, optimset('TolX', 1e-14));
sArr  = mod(sStar, 1);
xs    = L.slotState(sArr);
dKm   = L.lStar*sqrt(sum((xs(1:3) - rv(1:3)).^2));
dVms  = 1000*L.lStar/L.tStar*sqrt(sum((xs(4:6) - rv(4:6)).^2));
end

% ------------------------------------------------------------------------
function d2 = pos_dist2(x, rv)
% POS_DIST2  Squared position distance.  INPUTS: x, rv [>= 3 x 1].
%   OUTPUTS: d2 scalar ND^2.
d2 = sum((x(1:3) - rv(1:3)).^2);
end

% ------------------------------------------------------------------------
function S = second_order(L, P, o)
% SECOND_ORDER  The leg-2 second-order instruments on one transfer
%   (conjugate re-solve + mintime_hypothesis_gates), reported.
%   INPUTS: L; P; o knobs.  OUTPUTS: S struct (see the header).
S = struct('status', 'skipped', 'conjPass', NaN, 'conjCovered', NaN, 'conjVerdict', '', 'minLamV', NaN, ...
           'minQmt', NaN, 'dimS', NaN, 'reproduced', false, 'dz', NaN, 'resid', NaN, ...
           'reason', 'skipped', 'skipped', ~o.doSecond, 'wallS', 0, 'pass', NaN);
if ~o.doSecond
    return
end
t0 = tic;
Tnd = L.phys.Tnd;  cnd = L.phys.cnd;  mu = L.muStar;
try
    if fdef(o.inject, 'secondError', false)
        error('leg1_gates:inject', 'injected second-order error');
    end
    seed = seed_from_z8(P.z8, P.rv0, o.K, Tnd, cnd, mu);
    [zc, info] = ms_tfmin(P.rv0, P.rvf, seed, Tnd, cnd, mu, ...
                          struct('tolR', 3e-11, 'wallSec', 600, 'conjTest', true));
    S.dz = max(abs(zc - P.z8));
    S.resid = info.normR;
    S.reproduced = info.converged && S.dz < o.reproTol;
    if S.reproduced
        S.conjPass    = double(info.conj.pass);
        S.conjCovered = double(info.conj.covered);
        S.conjVerdict = info.conj.verdict;
        S.reason = sprintf('reproduced; %s: %s', info.conj.verdict, info.conj.reason);
    else
        S.reason = sprintf('not reproduced (converged %d, |dz| = %.1e)', info.converged, S.dz);
    end
    hy = mintime_hypothesis_gates(P.z8, P.rv0, Tnd, cnd, mu, struct('relTol', 1e-12));
    S.minLamV = hy.minLamV;  S.minQmt = hy.minQmt;  S.dimS = hy.dimS;
    S.pass = S.reproduced && S.conjPass == 1 && S.conjCovered == 1 && ...
             S.minLamV > 0 && S.minQmt > 0 && S.dimS == 1;
    if S.pass
        S.status = 'pass';
    elseif S.reproduced
        S.status = 'refuted';
    else
        S.status = 'notReproduced';
    end
catch err
    S.reason = ['error: ' err.message];
    S.pass = false;
    S.status = 'error';
end
S.wallS = toc(t0);
end

% ------------------------------------------------------------------------
function Ms = mass_block(mfFrac, L, dryFrac)
% MASS_BLOCK  Closed-form stage mass and power (see the header), g0 =
%   9.80665 m/s^2. Guarded: mfFrac not finite or not in (0, 1], or a
%   non-positive denominator 1 - f - dryFrac*f, gives NaN fields and
%   feasible = false with a reason.  INPUTS: mfFrac final/initial mass;
%   L (.m0kg .thrustN .ispS); dryFrac.  OUTPUTS: Ms struct.
g0 = 9.80665;
m0 = L.m0kg;
f  = 1 - mfFrac;
den = 1 - f - dryFrac*f;
Ms = struct('mfFrac', mfFrac, 'propKgPer150', NaN, 'stagePropKg', NaN, ...
            'stageDryKg', NaN, 'wetKg', NaN, 'thrustN', NaN, 'jetKW', NaN, ...
            'dvKms', NaN, 'dryFrac', dryFrac, 'feasible', false, 'reason', '');
if ~(isfinite(mfFrac) && mfFrac > 0 && mfFrac <= 1)
    Ms.reason = sprintf('mfFrac %g not in (0, 1]', mfFrac);
    return
end
if ~(den > 0)
    Ms.reason = sprintf(['stage cannot exist: 1 - f - dryFrac*f = %.4g <= 0 ' ...
                         '(f %.4g, dryFrac %.4g)'], den, f, dryFrac);
    return
end
mp  = m0*f/den;
d   = dryFrac*mp;
wet = m0 + d + mp;
thr = L.thrustN*wet/m0;
Ms.propKgPer150 = m0*f;
Ms.stagePropKg  = mp;
Ms.stageDryKg   = d;
Ms.wetKg        = wet;
Ms.thrustN      = thr;
Ms.jetKW        = thr*L.ispS*g0/2/1000;
Ms.dvKms        = -L.ispS*g0*log(mfFrac)/1000;
Ms.feasible     = true;
end

% ------------------------------------------------------------------------
function o = knobs(opts)
% KNOBS  Defaults; unknown names throw leg1_gates:opts.
%   INPUTS: opts struct.  OUTPUTS: o struct.
o = struct( ...
    'flyKmMax',        100, ...
    'flyVmsMax',       10, ...
    'oneGtoTolNd',     1e-9, ...
    'handoffTolPhase', 1e-5, ...
    'handoffKmMax',    100, ...
    'earthFloorKm',    200, ...
    'moonFloorKm',     500, ...
    'dryFrac',         0.15, ...
    'doSecond',        true, ...
    'K',               24, ...
    'reproTol',        1e-6, ...
    'nSub',            16, ...
    'nDro',            4000, ...
    'inject',          struct(), ...
    'verbose',         false);
fn = fieldnames(opts);
for kk = 1:numel(fn)
    if ~isfield(o, fn{kk})
        error('leg1_gates:opts', 'unknown option "%s"', fn{kk});
    end
    o.(fn{kk}) = opts.(fn{kk});
end
end

% ------------------------------------------------------------------------
function print_gates(G)
% PRINT_GATES  One block of every measured number.  INPUTS: G.
%   OUTPUTS: none (prints).
a = G.arrival;  g = G.oneGto;  h = G.handoff;  e = G.earth;  m = G.moon;
s = G.second;   w = G.mass;
fprintf('leg1_gates: pass %d  refuted %d  secondUnknown %d  (%.1f s)\n', G.pass, G.refuted, ...
        G.secondUnknown, G.wallS);
fprintf('  flight   admissible %d  %s\n', G.flight.ok, char(string(G.flight.reason)));
fprintf('  arrival  miss %.4g km %.4g m/s  pass %d\n', a.missKm, a.missMps, a.pass);
fprintf('  oneGto   dState %.3g ND (orient %.6f deg, phase %.9f)  pass %d\n', ...
        g.dStateNd, g.orientDeg, g.phase, g.pass);
fprintf(['  handoff  phaseArr %.9f  dist %.4g km  dVel %.4g m/s  tArrFlown %.12f ' ...
         '(dtArr %.2e)  phaseAtT0 %.9f  target %.9f  err %+.3e  pass %d %s\n'], h.phaseArr, ...
        h.distKm, h.dVelMps, h.tArrFlown, h.dtArr, h.phaseAtT0, h.target, h.errPhase, ...
        h.pass, h.reason);
fprintf('  earth    located %.4f km at t %.6f (sample %.4f km at %.6f; %d minima)  pass %d\n', ...
        e.hMinKm, e.tAtMin, e.sampleMinKm, e.sampleTAtMin, size(e.minima, 1), e.pass);
if ~isempty(e.reason), fprintf('           %s\n', e.reason); end
fprintf('  moon     located %.4f km at t %.6f (sample %.4f km at %.6f; %d minima)  pass %d\n', ...
        m.hMinKm, m.tAtMin, m.sampleMinKm, m.sampleTAtMin, size(m.minima, 1), m.pass);
if ~isempty(m.reason), fprintf('           %s\n', m.reason); end
fprintf(['  second   status %s  conjPass %g covered %g  minLamV %.4g  minQmt %.4g  dimS %g  ' ...
         '|dz| %.2e  resid %.2e  pass %g  (%.1f s)  %s\n'], s.status, s.conjPass, s.conjCovered, ...
        s.minLamV, s.minQmt, s.dimS, s.dz, s.resid, s.pass, s.wallS, s.reason);
fprintf(['  mass     mfFrac %.9f  prop/150 %.4f kg  stage prop %.4f kg  dry %.4f kg  ' ...
         'wet %.4f kg  thrust %.4f N  jet %.4f kW  dV %.5f km/s  feasible %d %s\n'], w.mfFrac, ...
        w.propKgPer150, w.stagePropKg, w.stageDryKg, w.wetKg, w.thrustN, w.jetKW, w.dvKms, ...
        w.feasible, w.reason);
end

% ------------------------------------------------------------------------
function s = ternary(c, a, b)
% TERNARY  a if c else b.  INPUTS: c logical; a, b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
