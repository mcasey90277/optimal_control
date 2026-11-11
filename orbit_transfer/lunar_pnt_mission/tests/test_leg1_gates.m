function ok = test_leg1_gates()
%% Purpose:
%
%   Verifies leg1_gates (leg 1, GTO -> DRO) on the transfer that
%   test_polish_gto_dro case (2) polishes (cell (1,1) of the first sheet,
%   +10 deg, a slot 0.03 of a period past the cell's arrival phase,
%   tDep 0.25 ND), with a ctx built here so that the inertial GTO at tDep
%   has that orientation and phase and the arrival slot is satellite kSat
%   at T0:
%   (1) every gate passes (arrival, oneGto, handoff, earth, moon) and the
%   second-order verdicts are reported (expected: not refuted); the
%   located minima never exceed the sample minima;
%   (2) the mass block: mfFrac is the flight's 7th state at the end, the
%   closed form holds to 1e-12 on the test transfer, and the hand-worked
%   case f = 0.42, dryFrac = 0.15 (exact fractions: mp = 63000/517 kg,
%   d = 9450/517 kg, wet = 150000/517 kg, thrust = 5000/517 N, jet power
%   33538743/413600 kW, dV = -1710*9.80665*log(0.58)/1000 km/s) to 1e-12;
%   (3) deliberate failures, one per gate: arrival (z8 tf perturbed),
%   oneGto (ctx.orient0 + 1 deg), handoff (ctx.T0 + 0.01 ND), Review
%   Focus 5 (kSat shifted by one), Moon (floor 1 km above the measured
%   located minimum fails, 1 km below passes);
%   (4) Review Focus 2: a synthetic ballistic CR3BP perigee pass at
%   150 km altitude, sampled so that both samples adjacent to perigee are
%   above 200 km: the sample minimum passes the 200 km floor, the located
%   minimum fails it and agrees with a fine-grid truth to 0.01 km;
%   (5) an unknown option name is refused;
%   (6, fix round 1) the hand-off default tolerance is 1e-5; the mass
%   guard (mfFrac 0.1, 0, NaN, 1.2 give NaN fields, feasible false) and
%   an infeasible stage (dryFrac 1.5) failing G.pass; a truncated flight
%   (inject.truncateFlight 0.5, re-validated) is not admissible and
%   fails G.pass with a NaN mass block; the flown arrival epoch must
%   agree with P.tArr to 1e-9; a fine re-integration that stops early
%   (inject.propStopFrac 0.5) and a NaN propagator fail the clearance
%   gates with a reason; an injected second-order error is status
%   'error', secondUnknown, NOT refuted.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 true if every assertion
%                                                   held (else it throws)
%
%% Revision History:
%  M. Casey                                                   (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_lunar_pnt_paths.m'));
assert(strcmp(fileparts(which('leg1_gates')), fullfile(fileparts(here), 'lib')), ...
       'leg1_gates not from this campaign');

L  = leg1_inputs(struct());
iS = 1;  iD = 1;  iA = 1;
tfCell = L.tf(iD, iA, iS);
cand = struct('iD', iD, 'iA', iA, 'iS', iS, 'tDep', 0.25, ...
              'orientDeg', L.sheetOrient(iS) + 10, ...
              'slotPhaseAtDep', mod(L.sA(iA) + 0.03 - tfCell/L.tauD, 1));
P = polish_gto_dro(L, cand, struct());
assert(P.ok, 'polish not ok: %s', P.reason);
ctx = good_ctx(L, P);
fprintf('  transfer: tf %.9f ND  tDep %.3f  tArr %.9f  arrPhase %.6f  ctx: orient0 %.6f deg  gtoPhase0 %.6f  T0 %.9f  kSat %d\n', ...
        P.tf, P.tDep, P.tArr, P.arrPhase, ctx.orient0, ctx.gtoPhase0, ctx.T0, ctx.kSat);

% ---- 1. every gate on the good transfer ----
t0 = tic;
G = leg1_gates(L, P, ctx, struct('verbose', true));
fprintf('  leg1_gates (with second order): %.1f s (second order %.1f s)\n', toc(t0), G.second.wallS);
assert(G.arrival.pass, 'arrival: %.3g km %.3g m/s', G.arrival.missKm, G.arrival.missMps);
assert(G.oneGto.pass, 'oneGto: %.3g ND', G.oneGto.dStateNd);
assert(G.handoff.pass, 'handoff: err %.3g, %.3g km', G.handoff.errPhase, G.handoff.distKm);
assert(G.earth.pass && G.moon.pass, 'clearance: earth %.1f moon %.1f', ...
       G.earth.hMinKm, G.moon.hMinKm);
assert(G.pass, 'G.pass false');
assert(G.earth.hMinKm <= G.earth.sampleMinKm + 1e-9 && G.moon.hMinKm <= G.moon.sampleMinKm + 1e-9, ...
       'a located minimum exceeds its sample minimum');
assert(~G.refuted && G.second.pass, 'second order refuted: %s', G.second.reason);
assert(islogical(G.refuted) && strcmp(G.second.status, 'pass') && ~G.secondUnknown, ...
       'second status %s', G.second.status);
assert(G.flight.ok && G.mass.feasible, 'flight admissible %d, mass feasible %d', ...
       G.flight.ok, G.mass.feasible);
assert(G.opts.handoffTolPhase == 1e-5, 'handoff tolerance default %g', G.opts.handoffTolPhase);
assert(abs(G.handoff.dtArr) < 1e-9 && ...
       abs(G.handoff.tArrFlown - (P.tDep + G.handoff.tfFlown)) == 0, 'flown arrival epoch');
fprintf('  flown arrival epoch %.12f vs P.tArr %.12f (%.2e)\n', G.handoff.tArrFlown, P.tArr, ...
        G.handoff.dtArr);

% ---- 2. mass block ----
Ms = G.mass;
F  = fly_transfer(P.z8, P.rv0, P.rvf, L.phys);
assert(abs(Ms.mfFrac - F.Y(end, 7)) < 1e-15 && abs(Ms.mfFrac - F.mf) < 1e-15, 'mfFrac');
f   = 1 - Ms.mfFrac;
mp  = L.m0kg*f/(1 - f - 0.15*f);
cf  = [L.m0kg*f, mp, 0.15*mp, L.m0kg + 1.15*mp, L.thrustN*(L.m0kg + 1.15*mp)/L.m0kg];
got = [Ms.propKgPer150, Ms.stagePropKg, Ms.stageDryKg, Ms.wetKg, Ms.thrustN];
assert(max(abs(got - cf)./abs(cf)) < 1e-12, 'closed form on the test transfer: %.3g', ...
       max(abs(got - cf)./abs(cf)));
assert(abs(Ms.jetKW - Ms.thrustN*L.ispS*9.80665/2/1000)/Ms.jetKW < 1e-12, 'jet power');
assert(abs(Ms.dvKms - F.dvKms)/F.dvKms < 1e-10, 'dV vs fly_transfer: %.3g', ...
       abs(Ms.dvKms - F.dvKms)/F.dvKms);
fprintf(['  mass (test transfer): mfFrac %.12f  f %.6g  prop/150 %.6f kg  stage prop %.6f kg  ' ...
         'dry %.6f kg  wet %.6f kg  thrust %.6f N  jet %.6f kW  dV %.6f km/s\n'], ...
        Ms.mfFrac, f, Ms.propKgPer150, Ms.stagePropKg, Ms.stageDryKg, Ms.wetKg, ...
        Ms.thrustN, Ms.jetKW, Ms.dvKms);
Mh  = leg1_gates('mass', 0.58, L, struct('dryFrac', 0.15));
hw  = [63, 63000/517, 9450/517, 150000/517, 5000/517, 33538743/413600, ...
       -1710*9.80665*log(0.58)/1000];
got = [Mh.propKgPer150, Mh.stagePropKg, Mh.stageDryKg, Mh.wetKg, Mh.thrustN, Mh.jetKW, Mh.dvKms];
fprintf('  mass (f = 0.42, dryFrac 0.15): %s\n', sprintf('%.12g ', got));
assert(max(abs(got - hw)./abs(hw)) < 1e-12, 'hand-worked mass block off by %.3g', ...
       max(abs(got - hw)./abs(hw)));
assert(Mh.feasible, 'hand-worked case infeasible');
% guard: mfFrac 0.1 (f 0.9: 1 - 0.9 - 0.135 < 0) and mfFrac 0 (no mass left)
for mfBad = [0.1, 0, NaN, 1.2]
    Mb = leg1_gates('mass', mfBad, L, struct('dryFrac', 0.15));
    nums = [Mb.propKgPer150, Mb.stagePropKg, Mb.stageDryKg, Mb.wetKg, Mb.thrustN, Mb.jetKW, Mb.dvKms];
    fprintf('  mass guard mfFrac %g: feasible %d (%s)\n', mfBad, Mb.feasible, Mb.reason);
    assert(~Mb.feasible && all(isnan(nums)) && ~isempty(Mb.reason), 'mass guard at %g', mfBad);
end

% ---- 3. deliberate failures, one per gate (second order skipped) ----
qs = struct('doSecond', false);
Pa = P;  Pa.z8(8) = Pa.z8(8) + 0.01;
Ga = leg1_gates(L, Pa, ctx, qs);
fprintf('  arrival, z8 tf + 0.01: pass %d  miss %.1f km %.1f m/s\n', Ga.arrival.pass, ...
        Ga.arrival.missKm, Ga.arrival.missMps);
assert(~Ga.arrival.pass && ~Ga.pass, 'perturbed z8 passed the arrival gate');

cb = ctx;  cb.orient0 = cb.orient0 + 1;
Gb = leg1_gates(L, P, cb, qs);
fprintf('  oneGto, orient0 + 1 deg: pass %d  dState %.3g ND\n', Gb.oneGto.pass, Gb.oneGto.dStateNd);
assert(~Gb.oneGto.pass && ~Gb.pass && Gb.arrival.pass && Gb.handoff.pass, 'oneGto failure');

cc = ctx;  cc.T0 = cc.T0 + 0.01;
Gc = leg1_gates(L, P, cc, qs);
fprintf('  handoff, T0 + 0.01: pass %d  errPhase %+.3g\n', Gc.handoff.pass, Gc.handoff.errPhase);
assert(~Gc.handoff.pass && ~Gc.pass && Gc.oneGto.pass, 'handoff T0 failure');

cr = ctx;  cr.kSat = mod(cr.kSat, 6) + 1;
Gr = leg1_gates(L, P, cr, qs);
fprintf('  Review Focus 5, kSat %d -> %d: pass %d  errPhase %+.6f\n', ctx.kSat, cr.kSat, ...
        Gr.handoff.pass, Gr.handoff.errPhase);
assert(~Gr.handoff.pass && abs(abs(Gr.handoff.errPhase) - 1/6) < 1e-5, 'kSat shift');

hM = G.moon.hMinKm;
Gm = leg1_gates(L, P, ctx, struct('doSecond', false, 'moonFloorKm', hM + 1));
Gp = leg1_gates(L, P, ctx, struct('doSecond', false, 'moonFloorKm', hM - 1));
fprintf('  Moon floor %.1f km: pass %d; floor %.1f km: pass %d\n', hM + 1, Gm.moon.pass, ...
        hM - 1, Gp.moon.pass);
assert(~Gm.moon.pass && ~Gm.pass && Gp.moon.pass && Gp.pass, 'Moon floor experiment');

% mass feasibility in G.pass: dryFrac 1.5 > (1 - f)/f = 1.41 makes the stage impossible
Gd = leg1_gates(L, P, ctx, struct('doSecond', false, 'dryFrac', 1.5));
fprintf('  dryFrac 1.5: mass feasible %d (%s)  pass %d\n', Gd.mass.feasible, Gd.mass.reason, Gd.pass);
assert(~Gd.mass.feasible && ~Gd.pass && Gd.arrival.pass && Gd.handoff.pass, 'infeasible stage passed');

% a truncated flight (first half, re-validated by validate_flight) is not admissible
Gt = leg1_gates(L, P, ctx, struct('doSecond', false, 'inject', struct('truncateFlight', 0.5)));
fprintf('  truncated flight: admissible %d (%s)  mass feasible %d  pass %d\n', Gt.flight.ok, ...
        Gt.flight.reason, Gt.mass.feasible, Gt.pass);
assert(~Gt.flight.ok && ~Gt.pass && ~Gt.mass.feasible && isnan(Gt.mass.wetKg), 'truncated flight');

% flown arrival epoch vs P.tArr
Pe = P;  Pe.tArr = Pe.tArr + 1e-8;  cE = ctx;  cE.T0 = cE.T0 + 1e-8;
Ge = leg1_gates(L, Pe, cE, qs);
fprintf('  P.tArr + 1e-8: handoff pass %d  dtArr %.2e  (%s)\n', Ge.handoff.pass, Ge.handoff.dtArr, ...
        Ge.handoff.reason);
assert(~Ge.handoff.pass && ~Ge.pass && abs(Ge.handoff.dtArr + 1e-8) < 1e-12, 'tArr disagreement');

% a fine re-integration that stops early fails the clearance gates
Gs = leg1_gates(L, P, ctx, struct('doSecond', false, 'inject', struct('propStopFrac', 0.5)));
fprintf('  re-integration stopped early: earth pass %d (%s)  moon pass %d  pass %d\n', ...
        Gs.earth.pass, Gs.earth.reason, Gs.moon.pass, Gs.pass);
assert(~Gs.earth.pass && ~Gs.moon.pass && ~Gs.pass && contains(Gs.earth.reason, 'stopped'), ...
       'early stop not failed');

% second order: an error inside the instruments is UNKNOWN, not refuted
Gx = leg1_gates(L, P, ctx, struct('inject', struct('secondError', true)));
fprintf('  second-order error: status %s  refuted %d  unknown %d  (%s)  pass %d\n', ...
        Gx.second.status, Gx.refuted, Gx.secondUnknown, Gx.second.reason, Gx.pass);
assert(strcmp(Gx.second.status, 'error') && ~Gx.refuted && Gx.secondUnknown && Gx.pass, ...
       'second error path');

% ---- 4. Review Focus 2: a perigee pass between two samples ----
[ts, Xs, prop, truthKm, rp] = perigee_arc(L);
pE = struct('center', [-L.muStar 0 0], 'radiusKm', 6378.137, 'lStar', L.lStar, ...
            'muStar', L.muStar, 'floorKm', 200);
hS = L.lStar*sqrt(sum((Xs(:, 1:3) - pE.center).^2, 2)) - pE.radiusKm;
[~, kS] = sort(abs(ts - rp.tPeri));
fprintf('  synthetic arc: perigee %.4f km at t %.3e ND; adjacent samples %.2f / %.2f km\n', ...
        truthKm, rp.tPeri, hS(kS(1)), hS(kS(2)));
assert(all(hS(kS(1:2)) > 200) && truthKm < 160, 'synthetic arc not as designed');
C = leg1_gates('located_min', ts, Xs, prop, pE);
fprintf('  located_min: sample %.3f km (pass %d)  located %.4f km at t %.3e (pass %d)  truth %.4f km\n', ...
        C.sampleMinKm, C.samplePass, C.hMinKm, C.tAtMin, C.pass, truthKm);
assert(C.samplePass && ~C.pass, 'Review Focus 2: located minimum did not fail the floor');
assert(abs(C.hMinKm - truthKm) < 0.01, 'located minimum %.4f vs truth %.4f', C.hMinKm, truthKm);
C0 = leg1_gates('located_min', ts, Xs, [], pE);
fprintf('  located_min without a propagator (Hermite only): %.4f km (pass %d)\n', C0.hMinKm, C0.pass);
assert(~C0.pass, 'Hermite-only located minimum passed');
CN = leg1_gates('located_min', ts, Xs, @(tA, xA, tB) NaN(6, 1), pE);
fprintf('  located_min with a NaN propagator: %.4g km  pass %d (%s)\n', CN.hMinKm, CN.pass, CN.reason);
assert(~CN.pass && ~isempty(CN.reason), 'NaN located minimum silently dropped');

% ---- 5. option names ----
threw = false;
try
    leg1_gates(L, P, ctx, struct('noSuchOpt', 1, 'doSecond', false));
catch err
    threw = strcmp(err.identifier, 'leg1_gates:opts');
end
assert(threw, 'unknown option name not refused');

fprintf('  PASS\n');
ok = true;
end

% ------------------------------------------------------------------------
function ctx = good_ctx(L, P)
% GOOD_CTX  A leg-1 context consistent with the polished transfer P: the
%   inertial GTO has P's orientation and departure phase at P.tDep; the
%   arrival slot is satellite kSat, the next 1/6 above P.arrPhase, at T0.
%   INPUTS: L leg1_inputs; P polish_gto_dro.  OUTPUTS: ctx struct.
ctx.orient0   = P.orientDeg - L.sgn*(180/pi)*P.tDep;
ctx.gtoPhase0 = mod(P.sD - P.tDep/L.tauG, 1);
ctx.kSat      = floor(6*P.arrPhase) + 2;
ctx.T0        = P.tArr + mod((ctx.kSat - 1)/6 - P.arrPhase, 1)*L.tauD;
end

% ------------------------------------------------------------------------
function [ts, Xs, prop, truthKm, rp] = perigee_arc(L)
% PERIGEE_ARC  A ballistic CR3BP arc through an Earth perigee at 150 km
%   altitude (inertial Kepler perigee, e = 0.7, converted with pumpkyn
%   fromPCI at t = 0), sampled every 8e-4 ND (306 s) with the perigee
%   midway between two samples.
%   INPUTS: L leg1_inputs.  OUTPUTS: ts [M x 1] ND; Xs [M x 6] ND; prop
%   handle (tA, xA, tB) -> [6 x 1]; truthKm fine-grid minimum altitude;
%   rp.tPeri its time.
muE = 6.67384e-20*(1 - L.muStar)*(5.9736E24 + 7.35E22);
rP  = 6378.137 + 150;   ecc = 0.7;
x0  = pumpkyn.cr3bp.fromPCI(0, [rP 0 0 0 sqrt(muE*(1 + ecc)/rP) 0], L.muStar, ...
                            L.tStar, L.lStar, 1).';
ode = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
prop = @(tA, xA, tB) fly6(tA, xA, tB, L.muStar, ode);
% fine truth around t = 0
tf_ = linspace(-2e-3, 2e-3, 40001).';
[~, Yf1] = ode113(@pumpkyn.cr3bp.eom, tf_(20001:-1:1), x0, ode, L.muStar, 1);
[~, Yf2] = ode113(@pumpkyn.cr3bp.eom, tf_(20001:end), x0, ode, L.muStar, 1);
Yf = [Yf1(end:-1:2, :); Yf2];
hF = L.lStar*sqrt(sum((Yf(:, 1:3) - [-L.muStar 0 0]).^2, 2)) - 6378.137;
[truthKm, kM] = min(hF);
rp.tPeri = tf_(kM);
dt = 8e-4;
ts = rp.tPeri + dt*((-6:5).' + 0.5);
Xs = zeros(numel(ts), 6);
for kk = 1:numel(ts)
    Xs(kk, :) = prop(0, x0, ts(kk)).';
end
end

% ------------------------------------------------------------------------
function x = fly6(tA, xA, tB, muStar, ode)
% FLY6  Ballistic CR3BP state at tB from (tA, xA).
%   INPUTS: tA, tB ND; xA [6 x 1]; muStar; ode odeset.  OUTPUTS: x [6 x 1].
if tB == tA
    x = xA(:);
    return
end
[~, Y] = ode113(@pumpkyn.cr3bp.eom, [tA, 0.5*(tA + tB), tB], xA(:), ode, muStar, 1);
x = Y(end, :).';
end
