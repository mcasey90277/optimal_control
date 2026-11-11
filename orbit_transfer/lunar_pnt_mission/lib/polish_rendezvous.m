function P = polish_rendezvous(M, k, jSlot, cand, opts)
%% Purpose:
%
%   EXACT-SLOT POLISH of one scheduled rendezvous. The scheduler picked a
%   catalog cell (iD, iA) whose arrival phase lies within a tolerance of
%   slot j's phase at the arrival epoch; this function turns it into a
%   minimum-time transfer that arrives ON slot j's exact state at the
%   PINNED arrival epoch tArr = cand.tArr.
%
%   For a trial duration tf the departure epoch is tDep = tArr - tf, so
%   the departure state slides along spacecraft k's DRO as tf changes:
%     rvf    = arrState(mod(slotPhase0(jSlot) + tArr/tauT, 1))      (fixed)
%     rv0(tf)= depState(mod(droPhase0(k) + tDep/tauD, 1))
%     T(tf)  = minimum time from rv0(tf) to rvf (eighth component of z8),
%              warm-started from the evaluated point nearest in tf (first
%              from the cell's z8); see INNER RULE for the solver
%     g(tf)  = tf - T(tf)
%   The root of g is found by secant from tf0 = M.tf(iD, iA) and T(tf0)
%   (the first step is the fixed-point step tf1 = T(tf0)); once a sign
%   change is bracketed a secant step that leaves the bracket is replaced
%   by bisection. Stop at |g| < tolG. An unbracketed step is capped at
%   familyFrac*tfCell (about 0.08 ND).
%
%   DISCONTINUITY: t_f jumps at transfer-family hand-overs, so g can be
%   discontinuous. A bracket narrower than min(|g_a|, |g_b|)/1e3 (a local
%   slope above 1e3 where the catalog's is about 0.4) is a jump, not a
%   root: ok = false, reason 'root not found'. maxEval exhausted is
%   'root not found' too.
%
%   INNER RULE: pumpkyn tfMin (fsolve table swallowed by evalc) is used
%   only when the warm z8's own departure state equals rv0, i.e. on the
%   first evaluation from the cell; every other evaluation (and a tfMin
%   miss) goes through ms_tfmin seeded by seed_from_z8 of the previous
%   solution flown from ITS OWN departure state, so only the first
%   junction is off. Measured 2026-11-10: single-shooting tfMin across a
%   0.0072 ND departure shift ran 241 s and missed by 9105 km, while
%   this seed converged in 7 iterations (1.5 s). A result counts only if
%   fly_transfer on it misses rvf by < 100 km and < 10 m/s. On a failed
%   trial the step is halved toward the last good point and retried up
%   to 3 times, then reason 'inner solve not converged'. This rule is
%   fixed and NOT subject to opts.gateOverride.
%
%   REFUSAL: an arrival earlier than the shortest catalog transfer
%   (tArr < min(M.tf(:))) is refused before any solve, and a root with
%   tDep < 0 is refused after the root-find; both give reason
%   'departure before epoch 0 ...' and run no gates.
%
%   ACCEPTANCE (opts.gates true), thresholds from DRO_tulip/indirect/
%   transfer_study.m held in ONE struct (gate_thresholds) that
%   opts.gateOverride may replace field by field:
%     resid  one ms_tfmin(tolR 3e-11, wallSec 600, conjTest) from
%            seed_from_z8(P.z8, rv0, 24) counted only if it REPRODUCES
%            the root (converged and max|z - P.z8| < 1e-6), normR < 1e-8
%     fly    flight admissible, flyKm < 100, flyVms < 10
%     pw     pmp_pointwise_checks: Hmax < 1e-6, lamMf < 1e-6, adjErr <
%            1e-7, max(fullGap, fieldGap) <= 1e-12, throttleErr < 1e-10,
%            minQmt >= 0
%     hyp    mintime_hypothesis_gates (relTol 1e-12): minLamV > 0,
%            minQmt > 0, dimS == 1
%     conj   reproduced and info.conj.pass == 1 and info.conj.covered
%     acc    pumpkyn tfMin from the final z8 with the final rv0/rvf
%            reproduces it: accDz = max|zAcc - z8| < 1e-6 (the catalog's
%            acceptance rule: the extremal is accepted by pumpkyn's own
%            solver)
%   With opts.gates false only the fly gate is evaluated (the final
%   flight is always flown); resid/pw/hyp/conj/acc are NaN / 'skipped'.
%   A throw in any certification call gives ok = false, reason
%   'gate error: <stage>: <message>', with the same struct shape. Test
%   seam: gateOverride.forceError = '<stage>' (fly|conj|pw|hyp|acc)
%   throws at that stage, checked in ONE place (force_check).
%
%   Self-demo (nargin == 0): spacecraft 3 of the real schedule (the
%   acceptance path), gates on, one table row printed.
%
%% Inputs:
%
%  M                        struct                  mission_inputs() output
%
%  k                        scalar                  spacecraft index (1..6)
%
%  jSlot                    scalar                  slot index (1..6)
%
%  cand                     struct                  .iD .iA .nRev .tDep .tArr
%                                                   (a schedule_rendezvous
%                                                   candidate; tDep, tArr ND
%                                                   clock times)
%
%  opts                     struct (optional)       .tolG [1e-9] ND
%                                                   .maxEval [25]
%                                                   .familyFrac [0.02]
%                                                   .gates [true]
%                                                   .gateOverride [struct()]
%                                                   fields of gate_thresholds
%                                                   to replace
%                                                   .verbose [false]
%
%% Outputs:
%
%  P                        struct                  .ok (logical) .reason
%                                                   (char) .z8 [8 x 1] .tf
%                                                   .tDep .tArr (ND) .rv0
%                                                   .rvf [6 x 1] .nEval
%                                                   .gResid (tf - z8(8), ND)
%                                                   .tfCell .familyChanged
%                                                   .flight (fly_transfer)
%                                                   .gates (.resid .flyKm
%                                                   .flyVms .passFly .pw
%                                                   .passPW .hyp .passHyp
%                                                   .conjPass .conjCovered
%                                                   .conjReason .passConj
%                                                   .passResid .reproduced
%                                                   .dz .accDz .passAcc
%                                                   .failed {names})
%                                                   .propKg .dvKms .wallS
%                                                   .wallRootS .wallCertS
%                                                   .history [nEval x 4]
%                                                   [tf T g innerRoute]
%                                                   (innerRoute 1 tfMin,
%                                                   2 ms_tfmin)
%                                                   .thresholds
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  M. Casey  slot index renamed jSlot (house rule)            (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    run(fullfile(here, 'setup_lunar_pnt_paths.m'));
    M = mission_inputs();
    S = schedule_rendezvous(M.tf, M.sD, M.sA, M.tauD, M.tauT, M.droPhase0, ...
                            M.slotPhase0, struct('propKgPerNd', M.propKgPerNd));
    k = 3;
    jSlot = S.assign(k);
    cand = struct('iD', S.iD(k), 'iA', S.iA(k), 'nRev', S.nRev(k), ...
                  'tDep', S.tDep(k), 'tArr', S.tArr(k));
    opts = struct('gates', true, 'verbose', true);
end
if nargin > 0 && nargin < 5, opts = struct(); end
P = polish_core(M, k, jSlot, cand, opts);
if nargin == 0 || fdef(opts, 'verbose', false)
    print_row(M, k, jSlot, cand, P);
end
end

% ------------------------------------------------------------------------
function P = polish_core(M, k, jSlot, cand, opts)
% POLISH_CORE  The root-find and certification (see the main header).
%   INPUTS: as polish_rendezvous.  OUTPUTS: P as polish_rendezvous.
tolG       = fdef(opts, 'tolG',       1e-9);
maxEval    = fdef(opts, 'maxEval',    25);
familyFrac = fdef(opts, 'familyFrac', 0.02);
doGates    = fdef(opts, 'gates',      true);
verbose    = fdef(opts, 'verbose',    false);
thr        = gate_thresholds(fdef(opts, 'gateOverride', struct()));

innerKm    = 100;            % inner-convergence rule (fixed, not overridable)
innerVms   = 10;
jumpRatio  = 1e3;            % bracket width vs min |g| at its ends
maxRetry   = 3;             % step halvings after a failed inner solve
Kseg       = 24;
tStart     = tic;

tArr   = cand.tArr;
tfCell = M.tf(cand.iD, cand.iA);
maxStep = familyFrac*tfCell;   % ND cap on an unbracketed secant step
z8Cell = M.sh.z8(:, M.sh.entry_index(cand.iD, cand.iA));
rvf    = M.arrState(mod(M.slotPhase0(jSlot) + tArr/M.tauT, 1));
depAt  = @(tf) M.depState(mod(M.droPhase0(k) + (tArr - tf)/M.tauD, 1));

P = struct('ok', false, 'reason', '', 'z8', nan(8, 1), 'tf', NaN, ...
           'tDep', NaN, 'tArr', tArr, 'rv0', nan(6, 1), 'rvf', rvf, ...
           'nEval', 0, 'gResid', NaN, 'tfCell', tfCell, ...
           'familyChanged', false, 'flight', struct(), ...
           'gates', skipped_gates(), 'propKg', NaN, 'dvKms', NaN, ...
           'wallS', NaN, 'wallRootS', NaN, 'wallCertS', NaN, ...
           'history', zeros(0, 4), 'thresholds', thr);

% ---- preflight refusal ----
if tArr < min(M.tf(:))
    P.reason = sprintf(['departure before epoch 0 (arrival %.4f ND is earlier ' ...
                        'than the shortest catalog transfer %.4f ND)'], ...
                       tArr, min(M.tf(:)));
    P.wallS = toc(tStart);  P.wallRootS = P.wallS;  P.wallCertS = 0;
    return
end

% ---- root-find on tf ----
Xa = nan(1, maxEval);  Ga = nan(1, maxEval);  Za = nan(8, maxEval);
H  = nan(maxEval, 4);
nE = 0;  nH = 0;  nRetry = 0;
x = tfCell;
found = false;  why = '';
while nH < maxEval
    X = Xa(1:nE);  Z = Za(:, 1:nE);
    % warm start: the evaluated point nearest in tf (first: the cell)
    if isempty(X)
        zw  = z8Cell;
        rvw = depAt(tfCell);
    else
        [~, in] = min(abs(X - x));
        zw  = Z(:, in);
        rvw = depAt(X(in));
    end
    rv0 = depAt(x);
    [z, route, inOk] = inner_solve(rv0, rvf, zw, rvw, M, Kseg, innerKm, innerVms);
    if ~inOk
        nH = nH + 1;  H(nH, :) = [x NaN NaN route];
        if verbose
            fprintf('  polish eval %2d: tf %.12f  inner solve FAILED (route %d)\n', ...
                    nH, x, route);
        end
        if nE == 0 || nRetry >= maxRetry
            why = sprintf('inner solve not converged (tf %.10f ND, %d halvings)', ...
                          x, nRetry);
            break
        end
        nRetry = nRetry + 1;
        x = Xa(nE) + 0.5*(x - Xa(nE));     % halve toward the last good point
        continue
    end
    nRetry = 0;
    gx = x - z(8);
    nE = nE + 1;  Xa(nE) = x;  Ga(nE) = gx;  Za(:, nE) = z;
    nH = nH + 1;  H(nH, :) = [x z(8) gx route];
    X = Xa(1:nE);  Gv = Ga(1:nE);
    if verbose
        fprintf('  polish eval %2d: tf %.12f  T %.12f  g %+.3e  route %d\n', ...
                nE, x, z(8), gx, route);
    end
    if abs(gx) < tolG
        found = true;
        break
    end
    % bracket from every evaluated point: tightest sign change
    [a, b, ga, gb] = tightest_bracket(X, Gv);
    if ~isnan(a) && (b - a) < min(abs(ga), abs(gb))/jumpRatio
        why = sprintf(['root not found (g discontinuous: bracket [%.10f, %.10f] ' ...
                       'width %.2e, g %+.2e / %+.2e)'], a, b, b - a, ga, gb);
        break
    end
    % next trial
    if nE == 1
        xn = z(8);                                   % fixed-point step
    else
        x1 = X(end-1);  g1 = Gv(end-1);  x2 = X(end);  g2 = Gv(end);
        if g2 ~= g1
            xn = x2 - g2*(x2 - x1)/(g2 - g1);
        else
            xn = NaN;
        end
    end
    if ~isnan(a)
        if ~(isfinite(xn) && xn > a && xn < b)
            xn = 0.5*(a + b);
        end
    else
        if ~isfinite(xn), xn = z(8); end
        xn = X(end) + max(-maxStep, min(maxStep, xn - X(end)));
    end
    x = xn;
end
X = Xa(1:nE);  Gv = Ga(1:nE);  Z = Za(:, 1:nE);
P.nEval   = nH;
P.history = H(1:nH, :);
if ~found && isempty(why)
    why = sprintf('root not found (maxEval %d reached, min |g| %.2e)', ...
                  maxEval, min(abs(Gv)));
end
P.wallRootS = toc(tStart);

if ~found
    P.reason = why;
    if ~isempty(X)
        [~, ib] = min(abs(Gv));
        P.tf = X(ib);  P.z8 = Z(:, ib);  P.gResid = Gv(ib);
        P.tDep = tArr - P.tf;  P.rv0 = depAt(P.tf);
    end
    P.wallS = toc(tStart);  P.wallCertS = 0;
    return
end

P.tf     = x;
P.z8     = z;
P.gResid = gx;
P.tDep   = tArr - x;
P.rv0    = depAt(x);
P.familyChanged = abs(P.tf - tfCell)/tfCell > familyFrac;

if P.tDep < 0
    P.reason = sprintf('departure before epoch 0 (root tDep = %.6f ND)', P.tDep);
    P.wallS = toc(tStart);  P.wallCertS = 0;
    return
end

% ---- certification ----
tCert = tic;
G = skipped_gates();
stage = 'fly';
try
    force_check(thr, stage);
    P.flight = fly_transfer(P.z8, P.rv0, P.rvf, M.phys);
    P.propKg = P.flight.propellantKg;
    P.dvKms  = P.flight.dvKms;
    G.flyKm  = P.flight.flyKm;
    G.flyVms = P.flight.flyVms;
    G.flyAdmissible = P.flight.admissibility.ok;
    if doGates
        stage = 'conj';
        force_check(thr, stage);
        seed = seed_from_z8(P.z8, P.rv0, Kseg, M.Tnd, M.cnd, M.muStar);
        [zc, info] = ms_tfmin(P.rv0, P.rvf, seed, M.Tnd, M.cnd, M.muStar, ...
                              struct('tolR', 3e-11, 'wallSec', 600, 'conjTest', true));
        G.dz = max(abs(zc - P.z8));
        G.reproduced = info.converged && G.dz < thr.reproTol;
        G.resid = info.normR;
        if G.reproduced
            G.conjPass    = double(info.conj.pass);
            G.conjCovered = double(info.conj.covered);
            G.conjReason  = sprintf('reproduced; %s: %s', info.conj.verdict, info.conj.reason);
        else
            G.conjReason  = sprintf('not reproduced (converged %d, |dz| = %.1e)', ...
                                    info.converged, G.dz);
        end
        stage = 'pw';
        force_check(thr, stage);
        G.pw  = pmp_pointwise_checks(P.flight.t, P.flight.Y, M.Tnd, M.cnd, M.muStar);
        stage = 'hyp';
        force_check(thr, stage);
        G.hyp = mintime_hypothesis_gates(P.z8, P.rv0, M.Tnd, M.cnd, M.muStar, ...
                                         struct('relTol', 1e-12));
        stage = 'acc';
        force_check(thr, stage);
        zAcc = quiet_tfmin(P.rv0, P.rvf, P.z8, M);
        G.accDz = max(abs(zAcc - P.z8));
    end
catch gerr
    G.failed = {stage};
    P.gates = G;
    P.wallCertS = toc(tCert);
    P.reason = sprintf('gate error: %s: %s', stage, gerr.message);
    P.wallS = toc(tStart);
    return
end
P.gates = evaluate_gates(G, thr, doGates);
P.wallCertS = toc(tCert);
if isempty(P.gates.failed)
    P.ok = true;
    P.reason = 'accepted';
else
    P.reason = ['gate failed: ' strjoin(P.gates.failed, ', ')];
end
P.wallS = toc(tStart);
end

% ------------------------------------------------------------------------
function print_row(M, k, jSlot, cand, P)
% PRINT_ROW  One summary block of a polish (any outcome).
%   INPUTS: M; k; jSlot; cand; P (polish_rendezvous output).  OUTPUTS: none.
dayNd = M.tStar/86400;
gapH  = mod(M.slotPhase0(jSlot) + cand.tArr/M.tauT, 1) - M.sA(cand.iA);
gapH  = (gapH - round(gapH))*M.tauT*M.tStar/3600;
fprintf(['sc %d slot %d cell (%d,%d): ok %d, %s\n  tf cell %.6f d -> polished %.6f d ' ...
         '(%+.3f h), arrival phase gap closed %.3f h, tDep %.6f d, nEval %d, ' ...
         'g %.2e\n  familyChanged %d, wall %.1f s (root %.1f, cert %.1f)\n'], ...
        k, jSlot, cand.iD, cand.iA, P.ok, P.reason, P.tfCell*dayNd, P.tf*dayNd, ...
        (P.tf - P.tfCell)*dayNd*24, gapH, P.tDep*dayNd, P.nEval, P.gResid, ...
        P.familyChanged, P.wallS, P.wallRootS, P.wallCertS);
if isfield(P.flight, 'flyKm')
    fprintf('  fly %.3g km %.3g m/s, prop %.3f kg, dV %.4f km/s\n', ...
            P.flight.flyKm, P.flight.flyVms, P.propKg, P.dvKms);
end
g = P.gates;
if isfield(g.pw, 'Hmax') && isfield(g.hyp, 'minLamV')
    pw = g.pw;  hy = g.hyp;
    fprintf(['  resid %.2e (reproduced %d, |dz| %.1e)  conj %g covered %g: %s\n' ...
             '  tfMin acceptance |dz| %.2e\n'], g.resid, g.reproduced, g.dz, ...
            g.conjPass, g.conjCovered, g.conjReason, g.accDz);
    fprintf(['  PW Hmax %.2e lamMf %.2e adjErr %.2e fullGap %.2e fieldGap %.2e ' ...
             'throttleErr %.2e minQmt %.3g\n  hyp minLamV %.3g minQmt %.3g dimS %d\n'], ...
            pw.Hmax, pw.lamMf, pw.adjErr, pw.fullGap, pw.fieldGap, ...
            pw.throttleErr, pw.minQmt, hy.minLamV, hy.minQmt, hy.dimS);
end
end

% ------------------------------------------------------------------------
function [z, route, ok] = inner_solve(rv0, rvf, zw, rvw, M, Kseg, innerKm, innerVms)
% INNER_SOLVE  T(tf) judged by the flown miss. tfMin (silenced) is tried
%   only when the warm z8 belongs to the SAME departure state (the cell
%   evaluation); otherwise, or when tfMin misses, ms_tfmin seeded by
%   seed_from_z8 of the warm z8 flown from ITS OWN departure state rvw
%   (only the first junction is then off; measured 2026-11-10: tfMin from
%   a 0.0072 ND departure shift exhausted fsolve's 1000 evaluations in
%   241 s and missed by 9105 km, and a seed flown from the NEW departure
%   missed by 2e4 km and did not converge, while this seed converged in
%   7 iterations, 1.5 s).
%   INPUTS: rv0, rvf [6 x 1]; zw [8 x 1] warm z8; rvw [6 x 1] the warm
%   z8's departure state; M mission inputs; Kseg segments; innerKm,
%   innerVms miss limits.
%   OUTPUTS: z [8 x 1]; route (1 tfMin, 2 ms_tfmin); ok logical.
z  = nan(8, 1);
if max(abs(rv0 - rvw)) < 1e-12
    route = 1;
    try
        z = quiet_tfmin(rv0, rvf, zw, M);
        ok = flies_onto(z, rv0, rvf, M, innerKm, innerVms);
    catch
        ok = false;
    end
    if ok, return; end
end
route = 2;
try
    seed = seed_from_z8(zw, rvw, Kseg, M.Tnd, M.cnd, M.muStar);
    [z, info] = ms_tfmin(rv0, rvf, seed, M.Tnd, M.cnd, M.muStar, ...
                         struct('tolR', 3e-11, 'wallSec', 300));
    ok = info.converged && flies_onto(z, rv0, rvf, M, innerKm, innerVms);
catch
    ok = false;
end
end

% ------------------------------------------------------------------------
function ok = flies_onto(z, rv0, rvf, M, innerKm, innerVms)
% FLIES_ONTO  Flown miss of z8 under the inner limits.
%   INPUTS: z [8 x 1]; rv0, rvf [6 x 1]; M; innerKm, innerVms.
%   OUTPUTS: ok logical.
ok = false;
if ~all(isfinite(z)) || z(8) <= 0, return; end
F = fly_transfer(z, rv0, rvf, M.phys);
ok = F.flyKm < innerKm && F.flyVms < innerVms;
end

% ------------------------------------------------------------------------
function [a, b, ga, gb] = tightest_bracket(X, Gv)
% TIGHTEST_BRACKET  Narrowest interval between adjacent (in tf) evaluated
%   points whose g values change sign.
%   INPUTS: X [1 x n] tf values; Gv [1 x n] g values.
%   OUTPUTS: a < b bracket ends and their g (NaN when none).
a = NaN;  b = NaN;  ga = NaN;  gb = NaN;
[xs, is] = sort(X);
gs = Gv(is);
w = inf;
for kk = 1:numel(xs) - 1
    if sign(gs(kk)) ~= sign(gs(kk+1)) && (xs(kk+1) - xs(kk)) < w
        w = xs(kk+1) - xs(kk);
        a = xs(kk);  b = xs(kk+1);  ga = gs(kk);  gb = gs(kk+1);
    end
end
end

% ------------------------------------------------------------------------
function thr = gate_thresholds(override)
% GATE_THRESHOLDS  Acceptance thresholds (transfer_study.m values); the
%   override replaces named fields and refuses unknown names.
%   INPUTS: override struct.  OUTPUTS: thr struct.
thr = struct( ...
    'residMax',       1e-8, ...
    'reproTol',       1e-6, ...
    'flyKmMax',       100, ...
    'flyVmsMax',      10, ...
    'HmaxMax',        1e-6, ...
    'lamMfMax',       1e-6, ...
    'adjErrMax',      1e-7, ...
    'gapMax',         1e-12, ...
    'throttleErrMax', 1e-10, ...
    'pwMinQmtMin',    0, ...
    'hypMinLamVMin',  0, ...
    'hypMinQmtMin',   0, ...
    'dimS',           1, ...
    'accDzMax',       1e-6, ...
    'forceError',     '');
fn = fieldnames(override);
for kk = 1:numel(fn)
    if ~isfield(thr, fn{kk})
        error('polish_rendezvous:gateOverride', 'unknown threshold "%s"', fn{kk});
    end
    thr.(fn{kk}) = override.(fn{kk});
end
end

% ------------------------------------------------------------------------
function G = skipped_gates()
% SKIPPED_GATES  Gate struct with every value NaN / 'skipped'.
%   INPUTS: none.  OUTPUTS: G struct.
G = struct('resid', NaN, 'flyKm', NaN, 'flyVms', NaN, 'flyAdmissible', NaN, ...
           'passFly', NaN, 'pw', struct(), 'passPW', NaN, 'hyp', struct(), ...
           'passHyp', NaN, 'conjPass', NaN, 'conjCovered', NaN, ...
           'conjReason', 'skipped', 'passConj', NaN, 'passResid', NaN, ...
           'reproduced', NaN, 'dz', NaN, 'accDz', NaN, 'passAcc', NaN, ...
           'failed', {{}});
end

% ------------------------------------------------------------------------
function G = evaluate_gates(G, thr, doGates)
% EVALUATE_GATES  Pass flags against thresholds and the list of failed
%   gate names (fly always; resid, pw, hyp, conj only when measured).
%   INPUTS: G gate values; thr thresholds; doGates logical.
%   OUTPUTS: G with .pass* and .failed.
failed = {};
G.passFly = logical(G.flyAdmissible) && G.flyKm < thr.flyKmMax && ...
            G.flyVms < thr.flyVmsMax;
if ~G.passFly, failed{end+1} = 'fly'; end
if doGates
    G.passResid = logical(G.reproduced) && G.resid < thr.residMax;
    pw = G.pw;
    G.passPW = pw.Hmax < thr.HmaxMax && pw.lamMf < thr.lamMfMax && ...
               pw.adjErr < thr.adjErrMax && ...
               max(pw.fullGap, pw.fieldGap) <= thr.gapMax && ...
               pw.throttleErr < thr.throttleErrMax && pw.minQmt >= thr.pwMinQmtMin;
    hy = G.hyp;
    G.passHyp = hy.minLamV > thr.hypMinLamVMin && hy.minQmt > thr.hypMinQmtMin && ...
                hy.dimS == thr.dimS;
    G.passConj = logical(G.reproduced) && G.conjPass == 1 && G.conjCovered == 1;
    G.passAcc  = G.accDz < thr.accDzMax;
    if ~G.passResid, failed{end+1} = 'resid'; end
    if ~G.passPW,    failed{end+1} = 'pw';    end
    if ~G.passHyp,   failed{end+1} = 'hyp';   end
    if ~G.passConj,  failed{end+1} = 'conj';  end
    if ~G.passAcc,   failed{end+1} = 'acc';   end
end
G.failed = failed;
end

% ------------------------------------------------------------------------
function z = quiet_tfmin(rv0, rvf, zw, M)
% QUIET_TFMIN  pumpkyn tfMin with its fsolve iteration table swallowed.
%   INPUTS: rv0, rvf [6 x 1]; zw [8 x 1] guess; M (.Tnd .cnd .muStar).
%   OUTPUTS: z [8 x 1].
z = quiet_call(@pumpkyn.cr3bp.tfMin, rv0, rvf, zw, M.Tnd, M.cnd, M.muStar);
end

% ------------------------------------------------------------------------
function out = quiet_call(varargin)
% QUIET_CALL  feval(varargin{:}) with its printed output captured by evalc.
%   INPUTS: varargin {fhandle, args...}.  OUTPUTS: out (first output).
[~, out] = evalc('feval(varargin{:})');
end

% ------------------------------------------------------------------------
function force_check(thr, stage)
% FORCE_CHECK  Test seam: throw at the named certification stage.
%   INPUTS: thr thresholds (.forceError char); stage char.  OUTPUTS: none.
if strcmp(thr.forceError, stage)
    error('polish_rendezvous:forcedError', 'forced error at stage %s (test seam)', stage);
end
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
