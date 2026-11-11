function P = polish_gto_dro(L, cand, opts)
%% Purpose:
%
%   POLISH OF ONE LEG-1 TRANSFER (GTO -> DRO, 5 N, minimum time). A
%   catalog cell (iD, iA) of sheet iS was solved at the sheet's GTO
%   orientation and onto the DRO grid phase sA(iA); this function turns
%   it into a minimum-time transfer that departs from the GTO at its TRUE
%   orientation cand.orientDeg (departure phase sD(iD), departure state
%   PINNED) and arrives ON a slot whose DRO phase at the departure epoch
%   is cand.slotPhaseAtDep, so the arrival target slides along the DRO
%   with t_f (the reverse of leg 2, where the arrival was pinned).
%
%   (a) WALK. From the sheet orientation to cand.orientDeg, arrival held
%   at the cell's own phase sA(iA), in steps of opts.stepDeg (the last
%   step lands exactly on the target). Each step is one INNER SOLVE. A
%   failed step is halved up to opts.maxHalvings times, then ok = false,
%   reason 'orientation walk stalled'. A step whose t_f differs from the
%   previous one by more than opts.jumpFrac of it is a BRANCH JUMP: same
%   halving, then reason 'branch jump'. After an accepted step the step
%   returns to opts.stepDeg. orientReached is the last accepted
%   orientation.
%
%   (b) ROOT. Departure fixed at rv0 = L.gtoState(cand.orientDeg, sD);
%     rvf(tf) = L.slotState(mod(cand.slotPhaseAtDep + tf/L.tauD, 1))
%     T(tf)   = minimum time from rv0 to rvf(tf) (inner solve, warm from
%               the evaluated point nearest in tf; first from the walk)
%     g(tf)   = tf - T(tf)
%   Secant from the walk's t_f (first step the fixed-point step T), with
%   bracketing: once a sign change is bracketed a secant step leaving the
%   bracket is replaced by bisection; an unbracketed step is capped at
%   opts.familyFrac of the cell t_f. DISCONTINUITY (polish_rendezvous
%   rule): a bracket narrower than min(|g_a|, |g_b|)/1e3 is a jump in
%   T, not a root: reason 'root not found'. Stop at |g| < opts.tolG. A
%   failed inner solve is halved toward the last good point up to 3
%   times, then reason 'inner solve not converged'.
%
%   INNER SOLVE (fixed, not subject to gateOverride; the leg-2 rule):
%   ms_tfmin(rv0, rvf, seed_from_z8(zPrev, rvPrev, K), tolR 3e-11,
%   wallSec 300) where zPrev is the previous solution flown from ITS OWN
%   departure state rvPrev; counted only if converged and fly_transfer
%   misses rvf by < 100 km and < 10 m/s. Measured 2026-11-11 (cell (1,1)
%   and (7,1) of the 0 deg sheet, K = 24): 2 deg steps converge in 4-10
%   iterations, about 0.6 s each; one 60 deg step from perigee does not
%   converge (100 iterations, residual 1.5).
%
%   ARRIVAL GATE (on the root, thresholds in ONE struct that
%   opts.gateOverride may replace field by field; unknown names throw
%   polish_gto_dro:gateOverride): fly_transfer flight admissible is NOT
%   required here (Earth clearance is gated by the caller), flyKm <
%   flyKmMax [100], flyVms < flyVmsMax [10], final residual res <
%   residMax [1e-8]. missKm / missMps report the achieved values.
%
%   CONTRACT: P.z8 is non-empty IFF P.ok. On any failure z8 = [] and
%   tf/tArr/arrPhase are NaN, so no transfer from a different departure
%   state (or from an unpolished cell) can be used by mistake; rv0 is
%   always the REQUESTED departure state.
%
%   Self-demo (nargin == 0): cell (1,1) of the first sheet, +10 deg, a
%   slot 0.03 past the arrival phase, verbose.
%
%% Inputs:
%
%  L                        struct                  leg1_inputs(struct())
%
%  cand                     struct                  .iD .iA .iS cell and
%                                                   sheet indices
%                                                   .orientDeg true GTO
%                                                   orientation at departure
%                                                   [deg]
%                                                   .tDep departure epoch,
%                                                   leg-1 clock [ND]
%                                                   .slotPhaseAtDep slot DRO
%                                                   phase at tDep [fraction]
%
%  opts                     struct (optional)       .stepDeg [2]
%                                                   .maxHalvings [3]
%                                                   .jumpFrac [0.05]
%                                                   .tolG [1e-10] ND
%                                                   .maxEval [25]
%                                                   .familyFrac [0.02]
%                                                   .K [24] ms segments
%                                                   .gateOverride [struct()]
%                                                   .verbose [false]
%
%% Outputs:
%
%  P                        struct                  .ok (logical) .reason
%                                                   (char) .z8 [8 x 1] or []
%                                                   .tf .tDep .tArr [ND]
%                                                   .rv0 .rvf [6 x 1]
%                                                   .orientDeg [deg] .sD
%                                                   .arrPhase [fraction]
%                                                   .missKm .missMps .res
%                                                   .nEval (root evaluations
%                                                   incl. failed) .walk
%                                                   (.orientDeg [1 x n]
%                                                   accepted, .tf [1 x n],
%                                                   .attempts [m x 6]
%                                                   [orient step ok tf iters
%                                                   seconds]) .orientReached
%                                                   .flight (fly_transfer
%                                                   of the root, or struct())
%                                                   .history [nEval x 4]
%                                                   [tf T g iters] .tfCell
%                                                   .familyChanged (logical:
%                                                   |tf - tfCell|/tfCell >
%                                                   familyFrac; false unless
%                                                   ok) .thresholds .wallS
%                                                   .wallWalkS .wallRootS
%
%% Revision History:
%  M. Casey                                                   (c) 11/11/2026
%  M. Casey  final fix wave: P.familyChanged                  (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    run(fullfile(here, 'setup_lunar_pnt_paths.m'));
    L = leg1_inputs(struct());
    cand = struct('iD', 1, 'iA', 1, 'iS', 1, 'orientDeg', L.sheetOrient(1) + 10, ...
                  'tDep', 0, 'slotPhaseAtDep', ...
                  mod(L.sA(1) + 0.03 - L.tf(1, 1, 1)/L.tauD, 1));
    opts = struct('verbose', true);
end
if nargin > 0 && nargin < 3, opts = struct(); end

stepDeg     = fdef(opts, 'stepDeg',     2);
maxHalvings = fdef(opts, 'maxHalvings', 3);
jumpFrac    = fdef(opts, 'jumpFrac',    0.05);
tolG        = fdef(opts, 'tolG',        1e-10);
maxEval     = fdef(opts, 'maxEval',     25);
familyFrac  = fdef(opts, 'familyFrac',  0.02);
Kseg        = fdef(opts, 'K',           24);
verbose     = fdef(opts, 'verbose',     false);
thr         = gate_thresholds(fdef(opts, 'gateOverride', struct()));

jumpRatio   = 1e3;          % root bracket width vs min |g| at its ends
maxRetry    = 3;            % root: halvings after a failed inner solve
tStart      = tic;

sD     = L.sD(cand.iD);
o0     = L.sheetOrient(cand.iS);
tfCell = L.tf(cand.iD, cand.iA, cand.iS);
rv0    = L.gtoState(cand.orientDeg, sD);
P = struct('ok', false, 'reason', '', 'z8', zeros(0, 1), 'tf', NaN, ...
           'tDep', cand.tDep, 'tArr', NaN, 'rv0', rv0, 'rvf', nan(6, 1), ...
           'orientDeg', cand.orientDeg, 'sD', sD, 'arrPhase', NaN, ...
           'missKm', NaN, 'missMps', NaN, 'res', NaN, 'nEval', 0, ...
           'walk', struct('orientDeg', o0, 'tf', tfCell, 'attempts', zeros(0, 6)), ...
           'orientReached', o0, 'flight', struct(), 'history', zeros(0, 4), ...
           'tfCell', tfCell, 'familyChanged', false, 'thresholds', thr, 'wallS', NaN, ...
           'wallWalkS', NaN, 'wallRootS', NaN);
if ~L.clean(cand.iD, cand.iA, cand.iS)
    P.reason = sprintf('cell (%d,%d) of sheet %d is not clean', cand.iD, cand.iA, cand.iS);
    P.wallS = toc(tStart);
    return
end

% ---- (a) orientation walk, arrival at the cell's own phase ----
rvfCell = L.slotState(L.sA(cand.iA));
z       = L.z8(:, cand.iD, cand.iA, cand.iS);
o       = o0;
rvPrev  = L.gtoState(o, sD);
% bound: every accepted step advances >= stepDeg/2^maxHalvings and costs
% at most maxHalvings + 1 attempts
nAccMax = ceil(abs(cand.orientDeg - o0)*2^maxHalvings/stepDeg) + 1;
wO  = [o0, nan(1, nAccMax)];  wT = [z(8), nan(1, nAccMax)];  nW = 1;
att = nan((nAccMax + 1)*(maxHalvings + 1), 6);                 nAt = 0;
nHalf = 0;  dO = stepDeg;  walkWhy = '';
while o ~= cand.orientDeg
    rem_ = cand.orientDeg - o;
    if dO >= abs(rem_)
        on = cand.orientDeg;                     % land exactly on the target
    else
        on = o + sign(rem_)*dO;
    end
    rvOn = L.gtoState(on, sD);
    t0 = tic;
    [zn, okIn, its] = inner_solve(rvOn, rvfCell, z, rvPrev, L, Kseg);
    jump = okIn && abs(zn(8) - z(8)) > jumpFrac*z(8);
    okStep = okIn && ~jump;
    nAt = nAt + 1;
    att(nAt, :) = [on, on - o, okStep, zn(8), its, toc(t0)];
    if verbose
        fprintf('  walk %+9.4f deg (step %+7.4f): ok %d  tf %.9f  iters %d  %.1f s%s\n', ...
                on, on - o, okStep, zn(8), its, att(nAt, 6), ...
                ternary(jump, '  BRANCH JUMP', ''));
    end
    if okStep
        o = on;  z = zn;  rvPrev = rvOn;
        nW = nW + 1;  wO(nW) = o;  wT(nW) = z(8);
        nHalf = 0;  dO = stepDeg;
        continue
    end
    if jump
        walkWhy = sprintf(['branch jump (step %+.4g deg to %.4g deg moved tf %.6f -> ' ...
                           '%.6f ND, more than %.3g of it)'], on - o, on, z(8), zn(8), jumpFrac);
    else
        walkWhy = sprintf('orientation walk stalled (step %+.4g deg to %.4g deg not converged)', ...
                          on - o, on);
    end
    if nHalf >= maxHalvings
        break
    end
    nHalf = nHalf + 1;
    dO = abs(on - o)/2;
end
P.walk = struct('orientDeg', wO(1:nW), 'tf', wT(1:nW), 'attempts', att(1:nAt, :));
P.orientReached = o;
P.wallWalkS = toc(tStart);
if o ~= cand.orientDeg
    P.reason = sprintf('%s after %d halvings; reached %.4g deg of %.4g', walkWhy, ...
                       nHalf, o, cand.orientDeg);
    P.wallS = toc(tStart);  P.wallRootS = 0;
    return
end

% ---- (b) root on tf, departure pinned, arrival slides ----
rvfAt = @(tf) L.slotState(mod(cand.slotPhaseAtDep + tf/L.tauD, 1));
maxStep = familyFrac*tfCell;
Xa = nan(1, maxEval);  Ga = nan(1, maxEval);  Za = nan(8, maxEval);  Ra = nan(1, maxEval);
H  = nan(maxEval, 4);
nE = 0;  nH = 0;  nRetry = 0;
x = z(8);  zWalk = z;
found = false;  why = '';
while nH < maxEval
    if nE == 0
        zw = zWalk;
    else
        [~, in] = min(abs(Xa(1:nE) - x));
        zw = Za(:, in);
    end
    [zr, okIn, its, nr] = inner_solve(rv0, rvfAt(x), zw, rv0, L, Kseg);
    nH = nH + 1;
    if ~okIn
        H(nH, :) = [x NaN NaN its];
        if verbose
            fprintf('  root eval %2d: tf %.12f  inner solve FAILED\n', nH, x);
        end
        if nE == 0 || nRetry >= maxRetry
            why = sprintf('inner solve not converged (tf %.10f ND, %d halvings)', x, nRetry);
            break
        end
        nRetry = nRetry + 1;
        x = Xa(nE) + 0.5*(x - Xa(nE));
        continue
    end
    nRetry = 0;
    gx = x - zr(8);
    nE = nE + 1;  Xa(nE) = x;  Ga(nE) = gx;  Za(:, nE) = zr;  Ra(nE) = nr;
    H(nH, :) = [x zr(8) gx its];
    if verbose
        fprintf('  root eval %2d: tf %.12f  T %.12f  g %+.3e  iters %d\n', nH, x, zr(8), gx, its);
    end
    if abs(gx) < tolG
        found = true;
        break
    end
    X = Xa(1:nE);  Gv = Ga(1:nE);
    [a, b, ga, gb] = tightest_bracket(X, Gv);
    if ~isnan(a) && (b - a) < min(abs(ga), abs(gb))/jumpRatio
        why = sprintf(['root not found (g discontinuous: bracket [%.10f, %.10f] ' ...
                       'width %.2e, g %+.2e / %+.2e)'], a, b, b - a, ga, gb);
        break
    end
    if nE == 1
        xn = zr(8);                                   % fixed-point step
    else
        x1 = X(end-1);  g1 = Gv(end-1);  x2 = X(end);  g2 = Gv(end);
        if g2 ~= g1, xn = x2 - g2*(x2 - x1)/(g2 - g1); else, xn = NaN; end
    end
    if ~isnan(a)
        if ~(isfinite(xn) && xn > a && xn < b), xn = 0.5*(a + b); end
    else
        if ~isfinite(xn), xn = zr(8); end
        xn = X(end) + max(-maxStep, min(maxStep, xn - X(end)));
    end
    x = xn;
end
P.nEval   = nH;
P.history = H(1:nH, :);
P.wallRootS = toc(tStart) - P.wallWalkS;
if ~found
    if isempty(why)
        why = sprintf('root not found (maxEval %d reached, min |g| %.2e)', ...
                      maxEval, min(abs(Ga(1:nE))));
    end
    P.reason = why;
    P.wallS = toc(tStart);
    return
end

% ---- arrival gate ----
F = fly_transfer(zr, rv0, rvfAt(x), L.phys);
P.flight  = F;
P.missKm  = F.flyKm;
P.missMps = F.flyVms;
P.res     = Ra(nE);
P.arrPhase = mod(cand.slotPhaseAtDep + x/L.tauD, 1);
P.rvf     = L.slotState(P.arrPhase);
failed = {};
if ~(F.flyKm < thr.flyKmMax && F.flyVms < thr.flyVmsMax), failed{end+1} = 'fly'; end
if ~(P.res < thr.residMax), failed{end+1} = 'resid'; end
if isempty(failed)
    P.ok   = true;
    P.reason = 'accepted';
    P.z8   = zr;
    P.tf   = x;
    P.tArr = cand.tDep + x;
    P.familyChanged = abs(x - tfCell)/tfCell > familyFrac;
else
    P.reason = ['gate failed: ' strjoin(failed, ', ')];
    P.arrPhase = NaN;
end
P.wallS = toc(tStart);
if verbose
    fprintf(['polish_gto_dro cell (%d,%d) sheet %d -> %.4g deg: ok %d (%s)\n' ...
             '  tf %.9f ND (cell %.9f)  arrPhase %.6f  miss %.3g km %.3g m/s  res %.2e\n' ...
             '  walk %d steps %.1f s, root %d evals %.1f s\n'], cand.iD, cand.iA, ...
            cand.iS, cand.orientDeg, P.ok, P.reason, x, tfCell, ...
            mod(cand.slotPhaseAtDep + x/L.tauD, 1), P.missKm, P.missMps, P.res, ...
            nW - 1, P.wallWalkS, P.nEval, P.wallRootS);
end
end

% ------------------------------------------------------------------------
function [z, ok, its, normR] = inner_solve(rv0, rvf, zw, rvw, L, Kseg)
% INNER_SOLVE  Minimum-time z8 from rv0 to rvf: ms_tfmin seeded by
%   seed_from_z8 of the warm z8 flown from ITS OWN departure state rvw;
%   counted only if converged and the flight misses rvf by < 100 km and
%   < 10 m/s (the leg-2 inner rule, fixed).
%   INPUTS: rv0, rvf [6 x 1]; zw [8 x 1]; rvw [6 x 1]; L leg1_inputs;
%   Kseg segments.  OUTPUTS: z [8 x 1] (NaN on throw); ok logical; its
%   iterations; normR final residual.
z = nan(8, 1);  ok = false;  its = NaN;  normR = NaN;
try
    seed = seed_from_z8(zw, rvw, Kseg, L.phys.Tnd, L.phys.cnd, L.muStar);
    [z, info] = ms_tfmin(rv0, rvf, seed, L.phys.Tnd, L.phys.cnd, L.muStar, ...
                         struct('tolR', 3e-11, 'wallSec', 300));
    its = info.iters;  normR = info.normR;
    if info.converged && all(isfinite(z)) && z(8) > 0
        F  = fly_transfer(z, rv0, rvf, L.phys);
        ok = F.flyKm < 100 && F.flyVms < 10;
    end
catch
    ok = false;
end
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
% GATE_THRESHOLDS  Arrival-gate thresholds; the override replaces named
%   fields and refuses unknown names.
%   INPUTS: override struct.  OUTPUTS: thr struct.
thr = struct('flyKmMax', 100, 'flyVmsMax', 10, 'residMax', 1e-8);
fn = fieldnames(override);
for kk = 1:numel(fn)
    if ~isfield(thr, fn{kk})
        error('polish_gto_dro:gateOverride', 'unknown threshold "%s"', fn{kk});
    end
    thr.(fn{kk}) = override.(fn{kk});
end
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
