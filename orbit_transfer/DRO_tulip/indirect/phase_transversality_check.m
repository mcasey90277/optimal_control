function X = phase_transversality_check(catMat, opts)
%% Purpose:
%
%   CROSS-CHECK X3 for a phase catalog: do the stored costates predict how
%   the minimum flight time changes with the two phases? Three things, from
%   the catalog file alone:
%
%     1. SENSITIVITIES at every entry: dT/ds_D = +lam(0) . x_D'(s_D) and
%        dT/ds_A = -lam(t_f) . x_A'(s_A)  (costate_common/phase_sensitivity).
%        One flight per entry; about ten seconds for 576 entries.
%     2. THE FINITE-DIFFERENCE TEST on a sample of entries: re-solve the
%        transfer at s +/- delta in each phase and difference the two flight
%        times. The OBSERVABLE does not use the returned costates; the re-solve
%        itself is an indirect solve seeded from the stored ones, and each
%        derivative tests ONE scalar projection of an endpoint costate. So it
%        is a consistency test, not an independent oracle (that is the
%        single-integrator test of phase_sensitivity). A converged re-solve is
%        assumed to have stayed on the branch; that is not checked. Gate:
%        |formula - FD| <= .absTol + .relTol max(|.|);
%        the verdict is PASS / FAIL / UNRESOLVED (a re-solve that does not
%        converge is not evidence about a costate).
%     3. THE EDGE MAP in departure phase: along every rib, the trapezoid
%        residual between neighbouring cells (phase_edge_residuals). Small is
%        CONSISTENT with one smooth branch of flight times, large suggests a
%        jump -- a heuristic in both directions (see section 3b in the code).
%        Reported, not gated. With it: time-consistent edges, the candidate
%        components they join (through s_D only), and the COSTATE jump across
%        each edge.
%
%   It also lists the cells nearest FIRST-ORDER STATIONARITY in both phases
%   (smallest |grad T|): the library's entries are minima at fixed phases,
%   and these are its candidates for the orbit-to-orbit minimum.
%
%  ASSUMPTIONS / NOTES:
%
% • THE ARRIVAL-PHASE EDGE MAP IS NOT JUDGED. At 24 arrival phases the flight
%   time is under-resolved along s_A (a 7-petal tulip puts structure at a
%   seventh of a period): measured median residual 105 min, against 0.5 min
%   along s_D. The exact test shows the arrival sensitivity itself is right
%   to 4-5 digits, so the residual is the grid's, not the costates'. It is
%   returned (.edgeA) for whoever refines the grid.
% • The exact test costs four polishes per sampled entry (about a minute).
% • One-sheet DRO -> tulip phase catalogs (schema v2), as the audit.
%
%% Inputs:
%
%  catMat                   char                    the catalog .mat
%  opts                     struct (optional)
%   .exactCells             [k x 2]                 (iD, iA) cells for the exact
%                                                   test; default a fixed-seed
%                                                   random sample of .nExact
%   .nExact                 int                     sample size [6]
%   .delta                  double                  phase offset of the re-solve [2e-4]
%   .relTol                 double                  exact-test gate [1e-3]
%   .consistentMinutes      double                  an edge within this is
%                                                   time-consistent [5]
%   .absTol                 double                  absolute part of the finite-
%                                                   difference gate [1e-5]
%   .jumpMinutes            double                  an s_D edge above this is a
%                                                   JUMP [10]
%   .print                  logical                 print the report [true]
%   .out                    char                    save X here ['']
%
%% Outputs:
%
%  X                        struct                  .dTf_dsD .dTf_dsA [nD x nA,
%                                                   ND time per unit phase]
%                                                   .edgeD .edgeA [minutes]
%                                                   .timeConsistentD/.A (per
%                                                   edge) .component .nComponent
%                                                   (phase_components)
%                                                   .costateJumpD
%                                                   .sameFamD .nEdgeD
%                                                   .nEdgeDJump .jumps (table)
%                                                   .exact (struct array)
%                                                   .verdictExact .okExact
%                                                   .stationary (table)
%
%% Revision History:
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

% TEST SEAM: handles to this file's local functions (tests/test_phase_transversality_check)
if ischar(catMat) && strcmp(catMat, 'localfunctions'), X = localHandles(localfunctions);  return, end
if nargin < 2, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
delta = pick(opts, 'delta', 2e-4);        relTol = pick(opts, 'relTol', 1e-3);
jumpMinutes = pick(opts, 'jumpMinutes', 10);   say = pick(opts, 'print', true);
consistentMinutes = pick(opts, 'consistentMinutes', 5);   absTol = pick(opts, 'absTol', 1e-5);

%% 1. The catalog, and the endpoint closures of ITS problem:
L = load(catMat);  fn = fieldnames(L);  c = L.(fn{1});
assert(isscalar(c.sheets), 'phase_transversality_check: a one-sheet phase catalog is expected');
sh = c.sheets(1);
sD = sh.sD_frac(:).';  sA = sh.sA_frac(:).';  nD = numel(sD);  nA = numel(sA);
has = logical(sh.has_solution(:, :, 1));   TF = sh.tf_nd(:, :, 1);   TF(~has) = NaN;
tStar = c.constants.tStar_s;
[B, ~] = arclength_arrival('setup', struct('thrustN', c.rungs_N(1), 'ispS', c.thruster.isp_s, 'm0kg', c.thruster.m0_kg, ...
             'tauDRO', sh.tauDRO, 'NpTulip', sh.Np, 'pmTulip', sh.pm, 'sD', sD(1), 'physicsOnly', true));

%% 2. The two sensitivities at every entry (one flight each):
GD = nan(nD, nA);  GA = nan(nD, nA);
for iD = 1:nD
    rv0 = B.stateD(sD(iD));  rv0 = rv0(1:6);
    dxD = departureTangent(B, sD(iD));
    for iA = find(has(iD, :))
        z8 = sh.z8(:, sh.entry_index(iD, iA, 1));
        [~, Y] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0; 1; z8(1:7)], B.Tnd, B.cnd, B.mu);
        dxA = B.dstateA(sA(iA));
        [GD(iD, iA), GA(iD, iA)] = phase_sensitivity(z8(1:6), Y(end, 8:13).', dxD, dxA(1:6));
    end
end

%% 3. The edge maps (minutes), and the jumps along departure phase:
toMin = tStar/60;
edgeD = phase_edge_residuals(TF, GD, sD, 1)*toMin;
edgeA = phase_edge_residuals(TF, GA, sA, 2)*toMin;
sameFamD = false(nD, nA);
if isfield(sh, 'family_index')
    FI = double(sh.family_index(:, :, 1));
    sameFamD = FI == circshift(FI, -1, 1);
end
[jD, jA] = find(abs(edgeD) > jumpMinutes);
nextD = mod(jD, nD) + 1;
jumps = table(jD, nextD, jA, sD(jD).', sD(nextD).', sA(jA).', edgeD(sub2ind([nD nA], jD, jA)), sameFamD(sub2ind([nD nA], jD, jA)), ...
              'VariableNames', {'iD', 'iDnext', 'iA', 'sD', 'sDnext', 'sA', 'residualMin', 'sameFamily'});
jumps = sortrows(jumps, 'residualMin', 'descend', 'ComparisonMethod', 'abs');

%% 3b. Time-consistent edges, candidate components, and the costate jump
% An edge is TIME-CONSISTENT when its trapezoid residual is within
% consistentMinutes: the two flight times and the two phase sensitivities fit
% one smooth branch of flight times. That is a HEURISTIC (review 2026-09-19):
% two different branches can give a small residual (equal times and slopes; a
% jump cancelled by the slope term), one branch can give a large one (large
% T''', a coarse grid, a fold), and it tests one scalar, not the costates. So:
%  - nothing here is called "safe"; whether a guess interpolated across an edge
%    is USABLE is an experiment (polish it), not a residual;
%  - the costate jump across each departure-phase edge is reported beside it
%    (.costateJumpD, relative): on the 70 mN library it is ~30% even on
%    time-consistent edges, so linear interpolation of costates is crude;
%  - the under-resolved arrival axis is NOT used to join cells: an accidentally
%    small residual there must not merge unrelated components.
timeConsistentD = abs(edgeD) <= consistentMinutes;   timeConsistentA = abs(edgeA) <= consistentMinutes;   % NaN is not consistent
[component, nComponent] = phase_components(timeConsistentD, false(nD, nA), has);
Z = nan(7, nD, nA);
for iD = 1:nD, for iA = find(has(iD, :)), Z(:, iD, iA) = sh.z8(1:7, sh.entry_index(iD, iA, 1)); end, end
costateJumpD = squeeze(sqrt(sum((circshift(Z, -1, 2) - Z).^2, 1))./sqrt(sum(Z.^2, 1)));

%% 4. The exact test: re-solve at s +/- delta, never reading a costate:
cells = pick(opts, 'exactCells', []);
if isempty(cells)
    [aD, aA] = find(has);  rs = RandStream('mt19937ar', 'Seed', 20260918);        % the same sample every run
    k = randperm(rs, numel(aD), min(pick(opts, 'nExact', 6), numel(aD)));
    cells = [aD(k), aA(k)];
end
exact = struct('iD', {}, 'iA', {}, 'formulaD', {}, 'exactD', {}, 'relErrD', {}, 'formulaA', {}, 'exactA', {}, 'relErrA', {}, 'resolved', {}, ...
               'stepD', {}, 'stepA', {}, 'absErrD', {}, 'absErrA', {});
for q = 1:size(cells, 1)
    iD = cells(q, 1);  iA = cells(q, 2);
    assert(has(iD, iA), 'phase_transversality_check: cell (%d,%d) holds no entry', iD, iA);
    z8 = sh.z8(:, sh.entry_index(iD, iA, 1));
    [exD, exA, stepD, stepA] = exactDerivatives(B, z8, sD(iD), sA(iA), delta);
    e = struct('iD', iD, 'iA', iA, 'formulaD', GD(iD, iA), 'exactD', exD, 'relErrD', relErr(GD(iD, iA), exD, absTol/relTol), ...
               'formulaA', GA(iD, iA), 'exactA', exA, 'relErrA', relErr(GA(iD, iA), exA, absTol/relTol), 'resolved', isfinite(exD) && isfinite(exA), ...
               'stepD', stepD, 'stepA', stepA, 'absErrD', abs(GD(iD, iA) - exD), 'absErrA', abs(GA(iD, iA) - exA));
    exact(end+1) = e; %#ok<AGROW>
end
verdictExact = exactVerdict(exact, relTol);
okExact = strcmp(verdictExact, 'PASS');

%% 5. The cells nearest first-order stationarity in BOTH phases:
gradNorm = hypot(GD, GA);
[gs, order] = sort(gradNorm(:));  order = order(isfinite(gs));  order = order(1:min(10, end));
[oD, oA] = ind2sub([nD nA], order);
stationary = table(oD, oA, sD(oD).', sA(oA).', TF(order)*tStar/86400, GD(order)*tStar/86400, GA(order)*tStar/86400, gradNorm(order), ...
                   'VariableNames', {'iD', 'iA', 'sD', 'sA', 'tfDays', 'dTfDays_dsD', 'dTfDays_dsA', 'gradNorm'});

X = struct('catMat', catMat, 'dTf_dsD', GD, 'dTf_dsA', GA, 'edgeD', edgeD, 'edgeA', edgeA, 'sameFamD', sameFamD, ...
           'nEdgeD', nnz(isfinite(edgeD)), 'nEdgeDJump', height(jumps), 'jumpMinutes', jumpMinutes, 'jumps', jumps, ...
           'consistentMinutes', consistentMinutes, 'timeConsistentD', timeConsistentD, 'timeConsistentA', timeConsistentA, ...
           'component', component, 'nComponent', nComponent, 'costateJumpD', costateJumpD, 'absTol', absTol, ...
           'exact', exact, 'relTol', relTol, 'delta', delta, 'verdictExact', verdictExact, 'okExact', okExact, 'stationary', stationary, 'when', char(datetime('now')));
if say, printReport(X, TF, tStar); end
if ~isempty(pick(opts, 'out', '')), save(opts.out, 'X'); end
end

% ==========================================================================
function dxD = departureTangent(B, s)
% DEPARTURETANGENT  d(departure state)/d(phase) of the closure the solver
% uses, by central difference (the setup exposes the arrival derivative
% B.dstateA but not the departure one; the closure is a smooth periodic
% spline, so a 1e-6 step is accurate to ~1e-10).  INPUTS: B; s.
% OUTPUTS: dxD [6 x 1].
e = 1e-6;
a = B.stateD(s + e);  b = B.stateD(s - e);
dxD = (a(1:6) - b(1:6))/(2*e);
dxD = dxD(:);
end

% ==========================================================================
function [exD, exA, stepD, stepA] = exactDerivatives(B, z8, sD, sA, delta)
% EXACTDERIVATIVES  dT/ds_D and dT/ds_A by RE-SOLVING the transfer with one
% phase moved to s - d and s + d (ms_tfmin, seeded from this entry) and
% differencing the two flight times. No costate is read. Each derivative is
% tried at d = delta with a 120 s budget, then at delta/2 with 300 s; NaN
% when neither attempt converges on both sides (UNRESOLVED, not a failure).
% INPUTS: B; z8; sD; sA; delta.  OUTPUTS: exD; exA.
rv0 = pickState(B.stateD(sD));   rvf = pickState(B.stateA(sA));
seed = seed_from_z8(z8, rv0, 24, B.Tnd, B.cnd, B.mu);
ladder = [delta, 120; delta/2, 300];                 % [phase offset, wall seconds]
exD = NaN;  exA = NaN;  stepD = NaN;  stepA = NaN;
for a = 1:size(ladder, 1)
    d = ladder(a, 1);  wall = ladder(a, 2);
    if ~isfinite(exA)
        tm = resolveTf(B, rv0, pickState(B.stateA(sA - d)), seed, wall);
        tp = resolveTf(B, rv0, pickState(B.stateA(sA + d)), seed, wall);
        exA = (tp - tm)/(2*d);                       % NaN if either side is NaN
        if isfinite(exA), stepA = d; end
    end
    if ~isfinite(exD)
        tm = resolveTf(B, pickState(B.stateD(sD - d)), rvf, seed, wall);
        tp = resolveTf(B, pickState(B.stateD(sD + d)), rvf, seed, wall);
        exD = (tp - tm)/(2*d);
        if isfinite(exD), stepD = d; end
    end
end
end

function tf = resolveTf(B, rv0, rvf, seed, wallSec)
% RESOLVETF  The minimum time between two given states, by polishing this
% entry's own multiple-shooting seed onto them; NaN if it does not converge.
% INPUTS: B; rv0; rvf [6 x 1]; seed; wallSec.  OUTPUTS: tf (ND) | NaN.
seed.Y(1:6, 1) = rv0;
[z, it] = ms_tfmin(rv0, rvf, seed, B.Tnd, B.cnd, B.mu, struct('tolR', 3e-11, 'wallSec', wallSec));
if it.converged, tf = z(8); else, tf = NaN; end
end

function x = pickState(x)
% PICKSTATE  Position and velocity of a closure's output, as a column.
% INPUTS: x.  OUTPUTS: x [6 x 1].
x = x(1:6);  x = x(:);
end

function v = exactVerdict(exact, relTol)
% EXACTVERDICT  The exact test's verdict, THREE-valued:
%   FAIL        some RESOLVED derivative disagrees with the formula
%   UNRESOLVED  none disagrees, but a re-solve did not converge (or nothing
%               was sampled): that says nothing about a costate, and it is
%               never counted as a pass either
%   PASS        every sampled derivative resolved and agrees
% INPUTS: exact (struct array: .relErrD .relErrA, NaN = unresolved); relTol.
% OUTPUTS: v char.
errs = [[exact.relErrD], [exact.relErrA]];
if any(errs > relTol),                    v = 'FAIL';
elseif isempty(errs) || any(isnan(errs)), v = 'UNRESOLVED';
else,                                     v = 'PASS';
end
end

% ==========================================================================
function printReport(X, TF, tStar)
% PRINTREPORT  The check as text, study-script style: value / threshold
% PASS|FAIL for the gate, plain numbers for what is reported.
% INPUTS: X; TF; tStar.  OUTPUTS: none.
fprintf('PHASE-TRANSVERSALITY CROSS-CHECK (X3)  %s\n', X.catMat);
fprintf('  sensitivities : %d entries; |dT/ds_D| median %.3f d per unit phase, |dT/ds_A| median %.3f\n', nnz(isfinite(X.dTf_dsD)), ...
        median(abs(X.dTf_dsD(:)), 'omitnan')*tStar/86400, median(abs(X.dTf_dsA(:)), 'omitnan')*tStar/86400);
fprintf('  X3 finite-difference test : formula against re-solves at s +/- %.0e (the observable uses no returned costate)\n', X.delta);
for e = X.exact
    fprintf('     cell (%2d,%2d)  dT/ds_D %+.6f vs %+.6f (rel %.1e)   dT/ds_A %+.6f vs %+.6f (rel %.1e)%s\n', e.iD, e.iA, ...
            e.formulaD, e.exactD, e.relErrD, e.formulaA, e.exactA, e.relErrA, passIf(e.resolved, '', '   UNRESOLVED'));
end
errs = [[X.exact.relErrD], [X.exact.relErrA]];
fprintf('     worst relative error %.1e / %.0e over %d resolved derivative(s), %d unresolved   %s\n', ...
        max([errs(isfinite(errs)), 0]), X.relTol, nnz(isfinite(errs)), nnz(~isfinite(errs)), X.verdictExact);
rD = abs(X.edgeD(isfinite(X.edgeD)));  rA = abs(X.edgeA(isfinite(X.edgeA)));
fprintf('  edges in s_D  : %d; residual median %.2f min, 90%% %.2f min; %d above %g min (JUMPS between neighbours; %d of them inside one family label)\n', ...
        numel(rD), median(rD), prctile(rD, 90), X.nEdgeDJump, X.jumpMinutes, nnz(X.jumps.sameFamily));
fprintf('  edges in s_A  : %d; residual median %.0f min -- NOT JUDGED: the sheet is under-resolved along s_A at this spacing\n', numel(rA), median(rA));
for q = 1:min(5, height(X.jumps))
    j = X.jumps(q, :);
    fprintf('     jump: column %2d (sA %.4f), sD %.4f -> %.4f: %.0f min%s\n', j.iA, j.sA, j.sD, j.sDnext, j.residualMin, passIf(j.sameFamily, '  [same family label]', ''));
end
sz = arrayfun(@(b) nnz(X.component == b), 1:X.nComponent);  jt = X.costateJumpD(X.timeConsistentD);
fprintf('  time-consistent edges (residual <= %g min; a heuristic, not an interpolation licence): %d of %d in s_D, %d of %d in s_A\n', ...
        X.consistentMinutes, nnz(X.timeConsistentD), numel(X.timeConsistentD), nnz(X.timeConsistentA), numel(X.timeConsistentA));
fprintf('  components    : %d, joined through s_D edges only; the largest hold %s cells\n', X.nComponent, mat2str(sz(1:min(5, end))));
fprintf('  costate jump  : across a time-consistent s_D edge the costates change by %.0f%% (median), %.0f%% (90%%)\n', 100*median(jt), 100*prctile(jt, 90));
[tfMin, kMin] = min(TF(:));  [mD, mA] = ind2sub(size(TF), kMin);
fprintf('  stationarity  : the fastest entry (%d,%d), %.3f d, has dT/ds = (%+.3f, %+.3f) d per unit phase\n', mD, mA, tfMin*tStar/86400, ...
        X.dTf_dsD(kMin)*tStar/86400, X.dTf_dsA(kMin)*tStar/86400);
s = X.stationary(1, :);
fprintf('                  the entry nearest stationarity in both phases: (%d,%d), %.3f d, dT/ds = (%+.3f, %+.3f)\n', s.iD, s.iA, s.tfDays, s.dTfDays_dsD, s.dTfDays_dsA);
end

% ==========================================================================
function r = relErr(a, b, floor_)
% RELERR  The error of a against b, scaled so that  r <= relTol  means
% |a - b| <= absTol + relTol * max(|a|, |b|)  when floor_ = absTol/relTol. A
% purely relative error is meaningless near a STATIONARY phase, where the
% derivative itself goes to zero -- the very points this tool looks for.
% NaN if either is not finite.  INPUTS: a; b; floor_.  OUTPUTS: r.
if isfinite(a) && isfinite(b), r = abs(a - b)/(max(abs(a), abs(b)) + floor_); else, r = NaN; end
end

function t = passIf(c, a, b)
% PASSIF  a when c, else b.  INPUTS: c; a; b.  OUTPUTS: t.
if c, t = a; else, t = b; end
end

function v = pick(s, f, d_)
% PICK  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end

function H = localHandles(fh)
% LOCALHANDLES  This file's local functions as a struct of handles keyed by
% name -- the TEST SEAM.  INPUTS: fh (cell of handles).  OUTPUTS: H struct.
H = struct();
for k = 1:numel(fh), H.(func2str(fh{k})) = fh{k}; end
end
