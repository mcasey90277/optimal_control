function ok = test_dro_elfo_library()
%% Purpose:
%
%   Test of the DRO -> ELFO library functions (dro_elfo_describe,
%   dro_elfo_transfers_at, dro_elfo_fly, dro_elfo_extremes,
%   dro_elfo_torus_plot, dro_elfo_movie) on the REAL pilot catalog
%   results/pilot_chain/costate_catalog_dro_elfo.mat (one sheet, read
%   only) and a synthetic 2-sheet catalog built IN MEMORY from it (sheet 2
%   = sheet 1 at raanDeg 30 with every t_f scaled by 1.1, so a wrong
%   routing shows). About 20 s. Without the pilot catalog (a fresh clone:
%   results/ is not tracked) the test prints SKIP and returns true.
%     1  describe: counts equal counts taken another way (z8 columns, a
%        loop over cells, the conjugate verdicts entry by entry, per-rung
%        t_f median from the entries' own z8)
%     2  transfers_at: exactly the rungs with has_solution at the cell, high
%        thrust first, Delta-V and propellant by the derive formulas
%        computed here; wrapped phases; an empty cell gives a 0 x 1 struct
%        with the same fields
%     3  refusals by name: an unknown sheet errors :noSheet in
%        transfers_at, fly and torus_plot; extremes :noMatch and :filter;
%        torus_plot :noRung
%     4  the 2-sheet catalog routes to the right sheet (raanDeg 30 and 390
%        give sheet 2 with t_f x 1.1; raanDeg 0 gives sheet 1)
%     5  fly: a 1 N entry arrives within 1 km, its located minimum Moon
%        altitude is at least 500 km, mass and Delta-V agree with the
%        catalog formulas; the plot is drawn invisibly and closed
%     6  extremes equal a brute-force min / max over every entry for each
%        filter field and both metrics; conjPassOnly skips a refuted entry
%        planted as the fastest
%     7  torus_plot: NaN exactly where no entry, values = t_f in days,
%        verdicts = conj_pass; the 'all' mode stacks both sheets
%     8  movie: 4 frames, a non-empty .mp4 of 1280 x 720 and a .gif, for one
%        and for two entries (SKIP, not FAIL, without SatelliteAnimator)
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 true iff every check
%                                                   passed (throws on a
%                                                   failure)
%
%% Revision History:
%  M. Casey                                                   (c) 11/21/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here    = fileparts(mfilename('fullpath'));
oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(fullfile(fileparts(here), 'catalog'));
catMat = fullfile(fileparts(here), 'results', 'pilot_chain', 'costate_catalog_dro_elfo.mat');
if ~isfile(catMat)
    fprintf('test_dro_elfo_library: SKIP -- no pilot catalog at %s\n', catMat);
    ok = true;  return
end
L  = load(catMat);
C  = L.costate_catalog_dro_elfo;
C2 = C;
C2.sheets(2) = C.sheets(1);
C2.sheets(2).sheet_key.raanDeg  = 30;
C2.sheets(2).arr_params.raanDeg = 30;
C2.sheets(2).z8(8, :) = 1.1*C2.sheets(2).z8(8, :);
C2.sheets(2).tf_nd    = 1.1*C2.sheets(2).tf_nd;
sh    = C.sheets(1);
rungs = C.rungs_N;  nR = numel(rungs);
[nD, nA, ~] = size(sh.has_solution);
tStar = C.constants.tStar_s;  lStar = C.constants.lStar_km;
m0kg  = C.thruster.m0_kg;     cnd   = C.thruster.c_nd;
dvOf  = @(TN, tf) cnd*log(1./(1 - ((TN/m0kg)*tStar^2/(lStar*1000))*tf/cnd))*lStar/tStar;

%% 1. describe
D = dro_elfo_describe(C);
nFull = 0;  nLow = 0;  cp = [0 0 0];
for iD = 1:nD
    for iA = 1:nA
        hk = squeeze(sh.has_solution(iD, iA, :));
        nFull = nFull + all(hk);
        nLow  = nLow + hk(end);
        for kr = find(hk).'
            v = double(sh.conj_pass(iD, iA, kr));
            if v == 1, cp(1) = cp(1) + 1; elseif v == 0, cp(2) = cp(2) + 1; else, cp(3) = cp(3) + 1; end
        end
    end
end
assert(D.nEntries == size(sh.z8, 2) && D.nEntries == C.n_entries, 'describe entries %d vs z8 columns %d', D.nEntries, size(sh.z8, 2));
assert(D.sheetFullLadders == nFull && D.sheetLowestRungCells == nLow, 'describe full ladders / lowest-rung cells');
assert(isequal(D.conj, cp), 'describe conjugate counts %s vs %s', mat2str(D.conj), mat2str(cp));
for kr = 1:nR
    cols = sh.entry_index(:, :, kr);  cols = cols(sh.has_solution(:, :, kr));
    tfd  = sh.z8(8, cols)*tStar/86400;
    assert(D.rungEntries(kr) == numel(cols), 'describe rung %d entries', kr);
    assert(abs(D.rungTfDays(kr, 2) - median(tfd)) < 1e-12 && abs(D.rungDvKms(kr, 2) - max(dvOf(rungs(kr), sh.z8(8, cols)))) < 1e-12, ...
           'describe rung %d t_f median / Delta-V max', kr);
end
D2 = dro_elfo_describe(C2);
assert(D2.nSheets == 2 && D2.nEntries == 2*D.nEntries && isequal(D2.raanDeg, [0 30]), 'describe on 2 sheets');
fprintf('  1 describe: %d entries, %d full ladders, conj %s -- PASS\n', D.nEntries, nFull, mat2str(cp));

%% 2. transfers_at
cellsTried = 0;
for iD = 1:3:nD
    for iA = 2:3:nA
        T = dro_elfo_transfers_at(C, 1, 0, sh.sD_frac(iD) + 1, sh.sA_frac(iA) - 1 + 1e-4);   % wrapped, slightly off-grid
        krWant = find(squeeze(sh.has_solution(iD, iA, :))).';
        assert(isequal(sort([T.kr]), krWant), 'cell (%d,%d): rungs %s vs %s', iD, iA, mat2str([T.kr]), mat2str(krWant));
        assert(issorted(-[T.thrustN]), 'rows high to low thrust');
        for kk = 1:numel(T)
            tf = sh.z8(8, sh.entry_index(iD, iA, T(kk).kr));
            mf = 1 - ((T(kk).thrustN/m0kg)*tStar^2/(lStar*1000))*tf/cnd;
            assert(T(kk).iD == iD && T(kk).iA == iA && T(kk).tf_nd == tf && T(kk).sheet == 1, 'row cell / tf');
            assert(abs(T(kk).dvKms - cnd*log(1/mf)*lStar/tStar) < 1e-12 && abs(T(kk).propellantKg - (1 - mf)*m0kg) < 1e-10, ...
                   'Delta-V / propellant formula');
            assert(abs(T(kk).tfDays - tf*tStar/86400) < 1e-12, 'tfDays');
        end
        cellsTried = cellsTried + 1;
    end
end
Ce = C;  Ce.sheets(1).has_solution(2, 2, :) = false;
Te = dro_elfo_transfers_at(Ce, 1, 0, sh.sD_frac(2), sh.sA_frac(2));
Tf = dro_elfo_transfers_at(C, 1, 0, sh.sD_frac(1), sh.sA_frac(1));
assert(isempty(Te) && isstruct(Te) && isequal(fieldnames(Te), fieldnames(Tf)), 'empty cell: 0 x 1 struct, same fields');
fprintf('  2 transfers_at: %d cells checked, empty cell gives %d x %d struct -- PASS\n', cellsTried, size(Te, 1), size(Te, 2));

%% 3. refusals by name
expect_err(@() dro_elfo_transfers_at(C2, 1, 45, 0, 0), 'dro_elfo_transfers_at:noSheet', 'raanDeg 30');
expect_err(@() dro_elfo_transfers_at(C, 1 + 1e-6, 0, 0, 0), 'dro_elfo_transfers_at:noSheet', 'depSize 1');
eBad = Tf(1);  eBad.sheet_key.raanDeg = 45;
expect_err(@() dro_elfo_fly(C, eBad), 'dro_elfo_fly:noSheet', 'raanDeg 0');
expect_err(@() dro_elfo_torus_plot(C, 45, 1, struct('visible', false)), 'dro_elfo_torus_plot:noSheet', 'raanDeg 0');
expect_err(@() dro_elfo_torus_plot(C, 0, 4, struct('visible', false)), 'dro_elfo_torus_plot:noRung', '0.5');
expect_err(@() dro_elfo_extremes(C, struct('raanDeg', 45)), 'dro_elfo_extremes:noMatch', 'raanDeg');
expect_err(@() dro_elfo_extremes(C, struct('thrust', 1)), 'dro_elfo_extremes:filter', 'thrust');
fprintf('  3 refusals -- PASS\n');

%% 4. routing on the 2-sheet catalog
T30  = dro_elfo_transfers_at(C2, 1, 30, 0.25, 0.5);
T390 = dro_elfo_transfers_at(C2, 1, 390, 0.25, 0.5);
T0   = dro_elfo_transfers_at(C2, 1, 0, 0.25, 0.5);
assert(all([T30.sheet] == 2) && all([T30.raanDeg] == 30) && isequal([T30.tf_nd], [T390.tf_nd]), '30 and 390 route to sheet 2');
assert(all([T0.sheet] == 1) && max(abs([T30.tf_nd] - 1.1*[T0.tf_nd])) < 1e-14, 'sheet 2 carries t_f x 1.1, sheet 1 does not');
assert(T30(1).arr_params.raanDeg == 30 && T0(1).arr_params.raanDeg == 0, 'each row carries its sheet''s arr_params');
fprintf('  4 routing: raanDeg 30 / 390 -> sheet 2, 0 -> sheet 1 -- PASS\n');

%% 5. fly a 1 N entry
e1 = dro_elfo_pick(C, 1, 0, 0.25, 0.5, 1);
[F, hF] = dro_elfo_fly(C, e1, struct('plot', true, 'visible', false));
assert(isgraphics(hF) && strcmp(get(hF, 'Visible'), 'off'), 'invisible figure');
close(hF);
assert(F.flyKm < 1, 'arrival miss %.3g km >= 1 km', F.flyKm);
assert(F.minMoonAltKm >= 500, 'minimum Moon altitude %.1f km < 500 km', F.minMoonAltKm);
assert(abs(F.dvKms - dvOf(1, e1.tf_nd)) < 1e-6 && abs(F.tfNd - e1.tf_nd) == 0, 'flown Delta-V vs formula');
assert(abs(F.massKg(1) - m0kg) < 1e-12 && all(abs(sqrt(sum(F.u.^2, 2)) - 1) < 1e-12), 'mass history, unit thrust direction');
T1 = dro_elfo_transfers_at(C2, 1, 0, e1.sD, e1.sA);
F2 = dro_elfo_fly(C2, T1(1));                          % a row from transfers_at flies too (top rung)
assert(F2.flyKm < 1 && F2.sheet == 1, 'a transfers_at row flies');
fprintf('  5 fly at 1 N: miss %.2e km, min Moon altitude %.0f km, %.2f revs -- PASS\n', F.flyKm, F.minMoonAltKm, F.moonRevs);

%% 6. extremes against brute force
Cx = C2;                                               % plant a refuted entry as the fastest of sheet 1
kRef = find(Cx.sheets(1).conj_pass == 0 & Cx.sheets(1).has_solution, 1);
[rD, rA, rR] = ind2sub(size(Cx.sheets(1).has_solution), kRef);
Cx.sheets(1).z8(8, Cx.sheets(1).entry_index(rD, rA, rR)) = 1e-3;
filters = {struct(), struct('thrustN', 1), struct('thrustN', [1 0.5]), struct('raanDeg', 30), ...
           struct('raanDeg', [0 30], 'thrustN', 2), struct('droTau', 1, 'raanDeg', 390), ...
           struct('conjPassOnly', true), struct('conjPassOnly', true, 'thrustN', rungs(rR), 'raanDeg', 0)};
for kf = 1:numel(filters)
    for metric = {'time', 'deltaV'}
        [eMin, eMax] = dro_elfo_extremes(Cx, filters{kf}, metric{1});
        [vLo, vHi]   = brute(Cx, filters{kf}, metric{1}, dvOf);
        if strcmp(metric{1}, 'time'), gLo = eMin.tf_nd; gHi = eMax.tf_nd; else, gLo = eMin.dvKms; gHi = eMax.dvKms; end
        assert(abs(gLo - vLo) < 1e-12 && abs(gHi - vHi) < 1e-12, 'filter %d %s: got [%g %g], brute [%g %g]', ...
               kf, metric{1}, gLo, gHi, vLo, vHi);
        for ee = [eMin eMax]
            f = filters{kf};
            if isfield(f, 'thrustN'), assert(any(abs(f.thrustN - ee.thrustN) < 1e-9), 'thrust filter honoured'); end
            if isfield(f, 'raanDeg'), assert(any(mod(f.raanDeg, 360) == ee.raanDeg), 'raanDeg filter honoured'); end
            if isfield(f, 'conjPassOnly'), assert(ee.conj_pass == 1, 'conjPassOnly honoured'); end
        end
    end
end
eAny = dro_elfo_extremes(Cx, struct(), 'time');
assert(eAny.conj_pass == 0 && eAny.tf_nd == 1e-3, 'the planted refuted entry is the fastest without conjPassOnly');
fprintf('  6 extremes: %d filters x 2 metrics equal brute force; planted refuted entry skipped -- PASS\n', numel(filters));

%% 7. torus plot
kr1 = find(rungs == 1);
H = dro_elfo_torus_plot(C, 0, 1, struct('visible', false));
close(H.fig);
has = sh.has_solution(:, :, kr1);
assert(isequal(isnan(H.tfDays), ~has) && isequal(isnan(H.verdict), ~has), 'NaN exactly where no entry');
cols = sh.entry_index(:, :, kr1);
assert(max(abs(H.tfDays(has) - sh.z8(8, cols(has)).'*tStar/86400)) < 1e-12, 'torus values = t_f [d]');
cpk = double(sh.conj_pass(:, :, kr1));
assert(isequal(H.verdict(has), cpk(has)), 'torus verdicts = conj_pass');
Ha = dro_elfo_torus_plot(C2, 'all', 0.5, struct('visible', false));
close(Ha.fig);
assert(isequal(size(Ha.tfDays), [nD nA 2]) && isequal(Ha.raanDeg, [0 30]) && ...
       max(abs(Ha.tfDays(:, :, 2) - 1.1*Ha.tfDays(:, :, 1)), [], 'all', 'omitnan') < 1e-12 && ...
       isequal(isnan(Ha.tfDays(:, :, 1)), ~sh.has_solution(:, :, end)), '''all'' mode stacks the sheets');
fprintf('  7 torus: %d empty cells NaN at 1 N, all-mode %s -- PASS\n', nnz(~has), mat2str(size(Ha.tfDays)));

%% 8. movie
if isempty(which('pumpkynPie.plot.SatelliteAnimator'))
    fprintf('  8 movie: SKIP -- pumpkynPie.plot.SatelliteAnimator not on the path\n');
else
    stem = fullfile(tempdir, sprintf('test_dro_elfo_movie_%d', round(1e6*rand)));
    ws = warning('off', 'dro_elfo_movie:noEphemData');
    [s1, i1] = dro_elfo_movie(C, e1, [stem '_1'], struct('nFrames', 4));
    [s2, i2] = dro_elfo_movie(C2, {T0(1), T30(1)}, [stem '_2'], struct('nFrames', 4, 'labels', {{'raan 0', 'raan 30'}}));
    warning(ws);
    for mv = {i1.mp4, i2.mp4}
        dd = dir(mv{1});
        assert(~isempty(dd) && dd.bytes > 0, 'movie %s written', mv{1});
        vr = VideoReader(mv{1});
        assert(vr.Width == 1280 && vr.Height == 720, 'frames %d x %d, want 1280 x 720', vr.Width, vr.Height);
    end
    assert(~isempty(s1) && ~isempty(s2) && isfile(i1.gif) && isfile(i2.gif) && i2.nFrames == 4, 'stems and gifs');
    delete([stem '_1.mp4'], [stem '_1.gif'], [stem '_2.mp4'], [stem '_2.gif']);
    fprintf('  8 movie: one- and two-panel 1280 x 720 .mp4 + .gif -- PASS\n');
end

ok = true;
fprintf('test_dro_elfo_library: PASS\n');
end

% ------------------------------------------------------------------------
function [vLo, vHi] = brute(C, f, metric, dvOf)
% BRUTE  Min and max of the metric over every entry, filter applied by
% plain loops.  INPUTS: C [catalog]; f [struct] filter; metric [char];
% dvOf [fhandle] (TN, tf) -> km/s.  OUTPUTS: vLo, vHi [1 x 1].
vLo = inf;  vHi = -inf;
for ks = 1:numel(C.sheets)
    sh = C.sheets(ks);
    if isfield(f, 'droTau') && sh.sheet_key.depSize ~= f.droTau, continue, end
    if isfield(f, 'raanDeg') && ~any(mod(f.raanDeg, 360) == sh.sheet_key.raanDeg), continue, end
    for kr = 1:numel(C.rungs_N)
        if isfield(f, 'thrustN') && ~any(f.thrustN == C.rungs_N(kr)), continue, end
        for iD = 1:size(sh.has_solution, 1)
            for iA = 1:size(sh.has_solution, 2)
                if ~sh.has_solution(iD, iA, kr), continue, end
                if isfield(f, 'conjPassOnly') && f.conjPassOnly && sh.conj_pass(iD, iA, kr) ~= 1, continue, end
                tf = sh.z8(8, sh.entry_index(iD, iA, kr));
                if strcmp(metric, 'time'), v = tf; else, v = dvOf(C.rungs_N(kr), tf); end
                vLo = min(vLo, v);  vHi = max(vHi, v);
            end
        end
    end
end
end

function expect_err(fh, id, msgPart)
% EXPECT_ERR  fh must throw identifier id, its message containing msgPart.
% INPUTS: fh [1 x 1] function handle; id [1 x n] char; msgPart [1 x m] char.
% OUTPUTS: none (throws on a mismatch).
try
    fh();
catch err
    assert(strcmp(err.identifier, id), 'expected %s, got %s (%s)', id, err.identifier, err.message);
    assert(contains(err.message, msgPart), 'message of %s lacks "%s": %s', id, msgPart, err.message);
    fprintf('    refused as expected: %s\n', err.identifier);
    return
end
error('test_dro_elfo_library:noThrow', 'expected %s, nothing thrown', id);
end
