function ok = test_mission_gdop()
%% Purpose:
%
%   Verifies mission_grid and mission_gdop: (1) the two grids; (2) the
%   stitch of per-interval engine calls when satellites come on station
%   (no duplicated or dropped epoch, nOn steps 4 -> 5 -> 6, three runs,
%   run 2 equal to a direct call on its satellite subset);
%   (3) equivalence with a direct gdop_grid_from_states call on the same
%   epochs (the wrapper adds nothing), for one-call, time-blocked and
%   one- and two-epoch runs; (4) fewer than four on station: NaN GDOP and
%   the visibility count of the three, no error; runs with 0, 1 and 2 on
%   station (NaN GDOP, the subset's visibility, R1's vis >= 2 flag); a
%   scored window shorter than a day is refused; a phase/tArr
%   inconsistency is refused; an H whose on-station position is moved
%   1e-6 ND off its slot is refused (mission_gdop:slot, the deliberate
%   failure of that check); (5) the validation seam on paper chunk 5
%   against the saved Table 1 column; (6) UNRESOLVED for a missing
%   reference file, one without S2, a reference phasing (refPhasing seam)
%   or an M.dtau6 that differs from the pinned identity; (7) deliberate failure: one
%   satellite moved 0.3 ND off its slot turns the validation to FAIL.
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
%  M. Casey                                                   (c) 11/10/2026
%  M. Casey  mission_gdop:slot deliberate failure              (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_lunar_pnt_paths.m'));

M   = mission_inputs();
day = 86400/M.tStar;                                   % one day [ND]

% ---- 1. grids ----
gp = mission_grid('paper');
want = {[-89.9 -87], [-86.5 -84], [-83.5 -81], [-80.5 -78], [-77.5 -75]};
assert(numel(gp.chunks) == 5 && isequal(gp.chunks, want), 'paper chunks');
nP = 0;
for cc = 1:5
    nP = nP + numel(gp.chunks{cc}(1):gp.lat_res:gp.chunks{cc}(2)) * ...
              numel(gp.lon_range(1):gp.lon_res:gp.lon_range(2));
end
assert(nP == 10770 && gp.lat_res == 0.5 && gp.lon_res == 1 && ...
       isequal(gp.lon_range, [0 358]) && gp.dtSec == 348.1 && gp.mask == 0, ...
       'paper grid values (%d points)', nP);
gc = mission_grid('coarse');
assert(gc.dtSec == 300 && gc.lat_res == 2 && gc.lon_res == 2 && gc.mask == 0, ...
       'coarse grid values');
fprintf('  grids: paper %d points, coarse cadence %g s\n', nP, gc.dtSec);

% ---- 2. stitch: satellites 5 and 6 arrive at 0.7 d and 1.3 d ----
H0 = assemble_histories(M, [], 1:6, struct('onSlots', true, 'tEnd', 2*day));
H  = H0;
tA = [0.7 1.3]*day;
for kk = 1:2
    sc = 4 + kk;
    H.tArr(sc) = tA(kk);
    H.phase(H.t < tA(kk), sc) = 2;
end
G = mission_gdop(H, M, [0 2*day], gc, struct('keepSeries', true, 'score', false));
dt = 300/M.tStar;
Nt = round(2*86400/300);
assert(numel(G.t) == Nt, 'epoch count %d vs %d', numel(G.t), Nt);
assert(all(abs(diff(G.t) - dt) < 1e-12), 'duplicated or dropped epoch');
wantOn = 4 + (G.t >= tA(1)) + (G.t >= tA(2));
assert(isequal(G.nOn, wantOn), 'nOn steps');
assert(size(G.runs, 1) == 3 && isequal(G.runs(:, 3).', [4 5 6]), 'three runs');
assert(G.runs(1, 1) == 1 && G.runs(end, 2) == Nt && ...
       all(G.runs(2:end, 1) == G.runs(1:end-1, 2) + 1), 'runs tile the epochs');
assert(G.t(G.runs(2, 1)) >= tA(1) && G.t(G.runs(2, 1) - 1) < tA(1), 'run 2 start');
assert(isequal(size(G.gdop), [Nt numel(G.lat)]), 'gdop size');
% the satellite SET inside a run: run 2 (spacecraft 1-5) against a direct
% engine call on that subset over the same epochs
r2r = G.runs(2, 1):G.runs(2, 2);
[gdD, vcD] = direct_engine(H, M, gc, r2r, 1:5);
relR2 = relgap(G.gdop(r2r, :), gdD);
assert(isequal(isnan(G.gdop(r2r, :)), isnan(gdD)) && relR2 < 1e-9 && ...
       isequal(G.gdop(r2r, :) < 6, gdD < 6) && isequal(double(G.vis(r2r, :)), vcD), ...
       'run 2 vs a direct call on spacecraft 1-5: rel %.3e', relR2);
[gd6, ~] = direct_engine(H, M, gc, r2r, 1:6);
assert(~isequal(gd6 < 6, gdD < 6), 'the subset test must be able to fail');
fprintf(['  stitch: Nt = %d, runs = %s, nOn steps at epochs %d and %d; ' ...
         'run 2 vs direct 5-satellite call rel %.2e\n'], Nt, ...
        mat2str(G.runs), G.runs(2, 1), G.runs(3, 1), relR2);

% ---- 3. equivalence with a direct engine call (all six on station) ----
G3 = mission_gdop(H0, M, [0 2*day], gc, struct('keepSeries', true, 'score', false));
gdC = cell(1, numel(gc.chunks));  vcC = gdC;
for cc = 1:numel(gc.chunks)
    p = struct('lat_range', gc.chunks{cc}, 'lon_range', gc.lon_range, ...
               'lat_res', gc.lat_res, 'lon_res', gc.lon_res, ...
               'time_samples', Nt, 'elevation_mask', gc.mask);
    R = gdop_grid_from_states(H0.t(1:Nt+1), H0.r(1:Nt+1, :, :), M.muStar, ...
                              M.lStar, M.tStar, p);
    gdC{cc} = R.gdop_time_series;
    vcC{cc} = R.visibility_count;
end
gd = [gdC{:}];  vc = [vcC{:}];
relGap = relgap(G3.gdop, gd);
assert(isequal(isnan(G3.gdop), isnan(gd)) && relGap < 1e-9, ...
       'wrapper vs direct engine: rel %.3e', relGap);
assert(isequal(double(G3.vis), vc), 'visibility vs direct engine');
% time-blocked engine calls are exact
G3b = mission_gdop(H0, M, [0 2*day], gc, struct('keepSeries', true, ...
                   'score', false, 'maxBaseBytes', 2e6));
% (node-time rounding of order 1e-16 ND is amplified only where the
% geometry is ill-conditioned: measured 1.7e-10 at GDOP 614, 2.4e-14 below
% GDOP 20; the coverage flags the scores use are identical)
relB  = relgap(G3b.gdop, G3.gdop);
wellC = G3.gdop < 20;
relB20 = relgap(G3b.gdop(wellC), G3.gdop(wellC));
assert(G3b.nCalls > G3.nCalls && relB < 1e-9 && relB20 < 1e-12, ...
       'time blocks: rel %.3e (GDOP < 20: %.3e)', relB, relB20);
assert(isequal(G3b.gdop < 6, G3.gdop < 6) && isequal(G3b.gdop < 7, G3.gdop < 7) ...
       && isequal(G3b.vis, G3.vis), 'time blocks change a coverage flag');
% one- and two-epoch runs
for nE = [1 2]
    i0 = 101;
    Gs = mission_gdop(H0, M, G3.t(i0) + [0 nE*dt], gc, ...
                      struct('keepSeries', true, 'score', false));
    assert(numel(Gs.t) == nE, 'short window length');
    relS = relgap(Gs.gdop, G3.gdop(i0:i0+nE-1, :));
    assert(relS < 1e-9, '%d-epoch run: rel %.3e', nE, relS);
end
fprintf(['  equivalence: rel gap %.2e (one call %d engine calls), ' ...
         'blocked %.2e / %.2e below GDOP 20 (%d calls), short runs ok\n'], ...
        relGap, G3.nCalls, relB, relB20, G3b.nCalls);

% ---- 4. fewer than four on station for the first half day ----
H4 = H0;
H4.tArr(4:6) = 0.5*day;
H4.phase(H4.t < 0.5*day, 4:6) = 2;
G4 = mission_gdop(H4, M, [0 1*day], gc, struct('keepSeries', true, 'score', false));
early = G4.t < 0.5*day;
assert(all(G4.nOn(early) == 3) && all(G4.nOn(~early) == 6), 'nOn with three');
assert(all(isnan(G4.gdop(early, :)), 'all'), 'GDOP must be NaN with three');
assert(any(isfinite(G4.gdop(~early, :)), 'all'), 'GDOP after all six arrive');
nE4 = nnz(early);
vc3C = cell(1, numel(gc.chunks));
for cc = 1:numel(gc.chunks)
    p = struct('lat_range', gc.chunks{cc}, 'lon_range', gc.lon_range, ...
               'lat_res', gc.lat_res, 'lon_res', gc.lon_res, ...
               'time_samples', nE4, 'elevation_mask', gc.mask);
    R = gdop_grid_from_states(H0.t(1:nE4+1), H0.r(1:nE4+1, :, 1:3), M.muStar, ...
                              M.lStar, M.tStar, p);
    vc3C{cc} = R.visibility_count;
end
vc3 = [vc3C{:}];
assert(isequal(double(G4.vis(early, :)), vc3), 'visibility count of the three');
fprintf('  three on station: %d epochs NaN, max visible %d\n', nE4, ...
        max(G4.vis(early, :), [], 'all'));
% 0, 1 and 2 on station through the real engine path (rT collapses to
% [n x 3] for one satellite): NaN GDOP, visibility of exactly that subset,
% and R1's vis >= 2 iff two on-station satellites are above the mask
H6 = H0;
H6.tArr = [0.1 0.25 0.5 0.5 0.5 0.5]*day;
for sc = 1:6, H6.phase(H6.t < H6.tArr(sc), sc) = 2; end
G6 = mission_gdop(H6, M, [0 0.75*day], gc, struct('keepSeries', true, 'score', false));
assert(isequal(G6.runs(:, 3).', [0 1 2 6]), 'runs 0/1/2/6: %s', mat2str(G6.runs));
for rr = 1:3
    rows = G6.runs(rr, 1):G6.runs(rr, 2);
    assert(all(isnan(G6.gdop(rows, :)), 'all'), 'GDOP NaN with %d', rr - 1);
    if rr == 1
        assert(all(G6.vis(rows, :) == 0, 'all'), 'zero visibility with none');
        continue;
    end
    [~, vcS] = direct_engine(H6, M, gc, rows, 1:rr-1);
    assert(isequal(double(G6.vis(rows, :)), vcS), 'visibility with %d', rr - 1);
    assert(isequal(G6.vis(rows, :) >= 2, vcS >= 2), 'R1 flag with %d', rr - 1);
end
rows2 = G6.runs(3, 1):G6.runs(3, 2);
n2 = nnz(G6.vis(rows2, :) >= 2);
assert(n2 > 0 && n2 < numel(G6.vis(rows2, :)), 'R1 flag must vary in the 2-run');
fprintf('  0/1/2 on station: runs %s, R1 flag true on %d of %d samples with two\n', ...
        mat2str(G6.runs(:, 3).'), n2, numel(G6.vis(rows2, :)));
% a scored window shorter than one day bin is refused
try
    mission_gdop(H0, M, [0 0.5*day], gc, struct());
    error('test:noRefusal', 'a scored half-day window was accepted');
catch err
    assert(strcmp(err.identifier, 'mission_gdop:window'), 'window id: %s', err.identifier);
end
% inconsistent phase (phase says transfer, tArr says on station) is refused
H5 = H0;  H5.phase(1:10, 2) = 2;
try
    mission_gdop(H5, M, [0 0.2*day], gc, struct('score', false));
    error('test:noRefusal', 'phase/tArr inconsistency accepted');
catch err
    assert(strcmp(err.identifier, 'mission_gdop:phase'), 'phase id: %s', err.identifier);
end
% deliberate failure of the slot check: one on-station position of H.r
% moved 1e-6 ND (0.39 km), 1000x the 1e-9 ND tolerance, is refused
H7 = H0;  H7.r(25, 1, 3) = H7.r(25, 1, 3) + 1e-6;
try
    mission_gdop(H7, M, [0 0.2*day], gc, struct('score', false, 'quiet', true));
    error('test:noRefusal', 'a perturbed on-station position was accepted');
catch err
    assert(strcmp(err.identifier, 'mission_gdop:slot'), 'slot id: %s', err.identifier);
    fprintf('  perturbed on-station position refused: %s\n', err.message);
end

% ---- 5. validation seam, paper chunk 5 ----
refFile = fullfile(getenv('HOME'), 'Desktop', 'proj7', 'figures', ...
                   'requirements_table_results.mat');
Lr = load(refFile, 'S2');
wpLat = Lr.S2(1).r.r2_worstpt(1);
t5 = tic;
V = mission_gdop('validate', M, struct('chunkSel', 5));
fprintf(['  validate chunk 5 (%.0f s): %s | ours %.6f %.6f %.6f h | ' ...
         'ref %.6f %.6f %.6f h\n'], toc(t5), V.status, V.ours.r1Floor, ...
        V.ours.r2Floor, V.ours.r3Floor, V.ref.r1Floor, V.ref.r2Floor, V.ref.r3Floor);
if isfield(V, 'perPointMaxDiff')
    fprintf('  per point max |diff|: R1 %.6f  R2 %.6f  R3 %.6f h (%d points)\n', ...
            V.perPointMaxDiff.r1, V.perPointMaxDiff.r2, V.perPointMaxDiff.r3, ...
            V.perPointMaxDiff.nPts);
end
fprintf('  reason: %s\n', V.reason);
assert(isequal(V.G.worst.r3g7, [NaN NaN]) || all(isnan(V.G.worst.r3g7)), ...
       'worst.r3g7 must be NaN when the row is off');
assert(strcmp(V.status, 'PASS'), 'chunk-5 validation: %s (%s)', V.status, V.reason);
if wpLat >= -77.5 - 1e-9 && wpLat <= -75 + 1e-9
    assert(abs(V.ours.r2Floor - 0.754032) < 0.005, 'R2 floor %.6f', V.ours.r2Floor);
else
    fprintf('  R2 worst point latitude %.2f is not in chunk 5: assertion skipped\n', wpLat);
end

% ---- 6. UNRESOLVED ----
V6 = mission_gdop('validate', M, struct('refFile', fullfile(tempdir, 'nope.mat')));
assert(strcmp(V6.status, 'UNRESOLVED'), 'missing file: %s', V6.status);
bad = [tempname '.mat'];
notS2 = 1;  save(bad, 'notS2');
cln = onCleanup(@() delete(bad));
V6b = mission_gdop('validate', M, struct('refFile', bad));
assert(strcmp(V6b.status, 'UNRESOLVED'), 'file without S2: %s', V6b.status);
fprintf('  UNRESOLVED: "%s" | "%s"\n', V6.reason, V6b.reason);
% reference identity: a different reference phasing, or a perturbed M.dtau6
V6c = mission_gdop('validate', M, struct('refPhasing', M.dtau6 + [0 0 0 0 0 1e-6]));
assert(strcmp(V6c.status, 'UNRESOLVED') && startsWith(V6c.reason, ...
       'reference predates current phasing'), 'refPhasing: %s', V6c.status);
Mp = M;  Mp.dtau6(6) = Mp.dtau6(6) + 1e-6;
V6d = mission_gdop('validate', Mp, struct());
assert(strcmp(V6d.status, 'UNRESOLVED') && startsWith(V6d.reason, ...
       'reference predates current phasing'), 'M.dtau6: %s', V6d.status);
fprintf('  UNRESOLVED (identity): "%s"\n', V6c.reason);

% ---- 7. deliberate failure: satellite 6 moved off its slot ----
V7 = mission_gdop('validate', M, struct('chunkSel', 5, ...
                  'dtauOverride', M.dtau6 + [0 0 0 0 0 0.3]));
fprintf('  moved slot 6: %s | ours %.6f %.6f %.6f h\n', V7.status, ...
        V7.ours.r1Floor, V7.ours.r2Floor, V7.ours.r3Floor);
assert(strcmp(V7.status, 'FAIL'), 'moved satellite must FAIL, got %s', V7.status);

ok = true;
end

function r = relgap(a, b)
f = isfinite(a) & isfinite(b);
r = max(abs(a(f) - b(f))./max(abs(b(f)), 1e-300), [], 'all');
if isempty(r), r = 0; end
end

function [gd, vc] = direct_engine(H, M, gc, rows, sats)
% gdop_grid_from_states on H's own samples rows(1):rows(end)+1 and the
% spacecraft subset sats, half-open onto exactly those rows' epochs.
n   = numel(rows);
idx = rows(1):rows(end) + 1;
gdC = cell(1, numel(gc.chunks));  vcC = gdC;
for cc = 1:numel(gc.chunks)
    p = struct('lat_range', gc.chunks{cc}, 'lon_range', gc.lon_range, ...
               'lat_res', gc.lat_res, 'lon_res', gc.lon_res, ...
               'time_samples', n, 'elevation_mask', gc.mask);
    R = gdop_grid_from_states(H.t(idx), H.r(idx, :, sats), M.muStar, M.lStar, ...
                              M.tStar, p);
    gdC{cc} = R.gdop_time_series;
    vcC{cc} = R.visibility_count;
end
gd = [gdC{:}];  vc = [vcC{:}];
end
