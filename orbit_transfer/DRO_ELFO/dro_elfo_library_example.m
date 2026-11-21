%% dro_elfo_library_example
%
%   Worked example of the DRO -> ELFO costate catalog's library functions:
%   what the catalog holds, how to pick a transfer, how to hand it to
%   pumpkyn's tfMin, what one phase cell holds at every thrust, how to fly
%   and plot an entry, the extremes at one thrust, the phase torus, and a
%   side-by-side movie.
%
%     1. describe the catalog (orbits, sheets, rungs, counts, verdicts)
%     2. pick one entry by (DRO size, ELFO plane orientation, phases, thrust)
%     3. hand the pick to pumpkyn.cr3bp.tfMin at its own endpoints
%     4. every thrust rung recorded at that phase cell
%     5. fly the pick and plot it in 3D (rotating frame, Moon-centred)
%     6. the fastest and slowest entries at one thrust
%     7. the phase torus at that thrust (one sheet, then every sheet)
%     8. a movie: fastest beside slowest, on one clock
%
%   Catalog: results/costate_catalog_dro_elfo.mat (12 sheets, raanDeg 0 to
%   330), or the one-sheet pilot results/pilot_chain/... when the full
%   catalog is not there yet. Outputs (pictures, movie) go to
%   results/library_demo/. Nothing is solved or re-certified except by
%   tfMin in section 3. Set showFigures = false before running to keep the
%   figures off screen. Run time about a minute (the movie takes most).
%
%  M. Casey                                                     11/21/2026

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
run(fullfile(here, 'setup_dro_elfo_paths.m'));
if ~exist('showFigures', 'var'), showFigures = usejava('desktop'); end
outDir = fullfile(here, 'results', 'library_demo');
if ~isfolder(outDir), mkdir(outDir); end
figsBefore = findall(groot, 'Type', 'figure');

%% 1. Load the catalog and describe it:
catFile = fullfile(here, 'results', 'costate_catalog_dro_elfo.mat');
if ~isfile(catFile)
    fprintf('NOTE: %s is not there yet; using the one-sheet pilot catalog instead\n', catFile);
    catFile = fullfile(here, 'results', 'pilot_chain', 'costate_catalog_dro_elfo.mat');
end
L    = load(catFile);
fn   = fieldnames(L);
cat_ = L.(fn{1});
D    = dro_elfo_describe(cat_);
fprintf('[1] verdict: %d sheet(s), %d entries, %d rungs from %g to %g N\n', ...
        D.nSheets, D.nEntries, numel(D.rungsN), D.rungsN(1), D.rungsN(end));

%% 2. Pick one entry:
droTau  = D.droTau(1);      % DRO size (period, ND)
raanDeg = D.raanDeg(1);     % ELFO plane orientation (deg), one sheet per value
sD      = 0.30;             % DRO departure phase (fraction of the DRO period)
sA      = 0.60;             % ELFO arrival phase (fraction of the ELFO period, from perilune)
thrustN = 1;                % a catalog rung
e = dro_elfo_pick(cat_, droTau, raanDeg, sD, sA, thrustN);
fprintf('[2] verdict: asked (%.3f, %.3f), got cell (%d,%d) at (%.4f, %.4f), t_f %.4f ND, conj_pass %d\n', ...
        sD, sA, e.iD, e.iA, e.sD, e.sA, e.tf_nd, e.conj_pass);

%% 3. Hand the pick to pumpkyn's tfMin at its own endpoints:
%  The endpoints follow the catalog's orbit_reconstruction rule. tfMin
%  starts from z8 and returns the root it converges to; |dz| near zero means
%  the stored costates already solve the shooting problem.
mu    = cat_.constants.muStar;
tStar = cat_.constants.tStar_s;
lStar = cat_.constants.lStar_km;
cnd   = cat_.thruster.c_nd;
Tnd   = (thrustN/cat_.thruster.m0_kg)*tStar^2/(lStar*1000);
pD = e.dep_params;  pD.muStar = mu;
pA = e.arr_params;  pA.muStar = mu;
[tD, rvD] = get_family_orbit(e.dep_family, pD);
[tA, rvA] = get_family_orbit(e.arr_family, pA);
xD  = phase_state(tD, rvD);  xA = phase_state(tA, rvA);
rv0 = reshape(xD(e.sD), 1, 6);
rvf = reshape(xA(e.sA), 1, 6);
if license('test', 'Optimization_Toolbox') && ~isempty(which('fsolve'))
    tic;  zTf = pumpkyn.cr3bp.tfMin(rv0, rvf, e.z8, Tnd, cnd, mu);  wall = toc;
    dz = sqrt(sum((zTf(:) - e.z8(:)).^2));
    fprintf('[3] verdict: tfMin returned in %.1f s with |dz| = %.2e (t_f %.5f ND -> %.5f ND)\n', ...
            wall, dz, e.z8(8), zTf(8));
else
    fprintf('[3] SKIPPED: tfMin needs fsolve (Optimization Toolbox), not available here\n');
end

%% 4. Every thrust rung recorded at that phase cell:
T = dro_elfo_transfers_at(cat_, droTau, raanDeg, e.sD, e.sA);
fprintf('  %8s %9s %10s %9s %5s\n', 'thrust N', 't_f [d]', 'dV [km/s]', 'prop [kg]', 'conj');
for k = 1:numel(T)
    fprintf('  %8g %9.3f %10.4f %9.3f %5d\n', T(k).thrustN, T(k).tfDays, T(k).dvKms, T(k).propellantKg, T(k).conj_pass);
end
fprintf('[4] verdict: %d of %d rungs recorded at cell (%d,%d)\n', numel(T), numel(cat_.rungs_N), e.iD, e.iA);

%% 5. Fly the pick and plot it:
%  F holds the history in ND and in days / km / kg (F.tDays, F.rMoonKm,
%  F.vKms, F.massKg, F.u thrust direction) for use elsewhere.
[F, hF] = dro_elfo_fly(cat_, e, struct('plot', true, 'visible', showFigures, ...
                                       'out', fullfile(outDir, 'example_fly.png')));
fprintf('[5] verdict: %d samples over %.3f d, arrival miss %.2e km, min Moon altitude %.0f km at day %.3f, %.2f revs\n', ...
        numel(F.t), F.tfDays, F.flyKm, F.minMoonAltKm, F.tMinMoonDays, F.moonRevs);

%% 6. The fastest and slowest entries at one thrust (every sheet):
[eMin, eMax] = dro_elfo_extremes(cat_, struct('thrustN', thrustN), 'time');
fprintf('[6] verdict: at %g N the flight time spans %.3f to %.3f d\n', thrustN, eMin.tfDays, eMax.tfDays);

%% 7. The phase torus at that thrust:
H = dro_elfo_torus_plot(cat_, raanDeg, thrustN, struct('visible', showFigures, ...
                        'out', fullfile(outDir, 'example_torus.png')));
fprintf('[7] verdict: %d of %d cells hold an entry at %g N on raanDeg %g\n', ...
        nnz(isfinite(H.tfDays)), numel(H.tfDays), thrustN, raanDeg);
if D.nSheets > 1
    Ha = dro_elfo_torus_plot(cat_, 'all', thrustN, struct('visible', showFigures, 'droTau', droTau, ...
                             'out', fullfile(outDir, 'example_torus_all.png')));
    fprintf('    every sheet: %d tiles on one scale, %.3f to %.3f d\n', size(Ha.tfDays, 3), Ha.clim(1), Ha.clim(2));
end

%% 8. A movie: the fastest beside the slowest at that thrust, one clock:
%  (120 frames here; omit nFrames for the 240-frame default)
[stem, info] = dro_elfo_movie(cat_, {eMin, eMax}, fullfile(outDir, 'example_extremes_movie'), ...
                              struct('labels', {{'FASTEST', 'SLOWEST'}}, 'nFrames', 120));
if isempty(stem)
    fprintf('[8] SKIPPED: pumpkynPie.plot.SatelliteAnimator is not on the path\n');
else
    fprintf('[8] verdict: %s.mp4 / .gif, %d frames, clock to %.2f d\n', stem, info.nFrames, info.tEndDays);
end

fprintf('ALL SECTIONS DONE; outputs in %s\n', outDir);
if ~showFigures, close(setdiff(findall(groot, 'Type', 'figure'), figsBefore)); end
