function render_phase1_movies()
% RENDER_PHASE1_MOVIES  The three Phase-1 movies: min-time / min-energy / min-fuel.
%
% Renders, in the 2-body campaign's house movie style (transfer_movie.m):
%   1. movie_mintime_2body_T10N   - the TWO-BODY certified min-time solution
%      (MEE_mintime_T10_npr15.mat; all-burn, c_tf=1). No CR3BP min-time
%      exists (not part of Phase 1) -- labeled honestly as 2-body.
%   2. movie_energy_cr3bp_T10N    - CR3BP min-ENERGY at gain=1 (bridge output).
%   3. movie_minfuel_cr3bp_T10N   - CR3BP certified min-FUEL (the +0.0545 kg
%      headline solution).
% Each is converted MEE->Cartesian via mee_res_to_cart_res (nDense densify --
% the polygonal-movie lesson) and written as MP4+GIF under results/movies/.
%
% INPUTS:  none (paths hardcoded to the Phase-1 artifacts)
% OUTPUTS: none (files written; paths printed)
%
% REFERENCES:
%   [1] ../../earth_elliptic_to_geo/direct/viz/transfer_movie.m (renderer).
%   [2] doc/cr3bp_geo_phase1_note.tex (the solutions being rendered).
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'setup_paths.m'));
addpath(fullfile(here, '..', '..', '..', 'earth_elliptic_to_geo', 'direct', 'viz'));
resDir = fullfile(here, '..', 'results');
movDir = fullfile(resDir, 'movies');
if ~exist(movDir, 'dir'), mkdir(movDir); end
e2bRes = fullfile(here, '..', '..', '..', 'earth_elliptic_to_geo', 'direct', 'results');
nDense = 8;                                % densify per segment (smooth curves)

% --- 1. min-time (2-body certified, all-burn) --------------------------------
S = load(fullfile(e2bRes, 'MEE_mintime_T10_npr15.mat'));
o = S.out.solverOut;                       % driver-level out wraps the solver struct
sg = linspace(0, 1, size(o.X,2)).';
res = mee_res_to_cart_res(o.X, o.U, o.dL, sg, 10, 1.0, 1, nDense);
res.cfg.label = 'min-time (2-body)';
transfer_movie(res, fullfile(movDir, 'movie_mintime_2body_T10N'));

% --- 2. min-energy (CR3BP, Moon on, gain=1) ----------------------------------
E = load(fullfile(resDir, 'energy_cr3bp_T10N_phi0.mat'));
res = mee_res_to_cart_res(E.X, E.U, E.dL, E.sigma, 10, 1.5, 1, nDense);
res.cfg.label = 'min-energy (CR3BP)';
transfer_movie(res, fullfile(movDir, 'movie_energy_cr3bp_T10N'));

% --- 3. min-fuel (CR3BP, certified) ------------------------------------------
F = load(fullfile(resDir, 'minfuel_cr3bp_T10N_phi0.mat'));
b = F.best;
sg = linspace(0, 1, size(b.X,2)).';
res = mee_res_to_cart_res(b.X, b.U, b.dL, sg, 10, 1.5, 1, nDense);
res.cfg.label = 'min-fuel (CR3BP)';
transfer_movie(res, fullfile(movDir, 'movie_minfuel_cr3bp_T10N'));

fprintf('RENDER_PHASE1_MOVIES done: 3 movies in %s\n', movDir);
end
