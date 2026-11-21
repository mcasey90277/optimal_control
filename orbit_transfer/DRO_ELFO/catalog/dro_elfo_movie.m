function [outStem, info] = dro_elfo_movie(catalog, entries, outStem, opts)
%% Purpose:
%
%   A MOVIE (.mp4 and .gif) of one or two DRO -> ELFO catalog transfers in
%   the pumpkynPie house style: the lit dark CR3BP scene (textured Moon,
%   Earth, coloured DRO and ELFO, coloured trail, satellite marker) built by
%   pumpkynPie.plot.SatelliteAnimator -- the scene of
%   costate_catalog_extremes_movies (DRO_tulip/indirect), whose recipe this
%   copies.
%
%   ONE entry: a single full-frame panel. TWO entries (a 2-element cell):
%   side by side, left and right separated by a vertical divider, BOTH ON
%   ONE CLOCK: frames advance on a shared time base, so the shorter
%   transfer arrives early and freezes on its final frame while the longer
%   one plays out; the clock (days) shows identical values in both panels.
%
%  ASSUMPTIONS / NOTES:
%
% - Trajectories and orbits come from dro_elfo_fly (pumpkyn's tfMinProp
%   from z8; orbits by the sheet's own recipe), so what you watch is the
%   catalog's PMP solution.
% - Aerospace-Toolbox guard (the ShowEarth trap): planetEphemeris without
%   its data pack is hidden from the path for this render only, and only
%   when it lives under matlabroot.
% - Every output frame is EXACTLY 1280 x 720 (one 1280-px panel, or two
%   637-px panels and a 6-px divider), enforced by index resampling -- the
%   cure for the H.264 diagonal-shear artifact. The .gif is written at the
%   same size, every third frame at 10 fps.
% - Needs pumpkynPie.plot.SatelliteAnimator; if it is not on the path the
%   function warns dro_elfo_movie:noAnimator and returns ''.
%
%% Inputs:
%
%  catalog                  struct | char           The catalog struct, or
%                                                   the path of its .mat
%
%  entries                  struct | cell           One entry [1 x 1]
%                                                   struct, or a {1 x 2}
%                                                   cell of entries (pick,
%                                                   transfers_at or
%                                                   extremes rows)
%
%  outStem                  char (optional)         Output stem for
%                                                   '<stem>.mp4' and
%                                                   '<stem>.gif'
%                                                   ['dro_elfo_movie'];
%                                                   its folder is created
%
%  opts                     struct (optional)       .labels {1 x nP} panel
%                                                   labels [{'A','B'}, or
%                                                   {''} for one panel]
%                                                   .nFrames [240] frames
%                                                   .visible [false] scene
%                                                   figures visible
%
%% Outputs:
%
%  outStem                  char                    [1 x n] the stem
%                                                   written, or '' when the
%                                                   animator is unavailable
%
%  info                     struct                  [1 x 1] .F {1 x nP}
%                                                   dro_elfo_fly structs
%                                                   .nFrames .mp4 .gif
%                                                   (paths) .tEndDays
%
%% Revision History:
%  M. Casey                                                   (c) 11/21/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    catalog = fullfile(here, 'results', 'costate_catalog_dro_elfo.mat');
    if ~isfile(catalog), catalog = fullfile(here, 'results', 'pilot_chain', 'costate_catalog_dro_elfo.mat'); end
    if ~isfile(catalog), fprintf('dro_elfo_movie: self-demo needs a catalog under results/\n');  outStem = '';  info = [];  return, end
    [eMin, eMax] = dro_elfo_extremes(catalog, struct('thrustN', 1), 'time');
    [outStem, info] = dro_elfo_movie(catalog, {eMin, eMax}, fullfile(here, 'results', 'library_demo', 'demo_movie'), ...
                                     struct('labels', {{'FASTEST', 'SLOWEST'}}, 'nFrames', 120));
    return
end
if nargin < 3 || isempty(outStem), outStem = 'dro_elfo_movie'; end
if nargin < 4 || isempty(opts), opts = struct(); end
outStem = char(outStem);
info = struct('F', {{}}, 'nFrames', 0, 'mp4', '', 'gif', '', 'tEndDays', NaN);
if ischar(catalog) || isstring(catalog)
    L = load(char(catalog));  fn = fieldnames(L);
    assert(isscalar(fn), 'dro_elfo_movie:catalog', 'expected one variable in %s', char(catalog));
    catalog = L.(fn{1});
end
if isstruct(entries), R = {entries}; elseif iscell(entries), R = entries(:).'; else, R = {}; end
if ~(isscalar(R) || numel(R) == 2) || ~all(cellfun(@isstruct, R))
    error('dro_elfo_movie:entries', 'entries must be one entry struct or a 2-element cell of entries');
end
nP = numel(R);
if isfield(opts, 'labels') && ~isempty(opts.labels)
    lab = cellstr(opts.labels);
    if numel(lab) ~= nP
        error('dro_elfo_movie:labels', 'opts.labels must hold %d label(s), got %d', nP, numel(lab));
    end
elseif nP == 2, lab = {'A', 'B'};
else,           lab = {''};
end
nFrames = 240;  if isfield(opts, 'nFrames') && ~isempty(opts.nFrames), nFrames = opts.nFrames; end
visible = isfield(opts, 'visible') && ~isempty(opts.visible) && logical(opts.visible);

if isempty(which('pumpkynPie.plot.SatelliteAnimator'))
    warning('dro_elfo_movie:noAnimator', 'pumpkynPie.plot.SatelliteAnimator not on the path -- no movie.');
    outStem = '';  return
end

%% Ephemeris guard: planetEphemeris on the path without its data package
%  is hidden for this render (restored on exit) so pumpkyn's analytic
%  Earth-Moon model takes over; only a MathWorks folder is ever hidden:
if exist('planetEphemeris', 'file')
    try
        planetEphemeris(juliandate(datetime(2030, 1, 1)), 'Earth', 'Moon');
    catch
        aeroDir = fileparts(which('planetEphemeris'));
        if startsWith(aeroDir, matlabroot)
            rmpath(aeroDir);
            cleanupPath = onCleanup(@() addpath(aeroDir));
            warning('dro_elfo_movie:noEphemData', ...
                'planetEphemeris found but its data package is missing; using pumpkyn''s analytic Earth-Moon model for this render.');
        end
    end
end

%% Flights:
tStar  = catalog.constants.tStar_s;
params = struct('muStar', catalog.constants.muStar, 'lStar', catalog.constants.lStar_km, 'tStar', tStar);
epoch  = juliandate(datetime(2030, 1, 1, 0, 0, 0));
Fl = cell(1, nP);
for k = 1:nP, Fl{k} = dro_elfo_fly(catalog, R{k}); end

%% Scenes and frames (figures closed on success and on error):
od = fileparts(outStem);
if ~isempty(od) && ~isfolder(od), mkdir(od); end
S = struct('fig', cell(1, nP), 'anim', [], 'clock', [], 'tf', []);
figs0 = findall(groot, 'Type', 'figure');
try
    S = build_scenes(S, Fl, lab, visible, params, epoch);
    tEnd = write_frames(S, outStem, nFrames, epoch, tStar);
catch err
    close(setdiff(findall(groot, 'Type', 'figure'), figs0));   % only the scenes opened here
    rethrow(err);
end
close_all(S);
info.F = Fl;  info.nFrames = nFrames;  info.mp4 = [outStem '.mp4'];  info.gif = [outStem '.gif'];
info.tEndDays = tEnd*tStar/86400;
fprintf('  movie -> %s.mp4 / .gif  (%d frame(s), 1280x720, %d panel(s), clock to %.2f d)\n', ...
        outStem, nFrames, nP, info.tEndDays);
end

% ---------------------------------------------------------------------------
function S = build_scenes(S, Fl, lab, visible, params, epoch)
% BUILD_SCENES  One pumpkyn-style lit scene per panel (the recipe of
% costate_catalog_extremes_movies).  INPUTS: S [1 x nP] struct (.fig .anim
% .clock .tf, filled here); Fl {1 x nP} dro_elfo_fly structs; lab {1 x nP}
% labels; visible [logical]; params [struct] .muStar .lStar .tStar; epoch
% [1 x 1] Julian date.  OUTPUTS: S [1 x nP] filled.
nP = numel(Fl);
tStar = params.tStar;
if nP == 1, pw = 1280; else, pw = 637; end          % 2*637 + 6 = 1280
if visible, vis = 'on'; else, vis = 'off'; end
for k = 1:nP
    F  = Fl{k};
    O  = F.orbits;
    [tu, iu] = unique(F.t);
    S(k).fig = figure('Color', 'k', 'Position', [60 60 pw 720], 'Visible', vis);
    S(k).anim = pumpkynPie.plot.SatelliteAnimator(S(k).fig, ...
        'Epoch',                 epoch, ...
        'CR3BPParameters',       params, ...
        'ShowEarth',             true, ...
        'FigureTheme',           'dark', ...
        'BodyLighting',          'lambertian', ...
        'FitMode',               'trajectory', ...
        'TrailMode',             'distance', ...
        'EnableTrailColor',      true, ...
        'EnableSparkles',        false, ...
        'EnableSatelliteMarker', true, ...
        'CameraUpVector',        [0 0 1], ...
        'View',                  [-35 22]);
    S(k).anim.addOrbit(epoch + O.tD*tStar/86400, O.rvD, ...
        'Name', 'DRO', 'Color', [0.25 0.60 0.35], 'LineWidth', 0.9);
    S(k).anim.addOrbit(epoch + O.tA*tStar/86400, O.rvA, ...
        'Name', 'ELFO', 'Color', [0.75 0.28 0.28], 'LineWidth', 0.9);
    S(k).anim.addSatellite(epoch + tu*tStar/86400, F.Y(iu, 1:6), ...
        'Name', strtrim(sprintf('%s transfer', lab{k})), ...
        'MarkerColor', [1.00 0.92 0.55], 'TrailColor', [0.35 0.72 1.00]);
    S(k).anim.initialize();
    camzoom(1.4);
    head = sprintf('DRO %.4g ND -> ELFO raan %g deg, %.3g N, s_D %.3f s_A %.3f', ...
                   F.sheet_key.depSize, F.sheet_key.raanDeg, F.thrustN, F.sD, F.sA);
    if ~isempty(lab{k}), head = sprintf('%s:  %s', lab{k}, head); end
    annotation(S(k).fig, 'textbox', [0.02 0.90 0.96 0.09], 'String', ...
        {head, sprintf('\\DeltaV %.3f km/s,  t_f %.3f d,  min Moon alt %.0f km', F.dvKms, F.tfDays, F.minMoonAltKm)}, ...
        'Color', 'w', 'EdgeColor', 'none', 'FontSize', 11, 'HorizontalAlignment', 'center');
    S(k).clock = annotation(S(k).fig, 'textbox', [0.02 0.02 0.96 0.05], ...
        'String', '', 'Color', [1 0.92 0.55], 'EdgeColor', 'none', ...
        'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    S(k).tf = F.tfNd;
end
end

% ---------------------------------------------------------------------------
function tEnd = write_frames(S, outStem, nFrames, epoch, tStar)
% WRITE_FRAMES  The frames on one shared clock (a shorter panel freezes at
% its final pose), each exactly 1280 x 720, to <stem>.mp4 at 30 fps and
% every third frame to <stem>.gif at 10 fps.  INPUTS: S [1 x nP] scenes;
% outStem [char]; nFrames [1 x 1]; epoch [1 x 1] Julian date; tStar
% [1 x 1] s.  OUTPUTS: tEnd [1 x 1] the clock's last tick, ND.
nP = numel(S);
if nP == 1, pw = 1280; else, pw = 637; end
tEnd = max([S.tf]);
divider = uint8(90*ones(720, 1280 - 2*pw, 3));
vw = VideoWriter([outStem '.mp4'], 'MPEG-4');
vw.FrameRate = 30;  vw.Quality = 95;  open(vw);
cleanupVW = onCleanup(@() close(vw));
gifFile = [outStem '.gif'];
for kf = 1:nFrames
    tNow = tEnd*(kf - 1)/max(nFrames - 1, 1);
    for k = 1:nP
        S(k).anim.update(epoch + min(tNow, S(k).tf)*tStar/86400);
        set(S(k).clock, 'String', sprintf('t = %.2f days', tNow*tStar/86400));
    end
    drawnow limitrate
    if nP == 1
        img = panel_frame(S(1).fig, pw);
    else
        img = [panel_frame(S(1).fig, pw), divider, panel_frame(S(2).fig, pw)];
    end
    writeVideo(vw, img);
    if mod(kf, 3) == 1
        [A, cmap] = rgb2ind(img, 256);
        if kf == 1
            imwrite(A, cmap, gifFile, 'gif', 'LoopCount', inf, 'DelayTime', 0.1);
        else
            imwrite(A, cmap, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end
end
clear cleanupVW
end

% ---------------------------------------------------------------------------
function img = panel_frame(fig, w)
% PANEL_FRAME  getframe forced to 720 x w by index resampling (the H.264
% shear cure).  INPUTS: fig [figure]; w [1 x 1] target width, px.
% OUTPUTS: img [720 x w x 3 uint8].
F = getframe(fig);
ri = max(1, min(size(F.cdata, 1), round(linspace(1, size(F.cdata, 1), 720))));
ci = max(1, min(size(F.cdata, 2), round(linspace(1, size(F.cdata, 2), w))));
img = F.cdata(ri, ci, :);
end

% ---------------------------------------------------------------------------
function close_all(S)
% CLOSE_ALL  Close the scene figures that are still open.  INPUTS: S
% [1 x nP] struct with .fig.  OUTPUTS: none.
for k = 1:numel(S)
    if ~isempty(S(k).fig) && isgraphics(S(k).fig), close(S(k).fig); end
end
end
