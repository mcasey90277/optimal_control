function outFile = record_tulip_animation(solFile, outFile, figSize, dtSec, viewAzEl)
% RECORD_TULIP_ANIMATION  GTO->tulip transfer in the pumpkynPie scene style.
%
% Renders our min-fuel transfer with the proj7 constellation-demo look --
% starfield background, lit Earth and Moon, satellite model, sparkle trail --
% instead of the campaign's 2-D line plot. The trajectory is drawn as
% alternating BURN (red) and COAST (blue) arcs, so the bang-bang control
% structure is visible as colour on the path itself.
%
% Mirrors the scene configuration of
% proj7/visualization/record_constellation_animation.m (which in turn mirrors
% pumpkynPie demos/constellation_Animation.m). Two things necessarily differ,
% and both are about our trajectory rather than the scene:
%
%   1. TIME RESAMPLING IS MANDATORY. The demo propagates a closed orbit on a
%      uniform time grid. Our solution lives on a Sundman mesh, which is
%      non-uniform in time by construction -- measured ratio between the largest
%      and smallest step is ~1.2e8. Played back node-by-node the spacecraft
%      would appear frozen for most of the movie and then jump. So the state is
%      interpolated onto a uniform grid in PHYSICAL time before rendering.
%
%   2. THE PATH IS AN OPEN SPIRAL, not a periodic orbit. The demo derives its
%      trail length from one orbital period; here it is taken as a fraction of
%      the transfer's own arc length.
%
% The throttle is resampled with NEAREST interpolation, never linear: linear
% would blur a bang-bang switch across an interval and paint a purple arc that
% does not exist in the solution.
%
% INPUTS:
%   solFile - solution .mat: needs out.X ([r;v;m;t] in rows 1-8, a 9th cScale
%             row is ignored) and out.U (row 4 = throttle). Accepts the tulip
%             seed layout (S.out) or a top-level S.X/S.U [char|struct]
%   outFile - output MP4 path [char]. Default: alongside the solution
%   figSize - figure position [x y w h]. Default [100 100 1440 1080]
%   dtSec   - animation time step in SECONDS of physical transfer time
%             [scalar]. Default sized to ~1200 frames over the transfer.
%
% OUTPUTS:
%   outFile - path to the written MP4 [char]
%
% SIDE EFFECTS: writes the MP4 plus three preview PNGs (_p15/_p50/_p85).
%
% REFERENCES:
%   [1] proj7/visualization/record_constellation_animation.m (the scene this
%       mirrors; re-diff against it if the demo changes).
%   [2] ../../doc/gto_tulip_guide.pdf sec 3.4 (why the mesh is non-uniform).

if nargin < 1 || isempty(solFile)
    solFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                       'lib', 'sundman_minfuel_certified.mat');
end
if isstruct(solFile), S = solFile; solPath = 'struct';
else,                 S = load(solFile); solPath = solFile; end
if isfield(S,'out'), X = S.out.X; U = S.out.U; else, X = S.X; U = S.U; end

if nargin < 2 || isempty(outFile)
    [d,n] = fileparts(solPath);
    if isempty(d), d = pwd; n = 'tulip_transfer'; end
    outFile = fullfile(d, [n '_scene.mp4']);
end
if nargin < 3 || isempty(figSize), figSize = [100 100 1440 1080]; end
% DEFAULT MATCHES THE DEMO: [90 0], as in
% proj7/visualization/record_constellation_animation.m. Note this is edge-on to
% the Earth-bound spiral, which dominates our transfer, so the spiral reads as a
% near-horizontal band rather than as loops -- that is inherent to the view, not
% a bug. It does put the camera in the Moon's plane, which is what makes the
% lunar arrival legible. Pass e.g. [-40 26] for an oblique view that opens the
% spiral out instead.
if nargin < 5 || isempty(viewAzEl), viewAzEl = [90 0]; end

FPS         = 20;
PUMPKYN_PIE = fullfile(getenv('HOME'),'Desktop','proj7','external','pumpkynPie');
assert(isfolder(PUMPKYN_PIE), 'record_tulip_animation:noPumpkyn', ...
    'pumpkynPie not found at %s', PUMPKYN_PIE);

here = pwd;
cleaner = onCleanup(@() cd(here));      %#ok<NASGU> startup() cd's; come back
cd(PUMPKYN_PIE);
startup();

params = pumpkynPie.cr3bp.getCR3BParams();

%% ---- our trajectory, resampled uniformly in physical time ---------------
tND = X(8,:);                            % carried time state (ND)
tS  = tND * params.tStar;                % seconds
if nargin < 4 || isempty(dtSec)
    dtSec = (tS(end) - tS(1)) / 1200;    % ~1200 frames
end
tU  = (tS(1):dtSec:tS(end)).';
if tU(end) < tS(end), tU(end+1) = tS(end); end

% state: linear in time is right here -- the uniform grid is FINER than the
% Sundman mesh almost everywhere, so this is interpolation, not smoothing
xU = interp1(tS(:), X(1:6,:).', tU, 'linear');

% throttle: NEAREST, to keep the bang-bang edges crisp (see header)
sU = interp1(tS(:), double(U(4,:) > 0.5).', tU, 'nearest') > 0.5;

jd0 = pumpkyn.util.juliandate('01/01/2000 00:00:00');
jd  = jd0 + (tU - tU(1)) / 86400;

% Resampling can DESTROY switches: a bang-bang arc shorter than dtSec simply
% vanishes from the nearest-interpolated throttle, and the movie would then show
% a control structure the solution does not have. Compare against the source and
% say so rather than rendering a quiet lie.
swSrc = sum(abs(diff(double(U(4,:) > 0.5))));
swRes = sum(abs(diff(double(sU(:).'))));
nBurn = sum(diff([false; sU(:)]) == 1);
fprintf('record_tulip_animation: %d nodes -> %d uniform frames (dt %.0f s), %d burn arcs\n', ...
        size(X,2), numel(tU), dtSec, nBurn);
if swRes < swSrc
    warning('record_tulip_animation:switchesLost', ...
        ['resampling at dt = %.0f s kept %d of the solution''s %d switches. Arcs ' ...
         'shorter than dt vanish, so the rendered control structure is COARSER ' ...
         'than the solution. Reduce dtSec (try %.0f s) for a faithful render.'], ...
        dtSec, swRes, swSrc, dtSec * swRes / max(swSrc,1));
end

%% ---- scene (mirrors the proj7 constellation demo) -----------------------
sa = pumpkynPie.plot.SatelliteAnimator( ...
    'Epoch', jd0, ...
    'CR3BPParameters', params, ...
    'FigurePosition', figSize, ...
    'Background', 'stars', ...
    'BackgroundColor', 'k', ...
    'Quality', 'interactive', ...
    'View', viewAzEl, ...
    'CameraUpVector', [0 0 1], ...
    'ViewingAngle', [], ...
    'CameraOrbit', [0 0], ...
    'DataAspectRatio', [1 1 1], ...
    'FitMode', 'trajectory', ...
    'FitPadding', 0.05, ...
    'ShowLogo', true, ...
    'TrailMode', 'distance', ...
    'TrailLength', 0.12 * pumpkyn.util.arclength(xU(:,1:3), 2), ...
    'TrailPointCount', 96, ...
    'TrailLineWidth', 2, ...
    'SparkleTrailWidth', 0.4, ...
    'EnableTrailColor', false, ...
    'EnableSparkles', true);

% --- the transfer path, one addOrbit call per constant-throttle arc -------
% Arcs are drawn separately because addOrbit takes ONE colour; splitting on
% throttle transitions is what puts the control structure into the picture.
% Consecutive arcs share their boundary sample so the path has no visible gaps.
BURN  = [0.90 0.20 0.15];
COAST = [0.25 0.55 0.95];
edges = [1; find(diff(sU(:)) ~= 0) + 1; numel(sU)];
for k = 1:numel(edges)-1
    i1 = edges(k);  i2 = min(edges(k+1) + 1, numel(sU));
    if i2 <= i1, continue; end
    if sU(i1), c = BURN; nm = 'burn'; else, c = COAST; nm = 'coast'; end
    sa.addOrbit(jd(i1:i2), xU(i1:i2,:), 'Name', nm, 'Color', c, 'LineWidth', 1.6);
end

% --- the spacecraft ------------------------------------------------------
% Nt-by-6 is the single-satellite form the animator accepts (the demo's
% permute to Nt-by-6-by-N is only needed for a constellation).
sa.addSatellite(jd, xU, 'Name', 'GTO->tulip', 'MarkerColor', [1 1 1], ...
                'TrailColor', [1 1 1], 'SparkleColor', [1.00 0.82 0.45]);
sa.initialize();

fig = sa.Figure;  ax = sa.Axes;
try, ax.Toolbar.Visible = 'off';        catch, end %#ok<CTCH>
try, axtoolbar(ax, 'none');             catch, end %#ok<CTCH>
try, disableDefaultInteractivity(ax);   catch, end %#ok<CTCH>
drawnow;

%% ---- record ------------------------------------------------------------
vw = VideoWriter(outFile, 'MPEG-4');
vw.FrameRate = FPS;  vw.Quality = 92;
open(vw);
closer = onCleanup(@() close(vw));      %#ok<NASGU> no truncated file on error

t0 = tic;
for k = 1:numel(jd)
    sa.update(jd(k));
    drawnow;
    writeVideo(vw, getframe(fig));
    if mod(k, 250) == 0
        fprintf('  %4d/%4d (%.0f s)\n', k, numel(jd), toc(t0));
    end
end
clear closer;                            % flush before reading back

%% ---- preview stills ----------------------------------------------------
v = VideoReader(outFile);
[outDir, outName] = fileparts(outFile);
for frac = [0.15 0.5 0.85]
    imwrite(read(v, max(1, round(v.NumFrames*frac))), ...
            fullfile(outDir, sprintf('%s_p%02d.png', outName, round(frac*100))));
end
fprintf('Done: %s (%.0f s render, %d frames, %.1f s runtime)\n', ...
        outFile, toc(t0), numel(jd), v.Duration);
end
