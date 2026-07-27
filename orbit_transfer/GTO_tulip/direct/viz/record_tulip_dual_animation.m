function outFile = record_tulip_dual_animation(solFile, outFile, panelSize, dtSec, opts)
% RECORD_TULIP_DUAL_ANIMATION  Two-panel GTO->tulip movie with a fading tail.
%
% Left panel  : looking DOWN on the Earth-Moon system (view [0 90]).
% Right panel : the constellation-demo view [90 0], in the Moon's plane, with
%               the full target tulip drawn in white so the spacecraft can be
%               seen entering it.
%
% Both panels show a FADING tail rather than the whole path: colour encodes the
% control (red = burn, blue = coast) and opacity encodes recency. This is the
% change that makes the edge-on right panel legible at all -- with the complete
% path drawn, ~40 stacked revolutions collapse into a solid band.
%
% WHY THE TAIL IS DRAWN BY HAND. The animator's built-in trail takes ONE colour
% per satellite (TrailColor is a single RGB), so it cannot express burn/coast.
% The tail here is a patch whose per-vertex CData carries the throttle state and
% whose per-vertex AlphaData carries recency, which gives both properties in one
% object and keeps the fade under our control.
%
% WHY TWO ANIMATORS. SatelliteAnimator owns its figure, so a side-by-side layout
% is two instances rendered independently and their frames concatenated. Both
% panels must therefore share a pixel height.
%
% Differences from record_tulip_animation (the single-panel version) are all
% requested styling: sparkles off, larger satellite marker, fading tail instead
% of a persistent path, and the white tulip reference on the right.
%
% INPUTS:
%   solFile   - solution .mat or struct (out.X rows 1-8, out.U row 4 = throttle)
%   outFile   - output MP4 path [char]
%   panelSize - per-panel figure size [x y w h]; the movie is 2*w wide
%               [default [100 100 640 720] -> a 1280x720 movie]
%   dtSec     - uniform physical-time step, seconds [default ~1200 frames]
%   opts      - (optional) struct:
%       .tailSeconds  fade length in transfer seconds        [default 4 days]
%       .markerSize   animator MarkerSize triple             [default [26 18 9]]
%       .tailWidth    tail line width                        [default 2.5]
%       .tulipColor   colour of the reference tulip           [default 0.38 grey]
%       .tulipWidth   its line width                          [default 0.7]
%       .zoomPad      right-panel padding, fraction of the tulip box [default 0.25]
%
% OUTPUTS: outFile - path to the written MP4 [char]
%
% REFERENCES:
%   [1] record_tulip_animation.m (single-panel; shares the resampling logic).
%   [2] proj7/visualization/record_constellation_animation.m (the scene).

if nargin < 1 || isempty(solFile)
    solFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                       'lib', 'sundman_minfuel_certified.mat');
end
if isstruct(solFile), S = solFile; solPath = 'struct';
else,                 S = load(solFile); solPath = solFile; end
if isfield(S,'out'), X = S.out.X; U = S.out.U; else, X = S.X; U = S.U; end
% 640x720 per panel -> a 1280x720 movie. BOTH dimensions divisible by 16: an
% H.264 frame size that is not produces the diagonal colour shear seen before.
if nargin < 3 || isempty(panelSize), panelSize = [100 100 640 720]; end
if nargin < 5, opts = struct(); end
gd = @(f,d) pumpkyn_local_default(opts, f, d);
markerSize  = gd('markerSize', [26 18 9]);
tailWidth   = gd('tailWidth',  2.5);
tulipColor  = gd('tulipColor', [0.38 0.38 0.38]);   % faint: a reference, not a subject
tulipWidth  = gd('tulipWidth', 0.7);
zoomPad     = gd('zoomPad',    0.25);               % right-panel framing margin
cBurn       = gd('tailColorBurn',  [0.95 0.25 0.15]);
cCoast      = gd('tailColorCoast', [0.25 0.60 1.00]);

FPS = 20;
PP  = fullfile(getenv('HOME'),'Desktop','proj7','external','pumpkynPie');
assert(isfolder(PP), 'record_tulip_dual_animation:noPumpkyn', 'no pumpkynPie at %s', PP);
here0 = pwd;  cleaner = onCleanup(@() cd(here0)); %#ok<NASGU>
cd(PP);  startup();
params = pumpkynPie.cr3bp.getCR3BParams();

if nargin < 2 || isempty(outFile)
    [d,n] = fileparts(solPath);  if isempty(d), d = here0; n = 'tulip'; end
    outFile = fullfile(d, [n '_dual.mp4']);
end

%% ---- resample uniformly in physical time --------------------------------
tS = X(8,:) * params.tStar;
if nargin < 4 || isempty(dtSec), dtSec = (tS(end)-tS(1))/1200; end
tU = (tS(1):dtSec:tS(end)).';
if tU(end) < tS(end), tU(end+1) = tS(end); end
xU = interp1(tS(:), X(1:6,:).', tU, 'linear');
sU = interp1(tS(:), double(U(4,:) > 0.5).', tU, 'nearest') > 0.5;

swSrc = sum(abs(diff(double(U(4,:) > 0.5))));
swRes = sum(abs(diff(double(sU(:).'))));
if swRes < swSrc
    warning('record_tulip_dual_animation:switchesLost', ...
        'dt = %.0f s kept %d of %d switches; reduce dtSec for a faithful tail.', ...
        dtSec, swRes, swSrc);
end
tailSeconds = gd('tailSeconds', 4*86400);
nTail = max(2, round(tailSeconds / dtSec));

jd0 = pumpkyn.util.juliandate('01/01/2000 00:00:00');
jd  = jd0 + (tU - tU(1))/86400;
fprintf('dual animation: %d frames (dt %.0f s), tail %d frames (%.1f d), %d switches\n', ...
        numel(tU), dtSec, nTail, tailSeconds/86400, swRes);

%% ---- the target tulip, for the white reference on the right -------------
% Same call as cr3bp_common/gto_tulip_endpoints -- this is OUR target orbit,
% not a generic one.
[tauT, x0T] = pumpkyn.cr3bp.getTulip((5/6)*2*pi, 7, -1, 1e-12);
x0T         = pumpkyn.cr3bp.cont_np(x0T, tauT, params.muStar, 1e-12);
[tT, xT]    = pumpkyn.cr3bp.prop(tauT, x0T, params.muStar);
jdT         = jd0 + tT*params.tStar/86400;

%% ---- build the two scenes ----------------------------------------------
% NOTE the up vector is set PER PANEL below, not here. The top-down view looks
% along -z, so an up vector of [0 0 1] is parallel to the view direction and the
% camera is degenerate -- pumpkyn warns "Camera up vector must be independent
% from camera direction" and the fit collapses (measured: the panel rendered as
% a single body filling the frame). [0 1 0] is the correct up for a top-down.
% RIGHT-PANEL FRAMING. The demo frames the Moon because the only thing in its
% scene is the tulip. Ours also holds a spiral reaching back to Earth, and the
% animator's 'trajectory' fit accumulates BOTH the orbit histories AND every
% satellite position history (SatelliteAnimator.getTrajectoryBoundsInCameraFrame)
% -- so no fit setting can frame the tulip while the spacecraft still flies the
% full transfer. The box is therefore pinned by hand to the tulip's own extent
% and the fit turned off, which reproduces the demo's zoom exactly. The
% spacecraft simply flies in from off-frame, which is what was asked for.
bb   = [min(xT(:,1:3),[],1); max(xT(:,1:3),[],1)].';    % 3x2
half = max(0.5*(bb(:,2)-bb(:,1))) * (1 + zoomPad);      % cube: aspect stays 1:1
ctr  = 0.5*(bb(:,1)+bb(:,2));
tulipBox = [ctr - half, ctr + half];

common = { 'Epoch', jd0, 'CR3BPParameters', params, 'FigurePosition', panelSize, ...
           'Background','stars', 'BackgroundColor','k', 'Quality','interactive', ...
           'ViewingAngle',[], 'CameraOrbit',[0 0], ...
           'DataAspectRatio',[1 1 1], 'FitPadding',0.05, ...
           'ShowLogo',true, 'MarkerSize', markerSize, ...
           'EnableSparkles', false, ...       % requested: no yellow sparkle
           'TrailMode','distance', 'TrailLength', eps, ... % kill the built-in
           'TrailPointCount', 4, 'TrailLineWidth', 0.5, 'EnableTrailColor', false };

saL = pumpkynPie.plot.SatelliteAnimator(common{:}, 'View',[0 90], ...
        'CameraUpVector',[0 1 0], 'FitMode','trajectory');   % top-down, whole system
saR = pumpkynPie.plot.SatelliteAnimator(common{:}, 'View',[90 0], ...
        'CameraUpVector',[0 0 1], 'FitMode','none', ...      % Moon plane, as the demo
        'AxesLimits', tulipBox, 'CameraTarget', ctr(:).');

% Invisible full path defines the LEFT fit box, so that camera does not jump as
% the tail moves. Colour 'none' keeps it out of the picture.
saL.addOrbit(jd, xU, 'Name','extent', 'Color','none');   % LEFT only -- on the
% right it would be ignored anyway (fit off), and it would cost memory.

% The target tulip, white, right panel only -- so the arrival is visible as the
% spacecraft flies into a drawn orbit rather than into empty space.
saR.addOrbit(jdT, xT(:,1:6), 'Name','target tulip', ...
             'Color', tulipColor, 'LineWidth', tulipWidth);

saL.addSatellite(jd, xU, 'Name','s/c', 'MarkerColor',[1 1 1], ...
                 'TrailColor',[1 1 1], 'SparkleColor',[1 1 1]);
saR.addSatellite(jd, xU, 'Name','s/c', 'MarkerColor',[1 1 1], ...
                 'TrailColor',[1 1 1], 'SparkleColor',[1 1 1]);
saL.initialize();  saR.initialize();

for sa = [saL saR]
    ax = sa.Axes;
    try, ax.Toolbar.Visible = 'off';      catch, end %#ok<CTCH>
    try, axtoolbar(ax,'none');            catch, end %#ok<CTCH>
    try, disableDefaultInteractivity(ax); catch, end %#ok<CTCH>
end

% The fading tails: one patch per panel. EdgeColor/EdgeAlpha 'interp' read the
% per-vertex CData/AlphaData we rewrite each frame.
tailL = local_make_tail(saL.Axes, tailWidth);
tailR = local_make_tail(saR.Axes, tailWidth);
drawnow;

%% ---- record ------------------------------------------------------------
vw = VideoWriter(outFile,'MPEG-4');  vw.FrameRate = FPS;  vw.Quality = 92;
open(vw);  closer = onCleanup(@() close(vw)); %#ok<NASGU>

t0 = tic;
for k = 1:numel(jd)
    saL.update(jd(k));  saR.update(jd(k));
    local_set_tail(tailL, xU, sU, k, nTail, cBurn, cCoast);
    local_set_tail(tailR, xU, sU, k, nTail, cBurn, cCoast);
    drawnow;
    fL = getframe(saL.Figure);  fR = getframe(saR.Figure);
    h  = min(size(fL.cdata,1), size(fR.cdata,1));
    writeVideo(vw, [fL.cdata(1:h,:,:), fR.cdata(1:h,:,:)]);
    if mod(k,250)==0, fprintf('  %4d/%4d (%.0f s)\n', k, numel(jd), toc(t0)); end
end
clear closer;

v = VideoReader(outFile);
[od,on] = fileparts(outFile);
for frac = [0.15 0.5 0.85]
    imwrite(read(v, max(1,round(v.NumFrames*frac))), ...
            fullfile(od, sprintf('%s_p%02d.png', on, round(frac*100))));
end
fprintf('Done: %s (%.0f s render, %d frames, %.1f s runtime)\n', ...
        outFile, toc(t0), numel(jd), v.Duration);
end

% ---------------------------------------------------------------------------
function h = local_make_tail(ax, lw)
% LOCAL_MAKE_TAIL  Empty polyline patch supporting per-vertex colour and alpha.
% INPUTS:  ax - axes; lw - line width      OUTPUTS: h - patch handle
h = patch('Parent', ax, 'Vertices', nan(2,3), 'Faces', [1 2], ...
          'FaceColor','none', 'EdgeColor','interp', 'EdgeAlpha','interp', ...
          'FaceVertexCData', [1 1 1; 1 1 1], ...
          'FaceVertexAlphaData', [0; 0], 'AlphaDataMapping','none', ...
          'LineWidth', lw, 'PickableParts','none', 'HitTest','off');
end

% ---------------------------------------------------------------------------
function local_set_tail(h, xU, sU, k, nTail, BURN, COAST)
% LOCAL_SET_TAIL  Point the tail at the last nTail samples ending at k.
% Colour carries the control state, alpha carries recency (oldest transparent).
% INPUTS: h - patch; xU - Nx6 states; sU - logical throttle; k - frame; nTail
% OUTPUTS: none
i0 = max(1, k - nTail + 1);
idx = i0:k;
n   = numel(idx);
if n < 2
    set(h, 'Vertices', nan(2,3), 'Faces', [1 2], ...
           'FaceVertexCData', [1 1 1; 1 1 1], 'FaceVertexAlphaData', [0;0]);
    return
end
C = repmat(COAST, n, 1);  C(sU(idx), :) = repmat(BURN, sum(sU(idx)), 1);
A = linspace(0, 1, n).';            % oldest fully transparent -> newest opaque
set(h, 'Vertices', xU(idx,1:3), 'Faces', [(1:n-1).' (2:n).'], ...
       'FaceVertexCData', C, 'FaceVertexAlphaData', A);
end

% ---------------------------------------------------------------------------
function v = pumpkyn_local_default(s, f, d)
% Local option default. INPUTS: s struct, f field, d default. OUTPUTS: v
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
