function outFile = dro_scene_movie(o, p, outFile, mopts)
% DRO_SCENE_MOVIE  pumpkynPie-style lit CR3BP scene of the direct transfer.
%
% Mirrors the scene pumpkynPie builds in demos/lowThrustDRO2Tulip.m and animates
% it the way Demos/lowThrust_constellation_animation.m does, but driven by OUR
% direct-collocation solution rather than the indirect one. Nothing is written
% into pumpkyn or pumpkynPie -- this only calls into them.
%
% GRACEFUL DEGRADATION. The lit scene needs pumpkynPie.plot.SatelliteAnimator
% (and, through it, pumpkyn's showMoon). If that class is not on the path the
% function returns '' without erroring, and the caller should fall back to
% dro_traj_movie, which uses base MATLAB only.
%
% INPUTS:
%   o       - solution struct from casadi_mintime_dro
%   p       - params from dro_tulip_endpoints (.muStar .lStar .tStar)
%   outFile - output path, .mp4 [char]
%   mopts   - struct (optional):
%             .nDense   spline samples per collocation interval [default 6]
%             .fps      frame rate                              [default 30]
%             .nFrames  frames                                  [default 360]
%             .showEarth                                        [default true]
%             .sparkles keep the sparkle particles              [default false]
%             .zoom     camzoom factor after initialize         [default 1.5]
%
% OUTPUTS:
%   outFile - the file written, or '' if the animator was unavailable [char]
%
% REFERENCES:
%   [1] pumpkynPie/src/+pumpkynPie/+plot/SatelliteAnimator.m
%   [2] pumpkynPie/Demos/lowThrust_constellation_animation.m

if nargin < 4, mopts = struct(); end
g = @(f,v) local_default(mopts, f, v);
nDense    = g('nDense', 6);
fps       = g('fps', 30);
nFrames   = g('nFrames', 360);
showEarth = g('showEarth', true);
sparkles  = g('sparkles', false);
zoomFac   = g('zoom', 1.5);

if isempty(which('pumpkynPie.plot.SatelliteAnimator'))
    warning('dro_scene_movie:noAnimator', ...
        'pumpkynPie.plot.SatelliteAnimator not on the path -- skipping the lit scene.');
    outFile = '';  return
end

mu = p.muStar;
params = struct('muStar', mu, 'lStar', p.lStar, 'tStar', p.tStar);
epoch  = juliandate(datetime(2030,1,1,0,0,0));

% --- densify the transfer, and build the two terminal orbits ---------------
tN = o.s(:).' * o.tf;
tD = linspace(tN(1), tN(end), nDense*(numel(tN)-1) + 1).';
XD = interp1(tN, o.X(1:6,:).', tD, 'spline');            % nD x 6
jdD = epoch + tD*p.tStar/86400;

[tauT, rvT] = local_orbit(mu, 'tulip');
[tauD, rvD] = local_orbit(mu, 'dro');

fig = figure('Color','k','Position',[80 80 1280 720],'Visible','off', ...
             'InvertHardcopy','off');
animator = pumpkynPie.plot.SatelliteAnimator(fig, ...
    'Epoch',              epoch, ...
    'CR3BPParameters',    params, ...
    'ShowEarth',          showEarth, ...
    'FigureTheme',        'dark', ...
    'BodyLighting',       'lambertian', ...
    'FitMode',            'trajectory', ...
    'TrailMode',          'distance', ...
    'EnableTrailColor',   true, ...
    'EnableSparkles',     sparkles, ...
    'EnableSatelliteMarker', true, ...
    'CameraUpVector',     [0 0 1], ...
    'View',               [-35 22]);

if ~isempty(rvD)
    animator.addOrbit(epoch + tauD*p.tStar/86400, rvD, ...
        'Name','Departure DRO', 'Color',[0.25 0.60 0.35], 'LineWidth',0.9);
end
if ~isempty(rvT)
    animator.addOrbit(epoch + tauT*p.tStar/86400, rvT, ...
        'Name','Arrival tulip orbit', 'Color',[0.75 0.28 0.28], 'LineWidth',0.9);
end
animator.addSatellite(jdD, XD, 'Name','DRO to tulip (direct)', ...
    'MarkerColor',[1.00 0.92 0.55], 'TrailColor',[0.35 0.72 1.00]);
animator.initialize();
camzoom(zoomFac);

vw = VideoWriter(outFile,'MPEG-4');  vw.FrameRate = fps;  vw.Quality = 95;
open(vw);
cleanupVW = onCleanup(@() close(vw)); %#ok<NASGU>

jdFrames = linspace(jdD(1), jdD(end), nFrames);
for kf = 1:nFrames
    animator.update(jdFrames(kf));
    drawnow limitrate
    writeVideo(vw, local_frame(fig));
end
close(fig);
fprintf('  movie -> %s  (%d frames @ %d fps)\n', outFile, nFrames, fps);
end

% ---------------------------------------------------------------------------
function [tau, rv] = local_orbit(mu, which)
% LOCAL_ORBIT  A terminal orbit, propagated one period. Same construction as
% pumpkynPie demos/lowThrustDRO2Tulip.m: catalog getter, then cont_np.
% INPUTS: mu; which ('tulip' | 'dro')
% OUTPUTS: tau [Mx1], rv [Mx6] (both empty if the getter is unavailable)
tau = [];  rv = [];
try
    switch lower(which)
        case 'tulip'
            tau0 = 5*2*pi/6;
            [~, rv0] = pumpkyn.cr3bp.getTulip(tau0, 7, -1);
        case 'dro'
            tau0 = 1.0;
            [~, rv0] = pumpkynPie.cr3bp.getDRO(tau0);
        otherwise
            return
    end
    rv0 = pumpkyn.cr3bp.cont_np(rv0, tau0, mu, 1e-12);
    [tau, rv] = pumpkyn.cr3bp.prop(tau0, rv0, mu);
catch, end
end

% ---------------------------------------------------------------------------
function F = local_frame(fig)
% LOCAL_FRAME  getframe forced to 1280x720 so H.264 does not shear the frame.
% INPUTS: fig   OUTPUTS: F (struct with .cdata [720x1280x3 uint8])
F = getframe(fig);
if size(F.cdata,1) ~= 720 || size(F.cdata,2) ~= 1280
    ri = max(1, min(size(F.cdata,1), round(linspace(1, size(F.cdata,1), 720))));
    ci = max(1, min(size(F.cdata,2), round(linspace(1, size(F.cdata,2), 1280))));
    F.cdata = F.cdata(ri, ci, :);
end
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
% INPUTS: s; f; dflt   OUTPUTS: v
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
