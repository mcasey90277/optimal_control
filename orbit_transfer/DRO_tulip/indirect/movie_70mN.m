function outStem = movie_70mN(outStem)
%% Purpose:
%
%   Render the CERTIFIED 70 mN minimum-time DRO -> tulip transfer (FINDINGS
%   section 32: t_f = 26.436 d, dV = 1.136 km/s, 18.12 kg) as a single-panel
%   pumpkyn-lit movie, for the cislunar poster and the IEEE paper.
%
%   Frame size is forced to exactly 1280 x 720. A frame size not divisible
%   by 16 makes H.264 shear the picture into diagonal colour streaks, which
%   reads as a render bug and is not one (memory: matlab-movie-diagonal-streaks).
%
%% Inputs:
%
%  outStem                  char (optional)         output stem; default
%                                                   direct/results/mintime_70mN_movie
%
%% Outputs:
%
%  outStem                  char                    stem actually written
%                                                   (.mp4 and .gif)
%
%% Revision History:
%  M. Casey                                                   (c) 09/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(outStem)
    outStem = fullfile(here, '..', 'direct', 'results', 'mintime_70mN_movie');
end
addpath(fullfile(fileparts(here), '..', 'costate_common'));

C = load(fullfile(here, 'results', 'mintime_70mN_certified.mat'));
D = load(fullfile(here, 'results', 'mintime_70mN_direct.mat'));
Q = load(fullfile(here, '..', 'direct', 'results', 'thrust_ladder_12x12.mat'));
ob = Q.meta;  mu = ob.muStar;  lStar = ob.lStar;  tStar = ob.tStar;
z = C.z;  Tnd = D.Tnd;  cnd = D.cnd;

% the two periodic orbits and the transfer, flown in the certified field
[tD, rvD, tT, rvT] = ladder_endpoints(ob);
[tj, yj] = pumpkyn.cr3bp.tfMinProp(z(8), [D.rv0(1:6); 1; z(1:7)], Tnd, cnd, mu);
[tu, iu] = unique(tj);  yj = yj(iu, :);
tfDays = z(8)*tStar/86400;
mf = yj(end, 7);
dV = cnd*log(1/mf)*lStar/tStar;
fprintf('70 mN transfer: t_f = %.3f d, dV = %.4f km/s, prop = %.2f kg\n', ...
        tfDays, dV, ob.m0kg*(1 - mf));

%% Ephemeris guard: the R2026a Aerospace Toolbox puts planetEphemeris on the
%  path without its data package, and pumpkyn's 'auto' selector picks it on
%  PRESENCE, blocking the analytic Earth-Moon fallback that actually works
%  here. Hide it for this render only; the path is restored on exit.
if exist('planetEphemeris', 'file')
    try
        planetEphemeris(juliandate(datetime(2030,1,1)), 'Earth', 'Moon');
    catch
        aeroDir = fileparts(which('planetEphemeris'));
        rmpath(aeroDir);
        cleanupPath = onCleanup(@() addpath(aeroDir)); %#ok<NASGU>
        warning('movie_70mN:noEphemData', ...
            ['planetEphemeris found but its data package is missing; ' ...
             'using pumpkyn''s analytic Earth-Moon model for this render.']);
    end
end

params = struct('muStar', mu, 'lStar', lStar, 'tStar', tStar);
epoch  = juliandate(datetime(2030,1,1,0,0,0));
fig = figure('Color','k', 'Position',[60 60 1280 720], ...
             'Visible','off', 'InvertHardcopy','off');
anim = pumpkynPie.plot.SatelliteAnimator(fig, ...
    'Epoch', epoch, 'CR3BPParameters', params, 'ShowEarth', true, ...
    'FigureTheme', 'dark', 'BodyLighting', 'lambertian', ...
    'FitMode', 'trajectory', 'TrailMode', 'distance', ...
    'EnableTrailColor', true, 'EnableSparkles', false, ...
    'EnableSatelliteMarker', true, 'CameraUpVector', [0 0 1], 'View', [-35 22]);
anim.addOrbit(epoch + tD*tStar/86400, rvD, ...
    'Name', sprintf('DRO (\\tau = %.2f)', ob.tauDRO), ...
    'Color', [0.25 0.60 0.35], 'LineWidth', 1.0);
anim.addOrbit(epoch + tT*tStar/86400, rvT, ...
    'Name', sprintf('%d-petal tulip', ob.NpTulip), ...
    'Color', [0.75 0.28 0.28], 'LineWidth', 1.0);
anim.addSatellite(epoch + tu*tStar/86400, yj(:,1:6), ...
    'Name', 'Minimum-time transfer', ...
    'MarkerColor', [1.00 0.92 0.55], 'TrailColor', [0.35 0.72 1.00]);
anim.initialize();
camzoom(1.35);
annotation(fig, 'textbox', [0.02 0.905 0.96 0.085], 'String', ...
    {sprintf(['Minimum-time DRO \\rightarrow %d-petal tulip   |   ' ...
              '%.0f mN, I_{sp} %g s, %g kg'], ob.NpTulip, Tnd*ob.m0kg* ...
              (lStar*1000)/tStar^2*1000, 900, ob.m0kg), ...
     sprintf(['t_f = %.2f days    \\DeltaV = %.3f km/s    propellant %.1f kg ' ...
              '(%.1f%%)'], tfDays, dV, ob.m0kg*(1-mf), 100*(1-mf))}, ...
    'Color','w', 'EdgeColor','none', 'FontSize',13, ...
    'HorizontalAlignment','center');
annotation(fig, 'textbox', [0.02 0.055 0.40 0.045], 'String', ...
    'certified: tfMin |\Deltaz| = 2e-10, conjugate PASS, gates pass', ...
    'Color',[0.55 0.85 0.55], 'EdgeColor','none', 'FontSize',10, ...
    'HorizontalAlignment','left');
clockBox = annotation(fig, 'textbox', [0.02 0.02 0.96 0.05], 'String','', ...
    'Color',[1 0.92 0.55], 'EdgeColor','none', 'FontSize',14, ...
    'FontWeight','bold', 'HorizontalAlignment','center');

nFrames = 300;  fps = 30;
vw = VideoWriter([outStem '.mp4'], 'MPEG-4');
vw.FrameRate = fps;  vw.Quality = 95;  open(vw);
gifFile = [outStem '.gif'];  wroteGif = false;
for kf = 1:nFrames
    tNow = z(8)*(kf-1)/(nFrames-1);
    anim.update(epoch + tNow*tStar/86400);
    clockBox.String = sprintf('t = %.2f / %.2f days', tNow*tStar/86400, tfDays);
    img = frame1280(fig);
    writeVideo(vw, img);
    if mod(kf, 3) == 1                                  % gif at 10 fps
        [A, map] = rgb2ind(img, 256);
        if ~wroteGif
            imwrite(A, map, gifFile, 'gif', 'LoopCount', inf, 'DelayTime', 0.1);
            wroteGif = true;
        else
            imwrite(A, map, gifFile, 'gif', 'WriteMode','append', 'DelayTime', 0.1);
        end
    end
end
close(vw);  close(fig);
fprintf('  movie -> %s.mp4 / .gif  (%d frames, 1280x720, %.2f d)\n', ...
        outStem, nFrames, tfDays);
end

% ------------------------------------------------------------------------
function img = frame1280(fig)
% FRAME1280  getframe forced to exactly 720 x 1280 by index resampling.
% A size not divisible by 16 makes H.264 shear the frame diagonally.
% INPUTS: fig.  OUTPUTS: img [720 x 1280 x 3 uint8].
F = getframe(fig);
ri = max(1, min(size(F.cdata,1), round(linspace(1, size(F.cdata,1), 720))));
ci = max(1, min(size(F.cdata,2), round(linspace(1, size(F.cdata,2), 1280))));
img = F.cdata(ri, ci, :);
end
