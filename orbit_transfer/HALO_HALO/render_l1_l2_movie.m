function outStem = render_l1_l2_movie(rungN, probeMat)
%% Purpose:
%
%   Renders the L1 -> L2 halo-to-halo probe transfer (probe_l1_l2_halo) as
%   a pumpkyn-style MOVIE (.mp4 and .gif): the lit dark CR3BP scene
%   (textured Moon, Earth, both halo orbits, colored trail, satellite
%   marker) built by pumpkynPie.plot.SatelliteAnimator. Each frame is
%   exactly 1280 x 720 via index resampling (the H.264 shear cure).
%
%% Inputs:
%
%  rungN                    double                  Thrust rung to render
%                                                   (N) [default 1, the
%                                                   longest arc]
%
%  probeMat                 char                    Probe sheet file under
%                                                   direct/results/l1l2/
%                                                   [default
%                                                   'probe_L1_L2.mat']
%
%% Outputs:
%
%  outStem                  char                    Stem of the .mp4/.gif
%                                                   written under
%                                                   direct/results/l1l2/
%
%% Revision History:
%  M. Casey                                                   (c) 08/07/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if ~exist('rungN','var'), rungN = 1; end
if ~exist('probeMat','var'), probeMat = 'probe_L1_L2.mat'; end

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'costate_common'));
outDir = fullfile(here, 'direct', 'results', 'l1l2');
P = load(fullfile(outDir, probeMat));

kr = find(P.rungs == rungN, 1);
assert(~isempty(kr) && P.OK(1,1,kr), 'rung %.3g N not solved in the probe', rungN);

   mu = P.meta.muStar;  tStar = P.meta.tStar;  lStar = P.meta.lStar;
   g0 = 9.80665*tStar^2/(1000*lStar);
  cnd = (P.meta.ispS/tStar)*g0;
  Tnd = (rungN/P.meta.m0kg)*tStar^2/(lStar*1000);
   z8 = P.Z8(:,1,1,kr);
   tf = P.TF(1,1,kr);

%% Endpoint halos and the flown transfer (the actual PMP solution):
[tD, rvD] = get_family_orbit(P.meta.depFamily, P.meta.depParams);
[tA, rvA] = get_family_orbit(P.meta.arrFamily, P.meta.arrParams);
rv0 = rvD(1,:);                                  % probe cell: phase 0 / 0
[tj, yj] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(1:6)'; 1; z8(1:7)], Tnd, cnd, mu);
[tu, iu] = unique(tj);
mf = 1 - Tnd*tf/cnd;
dV = cnd*log(1/mf)*lStar/tStar;

%% The pumpkyn-style lit scene:
epoch = juliandate(datetime(2030,1,1,0,0,0));
fig = figure('Color','k', 'Position',[60 60 1280 720], ...
             'Visible','off', 'InvertHardcopy','off');
anim = pumpkynPie.plot.SatelliteAnimator(fig, ...
    'Epoch',              epoch, ...
    'CR3BPParameters',    struct('muStar',mu,'lStar',lStar,'tStar',tStar), ...
    'ShowEarth',          true, ...
    'FigureTheme',        'dark', ...
    'BodyLighting',       'lambertian', ...
    'FitMode',            'trajectory', ...
    'TrailMode',          'distance', ...
    'EnableTrailColor',   true, ...
    'EnableSparkles',     false, ...
    'EnableSatelliteMarker', true, ...
    'CameraUpVector',     [0 0 1], ...
    'View',               [-35 22]);
anim.addOrbit(epoch + tD*tStar/86400, rvD, ...
    'Name','L1 southern halo', 'Color',[0.25 0.60 0.35], 'LineWidth',1.1);
anim.addOrbit(epoch + tA*tStar/86400, rvA, ...
    'Name','L2 southern halo', 'Color',[0.75 0.28 0.28], 'LineWidth',1.1);
anim.addSatellite(epoch + tu*tStar/86400, yj(iu,1:6), ...
    'Name','L1 -> L2 transfer', ...
    'MarkerColor',[1.00 0.92 0.55], 'TrailColor',[0.35 0.72 1.00]);
anim.initialize();
camzoom(1.35);
annotation(fig, 'textbox', [0.02 0.90 0.96 0.09], 'String', ...
    {sprintf('L1 halo (\\tau=%.2f) \\rightarrow L2 halo (\\tau=%.2f),  %.3g N', ...
             P.meta.depParams.tau, P.meta.arrParams.tau, rungN), ...
     sprintf('\\DeltaV %.3f km/s,  t_f %.2f d,  %.1f kg propellant', ...
             dV, tf*tStar/86400, (1-mf)*P.meta.m0kg)}, ...
    'Color','w', 'EdgeColor','none', 'FontSize',12, ...
    'HorizontalAlignment','center');
clockBox = annotation(fig, 'textbox', [0.02 0.02 0.96 0.05], 'String','', ...
    'Color',[1 0.92 0.55], 'EdgeColor','none', 'FontSize',13, ...
    'FontWeight','bold', 'HorizontalAlignment','center');

%% Frames (exact 1280 x 720 by index resampling):
nFrames = 240;  fps = 30;
[~, probeStem] = fileparts(probeMat);
tag = strrep(probeStem, 'probe_L1_L2', 'l1_l2_halo');
outStem = fullfile(outDir, sprintf('%s_%gN', tag, rungN));
vw = VideoWriter([outStem '.mp4'], 'MPEG-4');
vw.FrameRate = fps;  vw.Quality = 95;  open(vw);
cleanupVW = onCleanup(@() close(vw));
gifFile = [outStem '.gif'];
for kf = 1:nFrames
    tNow = tf*(kf-1)/(nFrames-1);
    anim.update(epoch + tNow*tStar/86400);
    set(clockBox, 'String', sprintf('t = %.2f days', tNow*tStar/86400));
    drawnow limitrate
    fr = getframe(fig);  img = fr.cdata;
    [h, w, ~] = size(img);
    img = img(round(linspace(1, h, 720)), round(linspace(1, w, 1280)), :);
    writeVideo(vw, img);
    if mod(kf,3) == 1
        [A, cmap] = rgb2ind(img, 256);
        if kf == 1
            imwrite(A, cmap, gifFile, 'gif', 'LoopCount', inf, 'DelayTime', 0.1);
        else
            imwrite(A, cmap, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end
end
close(fig);
fprintf('movie -> %s.mp4 / .gif  (%d frames, 1280x720, t_f %.2f d)\n', ...
        outStem, nFrames, tf*tStar/86400);
end
