function outStem = costate_catalog_extremes_movies(cat_, filt, metric, outStem)
%% Purpose:
%
%   Renders the two EXTREME transfers found by costate_catalog_extremes as a
%   side-by-side MOVIE (.mp4 and .gif) in the pumpkynPie house style: the
%   lit dark CR3BP scene (textured Moon, Earth, colored terminal orbits,
%   colored trail, satellite marker) built by
%   pumpkynPie.plot.SatelliteAnimator -- the same scene as
%   dro_scene_movie.m and the pumpkynPie demos.
%
%   Left panel = the minimum, right = the maximum, separated by a vertical
%   divider, and BOTH RUN ON ONE CLOCK: frames advance on a shared time
%   base, so the shorter transfer arrives early and freezes on its final
%   frame while the longer one plays out. The clock (days) is shown in both
%   panels with identical values.
%
%  ASSUMPTIONS / NOTES:
%
% • Each output frame is EXACTLY 1280 x 720 (two 637-px panels + a 6-px
%   divider), enforced by index resampling -- the known cure for the H.264
%   diagonal-shear artifact. The .gif is written at the same size.
% • The trajectories are FLOWN from the stored costates (tfMinProp): what
%   you watch is the actual PMP solution.
% • Needs pumpkynPie.plot.SatelliteAnimator; if it is not on the path the
%   function warns and returns '' (same graceful degradation as
%   dro_scene_movie).
%
%% Inputs:
%
%  cat_                     struct                  A compact costate catalog
%                                                   (any family: dro_tulip,
%                                                   halo_tulip, ...)
%
%  filt                     struct                  Same filter as
%                                                   costate_catalog_extremes
%                                                   (.tauDRO .Np .thrustN,
%                                                   any subset; [] = all)
%
%  metric                   char                    'deltaV' (default) or
%                                                   'time'
%
%  outStem                  char                    Output stem for
%                                                   '<stem>.mp4' and
%                                                   '<stem>.gif'  [default
%                                                   'catalog_extremes_movie']
%
%% Outputs:
%
%  outStem                  char                    The stem written, or ''
%                                                   if the animator was
%                                                   unavailable
%
%% Revision History:
%  M. Casey                                                   (c) 08/06/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: phasing extremes on the mid sheet at 5 N. Runs on whichever
   %compact catalog sits in the current folder:
     F = dir('costate_catalog_*.mat');
     assert(~isempty(F), ...
            'demo: no costate_catalog_*.mat in the current folder');
       L = load(F(1).name);
      fn = fieldnames(L);
    cat_ = L.(fn{1});
   fprintf('demo catalog: %s\n', F(1).name);
      sh = cat_.sheets(ceil(numel(cat_.sheets)/2));
 outStem = costate_catalog_extremes_movies(cat_, ...
               struct('tauDRO',sh.tauDRO,'Np',sh.Np,'thrustN',5));
     return;
end
if ~exist('filt','var'),    filt = struct(); end
if ~exist('metric','var') || isempty(metric), metric = 'deltaV'; end
if ~exist('outStem','var'), outStem = 'catalog_extremes_movie'; end

if isempty(which('pumpkynPie.plot.SatelliteAnimator'))
    warning('costate_catalog_extremes_movies:noAnimator', ...
        'pumpkynPie.plot.SatelliteAnimator not on the path -- no movie.');
    outStem = '';  return
end

%% The two extremes, and their flown trajectories:
[eMin, eMax] = costate_catalog_extremes(cat_, filt, false, metric);
      mu = cat_.constants.muStar;
   tStar = cat_.constants.tStar_s;
   lStar = cat_.constants.lStar_km;
     cnd = cat_.thruster.c_nd;
    m0kg = cat_.thruster.m0_kg;
     ndT = @(TN) (TN/m0kg)*tStar^2/(lStar*1000);
  params = struct('muStar',mu, 'lStar',lStar, 'tStar',tStar);
   epoch = juliandate(datetime(2030,1,1,0,0,0));
if strncmpi(metric,'t',1), lab = {'FASTEST','SLOWEST'};
else,                      lab = {'MIN \DeltaV','MAX \DeltaV'}; end

E = {eMin, eMax};
for k = 1:2
    e = E{k};
    [tauD, rvD] = get_family_orbit(e.dep_family, e.dep_params);
    tT0 = 2*pi*(e.Np-2)/(e.Np-1);
    [~, rvT0] = pumpkyn.cr3bp.getTulip(tT0, e.Np, -1);
    rvT0 = pumpkyn.cr3bp.cont_np(rvT0, tT0, mu, 1e-12);
    [tauT, rvT] = pumpkyn.cr3bp.prop(tT0, rvT0, mu);
    rv0 = interp1(tauD/tauD(end), rvD, e.dep_frac, 'spline');
    [tj, yj] = pumpkyn.cr3bp.tfMinProp(e.z8(8), [rv0(1:6)'; 1; e.z8(1:7)], ...
                   ndT(e.thrustN), cnd, mu);
    [tu, iu] = unique(tj);

    %% One pumpkyn-style lit scene per panel:
    S(k).fig = figure('Color','k', 'Position',[60 60 640 720], ...
                      'Visible','off', 'InvertHardcopy','off');
    S(k).anim = pumpkynPie.plot.SatelliteAnimator(S(k).fig, ...
        'Epoch',              epoch, ...
        'CR3BPParameters',    params, ...
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
    S(k).anim.addOrbit(epoch + tauD*tStar/86400, rvD, ...
        'Name','Departure orbit', 'Color',[0.25 0.60 0.35], 'LineWidth',0.9);
    S(k).anim.addOrbit(epoch + tauT*tStar/86400, rvT, ...
        'Name','Arrival tulip', 'Color',[0.75 0.28 0.28], 'LineWidth',0.9);
    S(k).anim.addSatellite(epoch + tu*tStar/86400, yj(iu,1:6), ...
        'Name', sprintf('%s transfer', lab{k}), ...
        'MarkerColor',[1.00 0.92 0.55], 'TrailColor',[0.35 0.72 1.00]);
    S(k).anim.initialize();
    camzoom(1.4);
    annotation(S(k).fig, 'textbox', [0.02 0.90 0.96 0.09], 'String', ...
        {sprintf('%s:  %s \\tau=%.2f, N_p=%d, %.2f N', ...
                 lab{k}, upper(e.dep_family), e.tauDRO, e.Np, e.thrustN), ...
         sprintf('\\DeltaV %.3f km/s,  t_f %.3f d', e.deltaV_kms, e.tf_days)}, ...
        'Color','w', 'EdgeColor','none', 'FontSize',11, ...
        'HorizontalAlignment','center');
    S(k).clock = annotation(S(k).fig, 'textbox', [0.02 0.02 0.96 0.05], ...
        'String','', 'Color',[1 0.92 0.55], 'EdgeColor','none', ...
        'FontSize',13, 'FontWeight','bold', 'HorizontalAlignment','center');
    S(k).tf = e.z8(8);
end
tEnd = max(S(1).tf, S(2).tf);            % the shared clock's last tick

%% Frames: shared clock; shorter panel freezes at its final pose:
 nFrames = 240;
     fps = 30;
      pw = 637;                          % panel width; 2*637 + 6 = 1280
 divider = uint8(90*ones(720, 1280-2*pw, 3));
vw = VideoWriter([outStem '.mp4'], 'MPEG-4');
vw.FrameRate = fps;  vw.Quality = 95;  open(vw);
cleanupVW = onCleanup(@() close(vw)); %#ok<NASGU>
gifFile = [outStem '.gif'];
for kf = 1:nFrames
    tNow = tEnd*(kf-1)/(nFrames-1);
    for k = 1:2
        S(k).anim.update(epoch + min(tNow, S(k).tf)*tStar/86400);
        set(S(k).clock, 'String', sprintf('t = %.2f days', tNow*tStar/86400));
    end
    drawnow limitrate
    img = [panel_frame(S(1).fig, pw), divider, panel_frame(S(2).fig, pw)];
    writeVideo(vw, img);
    if mod(kf,3) == 1                    % gif at 10 fps, same 1280x720 size
        [A, cmap] = rgb2ind(img, 256);
        if kf == 1
            imwrite(A, cmap, gifFile, 'gif', 'LoopCount', inf, 'DelayTime', 0.1);
        else
            imwrite(A, cmap, gifFile, 'gif', 'WriteMode','append', 'DelayTime', 0.1);
        end
    end
end
close(S(1).fig);  close(S(2).fig);
fprintf('  movies -> %s.mp4 / .gif  (%d frames, 1280x720, shared clock to %.2f d)\n', ...
        outStem, nFrames, tEnd*tStar/86400);
end

% ---------------------------------------------------------------------------
function img = panel_frame(fig, w)
% PANEL_FRAME  getframe forced to 720 x w by index resampling (the H.264
% shear cure from dro_scene_movie, generalized to a panel width).
% INPUTS: fig; w target width.  OUTPUTS: img [720 x w x 3 uint8].
F = getframe(fig);
ri = max(1, min(size(F.cdata,1), round(linspace(1, size(F.cdata,1), 720))));
ci = max(1, min(size(F.cdata,2), round(linspace(1, size(F.cdata,2), w))));
img = F.cdata(ri, ci, :);
end
