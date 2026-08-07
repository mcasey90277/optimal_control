function costate_catalog_extremes_movies(cat_, filt, metric, outStem)
%% Purpose:
%
%   Renders side-by-side MOVIES (both .mp4 and .gif) of the two EXTREME
%   transfers found by costate_catalog_extremes: the minimum and the maximum
%   of the chosen metric under the same filter. Left panel = minimum, right
%   panel = maximum, separated by a vertical dividing line, and BOTH RUN ON
%   THE SAME CLOCK: the shorter transfer finishes early and freezes on its
%   final frame while the longer one plays out, so the pacing difference is
%   visible, not just stated.
%
%   Each panel shows its sheet's DRO (black dashed) and tulip (red dashed),
%   the transfer trail growing in blue, the spacecraft as a moving dot, and
%   departure/arrival markers. A shared clock (days) runs across the top.
%
%  ASSUMPTIONS / NOTES:
%
% • Frames are locked to 1280 x 720 (dimensions divisible by 16). H.264
%   shears frames whose size is not; this is the known cure.
% • The transfers are flown by propagating each extreme's stored costates
%   (tfMinProp) -- what you watch is the actual PMP trajectory, not an
%   interpolation of saved states.
%
%% Inputs:
%
%  cat_                     struct                  costate_catalog_dro_tulip
%
%  filt                     struct                  Same filter as
%                                                   costate_catalog_extremes:
%                                                   any subset of .tauDRO,
%                                                   .Np, .thrustN; [] = the
%                                                   whole catalog
%
%  metric                   char                    'deltaV' (default) or
%                                                   'time'
%
%  outStem                  char                    Output path stem;
%                                                   '<stem>.mp4' and
%                                                   '<stem>.gif' are written
%                                                   [default
%                                                   'catalog_extremes_movie']
%
%% Outputs:
%
%  none (two movie files written)
%
%% Revision History:
%  M. Casey                                                   (c) 08/06/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: phasing extremes on one sheet at 5 N:
       L = load('costate_catalog_dro_tulip.mat');
    cat_ = L.costate_catalog_dro_tulip;
     costate_catalog_extremes_movies(cat_, ...
         struct('tauDRO',2.0,'Np',7,'thrustN',5));
     return;
end
if ~exist('filt','var'),    filt = struct(); end
if ~exist('metric','var') || isempty(metric), metric = 'deltaV'; end
if ~exist('outStem','var'), outStem = 'catalog_extremes_movie'; end

%% The two extremes (report prints once; no static plot):
[eMin, eMax] = costate_catalog_extremes(cat_, filt, false, metric);

%% Fly both and gather everything the frames need:
      mu = cat_.constants.muStar;
   tStar = cat_.constants.tStar_s;
   lStar = cat_.constants.lStar_km;
     cnd = cat_.thruster.c_nd;
    m0kg = cat_.thruster.m0_kg;
     ndT = @(TN) (TN/m0kg)*tStar^2/(lStar*1000);
E = {eMin, eMax};
for k = 1:2
    e = E{k};
    [~, rvD0] = pumpkynPie.cr3bp.getDRO(e.tauDRO);
    rvD0 = pumpkyn.cr3bp.cont_np(rvD0, e.tauDRO, mu, 1e-12);
    [~, P(k).rvD] = pumpkyn.cr3bp.prop(e.tauDRO, rvD0, mu);       %#ok<*AGROW>
    tauT = 2*pi*(e.Np-2)/(e.Np-1);
    [~, rvT0] = pumpkyn.cr3bp.getTulip(tauT, e.Np, -1);
    rvT0 = pumpkyn.cr3bp.cont_np(rvT0, tauT, mu, 1e-12);
    [tT, P(k).rvT] = pumpkyn.cr3bp.prop(tauT, rvT0, mu);
    rv0 = interp1(linspace(0,1,size(P(k).rvD,1)), P(k).rvD, e.dep_frac, 'spline');
    P(k).rv0 = rv0;
    P(k).rvf = interp1(tT/tT(end), P(k).rvT, e.arr_frac, 'spline');
    [tj, yj] = pumpkyn.cr3bp.tfMinProp(e.z8(8), [rv0(1:6)'; 1; e.z8(1:7)], ...
                   ndT(e.thrustN), cnd, mu);
    [P(k).t, iu] = unique(tj);
    P(k).y  = yj(iu, 1:3);
    P(k).tf = e.z8(8);
    P(k).e  = e;
end
tEnd = max(P(1).tf, P(2).tf);            % the SHARED clock's last tick

%% Render:
 nFrames = 150;
     fps = 25;
fig = figure('Color','w', 'Visible','off', 'Position',[50 50 1280 720]);
if strncmpi(metric,'t',1), ttl = {'FASTEST','SLOWEST'};
else,                      ttl = {'MINIMUM \DeltaV','MAXIMUM \DeltaV'}; end
for k = 1:2
    ax(k) = subplot(1,2,k);  hold(ax(k),'on');
    plot3(ax(k), P(k).rvD(:,1), P(k).rvD(:,2), P(k).rvD(:,3), 'k--');
    plot3(ax(k), P(k).rvT(:,1), P(k).rvT(:,2), P(k).rvT(:,3), 'r--');
    plot3(ax(k), P(k).rv0(1), P(k).rv0(2), P(k).rv0(3), '.g', 'MarkerSize',18);
    plot3(ax(k), P(k).rvf(1), P(k).rvf(2), P(k).rvf(3), '.r', 'MarkerSize',18);
    hTrail(k) = plot3(ax(k), NaN, NaN, NaN, 'b-', 'LineWidth', 1.5);
    hCraft(k) = plot3(ax(k), NaN, NaN, NaN, 'ob', 'MarkerFaceColor','b', ...
                      'MarkerSize', 7);
    axis(ax(k), 'equal');  grid(ax(k), 'on');  view(ax(k), -35, 25);
    xlabel(ax(k),'x [ND]'); ylabel(ax(k),'y [ND]'); zlabel(ax(k),'z [ND]');
    title(ax(k), {sprintf('%s:  DRO \\tau=%.2f, N_p=%d, %.2f N', ...
              ttl{k}, P(k).e.tauDRO, P(k).e.Np, P(k).e.thrustN), ...
          sprintf('\\DeltaV %.3f km/s,  t_f %.3f d', ...
              P(k).e.deltaV_kms, P(k).e.tf_days)});
end
% the vertical dividing line and the shared clock
annotation(fig, 'line', [0.5 0.5], [0.03 0.93], 'LineWidth', 1.2, ...
           'Color', [0.3 0.3 0.3]);
hClock = annotation(fig, 'textbox', [0.38 0.94 0.24 0.05], 'String', '', ...
    'HorizontalAlignment','center', 'EdgeColor','none', ...
    'FontSize', 14, 'FontWeight','bold');

vw = VideoWriter([outStem '.mp4'], 'MPEG-4');
vw.FrameRate = fps;  open(vw);
gifFile = [outStem '.gif'];
for kf = 1:nFrames
    tNow = tEnd * (kf-1)/(nFrames-1);
    for k = 1:2
        tk = min(tNow, P(k).tf);         % freeze after this panel finishes
        sel = P(k).t <= tk;
        pos = interp1(P(k).t, P(k).y, tk);
        set(hTrail(k), 'XData', [P(k).y(sel,1); pos(1)], ...
                       'YData', [P(k).y(sel,2); pos(2)], ...
                       'ZData', [P(k).y(sel,3); pos(3)]);
        set(hCraft(k), 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3));
    end
    set(hClock, 'String', sprintf('t = %.2f days', tNow*tStar/86400));
    fr = getframe(fig);
    img = fr.cdata;
    if size(img,1) ~= 720                % Retina renders at 2x: downsample
        img = img(1:2:end, 1:2:end, :);
    end
    img = img(1:min(end,720), 1:min(end,1280), :);
    writeVideo(vw, img);
    if mod(kf,2) == 1                    % gif at half rate to keep size sane
        [A, cmap] = rgb2ind(img, 128);
        if kf == 1
            imwrite(A, cmap, gifFile, 'gif', 'LoopCount', inf, ...
                    'DelayTime', 2/fps);
        else
            imwrite(A, cmap, gifFile, 'gif', 'WriteMode','append', ...
                    'DelayTime', 2/fps);
        end
    end
end
close(vw);  close(fig);
fprintf('movies -> %s.mp4 and %s.gif  (%d frames, shared clock to %.2f days)\n', ...
        outStem, outStem, nFrames, tEnd*tStar/86400);
end
