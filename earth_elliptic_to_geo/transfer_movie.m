function transfer_movie(matFile, outStem)
% TRANSFER_MOVIE  Animate a LEO-ellipse->GEO min-fuel transfer (MP4 + GIF).
%
% Top panel: inertial 3D trajectory colored by throttle (red burn / blue coast),
% GEO ring + target star; bottom: throttle strip with time cursor; text meter
% with running mass/DeltaV. Frames uniform in PHYSICAL time (X(8,:)).
%
% INPUTS:  matFile - run_transfer results .mat;  outStem - output basename
% OUTPUTS: none (writes <outStem>.mp4 and <outStem>.gif)
%
% REFERENCES: [1] PSR/psr_movie.m (layout). [2] memory: matlab-movie-diagonal-
%   streaks (fixed 1280x720 divisible-by-16 frame -> no H.264 shear).
S = load(matFile);  res = S.res;
p  = kepler_lt_params(res.cfg.thrustN, 1500, 2000);
X = res.fuel.X;  U = res.fuel.U;
r = X(1:3,:);  t = X(8,:);  m = X(7,:);  ss = U(4,:);
tD = t * p.TU_s/86400;  burn = ss > 0.5;
dV = p.c*log(1./m)*p.VU_kms;
th = linspace(0, 2*pi, 361);
fig = figure('Color','w','Position',[80 80 1000 750],'Visible','off');
axT = subplot('Position',[0.06 0.32 0.90 0.62]);  hold(axT,'on'); grid(axT,'on');
plot3(axT, r(1,:), r(2,:), r(3,:), '-', 'Color',[0.8 0.8 0.83], 'LineWidth',0.5);
plot3(axT, cos(th), sin(th), 0*th, 'g-', 'LineWidth', 1.0);
plot3(axT, 0,0,0, 'o', 'MarkerFaceColor',[0.1 0.35 0.8], 'MarkerSize',10);
hB = plot3(axT, nan,nan,nan, 'r-', 'LineWidth',1.8);
hC = plot3(axT, nan,nan,nan, 'b-', 'LineWidth',1.5);
hS = plot3(axT, nan,nan,nan, 'ko', 'MarkerFaceColor','k', 'MarkerSize',5);
hTx = text(axT, 0.02, 0.95, '', 'Units','normalized', 'FontName','Menlo', 'FontSize',10);
axis(axT, 'equal');  view(axT, -30, 25);
title(axT, sprintf('LEO ellipse %s GEO min-fuel  (T=%g N, c_{tf}=%.2f)', ...
      char(8594), res.cfg.thrustN, res.cfg.ctf));
axS = subplot('Position',[0.06 0.06 0.90 0.18]);  hold(axS,'on'); grid(axS,'on');
stairs(axS, tD, ss, '-', 'Color',[0.4 0.4 0.4]);  ylim(axS, [-0.05 1.08]);
xlabel(axS, 'time [days]');  ylabel(axS, 'throttle');
hCur = plot(axS, [0 0], [-0.05 1.08], 'k-');
vw = VideoWriter(outStem, 'MPEG-4');  vw.FrameRate = 24;  vw.Quality = 95;  open(vw);
gifFile = [outStem '.gif'];  gifMap = [];  tmp = [outStem '_tmp.png'];
vidHW = [720 1280];                                  % divisible by 16: no H.264 shear
tFr = linspace(t(1), t(end), 300);
for fc = 1:numel(tFr)
    k = find(t <= tFr(fc), 1, 'last');
    mask = @(vv, mm) subsasgn(vv, substruct('()', {~mm}), nan);
    xb = r(1,1:k);  yb = r(2,1:k);  zb = r(3,1:k);  bm = burn(1:k);
    set(hB, 'XData',mask(xb,bm),  'YData',mask(yb,bm),  'ZData',mask(zb,bm));
    set(hC, 'XData',mask(xb,~bm), 'YData',mask(yb,~bm), 'ZData',mask(zb,~bm));
    set(hS, 'XData',r(1,k), 'YData',r(2,k), 'ZData',r(3,k));
    set(hTx,'String', sprintf('t=%5.1f d  m=%7.1f kg  dV=%5.3f km/s', ...
        tD(k), 1500*m(k), dV(k)));
    set(hCur, 'XData', [tD(k) tD(k)]);
    drawnow;
    exportgraphics(fig, tmp, 'Resolution', 120);
    img = imresize(imread(tmp), vidHW);
    writeVideo(vw, img);
    if mod(fc-1, 2) == 0
        gi = imresize(img, [360 640]);
        if isempty(gifMap)
            [gInd, gifMap] = rgb2ind(gi, 256, 'nodither');
            imwrite(gInd, gifMap, gifFile, 'gif', 'LoopCount', Inf, 'DelayTime', 1/12);
        else
            gInd = rgb2ind(gi, gifMap, 'nodither');
            imwrite(gInd, gifMap, gifFile, 'gif', 'WriteMode','append', 'DelayTime', 1/12);
        end
    end
end
close(vw);  if isfile(tmp), delete(tmp); end
close(fig);
fprintf('WROTE %s.mp4 / .gif (%d frames)\n', outStem, numel(tFr));
end
