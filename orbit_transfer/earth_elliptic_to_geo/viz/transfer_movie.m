function transfer_movie(matFile, outStem)
% TRANSFER_MOVIE  Control movie for a LEO-ellipse->GEO min-fuel transfer.
%
% Renders the inertial-frame transfer with the control law: trajectory
% colored by throttle state (red = full-thrust burn, blue = coast), primer
% thrust-direction arrow on the spacecraft during burns, a synced bang-bang
% throttle strip, a running Delta-V accumulation curve, and a bold Delta-V
% meter. Playback is uniform in PHYSICAL time (X(8,:)), matching the
% campaign's Earth-Moon transfer movies (PSR/psr_movie.m house style).
%
% INPUTS:
%   matFile - EITHER a path to a run_transfer results .mat holding a `res`
%             variable [char/string], OR the `res` struct itself (e.g. from
%             mee_res_to_cart_res) [struct]. Either way res.cfg (.thrustN
%             .ctf), res.fuel.X (9xnN = [r;v;m;t;cScale], row 9 unused
%             here), res.fuel.U (4xnN = [alpha;s])
%   outStem - output basename WITHOUT extension; writes <outStem>.mp4 and
%             <outStem>.gif [char]
%
% OUTPUTS: none (files written; paths printed)
%
% REFERENCES:
%   [1] GTO_tulip/PSR/psr_movie.m (house-style movie this file
%       is restyled to match: layout, colors, Delta-V meter design).
%   [2] memory: matlab-movie-diagonal-streaks (fixed divisible-by-16 frame ->
%       no H.264 shear).
%
% NOTE: The Delta-V/mass meter assumes campaign Isp=2000 s and m0=1500 kg.
if isstruct(matFile)
    res = matFile;                       % already a Cartesian res struct
else
    S = load(matFile); res = S.res;      % path to a .mat holding `res`
end
p = kepler_lt_params(res.cfg.thrustN, 1500, 2000);
m0kg = p.m0kg;  c = p.c;

X = res.fuel.X;  U = res.fuel.U;
r = X(1:3,:);  m = X(7,:);  t = X(8,:);
s = U(4,:);   al = U(1:3,:);
burnTol = 0.05;                          % throttle above this = BURNING (red).
burn  = s > burnTol;                     % red whenever thrusting to ANY extent
nSw   = sum(abs(diff(burn)));            % on/off transition count (display only)
tDays = t * p.TU_s/86400;
dV    = c*log(1./m) * p.VU_kms;          % running Delta-V (km/s), per node
mKg   = m0kg*m;  propKg = m0kg*(1-m);
dVtot = dV(end);

% --- GEO ring backdrop (target orbit: radius-1 circle, equatorial plane) ----
thG  = linspace(0, 2*pi, 361);
geoRing = [cos(thG); sin(thG); 0*thG].';

% --- fixed axis limits -------------------------------------------------------
allP = [r.'; geoRing];
pad = 0.05*(max(allP)-min(allP)+eps);
xl = [min(allP(:,1))-pad(1), max(allP(:,1))+pad(1)];
yl = [min(allP(:,2))-pad(2), max(allP(:,2))+pad(2)];
zl = [min(allP(:,3))-pad(3), max(allP(:,3))+pad(3)];
if diff(zl) < eps, zl = [-0.1 0.1]; end  % near-coplanar case: avoid degenerate zlim

% --- uniform-TIME frame schedule ---------------------------------------------
tFrames = linspace(t(1), t(end), 300);

% --- figure scaffold ----------------------------------------------------------
fig = figure('Color','w','Position',[100 100 1000 860],'Visible','off');
try, theme(fig,'light'); catch, end
tl = tiledlayout(fig, 5, 1, 'TileSpacing','compact','Padding','compact');

axT = nexttile(tl, 1, [3 1]);   hold(axT,'on'); grid(axT,'on'); box(axT,'on');
plot3(axT, r(1,:), r(2,:), r(3,:), '-', 'Color',[0.74 0.74 0.78], 'LineWidth',0.5);
plot3(axT, geoRing(:,1), geoRing(:,2), geoRing(:,3), '-', 'Color',[0.45 0.72 0.45], 'LineWidth',1.0);
plot3(axT, 0,0,0,'o','MarkerFaceColor',[0.10 0.35 0.80],'MarkerEdgeColor','k','MarkerSize',11);
text(axT, 0,0,0.05,'Earth','FontSize',9);
xlim(axT,xl); ylim(axT,yl); zlim(axT,zl);
xlabel(axT,'x (ND)'); ylabel(axT,'y'); zlabel(axT,'z');
view(axT,-37,24); daspect(axT,[1 1 1]);
title(axT, sprintf('LEO ellipse \\rightarrow GEO min-fuel  (T=%g N, c_{tf}=%.2f)', ...
      res.cfg.thrustN, res.cfg.ctf));

hBurn = plot3(axT, nan,nan,nan,'-','Color',[0.85 0.15 0.15],'LineWidth',2.2);
hCoast= plot3(axT, nan,nan,nan,'-','Color',[0.15 0.35 0.85],'LineWidth',2.0);
hSC   = plot3(axT, nan,nan,nan,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6);
hArr  = quiver3(axT, nan,nan,nan,nan,nan,nan,0,'Color',[0.85 0.15 0.15],'LineWidth',1.6,'MaxHeadSize',2);
arrScale = 0.06*diff(xl);

% BOLD running Delta-V meter, axes-normalized so it never sits on the spiral
hDV = text(axT, 0.605, 0.42, '', 'Units','normalized', ...
    'FontSize',16,'FontWeight','bold','FontName','Menlo','Color',[0.60 0.10 0.10], ...
    'BackgroundColor',[1 1 1 0.80],'EdgeColor',[0.6 0.1 0.1],'Margin',5, ...
    'HorizontalAlignment','left','VerticalAlignment','middle');
hTxt = text(axT, 0.605, 0.24, '', 'Units','normalized', ...
    'FontSize',10,'FontName','Menlo','BackgroundColor',[1 1 1 0.6], ...
    'HorizontalAlignment','left','VerticalAlignment','top');

% throttle strip (tile 4) — the control law, synced to the animation
axS = nexttile(tl, 4);  hold(axS,'on'); grid(axS,'on'); box(axS,'on');
stairs(axS, tDays, s, '-','Color',[0.45 0.45 0.45],'LineWidth',0.9);
ylim(axS,[-0.05 1.08]); xlim(axS,[0 tDays(end)]);
ylabel(axS,'throttle s');
if mean(s>0.05 & s<0.95) > 0.05      % smooth control (many interior nodes)
    title(axS, sprintf('throttle: smooth control, %d on/off transitions (red burning / blue coasting)', nSw));
else
    title(axS, sprintf('throttle: %d-switch bang-bang (red burn / blue coast)', nSw));
end
set(axS,'XTickLabel',[]);
hCurS = plot(axS,[0 0],[-0.05 1.08],'k-','LineWidth',1.0);
hDotS = plot(axS, nan,nan,'o','MarkerFaceColor',[0.85 0.15 0.15],'MarkerEdgeColor','k','MarkerSize',6);

% Delta-V accumulation curve (tile 5)
axV = nexttile(tl, 5);  hold(axV,'on'); grid(axV,'on'); box(axV,'on');
plot(axV, tDays, dV, '-','Color',[0.80 0.80 0.82],'LineWidth',1.2);
ylim(axV,[0 1.06*dVtot]); xlim(axV,[0 tDays(end)]);
xlabel(axV,'time (days)'); ylabel(axV,'\DeltaV (km/s)');
title(axV, sprintf('running \\DeltaV  (total %.4f km/s,  prop %.4f kg,  %d switches)', ...
      dVtot, propKg(end), nSw));
hVfill = plot(axV, nan,nan,'-','Color',[0.60 0.10 0.10],'LineWidth',2.4);
hCurV  = plot(axV,[0 0],[0 1.06*dVtot],'k-','LineWidth',1.0);
hDotV  = plot(axV, nan,nan,'o','MarkerFaceColor',[0.60 0.10 0.10],'MarkerEdgeColor','k','MarkerSize',7);

% --- render --------------------------------------------------------------------
vw = VideoWriter(outStem,'MPEG-4'); vw.FrameRate = 24; vw.Quality = 95; open(vw);
gifFile = [outStem '.gif']; gifStride = 2; gifW = 640; gifDelay = 1/12; gifMap = [];
tmpPng = [outStem '_frame_tmp.png'];
vidHW = [720 1280];                                  % divisible by 16: no H.264 shear
for fc = 1:numel(tFrames)
    tF = tFrames(fc);
    f = find(t <= tF, 1, 'last'); if isempty(f), f = 1; end
    bMask = burn(1:f);  cMask = ~bMask;
    xb = r(1,1:f); yb = r(2,1:f); zb = r(3,1:f);
    set(hBurn, 'XData',maskv(xb,bMask),'YData',maskv(yb,bMask),'ZData',maskv(zb,bMask));
    set(hCoast,'XData',maskv(xb,cMask),'YData',maskv(yb,cMask),'ZData',maskv(zb,cMask));
    set(hSC,'XData',r(1,f),'YData',r(2,f),'ZData',r(3,f));
    if burn(f)
        % arrow length scales with throttle magnitude: full at s=1, half at s=0.5
        a = al(:,f); an = a/max(norm(a),eps)*arrScale*s(f);
        set(hArr,'XData',r(1,f),'YData',r(2,f),'ZData',r(3,f), ...
                 'UData',an(1),'VData',an(2),'WData',an(3),'Visible','on');
        state = sprintf('BURN  s=%.2f', s(f));
    else
        set(hArr,'Visible','off'); state = 'COAST';
    end
    set(hDV,'String',sprintf(' \\DeltaV = %6.3f km/s ', dV(f)));
    set(hTxt,'String',sprintf(' t = %5.2f d  (%.0f%%)\n mass = %6.3f kg\n prop = %5.3f kg\n %s', ...
        tDays(f), 100*t(f)/t(end), mKg(f), propKg(f), state));
    set(hCurS,'XData',[tDays(f) tDays(f)]);  set(hDotS,'XData',tDays(f),'YData',s(f));
    set(hVfill,'XData',tDays(1:f),'YData',dV(1:f));
    set(hCurV,'XData',[tDays(f) tDays(f)]);  set(hDotV,'XData',tDays(f),'YData',dV(f));

    exportgraphics(fig, tmpPng, 'Resolution',110);
    img = imread(tmpPng);
    if ~isequal([size(img,1) size(img,2)], vidHW), img = imresize(img, vidHW); end
    writeVideo(vw, img);
    if mod(fc-1,gifStride)==0
        gr = round(gifW*size(img,1)/size(img,2));  gimg = imresize(img,[gr gifW]);
        if isempty(gifMap)
            [gInd,gifMap] = rgb2ind(gimg,256,'nodither');
            imwrite(gInd,gifMap,gifFile,'gif','LoopCount',Inf,'DelayTime',gifDelay);
        else
            gInd = rgb2ind(gimg,gifMap,'nodither');
            imwrite(gInd,gifMap,gifFile,'gif','WriteMode','append','DelayTime',gifDelay);
        end
    end
end
close(vw);  if isfile(tmpPng), delete(tmpPng); end
close(fig);
fprintf('WROTE %s.mp4 (%d frames)\n', outStem, numel(tFrames));
fprintf('WROTE %s.gif (%d frames)\n', outStem, numel(1:gifStride:numel(tFrames)));
end

function v = maskv(x,m)
% MASKV  NaN-mask a vector so plotted segments break where the mask is false.
v = x; v(~m) = nan;
end
