function psr_movie(solFile, outStem, titleStr, mode, bgTrace)
% PSR_MOVIE  Control movie for a PSR (or any Sundman min-fuel) solution.
%
% Renders the GTO->tulip transfer in the rotating CR3BP frame with the control
% law: trajectory colored by throttle state (red = full-thrust burn, blue =
% coast), primer thrust-direction arrow on the spacecraft during burns, a
% synced bang-bang throttle strip, a running Delta-V accumulation curve, and a
% bold Delta-V meter. Playback is uniform in PHYSICAL time (the Sundman mesh
% is non-uniform in time by construction, so frames are scheduled on the
% carried time state X(8,:), not on node index).
%
% Generalized from movie/animate_sundman_minfuel.m (which hardcodes the
% certified 1.15x file) to take any solution .mat in the campaign seed layout.
%
% INPUTS:
%   solFile  - solution .mat in seed layout: out (with X [8xnN] = [r;v;m;t],
%              U [4xnN] = [alpha;s]), sigma, tauf0, rv0, rvf, factor
%   outStem  - output basename WITHOUT extension; writes <outStem>.mp4 and
%              <outStem>.gif ('movie') or <outStem>_{early,mid,late}.png
%              ('preview') [char]
%   titleStr - trajectory-panel title [char]
%   mode     - 'preview' (3 stills, seconds) | 'movie' (MP4+GIF, ~10-15 min
%              at 320 frames) [char, default 'preview']
%   bgTrace  - (optional) [Kx3+] backdrop orbit trace to draw instead of the
%              default tulip (e.g. the ELFO one-period trace). Omit for tulip.
%
% Accepts BOTH the tulip seed layout (S.out.X 8xnN) and the free-t_f ELFO layout
% (top-level S.X 9xnN); only rows 1-8 [r;v;m;t] are used, so a 9th cScale row is
% ignored. All panels (trajectory, throttle, running dV via the rocket equation)
% are model-independent, so this renders either transfer correctly.
%
% OUTPUTS: none (files written; paths printed)
%
% REFERENCES:
%   [1] movie/animate_sundman_minfuel.m (layout + Delta-V meter design).
%   [2] PSR/run_psr.m section 5 (pipeline caller).

if nargin < 4 || isempty(mode), mode = 'preview'; end

% --- constants (campaign values; match minfuel_config / cr3bp_lt_params) ----
p = cr3bp_lt_params(0.025, 15, 2100);
muStar = p.muStar;  lStar = p.lStar;  tStar = p.tStar;
m0kg = p.m0kg;  c = p.c;
earth = [-muStar 0 0];  moon = [1-muStar 0 0];

% --- load solution ----------------------------------------------------------
S  = load(solFile);
if isfield(S,'out') && isfield(S.out,'X')      % tulip seed layout
    X = S.out.X;  U = S.out.U;
else                                           % free-t_f ELFO layout (top-level)
    X = S.X;      U = S.U;
end
rvf = S.rvf;
r  = X(1:3,:);  m = X(7,:);  t = X(8,:);
s  = U(4,:);   al = U(1:3,:);
burnTol = 0.05;                          % throttle above this = BURNING (red).
burn  = s > burnTol;                     % red whenever thrusting to ANY extent
nSw   = sum(abs(diff(burn)));            % on/off transition count (display only)
tDays = t * tStar/86400;
dV    = c*log(1./m) * lStar/tStar;       % running Delta-V (km/s), per node
mKg   = m0kg*m;  propKg = m0kg*(1-m);
dVtot = dV(end);

% --- backdrop orbit trace: caller-supplied (e.g. ELFO) or tulip default -----
yTul = [];
if nargin >= 5 && ~isempty(bgTrace)
    yTul = bgTrace(:, 1:3);                    % caller's target orbit (e.g. ELFO)
else
    try                                        % default tulip (best-effort)
        [~, x0T] = pumpkyn.cr3bp.getTulip((5/6)*2*pi, 7, -1, 1e-12);
        [~, yT ] = pumpkyn.cr3bp.prop((5/6)*2*pi, x0T, muStar);
        yTul = yT(:,1:3);
    catch
    end
end

% --- fixed axis limits -------------------------------------------------------
allP = r.';  if ~isempty(yTul), allP = [allP; yTul]; end
pad = 0.05*(max(allP)-min(allP)+eps);
xl = [min(allP(:,1))-pad(1), max(allP(:,1))+pad(1)];
yl = [min(allP(:,2))-pad(2), max(allP(:,2))+pad(2)];
zl = [min(allP(:,3))-pad(3), max(allP(:,3))+pad(3)];

% --- uniform-TIME frame schedule ---------------------------------------------
if strcmp(mode,'preview')
    tFrames = t(1) + [0.30 0.60 0.92]*(t(end)-t(1));
    tags = {'early','mid','late'};
else
    tFrames = linspace(t(1), t(end), 320);
end

% --- figure scaffold ----------------------------------------------------------
fig = figure('Color','w','Position',[100 100 1000 860],'Visible','off');
try, theme(fig,'light'); catch, end
tl = tiledlayout(fig, 5, 1, 'TileSpacing','compact','Padding','compact');

axT = nexttile(tl, 1, [3 1]);   hold(axT,'on'); grid(axT,'on'); box(axT,'on');
plot3(axT, r(1,:), r(2,:), r(3,:), '-', 'Color',[0.74 0.74 0.78], 'LineWidth',0.5);
if ~isempty(yTul)
    plot3(axT, yTul(:,1), yTul(:,2), yTul(:,3), '-', 'Color',[0.45 0.72 0.45], 'LineWidth',1.0);
end
plot3(axT, earth(1),earth(2),earth(3),'o','MarkerFaceColor',[0.10 0.35 0.80],'MarkerEdgeColor','k','MarkerSize',11);
plot3(axT, moon(1),moon(2),moon(3),'o','MarkerFaceColor',[0.60 0.60 0.60],'MarkerEdgeColor','k','MarkerSize',8);
plot3(axT, rvf(1),rvf(2),rvf(3),'p','MarkerFaceColor',[0.20 0.65 0.20],'MarkerEdgeColor','k','MarkerSize',13);
text(axT, earth(1),earth(2),earth(3)+0.05,'Earth','FontSize',9);
text(axT, moon(1),moon(2),moon(3)+0.05,'Moon','FontSize',9);
xlim(axT,xl); ylim(axT,yl); zlim(axT,zl);
xlabel(axT,'x (rot, ND)'); ylabel(axT,'y'); zlabel(axT,'z');
view(axT,-37,24); daspect(axT,[1 1 1]);
title(axT, titleStr);

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
if ~strcmp(mode,'preview')
    vw = VideoWriter(outStem,'MPEG-4'); vw.FrameRate = 24; vw.Quality = 95; open(vw);
    gifFile = [outStem '.gif']; gifStride = 2; gifW = 640; gifDelay = 1/12; gifMap = [];
    tmpPng = [outStem '_frame_tmp.png'];
end
fc = 0; vidHW = [];
for tF = tFrames
    fc = fc + 1;
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

    if strcmp(mode,'preview')
        exportgraphics(fig, sprintf('%s_%s.png', outStem, tags{fc}), 'Resolution',130);
    else
        exportgraphics(fig, tmpPng, 'Resolution',110);
        img = imread(tmpPng);
        if isempty(vidHW), vidHW = 2*floor([size(img,1) size(img,2)]/2); end
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
end
if ~strcmp(mode,'preview')
    close(vw);
    fprintf('WROTE %s.mp4 (%d frames)\n', outStem, numel(tFrames));
    fprintf('WROTE %s.gif (%d frames)\n', outStem, numel(1:gifStride:numel(tFrames)));
    if isfile(tmpPng), delete(tmpPng); end
else
    fprintf('WROTE %s_{early,mid,late}.png\n', outStem);
end
close(fig);
end

function v = maskv(x,m)
% MASKV  NaN-mask a vector so plotted segments break where the mask is false.
v = x; v(~m) = nan;
end
