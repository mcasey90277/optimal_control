function animate_minfuel(mode)
% ANIMATE_MINFUEL  Render the min-fuel arrival-leg solution as a movie.
%
% Rotating CR3BP frame (Earth, Moon, tulip fixed). The min-time GTO->tulip
% spiral is drawn faint for context; the verified min-fuel NLP leg is
% animated in bold, colored by throttle state (red = burn, blue = coast),
% with primer-direction thrust arrows during the burn and a synced
% throttle subplot + mass-fraction readout.
%
% INPUTS:
%   mode - 'preview' saves three still PNGs (early / switch / late) and
%          exits; 'movie' renders BOTH the MP4 and a smaller looping GIF.
%          [default 'movie']

if nargin < 1, mode = 'movie'; end

sp = [fileparts(mfilename('fullpath')) filesep];   % this movie/ folder
S  = load([sp 'minfuel_movie_data.mat']);

muStar = S.muStar;  tStar = S.tStar;
earth  = [-muStar, 0, 0];
moon   = [1 - muStar, 0, 0];

% --- leg solution ----------------------------------------------------------
X  = S.out.X;                 % 7 x (N+1)
U  = S.out.U;                 % 4 x (N+1)  [w; s]
tM = S.out.tauMesh(:).';      % 1 x (N+1), leg-local ND time
rL = X(1:3, :);
mL = X(7, :);
sL = U(4, :);                 % throttle
wL = U(1:3, :);               % thrust direction (|w| = s)
burn = sL > 0.5;
tDaysLeg = tM * tStar/86400;
tfLeg    = tM(end);

% switching function proxy from the throttle plot (s already encodes it)
% mass in kg
mKg = S.m0kg * mL;
mProp = S.m0kg * (mL(1) - mL);

% --- backdrop spiral -------------------------------------------------------
rSpiral = S.yF(:, 1:3);
iLeg    = S.iLeg;

% --- fixed axis limits (whole transfer) ------------------------------------
allP = [rSpiral; rL.'];
pad  = 0.05*(max(allP) - min(allP) + eps);
xl = [min(allP(:,1))-pad(1), max(allP(:,1))+pad(1)];
yl = [min(allP(:,2))-pad(2), max(allP(:,2))+pad(2)];
zl = [min(allP(:,3))-pad(3), max(allP(:,3))+pad(3)];

% --- frame schedule --------------------------------------------------------
nNode = numel(tM);
if strcmp(mode, 'preview')
    iSw = find(diff(burn) ~= 0, 1, 'first'); if isempty(iSw), iSw = round(nNode/2); end
    frames = [round(0.25*iSw), iSw, round((iSw+nNode)/2)];
    tags   = {'early','switch','late'};
else
    frames = unique(round(linspace(2, nNode, 260)));
end

% --- figure scaffold -------------------------------------------------------
fig = figure('Color','w','Position',[100 100 1000 760],'Visible','off');
tl  = tiledlayout(fig, 4, 1, 'TileSpacing','compact','Padding','compact');

axT = nexttile(tl, 1, [3 1]);   % trajectory (3D)
hold(axT,'on'); grid(axT,'on'); box(axT,'on');
% backdrop
plot3(axT, rSpiral(:,1), rSpiral(:,2), rSpiral(:,3), '-', ...
      'Color',[0.72 0.72 0.75], 'LineWidth',0.6);
if ~isempty(S.yTul)
    plot3(axT, S.yTul(:,1), S.yTul(:,2), S.yTul(:,3), '-', ...
          'Color',[0.55 0.75 0.55], 'LineWidth',0.8);
end
% Earth + Moon
plot3(axT, earth(1),earth(2),earth(3),'o','MarkerFaceColor',[0.1 0.35 0.8], ...
      'MarkerEdgeColor','k','MarkerSize',11);
plot3(axT, moon(1),moon(2),moon(3),'o','MarkerFaceColor',[0.6 0.6 0.6], ...
      'MarkerEdgeColor','k','MarkerSize',8);
text(axT, earth(1),earth(2),earth(3)+0.04,'Earth','FontSize',9);
text(axT, moon(1),moon(2),moon(3)+0.04,'Moon','FontSize',9);
xlim(axT,xl); ylim(axT,yl); zlim(axT,zl);
xlabel(axT,'x (rot, ND)'); ylabel(axT,'y'); zlabel(axT,'z');
view(axT, -37, 22); daspect(axT,[1 1 1]);
title(axT,'Minimum-fuel arrival leg (direct-NLP solution), rotating CR3BP frame');

% animated graphics handles
hBurn = plot3(axT, nan, nan, nan, '-', 'Color',[0.85 0.15 0.15], 'LineWidth',2.2);
hCoast= plot3(axT, nan, nan, nan, '-', 'Color',[0.15 0.35 0.85], 'LineWidth',2.2);
hSC   = plot3(axT, nan, nan, nan, 'o', 'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6);
hArr  = quiver3(axT, nan,nan,nan, nan,nan,nan, 0, 'Color',[0.85 0.15 0.15],'LineWidth',1.6,'MaxHeadSize',2);
hTxt  = text(axT, xl(1)+0.03*diff(xl), yl(2)-0.05*diff(yl), zl(2), '', ...
             'FontSize',10,'FontName','Menlo','BackgroundColor',[1 1 1 0.6],'VerticalAlignment','top');

axS = nexttile(tl, 4);          % throttle subplot (row 4)
hold(axS,'on'); grid(axS,'on'); box(axS,'on');
stairs(axS, tDaysLeg, sL, '-', 'Color',[0.4 0.4 0.4], 'LineWidth',1.0);
ylim(axS,[-0.05 1.08]); xlim(axS,[0 tDaysLeg(end)]);
xlabel(axS,'leg time (days)'); ylabel(axS,'throttle s');
title(axS,'throttle: bang-bang (burn \rightarrow coast)');
hCur = plot(axS, [0 0], [-0.05 1.08], 'k-', 'LineWidth',1.0);
hDot = plot(axS, nan, nan, 'o','MarkerFaceColor',[0.85 0.15 0.15],'MarkerEdgeColor','k','MarkerSize',6);

arrScale = 0.06*diff(xl);   % thrust-arrow length in ND

% --- render ----------------------------------------------------------------
if ~strcmp(mode,'preview')
    vw = VideoWriter([sp 'minfuel_solution'], 'MPEG-4');
    vw.FrameRate = 24; vw.Quality = 95;
    open(vw);
    % companion GIF: downscaled, every gifStride-th frame, shared palette
    gifFile   = [sp 'minfuel_solution.gif'];
    gifStride = 2;        % keep ~half the frames
    gifW      = 640;      % target width (px)
    gifDelay  = 1/12;     % s per GIF frame (~12 fps after striding)
    gifMap    = [];
end

fc = 0;
vidHW = [];      % locked frame size (H W), set on first movie frame
for f = frames
    fc = fc + 1;
    isB = burn(1:f);
    % burn / coast segments up to node f (NaN-break at state changes)
    xb = rL(1,1:f); yb = rL(2,1:f); zb = rL(3,1:f);
    bMask = isB;  cMask = ~isB;
    set(hBurn, 'XData', maskv(xb,bMask), 'YData', maskv(yb,bMask), 'ZData', maskv(zb,bMask));
    set(hCoast,'XData', maskv(xb,cMask), 'YData', maskv(yb,cMask), 'ZData', maskv(zb,cMask));
    set(hSC, 'XData', rL(1,f), 'YData', rL(2,f), 'ZData', rL(3,f));

    if burn(f)
        a = wL(:,f); an = a/max(norm(a),eps)*arrScale;
        set(hArr, 'XData',rL(1,f),'YData',rL(2,f),'ZData',rL(3,f), ...
                  'UData',an(1),'VData',an(2),'WData',an(3),'Visible','on');
        state = 'BURN';
    else
        set(hArr,'Visible','off'); state = 'COAST';
    end
    set(hTxt,'String',sprintf(' t = %5.2f d  (%.0f%% of leg)\n mass  = %6.3f kg\n prop  = %5.3f kg\n %s', ...
        tDaysLeg(f), 100*tM(f)/tfLeg, mKg(f), mProp(f), state));

    set(hCur,'XData',[tDaysLeg(f) tDaysLeg(f)]);
    set(hDot,'XData',tDaysLeg(f),'YData',sL(f));

    if strcmp(mode,'preview')
        exportgraphics(fig, sprintf('%spreview_%s.png', sp, tags{fc}), 'Resolution',130);
    else
        % headless-safe: rasterize each frame deterministically, then assemble
        tmpPng = [sp 'frame_tmp.png'];
        exportgraphics(fig, tmpPng, 'Resolution',110);
        img = imread(tmpPng);
        if isempty(vidHW)
            vidHW = 2*floor([size(img,1) size(img,2)]/2);  % even dims for H.264
        end
        if ~isequal([size(img,1) size(img,2)], vidHW)
            img = imresize(img, vidHW);   % lock every frame to the first size
        end
        writeVideo(vw, img);

        % --- companion GIF frame (strided, downscaled, shared palette) -----
        if mod(fc-1, gifStride) == 0
            gr   = round(gifW*size(img,1)/size(img,2));
            gimg = imresize(img, [gr gifW]);
            if isempty(gifMap)
                [gInd, gifMap] = rgb2ind(gimg, 256, 'nodither');
                imwrite(gInd, gifMap, gifFile, 'gif', ...
                        'LoopCount', Inf, 'DelayTime', gifDelay);
            else
                gInd = rgb2ind(gimg, gifMap, 'nodither');
                imwrite(gInd, gifMap, gifFile, 'gif', ...
                        'WriteMode', 'append', 'DelayTime', gifDelay);
            end
        end
    end
end

if ~strcmp(mode,'preview')
    close(vw);
    nGif = numel(1:gifStride:numel(frames));
    fprintf('WROTE %sminfuel_solution.mp4 (%d frames)\n', sp, numel(frames));
    fprintf('WROTE %sminfuel_solution.gif (%d frames)\n', sp, nGif);
else
    fprintf('WROTE preview PNGs\n');
end
close(fig);
end

% ---------------------------------------------------------------------------
function v = maskv(x, m)
% keep entries where m is true, NaN elsewhere (breaks the polyline cleanly)
v = x; v(~m) = nan;
end
