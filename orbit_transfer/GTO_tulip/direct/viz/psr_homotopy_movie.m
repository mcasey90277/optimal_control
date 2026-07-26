function psr_homotopy_movie(framesFile, outStem, opts)
% PSR_HOMOTOPY_MOVIE  Animate the throttle sharpening as epsilon: 1 -> 0.
%
% Renders the energy->fuel homotopy captured by PSR_HOMOTOPY_FRAMES as a
% dark-theme movie: the throttle profile s(t) morphs from the smooth min-energy
% ramp (eps=1) into the bang-bang min-fuel control (eps=0), with a big easing
% epsilon readout, live metrics (dV, propellant, switches, bang-bang fraction),
% and a "bang-bang fraction vs epsilon" progress strip. The captured epsilon
% steps are the KEYFRAMES; in-between frames linearly tween the throttle profile
% (on a common physical-time grid) and epsilon, so the morph is smooth -- a
% visual interpolation between two genuinely-solved profiles, with epsilon
% labeled continuously. Reads only the cache, so re-rendering never re-solves.
%
% INPUTS:
%   framesFile - .mat from psr_homotopy_frames (vars `frames`, `meta`) [char]
%   outStem    - output basename WITHOUT extension; writes <outStem>.mp4 and
%                <outStem>.gif [char]
%   opts       - (optional) struct:
%                nBetween  - tween frames per keyframe segment [scalar, def 14]
%                holdEnd   - extra frames held on the eps=0 solution [def 30]
%                fps       - MP4 frame rate [def 24]
%                nTimeGrid - throttle-curve resolution [def 700]
%
% OUTPUTS: none (files written; paths printed)
%
% REFERENCES:
%   [1] PSR/psr_homotopy_frames.m (the capture that feeds this).
%   [2] ~/.claude/skills/matlab-polished-graphics (polish patterns applied).

if nargin < 3, opts = struct(); end
d = @(f,v) local_default(opts, f, v);
nBetween  = d('nBetween', 14);
holdEnd   = d('holdEnd', 30);
fps       = d('fps', 24);
Ng        = d('nTimeGrid', 700);

S = load(framesFile);  frames = S.frames;  meta = S.meta;
K = numel(frames);
assert(K >= 2, 'need >= 2 keyframes to animate');

tfDays = meta.tf_days;
tg = linspace(0, tfDays, Ng);              % common physical-time grid (days)

% --- resample every keyframe throttle onto the common time grid --------------
% (the Sundman mesh clusters nodes near perigee; interp onto uniform time keeps
% the burn/coast structure and makes keyframes blendable.)
Sgrid = zeros(K, Ng);
epsK  = zeros(1, K);   edgeK = zeros(1,K);  dVK = zeros(1,K);
propK = zeros(1, K);   swK   = zeros(1,K);  defK = zeros(1,K);
for k = 1:K
    tt = frames(k).tDays;  ss = frames(k).s;
    [tt, iu] = unique(tt, 'stable');  ss = ss(iu);
    Sgrid(k,:) = min(max(interp1(tt, ss, tg, 'linear', 'extrap'), 0), 1);
    epsK(k)  = frames(k).eps;   edgeK(k) = 100*frames(k).edge;
    dVK(k)   = frames(k).dV;    propK(k) = frames(k).prop_kg;
    swK(k)   = frames(k).switches;  defK(k) = frames(k).defect;
end

% --- build the tweened frame timeline ----------------------------------------
% Per segment k->k+1: nBetween frames, smoothstep ease so the morph lingers
% slightly at each solved keyframe.
segW = @(u) u.*u.*(3-2*u);                  % smoothstep
Sfr = {};  epsFr = [];  edgeFr = []; dVFr = []; propFr = []; swFr = []; defFr = []; keyFr = [];
for k = 1:K-1
    for j = 1:nBetween
        u = segW((j-1)/nBetween);
        Sfr{end+1}  = (1-u)*Sgrid(k,:) + u*Sgrid(k+1,:); %#ok<AGROW>
        epsFr(end+1)  = (1-u)*epsK(k)  + u*epsK(k+1);   %#ok<AGROW>
        edgeFr(end+1) = (1-u)*edgeK(k) + u*edgeK(k+1);  %#ok<AGROW>
        dVFr(end+1)   = (1-u)*dVK(k)   + u*dVK(k+1);    %#ok<AGROW>
        propFr(end+1) = (1-u)*propK(k) + u*propK(k+1);  %#ok<AGROW>
        % switches is an integer structural count: step to the nearer keyframe.
        swFr(end+1)   = swK(k + (u >= 0.5));            %#ok<AGROW>
        defFr(end+1)  = 10^((1-u)*log10(defK(k)) + u*log10(defK(k+1))); %#ok<AGROW>
        keyFr(end+1)  = (j==1);                         %#ok<AGROW> flash at a solved keyframe
    end
end
% land exactly on eps=0 and hold
for h = 1:holdEnd
    Sfr{end+1}=Sgrid(K,:); epsFr(end+1)=epsK(K); edgeFr(end+1)=edgeK(K);
    dVFr(end+1)=dVK(K); propFr(end+1)=propK(K); swFr(end+1)=swK(K);
    defFr(end+1)=defK(K); keyFr(end+1)=(h==1);           %#ok<AGROW>
end
nFr = numel(Sfr);

% --- colors ------------------------------------------------------------------
bg    = [0.06 0.07 0.10];  fgTxt = [0.90 0.92 0.95];  dim = [0.55 0.60 0.68];
% epsilon color: min-ENERGY (eps=1) cool blue  ->  min-FUEL (eps=0) hot red.
epsCol = @(e) (1-e)*[0.98 0.35 0.22] + e*[0.30 0.62 0.98];
burnCol = [0.98 0.42 0.28];  fillLo = [0.98 0.42 0.28];
% Headless MATLAB (-batch) uses SOFTWARE OpenGL, which renders FaceAlpha/RGBA
% transparency with diagonal banding artifacts. So use NO alpha anywhere -- fake
% the translucent fill by pre-blending the fill color into the background.
blendbg = @(c,a) a*c + (1-a)*bg;

% --- figure scaffold (1080p, dark) -------------------------------------------
fig = figure('Color', bg, 'Position', [80 80 1280 720], 'Visible', 'off', ...
             'InvertHardcopy', 'off');

% main throttle panel
axT = axes(fig, 'Position', [0.075 0.365 0.885 0.500], 'Color', bg);
hold(axT,'on'); box(axT,'on');
set(axT, 'XColor', dim, 'YColor', dim, 'GridColor', [0.30 0.34 0.40], ...
    'GridAlpha', 0.4, 'FontName', 'Helvetica', 'FontSize', 12, 'LineWidth', 1.0);
grid(axT,'on');
xlim(axT, [0 tfDays]);  ylim(axT, [-0.04 1.10]);
xlabel(axT, 'transfer time (days)', 'Color', fgTxt, 'FontSize', 13);
ylabel(axT, 'throttle  s', 'Color', fgTxt, 'FontSize', 14);
yline(axT, 1, ':', 'Color', [0.42 0.47 0.55], 'LineWidth', 1.0);
yline(axT, 0, ':', 'Color', [0.42 0.47 0.55], 'LineWidth', 1.0);
% ghost of the min-energy seed (eps=1) as a static reference under everything
% (solid dim blue -- no RGBA alpha, see blendbg note above)
plot(axT, tg, Sgrid(1,:), '-', 'Color', [0.24 0.38 0.55], 'LineWidth', 1.2);

hFill = fill(axT, [tg(1) tg tg(end)], [0 Sgrid(1,:) 0], blendbg(fillLo,0.30), ...
             'FaceAlpha', 1, 'EdgeColor', 'none');
hLine = plot(axT, tg, Sgrid(1,:), '-', 'Color', burnCol, 'LineWidth', 2.4);

% narrative title + subtitle
tt1 = title(axT, 'From smooth to bang-bang: the min-fuel throttle emerging', ...
    'Color', fgTxt, 'FontSize', 16, 'FontWeight', 'bold'); %#ok<NASGU>
annotation(fig, 'textbox', [0.075 0.905 0.885 0.05], 'EdgeColor','none', ...
    'Color', dim, 'FontSize', 12.5, 'HorizontalAlignment','left', ...
    'VerticalAlignment','middle', 'Interpreter','none', ...
    'String', sprintf('GTO %s tulip  -  %.2fx min-time (%.1f d)  -  energy->fuel homotopy (Bertrand-Epenoy)', ...
                      char(8594), meta.factor, tfDays));

% big epsilon money-metric readout (top-right, inside the panel)
hEps = text(axT, 0.985, 0.90, '', 'Units','normalized', 'FontName','Menlo', ...
    'FontSize', 40, 'FontWeight','bold', 'HorizontalAlignment','right', ...
    'VerticalAlignment','middle', 'Color', epsCol(1));
hEpsLbl = text(axT, 0.985, 0.70, '', 'Units','normalized', 'FontName','Menlo', ...
    'FontSize', 12, 'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
    'Color', dim); %#ok<NASGU>

% live metrics (upper-left, monospace, fixed width)
hMet = text(axT, 0.018, 0.94, '', 'Units','normalized', 'FontName','Menlo', ...
    'FontSize', 12.5, 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Color', fgTxt, 'BackgroundColor', [0.03 0.04 0.06], 'Margin', 7, ...
    'EdgeColor', [0.22 0.26 0.32]);

% --- bottom strip: bang-bang fraction vs epsilon -----------------------------
axB = axes(fig, 'Position', [0.075 0.115 0.885 0.165], 'Color', bg);
hold(axB,'on'); box(axB,'on');
set(axB, 'XColor', dim, 'YColor', dim, 'GridColor', [0.30 0.34 0.40], ...
    'GridAlpha', 0.4, 'FontName','Helvetica', 'FontSize', 11, 'LineWidth', 1.0, ...
    'XDir', 'reverse');                         % eps 1 (left) -> 0 (right): time flows L->R
grid(axB,'on');
xlim(axB, [0 1]);  ylim(axB, [0 105]);
xlabel(axB, 'homotopy parameter  \epsilon   (1 = min-energy    \rightarrow    0 = min-fuel)', ...
       'Color', fgTxt, 'FontSize', 12);
ylabel(axB, 'bang-bang %', 'Color', fgTxt, 'FontSize', 11);
plot(axB, epsK, edgeK, '-', 'Color', [0.45 0.50 0.58], 'LineWidth', 1.4);
plot(axB, epsK, edgeK, 'o', 'MarkerSize', 4, 'MarkerFaceColor',[0.45 0.50 0.58], ...
     'MarkerEdgeColor','none');
hBdot = plot(axB, epsK(1), edgeK(1), 'o', 'MarkerSize', 12, 'LineWidth', 2, ...
     'MarkerFaceColor', epsCol(1), 'MarkerEdgeColor', 'w');

% branding
annotation(fig, 'textbox', [0.60 0.005 0.395 0.035], 'EdgeColor','none', ...
    'Color', [0.42 0.46 0.52], 'FontSize', 10.5, 'HorizontalAlignment','right', ...
    'VerticalAlignment','middle', ...
    'String', 'Coorbital  |  Casey & Koblick  |  CR3BP low-thrust min-fuel');

% --- render ------------------------------------------------------------------
vw = VideoWriter(outStem, 'MPEG-4');  vw.FrameRate = fps;  vw.Quality = 96;  open(vw);
gifFile = [outStem '.gif'];  gifStride = 2;  gifW = 720;  gifMap = [];  gifDelay = 1/16;
tmpPng = [outStem '_tmp.png'];
% FIXED encoder-friendly output size (both divisible by 16): exportgraphics'
% raw size is not H.264-safe and a stride mismatch shears each frame into
% diagonal streaks -- resizing every frame to 1280x720 removes it.
vidHW = [720 1280];

for fi = 1:nFr
    s   = Sfr{fi};   e = epsFr(fi);   ecol = epsCol(max(e,0));
    set(hLine, 'YData', s, 'Color', ecol);
    set(hFill, 'YData', [0 s 0], 'FaceColor', blendbg(ecol, 0.30));
    set(hEps,  'String', sprintf('%s = %4.2f', char(949), e), 'Color', ecol);
    set(hEpsLbl,'String', eps_phase(e));
    set(hMet, 'String', sprintf([ ...
        '  \\DeltaV   = %6.3f km/s\n' ...
        '  prop    = %6.3f kg\n' ...
        '  switches= %6d\n' ...
        '  bang-bang %5.1f%%\n' ...
        '  defect   = %.0e'], ...
        dVFr(fi), propFr(fi), round(swFr(fi)), edgeFr(fi), defFr(fi)));
    set(hBdot, 'XData', e, 'YData', edgeFr(fi), 'MarkerFaceColor', ecol);

    drawnow;
    exportgraphics(fig, tmpPng, 'Resolution', 130, 'BackgroundColor', bg);
    img = imresize(imread(tmpPng), vidHW);          % fixed 720x1280 (H.264-safe)
    writeVideo(vw, img);
    if mod(fi-1, gifStride) == 0
        gr = round(gifW*vidHW(1)/vidHW(2));  gimg = imresize(img,[gr gifW]);
        if isempty(gifMap)
            [gI,gifMap] = rgb2ind(gimg,256,'nodither');
            imwrite(gI,gifMap,gifFile,'gif','LoopCount',Inf,'DelayTime',gifDelay);
        else
            gI = rgb2ind(gimg,gifMap,'nodither');
            imwrite(gI,gifMap,gifFile,'gif','WriteMode','append','DelayTime',gifDelay);
        end
    end
end
close(vw);
if isfile(tmpPng), delete(tmpPng); end
close(fig);
fprintf('WROTE %s.mp4 (%d frames)\n', outStem, nFr);
fprintf('WROTE %s.gif (%d frames)\n', outStem, numel(1:gifStride:nFr));
end

% ---------------------------------------------------------------------------
function str = eps_phase(e)
% Short phase label under the big epsilon readout.
if     e > 0.85, str = 'min-energy (smooth)';
elseif e > 0.10, str = 'sharpening...';
elseif e > 0.0,  str = 'nearly bang-bang';
else,            str = 'min-fuel (bang-bang)';
end
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
