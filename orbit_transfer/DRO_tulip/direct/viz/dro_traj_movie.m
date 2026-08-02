function outFile = dro_traj_movie(o, p, outFile, mopts)
% DRO_TRAJ_MOVIE  Moon-centred movie of a direct DRO->tulip transfer.
%
% Base MATLAB only -- no toolboxes beyond VideoWriter -- so this always runs.
% The companion dro_scene_movie renders the pumpkynPie-style lit scene, which
% needs that class to be on the path.
%
% The satellite path is densified by SPLINE interpolation of the collocation
% nodes (the same device pumpkynPie's animation demo uses), at nDense samples
% per interval. That is a rendering choice: it makes the curve smooth without
% claiming the sub-node path is dynamically exact. Compare with panel (d) of
% plot_dro_diagnostics, which measures how far the true dynamics depart from
% the collocation between nodes.
%
% What the frame shows, Moon-centred and rotating-frame:
%   - the lunar disc to scale, and the altitude floor as a dashed ring;
%   - the departure DRO (green) and the arrival tulip (red), both propagated;
%   - the transfer as a fading tail behind a bright satellite;
%   - a live readout of elapsed time and lunar altitude.
%
% INPUTS:
%   o       - solution struct from casadi_mintime_dro (.X .U .s .tf .minAltKm)
%   p       - params from dro_tulip_endpoints (.muStar .lStar .tStar .rv0DRO
%             .rvfTulip -- the last two optional; drawn only if present)
%   outFile - output path, .mp4 [char]
%   mopts   - struct (optional):
%             .nDense   spline samples per collocation interval [default 6]
%             .fps      frame rate                              [default 30]
%             .nFrames  frames in the movie                     [default 450]
%             .tailFrac tail length as a fraction of t_f        [default 0.18]
%             .spanKm   half-width of the view, km. [] = auto   [default []]
%
% OUTPUTS:
%   outFile - the file actually written [char]
%
% REFERENCES:
%   [1] pumpkynPie demos/lowThrustDRO2Tulip.m -- the scene this mirrors.
%   [2] MEMORY matlab-movie-diagonal-streaks -- frame size is forced to a
%       multiple of 16 (1280x720) because H.264 shears otherwise.

if nargin < 4, mopts = struct(); end
g = @(f,v) local_default(mopts, f, v);
nDense  = g('nDense', 6);
fps     = g('fps', 30);
nFrames = g('nFrames', 450);
tailFrac= g('tailFrac', 0.18);
spanKm  = g('spanKm', []);

mu = p.muStar;  lStar = p.lStar;  tStar = p.tStar;
rMoonKm = 1737.4;

% --- densify by spline, Moon-centred, in km --------------------------------
% Node times: uniform s*tf normally; the solver's own PHYSICAL times under
% Sundman, where the mesh is non-uniform in t by design.
if isfield(o,'tNodes') && ~isempty(o.tNodes), tN = o.tNodes(:).';
else, tN = o.s(:).' * o.tf; end
tD = linspace(tN(1), tN(end), nDense*(numel(tN)-1) + 1);
XD = interp1(tN, o.X(1:3,:).', tD, 'spline').';           % 3 x nD
xs = (XD(1,:) - (1-mu))*lStar;  ys = XD(2,:)*lStar;  zs = XD(3,:)*lStar;
altD = (vecnorm(XD - [1-mu;0;0], 2, 1))*lStar - rMoonKm;

% --- the two terminal orbits, propagated for context -----------------------
[droX, tulX] = local_orbits(p, mu, lStar);

if isempty(spanKm)
    spanKm = 1.08*max([abs(xs), abs(ys), abs(droX(1,:)), abs(droX(2,:))]);
end

% --- figure. The axes are square in DATA units (axis equal) but the frame is
%     16:9, so the box is given a square slice of the frame and the readouts
%     go in the title rather than beside the plot.
fig = figure('Color','k','Position',[80 80 1280 720],'Visible','off', ...
             'InvertHardcopy','off');
ax = axes('Parent',fig,'Color','k','XColor','w','YColor','w','ZColor','w', ...
          'Position',[0.24 0.05 0.52 0.86]);
hold(ax,'on');  grid(ax,'on');  set(ax,'GridColor',[0.3 0.3 0.3]);
axis(ax,'equal');  axis(ax,[-spanKm spanKm -spanKm spanKm]);

th = linspace(0,2*pi,240);
hLeg = gobjects(0);  lLeg = {};
hLeg(end+1) = fill(ax, rMoonKm*cos(th), rMoonKm*sin(th), [0.62 0.60 0.58], ...
    'EdgeColor','none');  lLeg{end+1} = 'Moon';
if isfinite(o.minAltKm)
    hLeg(end+1) = plot(ax, (rMoonKm+o.minAltKm)*cos(th), (rMoonKm+o.minAltKm)*sin(th), ...
        '--', 'Color',[1 0.35 0.35], 'LineWidth',1.0);
    lLeg{end+1} = sprintf('%.0f km floor', o.minAltKm);
end
if ~isempty(droX)
    hLeg(end+1) = plot(ax, droX(1,:), droX(2,:), '-', ...
        'Color',[0.25 0.75 0.35 0.55], 'LineWidth',1.0);  lLeg{end+1} = 'departure DRO';
end
if ~isempty(tulX)
    hLeg(end+1) = plot(ax, tulX(1,:), tulX(2,:), '-', ...
        'Color',[0.90 0.30 0.30 0.55], 'LineWidth',1.0);  lLeg{end+1} = 'arrival tulip';
end

hTail = gobjects(1,24);                       % fading tail, drawn in segments
for k = 1:numel(hTail)
    a = 0.06 + 0.94*(k/numel(hTail))^2;       % oldest faint, newest bright
    hTail(k) = plot(ax, nan, nan, '-', 'Color',[0.45 0.75 1 a], 'LineWidth',1.1+1.6*a);
end
hSat  = plot(ax, nan, nan, 'o', 'MarkerSize',9, 'MarkerFaceColor',[1 1 0.85], ...
             'MarkerEdgeColor',[1 0.9 0.4], 'LineWidth',1.2);
xlabel(ax,'x - x_{Moon}  [km]');  ylabel(ax,'y  [km]');
% The readout is a figure-level annotation, not an axes title: the axes are
% square in data units, so a title would collide with the exponent label.
hTitle = annotation(fig,'textbox',[0.02 0.93 0.96 0.06], 'String','', ...
    'Color','w','EdgeColor','none','Interpreter','tex','FontSize',13, ...
    'HorizontalAlignment','center','VerticalAlignment','middle');
legend(ax, hLeg, lLeg, 'TextColor','w','Color','k','EdgeColor',[0.4 0.4 0.4], ...
    'Location','northeastoutside','AutoUpdate','off','FontSize',9);

vw = VideoWriter(outFile,'MPEG-4');  vw.FrameRate = fps;  vw.Quality = 95;
open(vw);
cleanupVW = onCleanup(@() close(vw)); %#ok<NASGU>

tFrames = linspace(tD(1), tD(end), nFrames);
tailND  = tailFrac*o.tf;
for kf = 1:nFrames
    tNow = tFrames(kf);
    inTail = tD <= tNow & tD >= tNow - tailND;
    idx = find(inTail);
    if numel(idx) >= 2
        edges = round(linspace(idx(1), idx(end), numel(hTail)+1));
        for ks = 1:numel(hTail)
            seg = edges(ks):edges(ks+1);
            set(hTail(ks), 'XData', xs(seg), 'YData', ys(seg));
        end
    end
    xi = interp1(tD, xs, tNow);  yi = interp1(tD, ys, tNow);
    set(hSat, 'XData', xi, 'YData', yi);
    set(hTitle, 'String', sprintf(['DRO \\rightarrow tulip, min-time  |  ' ...
        't = %6.2f d of %6.2f  |  lunar altitude %7.0f km'], ...
        tNow*tStar/86400, o.tf*tStar/86400, interp1(tD, altD, tNow)));
    writeVideo(vw, local_frame(fig));
end
close(fig);
fprintf('  movie -> %s  (%d frames @ %d fps)\n', outFile, nFrames, fps);
end

% ---------------------------------------------------------------------------
function [droX, tulX] = local_orbits(p, mu, lStar)
% LOCAL_ORBITS  Propagate the departure DRO and arrival tulip, Moon-centred km.
% INPUTS: p; mu; lStar   OUTPUTS: droX [3xM], tulX [3xM] (empty if unavailable)
droX = zeros(3,0);  tulX = zeros(3,0);
try
    [~, rvDRO] = pumpkynPie.cr3bp.getDRO(1.0);
    rvDRO = pumpkyn.cr3bp.cont_np(rvDRO, 1.0, mu, 1e-12);
    [~, rv] = pumpkyn.cr3bp.prop(1.0, rvDRO, mu);
    droX = [(rv(:,1).'-(1-mu))*lStar; rv(:,2).'*lStar; rv(:,3).'*lStar];
catch, end
try
    tauf = 5*2*pi/6;
    [~, rvT] = pumpkyn.cr3bp.getTulip(tauf, 7, -1);
    rvT = pumpkyn.cr3bp.cont_np(rvT, tauf, mu, 1e-12);
    [~, rv] = pumpkyn.cr3bp.prop(tauf, rvT, mu);
    tulX = [(rv(:,1).'-(1-mu))*lStar; rv(:,2).'*lStar; rv(:,3).'*lStar];
catch, end
end

% ---------------------------------------------------------------------------
function F = local_frame(fig)
% LOCAL_FRAME  getframe forced to 1280x720 so H.264 does not shear the frame.
% INPUTS: fig   OUTPUTS: F (struct with .cdata [720x1280x3 uint8])
F = getframe(fig);
if ~isequal(size(F.cdata,1), 720) || ~isequal(size(F.cdata,2), 1280)
    F.cdata = imresize_nn(F.cdata, [720 1280]);
end
end

% ---------------------------------------------------------------------------
function B = imresize_nn(A, sz)
% IMRESIZE_NN  Nearest-neighbour resize with no Image Processing Toolbox.
% INPUTS: A [MxNx3 uint8]; sz [1x2]   OUTPUTS: B [sz(1) x sz(2) x 3 uint8]
ri = max(1, min(size(A,1), round(linspace(1, size(A,1), sz(1)))));
ci = max(1, min(size(A,2), round(linspace(1, size(A,2), sz(2)))));
B = A(ri, ci, :);
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
% INPUTS: s; f; dflt   OUTPUTS: v
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
