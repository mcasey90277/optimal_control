function movie_landing_cinematic(out, sol, P, outfile, opts)
% MOVIE_LANDING_CINEMATIC  Cinematic 1920x1080 landing movie: lit 3D
% booster, throttle-driven plume, animated camera, slow-motion finale and
% an instrument strip.
%
% This is a PRESENTATION product and is deliberately SEPARATE from
% viz/movie_landing.m, which stays the campaign's DIAGNOSTIC movie (2-D
% altitude view + throttle trace) and is what the front door
% run_booster_landing renders. Nothing here is part of the front-door
% pipeline: render it by hand from a saved results/booster_run.mat (see
% README, "Cinematic movie").
%
% WHAT IS DRAWN, AND FROM WHAT
%   Everything animated comes from the CLOSED-LOOP trace `out` (out.t,
%   out.X, out.Tcmd) -- i.e. what actually flew under the TVLQR tracker,
%   not the open-loop guidance solution. `sol` is used only for the
%   headline numbers on the title card (sol.tf, sol.mf), which are
%   guidance quantities and are labelled as such.
%
% BODY-AXIS CONVENTION (important, and stated on screen)
%   The campaign model is 3-DOF: a point mass with a thrust VECTOR and no
%   attitude state at all. A drawn booster nevertheless needs an axis, so
%   the body is drawn ALIGNED WITH THE CURRENT THRUST DIRECTION -- the
%   standard visual convention for 3-DOF powered-descent animations, and
%   the one attitude a real vehicle would have to hold to deliver that
%   thrust vector (gimbal angles aside). It is a rendering choice, not a
%   simulation output; no attitude dynamics, rate limit or gimbal limit is
%   implied or checked. The overlay says so in small print.
%
% SCALE
%   The booster is drawn at TRUE scale by default (opts.boosterScale = 1:
%   P.boosterDiam-ish 3.7 m x 45 m, the Falcon-9 first stage). At the top
%   of the trajectory (2 km) that is a small object, which is honest --
%   the camera zooms in as the vehicle descends rather than the vehicle
%   being inflated (house rule: honest scale, not cropped scale). If a
%   caller does set opts.boosterScale > 1, the exaggeration factor is
%   PRINTED ON SCREEN for the whole movie; it is never silently applied.
%
% FRAME SIZE is LOCKED to 1920x1080 (both divisible by 16) via imresize on
% every frame, exactly as movie_landing.m locks 1280x720: an unlocked
% figure size produces the diagonal-color-streak H.264 shear artifact in
% the exported MP4 (house lesson, see matlab-polished-graphics skill /
% MEMORY.md "MATLAB movie diagonal streaks").
%
% TIMELINE (playback seconds, at opts.fps)
%   [title card]  opts.titleSec
%   [fade-in]     opts.fadeSec (black overlay ramping off the live scene)
%   [flight]      real time until opts.slowLead seconds before touchdown,
%                 then a cosine ramp down to opts.slowRate playback speed
%                 for the finale (a "0.25x" tag is shown while slowed)
%   [hold]        opts.holdSec on the touchdown frame
%
% INPUTS:
%   out     - sim_closed_loop trace: .t [Mx1 s], .X [Mx7: r(3) v(3) m],
%             .Tcmd [Mx3 N], .td (uses .vtd, .miss, .landed)
%   sol     - guidance solution (solve_pdg_colloc): .tf, .mf -- title-card
%             headline numbers only
%   P       - booster_params: .Tmax .Tmin .etaT .m0 .mdry .pad_radius
%   outfile - output .mp4 path
%   opts    - (optional) all fields optional:
%             .fps          [def 60]    frame rate
%             .titleSec     [def 2.0]   title-card duration [s]
%             .fadeSec      [def 0.6]   card->scene fade [s]
%             .holdSec      [def 2.2]   touchdown hold [s]
%             .slowLead     [def 3.0]   seconds of FLIGHT time run slow
%             .slowRate     [def 0.25]  playback speed during the finale
%             .slowRamp     [def 0.9]   flight-seconds spent easing into it
%             .trailSec     [def 6.0]   fading trajectory ribbon length [s]
%             .boosterScale [def 1]     >1 exaggerates the vehicle AND
%                                       prints the factor on screen
%             .mcRate       [def []]    Monte-Carlo success rate, fraction
%                                       (0.995) or percent (99.5); printed
%                                       on the title card when given
%             .mcRuns       [def []]    run count for that rate
%             .title        [def 'Falcon-9-class Powered Descent']
%             .subtitle     [def 'Certified Min-Fuel Hoverslam']
%             .stills       [def []]    PLAYBACK times [s]; when non-empty
%                                       the function writes one PNG per
%                                       entry (see .stillDir) and writes
%                                       NO video -- the key-frame /
%                                       poster-frame path, used to judge
%                                       composition without paying for a
%                                       full render
%             .stillDir     [def fileparts(outfile)] where those PNGs go
% OUTPUTS: none (writes outfile, or the stills)
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-08-08-booster-landing-design.md
%   [2] matlab-polished-graphics skill -- divide-by-16 frame lock, dark
%       theme, fading trails, static path sketch, honest scale, branding
if nargin < 5, opts = struct(); end
def = struct('fps',60, 'titleSec',2.0, 'fadeSec',0.6, 'holdSec',2.2, ...
             'slowLead',3.0, 'slowRate',0.25, 'slowRamp',0.9, ...
             'trailSec',6.0, 'boosterScale',1, 'mcRate',[], 'mcRuns',[], ...
             'title','Falcon-9-class Powered Descent', ...
             'subtitle','Certified Min-Fuel Hoverslam', ...
             'stills',[], 'stillDir','');
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts, fn{k}), opts.(fn{k}) = def.(fn{k}); end
end
if isempty(opts.stillDir), opts.stillDir = fileparts(outfile); end

W = 1920;  H = 1080;                     % divide-by-16 lock (see header)

%% ---------------- Palette (dark theme) ----------------
col.bg      = [0.035 0.043 0.055];   % scene background (near-black blue)
col.hudBg   = [0.070 0.080 0.098];   % instrument strip
col.ground  = [0.105 0.115 0.135];   % ground plane
col.grid    = [0.190 0.215 0.255];   % ground grid lines
col.pad     = [0.180 0.850 0.450];   % aviation-safety green (same
                                     % convention as movie_landing.m)
col.body    = [0.760 0.790 0.830];   % booster skin
col.trail   = [0.450 0.760 1.000];   % trajectory ribbon (cool blue)
col.text    = [0.930 0.945 0.960];
col.dim     = [0.520 0.560 0.610];
col.accent  = [1.000 0.520 0.150];   % throttle / plume orange
col.hot     = [1.000 0.980 0.930];   % white-hot
col.warn    = [1.000 0.760 0.250];

%% ---------------- Flight data ----------------
t    = out.t(:);
R3   = out.X(:,1:3);
V3   = out.X(:,4:6);
Mkg  = out.X(:,7);
Tv   = out.Tcmd;
Tmag = sqrt(sum(Tv.^2, 2));
thr  = Tmag / P.Tmax;                     % throttle, fraction of Tmax
tEnd = t(end);
spd  = sqrt(sum(V3.^2, 2));

% Bang-bang switch instant: first upward crossing of the midpoint between
% the two certified levels (Tmin/Tmax and the guidance ceiling etaT).
thrSwitchLevel = 0.5 * (P.Tmin/P.Tmax + P.etaT);
kSw = find(thr > thrSwitchLevel, 1, 'first');
if isempty(kSw), tSw = NaN; else, tSw = t(kSw); end

fuelKg = P.m0 - sol.mf;
vtd    = out.td.vtd;

% Booster geometry [m] -- Falcon-9 first stage, times the (documented,
% on-screen) exaggeration factor.
sc    = opts.boosterScale;
bDiam = 3.7 * sc;
bLen  = 45.0 * sc;

%% ---------------- Playback <-> flight time map (slow-motion) ----------
% Playback speed as a function of FLIGHT time: 1x until slowLead+slowRamp
% before the end, cosine-eased down to slowRate, then held. Playback time
% is the integral of 1/speed, inverted onto a uniform frame grid.
tSlow0 = max(tEnd - opts.slowLead - opts.slowRamp, 0);
tSlow1 = max(tEnd - opts.slowLead, 0);
td_    = linspace(0, tEnd, 6000).';
sp_    = ones(size(td_));
inR    = td_ > tSlow0 & td_ < tSlow1;
u_     = (td_(inR) - tSlow0) / max(tSlow1 - tSlow0, eps);
sp_(inR) = 1 + (opts.slowRate - 1) * 0.5 .* (1 - cos(pi*u_));
sp_(td_ >= tSlow1) = opts.slowRate;
pb_    = cumtrapz(td_, 1 ./ sp_);         % playback time at each flight t
nFly   = max(round(pb_(end) * opts.fps), 2);
pbq    = linspace(0, pb_(end), nFly).';
tq     = interp1(pb_, td_, pbq, 'linear');
tq(1)  = 0;  tq(end) = tEnd;
spq    = interp1(td_, sp_, tq, 'linear'); % playback speed per frame

nCard = round(opts.titleSec * opts.fps);
nFade = round(opts.fadeSec  * opts.fps);
nHold = round(opts.holdSec  * opts.fps);
nTot  = nCard + nFly + nHold;             % fade overlays the first frames

%% ---------------- Per-frame state (interpolated ONCE) ----------------
Rq   = interp1(t, R3,  tq);
Vq   = interp1(t, V3,  tq);
Mq   = interp1(t, Mkg, tq);
Tq   = interp1(t, Tv,  tq);
thq  = sqrt(sum(Tq.^2, 2)) / P.Tmax;
spdq = sqrt(sum(Vq.^2, 2));

% Unit thrust direction per frame (vertical fallback if ever degenerate --
% unreachable here, |T| >= Tmin = 338 kN throughout).
Tn   = sqrt(sum(Tq.^2, 2));
Uq   = Tq ./ max(Tn, 1e-9);
Uq(Tn < 1e-9, :) = repmat([0 0 1], nnz(Tn < 1e-9), 1);

%% ---------------- Camera path (precomputed, no autoscale jitter) ------
% Wide establishing shot easing to a low pad-side angle. The target rides
% the vehicle (damped, so the ground stays in frame) and the camera
% distance shrinks with altitude, so the booster grows naturally instead
% of being scaled up.
zq   = Rq(:,3);
% Ease on FLIGHT-time progress, not frame index: the slow-motion finale
% owns ~40% of the frames, and easing per frame would spend 40% of the
% camera move on the last 3 seconds (a lurch exactly where the shot should
% settle).
easeS = 0.5 * (1 - cos(pi * min(max(tq/tEnd,0),1)));    % smooth 0 -> 1
camVA = 34;                                             % view angle [deg]
% Framing is driven by a desired VIEW HEIGHT in metres at the target
% plane, from which the camera distance follows.
%   WIDE  - the establishing shot: whole 2 km descent and the pad in one
%           frame. A true-scale 45 m booster is ~20 px tall here, which is
%           what the tracking reticle below is for.
%   TIGHT - the working shot: view height ~0.30 of the altitude, so the
%           vehicle holds ~50-180 px and reads as a rocket.
% The movie pushes in from one to the other over the first few seconds
% (a camera move, not a scale change -- the vehicle stays true-scale).
viewHwide  = 1.15*zq + 300;
viewHtight = max(0.30*zq + 165, 165);
ePush = 0.5 * (1 - cos(pi * min(max((tq - 0.4)/3.6, 0), 1)));
viewH = (1 - ePush).*viewHwide + ePush.*viewHtight;     % [m]
camD  = viewH / (2*tand(camVA/2));                      % [m]
% The target sits BELOW the vehicle by a fraction of the view height, so
% the booster rides high in frame with the ground underneath it. That
% offset goes slightly NEGATIVE by touchdown: with the camera nearly at
% ground level the horizon lands mid-frame and the vehicle-plus-pad sits
% centred, instead of being pinned to the top third over a huge empty
% apron.
offF  = 0.32 + (-0.10 - 0.32) * easeS;
azDeg = -52 + (-108 - -52) * easeS;
elDeg =  21 + (  6  -  21) * easeS;
camT  = [0.65*Rq(:,1), 0.65*Rq(:,2), zq - offF .* viewH];
azR   = deg2rad(azDeg);  elR = deg2rad(elDeg);
camP  = camT + camD .* [cos(elR).*cos(azR), cos(elR).*sin(azR), sin(elR)];

%% ---------------- Figure + static scene ----------------
fig = figure('Visible','off', 'Position',[60 60 W H], 'Color', col.bg, ...
             'InvertHardcopy','off');

ax = axes(fig, 'Position',[0 0.150 1 0.850], 'Color', col.bg);
hold(ax,'on');  axis(ax,'off');  axis(ax,'vis3d');
set(ax, 'DataAspectRatio',[1 1 1], 'Projection','perspective', ...
        'XLim',[-1.2e4 1.2e4], 'YLim',[-1.2e4 1.2e4], 'ZLim',[-60 2600], ...
        'CameraUpVector',[0 0 1], 'CameraPositionMode','manual', ...
        'CameraTargetMode','manual', 'CameraViewAngleMode','manual', ...
        'CameraViewAngle', camVA, 'Clipping','off');

% Ground: one big plate reaching well past the widest shot, shaded so it
% DARKENS toward the horizon into the background colour. That radial fade
% is doing the job atmospheric haze does in a real long-lens shot -- it
% keeps the plate's far edge from reading as a hard line across the sky in
% the establishing frame, without a fog model MATLAB does not have.
gRad = 1.1e4;
gv1  = linspace(-gRad, gRad, 81);
[GX, GY] = meshgrid(gv1, gv1);
wFade = exp(-sqrt(GX.^2 + GY.^2) / 2200);
gC = zeros(size(GX,1), size(GX,2), 3);
for kc = 1:3
    gC(:,:,kc) = col.bg(kc) + (col.ground(kc) - col.bg(kc)) * wFade;
end
surf(ax, GX, GY, zeros(size(GX)), gC, 'EdgeColor','none', 'FaceColor','interp', ...
     'FaceLighting','none');
% Grid lines only over the working area -- a grid out to 11 km would be a
% moire field, and the grid is there to give the near ground scale.
gx = [-1200 1400];  gy = [-1300 1300];
for gv = gx(1):100:gx(2)
    plot3(ax, [gv gv], gy, [0.05 0.05], '-', 'Color',[col.grid 0.5], 'LineWidth',0.5);
end
for gv = gy(1):100:gy(2)
    plot3(ax, gx, [gv gv], [0.05 0.05], '-', 'Color',[col.grid 0.5], 'LineWidth',0.5);
end

% Landing pad: filled disk + concentric rings, outer ring = P.pad_radius
thP = linspace(0, 2*pi, 121);
patch(ax, 'XData', P.pad_radius*cos(thP), 'YData', P.pad_radius*sin(thP), ...
      'ZData', 0.10*ones(size(thP)), 'FaceColor',[0.09 0.16 0.12], ...
      'EdgeColor','none');
for rr = [P.pad_radius, 2*P.pad_radius/3, P.pad_radius/3]
    aA = 0.95 * (rr == P.pad_radius) + 0.45 * (rr ~= P.pad_radius);
    plot3(ax, rr*cos(thP), rr*sin(thP), 0.20*ones(size(thP)), '-', ...
          'Color',[col.pad aA], 'LineWidth', 2.2 - 1.0*(rr ~= P.pad_radius));
end
for aM = 0:pi/2:3*pi/2                       % four cross ticks
    plot3(ax, [0.45 1.0]*P.pad_radius*cos(aM), [0.45 1.0]*P.pad_radius*sin(aM), ...
          [0.20 0.20], '-', 'Color',[col.pad 0.8], 'LineWidth',1.6);
end

% Static path sketch: the whole closed-loop track, thin and dim
plot3(ax, R3(:,1), R3(:,2), R3(:,3), '-', 'Color',[col.trail 0.22], 'LineWidth',0.9);
plot3(ax, R3(:,1), R3(:,2), zeros(size(R3,1),1), '-', ...
      'Color',[col.dim 0.30], 'LineWidth',0.8);        % ground track

% Lighting: fixed infinite lights (a camlight would have to be re-aimed
% every frame and flickers when the camera eases).
light(ax, 'Style','infinite', 'Position',[-0.7 -1.0 0.9], 'Color',[0.95 0.95 1.0]);
light(ax, 'Style','infinite', 'Position',[ 1.0  0.6 0.4], 'Color',[0.30 0.34 0.42]);

%% ---------------- Dynamic scene handles (updated, never recreated) ----
nTrail = 14;
hTrail = gobjects(nTrail,1);
for k = 1:nTrail
    hTrail(k) = plot3(ax, nan, nan, nan, '-', 'Color',[col.trail 0.2], 'LineWidth',2.2);
end

% Tracking reticle: a camera-facing ring around the vehicle whose radius
% is a fixed fraction of the VIEW HEIGHT, so it holds a constant size on
% screen. This is the honest answer to "a true-scale 45 m booster is 16
% pixels tall in the establishing shot": annotate the vehicle instead of
% inflating it. It fades out as the shot tightens and the real vehicle
% takes over.
retTh = linspace(0, 2*pi, 73);
hRet  = plot3(ax, nan, nan, nan, '-', 'Color',[col.text 0.5], 'LineWidth',1.4);

% Booster body: unit cylinder (radius 1, height 1 along +z) + top cap.
nSide = 32;
thB   = linspace(0, 2*pi, nSide+1);
cylX  = [cos(thB); cos(thB)];
cylY  = [sin(thB); sin(thB)];
cylZ  = [zeros(1,nSide+1); ones(1,nSide+1)];
hBody = surf(ax, nan(2,nSide+1), nan(2,nSide+1), nan(2,nSide+1), ...
             'FaceColor', col.body, 'EdgeColor','none', 'FaceLighting','gouraud', ...
             'AmbientStrength',0.30, 'DiffuseStrength',0.85, ...
             'SpecularStrength',0.85, 'SpecularExponent',22, 'BackFaceLighting','lit');
capX  = [zeros(1,nSide+1); cos(thB)];
capY  = [zeros(1,nSide+1); sin(thB)];
capZ  = ones(2,nSide+1);
hCap  = surf(ax, nan(2,nSide+1), nan(2,nSide+1), nan(2,nSide+1), ...
             'FaceColor', col.body*0.62, 'EdgeColor','none', 'FaceLighting','gouraud');
% Dark engine-bay band at the base, so the vehicle reads as a rocket and
% not as a featureless rod.
bndX  = [cos(thB); cos(thB)]*1.01;
bndY  = [sin(thB); sin(thB)]*1.01;
bndZ  = [zeros(1,nSide+1); 0.09*ones(1,nSide+1)];
hBand = surf(ax, nan(2,nSide+1), nan(2,nSide+1), nan(2,nSide+1), ...
             'FaceColor',[0.16 0.17 0.19], 'EdgeColor','none', 'FaceLighting','gouraud');

% Plume: two nested translucent cones, apex at the engine, pointing along
% -thrust. Alpha decays along the length so the tip dissolves.
nRing = 26;
sRing = linspace(0, 1, nRing).';
outR  = 0.55 + 3.85*sRing;             % x plume base radius (flares out)
inR   = 0.22 + 1.55*sRing;
outA  = 0.55 * (1 - sRing).^1.10;
inA   = 0.90 * (1 - sRing).^1.45;
hPlumeO = surf(ax, nan(nRing,nSide+1), nan(nRing,nSide+1), nan(nRing,nSide+1), ...
               'FaceColor', col.accent, 'EdgeColor','none', 'FaceLighting','none', ...
               'FaceAlpha','interp', 'AlphaDataMapping','none', ...
               'AlphaData', repmat(outA,1,nSide+1));
hPlumeI = surf(ax, nan(nRing,nSide+1), nan(nRing,nSide+1), nan(nRing,nSide+1), ...
               'FaceColor', col.hot, 'EdgeColor','none', 'FaceLighting','none', ...
               'FaceAlpha','interp', 'AlphaDataMapping','none', ...
               'AlphaData', repmat(inA,1,nSide+1));

% Ground wash: once the plume is long enough to reach the pad, the exhaust
% spreads instead of continuing into the dirt. Drawn as a flat glowing
% disc whose radius and opacity grow with the overlap -- the plume itself
% is CLIPPED at the ground in the same frame, so the two together read as
% one jet turning at the surface.
% A hard-edged disc reads as a SHADOW, not as light, so the wash is a
% polar mesh with the opacity falling off radially -- a soft glow.
washS = linspace(0, 1, 12).';
washA = 0.75 * (1 - washS).^1.8;
hWash = surf(ax, nan(12,nSide+1), nan(12,nSide+1), nan(12,nSide+1), ...
             'FaceColor', col.accent, 'EdgeColor','none', 'FaceLighting','none', ...
             'FaceAlpha','interp', 'AlphaDataMapping','none', ...
             'AlphaData', repmat(washA,1,nSide+1), 'Visible','off');

%% ---------------- Instrument strip ----------------
axH = axes(fig, 'Position',[0 0 1 0.150], 'Color', col.hudBg, ...
           'XLim',[0 1], 'YLim',[0 1], 'XTick',[], 'YTick',[], ...
           'XColor','none', 'YColor','none');
hold(axH,'on');
plot(axH, [0 1], [0.985 0.985], '-', 'Color',[col.dim 0.55], 'LineWidth',1.2);

% -- throttle gauge
gX0 = 0.035;  gX1 = 0.265;  gY0 = 0.34;  gY1 = 0.56;
text(axH, gX0, 0.80, 'THROTTLE   (T / Tmax)', 'Color', col.dim, ...
     'FontSize',13, 'FontWeight','bold', 'FontName','Helvetica');
rectangle(axH, 'Position',[gX0 gY0 gX1-gX0 gY1-gY0], 'FaceColor',[0.13 0.14 0.17], ...
          'EdgeColor',[col.dim 0.7], 'LineWidth',1.0);
hThrFill = rectangle(axH, 'Position',[gX0 gY0 1e-6 gY1-gY0], ...
                     'FaceColor', col.accent, 'EdgeColor','none');
gmk = [P.Tmin/P.Tmax, P.etaT, 1.0];
gml = {sprintf('%.2f', P.Tmin/P.Tmax), sprintf('%.2f', P.etaT), '1.00'};
gtx = {'Tmin', 'etaT', 'Tmax'};
for k = 1:3
    xk = gX0 + (gX1-gX0)*gmk(k);
    plot(axH, [xk xk], [gY0-0.06 gY1+0.06], '-', 'Color',[col.text 0.85], 'LineWidth',1.2);
    text(axH, xk, 0.20, sprintf('%s\n%s', gml{k}, gtx{k}), 'Color', col.dim, ...
         'FontSize',10, 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
         'FontName','Helvetica');
end
hThrVal = text(axH, gX1+0.012, 0.45, '0.00', 'Color', col.accent, 'FontSize',26, ...
               'FontWeight','bold', 'FontName','FixedWidth', 'VerticalAlignment','middle');

% -- digital readouts. THRUST TILT earns its place: the drawn booster
% visibly leans (up to 21 deg off vertical on this trajectory, at both
% ends of the burn), and without a number next to it a viewer reads that
% as a rendering error rather than as the primer direction the min-fuel
% solution actually commands.
hudReadout(axH, 0.325, 'ALTITUDE');
hAlt = text(axH, 0.325, 0.40, '', 'Color', col.text, 'FontSize',33, ...
            'FontWeight','bold', 'FontName','FixedWidth', 'VerticalAlignment','middle');
hudReadout(axH, 0.475, 'SPEED');
hSpd = text(axH, 0.475, 0.40, '', 'Color', col.text, 'FontSize',33, ...
            'FontWeight','bold', 'FontName','FixedWidth', 'VerticalAlignment','middle');
hudReadout(axH, 0.640, 'THRUST TILT');
hTlt = text(axH, 0.640, 0.40, '', 'Color', col.text, 'FontSize',33, ...
            'FontWeight','bold', 'FontName','FixedWidth', 'VerticalAlignment','middle');

% -- propellant bar
pX0 = 0.775;  pX1 = 0.965;  pY0 = 0.34;  pY1 = 0.56;
text(axH, pX0, 0.80, 'PROPELLANT', 'Color', col.dim, 'FontSize',13, ...
     'FontWeight','bold', 'FontName','Helvetica');
rectangle(axH, 'Position',[pX0 pY0 pX1-pX0 pY1-pY0], 'FaceColor',[0.13 0.14 0.17], ...
          'EdgeColor',[col.dim 0.7], 'LineWidth',1.0);
hPrpFill = rectangle(axH, 'Position',[pX0 pY0 pX1-pX0 pY1-pY0], ...
                     'FaceColor',[0.35 0.72 0.95], 'EdgeColor','none');
hPrpVal = text(axH, pX0, 0.20, '', 'Color', col.dim, 'FontSize',13, ...
               'FontName','FixedWidth', 'VerticalAlignment','top');

%% ---------------- Overlay (titles, callouts, branding) ----------------
axO = axes(fig, 'Position',[0 0 1 1], 'Color','none', 'XLim',[0 1], ...
           'YLim',[0 1], 'XTick',[], 'YTick',[], 'XColor','none', 'YColor','none');
hold(axO,'on');
text(axO, 0.022, 0.965, upper(opts.title), 'Color', col.text, 'FontSize',22, ...
     'FontWeight','bold', 'FontName','Helvetica');
text(axO, 0.022, 0.933, sprintf(['3-DOF min-fuel guidance  |  TVLQR tracker  |  ' ...
     'closed-loop truth sim'] ), 'Color', col.dim, 'FontSize',14, 'FontName','Helvetica');
hClock = text(axO, 0.978, 0.965, '', 'Color', col.text, 'FontSize',26, ...
              'FontWeight','bold', 'FontName','FixedWidth', 'HorizontalAlignment','right');
hSlow  = text(axO, 0.978, 0.925, '0.25x', 'Color', col.warn, 'FontSize',20, ...
              'FontWeight','bold', 'FontName','Helvetica', ...
              'HorizontalAlignment','right', 'Visible','off');
hEvent = text(axO, 0.270, 0.470, '', 'Color', col.accent, 'FontSize',64, ...
              'FontWeight','bold', 'FontName','Helvetica', ...
              'HorizontalAlignment','center', 'Visible','off');
noteStr = 'body axis drawn along the thrust vector (3-DOF model has no attitude state)';
if sc ~= 1
    noteStr = sprintf('booster drawn at %.3gx scale  |  %s', sc, noteStr);
end
text(axO, 0.022, 0.175, noteStr, 'Color', col.dim, 'FontSize',12, 'FontName','Helvetica');
text(axO, 0.978, 0.175, 'Coorbital  |  M. Casey', 'Color', col.dim*0.85, ...
     'FontSize',12, 'FontName','Helvetica', 'HorizontalAlignment','right');

% Fade / card layer: a full-frame black patch whose alpha is animated.
hFade = patch(axO, 'XData',[0 1 1 0], 'YData',[0 0 1 1], 'FaceColor', col.bg, ...
              'EdgeColor','none', 'FaceAlpha',1);

%% ---------------- Title card ----------------
cardH = build_card(axO, opts, col, sol, fuelKg, vtd, P);

%% ---------------- Render ----------------
doStills = ~isempty(opts.stills);
if doStills
    if ~exist(opts.stillDir,'dir'), mkdir(opts.stillDir); end
    frameIdx = unique(max(min(round(opts.stills(:).'*opts.fps)+1, nTot), 1));
    vw = [];
else
    frameIdx = 1:nTot;
    vw = VideoWriter(outfile, 'MPEG-4');  vw.FrameRate = opts.fps;
    vw.Quality = 97;  open(vw);
end

nStill = 0;
for kf = frameIdx
    if kf <= nCard                                   % ---- title card ----
        set(cardH, 'Visible','on');
        set(hFade, 'FaceAlpha', 1);       % opaque: the live scene is hidden
        kFly = 1;
    else                                             % ---- live scene ----
        set(cardH, 'Visible','off');
        kFly = min(kf - nCard, nFly);
        fa = max(1 - (kf - nCard) / max(nFade,1), 0);
        set(hFade, 'FaceAlpha', fa);
    end
    % Engine cutoff through the touchdown hold: the sim ends AT touchdown,
    % so replaying its last frame verbatim would leave a landed vehicle
    % firing its engine for two seconds. The plume is faded out instead.
    pMul = 1;
    if kf > nCard + nFly
        pMul = max(1 - (kf - nCard - nFly) / max(0.5*opts.fps, 1), 0);
    end
    update_scene(kFly, pMul);
    drawnow;
    frame = getframe(fig);
    frame.cdata = imresize(frame.cdata, [H W]);      % divide-by-16 lock
    if doStills
        nStill = nStill + 1;
        imwrite(frame.cdata, fullfile(opts.stillDir, sprintf('still_%02d.png', nStill)));
    else
        writeVideo(vw, frame);
    end
end
if ~doStills, close(vw); end
close(fig);

%% ================= nested updater =================
    function update_scene(kk, pMul)
        % kk   - frame index into the per-frame flight arrays
        % pMul - plume scale (1 in flight, ramped to 0 through the hold)
        rk = Rq(kk,:).';   uk = Uq(kk,:).';   thk = thq(kk);
        % --- body: rotate unit cylinder [0;0;1] -> uk, base at rk
        Rot = rot_to(uk);
        [bx, by, bz] = xform(Rot, rk, cylX*(bDiam/2), cylY*(bDiam/2), cylZ*bLen);
        set(hBody, 'XData',bx, 'YData',by, 'ZData',bz);
        [cx, cy, cz] = xform(Rot, rk, capX*(bDiam/2), capY*(bDiam/2), capZ*bLen);
        set(hCap,  'XData',cx, 'YData',cy, 'ZData',cz);
        [ex, ey, ez] = xform(Rot, rk, bndX*(bDiam/2), bndY*(bDiam/2), bndZ*bLen);
        set(hBand, 'XData',ex, 'YData',ey, 'ZData',ez);
        % --- plume: length and colour ramp with throttle
        % Plume length is strictly PROPORTIONAL to |T|/Tmax (the brief's
        % rule), scaled so full throttle is ~1.75 booster lengths: longer
        % than that and the jet simply runs off the bottom of the frame in
        % the terminal shots, which hides the very throttle change the
        % length is there to show.
        pLenFree = 2.0 * thk * bLen * pMul;
        dGnd = rk(3) / max(uk(3), 0.20);          % length along -u to z=0
        pLen = min(pLenFree, dGnd + 0.04*bLen);   % clip the jet at the pad
        pRad = 0.60 * bDiam/2;
        wHot = min(max((thk - P.Tmin/P.Tmax) / max(1 - P.Tmin/P.Tmax, eps), 0), 1);
        [ox, oy, oz] = xform(Rot, rk, (outR*pRad).*cos(thB), (outR*pRad).*sin(thB), ...
                             repmat(-sRing*pLen, 1, nSide+1));
        set(hPlumeO, 'XData',ox, 'YData',oy, 'ZData',oz, ...
                     'AlphaData', repmat(outA*pMul, 1, nSide+1), ...
                     'FaceColor', col.accent + ([1 0.80 0.45] - col.accent)*wHot);
        [ix, iy, iz] = xform(Rot, rk, (inR*pRad).*cos(thB), (inR*pRad).*sin(thB), ...
                             repmat(-sRing*pLen*0.72, 1, nSide+1));
        set(hPlumeI, 'XData',ix, 'YData',iy, 'ZData',iz, ...
                     'AlphaData', repmat(inA*pMul, 1, nSide+1), ...
                     'FaceColor', [1 0.72 0.32] + (col.hot - [1 0.72 0.32])*wHot);
        % --- ground wash where the jet is cut off by the pad
        ovl = pLenFree - dGnd;
        if ovl > 0 && pMul > 0.02
            wR  = 0.32*ovl + 2.0*pRad;
            wF  = min(ovl / max(pLenFree, eps), 1);
            set(hWash, 'XData', rk(1) - dGnd*uk(1) + (washS*wR).*cos(thB), ...
                       'YData', rk(2) - dGnd*uk(2) + (washS*wR).*sin(thB), ...
                       'ZData', 0.30*ones(12, nSide+1), 'Visible','on', ...
                       'AlphaData', repmat(washA*wF*pMul, 1, nSide+1), ...
                       'FaceColor', col.accent + (col.hot - col.accent)*wHot);
        else
            set(hWash, 'Visible','off');
        end
        % --- fading trajectory ribbon (last opts.trailSec of flight)
        tk = tq(kk);
        mask = t >= tk - opts.trailSec & t <= tk;
        trX = [R3(mask,:); rk.'];
        nP  = size(trX,1);
        if nP >= 2
            edges = round(linspace(1, nP, nTrail+1));
            for kt = 1:nTrail
                iA = edges(kt);  iB = edges(kt+1);
                aT = 0.85 * (kt/nTrail)^2;
                set(hTrail(kt), 'XData',trX(iA:iB,1), 'YData',trX(iA:iB,2), ...
                    'ZData',trX(iA:iB,3), 'Color',[col.trail aT]);
            end
        else
            set(hTrail, 'XData',nan, 'YData',nan, 'ZData',nan);
        end
        % --- camera
        set(ax, 'CameraPosition', camP(kk,:), 'CameraTarget', camT(kk,:));
        % --- tracking reticle (camera-facing ring, constant screen size)
        vh   = viewH(kk);
        aRet = min(max((vh - 200)/700, 0), 1) * 0.70;
        if aRet > 0.02
            cdir = camT(kk,:).' - camP(kk,:).';   cdir = cdir / norm(cdir);
            e1 = cross(cdir, [0;0;1]);            e1 = e1 / max(norm(e1), 1e-9);
            e2 = cross(e1, cdir);
            ctr = rk + 0.5*bLen*uk;
            rRet = 0.055 * vh;
            ring = ctr + rRet*(e1*cos(retTh) + e2*sin(retTh));
            tick = [];
            for kg = 0:3
                aG = kg*pi/2;
                pA = ctr + 1.02*rRet*(e1*cos(aG) + e2*sin(aG));
                pB = ctr + 1.35*rRet*(e1*cos(aG) + e2*sin(aG));
                tick = [tick, pA, pB, nan(3,1)];  %#ok<AGROW>
            end
            seg = [ring, nan(3,1), tick];
            set(hRet, 'XData',seg(1,:), 'YData',seg(2,:), 'ZData',seg(3,:), ...
                'Color', [col.text aRet], 'Visible','on');
        else
            set(hRet, 'Visible','off');
        end
        % --- instruments
        set(hThrFill, 'Position', [gX0 gY0 max((gX1-gX0)*thk,1e-6) gY1-gY0]);
        set(hThrVal, 'String', sprintf('%4.2f', thk), ...
                     'Color', col.accent + (col.hot - col.accent)*wHot);
        set(hAlt, 'String', sprintf('%6.0f m', max(Rq(kk,3),0)));
        set(hSpd, 'String', sprintf('%6.1f m/s', spdq(kk)));
        set(hTlt, 'String', sprintf('%5.1f deg', acosd(min(max(uk(3),-1),1))));
        pf = max(min((Mq(kk) - P.mdry) / (P.m0 - P.mdry), 1), 0);
        set(hPrpFill, 'Position', [pX0 pY0 max((pX1-pX0)*pf,1e-6) pY1-pY0]);
        set(hPrpVal, 'String', sprintf('%5.0f kg remaining   (%4.0f kg burned)', ...
                                       Mq(kk)-P.mdry, P.m0-Mq(kk)));
        set(hClock, 'String', sprintf('T+%5.2f s', tk));
        % --- slow-motion tag
        if spq(kk) < 0.95
            set(hSlow, 'Visible','on', 'String', sprintf('%.2gx', spq(kk)));
        else
            set(hSlow, 'Visible','off');
        end
        % --- event callouts
        if ~isnan(tSw) && tk >= tSw && tk < tSw + 1.1
            set(hEvent, 'Visible','on', 'String','SLAM', 'Position',[0.255 0.470 0], ...
                'Color', col.accent + (col.hot - col.accent)*max(1-(tk-tSw)/0.30, 0), ...
                'FontSize', 66 + 38*max(1 - (tk-tSw)/0.35, 0));
        elseif kk >= nFly
            set(hEvent, 'Visible','on', 'FontSize',54, 'Color', col.pad, ...
                'Position',[0.500 0.300 0], ...
                'String', sprintf('TOUCHDOWN   %.2f m/s', vtd));
        else
            set(hEvent, 'Visible','off');
        end
    end
end

%% ================= helpers =================
function hudReadout(axH, x0, label)
% HUDREADOUT  Small dim label above a digital readout in the strip.
% INPUTS: axH - HUD axes; x0 - left edge [normalized]; label - char
% OUTPUTS: none
text(axH, x0, 0.80, label, 'Color', [0.520 0.560 0.610], 'FontSize',13, ...
     'FontWeight','bold', 'FontName','Helvetica');
end

function Rot = rot_to(u)
% ROT_TO  Rotation taking +z onto the unit vector u (Rodrigues).
% INPUTS:  u   - unit vector [3x1]
% OUTPUTS: Rot - rotation matrix [3x3] with Rot*[0;0;1] = u
z = [0;0;1];
v = cross(z, u);  s = norm(v);  c = dot(z, u);
if s < 1e-12
    Rot = eye(3) * sign(c + (c==0));
    if c < 0, Rot = diag([1 -1 -1]); end
    return
end
K = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
Rot = eye(3) + K + K*K*((1-c)/s^2);
end

function [Xo, Yo, Zo] = xform(Rot, r0, Xi, Yi, Zi)
% XFORM  Rotate + translate a surface mesh.
% INPUTS:  Rot - rotation [3x3]; r0 - translation [3x1]; Xi,Yi,Zi - mesh
%          coordinate arrays [MxN each]
% OUTPUTS: Xo,Yo,Zo - transformed mesh arrays [MxN each]
sz = size(Xi);
Vv = Rot * [Xi(:).'; Yi(:).'; Zi(:).'] + r0;
Xo = reshape(Vv(1,:), sz);  Yo = reshape(Vv(2,:), sz);  Zo = reshape(Vv(3,:), sz);
end

function hs = build_card(axO, opts, col, sol, fuelKg, vtd, P)
% BUILD_CARD  Title-card text objects (created once, toggled by Visible).
%
% INPUTS:
%   axO    - overlay axes (normalized 0-1 in both directions)
%   opts   - movie options (.title .subtitle .mcRate .mcRuns)
%   col    - palette struct
%   sol    - guidance solution (.tf .mf)
%   fuelKg - propellant used [scalar, kg]
%   vtd    - touchdown speed [scalar, m/s]
%   P      - booster_params (.m0 .Tmax .etaT .pad_radius)
% OUTPUTS:
%   hs     - handle array of every card object [Nx1]
hs = gobjects(0,1);
hs(end+1) = text(axO, 0.5, 0.735, upper(opts.title), 'Color', col.text, ...
    'FontSize',58, 'FontWeight','bold', 'FontName','Helvetica', ...
    'HorizontalAlignment','center');
hs(end+1) = text(axO, 0.5, 0.660, opts.subtitle, 'Color', col.accent, ...
    'FontSize',30, 'FontName','Helvetica', 'HorizontalAlignment','center');
hs(end+1) = plot(axO, [0.30 0.70], [0.615 0.615], '-', 'Color',[col.dim 0.8], ...
    'LineWidth',1.2);

lab = {'FLIGHT TIME', 'PROPELLANT', 'TOUCHDOWN'};
val = {sprintf('%.2f s', sol.tf), sprintf('%.0f kg', fuelKg), sprintf('%.2f m/s', vtd)};
xs  = [0.22 0.50 0.78];
for k = 1:3
    hs(end+1) = text(axO, xs(k), 0.500, val{k}, 'Color', col.text, 'FontSize',44, ...
        'FontWeight','bold', 'FontName','FixedWidth', 'HorizontalAlignment','center'); %#ok<AGROW>
    hs(end+1) = text(axO, xs(k), 0.440, lab{k}, 'Color', col.dim, 'FontSize',17, ...
        'FontWeight','bold', 'FontName','Helvetica', 'HorizontalAlignment','center'); %#ok<AGROW>
end
if ~isempty(opts.mcRate)
    rate = opts.mcRate;  if rate <= 1, rate = 100*rate; end
    if isempty(opts.mcRuns)
        mcs = sprintf('Monte Carlo: %.1f%%', rate);
    else
        mcs = sprintf('Monte Carlo: %.1f%% (%d runs)', rate, opts.mcRuns);
    end
    hs(end+1) = text(axO, 0.5, 0.360, mcs, 'Color', col.pad, 'FontSize',24, ...
        'FontWeight','bold', 'FontName','Helvetica', 'HorizontalAlignment','center');
end
hs(end+1) = text(axO, 0.5, 0.285, sprintf(['%.0f t vehicle  |  %.0f kN engine, ' ...
    'guidance de-rated to %.0f%%  |  %.0f m pad'], P.m0/1000, P.Tmax/1000, ...
    100*P.etaT, P.pad_radius), 'Color', col.dim, 'FontSize',18, ...
    'FontName','Helvetica', 'HorizontalAlignment','center');
hs = hs(:);
end
