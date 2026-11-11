function info = mission_movie(resultsMat, opts)
%% Purpose:
%
%   Renders the public movie of the lunar PNT mission from the SAVED study
%   results (nothing is re-solved): six spacecraft parked in the DRO each
%   fly their polished 70 mN transfer onto a tulip slot, with a yellow-orange
%   thrust plume trailing along the exhaust direction while they thrust,
%   then orbit on station while the south-pole GDOP time series fills the
%   strip at the bottom of the frame.
%
%   Frame (1280 x 720, dark): title and story caption on top; the 3-D scene
%   (Earth-Moon rotating frame, Moon-centred km, orthographic camera) in the
%   middle with counters at the left and a state legend at the right; the
%   GDOP strip (mission day vs south-pole GDOP, satellites on station on a
%   right axis) at the bottom. Every axes has a fixed pixel position and
%   every camera property is manual, so nothing auto-rescales.
%
%   Time map: NON-UNIFORM and smooth. A hold at day 0, then a dwell density
%   (video seconds per mission day) that spends about 5 s on the first DRO
%   revolution, about 29.5 s on the transfers to the last arrival and about
%   30 s on the service window, Gaussian-smoothed in mission day (no speed
%   jumps) with an ease-in and an ease-out, then a hold on the last day.
%   The on-screen day counter and a days-per-second readout keep it honest.
%
%   Camera: fixed view direction (azimuth 12 deg, elevation 18 deg) and
%   fixed target; only the zoom (scene half-height Rv) changes. The first and
%   last frames are the END framing (centre and size from the tulip
%   outline's projected bounding box: DRO and tulip fill the scene). The
%   zoom is built from the data: per map frame, the half-height that keeps
%   all six spacecraft (12 % margin while transferring, tapering to none
%   over the last day before arrival) and their trails inside the scene
%   area, never below the end framing; then, in log space, a 2 s running
%   max each side, a unimodal hull (one zoom-out, one zoom-in; dips are
%   bridged), a raised-cosine smoothing of 6 s half-width, re-inflation
%   until the requirement is covered, and a cosine ramp over the opening
%   hold so frame 1 is exactly the end framing. Every map frame is checked:
%   all six markers at least 4 px inside the scene area, else
%   mission_movie:containment. Markers, trails and plumes are still clipped
%   in screen space to the scene area clipBox (clear of the title,
%   counters, legend and strip): trails by Liang-Barsky, plumes by polygon
%   clipping with colour and alpha carried; a marker outside the box is not
%   drawn, nor its plume. The DRO and tulip outlines (one CR3BP period each,
%   closure asserted) are drawn in every frame at constant brightness, with
%   labels. The scale bar shows the largest round length (1,000 ..
%   100,000 km) whose bar is at most 160 px, with its length recomputed
%   every frame.

%   Moon: pumpkyn.util.moon3D's textured globe (its mesh convention and its
%   Moon_2k.jpg map, loaded by path from the pumpkyn checkout; north +z,
%   prime meridian toward the Earth at -x, no rotation in this frame). A
%   missing texture file falls back to a shaded grey sphere with one line
%   printed.
%
%   Plume: exhaust direction e = +lam_v/|lam_v| (the thrust flown is
%   u = -lam_v/|lam_v|), lam_v linearly interpolated in flight time
%   (clock - tDep) on the de-duplicated flight grid, then normalised. On
%   only for tDep <= t < tArr (a spacecraft with NaN tDep is on station
%   from t = 0, from H.phase). Its length is a fixed number of PIXELS
%   (opts.plumePx; illustrative, a 70 mN thrust has no visual scale); it is
%   foreshortened by the view like any 3-D segment.
%
%   Modes: MP4 (default; frames 1..min(nFrames, maxFrames); written to
%   <name>.partial.mp4 and moved onto outFile only after VideoReader confirms
%   the frame count and size, so an error never truncates the movie; a
%   maxFrames render needs an explicit outFile other than the default,
%   mission_movie:partialOutFile), or PREVIEW
%   PNGs when frameList or dayList is non-empty or preview is true (with
%   both lists empty, preview renders 14 story frames). The layout is
%   designed at 1280 x 720; another frameSize is a resample of that render.
%
%% Inputs:
%
%  resultsMat               char                    lunar_pnt_mission_study.mat
%                                                   ('' = results/ beside
%                                                   this campaign); needs H,
%                                                   P, pole; reads
%                                                   fingerprint (muStar,
%                                                   lStar, tStar; else
%                                                   mission_inputs()) and
%                                                   run_meta when present
%  opts                     struct (optional)
%    .outFile               char                    MP4 path
%                                                   [results/lunar_pnt_mission.mp4]
%    .fps                   scalar [30]             Frame rate [1/s]
%    .durationSec           scalar [70]             Total length [s]
%    .frameSize             [1 x 2] [1280 720]      Output width, height [px]
%    .visible               logical [false]         Show the figure
%    .maxFrames             scalar [inf]            MP4 of the first n frames
%                                                   (explicit outFile required)
%    .frameList             [1 x m] []              Frame indices -> PNGs
%    .dayList               [1 x m] []              Mission days in [0, end]
%                                                   -> PNGs (the zoom of the
%                                                   next map frame, the exact
%                                                   day)
%    .preview               logical [false]         PNG mode; 14 story
%                                                   frames if both lists empty
%    .previewDir            char                    PNG folder
%                                                   [results/movie_preview]
%    .tailDays              scalar [2]              Trail length [d]
%    .plumePx               scalar [62]             On-screen plume length [px]
%    .moonTexture           char                    Moon texture image
%                                                   [~/Desktop/proj7/external/
%                                                   pumpkyn/src/+pumpkyn/+util/
%                                                   Moon_2k.jpg]
%    .probeKm               [m x 3] []              Test probes: magenta
%                                                   markers at these Moon-
%                                                   centred points [km]
%    .quality               scalar [95]             VideoWriter quality
%    .zoomScale             scalar [1]              TEST SEAM: multiplies the
%                                                   zoom envelope before the
%                                                   containment check (< 1
%                                                   forces it too small: the
%                                                   deliberate-failure
%                                                   experiment of
%                                                   mission_movie:containment)
%
%% Outputs:
%
%  info.outFile             char                    MP4 written ('' in PNG mode)
%  info.nFrames             scalar                  Frames rendered
%  info.nFramesMap          scalar                  Frames in the full time map
%  info.frameSize           [1 x 2]                 Width, height [px]
%  info.wallS               scalar                  Wall time [s]
%  info.previewPngs         {m x 1} cell            PNG paths (PNG mode)
%  info.frameIdx            [n x 1]                 Map frame of each render
%  info.frameDay            [n x 1]                 Mission day of each render
%  info.plumeCount          [n x 1]                 Visible plumes per render
%  info.plumeDir            [n x 3 x 6]             Unit exhaust direction of
%                                                   each thrusting spacecraft
%                                                   (NaN when not thrusting)
%  info.plumeDrawn          [n x 6] logical         Plume drawn (thrusting and
%                                                   the spacecraft inside the
%                                                   scene area)
%  info.timeMap             struct                  .tv [F x 1] video s,
%                                                   .day [F x 1], .rate
%                                                   [F x 1] days per second,
%                                                   .budgetS [1 x 5] hold,
%                                                   rev 1, transfers,
%                                                   service, hold
%  info.camera              struct                  .az .el [deg]; .RvEnd end
%                                                   half-height [km]; .Rv .req
%                                                   [F x 1] zoom and its
%                                                   requirement [km]; .hc .vc
%                                                   .target .position
%                                                   .viewAngleEnd .clipBox;
%                                                   .zoomMax .zoomMaxDay
%                                                   .endFramingFromDay
%                                                   .maxRatePctPerS
%                                                   .minMarginPx .minMarginDay
%                                                   .moonDiamPxEnd
%                                                   .moonDiamPxMin .droInside
%                                                   .tulInside; .toPx handle
%                                                   (X [m x 3] km, Rv) -> [m x 2]
%                                                   figure px
%  info.Rv                  [n x 1]                 Half-height per render [km]
%  info.scMarginPx          [n x 1]                 Smallest marker distance to
%                                                   the scene-area edge [px]
%  info.trailAtEdge         [n x 1] logical         A trail was cut at the edge
%  info.plumeScreenDir      [n x 2 x 6]             Drawn plume tip minus
%                                                   nozzle on screen [px]
%  info.probePx             [n x 2 x m]             toPx of the probes [px]
%  info.camPos, .camTarget  [n x 3]                 Camera as set, per render
%  info.camViewAngle        [n x 1]                 likewise [deg]
%  info.camModesManual      [n x 1] logical         All four camera modes manual
%  info.phase               [n x 6]                 Drawn state: 1 parked,
%                                                   2 thrusting, 3 on station
%  info.outlineVisible      [n x 2] logical         DRO, tulip outline visible
%  info.moonTextured        logical                 Moon surface carries the
%                                                   texture (read back)
%  info.moonTexture         char                    Texture file used
%  info.tauDND, .tauTND     scalar                  DRO, tulip periods [ND]
%  info.extentKm            struct                  .dro .transfer .tulip
%                                                   [3 x 2] Moon-centred
%                                                   min/max [km]
%  info.provenance          struct                  .fresh (logical: a fresh,
%                                                   full-knob run), .note
%                                                   (char: '' when fresh,
%                                                   else REUSED / reduced /
%                                                   no run_meta), .units
%                                                   ('fingerprint' |
%                                                   'mission_inputs'); NOT
%                                                   drawn on screen
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  Copyright Coorbital Inc.
%  Revision 2: screen-space clipping, textured Moon, constant outlines,
%  hard-cut captions, distinct spacecraft colours.
%  Revision 3: data-driven unimodal zoom (fixed direction and target),
%  containment assertion, safe partial-file MP4 writing, outlines from CR3BP
%  periods, plume screen direction and probes for tests, review fixes.
%  M. Casey  units from the results fingerprint, provenance in info,
%            zoomScale test seam                              (c) 11/10/2026
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(fileparts(mfilename('fullpath')));        % lunar_pnt_mission/
if nargin == 0
    info = mission_movie('', struct('dayList', [0 8 35]));
    fprintf('mission_movie self-demo: %d preview PNGs in %s (%.1f s)\n', ...
            numel(info.previewPngs), fileparts(info.previewPngs{1}), info.wallS);
    return
end
if nargin < 2 || isempty(opts), opts = struct(); end
if isempty(resultsMat)
    resultsMat = fullfile(here, 'results', 'lunar_pnt_mission_study.mat');
end
tWall = tic;

% ---- options ----
def = struct('outFile',     fullfile(here, 'results', 'lunar_pnt_mission.mp4'), ...
             'fps',         30, ...
             'durationSec', 70, ...
             'frameSize',   [1280 720], ...
             'visible',     false, ...
             'maxFrames',   inf, ...
             'frameList',   [], ...
             'dayList',     [], ...
             'preview',     false, ...
             'previewDir',  fullfile(here, 'results', 'movie_preview'), ...
             'tailDays',    2, ...
             'moonTexture', fullfile(getenv('HOME'), 'Desktop', 'proj7', 'external', 'pumpkyn', ...
                                     'src', '+pumpkyn', '+util', 'Moon_2k.jpg'), ...
             'plumePx',     62, ...
             'probeKm',     zeros(0, 3), ...
             'quality',     95, ...
             'zoomScale',   1);
for f = fieldnames(opts).'
    if ~isfield(def, f{1})
        error('mission_movie:opts', 'unknown option ''%s''', f{1});
    end
    def.(f{1}) = opts.(f{1});
end
o = def;
def0 = struct('outFile', fullfile(here, 'results', 'lunar_pnt_mission.mp4'));
if any(mod(o.frameSize, 2)) || numel(o.frameSize) ~= 2
    error('mission_movie:frameSize', 'frameSize must be two even integers');
end

% ---- layout constants (the CR3BP units muStar, lStar, tStar are read
%      with the data below: the results file's fingerprint, else mission_inputs) ----
rMoonKm = 1737.4;
W       = 1280;                          % layout width  [px]
Hh      = 720;                             % layout height [px]
scnPos  = [0 194 1280 450];                % scene axes [px]
stripPos = [96 62 1088 106];               % GDOP strip axes [px]
halfH   = scnPos(4)/2;                     % scene half-height [px]
vAllow  = 0.90*halfH;                      % content vertical allowance [px]
hAllow  = 370;                             % content horizontal allowance [px]
clipBox = [300 1030 200 640];              % visible scene area x0 x1 y0 y1 [px]:
                                           % clear of title, counters, legend, strip
gLim    = [2 8];                           % GDOP axis (values above clipped)
holdS   = [2.0 3.0];                       % start / end holds [s]
budget0 = [5 29.5 30];                     % rev 1 / transfers / service [s]

% ---- data ----
if ~isfile(resultsMat)
    error('mission_movie:noResults', 'results file not found: %s', resultsMat);
end
L = load(resultsMat, 'H', 'P', 'pole');
for nm = {'H', 'P', 'pole'}
    if ~isfield(L, nm{1}) || isempty(L.(nm{1}))
        error('mission_movie:results', '%s missing or empty in %s', nm{1}, resultsMat);
    end
end
Hs = L.H;  P = L.P;  pole = L.pole;
% CR3BP units from ONE place (the catalog, via the study's fingerprint;
% mission_inputs when the file has none; the DRO closure below checks
% muStar), and the run's provenance for info (not drawn on screen)
[muStar, lStar, tStar, prov] = results_units(resultsMat);
dayNd   = tStar/86400;                     % days per ND time unit
if ~isfield(pole, 'tDays') || ~isfield(pole, 'gdop')
    error('mission_movie:pole', 'pole needs tDays and gdop (%s)', resultsMat);
end
pole.gdop  = double(pole.gdop(:));                 % NaN = undefined (gap); Inf = clipped spike
pole.tDays = double(pole.tDays(:));
if ~iscell(P), P = num2cell(P); end
nSc   = size(Hs.r, 3);
hDay  = double(Hs.t(:))*dayNd;
R     = double(Hs.r);
R(:, 1, :) = R(:, 1, :) - (1 - muStar);
R     = R*lStar;                           % Moon-centred km [N x 3 x 6]
tDepD = double(Hs.tDep(:)).'*dayNd;             % NaN: on station from t = 0
tArrD = double(Hs.tArr(:)).'*dayNd;
flies = isfinite(tDepD) & isfinite(tArrD);
if ~any(flies)
    error('mission_movie:results', 'no spacecraft flies a transfer in %s', resultsMat);
end
dEnd  = hDay(end);
dLast = max(tArrD(flies));
dFirstDep = min(tDepD(flies));

% flights (only for spacecraft that fly): de-duplicated grids, lam_v
FL = cell(1, nSc);
for kk = find(flies)
    if numel(P) < kk || isempty(P{kk}) || ~isfield(P{kk}, 'flight')
        error('mission_movie:results', 'spacecraft %d flies but has no P.flight', kk);
    end
    [tf, ia] = unique(double(P{kk}.flight.t(:)), 'stable');
    if any(diff(tf) <= 0)
        error('mission_movie:flightGrid', 'spacecraft %d flight grid not increasing', kk);
    end
    FL{kk} = struct('t', tf, 'lamv', double(P{kk}.flight.Y(ia, 11:13)), ...
                    'tDep', double(P{kk}.tDep));
end

% outlines: one full period of each orbit, propagated in the CR3BP from a
% state on it (a departure state for the DRO, an arrival state for the
% tulip); closure to 10 km is asserted, so each outline is a closed curve
kF = find(flies);
[droKm, tauD] = periodic_curve(double(P{kF(1)}.rv0(:)), muStar, lStar, [0.5 1.4], 'DRO');
[tulKm, tauT] = periodic_curve(double(P{kF(1)}.rvf(:)), muStar, lStar, [3 8], 'tulip');

% extents
trM = false(size(Hs.phase));  trM(Hs.phase == 2) = true;
ext.dro      = [min(droKm).' max(droKm).'];
ext.transfer = page_extent(R, trM);
ext.tulip    = [min(tulKm).' max(tulKm).'];

% ---- time map ----
nF   = round(o.durationSec*o.fps);
tv   = (0:nF - 1).'/o.fps;
tMov = o.durationSec - sum(holdS);
dRev = tauD*dayNd;
budget = budget0*tMov/sum(budget0);
dd   = linspace(0, dEnd, 40001).';
w    = budget(1)/dRev*(dd < dRev) + budget(2)/(dLast - dRev)*(dd >= dRev & dd < dLast) + ...
       budget(3)/(dEnd - dLast)*(dd >= dLast);
sig  = 0.6;                                % smoothing width [d]
hd   = dd(2) - dd(1);
kx   = (-ceil(4*sig/hd):ceil(4*sig/hd)).'*hd;
ker  = exp(-0.5*(kx/sig).^2);  ker = ker/sum(ker);
nPad = numel(kx);
wp   = [repmat(w(1), nPad, 1); w; repmat(w(end), nPad, 1)];
ws   = conv(wp, ker, 'same');
ws   = ws(nPad + 1:nPad + numel(dd));
ws   = ws.*(1 + 1.2*exp(-(dd/0.25).^2) + 1.2*exp(-((dEnd - dd)/0.5).^2));
sCum = cumtrapz(dd, ws);
scl  = tMov/sCum(end);
sCum = sCum*scl;  ws = ws*scl;
tm   = min(max(tv - holdS(1), 0), tMov);
day  = interp1(sCum, dd, tm);
day(tv >= holdS(1) + tMov) = dEnd;
rate = 1./interp1(dd, ws, day);
rate(tv < holdS(1) | tv >= holdS(1) + tMov) = 0;
TM   = struct('tv', tv, 'day', day, 'rate', rate, ...
              'budgetS', [holdS(1) interp1(dd, sCum, dRev) ...
                          interp1(dd, sCum, dLast) - interp1(dd, sCum, dRev) ...
                          tMov - interp1(dd, sCum, dLast) holdS(2)]);

% ---- camera: fixed direction and target, zoom from the data ----
% The END framing (revision 2's fixed camera: azimuth 12 deg, elevation
% 18 deg, centre and size from the tulip outline's projected bounding box)
% is the first and the last frame. Only the zoom (scene half-height Rv)
% changes: it widens so all six spacecraft stay inside the scene area.
az   = 12;  el = 18;
Tv   = viewmtx(az, el);
eh   = Tv(1, 1:3);  ev = Tv(2, 1:3);
vdir = [cosd(el)*sind(az) -cosd(el)*cosd(az) sind(el)];
ht   = tulKm*eh.';  vt = tulKm*ev.';
hc   = (max(ht) + min(ht))/2;  vc = (max(vt) + min(vt))/2;
RvEnd = max(max(abs(vt - vc))/vAllow, max(abs(ht - hc))/hAllow)*halfH;
tgt  = hc*eh + vc*ev;
cxy  = [scnPos(1) + scnPos(3)/2, scnPos(2) + scnPos(4)/2];      % scene centre [px]
toPx = @(X, Rv) cxy + [X*eh.' - hc, X*ev.' - vc]*(halfH/Rv);   % km -> figure px
inBox = @(S) S(:, 1) >= clipBox(1) & S(:, 1) <= clipBox(2) & ...
             S(:, 2) >= clipBox(3) & S(:, 2) <= clipBox(4);
boxMargin = @(S) min([S(:, 1) - clipBox(1), clipBox(2) - S(:, 1), ...
                      S(:, 2) - clipBox(3), clipBox(4) - S(:, 2)], [], 2);
posF = reshape(interp1(hDay, reshape(R, [], 3*nSc), day), nF, 3, nSc);

% requirement per frame: the half-height that keeps every spacecraft inside
% the scene area, with a 12 % margin while it transfers (tapering to none
% over the last day before its arrival, so the zoom can come home), never
% below the end framing
lim  = [clipBox(1) - cxy(1), clipBox(2) - cxy(1), clipBox(3) - cxy(2), clipBox(4) - cxy(2)];
hS   = squeeze(sum(posF.*reshape(eh, 1, 3), 2)) - hc;          % [nF x 6] km
vS   = squeeze(sum(posF.*reshape(ev, 1, 3), 2)) - vc;
mg   = ones(nF, nSc);
for kk = find(flies)
    tr = day >= tDepD(kk) & day < tArrD(kk);
    mg(tr, kk) = 1 + 0.12*min(1, (tArrD(kk) - day(tr))/1.0);
end
need = max([max(max(hS, 0).*mg/lim(2), [], 2), max(max(-hS, 0).*mg/(-lim(1)), [], 2), ...
            max(max(vS, 0).*mg/lim(4), [], 2), max(max(-vS, 0).*mg/(-lim(3)), [], 2)], [], 2)*halfH;
% trails (tailDays behind each spacecraft, sampled every hour) must fit too,
% without extra margin, so no trail is cut at the scene edge
tStep = max(1, round(3600/(median(diff(hDay))*86400)));
for ff = 1:nF
    tw = find(hDay > day(ff) - o.tailDays & hDay < day(ff));
    if isempty(tw), continue; end
    tw = tw(1:tStep:end);
    q  = reshape(permute(R(tw, :, :), [1 3 2]), [], 3);
    hq = q*eh.' - hc;  vq = q*ev.' - vc;
    need(ff) = max([need(ff), max(hq)/lim(2)*halfH, max(-hq)/(-lim(1))*halfH, ...
                    max(vq)/lim(4)*halfH, max(-vq)/(-lim(3))*halfH]);
end
req  = max(need, RvEnd);

% smooth envelope in log space: running max (2 s each side), unimodal hull
% (bridges any dip: one zoom-out, one zoom-in), raised-cosine smoothing
% (6 s half-width), then inflate until it covers the requirement again
LEnd = log(RvEnd);
Lq   = log(req);
wMax = round(2*o.fps);  hw = round(6*o.fps);
Lz   = unimodal_hull(movmax(Lq, [wMax wMax]));
Lz   = hann_smooth(Lz, hw);
for it = 1:10
    gap = Lq - Lz;
    if all(gap <= 0), break; end
    Lz = unimodal_hull(Lz + hann_smooth(movmax(max(gap, 0) + 0.002*(gap > 0), [wMax wMax]), hw));
end
% the first frame IS the end framing: the excess is ramped in over the 2 s
% hold (the requirement there is the end framing itself)
ramp = (1 - cos(pi*min(tv/holdS(1), 1)))/2;
Lz   = LEnd + (Lz - LEnd).*ramp;
Lz(abs(Lz - LEnd) < 1e-12) = LEnd;
Lz   = max(Lz, LEnd);
RvF  = exp(Lz);
RvF(Lz == LEnd) = RvEnd;
RvF  = RvF*o.zoomScale;                   % TEST SEAM (1 in every real render)
% containment over EVERY map frame: all six markers inside by >= 4 px
mPx = zeros(nF, 1);
for ff = 1:nF
    mPx(ff) = min(boxMargin(toPx(reshape(posF(ff, :, :), 3, nSc).', RvF(ff))));
end
if min(mPx) < 4
    [mm, kw] = min(mPx);
    error('mission_movie:containment', 'a spacecraft is %.1f px from the scene edge at day %.3f', ...
          mm, day(kw));
end
backAt = find(RvF > RvEnd, 1, 'last') + 1;
[RvMax, kMax] = max(RvF);
CAM  = struct('az', az, 'el', el, 'RvEnd', RvEnd, 'Rv', RvF, 'req', req, 'hc', hc, 'vc', vc, ...
              'target', tgt, 'position', tgt + 1e6*vdir, ...
              'viewAngleEnd', 2*atand(RvEnd/1e6), 'clipBox', clipBox, ...
              'zoomMax', RvMax/RvEnd, 'zoomMaxDay', day(kMax), ...
              'endFramingFromDay', day(min(backAt, nF)), ...
              'maxRatePctPerS', 100*max(abs(diff(log(RvF))))*o.fps, ...
              'minMarginPx', min(mPx), 'minMarginDay', day(find(mPx == min(mPx), 1)), ...
              'moonDiamPxEnd', 2*rMoonKm*halfH/RvEnd, 'moonDiamPxMin', 2*rMoonKm*halfH/RvMax, ...
              'droInside', all(inBox(toPx(droKm, RvEnd))), 'tulInside', all(inBox(toPx(tulKm, RvEnd))), ...
              'toPx', toPx);

% ---- what to render ----
pngMode = ~isempty(o.frameList) || ~isempty(o.dayList) || o.preview;
if pngMode
    fIdx = o.frameList(:);  fDay = day(min(max(fIdx, 1), nF));
    if any(fIdx < 1 | fIdx > nF | fIdx ~= round(fIdx))
        error('mission_movie:frameList', 'frame indices must be integers in 1..%d', nF);
    end
    if o.preview && isempty(o.frameList) && isempty(o.dayList)
        o.dayList = [0 1 2.2 3.5 5 8 10 12 14 16 17.4 dLast 33 dEnd];
    end
    dq  = o.dayList(:);
    if any(~isfinite(dq) | dq < 0 | dq > dEnd + 1e-9)
        error('mission_movie:dayList', 'dayList entries must lie in [0, %.4f] (mission days)', dEnd);
    end
    kfs = nF*ones(numel(dq), 1);
    for kq = 1:numel(dq)
        kf = find(day >= dq(kq) - 1e-12, 1);
        if ~isempty(kf), kfs(kq) = kf; end
    end
    fIdx = [fIdx; kfs];
    fDay = [fDay; dq];
    if ~isfolder(o.previewDir), mkdir(o.previewDir); end
else
    nR   = min(nF, o.maxFrames);
    fIdx = (1:nR).';  fDay = day(fIdx);
    % a partial render must never replace the default movie
    if nR < nF && (~isfield(opts, 'outFile') || strcmp(full_path(o.outFile), full_path(def0.outFile)))
        error('mission_movie:partialOutFile', ...
              'maxFrames = %d < %d frames needs an explicit outFile other than the default %s', ...
              nR, nF, def0.outFile);
    end
end
nR = numel(fIdx);

% ---- figure ----
fig = figure('Visible', onoff(o.visible), 'Units', 'pixels', 'Position', [20 20 W Hh], ...
             'Color', [0 0 0], 'MenuBar', 'none', 'ToolBar', 'none', ...
             'Resize', 'off');
cleanFig = onCleanup(@() close_if(fig));
% outlines lie inside the scene area at the end framing, hence at every
% wider zoom: drawn unclipped
G = build_graphics(fig, scnPos, stripPos, droKm, tulKm, rMoonKm, pole, tArrD, dEnd, gLim, nSc, ...
                   o.moonTexture);
set(G.ax, 'CameraTarget', CAM.target, 'CameraPosition', CAM.position, ...
          'CameraViewAngle', CAM.viewAngleEnd);
pd = droKm*[eh; ev].';
[~, kd] = max(pd(:, 1) - 0.4*pd(:, 2));
pt = tulKm*[eh; ev].';
[~, kt] = max(pt(:, 1) - 0.5*pt(:, 2));
nProbe = size(o.probeKm, 1);
if nProbe > 0
    G.probe = plot3(G.ax, o.probeKm(:, 1), o.probeKm(:, 2), o.probeKm(:, 3), 's', ...
                    'MarkerSize', 9, 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', [1 0 1]);
end

% MP4: written to <name>.partial.mp4 and moved onto outFile only after a
% complete render that VideoReader confirms
vw = [];  partFile = '';
if ~pngMode
    [pDir, pName] = fileparts(o.outFile);
    partFile = fullfile(pDir, [pName '.partial.mp4']);
    vw = VideoWriter(partFile, 'MPEG-4');
    vw.FrameRate = o.fps;
    vw.Quality   = o.quality;
    open(vw);
end
cleanVw = onCleanup(@() close_vw(vw, partFile));

% caption schedule (mission-day boundaries -> video times)
trDays = tArrD - tDepD;
caps = {'Six spacecraft parked in a distant retrograde orbit', ...
        sprintf('Continuous-thrust transfers (%d-%d days each)', floor(min(trDays)), ceil(max(trDays))), ...
        sprintf('Constellation complete: day %.1f', dLast), ...
        sprintf('South-pole navigation geometry over %.0f days', dEnd - dLast)};
capD = [-inf dFirstDep dLast dLast + 4 inf];        % hard cuts at these days (the
                                                     % parked caption ends at the FIRST departure)

info = struct('outFile', '', 'nFrames', nR, 'nFramesMap', nF, 'frameSize', o.frameSize, ...
              'wallS', NaN, 'previewPngs', {cell(0, 1)}, 'frameIdx', fIdx, 'frameDay', fDay, ...
              'plumeCount', zeros(nR, 1), 'plumeDir', nan(nR, 3, nSc), 'timeMap', TM, ...
              'plumeDrawn', false(nR, nSc), 'plumeScreenDir', nan(nR, 2, nSc), ...
              'camera', CAM, 'extentKm', ext, 'Rv', nan(nR, 1), 'scMarginPx', nan(nR, 1), ...
              'trailAtEdge', false(nR, 1), 'probePx', nan(nR, 2, nProbe), ...
              'camPos', nan(nR, 3), 'camTarget', nan(nR, 3), 'camViewAngle', nan(nR, 1), ...
              'camModesManual', false(nR, 1), 'tauDND', tauD, 'tauTND', tauT, ...
              'phase', zeros(nR, nSc), ...
              'outlineVisible', false(nR, 2), 'moonTextured', G.moonTextured, ...
              'moonTexture', G.moonTexFile, 'provenance', prov);
gFirst = pole.tDays(find(isfinite(pole.gdop), 1));
if isempty(gFirst), gFirst = inf; end

for rr = 1:nR
    ff = fIdx(rr);  d0 = fDay(rr);
    tc = d0/dayNd;                                    % clock [ND]
    if pngMode && rr > numel(o.frameList)
        pos = reshape(interp1(hDay, reshape(R, [], 3*nSc), d0), 3, nSc).';
    else
        pos = reshape(posF(ff, :, :), 3, nSc).';
    end
    % camera: fixed direction and target, zoom of this map frame
    Rv = RvF(ff);
    kmPerPx = Rv/halfH;
    toPxR = @(X) toPx(X, Rv);
    set(G.ax, 'CameraViewAngle', 2*atand(Rv/1e6));
    info.Rv(rr) = Rv;
    info.scMarginPx(rr) = min(boxMargin(toPxR(pos)));
    set(G.droLbl, 'Position', droKm(kd, :) + 14*kmPerPx*eh);
    set(G.tulLbl, 'Position', tulKm(kt, :) + 14*kmPerPx*eh);
    update_overlay(G, Tv, kmPerPx, scnPos);
    if nProbe > 0
        info.probePx(rr, :, :) = reshape(toPxR(o.probeKm).', 1, 2, nProbe);
    end
    % spacecraft
    nOn = 0;  nPl = 0;
    iH = max(1, find(hDay <= d0 + 1e-9, 1, 'last'));
    for kk = 1:nSc
        if flies(kk)
            ph = 1 + (d0 >= tDepD(kk)) + (d0 >= tArrD(kk));
        else
            ph = double(Hs.phase(iH, kk));            % e.g. on station from t = 0
        end
        info.phase(rr, kk) = ph;
        p  = pos(kk, :);
        pIn = inBox(toPxR(p));
        set(G.sc(kk), 'XData', p(1), 'YData', p(2), 'ZData', p(3), 'Visible', onoff(pIn));
        switch ph
            case 1
                set(G.sc(kk), 'MarkerFaceColor', 'none', 'MarkerSize', 7, 'LineWidth', 1.6);
                set(G.ring(kk), 'Visible', 'off');
            case 2
                set(G.sc(kk), 'MarkerFaceColor', G.col(kk, :), 'MarkerSize', 8, 'LineWidth', 1);
                set(G.ring(kk), 'Visible', 'off');
            otherwise
                set(G.sc(kk), 'MarkerFaceColor', G.col(kk, :), 'MarkerSize', 7, 'LineWidth', 1);
                set(G.ring(kk), 'Visible', onoff(pIn), 'XData', p(1), 'YData', p(2), 'ZData', p(3));
                nOn = nOn + 1;
        end
        % trail
        tw = find(hDay > d0 - o.tailDays & hDay < d0);
        if ~isempty(tw), tw = tw(mod(tw - tw(end), 2) == 0); end
        tp = [R(tw, :, kk); p];
        nt = size(tp, 1);
        al = 0.85*((1:nt).'/nt).^1.6;
        info.trailAtEdge(rr) = info.trailAtEdge(rr) || ~all(inBox(toPxR(tp)));
        [tv3, ta] = clip_poly(tp, al, toPxR, clipBox);  % clipped to the scene area
        if size(tv3, 1) < 2, tv3 = nan(2, 3);  ta = [0; 0]; end
        set(G.trail(kk), 'Vertices', tv3, 'Faces', 1:size(tv3, 1), 'FaceVertexAlphaData', ta);
        % plume
        if ph == 2
            lam = interp1(FL{kk}.t, FL{kk}.lamv, min(max(tc - FL{kk}.tDep, 0), FL{kk}.t(end)), 'linear');
            e   = lam/sqrt(sum(lam.^2));
            info.plumeDir(rr, :, kk) = e;
            nv  = cross(e, vdir);                     % across the plume, in the screen plane
            if sqrt(sum(nv.^2)) < 0.05, nv = eh; end
            nv  = nv/sqrt(sum(nv.^2));
            c0  = p + 4*kmPerPx*e;                    % nozzle, just behind the marker
            Lkm = o.plumePx*kmPerPx;
            S   = G.plumeShape;
            Vg  = plume_poly(c0, e, nv, S.glowL*Lkm, S.waist, S.glowW*kmPerPx);
            Vc  = plume_poly(c0, e, nv, Lkm, S.waist, S.coreW*kmPerPx);
            % clipped to the scene area polygon-wise, colour and alpha carried
            sp  = toPxR([Vc(3, :); (Vc(1, :) + Vc(5, :))/2]);   % tip, nozzle [px]
            info.plumeScreenDir(rr, :, kk) = sp(1, :) - sp(2, :);
            [Vg, Cg, Ag] = clip_polygon(Vg, G.glowC, G.glowA, toPxR, clipBox);
            [Vc, Cc, Ac] = clip_polygon(Vc, G.coreC, G.coreA, toPxR, clipBox);
            if ~pIn                                   % no plume without its spacecraft
                Vg = zeros(0, 3);  Vc = zeros(0, 3);
            end
            drawn = size(Vc, 1) >= 3;
            set_plume(G.plume{kk}(1), Vg, Cg, Ag);
            set_plume(G.plume{kk}(2), Vc, Cc, Ac);
            info.plumeDrawn(rr, kk) = drawn;
            nPl = nPl + drawn;
        else
            set(G.plume{kk}, 'Visible', 'off');
        end
    end
    info.plumeCount(rr) = nPl;
    info.outlineVisible(rr, :) = [strcmp(G.droLine.Visible, 'on') && any(isfinite(G.droLine.XData)), ...
                                  strcmp(G.tulLine.Visible, 'on') && any(isfinite(G.tulLine.XData))];
    info.camPos(rr, :) = G.ax.CameraPosition;
    info.camTarget(rr, :) = G.ax.CameraTarget;
    info.camViewAngle(rr) = G.ax.CameraViewAngle;
    info.camModesManual(rr) = all(strcmp({G.ax.CameraPositionMode, G.ax.CameraTargetMode, ...
                                          G.ax.CameraViewAngleMode, G.ax.CameraUpVectorMode}, 'manual'));
    % text
    gNow = NaN;
    kp = find(pole.tDays <= d0 + 1e-9, 1, 'last');
    if ~isempty(kp), gNow = pole.gdop(kp); end
    set(G.txtDay, 'String', sprintf('Day %4.1f', d0));
    set(G.txtOn,  'String', sprintf('On station %d / %d', nOn, nSc));
    if isinf(gNow)
        set(G.txtG, 'String', sprintf('South-pole GDOP %5s', 'inf'), 'Color', [1 0.45 0.45]);
    elseif isfinite(gNow)
        if gNow < 6, gc = [0.45 0.95 0.55]; elseif gNow < 10, gc = [1 0.78 0.35]; else, gc = [1 0.45 0.45]; end
        set(G.txtG, 'String', sprintf('South-pole GDOP %5.2f', gNow), 'Color', gc);
    else
        set(G.txtG, 'String', sprintf('South-pole GDOP %5s', '--'), 'Color', [0.65 0.65 0.7]);
    end
    rNow = rate(ff);
    if rNow > 0
        set(G.txtRate, 'String', sprintf('%4.2f days per second', rNow));
    else
        set(G.txtRate, 'String', '');
    end
    kc = find(capD(2:end) > d0, 1);
    set(G.cap, 'String', caps{kc});
    % strip
    kg = pole.tDays <= d0 + 1e-9;
    set(G.gLine, 'XData', pole.tDays(kg), 'YData', clip_top(pole.gdop(kg), gLim(2)));
    kx2 = kg & pole.gdop > gLim(2);
    set(G.gClip, 'XData', pole.tDays(kx2), 'YData', gLim(2)*ones(nnz(kx2), 1));
    [sx, sy] = on_station_steps(tArrD, d0);
    set(G.nLine, 'XData', sx, 'YData', sy);
    set(G.cursor, 'XData', [d0 d0]);
    for kk = 1:nSc
        set(G.arrTick(kk), 'Visible', onoff(d0 >= tArrD(kk)));
    end
    set(G.needTxt, 'Visible', onoff(d0 < gFirst));
    drawnow;
    F = getframe(fig);
    img = F.cdata;
    if size(img, 1) ~= o.frameSize(2) || size(img, 2) ~= o.frameSize(1)
        img = imresize(img, [o.frameSize(2) o.frameSize(1)]);
    end
    if pngMode
        pf = fullfile(o.previewDir, sprintf('frame_%04d_day%05.2f.png', ff, d0));
        imwrite(img, pf);
        info.previewPngs{end + 1, 1} = pf;
    else
        writeVideo(vw, img);
    end
end
if ~pngMode
    close(vw);
    vr = VideoReader(partFile);
    okV = vr.NumFrames == nR && vr.Width == o.frameSize(1) && vr.Height == o.frameSize(2);
    clear vr
    if ~okV
        error('mission_movie:verify', 'written video failed its check: %s', partFile);
    end
    movefile(partFile, o.outFile, 'f');               % only now is outFile replaced
    info.outFile = o.outFile;
end
clear cleanVw cleanFig
info.wallS = toc(tWall);
end

%% ------------------------------------------------------------------------
function G = build_graphics(fig, scnPos, stripPos, droKm, tulKm, rMoonKm, pole, tArrD, dEnd, ...
                           gLim, nSc, texFile)
% BUILD_GRAPHICS  Create every graphics object once; the loop only updates.
%   INPUTS:  fig figure; scnPos, stripPos [1 x 4] px; droKm [m x 3] and
%            tulKm [q x 3] outline vertices (NaN rows break the line);
%            rMoonKm scalar; pole struct; tArrD [1 x 6] d; dEnd scalar d;
%            gLim [1 x 2]; nSc scalar; texFile char Moon texture image.
%   OUTPUTS: G struct of handles and styles.
fg   = [0.93 0.94 0.97];
fg2  = [0.70 0.73 0.80];
mono = 'Menlo';
if ~any(strcmp(listfonts, mono)), mono = 'FixedWidth'; end
sans = 'Helvetica';
G.col = [0.40 0.80 0.93;      % cyan
         0.93 0.30 0.33;      % red
         0.30 0.78 0.40;      % green
         0.36 0.58 1.00;      % blue
         0.72 0.50 1.00;      % violet
         1.00 0.55 0.80];     % pink
% scene
ax = axes(fig, 'Units', 'pixels', 'Position', scnPos, 'Color', 'none', ...
          'XLim', [-1.5e5 1.5e5], 'YLim', [-1.5e5 1.5e5], 'ZLim', [-1.5e5 1.5e5], ...
          'DataAspectRatio', [1 1 1], 'Projection', 'orthographic', 'Clipping', 'off', ...
          'CameraUpVector', [0 0 1], 'Visible', 'off', 'SortMethod', 'depth');
hold(ax, 'on');
% Moon: pumpkyn.util.moon3D's construction (ellipsoid mesh with x, y, z all
% negated, then the equirectangular map as a texture), which puts north
% along +z and the prime meridian toward -x (the Earth) in the rotating
% frame; the Moon is tidally locked, so it does not turn in this frame
[xs, ys, zs] = ellipsoid(0, 0, 0, rMoonKm, rMoonKm, rMoonKm, 100);
G.moon = surf(ax, -xs, -ys, -zs, 'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
              'AmbientStrength', 0.65, 'DiffuseStrength', 0.75, 'SpecularStrength', 0, ...
              'BackFaceLighting', 'unlit');
G.moonTexFile  = texFile;
G.moonTextured = isfile(texFile);
if G.moonTextured
    tex = imread(texFile);
    if size(tex, 3) == 1, tex = repmat(tex, 1, 1, 3); end
    set(G.moon, 'FaceColor', 'texturemap', 'CData', tex);
else
    set(G.moon, 'FaceColor', [0.66 0.66 0.68]);
    fprintf('mission_movie: Moon texture not found (%s); drawing a shaded grey sphere\n', texFile);
end
G.moonTextured = ischar(G.moon.FaceColor) && strcmp(G.moon.FaceColor, 'texturemap') && ...
                 size(G.moon.CData, 3) == 3;              % read back from the surface
light(ax, 'Style', 'infinite', 'Position', [0.35 -1 0.45]);
G.droRGB = [0.45 0.62 0.95];
G.tulRGB = [0.86 0.86 0.92];
G.droLine = plot3(ax, droKm(:, 1), droKm(:, 2), droKm(:, 3), '-', 'Color', [G.droRGB 0.70], 'LineWidth', 1.1);
G.tulLine = plot3(ax, tulKm(:, 1), tulKm(:, 2), tulKm(:, 3), '-', 'Color', [G.tulRGB 0.55], 'LineWidth', 1.0);
G.droLbl = text(ax, 0, 0, 0, 'DRO', 'Color', [0.55 0.70 1.0], 'FontName', sans, 'FontSize', 12, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
G.tulLbl = text(ax, 0, 0, 0, 'tulip orbit', 'Color', [0.85 0.85 0.9], 'FontName', sans, ...
                'FontSize', 12, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
% plume: two tapered patches (soft glow under a bright core), vertex order
% nozzle+, waist+, tip, waist-, nozzle-; colour and alpha interpolated
G.plumeShape = struct('waist', 0.30, ...                      % waist at 30 % of the length
                      'coreW', [9 6 0], 'glowW', [20 14 0], ...   % widths [px] nozzle, waist, tip
                      'glowL', 0.80);                         % glow length / core length
coreC = [1.00 1.00 0.80; 1.00 0.86 0.30; 1.00 0.86 0.30; 1.00 0.86 0.30; 1.00 1.00 0.80];
coreC(3, :) = [1.00 0.42 0.05];
coreA = [1; 0.85; 0; 0.85; 1];
glowC = repmat([1.00 0.78 0.25], 5, 1);
glowA = [0.30; 0.18; 0; 0.18; 0.30];
G.coreC = coreC;  G.coreA = coreA;  G.glowC = glowC;  G.glowA = glowA;
G.plume = cell(1, nSc);  G.trail = gobjects(1, nSc);
G.sc = gobjects(1, nSc);  G.ring = gobjects(1, nSc);
for kk = 1:nSc
    G.trail(kk) = patch(ax, 'Vertices', nan(2, 3), 'Faces', [1 2], 'FaceColor', 'none', ...
                        'EdgeColor', G.col(kk, :), 'EdgeAlpha', 'interp', ...
                        'FaceVertexAlphaData', [0; 0], 'AlphaDataMapping', 'none', 'LineWidth', 1.5);
end
for kk = 1:nSc
    G.plume{kk}(1) = patch(ax, 'Vertices', zeros(5, 3), 'Faces', 1:5, 'FaceColor', 'interp', ...
                           'FaceVertexCData', glowC, 'FaceAlpha', 'interp', ...
                           'FaceVertexAlphaData', glowA, 'AlphaDataMapping', 'none', ...
                           'EdgeColor', 'none', 'FaceLighting', 'none', 'Visible', 'off');
    G.plume{kk}(2) = patch(ax, 'Vertices', zeros(5, 3), 'Faces', 1:5, 'FaceColor', 'interp', ...
                           'FaceVertexCData', coreC, 'FaceAlpha', 'interp', ...
                           'FaceVertexAlphaData', coreA, 'AlphaDataMapping', 'none', ...
                           'EdgeColor', 'none', 'FaceLighting', 'none', 'Visible', 'off');
end
for kk = 1:nSc
    G.ring(kk) = plot3(ax, 0, 0, 0, 'o', 'MarkerSize', 14, 'LineWidth', 0.8, ...
                       'Color', G.col(kk, :), 'MarkerFaceColor', 'none', 'Visible', 'off');
    G.sc(kk) = plot3(ax, 0, 0, 0, 'o', 'Color', G.col(kk, :), 'MarkerSize', 7, 'LineWidth', 1.6);
end
G.ax = ax;

% GDOP strip: left axes GDOP, right axes satellites on station
axG = axes(fig, 'Units', 'pixels', 'Position', stripPos, 'Color', [0.05 0.06 0.09], ...
           'XLim', [0 dEnd], 'YLim', gLim, 'XColor', fg2, 'YColor', fg2, 'FontName', sans, ...
           'FontSize', 11, 'Box', 'on', 'Layer', 'top', 'XTick', 0:5:45, 'YTick', 2:2:8, ...
           'TickDir', 'out', 'GridColor', [1 1 1], 'GridAlpha', 0.10, 'XGrid', 'on', 'YGrid', 'on');
hold(axG, 'on');
xlabel(axG, 'mission day', 'Color', fg2, 'FontSize', 12);
ylabel(axG, 'south-pole GDOP', 'Color', fg2, 'FontSize', 12);
for kk = 1:nSc
    G.arrTick(kk) = plot(axG, [tArrD(kk) tArrD(kk)], gLim(1) + [0 0.16]*diff(gLim), '-', ...
                         'Color', G.col(kk, :), 'LineWidth', 2, 'Visible', 'off');
end
h6 = plot(axG, [0 dEnd], [6 6], '--', 'Color', [1 0.45 0.45 0.85], 'LineWidth', 1.1);
text(axG, dEnd - 0.3, 6, 'GDOP = 6', 'Color', [1 0.55 0.55], 'FontSize', 10, 'FontName', sans, ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
G.gLine = plot(axG, NaN, NaN, '-', 'Color', [0.97 0.97 1.0], 'LineWidth', 1.3);
G.gClip = plot(axG, NaN, NaN, '^', 'Color', [1 0.45 0.45], 'MarkerFaceColor', [1 0.45 0.45], ...
               'MarkerSize', 4, 'LineStyle', 'none');
G.cursor = plot(axG, [0 0], gLim, '-', 'Color', [1 1 1 0.55], 'LineWidth', 1);
hN = plot(axG, NaN, NaN, '-', 'Color', [0.55 0.72 1.0], 'LineWidth', 1.6);
G.needTxt = text(axG, 9, 3.9, 'GDOP needs 4 satellites on station', 'Color', fg2, ...
                 'FontName', sans, 'FontSize', 12, 'FontAngle', 'italic', ...
                 'HorizontalAlignment', 'center');
lg = legend(axG, [G.gLine hN h6 G.gClip], {'south-pole GDOP', 'satellites on station', ...
            'GDOP = 6', sprintf('GDOP > %d (clipped)', gLim(2))}, 'NumColumns', 2, ...
            'TextColor', fg, 'Color', [0.08 0.09 0.13], 'EdgeColor', [0.3 0.32 0.38], ...
            'FontSize', 10, 'FontName', sans, 'Location', 'northwest', 'AutoUpdate', 'off');
lg.ItemTokenSize = [18 10];
axN = axes(fig, 'Units', 'pixels', 'Position', stripPos, 'Color', 'none', 'XLim', [0 dEnd], ...
           'YLim', [0 6.6], 'YTick', 0:2:6, 'YAxisLocation', 'right', 'XColor', 'none', ...
           'YColor', [0.55 0.72 1.0], 'FontName', sans, 'FontSize', 11, 'TickDir', 'out', ...
           'Box', 'off');
hold(axN, 'on');
ylabel(axN, 'satellites on station', 'Color', [0.55 0.72 1.0], 'FontSize', 12);
G.nLine = plot(axN, NaN, NaN, '-', 'Color', [0.55 0.72 1.0], 'LineWidth', 1.6);
if ~isfield(pole, 'tDays') || ~isfield(pole, 'gdop')
    error('mission_movie:pole', 'pole needs tDays and gdop');
end

% overlay in figure pixels
ov = axes(fig, 'Units', 'pixels', 'Position', [0 0 1280 720], 'XLim', [0 1280], 'YLim', [0 720], ...
          'Visible', 'off', 'Color', 'none');
hold(ov, 'on');
text(ov, 640, 700, 'Lunar PNT constellation deployment: DRO to tulip, 70 mN low thrust', ...
     'Color', fg, 'FontName', sans, 'FontSize', 20, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
G.cap = text(ov, 640, 668, '', 'Color', fg2, 'FontName', sans, 'FontSize', 15, ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
G.txtDay  = text(ov, 28, 612, '', 'Color', fg, 'FontName', mono, 'FontSize', 20, 'FontWeight', 'bold');
G.txtOn   = text(ov, 28, 580, '', 'Color', fg, 'FontName', mono, 'FontSize', 15);
G.txtG    = text(ov, 28, 554, '', 'Color', fg, 'FontName', mono, 'FontSize', 15);
G.txtRate = text(ov, 28, 528, '', 'Color', [0.55 0.58 0.65], 'FontName', mono, 'FontSize', 11);
% state legend (right)
x0 = 1062;  y0 = 612;  dy = 25;
text(ov, x0 - 8, y0 + 4, 'spacecraft state', 'Color', fg2, 'FontName', sans, 'FontSize', 11);
yy = y0 - dy;
plot(ov, x0 + 8, yy, 'o', 'Color', fg, 'MarkerSize', 7, 'LineWidth', 1.6);
text(ov, x0 + 26, yy, 'parked in DRO', 'Color', fg, 'FontName', sans, 'FontSize', 12);
yy = yy - dy;
plot(ov, [x0 - 14 x0 + 2], [yy yy], '-', 'Color', [1 0.66 0.15 0.7], 'LineWidth', 3);
plot(ov, [x0 - 6 x0 + 4], [yy yy], '-', 'Color', [1 0.95 0.6], 'LineWidth', 5);
plot(ov, x0 + 8, yy, 'o', 'Color', fg, 'MarkerFaceColor', fg, 'MarkerSize', 8);
text(ov, x0 + 26, yy, 'thrusting (plume)', 'Color', fg, 'FontName', sans, 'FontSize', 12);
yy = yy - dy;
plot(ov, x0 + 8, yy, 'o', 'Color', fg, 'MarkerSize', 14, 'LineWidth', 0.8);
plot(ov, x0 + 8, yy, 'o', 'Color', fg, 'MarkerFaceColor', fg, 'MarkerSize', 7);
text(ov, x0 + 26, yy, 'on station', 'Color', fg, 'FontName', sans, 'FontSize', 12);
yy = yy - dy - 6;
text(ov, x0 - 8, yy, 'spacecraft', 'Color', fg2, 'FontName', sans, 'FontSize', 11);
for kk = 1:nSc
    text(ov, x0 + 84 + 20*(kk - 1), yy, sprintf('%d', kk), 'Color', G.col(kk, :), ...
         'FontName', sans, 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
% scale bar, frame label, Earth arrow
G.sbLine = plot(ov, [0 1], [212 212], '-', 'Color', fg2, 'LineWidth', 2);
G.sbEnds = plot(ov, [0 0 NaN 1 1], [207 217 NaN 207 217], '-', 'Color', fg2, 'LineWidth', 1.2);
G.sbTxt  = text(ov, 28, 228, '', 'Color', fg2, 'FontName', mono, 'FontSize', 11);
text(ov, 1252, 212, 'Earth-Moon rotating frame', 'Color', [0.55 0.58 0.65], 'FontName', sans, ...
     'FontSize', 11, 'HorizontalAlignment', 'right');
G.eaLine = plot(ov, [0 1], [0 1], '-', 'Color', [0.55 0.75 1.0], 'LineWidth', 1.6);
G.eaHead = patch(ov, 'XData', [0 0 0], 'YData', [0 0 0], 'FaceColor', [0.55 0.75 1.0], 'EdgeColor', 'none');
G.eaTxt  = text(ov, 0, 0, 'to Earth', 'Color', [0.55 0.75 1.0], 'FontName', sans, 'FontSize', 11, ...
                'HorizontalAlignment', 'center');
G.ov = ov;
end

%% ------------------------------------------------------------------------
function update_overlay(G, Tv, kmPerPx, scnPos)
% UPDATE_OVERLAY  Earth-direction arrow at the scene edge and the scale bar.
%   INPUTS:  G handles; Tv [4 x 4] viewmtx; kmPerPx scalar; scnPos [1 x 4] px.
%   OUTPUTS: none (updates G).
% Earth direction (-x of the rotating frame) as seen on screen, drawn as a
% fixed-anchor pointer at the left edge of the scene (Earth is 389,703 km
% away, far outside the frame)
dE = (Tv(1:2, 1:3)*[-1; 0; 0]).';
if sqrt(sum(dE.^2)) < 1e-6, dE = [-1 0]; end
dE = dE/sqrt(sum(dE.^2));
anc  = [scnPos(1) + 78, scnPos(2) + 120];
tip  = anc + 30*dE;
tail = anc - 30*dE;
nrm  = [-dE(2) dE(1)];
base = tip - 12*dE;
set(G.eaLine, 'XData', [tail(1) base(1)], 'YData', [tail(2) base(2)]);
set(G.eaHead, 'XData', [tip(1) base(1) + 5*nrm(1) base(1) - 5*nrm(1)], ...
              'YData', [tip(2) base(2) + 5*nrm(2) base(2) - 5*nrm(2)]);
set(G.eaTxt, 'Position', [anc(1), anc(2) - 34, 0]);
% scale bar: the largest round length whose bar is at most 160 px (under a
% quarter of the scene area's width); its on-screen length follows the zoom
nice = [1000 2000 5000 10000 20000 50000 100000];
kb  = max([1, find(nice/kmPerPx <= 160, 1, 'last')]);
Lpx = nice(kb)/kmPerPx;
x0 = 28;
set(G.sbLine, 'XData', [x0 x0 + Lpx]);
set(G.sbEnds, 'XData', [x0 x0 NaN x0 + Lpx x0 + Lpx]);
set(G.sbTxt, 'String', sprintf('%-10s', [thousands(nice(kb)) ' km']));   % left-anchored, fixed width
end

%% ------------------------------------------------------------------------
function [km, tau] = periodic_curve(rv0, muStar, lStar, tWin, name)
% PERIODIC_CURVE  One period of a periodic CR3BP orbit propagated from a
%   state on it. The period is the closest full-state return to rv0 inside
%   tWin; closure (position back within 10 km) is asserted, so the outline is
%   a closed curve, which also checks muStar against the saved states.
%   INPUTS:  rv0 [6 x 1] ND barycentric rotating state; muStar, lStar
%            scalar; tWin [1 x 2] search window for the period [ND]; name
%            char (for the error).
%   OUTPUTS: km [m x 3] Moon-centred positions [km], closed; tau scalar [ND].
f   = @(~, y) cr3bp(y, muStar);
op  = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
sol = ode113(f, [0 tWin(2)], rv0, op);
tq  = linspace(tWin(1), tWin(2), 40001);
yq  = deval(sol, tq);
dq  = sqrt(sum((yq - rv0).^2, 1));
[~, kq] = min(dq);
tau = fminbnd(@(t) sqrt(sum((deval(sol, t) - rv0).^2)), tq(max(1, kq - 2)), tq(min(end, kq + 2)), ...
              optimset('TolX', 1e-12));
dKm = sqrt(sum((deval(sol, tau, 1:3) - rv0(1:3)).^2))*lStar;
if dKm > 10
    error('mission_movie:outline', '%s outline does not close: %.1f km after one period', name, dKm);
end
y  = deval(sol, linspace(0, tau, 6000));
km = (y(1:3, :).' - [1 - muStar 0 0])*lStar;
km(end, :) = km(1, :);
end

function dy = cr3bp(y, mu)
% CR3BP  Circular restricted three-body equations, rotating frame.
%   INPUTS: y [6 x 1] ND state; mu scalar.  OUTPUTS: dy [6 x 1].
r1 = sqrt((y(1) + mu)^2 + y(2)^2 + y(3)^2);
r2 = sqrt((y(1) - 1 + mu)^2 + y(2)^2 + y(3)^2);
dy = [y(4:6);
      y(1) + 2*y(5) - (1 - mu)*(y(1) + mu)/r1^3 - mu*(y(1) - 1 + mu)/r2^3;
      y(2) - 2*y(4) - (1 - mu)*y(2)/r1^3 - mu*y(2)/r2^3;
      -(1 - mu)*y(3)/r1^3 - mu*y(3)/r2^3];
end

function y = clip_top(y, top)
% CLIP_TOP  Values above top set to top; NaN stays NaN (min would fill it).
%   INPUTS: y [n x 1]; top scalar.  OUTPUTS: y [n x 1].
y(y > top) = top;
end

function [V, A] = clip_poly(X, a, toPx, box)
% CLIP_POLY  Clip a 3-D polyline to a screen rectangle (Liang-Barsky per
%   segment in screen coordinates; the projection is orthographic, hence
%   linear, so the screen parameter is also the 3-D parameter). Each visible
%   piece is emitted as [start; end; NaN] so no segment runs along the box
%   edge and nothing outside the box is drawn.
%   INPUTS:  X [n x 3] vertices [km]; a [n x 1] per-vertex alpha; toPx
%            handle [m x 3] -> [m x 2] px; box [1 x 4] x0 x1 y0 y1 [px].
%   OUTPUTS: V [3k x 3] vertices with NaN separators; A [3k x 1] alpha.
n = size(X, 1);
V = zeros(0, 3);  A = zeros(0, 1);
if n < 2, return; end
S  = toPx(X);
S0 = S(1:end - 1, :);  dS = diff(S);
t0 = zeros(n - 1, 1);  t1 = ones(n - 1, 1);  ok = true(n - 1, 1);
pp = [-dS(:, 1), dS(:, 1), -dS(:, 2), dS(:, 2)];
qq = [S0(:, 1) - box(1), box(2) - S0(:, 1), S0(:, 2) - box(3), box(4) - S0(:, 2)];
for cc = 1:4
    p = pp(:, cc);  q = qq(:, cc);
    ok(p == 0 & q < 0) = false;
    t = q./p;
    ng = p < 0;  t0(ng) = max(t0(ng), t(ng));
    ps = p > 0;  t1(ps) = min(t1(ps), t(ps));
end
k = find(ok & t0 <= t1);
m = numel(k);
if m == 0, return; end
X0 = X(k, :);  dX = X(k + 1, :) - X0;
a0 = a(k);     da = a(k + 1) - a0;
Pa = X0 + t0(k).*dX;  Pb = X0 + t1(k).*dX;
V  = reshape(permute(cat(3, Pa, Pb, nan(m, 3)), [3 1 2]), 3*m, 3);
A  = reshape([a0 + t0(k).*da, a0 + t1(k).*da, zeros(m, 1)].', [], 1);
end

function [V, C, A] = clip_polygon(V, C, A, toPx, box)
% CLIP_POLYGON  Sutherland-Hodgman clip of a planar polygon to a screen
%   rectangle, interpolating per-vertex colour and alpha (the projection is
%   orthographic, so screen and 3-D parameters agree).
%   INPUTS:  V [m x 3] vertices [km]; C [m x 3] colours; A [m x 1] alpha;
%            toPx handle; box [1 x 4] x0 x1 y0 y1 [px].
%   OUTPUTS: V, C, A clipped (fewer than 3 rows = nothing visible).
for ed = 1:4
    m = size(V, 1);
    if m == 0, break; end
    S = toPx(V);
    switch ed
        case 1, dist = S(:, 1) - box(1);
        case 2, dist = box(2) - S(:, 1);
        case 3, dist = S(:, 2) - box(3);
        otherwise, dist = box(4) - S(:, 2);
    end
    Vo = zeros(2*m, 3);  Co = zeros(2*m, 3);  Ao = zeros(2*m, 1);  no = 0;
    for kk = 1:m
        nx = mod(kk, m) + 1;
        inA = dist(kk) >= 0;  inB = dist(nx) >= 0;
        if inA
            no = no + 1;
            Vo(no, :) = V(kk, :);  Co(no, :) = C(kk, :);  Ao(no) = A(kk);
        end
        if inA ~= inB
            t  = dist(kk)/(dist(kk) - dist(nx));
            no = no + 1;
            Vo(no, :) = V(kk, :) + t*(V(nx, :) - V(kk, :));
            Co(no, :) = C(kk, :) + t*(C(nx, :) - C(kk, :));
            Ao(no) = A(kk) + t*(A(nx) - A(kk));
        end
    end
    V = Vo(1:no, :);  C = Co(1:no, :);  A = Ao(1:no);
end
end

function set_plume(h, V, C, A)
% SET_PLUME  Update one plume patch; hide it when nothing is left.
%   INPUTS: h patch; V [m x 3]; C [m x 3]; A [m x 1].  OUTPUTS: none.
if size(V, 1) < 3
    set(h, 'Visible', 'off');
else
    set(h, 'Vertices', V, 'Faces', 1:size(V, 1), 'FaceVertexCData', C, ...
           'FaceVertexAlphaData', A, 'Visible', 'on');
end
end

function V = plume_poly(c0, e, nv, Lkm, waist, wKm)
% PLUME_POLY  Tapered plume polygon: nozzle+, waist+, tip, waist-, nozzle-.
%   INPUTS:  c0 [1 x 3] nozzle [km]; e [1 x 3] unit exhaust; nv [1 x 3] unit
%            across; Lkm scalar length; waist scalar fraction; wKm [1 x 3]
%            full widths at nozzle, waist, tip [km].
%   OUTPUTS: V [5 x 3] vertices [km].
cw = c0 + waist*Lkm*e;
V  = [c0 + wKm(1)/2*nv; cw + wKm(2)/2*nv; c0 + Lkm*e; cw - wKm(2)/2*nv; c0 - wKm(1)/2*nv];
end

function ex = page_extent(R, mask)
% PAGE_EXTENT  Min/max of R over the samples where mask is true.
%   INPUTS: R [N x 3 x K]; mask [N x K] logical.  OUTPUTS: ex [3 x 2].
ex = [inf(3, 1) -inf(3, 1)];
for kk = 1:size(R, 3)
    q = R(mask(:, kk), :, kk);
    if ~isempty(q)
        ex = [min(ex(:, 1), min(q, [], 1).') max(ex(:, 2), max(q, [], 1).')];
    end
end
end

function [x, y] = on_station_steps(tArrD, d0)
% ON_STATION_STEPS  Step curve of satellites on station from day 0 to d0.
%   INPUTS: tArrD [1 x K] arrival days; d0 scalar day.  OUTPUTS: x, y [1 x m].
ta = sort(tArrD(tArrD <= d0));
n  = numel(ta);
x  = [0, reshape([ta; ta], 1, []), d0];
y  = [0, reshape([0:n - 1; 1:n], 1, []), n];
end

function y = unimodal_hull(x)
% UNIMODAL_HULL  Smallest single-peaked sequence above x: rises to the peak,
%   then falls (dips are bridged).  INPUTS: x [n x 1].  OUTPUTS: y [n x 1].
y = min(cummax(x), flipud(cummax(flipud(x))));
end

function y = hann_smooth(x, h)
% HANN_SMOOTH  Raised-cosine moving average of half-width h samples, edges
%   padded with the end values (compact support: a constant stretch longer
%   than the window stays exactly constant; log-concave kernel: a unimodal
%   input stays unimodal).  INPUTS: x [n x 1]; h integer.  OUTPUTS: y [n x 1].
k  = 0.5*(1 + cos(pi*(-h:h).'/(h + 1)));
k  = k/sum(k);
xp = [repmat(x(1), h, 1); x(:); repmat(x(end), h, 1)];
y  = conv(xp, k, 'valid');
end

function f = full_path(f)
% FULL_PATH  Absolute, normalised path for comparisons.  INPUTS/OUTPUTS: char.
jf = java.io.File(f);
f  = char(jf.getCanonicalPath());
end

function s = thousands(n)
% THOUSANDS  Integer with comma separators.  INPUTS: n scalar.  OUTPUTS: s char.
s = sprintf('%d', round(n));
for kk = numel(s) - 3:-3:1
    s = [s(1:kk) ',' s(kk + 1:end)];
end
end

function s = onoff(c)
% ONOFF  'on' / 'off' from a logical.  INPUTS: c logical.  OUTPUTS: s char.
if c, s = 'on'; else, s = 'off'; end
end

function close_if(fig)
% CLOSE_IF  Close the figure if it still exists.  INPUTS: fig handle.
if isgraphics(fig), close(fig); end
end

function close_vw(vw, partFile)
% CLOSE_VW  Close a VideoWriter if open and delete an unfinished partial
%   file (a render that did not complete never touches outFile).
%   INPUTS: vw VideoWriter or []; partFile char ('' = nothing to delete).
if ~isempty(vw) && isvalid(vw), close(vw); end
if ~isempty(partFile) && isfile(partFile), delete(partFile); end
end

function [muStar, lStar, tStar, prov] = results_units(resultsMat)
% RESULTS_UNITS  The CR3BP units and the run's provenance of a results file.
%   Units come from its fingerprint (written by the study from the
%   catalog); a file without one falls back to mission_inputs(). The
%   provenance note is run_meta.verdictSuffix when present (the text the
%   study appended to its VERDICT), else rebuilt from run_meta.reduced and
%   run_meta.polishSource; '' marks a fresh, full-knob run.
%   INPUTS:  resultsMat char.
%   OUTPUTS: muStar [ND], lStar [km], tStar [s] scalars; prov struct
%            (.fresh logical, .note char, .units char).
w  = whos('-file', resultsMat);
nm = {w.name};
if any(strcmp(nm, 'fingerprint'))
    F = load(resultsMat, 'fingerprint');
    muStar = F.fingerprint.muStar;  lStar = F.fingerprint.lStar;  tStar = F.fingerprint.tStar;
    src = 'fingerprint';
else
    Mi = mission_inputs();
    muStar = Mi.muStar;  lStar = Mi.lStar;  tStar = Mi.tStar;
    src = 'mission_inputs';
end
if ~any(strcmp(nm, 'run_meta'))
    prov = struct('fresh', false, 'note', 'no run_meta in the results file (provenance unknown)', ...
                  'units', src);
    return
end
R  = load(resultsMat, 'run_meta');
rm = R.run_meta;
if isfield(rm, 'verdictSuffix')
    note = strtrim(rm.verdictSuffix);
else
    parts = {};
    if isfield(rm, 'polishSource') && any(strcmp({rm.polishSource.origin}, 'reused'))
        parts{end+1} = 'REUSED polishes (not re-certified in that run)';
    end
    if isfield(rm, 'reduced') && ~isempty(rm.reduced)
        parts{end+1} = ['reduced: ' strjoin(rm.reduced(:).', ', ')];
    end
    note = strjoin(parts, '; ');
end
prov = struct('fresh', isempty(note), 'note', note, 'units', src);
end
