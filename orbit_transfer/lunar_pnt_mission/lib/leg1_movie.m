function info = leg1_movie(leg1Mat, opts)
%% Purpose:
%
%   Renders the movie of LEG 1 of the lunar PNT mission from the SAVED
%   leg-1 study results (nothing is re-solved): six 5 N transfer stages
%   stacked on one inertially fixed GTO near the Earth leave it one by one,
%   fly their polished minimum-time transfers with a yellow-orange thrust
%   plume along the exhaust direction, and coast on six slots of the DRO
%   around the Moon until the hand-off epoch T0, where leg 2 (mission_movie)
%   begins. The last frame is mission_movie's FIRST frame (same view
%   direction, target and scene size), so the two movies join there.
%
%   Frame (1280 x 720, dark), the house layout of mission_movie: title and
%   story caption on top; the scene (Earth-Moon rotating frame, Moon-centred
%   km, orthographic camera) in the middle with counters and the regime
%   statement at the left and a state legend at the right; the strip
%   (satellites in the DRO versus mission day, time cursor, T0 marked) at
%   the bottom. Every axes has a fixed pixel position; every camera property
%   is manual. A black mask axes above the scene and below the strip and
%   overlay limits everything drawn in the scene (Earth, Moon, outlines,
%   trails, plumes) to the scene area clipBox.
%
%   Positions of stage k (leg-1 time t in ND, P1(k) of the results):
%     t <  tDep : on the GTO, phase mod(gtoPhase0 + t/tauG, 1) at
%                 orientation orient0 + sgn*(180/pi)*t deg. The shape is
%                 L.gtoState(oRef, phase) rotated about the EARTH's z-axis
%                 by (orientation - oRef); this rotation is verified against
%                 L.gtoState at another orientation to 1e-9 ND at load
%                 (leg1_movie:gtoRotation; a rotation about the barycentre
%                 misses by 6e-3 ND).
%     tDep..tArr: the flown history P1(k).t, P1(k).X, cubic Hermite in
%                 position with the stored velocities.
%     t >  tArr : on the DRO, phase mod(droPhaseArr + (t - tArr)/tauD, 1).
%   Continuity at both junctions is asserted (< 1 km, leg1_movie:consistency).
%
%   Stack: the six stages start as one stack on the GTO; markers closer
%   than 7 px carry a count label ("6 stacked") until 0.1 d after the last
%   departure.
%
%   Plume: exhaust e = -alpha/|alpha| (alpha is the flown thrust direction,
%   -lambda_v/|lambda_v|), alpha linearly interpolated in leg-1 time on the
%   flight grid; drawn only while tDep <= t < tArr. Its length is a fixed
%   number of pixels (opts.plumePx; illustrative).
%
%   Time map (NON-UNIFORM, smooth): a hold on day 0, then a dwell density
%   that spends budgetS(1) video seconds on [0, splitDays(1)] (departures,
%   GTO revolutions), budgetS(2) on the cruise to splitDays(2) (default:
%   the first arrival minus 0.25 d) and budgetS(3) on the arrivals and the
%   coast to T0, Gaussian-smoothed in mission day with an ease-in and an
%   ease-out, then a hold on T0. Day counter and days-per-second readout.
%
%   Camera: azimuth fixed at mission_movie's; elevation elTop (a nearly
%   top-down view of the planar geometry) blended by a raised cosine over
%   the arrival segment to mission_movie's elevation. Target and zoom from
%   the data, per map frame: the world bounding box of the six spacecraft
%   plus the Earth and the GTO while any stage is within earthKeepKm of the
%   Earth, plus the Moon and the DRO once any stage is within moonKeepKm of
%   the Moon. Each box edge is smoothed (1 s running max/min, then a
%   raised cosine of targetSmoothS half-width); the target is the smoothed
%   box's centre (shifted for the off-centre scene area). Zoom requirement:
%   the smoothed box with 6 % margin, and every point with its own margin
%   (12 % on a transferring stage, tapering to none over the last 0.3 d
%   before its arrival); then running max, unimodal hull (one zoom-out, one
%   zoom-in), raised-cosine smoothing (zoomSmoothS), re-inflation until the
%   requirement is covered. Nothing auto-rescales. Over the last endRampS
%   before the closing hold, target and zoom are ramped (raised cosine) onto
%   mission_movie's end framing, which every frame of the hold uses
%   exactly. Every map frame is checked: all six markers >= 4 px inside
%   clipBox, else leg1_movie:containment.
%
%   Earth: a sphere of radius 6378.137 km at (-lStar, 0, 0) Moon-centred,
%   textured with pumpkyn's earth-clouds-4k.jpg (subsampled by 4) when the
%   file exists (else a shaded blue sphere), turning at the Earth's
%   rotation rate relative to the rotating frame. Moon: mission_movie's
%   textured globe (pumpkyn Moon_2k.jpg).
%
%   VETTING: the results are rendered only if run_meta.verdict is exactly
%   'PASS', six transfers are accepted with G.pass true and every test seam
%   recorded in run_meta.knobs (seamRefuse, seamRefute, seamFinalT0) is
%   empty / zero; else leg1_movie:verdict, unless opts.allowUnvetted (then
%   "UNVETTED RESULTS" is stamped above the regime statement).
%
%   PROTECTED FILE: results/lunar_pnt_mission.mp4 (the owner's leg-2 movie)
%   is never written: an outFile (or its partial file) that resolves to it,
%   relative paths taken against MATLAB's pwd and '..' resolved
%   (canon_path), throws leg1_movie:protected before anything is loaded.
%
%   Modes: MP4 (default; written to <name>.partial.mp4, moved onto outFile
%   after VideoReader confirms frame count and size); PREVIEW MP4
%   (preview = true: the whole time map subsampled to previewSec*fps
%   frames, default outFile results/lunar_pnt_leg1_preview.mp4; refused onto
%   the full movie's file, leg1_movie:partialOutFile); PNG stills when
%   frameList or dayList is non-empty (<previewDir>/<pngPrefix>_<frame>_
%   day<d>.png).
%
%% Inputs:
%
%  leg1Mat                  char                    lunar_pnt_leg1_study.mat
%                                                   ('' = results/ beside
%                                                   this campaign); needs P1,
%                                                   T0, run_meta
%  opts                     struct (optional)
%    .outFile               char                    MP4 path
%                                                   [results/lunar_pnt_leg1.mp4]
%    .leg2Mat               char                    leg-2 results whose
%                                                   mission_movie first frame
%                                                   is this movie's last
%                                                   [results/lunar_pnt_mission_study.mat]
%    .catalogMat            char ['']               GTO -> DRO package for
%                                                   leg1_inputs ('' = the file
%                                                   named in run_meta); found
%                                                   by SHA-256 BEFORE loading
%                                                   (leg1_package_for), else
%                                                   run_meta.packageSnapshot,
%                                                   else a snapshot in inputs/
%                                                   beside the results or in
%                                                   results/inputs; none ->
%                                                   leg1_movie:identity
%    .fps                   scalar [30]             Frame rate [1/s]
%    .durationSec           scalar [30]             Total length [s]
%    .holdS                 [1 x 2] [1.5 2.5]       Opening / closing hold [s]
%    .budgetS               [1 x 3] [6 14 6]        Departures / cruise /
%                                                   arrivals [s] (rescaled to
%                                                   the time between holds)
%    .splitDays             [1 x 2] [0.3 NaN]       Segment boundaries [d]
%                                                   (NaN = first arrival
%                                                   - 0.25 d)
%    .elTop                 scalar [70]             Elevation before the
%                                                   arrival blend [deg]
%    .earthKeepKm           scalar [150000]         Earth + GTO framed while a
%                                                   stage is this close [km]
%    .moonKeepKm            scalar [330000]         Moon + DRO framed once a
%                                                   stage is this close [km]
%    .targetSmoothS         scalar [2.5]            Target smoothing [s]
%    .zoomSmoothS           scalar [3]              Zoom smoothing [s]
%    .endRampS              scalar [3]              Ramp onto the end framing [s]
%    .frameSize             [1 x 2] [1280 720]      Output width, height [px]
%    .visible               logical [false]         Show the figure
%    .preview               logical [false]         Preview MP4 (see Modes)
%    .previewSec            scalar [2]              Preview length [s]
%    .frameList             [1 x m] []              Frame indices -> PNGs
%    .dayList               [1 x m] []              Mission days in [0, T0]
%                                                   -> PNGs
%    .previewDir            char                    PNG folder
%                                                   [results/movie_preview]
%    .pngPrefix             char ['leg1']           PNG file name prefix
%    .tailDays              scalar [0.3]            Trail length [d]
%    .plumePx               scalar [44]             On-screen plume length [px]
%    .moonTexture           char                    Moon texture image
%                                                   [pumpkyn +util/Moon_2k.jpg]
%    .earthTexture          char                    Earth texture image
%                                                   [pumpkyn +util/
%                                                   earth-clouds-4k.jpg]
%    .quality               scalar [95]             VideoWriter quality
%    .zoomScale             scalar [1]              TEST SEAM: multiplies the
%                                                   zoom before the
%                                                   containment check (< 1
%                                                   forces
%                                                   leg1_movie:containment)
%    .allowUnvetted         logical [false]         Render results that fail
%                                                   the vetting below, with
%                                                   UNVETTED RESULTS stamped
%                                                   in the regime text
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
%  info.timeMap             struct                  .tv .day .rate [F x 1]
%                                                   (video s, day, days per
%                                                   s); .budgetS [1 x 5] hold,
%                                                   departures, cruise,
%                                                   arrivals, hold [s];
%                                                   .splitDays [1 x 2]
%  info.camera              struct                  .az [deg]; .el [F x 1]
%                                                   [deg]; .target [F x 3];
%                                                   .Rv .req [F x 1] [km];
%                                                   .RvEnd .targetEnd
%                                                   .positionEnd
%                                                   .viewAngleEnd .elEnd
%                                                   (mission_movie's first
%                                                   frame); .clipBox;
%                                                   .minMarginPx .minMarginDay;
%                                                   .zoomMax (max Rv / RvEnd);
%                                                   .toPx handle (X [m x 3]
%                                                   km, f map frame) ->
%                                                   [m x 2] figure px
%  info.posKm               [n x 3 x 6]             Drawn spacecraft positions
%                                                   (Moon-centred km)
%  info.scMarginPx          [n x 1]                 Smallest marker distance to
%                                                   the clipBox edge [px]
%  info.Rv                  [n x 1]                 Half-height per render [km]
%  info.camPos, .camTarget  [n x 3]                 Camera as set, per render
%  info.camViewAngle        [n x 1]                 likewise [deg]
%  info.camModesManual      [n x 1] logical         All four camera modes manual
%  info.phase               [n x 6]                 1 on the GTO, 2 thrusting,
%                                                   3 in the DRO
%  info.plumeDir            [n x 3 x 6]             Unit exhaust direction
%                                                   (NaN when not thrusting)
%  info.plumeDrawn          [n x 6] logical         Plume drawn
%  info.plumeScreenDir      [n x 2 x 6]             Drawn plume tip minus
%                                                   nozzle on screen [px]
%  info.stackMax            [n x 1]                 Largest group of markers
%                                                   within 7 px (a count label
%                                                   is drawn on each group up
%                                                   to 0.1 d after the last
%                                                   departure: the stack)
%  info.outlineDrawn        [n x 2] logical         GTO, DRO outline drawn
%  info.outlineInBox        [n x 2] logical         ... with a vertex inside
%                                                   clipBox
%  info.moonTextured        logical                 Moon carries its texture
%  info.earthTextured       logical                 Earth carries its texture
%  info.T0Days              scalar                  Hand-off epoch [d]
%  info.tDepDays, .tArrDays [1 x 6]                 Departures, arrivals [d]
%  info.sgn, .orient0       scalars                 GTO drift sign, orientation
%                                                   at t = 0 [deg]
%  info.gtoRotationErr      scalar                  Rotation check [ND]
%  info.junctionKm          [6 x 2]                 Departure, arrival
%                                                   continuity [km]
%  info.dayText             {n x 1} cell            The day counter's String
%                                                   as DISPLAYED (read back
%                                                   from the text object)
%  info.regimeText          {1 x m} cell            Regime statement lines
%                                                   (first 'UNVETTED RESULTS'
%                                                   when unvetted)
%  info.unvetted            logical                 Rendered under
%                                                   allowUnvetted
%  info.vetting             char                    '' or why the results are
%                                                   not vetted
%  info.packageFile         char                    Package file used
%
%% Revision History:
%  M. Casey                                                   (c) 11/11/2026
%  M. Casey  final fix wave: vetting (allowUnvetted), protected
%            owner's movie, package by SHA-256 before loading,
%            info.dayText                                     (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(fileparts(mfilename('fullpath')));        % lunar_pnt_mission/
if nargin == 0
    run(fullfile(here, 'setup_lunar_pnt_paths.m'));
    info = leg1_movie('', struct('dayList', [0 0.15 1.2 2.6], 'pngPrefix', 'leg1_demo'));
    fprintf('leg1_movie self-demo: %d preview PNGs in %s (%.1f s)\n', ...
            numel(info.previewPngs), fileparts(info.previewPngs{1}), info.wallS);
    return
end
if nargin < 2 || isempty(opts), opts = struct(); end
if isempty(leg1Mat)
    leg1Mat = fullfile(here, 'results', 'lunar_pnt_leg1_study.mat');
end
tWall = tic;

% ---- options ----
utilDir = fullfile(getenv('HOME'), 'Desktop', 'proj7', 'external', 'pumpkyn', 'src', '+pumpkyn', '+util');
def = struct('outFile',       fullfile(here, 'results', 'lunar_pnt_leg1.mp4'), ...
             'leg2Mat',       fullfile(here, 'results', 'lunar_pnt_mission_study.mat'), ...
             'catalogMat',    '', ...
             'fps',           30, ...
             'durationSec',   30, ...
             'holdS',         [1.5 2.5], ...
             'budgetS',       [6 14 6], ...
             'splitDays',     [0.3 NaN], ...
             'elTop',         70, ...
             'earthKeepKm',   150000, ...
             'moonKeepKm',    330000, ...
             'targetSmoothS', 2.5, ...
             'zoomSmoothS',   3, ...
             'endRampS',      3, ...
             'frameSize',     [1280 720], ...
             'visible',       false, ...
             'preview',       false, ...
             'previewSec',    2, ...
             'frameList',     [], ...
             'dayList',       [], ...
             'previewDir',    fullfile(here, 'results', 'movie_preview'), ...
             'pngPrefix',     'leg1', ...
             'tailDays',      0.3, ...
             'plumePx',       44, ...
             'moonTexture',   fullfile(utilDir, 'Moon_2k.jpg'), ...
             'earthTexture',  fullfile(utilDir, 'earth-clouds-4k.jpg'), ...
             'quality',       95, ...
             'zoomScale',     1, ...
             'allowUnvetted', false);
for f = fieldnames(opts).'
    if ~isfield(def, f{1})
        error('leg1_movie:opts', 'unknown option ''%s''', f{1});
    end
    def.(f{1}) = opts.(f{1});
end
o = def;
fullOut    = fullfile(here, 'results', 'lunar_pnt_leg1.mp4');
previewOut = fullfile(here, 'results', 'lunar_pnt_leg1_preview.mp4');
if o.preview && ~isfield(opts, 'outFile')
    o.outFile = previewOut;
end
% the owner's leg-2 movie is never an output (paths resolved against MATLAB's pwd)
protectedMp4 = canon_path(fullfile(here, 'results', 'lunar_pnt_mission.mp4'));
[pDir0, pName0] = fileparts(canon_path(o.outFile));
if any(strcmp(protectedMp4, {canon_path(o.outFile), fullfile(pDir0, [pName0 '.partial.mp4'])}))
    error('leg1_movie:protected', 'refusing to write %s (the owner''s movie): outFile resolves to it', protectedMp4);
end
if numel(o.frameSize) ~= 2 || any(mod(o.frameSize, 2))
    error('leg1_movie:frameSize', 'frameSize must be two even integers');
end

% ---- layout constants (mission_movie's) ----
rMoonKm  = 1737.4;
rEarthKm = 6378.137;
W        = 1280;                           % layout width  [px]
Hh       = 720;                            % layout height [px]
scnPos   = [0 194 1280 450];               % scene axes [px] (as mission_movie)
stripPos = [96 62 1088 106];               % strip axes [px]
halfH    = scnPos(4)/2;                    % scene half-height [px]
clipBox  = [250 1040 200 640];             % visible scene area x0 x1 y0 y1 [px]
stackPx  = 7;                              % markers closer than this form a stack [px]
omegaE   = 7.2921159e-5;                   % Earth rotation rate, inertial [rad/s]

% ---- data ----
if ~isfile(leg1Mat)
    error('leg1_movie:noResults', 'results file not found: %s', leg1Mat);
end
D = load(leg1Mat, 'P1', 'T0', 'run_meta');
for nm = {'P1', 'T0', 'run_meta'}
    if ~isfield(D, nm{1}) || isempty(D.(nm{1}))
        error('leg1_movie:results', '%s missing or empty in %s', nm{1}, leg1Mat);
    end
end
P1 = D.P1;  T0 = double(D.T0);  rm = D.run_meta;
vetWhy = vet_results(P1, rm);
unvetted = ~isempty(vetWhy);
if unvetted && ~o.allowUnvetted
    error('leg1_movie:verdict', 'refusing to render unvetted leg-1 results (%s); opts.allowUnvetted overrides', vetWhy);
end
nSc = numel(P1);
if nSc ~= 6 || ~all([P1.accepted])
    error('leg1_movie:results', 'leg-1 results need six accepted transfers (%d, %d accepted)', ...
          nSc, nnz([P1.accepted]));
end
if ~isequal(sort([P1.k]), 1:6)
    error('leg1_movie:results', 'satellite numbers k are not 1..6');
end
[~, ordK] = sort([P1.k]);
P1 = P1(ordK);                                           % index = leg-2 number k
L  = load_inputs(o.catalogMat, rm, here, leg1Mat);
mu = L.muStar;  lStar = L.lStar;  tStar = L.tStar;
dayNd = tStar/86400;
orient0 = P1(1).orient0;  ph0 = P1(1).gtoPhase0;
if any(abs([P1.orient0] - orient0) > 0) || any(abs([P1.gtoPhase0] - ph0) > 0)
    error('leg1_movie:results', 'the six stages do not share one GTO (orient0, gtoPhase0)');
end
if abs(rm.sgn - L.sgn) > 0 || abs(rm.tauG - L.tauG) > 1e-12 || abs(rm.tauD - L.tauD) > 1e-12
    error('leg1_movie:results', 'run_meta sgn/tauG/tauD differ from leg1_inputs');
end

% S: everything the position function needs
S = struct('mu', mu, 'lStar', lStar, 'sgn', L.sgn, 'tauG', L.tauG, 'tauD', L.tauD, ...
           'oRef', L.sheetOrient(1), 'orient0', orient0, 'ph0', ph0, ...
           'gtoFun', L.gtoState, 'slotFun', L.slotState);
% the rotation of the GTO shape: about the Earth, checked against gtoState
oChk = S.oRef - 30;  phC = (0:11)/12;
xA = L.gtoState(oChk, phC);
xR = rotate_earth(L.gtoState(S.oRef, phC), oChk - S.oRef, mu);
gtoRotErr = max(abs(xA(1:3, :) - xR), [], 'all');
if gtoRotErr > 1e-9
    error('leg1_movie:gtoRotation', 'GTO rotation about the Earth misses gtoState by %.3g ND', gtoRotErr);
end
sat = struct('tDep', cell(1, nSc), 'tArr', [], 't', [], 'X', [], 'alpha', [], 'phArr', []);
for kk = 1:nSc
    [tg, ia] = unique(double(P1(kk).t(:)), 'stable');
    if any(diff(tg) <= 0)
        error('leg1_movie:flightGrid', 'stage %d flight grid not increasing', kk);
    end
    sat(kk).t = tg;
    sat(kk).X = double(P1(kk).X(1:6, ia));
    sat(kk).alpha = double(P1(kk).alpha(:, ia)).';
    sat(kk).tDep = double(P1(kk).tDep);
    sat(kk).tArr = double(P1(kk).tArr);
    sat(kk).phArr = double(P1(kk).droPhaseArr);
end
S.sat = sat;
% continuity at the two junctions (GTO -> flight, flight -> DRO)
jKm = zeros(nSc, 2);
for kk = 1:nSc
    a = sat(kk);
    pg = gto_pos(S, a.tDep - 1e-12);
    jKm(kk, 1) = sqrt(sum((pg - a.X(1:3, 1).').^2))*lStar;
    xs = L.slotState(a.phArr);
    jKm(kk, 2) = sqrt(sum((xs(1:3) - a.X(1:3, end)).^2))*lStar;
end
if any(jKm(:) > 1)
    error('leg1_movie:consistency', 'stage path discontinuous: max %.3g km at a junction', max(jKm(:)));
end
tDepD = [sat.tDep]*dayNd;  tArrD = [sat.tArr]*dayNd;  T0d = T0*dayNd;
if T0d < max(tArrD)
    error('leg1_movie:results', 'T0 %.4f d precedes the last arrival %.4f d', T0d, max(tArrD));
end
dEnd = T0d;

% outlines
phO    = linspace(0, 1, 721);
gtoRef = L.gtoState(S.oRef, phO);                         % [6 x m] at orientation oRef
xd     = L.slotState(phO);
droKm  = (xd(1:3, :).' - [1 - mu 0 0])*lStar;
droKm(end, :) = droKm(1, :);
xEarth = [-lStar 0 0];                                    % Earth, Moon-centred [km]
gtoAt  = @(tc) (rotate_earth(gtoRef, S.orient0 + S.sgn*(180/pi)*tc - S.oRef, mu).' - ...
                [1 - mu 0 0])*lStar;

% ---- end framing: mission_movie's first frame ----
EF = end_framing(o.leg2Mat);

% ---- time map ----
nF   = round(o.durationSec*o.fps);
tv   = (0:nF - 1).'/o.fps;
holdS = o.holdS;
tMov = o.durationSec - sum(holdS);
split = o.splitDays;
if isnan(split(2)), split(2) = min(tArrD) - 0.25; end
if ~(split(1) > 0 && split(2) > split(1) && split(2) < dEnd)
    error('leg1_movie:splitDays', 'splitDays [%g %g] must satisfy 0 < s1 < s2 < T0 = %.4f d', ...
          split, dEnd);
end
budget = o.budgetS*tMov/sum(o.budgetS);
dd   = linspace(0, dEnd, 20001).';
w    = budget(1)/split(1)*(dd < split(1)) + budget(2)/(split(2) - split(1))*(dd >= split(1) & dd < split(2)) + ...
       budget(3)/(dEnd - split(2))*(dd >= split(2));
sig  = 0.06;                                % smoothing width [d]
hd   = dd(2) - dd(1);
kx   = (-ceil(4*sig/hd):ceil(4*sig/hd)).'*hd;
ker  = exp(-0.5*(kx/sig).^2);  ker = ker/sum(ker);
nPad = numel(kx);
wp   = [repmat(w(1), nPad, 1); w; repmat(w(end), nPad, 1)];
ws   = conv(wp, ker, 'same');
ws   = ws(nPad + 1:nPad + numel(dd));
ws   = ws.*(1 + 1.0*exp(-(dd/0.03).^2) + 1.0*exp(-((dEnd - dd)/0.03).^2));
sCum = cumtrapz(dd, ws);
scl  = tMov/sCum(end);
sCum = sCum*scl;  ws = ws*scl;
tm   = min(max(tv - holdS(1), 0), tMov);
day  = interp1(sCum, dd, tm);
day(tv >= holdS(1) + tMov) = dEnd;
day(1) = 0;
rate = 1./interp1(dd, ws, day);
rate(tv < holdS(1) | tv >= holdS(1) + tMov) = 0;
if any(diff(day) < 0) || day(end) ~= dEnd
    error('leg1_movie:timeMap', 'time map not monotone or does not reach T0');
end
TM = struct('tv', tv, 'day', day, 'rate', rate, 'splitDays', split, ...
            'budgetS', [holdS(1) interp1(dd, sCum, split(1)) ...
                        interp1(dd, sCum, split(2)) - interp1(dd, sCum, split(1)) ...
                        tMov - interp1(dd, sCum, split(2)) holdS(2)]);
tcF = day/dayNd;                            % clock per map frame [ND]

% ---- positions per map frame ----
posF = zeros(nF, 3, nSc);
for kk = 1:nSc
    posF(:, :, kk) = sat_pos(S, kk, tcF);
end

% ---- camera ----
az    = EF.az;
tvHold = holdS(1) + tMov;                   % first frame of the closing hold [s]
tvC0  = holdS(1) + TM.budgetS(2) + TM.budgetS(3);          % start of the arrival segment [s]
wEl   = raised_ramp(tv, tvC0, tvHold);
el    = o.elTop + (EF.el - o.elTop)*wEl;
wEnd  = raised_ramp(tv, tvHold - o.endRampS, tvHold);
EH = zeros(nF, 3);  EV = zeros(nF, 3);  VD = zeros(nF, 3);
for ff = 1:nF
    Tv = viewmtx(az, el(ff));
    EH(ff, :) = Tv(1, 1:3);  EV(ff, :) = Tv(2, 1:3);
    VD(ff, :) = [cosd(el(ff))*sind(az) -cosd(el(ff))*cosd(az) sind(el(ff))];
end
cxy = [scnPos(1) + scnPos(3)/2, scnPos(2) + scnPos(4)/2];   % scene centre [px]
ccb = [(clipBox(1) + clipBox(2))/2, (clipBox(3) + clipBox(4))/2];
lim = [clipBox(1) - cxy(1), clipBox(2) - cxy(1), clipBox(3) - cxy(2), clipBox(4) - cxy(2)];
% what must be in view, per map frame
dEk = squeeze(sqrt(sum((posF - xEarth).^2, 2)));           % [nF x 6] km
dMk = squeeze(sqrt(sum(posF.^2, 2)));
needEarth = any(dEk < o.earthKeepKm, 2);
moonOn = find(any(dMk < o.moonKeepKm, 2), 1);
needMoon = false(nF, 1);
if ~isempty(moonOn), needMoon(moonOn:end) = true; end
mg = ones(nF, nSc);
for kk = 1:nSc
    tr = day >= tDepD(kk) & day < tArrD(kk);
    mg(tr, kk) = 1 + 0.12*min(1, (tArrD(kk) - day(tr))/0.3);
end
gtoSub = 1:10:numel(phO);
droSub = 1:10:size(droKm, 1);
% world bounding box of what must be seen, per map frame; its edges are
% smoothed (running max/min, then raised cosine) so target and zoom glide
bx = zeros(nF, 6);                          % xmin xmax ymin ymax zmin zmax [km]
pts = cell(nF, 1);  pmg = cell(nF, 1);
for ff = 1:nF
    qE = zeros(0, 3);  qM = zeros(0, 3);
    if needEarth(ff)
        gq = gtoAt(tcF(ff));
        qE = [gq(gtoSub, :); xEarth + rEarthKm*[eye(3); -eye(3)]];
    end
    if needMoon(ff)
        qM = [droKm(droSub, :); rMoonKm*[eye(3); -eye(3)]];
    end
    q = [reshape(posF(ff, :, :), 3, nSc).'; qE; qM];
    pts{ff} = q;
    pmg{ff} = [mg(ff, :).'; 1.05*ones(size(qE, 1) + size(qM, 1), 1)];
    bx(ff, :) = reshape([min(q, [], 1); max(q, [], 1)], 1, 6);
end
hwT = max(1, round(o.targetSmoothS*o.fps));
wMax = round(1*o.fps);
for cc = 1:6
    if mod(cc, 2) == 1
        bx(:, cc) = hann_smooth(-movmax(-bx(:, cc), [wMax wMax]), hwT);
    else
        bx(:, cc) = hann_smooth(movmax(bx(:, cc), [wMax wMax]), hwT);
    end
end
tgS = zeros(nF, 3);  needB = zeros(nF, 1);
for ff = 1:nF
    cr = [bx(ff, [1 3 5]); bx(ff, [2 3 5]); bx(ff, [1 4 5]); bx(ff, [2 4 5]); ...
          bx(ff, [1 3 6]); bx(ff, [2 3 6]); bx(ff, [1 4 6]); bx(ff, [2 4 6])];
    hq = cr*EH(ff, :).';  vq = cr*EV(ff, :).';
    hmn = min(hq);  hmx = max(hq);  vmn = min(vq);  vmx = max(vq);
    kpp = 1.06*max((hmx - hmn)/(clipBox(2) - clipBox(1)), (vmx - vmn)/(clipBox(4) - clipBox(3)));
    hcn = (hmn + hmx)/2 - (ccb(1) - cxy(1))*kpp;
    vcn = (vmn + vmx)/2 - (ccb(2) - cxy(2))*kpp;
    tgS(ff, :) = hcn*EH(ff, :) + vcn*EV(ff, :);
    needB(ff) = kpp*halfH;
end
tg  = tgS.*(1 - wEnd) + EF.target.*wEnd;
% zoom requirement against this target: the smoothed box, and every point
% with its margin (12 % on a transferring stage)
need = zeros(nF, 1);
for ff = 1:nF
    q = pts{ff};  m = pmg{ff};
    hq = (q - tg(ff, :))*EH(ff, :).';  vq = (q - tg(ff, :))*EV(ff, :).';
    need(ff) = max([max(max(hq, 0).*m/lim(2)), max(max(-hq, 0).*m/(-lim(1))), ...
                    max(max(vq, 0).*m/lim(4)), max(max(-vq, 0).*m/(-lim(3)))])*halfH;
end
need = max(need, needB.*(1 - wEnd));
LEnd = log(EF.RvEnd);
Lq   = log(max(need, 1));
hw   = max(1, round(o.zoomSmoothS*o.fps));
Lz   = unimodal_hull(movmax(Lq, [wMax wMax]));
Lz   = hann_smooth(Lz, hw);
for it = 1:20
    gap = Lq - Lz;
    if all(gap <= 0), break; end
    Lz = unimodal_hull(Lz + hann_smooth(movmax(max(gap, 0) + 0.002*(gap > 0), [wMax wMax]), hw));
end
Lz  = LEnd + (Lz - LEnd).*(1 - wEnd);       % onto the end framing over endRampS
Lz(wEnd == 1) = LEnd;
RvF = exp(Lz);
RvF(wEnd == 1) = EF.RvEnd;
RvF = RvF*o.zoomScale;                      % TEST SEAM (1 in every real render)
toPx = @(X, f) cxy + [(X - tg(f, :))*EH(f, :).', (X - tg(f, :))*EV(f, :).']*(halfH/RvF(f));
boxMargin = @(Sx) min([Sx(:, 1) - clipBox(1), clipBox(2) - Sx(:, 1), ...
                       Sx(:, 2) - clipBox(3), clipBox(4) - Sx(:, 2)], [], 2);
mPx = zeros(nF, 1);
for ff = 1:nF
    mPx(ff) = min(boxMargin(toPx(reshape(posF(ff, :, :), 3, nSc).', ff)));
end
if min(mPx) < 4
    [mm, kw] = min(mPx);
    error('leg1_movie:containment', 'a spacecraft is %.1f px from the scene edge at day %.3f', ...
          mm, day(kw));
end
CAM = struct('az', az, 'el', el, 'target', tg, 'Rv', RvF, 'req', need, 'RvEnd', EF.RvEnd, ...
             'targetEnd', EF.target, 'positionEnd', EF.position, 'viewAngleEnd', EF.viewAngleEnd, ...
             'elEnd', EF.el, 'clipBox', clipBox, 'minMarginPx', min(mPx), ...
             'minMarginDay', day(find(mPx == min(mPx), 1)), 'zoomMax', max(RvF)/EF.RvEnd, ...
             'toPx', toPx);

% ---- what to render ----
pngMode = ~isempty(o.frameList) || ~isempty(o.dayList);
if pngMode
    fIdx = o.frameList(:);
    if any(fIdx < 1 | fIdx > nF | fIdx ~= round(fIdx))
        error('leg1_movie:frameList', 'frame indices must be integers in 1..%d', nF);
    end
    fDay = day(fIdx);
    dq = o.dayList(:);
    if any(~isfinite(dq) | dq < 0 | dq > dEnd + 1e-9)
        error('leg1_movie:dayList', 'dayList entries must lie in [0, %.4f] (mission days)', dEnd);
    end
    kfs = nF*ones(numel(dq), 1);
    for kq = 1:numel(dq)
        kf = find(day >= dq(kq) - 1e-12, 1);
        if ~isempty(kf), kfs(kq) = kf; end
    end
    fIdx = [fIdx; kfs];
    fDay = [fDay; dq];
    if ~isfolder(o.previewDir), mkdir(o.previewDir); end
elseif o.preview
    if strcmp(canon_path(o.outFile), canon_path(fullOut))
        error('leg1_movie:partialOutFile', 'a preview may not be written onto the full movie %s', fullOut);
    end
    nP   = round(o.previewSec*o.fps);
    fIdx = unique(round(linspace(1, nF, nP))).';
    fDay = day(fIdx);
else
    fIdx = (1:nF).';  fDay = day;
end
nR = numel(fIdx);

% ---- figure ----
fig = figure('Visible', onoff(o.visible), 'Units', 'pixels', 'Position', [20 20 W Hh], ...
             'Color', [0 0 0], 'MenuBar', 'none', 'ToolBar', 'none', 'Resize', 'off');
cleanFig = onCleanup(@() close_if(fig));
G = build_graphics(fig, scnPos, stripPos, clipBox, droKm, rMoonKm, rEarthKm, xEarth, ...
                   tArrD, tDepD, dEnd, nSc, o.moonTexture, o.earthTexture, rm.regime, unvetted);

vw = [];  partFile = '';
if ~pngMode
    [pDir, pName] = fileparts(o.outFile);
    if ~isfolder(pDir), mkdir(pDir); end
    partFile = fullfile(pDir, [pName '.partial.mp4']);
    vw = VideoWriter(partFile, 'MPEG-4');
    vw.FrameRate = o.fps;
    vw.Quality   = o.quality;
    open(vw);
end
cleanVw = onCleanup(@() close_vw(vw, partFile));

trD  = tArrD - tDepD;
caps = {'Six transfer stages leave one GTO', ...
        sprintf('%g N minimum-time transfers (%.1f-%.1f days each)', rm.regime.thrustN, ...
                min(trD), max(trD)), ...
        'Arrivals on the DRO slots', ...
        sprintf('All six on their DRO slots: hand-off at day %.2f', T0d)};
capD = [-inf max(tDepD) + 0.05 min(tArrD) dEnd - 1e-9 inf];

info = struct('outFile', '', 'nFrames', nR, 'nFramesMap', nF, 'frameSize', o.frameSize, ...
              'wallS', NaN, 'previewPngs', {cell(0, 1)}, 'frameIdx', fIdx, 'frameDay', fDay, ...
              'timeMap', TM, 'camera', CAM, 'posKm', nan(nR, 3, nSc), 'scMarginPx', nan(nR, 1), ...
              'Rv', nan(nR, 1), 'camPos', nan(nR, 3), 'camTarget', nan(nR, 3), ...
              'camViewAngle', nan(nR, 1), 'camModesManual', false(nR, 1), 'phase', zeros(nR, nSc), ...
              'plumeDir', nan(nR, 3, nSc), 'plumeDrawn', false(nR, nSc), ...
              'plumeScreenDir', nan(nR, 2, nSc), 'stackMax', zeros(nR, 1), ...
              'outlineDrawn', false(nR, 2), 'outlineInBox', false(nR, 2), ...
              'moonTextured', G.moonTextured, 'earthTextured', G.earthTextured, ...
              'T0Days', T0d, 'tDepDays', tDepD, 'tArrDays', tArrD, 'sgn', S.sgn, ...
              'orient0', orient0, 'gtoRotationErr', gtoRotErr, 'junctionKm', jKm, ...
              'dayText', {cell(nR, 1)}, 'regimeText', {G.regimeText}, 'unvetted', unvetted, ...
              'vetting', vetWhy, 'packageFile', L.identity.file);
inBox = @(Sx) Sx(:, 1) >= clipBox(1) & Sx(:, 1) <= clipBox(2) & ...
              Sx(:, 2) >= clipBox(3) & Sx(:, 2) <= clipBox(4);
spin = omegaE*tStar - 1;                    % Earth rotation in this frame [rad per ND]
stackEndD = max(tDepD) + 0.1;               % stack labels while the stack disperses [d]

for rr = 1:nR
    ff = fIdx(rr);  d0 = fDay(rr);
    tc = d0/dayNd;
    if pngMode && rr > numel(o.frameList)
        pos = zeros(nSc, 3);
        for kk = 1:nSc, pos(kk, :) = sat_pos(S, kk, tc); end
    else
        pos = reshape(posF(ff, :, :), 3, nSc).';
    end
    info.posKm(rr, :, :) = reshape(pos.', 1, 3, nSc);
    Rv = RvF(ff);
    kmPerPx = Rv/halfH;
    eh = EH(ff, :);  ev = EV(ff, :);  vdir = VD(ff, :);
    toPxR = @(X) toPx(X, ff);
    set(G.ax, 'CameraTarget', tg(ff, :), 'CameraPosition', tg(ff, :) + 1e6*vdir, ...
              'CameraUpVector', [0 0 1], 'CameraViewAngle', 2*atand(Rv/1e6));
    info.Rv(rr) = Rv;
    sp = toPxR(pos);
    info.scMarginPx(rr) = min(boxMargin(sp));
    % bodies and outlines
    set(G.earthT, 'Matrix', makehgtform('translate', xEarth, 'zrotate', spin*tc));
    gq = gtoAt(tc);
    set(G.gtoLine, 'XData', gq(:, 1), 'YData', gq(:, 2), 'ZData', gq(:, 3));
    place_label(G.gtoLbl, gq, 0.3, eh, ev, kmPerPx, toPxR, clipBox);
    place_label(G.droLbl, droKm, -0.4, eh, ev, kmPerPx, toPxR, clipBox);
    info.outlineDrawn(rr, :) = [strcmp(G.gtoLine.Visible, 'on'), strcmp(G.droLine.Visible, 'on')];
    info.outlineInBox(rr, :) = [any(inBox(toPxR(gq))), any(inBox(toPxR(droKm)))];
    update_scalebar(G, kmPerPx);
    % spacecraft
    nG = 0;  nT = 0;  nD = 0;
    for kk = 1:nSc
        ph = 1 + (d0 >= tDepD(kk)) + (d0 >= tArrD(kk));
        info.phase(rr, kk) = ph;
        p = pos(kk, :);
        set(G.sc(kk), 'XData', p(1), 'YData', p(2), 'ZData', p(3));
        switch ph
            case 1
                set(G.sc(kk), 'MarkerFaceColor', 'none', 'MarkerSize', 7, 'LineWidth', 1.6);
                set(G.ring(kk), 'Visible', 'off');
                nG = nG + 1;
            case 2
                set(G.sc(kk), 'MarkerFaceColor', G.col(kk, :), 'MarkerSize', 8, 'LineWidth', 1);
                set(G.ring(kk), 'Visible', 'off');
                nT = nT + 1;
            otherwise
                set(G.sc(kk), 'MarkerFaceColor', G.col(kk, :), 'MarkerSize', 7, 'LineWidth', 1);
                set(G.ring(kk), 'Visible', 'on', 'XData', p(1), 'YData', p(2), 'ZData', p(3));
                nD = nD + 1;
        end
        % trail
        tq = linspace(max(0, tc - o.tailDays/dayNd), tc, 60).';
        tp = sat_pos(S, kk, tq);
        al = 0.85*((1:size(tp, 1)).'/size(tp, 1)).^1.6;
        % a NaN vertex ends the polyline, so the patch does not close it
        set(G.trail(kk), 'Vertices', [tp; nan(1, 3)], 'Faces', 1:size(tp, 1) + 1, ...
                         'FaceVertexAlphaData', [al; 0]);
        % plume
        if ph == 2
            a = sat(kk);
            al3 = interp1(a.t, a.alpha, min(max(tc, a.t(1)), a.t(end)), 'linear');
            e = -al3/sqrt(sum(al3.^2));               % exhaust, opposite the thrust
            info.plumeDir(rr, :, kk) = e;
            nv = cross(e, vdir);
            if sqrt(sum(nv.^2)) < 0.05, nv = eh; end
            nv = nv/sqrt(sum(nv.^2));
            c0 = p + 4*kmPerPx*e;
            Lkm = o.plumePx*kmPerPx;
            Sh = G.plumeShape;
            Vg = plume_poly(c0, e, nv, Sh.glowL*Lkm, Sh.waist, Sh.glowW*kmPerPx);
            Vc = plume_poly(c0, e, nv, Lkm, Sh.waist, Sh.coreW*kmPerPx);
            spl = toPxR([Vc(3, :); (Vc(1, :) + Vc(5, :))/2]);
            info.plumeScreenDir(rr, :, kk) = spl(1, :) - spl(2, :);
            set(G.plume{kk}(1), 'Vertices', Vg, 'Visible', 'on');
            set(G.plume{kk}(2), 'Vertices', Vc, 'Visible', 'on');
            info.plumeDrawn(rr, kk) = true;
        else
            set(G.plume{kk}, 'Visible', 'off');
        end
    end
    % stacks: markers within stackPx of one another carry a count label
    grp = zeros(nSc, 1);  ng = 0;
    for kk = 1:nSc
        if grp(kk) > 0, continue; end
        ng = ng + 1;
        near = grp == 0 & sqrt(sum((sp - sp(kk, :)).^2, 2)) < stackPx;
        grp(near) = ng;
    end
    set(G.stackTxt, 'Visible', 'off');
    placed = nan(nSc, 2);  nPl = 0;             % label anchors already used [px]
    for gg = 1:ng
        mem = find(grp == gg);
        info.stackMax(rr) = max(info.stackMax(rr), numel(mem));
        if numel(mem) > 1 && d0 < stackEndD
            lp = [sp(mem(1), 1) + 9, sp(mem(1), 2) + 11];
            while any(abs(placed(:, 1) - lp(1)) < 70 & abs(placed(:, 2) - lp(2)) < 16)
                lp(2) = lp(2) - 16;             % step down past a label already there
            end
            nPl = nPl + 1;  placed(nPl, :) = lp;
            set(G.stackTxt(gg), 'Position', [lp 0], ...
                'String', sprintf('%d stacked', numel(mem)), 'Visible', 'on');
        end
    end
    info.camPos(rr, :) = G.ax.CameraPosition;
    info.camTarget(rr, :) = G.ax.CameraTarget;
    info.camViewAngle(rr) = G.ax.CameraViewAngle;
    info.camModesManual(rr) = all(strcmp({G.ax.CameraPositionMode, G.ax.CameraTargetMode, ...
                                          G.ax.CameraViewAngleMode, G.ax.CameraUpVectorMode}, 'manual'));
    % text
    set(G.txtDay, 'String', sprintf('Day %5.2f', d0));
    info.dayText{rr} = char(G.txtDay.String);          % what the counter DISPLAYS
    set(G.txtGto, 'String', sprintf('On the GTO  %d', nG));
    set(G.txtThr, 'String', sprintf('Thrusting   %d', nT));
    set(G.txtDro, 'String', sprintf('In the DRO  %d / %d', nD, nSc));
    rNow = rate(ff);
    if rNow > 0 && ~(pngMode && rr > numel(o.frameList) && d0 >= dEnd)
        set(G.txtRate, 'String', sprintf('%4.3f days per second', rNow));
    else
        set(G.txtRate, 'String', '');
    end
    kc = find(capD(2:end) > d0, 1);
    set(G.cap, 'String', caps{kc});
    % strip
    [sx, sy] = on_station_steps(tArrD, d0);
    set(G.nLine, 'XData', sx, 'YData', sy);
    set(G.cursor, 'XData', [d0 d0]);
    for kk = 1:nSc
        set(G.arrTick(kk), 'Visible', onoff(d0 >= tArrD(kk)));
        set(G.depTick(kk), 'Visible', onoff(d0 >= tDepD(kk)));
    end
    drawnow;
    F = getframe(fig);
    img = F.cdata;
    if size(img, 1) ~= o.frameSize(2) || size(img, 2) ~= o.frameSize(1)
        img = imresize(img, [o.frameSize(2) o.frameSize(1)]);
    end
    if pngMode
        pf = fullfile(o.previewDir, sprintf('%s_%04d_day%05.3f.png', o.pngPrefix, ff, d0));
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
        error('leg1_movie:verify', 'written video failed its check: %s', partFile);
    end
    movefile(partFile, o.outFile, 'f');
    info.outFile = o.outFile;
end
clear cleanVw cleanFig
info.wallS = toc(tWall);
end

%% ------------------------------------------------------------------------
function L = load_inputs(catalogMat, rm, here, leg1Mat)
% LOAD_INPUTS  leg1_inputs on the package the leg-1 run used, found by its
%   SHA-256 BEFORE anything is loaded (leg1_package_for): the opts value
%   (else the file named in run_meta), else run_meta.packageSnapshot, else a
%   snapshot in inputs/ beside the results or in results/inputs
%   (leg1_movie:identity when none matches).
%   INPUTS:  catalogMat char ('' = run_meta's file); rm run_meta struct;
%            here char; leg1Mat char.
%   OUTPUTS: L leg1_inputs struct.
snapX = {};
if isfield(rm, 'packageSnapshot') && ~isempty(rm.packageSnapshot), snapX = {char(rm.packageSnapshot)}; end
inDirs = unique({fullfile(fileparts(canon_path(leg1Mat)), 'inputs'), fullfile(here, 'results', 'inputs')}, 'stable');
[pkg, why] = leg1_package_for(rm.identity, struct('primary', catalogMat, 'extra', {snapX}, 'inputsDirs', {inDirs}));
if isempty(pkg)
    error('leg1_movie:identity', '%s', why);
end
L = leg1_inputs(struct('catalogMat', pkg));
if ~strcmp(L.identity.sha256, rm.identity.sha256)
    error('leg1_movie:identity', 'package %s changed while it was being loaded', pkg);
end
end

function why = vet_results(P1, rm)
% VET_RESULTS  '' when the leg-1 results may be rendered: run_meta.verdict
%   exactly 'PASS', six transfers all accepted with G.pass true, and every
%   test seam recorded in run_meta.knobs empty / zero; else the reasons.
%   INPUTS:  P1 struct array; rm run_meta struct.
%   OUTPUTS: why char.
r = {};
v = '';
if isfield(rm, 'verdict'), v = char(rm.verdict); end
if ~strcmp(v, 'PASS'), r{end+1} = sprintf('verdict "%s", not PASS', v); end
nAcc = 0;  nPass = 0;
if isfield(P1, 'accepted'), nAcc = nnz([P1.accepted]); end
if isfield(P1, 'G')
    nPass = nnz(arrayfun(@(q) isstruct(q.G) && isfield(q.G, 'pass') && isequal(q.G.pass, true), P1));
end
if ~(numel(P1) == 6 && nAcc == 6 && nPass == 6)
    r{end+1} = sprintf('%d transfers, %d accepted, %d with G.pass', numel(P1), nAcc, nPass);
end
if ~isfield(rm, 'knobs')
    r{end+1} = 'no knobs recorded in run_meta (test seams unknown)';
else
    kn = rm.knobs;
    seamNm  = {'seamRefuse', 'seamRefute', 'seamFinalT0'};
    seamOn  = [isfield(kn, 'seamRefuse') && ~isempty(kn.seamRefuse), ...
               isfield(kn, 'seamRefute') && ~isempty(kn.seamRefute), ...
               isfield(kn, 'seamFinalT0') && ~isequal(kn.seamFinalT0, 0)];
    if any(seamOn), r{end+1} = sprintf('test seams set: %s', strjoin(seamNm(seamOn), ', ')); end
end
why = strjoin(r, '; ');
end

function EF = end_framing(leg2Mat)
% END_FRAMING  mission_movie's first frame (view, target, scene size), read
%   from its own info by a one-still render into a temporary folder.
%   INPUTS:  leg2Mat char leg-2 results.
%   OUTPUTS: EF struct .az .el .target .position .viewAngleEnd .RvEnd.
if ~isfile(leg2Mat)
    error('leg1_movie:noLeg2', 'leg-2 results not found (needed for the end framing): %s', leg2Mat);
end
tmp = tempname;
mkdir(tmp);
cleanTmp = onCleanup(@() rmdir(tmp, 's'));
I = mission_movie(leg2Mat, struct('dayList', 0, 'previewDir', tmp));
C = I.camera;
EF = struct('az', C.az, 'el', C.el, 'target', C.target, 'position', C.position, ...
            'viewAngleEnd', C.viewAngleEnd, 'RvEnd', C.RvEnd);
clear cleanTmp
end

function w = raised_ramp(t, t0, t1)
% RAISED_RAMP  0 before t0, 1 after t1, raised cosine between.
%   INPUTS: t [n x 1]; t0, t1 scalars.  OUTPUTS: w [n x 1].
w = (1 - cos(pi*min(max((t - t0)/max(t1 - t0, eps), 0), 1)))/2;
end

function r = rotate_earth(X, dDeg, mu)
% ROTATE_EARTH  Positions of X rotated about the Earth's z-axis.
%   INPUTS:  X [>=3 x n] ND barycentric rotating states; dDeg scalar or
%            [1 x n] angle [deg]; mu scalar.
%   OUTPUTS: r [3 x n] rotated positions [ND].
c = cosd(dDeg);  s = sind(dDeg);
rx = X(1, :) + mu;  ry = X(2, :);
r = [c.*rx - s.*ry - mu; s.*rx + c.*ry; X(3, :)];
end

function R = gto_pos(S, tq)
% GTO_POS  Positions on the one inertial GTO (shared stack phase) at leg-1
%   times tq.  INPUTS: S struct; tq [n x 1] ND.  OUTPUTS: R [n x 3] ND.
tq = tq(:).';
ph = mod(S.ph0 + tq/S.tauG, 1);
x  = S.gtoFun(S.oRef, ph);
R  = rotate_earth(x, S.orient0 + S.sgn*(180/pi)*tq - S.oRef, S.mu).';
end

function R = sat_pos(S, kk, tq)
% SAT_POS  Moon-centred position of stage kk at leg-1 times tq: GTO before
%   tDep, flown history (cubic Hermite) to tArr, DRO after.
%   INPUTS:  S struct; kk index; tq [n x 1] ND.
%   OUTPUTS: R [n x 3] km.
a  = S.sat(kk);
tq = tq(:);
R  = zeros(numel(tq), 3);
b1 = tq < a.tDep;  b2 = tq >= a.tDep & tq <= a.tArr;  b3 = tq > a.tArr;
if any(b1), R(b1, :) = gto_pos(S, tq(b1)); end
if any(b2), R(b2, :) = hermite3(a.t, a.X, tq(b2)); end
if any(b3)
    x = S.slotFun(mod(a.phArr + (tq(b3).' - a.tArr)/S.tauD, 1));
    R(b3, :) = x(1:3, :).';
end
R = (R - [1 - S.mu 0 0])*S.lStar;
end

function P = hermite3(t, X, tq)
% HERMITE3  Cubic Hermite interpolation of position with velocity.
%   INPUTS:  t [n x 1] grid; X [6 x n] states; tq [m x 1] in [t(1), t(end)].
%   OUTPUTS: P [m x 3].
k  = discretize(tq, t);
k(isnan(k)) = numel(t) - 1;
k  = min(k, numel(t) - 1);
h  = t(k + 1) - t(k);
s  = (tq - t(k))./h;
p0 = X(1:3, k).';  p1 = X(1:3, k + 1).';
v0 = X(4:6, k).'.*h;  v1 = X(4:6, k + 1).'.*h;
P  = (2*s.^3 - 3*s.^2 + 1).*p0 + (s.^3 - 2*s.^2 + s).*v0 + (-2*s.^3 + 3*s.^2).*p1 + (s.^3 - s.^2).*v1;
end

%% ------------------------------------------------------------------------
function G = build_graphics(fig, scnPos, stripPos, clipBox, droKm, rMoonKm, rEarthKm, xEarth, ...
                            tArrD, tDepD, dEnd, nSc, moonTex, earthTex, regime, unvetted)
% BUILD_GRAPHICS  Create every graphics object once; the loop only updates.
%   INPUTS:  fig figure; scnPos, stripPos [1 x 4] px; clipBox [1 x 4] px;
%            droKm [m x 3] km; rMoonKm, rEarthKm scalars; xEarth [1 x 3] km;
%            tArrD, tDepD [1 x 6] d; dEnd scalar d; nSc scalar; moonTex,
%            earthTex char texture files; regime struct (thrustN ispS m0kg);
%            unvetted logical (stamp UNVETTED RESULTS in the regime text).
%   OUTPUTS: G struct of handles and styles.
fg   = [0.93 0.94 0.97];
fg2  = [0.70 0.73 0.80];
mono = 'Menlo';
if ~any(strcmp(listfonts, mono)), mono = 'FixedWidth'; end
sans = 'Helvetica';
G.col = [0.40 0.80 0.93;      % cyan     (mission_movie's spacecraft colours,
         0.93 0.30 0.33;      % red       same numbering k as leg 2)
         0.30 0.78 0.40;      % green
         0.36 0.58 1.00;      % blue
         0.72 0.50 1.00;      % violet
         1.00 0.55 0.80];     % pink
% scene
ax = axes(fig, 'Units', 'pixels', 'Position', scnPos, 'Color', 'none', ...
          'XLim', [-5e5 5e5], 'YLim', [-5e5 5e5], 'ZLim', [-5e5 5e5], ...
          'DataAspectRatio', [1 1 1], 'Projection', 'orthographic', 'Clipping', 'off', ...
          'CameraUpVector', [0 0 1], 'Visible', 'off', 'SortMethod', 'depth');
hold(ax, 'on');
set(ax, 'CameraPositionMode', 'manual', 'CameraTargetMode', 'manual', ...
        'CameraViewAngleMode', 'manual', 'CameraUpVectorMode', 'manual');
% Moon: mission_movie's construction (pumpkyn.util.moon3D's mesh and map)
[xs, ys, zs] = ellipsoid(0, 0, 0, rMoonKm, rMoonKm, rMoonKm, 100);
G.moon = surf(ax, -xs, -ys, -zs, 'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
              'AmbientStrength', 0.65, 'DiffuseStrength', 0.75, 'SpecularStrength', 0, ...
              'BackFaceLighting', 'unlit');
if isfile(moonTex)
    tex = imread(moonTex);
    if size(tex, 3) == 1, tex = repmat(tex, 1, 1, 3); end
    set(G.moon, 'FaceColor', 'texturemap', 'CData', tex);
else
    set(G.moon, 'FaceColor', [0.66 0.66 0.68]);
    fprintf('leg1_movie: Moon texture not found (%s); drawing a shaded grey sphere\n', moonTex);
end
G.moonTextured = ischar(G.moon.FaceColor) && strcmp(G.moon.FaceColor, 'texturemap') && ...
                 size(G.moon.CData, 3) == 3;
% Earth: same mesh convention, centred by an hgtransform that also spins it
G.earthT = hgtransform(ax, 'Matrix', makehgtform('translate', xEarth));
[xs, ys, zs] = ellipsoid(0, 0, 0, rEarthKm, rEarthKm, rEarthKm, 100);
G.earth = surf(ax, -xs, -ys, -zs, 'Parent', G.earthT, 'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
               'AmbientStrength', 0.70, 'DiffuseStrength', 0.75, 'SpecularStrength', 0.1, ...
               'BackFaceLighting', 'unlit');
if isfile(earthTex)
    tex = imread(earthTex);
    tex = tex(1:4:end, 1:4:end, :);
    if size(tex, 3) == 1, tex = repmat(tex, 1, 1, 3); end
    set(G.earth, 'FaceColor', 'texturemap', 'CData', tex);
else
    set(G.earth, 'FaceColor', [0.25 0.45 0.85]);
    fprintf('leg1_movie: Earth texture not found (%s); drawing a shaded blue sphere\n', earthTex);
end
G.earthTextured = ischar(G.earth.FaceColor) && strcmp(G.earth.FaceColor, 'texturemap') && ...
                  size(G.earth.CData, 3) == 3;
light(ax, 'Style', 'infinite', 'Position', [0.35 -1 0.45]);
% outlines (always on)
G.droRGB = [0.45 0.62 0.95];
G.gtoRGB = [0.98 0.80 0.45];
G.droLine = plot3(ax, droKm(:, 1), droKm(:, 2), droKm(:, 3), '-', 'Color', [G.droRGB 0.70], 'LineWidth', 1.1);
G.gtoLine = plot3(ax, NaN, NaN, NaN, '-', 'Color', [G.gtoRGB 0.60], 'LineWidth', 1.1);
G.droLbl = text(ax, 0, 0, 0, 'DRO', 'Color', [0.55 0.70 1.0], 'FontName', sans, 'FontSize', 12, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
G.gtoLbl = text(ax, 0, 0, 0, 'GTO', 'Color', G.gtoRGB, 'FontName', sans, 'FontSize', 12, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
% plume: mission_movie's two tapered patches
G.plumeShape = struct('waist', 0.30, 'coreW', [9 6 0], 'glowW', [20 14 0], 'glowL', 0.80);
coreC = [1.00 1.00 0.80; 1.00 0.86 0.30; 1.00 0.42 0.05; 1.00 0.86 0.30; 1.00 1.00 0.80];
coreA = [1; 0.85; 0; 0.85; 1];
glowC = repmat([1.00 0.78 0.25], 5, 1);
glowA = [0.30; 0.18; 0; 0.18; 0.30];
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

% mask: black everywhere outside clipBox, above the scene, below the rest
mk = axes(fig, 'Units', 'pixels', 'Position', [0 0 1280 720], 'XLim', [0 1280], 'YLim', [0 720], ...
          'Visible', 'off', 'Color', 'none');
hold(mk, 'on');
x0 = clipBox(1);  x1 = clipBox(2);  y0 = clipBox(3);  y1 = clipBox(4);
rects = [0 1280 0 y0; 0 1280 y1 720; 0 x0 y0 y1; x1 1280 y0 y1];
for kr = 1:4
    patch(mk, rects(kr, [1 2 2 1]), rects(kr, [3 3 4 4]), [0 0 0], 'EdgeColor', 'none');
end

% strip: satellites in the DRO versus mission day
xEndS = dEnd*1.04;
axN = axes(fig, 'Units', 'pixels', 'Position', stripPos, 'Color', [0.05 0.06 0.09], ...
           'XLim', [0 xEndS], 'YLim', [0 6.6], 'YTick', 0:2:6, 'XTick', 0:0.5:floor(xEndS*2)/2, ...
           'XColor', fg2, 'YColor', [0.55 0.72 1.0], 'FontName', sans, 'FontSize', 11, ...
           'Box', 'on', 'Layer', 'top', 'TickDir', 'out', 'GridColor', [1 1 1], 'GridAlpha', 0.10, ...
           'XGrid', 'on', 'YGrid', 'on');
hold(axN, 'on');
xlabel(axN, 'mission day', 'Color', fg2, 'FontSize', 12);
ylabel(axN, 'satellites in the DRO', 'Color', [0.55 0.72 1.0], 'FontSize', 12);
for kk = 1:nSc
    G.arrTick(kk) = plot(axN, [tArrD(kk) tArrD(kk)], [0 1.0], '-', 'Color', G.col(kk, :), ...
                         'LineWidth', 2, 'Visible', 'off');
    G.depTick(kk) = plot(axN, tDepD(kk), 6.25, 'v', 'Color', G.col(kk, :), ...
                         'MarkerFaceColor', G.col(kk, :), 'MarkerSize', 5, 'Visible', 'off');
end
plot(axN, [dEnd dEnd], [0 6.6], '--', 'Color', [1 0.80 0.45 0.9], 'LineWidth', 1.2);
text(axN, min(tArrD) - 0.04, 4.6, sprintf('hand-off T0 = day %.2f', dEnd), 'Color', [1 0.80 0.45], ...
     'FontName', sans, 'FontSize', 10, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
text(axN, max(tDepD) + 0.04, 6.25, 'departures', 'Color', fg2, 'FontName', sans, 'FontSize', 10, ...
     'VerticalAlignment', 'middle');
G.nLine = plot(axN, NaN, NaN, '-', 'Color', [0.55 0.72 1.0], 'LineWidth', 1.8);
G.cursor = plot(axN, [0 0], [0 6.6], '-', 'Color', [1 1 1 0.55], 'LineWidth', 1);

% overlay in figure pixels
ov = axes(fig, 'Units', 'pixels', 'Position', [0 0 1280 720], 'XLim', [0 1280], 'YLim', [0 720], ...
          'Visible', 'off', 'Color', 'none');
hold(ov, 'on');
text(ov, 640, 700, sprintf('Lunar PNT constellation deployment, leg 1: GTO to DRO, %g N transfer stages', ...
     regime.thrustN), 'Color', fg, 'FontName', sans, 'FontSize', 20, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
G.cap = text(ov, 640, 668, '', 'Color', fg2, 'FontName', sans, 'FontSize', 15, ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
G.txtDay  = text(ov, 28, 612, '', 'Color', fg, 'FontName', mono, 'FontSize', 20, 'FontWeight', 'bold');
G.txtGto  = text(ov, 28, 580, '', 'Color', fg, 'FontName', mono, 'FontSize', 14);
G.txtThr  = text(ov, 28, 556, '', 'Color', fg, 'FontName', mono, 'FontSize', 14);
G.txtDro  = text(ov, 28, 532, '', 'Color', fg, 'FontName', mono, 'FontSize', 14);
G.txtRate = text(ov, 28, 506, '', 'Color', [0.55 0.58 0.65], 'FontName', mono, 'FontSize', 11);
G.regimeText = {sprintf('%g N per %g kg, Isp %g s,', regime.thrustN, regime.m0kg, regime.ispS), ...
                'engine on throughout:', 'catalog regime,', 'not a vehicle design'};
if unvetted
    G.regimeText = [{'UNVETTED RESULTS'}, G.regimeText];
    text(ov, 28, 438, 'UNVETTED RESULTS', 'Color', [1 0.35 0.35], 'FontName', sans, 'FontSize', 13, ...
         'FontWeight', 'bold', 'VerticalAlignment', 'top');
end
text(ov, 28, 420, G.regimeText(1 + unvetted:end), ...
     'Color', [0.62 0.64 0.70], 'FontName', sans, 'FontSize', 11, 'VerticalAlignment', 'top');
% state legend (right)
x0 = 1062;  y0 = 612;  dy = 25;
text(ov, x0 - 8, y0 + 4, 'spacecraft state', 'Color', fg2, 'FontName', sans, 'FontSize', 11);
yy = y0 - dy;
plot(ov, x0 + 8, yy, 'o', 'Color', fg, 'MarkerSize', 7, 'LineWidth', 1.6);
text(ov, x0 + 26, yy, 'on the GTO', 'Color', fg, 'FontName', sans, 'FontSize', 12);
yy = yy - dy;
plot(ov, [x0 - 14 x0 + 2], [yy yy], '-', 'Color', [1 0.66 0.15 0.7], 'LineWidth', 3);
plot(ov, [x0 - 6 x0 + 4], [yy yy], '-', 'Color', [1 0.95 0.6], 'LineWidth', 5);
plot(ov, x0 + 8, yy, 'o', 'Color', fg, 'MarkerFaceColor', fg, 'MarkerSize', 8);
text(ov, x0 + 26, yy, 'thrusting (plume)', 'Color', fg, 'FontName', sans, 'FontSize', 12);
yy = yy - dy;
plot(ov, x0 + 8, yy, 'o', 'Color', fg, 'MarkerSize', 14, 'LineWidth', 0.8);
plot(ov, x0 + 8, yy, 'o', 'Color', fg, 'MarkerFaceColor', fg, 'MarkerSize', 7);
text(ov, x0 + 26, yy, 'in the DRO', 'Color', fg, 'FontName', sans, 'FontSize', 12);
yy = yy - dy - 6;
text(ov, x0 - 8, yy, 'spacecraft', 'Color', fg2, 'FontName', sans, 'FontSize', 11);
for kk = 1:nSc
    text(ov, x0 + 84 + 20*(kk - 1), yy, sprintf('%d', kk), 'Color', G.col(kk, :), ...
         'FontName', sans, 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
yy = yy - dy - 6;
plot(ov, [x0 - 6 x0 + 20], [yy yy], '-', 'Color', [G.gtoRGB 0.8], 'LineWidth', 1.4);
text(ov, x0 + 26, yy, 'GTO (inertial)', 'Color', fg, 'FontName', sans, 'FontSize', 12);
yy = yy - dy;
plot(ov, [x0 - 6 x0 + 20], [yy yy], '-', 'Color', [G.droRGB 0.9], 'LineWidth', 1.4);
text(ov, x0 + 26, yy, 'DRO', 'Color', fg, 'FontName', sans, 'FontSize', 12);
G.stackTxt = gobjects(1, nSc);
for kk = 1:nSc
    G.stackTxt(kk) = text(ov, 0, 0, '', 'Color', fg, 'FontName', sans, 'FontSize', 11, ...
                          'Visible', 'off');
end
% scale bar and frame label
G.sbLine = plot(ov, [0 1], [212 212], '-', 'Color', fg2, 'LineWidth', 2);
G.sbEnds = plot(ov, [0 0 NaN 1 1], [207 217 NaN 207 217], '-', 'Color', fg2, 'LineWidth', 1.2);
G.sbTxt  = text(ov, 28, 228, '', 'Color', fg2, 'FontName', mono, 'FontSize', 11);
text(ov, 1252, 212, 'Earth-Moon rotating frame', 'Color', [0.55 0.58 0.65], 'FontName', sans, ...
     'FontSize', 11, 'HorizontalAlignment', 'right');
G.ov = ov;
end

function place_label(h, X, tilt, eh, ev, kmPerPx, toPxR, clipBox)
% PLACE_LABEL  Put an outline label beside the outline's right-most point
%   (tilted by tilt toward screen up), or beside its left-most point,
%   right-aligned, when the right one would leave the scene area.
%   INPUTS:  h text handle; X [m x 3] outline [km]; tilt scalar; eh, ev
%            [1 x 3] screen axes; kmPerPx scalar; toPxR handle km -> px;
%            clipBox [1 x 4] px.
%   OUTPUTS: none (updates h).
[~, kr] = max(X*eh.' + tilt*(X*ev.'));
pr = toPxR(X(kr, :));
if pr(1) + 12 + 50 <= clipBox(2)
    set(h, 'Position', X(kr, :) + 12*kmPerPx*eh, 'HorizontalAlignment', 'left');
else
    [~, kl] = min(X*eh.' - tilt*(X*ev.'));
    set(h, 'Position', X(kl, :) - 12*kmPerPx*eh, 'HorizontalAlignment', 'right');
end
end

function update_scalebar(G, kmPerPx)
% UPDATE_SCALEBAR  The largest round length whose bar is at most 160 px.
%   INPUTS: G handles; kmPerPx scalar.  OUTPUTS: none (updates G).
nice = [1000 2000 5000 10000 20000 50000 100000 200000];
kb  = max([1, find(nice/kmPerPx <= 160, 1, 'last')]);
Lpx = nice(kb)/kmPerPx;
x0 = 28;
set(G.sbLine, 'XData', [x0 x0 + Lpx]);
set(G.sbEnds, 'XData', [x0 x0 NaN x0 + Lpx x0 + Lpx]);
set(G.sbTxt, 'String', sprintf('%-10s', [thousands(nice(kb)) ' km']));
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

function [x, y] = on_station_steps(tArrD, d0)
% ON_STATION_STEPS  Step curve of satellites in the DRO from day 0 to d0.
%   INPUTS: tArrD [1 x K] arrival days; d0 scalar day.  OUTPUTS: x, y [1 x m].
ta = sort(tArrD(tArrD <= d0));
n  = numel(ta);
x  = [0, reshape([ta; ta], 1, []), d0];
y  = [0, reshape([0:n - 1; 1:n], 1, []), n];
end

function y = unimodal_hull(x)
% UNIMODAL_HULL  Smallest single-peaked sequence above x.
%   INPUTS: x [n x 1].  OUTPUTS: y [n x 1].
y = min(cummax(x), flipud(cummax(flipud(x))));
end

function y = hann_smooth(x, h)
% HANN_SMOOTH  Raised-cosine moving average of half-width h samples, edges
%   padded with the end values.  INPUTS: x [n x 1]; h integer.
%   OUTPUTS: y [n x 1].
k  = 0.5*(1 + cos(pi*(-h:h).'/(h + 1)));
k  = k/sum(k);
xp = [repmat(x(1), h, 1); x(:); repmat(x(end), h, 1)];
y  = conv(xp, k, 'valid');
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
% CLOSE_VW  Close a VideoWriter if open and delete an unfinished partial file.
%   INPUTS: vw VideoWriter or []; partFile char ('' = nothing to delete).
if ~isempty(vw) && isvalid(vw), close(vw); end
if ~isempty(partFile) && isfile(partFile), delete(partFile); end
end
