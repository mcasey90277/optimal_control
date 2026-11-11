function [ok, skipped] = test_leg1_movie()
%% Purpose:
%
%   Verifies leg1_movie, mission_movie's dayOffset and full_mission_movie on
%   the SAVED leg-1 and leg-2 study results (nothing re-solved):
%   (1) PREVIEW: a 2 s preview MP4 (opts.preview) into a temporary folder
%   exists, is 1280 x 720 and has exactly 2*fps frames (VideoReader); a
%   preview aimed at the full movie's file is refused
%   (leg1_movie:partialOutFile); (2) TIME MAP: monotone, starts at day 0,
%   reaches the last arrival and ends at T0; (3) IN VIEW: at 20 map frames
%   spread over the movie, every spacecraft's projected position (computed
%   here from info.camera.toPx and info.posKm) lies inside the scene box
%   info.camera.clipBox, and the map-wide minimum margin is >= 4 px;
%   (4) POSITIONS, independently: mid-flight drawn positions equal the saved
%   flown history (interp1 spline of P1.X, < 1 km), a pre-departure position
%   equals L.gtoState at the drifted orientation (< 1e-6 km), a post-arrival
%   one equals L.slotState at the advanced phase; (5) PLUME: only while
%   thrusting; its direction is -alpha/|alpha| of the saved flight at a node
%   (1e-12) and the drawn plume points along the projected exhaust;
%   (6) CAMERA: all modes manual; the LAST frame equals mission_movie's
%   FIRST frame (target, position, view angle to 1e-9); (7) dayOffset: a
%   mission_movie still with dayOffset = 3 reports a first displayed day 3
%   greater than its mission day, the default reports it unchanged;
%   (8) JOIN: full_mission_movie on two short movies writes n1 + n2 - fade
%   frames at 1280 x 720, and refuses a second movie of another frame size
%   (full_mission_movie:mismatch); (9) DELIBERATE FAILURES: zoom forced to
%   half (zoomScale seam) -> leg1_movie:containment; GTO turned the wrong
%   way would put a pre-departure stage off its GTO (checked as a number);
%   a missing results file -> leg1_movie:noResults; dayList beyond T0 ->
%   leg1_movie:dayList. Returns skipped = true (ok = true) when a results
%   file is missing.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 true if every assertion
%                                                   held (else it throws)
%  skipped                  logical                 true if a results file
%                                                   is missing
%
%% Revision History:
%  M. Casey                                                   (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_lunar_pnt_paths.m'));
resDir  = fullfile(fileparts(here), 'results');
leg1Mat = fullfile(resDir, 'lunar_pnt_leg1_study.mat');
leg2Mat = fullfile(resDir, 'lunar_pnt_mission_study.mat');
skipped = false;
if ~isfile(leg1Mat) || ~isfile(leg2Mat)
    fprintf('  SKIPPED: leg-1 or leg-2 results missing; run both studies first\n');
    ok = true;  skipped = true;
    return
end
tmp = tempname;
mkdir(tmp);
cleanTmp = onCleanup(@() rmdir(tmp, 's'));
fullMp4 = fullfile(resDir, 'lunar_pnt_leg1.mp4');
fullInfo = dir(fullMp4);

% ---- 1. preview MP4 ----
mp4 = fullfile(tmp, 'leg1_preview.mp4');
I = leg1_movie(leg1Mat, struct('preview', true, 'outFile', mp4));
assert(isfile(mp4) && ~isfile(fullfile(tmp, 'leg1_preview.partial.mp4')), 'preview MP4 not written cleanly');
vr = VideoReader(mp4);
fprintf('  preview MP4: %d x %d, %d frames at %g fps, wall %.1f s\n', vr.Width, vr.Height, ...
        vr.NumFrames, vr.FrameRate, I.wallS);
assert(vr.Width == 1280 && vr.Height == 720, 'preview size %d x %d', vr.Width, vr.Height);
assert(mod(vr.Width, 16) == 0 && mod(vr.Height, 16) == 0, 'frame size not a multiple of 16');
assert(vr.NumFrames == 60 && I.nFrames == 60, 'preview frames %d (info %d), expected 60', ...
       vr.NumFrames, I.nFrames);
n1 = vr.NumFrames;
clear vr
thrown = '';
try
    leg1_movie(leg1Mat, struct('preview', true, 'outFile', fullMp4));
catch err
    thrown = err.identifier;
end
assert(strcmp(thrown, 'leg1_movie:partialOutFile'), 'preview onto the full movie not refused (%s)', thrown);
if ~isempty(fullInfo)
    nowInfo = dir(fullMp4);
    assert(nowInfo.datenum == fullInfo.datenum && nowInfo.bytes == fullInfo.bytes, 'full movie touched');
end
fprintf('  preview onto the full movie refused (leg1_movie:partialOutFile)\n');

% ---- 2. time map ----
TM = I.timeMap;
assert(all(diff(TM.day) >= 0), 'time map not monotone');
assert(TM.day(1) == 0 && abs(TM.day(end) - I.T0Days) < 1e-12, 'time map ends %.6f, T0 %.6f', ...
       TM.day(end), I.T0Days);
assert(max(TM.day) >= max(I.tArrDays), 'time map does not reach the last arrival');
fprintf(['  time map: %d frames, monotone, day 0 -> T0 %.4f d (last arrival %.4f d); budget [s] %s; ', ...
         'rate %.3f..%.3f d/s\n'], numel(TM.day), I.T0Days, max(I.tArrDays), mat2str(TM.budgetS, 3), ...
        min(TM.rate(TM.rate > 0)), max(TM.rate));

% ---- 3. in view at 20 map frames (positions recomputed below) ----
C  = I.camera;
nF = I.nFramesMap;
D  = load(leg1Mat, 'P1', 'T0');
P1 = D.P1;  [~, ord] = sort([P1.k]);  P1 = P1(ord);
L  = leg1_inputs(struct('catalogMat', fullfile(resDir, 'inputs', ...
                 'costate_catalog_gto_dro_pilot_e57d9236.mat')));
mu = L.muStar;  lStar = L.lStar;  dayNd = L.tStar/86400;
ff20 = round(linspace(1, nF, 20));
J = leg1_movie(leg1Mat, struct('frameList', ff20, 'previewDir', tmp, 'pngPrefix', 't'));
minIn = inf;
for rr = 1:20
    X  = reshape(J.posKm(rr, :, :), 3, 6).';
    sp = C.toPx(X, ff20(rr));
    m  = min([sp(:, 1) - C.clipBox(1), C.clipBox(2) - sp(:, 1), sp(:, 2) - C.clipBox(3), ...
              C.clipBox(4) - sp(:, 2)], [], 'all');
    minIn = min(minIn, m);
end
assert(minIn > 0 && C.minMarginPx >= 4, 'a spacecraft outside the scene box (%.1f px)', minIn);
assert(all(J.camModesManual), 'a camera mode is not manual');
fprintf('  in view: 20 sampled frames, min margin %.1f px (map-wide %.1f px at day %.3f); zoom max x%.1f\n', ...
        minIn, C.minMarginPx, C.minMarginDay, C.zoomMax);

% ---- 4. drawn positions vs the saved data, independently ----
kS = 2;  p = P1(kS);
tMid = p.tDep + 0.4*(p.tArr - p.tDep);
tPre = 0.5*p.tDep;                                    % stage 2 still on the GTO
tPost = p.tArr + 0.5*(D.T0 - p.tArr);
K = leg1_movie(leg1Mat, struct('dayList', [tPre tMid tPost]*dayNd, 'previewDir', tmp, 'pngPrefix', 't'));
xMid = interp1(p.t(:), p.X(1:3, :).', tMid, 'spline');
eMid = sqrt(sum(((xMid - [1 - mu 0 0])*lStar - reshape(K.posKm(2, :, kS), 1, 3)).^2));
oPre = p.orient0 + L.sgn*(180/pi)*tPre;
xPre = L.gtoState(oPre, mod(p.gtoPhase0 + tPre/L.tauG, 1));
ePre = sqrt(sum(((xPre(1:3).' - [1 - mu 0 0])*lStar - reshape(K.posKm(1, :, kS), 1, 3)).^2));
xWrong = L.gtoState(p.orient0 - L.sgn*(180/pi)*tPre, mod(p.gtoPhase0 + tPre/L.tauG, 1));
eWrong = sqrt(sum(((xWrong(1:3).' - [1 - mu 0 0])*lStar - reshape(K.posKm(1, :, kS), 1, 3)).^2));
xPost = L.slotState(mod(p.droPhaseArr + (tPost - p.tArr)/L.tauD, 1));
ePost = sqrt(sum(((xPost(1:3).' - [1 - mu 0 0])*lStar - reshape(K.posKm(3, :, kS), 1, 3)).^2));
fprintf(['  positions sat %d: mid-flight %.3g km, on the GTO %.3g km (the GTO turned the wrong way: ', ...
         '%.0f km), on the DRO %.3g km\n'], kS, eMid, ePre, eWrong, ePost);
assert(eMid < 1 && ePre < 1e-6 && ePost < 1e-6, 'drawn positions differ from the saved data');
assert(eWrong > 100, 'deliberate failure: the wrong drift sign is not detectable (%.3g km)', eWrong);
assert(isequal(K.phase(:, kS).', [1 2 3]), 'phases of sat %d: %s', kS, mat2str(K.phase(:, kS).'));

% ---- 5. plume ----
[tg, ia] = unique(double(p.t(:)), 'stable');
al = double(p.alpha(:, ia)).';
nNode = round(numel(tg)/2);
Q = leg1_movie(leg1Mat, struct('dayList', [tg(nNode)*dayNd 0 I.T0Days], 'previewDir', tmp, 'pngPrefix', 't'));
eNode = -al(nNode, :)/sqrt(sum(al(nNode, :).^2));
errN = max(abs(Q.plumeDir(1, :, kS) - eNode));
Tv = viewmtx(Q.camera.az, Q.camera.el(Q.frameIdx(1)));
pe = (Tv(1:2, 1:3)*eNode.').';
sd = Q.plumeScreenDir(1, :, kS);
fprintf('  plume sat %d at node %d: direction err %.2e; drawn [%.1f %.1f] px vs exhaust [%.3f %.3f]\n', ...
        kS, nNode, errN, sd, pe);
assert(errN < 1e-12, 'plume direction differs from -alpha/|alpha|');
assert(dot(sd, pe) > 0, 'drawn plume points against the exhaust');
assert(Q.plumeDrawn(1, kS) && ~any(Q.plumeDrawn(3, :)), 'plume on/off wrong (mid-flight / T0)');
assert(isequal(Q.plumeDrawn(2, :), Q.phase(2, :) == 2), 'day 0: plume on a stage that is not thrusting');
assert(Q.stackMax(2) == 6, 'day 0: the six-stage stack is not labelled (%d)', Q.stackMax(2));
fprintf('  plumes at day 0: %s (thrusting %s); stack of %d labelled; none at T0\n', ...
        mat2str(Q.plumeDrawn(2, :)), mat2str(Q.phase(2, :) == 2), Q.stackMax(2));
assert(Q.earthTextured || ~isfile(fullfile(getenv('HOME'), 'Desktop', 'proj7', 'external', ...
       'pumpkyn', 'src', '+pumpkyn', '+util', 'earth-clouds-4k.jpg')), 'Earth not textured');
assert(Q.moonTextured, 'Moon not textured');
assert(all(Q.outlineDrawn, 'all'), 'an outline not drawn');

% ---- 6. last frame = mission_movie's first frame ----
M = mission_movie(leg2Mat, struct('dayList', 0, 'previewDir', tmp));
dT = max(abs(Q.camTarget(3, :) - M.camTarget(1, :)));
dP = max(abs(Q.camPos(3, :) - M.camPos(1, :)));
dV = abs(Q.camViewAngle(3) - M.camViewAngle(1));
fprintf('  last frame vs mission_movie first frame: target %.2e km, position %.2e km, view angle %.2e deg\n', ...
        dT, dP, dV);
assert(dT < 1e-9 && dP < 1e-6 && dV < 1e-9, 'last frame does not match mission_movie''s first');

% ---- 7. dayOffset ----
assert(isequal(M.displayDay, M.frameDay), 'default dayOffset changes the displayed day');
M3 = mission_movie(leg2Mat, struct('dayList', 0, 'previewDir', tmp, 'dayOffset', 3));
fprintf('  mission_movie displayed first day: default %.3f, dayOffset 3 -> %.3f\n', M.displayDay(1), ...
        M3.displayDay(1));
assert(abs(M3.displayDay(1) - M.displayDay(1) - 3) < 1e-12 && M3.dayOffset == 3, 'dayOffset not applied');

% ---- 8. join ----
mp42 = fullfile(tmp, 'leg2_short.mp4');
mission_movie(leg2Mat, struct('outFile', mp42, 'maxFrames', 30, 'dayOffset', I.T0Days));
joined = fullfile(tmp, 'joined.mp4');
Jn = full_mission_movie(leg1Mat, leg2Mat, struct('render', false, 'leg1Mp4', mp4, 'leg2Mp4', mp42, ...
                                                  'outFile', joined, 'fadeFrames', 10));
vr = VideoReader(joined);
fprintf('  join: %d + %d - %d = %d frames (read %d), %d x %d\n', n1, 30, 10, Jn.nFrames, ...
        vr.NumFrames, vr.Width, vr.Height);
assert(vr.NumFrames == n1 + 30 - 10 && Jn.nFrames == vr.NumFrames, 'joined frame count');
assert(vr.Width == 1280 && vr.Height == 720, 'joined size');
clear vr
small = fullfile(tmp, 'small.mp4');
vw = VideoWriter(small, 'MPEG-4');  vw.FrameRate = 30;  open(vw);
for kk = 1:5, writeVideo(vw, zeros(360, 640, 3, 'uint8')); end
close(vw);
thrown = '';
try
    full_mission_movie(leg1Mat, leg2Mat, struct('render', false, 'leg1Mp4', mp4, 'leg2Mp4', small, ...
                                                 'outFile', fullfile(tmp, 'bad.mp4')));
catch err
    thrown = err.identifier;
end
assert(strcmp(thrown, 'full_mission_movie:mismatch') && ~isfile(fullfile(tmp, 'bad.mp4')), ...
       'frame-size mismatch not refused (%s)', thrown);
fprintf('  join of a 640 x 360 movie refused (full_mission_movie:mismatch), nothing written\n');

% ---- 9. deliberate failures and refusals ----
ids = cell(1, 3);
calls = {{leg1Mat, struct('dayList', 0, 'previewDir', tmp, 'zoomScale', 0.5)}, ...
         {fullfile(tmp, 'nope.mat'), struct()}, ...
         {leg1Mat, struct('dayList', I.T0Days + 1, 'previewDir', tmp)}};
for kb = 1:3
    try
        leg1_movie(calls{kb}{:});
    catch err
        ids{kb} = err.identifier;
    end
end
fprintf('  refusals: %s\n', strjoin(ids, ', '));
assert(isequal(ids, {'leg1_movie:containment', 'leg1_movie:noResults', 'leg1_movie:dayList'}), ...
       'refusal ids: %s', strjoin(ids, ' '));
ok = true;
fprintf('  PASS\n');
end
