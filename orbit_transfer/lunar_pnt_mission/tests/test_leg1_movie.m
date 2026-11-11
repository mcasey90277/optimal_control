function [ok, skipped, note] = test_leg1_movie()
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
%   FIRST frame (target, position, view angle to 1e-9); (7) dayOffset, on
%   what is DISPLAYED: the day counter's String (info.dayText, read back
%   from the text object) of the same map frame reads the mission day with
%   the default and the mission day + 3 with dayOffset = 3; the default
%   dayFormat leaves the frame pixel-identical to the explicit 'Day %4.1f';
%   (7b) the counter does not step back at the join: leg 1's T0 still vs
%   leg 2's first still with dayOffset = T0 and the two-decimal format
%   (the default one-decimal format WOULD step back: asserted, so the
%   check can fail);
%   (8) JOIN: full_mission_movie on two short movies writes n1 + n2 - fade
%   frames at 1280 x 720, and refuses a second movie of another frame size
%   (full_mission_movie:mismatch); (8b) full_mission_movie rendering both
%   legs short (leg-1 preview, leg-2 20 frames): leg 2 in 'Day %5.2f', its
%   first displayed day >= leg 1's last, leg 1's end framing from the
%   joined leg2Mat; (8c) full_mission_movie refusals before any render:
%   a leg-1 MP4 older than its results with render = false (stale), leg-1
%   results not PASS (combined_timeline:verdict), leg1Opts.leg2Mat not the
%   joined one (leg2Mat), the owner's results/lunar_pnt_mission.mp4 as
%   outFile by a relative path and as leg2Mp4 by a '..' path (protected);
%   (8d) leg1_movie: the owner's movie by a relative and a '..' path
%   (leg1_movie:protected, the file untouched); results with verdict not
%   PASS, a G.pass false, a seam knob set -> leg1_movie:verdict; with
%   allowUnvetted the frame's regime text is stamped UNVETTED RESULTS (and
%   the vetted results' is not). The package is found by SHA-256
%   (leg1_package_for), never by file name. (9) DELIBERATE FAILURES: zoom forced to
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
%  note                     char                    '' or which package
%                                                   snapshot was used
%
%% Revision History:
%  M. Casey                                                   (c) 11/11/2026
%  M. Casey  final fix wave: displayed-day checks, join counter,
%            vetting, protected paths, full_mission_movie
%            refusals, package by SHA-256                     (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_lunar_pnt_paths.m'));
resDir  = fullfile(fileparts(here), 'results');
leg1Mat = fullfile(resDir, 'lunar_pnt_leg1_study.mat');
leg2Mat = fullfile(resDir, 'lunar_pnt_mission_study.mat');
skipped = false;
note = '';
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
Rm = load(leg1Mat, 'run_meta');
[pkg, whyPkg, pkgInfo] = leg1_package_for(Rm.run_meta.identity);      % by SHA-256, never by file name
assert(~isempty(pkg), 'no package matches the saved leg-1 results: %s', whyPkg);
assert(strcmp(I.packageFile, pkg) || strcmp(file_sha256(I.packageFile), Rm.run_meta.identity.sha256), ...
       'leg1_movie used a package that is not the results''');
fprintf('  package %s... found by hash: %s (%s)\n', Rm.run_meta.identity.sha256(1:8), pkg, pkgInfo.source);
if ~strcmp(pkgInfo.source, 'primary')
    note = sprintf('the default package is not the saved results''; used %s', pkg);
end
L  = leg1_inputs(struct('catalogMat', pkg));
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

% ---- 7. dayOffset and dayFormat, on what is DISPLAYED (the counter's String) ----
M3 = mission_movie(leg2Mat, struct('dayList', 0, 'previewDir', fullfile(tmp, 'off3'), 'dayOffset', 3));
fprintf('  mission_movie day counter at map frame %d (day %.3f): default "%s", dayOffset 3 "%s"\n', ...
        M.frameIdx(1), M.frameDay(1), M.dayText{1}, M3.dayText{1});
assert(strcmp(M.dayText{1}, sprintf('Day %4.1f', M.frameDay(1))), 'default counter shows "%s"', M.dayText{1});
assert(strcmp(M3.dayText{1}, sprintf('Day %4.1f', M3.frameDay(1) + 3)) && M3.frameIdx(1) == M.frameIdx(1), ...
       'dayOffset 3: counter shows "%s"', M3.dayText{1});
assert(abs(day_of(M3.dayText{1}) - day_of(M.dayText{1}) - 3) < 1e-9, 'dayOffset not in the DISPLAYED day');
% the default format leaves frames pixel-identical: default vs the default
% format passed explicitly (and, measured once at the fix, vs the pre-fix code)
Mx = mission_movie(leg2Mat, struct('dayList', 0, 'previewDir', fullfile(tmp, 'fmt'), 'dayFormat', 'Day %4.1f'));
assert(isequal(imread(M.previewPngs{1}), imread(Mx.previewPngs{1})), 'default dayFormat changes the frame');
fprintf('  default dayFormat: frame pixel-identical to the explicit ''Day %%4.1f''\n');

% ---- 7b. the counter does not step back at the join ----
dLeg1 = day_of(Q.dayText{3});                          % leg 1's last frame (T0)
M2f = mission_movie(leg2Mat, struct('dayList', 0, 'previewDir', fullfile(tmp, 'join2'), 'dayOffset', I.T0Days, ...
                                    'dayFormat', 'Day %5.2f'));
M1f = mission_movie(leg2Mat, struct('dayList', 0, 'previewDir', fullfile(tmp, 'join1'), 'dayOffset', I.T0Days));
fprintf(['  join stills: leg 1 last "%s"; leg 2 first with the two-decimal format "%s" (with the default ' ...
         'format it would read "%s")\n'], Q.dayText{3}, M2f.dayText{1}, M1f.dayText{1});
assert(day_of(M2f.dayText{1}) >= dLeg1, 'the day counter steps back at the join: %s -> %s', ...
       Q.dayText{3}, M2f.dayText{1});
assert(day_of(M1f.dayText{1}) < dLeg1, ['deliberate failure: the default one-decimal format no longer ' ...
       'steps back (%s -> %s); this check would not catch a regression'], Q.dayText{3}, M1f.dayText{1});

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

% ---- 8b. full_mission_movie renders both legs (short): end framing, counter ----
jr = fullfile(tmp, 'joined_render.mp4');
Jr = full_mission_movie(leg1Mat, leg2Mat, struct('outFile', jr, 'leg1Mp4', fullfile(tmp, 'r1.mp4'), ...
                        'leg2Mp4', fullfile(tmp, 'r2.mp4'), 'fadeFrames', 5, ...
                        'leg1Opts', struct('preview', true), 'leg2Opts', struct('maxFrames', 20)));
dEnd1 = day_of(Jr.leg1.dayText{end});  dBeg2 = day_of(Jr.leg2.dayText{1});
fprintf('  rendered join: leg 1 ends "%s", leg 2 starts "%s" (format %s); %d frames\n', Jr.leg1.dayText{end}, ...
        Jr.leg2.dayText{1}, Jr.leg2.dayFormat, Jr.nFrames);
assert(dBeg2 >= dEnd1 && strcmp(Jr.leg2.dayFormat, 'Day %5.2f'), 'rendered join: counter steps back');
assert(max(abs(Jr.leg1.camera.targetEnd - M.camTarget(1, :))) < 1e-9, 'rendered join: end framing not from leg2Mat');
assert(Jr.nFrames == 60 + 20 - 5 && ~isempty(Jr.checks), 'rendered join: frames %d', Jr.nFrames);

% ---- 8c. full_mission_movie refusals (cheap: all before any rendering) ----
protMp4 = fullfile(resDir, 'lunar_pnt_mission.mp4');
protInfo = dir(protMp4);
oldMp4 = fullfile(tmp, 'old_leg1.mp4');
copyfile(mp4, oldMp4);
java.io.File(oldMp4).setLastModified(int64(946684800000));             % 2000-01-01: older than the results
T = load(leg1Mat);
T.run_meta.verdict = 'not passed: GH';
badLeg1 = fullfile(tmp, 'leg1_not_passed.mat');
save(badLeg1, '-struct', 'T', '-v7.3');
% relative paths are resolved against MATLAB's pwd; the JVM's user.dir
% does not follow cd, so work from a folder that is NOT user.dir (the old
% java-only resolution missed the protected file from there)
here0 = pwd;
cleanCd = onCleanup(@() cd(here0));
jvmDir = char(java.lang.System.getProperty('user.dir'));
if ~strcmp(canon_path(jvmDir), canon_path(resDir))
    cd(resDir);
    relProt = 'lunar_pnt_mission.mp4';
    dotProt = fullfile('inputs', '..', 'lunar_pnt_mission.mp4');
else
    cd(fullfile(resDir, 'inputs'));
    relProt = fullfile('..', 'lunar_pnt_mission.mp4');
    dotProt = fullfile('..', 'inputs', '..', 'lunar_pnt_mission.mp4');
end
fprintf('  protected-path checks from pwd %s (JVM user.dir %s): "%s", "%s"\n', pwd, jvmDir, relProt, dotProt);
fm = {{leg1Mat, leg2Mat, struct('render', false, 'leg1Mp4', oldMp4, 'leg2Mp4', mp42, ...
                                'outFile', fullfile(tmp, 'x1.mp4'))}, 'full_mission_movie:stale'; ...
      {badLeg1, leg2Mat, struct('render', false, 'leg1Mp4', mp4, 'leg2Mp4', mp42, ...
                                'outFile', fullfile(tmp, 'x2.mp4'))}, 'combined_timeline:verdict'; ...
      {leg1Mat, leg2Mat, struct('leg1Opts', struct('leg2Mat', badLeg1), ...
                                'outFile', fullfile(tmp, 'x3.mp4'))}, 'full_mission_movie:leg2Mat'; ...
      {leg1Mat, leg2Mat, struct('render', false, 'leg1Mp4', mp4, 'leg2Mp4', mp42, ...
                                'outFile', relProt)}, 'full_mission_movie:protected'; ...
      {leg1Mat, leg2Mat, struct('render', false, 'leg1Mp4', mp4, ...
                                'leg2Mp4', dotProt, ...
                                'outFile', fullfile(tmp, 'x5.mp4'))}, 'full_mission_movie:protected'};
for kf = 1:size(fm, 1)
    id = '';
    try
        full_mission_movie(fm{kf, 1}{:});
    catch err
        id = err.identifier;
        msg = err.message;
    end
    assert(strcmp(id, fm{kf, 2}), 'full_mission_movie refusal %d: expected %s, got "%s"', kf, fm{kf, 2}, id);
    fprintf('  full_mission_movie refused (%s): %s\n', id, msg(1:min(150, end)));
end
assert(~isfile(fullfile(tmp, 'x1.mp4')) && ~isfile(fullfile(tmp, 'x2.mp4')) && ~isfile(fullfile(tmp, 'x5.mp4')), ...
       'a refused join wrote its output');

% ---- 8d. leg1_movie: the protected file (relative and '..' paths) and vetting ----
lm = {struct('dayList', 0, 'previewDir', tmp, 'outFile', relProt), ...
      struct('preview', true, 'outFile', dotProt)};
for kl = 1:numel(lm)
    id = '';
    try
        leg1_movie(leg1Mat, lm{kl});
    catch err
        id = err.identifier;
    end
    assert(strcmp(id, 'leg1_movie:protected'), 'leg1_movie protected path %d: got "%s"', kl, id);
end
fprintf('  leg1_movie refused the owner''s movie by a relative and a ''..'' path (leg1_movie:protected)\n');
if ~isempty(protInfo)
    nowP = dir(protMp4);
    assert(nowP.datenum == protInfo.datenum && nowP.bytes == protInfo.bytes, 'the owner''s movie was touched');
end
clear cleanCd
T2 = load(leg1Mat);  T2.P1(3).G.pass = false;
badG = fullfile(tmp, 'leg1_gpass.mat');  save(badG, '-struct', 'T2', '-v7.3');
T3 = load(leg1Mat);  T3.run_meta.knobs.seamFinalT0 = 0.01;
badS = fullfile(tmp, 'leg1_seam.mat');  save(badS, '-struct', 'T3', '-v7.3');
bads = {badLeg1, badG, badS};
for kv = 1:numel(bads)
    id = '';
    try
        leg1_movie(bads{kv}, struct('dayList', 0, 'previewDir', tmp));
    catch err
        id = err.identifier;
        msg = err.message;
    end
    assert(strcmp(id, 'leg1_movie:verdict'), 'unvetted results %d rendered (got "%s")', kv, id);
    fprintf('  leg1_movie refused: %s\n', msg);
end
U = leg1_movie(badLeg1, struct('dayList', 0, 'previewDir', tmp, 'pngPrefix', 'unvetted', 'allowUnvetted', true));
assert(U.unvetted && strcmp(U.regimeText{1}, 'UNVETTED RESULTS') && ~Q.unvetted && ...
       ~any(contains(Q.regimeText, 'UNVETTED')), 'allowUnvetted: stamp missing (or on vetted results)');
fprintf('  allowUnvetted: rendered with "%s" stamped (%s); vetted results carry no stamp\n', U.regimeText{1}, U.vetting);

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

% ------------------------------------------------------------------------
function d = day_of(txt)
% DAY_OF  The number a day counter displays ('Day  2.65' -> 2.65).
%   INPUTS: txt char.  OUTPUTS: d scalar.
d = sscanf(regexprep(txt, '^\D*', ''), '%f', 1);
assert(isscalar(d) && isfinite(d), 'no day number in "%s"', txt);
end
