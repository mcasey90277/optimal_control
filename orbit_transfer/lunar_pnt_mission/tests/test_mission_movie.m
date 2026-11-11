function [ok, skipped] = test_mission_movie()
%% Purpose:
%
%   Verifies mission_movie on the SAVED study results (nothing re-solved):
%   (1) MP4 mode with maxFrames = 12 and an explicit outFile writes a
%   1280 x 720 MPEG-4 of exactly 12 frames (VideoReader) and leaves no
%   .partial.mp4 behind; a maxFrames render without an explicit outFile, or
%   aimed at the default movie, is refused (mission_movie:partialOutFile)
%   and the default movie is untouched; (2) CAMERA: every camera mode is
%   manual; the first and last frames sit exactly at the end framing; the
%   widest frame is wider than it; the view direction and target never
%   change; all six spacecraft lie inside the scene area by >= 4 px in every
%   sampled frame (and over the whole map, camera.minMarginPx); the DRO and
%   tulip outlines lie inside the scene area at the end framing and are
%   visible in every sampled frame; (3) the km -> pixel mapping the
%   containment check, clip box and scale bar rely on matches where MATLAB
%   actually draws: two probe markers, found by colour in the rendered PNG,
%   land within 2 px of camera.toPx at the end framing and at the widest
%   zoom; (4) PLUME: visible mid-transfer, none after the last arrival; its
%   direction equals +lam_v/|lam_v| of the saved flight to 1e-12 at a flight
%   node AND at a time between nodes (an independent interp1 of the saved
%   flight); the DRAWN plume (tip minus nozzle on screen) points along the
%   projected exhaust and against the projected thrust; (5) the Moon carries
%   the texture when the texture file exists; (6) refusals: a missing
%   results file and dayList entries outside [0, end] throw their ids; a
%   spacecraft with NaN tDep (on station from t = 0, no P entry) is drawn on
%   station with no plume. Returns skipped = true (and
%   ok = true) with a message when the results file does not exist, so the
%   runner can show SKIP.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 true if every assertion
%                                                   held (else it throws)
%  skipped                  logical                 true if the results file
%                                                   is missing (nothing run)
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_lunar_pnt_paths.m'));
resultsMat = fullfile(fileparts(here), 'results', 'lunar_pnt_mission_study.mat');
skipped = false;
if ~isfile(resultsMat)
    fprintf('  SKIPPED: no results file (%s); run lunar_pnt_mission_study first\n', resultsMat);
    ok = true;  skipped = true;
    return
end
tmp = tempname;
mkdir(tmp);
cleanTmp = onCleanup(@() rmdir(tmp, 's'));
dayNd = 382981.289129055/86400;
defMp4 = fullfile(fileparts(here), 'results', 'lunar_pnt_mission.mp4');
defInfo = dir(defMp4);

% ---- 1. MP4, 12 frames, explicit outFile; partial renders refused ----
mp4  = fullfile(tmp, 'movie_test.mp4');
info = mission_movie(resultsMat, struct('outFile', mp4, 'maxFrames', 12));
assert(isfile(mp4), 'MP4 not written');
assert(~isfile(fullfile(tmp, 'movie_test.partial.mp4')), 'partial file left behind');
vr = VideoReader(mp4);
fprintf('  MP4: %d x %d, %d frames, %.3f s, wall %.1f s\n', vr.Width, vr.Height, ...
        vr.NumFrames, vr.Duration, info.wallS);
assert(vr.Width == 1280 && vr.Height == 720, 'MP4 size %d x %d', vr.Width, vr.Height);
assert(vr.NumFrames == 12 && info.nFrames == 12, 'MP4 frames %d (info %d)', vr.NumFrames, info.nFrames);
clear vr
for tryOpts = {struct('maxFrames', 3), struct('maxFrames', 3, 'outFile', defMp4)}
    thrown = '';
    try
        mission_movie(resultsMat, tryOpts{1});
    catch err
        thrown = err.identifier;
    end
    assert(strcmp(thrown, 'mission_movie:partialOutFile'), 'partial render not refused (%s)', thrown);
end
if ~isempty(defInfo)
    nowInfo = dir(defMp4);
    assert(nowInfo.datenum == defInfo.datenum && nowInfo.bytes == defInfo.bytes, 'default movie touched');
end
fprintf('  partial renders refused (mission_movie:partialOutFile); default movie untouched\n');
C = info.camera;
assert(all(info.camModesManual), 'a camera mode is not manual');
assert(isequal(info.camPos(1, :), C.position) && isequal(info.camTarget(1, :), C.target) && ...
       abs(info.camViewAngle(1) - C.viewAngleEnd) < 1e-12, 'first frame is not the end framing');
assert(all(info.outlineVisible, 'all') && C.droInside && C.tulInside, 'outline check (MP4 frames)');

% ---- 2-5. previews: frame 1, a node, between nodes, widest, after arrival, end ----
L  = load(resultsMat, 'H', 'P');
P  = L.P;  if iscell(P), P = [P{:}]; end
kS = 5;
[tf, ia] = unique(double(P(kS).flight.t(:)), 'stable');
Y    = double(P(kS).flight.Y(ia, :));
nNode = round(numel(tf)/2);
dNode = (double(P(kS).tDep) + tf(nNode))*dayNd;
tauOff = 0.37*tf(nNode) + 0.63*tf(nNode + 1);
dOff  = (double(P(kS).tDep) + tauOff)*dayNd;
dPost = max(double(L.H.tArr))*dayNd + 5;
dEnd  = double(L.H.t(end))*dayNd;
days  = [dNode dOff C.zoomMaxDay dPost dEnd];
I = mission_movie(resultsMat, struct('frameList', 1, 'dayList', days, 'previewDir', tmp));
nP = numel(I.previewPngs);
assert(nP == 6, 'expected 6 PNGs, got %d', nP);
for kk = 1:nP
    sz = size(imread(I.previewPngs{kk}));
    assert(isequal(sz, [720 1280 3]), 'PNG %d size %s', kk, mat2str(sz));
end
% camera
vd = (I.camPos - I.camTarget)./sqrt(sum((I.camPos - I.camTarget).^2, 2));
assert(all(abs(vd - vd(1, :)) < 1e-12, 'all') && all(abs(I.camTarget - C.target) < 1e-9, 'all'), ...
       'view direction or target changed');
assert(all(I.camModesManual), 'a camera mode is not manual');
assert(abs(I.camViewAngle(1) - C.viewAngleEnd) < 1e-12 && abs(I.camViewAngle(end) - C.viewAngleEnd) < 1e-12, ...
       'first/last frame not at the end framing');
assert(I.camViewAngle(4) > C.viewAngleEnd*1.2, 'widest frame not wider than the end framing');
assert(all(I.scMarginPx >= 4) && C.minMarginPx >= 4, 'a spacecraft within 4 px of the scene edge');
assert(all(I.outlineVisible, 'all'), 'an outline missing in a preview frame');
fprintf(['  camera: modes manual, direction/target fixed, end framing at first/last, widest x%.3f ', ...
         '(day %.2f); min margin %.1f px (map) / %.1f px (sampled)\n'], C.zoomMax, C.zoomMaxDay, ...
        C.minMarginPx, min(I.scMarginPx));
% plume on / off
assert(I.plumeCount(2) >= 1, 'no plume mid-transfer');
assert(I.plumeCount(5) == 0 && ~any(I.plumeDrawn(5, :)) && all(isnan(I.plumeDir(5, :, :)), 'all'), ...
       'plume after the last arrival');
% plume direction: node and between nodes (independent interp1 of the saved flight)
eNode = Y(nNode, 11:13)/sqrt(sum(Y(nNode, 11:13).^2));
lamOff = interp1(tf, Y(:, 11:13), tauOff);
eOff  = lamOff/sqrt(sum(lamOff.^2));
errN  = max(abs(I.plumeDir(2, :, kS) - eNode));
errO  = max(abs(I.plumeDir(3, :, kS) - eOff));
fprintf('  plume direction sc %d: node %d err %.2e, between nodes err %.2e\n', kS, nNode, errN, errO);
assert(errN < 1e-12 && errO < 1e-12, 'plume direction differs from +lam_v/|lam_v|');
% the DRAWN plume points along the projected exhaust, against the projected thrust
Tv = viewmtx(C.az, C.el);
pe = (Tv(1:2, 1:3)*eNode.').';
sd = I.plumeScreenDir(2, :, kS);
fprintf('  drawn plume on screen [%.1f %.1f] px; projected exhaust [%.3f %.3f]; dot %.2f / thrust %.2f\n', ...
        sd, pe, dot(sd, pe), dot(sd, -pe));
assert(dot(sd, pe) > 0 && dot(sd, -pe) < 0, 'drawn plume points the wrong way');
% Moon
if isfile(I.moonTexture)
    assert(I.moonTextured, 'texture file present but the Moon is not textured');
    fprintf('  Moon textured from %s\n', I.moonTexture);
else
    fprintf('  Moon texture file absent (%s): grey fallback, texture check skipped\n', I.moonTexture);
end

% ---- km -> pixel mapping vs where MATLAB draws, at two zoom levels ----
probe = [9000 -6000 5000; -14000 4000 -16000];
J = mission_movie(resultsMat, struct('dayList', [0 C.zoomMaxDay], 'previewDir', tmp, 'probeKm', probe));
for kk = 1:2
    img = double(imread(J.previewPngs{kk}));
    mag = img(:, :, 1) > 200 & img(:, :, 2) < 70 & img(:, :, 3) > 200;
    st  = regionprops(bwconncomp(mag), 'Centroid', 'Area');   % one region per blob
    st  = st([st.Area] >= 6);
    assert(numel(st) == 2, 'frame %d: expected 2 probe blobs, found %d', kk, numel(st));
    found = cat(1, st.Centroid);
    found = [found(:, 1) - 0.5, 720 - (found(:, 2) - 0.5)];      % figure px, origin bottom-left
    want  = reshape(J.probePx(kk, :, :), 2, []).';
    dErr  = zeros(2, 1);
    for qq = 1:2
        dErr(qq) = min(sqrt(sum((found - want(qq, :)).^2, 2)));
    end
    fprintf('  toPx vs drawn (Rv %.0f km): probe errors %.2f %.2f px\n', J.Rv(kk), dErr);
    assert(all(dErr < 2), 'toPx misses the drawn probe by %.2f px', max(dErr));
end
assert(J.Rv(2) > J.Rv(1)*1.2, 'probe frames not at two different zoom levels');
% ---- refusals and the on-station-from-t0 path ----
ids = cell(1, 3);
try, mission_movie(fullfile(tmp, 'nope.mat')); catch err, ids{1} = err.identifier; end
try, mission_movie(resultsMat, struct('dayList', -1, 'previewDir', tmp)); catch err, ids{2} = err.identifier; end
try, mission_movie(resultsMat, struct('dayList', dEnd + 1, 'previewDir', tmp)); catch err, ids{3} = err.identifier; end
assert(isequal(ids, {'mission_movie:noResults', 'mission_movie:dayList', 'mission_movie:dayList'}), ...
       'refusal ids: %s', strjoin(ids, ' '));
M = load(resultsMat, 'H', 'P', 'pole');
H = M.H;  Pc = num2cell(M.P);
H.tDep(6) = NaN;  H.tArr(6) = NaN;  H.phase(:, 6) = 3;  Pc{6} = [];
Sv = struct('H', H, 'P', {Pc}, 'pole', M.pole);
save(fullfile(tmp, 'onstation.mat'), '-struct', 'Sv');
K = mission_movie(fullfile(tmp, 'onstation.mat'), struct('dayList', [0 5], 'previewDir', tmp));
assert(all(K.phase(:, 6) == 3) && all(isnan(K.plumeDir(:, :, 6)), 'all'), ...
       'spacecraft with NaN tDep not drawn on station');
fprintf('  refusals: %s; NaN-tDep spacecraft drawn on station, no plume\n', strjoin(unique(ids), ', '));
ok = true;
fprintf('  PASS\n');
end
