function info = full_mission_movie(leg1Mat, leg2Mat, opts)
%% Purpose:
%
%   The FULL lunar PNT mission movie: leg 1 (six 5 N stages, one GTO -> the
%   DRO slots, leg1_movie) followed by leg 2 (DRO -> tulip constellation with
%   the south-pole GDOP strip, mission_movie rendered with dayOffset = T0 in
%   days so its day counter continues leg 1's clock), joined into ONE MP4.
%
%   Join: one VideoWriter reads both rendered files frame by frame
%   (VideoReader; no ffmpeg). Frame size and frame rate must match
%   (full_mission_movie:mismatch, nothing written). The scene of leg 1's
%   last frame IS leg 2's first (same camera; leg1_movie's contract), but
%   the overlays differ (strip, captions, legend, tulip outline), so the
%   last fadeFrames of leg 1 are cross-faded linearly into the first
%   fadeFrames of leg 2 (both are inside their movies' holds, so nothing
%   moves during the fade): frames = n1 + n2 - fadeFrames. Written to
%   <name>.partial.mp4 and moved onto outFile after VideoReader confirms
%   the frame count and size.
%
%   The leg-2 render goes to its OWN file (default
%   results/lunar_pnt_mission_dayoffset.mp4) with the day counter in leg
%   1's two-decimal format (mission_movie dayFormat 'Day %5.2f'), so the
%   displayed day never steps back at the join. PROTECTED: no output path
%   (outFile, leg1Mp4, leg2Mp4, or their partial files) may resolve to
%   results/lunar_pnt_mission.mp4, relative paths taken against MATLAB's
%   pwd and '..' resolved (full_mission_movie:protected).
%
%   BEFORE anything is rendered or joined: combined_timeline's
%   compatibility checks on the two results (package identity, DRO,
%   phases, units, both verdicts, epoch; no figure), so unvetted or
%   mismatched results throw combined_timeline:<check>. The leg-1 render
%   gets leg2Mat, so its end framing comes from the SAME leg-2 results
%   that are joined (a different leg1Opts.leg2Mat: full_mission_movie:
%   leg2Mat). render = false joins only MP4s that exist and are newer than
%   their results (the leg-2 one also newer than the leg-1 results), else
%   full_mission_movie:stale; frame size and rate must match (mismatch).
%
%% Inputs:
%
%  leg1Mat                  char                    leg-1 study results ('' =
%                                                   results/lunar_pnt_leg1_study.mat)
%  leg2Mat                  char                    leg-2 study results ('' =
%                                                   results/lunar_pnt_mission_study.mat)
%  opts                     struct (optional)
%    .outFile               char                    joined MP4
%                                                   [results/lunar_pnt_full_mission.mp4]
%    .render                logical [true]          Render both legs first
%                                                   (false: join the files
%                                                   named below as they are)
%    .leg1Mp4               char                    leg-1 MP4
%                                                   [results/lunar_pnt_leg1.mp4]
%    .leg2Mp4               char                    leg-2 MP4 (day offset)
%                                                   [results/lunar_pnt_mission_dayoffset.mp4]
%    .leg1Opts, .leg2Opts   struct [struct()]       Extra options for the
%                                                   renders
%    .fadeFrames            scalar [15]             Cross-fade length [frames]
%    .quality               scalar [95]             VideoWriter quality
%
%% Outputs:
%
%  info.outFile             char                    Joined MP4
%  info.nFrames             scalar                  Frames written
%  info.nLeg1, .nLeg2       scalars                 Frames read from each leg
%  info.fps                 scalar                  Frame rate [1/s]
%  info.frameSize           [1 x 2]                 Width, height [px]
%  info.dayOffset           scalar                  T0 [d] given to leg 2
%  info.joinFrame           scalar                  First frame of the fade
%  info.leg1, .leg2         struct                  Render infos (empty when
%                                                   render = false)
%  info.checks              {1 x m} cell            combined_timeline checks
%  info.wallS               scalar                  Wall time [s]
%
%% Revision History:
%  M. Casey                                                   (c) 11/11/2026
%  M. Casey  final fix wave: compatibility first, leg2Mat passed
%            through, stale-file and protected refusals, two-
%            decimal day counter on leg 2                     (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(fileparts(mfilename('fullpath')));        % lunar_pnt_mission/
if nargin == 0
    run(fullfile(here, 'setup_lunar_pnt_paths.m'));
    info = full_mission_movie('', '', struct());
    fprintf('full_mission_movie self-demo: %s, %d frames (%.0f s)\n', info.outFile, ...
            info.nFrames, info.wallS);
    return
end
if nargin < 3 || isempty(opts), opts = struct(); end
if isempty(leg1Mat), leg1Mat = fullfile(here, 'results', 'lunar_pnt_leg1_study.mat'); end
if isempty(leg2Mat), leg2Mat = fullfile(here, 'results', 'lunar_pnt_mission_study.mat'); end
tWall = tic;

def = struct('outFile',    fullfile(here, 'results', 'lunar_pnt_full_mission.mp4'), ...
             'render',     true, ...
             'leg1Mp4',    fullfile(here, 'results', 'lunar_pnt_leg1.mp4'), ...
             'leg2Mp4',    fullfile(here, 'results', 'lunar_pnt_mission_dayoffset.mp4'), ...
             'leg1Opts',   struct(), ...
             'leg2Opts',   struct(), ...
             'fadeFrames', 15, ...
             'quality',    95);
for f = fieldnames(opts).'
    if ~isfield(def, f{1})
        error('full_mission_movie:opts', 'unknown option ''%s''', f{1});
    end
    def.(f{1}) = opts.(f{1});
end
o = def;
% the owner's leg-2 movie is never an output (paths resolved against MATLAB's pwd)
protectedMp4 = canon_path(fullfile(here, 'results', 'lunar_pnt_mission.mp4'));
outs = {o.outFile, o.leg1Mp4, o.leg2Mp4};
outNames = {'outFile', 'leg1Mp4', 'leg2Mp4'};
for kk = 1:numel(outs)
    [pd, pn] = fileparts(canon_path(outs{kk}));
    if any(strcmp(protectedMp4, {canon_path(outs{kk}), fullfile(pd, [pn '.partial.mp4'])}))
        error('full_mission_movie:protected', 'refusing to write %s (the owner''s movie): %s resolves to it', ...
              protectedMp4, outNames{kk});
    end
end
if isfield(o.leg1Opts, 'leg2Mat') && ~strcmp(canon_path(o.leg1Opts.leg2Mat), canon_path(leg2Mat))
    error('full_mission_movie:leg2Mat', ['leg1Opts.leg2Mat (%s) is not the leg-2 results being joined (%s): ' ...
          'the end framing must come from the same results'], o.leg1Opts.leg2Mat, leg2Mat);
end

% compatibility of the two results BEFORE anything is rendered or joined
% (combined_timeline's checks: package identity, DRO, phases, units,
% verdicts, epoch; no figure)
[~, tl] = combined_timeline(leg1Mat, leg2Mat, struct());

% day offset: T0 in days, from the leg-1 results
D = load(leg1Mat, 'T0', 'run_meta');
dayOffset = double(D.T0)*D.run_meta.regime.tStar/86400;

I1 = [];  I2 = [];
if o.render
    o1 = o.leg1Opts;  o1.outFile = o.leg1Mp4;  o1.leg2Mat = leg2Mat;
    I1 = leg1_movie(leg1Mat, o1);
    o2 = o.leg2Opts;  o2.outFile = o.leg2Mp4;  o2.dayOffset = dayOffset;
    if ~isfield(o2, 'dayFormat'), o2.dayFormat = 'Day %5.2f'; end      % leg 1's: no step back at the join
    I2 = mission_movie(leg2Mat, o2);
else
    % joining existing files: each must exist and be newer than the
    % results it was rendered from (the leg-2 file also than the leg-1
    % results, which set its day offset)
    t1 = file_time(leg1Mat);  t2 = file_time(leg2Mat);
    m1 = file_time(o.leg1Mp4);  m2 = file_time(o.leg2Mp4);
    if ~(m1 > t1 && m2 > max(t1, t2))
        error('full_mission_movie:stale', ['render = false joins only movies newer than their results: ' ...
              'leg-1 MP4 %s (%s) vs leg-1 results (%s); leg-2 MP4 %s (%s) vs leg-2 results (%s) and ' ...
              'leg-1 results; render them again (render = true)'], o.leg1Mp4, stamp(m1), stamp(t1), ...
              o.leg2Mp4, stamp(m2), stamp(t2));
    end
end

% ---- join ----
r1 = VideoReader(o.leg1Mp4);
r2 = VideoReader(o.leg2Mp4);
n1 = r1.NumFrames;  n2 = r2.NumFrames;
% VideoReader reports the rate of a short MP4 in single precision
% (30 fps reads 29.9999981 for 5 or 20 frames): compare to 1e-6 relative
if r1.Width ~= r2.Width || r1.Height ~= r2.Height || abs(r1.FrameRate - r2.FrameRate) > 1e-6*r1.FrameRate
    error('full_mission_movie:mismatch', 'leg 1 %d x %d at %g fps, leg 2 %d x %d at %g fps', ...
          r1.Width, r1.Height, r1.FrameRate, r2.Width, r2.Height, r2.FrameRate);
end
K = o.fadeFrames;
if K < 0 || K ~= round(K) || K > min(n1, n2)
    error('full_mission_movie:fade', 'fadeFrames must be an integer in 0..%d', min(n1, n2));
end
[pDir, pName] = fileparts(o.outFile);
if ~isfolder(pDir), mkdir(pDir); end
partFile = fullfile(pDir, [pName '.partial.mp4']);
vw = VideoWriter(partFile, 'MPEG-4');
vw.FrameRate = r1.FrameRate;
vw.Quality   = o.quality;
open(vw);
cleanVw = onCleanup(@() close_vw(vw, partFile));
for kk = 1:n1 - K
    writeVideo(vw, readFrame(r1));
end
for kk = 1:K
    a = double(readFrame(r1));
    b = double(readFrame(r2));
    w = kk/(K + 1);
    writeVideo(vw, uint8(round((1 - w)*a + w*b)));
end
for kk = K + 1:n2
    writeVideo(vw, readFrame(r2));
end
close(vw);
fps = r1.FrameRate;  sz = [r1.Width r1.Height];
clear r1 r2
vr = VideoReader(partFile);
nOut = vr.NumFrames;
okV = nOut == n1 + n2 - K && vr.Width == sz(1) && vr.Height == sz(2);
clear vr
if ~okV
    error('full_mission_movie:verify', 'joined video failed its check: %s (%d frames)', partFile, nOut);
end
movefile(partFile, o.outFile, 'f');
clear cleanVw
info = struct('outFile', o.outFile, 'nFrames', nOut, 'nLeg1', n1, 'nLeg2', n2, 'fps', fps, ...
              'frameSize', sz, 'dayOffset', dayOffset, 'joinFrame', n1 - K + 1, ...
              'leg1', I1, 'leg2', I2, 'checks', {tl.checks}, 'wallS', toc(tWall));
end

%% ------------------------------------------------------------------------
function t = file_time(f)
% FILE_TIME  Modification time of a file as a datenum, -Inf when missing.
%   INPUTS: f char.  OUTPUTS: t scalar.
d = dir(f);
if isscalar(d) && ~d.isdir, t = d.datenum; else, t = -Inf; end
end

function s = stamp(t)
% STAMP  A datenum as text, 'missing' for -Inf.  INPUTS: t scalar.
%   OUTPUTS: s char.
if isinf(t), s = 'missing'; else, s = char(datetime(t, 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm:ss')); end
end

function close_vw(vw, partFile)
% CLOSE_VW  Close a VideoWriter if open and delete an unfinished partial file.
%   INPUTS: vw VideoWriter; partFile char.
if isvalid(vw), close(vw); end
if isfile(partFile), delete(partFile); end
end
