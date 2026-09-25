function out = make_conjugate_video(opts)
%% Purpose:
%
%   The DIRECTOR: build the conjugate-point video from the explorer, the
%   shot list (video_scenes.m) and one narration audio file per scene.
%
%     1. For each scene, read its audio (audio/sceneNN.mp3 if present -- the
%        final voice -- else audio/sceneNN.wav, the draft from
%        make_draft_voice.sh) and measure its length.
%     2. Place every beat at the time the narrator reaches its anchor
%        phrase: the phrase's character position in the text, as a
%        fraction of the text, times the audio length.
%     3. Drive the explorer (headless, 1920 x 1080) through the beats and
%        render each shot: the full window, a close-up of one panel, or a
%        title card, with a caption bar; cross-fade between shots.
%     4. Write the silent video (out/conjugate_points_silent.mp4) and the
%        narration track laid out on the same clock (out/narration.wav),
%        then join them with ffmpeg if it is installed
%        (out/conjugate_points_video.mp4).
%
%   Usage:
%     make_conjugate_video                         % full quality, 30 fps
%     make_conjugate_video(struct('preview', true)) % 12 fps, quick look
%     make_conjugate_video(struct('scenes', 3))     % one scene only
%
%% Inputs:
%
%  opts                     struct (optional)       .preview [false]
%                                                   .fps [30; 12 in preview]
%                                                   .scenes [all] indices
%                                                   .lead [0.4 s] silence
%                                                   before each scene
%                                                   .tail [0.6 s] after
%                                                   .fade [0.35 s]
%                                                   .endCard [3 s]
%
%% Outputs:
%
%  out                      struct                  .video .audio .final
%                                                   (paths; final = '' if
%                                                   ffmpeg is missing)
%                                                   .duration [s] .beats
%                                                   (table of beat times)
%
%% Revision History:
%  M. Casey                                                   (c) 09/24/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, opts = struct(); end
d = @(f, v) getfielddef(opts, f, v);
preview = d('preview', false);
fps     = d('fps', 30 - 18*preview);
lead    = d('lead', 0.4);
tail    = d('tail', 0.6);
fadeS   = d('fade', 0.35);
endCard = d('endCard', 3);
W = 1920;  H = 1080;  capH = 120;

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));                              % the explorer and cov_* files
outDir = fullfile(here, 'out');
if ~isfolder(outDir), mkdir(outDir); end
tmpPng = [tempname '.png'];
S = video_scenes();
scenes = d('scenes', 1:numel(S));

app = conjugate_point_explorer('Visible', 'off', 'Preset', 1, 'Size', [W H]);
closer = onCleanup(@() delete(app.fig));
pinLight(app.fig);          % the Mac's evening dark mode must not reach the video

tag = '';  if preview, tag = '_preview'; end
vidFile = fullfile(outDir, ['conjugate_points_silent' tag '.mp4']);
vw = VideoWriter(vidFile, 'MPEG-4');
vw.FrameRate = fps;  vw.Quality = 92;
open(vw);
fsA = 44100;
track = zeros(0, 1);
tNow = 0;                                              % the video clock [s]
nWritten = 0;
prev = [];                                             % last frame written
beatLog = cell(0, 4);
srcExt = {};                                           % audio type per scene
view = 'full';  caption = '';  card = {};

for si = scenes
    sc = S(si);
    txt = strtrim(fileread(fullfile(here, 'narration', sc.file)));
    [y, dur, src] = sceneAudio(fullfile(here, 'audio', erase(sc.file, '.txt')), fsA);
    srcExt{end+1} = src(find(src == '.', 1, 'last'):end); %#ok<AGROW>
    fprintf('scene %d  %-12s  %5.1f s  (%s)\n', si, sc.file, dur, src);
    track = [track; zeros(round(lead*fsA), 1); y; zeros(round(tail*fsA), 1)]; %#ok<AGROW>
    T = lead + dur + tail;

    % beat times from anchor positions in the text
    nb = numel(sc.beats);
    tb = zeros(1, nb);
    for k = 1:nb
        tb(k) = lead*(k > 1) + dur*anchorFrac(txt, sc.beats(k).anchor);
    end
    assert(all(diff(tb) > 0), 'make_conjugate_video:order', ...
           'scene %d: beats are not in narration order', si);

    for k = 1:nb
        tEnd = T;  if k < nb, tEnd = tb(k+1); end
        ops = sc.beats(k).ops;
        shrink = false;
        for m = 1:numel(ops)
            [view, caption, card, shrink] = applyOp(app, ops{m}, view, caption, card, shrink);
        end
        beatLog(end+1, :) = {si, tNow + tb(k), sc.beats(k).anchor, view}; %#ok<AGROW>
        nF = round((tNow + tEnd)*fps) - nWritten;      % frames owed to this beat
        if shrink
            frames = shrinkFrames(app, nF, view, caption, card);
        else
            frames = {render(app, view, caption, card)};
        end
        writeShot(frames, nF);
    end
    tNow = tNow + T;
end

% end card
view = 'full';  caption = '';
card = {'Conjugate Point Explorer', 'MATLAB  \cdot  optimal\_control\_examples/conjugate\_points'};
nF = round((tNow + endCard)*fps) - nWritten;
writeShot({render(app, view, caption, card)}, nF);
track = [track; zeros(round(endCard*fsA), 1)];
tNow = tNow + endCard;
close(vw);
if numel(unique(srcExt)) > 1
    % a misnamed final-voice file silently falls back to the draft: say so
    warning('make_conjugate_video:mixedVoices', ['scenes use different audio types (%s): ' ...
            'check audio/ for a misnamed file (must be sceneNN.mp3)'], strjoin(srcExt, ' '));
end

audFile = fullfile(outDir, ['narration' tag '.wav']);
peak = max(abs(track));
if peak > 0, track = 0.9*track/peak; end
audiowrite(audFile, track, fsA);

finalFile = fullfile(outDir, ['conjugate_points_video' tag '.mp4']);
ff = findFfmpeg();
if isempty(ff)
    finalFile = '';
    fprintf(['\nffmpeg not found: install it (brew install ffmpeg) and run mux_video.sh,\n' ...
             'or join %s and %s in iMovie.\n'], vidFile, audFile);
else
    cmd = sprintf('"%s" -y -loglevel error -i "%s" -i "%s" -c:v copy -c:a aac -b:a 192k -shortest "%s"', ...
                  ff, vidFile, audFile, finalFile);
    [st, msg] = system(cmd);
    if st ~= 0, warning('make_conjugate_video:ffmpeg', 'ffmpeg failed: %s', msg); finalFile = ''; end
end
out = struct('video', vidFile, 'audio', audFile, 'final', finalFile, 'duration', tNow, ...
             'beats', cell2table(beatLog, 'VariableNames', {'scene', 't', 'anchor', 'view'}));
fprintf('\nvideo %s\naudio %s\nfinal %s\nduration %.1f s (%d frames at %d fps)\n', ...
        vidFile, audFile, finalFile, tNow, nWritten, fps);

    % -------------------------------------------------------- nested
    function writeShot(frames, nF)
        % write nF frames: cross-fade from the previous shot, then the new
        % frame(s) spread evenly over the remainder
        if nF <= 0, return, end
        nFade = 0;
        if ~isempty(prev), nFade = min(round(fadeS*fps), nF - 1); end
        for q = 1:nF
            fr = frames{min(numel(frames), 1 + floor((q-1)*numel(frames)/nF))};
            if q <= nFade
                a = q/(nFade + 1);
                fr = uint8((1 - a)*double(prev) + a*double(fr));
            end
            writeVideo(vw, fr);
            nWritten = nWritten + 1;
        end
        prev = frames{end};
    end

    function img = render(app, view, caption, card)
        % one frame: title card, or the window / a panel plus a caption bar
        if ~isempty(card)
            img = cardImage(card{1}, card{2}, W, H);
            return
        end
        drawnow;                                       % never capture mid-redraw
        if strcmp(view, 'full')
            exportapp(app.fig, tmpPng);
            img = fitCanvas(imread(tmpPng), W, H, [245 245 245]);
        else
            img = closeUp(app.axes.(view), W, H - capH);
            img = [img; repmat(uint8(255), capH, W, 3)];
        end
        if ~isempty(caption)
            bar = captionImage(caption, W, capH);
            rows = H - capH + 1:H;
            img(rows, :, :) = uint8(0.12*double(img(rows, :, :)) + 0.88*double(bar));
        end
    end

    function frames = shrinkFrames(app, nF, view, caption, card)
        % animate delta -> 0: 20 captures over the first 75% of the beat,
        % then hold the last
        st0 = app.state();
        d0 = 0.6;  if ~isempty(st0.fan), d0 = st0.fan(1).delta; end
        nSteps = 20;
        frames = cell(1, 0);
        for q = 1:nSteps
            app.setDelta(d0*0.82^q);
            frames{end+1} = render(app, view, caption, card); %#ok<AGROW>
        end
        nHold = max(1, round(nSteps/3));
        frames = [frames, repmat(frames(end), 1, nHold)];
        frames = frames(1:min(numel(frames), max(1, nF)));
    end
end

% ---------------------------------------------------------------------------
function [view, caption, card, shrink] = applyOp(app, op, view, caption, card, shrink)
% APPLYOP  Apply one beat operation to the explorer or the shot state.
% INPUTS: app handle struct; op {name, arg}; current view/caption/card/shrink.
% OUTPUTS: the updated view, caption, card, shrink.
switch op{1}
    case 'preset',   app.setPreset(op{2});  card = {};
    case 'extremal'
        st = app.state();
        if numel(st.E) >= op{2}, app.selectExtremal(op{2}); end
    case 'delta',    app.setDelta(op{2});
    case 'scaled',   app.setScaled(op{2});
    case 'problem',  p = op{2};  app.setProblem(p{:});  card = {};
    case 'view',     view = op{2};  card = {};
    case 'caption',  caption = op{2};
    case 'card',     card = op{2};
    case 'shrink',   shrink = true;
    otherwise,       error('make_conjugate_video:op', 'unknown beat operation "%s"', op{1});
end
end

function f = anchorFrac(txt, anchor)
% ANCHORFRAC  Where the anchor phrase starts, as a fraction of the spoken
% text by characters (speech time is close to proportional to characters).
% INPUTS: txt narration; anchor phrase ('' = start). OUTPUTS: f in [0,1).
if isempty(anchor), f = 0; return, end
norm = @(s) regexprep(lower(regexprep(s, '''', '')), '[^a-z0-9]+', ' ');
nt = norm(txt);  na = strtrim(norm(anchor));
k = strfind(nt, na);
assert(~isempty(k), 'make_conjugate_video:anchor', 'anchor "%s" not found in the narration', anchor);
f = (k(1) - 1)/numel(nt);
end

function [y, dur, src] = sceneAudio(stem, fs)
% SCENEAUDIO  The scene's narration, mono at fs: .mp3 (final) preferred over
% .wav (draft). INPUTS: stem path without extension; fs. OUTPUTS: y column,
% dur seconds, src file name.
cand = {[stem '.mp3'], [stem '.m4a'], [stem '.wav'], [stem '.aiff']};
k = find(cellfun(@isfile, cand), 1);
assert(~isempty(k), 'make_conjugate_video:audio', ...
       'no audio for %s: run make_draft_voice.sh or add an mp3', stem);
[y, f0] = audioread(cand{k});
y = mean(y, 2);
if f0 ~= fs, y = resample(y, fs, f0); end
dur = numel(y)/fs;
[~, n, e] = fileparts(cand{k});  src = [n e];
end

function img = closeUp(src, W, H)
% CLOSEUP  Re-draw one explorer panel full-screen: copy its graphics into an
% offscreen figure of the frame's size, with presentation fonts and thicker
% lines (exporting the panel itself gives a thin strip for the wide panels).
% INPUTS: src uiaxes; W, H pixels. OUTPUTS: img [H x W x 3] uint8.
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', 'Position', [0 0 W H], ...
           'InvertHardcopy', 'off');
pinLight(f);
ax = axes(f, 'Position', [0.07 0.14 0.90 0.76], 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
          'GridColor', [0.15 0.15 0.15]);
copyobj(allchild(src), ax);
set(ax, 'XLim', src.XLim, 'YLim', src.YLim, 'YDir', src.YDir, 'XGrid', 'on', 'YGrid', 'on', ...
        'Box', 'on', 'FontSize', 20, 'FontWeight', 'bold', 'LineWidth', 1.5, 'Layer', 'top');
for L = findall(ax, 'Type', 'line').'
    L.LineWidth = 1.6*L.LineWidth;
    if ~strcmp(L.Marker, 'none'), L.MarkerSize = 1.6*L.MarkerSize; end
end
for C = findall(ax, 'Type', 'constantline').'
    C.LineWidth = 1.6*C.LineWidth;  C.FontSize = 18;
end
title(ax, src.Title.String, 'Interpreter', src.Title.Interpreter, 'FontSize', 22, 'Color', 'k');
xlabel(ax, src.XLabel.String, 'FontSize', 20, 'Color', 'k');
ylabel(ax, src.YLabel.String, 'FontSize', 20, 'Color', 'k');
img = imresize(print(f, '-RGBImage', '-r0'), [H W]);
close(f);
end

function pinLight(f)
% PINLIGHT  Force the light theme (R2025a+ figures follow the OS appearance,
% so a render after sunset came out with black plot areas). No-op on
% releases without figure themes. INPUTS: f figure or uifigure.
try
    f.Theme = 'light';
catch
end
end

function img = fitCanvas(I, W, H, bg)
% FITCANVAS  Scale I to fit W x H (keeping aspect) and centre it on a
% background colour. INPUTS: I image, W, H, bg [1x3]. OUTPUTS: img uint8.
if size(I, 3) == 1, I = repmat(I, 1, 1, 3); end
s = min(W/size(I, 2), H/size(I, 1));
I = imresize(I, s);
I = I(1:min(end, H), 1:min(end, W), :);
img = repmat(reshape(uint8(bg), 1, 1, 3), H, W);
r0 = floor((H - size(I, 1))/2);  c0 = floor((W - size(I, 2))/2);
img(r0 + (1:size(I, 1)), c0 + (1:size(I, 2)), :) = I;
end

function img = captionImage(txt, W, H)
% CAPTIONIMAGE  A dark caption bar with centred white TeX text.
% INPUTS: txt, W, H pixels. OUTPUTS: img [H x W x 3] uint8.
f = figure('Visible', 'off', 'Color', [0.08 0.10 0.16], 'Units', 'pixels', ...
           'Position', [0 0 W H], 'InvertHardcopy', 'off');
ax = axes(f, 'Position', [0 0 1 1], 'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1]);
text(ax, 0.5, 0.5, txt, 'Color', 'w', 'FontSize', 30, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Interpreter', 'tex');
img = imresize(print(f, '-RGBImage', '-r0'), [H W]);
close(f);
end

function img = cardImage(title, sub, W, H)
% CARDIMAGE  A full-screen title card. INPUTS: title, sub (TeX), W, H.
% OUTPUTS: img [H x W x 3] uint8.
f = figure('Visible', 'off', 'Color', [0.06 0.09 0.18], 'Units', 'pixels', ...
           'Position', [0 0 W H], 'InvertHardcopy', 'off');
ax = axes(f, 'Position', [0 0 1 1], 'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1]);
text(ax, 0.5, 0.56, title, 'Color', 'w', 'FontSize', 64, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Interpreter', 'tex');
line(ax, [0.4 0.6], [0.47 0.47], 'Color', [1 0.82 0.3], 'LineWidth', 5);
text(ax, 0.5, 0.39, sub, 'Color', [0.8 0.85 0.95], 'FontSize', 32, ...
     'HorizontalAlignment', 'center', 'Interpreter', 'tex');
img = imresize(print(f, '-RGBImage', '-r0'), [H W]);
close(f);
end

function ff = findFfmpeg()
% FINDFFMPEG  Path of ffmpeg, or '' (MATLAB's PATH may lack Homebrew).
% INPUTS: none. OUTPUTS: ff char.
ff = '';
for c = {'/opt/homebrew/bin/ffmpeg', '/usr/local/bin/ffmpeg'}
    if isfile(c{1}), ff = c{1}; return, end
end
[st, p] = system('command -v ffmpeg');
if st == 0, ff = strtrim(p); end
end

function v = getfielddef(s, f, v)
% GETFIELDDEF  s.f if present, else the default v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
end
