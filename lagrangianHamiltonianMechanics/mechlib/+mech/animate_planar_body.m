function file = animate_planar_body(t, pose, body, stem, opts)
%% Purpose:
%
%   A 1280 x 720 top-view MP4 of a planar rigid body: its outline, markers
%   and decorative segments given once in body coordinates and carried by
%   the pose (x, y, theta) of a reference point; arrows and a subtitle
%   supplied by the caller as functions of time; traces that grow with
%   time (drawn from the caller's own samples plus the current
%   interpolated point, so they stay smooth at any fps); a fixed camera
%   and fixed, equal-aspect limits covering the whole motion (a frame that
%   rescales reads as motion that is not there).
%   Frames are forced to an exact 1280 x 720: H.264 shears a frame whose
%   size is not a multiple of 16 into diagonal streaks. Knows no system.
%
%% Inputs:
%
%  t                        [1 x N]                 Times (s), strictly
%                                                   increasing, N >= 2
%
%  pose                     [3 x N]                 (x, y, theta) of the
%                                                   reference point in the
%                                                   world frame; theta
%                                                   counterclockwise from x
%
%  body                     struct                  .outline [2 x P] body
%                                                   coordinates (closed by
%                                                   this function);
%                                                   .markers struct array
%                                                   (.name, .pos [2 x 1],
%                                                   .color); .segments
%                                                   struct array (.pts [2
%                                                   x 2] endpoints as
%                                                   columns, .color,
%                                                   .width); .arrows struct
%                                                   array (.name, .fun =
%                                                   @(tk) -> [base; vec]
%                                                   [4 x 1] WORLD, .color);
%                                                   .subtitle = @(tk) ->
%                                                   char. All but .outline
%                                                   optional
%
%  stem                     char                    Output path without
%                                                   extension
%
%  opts                     struct (optional)       .fps [30]; .seconds
%                                                   [t(end) - t(1)];
%                                                   .traces [2 x N x T]
%                                                   world paths on t
%                                                   [default: the
%                                                   reference point];
%                                                   .traceColors [T x 3];
%                                                   .title ['']
%
%% Outputs:
%
%  file                     char                    The MP4 written
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: a rectangle sliding on a circle, two seconds, to tempdir
     t = linspace(0, 2, 61);
     pose = [cos(pi*t); sin(pi*t); pi*t + pi/2];
     body.outline = [-0.15 0.15 0.15 -0.15; -0.08 -0.08 0.08 0.08];
     body.markers = struct('name', 'nose', 'pos', [0.15; 0], 'color', [1 0.3 0.3]);
     body.subtitle = @(tk) sprintf('t = %.2f s', tk);
     file = mech.animate_planar_body(t, pose, body, fullfile(tempdir, 'demo_planar_body'), struct('title', 'demo'));
     fprintf('wrote %s\n', file);
     return
end

if nargin < 5 || isempty(opts), opts = struct(); end
t = t(:).';
N = numel(t);
if N < 2 || any(~isfinite(t)) || any(diff(t) <= 0)
    error('mech:animate_planar_body:t', 't must be finite and strictly increasing with at least two entries (got %d)', N);
end
if ~isequal(size(pose), [3 N])
    error('mech:animate_planar_body:shape', 'pose must be [3 x N] with N = numel(t) = %d (got %d x %d)', N, size(pose, 1), size(pose, 2));
end
if ~isfield(body, 'outline') || size(body.outline, 1) ~= 2
    error('mech:animate_planar_body:shape', 'body.outline must be [2 x P]');
end
markers  = getf(body, 'markers',  struct('name', {}, 'pos', {}, 'color', {}));
segments = getf(body, 'segments', struct('pts', {}, 'color', {}, 'width', {}));
arrows   = getf(body, 'arrows',   struct('name', {}, 'fun', {}, 'color', {}));
subtitle = getf(body, 'subtitle', @(tk) sprintf('t = %.2f s', tk));
fps      = getf(opts, 'fps', 30);
seconds  = getf(opts, 'seconds', t(end) - t(1));
ttl      = getf(opts, 'title', '');
traces   = getf(opts, 'traces', reshape(pose(1:2, :), 2, N, 1));
if ~isequal(size(traces, 1:2), [2 N])
    error('mech:animate_planar_body:shape', 'opts.traces must be [2 x N x T] with N = %d (got %s)', N, mat2str(size(traces)));
end
nTr      = size(traces, 3);
trColors = getf(opts, 'traceColors', repmat([0.45 0.70 1.00], nTr, 1));

nFrame = max(2, round(fps*seconds));
tMov   = linspace(t(1), t(end), nFrame);
P      = interp1(t, pose.', tMov, 'pchip').';                 % [3 x nFrame]
Tr     = zeros(2, nFrame, nTr);
for tt = 1:nTr
    Tr(:, :, tt) = interp1(t, traces(:, :, tt).', tMov, 'pchip').';
end

% body size and fixed limits: every trace, every reference point, plus the
% body's own reach, with a margin; equal aspect at 16:9 so the frame fills
pts   = [body.outline, [markers.pos], reshape([segments.pts], 2, [])];
if isempty(pts), pts = zeros(2, 1); end
reach = max([sqrt(sum(pts.^2, 1)), eps]);
allXY = [reshape(Tr, 2, []), P(1:2, :)];
lo    = min(allXY, [], 2) - 1.2*reach;
hi    = max(allXY, [], 2) + 1.2*reach;
span  = max(hi - lo);
lo    = lo - 0.05*span;  hi = hi + 0.05*span;
mid   = (lo + hi)/2;
half  = max((hi - lo)/2, 1.5*reach);
if half(1)/half(2) < 16/9, half(1) = half(2)*16/9; else, half(2) = half(1)*9/16; end
lim   = [mid - half, mid + half];                             % [2 x 2] rows x, y

fig = figure('Color', 'k', 'Position', [100 100 1280 720], 'InvertHardcopy', 'off', 'Visible', 'off');
ax  = axes(fig, 'Position', [0.06 0.08 0.88 0.80]);
frames(nFrame) = struct('cdata', [], 'colormap', []);
closedOutline = [body.outline, body.outline(:, 1)];
for kk = 1:nFrame
    cla(ax);  hold(ax, 'on');
    for tt = 1:nTr
        upTo = t <= tMov(kk);
        plot(ax, [traces(1, upTo, tt), Tr(1, kk, tt)], [traces(2, upTo, tt), Tr(2, kk, tt)], '-', 'Color', [trColors(tt, :) 0.8], 'LineWidth', 1.6);
    end
    c  = P(1:2, kk);  th = P(3, kk);
    R  = [cos(th), -sin(th); sin(th), cos(th)];
    O  = c + R*closedOutline;
    patch(ax, O(1, :), O(2, :), [1 0.75 0.3], 'FaceAlpha', 0.35, 'EdgeColor', [1 0.75 0.3], 'LineWidth', 2);
    for ss = 1:numel(segments)
        S = c + R*segments(ss).pts;
        plot(ax, S(1, :), S(2, :), '-', 'Color', segments(ss).color, 'LineWidth', segments(ss).width);
    end
    for mm = 1:numel(markers)
        Mk = c + R*markers(mm).pos;
        plot(ax, Mk(1), Mk(2), 'o', 'MarkerFaceColor', markers(mm).color, 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
    end
    for aa = 1:numel(arrows)
        bv = arrows(aa).fun(tMov(kk));
        quiver(ax, bv(1), bv(2), bv(3), bv(4), 0, 'Color', arrows(aa).color, 'LineWidth', 2, 'MaxHeadSize', 0.6);
    end
    set(ax, 'Color', 'k', 'XColor', [0.75 0.75 0.75], 'YColor', [0.75 0.75 0.75], ...
            'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.35, 'XLim', lim(1, :), 'YLim', lim(2, :), ...
            'DataAspectRatio', [1 1 1]);
    grid(ax, 'on');
    xlabel(ax, 'x  [m]', 'Color', [0.8 0.8 0.8]);
    ylabel(ax, 'y  [m]', 'Color', [0.8 0.8 0.8]);
    annotation(fig, 'textbox', [0 0.935 1 0.05], 'String', ttl, 'Color', 'w', 'FontSize', 18, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'Interpreter', 'none');
    annotation(fig, 'textbox', [0 0.895 1 0.04], 'String', subtitle(tMov(kk)), 'Color', [0.72 0.72 0.72], ...
        'FontSize', 13, 'FontName', 'FixedWidth', 'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'Interpreter', 'none');
    drawnow;
    frames(kk) = getframe(fig);
    delete(findall(fig, 'Type', 'textboxshape'));
end
close(fig);

for kk = 1:nFrame
    if ~isequal(size(frames(kk).cdata, 1:2), [720 1280])
        frames(kk).cdata = imresize(frames(kk).cdata, [720 1280]);
    end
end
file = [stem '.mp4'];
v = VideoWriter(file, 'MPEG-4');
v.FrameRate = fps;  v.Quality = 95;
open(v);
for kk = 1:nFrame, writeVideo(v, frames(kk)); end
close(v);
end

% ------------------------------------------------------------------------
function v = getf(s, name, default)
%% Purpose:
%
%   Field with a default.
%
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
