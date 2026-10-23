function file = movie_quad(t, X, info, stem, opts)
%% Purpose:
%
%   A 1280 x 720 MP4 of a quadrotor flight in the user's frame: the
%   vehicle as a "+" of two body arms (rotor 1 marked, arms scaled so they
%   read at the path's scale) with its thrust axis, over the path flown so
%   far, a fixed camera and fixed axis limits (a frame that rescales reads
%   as motion that is not there), with time, position and attitude in the
%   subtitle. Frames are forced to an exact 1280 x 720: H.264 shears a frame
%   whose size is not a multiple of 16 into diagonal streaks.
%
%% Inputs:
%
%  t                        [1 x N]                 Times (s)
%
%  X                        [12 x N]                States, LIBRARY frame
%
%  info                     struct                  .frame 'ENU' | 'NED';
%                                                   .lla0; .label
%
%  stem                     char                    Output path without
%                                                   extension
%
%  opts                     struct (optional)       .fps [30]; .seconds
%                                                   [t(end)] movie length
%                                                   (real time by default)
%
%% Outputs:
%
%  file                     char                    The MP4 written
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: two seconds of hover-and-climb:
     p = drone_params('rigid6dof');  h = p.m*p.g/4;
     [t, X] = simulate_quad([0; 0; 1; zeros(9, 1)], struct('t', [0 1], 'f', h*[1 1 1 1; 1.1 1.1 1.1 1.1].'), 2, p);
     file = movie_quad(t, X, struct('frame', 'ENU', 'lla0', [33.77 -118.39 50], 'label', 'demo'), fullfile(tempdir, 'demo'));
     fprintf('wrote %s\n', file);
     return
end

if nargin < 5 || isempty(opts), opts = struct(); end
fps     = getf(opts, 'fps', 30);
seconds = getf(opts, 'seconds', t(end));
frame   = upper(info.frame);
p       = drone_params('rigid6dof');

nFrame = max(2, round(fps*seconds));
tMov   = linspace(t(1), t(end), nFrame);
Xi     = interp1(t, X.', tMov, 'pchip').';
Xu     = quad_frames('from_lib', Xi, frame);
switch frame
    case 'ENU', axName = {'east', 'north', 'up'};    dirs = {'normal', 'normal'};  Pf = eye(3);
    case 'NED', axName = {'north', 'east', 'down'};  dirs = {'reverse', 'reverse'}; Pf = [0 1 0; 1 0 0; 0 0 -1];
end

% fixed limits with a margin that also fits the drawn arms
span   = max(max(max(Xu(1:3, :), [], 2) - min(Xu(1:3, :), [], 2)), 1);
armLen = max(0.12*span, 4*p.l);
% per-axis limits at true scale (DataAspectRatio 1:1:1), each axis at least
% wide enough to show the vehicle, so a level flight is not drawn in a cube
% that is mostly empty sky
lo = min(Xu(1:3, :), [], 2) - armLen - 0.05*span;
hi = max(Xu(1:3, :), [], 2) + armLen + 0.05*span;
mid = (lo + hi)/2;  half = max((hi - lo)/2, 1.5*armLen);
lim = [mid - half, mid + half];

fig = figure('Color', 'k', 'Position', [100 100 1280 720], 'InvertHardcopy', 'off');
ax  = axes(fig, 'Position', [0.06 0.08 0.88 0.80]);
frames(nFrame) = struct('cdata', [], 'colormap', []);
for k = 1:nFrame
    cla(ax);  hold(ax, 'on');
    plot3(ax, Xu(1, 1:k), Xu(2, 1:k), Xu(3, 1:k), '-', 'Color', [0.45 0.70 1.00 0.8], 'LineWidth', 1.6);
    c = Xu(1:3, k);
    R = Pf*euler_zyx(Xi(7:9, k));                    % library body axes in the display frame
    a1 = c + armLen*R(:, 1);  a3 = c - armLen*R(:, 1);
    a2 = c + armLen*R(:, 2);  a4 = c - armLen*R(:, 2);
    tz = c + 0.8*armLen*R(:, 3);
    plot3(ax, [a1(1) a3(1)], [a1(2) a3(2)], [a1(3) a3(3)], '-', 'Color', [1 0.75 0.3], 'LineWidth', 3);
    plot3(ax, [a2(1) a4(1)], [a2(2) a4(2)], [a2(3) a4(3)], '-', 'Color', [1 0.75 0.3], 'LineWidth', 3);
    plot3(ax, [c(1) tz(1)], [c(2) tz(2)], [c(3) tz(3)], '-', 'Color', [0.5 0.9 1], 'LineWidth', 2);
    plot3(ax, a1(1), a1(2), a1(3), 'o', 'MarkerFaceColor', [1 0.3 0.3], 'MarkerEdgeColor', 'none', 'MarkerSize', 7);
    for a = {a2, a3, a4}
        plot3(ax, a{1}(1), a{1}(2), a{1}(3), 'o', 'MarkerFaceColor', [1 0.75 0.3], 'MarkerEdgeColor', 'none', 'MarkerSize', 6);
    end
    set(ax, 'Color', 'k', 'XColor', [0.75 0.75 0.75], 'YColor', [0.75 0.75 0.75], 'ZColor', [0.75 0.75 0.75], ...
            'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.35, 'XLim', lim(1, :), 'YLim', lim(2, :), 'ZLim', lim(3, :), ...
            'ZDir', dirs{1}, 'YDir', dirs{2}, 'DataAspectRatio', [1 1 1]);
    grid(ax, 'on');  view(ax, 35, 22);
    xlabel(ax, [axName{1} '  [m]'], 'Color', [0.8 0.8 0.8]);
    ylabel(ax, [axName{2} '  [m]'], 'Color', [0.8 0.8 0.8]);
    zlabel(ax, [axName{3} '  [m]'], 'Color', [0.8 0.8 0.8]);

    annotation(fig, 'textbox', [0 0.935 1 0.05], 'String', sprintf('%s   (%s frame, origin %.5f, %.5f, %.1f m)', ...
        info.label, frame, info.lla0(1), info.lla0(2), info.lla0(3)), 'Color', 'w', 'FontSize', 18, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'Interpreter', 'none', 'Tag', 'ttl');
    annotation(fig, 'textbox', [0 0.895 1 0.04], 'String', sprintf( ...
        't =%6.2f s   pos [%7.2f %7.2f %7.2f] m   roll/pitch/yaw [%6.1f %6.1f %6.1f] deg', ...
        tMov(k), Xu(1, k), Xu(2, k), Xu(3, k), Xu(7:9, k)*180/pi), 'Color', [0.72 0.72 0.72], 'FontSize', 13, ...
        'FontName', 'FixedWidth', 'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'Tag', 'sub');
    drawnow;
    frames(k) = getframe(fig);
    delete(findall(fig, 'Type', 'textboxshape'));
end
close(fig);

for k = 1:nFrame
    if ~isequal(size(frames(k).cdata, 1:2), [720 1280])
        frames(k).cdata = imresize(frames(k).cdata, [720 1280]);
    end
end
file = [stem '.mp4'];
v = VideoWriter(file, 'MPEG-4');
v.FrameRate = fps;  v.Quality = 95;
open(v);
for k = 1:nFrame, writeVideo(v, frames(k)); end
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
