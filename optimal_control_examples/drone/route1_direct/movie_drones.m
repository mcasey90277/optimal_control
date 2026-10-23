function movie_drones(out, file, opts)
%% Purpose:
%
%   An MP4 of one solve: agents as circles moving along faint paths. Every
%   frame is 1280 x 720 (a multiple of 16; H.264 shears any other size
%   into diagonal streaks -- memory matlab-movie-diagonal-streaks).
%
%% Inputs:
%
%  out                      struct                  From solve_drones_direct
%
%  file                     char                    Output .mp4
%
%  opts                     struct (optional)       .fps [30] .nFrames
%                                                   [300]
%
%% Outputs:
%
%  none (writes file)
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the M = 4 fixture, a short clip in the temp folder:
     here = fileparts(mfilename('fullpath'));
     S = load(fullfile(here, 'data', 'drones_direct_M4.mat'));
     movie_drones(S.out, fullfile(tempdir, 'drones_demo'), struct('nFrames', 60));
     return
end
if nargin < 3, opts = struct(); end
fps = 30;  if isfield(opts, 'fps'), fps = opts.fps; end
nF = 300;  if isfield(opts, 'nFrames'), nF = opts.nFrames; end
geom = out.geom;  M = geom.M;
[Xf, ~, tf] = hermite_resample(out.X, out.V, out.tN, 5);
tt = linspace(0, out.T, nF);
cols = lines(M);  th = linspace(0, 2*pi, 40);
lim = max(abs(Xf(:))) + 1.5*geom.r;
fig = figure('color', [1 1 1], 'Visible', 'off', 'Position', [100 100 1280 720]);
ax = axes(fig);  hold(ax, 'on');  axis(ax, 'equal');  grid(ax, 'on');
if geom.d == 2
    xlim(ax, [-lim lim]);  ylim(ax, [-lim lim]);
else
    view(ax, 3);
    if ~isempty(geom.ball)
        [sx, sy, sz] = sphere(24);
        surf(ax, geom.ball.c(1) + geom.ball.R*sx, geom.ball.c(2) + geom.ball.R*sy, geom.ball.c(3) + geom.ball.R*sz, ...
             'FaceColor', [0.6 0.6 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    end
end
for a = 1:M
    ia = 3*(a-1) + (1:3);
    if geom.d == 2, plot(ax, Xf(ia(1), :), Xf(ia(2), :), '-', 'Color', [0.85 0.85 0.85]);
    else, plot3(ax, Xf(ia(1), :), Xf(ia(2), :), Xf(ia(3), :), '-', 'Color', [0.85 0.85 0.85]);
    end
end
hA = gobjects(1, M);
for a = 1:M
    if geom.d == 2, hA(a) = plot(ax, NaN, NaN, '-', 'Color', cols(a, :), 'LineWidth', 1.5);
    else, hA(a) = plot3(ax, NaN, NaN, NaN, 'o', 'Color', cols(a, :), 'MarkerFaceColor', cols(a, :));
    end
end
vw = VideoWriter(file, 'MPEG-4');  vw.FrameRate = fps;  open(vw);
for f = 1:nF
    Xt = interp1(tf, Xf.', tt(f)).';
    for a = 1:M
        ia = 3*(a-1) + (1:3);
        if geom.d == 2, set(hA(a), 'XData', Xt(ia(1)) + geom.r*cos(th), 'YData', Xt(ia(2)) + geom.r*sin(th));
        else, set(hA(a), 'XData', Xt(ia(1)), 'YData', Xt(ia(2)), 'ZData', Xt(ia(3)));
        end
    end
    title(ax, sprintf('t = %.2f / %.0f', tt(f), out.T));
    fr = getframe(fig);
    im = imresize(fr.cdata, [720 1280]);
    writeVideo(vw, im);
end
close(vw);  close(fig);
end
