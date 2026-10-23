function fig = plot_drones(out, opts)
%% Purpose:
%
%   The paper's Figure 4 view of one solve: four panels at t = 0, T/3,
%   2T/3, T, paths as faint grey curves, agents as circles of radius r.
%   2D layouts are drawn in the plane; 3D ones in a 3D axes with the ball
%   when present.
%
%% Inputs:
%
%  out                      struct                  From solve_drones_direct
%
%  opts                     struct (optional)       .visible [true] .file
%                                                   [''] PNG to export
%
%% Outputs:
%
%  fig                      handle                  The figure
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the M = 4 fixture:
     here = fileparts(mfilename('fullpath'));
     S = load(fullfile(here, 'data', 'drones_direct_M4.mat'));
     plot_drones(S.out);
     return
end
if nargin < 2, opts = struct(); end
vis = true;  if isfield(opts, 'visible'), vis = opts.visible; end
file = '';   if isfield(opts, 'file'), file = opts.file; end
geom = out.geom;  M = geom.M;
[Xf, ~, tf] = hermite_resample(out.X, out.V, out.tN, 5);
snaps = [0, 1/3, 2/3, 1]*out.T;
cols = lines(M);
th = linspace(0, 2*pi, 40);
lim = max(abs(Xf(:))) + 1.5*geom.r;
fig = figure('color', [1 1 1], 'Position', [100 100 1600 420], 'Visible', tern(vis, 'on', 'off'));
for s = 1:4
    ax = subplot(1, 4, s);  hold(ax, 'on');
    [~, k] = min(abs(tf - snaps(s)));
    for a = 1:M
        ia = 3*(a-1) + (1:3);
        if geom.d == 2
            plot(ax, Xf(ia(1), :), Xf(ia(2), :), '-', 'Color', [0.8 0.8 0.8]);
            plot(ax, Xf(ia(1), k) + geom.r*cos(th), Xf(ia(2), k) + geom.r*sin(th), '-', 'Color', cols(a, :), 'LineWidth', 1.2);
        else
            plot3(ax, Xf(ia(1), :), Xf(ia(2), :), Xf(ia(3), :), '-', 'Color', [0.8 0.8 0.8]);
            plot3(ax, Xf(ia(1), k), Xf(ia(2), k), Xf(ia(3), k), 'o', 'Color', cols(a, :), 'MarkerFaceColor', cols(a, :));
        end
    end
    if geom.d == 3
        view(ax, 3);
        if ~isempty(geom.ball)
            [sx, sy, sz] = sphere(24);
            surf(ax, geom.ball.c(1) + geom.ball.R*sx, geom.ball.c(2) + geom.ball.R*sy, geom.ball.c(3) + geom.ball.R*sz, ...
                 'FaceColor', [0.6 0.6 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        end
    end
    axis(ax, 'equal');  grid(ax, 'on');
    title(ax, sprintf('t = %.2f', snaps(s)));
    if geom.d == 2, xlim(ax, [-lim lim]); ylim(ax, [-lim lim]); end
end
sgtitle(fig, sprintf('M = %d, %s, J = %.4f (J_{LB} %.4f)', M, out.mode, out.J, out.JLB));
if ~isempty(file), exportgraphics(fig, file, 'Resolution', 120); end
end

% ------------------------------------------------------------------------
function s = tern(c, a, b)
%% Purpose:
%
%   a if c, else b.
%
if c, s = a; else, s = b; end
end
