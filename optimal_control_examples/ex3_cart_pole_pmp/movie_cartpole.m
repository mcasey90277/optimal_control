function movie_cartpole(opts)
%% Purpose:
%
%   Animate the swing-up both ways at once: the DIRECT collocation fixture
%   and the INDIRECT (Pontryagin) solve, side by side on a shared clock,
%   with the two control histories traced underneath. Built on the animation
%   in ex2_cart_pole_swing_up/try2 (same cart, rod and bob drawing, same
%   geometry: bob at cart + L*[sin q2; -cos q2]).
%
%   The point of the pairing: these are two different methods -- discretize
%   then optimize, versus write the optimality conditions then solve them --
%   landing on the same motion. They agree to 0.062% in cost and 0.42% RMS
%   in control, and the frame shows that as one cart, not two.
%
%  ASSUMPTIONS / NOTES:
%
% • Frames are forced to an EXACT 1280 x 720: H.264 shears a frame whose
%   dimensions are not multiples of 16, which shows up as diagonal coloured
%   streaks in the .mp4 (and is a frame-size artefact, not a render bug).
% • Both .mp4 and .gif are written; the .gif usually looks better.
% • The indirect solve runs here (a few seconds) unless .out is supplied.
%
%% Inputs:
%
%  opts                     struct (optional)       .fps [30], .seconds [10]
%                                                   playback length, .out
%                                                   (a run_cartpole_pmp
%                                                   output struct, to skip
%                                                   re-solving), .file
%                                                   [cartpole_direct_vs_indirect]
%
%% Outputs:
%
%  none (writes <file>.mp4 and <file>.gif beside this function)
%
%% Revision History:
%  M. Casey                                                   (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'cartpole_common'));
d = @(f, v) fieldd(opts, f, v);
fps     = d('fps', 30);
seconds = d('seconds', 10);
stem    = d('file', fullfile(here, 'cartpole_direct_vs_indirect'));

R = load(fullfile(here, 'data', 'cartpole_direct_ref.mat'));
p = R.p;  tf = R.tf;  L = p.L;
out = d('out', []);
if isempty(out), out = run_cartpole_pmp(struct('plot', false)); end
assert(out.ok, 'movie_cartpole:notSolved', 'the indirect solve is not a solution: %s', out.why);

%% One clock for both: sample each solution at the same instants
nFrame = round(fps*seconds);
tMov   = linspace(0, tf, nFrame);
Xd = interp1(R.tN,  R.X.',  tMov, 'pchip').';          % direct
Ud = interp1(R.tN,  R.U,    tMov, 'pchip');
Xi = interp1(out.t, out.X.', tMov, 'pchip').';         % indirect
Ui = interp1(out.t, out.U,   tMov, 'pchip');

%% Fixed limits for every frame: a cart that rescales its own axes reads as
%  motion that is not there
cartW = 0.8;  cartH = 0.4;
xAll = [Xd(1,:), Xi(1,:)];
xPad = max(L, 0.6) + cartW;
xLim = [min(xAll) - xPad, max(xAll) + xPad];
yLim = [-(L + 0.6), L + 0.6];
uLim = [min([Ud, Ui]), max([Ud, Ui])] + [-5, 5];

fig = figure('Color', 'k', 'Position', [100 100 1280 720], 'InvertHardcopy', 'off');
axD = axes(fig, 'Position', [0.045 0.42 0.44 0.48]);
axI = axes(fig, 'Position', [0.525 0.42 0.44 0.48]);
axU = axes(fig, 'Position', [0.075 0.09 0.885 0.24]);

frames(nFrame) = struct('cdata', [], 'colormap', []);
for k = 1:nFrame
    drawCart(axD, Xd(:,k), p, xLim, yLim, cartW, cartH, [0.45 0.70 1.00], ...
             sprintf('DIRECT  collocation   ·   J = %.1f', R.J));
    drawCart(axI, Xi(:,k), p, xLim, yLim, cartW, cartH, [1.00 0.72 0.30], ...
             sprintf('INDIRECT  Pontryagin BVP   ·   J = %.1f', out.J));

    cla(axU);  hold(axU, 'on');
    plot(axU, tMov, Ud, '-',  'Color', [0.45 0.70 1.00 0.55], 'LineWidth', 1.2);
    plot(axU, tMov, Ui, '--', 'Color', [1.00 0.72 0.30 0.55], 'LineWidth', 1.2);
    plot(axU, tMov(1:k), Ud(1:k), '-',  'Color', [0.45 0.70 1.00], 'LineWidth', 2.0);
    plot(axU, tMov(1:k), Ui(1:k), '--', 'Color', [1.00 0.72 0.30], 'LineWidth', 2.0);
    plot(axU, tMov(k), Ud(k), 'o', 'MarkerFaceColor', [0.45 0.70 1.00], ...
         'MarkerEdgeColor', 'none', 'MarkerSize', 7);
    set(axU, 'Color', 'k', 'XColor', [0.75 0.75 0.75], 'YColor', [0.75 0.75 0.75], ...
             'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.3, 'XLim', [0 tf], 'YLim', uLim);
    grid(axU, 'on');
    xlabel(axU, 'time  [s]', 'Color', [0.8 0.8 0.8]);
    ylabel(axU, 'u  [N]',    'Color', [0.8 0.8 0.8]);
    legend(axU, {'direct', 'indirect'}, 'TextColor', 'w', 'Color', [0.15 0.15 0.15], ...
           'EdgeColor', [0.3 0.3 0.3], 'Location', 'northeast', 'Orientation', 'horizontal');

    % FixedWidth so the numbers do not jitter the title from frame to frame
    annotation(fig, 'textbox', [0 0.945 1 0.05], 'String', ...
        'Two ways to swing up a pendulum, landing on the same motion', ...
        'Color', 'w', 'FontSize', 19, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'Tag', 'ttl');
    annotation(fig, 'textbox', [0 0.905 1 0.04], 'String', sprintf(...
        ['t =%5.2f s of %.0f   ·   cost within 0.062%%   ·   control within 0.42%% RMS   ·   ', ...
         'terminal miss %.1e'], tMov(k), tf, out.missTerminal), ...
        'Color', [0.72 0.72 0.72], 'FontSize', 13, 'FontName', 'FixedWidth', ...
        'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'Tag', 'sub');
    annotation(fig, 'textbox', [0.70 0.005 0.29 0.028], 'String', ...
        'Coorbital  |  collocation\_examples/ex3\_cart\_pole\_pmp', ...
        'Color', [0.45 0.45 0.45], 'FontSize', 10, 'HorizontalAlignment', 'right', ...
        'EdgeColor', 'none', 'Tag', 'brand');

    drawnow;
    frames(k) = getframe(fig);
    delete(findall(fig, 'Type', 'textboxshape'));
end
close(fig);

%% EXACT 1280 x 720 on every frame, or H.264 shears them into diagonal streaks
for k = 1:nFrame
    if ~isequal(size(frames(k).cdata, 1:2), [720 1280])
        frames(k).cdata = imresize(frames(k).cdata, [720 1280]);
    end
end

v = VideoWriter([stem '.mp4'], 'MPEG-4');
v.FrameRate = fps;  v.Quality = 98;
open(v);
for k = 1:nFrame, writeVideo(v, frames(k)); end
close(v);

for k = 1:nFrame
    [A, map] = rgb2ind(frames(k).cdata, 256);
    if k == 1
        imwrite(A, map, [stem '.gif'], 'gif', 'LoopCount', inf, 'DelayTime', 1/fps);
    else
        imwrite(A, map, [stem '.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', 1/fps);
    end
end

fprintf('movie_cartpole: %d frames at %d fps -> %s.mp4 and %s.gif\n', nFrame, fps, stem, stem);
end

% ------------------------------------------------------------------------
function drawCart(ax, x, p, xLim, yLim, cartW, cartH, col, ttl)
%% Purpose:
%
%   One cart-pole pose: rail, cart, rod, bob. Geometry as in ex2 -- the bob
%   sits at cart + L*[sin q2; -cos q2], so q2 = 0 hangs and q2 = pi is up.
%
cla(ax);  hold(ax, 'on');
set(ax, 'Color', 'k', 'XColor', [0.6 0.6 0.6], 'YColor', [0.6 0.6 0.6], ...
        'GridColor', [0.4 0.4 0.4], 'GridAlpha', 0.25);
cart = [x(1); 0];
bob  = cart + p.L*[sin(x(2)); -cos(x(2))];
plot(ax, xLim, [-cartH/2 -cartH/2], '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 2);
rectangle(ax, 'Position', [cart(1)-cartW/2, cart(2)-cartH/2, cartW, cartH], ...
          'FaceColor', col*0.65, 'EdgeColor', col, 'LineWidth', 1.5, 'Curvature', 0.15);
plot(ax, [cart(1) bob(1)], [cart(2) bob(2)], '-', 'Color', col, 'LineWidth', 3);
plot(ax, bob(1), bob(2), 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', ...
     'MarkerSize', 14, 'LineWidth', 1);
set(ax, 'XLim', xLim, 'YLim', yLim, 'DataAspectRatio', [1 1 1]);
grid(ax, 'on');
title(ax, ttl, 'Color', col, 'FontSize', 13, 'FontWeight', 'bold');
end

function v = fieldd(s, f, v0)
%% Purpose:
%
%   s.(f) if present, else the default.
%
if isfield(s, f), v = s.(f); else, v = v0; end
end
