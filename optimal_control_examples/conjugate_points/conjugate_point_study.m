%% conjugate_point_study.m
%
%   The guide's numbers and figures, computed in one place, in order. Each
%   section states one fact about conjugate points, computes it here, and
%   prints PASS/FAIL against an answer that does not come from the code
%   being checked (a closed form, or a second instrument).
%
%     0  settings
%     1  the oscillator past pi: extremal, Jacobi field, the neighbours meet at pi
%     2  the pendulum: a finite slope change crosses NEAR t_c; delta -> 0 converges
%     3  the Morse count: eigenvalues of the second variation against b
%     4  the catenary: two extremals, one minimiser, and Lindelof's tangents
%     5  the catenary fold: at the critical half-width the conjugate point reaches b
%     6  Delta J: the negative mode lowers the cost; the positive one does not
%     7  figures for the guide (doc/figures/)
%
%   Run:  conjugate_point_study          (about 30 s; writes doc/figures/*.png)
%
%   References: Gelfand & Fomin, Calculus of Variations (1963), Ch. 5;
%   Bliss, Calculus of Variations (1925).
%
%   M. Casey  (c) 09/24/2026   Copyright Coorbital Inc.

%% 0. Settings
here = fileparts(mfilename('fullpath'));
addpath(here);
figDir = fullfile(here, 'doc', 'figures');
if ~isfolder(figDir), mkdir(figDir); end
makeFigures = true;
nFail = 0;

%% 1. The oscillator past pi: F = y'^2 - y^2 on [0, 4], y(0) = 0, y(4) = 1
osc = cov_problem('yp^2 - y^2', 0, 4, 0, 1);
Eo  = cov_extremals(osc, [-5 5]);
So  = cov_shoot(osc, Eo(1).p, 5);
c1 = numel(Eo) == 1 && abs(Eo(1).p - 1/sin(4)) < 1e-8;
fprintf('1a  one extremal, y''(0) = %.8f  (1/sin 4 = %.8f)          %s\n', Eo(1).p, 1/sin(4), pf(c1));
c2 = abs(So.tConj(1) - pi) < 1e-8;
fprintf('1b  Jacobi field h = sin t: first zero %.10f  (pi)           %s\n', So.tConj(1), pf(c2));
dOsc = [1 0.5 0.25];
tcx = arrayfun(@(d) firstOr(cov_perturbed(osc, Eo(1), d, 5).tCross), dOsc);
c3 = all(abs(tcx - pi) < 1e-7);
fprintf('1c  neighbours at delta = %s cross at %s   %s\n', mat2str(dOsc), mat2str(tcx, 8), pf(c3));
nFail = nFail + ~c1 + ~c2 + ~c3;

%% 2. The pendulum: F = y'^2/2 + cos y on [0, 5], y(0) = 0, y(5) = 0.5
pen = cov_problem('yp^2/2 + cos(y)', 0, 5, 0, 0.5);
Ep  = cov_extremals(pen, [-1.9 1.9]);
[~, iw] = max(arrayfun(@(e) numel(e.S.tConj), Ep));   % the extremal with t_c < b
tcP = Ep(iw).S.tConj(1);
dPen = 0.6*2.^-(0:5);
crossP = arrayfun(@(d) firstOr(cov_perturbed(pen, Ep(iw), d, 5).tCross), dPen);
errP = abs(crossP - tcP);
rat = errP(1:end-1)./errP(2:end);
c4 = numel(Ep) == 2 && tcP < pen.b;
fprintf('2a  two extremals; the second has t_c = %.6f < b = 5            %s\n', tcP, pf(c4));
c5 = all(diff(errP) < 0) && abs(rat(end) - 2) < 0.1;
fprintf('2b  crossing error vs delta %s\n    halving ratios %s -> 2 (first order)     %s\n', ...
        mat2str(errP, 3), mat2str(rat, 3), pf(c5));
nFail = nFail + ~c4 + ~c5;

%% 3. The Morse count: lambda_k(b) for the oscillator, b from 2 to 8
bGrid = linspace(2, 8, 61);
lam = nan(3, numel(bGrid));  nNeg = nan(1, numel(bGrid));
for k = 1:numel(bGrid)
    pr = cov_problem('yp^2 - y^2', 0, bGrid(k), 0, 1);
    if abs(sin(bGrid(k))) < 1e-3, continue, end           % y(b) = 1 unreachable at k pi
    E = cov_extremals(pr, [-60 60], 600);
    if isempty(E), continue, end
    V = cov_second_variation(pr, E(1), 300);
    lam(:,k) = V.lambda(1:3);  nNeg(k) = V.nNeg;
end
lamTrue = 2*(((1:3).'*pi./bGrid).^2 - 1);
ok = isfinite(nNeg);
c6 = max(abs(lam(:,ok) - lamTrue(:,ok))./max(1, abs(lamTrue(:,ok))), [], 'all') < 1e-3;
c7 = all(nNeg(ok) == floor(bGrid(ok)/pi));
fprintf('3a  lambda_1..3(b) vs 2((k pi/b)^2 - 1), %d values of b          %s\n', nnz(ok), pf(c6));
fprintf('3b  negative modes = floor(b/pi) = number of conjugate points     %s\n', pf(c7));
nFail = nFail + ~c6 + ~c7;

%% 4. The catenary: F = y sqrt(1 + y'^2) on [-1/2, 1/2], y = 1 at both ends
cat_ = cov_problem('y*sqrt(1+yp^2)', -0.5, 0.5, 1, 1);
[Ec, scanC] = cov_extremals(cat_, [-8 2]);
fc = @(c) c.*cosh(0.5./c) - 1;
cc = [fzero(fc, [0.5 1.5]), fzero(fc, [0.1 0.4])];     % shallow, deep
c8 = numel(Ec) == 2 && all(abs([Ec.p] + sinh(0.5./cc)) < 1e-8);
fprintf('4a  two catenaries, c = %.6f (shallow), %.6f (deep)            %s\n', cc, pf(c8));
c = cc(2);  yC = @(t) c*cosh(t/c);  ypC = @(t) sinh(t/c);
tStar = -0.5 - yC(-0.5)/ypC(-0.5);                     % tangent at a meets y = 0
tcL = fzero(@(t) t - yC(t)/ypC(t) - tStar, [0.01 0.49]);
Sd = cov_shoot(cat_, Ec(2).p, 0.5);
c9 = isempty(Ec(1).S.tConj) && abs(Sd.tConj(1) - tcL) < 1e-7;
fprintf('4b  deep: t_c = %.8f, Lindelof tangents %.8f; shallow: none  %s\n', Sd.tConj(1), tcL, pf(c9));
c10 = Ec(1).J < Ec(2).J;
fprintf('4c  area/2pi: shallow %.6f < deep %.6f                        %s\n', Ec(1).J, Ec(2).J, pf(c10));
nFail = nFail + ~c8 + ~c9 + ~c10;

%% 5. The catenary fold: widen the gap, L = half-width
% two catenaries exist for L < L* where x tanh x = 1, c* = 1/cosh x, L* = x c*
xs = fzero(@(x) x*tanh(x) - 1, [0.5 2]);
Lstar = xs/cosh(xs);
Ls = [0.50 0.60 0.65 0.66];
tcMinusB = nan(size(Ls));
for k = 1:numel(Ls)
    pr = cov_problem('y*sqrt(1+yp^2)', -Ls(k), Ls(k), 1, 1);
    E = cov_extremals(pr, [-8 2], 800);
    if numel(E) == 2
        S = cov_shoot(pr, E(2).p, Ls(k));
        tcMinusB(k) = Ls(k) - S.tConj(1);
    end
end
prX = cov_problem('y*sqrt(1+yp^2)', -0.67, 0.67, 1, 1);
c11 = isempty(cov_extremals(prX, [-8 2], 800)) && all(diff(tcMinusB) < 0);
fprintf('5   L* = %.6f; b - t_c on the deep catenary at L = %s: %s\n', Lstar, mat2str(Ls), mat2str(tcMinusB, 3));
fprintf('    -> shrinks to 0 at the fold; no catenary at L = 0.67 > L*     %s\n', pf(c11));
nFail = nFail + ~c11;

%% 6. Delta J along the lowest mode
Vs = cov_second_variation(cat_, Ec(1));
Vd = cov_second_variation(cat_, Ec(2));
Vo = cov_second_variation(osc, Eo(1));
c12 = all(Vs.dJ > 0) && all(Vd.dJ < 0) && all(Vo.dJ < 0);
fprintf('6   Delta J: shallow %s, deep %s, oscillator %s   %s\n', ...
        mat2str(Vs.dJ, 3), mat2str(Vd.dJ, 3), mat2str(Vo.dJ, 3), pf(c12));
nFail = nFail + ~c12;

fprintf('\n%d check(s) failed\n', nFail);

%% 7. Figures for the guide
if makeFigures
    col = struct('ext', [0.10 0.35 0.75], 'nb', [0.95 0.45 0.10], 'tc', [0 0.6 0], 'x', [0.8 0.1 0.1]);

    % (a) oscillator: the fan meets at pi
    f = figure('Color', 'w', 'Position', [60 60 1100 420], 'Visible', 'off');
    tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
    ax = nexttile(tl);  hold(ax, 'on');
    tt = linspace(0, 5, 500);  z = deval(So.sol, tt);
    plot(ax, tt, z(1,:), 'LineWidth', 3, 'Color', col.ext);
    for d = dOsc
        P = cov_perturbed(osc, Eo(1), d, 5);  z1 = deval(P.S.sol, tt);
        plot(ax, tt, z1(1,:), 'LineWidth', 1.6, 'Color', col.nb);
    end
    xline(ax, pi, ':', 't_c = \pi', 'Color', col.tc, 'LineWidth', 2);
    xline(ax, 4, 'k--', 'b');
    plot(ax, pi, 0, 'o', 'MarkerSize', 10, 'LineWidth', 2, 'Color', col.x);
    title(ax, 'Neighbours of the extremal meet it at t = \pi');  xlabel(ax, 't');  ylabel(ax, 'y');
    styl(ax);
    ax = nexttile(tl);  hold(ax, 'on');
    plot(ax, tt, z(3,:), 'LineWidth', 2.5, 'Color', col.ext);  yline(ax, 0, 'k-');
    xline(ax, pi, ':', '\pi', 'Color', col.tc, 'LineWidth', 2);  xline(ax, 4, 'k--', 'b');
    title(ax, 'Jacobi field h(t) = \partial y/\partial p = sin t');  xlabel(ax, 't');
    styl(ax);
    exportgraphics(f, fullfile(figDir, 'oscillator_fan.png'), 'Resolution', 150);  close(f);

    % (b) pendulum: finite delta vs the limit
    f = figure('Color', 'w', 'Position', [60 60 1100 420], 'Visible', 'off');
    tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
    ax = nexttile(tl);  hold(ax, 'on');
    tt = linspace(0, 5, 500);  z = deval(Ep(iw).S.sol, tt);
    plot(ax, tt, z(1,:), 'LineWidth', 3, 'Color', col.ext);
    for d = dPen(1:3)
        P = cov_perturbed(pen, Ep(iw), d, 5);  z1 = deval(P.S.sol, tt);
        plot(ax, tt, z1(1,:), 'LineWidth', 1.6, 'Color', col.nb);
        zc = deval(Ep(iw).S.sol, P.tCross(1));
        plot(ax, P.tCross(1), zc(1), 'o', 'MarkerSize', 9, 'LineWidth', 2, 'Color', col.x);
    end
    xline(ax, tcP, ':', 't_c', 'Color', col.tc, 'LineWidth', 2);
    title(ax, 'Pendulum: finite \delta crosses near, not at, t_c');  xlabel(ax, 't');  ylabel(ax, 'y');
    styl(ax);
    ax = nexttile(tl);
    loglog(ax, dPen, errP, 'o-', 'LineWidth', 2, 'Color', col.x);  hold(ax, 'on');
    loglog(ax, dPen, errP(end)*dPen/dPen(end), 'k--');
    title(ax, '|crossing - t_c| vs \delta  (dashed: slope 1)');  xlabel(ax, '\delta');
    styl(ax);
    exportgraphics(f, fullfile(figDir, 'pendulum_limit.png'), 'Resolution', 150);  close(f);

    % (c) Morse: eigenvalues vs b
    f = figure('Color', 'w', 'Position', [60 60 700 420], 'Visible', 'off');
    ax = axes(f);  hold(ax, 'on');
    plot(ax, bGrid, lamTrue.', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    plot(ax, bGrid, lam.', 'o', 'MarkerSize', 4, 'LineWidth', 1.2);
    yline(ax, 0, 'k-');
    for k = 1:2, xline(ax, k*pi, ':', sprintf('%d\\pi', k), 'Color', col.tc, 'LineWidth', 2); end
    ylim(ax, [-3 3]);
    title(ax, 'Second-variation eigenvalues \lambda_1, \lambda_2, \lambda_3 vs b (oscillator)');
    xlabel(ax, 'b');  ylabel(ax, '\lambda');
    styl(ax);
    exportgraphics(f, fullfile(figDir, 'morse_eigenvalues.png'), 'Resolution', 150);  close(f);

    % (d) catenary: two extremals, their shooting function, and Delta J
    f = figure('Color', 'w', 'Position', [60 60 1300 400], 'Visible', 'off');
    tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact');
    ax = nexttile(tl);  hold(ax, 'on');
    tt = linspace(-0.5, 0.5, 300);
    zs = deval(Ec(1).S.sol, tt);  zd = deval(Ec(2).S.sol, tt);
    plot(ax, tt, zs(1,:), 'LineWidth', 3, 'Color', col.ext);
    plot(ax, tt, zd(1,:), 'LineWidth', 3, 'Color', col.x);
    xline(ax, Sd.tConj(1), ':', 't_c (deep)', 'Color', col.tc, 'LineWidth', 2);
    text(ax, 0, 0.80, 'shallow: minimiser', 'Color', col.ext, 'FontSize', 13, ...
         'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(ax, -0.02, 0.17, 'deep: not', 'Color', col.x, 'FontSize', 13, ...
         'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    ylim(ax, [0.1 1.05]);
    title(ax, 'Two catenaries, same ends');  xlabel(ax, 't');  ylabel(ax, 'y');
    styl(ax);
    ax = nexttile(tl);  hold(ax, 'on');
    plot(ax, scanC.p, scanC.r, 'k', 'LineWidth', 2);  yline(ax, 0, 'k-');
    plot(ax, [Ec.p], [0 0], 'o', 'MarkerSize', 9, 'LineWidth', 2, 'Color', col.x);
    ylim(ax, [-1 3]);
    title(ax, 'Shooting function r(p)');  xlabel(ax, 'p = y''(a)');
    styl(ax);
    ax = nexttile(tl);  hold(ax, 'on');
    e1 = linspace(-0.15, 0.15, 21);
    J0s = cov_functional(cat_, Ec(1).S.sol, zeros(size(Vs.eta)));
    J0d = cov_functional(cat_, Ec(2).S.sol, zeros(size(Vd.eta)));
    plot(ax, e1, arrayfun(@(e) cov_functional(cat_, Ec(1).S.sol, e*Vs.eta) - J0s, e1), ...
         'o-', 'LineWidth', 2, 'Color', col.ext);
    plot(ax, e1, arrayfun(@(e) cov_functional(cat_, Ec(2).S.sol, e*Vd.eta) - J0d, e1), ...
         'o-', 'LineWidth', 2, 'Color', col.x);
    yline(ax, 0, 'k-');
    title(ax, '\DeltaJ along \eta^*  (blue: shallow, red: deep)');  xlabel(ax, '\epsilon');
    styl(ax);
    exportgraphics(f, fullfile(figDir, 'catenary.png'), 'Resolution', 150);  close(f);

    % (e) the explorer itself (uifigure, invisible)
    app = conjugate_point_explorer('Visible', 'off', 'Preset', 5);
    app.selectExtremal(2);  drawnow;
    exportapp(app.fig, fullfile(figDir, 'explorer_pendulum.png'));
    delete(app.fig);
    fprintf('figures written to %s\n', figDir);
end

% ---------------------------------------------------------------------------
function t = firstOr(v)
% FIRSTOR  First element, or NaN when empty. INPUTS: v. OUTPUTS: t scalar.
if isempty(v), t = NaN; else, t = v(1); end
end

function styl(ax)
% STYL  Presentation-weight axes. INPUTS: ax axes handle.
ax.FontSize = 13;  ax.FontWeight = 'bold';  ax.LineWidth = 1.2;
ax.Box = 'on';  grid(ax, 'on');
end

function s = pf(ok)
% PF  PASS/FAIL text. INPUTS: ok logical. OUTPUTS: s char.
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
