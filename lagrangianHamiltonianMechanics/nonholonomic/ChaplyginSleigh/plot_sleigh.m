function files = plot_sleigh(t, Y, out, p, info, outDir)
%% Purpose:
%
%   The static diagnostics: a path figure (center-of-mass and blade
%   traces, heading ticks every 0.5 s, start and end markers, equal
%   aspect) beside a zoom, and a 2 x 3 panel (v(t) with the exact solution, the
%   asymptote and lambda(t) on a right axis; w(t) with the exact
%   solution; th(t); energy drift; the energy split into rotational and
%   translational parts with their (flat) total, the constraint residual
%   in the title; the (v, w) phase portrait on its energy ellipse). The
%   path figure is two panels: the whole run and a zoom on the transient
%   (bodies and heading every 0.5 s). Two PNGs; headless by default.
%
%% Inputs:
%
%  t                        [1 x N]                 Times (s)
%
%  Y                        [5 x N]                 Reduced solution
%
%  out                      struct                  From sleigh_quantities
%
%  p                        struct                  From sleigh_params
%
%  info                     struct                  .label; .visible
%                                                   [false]; .resid [1 x
%                                                   N] constraint residual
%                                                   (optional; computed
%                                                   from out.Yfull when
%                                                   absent)
%
%  outDir                   char                    Folder for the PNGs
%
%% Outputs:
%
%  files                    cell [1 x 2]            The PNG paths
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the defaults, three seconds, into tempdir
     p = sleigh_params();  t = 0:0.01:3;
     [~, yy] = ode113(@(a, b) sleigh_rhs(a, b, p), t, [0; 0; 0; 1; 2], odeset('RelTol', 1e-12, 'AbsTol', 1e-12));
     out = sleigh_quantities(t, yy.', p);
     files = plot_sleigh(t, yy.', out, p, struct('label', 'demo', 'visible', true), tempdir);
     fprintf('wrote %s and %s\n', files{1}, files{2});
     return
end
if ~isfield(info, 'visible') || isempty(info.visible), info.visible = false; end
if info.visible, vis = 'on'; else, vis = 'off'; end
th = Y(3, :);  v = Y(4, :);  w = Y(5, :);
if isfield(info, 'resid') && ~isempty(info.resid)
    resid = info.resid;
else
    vc = out.Yfull(4:5, :);
    resid = abs(-vc(1, :).*sin(th) + vc(2, :).*cos(th) - p.a*w);
end
E0 = out.E(1);
drift = (out.E - E0)/max(abs(E0), eps);

% 1. the path: the whole run on the left, a zoom on the transient on the right
tEnd = t(end);
if p.a ~= 0 && w(1) ~= 0 && out.C > 0
    tZoom = min(tEnd, 4/(abs(p.k)*sqrt(out.C)));
else
    tZoom = tEnd;
end
kz = t <= tZoom;
Lb = max(abs(p.a), 0.1);
f1 = figure('Visible', vis, 'Color', 'w', 'Position', [100 100 1400 700]);
tl1 = tiledlayout(f1, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile(tl1);  hold(ax, 'on');
plot(ax, out.rc(1, :), out.rc(2, :), '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.8);
plot(ax, out.rb(1, :), out.rb(2, :), '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.2);
tick = 0:0.5:t(end);
kt = arrayfun(@(tt) find(t >= tt, 1), tick);
quiver(ax, out.rb(1, kt), out.rb(2, kt), Lb*cos(th(kt)), Lb*sin(th(kt)), 0, 'Color', [0.4 0.4 0.4], 'MaxHeadSize', 0.8);
plot(ax, out.rb(1, 1), out.rb(2, 1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(ax, out.rb(1, end), out.rb(2, end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
axis(ax, 'equal');  grid(ax, 'on');
xlabel(ax, 'x  [m]');  ylabel(ax, 'y  [m]');
title(ax, sprintf('%s: a = %+.3g m, v0 = %.3g m/s, w0 = %.3g rad/s (ticks every 0.5 s)', info.label, p.a, v(1), w(1)), 'Interpreter', 'none');
legend(ax, {'center of mass', 'blade', 'heading', 'start', 'end'}, 'Location', 'best');

ax = nexttile(tl1);  hold(ax, 'on');
plot(ax, out.rc(1, kz), out.rc(2, kz), '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.8);
plot(ax, out.rb(1, kz), out.rb(2, kz), '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.2);
xyz  = [out.rc(:, kz), out.rb(:, kz)];
span = max(max(xyz, [], 2) - min(xyz, [], 2));
span = max(span, Lb);
kz5  = kt(tick <= tZoom);
quiver(ax, out.rb(1, kz5), out.rb(2, kz5), 0.06*span*cos(th(kz5)), 0.06*span*sin(th(kz5)), 0, 'Color', [0.4 0.4 0.4], 'MaxHeadSize', 0.8);
outl = sleigh_body(p).outline;  outl = [outl, outl(:, 1)];
for kk = kz5
    Rk = [cos(th(kk)), -sin(th(kk)); sin(th(kk)), cos(th(kk))];
    Ok = out.rb(:, kk) + Rk*outl;
    plot(ax, Ok(1, :), Ok(2, :), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
end
plot(ax, out.rb(1, 1), out.rb(2, 1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
axis(ax, 'equal');  grid(ax, 'on');
xlabel(ax, 'x  [m]');  ylabel(ax, 'y  [m]');
title(ax, sprintf('zoom: 0 <= t <= %.3g s (bodies and heading every 0.5 s)', tZoom));
files{1} = fullfile(outDir, [info.label '_path.png']);
exportgraphics(f1, files{1}, 'Resolution', 150);
close(f1);

% 2. the diagnostics
f2 = figure('Visible', vis, 'Color', 'w', 'Position', [100 100 1400 800]);
tl = tiledlayout(f2, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile(tl);  hold(ax, 'on');
yyaxis(ax, 'left');
plot(ax, t, v, 'b-', 'LineWidth', 1.6);  plot(ax, t, out.vExact, 'k--');
yline(ax, out.asymptote, ':', 'Color', [0.3 0.3 0.3]);
ylabel(ax, 'v  [m/s]');
yyaxis(ax, 'right');  plot(ax, t, out.lambda, '-', 'Color', [0.9 0.4 0.9]);  ylabel(ax, '\lambda  [N]');
ax.YAxis(2).Color = [0.9 0.4 0.9];
xlabel(ax, 't  [s]');  grid(ax, 'on');  title(ax, 'forward speed (exact dashed, asymptote dotted), blade force');
ax = nexttile(tl);  hold(ax, 'on');
plot(ax, t, w, 'b-', 'LineWidth', 1.6);  plot(ax, t, out.wExact, 'k--');
xlabel(ax, 't  [s]');  ylabel(ax, '\omega  [rad/s]');  grid(ax, 'on');  title(ax, 'yaw rate (exact dashed)');
ax = nexttile(tl);
plot(ax, t, th, 'b-', 'LineWidth', 1.6);
xlabel(ax, 't  [s]');  ylabel(ax, '\theta  [rad]');  grid(ax, 'on');  title(ax, 'heading');
ax = nexttile(tl);
plot(ax, t, drift, 'b-');
xlabel(ax, 't  [s]');  ylabel(ax, '(E - E_0)/E_0');  grid(ax, 'on');  title(ax, sprintf('energy drift, max %.1e', max(abs(drift))));
ax = nexttile(tl);  hold(ax, 'on');
vc = out.Yfull(4:5, :);
E0s = max(abs(E0), eps);
rot = 0.5*p.J*w.^2/E0s;
tra = 0.5*p.m*(vc(1, :).^2 + vc(2, :).^2)/E0s;
plot(ax, t, rot, '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.4);
plot(ax, t, tra, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.4);
plot(ax, t, rot + tra, 'k-', 'LineWidth', 1.2);
xlabel(ax, 't  [s]');  ylabel(ax, 'E_i / E_0');  grid(ax, 'on');
legend(ax, {'rotational', 'translational', 'total'}, 'Location', 'best');
title(ax, sprintf('energy split (max |A qd| = %.1e)', max(resid)));
ax = nexttile(tl);  hold(ax, 'on');
sC = sqrt(max(out.C, 0));  s = linspace(0, 2*pi, 361);
plot(ax, sC*cos(s), sC*sqrt(p.m/p.Jb)*sin(s), '-', 'Color', [0.6 0.6 0.6]);
plot(ax, v, w, 'b-', 'LineWidth', 1.6);
plot(ax, v(1), w(1), 'go', 'MarkerFaceColor', 'g');  plot(ax, v(end), w(end), 'rs', 'MarkerFaceColor', 'r');
axis(ax, 'equal');  grid(ax, 'on');
xlabel(ax, 'v  [m/s]');  ylabel(ax, '\omega  [rad/s]');  title(ax, 'phase portrait on the energy ellipse');
title(tl, sprintf('%s: m = %.3g kg, J = %.3g kg m^2, a = %+.3g m', info.label, p.m, p.J, p.a), 'Interpreter', 'none');
files{2} = fullfile(outDir, [info.label '_diag.png']);
exportgraphics(f2, files{2}, 'Resolution', 150);
close(f2);
end
