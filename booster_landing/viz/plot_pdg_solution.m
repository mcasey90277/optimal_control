function fig = plot_pdg_solution(solC, solV, outfile)
% PLOT_PDG_SOLUTION  2x2 comparison figure: collocation vs convex PDG.
%
% Panels: (a) 3D trajectory + glideslope cone + thrust arrows every 5th
% node, both solvers overlaid; (b) throttle ||T||/Tmax vs t (the MIN-MAX
% "money plot"), both solvers, WITH the two guidance bounds
% (Tmin/Tmax and the etaT de-rated ceiling) labeled; (c) mass vs t;
% (d) speed vs t (terminal target is P.vf, NOT zero -- see booster_params).
%
% The collocation throttle curve is sampled via lib/hs_quad_ctrl.m (the
% flyable per-segment reconstruction Hermite-Simpson's own defects are
% built against), NOT a global pchip spline through the nodes -- a pchip
% smears the bang-bang switch into a ramp (see hs_quad_ctrl.m's own
% header for the measured consequence of getting this wrong). The convex
% solver has no midpoint samples (no .Um field) and is plotted directly
% off its own dense node grid.
%
% "MIN-MAX", corrected 2026-08-09: this header called panel (b) the
% "max-min-max money plot" from the first draft onward. The certified
% structure for these boundary conditions is min-max -- a single
% Tmin->Tmax switch held to touchdown -- because P.Tmin already exceeds
% vehicle weight at every mass on this trajectory, so there is no
% hover-capable arc to open with. The finding is significant enough that
% the theory note gives it a section and G5's structure check was written
% around it; only this comment was left behind.
%
% INPUTS:
%   solC    - solve_pdg_colloc solution (.t .X .U .Um .tf .mf .P)
%   solV    - solve_pdg_convex solution (.t .X .U .tf .mf .P), same tf
%   outfile - (optional) PNG path; exportgraphics at 200 dpi if given
% OUTPUTS:
%   fig - figure handle (Visible off). OWNERSHIP CONTRACT (fix report,
%         2026-08-09): called with an output requested (fig = ...), the
%         CALLER owns the handle and is responsible for closing it.
%         Called as a bare statement (no output captured, nargout==0 --
%         e.g. a sweep or a smoke test that only wants the exported PNG),
%         this function closes the figure itself before returning, so
%         repeated calls in one long -batch process (Task 10's flagship,
%         Task 11's sweeps) do not accumulate invisible figure handles.
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-08-08-booster-landing-design.md
%   [2] lib/hs_quad_ctrl.m -- flyable control reconstruction + switch-step
if nargin < 3, outfile = ''; end
P = solC.P;

colC = [0 0.447 0.741];      % colorblind-safe blue  -- collocation
colV = [0.85 0.325 0.098];   % colorblind-safe orange -- convex

fig = figure('Visible', 'off', 'Position', [100 100 1400 1100], 'Color', 'w');
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% (a) 3D trajectory + glideslope cone + thrust arrows:
ax1 = nexttile(tl);
hold(ax1, 'on');  grid(ax1, 'on');  box(ax1, 'on');
cotg  = 1 / tand(P.gs_deg);
zmax  = max([solC.X(3,:), solV.X(3,:)]) * 1.05;
zcone = linspace(0, zmax, 20);
th    = linspace(0, 2*pi, 40);
[Zc, Th] = meshgrid(zcone, th);
Rc = cotg * Zc;
Xc = Rc .* cos(Th);  Yc = Rc .* sin(Th);
surf(ax1, Xc, Yc, Zc, 'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 0.12, ...
     'EdgeColor', 'none', 'HandleVisibility', 'off');

plot3(ax1, solC.X(1,:), solC.X(2,:), solC.X(3,:), '-', 'Color', colC, ...
      'LineWidth', 2, 'DisplayName', 'collocation');
plot3(ax1, solV.X(1,:), solV.X(2,:), solV.X(3,:), '--', 'Color', colV, ...
      'LineWidth', 2, 'DisplayName', 'convex');

quiverScale = 0.25 * zmax / P.Tmax;
plot_thrust_arrows(ax1, solC.X, solC.U, quiverScale, colC);
plot_thrust_arrows(ax1, solV.X, solV.U, quiverScale, colV);

plot3(ax1, 0, 0, 0, 'kp', 'MarkerFaceColor', 'k', 'MarkerSize', 12, ...
      'HandleVisibility', 'off');
xlabel(ax1, 'x [m]');  ylabel(ax1, 'y [m]');  zlabel(ax1, 'z [m]');
xmax = max(abs([solC.X(1:2,:), solV.X(1:2,:)]), [], 'all') * 1.1;
xmax = max(xmax, 50);
xlim(ax1, [-xmax xmax]);  ylim(ax1, [-xmax xmax]);  zlim(ax1, [0 zmax]);
view(ax1, -35, 20);  axis(ax1, 'vis3d');
title(ax1, 'trajectory + glideslope cone + thrust');
legend(ax1, 'Location', 'northeast');

%% (b) throttle ||T||/Tmax -- the money plot:
ax2 = nexttile(tl);
hold(ax2, 'on');  grid(ax2, 'on');  box(ax2, 'on');
tqC = linspace(0, solC.tf, 400);
Nseg = numel(solC.t) - 1;  h = solC.tf / Nseg;
TmagC = zeros(size(tqC));
for k = 1:numel(tqC)
    Tv = hs_quad_ctrl(tqC(k), solC.U, solC.Um, h, Nseg, P.Tmin, P.etaT*P.Tmax);
    TmagC(k) = sqrt(sum(Tv.^2));
end
plot(ax2, tqC, TmagC/P.Tmax, '-', 'Color', colC, 'LineWidth', 2, ...
     'DisplayName', 'collocation');
plot(ax2, solV.t, sqrt(sum(solV.U.^2,1))/P.Tmax, '--', 'Color', colV, ...
     'LineWidth', 2, 'DisplayName', 'convex');
yline(ax2, P.Tmin/P.Tmax, 'k:', 'LineWidth', 1.25, ...
      'Label', sprintf('Tmin/Tmax = %.2f', P.Tmin/P.Tmax), ...
      'LabelHorizontalAlignment', 'left');
yline(ax2, P.etaT, 'k--', 'LineWidth', 1.25, ...
      'Label', sprintf('etaT = %.2f (guidance ceiling)', P.etaT), ...
      'LabelHorizontalAlignment', 'left');
xlabel(ax2, 'time [s]');  ylabel(ax2, '||T|| / Tmax [-]');
xlim(ax2, [0 max(solC.tf, solV.tf)]);  ylim(ax2, [0 1.05]);
title(ax2, 'throttle (bang-bang, both solvers)');
legend(ax2, 'Location', 'east');

%% (c) mass vs time:
ax3 = nexttile(tl);
hold(ax3, 'on');  grid(ax3, 'on');  box(ax3, 'on');
plot(ax3, solC.t, solC.X(7,:), '-', 'Color', colC, 'LineWidth', 2, ...
     'DisplayName', 'collocation');
plot(ax3, solV.t, solV.X(7,:), '--', 'Color', colV, 'LineWidth', 2, ...
     'DisplayName', 'convex');
yline(ax3, P.mdry, 'k:', 'LineWidth', 1.25, 'Label', 'mdry');
xlabel(ax3, 'time [s]');  ylabel(ax3, 'mass [kg]');
xlim(ax3, [0 max(solC.tf, solV.tf)]);
ylim(ax3, [P.mdry*0.98, P.m0*1.02]);
title(ax3, 'mass depletion');
legend(ax3, 'Location', 'northeast');

%% (d) speed vs time:
ax4 = nexttile(tl);
hold(ax4, 'on');  grid(ax4, 'on');  box(ax4, 'on');
speedC = sqrt(sum(solC.X(4:6,:).^2, 1));
speedV = sqrt(sum(solV.X(4:6,:).^2, 1));
plot(ax4, solC.t, speedC, '-', 'Color', colC, 'LineWidth', 2, ...
     'DisplayName', 'collocation');
plot(ax4, solV.t, speedV, '--', 'Color', colV, 'LineWidth', 2, ...
     'DisplayName', 'convex');
vfMag = sqrt(sum(P.vf.^2));
yline(ax4, vfMag, 'k:', 'LineWidth', 1.25, ...
      'Label', sprintf('|v(tf)| target = %.1f m/s', vfMag));
xlabel(ax4, 'time [s]');  ylabel(ax4, 'speed [m/s]');
xlim(ax4, [0 max(solC.tf, solV.tf)]);
ylim(ax4, [0, max([speedC, speedV])*1.05]);
title(ax4, 'speed profile');
legend(ax4, 'Location', 'northeast');

title(tl, sprintf('PDG guidance: collocation vs convex   (tf=%6.2f s, mf_C=%7.1f kg, mf_V=%7.1f kg)', ...
      solC.tf, solC.mf, solV.mf));

if ~isempty(outfile)
    exportgraphics(fig, outfile, 'Resolution', 200);
end
if nargout == 0
    close(fig);        % no caller owns the handle -- close it here (see
                        % OUTPUTS' ownership-contract note above)
end
end

function plot_thrust_arrows(ax, X, U, scale, col)
% Quiver3 the thrust vector at every 5th node, scaled for visibility.
idx = 1:5:size(X,2);
Xp = X(1,idx);  Yp = X(2,idx);  Zp = X(3,idx);
Up = U(1,idx)*scale;  Vp = U(2,idx)*scale;  Wp = U(3,idx)*scale;
quiver3(ax, Xp, Yp, Zp, Up, Vp, Wp, 0, 'Color', col, 'LineWidth', 1, ...
        'MaxHeadSize', 0.4, 'HandleVisibility', 'off');
end
