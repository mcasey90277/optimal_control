function fig = plot_vacuum_vs_drag(solC, solD, outfile)
% PLOT_VACUUM_VS_DRAG  1x3 comparison figure: vacuum vs drag-on PDG.
%
% Panels: (a) trajectory downrange-altitude (planar range sqrt(x^2+y^2)
% vs z, both solutions overlaid); (b) throttle ||T||/Tmax vs t, both
% solutions sampled via lib/hs_quad_ctrl.m (the flyable per-segment
% reconstruction, same representation certify_pdg's G2 gate certifies --
% see that function's own header for why a global pchip is wrong here);
% (c) mass vs t. An annotation box reports the headline Phase-2 numbers:
% Delta-fuel (vacuum minus drag, i.e. what drag SAVES) and Delta-tf
% (drag minus vacuum) -- "the what drag-free guidance misses" number.
%
% Both solC and solD are solve_pdg_colloc solutions (each carries its own
% .P: solC.P.drag.on=false, solD.P.drag.on=true, identical otherwise --
% see run_booster_landing.m's Phase-2 branch, Pd = P with only .drag.on
% flipped), so Tmin/TmaxG are read off solD.P (equal to solC.P's) rather
% than assumed equal by construction.
%
% INPUTS:
%   solC    - solve_pdg_colloc solution, P.drag.on=false (vacuum)
%   solD    - solve_pdg_colloc solution, P.drag.on=true  (drag), same P
%             otherwise (warm-started from solC -- see run_booster_landing)
%   outfile - (optional) PNG path; exportgraphics at 200 dpi if given
% OUTPUTS:
%   fig - figure handle (Visible off). OWNERSHIP CONTRACT (same as
%         plot_pdg_solution.m / plot_footprint.m, task-9/10 fix reports):
%         with an output requested (fig = ...), the CALLER owns and must
%         close it; called as a bare statement (nargout==0), this
%         function closes the figure itself.
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-08-08-booster-landing-design.md
%   [2] lib/hs_quad_ctrl.m -- flyable control reconstruction + switch-step
if nargin < 3, outfile = ''; end
P = solD.P;   % Tmin/Tmax/etaT identical between solC.P and solD.P by
              % construction (Pd = P; Pd.drag.on = true) -- see header

colV = [0 0.447 0.741];      % colorblind-safe blue   -- vacuum
colD = [0.85 0.325 0.098];   % colorblind-safe orange -- drag

fig = figure('Visible', 'off', 'Position', [100 100 1500 520], 'Color', 'w');
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

%% (a) trajectory: downrange (planar range) vs altitude:
ax1 = nexttile(tl);
hold(ax1, 'on');  grid(ax1, 'on');  box(ax1, 'on');
rangeC = sqrt(solC.X(1,:).^2 + solC.X(2,:).^2);
rangeD = sqrt(solD.X(1,:).^2 + solD.X(2,:).^2);
plot(ax1, rangeC, solC.X(3,:), '-', 'Color', colV, 'LineWidth', 2, ...
     'DisplayName', 'vacuum');
plot(ax1, rangeD, solD.X(3,:), '--', 'Color', colD, 'LineWidth', 2, ...
     'DisplayName', 'drag');
plot(ax1, 0, 0, 'kp', 'MarkerFaceColor', 'k', 'MarkerSize', 12, ...
     'HandleVisibility', 'off');
xlabel(ax1, 'downrange [m]');  ylabel(ax1, 'altitude z [m]');
title(ax1, 'trajectory: downrange-altitude');
legend(ax1, 'Location', 'northeast');

%% (b) throttle ||T||/Tmax -- both solutions via hs_quad_ctrl:
ax2 = nexttile(tl);
hold(ax2, 'on');  grid(ax2, 'on');  box(ax2, 'on');
TmaxG = P.etaT * P.Tmax;
tqC = linspace(0, solC.tf, 400);
NsegC = numel(solC.t) - 1;  hC = solC.tf / NsegC;
TmagC = zeros(size(tqC));
for k = 1:numel(tqC)
    Tv = hs_quad_ctrl(tqC(k), solC.U, solC.Um, hC, NsegC, P.Tmin, TmaxG);
    TmagC(k) = sqrt(sum(Tv.^2));
end
tqD = linspace(0, solD.tf, 400);
NsegD = numel(solD.t) - 1;  hD = solD.tf / NsegD;
TmagD = zeros(size(tqD));
for k = 1:numel(tqD)
    Tv = hs_quad_ctrl(tqD(k), solD.U, solD.Um, hD, NsegD, P.Tmin, TmaxG);
    TmagD(k) = sqrt(sum(Tv.^2));
end
plot(ax2, tqC, TmagC/P.Tmax, '-', 'Color', colV, 'LineWidth', 2, ...
     'DisplayName', 'vacuum');
plot(ax2, tqD, TmagD/P.Tmax, '--', 'Color', colD, 'LineWidth', 2, ...
     'DisplayName', 'drag');
yline(ax2, P.Tmin/P.Tmax, 'k:', 'LineWidth', 1.25, ...
      'Label', sprintf('Tmin/Tmax = %.2f', P.Tmin/P.Tmax), ...
      'LabelHorizontalAlignment', 'left');
yline(ax2, P.etaT, 'k--', 'LineWidth', 1.25, ...
      'Label', sprintf('etaT = %.2f (guidance ceiling)', P.etaT), ...
      'LabelHorizontalAlignment', 'left');
xlabel(ax2, 'time [s]');  ylabel(ax2, '||T|| / Tmax [-]');
xlim(ax2, [0 max(solC.tf, solD.tf)]);  ylim(ax2, [0 1.05]);
title(ax2, 'throttle');
legend(ax2, 'Location', 'east');

%% (c) mass vs time:
ax3 = nexttile(tl);
hold(ax3, 'on');  grid(ax3, 'on');  box(ax3, 'on');
plot(ax3, solC.t, solC.X(7,:), '-', 'Color', colV, 'LineWidth', 2, ...
     'DisplayName', 'vacuum');
plot(ax3, solD.t, solD.X(7,:), '--', 'Color', colD, 'LineWidth', 2, ...
     'DisplayName', 'drag');
yline(ax3, P.mdry, 'k:', 'LineWidth', 1.25, 'Label', 'mdry');
xlabel(ax3, 'time [s]');  ylabel(ax3, 'mass [kg]');
xlim(ax3, [0 max(solC.tf, solD.tf)]);
ylim(ax3, [P.mdry*0.98, P.m0*1.02]);
title(ax3, 'mass depletion');
legend(ax3, 'Location', 'northeast');

%% Headline numbers -- "what drag-free guidance misses":
fuelV  = P.m0 - solC.mf;   fuelD  = P.m0 - solD.mf;
dfuel  = fuelV - fuelD;                 % >0 means drag SAVES fuel
dtf    = solD.tf - solC.tf;             % drag minus vacuum
txt = sprintf(['fuel: vac %.1f kg  drag %.1f kg\n' ...
               'dfuel (vac-drag) = %+.1f kg\n' ...
               'tf:   vac %.2f s   drag %.2f s\n' ...
               'dtf (drag-vac)   = %+.2f s'], ...
    fuelV, fuelD, dfuel, solC.tf, solD.tf, dtf);
annotation(fig, 'textbox', [0.005 0.72 0.20 0.24], 'String', txt, ...
    'FontName', 'FixedWidth', 'FontSize', 9, 'EdgeColor', [0.7 0.7 0.7], ...
    'BackgroundColor', [1 1 1], 'FitBoxToText', 'on');

title(tl, sprintf('Phase 2: vacuum vs drag   (dfuel=%+.1f kg, dtf=%+.2f s)', ...
      dfuel, dtf));

if ~isempty(outfile)
    exportgraphics(fig, outfile, 'Resolution', 200);
end
if nargout == 0
    close(fig);        % no caller owns the handle -- close it here (see
                        % OUTPUTS' ownership-contract note above)
end
end
