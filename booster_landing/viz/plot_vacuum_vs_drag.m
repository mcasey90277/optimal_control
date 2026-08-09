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
% flipped), so Tmin/TmaxG are read off solD.P; the equality (up to
% .drag.on) with solC.P is CHECKED below (task-11 close-out review: an
% earlier version of this file claimed this in prose without asserting
% it), not merely assumed by construction.
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
% Bounds-equality ASSERT (task-11 close-out review): make the header's
% "identical otherwise" claim checked, not just asserted in prose. Only
% .drag itself is allowed to differ between the two params structs.
assert(isequal(P.Tmin, solC.P.Tmin) && isequal(P.Tmax, solC.P.Tmax) && ...
       isequal(P.etaT, solC.P.etaT), ...
    'plot_vacuum_vs_drag:paramMismatch', ...
    ['solD.P and solC.P disagree on Tmin/Tmax/etaT -- this function reads ' ...
     'Tmin/TmaxG off solD.P only and assumes solC.P matches (see header).']);

colV = [0 0.447 0.741];      % colorblind-safe blue   -- vacuum
colD = [0.85 0.325 0.098];   % colorblind-safe orange -- drag

% Figure is taller than the original 520 px, and the layout is 2 ROWS not
% 1 (task-11 close-out review, Important-3a): the annotation box used to
% float at [0.005 0.72 ...] in hand-picked figure-normalized coords, which
% occluded panel (a)'s title and the top of its trajectory line (that
% panel's data starts at HIGH altitude/large downrange -- exactly the
% top-left corner the box sat in). Two things this went through before
% landing here (both caught by re-rendering and inspecting the actual PNG,
% not just re-running the tests -- a passing assert only proves the file
% exists, not that it looks right): (1) `tl.RowHeight` (the natural way to
% make row 2 short) does not exist on this MATLAB's TiledChartLayout --
% threw "Unrecognized property" outright; (2) a first Position-based
% attempt reserved too little margin above the tiles for the layout's OWN
% auto-placed title, which then got clipped by exportgraphics' tight crop
% together with the top pixel row of panel (a). Fix: 2 rows via GridSize
% at construction (row 2 stays whatever tiledlayout's default equal split
% gives it -- see the note by axFoot below for why that's fine here), and
% tl.Position left at its AUTOMATIC default (not hand-set) so the layout
% keeps doing its own title-space bookkeeping.
fig = figure('Visible', 'off', 'Position', [100 100 1500 680], 'Color', 'w');
tl  = tiledlayout(fig, 8, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

%% (a) trajectory: downrange (planar range) vs altitude:
ax1 = nexttile(tl, 1, [7 1]);
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
ax2 = nexttile(tl, 2, [7 1]);
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
      'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
yline(ax2, P.etaT, 'k--', 'LineWidth', 1.25, ...
      'Label', sprintf('etaT = %.2f (guidance ceiling)', P.etaT), ...
      'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xlabel(ax2, 'time [s]');  ylabel(ax2, '||T|| / Tmax [-]');
xlim(ax2, [0 max(solC.tf, solD.tf)]);  ylim(ax2, [0 1.05]);
title(ax2, 'throttle');
legend(ax2, 'Location', 'east');

%% (c) mass vs time:
ax3 = nexttile(tl, 3, [7 1]);
hold(ax3, 'on');  grid(ax3, 'on');  box(ax3, 'on');
plot(ax3, solC.t, solC.X(7,:), '-', 'Color', colV, 'LineWidth', 2, ...
     'DisplayName', 'vacuum');
plot(ax3, solD.t, solD.X(7,:), '--', 'Color', colD, 'LineWidth', 2, ...
     'DisplayName', 'drag');
yline(ax3, P.mdry, 'k:', 'LineWidth', 1.25, 'Label', 'mdry', ...
      'HandleVisibility', 'off');
xlabel(ax3, 'time [s]');  ylabel(ax3, 'mass [kg]');
xlim(ax3, [0 max(solC.tf, solD.tf)]);
ylim(ax3, [P.mdry*0.98, P.m0*1.02]);
title(ax3, 'mass depletion');
legend(ax3, 'Location', 'northeast');

%% Headline numbers -- "what drag-free guidance misses":
fuelV  = P.m0 - solC.mf;   fuelD  = P.m0 - solD.mf;
dfuel  = fuelV - fuelD;                 % >0 means drag SAVES fuel
dtf    = solD.tf - solC.tf;             % drag minus vacuum
txt = sprintf(['fuel: vac %.1f kg  drag %.1f kg   dfuel (vac-drag) = %+.1f kg   |   ' ...
               'tf: vac %.2f s  drag %.2f s   dtf (drag-vac) = %+.2f s'], ...
    fuelV, fuelD, dfuel, solC.tf, solD.tf, dtf);
% Footer TILE (task-11 close-out review, Important-3a), not a
% hand-positioned `annotation` box -- row 8 of an 8-row GridSize (the 3
% panels above each span rows 1-7 via nexttile(tl, col, [7 1]), so this is
% a genuine 1/8-height strip, not a manually-guessed fraction). A real
% tile is guaranteed by tiledlayout's own bookkeeping to sit in empty
% figure real estate below the panels, with the layout still handling its
% own title placement automatically above row 1 -- so this never overlaps
% a panel's title or data the way the original figure-normalized-
% coordinate box did.
axFoot = nexttile(tl, 22, [1 3]);   % tile 22 = row 8, col 1 (row-major, 3 cols/row)
axis(axFoot, 'off');
text(axFoot, 0.5, 0.5, txt, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontName', 'FixedWidth', 'FontSize', 9);

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
