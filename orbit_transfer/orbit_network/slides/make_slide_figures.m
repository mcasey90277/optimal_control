%% MAKE_SLIDE_FIGURES  Presentation-grade figures for the orbit-network slide deck.
%
%   Runs the two study scripts headless (their outputs go to a scratch
%   folder, never to results/), then redraws their figures from the
%   workspaces they leave behind at slide size and slide fonts, and adds
%   one figure the note does not have: the 18 node orbits in the rotating
%   frame. Every number on a slide comes from the same variables the
%   note's figures came from; nothing is recomputed here.
%
%   Output: slides/figures/*.png (committed with the deck).
%
%   Sections
%     0  settings
%     1  run the studies headless
%     2  the 18 orbits (new)
%     3  the graph at 1 N
%     4  the accessibility plane (raw coverage)
%     5  catalog closure down the ladder
%     6  relay closure
%     7  stability versus cost
%     8  census
%
%  M. Casey                                                   (c) 11/08/2026
%  Copyright Coorbital Inc.

%% 0  settings
clearvars; clc
here    = fileparts(mfilename('fullpath'));
root    = fileparts(here);                                  % orbit_network/
scratch = fullfile(tempdir, 'orbit_network_slides_run');    % the studies' own outputs (never results/)
if ~exist(scratch, 'dir'), mkdir(scratch); end
addpath(root, fullfile(root, 'lib'), fullfile(fileparts(root), 'costate_common'));

%% 1  run the studies headless; keep their workspaces
%   Each study begins with `clearvars -except study_opts`, so the first
%   study's workspace is carried through the second inside study_opts.
study_opts = struct('visible', 'off', 'outDir', scratch, 'root', root);   % the studies clear everything else
run(fullfile(study_opts.root, 'stability_cost_study.m'));  % leaves dCJ, sDep, sArr, tf1, surv, stDep, fits, resp, verdict
study_opts.keep = struct('dCJ', dCJ, 'sDep', sDep, 'tf1', tf1, 'surv', surv, 'stDep', stDep, 'fits', {fits}, 'resp', resp, 'verdict', verdict);
run(fullfile(study_opts.root, 'orbit_network_study.m'));   % leaves T, N, E, G, P, caps, rungs, rC, rD, nR, nS, best, R, ...
stabWS = study_opts.keep;
clearvars -except stabWS T N E G P caps rungs rC rD nR nS best R
here    = fileparts(mfilename('fullpath'));
figDir  = fullfile(here, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
FS   = 16;                                                  % axis font size on slides
FSt  = 17;                                                  % title font size
LW   = 2.5;
res  = 200;
axStyle = @(ax) set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', 'GridColor', [0.15 0.15 0.15], 'FontSize', FS, 'LineWidth', 1);
famColor = containers.Map({'dro', 'dpo', 'halo_L1', 'halo_L2', 'tulip'}, ...
                          {[0.00 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19], [0.49 0.18 0.56], [0.93 0.69 0.13]});
famName  = containers.Map({'dro', 'dpo', 'halo_L1', 'halo_L2', 'tulip'}, ...
                          {'DRO (4)', 'DPO (4)', 'L1 halo (2)', 'L2 halo (4)', 'tulip (4)'});
famOf   = @(k) regexprep(k, '_tau.*|_Np.*', '');
rungLab = arrayfun(@(r) sprintf('%g', r), rungs, 'UniformOutput', false);
mu      = T(1).constants.muStar;

%% 2  the 18 orbits in the rotating frame (the nodes), colored by family
fig = figure('Color', 'w', 'Position', [60 60 1280 720], 'Visible', 'off');
ax  = axes(fig);  axStyle(ax);  hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
[sx, sy, sz] = sphere(40);  rMoon = 1737.4/384400;
surf(ax, 1 - mu + rMoon*sx, rMoon*sy, rMoon*sz, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'HandleVisibility', 'off');
hFam = containers.Map();
for kk = 1:numel(N)
    [~, ~, rv] = evalc('get_family_orbit(N(kk).family, N(kk).params)');   % evalc's first output is the loader's log
    fam = famOf(N(kk).key);
    vis = 'off';
    if ~isKey(hFam, fam), hFam(fam) = true;  vis = 'on'; end
    plot3(ax, rv(:, 1), rv(:, 2), rv(:, 3), '-', 'Color', famColor(fam), 'LineWidth', 1.8, 'DisplayName', famName(fam), 'HandleVisibility', vis);
end
Lpts = pumpkyn.cr3bp.lagrangePts(mu);
plot3(ax, Lpts(1:2, 1), Lpts(1:2, 2), Lpts(1:2, 3), 'k^', 'MarkerFaceColor', 'k', 'MarkerSize', 7, 'DisplayName', 'L1, L2');
axis(ax, 'equal');  view(ax, [-35 22]);
xlabel(ax, 'x [ND]');  ylabel(ax, 'y [ND]');  zlabel(ax, 'z [ND]');     % ND = Earth-Moon distance, said on the slide
legend(ax, 'Location', 'northeast', 'FontSize', FS - 1, 'TextColor', 'k', 'Color', 'w');
title(ax, '18 periodic orbits in the Earth-Moon rotating frame (the nodes)', 'FontSize', FSt, 'Color', 'k');
exportgraphics(fig, fullfile(figDir, 'nodes_3d.png'), 'Resolution', res);  close(fig);

%% 3  the graph at 1 N (one panel, big labels)
fig = figure('Color', 'w', 'Position', [60 60 1280 720], 'Visible', 'off');
ax  = axes(fig);  axStyle(ax);
Er  = G.Edges(G.Edges.closed_frac(:, rC) > 0, :);
Gr  = digraph(Er, G.Nodes);
h   = plot(ax, Gr, 'Layout', 'layered', 'Direction', 'right', 'NodeLabel', Gr.Nodes.Name, 'Interpreter', 'none', ...
           'MarkerSize', 9, 'NodeColor', [0.1 0.3 0.7], 'ArrowSize', 10, 'NodeFontSize', FS - 2, 'NodeLabelColor', 'k');
h.LineWidth = 0.6 + 5*Gr.Edges.closed_frac(:, rC);
h.EdgeCData = Gr.Edges.tf_min_days(:, rC);
axStyle(ax);  set(ax, 'XTick', [], 'YTick', []);  xl = xlim(ax);  xlim(ax, [xl(1) - 0.05*(xl(2) - xl(1)), xl(2) + 0.25*(xl(2) - xl(1))]);
colormap(ax, 'parula');  cb = colorbar(ax);  cb.Color = 'k';  cb.FontSize = FS - 2;  cb.Label.String = 'min t_f over accepted phases [days]';  cb.Label.FontSize = FS - 1;
title(ax, sprintf('%g N: %d edges; width = accepted fraction of the 6 x 6 phasing grid, color = fastest accepted time', rungs(rC), numedges(Gr)), 'FontSize', FSt - 1, 'Color', 'k', 'Interpreter', 'none');
exportgraphics(fig, fullfile(figDir, 'network_graph_1N.png'), 'Resolution', res);  close(fig);

%% 4  the accessibility plane: raw coverage with the two diagnostics
fig = figure('Color', 'w', 'Position', [60 60 1150 720], 'Visible', 'off');
ax  = axes(fig);  axStyle(ax);  hold(ax, 'on');
imagesc(ax, caps, 1:nR, P.pairs);  axis(ax, 'tight');  set(ax, 'YDir', 'normal', 'YTick', 1:nR, 'YTickLabel', rungLab);
colormap(ax, 'hot');  cb = colorbar(ax);  cb.Color = 'k';  cb.FontSize = FS - 2;  cb.Label.String = sprintf('direct pairs covered (of %d)', nS);  cb.Label.FontSize = FS - 1;
mAll = double(P.pairs == nS);
if any(mAll(:)) && ~all(mAll(:)), contour(ax, caps, 1:nR, mAll, [0.5 0.5], 'c-', 'LineWidth', LW + 0.5, 'DisplayName', 'all 56 direct pairs covered'); end
mRel = double(P.multileg >= 1);
if any(mRel(:)) && ~all(mRel(:)), contour(ax, caps, 1:nR, mRel, [0.5 0.5], 'g--', 'LineWidth', LW + 0.5, 'DisplayName', 'first relay route closes'); end
xlabel(ax, 'time cap [days]');  ylabel(ax, 'thrust bound [N]');
legend(ax, 'Location', 'southeast', 'FontSize', FS - 1, 'TextColor', 'k', 'Color', 'w');
title(ax, 'Raw coverage: pairs with an accepted transfer at this thrust under the time cap', 'FontSize', FSt, 'Color', 'k');
exportgraphics(fig, fullfile(figDir, 'accessibility_plane_raw.png'), 'Resolution', res);  close(fig);

%% 5  catalog closure down the ladder, by departure family
fig = figure('Color', 'w', 'Position', [60 60 1150 680], 'Visible', 'off');
ax  = axes(fig);  axStyle(ax);  hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
cf  = vertcat(E.closed_frac);  depFam = cellfun(famOf, {E.dep_key}, 'UniformOutput', false);
for ff = {'dro', 'dpo', 'halo_L1', 'halo_L2'}
    sel = strcmp(depFam, ff{1});
    plot(ax, 1:nR, mean(cf(sel, :), 1), '-o', 'Color', famColor(ff{1}), 'LineWidth', LW, 'MarkerSize', 8, 'MarkerFaceColor', famColor(ff{1}), ...
         'DisplayName', sprintf('%s departures (%d edges)', strtok(famName(ff{1}), ' ('), nnz(sel)));
end
set(ax, 'XTick', 1:nR, 'XTickLabel', rungLab);  ylim(ax, [0.3 0.9]);
xlabel(ax, 'thrust rung [N], shallow to deep');  ylabel(ax, 'mean accepted fraction of the 6 x 6 grid');
legend(ax, 'Location', 'southwest', 'FontSize', FS - 1, 'TextColor', 'k', 'Color', 'w');
title(ax, 'Catalog closure down the ladder (an unaccepted cell is unresolved, not infeasible)', 'FontSize', FSt, 'Color', 'k');
exportgraphics(fig, fullfile(figDir, 'ladder_closure.png'), 'Resolution', res);  close(fig);

%% 6  relay closure: assembled time per (L1 halo, tulip) pair; solid = assembled, dotted = optimistic sum
fig = figure('Color', 'w', 'Position', [60 60 1150 680], 'Visible', 'off');
ax  = axes(fig);  axStyle(ax);  hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
cols = parula(numel(best) + 1);
for pp = 1:numel(best)
    lab = regexprep(best(pp).pair, 'halo_L1_tau', 'L1 \\tau=');  lab = regexprep(lab, ' ->tulip_Np', ' to tulip N_p=');
    plot(ax, 1:nR, best(pp).t_route_days, '-o', 'Color', cols(pp, :), 'LineWidth', LW, 'MarkerSize', 7, 'MarkerFaceColor', cols(pp, :), 'DisplayName', lab);
    plot(ax, 1:nR, best(pp).t_sum_days,   ':',  'Color', cols(pp, :), 'LineWidth', 1.5, 'HandleVisibility', 'off');
end
set(ax, 'XTick', 1:nR, 'XTickLabel', rungLab);
xlabel(ax, 'thrust rung [N], shallow to deep');  ylabel(ax, 'two-leg time [days]');
legend(ax, 'Location', 'northwest', 'FontSize', FS - 3, 'TextColor', 'k', 'Color', 'w', 'Interpreter', 'tex');
title(ax, 'L1 halo to L2 halo to tulip: best assembled relay time (solid) and optimistic sum (dotted)', 'FontSize', FSt - 1, 'Color', 'k');
exportgraphics(fig, fullfile(figDir, 'relay_closure.png'), 'Resolution', res);  close(fig);

%% 7  stability versus cost: two panels at slide size
fig = figure('Color', 'w', 'Position', [60 60 1400 600], 'Visible', 'off');
dCJ = stabWS.dCJ;  sDep = stabWS.sDep;  stDep = stabWS.stDep;  fits = stabWS.fits;  resp = stabWS.resp;
cS = [0.15 0.45 0.85];  cU = [0.85 0.30 0.10];
yl = {'min t_f at 1 N [days]', 'retention to 0.5 N'};
for rr = 1:2
    ax = subplot(1, 2, rr);  axStyle(ax);  hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
    y  = resp(rr).y;
    scatter(ax, dCJ(stDep),  y(stDep),  90, cS, 'filled', 'DisplayName', 'stable departure');
    scatter(ax, dCJ(~stDep), y(~stDep), 90, sDep(~stDep), 'filled', 'MarkerEdgeColor', cU, 'LineWidth', 1.5, 'DisplayName', 'unstable departure (color = \sigma_{dep})');
    colormap(ax, 'autumn');
    xx = linspace(min(dCJ), max(dCJ), 50).';
    plot(ax, xx, [ones(50, 1), xx]*fits{rr}.b0, 'k-', 'LineWidth', LW, 'DisplayName', 'fit on \DeltaC_J alone');
    xlabel(ax, '\Delta C_J = C_J(arrival) - C_J(departure)');  ylabel(ax, yl{rr});
    title(ax, sprintf('%s: \\sigma_{dep} coefficient p = %.2f (nominal OLS)', resp(rr).id, fits{rr}.pt1(3)), 'FontSize', FSt, 'Color', 'k');
    if rr == 1, loc = 'northeast'; else, loc = 'southwest'; end
    legend(ax, 'Location', loc, 'FontSize', FS - 3, 'TextColor', 'k', 'Color', 'w');
end
exportgraphics(fig, fullfile(figDir, 'stability_cost.png'), 'Resolution', res);  close(fig);

%% 8  census: the numbers the deck quotes, from the same workspaces
fprintf('\n== slide figures written to %s\n', figDir);
fprintf('sheets %d, nodes %d, rungs %s N\n', nS, numel(N), mat2str(rungs));
fprintf('accepted cells per rung: %s\n', mat2str(sum(cell2mat(arrayfun(@(r) reshape(sum(sum(r.accepted, 1), 2), 1, []), T, 'UniformOutput', false)), 1)));
fprintf('edges accepted at %g N: %d of %d; at %g N: %d of %d\n', rungs(rC), nnz(arrayfun(@(e) e.n_accepted(rC) > 0, E)), nS, rungs(rD), nnz(arrayfun(@(e) e.n_accepted(rD) > 0, E)), nS);
fprintf('first cap covering all %d pairs at %g N: %.3f d\n', nS, rungs(rC), caps(find(P.pairs(rC, :) == nS, 1, 'first')));
fprintf('relay pairs reaching %g N: %d of %d; routes with zero wait: %d of %d closed\n', rungs(rD), nnz([best.deepest_N] == rungs(rD)), numel(best), ...
        nnz(vertcat(R.closed) & vertcat(R.wait_days) < 1e-9), nnz(vertcat(R.closed)));
for rr = 1:2
    fprintf('%s: sigma_dep b = %.4f, p = %.4f, CI [%.4f, %.4f]; joint F p = %.4f\n', resp(rr).id, fits{rr}.b1(3), fits{rr}.pt1(3), fits{rr}.ci1(3, :), fits{rr}.p);
end
for kk = 1:numel(stabWS.verdict)
    if stabWS.verdict(kk).supported, s = 'DETECTED'; else, s = 'NOT DETECTED'; end
    fprintf('%s %s\n', stabWS.verdict(kk).id, s);
end
