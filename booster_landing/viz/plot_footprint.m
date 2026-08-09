function fig = plot_footprint(mc, P, outfile)
% PLOT_FOOTPRINT  Monte-Carlo landing-footprint scatter.
%
% Landing/terminal xy scatter colored by success (mc.ok), pad-radius
% circle, 3-sigma dispersion ellipse from cov(mc.land), and an annotation
% box with the success rate, the failure-mode breakdown
% (n_landed/n_arrest/n_horizon -- an arrest or horizon run never reaches
% a real touchdown, see sim_closed_loop.m's ADAPTATION note), and
% touchdown-speed statistics.
%
% INPUTS:
%   mc      - run_monte_carlo output (.land .ok .vtd .n_landed .n_arrest
%             .n_horizon .success_rate)
%   P       - booster_params (.pad_radius)
%   outfile - (optional) PNG path; exportgraphics at 200 dpi if given
% OUTPUTS:
%   fig - figure handle (Visible off). OWNERSHIP CONTRACT (fix report,
%         2026-08-09): same as plot_pdg_solution.m -- with an output
%         requested (fig = ...), the CALLER owns and must close it; called
%         as a bare statement (nargout==0), this function closes the
%         figure itself so repeated calls across a long -batch process
%         (sweeps, campaigns) do not leak invisible handles.
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-08-08-booster-landing-design.md
if nargin < 3, outfile = ''; end

colOk   = [0 0.447 0.741];      % colorblind-safe blue   -- success
colFail = [0.85 0.325 0.098];   % colorblind-safe orange -- failure

fig = figure('Visible', 'off', 'Position', [100 100 900 850], 'Color', 'w');
ax  = axes(fig);
hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');  axis(ax, 'equal');

ok = logical(mc.ok);
scatter(ax, mc.land(ok,1),  mc.land(ok,2),  55, colOk,   'filled', ...
        'MarkerFaceAlpha', 0.85, 'DisplayName', sprintf('success (n=%d)', nnz(ok)));
scatter(ax, mc.land(~ok,1), mc.land(~ok,2), 55, colFail, 'filled', ...
        'MarkerFaceAlpha', 0.85, 'DisplayName', sprintf('failure (n=%d)', nnz(~ok)));

%% Pad-radius circle:
th = linspace(0, 2*pi, 100);
plot(ax, P.pad_radius*cos(th), P.pad_radius*sin(th), 'k-', 'LineWidth', 2, ...
     'DisplayName', sprintf('pad radius = %.0f m', P.pad_radius));

%% 3-sigma dispersion ellipse from cov(mc.land), centered at the sample mean:
mu = mean(mc.land, 1);
C  = cov(mc.land);
[V, D] = eig(C);
[dvals, ord] = sort(diag(D), 'descend');
V = V(:, ord);
tvec = linspace(0, 2*pi, 100);
ell  = 3 * V * sqrt(diag(max(dvals,0))) * [cos(tvec); sin(tvec)];
plot(ax, mu(1) + ell(1,:), mu(2) + ell(2,:), '--', 'Color', [0.3 0.3 0.3], ...
     'LineWidth', 1.75, 'DisplayName', '3-sigma ellipse');
plot(ax, mu(1), mu(2), 'k+', 'MarkerSize', 10, 'LineWidth', 1.5, ...
     'HandleVisibility', 'off');
plot(ax, 0, 0, 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 14, ...
     'DisplayName', 'pad center');

%% Fixed limits (no autoscale jitter -- pad centered, generous margin):
lim = max([P.pad_radius*1.5; abs(mc.land(:)); 3*sqrt(max(dvals,0))+norm(mu)]);
lim = max(lim * 1.15, P.pad_radius*1.5);
xlim(ax, [-lim lim]);  ylim(ax, [-lim lim]);
xlabel(ax, 'downrange x [m]');  ylabel(ax, 'crossrange y [m]');
title(ax, sprintf('Monte Carlo landing footprint (N=%d runs)', numel(mc.ok)));
legend(ax, 'Location', 'southoutside', 'NumColumns', 2);

%% Annotation box: success rate, failure-mode breakdown, vtd stats:
txt = sprintf(['success rate: %5.1f%% (%d/%d)\n' ...
               'landed=%d  arrest=%d  horizon=%d\n' ...
               'vtd: mean=%.2f  max=%.2f  gate=%.1f m/s'], ...
    100*mc.success_rate, nnz(ok), numel(mc.ok), ...
    mc.n_landed, mc.n_arrest, mc.n_horizon, ...
    mean(mc.vtd), max(mc.vtd), P.vtd_max);
annotation(fig, 'textbox', [0.14 0.70 0.32 0.16], 'String', txt, ...
    'FontName', 'FixedWidth', 'FontSize', 10, 'EdgeColor', [0.7 0.7 0.7], ...
    'BackgroundColor', [1 1 1], 'FitBoxToText', 'on');

if ~isempty(outfile)
    exportgraphics(fig, outfile, 'Resolution', 200);
end
if nargout == 0
    close(fig);        % no caller owns the handle -- close it here (see
                        % OUTPUTS' ownership-contract note above)
end
end
