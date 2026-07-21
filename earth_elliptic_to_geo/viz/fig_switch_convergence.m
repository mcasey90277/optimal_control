function fig_switch_convergence(rows, outPng)
% FIG_SWITCH_CONVERGENCE  Two-panel figure for the P0 switch-count
% mesh-convergence certification (process/P0_SWITCH_MESH_CONVERGENCE.md):
%   LEFT  - switch count vs nodes/rev, with the converged fine-mesh band shaded
%           and the under-resolved coarse point flagged.
%   RIGHT - final mass vs nodes/rev, with the converged value drawn as an
%           asymptote (shows mass is mesh-converged while the count is not).
%
% INPUTS:
%   rows   - (optional) struct array from verify/meshstudy_switch.m
%            (fields .npr .nSw .m_f_kg). If omitted or [], the CERTIFIED 0.2 N
%            result (this repo, 2026-07-21, 8/16/24/40 nodes/rev, all
%            eps=0-certified, maxDefect ~2-3e-13) is used -- so the figure
%            regenerates without re-running the 3.6 h study.
%   outPng - (optional) output path [default results/p0_switch_mesh_convergence.png]
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/process/P0_SWITCH_MESH_CONVERGENCE.md.
%   [2] earth_elliptic_to_geo/verify/meshstudy_switch.m.
if nargin < 1 || isempty(rows)
    % Certified 0.2 N (T=0.2 N, c_tf=1.5, ~346.7 rev), all rows eps=0-certified:
    npr = [8      16      24      40    ];
    nSw = [823    865     871     863   ];
    mf  = [1377.287 1375.918 1375.836 1375.819];
    Tlab = '0.2 N';
else
    npr = [rows.npr];  nSw = [rows.nSw];  mf = [rows.m_f_kg];  Tlab = '';
end

fine   = npr > min(npr);                 % the refined meshes (exclude the base)
bandLo = min(nSw(fine));  bandHi = max(nSw(fine));
mfConv = mf(end);                        % finest mesh = converged mass

fig = figure('Color','w','Position',[100 100 1000 420]);

% ---- LEFT: switch count ----
ax1 = subplot(1,2,1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
xl = [min(npr)-3, max(npr)+4];
patch(ax1, [xl(1) xl(2) xl(2) xl(1)], [bandLo bandLo bandHi bandHi], ...
    [0.85 0.92 0.85], 'EdgeColor','none', 'FaceAlpha',0.8);
plot(ax1, xl, [bandLo bandLo], '--', 'Color',[0.4 0.6 0.4]);
plot(ax1, xl, [bandHi bandHi], '--', 'Color',[0.4 0.6 0.4]);
plot(ax1, npr(fine), nSw(fine), 'o-', 'Color',[0.10 0.45 0.10], ...
    'MarkerFaceColor',[0.10 0.45 0.10], 'LineWidth',1.6, 'MarkerSize',7);
plot(ax1, npr(~fine), nSw(~fine), 's', 'Color',[0.80 0.25 0.15], ...
    'MarkerFaceColor',[0.95 0.55 0.35], 'LineWidth',1.6, 'MarkerSize',10);
text(ax1, npr(1)+0.6, nSw(1), sprintf('  %d (8/rev:\n  ~%.0f%% undercount)', ...
    nSw(1), 100*(mean([bandLo bandHi])-nSw(1))/mean([bandLo bandHi])), ...
    'Color',[0.6 0.15 0.05], 'FontSize',9, 'VerticalAlignment','middle');
text(ax1, mean(xl), bandHi+2, sprintf('converged band %d-%d', bandLo, bandHi), ...
    'Color',[0.2 0.4 0.2], 'FontSize',9, 'HorizontalAlignment','center');
xlim(ax1, xl); xlabel(ax1,'nodes per revolution'); ylabel(ax1,'switch count');
title(ax1, 'Switch count: NOT mesh-invariant');

% ---- RIGHT: final mass ----
ax2 = subplot(1,2,2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(ax2, xl, [mfConv mfConv], '--', 'Color',[0.3 0.3 0.75], 'LineWidth',1.1);
plot(ax2, npr(fine), mf(fine), 'o-', 'Color',[0.15 0.25 0.6], ...
    'MarkerFaceColor',[0.15 0.25 0.6], 'LineWidth',1.6, 'MarkerSize',7);
plot(ax2, npr(~fine), mf(~fine), 's', 'Color',[0.80 0.25 0.15], ...
    'MarkerFaceColor',[0.95 0.55 0.35], 'LineWidth',1.6, 'MarkerSize',10);
text(ax2, npr(1)+0.6, mf(1), sprintf('  %.2f kg\n  (+%.2f kg high)', ...
    mf(1), mf(1)-mfConv), 'Color',[0.6 0.15 0.05], 'FontSize',9, ...
    'VerticalAlignment','middle');
text(ax2, max(npr), mfConv-0.03, sprintf('converged %.2f kg', mfConv), ...
    'Color',[0.3 0.3 0.75], 'FontSize',9, 'HorizontalAlignment','right', ...
    'VerticalAlignment','top');
xlim(ax2, xl); xlabel(ax2,'nodes per revolution'); ylabel(ax2,'final mass m_f [kg]');
title(ax2, 'Final mass: converged');

sgtitle(fig, sprintf('%s deep rung: primal mesh-convergence (all points \\epsilon=0-certified)', Tlab), ...
    'FontWeight','bold');

if nargin < 2 || isempty(outPng)
    here = fileparts(fileparts(mfilename('fullpath')));   % module root
    outPng = fullfile(here, 'results', 'p0_switch_mesh_convergence.png');
end
exportgraphics(fig, outPng, 'Resolution', 150);
fprintf('wrote %s\n', outPng);
close(fig);
end
