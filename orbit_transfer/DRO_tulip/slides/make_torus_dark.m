function make_torus_dark()
%% Purpose:
%
%   Redraw the 24 x 24 minimum-time phase torus of the 70 mN library of
%   record in the deck's dark style (black background, white bold labels),
%   marking the fastest and slowest certified cells.
%
%% Inputs:
%
%  (none)
%
%% Outputs:
%
%  (file)                   slides/assets/phase_torus_dark.png
%
%% Revision History:
%  M. Casey                                                   (c) 09/22/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here   = fileparts(mfilename('fullpath'));
catMat = fullfile(here, '..', 'indirect', 'results', 'library_70mN_24x24_final', ...
                  'costate_catalog_dro_tulip_70mN.mat');
L = load(catMat);  fn = fieldnames(L);  c = L.(fn{1});  s = c.sheets(1);
tStar = c.constants.tStar_s;

tf = s.tf_nd(:, :, 1)*tStar/86400;               % [nD x nA] days
tf(~s.has_solution(:, :, 1)) = NaN;
sD = s.sD_frac(:);  sA = s.sA_frac(:);
[~, kMin] = min(tf(:));  [~, kMax] = max(tf(:));
[dMin, aMin] = ind2sub(size(tf), kMin);  [dMax, aMax] = ind2sub(size(tf), kMax);

fig = figure('Color', 'k', 'Position', [60 60 1100 950], 'Visible', 'off', ...
             'InvertHardcopy', 'off');
ax = axes(fig, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'FontSize', 15, ...
          'FontWeight', 'bold', 'LineWidth', 1.2, 'Layer', 'top');
hold(ax, 'on');
imagesc(ax, sA, sD, tf, 'AlphaData', ~isnan(tf));
axis(ax, 'xy');  axis(ax, 'tight');
colormap(ax, parula(256));
cb = colorbar(ax, 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.2);
cb.Label.String = 'minimum transfer time  t_f  [days]';
cb.Label.FontSize = 16;
plot(ax, sA(aMin), sD(dMin), 'o', 'MarkerSize', 22, 'LineWidth', 3, 'Color', 'w');
plot(ax, sA(aMax), sD(dMax), 's', 'MarkerSize', 22, 'LineWidth', 3, 'Color', 'w');
text(ax, sA(aMin) + 0.03, sD(dMin) + 0.05, sprintf('fastest %.2f d', tf(kMin)), ...
     'Color', 'w', 'FontSize', 15, 'FontWeight', 'bold');
text(ax, sA(aMax) + 0.03, sD(dMax) + 0.05, sprintf('slowest %.2f d', tf(kMax)), ...
     'Color', 'w', 'FontSize', 15, 'FontWeight', 'bold');
xlabel(ax, 'arrival phase  s_A  (fraction of tulip period)', 'FontSize', 17);
ylabel(ax, 'departure phase  s_D  (fraction of DRO period)', 'FontSize', 17);
xticks(ax, 0:0.1:1);  yticks(ax, 0:0.1:1);
hA = sA(2) - sA(1);  hD = sD(2) - sD(1);
xlim(ax, [sA(1) sA(end)] + [-1 1]*hA/2);  ylim(ax, [sD(1) sD(end)] + [-1 1]*hD/2);
exportgraphics(fig, fullfile(here, 'assets', 'phase_torus_dark.png'), ...
               'Resolution', 150, 'BackgroundColor', 'k');
close(fig);
fprintf('fastest %.3f d at (%.4f, %.4f); slowest %.3f d at (%.4f, %.4f)\n', ...
        tf(kMin), sD(dMin), sA(aMin), tf(kMax), sD(dMax), sA(aMax));
end
