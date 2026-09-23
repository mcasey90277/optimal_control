function make_provenance_torus()
%% Purpose:
%
%   How the 70 mN library of record was BUILT, drawn on its minimum-time
%   phase torus: seed roots, the arrival-phase (spine) arclength walk, the
%   departure ribs, and the cells filled by direct solves. Provenance is
%   read from each entry's own entry_notes in the catalog. Rows are shown
%   in rib-walk order: the spine s_D = 0 on top, then s_D = 0.958 down to
%   0.042 (the ribs step in -s_D through the wrap at s_D = 1).
%
%% Inputs:
%
%  (none)
%
%% Outputs:
%
%  (files)                  slides/assets/phase_torus_provenance.png,
%                           slides/assets/provenance_counts.txt
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
[nD, nA] = size(s.has_solution(:, :, 1));
tf  = s.tf_nd(:, :, 1)*tStar/86400;
idx = s.entry_index(:, :, 1);

%% Provenance of every cell, from its entry note:
kind = zeros(nD, nA);           % 1 seed, 2 spine arclength, 3 rib, 4 direct fill
for iD = 1:nD
    for iA = 1:nA
        t = s.entry_notes{idx(iD, iA)};
        if     contains(t, '| seed:'),              kind(iD, iA) = 1;
        elseif contains(t, 'arc crossing'),         kind(iD, iA) = 2;
        elseif contains(t, 'rib step'),             kind(iD, iA) = 3;
        elseif contains(t, 'direct cell solve'),    kind(iD, iA) = 4;
        end
    end
end
assert(all(kind(:) > 0), 'unclassified entry note');
assert(all(kind(1, :) <= 2) && all(kind(2:end, :) >= 3, 'all'), ...
       'seeds/arclength expected on the spine row only');

%% Rows in walk order: sD = 0, then 0.958 ... 0.042:
ord = [1, nD:-1:2];
T = tf(ord, :);  K = kind(ord, :);
sA = s.sA_frac(:);  sDo = s.sD_frac(ord);

fig = figure('Color', 'k', 'Position', [60 60 1250 950], 'Visible', 'off', ...
             'InvertHardcopy', 'off');
ax = axes(fig, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'FontSize', 14, ...
          'FontWeight', 'bold', 'LineWidth', 1.2, 'Layer', 'top', ...
          'Position', [0.09 0.20 0.74 0.74]);
hold(ax, 'on');
imagesc(ax, 1:nA, 1:nD, T);
set(ax, 'YDir', 'reverse');
colormap(ax, parula(256));
cb = colorbar(ax, 'Color', 'w', 'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.2, ...
              'Position', [0.855 0.20 0.025 0.74]);
cb.Label.String = 'minimum transfer time  t_f  [days]';  cb.Label.FontSize = 15;
xlim(ax, [0.5 nA + 0.5]);  ylim(ax, [0.5 nD + 0.5]);

% the spine band
rectangle(ax, 'Position', [0.5 0.5 nA 1], 'EdgeColor', 'w', 'LineWidth', 2.5);

% ribs: a line from the spine down each contiguous run of rib cells
for iA = 1:nA
    r = 1;
    while r < nD && K(r + 1, iA) == 3, r = r + 1; end
    if r > 1
        plot(ax, [iA iA], [1 r], '-', 'Color', [0 0 0], 'LineWidth', 4.5);
        plot(ax, [iA iA], [1 r], '-', 'Color', [1 1 1], 'LineWidth', 2.0);
        plot(ax, iA, r, 'v', 'MarkerSize', 9, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
    end
    % rib cells below a direct-filled gap (a rib that resumed) get a dot
    rr = find(K(:, iA) == 3);  rr = rr(rr > r);
    if ~isempty(rr)
        plot(ax, iA*ones(size(rr)), rr, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w', ...
             'MarkerEdgeColor', 'k');
    end
end
[dr, dc] = find(K == 4);
plot(ax, dc, dr, 'x', 'MarkerSize', 13, 'LineWidth', 4.5, 'Color', 'k');
plot(ax, dc, dr, 'x', 'MarkerSize', 13, 'LineWidth', 2.0, 'Color', [1 0.35 0.35]);
[ar, ac] = find(K == 2);
plot(ax, ac, ar, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
[sr, sc] = find(K == 1);
plot(ax, sc, sr, 'p', 'MarkerSize', 26, 'MarkerFaceColor', [1 0.85 0.3], ...
     'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

xticks(ax, 1:2:nA);  xticklabels(ax, compose('%.2f', sA(1:2:nA)));
yticks(ax, 1:2:nD);  yticklabels(ax, compose('%.3f', sDo(1:2:nD)));
xtickangle(ax, 45);
xlabel(ax, 'arrival phase  s_A', 'FontSize', 16);
ylabel(ax, 'departure phase  s_D   (rows in rib-walk order)', 'FontSize', 16);

% legend, drawn as proxies below the axes
lg = axes(fig, 'Position', [0.09 0.02 0.74 0.08], 'Visible', 'off');  hold(lg, 'on');
xlim(lg, [0 4]);  ylim(lg, [0 1]);
plot(lg, 0.08, 0.5, 'p', 'MarkerSize', 22, 'MarkerFaceColor', [1 0.85 0.3], 'MarkerEdgeColor', 'k');
text(lg, 0.18, 0.5, sprintf('seed root (%d)', nnz(kind == 1)), 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');
plot(lg, 1.05, 0.5, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
text(lg, 1.15, 0.5, sprintf('arclength along s_A (%d)', nnz(kind == 2)), 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');
plot(lg, [2.15 2.15], [0.15 0.85], '-', 'Color', 'w', 'LineWidth', 2.5);
text(lg, 2.25, 0.5, sprintf('rib along s_D (%d)', nnz(kind == 3)), 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');
plot(lg, 3.05, 0.5, 'x', 'MarkerSize', 13, 'LineWidth', 2.5, 'Color', [1 0.35 0.35]);
text(lg, 3.15, 0.5, sprintf('direct fill (%d)', nnz(kind == 4)), 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');

exportgraphics(fig, fullfile(here, 'assets', 'phase_torus_provenance.png'), ...
               'Resolution', 150, 'BackgroundColor', 'k');
close(fig);
fid = fopen(fullfile(here, 'assets', 'provenance_counts.txt'), 'w');
fprintf(fid, 'seed=%d arc=%d rib=%d direct=%d\n', nnz(kind == 1), nnz(kind == 2), ...
        nnz(kind == 3), nnz(kind == 4));
fclose(fid);
fprintf('seed %d, arclength %d, rib %d, direct %d\n', nnz(kind == 1), nnz(kind == 2), ...
        nnz(kind == 3), nnz(kind == 4));
end
