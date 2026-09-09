function I = plot_phase_sheet(Q, outPng, ttl)
%% Purpose:
%
%   The phase-torus figure: minimum transfer time over the (departure phase
%   x arrival phase) grid of a sheet file. Certified cells are coloured by
%   t_f; cells with nothing certified are drawn as GAPS, not as zeros and
%   not as interpolated colour, because "we have no certified transfer here"
%   is the honest reading and an interpolated patch would invent one.
%
%   The fastest and slowest certified cells are marked, since the spread
%   across the torus is the headline: at 70 mN it runs from about 16 d to
%   about 26 d, and arrival phase carries far more of that than departure.
%
%% Inputs:
%
%  Q                        struct                  sheet file
%                                                   (sheet_to_catalog_file):
%                                                   .OK .TF (ND) .sD .sA
%                                                   .rungs .meta
%  outPng                   char (optional)         save path
%  ttl                      char (optional)         title override
%
%% Outputs:
%
%  I                        struct                  .fig .C (the drawn
%                                                   matrix, NaN where
%                                                   uncertified) .nCert
%                                                   .nCells .gridSize
%                                                   .tfMinDays .tfMaxDays
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, outPng = ''; end
if nargin < 3, ttl = ''; end
tStar = Q.meta.tStar;
TF = Q.TF(:,:,1);  OK = Q.OK(:,:,1);
C = nan(size(TF));  C(OK) = TF(OK)*tStar/86400;
nD = numel(Q.sD);  nA = numel(Q.sA);

I = struct('fig', [], 'C', C, 'nCert', nnz(OK), 'nCells', numel(OK), ...
           'gridSize', [nD nA], 'tfMinDays', min(C(:)), 'tfMaxDays', max(C(:)));

fig = figure('Color', 'w', 'Position', [100 100 900 720]);
I.fig = fig;
h = imagesc(1:nA, 1:nD, C);
set(h, 'AlphaData', ~isnan(C));            % gaps stay the axes background
set(gca, 'Color', [0.92 0.92 0.94], 'YDir', 'normal', 'Layer', 'top');
colormap(gca, parula);
cb = colorbar;  cb.Label.String = 'minimum transfer time t_f [days]';
xticks(1:nA);  yticks(1:nD);
xticklabels(arrayfun(@(s) sprintf('%.3f', s), Q.sA, 'UniformOutput', false));
yticklabels(arrayfun(@(s) sprintf('%.3f', s), Q.sD, 'UniformOutput', false));
xtickangle(60);
xlabel('arrival phase s_A  (fraction of the tulip period)');
ylabel('departure phase s_D  (fraction of the DRO period)');
if isempty(ttl)
    ttl = sprintf(['DRO \\rightarrow tulip minimum time over the phase torus\n' ...
                   '%.0f mN, I_{sp} %g s, m_0 %g kg  --  %d of %d phase pairs certified'], ...
                  Q.rungs(1)*1000, Q.meta.ispS, Q.meta.m0kg, I.nCert, I.nCells);
end
title(ttl, 'FontWeight', 'normal');
grid on;  set(gca, 'GridColor', 'w', 'GridAlpha', 0.35);

if I.nCert > 0
    [~, kMin] = min(C(:));  [~, kMax] = max(C(:));
    [iMin, jMin] = ind2sub(size(C), kMin);  [iMax, jMax] = ind2sub(size(C), kMax);
    hold on
    plot(jMin, iMin, 'o', 'MarkerSize', 15, 'LineWidth', 2, 'Color', [0 0.5 0]);
    plot(jMax, iMax, 's', 'MarkerSize', 15, 'LineWidth', 2, 'Color', [0.7 0 0]);
    text(jMin, iMin - 0.55, sprintf('%.2f d', I.tfMinDays), 'HorizontalAlignment', 'center', ...
         'Color', [0 0.4 0], 'FontWeight', 'bold');
    text(jMax, iMax - 0.55, sprintf('%.2f d', I.tfMaxDays), 'HorizontalAlignment', 'center', ...
         'Color', [0.6 0 0], 'FontWeight', 'bold');
end

if ~isempty(outPng), exportgraphics(fig, outPng, 'Resolution', 150); end
end
