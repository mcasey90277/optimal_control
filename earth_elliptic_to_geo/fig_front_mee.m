function fig_front_mee()
% FIG_FRONT_MEE  Fig-23-adjacent overlay: m_f vs thrust at fixed c_tf=1.5,
% the MEE thrust ladder's answer to the paper's own near-thrust-independence
% observation (Figs 21-23: minimum fuel consumption is nearly insensitive to
% T_max over their tested range).
%
% HONESTY NOTE (per Task-11 brief -- do not fabricate): the paper's Fig 23
% is a multi-c_tf overlay (several c_tf curves, each spanning several
% thrusts). This campaign only ever solved ONE c_tf (1.5) per thrust level
% -- there is no second c_tf column at 5/2.5/1/0.5 N to build a true
% multi-curve overlay from. This figure is therefore the honest
% SINGLE-c_tf version: our 5 certified c_tf=1.5 rungs (10/5/2.5/1/0.5 N)
% plotted against a horizontal band representing the paper's implied
% near-independence range (the 1370-1375 kg band this campaign's own M2
% cross-formulation gate landed inside, per README.md). See README.md's
% MEE-campaign section for the same caveat in prose.
%
% INPUTS:  none
% OUTPUTS: none (writes results/fig_front_mee.png, 300 dpi)
%
% REFERENCES:
%   [1] Haberkorn, Martinon, Gergaud, JGCD 27(6), 2004, Figs 21-23.
%   [2] README.md (M2 cross-formulation gate, 1370-1375 kg paper band).
%   [3] fig_basin_scatter.m / fig_switching.m (house figure-script style).

here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');

thrustN = [10, 5, 2.5, 1, 0.5];

S10  = load(fullfile(resDir, 'MEE_M2_10N.mat'));
S5   = load(fullfile(resDir, 'MEE_M2_5N.mat'));
S2p5 = load(fullfile(resDir, 'MEE_M2_2p5N.mat'));
S1   = load(fullfile(resDir, 'MEE_M2_1N_PSR_psr_final.mat'));
S0p5 = load(fullfile(resDir, 'MEE_M2_0p5N_PSR_psr_final.mat'));

mf = [S10.res.report.m_f_kg, S5.res.report.m_f_kg, S2p5.res.report.m_f_kg, ...
      S1.out.finalOut.m_f_kg, S0p5.out.finalOut.m_f_kg];

% ruling: 0.5 N row is PSR round-4-of-4, budget-limited (not confirmed
% mesh-stable) -- and its anchor tf/tfmin is itself an R0-law estimate
% (ruling 1), not an independently certified min-time solve. Flag it
% visually (hollow marker) rather than presenting it at equal confidence.
isProvisional = [false false false false true];

fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 800 560]);
try, theme(fig, 'light'); catch, end
hold on; grid on; box on

% paper's implied near-independence band (this campaign's own M2
% cross-formulation gate landed inside 1370-1375 kg, per README.md)
bandLo = 1370; bandHi = 1375;
xBand = [0.3 20];
hBand = fill([xBand fliplr(xBand)], [bandLo bandLo bandHi bandHi], ...
    [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

hLine = plot(thrustN, mf, '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.2);
hCert = plot(thrustN(~isProvisional), mf(~isProvisional), 'o', ...
    'MarkerFaceColor', [0.10 0.35 0.85], 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
hProv = plot(thrustN(isProvisional), mf(isProvisional), 'o', ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.10 0.35 0.85], 'MarkerSize', 10, 'LineWidth', 1.8);

for k = 1:numel(thrustN)
    text(thrustN(k)*1.08, mf(k), sprintf('%.2f kg', mf(k)), 'FontSize', 8);
end

set(gca, 'XScale', 'log', 'XDir', 'reverse');
xlim(xBand);
ylim([1360 1380]);
xlabel('T_{max} [N]'); ylabel('m_f [kg]');
title(sprintf('m_f vs thrust at c_{tf}=1.5 (single-c_{tf} overlay -- ours spans %.2f-%.2f kg)', min(mf), max(mf)));
legend([hBand, hLine, hCert, hProv], ...
    {'paper-implied near-independence band (1370-1375 kg)', 'our c_{tf}=1.5 front', ...
     'certified', '0.5 N (PSR round-4-of-4, budget-limited)'}, ...
    'Location', 'southwest', 'FontSize', 8);

fn = fullfile(resDir, 'fig_front_mee.png');
exportgraphics(fig, fn, 'Resolution', 300);
close(fig);
fprintf('WROTE %s\n', fn);
fprintf('m_f range: %.2f - %.2f kg (span %.2f kg) vs paper-implied band %d-%d kg\n', ...
    min(mf), max(mf), max(mf)-min(mf), bandLo, bandHi);
end
