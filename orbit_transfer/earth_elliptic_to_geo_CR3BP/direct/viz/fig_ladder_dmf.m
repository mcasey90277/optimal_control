function fig_ladder_dmf()
% FIG_LADDER_DMF  Moon effect down the thrust ladder: the two-regime figure.
%
% Panel (a): dmf(T) [g] vs transfer time in lunar months (log x) -- the
% phase-coherent plateau (~50 g, t_f << lunar month) rolling into the
% phase-averaged decay toward a nonzero secular floor (~30 g at 19 months).
% Panel (b): switch count vs thrust (log-log) with a 1/T reference -- the
% structural growth 19 -> 1724 switches across two thrust decades.
% All dmf values use SAME-CHAIN gain=0 controls. Writes
% doc/figs/ladder_dmf.png (150 dpi) for the campaign note.
%
% INPUTS:  none (loads the 7 rung products + 7 controls from results/)
% OUTPUTS: none (figure file written; values printed)
%
% REFERENCES:
%   [1] doc/cr3bp_geo_phase1_note.tex (thrust-ladder section; consumes this).
%   [2] results/compare_vs_2body.md (same numbers, table form).
here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, '..', 'results');
figDir = fullfile(here, '..', '..', 'doc', 'figs');
if ~exist(figDir, 'dir'), mkdir(figDir); end
rungs = {10,'T10N'; 5,'T5N'; 2.5,'T2p5N'; 1,'T1N'; 0.5,'T0p5N'; 0.2,'T0p2N'; 0.1,'T0p1N'};
nR = size(rungs,1);
T = zeros(nR,1); dmfG = zeros(nR,1); tfMo = zeros(nR,1); sw = zeros(nR,1);
monthDays = 27.32;
for k = 1:nR
    F = load(fullfile(resDir, sprintf('cr3bp_%s_phi0_fuel.mat',      rungs{k,2})), 'm_f_kg','switches','t_days','fp');
    C = load(fullfile(resDir, sprintf('cr3bp_%s_2body_control.mat',  rungs{k,2})), 'm_f_kg');
    T(k)    = rungs{k,1};
    dmfG(k) = 1000*(F.m_f_kg - C.m_f_kg);
    tfMo(k) = F.t_days(end)/monthDays;
    sw(k)   = F.switches;
    fprintf('T=%4.1f N: dmf=%+6.1f g  tf=%6.2f mo  sw=%d\n', T(k), dmfG(k), tfMo(k), sw(k));
end
fig = figure('Visible','off','Position',[60 60 1240 460]);
ax1 = subplot(1,2,1);  hold(ax1,'on');  grid(ax1,'on');  set(ax1,'XScale','log');
xline(ax1, 1, 'k--');  text(ax1, 1.06, 60, 'one lunar month', 'Rotation', 90, 'FontSize', 8);
plot(ax1, tfMo, dmfG, 'o-', 'Color',[0.3 0.3 0.9], 'LineWidth',1.3, ...
     'MarkerFaceColor',[0.15 0.6 0.2], 'MarkerSize',8);
for k = 1:nR
    text(ax1, tfMo(k), dmfG(k)+2.5, sprintf('%g N', T(k)), 'HorizontalAlignment','center','FontSize',8);
end
yline(ax1, 0, 'k-');
xlabel(ax1, 't_f  [lunar months]  (log)'); ylabel(ax1, '\Delta m_f  [g]');
ylim(ax1, [0 65]);
title(ax1, {'Moon effect down the ladder (\phi_0=0, same-chain baselines)', ...
            'phase-coherent plateau \rightarrow phase-averaged decay \rightarrow secular floor'});
ax2 = subplot(1,2,2);  hold(ax2,'on');  grid(ax2,'on');
loglog(ax2, T, sw, 's-', 'Color',[0.85 0.4 0.1], 'LineWidth',1.3, ...
       'MarkerFaceColor',[0.95 0.6 0.2], 'MarkerSize',8);
loglog(ax2, T, sw(1)*T(1)./T, 'k:', 'LineWidth',1);
set(ax2,'XScale','log','YScale','log','XDir','reverse');
xlabel(ax2, 'thrust T  [N]  (log, decreasing)'); ylabel(ax2, 'switch count (nodal)');
legend(ax2, {'certified solves','\propto 1/T reference'}, 'Location','northwest');
title(ax2, 'bang-bang structure growth: 19 \rightarrow 1724 switches');
outPng = fullfile(figDir, 'ladder_dmf.png');
exportgraphics(fig, outPng, 'Resolution', 150);
close(fig);
fprintf('WROTE %s\n', outPng);
end
