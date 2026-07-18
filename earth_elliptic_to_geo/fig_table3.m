function fig_table3()
% FIG_TABLE3  Table-3 analog (Haberkorn-Martinon-Gergaud 2004): switches and
% revolutions vs thrust (log-x), ours vs the paper's own printed Table 3
% counts, plus the empirical thrust law R0 = T_max*t_f,min panel across the
% four independently CERTIFIED min-time anchors (10/5/2.5/1 N) with the
% anchor-free 0.5 N point shown hollow (excluded from the fit -- see
% ruling below).
%
% All "ours" numbers are read live from the certified .mat deliverables in
% results/ (never hardcoded) except the paper's own published Table 3
% entries, which are literal transcriptions of the reference and have no
% .mat source.
%
% RULINGS carried verbatim from the Task-11 brief / SDD ledger (Tasks 7-9):
%   (1) 0.5 N row: t_f/t_f,min is an R0-LAW ESTIMATE (anchorSource='R0law'),
%       NOT an independently certified min-time solve (Task 9: min-time
%       anchor wall at 0.5 N is an OPEN finding). Its R0 value is therefore
%       exactly the mean of the 4 certified anchors BY CONSTRUCTION
%       (circular) -- shown hollow/gray, excluded from any R0 fit/spread.
%   (2) 0.5 N m_f/switches: PSR round-4-of-4, budget-limited (stopReason=
%       'maxRounds'), not confirmed mesh-stable.
%   (3) The ours-vs-paper revs gap is LADDER-WIDE and systematic
%       (approximately -5.6/-7.2/-7.2/-7.0% at 5/2.5/1/0.5 N) -- one
%       ladder-level footnote (inherited model/paper discrepancy), not a
%       per-rung caveat.
%   (4) 1 N number provenance: PSR-refined m_f=1371.44 kg (round 2 of the
%       1 N PSR port) SUPERSEDES the earlier uniform-mesh value 1370.36 kg;
%       the uniform value is retained only as provenance context, not on
%       the figure.
%
% INPUTS:  none
% OUTPUTS: none (writes results/fig_table3.png, 300 dpi)
%
% REFERENCES:
%   [1] Haberkorn, Martinon, Gergaud, JGCD 27(6), 2004, Table 3.
%   [2] .superpowers/sdd/progress.md, Tasks 6-9 (R0 law + PSR-refined 1N/0.5N
%       certification record).
%   [3] fig_switching.m / fig_basin_scatter.m (house figure-script style:
%       theme(fig,'light') try/catch, exportgraphics, 'Visible','off').

here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');

% =========================================================================
% (1) load "ours" -- fuel-solve structure counts (c_tf=1.5) --------------
% =========================================================================
thrustN = [10, 5, 2.5, 1, 0.5];

S10  = load(fullfile(resDir, 'MEE_M2_10N.mat'));
S5   = load(fullfile(resDir, 'MEE_M2_5N.mat'));
S2p5 = load(fullfile(resDir, 'MEE_M2_2p5N.mat'));
S1   = load(fullfile(resDir, 'MEE_M2_1N_PSR_psr_final.mat'));     % PSR r2, supersedes uniform 1370.36
S0p5 = load(fullfile(resDir, 'MEE_M2_0p5N_PSR_psr_final.mat'));   % PSR r4-of-4, budget-limited

ours_sw   = [S10.res.report.switches, S5.res.report.switches, S2p5.res.report.switches, ...
             S1.out.finalOut.switches, S0p5.out.finalOut.switches];
ours_revs = [S10.res.report.revs, S5.res.report.revs, S2p5.res.report.revs, ...
             S1.out.finalOut.dL/(2*pi), S0p5.out.finalOut.dL/(2*pi)];
ours_mf   = [S10.res.report.m_f_kg, S5.res.report.m_f_kg, S2p5.res.report.m_f_kg, ...
             S1.out.finalOut.m_f_kg, S0p5.out.finalOut.m_f_kg];

% paper Table 3 (literal transcription -- no .mat source)
paper_sw   = [18, 36, 73, 179, 360];
paper_revs = [7.5, 15, 30, 74.5, 149];

% =========================================================================
% (2) load "ours" -- min-time anchors + R0 law ----------------------------
% =========================================================================
Amt = {load(fullfile(resDir,'MEE_mintime_T100.mat')), ...   % 10 N
       load(fullfile(resDir,'MEE_mintime_T50.mat')),  ...   % 5 N
       load(fullfile(resDir,'MEE_mintime_T25.mat')),  ...   % 2.5 N
       load(fullfile(resDir,'MEE_mintime_T10.mat'))};       % 1 N

T4      = cellfun(@(s) s.out.thrustN, Amt);
tfmin4  = cellfun(@(s) s.out.tfmin,   Amt);
tfminH4 = cellfun(@(s) s.out.tfmin_h, Amt);
cert4   = cellfun(@(s) s.out.certified, Amt);
assert(all(T4 == [10 5 2.5 1]), 'fig_table3: mintime anchor thrust order unexpected');
assert(all(cert4 == 1), 'fig_table3: all 4 anchors must be certified');

R0_h    = T4 .* tfminH4;                 % N.h "family" units
R0mean  = mean(R0_h);
R0spreadPct = (max(R0_h) - min(R0_h)) / R0mean * 100;

% 0.5 N: anchor-free by construction (R0-law estimate, ruling 1) -- its R0
% is EXACTLY the mean of the 4 certified anchors, shown hollow/excluded.
T5th   = 0.5;
R0_05  = R0mean;

% =========================================================================
% (3) figure ---------------------------------------------------------------
% =========================================================================
fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1500 460]);
try, theme(fig, 'light'); catch, end

% --- panel 1: switches vs thrust ------------------------------------------
ax1 = subplot(1,3,1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
plot(ax1, thrustN, ours_sw,  '-o', 'Color',[0.10 0.35 0.85], 'MarkerFaceColor',[0.10 0.35 0.85], 'LineWidth',1.4, 'MarkerSize',7);
plot(ax1, thrustN, paper_sw, '--s', 'Color',[0.85 0.35 0.10], 'MarkerFaceColor',[0.85 0.35 0.10], 'LineWidth',1.4, 'MarkerSize',7);
set(ax1, 'XScale', 'log', 'XDir', 'reverse');
xlabel(ax1, 'T_{max} [N]'); ylabel(ax1, 'switches');
title(ax1, 'switches vs thrust');
legend(ax1, {'ours','paper Table 3'}, 'Location','northwest');

% --- panel 2: revs vs thrust -----------------------------------------------
ax2 = subplot(1,3,2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
plot(ax2, thrustN, ours_revs,  '-o', 'Color',[0.10 0.35 0.85], 'MarkerFaceColor',[0.10 0.35 0.85], 'LineWidth',1.4, 'MarkerSize',7);
plot(ax2, thrustN, paper_revs, '--s', 'Color',[0.85 0.35 0.10], 'MarkerFaceColor',[0.85 0.35 0.10], 'LineWidth',1.4, 'MarkerSize',7);
set(ax2, 'XScale', 'log', 'XDir', 'reverse');
xlabel(ax2, 'T_{max} [N]'); ylabel(ax2, 'revolutions');
title(ax2, 'revs vs thrust (paper gap: -5.6 to -7.2%, ladder-wide)');
legend(ax2, {'ours','paper Table 3'}, 'Location','northwest');

% --- panel 3: R0 law T_max*t_f,min --------------------------------------
ax3 = subplot(1,3,3); hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
yline(ax3, R0mean, ':', 'Color',[0.3 0.3 0.3], 'LineWidth',1.2);
yline(ax3, 850, '-.', 'Color',[0.85 0.35 0.10], 'LineWidth',1.0);
plot(ax3, T4, R0_h, 'o', 'MarkerFaceColor',[0.10 0.35 0.85], 'MarkerEdgeColor','k', 'MarkerSize',10);
plot(ax3, T5th, R0_05, 'o', 'MarkerFaceColor','none', 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize',10, 'LineWidth',1.8);
text(ax3, T5th*1.06, R0_05, sprintf('R0-law estimate\n(circular, excluded)'), 'FontSize',8, 'Color',[0.4 0.4 0.4]);
set(ax3, 'XScale', 'log', 'XDir', 'reverse');
xlabel(ax3, 'T_{max} [N]'); ylabel(ax3, 'T_{max}\cdott_{f,min}  [N\cdoth]');
title(ax3, sprintf('R0 law: mean=%.1f N.h, spread=%.2f%% (4 certified anchors)', R0mean, R0spreadPct));
legend(ax3, {sprintf('mean R0 = %.1f N.h', R0mean), 'paper \approx 850 N.h', 'certified anchors', '0.5 N (excluded)'}, ...
    'Location', 'southwest', 'FontSize', 7);

sgtitle(fig, 'Table-3 analog: structure counts + R0 law vs thrust (c_{tf}=1.5)', 'FontWeight','bold');

fn = fullfile(resDir, 'fig_table3.png');
exportgraphics(fig, fn, 'Resolution', 300);
close(fig);
fprintf('WROTE %s\n', fn);

% --- print the rendered table to the console (for the task report) --------
fprintf('\n%-8s %10s %10s %10s %10s %10s %10s\n', 'T[N]','ours_sw','paper_sw','ours_rev','paper_rev','ours_mf','R0[N.h]');
for k = 1:5
    if k <= 4
        r0str = sprintf('%.2f', R0_h(k));
    else
        r0str = sprintf('%.2f*', R0_05);
    end
    fprintf('%-8.1f %10d %10d %10.3f %10.1f %10.2f %10s\n', ...
        thrustN(k), ours_sw(k), paper_sw(k), ours_revs(k), paper_revs(k), ours_mf(k), r0str);
end
fprintf('(* = R0-law estimate, circular by construction, excluded from the fit)\n');
end
