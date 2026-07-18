function fig_switching(ver, out, tag, outDir)
% FIG_SWITCHING  Fig-16 analog (Haberkorn-Martinon-Gergaud 2004): switching
% function, thrust norm, and primer misalignment vs. time -- verifies H1
% (B'p, here the MEE primer, never zero on a burn) and H2 (psi, here S, has
% only PINPOINT zeros -- no singular arc) for a verify_pmp_mee.m result.
%
% Three stacked panels vs. physical time t (paper's own x-axis for its Fig
% 16, not true longitude L -- L would span the same range but be much less
% readable over the tens-to-hundreds of revolutions these MEE solutions
% cover):
%   (1) throttle thr (0/1 step) with true burn arcs shaded, overlaid with
%       vertical markers at the solver's own thr=0.5 switch times (solid)
%       and the interpolated S=0 crossing times (dashed) -- their
%       coincidence (or lack of it) is the direct visual H1/H2 check.
%   (2) switching function S (symlog-scaled for readability -- S's natural
%       scale is NOT O(1) like the paper's psi since this transcription has
%       no built-in normalization) with a zero line; sign(S) should track
%       -sign(thr-0.5).
%   (3) primer misalignment angle (deg, burn nodes only) with the 1 deg gate
%       line, so any breakdown of the primer-alignment gate is visible
%       directly against the campaign's own acceptance threshold.
%
% INPUTS:
%   ver    - verify_pmp_mee.m output struct (needs .t .thr .burn .S .tCross
%            .tSwitch .primerDeg .primerMedianDeg .overallSignPct .pass)
%   out    - casadi_lt_mee result struct (used only for .mf, title annotation)
%   tag    - string label for the title/filename, e.g. 'MEE\_M2\_10N'
%   outDir - directory to write '<tag>_fig_switching.png' into
%
% OUTPUTS: none (writes <outDir>/<tag>_fig_switching.png at 300 dpi)
%
% REFERENCES:
%   [1] Haberkorn, Martinon, Gergaud, JGCD 27(6), 2004, Fig. 16.
%   [2] earth_elliptic_to_geo/verify_pmp_mee.m (the ver struct this consumes).
%   [3] earth_elliptic_to_geo/fig_basin_scatter.m (house figure-script style:
%       theme(fig,'light') try/catch, exportgraphics, 'Visible','off').
if nargin < 4 || isempty(outDir)
    outDir = fullfile(fileparts(mfilename('fullpath')), 'results');
end

t    = ver.t;
thr  = ver.thr;
burn = ver.burn;
S    = ver.S;

fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 900 900]);
try, theme(fig, 'light'); catch, end

% --- panel 1: throttle + switch/crossing alignment --------------------------
ax1 = subplot(3,1,1); hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
% shade true burn arcs
d = diff([0, burn, 0]);
arcStart = find(d == 1);  arcEnd = find(d == -1) - 1;
for a = 1:numel(arcStart)
    xx = [t(arcStart(a)), t(min(arcEnd(a)+1, numel(t))), ...
          t(min(arcEnd(a)+1, numel(t))), t(arcStart(a))];
    yy = [0 0 1 1];
    fill(xx, yy, [1.0 0.85 0.55], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
end
plot(ax1, t, thr, '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.2);
for q = 1:numel(ver.tSwitch)
    xline(ax1, ver.tSwitch(q), '-', 'Color', [0.1 0.35 0.85], 'LineWidth', 1.0);
end
for q = 1:numel(ver.tCross)
    xline(ax1, ver.tCross(q), '--', 'Color', [0.85 0.15 0.15], 'LineWidth', 1.0);
end
ylim(ax1, [-0.1 1.1]);
ylabel(ax1, 'thr');
title(ax1, sprintf('%s -- switching structure (m_f=%.4f kg) -- solver switch (blue solid) vs S=0 (red dashed)', ...
    tag, out.m_f_kg), 'Interpreter', 'none');

% --- panel 2: switching function S (symlog) ---------------------------------
ax2 = subplot(3,1,2); hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
Ssym = sign(S) .* log10(1 + abs(S));           % symlog compression for plotting
plot(ax2, t, Ssym, '-', 'Color', [0.6 0.1 0.6], 'LineWidth', 1.0);
yline(ax2, 0, 'k-', 'LineWidth', 0.8);
ylabel(ax2, 'sign(S)\cdotlog_{10}(1+|S|)');
title(ax2, 'switching function S (symlog-compressed; S<0 should mean thr=1)');

% --- panel 3: primer misalignment angle (burn nodes only) -------------------
ax3 = subplot(3,1,3); hold(ax3, 'on'); grid(ax3, 'on'); box(ax3, 'on');
plot(ax3, t(burn), ver.primerDeg(burn), '.', 'Color', [0.1 0.5 0.2], 'MarkerSize', 6);
yline(ax3, 1.0, 'r--', 'LineWidth', 1.2);
ylabel(ax3, 'primer misalign [deg]');
xlabel(ax3, 't [ND]');
title(ax3, sprintf('primer misalignment on burns (median %.2f deg; 1 deg gate line; sign-agree %.1f%%; pass=%d)', ...
    ver.primerMedianDeg, ver.overallSignPct, ver.pass));
ymax = max(5, min(185, prctile(ver.primerDeg(burn), 99)));
ylim(ax3, [0, ymax]);

linkaxes([ax1, ax2, ax3], 'x');

fn = fullfile(outDir, sprintf('%s_fig_switching.png', tag));
exportgraphics(fig, fn, 'Resolution', 300);
close(fig);
fprintf('WROTE %s\n', fn);
end
