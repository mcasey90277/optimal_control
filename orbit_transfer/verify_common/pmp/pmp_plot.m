function pmp_plot(R, tag, outPng)
% PMP_PLOT  Six-panel diagnostic for a continuous PMP residual run.
%
% The panels are ordered by what they answer, not by convenience:
%
%   (a) RESIDUAL MAP. R_x and R_lam per interval against physical time, log
%       scale, switch-containing cells marked. Shows WHERE the discretization
%       struggles. The pre-registered expectation was that switch cells would
%       dominate; this panel is what tests it.
%   (b) THE COMPARISON THAT MATTERS. The per-interval DISCRETE defect and the
%       CONTINUOUS propagation residual on the same log axis. The gap between
%       them is discretization error made visible -- how much better the
%       solution satisfies the equations it was solved against than the
%       equations we actually care about. No ladder, no branch argument, no
%       extrapolation.
%   (c) CONSERVED QUANTITY. In this L-domain transcription the Hamiltonian in
%       sigma is NOT constant (the dynamics carry sigma explicitly through
%       L = pi + sigma*DeltaL). The conserved quantity is the TIME COSTATE:
%       t enters only through its own trivial equation, so lam_t is constant on
%       an extremal and equals -H_t. A drifting lam_t means the extremal
%       structure itself is off, independently of any mesh question.
%   (d) RESIDUAL vs LOCAL STEP, log-log, coloured by switch/interior. Because
%       the physical step varies by more than an order of magnitude across one
%       transfer, a SINGLE solution already samples many step sizes. The slope
%       is suggestive of a local order -- but see the caveat in the panel
%       title: step size is confounded with orbital geometry (large steps sit
%       near apogee, small near perigee), so this is a diagnostic, NOT an
%       observed order in the verification-and-validation sense.
%   (e) CONTROL CONTEXT. Throttle against time with switch times marked, so the
%       residual map above can be read against the bang-bang structure.
%   (f) RESIDUAL DISTRIBUTION, switch cells versus interior, as histograms.
%       Turns "do switches dominate?" from an eyeball question into a
%       comparison of two distributions.
%
% INPUTS:
%   R      - output of pmp_residual_mee [struct]
%   tag    - label for the figure title [char]
%   outPng - (optional) write the figure here [char]; omit for screen only
%
% OUTPUTS: none (draws, and optionally writes, a figure)
%
% REFERENCES:
%   [1] pmp_residual.m, pmp_residual_mee.m.
%   [2] docs/superpowers/plans/2026-07-25-mesh-convergence-study.md, PHASE 3.

if nargin < 3, outPng = ''; end
t   = R.X(7,:);                 % physical time at nodes
tm  = 0.5*(t(1:end-1) + t(2:end));
sw  = R.isSwitchCell;
rx  = R.Rx;  rl = R.Rlam;  dc = R.defectCell;
h   = R.hPhys;

fig = figure('Color','w','Position',[60 60 1500 900]);

% ---- (a) residual map ----------------------------------------------------
subplot(2,3,1);
semilogy(tm, rx, '.', 'MarkerSize', 5, 'Color', [0.15 0.35 0.75]); hold on; grid on;
semilogy(tm, rl, '.', 'MarkerSize', 5, 'Color', [0.85 0.35 0.10]);
if any(sw), semilogy(tm(sw), rx(sw), 'ko', 'MarkerSize', 5, 'LineWidth', 0.8); end
xlabel('time [ND]'); ylabel('per-interval residual');
title('(a) BVP residual map');
legend({'R_x','R_\lambda','switch cells'}, 'Location','southeast', 'FontSize', 7);

% ---- (b) discrete vs continuous -- the money panel -----------------------
subplot(2,3,2);
semilogy(tm, dc, '.', 'MarkerSize', 5, 'Color', [0.2 0.6 0.2]); hold on; grid on;
semilogy(tm, rx, '.', 'MarkerSize', 5, 'Color', [0.15 0.35 0.75]);
xlabel('time [ND]'); ylabel('norm');
rat = median(rx,'omitnan') / max(median(dc,'omitnan'), realmin);
title(sprintf('(b) DISCRETE defect vs CONTINUOUS residual\nmedian ratio %.1e', rat));
legend({'discrete defect (what was solved)','continuous residual (what we want)'}, ...
       'Location','east', 'FontSize', 7);

% ---- (c) conserved quantity ----------------------------------------------
subplot(2,3,3);
lt = R.lam(7,:);
plot(t, lt, '-', 'LineWidth', 1.0, 'Color', [0.4 0.2 0.6]); grid on;
xlabel('time [ND]'); ylabel('\lambda_t');
title(sprintf('(c) time costate (= -H_t), CoV %.2e\nconstant on an extremal', R.lamTimeCoV));

% ---- (d) residual vs local step ------------------------------------------
subplot(2,3,4);
good = isfinite(rx) & h > 0;
loglog(h(good & ~sw), rx(good & ~sw), '.', 'MarkerSize', 6, ...
       'Color', [0.15 0.35 0.75]); hold on; grid on;
if any(good & sw)
    loglog(h(good & sw), rx(good & sw), 'o', 'MarkerSize', 5, ...
           'Color', [0.85 0.2 0.1], 'LineWidth', 0.9);
end
cf = polyfit(log(h(good)), log(rx(good)), 1);
hh = linspace(min(h(good)), max(h(good)), 20);
loglog(hh, exp(cf(2))*hh.^cf(1), 'k--', 'LineWidth', 1.2);
xlabel('local physical step h'); ylabel('R_x');
title(sprintf('(d) residual vs step, slope %.2f\nDIAGNOSTIC ONLY: h confounded with geometry', cf(1)));
legend({'interior','switch','fit'}, 'Location','southeast', 'FontSize', 7);

% ---- (e) control context --------------------------------------------------
subplot(2,3,5);
stairs(t, R.U(4,:), '-', 'LineWidth', 0.9, 'Color', [0.2 0.2 0.2]); hold on; grid on;
yl = ylim;
for q = 1:numel(R.switchTimes)
    plot([R.switchTimes(q) R.switchTimes(q)], yl, '-', 'Color', [0.85 0.4 0.4 0.5]);
end
xlabel('time [ND]'); ylabel('throttle');
title(sprintf('(e) bang-bang structure, %d switches', numel(R.switchTimes)));
ylim([-0.1 1.1]);

% ---- (f) distributions ----------------------------------------------------
subplot(2,3,6);
lg = @(v) log10(v(isfinite(v) & v > 0));
histogram(lg(rx(~sw)), 20, 'FaceColor', [0.15 0.35 0.75], 'FaceAlpha', 0.6, ...
          'Normalization', 'probability'); hold on; grid on;
if any(sw)
    histogram(lg(rx(sw)), 10, 'FaceColor', [0.85 0.2 0.1], 'FaceAlpha', 0.6, ...
              'Normalization', 'probability');
end
xlabel('log_{10} R_x'); ylabel('fraction');
title(sprintf('(f) switch vs interior\nmedians %.2e / %.2e', ...
      median(rx(sw),'omitnan'), median(rx(~sw),'omitnan')));
legend({'interior','switch'}, 'Location','northwest', 'FontSize', 7);

sgtitle(sprintf('Continuous PMP residuals -- %s  (self-gate %.2e, %s)', ...
    tag, R.selfDefect, local_tf(R.gatePass,'PASS','FAIL')), 'FontWeight','bold');

if ~isempty(outPng)
    exportgraphics(fig, outPng, 'Resolution', 140);
    fprintf('figure -> %s\n', outPng);
end
end

% ---------------------------------------------------------------------------
function s = local_tf(c, a, b)
% LOCAL_TF  Pick a label by condition. INPUTS: c,a,b  OUTPUTS: s [char]
if c, s = a; else, s = b; end
end
