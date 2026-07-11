function history = refine_loop(seedFile, opts)
% REFINE_LOOP  PMP-residual-driven adaptive mesh refinement (prototype driver).
%
% Each round: measure the current solution's PMP switch-localization score
% (pmp_refine_indicator), refine the sigma mesh where switches are worst
% localized (refine_sigma), build a no-resample warm start (warmstart_on_mesh),
% and re-solve the direct Sundman solver at eps=0 warmTight. Stops when the
% switch times stabilize (max move < local interval width AND |dProp| < propTol
% AND switch count unchanged) or at maxRounds, or if a re-solve fails to
% converge tight. History is persisted every round (crash-recoverable) and a
% summary figure is written. No LM/shooting anywhere -- the indirect machinery
% (via pmp_refine_indicator) is a measurement tool only.
%
% INPUTS:
%   seedFile - prepared seed .mat carrying out.lamDef, factor, tauf0, sigma,
%              rv0, rvf (see prep_refine_seed)
%   opts     - struct: maxRounds [default 4], tag [default seed basename],
%              propTol [default 1e-4 kg], and pass-through indicator/refiner
%              opts M/epsEval/mode/nbr/K/hFloor/maxAdd
%
% OUTPUTS:
%   history - [1xR] struct array, row per measured round (row 1 = seed):
%             nNodes, switches, tauSwitch, maxSwitchMove, prop_kg, dProp,
%             nViol, HresMax, maxDefect, betaSpread, converged, ipoptStatus.
%             Saved to refine_history_<tag>.mat; figure refine_<tag>.png.
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-11-pmp-mesh-refine-design.md

here = fileparts(mfilename('fullpath'));
if nargin < 2, opts = struct(); end
[~, base] = fileparts(seedFile);
if ~isfield(opts, 'maxRounds'), opts.maxRounds = 4;      end
if ~isfield(opts, 'tag'),       opts.tag = base;         end
if ~isfield(opts, 'propTol'),   opts.propTol = 1e-4;     end

p = cr3bp_lt_params(0.025, 15, 2100);
S = load(seedFile);
out = S.out;  sigma = S.sigma;  tauf0 = S.tauf0;
rv0 = S.rv0;  rvf = S.rvf;  factor = S.factor;
tf  = out.X(8, end);
tmpFile = fullfile(here, sprintf('.refine_tmp_%s.mat', opts.tag));

history = struct([]);  prevSwitch = [];  prevProp = NaN;
for r = 1:(opts.maxRounds + 1)
    % --- measure current solution ---
    write_seed(tmpFile, out, factor, tauf0, sigma, rv0, rvf);
    [score, tauSwitch, diag] = pmp_refine_indicator(tmpFile, opts);
    prop = p.m0kg*(1 - out.mf);

    [maxMove, converged, dProp] = deal(NaN, false, NaN);
    if ~isempty(prevSwitch)
        [maxMove, localH] = switch_move(prevSwitch, tauSwitch, diag.tauN);
        dProp = prop - prevProp;
        converged = numel(tauSwitch) == numel(prevSwitch) ...
                    && maxMove < localH && abs(dProp) < opts.propTol;
    end
    history(r).nNodes = numel(sigma); %#ok<*AGROW>
    history(r).switches = numel(tauSwitch);
    history(r).tauSwitch = tauSwitch;
    history(r).maxSwitchMove = maxMove;
    history(r).prop_kg = prop;
    history(r).dProp = dProp;
    history(r).nViol = diag.nViol;
    history(r).HresMax = diag.HresMax;
    history(r).maxDefect = out.maxDefect;
    history(r).betaSpread = diag.betaSpread;
    history(r).converged = converged;
    history(r).ipoptStatus = out.ipoptStatus;
    save(fullfile(here, sprintf('refine_history_%s.mat', opts.tag)), 'history');
    fprintf(['[round %d] nodes=%d sw=%d maxMove=%.2e dProp=%.2e nViol=%d ' ...
             'HresMax=%.2e defect=%.2e conv=%d\n'], r-1, numel(sigma), ...
            numel(tauSwitch), maxMove, dProp, diag.nViol, diag.HresMax, ...
            out.maxDefect, converged);
    if converged || r > opts.maxRounds, break; end

    % --- refine + re-solve ---
    [sigmaNew, isNew, nDropped] = refine_sigma(sigma, score, opts);
    if nDropped > 0
        fprintf('  refine_sigma dropped %d sub-hFloor interval(s)\n', nDropped);
    end
    if nnz(isNew) == 0
        fprintf('  no intervals refinable (all sub-hFloor); stopping.\n');  break;
    end
    [X0, U0] = warmstart_on_mesh(out, sigma, sigmaNew, isNew);
    o = casadi_minfuel_sundman(sigmaNew, tf, rv0, rvf, p.Tmax, p.c, p.muStar, ...
                               X0, U0, tauf0, 1.5, 3000, 0, true);
    if ~(o.success && o.maxDefect < 1e-6)
        fprintf('  re-solve did NOT converge tight (defect=%.2e, %s); stopping.\n', ...
                o.maxDefect, o.ipoptStatus);
        break;
    end
    prevSwitch = tauSwitch;  prevProp = prop;
    sigma = sigmaNew;  out = o;
end
if isfile(tmpFile), delete(tmpFile); end

make_figure(history, diag, here, opts.tag);
end

% -------------------------------------------------------------------------
function write_seed(f, out, factor, tauf0, sigma, rv0, rvf)
% Persist the current solution in sms_seed_duals input layout.
save(f, 'out', 'factor', 'tauf0', 'sigma', 'rv0', 'rvf');
end

function [maxMove, localH] = switch_move(prev, curr, tauN)
% Nearest-neighbor match of switch times between rounds; max move + the
% local interval width at the worst-moving switch (acceptance scale).
maxMove = 0;  localH = Inf;
n = min(numel(prev), numel(curr));
if numel(prev) ~= numel(curr)
    maxMove = Inf;  localH = 1;  return;   % count changed -> not converged
end
sp = sort(prev);  sc = sort(curr);
mv = abs(sp - sc);
[maxMove, w] = max(mv);
% local width near this switch time on the CURRENT node grid
[~, kk] = min(abs(tauN - sc(w)));
kk = min(max(kk, 1), numel(tauN) - 1);
localH = tauN(kk+1) - tauN(kk);
end

function make_figure(history, diag, here, tag)
% Three-panel summary: switch times vs round, S(tau)+throttle, escalation row.
fig = figure('Visible', 'off', 'Position', [50 50 1000 780]);
subplot(3, 1, 1);
hold on;
for r = 1:numel(history)
    ts = history(r).tauSwitch;
    plot(r*ones(size(ts)), ts, 'o');
end
grid on; xlabel('round'); ylabel('switch time \tau');
title(sprintf('%s: switch times vs refinement round', tag), 'Interpreter', 'none');
subplot(3, 1, 2);
yyaxis left;  plot(diag.tauN, diag.Snode, 'LineWidth', 1); ylabel('S (duals)');
hold on; yline(0, 'k-');
plot(diag.tauCr, zeros(size(diag.tauCr)), 'go');
yyaxis right; ylabel('|H_\sigma|'); plot(diag.tauN, diag.Hres, '-');
xlabel('\tau'); grid on; title('final: switching function, S=0 crossings, |H_\sigma|');
subplot(3, 1, 3);
rr = 0:(numel(history)-1);
yyaxis left;  plot(rr, [history.nViol], 's-'); ylabel('nViol');
yyaxis right; plot(rr, [history.HresMax], 'o-'); ylabel('HresMax');
xlabel('round'); grid on; title('escalation dashboard (drive option-2 decision)');
exportgraphics(fig, fullfile(here, sprintf('refine_%s.png', tag)), 'Resolution', 140);
close(fig);
end
