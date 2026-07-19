% RUN_TASK9_DEEP  Top-level driver for Task 9: the deep thrust ladder
% 0.5 -> 0.2 -> 0.1 N. Descends from the certified 1 N rung (Task 7c,
% results/MEE_ladder_T10.mat), warm-chaining the anchor and fuel state
% exactly like run_ladder.m does across rungs, but hand-driving each rung
% via run_task9_rung.m (Task 9's own per-rung node-density strategy) instead
% of calling run_ladder.m itself. Resume-safe end to end: re-running this
% whole script after a crash reloads every already-completed rung's
% DEEP_rung_T<tag>.mat cache (run_task9_rung.m's own guard) and continues
% at the first uncached rung.
%
% BASE OF THE CHAIN: prefers the PSR-refined 1 N fuel result
% (results/MEE_M2_1N_PSR_psr_final.mat, Task 8) over the plain certified 1 N
% fuel solve (results/MEE_M2_1N.mat) as the fuel warm-start root, since it
% is a strictly more mesh-resolved (higher m_f, same switch family)
% converged state of the SAME 1 N problem -- a better starting shape to
% rescale down from. Falls back to the plain certified fuel solve if the
% PSR artifact is not present.
%
% REFERENCES: [1] .superpowers/sdd/task-9-brief.md (this task's spec).
%   [2] run_task9_rung.m (per-rung body). [3] run_ladder.m (the pattern this
%   hand-assembles a variant of).

resDir = fullfile(module_root(), 'results');
addpath(module_root());

S = load(fullfile(resDir, 'MEE_ladder_T10.mat'));   % the certified 1 N rung
rung1N = S.rung;
prevAnchor = rung1N.anchor;
prevThrust = 1.0;

psrFinalFile = fullfile(resDir, 'MEE_M2_1N_PSR_psr_final.mat');
if isfile(psrFinalFile)
    Sp = load(psrFinalFile);
    prevFuelSigma = Sp.out.finalSigma;
    prevFuelX     = Sp.out.finalOut.X;
    prevFuelU     = Sp.out.finalOut.U;
    prevFuelDL    = Sp.out.finalOut.dL;
    fprintf('BASE: using PSR-refined 1 N fuel as warm-start chain root (sw=%d mf=%.4f kg, N=%d)\n', ...
        Sp.out.finalOut.switches, Sp.out.finalOut.m_f_kg, numel(Sp.out.finalSigma) - 1);
else
    Sf = load(fullfile(resDir, [rung1N.fuelTag '.mat']));
    prevFuelSigma = Sf.res.sigma;
    prevFuelX     = Sf.res.fuel.X;
    prevFuelU     = Sf.res.fuel.U;
    prevFuelDL    = Sf.res.fuel.dL;
    fprintf('BASE: using plain-certified 1 N fuel as warm-start chain root (sw=%d mf=%.4f kg, N=%d)\n', ...
        Sf.res.report.switches, Sf.res.report.m_f_kg, numel(Sf.res.sigma) - 1);
end

thrustList          = [0.5 0.2 0.1];
% mtNodesPerRev/mtMaxIter for the 0.5 N rung, as of the LAST attempt this
% session (see task-9-report.md "0.5 N anchor" section for the full,
% honest 7-configuration tuning history -- none certified): 12/rev is the
% mesh that made the most real progress (3 clean rounds down to defect
% 0.0545 at cap=75); 15/rev was tried and was WORSE (stalled earlier, at
% defect ~0.40); cap=150 was tried and abandoned (30+ min on round 0 alone,
% no completion); cap=100 (the value below) was the last live attempt --
% it reached round 1 (defect 0.614 -> 0.474, sub-floor) before being killed
% mid-retry for wall-clock-budget reasons, so it is UNRESOLVED, not
% disproven. The 0.5 N anchor (Stage B continuation) is NOT certified as of
% this commit -- the fuel/PSR stages below were never reached for any of
% the three rungs. Whoever picks this up next should either (a) resume this
% exact cap=100/12-per-rev config from its round00/round01 cache, or
% (b) pursue the R0-law shortcut documented in task-9-report.md (tfMinAnchor
% ~= 223.14/thrustN ND, extrapolated from the 4 already-certified rungs'
% tight T*tfmin spread) to unblock the fuel+PSR stages without a fully
% certified anchor.
mtNodesPerRevList   = [12  12  12];
fuelNodesPerRevList = [12  9   9];
mtMaxIterList       = [100 100 150];
psrMaxRoundsList    = [4   4   4];

deepResults = cell(1, numel(thrustList));
for k = 1:numel(thrustList)
    thrustN = thrustList(k);
    fprintf('\n\n========== TASK 9 DEEP RUNG %d/%d: T=%g N ==========\n', k, numel(thrustList), thrustN);
    opts = struct('mtNodesPerRev', mtNodesPerRevList(k), 'fuelNodesPerRev', fuelNodesPerRevList(k), ...
        'mtMaxIter', mtMaxIterList(k), 'psrMaxRounds', psrMaxRoundsList(k), ...
        'psrGlobalEvery', 3, 'psrGlobalFactor', 1.3);
    deep = run_task9_rung(thrustN, prevThrust, prevAnchor, prevFuelSigma, prevFuelX, ...
        prevFuelU, prevFuelDL, opts);
    deepResults{k} = deep;

    prevAnchor    = deep.anchor;
    prevThrust    = thrustN;
    prevFuelSigma = deep.psr.finalSigma;
    prevFuelX     = deep.psr.finalOut.X;
    prevFuelU     = deep.psr.finalOut.U;
    prevFuelDL    = deep.psr.finalOut.dL;
end

save(fullfile(resDir, 'DEEP_ladder_task9.mat'), 'deepResults');

fprintf('\n\n=== TASK 9 DEEP LADDER COMPLETE ===\n');
fprintf('%-6s %-10s %-8s %-10s %-6s %-8s | %-10s %-6s %-6s %-12s\n', ...
    'T[N]', 'tfmin[ND]', 'revsMT', 'mf_coarse', 'sw_c', 'revs_c', 'mf_PSR', 'sw_P', 'N_P', 'stopReason');
for k = 1:numel(deepResults)
    dd = deepResults{k};
    fprintf('%-6g %-10.4f %-8.3f %-10.2f %-6d %-8.3f | %-10.4f %-6d %-6d %-12s\n', ...
        dd.thrustN, dd.anchor.tfmin, dd.anchor.revs, dd.fuelCoarse.m_f_kg, dd.fuelCoarse.switches, ...
        dd.fuelCoarse.revs, dd.psr.finalOut.m_f_kg, dd.psr.finalOut.switches, ...
        numel(dd.psr.finalSigma) - 1, dd.psr.stopReason);
end
