function deep = run_task9_rung(thrustN, prevThrust, prevAnchor, prevFuelSigma, prevFuelX, ...
    prevFuelU, prevFuelDL, opts)
% RUN_TASK9_RUNG  One rung of the Task 9 deep thrust ladder (0.5 -> 0.2 ->
% 0.1 N): a min-time anchor (small-N-first, C-law warm-hinted from the
% PREVIOUS rung's own converged anchor, exactly run_ladder.m's mechanism) +
% a COARSE-BASE fixed-tf fuel solve (interp_warmstart-ed from the previous
% rung's own converged fuel trajectory, C-law dL rescaled, direct eps=0
% entry via run_transfer_mee.m's cfg.warmStart path) + PSR refinement
% (psr_mee_refine.m, Task 8/9) to stabilization or budget.
%
% This is a hand-assembled analog of run_ladder.m's per-rung body, NOT a
% call into run_ladder.m itself: Task 9's brief strategy uses DIFFERENT
% node densities for the anchor (12-15/rev) vs the fuel coarse base
% (8-12/rev, PSR fills in) than run_ladder.m's single uniform nodesPerRev
% (25/rev, both stages) used for the already-certified 10/5/2.5/1 N rungs --
% reusing run_ladder.m directly across a MIXED nodesPerRev list would also
% throw its own cache-fingerprint guard on every already-certified rung.
%
% RESUME-SAFE: caches the WHOLE rung (anchor + coarse fuel + PSR final) to
% resDir/DEEP_rung_T<thrustTag>.mat; if that file exists, it is loaded
% verbatim and NONE of the three stages below re-run (each stage's OWN
% internal caching -- run_mintime_mee's per-round files, run_transfer_mee's
% probe/seed/warmdirect files, psr_mee_refine's per-round files -- already
% makes a mid-stage crash resumable at finer grain; this outer cache just
% skips re-driving an already-fully-certified rung entirely).
%
% INPUTS:
%   thrustN        - this rung's thrust [N]
%   prevThrust     - the PREVIOUS rung's thrust [N] (for the C-law rescale)
%   prevAnchor     - the PREVIOUS rung's run_mintime_mee output struct
%                     (.dL_mt, .solverOut.X/.U, .N)
%   prevFuelSigma/X/U/DL - the PREVIOUS rung's converged fuel trajectory
%                     (sigma grid, X [7x(Np+1)], U [4x(Np+1)], dL scalar) --
%                     the caller decides whether this is the plain
%                     certified fuel solve or a PSR-refined one; this
%                     function just consumes whatever it is handed
%   opts           - struct, all optional:
%                    .mtNodesPerRev [12], .fuelNodesPerRev [12],
%                    .mtMaxIter [100] (run_mintime_mee's per-round IPOPT
%                    cap, kept small since small-N anchors need far fewer
%                    iterations per round), .psrMaxRounds [4],
%                    .psrGlobalEvery [3], .psrGlobalFactor [1.3],
%                    .fuelMaxIter [1500], .m0kg [1500], .ispS [2000]
%
% OUTPUTS:
%   deep - struct: .thrustN, .anchor (run_mintime_mee out), .fuelTag,
%          .fuelCoarse (run_transfer_mee .report, the PRE-PSR coarse-base
%          certified solve), .fuelCoarseN (its node count), .psr
%          (psr_mee_refine out -- .finalSigma/.finalOut are the BEST
%          CERTIFIED state reached, guaranteed certified per
%          psr_mee_refine's own contract), .tf (fixedtf target, ND),
%          .wallAnchor/.wallFuel/.wallPsr [s], .opts (resolved)
%
% REFERENCES: [1] run_ladder.m (the per-rung body this hand-assembles a
%   variant of, with per-stage node densities the brief calls for).
%   [2] run_mintime_mee.m, run_transfer_mee.m, psr_mee_refine.m (the three
%   stages called here). [3] .superpowers/sdd/task-9-brief.md (this task's
%   STRATEGY section, the source of the node-density choices).
here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
if nargin < 8, opts = struct(); end
d = @(f, v) getdef_t9(opts, f, v);

opts.mtNodesPerRev   = d('mtNodesPerRev', 12);
opts.fuelNodesPerRev = d('fuelNodesPerRev', 12);
opts.mtMaxIter       = d('mtMaxIter', 100);
opts.psrMaxRounds    = d('psrMaxRounds', 4);
opts.psrGlobalEvery  = d('psrGlobalEvery', 3);
opts.psrGlobalFactor = d('psrGlobalFactor', 1.3);
opts.fuelMaxIter     = d('fuelMaxIter', 1500);
opts.m0kg            = d('m0kg', 1500);
opts.ispS            = d('ispS', 2000);

tagStr = strrep(sprintf('%g', thrustN), '.', 'p');
deepRungFile = fullfile(resDir, sprintf('DEEP_rung_T%s.mat', tagStr));
if isfile(deepRungFile)
    Sd = load(deepRungFile);
    deep = Sd.deep;
    fprintf('[cached deep rung] loaded %s (T=%g N, mf=%.4f kg sw=%d N=%d stopReason=%s)\n', ...
        deepRungFile, deep.thrustN, deep.psr.finalOut.m_f_kg, deep.psr.finalOut.switches, ...
        numel(deep.psr.finalSigma) - 1, deep.psr.stopReason);
    return;
end

% --- ANCHOR: min-time, small-N-first, C-law warm-hinted -------------------
dLGuess   = prevAnchor.dL_mt * (prevThrust / thrustN);
nRevGuess = max(1, round(dLGuess / (2*pi)));
mtCfg = struct('m0kg', opts.m0kg, 'ispS', opts.ispS, 'maxIter', opts.mtMaxIter, ...
    'nRevSeed', nRevGuess, 'warmStartAnchor', struct('X', prevAnchor.solverOut.X, ...
        'U', prevAnchor.solverOut.U, 'dL', dLGuess, 'N', prevAnchor.N));
fprintf(['  [anchor] T=%g N: warm hint dL_guess=%.4f rad -> nRevSeed=%d ' ...
         '(mtNodesPerRev=%d, mtMaxIter=%d)\n'], thrustN, dLGuess, nRevGuess, ...
        opts.mtNodesPerRev, opts.mtMaxIter);
tA = tic;
anchorOut = run_mintime_mee(thrustN, opts.mtNodesPerRev, mtCfg);
wallAnchor = toc(tA);
assert(anchorOut.certified, 'run_task9_rung:anchorUncertified', ...
    'T=%g N anchor did NOT certify -- rung BLOCKED at the anchor stage', thrustN);
fprintf('  [anchor DONE] tfmin=%.4f ND (%.2f h) revs=%.3f wall=%.1fs\n', ...
    anchorOut.tfmin, anchorOut.tfmin_h, anchorOut.revs, wallAnchor);

% --- FUEL: coarse base, warm-started from the prior rung's converged fuel,
% direct eps=0 entry (run_transfer_mee's cfg.warmStart path) -------------
fuelTag = mee_fuel_tag(thrustN);
dLGuessFuel = prevFuelDL * (prevThrust / thrustN);
fuelCfg = struct('thrustN', thrustN, 'ctf', 1.5, 'tfMinAnchor', anchorOut.tfmin, ...
    'tag', fuelTag, 'nodesPerRev', opts.fuelNodesPerRev, 'maxIter', opts.fuelMaxIter, ...
    'm0kg', opts.m0kg, 'ispS', opts.ispS, 'warmStart', struct('sigma', prevFuelSigma, ...
        'X', prevFuelX, 'U', prevFuelU, 'dL', dLGuessFuel));
fprintf('  [fuel coarse] T=%g N: fuelNodesPerRev=%d, dLGuessFuel=%.4f rad (revsGuess=%.3f)\n', ...
    thrustN, opts.fuelNodesPerRev, dLGuessFuel, dLGuessFuel / (2*pi));
tF = tic;
fuelRes = run_transfer_mee(fuelCfg);
wallFuel = toc(tF);
assert(fuelRes.report.certified, 'run_task9_rung:fuelUncertified', ...
    'T=%g N coarse fuel solve did NOT certify (defect=%.2e) -- rung BLOCKED at the fuel stage', ...
    thrustN, fuelRes.report.defect);
fprintf('  [fuel coarse DONE] N=%d mf=%.4f kg sw=%d revs=%.3f wall=%.1fs\n', ...
    numel(fuelRes.sigma) - 1, fuelRes.report.m_f_kg, fuelRes.report.switches, ...
    fuelRes.report.revs, wallFuel);

% --- PSR refinement to stabilization (or budget/resolveFailed) -----------
psrOpts = struct('tag', [fuelTag '_PSR'], 'maxRounds', opts.psrMaxRounds, ...
    'maxIter', opts.fuelMaxIter, 'nbr', 2, 'globalEvery', opts.psrGlobalEvery, ...
    'globalFactor', opts.psrGlobalFactor);
tP = tic;
psrOut = psr_mee_refine(fuelRes, psrOpts);
wallPsr = toc(tP);
fprintf(['  [PSR DONE] stopReason=%s certified=%d Nfinal=%d mf=%.4f kg sw=%d ' ...
         'wall=%.1fs\n'], psrOut.stopReason, psrOut.certified, numel(psrOut.finalSigma) - 1, ...
        psrOut.finalOut.m_f_kg, psrOut.finalOut.switches, wallPsr);

deep = struct('thrustN', thrustN, 'anchor', anchorOut, 'fuelTag', fuelTag, ...
    'fuelCoarse', fuelRes.report, 'fuelCoarseN', numel(fuelRes.sigma) - 1, ...
    'psr', psrOut, 'tf', fuelRes.tf, 'wallAnchor', wallAnchor, 'wallFuel', wallFuel, ...
    'wallPsr', wallPsr, 'opts', opts);
save(deepRungFile, 'deep');
fprintf('[deep rung DONE] saved %s\n', deepRungFile);
end

% ---------------------------------------------------------------------------
function v = getdef_t9(s, f, dflt)
% GETDEF_T9  Optional-field default (mirrors run_ladder.m's local helper).
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
