function deep = run_task9_rung(thrustN, prevThrust, prevAnchor, prevFuelSigma, prevFuelX, ...
    prevFuelU, prevFuelDL, opts)
% RUN_TASK9_RUNG  *** DEPRECATED (2026-07-18) ***
% Superseded by the Table-3 reproducer ENGINE: reproduce_row.m + the
% per-rung recipe registry table3_recipes.m (see README.md "Reproducing
% from scratch (best-found)"). table3_recipes.m now carries 'chain'-anchor
% recipes for 0.2 N (warmFrom 0.5) and 0.1 N (warmFrom 0.2) that cover the
% same anchor + coarse-fuel + PSR ground this function hand-assembled, PLUS
% a keep-best-mass fuel multi-start (this function ran a single coarse fuel
% solve) and from-scratch REPRO_-tag isolation (this function shares the
% campaign's own production tags/cache directory). Prefer
% `reproduce_row(0.2)` / `reproduce_row(0.1)` (after `reproduce_row(0.5)`
% has produced results/repro/REPRO_row_T5.mat) for any new work.
%
% KEPT CALLABLE, NOT DELETED: run_task9_deep.m still calls this function
% with this exact 8-argument signature (prevAnchor/prevFuelSigma/prevFuelX/
% prevFuelU/prevFuelDL threaded by hand from its own caller-side state, not
% from a results/repro/REPRO_row_*.mat file) -- reproduce_row.m's `chain`
% strategy instead loads the previous rung's state itself via its internal
% load_prev helper, reading a REPRO_row_T*.mat file. Those two calling
% conventions do not line up 1:1 (bridging them would mean either changing
% run_task9_deep.m's call site or teaching reproduce_row.m to accept an
% in-memory previous-rung state -- both out of scope for a docs/shim-only
% change), so THIS function's original three-stage body (below) is left
% intact rather than gutted, and a deprecation warning is emitted on every
% call instead. Do not add new callers of this function; drive new deep-rung
% work through reproduce_row.m/table3_recipes.m instead.
%
% ORIGINAL DESCRIPTION (still accurate for the body below): one rung of the
% Task 9 deep thrust ladder (0.5 -> 0.2 -> 0.1 N): a min-time anchor
% (small-N-first, C-law warm-hinted from the PREVIOUS rung's own converged
% anchor, exactly run_ladder.m's mechanism) + a COARSE-BASE fixed-tf fuel
% solve (interp_warmstart-ed from the previous rung's own converged fuel
% trajectory, C-law dL rescaled, direct eps=0 entry via run_transfer_mee.m's
% cfg.warmStart path) + PSR refinement (psr_mee_refine.m, Task 8/9) to
% stabilization or budget.
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
%   [4] reproduce_row.m / table3_recipes.m (Task-4 reproducer engine that
%   supersedes this function; see the DEPRECATED banner above).
%
% STATUS (deprecated 2026-07-18; last live-use status unchanged from the
% prior review): committed but not yet exercised live -- the 0.5 N rung was
% certified via the anchor-free R0-law-estimate path (process/DESIGN_thrust_ladder.md
% footnote 1), not through this function; run_task9_deep.m (its one live
% caller) has not yet been run past 0.5 N, so the first live use of this
% file's body would be the 0.2 N rung -- but new work should reach 0.2/0.1 N
% via reproduce_row.m instead (table3_recipes.m already carries 'chain'
% recipes for both, seeded/not-yet-run).
resDir = fullfile(module_root(), 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
if nargin < 8, opts = struct(); end
warning('run_task9_rung:deprecated', ['run_task9_rung is DEPRECATED -- ' ...
    'superseded by reproduce_row.m + table3_recipes.m (the Table-3 ' ...
    'reproducer engine; see README.md "Reproducing from scratch ' ...
    '(best-found)"). Kept callable only because run_task9_deep.m still ' ...
    'calls this exact signature; do not add new callers.']);
d = @(f, v) optdef(opts, f, v);

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
