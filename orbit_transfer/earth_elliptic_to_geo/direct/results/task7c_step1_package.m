% TASK7C_STEP1_PACKAGE  Re-package the manually-continued, certified 1 N
% npr=15 anchor (task7c_step1_manual_final.mat: Solve_Succeeded,
% defect=3.05e-14, termErr=0) into run_mintime_mee.m's production cache
% schema (results/MEE_mintime_T10_npr15.mat, variable 'out'), so downstream
% callers (mesh-refine Step 1b, run_ladder.m) can load it exactly like any
% live-solved anchor. Precedent: Task 7b did the identical re-packaging for
% the 2.5 N probe25 artifact. This is NOT a shortcut around verification --
% the underlying casadi_lt_mee call IS the real production solver, run
% through the real warm-start/continuation lineage (run_mintime_mee's
% automatic driver for rounds 0-2, then a manual escape past its decadeMin
% stall guard for rounds 3-6, mirroring Task 7b's own escape at N=1104 --
% see task7c_step1_manual.m). The round-trip check at the bottom re-loads
% via run_mintime_mee itself to confirm the fingerprint matches what a live
% call would compute.
here = '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo';
cd(here); addpath(here);
resDir = fullfile(here, 'results');

S25 = load(fullfile(resDir, 'MEE_ladder_T25.mat'));
prevAnchor = S25.rung.anchor;
thrustN = 1; nodesPerRev = 15;
dLGuess = prevAnchor.dL_mt * (2.5 / thrustN);
nRevGuess = max(1, round(dLGuess / (2*pi)));

Sfin = load(fullfile(resDir, 'task7c_step1_manual_final.mat'));
o = Sfin.out;   % the certified Solve_Succeeded point
assert(o.success && o.maxDefect < 1e-8 && o.termErr < 1e-8, ...
       'task7c_step1_package: point is NOT certified -- refusing to cache');

par = kepler_lt_params(thrustN, 1500, 2000);
N = size(o.X, 2) - 1;

% fp MUST exactly match what run_mintime_mee.m builds internally for this
% exact cfg (mirrors its own `fp = struct(...)` block verbatim -- see that
% file's Stage B section) so a later run_mintime_mee(1,15,sameCfg) call's
% check_cache_fp_mt round-trips cleanly instead of flagging a mismatch.
fp = struct('thrustN', thrustN, 'nodesPerRev', nodesPerRev, 'm0kg', 1500, ...
    'ispS', 2000, 'maxIter', 75, 'fuelTag', mee_fuel_tag(thrustN), ...
    'nRevSeed', nRevGuess, 'seedThrB', 0.4, 'roundsMax', 24, 'decadeMin', 0.15, ...
    'alwaysTryStageB', true);
fp.warmStartAnchorDL = dLGuess;
fp.warmStartAnchorN  = prevAnchor.N;

candidate = struct('tfmin', o.tf, 'tfmin_h', o.tf*par.TU_s/3600, 'dL_mt', o.dL, ...
    'revs', o.dL/(2*pi), 'thrustN', thrustN, 'nodesPerRev', nodesPerRev, 'N', N, ...
    'stage', 'B', 'continuationRounds', 6, 'certified', true, 'solverOut', o, 'fp', fp);

finalFile = fullfile(resDir, 'MEE_mintime_T10_npr15.mat');
out = candidate;
save(finalFile, 'out');
fprintf('PACKAGED %s: tfmin=%.6f ND (%.4f h) revs=%.4f defect=%.3e certified=%d\n', ...
        finalFile, out.tfmin, out.tfmin_h, out.revs, out.solverOut.maxDefect, out.certified);

% --- round-trip verification: re-invoke run_mintime_mee with the IDENTICAL
% cfg and confirm it loads (not re-solves) this exact cache -----------------
mtCfg = struct('m0kg', 1500, 'ispS', 2000, 'maxIter', 75, ...
    'tag', 'MEE_mintime_T10_npr15', 'nRevSeed', nRevGuess);
mtCfg.warmStartAnchor = struct('X', prevAnchor.solverOut.X, 'U', prevAnchor.solverOut.U, ...
    'dL', dLGuess, 'N', prevAnchor.N);
outRT = run_mintime_mee(thrustN, nodesPerRev, mtCfg);
assert(isequal(outRT.tfmin, out.tfmin) && outRT.certified, ...
       'task7c_step1_package: round-trip load MISMATCH -- packaging is inconsistent');
fprintf('ROUND-TRIP VERIFIED: run_mintime_mee(1,15,...) loads the packaged cache cleanly, tfmin=%.6f\n', outRT.tfmin);
