% TASK7C_STEP1B_PACKAGE  Re-package the mesh-refined, production-density
% (25 nodes/rev, N=1104) certified 1 N anchor (task7c_step1b_out.mat:
% Solve_Succeeded, defect=2.68e-14, termErr=0) into run_mintime_mee.m's
% DEFAULT production tag (MEE_mintime_T10.mat), fingerprinted EXACTLY as a
% live run_ladder([... 1], ...) invocation would compute it -- i.e. warm-
% hinted from the 2.5 N rung's own ORIGINAL npr=25 anchor (dL_mt=110.9571,
% N=434), NOT from the npr=15 intermediate used to reach this point. This
% makes the artifact a drop-in cache hit for run_ladder.m's normal per-rung
% flow, so a future `run_ladder([10 5 2.5 1])` call reuses it rather than
% re-solving (same reuse precedent as Task 7b's 2.5 N probe25 packaging).
here = '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo';
cd(here); addpath(here);
resDir = fullfile(here, 'results');

S25anchor = load(fullfile(resDir, 'MEE_ladder_T25.mat'));
prevAnchor = S25anchor.rung.anchor;   % 2.5 N's own npr=25 anchor (the one
                                       % run_ladder.m would actually chain from)
thrustN = 1; nodesPerRev = 25;
dLGuess = prevAnchor.dL_mt * (2.5 / thrustN);
nRevGuess = max(1, round(dLGuess / (2*pi)));

Sref = load(fullfile(resDir, 'task7c_step1b_out.mat'));
o = Sref.o;
assert(o.success && o.maxDefect < 1e-8 && o.termErr < 1e-8, ...
       'task7c_step1b_package: refined point is NOT certified -- refusing to cache');

par = kepler_lt_params(thrustN, 1500, 2000);
N = size(o.X, 2) - 1;

fp = struct('thrustN', thrustN, 'nodesPerRev', nodesPerRev, 'm0kg', 1500, ...
    'ispS', 2000, 'maxIter', 75, 'fuelTag', mee_fuel_tag(thrustN), ...
    'nRevSeed', nRevGuess, 'seedThrB', 0.4, 'roundsMax', 24, 'decadeMin', 0.15, ...
    'alwaysTryStageB', true);
fp.warmStartAnchorDL = dLGuess;
fp.warmStartAnchorN  = prevAnchor.N;

candidate = struct('tfmin', o.tf, 'tfmin_h', o.tf*par.TU_s/3600, 'dL_mt', o.dL, ...
    'revs', o.dL/(2*pi), 'thrustN', thrustN, 'nodesPerRev', nodesPerRev, 'N', N, ...
    'stage', 'B', 'continuationRounds', 7, 'certified', true, 'solverOut', o, 'fp', fp);

finalFile = fullfile(resDir, 'MEE_mintime_T10.mat');
out = candidate;
save(finalFile, 'out');
fprintf('PACKAGED %s: tfmin=%.6f ND (%.4f h) revs=%.4f defect=%.3e certified=%d\n', ...
        finalFile, out.tfmin, out.tfmin_h, out.revs, out.solverOut.maxDefect, out.certified);
fprintf('thr range: [%.6f, %.6f]\n', min(o.U(4,:)), max(o.U(4,:)));
R0 = thrustN * out.tfmin_h;
fprintf('R0 = T*tfmin = %.4f N.h\n', R0);

% --- round-trip verification: EXACTLY the cfg run_ladder.m would build for
% this rung when chaining from the 2.5 N rung (mtCfg.maxIter is irrelevant
% to the fingerprint, per run_ladder's own documented exclusion) -----------
% NOTE: maxIter=75 here MUST match fp.maxIter above -- unlike run_ladder.m's
% OWN separate rung-level fingerprint (which deliberately excludes mtMaxIter
% as checkpoint-granularity-only), run_mintime_mee.m's INTERNAL fp DOES
% fingerprint maxIter, so a future run_ladder call chaining into this cache
% must also pass mtMaxIter=75 (matching this task's own finding/recommendation
% that 75 is the safe/productive sweet spot for this rung) or it will
% correctly fail loud on a fingerprint mismatch rather than silently reuse.
mtCfg = struct('m0kg', 1500, 'ispS', 2000, 'maxIter', 75, 'nRevSeed', nRevGuess);
mtCfg.warmStartAnchor = struct('X', prevAnchor.solverOut.X, 'U', prevAnchor.solverOut.U, ...
    'dL', dLGuess, 'N', prevAnchor.N);
outRT = run_mintime_mee(thrustN, nodesPerRev, mtCfg);
assert(isequal(outRT.tfmin, out.tfmin) && outRT.certified, ...
       'task7c_step1b_package: round-trip load MISMATCH -- packaging is inconsistent');
fprintf('ROUND-TRIP VERIFIED: run_mintime_mee(1,25,...) [default tag, matches run_ladder''s own call] loads the packaged cache cleanly, tfmin=%.6f\n', outRT.tfmin);
