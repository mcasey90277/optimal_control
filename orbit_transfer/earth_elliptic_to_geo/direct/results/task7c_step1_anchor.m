% TASK7C_STEP1_ANCHOR  1 N min-time anchor, small-N-first (nodesPerRev=15),
% warm-started from the certified 2.5 N anchor via the C-law dL rescale.
% Driver script for Task 7c Step 1 (see .superpowers/sdd/task-7c-brief and
% task-7-report.md's Task 7b section for the recipe this follows).
here = '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo';
cd(here);
addpath(here);

S25 = load(fullfile(here, 'results', 'MEE_ladder_T25.mat'));
prevAnchor = S25.rung.anchor;
prevThrust = 2.5;
thrustN = 1;

dLGuess = prevAnchor.dL_mt * (prevThrust / thrustN);
nRevGuess = max(1, round(dLGuess / (2*pi)));
fprintf('C-law dL rescale: dL_mt(2.5N)=%.6f rad -> dLGuess(1N)=%.6f rad -> nRevGuess=%d\n', ...
        prevAnchor.dL_mt, dLGuess, nRevGuess);

mtCfg = struct('m0kg', 1500, 'ispS', 2000, 'maxIter', 75, ...
    'tag', 'MEE_mintime_T10_npr15', 'nRevSeed', nRevGuess);
mtCfg.warmStartAnchor = struct('X', prevAnchor.solverOut.X, 'U', prevAnchor.solverOut.U, ...
    'dL', dLGuess, 'N', prevAnchor.N);

fprintf('Launching run_mintime_mee(1, 15, ...) with mtMaxIter=75, tag=MEE_mintime_T10_npr15\n');
t0 = tic;
out = run_mintime_mee(thrustN, 15, mtCfg);
fprintf('STEP1 DONE in %.1f s: tfmin=%.6f ND (%.2f h) revs=%.4f defect=%.3e termErr=%.3e certified=%d\n', ...
        toc(t0), out.tfmin, out.tfmin_h, out.revs, out.solverOut.maxDefect, out.solverOut.termErr, out.certified);
save(fullfile(here, 'results', 'task7c_step1_out.mat'), 'out');
