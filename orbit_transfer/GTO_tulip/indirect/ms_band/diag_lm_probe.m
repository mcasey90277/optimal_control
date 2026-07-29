% DIAG_LM_PROBE  Bounded LM probe (30 iterations) at the seed_from_duals
% output to check whether the large seed residual (||R||~290) is converging
% geometrically or is genuinely stuck, before committing to the full
% 150-iteration x 3-sharpening-stage solve. Guard discipline: this result is
% NEVER treated as a saved/certified solution regardless of outcome.
setup_paths;
M = 24;
matFile = '/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/results/minfuel/legacy_ms_f1120.mat';
[Zseed, tJ, info] = seed_from_duals(matFile, 1.12, M);
fprintf('beta=%.4g spread=%.1f%% burnAgree=%.1f%% arcCheckErr=%.2e\n', ...
        info.beta, info.spreadPct, 100*info.burnAgree, info.arcCheckErr);
prob = ms_problem(1.12, 1e-2);
prob.tJ = tJ;
tic;
out1 = ms_solve(Zseed, prob, 1e-9, 30);
toc;
fprintf('PROBE DONE: ||R||=%.4e flag=%d iters=%d success=%d\n', ...
    out1.resNorm, out1.flag, out1.iterations, out1.success);
