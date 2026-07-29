% SMOKE_SEED  Quick sanity check of seed_from_duals before the long LM solves.
setup_paths;
M = 24;
matFile = '/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/results/minfuel/legacy_ms_f1120.mat';
tic;
[Zseed, tJ, info] = seed_from_duals(matFile, 1.12, M);
toc;
fprintf('size(Zseed)=%s  numel expected=%d\n', mat2str(size(Zseed)), 14*M-7);
fprintf('tJ(1)=%.6f tJ(end)=%.6f  M+1=%d numel(tJ)=%d\n', tJ(1), tJ(end), M+1, numel(tJ));
fprintf('beta=%.6g spread=%.2f%% burnAgree=%.2f%% coastAgree=%.2f%% arcCheckErr=%.4e factorSrc=%.3f\n', ...
    info.beta, info.spreadPct, 100*info.burnAgree, 100*info.coastAgree, info.arcCheckErr, info.factorSrc);
