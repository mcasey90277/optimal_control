% BRANCHA_112  Branch-A experiment: M=48, eps=1, maxIter=600 from the dual seed.
% If eps=1 converges (||R||<=1e-9), cascade the warm-started sharpen schedule
% [0.3 ... 1e-4] with maxIter 200/step (guard discipline: warm start advances
% only on success). Full gates evaluated at the sharpest converged eps.
setup_paths;
M = 48;
matFile = '/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/results/minfuel/legacy_ms_f1120.mat';
[Zseed, tJ, info] = seed_from_duals(matFile, 1.12, M);
fprintf('beta=%.4g spread=%.1f%% burnAgree=%.1f%% arcCheckErr=%.2e\n', ...
        info.beta, info.spreadPct, 100*info.burnAgree, info.arcCheckErr);

prob    = ms_problem(1.12, 1);
prob.tJ = tJ;
tA = tic;
out1 = ms_solve(Zseed, prob, 1e-9, 600);
fprintf('BRANCH-A eps=1: ||R||=%.3e iters=%d flag=%d success=%d wall=%.1fs\n', ...
        out1.resNorm, out1.iterations, out1.flag, out1.success, toc(tA));
if ~out1.success
    fprintf('BRANCH_A_STALLED\n');
    return;
end

Zwarm = out1.Z;
best  = out1;  best.eps = 1;
for epsS = [0.3 0.1 0.03 0.01 3e-3 1e-3 3e-4 1e-4]
    prob.epsSmooth = epsS;
    tStep = tic;
    outK  = ms_solve(Zwarm, prob, 1e-9, 200);
    fprintf('sharpen: eps=%.0e  ||R||=%.3e  iters=%d  flag=%d  success=%d  wall=%.1fs\n', ...
            epsS, outK.resNorm, outK.iterations, outK.flag, outK.success, toc(tStep));
    if outK.success
        Zwarm = outK.Z;
        best  = outK;  best.eps = epsS;
    end
end

prob.epsSmooth = best.eps;
traj = ms_traj(best.Z, prob);
fprintf('FINAL: eps=%.0e  ||R||=%.3e  dV=%.4f  switches=%d  bang=%.1f%%\n', ...
        best.eps, best.resNorm, traj.dV_kms, traj.switches, 100*traj.bangFrac);
ok = best.resNorm <= 1e-9 && best.eps <= 1e-3 ...
     && abs(traj.dV_kms - 3.8278) < 0.005 && traj.switches == 12;
save('branchA_result.mat', 'best', 'tJ', 'info', 'ok');
if ok
    fprintf('BRANCH_A_PASS\n');
else
    fprintf('BRANCH_A_GATES_FAIL: eps=%.0e ||R||=%.3e dV=%.4f switches=%d\n', ...
            best.eps, best.resNorm, traj.dV_kms, traj.switches);
end
