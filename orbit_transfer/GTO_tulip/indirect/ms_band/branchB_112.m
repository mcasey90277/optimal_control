% BRANCHB_112  Branch-B: eps=1 MS seed from the ENERGY solution's states
% (energy_f1120.mat, smooth-throttle trajectory — near the smoothed extremal)
% + the bang-bang solution's beta-scaled dual costates (legacy_ms_f1120.mat).
% Rationale: at eps=1 the smoothed BVP extremal lives near the energy
% solution, not the 12-switch one; LM repairs costates against pinned states.
% On eps=1 success: warm-started sharpen [0.3 ... 1e-4], maxIter 200/step,
% guard discipline; full gates at the sharpest converged eps.
setup_paths;
M = 48;
bangFile   = '/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/results/minfuel/legacy_ms_f1120.mat';
energyFile = '/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/results/energy/energy_f1120.mat';

% --- costates from the bang-bang duals (validated beta fit) ---
B = load(bangFile);
p = cr3bp_lt_params(0.025, 15, 2100);
[beta, bInfo] = beta_from_duals(B.out.X, B.out.U, B.out.lamDef, p.c);
fprintf('beta=%.4g spread=%.1f%% burnAgree=%.1f%%\n', beta, bInfo.spreadPct, 100*bInfo.burnAgree);
lamB = -beta*[B.out.lamDef(1:7, :), B.out.lamDef(1:7, end)];   % [7 x nN] on bang grid
tB   = B.out.X(8, :);

% --- states + joints from the energy solution's own grid ---
E   = load(energyFile);
tE  = E.X(8, :);
% same target tf for both (factor 1.12); rescale each time axis exactly to it
tfT = 1.12*6.290694;
tBs = tB*(tfT/tB(end));
tEs = tE*(tfT/tE(end));

tJ = arc_boundaries_tau(tEs, E.X(1:3, :), M, p.muStar);
[tuE, iuE] = unique(tEs);
xJ = interp1(tuE.', E.X(1:7, iuE).', tJ.', 'linear').';       % [7 x (M+1)] energy states
[tuB, iuB] = unique(tBs);
lJ = interp1(tuB.', lamB(:, iuB).', tJ.', 'linear').';        % [7 x (M+1)] bang costates
yJ = [xJ; lJ];
Zseed = ms_pack(yJ(8:14, 1), yJ(:, 2:M));

prob    = ms_problem(1.12, 1);
prob.tJ = tJ;
tA = tic;
out1 = ms_solve(Zseed, prob, 1e-9, 600);
fprintf('BRANCH-B eps=1: ||R||=%.3e iters=%d flag=%d success=%d wall=%.1fs\n', ...
        out1.resNorm, out1.iterations, out1.flag, out1.success, toc(tA));
if ~out1.success
    fprintf('BRANCH_B_STALLED\n');
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
save('branchB_result.mat', 'best', 'tJ', 'ok');
if ok
    fprintf('BRANCH_B_PASS\n');
else
    fprintf('BRANCH_B_GATES_FAIL: eps=%.0e ||R||=%.3e dV=%.4f switches=%d\n', ...
            best.eps, best.resNorm, traj.dV_kms, traj.switches);
end
