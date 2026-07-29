% TEST_MS_REPRODUCE_112  M1(b): MS reproduces the 1.12x direct solution from duals.
%
% Seeds multiple shooting from the direct solution's KKT-dual costates
% (seed_from_duals) and converges via a GUARDED eps-march over the full
% smoothing schedule 1 -> 1e-4 (campaign lesson: interval-averaged duals are
% too noisy as pointwise costate ICs at sharp eps — a direct eps=1e-2 start
% stalled at ||R||~8 with ode113 min-step failures; at eps=1 the throttle
% ramp is gentle and LM repairs the costates against the pinned states,
% then each sharper step warm-starts closer). Warm start advances ONLY on
% out.success (guard discipline).
%
% Gates: reached eps<=1e-3 with ||R||<=1e-9; |dV-3.8278|<0.005 km/s;
% switches == 12 (both measured at the final converged eps).
setup_paths;
M = 48;   % 24 stalls: LM hits nonzero local minima at every eps (see report);
          % halving arc length tames per-arc error amplification of the
          % dual-seeded costates.
matFile = '/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/results/minfuel/legacy_ms_f1120.mat';
[Zseed, tJ, info] = seed_from_duals(matFile, 1.12, M);
fprintf('beta=%.4g spread=%.1f%% burnAgree=%.1f%% arcCheckErr=%.2e\n', ...
        info.beta, info.spreadPct, 100*info.burnAgree, info.arcCheckErr);

epsSchedule = [1 0.3 0.1 0.03 0.01 3e-3 1e-3 3e-4 1e-4];
prob  = ms_problem(1.12, epsSchedule(1));
prob.tJ = tJ;
Zwarm = Zseed;                          % advances only on success (guard)
best  = [];                             % best converged (Z, eps, resNorm)
for epsS = epsSchedule
    prob.epsSmooth = epsS;
    tStep = tic;
    outK  = ms_solve(Zwarm, prob, 1e-9, 150);
    fprintf('eps-march: eps=%.0e  ||R||=%.3e  iters=%d  flag=%d  success=%d  wall=%.1fs\n', ...
            epsS, outK.resNorm, outK.iterations, outK.flag, outK.success, toc(tStep));
    if outK.success
        Zwarm = outK.Z;
        best  = outK;  best.eps = epsS;
    end
end
if isempty(best), error('FAIL test_ms_reproduce_112: no eps step converged'); end

prob.epsSmooth = best.eps;
traj = ms_traj(best.Z, prob);
fprintf('eps=%.0e  ||R||=%.3e  dV=%.4f  switches=%d  bang=%.1f%%\n', ...
        best.eps, best.resNorm, traj.dV_kms, traj.switches, 100*traj.bangFrac);
ok = best.resNorm <= 1e-9 && best.eps <= 1e-3 ...
     && abs(traj.dV_kms - 3.8278) < 0.005 && traj.switches == 12;
if ok
    fprintf('PASS test_ms_reproduce_112\n');
else
    error('FAIL test_ms_reproduce_112: eps=%.0e ||R||=%.3e dV=%.4f switches=%d', ...
          best.eps, best.resNorm, traj.dV_kms, traj.switches);
end
