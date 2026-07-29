% RUN_S1_GATED_DUALMAP  Gate D with the winning dual map + adaptive eps march.
%
% Coordinator directive 2026-07-10 (post GPT-5.6 review): seed from
% legacy_ms_f1120 with the winning midpoint-principled dual map (mode 'd',
% adjudicated by test_sms_dualmap: one-arc costate error 21x below the
% baseline), M = 40, and run the sharp-start ADAPTIVE eps march
% (eps_march_adaptive: start 1e-2, geometric steps with switch
% -displacement feedback, eps_march-style relays, kill-robust state file).
% Gate unchanged: eps <= 1e-3, ||R|| <= 1e-9, |dV - 3.8278| < 0.005,
% switches == 12. error() on fail.
setup_paths;
tD = tic;
matFile = '../sundman_minfuel/results/minfuel/legacy_ms_f1120.mat';
[Zseed, prob, info] = sms_seed_duals(matFile, 40, 1e-2, 'd');
fprintf(['GateD-d seed: mode %s  beta = %.5f  node1Err %.3e  ' ...
         'arcCheckErr %.3e\n'], info.mode, info.beta, info.node1Err, ...
        info.arcCheckErr);
fprintf('GateD-d seed ||R(eps=%.3g)|| = %.3e\n', prob.epsSmooth, ...
        norm(sms_residual(Zseed, prob)));

best = eps_march_adaptive(Zseed, prob, 1e-2, 1e-4, 1e-9, ...
                          'sms_gateD_dualmap_state.mat');

if isempty(best.Z)
    error('FAIL run_s1_gateD_dualmap: no accepted step (wall %.1f min)', ...
          toc(tD)/60);
end
probB = prob;  probB.epsSmooth = best.eps;
traj  = sms_traj(best.Z, probB);
fprintf(['GateD-d: eps = %.3g   ||R|| = %.3e   dV = %.4f km/s   ' ...
         'switches = %d   bang %.1f%%   maxJD = %.3e   wall %.1f min\n'], ...
        best.eps, best.resNorm, traj.dV_kms, traj.switches, ...
        100*traj.bangFrac, traj.maxJointDefect, toc(tD)/60);
save('sms_gateD_dualmap.mat', 'best', 'probB', 'info');
okD = best.eps <= 1e-3 && best.resNorm <= 1e-9 ...
      && abs(traj.dV_kms - 3.8278) < 0.005 && traj.switches == 12;
if okD
    fprintf('PASS run_s1_gateD_dualmap (1.12x from native duals, mode d, M = 40)\n');
else
    error('FAIL run_s1_gateD_dualmap: eps=%.3g ||R||=%.3e dV=%.4f switches=%d', ...
          best.eps, best.resNorm, traj.dV_kms, traj.switches);
end
