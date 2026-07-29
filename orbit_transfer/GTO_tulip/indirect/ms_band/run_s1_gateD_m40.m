% RUN_S1_GATED_M40  Gate D sharp-start eps-march at the M = 40 escalation.
%
% Coordinator-authorized escalation (2026-07-10 ruling, item 1): the
% M = 50 sharp-start first step plateaued with flat relays (99.7% ->
% 9.3% < 10% guard), so re-run the identical recipe at M = 40 (fewer,
% longer arcs; joints land on different node interpolants). Same gate:
% eps <= 1e-3, ||R|| <= 1e-9, |dV - 3.8278| < 0.005, switches == 12.
% error() on fail.
setup_paths;
tD = tic;
matFile = '../sundman_minfuel/results/minfuel/legacy_ms_f1120.mat';
[Zseed, prob, info] = sms_seed_duals(matFile, 40, 1e-2);
fprintf(['M40 seed: beta = %.5f  spread %.2f%%  burnAgree %.1f%%  ' ...
         'coastAgree %.1f%%  node1Err %.3e  arcCheckErr %.3e\n'], ...
        info.beta, info.spreadPct, 100*info.burnAgree, ...
        100*info.coastAgree, info.node1Err, info.arcCheckErr);
fprintf('M40 seed ||R(eps=%.3g)|| = %.3e\n', prob.epsSmooth, ...
        norm(sms_residual(Zseed, prob)));

best = eps_march(Zseed, prob, [1e-2 3e-3 1e-3 3e-4 1e-4], 1e-9);

if isempty(best.Z)
    error('FAIL run_s1_gateD_m40: eps_march produced no converged step');
end
probB = prob;  probB.epsSmooth = best.eps;
traj  = sms_traj(best.Z, probB);
fprintf(['M40 Gate D: eps = %.3g   ||R|| = %.3e   dV = %.4f km/s   ' ...
         'switches = %d   bang %.1f%%   maxJD = %.3e   wall %.1f min\n'], ...
        best.eps, best.resNorm, traj.dV_kms, traj.switches, ...
        100*traj.bangFrac, traj.maxJointDefect, toc(tD)/60);
save('sms_gateD_m40.mat', 'best', 'probB', 'info');
okD = best.eps <= 1e-3 && best.resNorm <= 1e-9 ...
      && abs(traj.dV_kms - 3.8278) < 0.005 && traj.switches == 12;
if okD
    fprintf('PASS run_s1_gateD_m40 (1.12x from native duals, M = 40)\n');
else
    error('FAIL run_s1_gateD_m40: eps=%.3g ||R||=%.3e dV=%.4f switches=%d', ...
          best.eps, best.resNorm, traj.dV_kms, traj.switches);
end
