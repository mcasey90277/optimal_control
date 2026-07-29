% DIAG_S1_GATED2  Trace capture for the SHARP-START Gate D plateau (Task S1).
%
% Same analysis as DIAG_S1_GATED but at the sharp end: reproduces the
% capped eps = 1e-2 solve (+ relay-1 equivalent) from the native dual
% seed (M = 50) and analyzes the returned iterate: residual row classes,
% first-order optimality, Gauss-Newton consistency, SVD of J, trajectory
% diagnostics, distance moved from the seed. Saves diag_s1_gateD2.mat.
% Analysis only (coordinator ruling 2026-07-10: on a sharp-start plateau,
% capture traces and report BLOCKED — no further strategies).
setup_paths;
tD = tic;
matFile = '../sundman_minfuel/results/minfuel/legacy_ms_f1120.mat';
[Zseed, prob, info] = sms_seed_duals(matFile, 50, 1e-2);
out1 = ms_solve(Zseed, prob, 1e-9, 200);
out2 = ms_solve(out1.Z, prob, 1e-9, 200);          % relay-1 equivalent
Z = out2.Z;

[R, J] = sms_residual(Z, prob);
M  = numel(prob.sJ) - 1;
Rm = reshape(R(1:16*(M-1)), 16, M-1);
Rt = R(16*(M-1)+1:end);
fprintf('\n=== sharp-start plateau analysis (eps=1e-2, ||R|| = %.4e) ===\n', norm(R));
fprintf('rowclass: r,v %.3e  m %.3e  t %.3e  lamR %.3e  lamV %.3e  lamM %.3e  lamT %.3e\n', ...
        max(max(abs(Rm(1:6,:)))), max(abs(Rm(7,:))), max(abs(Rm(8,:))), ...
        max(max(abs(Rm(9:11,:)))), max(max(abs(Rm(12:14,:)))), ...
        max(abs(Rm(15,:))), max(abs(Rm(16,:))));
[~, wj] = max(max(abs(Rm), [], 1));
fprintf('worst joint: %d of %d (sigma = %.2f)\n', wj, M-1, prob.sJ(wj+1));
fprintf('terminal rows: rv %.3e  lamM %.3e  t %.3e\n', ...
        max(abs(Rt(1:6))), abs(Rt(7)), abs(Rt(8)));
fprintf('optimality ||J''R||_inf = %.4e\n', norm(J.'*R, inf));

dGN    = -(J\R);
relres = norm(J*dGN + R)/norm(R);
fprintf('GN consistency: ||J dGN + R||/||R|| = %.4e   ||dGN|| = %.3e\n', ...
        relres, norm(dGN));

sv = svd(full(J));
fprintf('sv(J): max %.3e  min %.3e  cond %.3e  (5 smallest: %s)\n', ...
        sv(1), sv(end), sv(1)/sv(end), sprintf('%.2e ', sv(end-4:end)));

fprintf('distance moved: ||Z - Zseed|| = %.3e  (seed norm %.3e)\n', ...
        norm(Z - Zseed), norm(Zseed));

traj = sms_traj(Z, prob);
fprintf('traj @ iterate: dV = %.4f km/s  mf = %.5f  switches = %d  bang %.1f%%  maxJD = %.3e\n', ...
        traj.dV_kms, traj.mf, traj.switches, 100*traj.bangFrac, traj.maxJointDefect);
save('diag_s1_gateD2.mat', 'Z', 'Zseed', 'out1', 'out2', 'prob', 'info', 'sv', 'relres');
fprintf('diag wall: %.1f min\n', toc(tD)/60);
