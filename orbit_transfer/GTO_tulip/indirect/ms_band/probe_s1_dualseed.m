% PROBE_S1_DUALSEED  One-off: smoke the sms dual seed before the Gate D run.
%
% Builds the 1.12x seed from legacy_ms_f1120.mat at M = 50, prints the
% beta-fit quality, node-1 consistency, one-arc propagation check, and the
% seed residual at eps = 1 (the eps-march start) split into continuity vs
% terminal rows. Cheap early-warning before committing 1-3 h to Gate D.
setup_paths;
tP = tic;
matFile = '../sundman_minfuel/results/minfuel/legacy_ms_f1120.mat';
[Zseed, prob, info] = sms_seed_duals(matFile, 50, 1);
fprintf(['seed: beta = %.5f  spread %.2f%%  burnAgree %.1f%%  coastAgree %.1f%%\n' ...
         '      lamT relStd %.3e  node1Err %.3e  arcCheckErr %.3e\n'], ...
        info.beta, info.spreadPct, 100*info.burnAgree, 100*info.coastAgree, ...
        info.lamTrelStd, info.node1Err, info.arcCheckErr);
R = sms_residual(Zseed, prob);
M = numel(prob.sJ) - 1;
Rc = R(1:16*(M-1));  Rt = R(16*(M-1)+1:end);
fprintf('||R(eps=1)|| = %.3e   (continuity %.3e, terminal %.3e)\n', ...
        norm(R), norm(Rc), norm(Rt));
% worst continuity component classes
Rm = reshape(Rc, 16, M-1);
lbl = {'r','r','r','v','v','v','m','t','lR','lR','lR','lV','lV','lV','lM','lT'};
[wv, wi] = max(abs(Rm), [], 1);
[~, wk] = max(wv);
fprintf('worst joint defect: joint %d, comp %s, |d| = %.3e\n', wk, lbl{wi(wk)}, wv(wk));
fprintf('rowclass max: r,v %.2e  m %.2e  t %.2e  lam %.2e\n', ...
        max(max(abs(Rm(1:6,:)))), max(abs(Rm(7,:))), max(abs(Rm(8,:))), ...
        max(max(abs(Rm(9:16,:)))));
fprintf('probe wall: %.1f min\n', toc(tP)/60);
