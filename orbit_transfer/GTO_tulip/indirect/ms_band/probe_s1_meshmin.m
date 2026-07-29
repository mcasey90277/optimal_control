% PROBE_S1_MESHMIN  One-off: interval-width distribution of the legacy mesh.
setup_paths;
S = load('../sundman_minfuel/results/minfuel/legacy_ms_f1120.mat');
h = diff(S.sigma(:));
fprintf('h: min %.3e  max %.3e  mean %.3e  n(h < mean/10) = %d  n(h < mean/100) = %d\n', ...
        min(h), max(h), mean(h), nnz(h < mean(h)/10), nnz(h < mean(h)/100));
lam8 = S.out.lamDef(8, :);
fprintf('lamT dual adjacent-interval jump: max %.3e  rms %.3e\n', ...
        max(abs(diff(lam8))), rms(diff(lam8)));
