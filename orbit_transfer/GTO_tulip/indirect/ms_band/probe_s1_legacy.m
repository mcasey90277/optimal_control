% PROBE_S1_LEGACY  One-off: inspect legacy_ms_f1120.mat for the sms seed builder.
setup_paths;
S = load('../sundman_minfuel/results/minfuel/legacy_ms_f1120.mat');
fprintf('fields: %s\n', strjoin(fieldnames(S), ', '));
if isfield(S, 'sigma')
    sg = S.sigma(:);
    fprintf('sigma: N=%d  range [%.3g %.3g]  uniformity err=%.3e\n', ...
            numel(sg)-1, sg(1), sg(end), max(abs(diff(sg) - mean(diff(sg)))));
end
if isfield(S, 'tauf0'), fprintf('tauf0 = %.8f\n', S.tauf0); end
if isfield(S, 'pSund'), fprintf('pSund = %.3f\n', S.pSund); end
fprintf('factor = %.4f\n', S.factor);
fprintf('size X = %s, size U = %s, size lamDef = %s\n', mat2str(size(S.out.X)), ...
        mat2str(size(S.out.U)), mat2str(size(S.out.lamDef)));
fprintf('t(end) = %.8f  (factor*tfMin = %.8f)\n', S.out.X(8,end), S.factor*6.290694);
lam8 = S.out.lamDef(8, :);
fprintf('lamDef row 8 (lamT dual): mean %.6e  std %.6e  (rel %.3e)\n', ...
        mean(lam8), std(lam8), std(lam8)/abs(mean(lam8)));
fprintf('mf = %.6f  switches = %d\n', S.out.X(7,end), S.out.switches);
% sigf reconstruction cross-check: cumtrapz of 1/kappa along the solution
r1 = sqrt(sum((S.out.X(1:3,:) - [-0.012150585609624; 0; 0]).^2, 1));
kap = r1.^1.5;
sigRec = cumtrapz(S.out.X(8,:).', 1./kap(:));
fprintf('sigf reconstructed = %.6f\n', sigRec(end));
