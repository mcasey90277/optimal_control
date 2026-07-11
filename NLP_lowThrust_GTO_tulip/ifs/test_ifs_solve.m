function test_ifs_solve()
% TEST_IFS_SOLVE  Make-or-break gate: full 1.12x IFS solve converges + certifies.
%
% The interior-window gate was dropped (rank-deficient in the lambda_m gauge; see
% RESULTS.md). The full problem's rendezvous transversality pins lambda_m -> full rank.
%
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths();
seed = fullfile(here,'..','sundman_minfuel','results','minfuel','legacy_ms_f1120.mat');
[Z0, prob, meta] = ifs_seed(seed, struct('mode','full'));
out = ifs_solve(Z0, prob, struct('tolR',1e-8,'maxIter',400));
assert(out.success, 'full 1.12x IFS solve must converge, ||R||=%.2e (seed %.2e)', out.resNorm, out.seedResNorm);
cert = ifs_certify(out.Z, prob, meta);
assert(cert.ok, 'converged solution must certify (%s)', cert.text);
fprintf('ALL PASS (k=%d seed||R||=%.2e -> ||R||=%.2e, %d iters, maxSwitchMove=%.2e)\n', ...
        prob.k, out.seedResNorm, out.resNorm, out.iterations, cert.switchMoveFromSeed);
end
