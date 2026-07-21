function test_ifs_solve()
% TEST_IFS_SOLVE  End-to-end smoke: ifs_solve + ifs_certify run on the full 1.12x
% seed for a few LM iterations without crashing and REDUCE the residual.
%
% This is NOT a convergence gate. The full 1.12x solve does NOT converge (a
% conditioning-limited crawl; see RESULTS.md "Post-merge diagnostic
% investigation"). This smoke only verifies the solve + certify pipeline is
% intact end-to-end (catches integration regressions cheaply) and that a few LM
% steps decrease the residual.
%
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths();
seed = fullfile(here, '..', 'sundman_minfuel', 'results', 'minfuel', 'legacy_ms_f1120.mat');
[Z0, prob, meta] = ifs_seed(seed, struct('mode', 'full'));
out = ifs_solve(Z0, prob, struct('tolR', 1e-8, 'maxIter', 8));   % a few iters only
assert(out.resNorm < out.seedResNorm, ...
       'a few LM steps must reduce the residual, %.3e -> %.3e', out.seedResNorm, out.resNorm);
cert = ifs_certify(out.Z, prob, meta);                            % must run without crashing
assert(isstruct(cert) && isfield(cert, 'ok'), 'certify must return a verdict struct');
fprintf('ALL PASS (smoke: ||R|| %.3e -> %.3e in %d iters; certify ran, ok=%d)\n', ...
        out.seedResNorm, out.resNorm, out.iterations, cert.ok);
end
