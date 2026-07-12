% TEST_IFS_SEED_ADJOINT  Rung-A gates GA1/GA2 for the adjoint-smoother seed.
%
% Method under test: 'smooth' (adjoint smoother -- terminal costate fitted to
% the whole dual history with truncated-SVD GN; blended node costates). The
% pure 'sweep' method was falsified 2026-07-12 (backward amplification 1.2e12;
% see RESULTS_RUNG01_RUNG2.md) and is kept only as a diagnostic.
%
% GA1 (fit sanity): fit descends; S sign pattern reproduces the dual-map arcs;
%   q ~= 1 at switches; blend fraction reported.
% GA2 (seed residual): ifs_residual at the smoother seed vs the dual-map
%   baseline (1.96 at 1.12x). GA3 (full ifs_solve2) is a separate long job.
%
% REFERENCES: PLAN_OF_ATTACK_2.md Rung A (gates), ifs_seed_adjoint.m

here = fileparts(mfilename('fullpath'));
cd(here); setup_paths();
matFile = fullfile(here, '..', 'sundman_minfuel', 'results', 'minfuel', ...
                   'legacy_ms_f1120.mat');

fprintf('=== Rung A smoke: adjoint-SMOOTHER seed on 1.12x ===\n');

t0 = tic;
[Za, prA, ma] = ifs_seed_adjoint(matFile, struct('method','smooth','verbose',true));
tA = toc(t0);

fprintf('\n--- GA1: fit sanity ---\n');
fprintf('k (dual-S crossings)    : %d  (raw throttle count: %d)\n', ma.k, ma.kThrottle);
fprintf('fit ||R||               : %.4e -> %.4e (%d evals/sweeps)\n', ...
        ma.fitResHist(1), ma.fitResHist(end), ma.nSweeps);
fprintf('data-scale s            : %.6e (beta0 %.6e)\n', ma.sFit, ma.beta0);
fprintf('S-sign agreement        : %.2f %%\n', ma.signAgree);
fprintf('S zero crossings (seed) : %d (vs k=%d)\n', ma.nCrossSweep, ma.k);
fprintf('q at switches           : mean=%.4f  spread=[%.4f, %.4f]\n', ...
        mean(ma.qSw), min(ma.qSw), max(ma.qSw));
fprintf('misfit profile          : median=%.3e  max=%.3e\n', ...
        median(ma.misfit), max(ma.misfit));
fprintf('blend                   : %.1f %% of nodes trust the sweep\n', 100*ma.fracSweep);
fprintf('growth (post-blend)     : max||lam||/||lam(tauf)|| = %.3e at tau=%.3f\n', ...
        ma.growth, ma.growthTau);
fprintf('terminal-dual source    : %s\n', ma.dualSrc);

fprintf('\n--- GA2: seed residual (target: beat dual-map 1.96) ---\n');
fprintf('smoother seedRes        : %.6e\n', ma.seedResNorm);
fprintf('build time              : %.1f s\n', tA);

% dual-map baseline on the same file, same tauParam
[~, ~, md] = ifs_seed(matFile, struct('mode', 'full'));
fprintf('dual-map seedRes        : %.6e\n', md.seedResNorm);
fprintf('improvement factor      : %.2fx\n', md.seedResNorm/ma.seedResNorm);

assert(all(isfinite(Za)), 'GA1 FAIL: non-finite seed vector');
assert(ma.signAgree > 90, 'GA1 FAIL: S-sign agreement %.1f%% < 90%%', ma.signAgree);
assert(ma.seedResNorm < md.seedResNorm, ...
       'GA2 FAIL: smoother seed (%.3e) does not beat dual map (%.3e)', ...
       ma.seedResNorm, md.seedResNorm);
fprintf('\nPASS: GA1 sane, GA2 beats the dual map.\n');
