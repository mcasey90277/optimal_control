% P0B_MINTIME_LADDER  Preflight: reproduce the min-time thrust ladder (25 ->
% 200 mN) and SAVE the costate backbone this time (PLAN_PRONG_Z.md P0b).
%
% Re-runs Phase 1 of thrust_continuation_minfuel_indirect.m: march the
% min-time single-shooting solve (solve_tfmin_indirect, complex-step
% Jacobian) UP in thrust from the converged 25 mN reference, warm-starting
% each rung. The old run converged at all 7 rungs but only saved tfMin --
% not the costates. This run saves both: tfMin(T) and costMT(:,T) are the
% Z2/Z3 backbone (top-rung costates seed the P0c energy solve).
%
% Gate: converged (||R|| < 1e-9) at every rung; tfMin matches the old table
%   [6.2907 6.1073 6.0842 6.0842 6.0794 6.0651 6.0651] for f=[1 1.2 1.5 2.1
%   3.2 5 8] to ~1e-3 ND.
%
% Output: results/p0b_mintime_backbone.mat (factors, tfMin, costMT, resNorms)

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

[rv0, rvf, P] = ztl_endpoints();

% converged 25 mN min-time reference (run_gto_tulip_indirect solution)
zGuess25 = [190.476497248065; -79.7064866984696; -0.430399154713168; ...
             0.301159446575878; 0.586671892449694; -0.00711582435720301; ...
             4.32931089137559; 6.29081541876621];

factors = [1, 1.2, 1.5, 2.1, 3.2, 5, 8];          % 25 -> 200 mN (old rung set)
tfMinOld = [6.2907, 6.1073, 6.0842, 6.0842, 6.0794, 6.0651, 6.0651];

nF = numel(factors);
tfMin = nan(1, nF);  costMT = nan(7, nF);  resNorms = nan(1, nF);
zmt = zGuess25;

fprintf('=== P0b: min-time ladder 25 -> 200 mN (saving the backbone) ===\n');
for k = 1:nF
    f = factors(k);  Tmax = f*P.Tmax25;
    [zmt, rn, flag] = solve_tfmin_indirect(rv0, rvf, zmt, Tmax, P.c, P.muStar);
    tfMin(k) = zmt(8);  costMT(:,k) = zmt(1:7);  resNorms(k) = rn;
    fprintf('  f=%4.1f (%5.1f mN): tf_min=%.6f ND (old %.4f, d=%+.1e)  ||R||=%.2e  flag=%d\n', ...
        f, 25*f, zmt(8), tfMinOld(k), zmt(8)-tfMinOld(k), rn, flag);
end

save(fullfile(resDir, 'p0b_mintime_backbone.mat'), ...
     'factors', 'tfMin', 'costMT', 'resNorms', 'zGuess25');
fprintf('saved %s\n', fullfile(resDir, 'p0b_mintime_backbone.mat'));

ok = all(resNorms < 1e-9) && all(abs(tfMin - tfMinOld) < 5e-3);
if ok
    fprintf('GATE P0b: PASS (all rungs converged, tfMin matches old table)\n');
else
    fprintf('GATE P0b: FAIL (res: %s | dtf: %s)\n', ...
        mat2str(resNorms, 2), mat2str(tfMin - tfMinOld, 2));
end
