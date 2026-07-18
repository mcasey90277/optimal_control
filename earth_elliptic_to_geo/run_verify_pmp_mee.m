% RUN_VERIFY_PMP_MEE  Reproduction driver for the Task-10 verify_pmp_mee.m
% gate numbers (MEE/sigma-domain PMP verifier).
%
% Thin script, no new logic: loads the two certified MEE fuel solutions
% (10 N M2 and 1 N PSR), reconstructs the (out, par, sigma) triple each needs
% from what was actually saved (res.fuel/res.sigma/res.fp for the 10 N case;
% out.finalOut/out.finalSigma/fpFinal for the PSR case -- see
% results/tmp_inspect_mats*.m probes, not committed), calls verify_pmp_mee.m
% on each, prints both gate tables, saves the returned ver structs +
% fig_switching.m figures to results/ -- so the numbers quoted in
% .superpowers/sdd/task-10-report.md are reproducible from a committed
% script, per this campaign's every-experiment-writes-a-script convention.
%
% INPUTS:  none (loads results/MEE_M2_10N.mat and
%          results/MEE_M2_1N_PSR_psr_final.mat, both already committed
%          certified caches -- no NLP solve happens here)
% OUTPUTS: none (prints both verify_pmp_mee gate tables to stdout; writes
%          results/verify_pmp_mee_10N.mat, results/verify_pmp_mee_1N.mat
%          [ver structs] and results/MEE_M2_10N_fig_switching.png,
%          results/MEE_M2_1N_PSR_fig_switching.png)
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/verify_pmp_mee.m (verifier this drives).
%   [2] earth_elliptic_to_geo/fig_switching.m (Fig-16 analog figure).
%   [3] earth_elliptic_to_geo/kepler_lt_params.m (par reconstruction from
%       thrustN/m0kg/ispS, the only inputs verify_pmp_mee.m needs beyond
%       what's in the .mat files).
%   [4] .superpowers/sdd/task-10-report.md (gate numbers this reproduces).
here    = fileparts(mfilename('fullpath'));
resDir  = fullfile(here, 'results');

%% --- 10 N M2 case (results/MEE_M2_10N.mat: struct 'res') -------------------
S1  = load(fullfile(resDir, 'MEE_M2_10N.mat'));
res = S1.res;
out10   = res.fuel;
sigma10 = res.sigma;
par10   = kepler_lt_params(res.fp.thrustN, res.fp.m0kg, res.fp.ispS);

fprintf('\n=== verify_pmp_mee: 10 N M2 (results/MEE_M2_10N.mat) ===\n');
ver10 = verify_pmp_mee(out10, par10, sigma10, struct('eps', 0));
save(fullfile(resDir, 'verify_pmp_mee_10N.mat'), 'ver10');
fig_switching(ver10, out10, 'MEE_M2_10N', resDir);

%% --- 1 N PSR case (results/MEE_M2_1N_PSR_psr_final.mat: struct 'out') ------
S2  = load(fullfile(resDir, 'MEE_M2_1N_PSR_psr_final.mat'));
out1    = S2.out.finalOut;
sigma1  = S2.out.finalSigma;
par1    = kepler_lt_params(S2.fpFinal.thrustN, S2.fpFinal.m0kg, S2.fpFinal.ispS);

fprintf('\n=== verify_pmp_mee: 1 N PSR (results/MEE_M2_1N_PSR_psr_final.mat) ===\n');
ver1 = verify_pmp_mee(out1, par1, sigma1, struct('eps', 0));
save(fullfile(resDir, 'verify_pmp_mee_1N.mat'), 'ver1');
fig_switching(ver1, out1, 'MEE_M2_1N_PSR', resDir);

fprintf('\nrun_verify_pmp_mee: done.\n');
