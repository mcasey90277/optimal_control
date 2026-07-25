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
resDir  = fullfile(module_root(), 'results');

% DUAL REFRESH (2026-07-25). The caches store lamDef as it was extracted at
% solve time, via opti.dual(conDef{k}) -- the entry-wise sign-corrupted path
% that WAS the Campaign-B primer anomaly (root cause + minimal reproduction:
% results/dual_anomaly/diag_optidual_minimal.m). casadi_lt_mee.m now sources
% duals from opti.lam_g, but that fixes fresh solves only, so each row is
% warm re-solved here at its saved primal to re-derive correct duals
% (refresh_duals_mee.m, which refuses to proceed if the primal drifts). The
% certified PRIMAL numbers -- mass, switches, defect -- are untouched by any
% of this; only the dual-dependent PMP gates change.
rows = { 'MEE_M2_10N.mat',              'MEE_M2_10N',     '10 N M2'
         'MEE_M2_1N_PSR_psr_final.mat', 'MEE_M2_1N_PSR',  '1 N PSR' };
verOut = cell(size(rows,1), 1);

for q = 1:size(rows,1)
    matPath = fullfile(resDir, rows{q,1});
    tag     = rows{q,2};
    fprintf('\n=== verify_pmp_mee: %s (results/%s) ===\n', rows{q,3}, rows{q,1});
    try
        [outq, parq, sigmaq, infoq] = refresh_duals_mee(matPath);
    catch ME
        fprintf(['  DUAL REFRESH FAILED (%s): %s\n  Row NOT verified -- the ' ...
                 'certified primal is unaffected.\n'], ME.identifier, ME.message);
        continue
    end
    fprintf(['  dual refresh: %s | defect %.2e | m_f rel change %.2e | node ' ...
             'drift %.2e | dual signs flipped %.1f%% | |dual| rel diff %.2e\n'], ...
        infoq.ipoptStatus, infoq.maxDefect, infoq.mfRelChange, infoq.drift, ...
        100*infoq.maxSignFlipFrac, infoq.magRelDiff);
    verq = verify_pmp_mee(outq, parq, sigmaq, struct('eps', 0));
    verOut{q} = verq;
    save(fullfile(resDir, sprintf('verify_pmp_mee_%s.mat', tag)), 'verq', 'infoq');
    fig_switching(verq, outq, tag, resDir);
end

fprintf('\nrun_verify_pmp_mee: done.\n');
