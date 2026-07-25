function T = run_verify_pmp_all(rowList)
% RUN_VERIFY_PMP_ALL  Sweep first-order PMP verification across every certified
% MEE fuel rung, with refreshed (corrected) duals.
%
% Companion to run_verify_pmp_mee.m, which covers only the two headline rows.
% This driver walks the whole 10 -> 0.1 N ladder so the Campaign-B dual fix
% (2026-07-25) can be judged on BREADTH, not on a single rung: the corrupted
% opti.dual path failed every row it was ever run on, so a credible fix has to
% pass every row it can reach and report honestly on the ones it cannot.
%
% Each row is warm re-solved at its saved primal to re-derive duals
% (refresh_duals_mee.m), then passed to verify_pmp_mee.m. The refresh is gated
% on the CERTIFIED quantities -- Solve_Succeeded, machine-tight defect, and
% final mass unchanged -- not on node-wise drift, because these extremals are
% weak minima and the solver may legitimately slide within the flat optimal
% set; see that file's header. Node drift is reported anyway. Rows whose
% refresh fails are reported as REFRESH-FAIL and excluded -- they are NOT
% counted as verification failures, because their certified primal is
% untouched and only the dual-dependent gates are unavailable.
%
% Rungs are swept shallow-first so the cheap rows report before the deep ones
% (0.2 / 0.1 N carry ~10^4 nodes and can take a long time or hit the known
% MUMPS scale walls).
%
% Each refresh additionally requests the CasADi model (returnModel=true) so
% the generic AD-based first-order gate (verify_common/foc_check.m) can run
% alongside the physical verify_pmp_mee.m check on the SAME refreshed point;
% its standard report prints per row (foc_report.m, with the IPOPT-inertia
% 2nd-order verdict attached) and a foc_<row>.mat sidecar is saved to
% results/. The table's focPass column carries foc_check's advisory verdict
% (report-only burn-in; does not alter certified status).
%
% INPUTS:
%   rowList - optional cellstr of .mat basenames to override the default sweep
% OUTPUTS:
%   T - table: one row per rung with primer/sign/switch-alignment gates, the
%       generic FOC gate's advisory focPass, and the refresh forensics (also
%       saved to results/verify_pmp_all.mat)
%
% REFERENCES:
%   [1] verify/refresh_duals_mee.m (dual refresh + certified-quantity guard).
%   [2] core/casadi_lt_mee.m (the corrected extraction and its evidence).
%   [3] results/dual_anomaly/diag_optidual_minimal.m (minimal reproduction).
resDir = fullfile(module_root(), 'results');
if nargin < 1 || isempty(rowList)
    rowList = {'MEE_M2_10N.mat', 'MEE_M2_5N.mat', 'MEE_M2_2p5N.mat', ...
               'MEE_M2_1N.mat',  'MEE_M2_0p5N.mat', ...
               'MEE_M2_0p2N.mat','MEE_M2_0p1N.mat', ...
               'MEE_M2_1N_PSR_psr_final.mat', 'MEE_M2_0p5N_PSR_psr_final.mat'};
end

vcDir = fullfile(fileparts(fileparts(module_root())), 'verify_common');
addpath(vcDir);
setup_verify_common();
manEarthMee = foc_manifest('earth_mee');

tag = {};  status = {};  primer = [];  signPct = [];  swAlign = [];
flipPct = [];  magDiff = [];  drift = [];  passed = [];  focPass = [];

for q = 1:numel(rowList)
    f = fullfile(resDir, rowList{q});
    [~, base] = fileparts(rowList{q});
    fprintf('\n=== %s ===\n', base);
    if ~isfile(f)
        fprintf('  MISSING -- skipped\n');
        continue
    end
    tag{end+1} = base; %#ok<AGROW>
    try
        [outq, parq, sigq, infoq] = refresh_duals_mee(f, struct('returnModel', true));
    catch ME
        fprintf('  REFRESH-FAIL (%s): %s\n', ME.identifier, ME.message);
        status{end+1} = 'REFRESH-FAIL'; %#ok<AGROW>
        primer(end+1)  = nan;  signPct(end+1) = nan;  swAlign(end+1) = nan; %#ok<AGROW>
        flipPct(end+1) = nan;  magDiff(end+1) = nan;  drift(end+1)   = nan; %#ok<AGROW>
        passed(end+1)  = false; focPass(end+1) = false; %#ok<AGROW>
        continue
    end
    fprintf(['  refresh: %s | defect %.2e | m_f rel change %.2e | node drift ' ...
             '%.2e | dual signs flipped %.1f%% | |dual| rel diff %.2e\n'], ...
        infoq.ipoptStatus, infoq.maxDefect, infoq.mfRelChange, infoq.drift, ...
        100*infoq.maxSignFlipFrac, infoq.magRelDiff);
    v = verify_pmp_mee(outq, parq, sigq, struct('eps', 0));
    if infoq.improved, status{end+1} = 'OK-improved'; else, status{end+1} = 'OK'; end %#ok<AGROW>
    primer(end+1)  = v.primerMedianDeg;        %#ok<AGROW>
    signPct(end+1) = v.overallSignPct;         %#ok<AGROW>
    swAlign(end+1) = v.maxSwitchAlignErr;      %#ok<AGROW>
    flipPct(end+1) = 100*infoq.maxSignFlipFrac;%#ok<AGROW>
    magDiff(end+1) = infoq.magRelDiff;         %#ok<AGROW>
    drift(end+1)   = infoq.mfRelChange;        %#ok<AGROW>
    passed(end+1)  = v.pass;                   %#ok<AGROW>

    repq = foc_check(outq, sigq, manEarthMee, struct());
    if isfield(outq, 'regHistory') && ~isempty(outq.regHistory)
        repq.ipopt = foc_ipopt_inertia(outq.regHistory);
    else
        repq.ipopt = foc_ipopt_inertia([]);
    end
    foc_report(repq, base, resDir);
    focPass(end+1) = repq.pass; %#ok<AGROW>
end

T = table(tag(:), status(:), primer(:), signPct(:), swAlign(:), ...
          flipPct(:), magDiff(:), drift(:), passed(:), focPass(:), ...
    'VariableNames', {'row','status','primerMedDeg','signPct','swAlignErr', ...
                      'dualSignFlipPct','magRelDiff','mfRelChange','pass','focPass'});
fprintf('\n===== PMP verification sweep =====\n');
disp(T);
save(fullfile(resDir, 'verify_pmp_all.mat'), 'T');
fprintf('saved %s\n', fullfile(resDir, 'verify_pmp_all.mat'));
end
