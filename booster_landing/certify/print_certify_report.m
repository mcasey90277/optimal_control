function print_certify_report(rep)
% PRINT_CERTIFY_REPORT  Plain fprintf gate table (value / EFFECTIVE
% threshold / verdict) for a certify_pdg report struct, with a per-row
% verdict on every row (not just the last row of a multi-row gate) and an
% explicit gate-summary line for every multi-row gate.
%
% Standalone file (not nested in certify_pdg.m) so it is callable both
% from tests/test_certify_nominal.m and the Task-10 front door -- see the
% note at the bottom of certify_pdg.m for why.
%
% Thresholds shown are EFFECTIVE (i.e. nominal*tolScale for the gates
% tolScale actually scales -- currently only G3), read from rep.tolScale
% so the printed table can never say something false like "value < 0.1
% PASS" when the enforced threshold was actually 0.1*tolScale. The header
% echoes tolScale whenever it is not 1, and every row's own individual
% verdict is computed here from (value, effective threshold), not just
% copied from a single stored .*_pass flag -- so a table with several rows
% under one gate (G2, G3, G5) never shows a bare data row next to another
% row's verdict.
%
% INPUTS:  rep - report struct from certify_pdg (must include .tolScale;
%                falls back to 1 if an older-format struct lacks it)
% OUTPUTS: (none) -- prints to stdout
if isfield(rep, 'tolScale'), tolScale = rep.tolScale; else, tolScale = 1; end

fprintf('\n%-30s %14s   %20s   %s\n', 'Gate', 'value', 'threshold', 'verdict');
if tolScale ~= 1
    fprintf('(tolScale = %g, applied to G3 only -- G2/G5 thresholds below are the real, unscaled nominal gates)\n', tolScale);
end
fprintf('%s\n', repmat('-', 1, 84));

% G1 (single row, no tolScale):
prow('G1 max HS defect', rep.G1_defect, 1e-6, '<', 1);

% G2 (three rows, NOT scaled by tolScale -- see certify_pdg.m header):
prow('G2 pos residual [m]',   rep.G2_pos, 1,   '<', 1);
prow('G2 vel residual [m/s]', rep.G2_vel, 0.1, '<', 1);
prow('G2 mass residual [kg]', rep.G2_dm,  0.5, '<', 1);
gsummary('G2', rep.G2_pass);

% G3 (two data rows + one info row, dmf/dtf scaled by tolScale):
if isequal(rep.G3_pass, 'skipped')
    fprintf('%-30s %14s   %20s   %s\n', 'G3 cross-method (dmf, dtf)', '--', '--', 'skipped');
else
    prow('G3 |dmf| [kg]', rep.G3_dmf, 0.1, '<', tolScale);
    prow('G3 |dtf| [s]',  rep.G3_dtf, 0.2, '<', tolScale);
    fprintf('%-30s %14.6g   %20s   %s\n', 'G3 traj Linf [m]', rep.G3_traj_Linf, '(info)', '');
    gsummary('G3', rep.G3_pass);
end

% G4 (single row, no tolScale; threshold is a formula of P, not a bare
% number known to this printer, so its verdict comes straight from
% rep.G4_pass -- computed by certify_pdg against the true 1e-4*Tmax/m0
% value -- rather than being recomputed here from a guessed number):
if isequal(rep.G4_pass, 'skipped')
    fprintf('%-30s %14s   %20s   %s\n', 'G4 lossless gap', '--', '--', 'skipped');
else
    fprintf('%-30s %14.6g   %20s   %s\n', 'G4 lossless gap [m/s^2]', rep.G4_gap, ...
        '< 1e-4*Tmax/m0', bool2str(rep.G4_pass));
end

% G5 (three data rows + structure/primer verdicts, no tolScale):
prow('G5 bound fraction',    rep.G5_bound_frac, 0.95, '>=', 1);
prow('G5 interior switches', rep.G5_switches,   2,    '<=', 1);
fprintf('%-30s %14s   %20s   %s\n', 'G5 structure (bound+switch+max-last)', '', '', bool2str(rep.G5_structOk));
prow('G5 primer angle [deg]', rep.G5_primer_deg, 1, '<', 1);
gsummary('G5', rep.G5_pass);

fprintf('%s\n', repmat('-', 1, 84));
fprintf('%-30s %14s   %20s   %s\n\n', 'ALL GATES', '', '', bool2str(rep.all_pass));
end

function prow(name, val, nominalThresh, op, scale)
% Print one data row with its OWN verdict, computed from (val, nominal
% threshold * scale) -- never borrowed from a neighboring row.
effThresh = nominalThresh * scale;
if scale ~= 1
    threshStr = sprintf('%s %g (x%g=%.4g)', op, nominalThresh, scale, effThresh);
else
    threshStr = sprintf('%s %g', op, nominalThresh);
end
switch op
    case '<',  pass = val < effThresh;
    case '>=', pass = val >= effThresh;
    case '<=', pass = val <= effThresh;
end
fprintf('%-30s %14.6g   %20s   %s\n', name, val, threshStr, bool2str(pass));
end

function gsummary(gname, gatePass)
fprintf('%-30s %14s   %20s   %s\n', sprintf('  -> %s gate', gname), '', '', bool2str(gatePass));
end

function s = bool2str(b)
if isempty(b),          s = '';
elseif isequal(b,true), s = 'PASS';
else,                   s = 'FAIL';
end
end
