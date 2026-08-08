function print_certify_report(rep)
% PRINT_CERTIFY_REPORT  Plain fprintf gate table (value / threshold /
% verdict), one row per gate G1-G5, for a certify_pdg report struct.
%
% Standalone file (not nested in certify_pdg.m) so it is callable both
% from tests/test_certify_nominal.m and the Task-10 front door -- see the
% note at the bottom of certify_pdg.m for why.
%
% INPUTS:  rep - report struct from certify_pdg
% OUTPUTS: (none) -- prints to stdout
fprintf('\n%-30s %14s   %16s   %s\n', 'Gate', 'value', 'threshold', 'verdict');
fprintf('%s\n', repmat('-', 1, 78));
prow('G1 max HS defect',        rep.G1_defect, '< 1e-6',   rep.G1_pass);
prow('G2 pos residual [m]',     rep.G2_pos,    '< 1',      []);
prow('G2 vel residual [m/s]',   rep.G2_vel,    '< 0.1',    []);
prow('G2 mass residual [kg]',   rep.G2_dm,     '< 0.5',    rep.G2_pass);
if isequal(rep.G3_pass, 'skipped')
    fprintf('%-30s %14s   %16s   %s\n', 'G3 cross-method (dmf, dtf)', '--', '--', 'skipped');
else
    prow('G3 |dmf| [kg]',        rep.G3_dmf,       '< 0.1',  []);
    prow('G3 |dtf| [s]',         rep.G3_dtf,       '< 0.2',  rep.G3_pass);
    prow('G3 traj Linf [m]',     rep.G3_traj_Linf, '(info)', []);
end
if isequal(rep.G4_pass, 'skipped')
    fprintf('%-30s %14s   %16s   %s\n', 'G4 lossless gap', '--', '--', 'skipped');
else
    prow('G4 lossless gap [m/s^2]', rep.G4_gap, '< 1e-4*Tmax/m0', rep.G4_pass);
end
prow('G5 bound fraction',       rep.G5_bound_frac, '>= 0.95', []);
prow('G5 interior switches',    rep.G5_switches,   '<= 2',    []);
prow('G5 primer angle [deg]',   rep.G5_primer_deg, '< 1',     rep.G5_pass);
fprintf('%s\n', repmat('-', 1, 78));
fprintf('%-30s %14s   %16s   %s\n\n', 'ALL GATES', '', '', bool2str(rep.all_pass));
end

function prow(name, val, thresh, verdict)
fprintf('%-30s %14.6g   %16s   %s\n', name, val, thresh, bool2str(verdict));
end

function s = bool2str(b)
if isempty(b),          s = '';
elseif isequal(b,true), s = 'PASS';
else,                   s = 'FAIL';
end
end
