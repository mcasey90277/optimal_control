function [rep, ver, infoq] = run_foc_mee(matPath)
% RUN_FOC_MEE  Standard first-order optimality report for a certified MEE row.
%
% Thin wiring driver: warm re-solves a saved certified MEE fuel row with the
% CasADi model attached (refresh_duals_mee.m, returnModel=true), runs the
% generic AD-based FOC/KKT gate (foc_check.m) and the physical MEE-specific
% verifier (verify_pmp_mee.m) side by side on the SAME refreshed point, prints
% the standard FOC report (foc_report.m, with the IPOPT inertia 2nd-order
% verdict attached), and saves a foc_<tag>.mat sidecar next to the row.
%
% INPUTS:
%   matPath - certified row (.mat), e.g. 'results/MEE_M2_10N.mat' [char]
%
% OUTPUTS:
%   rep   - foc_check report struct, with .ipopt attached [struct]
%   ver   - verify_pmp_mee physical-layer report struct [struct]
%   infoq - refresh_duals_mee info struct (drift, ipoptStatus, etc.) [struct]
%
% REFERENCES:
%   [1] verify_common/foc_check.m, foc_report.m, foc_ipopt_inertia.m.
%   [2] verify/refresh_duals_mee.m (dual refresh + certified-quantity guard).
%   [3] verify_common/OPTIMALITY_CERTIFICATION.md Part A (generic FOC gate design).

vcDir = fullfile(fileparts(fileparts(module_root())), 'verify_common');
addpath(vcDir);
setup_verify_common();

[outq, parq, sigq, infoq] = refresh_duals_mee(matPath, struct('returnModel', true));
rep = foc_check(outq, sigq, foc_manifest('earth_mee'), struct());
rep.ipopt = foc_ipopt_inertia(getfield_default(outq, 'regHistory', []));
ver = verify_pmp_mee(outq, parq, sigq, struct('eps', 0));   % physical layer alongside

[~, tag] = fileparts(matPath);
foc_report(rep, tag, fullfile(module_root(), 'results'));
end

% =============================================================================
function v = getfield_default(s, f, dflt)
% GETFIELD_DEFAULT  s.(f) if present and nonempty, else dflt.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
