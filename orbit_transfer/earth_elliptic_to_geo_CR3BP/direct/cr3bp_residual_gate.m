function R = cr3bp_residual_gate(matFile, opts)
% CR3BP_RESIDUAL_GATE  Continuous-residual (G1) gate for a certified
%   CR3BP-GEO front-door artifact -- the campaign's first rung-2 check.
%
%   Thin wiring driver: load the artifact, rebuild par + par.pert from its
%   OWN fingerprint (the verify_cr3bp_pmp pattern -- read, never
%   recomputed, so the gate cannot silently diverge from the solved
%   physics), warm re-solve at the saved primal for the trajectory
%   (refresh_duals_cr3bp), then measure the true per-interval error with
%   the shared gate verify_common/mee_residual. The lunar third-body term
%   rides in through par.pert; the gate itself never looks.
%
% INPUTS:
%   matFile - front-door artifact, e.g. 'results/cr3bp_T10N_phi0_fuel.mat'
%             (default). Must be lunar-aware (fp.gain > 0). [char]
%   opts    - (optional) forwarded to mee_residual (.relTol .absTol)
%
% OUTPUTS:
%   R - mee_residual output struct, plus .matFile .fp .refreshInfo
%
% REFERENCES:
%   [1] verify_common/mee_residual.m (the shared gate; oc.local_residual)
%   [2] verify_cr3bp_pmp.m (the fingerprint-rebuild pattern copied here)
%   [3] DRO_tulip/FINDINGS.md (why defects are not accuracy: the 1e7 gap)

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(matFile)
    matFile = fullfile(here, 'results', 'cr3bp_T10N_phi0_fuel.mat');
end
if nargin < 2, opts = struct(); end
run(fullfile(here, 'setup_paths.m'));
if isempty(which('mee_residual'))
    addpath(fullfile(fileparts(fileparts(here)), 'verify_common'));
end
if isempty(which('oc.local_residual'))
    addpath(fullfile(fileparts(fileparts(fileparts(here))), 'oclib'));
end

S  = load(matFile);
fp = S.fp;
assert(isfield(fp, 'gain') && fp.gain > 0, 'cr3bp_residual_gate:notLunar', ...
    'artifact %s is not lunar-aware (gain=0): the 2-body campaign covers it', matFile);

par = kepler_lt_params(fp.thrustN, fp.m0kg, fp.ispS);
par.pert = struct('muM', fp.muM, 'DM', fp.DM, 'nM', fp.nM, ...
                  'phi0', fp.phi0, 'gain', fp.gain);

[out, info] = refresh_duals_cr3bp(S, par, struct());

R = mee_residual(out, par, S.sigma, opts);
R.matFile = matFile;  R.fp = fp;  R.refreshInfo = info;

LU = par.LU_km;  TU = par.TU_s;
fprintf('=== CR3BP-GEO continuous-residual gate: %g N, phi0 = %.3f ===\n', ...
        fp.thrustN, fp.phi0);
fprintf('re-solve: drift %.2e, defect %.2e, %s\n', info.drift, info.maxDefect, info.ipoptStatus);
fprintf('P    : max %.3e ND = %.4f km   (med %.3e ND)\n', R.RPMax, R.RPMax*LU, R.RPMed);
fprintf('shape: max %.3e   mass: max %.3e (%.4f kg)   time: max %.3e ND = %.2f s\n', ...
        R.ReMax, R.RmMax, R.RmMax*par.m0kg, R.RtMax, R.RtMax*TU);
fprintf('worst interval #%d of %d at L = %.2f rad (rev %.1f)\n', ...
        R.kWorst, numel(R.RP), R.LMid(R.kWorst), (R.LMid(R.kWorst)-pi)/(2*pi));
end
