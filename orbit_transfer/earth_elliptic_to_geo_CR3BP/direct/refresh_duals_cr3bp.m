function [out, info] = refresh_duals_cr3bp(S, par, opts)
% REFRESH_DUALS_CR3BP  Re-derive a CR3BP front-door artifact's defect duals
% under the corrected extraction, by warm re-solving at its saved primal.
%
% WHY. run_cr3bp_geo.m artifacts store defectDuals = best.lamDef as extracted
% at solve time, which before 2026-07-25 came from opti.dual(conDef{k}). That
% accessor returns multipliers in CasADi's CANONICALIZED constraint
% orientation rather than the orientation of the opti.g row they pair with in
% grad_f + A'*lam, so the stored duals carry entry-wise sign corruption
% (minimal reproduction: ../../earth_elliptic_to_geo/direct/results/
% dual_anomaly/diag_optidual_minimal.m). This is the single root cause of the
% primer/sign gate failures recorded for BOTH campaigns -- including the
% 10 N / 5 N CR3BP numbers (32.4 deg / 78.4%, 26.7 deg / 76.6%) that task B
% correctly identified as "the same pre-existing anomaly as the 2-body rows,
% not a lunar-specific defect". It was never lunar-specific because it was
% never physics: it was the dual accessor.
%
% core/casadi_lt_mee.m (shared with the 2-body campaign) now sources duals
% from opti.lam_g by row range, so FRESH solves are correct; this helper
% repairs already-banked artifacts without re-running a campaign.
%
% Guarding follows the 2-body sibling refresh_duals_mee.m: gate on the
% CERTIFIED quantities (Solve_Succeeded, machine-tight defect, final mass not
% degraded -- one-sided, since minimum fuel means maximize m_f) rather than on
% node-wise drift, which for these weak minima measures which member of a flat
% optimal set the solver stopped at, not whether the physics changed.
%
% INPUTS:
%   S    - loaded front-door artifact struct (fields X, U, dL, sigma, fp,
%          defectDuals) [struct]
%   par  - kepler_lt_params struct WITH par.pert already attached, exactly as
%          verify_cr3bp_pmp.m rebuilds it from the artifact fingerprint
%   opts - struct (optional): .massTol [1e-6] .feasTol [1e-8]
%          .liftDL [true] .maxIter [fp.maxIter] .printLevel [0]
%          .ipoptExtra [struct()] (MUMPS workspace overrides for deep rungs)
%          .returnModel [false] (true -> casadi_lt_mee ALSO attaches
%          out.model, the live opti/X/U/dL/creg the generic AD-based FOC gate
%          (verify_common/foc_check.m) needs to differentiate the Lagrangian
%          at this SAME refreshed point; false preserves the byte-identical
%          legacy out shape for callers that only want the corrected duals)
%
% OUTPUTS:
%   out  - casadi_lt_mee out-struct with corrected .lamDef [struct]
%   info - struct: .ipoptStatus .maxDefect .drift .mfRelChange .improved
%          .signFlipFrac (fraction of dual entries whose sign changed vs the
%          artifact) .magRelDiff (magnitude agreement -- expected ~1e-16, the
%          signature that ONLY signs were wrong)
%
% REFERENCES:
%   [1] ../../earth_elliptic_to_geo/direct/verify/refresh_duals_mee.m (sibling).
%   [2] ../../earth_elliptic_to_geo/direct/core/casadi_lt_mee.m (the fix).
%   [3] run_cr3bp_geo.m (artifact producer; fp field definitions).
if nargin < 3, opts = struct(); end
d = @(f,v) optdef(opts, f, v);
massTol = d('massTol', 1e-6);
feasTol = d('feasTol', 1e-8);

fp = S.fp;
sigma = S.sigma(:);
sopts = struct('par', par, 'mode', 'fixedtf', 'eps', 0, ...
    'tfTarget', fp.tfTarget, 'x0', S.X(:,1), 'xf', fp.xf(:), ...
    'maxIter', d('maxIter', fp.maxIter), 'warmTight', true, ...
    'printLevel', d('printLevel', 0), 'liftDL', d('liftDL', true), ...
    'ipoptExtra', d('ipoptExtra', struct()), ...
    'returnModel', d('returnModel', false));
out = casadi_lt_mee(sigma, S.X, S.U, S.dL, sopts);

% verify_common bootstrap: certified_guard lives there, and three of this
% function's callers (run_foc_mee, run_verify_pmp_mee, verify_pmp_mee) do not
% put verify_common on the path. Self-bootstrapping here -- the same pattern
% run_foc_tulip / run_foc_elfo already use -- keeps this function callable from
% any of them, and costs a caller that already did it nothing.
if isempty(which('certified_guard'))
    addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'verify_common'));
end

% --- certified-quantity guard (shared: verify_common/certified_guard.m) ------
% Success class + machine-tight defect + ONE-SIDED mass check. Extracted
% 2026-07-26 from four identical copies; behaviour unchanged.
mfSaved = S.X(6, end);
mfNew   = out.X(6, end);
info = certified_guard( ...
    struct('success', out.success, 'ipoptStatus', out.ipoptStatus, ...
           'maxDefect', out.maxDefect, 'value', mfNew), ...
    struct('caller','refresh_duals_cr3bp', 'label', '', 'saved', mfSaved, ...
           'name','m_f', 'errName','mass', 'better','higher', ...
           'feasTol', feasTol, 'tol', massTol));
% Field names this function has always exported (verify_cr3bp_pmp reads
% .mfRelChange by name).
info.drift       = max(abs(out.X(:) - S.X(:)));
info.mfRelChange = info.relChange;

% forensics: signs should change, magnitudes should not
if isequal(size(S.defectDuals), size(out.lamDef))
    info.signFlipFrac = mean(sign(S.defectDuals(:)) ~= sign(out.lamDef(:)));
    info.magRelDiff = max(abs(abs(S.defectDuals(:)) - abs(out.lamDef(:)))) / ...
                      max(max(abs(out.lamDef(:))), 1e-30);
else
    info.signFlipFrac = nan;  info.magRelDiff = nan;
end
end
