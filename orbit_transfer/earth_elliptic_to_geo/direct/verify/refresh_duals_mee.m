function [out, par, sigma, info] = refresh_duals_mee(matPath, opts)
% REFRESH_DUALS_MEE  Warm re-solve a saved certified MEE fuel row to recover
% its defect duals under the CORRECTED extraction, and return the (out, par,
% sigma) triple verify_pmp_mee.m consumes.
%
% WHY THIS EXISTS. The certified .mat caches were written before the
% Campaign-B dual bug was found (2026-07-25): their stored res.fuel.lamDef
% came from opti.dual(conDef{k}), which returns multipliers in CasADi's
% CANONICALIZED constraint orientation rather than the orientation of the
% opti.g row they pair with in grad_f + A'*lam -- an entry-wise sign
% corruption (minimal reproduction: results/dual_anomaly/
% diag_optidual_minimal.m). casadi_lt_mee.m now sources lamDef from
% opti.lam_g by row range, but that only fixes FRESH solves; the caches still
% carry the stale duals. This helper re-derives them without re-running any
% campaign: it rebuilds the NLP at the saved primal, warm re-solves (~0
% iterations at a converged point), and hands back the corrected duals.
%
% WHAT "SAME TRAJECTORY" MEANS HERE, AND WHY IT IS NOT NODE-WISE DRIFT. The
% obvious guard -- demand max|X_resolved - X_saved| be tiny -- is the WRONG
% test for this problem, and empirically it rejects almost every rung (5 N
% drifts 6.7e-2, the 1 N PSR row 5.8e-1) while 10 N reproduces to 3.6e-15.
% That is not instability: the campaign's own SOSC study found these min-fuel
% bang-bang extremals are WEAK, non-strict minima (270 numerically-flat
% reduced-Hessian directions at 10 N, process/DESIGN_sosc.md sec 11-12). Along
% a flat direction the solver may slide to a different point of the same
% optimal set at identical objective, so node-wise drift measures which member
% of the flat set IPOPT happened to stop at -- not whether the physics
% changed. Gating on it would make dual refresh impossible for exactly the
% rows that need it.
%
% The guard is therefore on the CERTIFIED QUANTITIES instead: the re-solve
% must report Solve_Succeeded, land machine-tight (maxDefect < feasTol), and
% reproduce the row's certified final mass to massTol. Node-wise drift is
% always measured and returned in info.drift, and the caller/driver prints it,
% so a large drift is visible and reportable -- it is just not, by itself,
% disqualifying. The duals returned are the multipliers of the point actually
% re-solved, and verify_pmp_mee is handed THAT point's X/U alongside them, so
% the verification is self-consistent by construction.
%
% INPUTS:
%   matPath - certified row (.mat): results/MEE_M2_*.mat or *_PSR_psr_final.mat
%   opts    - struct (optional):
%             .massTol   max relative |m_f change| accepted   [default 1e-6]
%             .feasTol   max accepted re-solve maxDefect      [default 1e-8]
%             .maxIter   IPOPT cap override      [default: the row's own]
%             .printLevel                         [default 0]
%             .returnModel  pass-through to casadi_lt_mee's returnModel opt --
%                           when true, out.model.{opti,creg} is attached so
%                           callers (run_foc_mee, run_verify_pmp_all) can feed
%                           out straight into foc_check          [default false]
%
% OUTPUTS:
%   out   - casadi_lt_mee result struct with the CORRECTED .lamDef [struct]
%   par   - kepler_lt_params for the row [struct]
%   sigma - node grid [(N+1)x1]
%   info  - struct: .drift .ipoptStatus .maxDefect .tag .lamDefStale (the
%           cache's original duals, for forensics) .maxSignFlipFrac (fraction
%           of defect-dual entries whose sign changed vs the cache)
%
% REFERENCES:
%   [1] core/casadi_lt_mee.m (the corrected extraction; its comment carries
%       the root-cause evidence).
%   [2] verify/sosc/sosc_load_row.m (the row-shape normalizer reused here).
%   [3] results/dual_anomaly/diag_rawdual.m (the fix test: 32.370 -> 0.000 deg).
if nargin < 2, opts = struct(); end
d = @(f,v) optdef(opts, f, v);
massTol = d('massTol', 1e-6);
feasTol = d('feasTol', 1e-8);
prnt    = d('printLevel', 0);
retModel = d('returnModel', false);

saved = sosc_load_row(matPath);
par   = kepler_lt_params(saved.thrustN, saved.m0kg, saved.ispS);
sigma = saved.sigma(:);
mIter = d('maxIter', saved.maxIter);

sopts = struct('par', par, 'mode', 'fixedtf', 'eps', 0, ...
    'tfTarget', saved.tfTarget, 'x0', saved.X(:,1), 'xf', saved.xf, ...
    'maxIter', mIter, 'warmTight', true, 'printLevel', prnt, ...
    'returnModel', retModel);
out = casadi_lt_mee(sigma, saved.X, saved.U, saved.dL, sopts);

% verify_common bootstrap: certified_guard lives there, and three of this
% function's callers (run_foc_mee, run_verify_pmp_mee, verify_pmp_mee) do not
% put verify_common on the path. Self-bootstrapping here -- the same pattern
% run_foc_tulip / run_foc_elfo already use -- keeps this function callable from
% any of them, and costs a caller that already did it nothing.
if isempty(which('certified_guard'))
    addpath(fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))), 'verify_common'));
end

% --- certified-quantity guard (shared: verify_common/certified_guard.m) ------
% Success class + machine-tight defect + the ONE-SIDED mass check whose
% convention originated in this campaign (reproduce/verify_row.m: minimum fuel
% means MAXIMIZE final mass, so a higher mass is an improvement, never a
% failure). Extracted 2026-07-26; behaviour unchanged -- casadi_lt_mee already
% folds the IPOPT status whitelist into out.success, so applying it again in
% the shared guard is a no-op here.
mfSaved = saved.X(6, end);
info = certified_guard( ...
    struct('success', out.success, 'ipoptStatus', out.ipoptStatus, ...
           'maxDefect', out.maxDefect, 'value', out.X(6,end)), ...
    struct('caller','refresh_duals_mee', 'label', saved.tag, 'saved', mfSaved, ...
           'name','m_f', 'errName','mass', 'better','higher', ...
           'feasTol', feasTol, 'tol', massTol));
% Field names this function has always exported (run_verify_pmp_all and
% run_verify_pmp_mee read .mfRelChange and .drift by name).
info.drift       = max(abs(out.X(:) - saved.X(:)));
info.tag         = saved.tag;
info.mfRelChange = info.relChange;
info.mfSaved     = mfSaved;
info.mfResolved  = out.X(6,end);

% forensics: how much did the corrected extraction actually change? Pull the
% cache's own (stale, opti.dual-sourced) duals straight from the .mat -- they
% are not part of sosc_load_row's normalized shape.
Sraw = load(matPath);
if isfield(Sraw, 'res') && isfield(Sraw.res.fuel, 'lamDef')
    info.lamDefStale = Sraw.res.fuel.lamDef;
elseif isfield(Sraw, 'out') && isfield(Sraw.out.finalOut, 'lamDef')
    info.lamDefStale = Sraw.out.finalOut.lamDef;
else
    info.lamDefStale = [];
end
if ~isempty(info.lamDefStale) && isequal(size(info.lamDefStale), size(out.lamDef))
    sgFlip = sign(info.lamDefStale) ~= sign(out.lamDef);
    info.maxSignFlipFrac = mean(sgFlip(:));
    info.magRelDiff = max(abs(abs(info.lamDefStale(:)) - abs(out.lamDef(:)))) / ...
                      max(max(abs(out.lamDef(:))), 1e-30);
else
    info.maxSignFlipFrac = nan;  info.magRelDiff = nan;
end
end
