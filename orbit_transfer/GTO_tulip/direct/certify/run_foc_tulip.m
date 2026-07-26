function [rep, cert, info] = run_foc_tulip(matPath)
% RUN_FOC_TULIP  Standard first-order optimality report for a certified tulip
% min-fuel row, with an independent LS-vs-dual costate cross-check.
%
% Thin wiring driver (Task 8, tulip precinct of the FOC-gate campaign):
% loads a certified Sundman min-fuel row, warm re-solves it at eps=0 with the
% CasADi model attached (casadi_minfuel_sundman.m, opts.returnModel=true),
% guards on certified quantities, runs the generic AD-based FOC/KKT gate
% (foc_check.m) with the IPOPT-inertia 2nd-order verdict attached
% (foc_ipopt_inertia.m), prints the standard report (foc_report.m), and THEN
% runs the pre-existing physical-layer certifier (certify_minfuel_pmp.m) on
% the SAME re-solved row so the two costate sources -- foc_check's raw KKT
% duals (opti.lam_g) vs certify_minfuel_pmp's least-squares-reconstructed
% continuous adjoint -- can be compared side by side. Two independent costate
% sources agreeing is the point of the cross-check; both numbers are recorded
% (and printed) whether they agree or not.
%
% INPUTS:
%   matPath - certified row (.mat), e.g. 'sundman_minfuel_certified.mat'
%             [char, default: that file beside this one -- the certified
%             25-switch, 1.15x-min-time flagship, README.md "Certified
%             result"]. Expected fields: out.X [8xnN], out.U [4xnN], sigma
%             [nNx1], tauf0, pSund, rv0, rvf (sundman_homotopy.m save list).
%
% OUTPUTS:
%   rep  - foc_check report struct, with .ipopt (2nd-order verdict) and
%          .crossCheck (LS-reconstructed vs raw-dual costate comparison)
%          attached [struct]
%   cert - certify_minfuel_pmp report struct on the SAME re-solved row
%          (the LS-reconstructed-costate side of the cross-check) [struct]
%   info - warm re-solve bookkeeping: .drift .ipoptStatus .maxDefect
%          .mfSaved .mfResolved .improved .crossCheck [struct]
%
% REFERENCES:
%   [1] verify_common/foc_check.m, foc_report.m, foc_ipopt_inertia.m.
%   [2] casadi_minfuel_sundman.m (Task 8 returnModel/creg + regHistory hook).
%   [3] certify_minfuel_pmp.m (independent LS-reconstructed-costate PMP
%       certifier -- the cross-check partner).
%   [4] earth_elliptic_to_geo/direct/verify/refresh_duals_mee.m header (the
%       certified-quantity-not-node-drift guard rationale this mirrors: these
%       min-fuel bang-bang extremals are weak/non-strict minima, so gating on
%       node-wise drift would wrongly reject a legitimate re-solve).
%   [5] sundman_homotopy.m (the artifact's own save-field list this loader
%       assumes: out,sigma,rv0,rvf,tauf0,pSund,eps,tbl,fp).

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(matPath)
    matPath = fullfile(here, '..', 'lib', 'sundman_minfuel_certified.mat');
end
addpath(here);
% cr3bp_common (cr3bp_lt_params) only -- pumpkyn/gto_tulip_endpoints not
% needed here since rv0/rvf come straight from the artifact (same addpath
% this folder's own sundman_homotopy.m uses).
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
vcDir = fullfile(here, '..', '..', '..', 'verify_common');
addpath(vcDir);
setup_verify_common();

S = load(matPath);
p = cr3bp_lt_params(0.025, 15, 2100);   % matches the certified campaign (README.md)
pSund = 1.5;   if isfield(S,'pSund') && ~isempty(S.pSund), pSund = S.pSund; end
tauf0 = S.tauf0;
tf    = S.out.X(8,end);

% --- warm re-solve AT the saved primal, with the FOC-gate model attached ----
maxIter = 800;
sopts = struct('returnModel', true);
out = casadi_minfuel_sundman(S.sigma, tf, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, ...
    S.out.X, S.out.U, tauf0, pSund, maxIter, 0, true, sopts);

% --- certified-quantity guard (shared: verify_common/certified_guard.m) ------
% Success class + machine-tight defect + a ONE-SIDED final-mass check. Was ~25
% lines here; the logic and thresholds were identical in four drivers across
% four campaigns, so it now lives in one place with a branch-level test.
mfSaved = S.out.mf;
info = certified_guard( ...
    struct('success', out.success, 'ipoptStatus', out.ipoptStatus, ...
           'maxDefect', out.maxDefect, 'value', out.mf), ...
    struct('caller','run_foc_tulip', 'label', matPath, 'saved', mfSaved, ...
           'name','m_f', 'errName','mass', 'better','higher', ...
           'feasTol', 1e-8, 'tol', 1e-6));
% Field names this driver has always exported, kept for its callers.
info.drift      = max(abs(out.X(:) - S.out.X(:)));
info.mfSaved    = mfSaved;
info.mfResolved = out.mf;

% --- generic AD-based FOC/KKT gate (raw-dual costate source) ----------------
man = foc_manifest('tulip');
rep = foc_check(out, S.sigma, man, struct());
rep.ipopt = foc_ipopt_inertia(getfield_default(out, 'regHistory', []));

% --- independent cross-check: LS-reconstructed continuous costate ----------
% Same re-solved row (out), same NLP solution -- an independent recovery
% pipeline (block-bidiagonal adjoint recursion + primer-direction LS, NOT the
% NLP's own KKT duals) rather than a duplicate of foc_check's raw-dual source.
solStruct = struct('out', out, 'sigma', S.sigma, 'tauf0', tauf0, 'pSund', pSund, ...
                    'eps', 0, 'rv0', S.rv0, 'rvf', S.rvf);
cert = certify_minfuel_pmp(solStruct, false);

% Compare the two sources on their OWN PMP-direction/sign-law terms, not the
% two drivers' full advisory verdicts (rep.pass also folds in unrelated
% checks -- e.g. sdotMinRel -- that have nothing to do with costate quality
% and would make an apples-to-oranges "agreement"). focCostateOK mirrors
% cert.passed's own primer/sign-law thresholds (certify_minfuel_pmp.m: primer
% dir err < 1e-2, sign match > 99%) applied to the raw-dual side.
focCostateOK  = (out.primerAlignDeg < 1) && (rep.signPct >= 99);
certCostateOK = cert.passed;
if focCostateOK, focVerdict = 'PASS'; else, focVerdict = 'FAIL'; end
if certCostateOK, certVerdict = 'PASS'; else, certVerdict = 'FAIL'; end
agree = (focCostateOK == certCostateOK);
if agree, agreeStr = 'AGREE'; else, agreeStr = 'DISAGREE'; end
rep.crossCheck = struct( ...
    'focPrimerAlignDeg', out.primerAlignDeg, ...      % raw-dual source (foc, solver's own opti.lam_g)
    'focSignPct',        rep.signPct, ...
    'certPrimerDirErr',  cert.primerDirErr, ...        % LS-reconstructed source (certify_minfuel_pmp)
    'certSignMatchFrac', cert.signMatchFrac, ...
    'focVerdict',        focVerdict, ...
    'certVerdict',       certVerdict, ...
    'agree',             agree);
info.crossCheck = rep.crossCheck;

fprintf(['\n---- Costate cross-check (two independent sources, same re-solved row) ----\n' ...
    ' raw-dual   (foc,     opti.lam_g)      : primerAlign %8.4f deg   signPct   %6.2f%%   %s\n' ...
    ' LS-recon.  (certify, adjoint+primer)  : primerDirErr %7.3e   signMatch %6.2f%%   %s\n' ...
    ' -> %s\n'], ...
    out.primerAlignDeg, rep.signPct, focVerdict, ...
    cert.primerDirErr, 100*cert.signMatchFrac, certVerdict, agreeStr);

% --- standard report (prints + saves the sidecar, WITH crossCheck attached) -
[~, tag] = fileparts(matPath);
resDir = fullfile(here, 'results');
foc_report(rep, tag, resDir);
end

% =============================================================================
function v = getfield_default(s, f, dflt)
% GETFIELD_DEFAULT  s.(f) if present and nonempty, else dflt.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
