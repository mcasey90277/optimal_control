function v = sosc_decide(K, AS, IN)
% SOSC_DECIDE  Map KKT/active-set/reduced-inertia results to a verdict
% (process/DESIGN_sosc.md sec 12.2, ordered logic; supersedes sec 11.5 and 5).
% Driven by the DIRECT reduced-Hessian inertia IN.red plus the zt-sensitivity
% flag IN.sensStable: a negative reduced eigenvalue is trusted only when its
% count is threshold-stable across the zt-band (else INCONCLUSIVE). LICQ and
% redConsistent are retired (Z gives the exact rank; no sprank).
%
% INPUTS:
%   K  - kkt_residual struct; needs .pass .signOK .sign
%   AS - active_set struct;   needs .nWeak (.licq REPORTED, never gates)
%   IN - inertia struct;      needs .robust .sensStable .red (.nneg .nzero)
%
% OUTPUTS:
%   v - struct .verdict .reason .status
%       .verdict in {PASS, WEAK_MIN, FAIL, INCONCLUSIVE, ERROR}
%       .status  in {certified-sosc, certified-weak-min, feasible-only,
%                    certified-feasibility+sosc-inconclusive}
%
% REFERENCES: process/DESIGN_sosc.md sec 12.2 (verdict), sec 11.6 (tiered
%             gate). Only FAIL demotes; WEAK_MIN is a POSITIVE certificate.
red = IN.red;
if ~K.pass || ~K.signOK
    v.verdict = 'ERROR';
    v.reason  = 'KKT residual/sign check failed (no trustworthy KKT point)';
elseif ~IN.robust
    v.verdict = 'INCONCLUSIVE';
    v.reason  = sprintf(['reduced inertia not computable at this scale ' ...
        '(dense null-space intractable, method=%s)'], IN.method);
elseif ~IN.sensStable
    v.verdict = 'INCONCLUSIVE';
    v.reason  = sprintf(['reduced-Hessian negative count is zt-sensitive across ' ...
        'the band [%s] -> near-flat directions of unresolvable sign'], ...
        num2str(IN.nnegBand));
elseif red.nneg > 0 && AS.nWeak == 0
    v.verdict = 'FAIL';
    v.reason  = sprintf(['reduced Hessian has %d stably-negative curvature ' ...
        'direction(s) with no weakly-active junctions -> genuine descent ' ...
        'direction on the critical cone (not a local min)'], red.nneg);
elseif red.nneg > 0 && AS.nWeak > 0
    v.verdict = 'INCONCLUSIVE';
    v.reason  = sprintf(['reduced Hessian has %d negative curvature direction(s) but ' ...
        '%d weakly-active junction(s) make the critical cone a strict subset -> ' ...
        'cannot conclude a saddle'], red.nneg, AS.nWeak);
elseif red.nneg == 0 && red.nzero == 0
    v.verdict = 'PASS';
    v.reason  = 'reduced Hessian PD on the critical cone (strict local min)';
else   % red.nneg == 0 && red.nzero > 0
    v.verdict = 'WEAK_MIN';
    v.reason  = sprintf(['reduced Hessian PSD (no negative curvature), %d flat ' ...
        'direction(s) -> weak local minimum (bang-bang signature)'], red.nzero);
end

switch v.verdict
    case 'PASS',     v.status = 'certified-sosc';
    case 'WEAK_MIN', v.status = 'certified-weak-min';
    case 'FAIL',     v.status = 'feasible-only';
    otherwise,       v.status = 'certified-feasibility+sosc-inconclusive';
end
end
