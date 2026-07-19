function v = sosc_decide(K, AS, IN)
% SOSC_DECIDE  Map KKT/active-set/reduced-inertia results to a verdict
% (process/DESIGN_sosc.md sec 11.5, ordered logic; supersedes sec 5).
%
% INPUTS:
%   K  - kkt_residual struct; needs .pass .signOK .sign
%   AS - active_set struct;   needs .nWeak .licq (licq REPORTED, never gates)
%   IN - inertia struct;      needs .red (.npos .nneg .nzero) .redConsistent
%
% OUTPUTS:
%   v - struct .verdict .reason .status
%       .verdict in {PASS, WEAK_MIN, FAIL, INCONCLUSIVE, ERROR}
%       .status  in {certified-sosc, certified-weak-min, feasible-only,
%                    certified-feasibility+sosc-inconclusive}
%
% REFERENCES: process/DESIGN_sosc.md sec 11.5 (verdict taxonomy),
%             sec 11.6 (tiered gate). Only FAIL demotes; WEAK_MIN is a
%             POSITIVE certificate. LICQ is reported but does NOT gate.
red = IN.red;
if ~K.pass || ~K.signOK
    v.verdict = 'ERROR';
    v.reason  = 'KKT residual/sign check failed (no trustworthy KKT point)';
elseif ~IN.redConsistent
    v.verdict = 'INCONCLUSIVE';
    v.reason  = sprintf(['reduced-inertia consistency failed (sprank untrustworthy): ' ...
        'red=[%d %d %d], rankA=%d'], red.npos, red.nneg, red.nzero, IN.rankA);
elseif red.nneg == 0 && red.nzero == 0
    v.verdict = 'PASS';
    v.reason  = 'reduced Hessian PD on the critical cone (strict local min)';
elseif red.nneg == 0 && red.nzero > 0
    v.verdict = 'WEAK_MIN';
    v.reason  = sprintf(['reduced Hessian PSD (no negative curvature), %d flat ' ...
        'direction(s) -> weak local minimum (bang-bang signature)'], red.nzero);
elseif red.nneg > 0 && AS.nWeak == 0
    v.verdict = 'FAIL';
    v.reason  = sprintf(['reduced Hessian has %d negative curvature direction(s) ' ...
        'with no weakly-active junctions -> genuine descent direction (not a local min)'], ...
        red.nneg);
else   % red.nneg > 0 && AS.nWeak > 0
    v.verdict = 'INCONCLUSIVE';
    v.reason  = sprintf(['reduced Hessian has %d negative curvature direction(s) but ' ...
        '%d weakly-active junction(s) make the critical cone a strict subset -> ' ...
        'cannot conclude a saddle'], red.nneg, AS.nWeak);
end

switch v.verdict
    case 'PASS',     v.status = 'certified-sosc';
    case 'WEAK_MIN', v.status = 'certified-weak-min';
    case 'FAIL',     v.status = 'feasible-only';
    otherwise,       v.status = 'certified-feasibility+sosc-inconclusive';
end
end
