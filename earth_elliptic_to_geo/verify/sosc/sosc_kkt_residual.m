function K = sosc_kkt_residual(R, tol)
% SOSC_KKT_RESIDUAL  Resolve the one global Lagrangian-sign ambiguity, then
% re-check first-order KKT residuals at the recovered point. Bounds and per-row
% constraint KIND are sourced from the canonical NLP bounds R.lbg/R.ubg (NOT
% from creg.bound); dual feasibility is checked PER KIND.
%
% INPUTS: R - sosc_recover_kkt struct (needs .lbg .ubg); tol - sosc_defaults
% OUTPUTS: K - struct .sign .signOK .stat .primalEq .primalIneq .dualFeas
%              .comp .pass
% REFERENCES: process/DESIGN_sosc.md sec 11.2 (bound sourcing), 11.3 (per-kind
%   dual feasibility); verify_pmp_mee.m:112-121 (sign trick).
lbg = R.lbg(:);  ubg = R.ubg(:);

% Per-row KIND from bounds (§11.2)
isEq     = (lbg == ubg);
hasUpper = isfinite(ubg) & ~isEq;
hasLower = isfinite(lbg) & ~isEq;

% (1) global sign: choose s minimizing ||grad_f + s*A'*lam||_inf
rP = R.grad_f + (R.A_all.' * R.lam_g);
rM = R.grad_f - (R.A_all.' * R.lam_g);
if norm(rP,inf) <= norm(rM,inf), K.sign=+1; stat=norm(rP,inf);
else,                            K.sign=-1; stat=norm(rM,inf); end
K.stat   = stat;
K.signOK = stat < tol.stat;

% (2) primal feasibility -- uniform bound violation for ALL rows (§11.2)
viol = max(lbg - R.gval, 0) + max(R.gval - ubg, 0);
K.primalEq   = safemax(viol(isEq));
K.primalIneq = safemax(viol(~isEq));

% (3) dual feasibility -- per kind (§11.3): with ls = s*lam_g,
%     upper-active rows require ls >= -tol.dual; lower-active rows require
%     ls <= +tol.dual; equality rows unconstrained in sign.
ls = K.sign * R.lam_g;
dvU = safemax(max(0, -ls(hasUpper)));   % upper: violation when ls < 0
dvL = safemax(max(0,  ls(hasLower)));   % lower: violation when ls > 0
K.dualFeas = max(dvU, dvL);

% (4) complementarity: slack to the NEAREST active bound (finite bounds only),
%     over inequality rows.
distLo = R.gval - lbg;  distLo(~isfinite(lbg)) = inf;
distHi = ubg - R.gval;  distHi(~isfinite(ubg)) = inf;
slack  = min(distLo, distHi);
ineqRow = ~isEq;
K.comp  = safemax(abs(R.lam_g(ineqRow) .* slack(ineqRow)));

K.pass = K.signOK && K.primalEq<tol.feas && K.primalIneq<tol.feas && ...
         K.dualFeas<tol.dual && K.comp<tol.comp;
end

function v = safemax(x)
% max over possibly-empty vector, ignoring NaN, defaulting to 0.
x = x(~isnan(x));
if isempty(x), v = 0; else, v = max(x); end
end
