function K = sosc_kkt_residual(R, tol)
% SOSC_KKT_RESIDUAL  Resolve the one global Lagrangian-sign ambiguity, then
% re-check first-order KKT residuals at the recovered point.
%
% INPUTS: R - sosc_recover_kkt struct; tol - sosc_defaults struct
% OUTPUTS: K - struct .sign .signOK .stat .primalEq .primalIneq .dualFeas
%              .comp .pass
% REFERENCES: process/DESIGN_sosc.md sec 4.3; verify_pmp_mee.m:112-121 (sign trick).
isEq   = strcmp({R.creg.kind},'eq'); %#ok<NASGU>
% expand per-group kind to per-row masks
kindRow = strings(R.m,1);
for i=1:numel(R.creg), kindRow(R.creg(i).rows) = R.creg(i).kind; end
eqRow   = kindRow=="eq";
ineqRow = ~eqRow;

% (1) global sign: choose s minimizing ||grad_f + s*A'*lam||_inf
rP = R.grad_f + (R.A_all.' * R.lam_g);
rM = R.grad_f - (R.A_all.' * R.lam_g);
if norm(rP,inf) <= norm(rM,inf), K.sign=+1; stat=norm(rP,inf);
else,                            K.sign=-1; stat=norm(rM,inf); end
K.stat   = stat;
K.signOK = stat < tol.stat;

% (2) primal feasibility (eq residual; ineq violation for g<=0 canonical form)
K.primalEq   = max(abs(R.gval(eqRow)), [], 'omitnan');
if isempty(K.primalEq), K.primalEq = 0; end
% ineqHi: g<=bound -> viol = max(0, g-bound); ineqLo: g>=bound -> viol=max(0,bound-g)
viol = zeros(R.m,1);
for i=1:numel(R.creg)
    c=R.creg(i); if strcmp(c.kind,'eq'), continue; end
    gv=R.gval(c.rows);
    if strcmp(c.kind,'ineqHi'), viol(c.rows)=max(0, gv - c.bound);
    else,                        viol(c.rows)=max(0, c.bound - gv); end
end
K.primalIneq = max(viol);

% (3) dual feasibility: s*lam_g >= 0 for all inequality rows (g<=0 convention)
lamSigned = K.sign * R.lam_g;
K.dualFeas = max(-lamSigned(ineqRow), [], 'omitnan');   % worst negative
if isempty(K.dualFeas), K.dualFeas = 0; end

% (4) complementarity: |lam * slack| over inequalities
slack = zeros(R.m,1);
for i=1:numel(R.creg)
    c=R.creg(i); if strcmp(c.kind,'eq'), continue; end
    gv=R.gval(c.rows);
    if strcmp(c.kind,'ineqHi'), slack(c.rows)=c.bound-gv; else, slack(c.rows)=gv-c.bound; end
end
K.comp = max(abs(R.lam_g(ineqRow).*slack(ineqRow)), [], 'omitnan');
if isempty(K.comp), K.comp = 0; end

K.pass = K.signOK && K.primalEq<tol.feas && K.primalIneq<tol.feas && ...
         K.dualFeas<tol.dual && K.comp<tol.comp;
end
