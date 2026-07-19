function R = sosc_recover_kkt(saved, tol)
% SOSC_RECOVER_KKT  Rebuild the NLP at a saved primal, warm re-solve to recover
% the full multiplier set, and assemble the KKT objects in Opti's native
% (unscaled) symbols.
%
% INPUTS:
%   saved - struct from sosc_load_row (primal + config)
%   tol   - struct from sosc_defaults
% OUTPUTS:
%   R - struct: .recoverOK .x[nx1] .lam_g[mx1] .gval[mx1] .grad_f[nx1]
%       .H[nxn sparse] .A_all[mxn sparse] .creg .vreg .drift .n .m .ipoptStatus
% REFERENCES:
%   [1] process/DESIGN_sosc.md sec 4.2. [2] casadi_lt_mee.m (returnModel hook).
import casadi.*
par  = kepler_lt_params(saved.thrustN, saved.m0kg, saved.ispS);
opts = struct('par',par,'mode','fixedtf','eps',0,'tfTarget',saved.tfTarget, ...
    'x0',saved.X(:,1),'xf',saved.xf,'maxIter',saved.maxIter, ...
    'warmTight',true,'printLevel',0,'returnModel',true);
o = casadi_lt_mee(saved.sigma, saved.X, saved.U, saved.dL, opts);
R.ipoptStatus = o.ipoptStatus;
R.recoverOK   = o.success && o.maxDefect < tol.feas;
R.drift = max(abs(o.X(:) - saved.X(:)));
if ~R.recoverOK, R.x=[]; R.lam_g=[]; R.gval=[]; R.grad_f=[]; R.H=[]; R.A_all=[];
    R.creg=[]; R.vreg=[]; R.n=0; R.m=0; return; end

opti = o.model.opti;  sol = opti.debug;   % the solved Opti (sol from last solve)
% Native symbols:
x   = opti.x;   g = opti.g;   f = opti.f;   lam = opti.lam_g;
R.x     = full(sol.value(x));
R.lam_g = full(sol.value(lam));
R.gval  = full(sol.value(g));
% Gradient, Jacobian, Lagrangian Hessian as CasADi Functions, evaluated at soln:
gradF = gradient(f, x);
Jg    = jacobian(g, x);
Hlag  = hessian(f + lam.'*g, x);           % returns [H, grad] -> take H
Fkkt  = Function('Fkkt', {x, lam}, {gradF, Jg, Hlag});
[gf, A, H] = Fkkt(R.x, R.lam_g);
R.grad_f = full(gf);
R.A_all  = sparse(A);
R.H      = sparse(H);
R.creg = o.model.creg;  R.vreg = o.model.vreg;
R.n = numel(R.x);  R.m = numel(R.gval);
end
