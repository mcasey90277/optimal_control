% TASK5_DIAG_COL49  h-sweep of the FD error on Jacobian-test column 49.
%
% Reproduces test_ms_jacobian's setup (RelTol 1e-13 now in ms_problem) and
% sweeps the central-FD step for column 49 (arc-4 lambda_m), comparing each
% FD column against the CS Jacobian column (certified h-independent to
% ~2e-13 by Gate 2). Shows whether the FD error is noise-dominated (~1/h),
% i.e. whether the h-vs-h/2 self-consistency tolerance underestimates it.
% Diagnostic only; not part of the committed test suite.
setup_paths;
rng(7);
epsS = 0.3;
prob = ms_problem(1.03, epsS);
ref  = run_gto_tulip_indirect(false);
lam0 = ref.zSol(1:7);
y0   = [prob.rv0; prob.m0; lam0];
sol  = ode113(@(t,y) lt_pmp_eom_minfuel(t, y, prob.Tmax, prob.c, ...
              prob.muStar, epsS), [0 prob.tf], y0, prob.odeOpts);
M    = 4;
prob.tJ = arc_boundaries_tau(sol.x, sol.y(1:3,:), M, prob.muStar);
yJ   = deval(sol, prob.tJ);
Z    = ms_pack(lam0, yJ(:, 2:M));

[~, J] = ms_residual(Z, prob);
cIdx  = 49;
scale = max(1, abs(Z(cIdx)));
jCol  = full(J(:, cIdx));
fprintf('col %d: |Z| = %.3e, scale = %.3g, max|Jcol| = %.3e\n', ...
        cIdx, abs(Z(cIdx)), scale, max(abs(jCol)));
for hFD = [1e-8 1e-7 1e-6 1e-5 1e-4]
    h  = hFD*scale;
    Zp = Z;  Zp(cIdx) = Zp(cIdx) + h;
    Zm = Z;  Zm(cIdx) = Zm(cIdx) - h;
    colFD   = (ms_residual(Zp, prob) - ms_residual(Zm, prob))/(2*h);
    Zp2 = Z; Zp2(cIdx) = Zp2(cIdx) + h/2;
    Zm2 = Z; Zm2(cIdx) = Zm2(cIdx) - h/2;
    colFD2  = (ms_residual(Zp2, prob) - ms_residual(Zm2, prob))/h;
    fprintf('hFD %.0e: |J-FD| %.3e   selfErr(h vs h/2) %.3e\n', ...
            hFD, max(abs(jCol - colFD)), max(abs(colFD - colFD2)));
end
