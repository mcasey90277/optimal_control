% TEST_MS_JACOBIAN  Task-4 gate: complex-step Jacobian vs finite differences.
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
n = numel(Z);
colSel = unique([1 3 7, 8 15 22, n-13, n, randi(n, 1, 4)]);   % lam0 + joint cols
hFD = 1e-7;  errMax = 0;
for cIdx = colSel
    scale = max(1, abs(Z(cIdx)));
    Zp = Z;  Zp(cIdx) = Zp(cIdx) + hFD*scale;
    Rp = ms_residual(Zp, prob);
    Zm = Z;  Zm(cIdx) = Zm(cIdx) - hFD*scale;
    Rm = ms_residual(Zm, prob);
    colFD  = (Rp - Rm)/(2*hFD*scale);
    denom  = max(1, max(abs(colFD)));
    errMax = max(errMax, max(abs(J(:, cIdx) - colFD))/denom);
end
% structure: rows of continuity block k depend only on y_k and y_{k+1}
bw = 0;
for k = 1:M-1
    rowsK = 14*(k-1)+(1:14);
    colsAllowed = false(1, n);
    if k == 1, colsAllowed(1:7) = true; else, colsAllowed(7+14*(k-2)+(1:14)) = true; end
    colsAllowed(7+14*(k-1)+(1:14)) = true;
    bw = max(bw, nnz(any(J(rowsK, ~colsAllowed), 1)));
end
fprintf('max rel col err vs FD: %.2e   off-structure cols: %d\n', errMax, bw);
ok = errMax < 1e-6 && bw == 0;
if ok, fprintf('PASS test_ms_jacobian\n'); else, fprintf('FAIL test_ms_jacobian\n'); end
