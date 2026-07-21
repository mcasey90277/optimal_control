% TEST_MS_REPRODUCE_MINTIME  M1(a): MS reproduces the min-time anchor at 1.00x.
setup_paths;
epsS = 1e-3;                     % S << 0 everywhere -> throttle saturates at 1
prob = ms_problem(1.00, epsS);
ref  = run_gto_tulip_indirect(false);
prob.tf = ref.zSol(8);           % use the CONVERGED min-time tf exactly
lam0 = ref.zSol(1:7);
y0   = [prob.rv0; prob.m0; lam0];
sol  = ode113(@(t,y) lt_pmp_eom_minfuel(t, y, prob.Tmax, prob.c, ...
              prob.muStar, epsS), [0 prob.tf], y0, prob.odeOpts);
M    = 24;
prob.tJ = arc_boundaries_tau(sol.x, sol.y(1:3,:), M, prob.muStar);
yJ   = deval(sol, prob.tJ);
Zseed = ms_pack(lam0, yJ(:, 2:M));

% At RelTol 1e-13 the seed residual rises to ~2e-5 (tighter integrator
% exposes the min-time reference's own ~1e-6 gate); LM absorbs it in one step.
out  = ms_solve(Zseed, prob, 1e-9, 100);
traj = ms_traj(out.Z, prob);
fprintf('||R|| = %.3e   dV = %.4f km/s   prop = %.4f kg   bang %.1f%%\n', ...
        out.resNorm, traj.dV_kms, traj.prop_kg, 100*traj.bangFrac);
ok = out.success && abs(traj.dV_kms - 4.4665) < 0.002 ...
     && abs(traj.prop_kg - 2.9247) < 0.002;
if ok
    fprintf('PASS test_ms_reproduce_mintime\n');
else
    error('FAIL test_ms_reproduce_mintime: ||R||=%.3e dV=%.4f prop=%.4f', ...
          out.resNorm, traj.dV_kms, traj.prop_kg);
end
