% TEST_MS_RESIDUAL  Task-3 gate: MS residual consistency vs single shooting.
setup_paths;
epsS = 1.0;                       % smooth throttle: benign dynamics
prob = ms_problem(1.05, epsS);
ref  = run_gto_tulip_indirect(false);          % min-time costates as lam0
lam0 = ref.zSol(1:7);
y0   = [prob.rv0; prob.m0; lam0];
sol  = ode113(@(t,y) lt_pmp_eom_minfuel(t, y, prob.Tmax, prob.c, ...
              prob.muStar, epsS), [0 prob.tf], y0, prob.odeOpts);
M    = 12;
prob.tJ = arc_boundaries_tau(sol.x, sol.y(1:3,:), M, prob.muStar);
% Snap interior joints onto sol's own step mesh so deval returns EXACT
% states of the one continuous integration. Raw deval interpolation error
% (~1e-10 abs on lamR ~ 2e2 at RelTol 1e-12) is amplified ~1e2 by
% perigee-arc STMs to ~3e-8, breaking the 1e-8 gate for correct code;
% mesh-snapped M=12 measures cont 3e-9 / term 6e-10 (see task-3 report).
for k = 2:M
    [~, idx] = min(abs(sol.x - prob.tJ(k)));
    prob.tJ(k) = sol.x(idx);
end
yJ   = deval(sol, prob.tJ);                    % 14 x (M+1)
Z    = ms_pack(lam0, yJ(:, 2:M));
R    = ms_residual(Z, prob);
contDef = max(abs(R(1:14*(M-1))));
Rss  = shoot_residual_minfuel(lam0, prob.tf, prob.rv0, prob.m0, prob.rvf, ...
                              prob.Tmax, prob.c, prob.muStar, epsS);
termErr = max(abs(R(14*(M-1)+1:end) - Rss));
traj = ms_traj(Z, prob);
fprintf('cont defect %.2e   term-vs-SS %.2e   traj joint defect %.2e\n', ...
        contDef, termErr, traj.maxJointDefect);
ok = contDef < 1e-8 && termErr < 1e-8 && traj.maxJointDefect < 1e-8 ...
     && isfinite(traj.dV_kms) && numel(traj.S) == numel(traj.t);
if ok, fprintf('PASS test_ms_residual\n'); else, fprintf('FAIL test_ms_residual\n'); end
