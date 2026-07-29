% TASK5_DIAG_ROWCLASS  Row-class breakdown of the Task-5 stalled LM solution.
%
% Re-runs the deterministic M1(a) solve (same seed, same options -> same
% converged Z), then evaluates ms_residual at that Z and reports max|R| and
% block 2-norms per row class: joint state rows (1:7 per joint), joint
% costate rows (8:14 per joint), terminal rows. Diagnostic only; not part
% of the committed test suite.
setup_paths;
epsS = 1e-3;
prob = ms_problem(1.00, epsS);
ref  = run_gto_tulip_indirect(false);
prob.tf = ref.zSol(8);
lam0 = ref.zSol(1:7);
y0   = [prob.rv0; prob.m0; lam0];
sol  = ode113(@(t,y) lt_pmp_eom_minfuel(t, y, prob.Tmax, prob.c, ...
              prob.muStar, epsS), [0 prob.tf], y0, prob.odeOpts);
M    = 24;
prob.tJ = arc_boundaries_tau(sol.x, sol.y(1:3,:), M, prob.muStar);
yJ   = deval(sol, prob.tJ);
Zseed = ms_pack(lam0, yJ(:, 2:M));

out = ms_solve(Zseed, prob, 1e-9, 100);
save('task5_converged.mat', 'out', 'prob');

R   = ms_residual(out.Z, prob);
Mm1 = M - 1;
Rj  = reshape(R(1:14*Mm1), 14, Mm1);
Rt  = R(14*Mm1+1:end);
fprintf('re-eval ||R|| = %.6e (solver reported %.6e)\n', norm(R), out.resNorm);
fprintf('max|R| joint STATE rows (1:7):    %.3e\n', max(abs(Rj(1:7,:)), [], 'all'));
fprintf('max|R| joint COSTATE rows (8:14): %.3e\n', max(abs(Rj(8:14,:)), [], 'all'));
fprintf('max|R| TERMINAL rows:             %.3e\n', max(abs(Rt)));
fprintf('block 2-norms: state %.3e  costate %.3e  terminal %.3e  total %.3e\n', ...
        norm(Rj(1:7,:), 'fro'), norm(Rj(8:14,:), 'fro'), norm(Rt), norm(R));
% top 10 residual rows, mapped to joint + component
[~, ord] = sort(abs(R), 'descend');
comp = {'rx','ry','rz','vx','vy','vz','m', ...
        'lrx','lry','lrz','lvx','lvy','lvz','lm'};
for d = 1:10
    rr = ord(d);
    if rr <= 14*Mm1
        kJ = ceil(rr/14); cc = rr - 14*(kJ-1);
        fprintf('row %3d |R|=%.3e  joint %2d  %s  (tJ=%.4f)\n', ...
                rr, abs(R(rr)), kJ, comp{cc}, prob.tJ(kJ+1));
    else
        fprintf('row %3d |R|=%.3e  TERMINAL comp %d\n', rr, abs(R(rr)), rr-14*Mm1);
    end
end
