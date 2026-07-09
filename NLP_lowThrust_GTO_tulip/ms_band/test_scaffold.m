% TEST_SCAFFOLD  Task-1 gate: paths resolve, problem factory sane.
setup_paths;
prob = ms_problem(1.05, 1.0);
ok = true;
ok = ok && abs(prob.c - 20.24)   < 0.05;          % ND exhaust velocity
ok = ok && abs(prob.Tmax - 0.627) < 0.005;        % ND max accel at m=1
ok = ok && numel(prob.rv0) == 6 && numel(prob.rvf) == 6;
ok = ok && abs(prob.tf - 1.05*6.290694) < 1e-12;
yDot = lt_pmp_eom_minfuel(0, [prob.rv0; 1; ones(7,1)], prob.Tmax, ...
                          prob.c, prob.muStar, prob.epsSmooth);
ok = ok && numel(yDot) == 14 && all(isfinite(yDot));
if ok, fprintf('PASS test_scaffold\n'); else, fprintf('FAIL test_scaffold\n'); end
