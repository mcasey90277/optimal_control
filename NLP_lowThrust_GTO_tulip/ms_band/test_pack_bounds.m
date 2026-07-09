% TEST_PACK_BOUNDS  Task-2 gate: pack/unpack round-trip + tau joint placement.
setup_paths;
prob = ms_problem(1.05, 1.0);
M = 6;
lam0 = (1:7).';  yInt = reshape(1:14*(M-1), 14, M-1) + 100;
Z = ms_pack(lam0, yInt);
prob.tJ = linspace(0, prob.tf, M+1);
[lam0b, yJ] = ms_unpack(Z, prob);
ok = isequal(lam0b, lam0) && isequal(yJ(:,2:end), yInt) ...
     && isequal(yJ(1:6,1), prob.rv0) && yJ(7,1) == 1 && isequal(yJ(8:14,1), lam0);

% tau joints: constant kappa (circular r1) -> uniform in time
K = 1001; t = linspace(0, 10, K);
rCirc = [cos(t) - prob.muStar; sin(t); zeros(1,K)];       % r1 = 1 everywhere
tJc = arc_boundaries_tau(t, rCirc, 5, prob.muStar);
ok = ok && max(abs(tJc - linspace(0, 10, 6))) < 1e-6;

% varying kappa: small r1 in first half -> more joints early
rHalf = rCirc;  rHalf(1, t < 5) = 0.1 - prob.muStar;  rHalf(2, t < 5) = 0;
tJv = arc_boundaries_tau(t, rHalf, 5, prob.muStar);
ok = ok && sum(tJv(2:end-1) < 5) > sum(tJv(2:end-1) >= 5) ...
     && all(diff(tJv) > 0) && tJv(1) == 0 && abs(tJv(end) - 10) < 1e-12;
if ok, fprintf('PASS test_pack_bounds\n'); else, fprintf('FAIL test_pack_bounds\n'); end
