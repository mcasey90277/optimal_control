% SMOKE_SEED_MINTIME  Quick sanity check of seed_from_mintime before the
% long eps-march: shape, finiteness, and the seed's residual norm at
% eps=1 (should be large but finite -- this is a wide-basin cold seed).
setup_paths;
factor = 1.01;  M = 24;
[Zseed, tJ] = seed_from_mintime(factor, M, 1.0);
fprintf('numel(Zseed) = %d (expect %d)\n', numel(Zseed), 14*M-7);
fprintf('numel(tJ)    = %d (expect %d)\n', numel(tJ), M+1);
fprintf('any(~isfinite(Zseed)) = %d\n', any(~isfinite(Zseed)));
prob = ms_problem(factor, 1.0);  prob.tJ = tJ;
R = ms_residual(Zseed, prob);
fprintf('||R(seed)|| at eps=1 = %.4e\n', norm(R));
fprintf('tJ = '); disp(tJ);
