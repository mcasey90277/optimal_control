function [C, t] = recert_one(S, B, pool, rv0, j, m, extra)
% RECERT_ONE  Re-certify sheet candidate (j, m) from its stored seed.
% INPUTS: S sheet; B setup; pool; rv0; j; m; extra copts.  OUTPUTS: C; t sec.
c = S.cand{j}(m);  z = c.z;  Y = c.Y;  K = size(Y, 2);
rvf = B.stateA(S.sA(j));
seed = struct('tf', z(8), 'tGrid', linspace(0, z(8), K+1), 'Y', [Y, Y(:, end)]);
seed.Y(8:14, 1) = z(1:7);  seed.Y(1:7, 1) = [rv0(1:6); 1];
copts = S.policy;  copts.pool = pool;  copts.sA = S.sA(j);  copts.sD = 0;
for f = fieldnames(extra)', copts.(f{1}) = extra.(f{1}); end
t0 = tic;  C = certify_root(seed, rv0, rvf, B, copts);  t = toc(t0);
end
