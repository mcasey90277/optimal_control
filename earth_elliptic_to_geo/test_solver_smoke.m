% TEST_SOLVER_SMOKE  Construction + short-iteration smoke (no convergence gate).
p  = kepler_lt_params(10, 1500, 2000);
P0 = 11625/p.LU_km;
[r0, v0] = elements_to_cart(P0, 0.75, 0, 0, 0, pi, p.mu);
rv0 = [r0; v0];
[sg, X0, U0, tauf0] = seed_2body(p, rv0, struct('sbar',1,'tDur',5,'N',80));
term = geo_terminal('manifold', p, []);
out = casadi_lt_2body(sg, X0, U0, tauf0, term, struct('par',p,'mode','mintime', ...
        'rv0',rv0,'maxIter',5,'printLevel',0));
assert(isstruct(out) && isfield(out,'maxDefect') && isfield(out,'lamDef'));
assert(size(out.X,1) == 9 && size(out.U,1) == 4);
% fixedtf construction path too
out2 = casadi_lt_2body(sg, X0, U0, tauf0, geo_terminal('fixed',p,pi+8), ...
        struct('par',p,'mode','fixedtf','eps',1,'tfTarget',6,'rv0',rv0, ...
               'maxIter',5,'printLevel',0));
assert(isfield(out2,'tf') && ~out2.success);   % 5 iters won't converge; must not error
fprintf('test_solver_smoke: ALL PASS\n');
