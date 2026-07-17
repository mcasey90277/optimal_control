% TEST_MEE_SOLVER_SMOKE  Construction + short-iteration smoke for the L-domain
% MEE solver core (Task 3), mirroring test_solver_smoke.m's non-convergence-
% gated pattern: assert the full out struct returns without error at maxIter=5
% in both modes, sizes are correct, ΔL is a scalar variable, and the single
% ΔL column does not create a dense-KKT pathology (wall-time budget).
p = kepler_lt_params(10, 1500, 2000);
seedOpts = struct('thr', 0.5, 'betaMode', 'transverse', 'N', 75, 'nRev', 3);
[sg, X0, U0, dL0, sinfo] = mee_seed(p, seedOpts);
x0 = X0(:,1);

fields_needed = {'X','U','dL','success','ipoptStatus','maxDefect','maxUnit', ...
    'termErr','mf','m_f_kg','dV_kms','tf','switches','edge','lamDef', ...
    'LdotMin','incDeg'};

% (a) fixedtf mode, eps=1, tfTarget = 1.3x seed tEnd, maxIter=5.
tA = tic;
outA = casadi_lt_mee(sg, X0, U0, dL0, struct('par',p,'mode','fixedtf', ...
        'eps',1,'tfTarget',1.3*sinfo.tEnd,'x0',x0,'maxIter',5,'printLevel',0));
wallA = toc(tA);
for kf = 1:numel(fields_needed)
    assert(isfield(outA, fields_needed{kf}), 'fixedtf out missing field %s', fields_needed{kf});
end
assert(isequal(size(outA.X), [7 76]), 'fixedtf out.X size');
assert(isequal(size(outA.U), [4 76]), 'fixedtf out.U size');
assert(isscalar(outA.dL) && isfinite(outA.dL), 'fixedtf out.dL must be a finite scalar');
assert(islogical(outA.success) || isnumeric(outA.success), 'fixedtf out.success must exist');

% (b) mintime mode, maxIter=5: same completeness + thr pin held.
tB = tic;
outB = casadi_lt_mee(sg, X0, U0, dL0, struct('par',p,'mode','mintime', ...
        'x0',x0,'maxIter',5,'printLevel',0));
wallB = toc(tB);
for kf = 1:numel(fields_needed)
    assert(isfield(outB, fields_needed{kf}), 'mintime out missing field %s', fields_needed{kf});
end
assert(isequal(size(outB.X), [7 76]), 'mintime out.X size');
assert(isequal(size(outB.U), [4 76]), 'mintime out.U size');
assert(isscalar(outB.dL) && isfinite(outB.dL), 'mintime out.dL must be a finite scalar');
assert(all(abs(outB.U(4,:) - 1) < 1e-9), 'mintime thr pin must hold at every node');

% (c) sparsity/wall-time budget: a dense ΔL column would blow this at N=75.
budget_s = 120;
assert(wallA < budget_s, 'fixedtf 5-iter solve took %.1fs, budget %.0fs (dense-column suspect)', wallA, budget_s);
assert(wallB < budget_s, 'mintime 5-iter solve took %.1fs, budget %.0fs (dense-column suspect)', wallB, budget_s);

fprintf('test_mee_solver_smoke: ALL PASS (wallA=%.2fs wallB=%.2fs, dL0=%.4f outA.dL=%.4f outB.dL=%.4f)\n', ...
    wallA, wallB, dL0, outA.dL, outB.dL);
