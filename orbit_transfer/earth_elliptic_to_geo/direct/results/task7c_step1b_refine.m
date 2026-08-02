% TASK7C_STEP1B_REFINE  Mesh-refine the certified 1 N npr=15 min-time anchor
% (MEE_mintime_T10_npr15.mat, defect=3.05e-14) UP to production density
% (25 nodes/rev), single warmTight=true eps=0-style solve -- same pattern as
% nodestudy_mee.m>solve_warm_node's fixedtf refinement, adapted to 'mintime'
% mode. No homotopy re-sweep: the npr=15 point is already essentially exact,
% so the refined mesh only needs to re-satisfy the (denser) defect/terminal
% conditions from an excellent starting shape.
here = '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo';
cd(here); addpath(here);
resDir = fullfile(here, 'results');

S15 = load(fullfile(resDir, 'MEE_mintime_T10_npr15.mat'));
base = S15.out;
thrustN = 1;
par = kepler_lt_params(thrustN, 1500, 2000);

nprTarget = 25;
revs = base.dL_mt / (2*pi);
N = round(nprTarget * revs);
fprintf('Refining npr=15 (N=%d) -> npr=25 (N=%d), revs=%.4f\n', base.N, N, revs);

sigmaSrc = linspace(0, 1, base.N + 1).';
sigmaDst = linspace(0, 1, N + 1).';
W = interp_warmstart(base.solverOut.X, base.solverOut.U, base.solverOut.dL, sigmaSrc, sigmaDst);
x0 = W.X(:,1);

t0 = tic;
o = casadi_lt_mee(sigmaDst, W.X, W.U, W.dL, struct('par', par, 'mode', 'mintime', ...
    'x0', x0, 'maxIter', 1500, 'warmTight', true, 'printLevel', 3));
fprintf('REFINE solve done in %.1f s: status=%s defect=%.4e termErr=%.4e tf=%.6f revs=%.4f\n', ...
        toc(t0), o.ipoptStatus, o.maxDefect, o.termErr, o.tf, o.dL/(2*pi));

certified = o.success && o.maxDefect < 1e-8 && o.termErr < 1e-8;
fprintf('certified=%d\n', certified);
save(fullfile(resDir, 'task7c_step1b_out.mat'), 'o', 'certified', 'N');
