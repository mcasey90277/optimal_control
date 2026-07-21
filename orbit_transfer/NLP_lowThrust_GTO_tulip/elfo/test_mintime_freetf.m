% TEST_MINTIME_FREETF  Smoke test: casadi_mintime_freetf constructs, runs, and
% returns the all-burn (s==1) min-time LAYOUT on the ELFO energy seed. These
% assertions are convergence-INDEPENDENT (they hold after a short 60-iter solve).
% The numeric all-burn mass identity mf=1-(Tmax/c)*tf holds only once the defects
% are converged, so it is checked at the converged anchor in Task 2, not here.
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
cfg = minfuel_config();  p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
S = load(fullfile(here,'results','energy_elfo_f0990.mat'));
o = struct('pSund',S.pSund,'qSund',S.qSund,'moonZone',S.moonZone, ...
           'maxIter',60,'warmTight',false);
out = casadi_mintime_freetf(S.sigma,S.rv0,S.rvf,p.Tmax,p.c,p.muStar, ...
                            S.X,S.U(1:3,:),S.tauf0,o);
nN = size(S.X,2);
assert(isequal(size(out.X),[9 nN]), 'X layout wrong: %dx%d', size(out.X,1),size(out.X,2));
% THE wiring proof: control is 3-row (throttle is NOT a decision variable => s==1
% is structural, thrust/mdot cannot carry a hidden s factor).
assert(isequal(size(out.U),[3 nN]), 'U must be 3-row steering, got %dx%d', size(out.U,1),size(out.U,2));
assert(all(isfield(out,{'tf','minR1','tMonotone','cScale','mf'})), 'missing min-time fields');
assert(isfinite(out.tf) && out.tf > 0, 'tf not a positive finite number: %g', out.tf);
assert(out.maxUnit < 1e-1, 'unit-steering not tracking (loose 60-iter bound): %.2e', out.maxUnit);
% informational (NOT asserted -- exact only at convergence, Task 2):
mf_pred = 1 - (p.Tmax/p.c)*out.tf;
fprintf('TEST_MINTIME_FREETF: PASS (tf=%.4f ND, maxUnit=%.2e; mf=%.4f vs all-burn pred %.4f)\n', ...
        out.tf, out.maxUnit, out.mf, mf_pred);
