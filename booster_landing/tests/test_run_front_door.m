% TEST_RUN_FRONT_DOOR  Front-door contract on a FAST config: runs end to
% end from only setup_paths, honors cfg overrides (the advertised-but-
% ignored-options bug is the canonical failure -- verify N actually
% reached the solver via solution size), writes its products.
%
% GRID ADAPTATION from the brief's literal N=20/Nconv=60 (documented per
% this task's own supersession note: "the brief predates tasks 5-9"):
% task-5's report established that certify_pdg's G5 primer-alignment gate
% needed N>=40 to clear the then-<1 deg threshold (it measured ~1.7 deg at
% N=20). That grid dependence turned out to be an O(h) time-base bug in
% the gate itself, fixed 2026-08-09 -- the primer angle is now at the acos
% machine-precision floor on every grid (see certify_pdg.m's "G5 primer
% TIME BASE" note), so this is no longer WHY N=40. It stays because it is
% the coarse pair test_certify_nominal.m characterizes and asserts
% all_pass on. N=40/Nconv=90 is the exact coarse pair
% test_certify_nominal.m already uses and already knows clears ALL FIVE
% gates outright, so this test reuses it rather than picking a new grid
% to characterize: fast (measured ~23 s -- dominated by the TVLQR Riccati
% integration + 6-run MC, not the coarse NLP solves themselves, which
% take well under 1 s each), and both a real all_pass AND a genuine
% cfg.P override check (N=40 is not any function's internal default).
%
% FIX (2026-08-09 review): dropped the `|| ischar(R.rep.G3_pass)` OR
% clause from the gate assertion below -- the front door ALWAYS supplies
% solV (see run_booster_landing.m step [2/5]), so G3 is never 'skipped'
% here and that clause was dead code that could only ever mask a real
% all_pass=false, never legitimately fire.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

od  = fullfile(tempdir, 'bl_front_smoke');
cfg = struct('doMovie', false, 'Nrun', 6, 'outdir', od, ...
             'P', struct('N', 40, 'Nconv', 90));
R   = run_booster_landing(cfg);
assert(size(R.solC.X, 2) == 41, 'cfg.P.N did not reach the solver');
assert(R.rep.all_pass, 'gates failed on fast config');
assert(isfile(fullfile(od, 'booster_run.mat')), 'products not written');
assert(isfile(fullfile(od, 'pdg_solution.png')), 'solution plot missing');
assert(isfile(fullfile(od, 'footprint.png')), 'footprint plot missing');
assert(numel(R.mc.ok) == cfg.Nrun, 'cfg.Nrun did not reach run_monte_carlo');
assert(~isfile(fullfile(od, 'landing.mp4')), ...
       'cfg.doMovie=false was ignored -- landing.mp4 should not exist');
fprintf('test_run_front_door PASS\n');
