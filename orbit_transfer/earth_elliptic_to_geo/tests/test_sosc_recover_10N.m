% TEST_SOSC_RECOVER_10N  Task-4 gate: sosc_recover_kkt rebuilds+warm-resolves
% the 10 N NLP at the saved primal and assembles the KKT objects (x, lam_g,
% gval, grad_f, sparse H, sparse A_all) in Opti's native symbols.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));
saved = sosc_load_row(fullfile(module_root(),'results','MEE_M2_10N.mat'));
tol = sosc_defaults();
R = sosc_recover_kkt(saved, tol);
assert(R.recoverOK, 'recovery failed: %s', R.ipoptStatus);
assert(numel(R.x)==R.n && numel(R.lam_g)==R.m, 'x/lam_g dims');
assert(isequal(size(R.H),[R.n R.n]) && issparse(R.H), 'H shape/sparse');
assert(size(R.A_all,1)==R.m && size(R.A_all,2)==R.n && issparse(R.A_all), 'A_all shape/sparse');
assert(R.drift < tol.drift, sprintf('drift %.3e >= %.1e', R.drift, tol.drift));
fprintf('test_sosc_recover_10N PASSED (n=%d, m=%d, drift=%.2e)\n', R.n, R.m, R.drift);
