% TEST_FOC_IPOPT_INERTIA  Unit test for the delta_w tail-verdict logic.
% Pure logic, no solve: exercises the four regHistory shapes from the brief
% (all-zero tail, nonzero tail, early-only regularization, empty history).
root = fileparts(fileparts(mfilename('fullpath'))); cd(root);
setup_verify_common;

ic = foc_ipopt_inertia(zeros(1,50));
assert(ic.certLocalMin, 'all-zero delta_w history must certify');
assert(contains(ic.verdict, 'LOCAL MIN'), 'certified verdict must read LOCAL MIN');

ic = foc_ipopt_inertia([zeros(1,45), 1e-3*ones(1,5)]);
assert(~ic.certLocalMin, 'nonzero tail delta_w must NOT certify');
assert(contains(ic.verdict, 'NOT CERTIFIED'), 'uncertified verdict must read NOT CERTIFIED');

ic = foc_ipopt_inertia([1e-2*ones(1,10), zeros(1,40)]);
assert(ic.certLocalMin, 'early regularization with a clean tail must certify');

ic = foc_ipopt_inertia([]);
assert(~ic.certLocalMin, 'empty history must NOT certify');
assert(contains(ic.verdict, 'NO-DATA'), 'empty-history verdict must read NO-DATA');

% opts.tailN / opts.tol are respected
ic = foc_ipopt_inertia([1e-9*ones(1,3), 0, 0], struct('tailN',2,'tol',1e-8));
assert(ic.certLocalMin, 'tailN=2 window with clean values must certify');
assert(ic.maxTailDw == 0, 'maxTailDw should be the tail max (0 here)');

ic = foc_ipopt_inertia([1e-6, 1e-6], struct('tailN',5,'tol',1e-8));
assert(~ic.certLocalMin, 'tailN longer than history clamps to available length, still fails on nonzero');

fprintf('test_foc_ipopt_inertia: ALL PASS\n');
