% TEST_ANCHOR_SMALLN_FIRST — unit-level check of the two PURE helpers that
% drive anchor_smallN_first.m's Stage-1 manual relaxed-stall continuation
% (smallN_warmtight, smallN_maxiter), harvested verbatim from
% results/task7c_step1_manual.m. This test does NOT run a solve (no CasADi
% needed, no cd into a solver mode) -- it is a fast logic-only gate; the
% full anchor_smallN_first(T,par,warmAnchor,aopts) live solve (~6+ min) is
% exercised separately by Task 6's background validation, per the brief.
here = fileparts(mfilename('fullpath')); cd(here);

% (a) the continuation's warmTight predicate (status-class allowlist AND
% defect < 1e-4):
assert(smallN_warmtight('Maximum_Iterations_Exceeded', 1e-5) == true, ...
    'benign status + tight defect must warmTight');
assert(smallN_warmtight('Restoration_Failed', 1e-9) == false, ...
    'bad status must never warmTight even with a tiny defect');
assert(smallN_warmtight('Solve_Succeeded', 1e-3) == false, ...
    'benign status but defect not yet < 1e-4 must not warmTight');

% (b) the maxIter ramp: 75 default, 150 once defect < 1e-3
assert(smallN_maxiter(1e-2) == 75, 'maxIter should be 75 above the 1e-3 floor');
assert(smallN_maxiter(5e-4) == 150, 'maxIter should ramp to 150 once defect < 1e-3');

fprintf('test_anchor_smallN_first PASSED\n');
