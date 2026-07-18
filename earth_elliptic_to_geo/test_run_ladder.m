% TEST_RUN_LADDER  No-solve regression test for run_ladder.m (Fix 4, review
% finding):
%   (a) a non-strictly-descending thrustList throws run_ladder:notDescending
%   (b) run_ladder([10]) with the existing certified 10 N artifacts REUSES
%       them with NO solve -- reused=true, wall time < 5 s
%   (c) the rung struct has the expected field set
%
% Deliberately single-path / no-solve, matching the sibling guard tests'
% style (test_mintime_mee_guard.m, test_stall_guard.m) -- this test must
% never trigger an actual NLP solve; if MEE_ladder_T100.mat (or its
% MEE_mintime_T100.mat / MEE_M2_10N.mat prerequisites) is ever deleted, this
% test should fail loudly (via the wall-time assertion) rather than
% silently solving for minutes.
%
% REFERENCES: [1] run_ladder.m. [2] task-6-report.md (10 N leg reuse,
%   the artifacts this test exercises).

% --- (a) non-descending thrustList must throw run_ladder:notDescending ----
threw = false;
try
    run_ladder([5 10]);
catch ME
    threw = true;
    assert(strcmp(ME.identifier, 'run_ladder:notDescending'), ...
        'test_run_ladder: expected identifier run_ladder:notDescending, got %s', ME.identifier);
end
assert(threw, 'test_run_ladder: run_ladder([5 10]) should have thrown run_ladder:notDescending');

% --- (b) run_ladder([10]) with the existing certified 10 N artifacts must
% REUSE (no solve) -- fast and reused=true ----------------------------------
% cfg.summaryFile is namespaced TEST_-prefixed (final-review Fix 1): the
% default summaryFile='MEE_ladder.mat' is the real 4-rung [10 5 2.5 1]
% campaign summary -- without this override, every test run would clobber
% it with a 1-rung summary.
tStart = tic;
results = run_ladder([10], struct('summaryFile', 'TEST_MEE_ladder.mat'));
wallSec = toc(tStart);
assert(isscalar(results), 'test_run_ladder: expected a 1x1 results struct for a single-rung call');
assert(results.reused == true, ['test_run_ladder: expected results.reused=true (no-solve ' ...
    'reuse of the 10 N cache), got false -- did the cache get invalidated?']);
assert(wallSec < 5, ['test_run_ladder: run_ladder([10]) took %.2f s -- expected a cache ' ...
    'HIT (<5 s); a solve appears to have run'], wallSec);

% --- (c) rung struct field set ----------------------------------------------
expectedFields = {'thrustN', 'anchor', 'fuelTag', 'fuel', 'tf', 'certified', 'reused', 'fp'};
for k = 1:numel(expectedFields)
    assert(isfield(results, expectedFields{k}), ...
        'test_run_ladder: rung struct missing expected field ''%s''', expectedFields{k});
end
assert(results.thrustN == 10, 'test_run_ladder: expected thrustN=10, got %g', results.thrustN);
assert(results.certified == true, 'test_run_ladder: expected certified=true for the cached 10 N rung');
assert(isfield(results.anchor, 'tfmin'), 'test_run_ladder: rung.anchor missing .tfmin');
assert(isfield(results.fuel, 'm_f_kg'), 'test_run_ladder: rung.fuel missing .m_f_kg');

fprintf(['test_run_ladder: ALL PASS (notDescending guard fires; run_ladder([10]) reused the ' ...
    'certified artifacts with no solve in %.3f s; rung struct has the expected fields)\n'], wallSec);
