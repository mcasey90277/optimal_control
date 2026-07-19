% TEST_REPRODUCE_COLLECT  Unit test for reproduce_table3_collect.m against
% the already-cached 10 N REPRO row (results/repro/REPRO_row_T100.mat, from
% Task 2's from-scratch validation smoke test). Pure reader test: does NOT
% solve anything (no reproduce_row call here).
%
% Requires results/repro/REPRO_row_T100.mat to exist; if it does not (e.g. a
% fresh checkout that has never run reproduce_row(10)), SKIPS with a clear
% message rather than failing -- this test is about reproduce_table3_collect.m's
% own reading/printing logic, not about re-earning that cache.
%
% Run:
%   matlab -batch "run('/abs/path/earth_elliptic_to_geo/test_reproduce_collect.m')"
%
% REFERENCES:
%   [1] .superpowers/sdd/task-4-brief.md (Step 1's exact test).
here = fileparts(mfilename('fullpath')); cd(here);

rowfile = fullfile(here, 'results', 'repro', 'REPRO_row_T100.mat');
if isfile(rowfile)
    tbl = reproduce_table3_collect([10]);
    assert(isscalar(tbl), 'reproduce_table3_collect([10]) must return exactly one row');
    assert(tbl(1).thrustN == 10, 'returned row must be the T=10 N rung');
    cert = table3_certified(10);
    assert(tbl(1).m_f_kg >= cert.m_f_kg - 0.5, ...
        sprintf('T=10 N: m_f=%.4f must be >= campaign floor %.4f - 0.5', ...
        tbl(1).m_f_kg, cert.m_f_kg));

    % MISSING-rung path: a thrust level with no cached REPRO row must print
    % MISSING and be OMITTED from the return (not error, not fabricate a row).
    tblMissing = reproduce_table3_collect([10, 5]);
    assert(isscalar(tblMissing) && tblMissing(1).thrustN == 10, ...
        'a rung with no cached REPRO_row file must be omitted from the return, not error');

    fprintf('test_reproduce_collect PASSED\n');
else
    fprintf('test_reproduce_collect SKIPPED (no results/repro/REPRO_row_T100.mat yet -- run reproduce_row(10) first)\n');
end
