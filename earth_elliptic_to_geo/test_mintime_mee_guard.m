% TEST_MINTIME_MEE_GUARD  No-solve pin: run_mintime_mee.m's continuation
% guard reads the SAME shared production constants as run_mintime.m
% (mintime_guard_constants.m), not a hardcoded local copy, and the guard
% arithmetic itself (decadeImprove = log10(prevDefect) - log10(newDefect),
% continue if decadeImprove >= decadeMin else stall) still holds against the
% same two observed data points test_stall_guard.m pins for the Cartesian
% driver -- the MEE continuation loop (mintime_mee_continue, local to
% run_mintime_mee.m) is a line-for-line port of the same arithmetic, so a
% regression in either file's guard math would show up here too.
%
% Task 6 review requirement: "unit-test the guard arithmetic reuse ...
% asserting the MEE driver reads the SAME production constants from
% mintime_guard_constants.m (no hardcoded copies) -- a no-solve test."

% --- (1) source-level check: run_mintime_mee.m calls the shared getter, and
% does not shadow it with a local hardcoded roundsMax/decadeMin pair --------
here = fileparts(mfilename('fullpath'));
srcFile = fullfile(here, 'run_mintime_mee.m');
assert(isfile(srcFile), 'test_mintime_mee_guard: run_mintime_mee.m not found at %s', srcFile);
src = fileread(srcFile);

assert(contains(src, '[roundsMax, decadeMin] = mintime_guard_constants()'), ...
    'test_mintime_mee_guard: run_mintime_mee.m must call mintime_guard_constants() to obtain roundsMax/decadeMin');

% Guard against a future edit reintroducing a hardcoded copy of the current
% production values (24 rounds / 0.15 decades) as a literal assignment --
% e.g. "roundsMax = 24;" -- which would silently drift from the shared file
% if mintime_guard_constants.m is ever recalibrated again.
badPatterns = {'roundsMax\s*=\s*24\s*;', 'decadeMin\s*=\s*0\.15\s*;'};
for k = 1:numel(badPatterns)
    assert(isempty(regexp(src, badPatterns{k}, 'once')), ...
        'test_mintime_mee_guard: run_mintime_mee.m appears to hardcode a local copy of the guard constants (pattern ''%s'' matched) -- it must only read mintime_guard_constants()', ...
        badPatterns{k});
end

% --- (2) production values are exactly what test_stall_guard.m pins -------
[roundsMax, decadeMin] = mintime_guard_constants();
assert(roundsMax == 24, 'test_mintime_mee_guard: expected roundsMax=24, got %d', roundsMax);
assert(abs(decadeMin - 0.15) < 1e-12, 'test_mintime_mee_guard: expected decadeMin=0.15, got %.4f', decadeMin);

% --- (3) guard arithmetic itself, same two observed data points as
% test_stall_guard.m (the MEE continuation loop's decadeImprove formula is
% identical: log10(prevDefect) - log10(newDefect)) --------------------------
d1 = log10(9.353e-03) - log10(5.336e-03);     % must NOT stall
assert(d1 >= decadeMin, 'case 1 (%.4f decades) should NOT stall at floor %.2f', d1, decadeMin);

d2 = log10(5.141e-03) - log10(3.912e-03);     % must stall
assert(d2 < decadeMin, 'case 2 (%.4f decades) SHOULD stall at floor %.2f', d2, decadeMin);

fprintf(['test_mintime_mee_guard: ALL PASS (run_mintime_mee.m reads the shared ' ...
    'mintime_guard_constants() -- roundsMax=%d, decadeMin=%.2f; case1=%.4f decades ' ...
    'continue, case2=%.4f decades stall)\n'], roundsMax, decadeMin, d1, d2);
