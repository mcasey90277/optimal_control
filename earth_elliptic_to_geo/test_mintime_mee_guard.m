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

assert(contains(src, '[roundsMax, decadeMin, maxLooseRetries] = mintime_guard_constants()'), ...
    'test_mintime_mee_guard: run_mintime_mee.m must call mintime_guard_constants() to obtain roundsMax/decadeMin/maxLooseRetries');

% Guard against a future edit reintroducing a hardcoded copy of the current
% production values (24 rounds / 0.15 decades / 1 loose retry) as a literal
% assignment -- e.g. "roundsMax = 24;" -- which would silently drift from
% the shared file if mintime_guard_constants.m is ever recalibrated again.
% (Fix 5, review finding: also cover the retry-loose constant added for the
% per-round keep-if-improved / retry-loose mechanism, Phase 0 item 1.)
badPatterns = {'roundsMax\s*=\s*24\s*;', 'decadeMin\s*=\s*0\.15\s*;', 'maxLooseRetries\s*=\s*1\s*;'};
for k = 1:numel(badPatterns)
    assert(isempty(regexp(src, badPatterns{k}, 'once')), ...
        'test_mintime_mee_guard: run_mintime_mee.m appears to hardcode a local copy of the guard constants (pattern ''%s'' matched) -- it must only read mintime_guard_constants()', ...
        badPatterns{k});
end

% --- (2) production values are exactly what test_stall_guard.m pins -------
[roundsMax, decadeMin, maxLooseRetries] = mintime_guard_constants();
assert(roundsMax == 24, 'test_mintime_mee_guard: expected roundsMax=24, got %d', roundsMax);
assert(abs(decadeMin - 0.15) < 1e-12, 'test_mintime_mee_guard: expected decadeMin=0.15, got %.4f', decadeMin);
assert(maxLooseRetries == 1, 'test_mintime_mee_guard: expected maxLooseRetries=1, got %d', maxLooseRetries);

% --- (3) guard arithmetic itself, same two observed data points as
% test_stall_guard.m (the MEE continuation loop's decadeImprove formula is
% identical: log10(prevDefect) - log10(newDefect)) --------------------------
d1 = log10(9.353e-03) - log10(5.336e-03);     % must NOT stall
assert(d1 >= decadeMin, 'case 1 (%.4f decades) should NOT stall at floor %.2f', d1, decadeMin);

d2 = log10(5.141e-03) - log10(3.912e-03);     % must stall
assert(d2 < decadeMin, 'case 2 (%.4f decades) SHOULD stall at floor %.2f', d2, decadeMin);

% --- (4) round_advance_decision.m: the retain-if-improved / retry-loose
% decision (DESIGN_thrust_ladder.md Phase 0 item 1, Fix 1) unit-tested with
% synthetic defect numbers -- no solve required. round_advance_decision.m
% is kept as its own file (not a nested local function of run_mintime_mee.m)
% specifically so it can be called directly here. -----------------------
guardC = struct('maxLooseRetries', maxLooseRetries);

% Case A: improved (newDefect < prevDefect) -> keep the new iterate, no retry.
[keepNew, retryLoose] = round_advance_decision(5.0e-03, 2.0e-03, false, guardC);
assert(keepNew == true && retryLoose == false, ...
    'round_advance_decision case A (improved): expected keepNew=true, retryLoose=false, got keepNew=%d retryLoose=%d', ...
    keepNew, retryLoose);

% Case B: regressed-first-time (newDefect >= prevDefect, no retry consumed
% yet) -> reject the new iterate, spend the one loose-retry budget.
[keepNew, retryLoose] = round_advance_decision(2.0e-03, 3.5e-03, false, guardC);
assert(keepNew == false && retryLoose == true, ...
    'round_advance_decision case B (regressed-first-time): expected keepNew=false, retryLoose=true, got keepNew=%d retryLoose=%d', ...
    keepNew, retryLoose);

% Case B': regressed-first-time at exact equality (newDefect == prevDefect,
% i.e. no improvement) is also a regression, not an improvement.
[keepNew, retryLoose] = round_advance_decision(2.0e-03, 2.0e-03, false, guardC);
assert(keepNew == false && retryLoose == true, ...
    'round_advance_decision case B'' (no-change treated as regression): expected keepNew=false, retryLoose=true, got keepNew=%d retryLoose=%d', ...
    keepNew, retryLoose);

% Case C: regressed-after-retry (newDefect >= prevDefect, the one loose
% retry has already been consumed) -> reject, no more retries: the caller's
% existing stall guard is what fires next.
[keepNew, retryLoose] = round_advance_decision(2.0e-03, 3.5e-03, true, guardC);
assert(keepNew == false && retryLoose == false, ...
    'round_advance_decision case C (regressed-after-retry): expected keepNew=false, retryLoose=false, got keepNew=%d retryLoose=%d', ...
    keepNew, retryLoose);

fprintf(['test_mintime_mee_guard: ALL PASS (run_mintime_mee.m reads the shared ' ...
    'mintime_guard_constants() -- roundsMax=%d, decadeMin=%.2f, maxLooseRetries=%d; ' ...
    'case1=%.4f decades continue, case2=%.4f decades stall; round_advance_decision ' ...
    'improved/regressed-first-time/no-change/regressed-after-retry all correct)\n'], ...
    roundsMax, decadeMin, maxLooseRetries, d1, d2);
