% TEST_STALL_GUARD  Unit check for the continuation stall-guard arithmetic in
% run_mintime.m>mintime_stage12_continue: decadeImprove = log10(prevDefect) -
% log10(newDefect), continue if decadeImprove >= decadeMin else stall/error.
% Task 14 controller triage round 4: a round-3 bug routed the 5 N cold-seed
% continuation through the WRONG roundsMax/decadeMin pair (3 rounds / 1.0
% decades, the 10 N-calibrated constants) instead of the intended 24 / 0.15,
% so a round that legitimately cleared the 0.15 floor was rejected against a
% floor of 1.0 instead. This pins the guard arithmetic itself against the two
% observed data points so it can be trusted independent of which
% roundsMax/decadeMin pair gets plumbed into which call site.
[~, decadeMin] = mintime_guard_constants();

% Case 1: observed 5 N cold-seed round 1 (9.353e-03 -> 5.336e-03) -- must NOT stall
d1 = log10(9.353e-03) - log10(5.336e-03);
assert(d1 >= decadeMin, 'case 1 (%.4f decades) should NOT stall at floor %.2f', d1, decadeMin);

% Case 2: observed 5 N anchor-seeded round 1 (5.141e-03 -> 3.912e-03) -- must stall
d2 = log10(5.141e-03) - log10(3.912e-03);
assert(d2 < decadeMin, 'case 2 (%.4f decades) SHOULD stall at floor %.2f', d2, decadeMin);

fprintf('test_stall_guard: ALL PASS (case1=%.4f decades continue, case2=%.4f decades stall)\n', d1, d2);
