function [keepNew, retryLoose] = round_advance_decision(prevDefect, newDefect, alreadyRetried, C)
% ROUND_ADVANCE_DECISION  Pure retain-if-improved / retry-loose-once decision
% for a single mintime continuation round (process/DESIGN_thrust_ladder.md Phase 0
% item 1: "retain the new primal iterate only if feasibility improved over
% the incoming one ... otherwise keep the prior iterate and retry with
% adaptive mu"). Factored out of run_mintime_mee.m's mintime_mee_continue so
% the retain/retry decision is unit-testable with synthetic numbers, without
% a solve (test_mintime_mee_guard.m).
%
% "IMPROVED" IS FLOOR-GATED, NOT A BARE INEQUALITY (re-review finding,
% 2026-07-17): the first cut of this function used improved = newDefect <
% prevDefect, i.e. ANY positive progress, however small, counted as
% "improved" and reset the retry budget. That silently defeated decadeMin:
% a sub-floor crawl (e.g. the guard's own case-2 data point, 0.1186 decades
% -- positive, but below the 0.15 floor) would return keepNew=true every
% single round and run all the way to roundsMax instead of ever reaching
% the stall guard. "Improved" now means the SAME decadeImprove computation
% the stall guard has always used clears C.decadeMin:
%   decadeImprove = log10(prevDefect) - log10(newDefect) >= C.decadeMin
% Anything below that floor -- a sub-floor-but-positive crawl OR an actual
% regression (newDefect >= prevDefect) -- is treated identically: NOT
% improved, spend the one loose retry, then let the caller's stall guard
% fire. (The caller separately decides, on a NOT-improved round, whether to
% retain outNew as the new baseline when it was at least numerically better
% than the incoming point -- "the better of prev/new by defect" -- even
% though sub-floor progress does not reset the retry budget; that retained-
% point bookkeeping lives in run_mintime_mee.m>mintime_mee_continue, not
% here, since this function only returns the two booleans.)
%
% Kept as its own single-purpose file rather than a nested subfunction of
% run_mintime_mee.m: MATLAB has no supported way to call a true local
% function from an external test script, and the review's own ask ("unit-
% test THAT in the guard test") requires exactly this. It is still "local"
% in the sense that matters here -- small, pure, single caller
% (mintime_mee_continue), no state of its own.
%
% INPUTS:  prevDefect     - incoming (currently retained) iterate's
%                           maxDefect [scalar]
%          newDefect      - candidate round's maxDefect [scalar]
%          alreadyRetried - true if the CURRENT retained iterate has already
%                           consumed a loose-regime retry without clearing
%                           the decadeMin floor [logical]
%          C              - struct with .decadeMin (minimum log10 defect
%                           improvement to count as "improved" -- the SAME
%                           value/formula the stall guard uses, so this
%                           function and the stall guard can never
%                           disagree about what counts as progress) and
%                           .maxLooseRetries (loose retries allowed per
%                           stuck iterate before the stall guard is allowed
%                           to fire), both from mintime_guard_constants
%                           [struct]
%
% OUTPUTS: keepNew    - true: the candidate cleared the decadeMin floor:
%                       retain it, reset the retry budget
%          retryLoose - true: the candidate did NOT clear the floor (sub-
%                       floor crawl OR true regression), but a loose-retry
%                       budget remains: retry once more (forced
%                       warmTight=false) before the caller's existing stall
%                       guard is allowed to fire
%
% REFERENCES: [1] process/DESIGN_thrust_ladder.md Phase 0 item 1 (the mandate).
%   [2] mintime_guard_constants.m (decadeMin/maxLooseRetries, shared).
%   [3] test_stall_guard.m / test_mintime_mee_guard.m (decadeImprove
%   arithmetic and floor value this must match exactly).
%   [4] run_mintime_mee.m>mintime_mee_continue (sole caller).
decadeImprove = log10(max(prevDefect, realmin)) - log10(max(newDefect, realmin));
if decadeImprove >= C.decadeMin
    keepNew    = true;
    retryLoose = false;
else
    keepNew    = false;
    retryLoose = ~alreadyRetried && C.maxLooseRetries >= 1;
end
end
