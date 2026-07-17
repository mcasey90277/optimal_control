function [keepNew, retryLoose] = round_advance_decision(prevDefect, newDefect, alreadyRetried, C)
% ROUND_ADVANCE_DECISION  Pure retain-if-improved / retry-loose-once decision
% for a single mintime continuation round (DESIGN_thrust_ladder.md Phase 0
% item 1: "retain the new primal iterate only if feasibility improved over
% the incoming one ... otherwise keep the prior iterate and retry with
% adaptive mu"). Factored out of run_mintime_mee.m's mintime_mee_continue so
% the retain/retry decision is unit-testable with synthetic numbers, without
% a solve (test_mintime_mee_guard.m).
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
%                           consumed a loose-regime retry without improving
%                           [logical]
%          C              - struct with .maxLooseRetries (from
%                           mintime_guard_constants), the number of loose
%                           retries allowed per stuck iterate before the
%                           stall guard is allowed to fire [struct]
%
% OUTPUTS: keepNew    - true: the candidate improved feasibility over the
%                       retained iterate: retain it, reset the retry budget
%          retryLoose - true: the candidate did NOT improve, but a loose-
%                       retry budget remains: keep the prior iterate and
%                       retry once more (forced warmTight=false) before the
%                       caller's existing stall guard is allowed to fire
%
% REFERENCES: [1] DESIGN_thrust_ladder.md Phase 0 item 1 (the mandate).
%   [2] mintime_guard_constants.m (maxLooseRetries, the shared constant).
%   [3] run_mintime_mee.m>mintime_mee_continue (sole caller).
improved = newDefect < prevDefect;
if improved
    keepNew    = true;
    retryLoose = false;
else
    keepNew    = false;
    retryLoose = ~alreadyRetried && C.maxLooseRetries >= 1;
end
end
