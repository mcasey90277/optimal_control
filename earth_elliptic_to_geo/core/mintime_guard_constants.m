function [roundsMax, decadeMin, maxLooseRetries] = mintime_guard_constants()
% MINTIME_GUARD_CONSTANTS  Single source of the anchor continuation guard triple.
%
% Shared by run_mintime (all continuation branches), run_mintime_mee, and
% test_stall_guard/test_mintime_mee_guard so the tests pin the PRODUCTION
% constants (review finding: a hardcoded copy in the test, or in a solver
% driver, could drift from the values actually wired into the solver).
%
% INPUTS:  none
% OUTPUTS: roundsMax       - max continuation rounds [scalar]
%          decadeMin       - minimum log10 defect improvement per round [scalar]
%          maxLooseRetries - per-round-regression retry budget [scalar]: when
%                            a continuation round's new iterate does NOT
%                            improve feasibility over the incoming one, the
%                            prior iterate is retained and retried this many
%                            times with the loose (adaptive-mu) regime before
%                            the stall guard is allowed to fire (review
%                            finding, process/DESIGN_thrust_ladder.md Phase 0 item 1:
%                            "retain the new primal iterate only if
%                            feasibility improved ... otherwise keep the
%                            prior iterate and retry with adaptive mu").
%
% REFERENCES: [1] task-14-report.md (guard recalibration record).
%   [2] process/DESIGN_thrust_ladder.md Phase 0 item 1 (retain-if-improved / retry-
%       loose mandate, origin of maxLooseRetries).
roundsMax = 24;
decadeMin = 0.15;
maxLooseRetries = 1;
end
