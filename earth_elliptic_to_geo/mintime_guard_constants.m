function [roundsMax, decadeMin] = mintime_guard_constants()
% MINTIME_GUARD_CONSTANTS  Single source of the anchor continuation guard pair.
%
% Shared by run_mintime (all continuation branches) and test_stall_guard so the
% test pins the PRODUCTION constants (review finding: a hardcoded copy in the
% test could drift from the values actually wired into the solver).
%
% INPUTS:  none
% OUTPUTS: roundsMax - max continuation rounds [scalar]
%          decadeMin - minimum log10 defect improvement per round [scalar]
%
% REFERENCES: [1] task-14-report.md (guard recalibration record).
roundsMax = 24;
decadeMin = 0.15;
end
