function tf = smallN_warmtight(ipoptStatus, maxDefect)
% SMALLN_WARMTIGHT  Feasibility-selected warmTight predicate for the
% smallN_first (1 N) manual relaxed-stall continuation (anchor_smallN_first.m,
% Task 3), harvested VERBATIM from results/task7c_step1_manual.m: the next
% continuation round is allowed to run warmTight=true iff the INCOMING
% warm-start point (a) exited via a benign IPOPT status class -- Solve_
% Succeeded, Solved_To_Acceptable_Level, or Maximum_Iterations_Exceeded --
% never a restoration-phase or infeasibility abort -- AND (b) its measured
% defect is already tight (< 1e-4). A point exiting via a bad status (e.g.
% Restoration_Failed_!, Infeasible_Problem_Detected) always gets loose/
% adaptive-mu treatment next round regardless of its self-reported defect,
% since a restoration-phase exit's duals cannot be trusted (casadi_lt_mee.m
% builds a fresh casadi.Opti() every call, so only this status-class gate --
% not warm dual state -- prevents trusting them).
%
% INPUTS:
%   ipoptStatus - IPOPT return_status of the incoming warm-start point [char]
%   maxDefect   - incoming point's measured max defect [scalar]
%
% OUTPUTS:
%   tf - true iff the NEXT continuation round should run warmTight=true [logical]
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/results/task7c_step1_manual.m (source line:
%       warmTight = any(strcmp(out.ipoptStatus, {'Solve_Succeeded', ...
%       'Solved_To_Acceptable_Level', 'Maximum_Iterations_Exceeded'})) && ...
%       out.maxDefect < 1e-4).
%   [2] earth_elliptic_to_geo/run_mintime_mee.m (benignStatus, the same
%       status-class allowlist shared with the automatic driver).
%   [3] .superpowers/sdd/task-3-brief.md (this task's spec).
tf = any(strcmp(ipoptStatus, {'Solve_Succeeded', 'Solved_To_Acceptable_Level', ...
    'Maximum_Iterations_Exceeded'})) && maxDefect < 1e-4;
end
