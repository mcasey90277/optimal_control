function reproduce_table3(thrustList)
% REPRODUCE_TABLE3  Thin IN-PROCESS convenience wrapper: reproduce every
% rung in thrustList from scratch (reproduce_row.m) then print the updated
% Table 3 (reproduce_table3_collect.m).
%
% *** RUNS IN ONE MATLAB PROCESS. *** This is convenient for the crash-free
% TOP rungs (10/5/2.5 N have not been observed to crash MATLAB in this
% campaign) or for quick iteration, where a single in-process loop is
% simplest. The DEEP/hang-prone rungs (1 N's smallN_first long multi-round
% continuation, 0.5 N's PSR-refinement budget, and any future seeded 0.2/
% 0.1 N rung) are known to be susceptible to the recurring UNCATCHABLE
% MEX/CasADi fatal crash documented throughout this campaign (e.g.
% run_task9_watchdog.sh's header comment, process/CAMPAIGN.md) -- a crash there
% takes this ENTIRE MATLAB process down, and this in-process loop does NOT
% auto-relaunch (rows already reproduced earlier in the SAME call are safe
% -- they were already written to results/repro/ -- but the loop itself
% stops dead and must be re-invoked by hand). For those rungs, use the
% per-process watchdog orchestrator reproduce_table3.sh instead, which
% isolates each rung in its own OS process and auto-relaunches across the
% crash (reproduce_row.m's own REPRO_-tagged driver caches make a relaunch
% resume close to where the dead process left off, not re-solve from
% absolute scratch).
%
% Every rung's reproduce_row(T) is called unconditionally (not skipped even
% if a results/repro/REPRO_row_T*.mat cache already exists) -- this
% function does not add its own resume-skip logic on top of
% reproduce_row.m's from-scratch-by-design contract; reproduce_row.m's
% internal REPRO_-tagged driver caches already make a repeat call cheap
% once a rung has converged.
%
% INPUTS:
%   thrustList - vector of thrust levels [N] to reproduce + print, in
%                order (e.g. [10 5 2.5 1 0.5])                           [1xK]
%
% OUTPUTS: none. Side effects: writes results/repro/REPRO_row_T*.mat for
%          each T (via reproduce_row.m) and prints the collected table
%          (via reproduce_table3_collect.m).
%
% REFERENCES:
%   [1] .superpowers/sdd/task-4-brief.md (this function's spec).
%   [2] earth_elliptic_to_geo/reproduce_row.m (the per-rung engine this
%       wrapper loops over).
%   [3] earth_elliptic_to_geo/reproduce_table3.sh (the crash-survivable,
%       per-process sibling of this convenience wrapper -- prefer it for
%       the deep rungs).
%   [4] earth_elliptic_to_geo/reproduce_table3_collect.m (the table
%       printer this wrapper calls at the end).

if nargin < 1 || isempty(thrustList)
    error('reproduce_table3:badInput', 'thrustList is required');
end

for T = thrustList
    reproduce_row(T);
end

reproduce_table3_collect(thrustList);

end
