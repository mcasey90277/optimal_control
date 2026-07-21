function maxIterR = smallN_maxiter(maxDefect)
% SMALLN_MAXITER  IPOPT maxIter ramp for the smallN_first (1 N) manual
% relaxed-stall continuation (anchor_smallN_first.m, Task 3), harvested
% VERBATIM from results/task7c_step1_manual.m: 75 iterations per round by
% default (cheap at the low nodes-per-rev density this stage runs at),
% raised to 150 once the incoming warm-start point's defect is already
% below 1e-3 (close enough to certification that a larger per-round budget
% is worth the cost).
%
% INPUTS:
%   maxDefect - incoming warm-start point's measured max defect [scalar]
%
% OUTPUTS:
%   maxIterR - IPOPT maxIter for the next continuation round [scalar]
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/results/task7c_step1_manual.m (source lines:
%       maxIterR = 75; if out.maxDefect < 1e-3, maxIterR = 150; end).
%   [2] .superpowers/sdd/task-3-brief.md (this task's spec).
maxIterR = 75;
if maxDefect < 1e-3
    maxIterR = 150;
end
end
