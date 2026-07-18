function iterCap = psr_iter_cap(floorIter, N)
% PSR_ITER_CAP  Pure per-round IPOPT max-iteration cap for psr_mee_refine.m
% (Task 9 Step 0 fix): the flat opts.maxIter=1500 default (Task 8) was seen
% hitting Maximum_Iterations_Exceeded at N~3000-3365 nodes, at defects only
% 1-2 orders of magnitude above the certify gate (5.7e-7/7.95e-7 vs <1e-8) --
% i.e. plausibly just short of budget, not diverging. Scaling the cap with
% the round's own node count (the brief's own suggested formula) gives later,
% bigger rounds proportionally more iterations instead of the same fixed
% budget every round.
%
% Factored into its own pure function (mirrors round_advance_decision.m /
% psr_switch_score_mee.m's testability rationale) so the scaling formula is
% unit-testable without a solve.
%
% INPUTS:
%   floorIter - opts.maxIter, the user-set FLOOR (never solve with fewer
%               iterations than this even at tiny N) [scalar, positive]
%   N         - node count (numel(sigma)-1) of the round about to be solved
%               [scalar, positive]
%
% OUTPUTS:
%   iterCap - max(floorIter, ceil(N)) [scalar]
%
% REFERENCES: [1] task-9-brief.md Step 0 item 1 ("scale opts.maxIter with N
%   (e.g. max(1500, ceil(N))"). [2] psr_mee_refine.m (sole caller).
iterCap = max(floorIter, ceil(N));
end
