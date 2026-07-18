function tf = psr_is_global_round(measuredRound, globalEvery)
% PSR_IS_GLOBAL_ROUND  Pure decision: should the refinement step LEAVING this
% measured round use a periodic GLOBAL (uniform, switch-blind) mesh densify
% instead of the default switch-local refinement (Task 9 Step 0 hybrid-
% refinement fix)?
%
% MOTIVATION: psr_switch_score_mee.m only scores collocation intervals near
% a switch the CURRENT mesh has already resolved well enough to see (a
% throttle crossing on the primal U(4,:)>0.5 test). It is structurally BLIND
% to a switch pair hiding inside a still-too-coarse region with no visible
% crossing yet -- Task 8/9's own comparison found a uniform 1.5x refine at
% 1 N discovering sw=177 vs windowed PSR's flat 171 across multiple rounds.
% A periodic fully-global round (every opts.globalEvery-th measured round)
% densifies EVERYWHERE, giving any as-yet-invisible switch pair a chance to
% resolve into a visible crossing on the next resolve, while the intervening
% rounds stay cheap and switch-local as before.
%
% Factored into its own pure function (mirrors round_advance_decision.m's
% testability rationale) so the periodicity rule is unit-testable without a
% solve.
%
% INPUTS:
%   measuredRound - 0-based round index just measured (history(row).round)
%                   [scalar, integer >= 0]
%   globalEvery   - period in measured rounds [scalar; <=0 or empty disables
%                   global rounds entirely -- the DEFAULT, since Task 8's
%                   evidence for this gap is real but this port's baseline
%                   behavior should not change unless a caller opts in]
%
% OUTPUTS:
%   tf - true iff measuredRound > 0 AND mod(measuredRound, globalEvery) == 0
%        (round 0, the seed, is NEVER a global round -- there is nothing to
%        refine away from yet)
%
% REFERENCES: [1] task-9-brief.md Step 0 item 2 (hybrid periodic global
%   refinement mandate). [2] psr_global_score_mee.m (the scorer this
%   triggers). [3] psr_mee_refine.m (sole caller).
if nargin < 2 || isempty(globalEvery) || globalEvery <= 0
    tf = false;
    return;
end
tf = measuredRound > 0 && mod(measuredRound, globalEvery) == 0;
end
