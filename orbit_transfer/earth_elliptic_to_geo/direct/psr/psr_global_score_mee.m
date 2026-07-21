function score = psr_global_score_mee(N, factor)
% PSR_GLOBAL_SCORE_MEE  Pure per-interval score for a periodic GLOBAL
% (uniform, switch-blind) PSR densify round (Task 9 Step 0 hybrid-refinement
% fix -- see psr_is_global_round.m for the round-selection rule and
% motivation). Unlike psr_switch_score_mee.m (local, peaked at KNOWN
% switches), this scorer is switch-agnostic: it selects an EVENLY SPACED
% subset of collocation intervals spanning the WHOLE mesh, sized so that
% roughly (factor-1)*N new nodes get inserted -- e.g. factor=1.3 densifies
% the full mesh by ~30%, spread uniformly, giving any switch pair the
% windowed scorer hasn't discovered yet a chance to resolve into a visible
% throttle crossing on the next resolve.
%
% Factored into its own pure function (mirrors psr_switch_score_mee.m's
% testability rationale) so the selection is unit-testable without a solve.
%
% INPUTS:
%   N      - number of collocation intervals (numel(sigma)-1) [scalar, >=1]
%   factor - global densify factor, e.g. 1.3 -> ~30% more nodes [scalar > 1]
%
% OUTPUTS:
%   score - per-interval score [1xN], 1 at round((factor-1)*N) EVENLY SPACED
%           selected intervals (via linspace(1,N,...), deduplicated), 0
%           elsewhere -- directly consumable by psr_refine_sigma_mee.m
%           exactly like psr_switch_score_mee.m's output (that function's
%           K=Inf default selects every score>0 interval, so this global
%           score plugs into the SAME bisection call unchanged).
%
% REFERENCES: [1] task-9-brief.md Step 0 item 2 ("uniform densify by 1.3x").
%   [2] psr_is_global_round.m (the round-selection rule that picks this
%   scorer over psr_switch_score_mee.m). [3] psr_refine_sigma_mee.m (the
%   consumer). [4] psr_mee_refine.m (sole caller).
assert(N >= 1, 'psr_global_score_mee:badN', 'N must be >= 1, got %g', N);
assert(factor > 1, 'psr_global_score_mee:badFactor', ...
    'factor must be > 1 (a densify factor), got %g', factor);

targetAdd = min(N, max(1, round((factor - 1) * N)));
idx = unique(round(linspace(1, N, targetAdd)));

score = zeros(1, N);
score(idx) = 1;
end
