function [swIdx, score] = psr_switch_score_mee(sigma, U, opts)
% PSR_SWITCH_SCORE_MEE  Primal throttle-crossing switch detector + per-
% interval mesh-refinement score for the MEE/L-domain PSR port. ADAPTED
% (not ported) from NLP_lowThrust_GTO_tulip/PSR/lib/pmp_refine_indicator.m's
% dual-based switching-function S(tau) scorer: that scorer recovers costates
% from the direct solve's own KKT defect duals and localizes a switch to a
% sub-interval S=0 root. This file instead scores directly off the PRIMAL
% throttle U(4,:) crossing 0.5 -- the same crossings casadi_lt_mee.m's own
% out.switches already counts (sum(abs(diff(thr>0.5)))). Porting the
% dual-recovery pipeline (sms_seed_duals's mode-'d' dual->costate map) to the
% MEE formulation is out of this task's scope per the brief, which explicitly
% accepts the primal-throttle version for this port. No sub-interval root
% localization is attempted here (a bang-bang thr in {0,1} carries no
% between-node information to interpolate against, unlike a continuous S).
%
% Factored out as its own pure function (mirrors interp_warmstart.m's
% rationale) so the mesh-insertion logic is unit-testable without a solve.
%
% INPUTS:
%   sigma - current node grid, 0->1 [(N+1)x1]
%   U     - controls [4x(N+1)], row 4 = throttle in [0,1]
%   opts  - struct: .nbr switch-window half-width in COLLOCATION INTERVALS
%           [default 2]
%
% OUTPUTS:
%   swIdx - collocation-interval indices (1..N) whose node pair straddles a
%           thr=0.5 crossing [1 x nsw], identical bracketing to
%           casadi_lt_mee.m's out.switches
%   score - per-interval refinement score [1xN], >= 0: a linear taper peaked
%           at each switch's own interval (weight 1) decaying to 0 at
%           +-(nbr+1) intervals away, weighted by the interval's own
%           normalized sigma-width (a wide interval near a switch is
%           worse-localized than a narrow one -- same intuition as the
%           tulip scorer's (hInt/sigf) factor)
%
% REFERENCES:
%   [1] NLP_lowThrust_GTO_tulip/PSR/lib/pmp_refine_indicator.m (dual-based
%       analog this file stands in for -- read-only source, not modified).
%   [2] earth_elliptic_to_geo/casadi_lt_mee.m (out.switches, same bracketing
%       convention).
%   [3] .superpowers/sdd/task-8-brief.md (this task's spec; ADAPTED choice
%       documented in psr_mee_refine.m's header).
if nargin < 3, opts = struct(); end
if ~isfield(opts, 'nbr') || isempty(opts.nbr), opts.nbr = 2; end

sigma = sigma(:);  N = numel(sigma) - 1;
thr  = U(4, :);
burn = thr > 0.5;
swIdx = find(diff(double(burn)) ~= 0);   % [1 x nsw], values in 1..N

score = zeros(1, N);
h    = diff(sigma).';                    % [1xN]
sigf = sigma(end) - sigma(1);
nbr  = opts.nbr;
for q = 1:numel(swIdx)
    kc = swIdx(q);
    for kk = max(1, kc - nbr):min(N, kc + nbr)
        w = 1 - abs(kk - kc) / (nbr + 1);         % linear taper, 1 at kc, 0 at +-(nbr+1)
        score(kk) = score(kk) + (h(kk) / sigf) * w;
    end
end
end
