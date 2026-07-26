function lam = foc_dual_to_costate(LamDef, sigma)
% FOC_DUAL_TO_COSTATE  Interval defect duals -> nodal costate, step-weighted.
%
% Trapezoid defect constraints couple node k and node k+1:
%   X(:,k+1) - X(:,k) - (dsig(k)/2)*dL*(dXdL(:,k)+dXdL(:,k+1)) == 0,  k=1..N
% so out.lamDef(:,k) is the multiplier of INTERVAL k (living between nodes k
% and k+1), not a per-node quantity. The continuous costate lambda(sigma) is
% sampled AT the nodes; the standard (alpha-stationarity) discrete-adjoint
% result -- adopted campaign-wide per process/DESIGN_dual_map.md's "[CORRECTNESS]"
% clause and the ms_band precedent (MS_BAND_CAMPAIGN.md) -- is that the
% interior nodal costate is the STEP-WEIGHTED AVERAGE of its two adjacent
% interval duals (reduces to the plain average on a uniform mesh), one-sided
% at the two endpoints (only one adjacent interval exists there):
%
%   lam(:,1)   = Lam(:,1)                                    (one-sided)
%   lam(:,k)   = (h(k-1)*Lam(:,k-1) + h(k)*Lam(:,k)) / (h(k-1)+h(k))   1<k<N+1
%   lam(:,N+1) = Lam(:,N)                                    (one-sided)
%
% with h(k) = sigma(k+1)-sigma(k) the interval's own sigma-width. This map is
% NOT in dispute for transcriptions without per-node slack states. For solvers
% where the row-scaling suspect exists (e.g. cScale-augmented 9-state solvers),
% see the referenced DESIGN_dual_map triage.
%
% INPUTS:
%   LamDef - interval defect-constraint duals [nx x N] (nx = state dimension;
%            N = number of collocation intervals)
%   sigma  - node grid, monotonic increasing [(N+1)x1] or [1x(N+1)]
%
% OUTPUTS:
%   lam  - nodal costate, step-weighted adjacent-interval average, ONE-SIDED
%          at the two endpoint columns [nx x (N+1)]. This is the STATIONARITY
%          combination: the object the discrete adjoint recursion pairs with
%          the node-k Jacobian. Use it for anything derived from stationarity.
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/process/DESIGN_dual_map.md sec "[CORRECTNESS]" (the
%       mandatory weighted-average formula, campaign-wide).
%   [2] GTO_tulip/indirect/ms_band/MS_BAND_CAMPAIGN.md (adjudicated
%       midpoint dual map, the method precedent).
%   [3] verify_common/doc/review_2026-07-25_{gpt56terra,gemini31pro}.md
%       (external review: both reviewers independently derived the
%       stationarity map above. NOTE: they also flagged that the ONE-SIDED
%       endpoint columns carry an O(h) representation offset. Do NOT patch
%       that here -- endpoint-VALUED conditions (free-mass transversality)
%       are handled correctly in foc_check.m section (6) via the mapped
%       terminal covector, read off the Lagrangian gradient. An earlier
%       extrapolation fix in this file was superseded by that and removed.)

sigma = sigma(:).';
N  = size(LamDef, 2);
assert(numel(sigma) == N + 1, 'foc_dual_to_costate: need N+1 sigma nodes');
h  = diff(sigma);
assert(all(h > 0), 'foc_dual_to_costate: sigma must be strictly increasing');
nx = size(LamDef, 1);
lam = zeros(nx, N + 1);
lam(:, 1) = LamDef(:, 1);   lam(:, N + 1) = LamDef(:, N);
for k = 2:N
    lam(:, k) = (h(k-1)*LamDef(:, k-1) + h(k)*LamDef(:, k)) / (h(k-1) + h(k));
end
end
