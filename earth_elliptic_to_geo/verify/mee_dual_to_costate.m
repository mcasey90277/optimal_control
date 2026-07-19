function lam = mee_dual_to_costate(LamDef, sigma)
% MEE_DUAL_TO_COSTATE  Interval defect duals -> nodal costate, step-weighted.
%
% casadi_lt_mee's trapezoid defect conDef{k} couples node k and node k+1:
%   X(:,k+1) - X(:,k) - (dsig(k)/2)*dL*(dXdL(:,k)+dXdL(:,k+1)) == 0,  k=1..N
% so out.lamDef(:,k) is the multiplier of INTERVAL k (living between nodes k
% and k+1), not a per-node quantity. The continuous costate lambda(sigma) is
% sampled AT the nodes; the standard (alpha-stationarity) discrete-adjoint
% result -- adopted campaign-wide per DESIGN_dual_map.md's "[CORRECTNESS]"
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
% NOT in dispute for this transcription: unlike casadi_lt_2body's 9-state
% cScale-augmented Sundman solver (open ~20 deg primer anomaly, believed
% extraction-side per DESIGN_dual_map.md's triage, NOT this averaging step),
% casadi_lt_mee.m has no per-node slack state at all -- DeltaL is a single
% free SCALAR (opti.variable(), one column, not a per-node row multiplying
% every dynamics row) -- so there is no analogous row-scaling suspect here.
% This file supplies the map anyway, correctly, from first principles, so
% Campaign-B's anomaly class cannot be silently imported into the MEE
% verifier by omission.
%
% INPUTS:
%   LamDef - interval defect-constraint duals [nx x N] (nx = state dim, 7 for
%            the MEE solver; N = number of collocation intervals)
%   sigma  - node grid, monotonic increasing [(N+1)x1] or [1x(N+1)]
%
% OUTPUTS:
%   lam - nodal costate, step-weighted adjacent-interval average [nx x (N+1)]
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/DESIGN_dual_map.md sec "[CORRECTNESS]" (the
%       mandatory weighted-average formula, campaign-wide).
%   [2] NLP_lowThrust_GTO_tulip/ms_band/MS_BAND_CAMPAIGN.md (adjudicated
%       midpoint dual map, the method precedent).
sigma = sigma(:).';                      % [1x(N+1)]
N  = size(LamDef, 2);
assert(numel(sigma) == N + 1, ...
    'mee_dual_to_costate: sigma must have N+1 = %d nodes for N = %d intervals', N+1, N);
h  = diff(sigma);                        % [1xN]

nx  = size(LamDef, 1);
lam = zeros(nx, N + 1);
lam(:, 1)     = LamDef(:, 1);
lam(:, N + 1) = LamDef(:, N);
for k = 2:N
    lam(:, k) = (h(k-1) * LamDef(:, k-1) + h(k) * LamDef(:, k)) / (h(k-1) + h(k));
end
end
