function Tv = hs_quad_ctrl(tt, U, Um, h, N)
% HS_QUAD_CTRL  Exact per-segment quadratic Lagrange control
% reconstruction through (U_k, Um_k, U_{k+1}) for segment k -- the
% control representation Hermite-Simpson's own defect equations (Simpson's
% rule) are actually built against within each segment.
%
% PROMOTED from certify/certify_pdg.m (task-7 fix report, 2026-08-08):
% originally a local function used only by certify_pdg's G2 gate. Task 6's
% tvlqr_design.m built ctrl.Tnom from a GLOBAL pchip spline across all
% nodes+midpoints instead -- exactly the representation certify_pdg's own
% G2 decisive experiment (task-5 fix report) disproved: a global spline
% does not reconstruct the same control HS's defect equations are built
% against, and measured a non-shrinking ~0.6-0.8 kg mass residual "floor"
% that the per-segment quadratic below collapsed to ~0.0001-0.0005 kg (see
% "G2 control reconstruction" note in certify_pdg.m). Sharing this single
% implementation (byte-identical logic to the former local copy, proven by
% test_certify_nominal.m still passing) means Task 6's TVLQR feedforward
% and Task 7's closed-loop truth sim now fly the SAME control
% reconstruction certify_pdg's G2 gate already certified as accurate to
% 1.6e-4 m, instead of a spline known to smear the bang-bang switch and
% front-load braking.
%
% INPUTS:
%   tt - query time [scalar, s] (ode45 calls this with scalar t)
%   U  - node control samples [3 x (N+1)]
%   Um - midpoint control samples [3 x N]
%   h  - segment duration [scalar, s] (= tf/N, uniform grid)
%   N  - number of segments [scalar]
% OUTPUTS:
%   Tv - reconstructed control at tt [3x1]
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-08-08-booster-landing-design.md
ttc = min(max(tt, 0), N*h);
k   = min(max(floor(ttc/h) + 1, 1), N);        % segment index 1..N
tau = (ttc - (k-1)*h) / h;                     % local var in [0,1]
L0  = 2*tau^2 - 3*tau + 1;                     % Lagrange basis at tau=0
L1  = -4*tau^2 + 4*tau;                        % at tau=0.5 (midpoint)
L2  = 2*tau^2 - tau;                           % at tau=1
Tv  = L0*U(:,k) + L1*Um(:,k) + L2*U(:,k+1);
end
