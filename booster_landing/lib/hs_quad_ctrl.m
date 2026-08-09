function Tv = hs_quad_ctrl(tt, U, Um, h, N, Tmin, Tmax)
% HS_QUAD_CTRL  Annulus-feasible per-segment quadratic control
% reconstruction through (U_k, Um_k, U_{k+1}) for segment k -- the control
% representation Hermite-Simpson's own defect equations (Simpson's rule)
% are actually built against within each segment.
%
% PROMOTED from certify/certify_pdg.m (task-7 fix report, 2026-08-08):
% originally a local function used only by certify_pdg's G2 gate. Task 6's
% tvlqr_design.m built ctrl.Tnom from a GLOBAL pchip spline across all
% nodes+midpoints instead -- exactly the representation certify_pdg's own
% G2 decisive experiment (task-5 fix report) disproved: a global spline
% does not reconstruct the same control HS's defect equations are built
% against, and measured a non-shrinking ~0.6-0.8 kg mass residual "floor"
% that the per-segment quadratic collapsed to ~0.0001-0.0005 kg (see "G2
% control reconstruction" note in certify_pdg.m). Sharing this single
% implementation means Task 6's TVLQR feedforward and Task 7's closed-loop
% truth sim fly the SAME control reconstruction certify_pdg's G2 gate
% certifies, instead of a spline known to smear the bang-bang switch and
% front-load braking.
%
% ADAPTATION (task-7 fix report round 4, 2026-08-08): direction and
% magnitude are now reconstructed SEPARATELY, magnitude clipped to
% [Tmin,Tmax] -- was a single quadratic Lagrange interpolant of the
% vector T directly (L0*U_k + L1*Um_k + L2*U_{k+1}). That single-vector
% form is NOT guaranteed annulus-feasible between nodes even when all
% three interpolated points sit exactly ON the annulus: interpolating a
% vector whose DIRECTION also turns across the segment shrinks the
% interpolated magnitude below the true |T| (a chord-vs-arc effect, since
% |a*p+b*q| <= a*|p|+b*|q| for unit-sum weights a+b<=... here L0+L1+L2=1
% but the L_i are not all nonnegative -- L1 can exceed 1 near tau~0.5,
% making this a genuine extrapolation in the vector's own component
% directions, not a plain convex combination). Measured at the N=30 test
% grid: 11.7% of the flight sat below Tmin (worst dip -18%), the
% closed-loop sim's saturating clamp turned that into a persistent
% +0.6-1+ m altitude bias baked into the FEEDFORWARD itself (not a
% tracking error TVLQR could ever correct, since it is present even with
% zero feedback), and two rounds of TVLQR weight tuning were, in
% hindsight, buying feedback gain to cancel this feedforward defect
% rather than fixing it. Reconstructing |T| as its OWN per-segment
% quadratic (same L0,L1,L2 basis, applied to the three scalar magnitudes
% instead of the vector) and clamping that scalar to [Tmin,Tmax] before
% re-applying the (separately reconstructed, unit-norm) direction cannot
% produce an infeasible magnitude by construction. Reduces to
% (numerically very close to) the prior single-vector formula wherever
% the thrust direction is locally near-constant across a segment (the
% common case away from the bang-bang switch) -- verified unchanged
% G2 residual at the N=60 production grid (test-7 fix report round 4).
%
% INPUTS:
%   tt        - query time [scalar, s] (ode45 calls this with scalar t)
%   U         - node control samples [3 x (N+1)]
%   Um        - midpoint control samples [3 x N]
%   h         - segment duration [scalar, s] (= tf/N, uniform grid)
%   N         - number of segments [scalar]
%   Tmin,Tmax - thrust annulus bounds [scalar, N] (booster_params P.Tmin,
%               P.Tmax)
% OUTPUTS:
%   Tv - reconstructed control at tt [3x1], |Tv| in [Tmin,Tmax] always
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-08-08-booster-landing-design.md
ttc = min(max(tt, 0), N*h);
k   = min(max(floor(ttc/h) + 1, 1), N);        % segment index 1..N
tau = (ttc - (k-1)*h) / h;                     % local var in [0,1]
L0  = 2*tau^2 - 3*tau + 1;                     % Lagrange basis at tau=0
L1  = -4*tau^2 + 4*tau;                        % at tau=0.5 (midpoint)
L2  = 2*tau^2 - tau;                           % at tau=1

%% Direction: quadratic vector interpolant, then normalized:
Traw = L0*U(:,k) + L1*Um(:,k) + L2*U(:,k+1);
Traw_mag = sqrt(sum(Traw.^2));
dirT = Traw / max(Traw_mag, 1e-9);

%% Magnitude: INDEPENDENT quadratic of |T| at the same three samples,
%% clipped to the annulus -- this is what guarantees feasibility:
m0 = sqrt(sum(U(:,k).^2));
mm = sqrt(sum(Um(:,k).^2));
m1 = sqrt(sum(U(:,k+1).^2));
Tmag = L0*m0 + L1*mm + L2*m1;
Tmag = min(max(Tmag, Tmin), Tmax);

Tv = dirT * Tmag;
end
