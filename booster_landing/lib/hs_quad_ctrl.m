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
%% clipped to the annulus -- this is what guarantees feasibility. EXCEPT
%% on the one bang-bang SWITCH segment, where it is a STEP (see below).
m0 = sqrt(sum(U(:,k).^2));
mm = sqrt(sum(Um(:,k).^2));
m1 = sqrt(sum(U(:,k+1).^2));

%% SWITCH-SEGMENT STEP (task-7b): a min-fuel |T*| is bang-bang -- G5
%% certifies exactly one interior switch -- and a quadratic through a
%% STEP's three samples is simply the wrong basis function. Measured
%% consequence: the ENTIRE G2 continuous residual is generated inside this
%% single segment (mass error 0.000 kg before t=5 s, +1.45 kg at t=6.5 s,
%% -1.34 kg from t=7.0 s onward and flat to touchdown), because the
%% smeared ramp burns a different amount of propellant across the switch
%% than the bang-bang the NLP actually solved for. This was NOT introduced
%% by the P.etaT de-rate -- it is a pre-existing grid lottery, and the
%% shipped (etaT=1, N=60) combination was simply LUCKY: identical physics
%% at N=80 gives G2_pos = 0.93 m against a 1 m gate, and the de-rated
%% tf merely rolled the dice differently (1.29 m, FAIL).
%%
%% Fix: on a TRANSITION segment (see the test below), reconstruct |T| as a
%% step between the certified bang-bang LEVELS A and B -- Tmin and Tmax,
%% ordered by the sample trend -- at the instant s that makes the step's
%% integral EQUAL the segment's own Simpson quadrature:
%%   s*A + (h-s)*B = (h/6)(m0 + 4*mm + m1) =: Isimp
%%   =>  s = (Isimp - h*B) / (A - B),
%% i.e. the reconstruction consumes exactly the propellant the NLP's own
%% quadrature charged it. Falls back to the quadratic if s leaves [0,h].
%% Feasibility is preserved by construction: both levels ARE the annulus
%% bounds, so G2ff cannot regress. Direction is untouched -- the primer
%% direction is continuous through a switch; only the magnitude steps.
tolB  = 1e-3;
onB   = @(v) v <= Tmin*(1+tolB) || v >= Tmax*(1-tolB);
anyOn = onB(m0) || onB(mm) || onB(m1);
%% "Transition segment" = not all three samples pinned to ONE bound. The
%% off-bound sample the NLP uses to encode the switch can land on either a
%% MIDPOINT or a NODE; when it lands on a node the switch straddles TWO
%% segments, which is why testing only for "endpoints on opposite bounds"
%% is not enough (it missed exactly the etaT=0.93/N=60 and etaT=1.00/N=80
%% cases, leaving them at G2_pos 1.29 m and 0.93 m). So the levels are
%% taken to be the certified bang-bang levels themselves, Tmin and Tmax,
%% with the ordering set by the sample trend:
allLo = m0 <= Tmin*(1+tolB) && mm <= Tmin*(1+tolB) && m1 <= Tmin*(1+tolB);
allHi = m0 >= Tmax*(1-tolB) && mm >= Tmax*(1-tolB) && m1 >= Tmax*(1-tolB);
%% A transition segment must (i) not be pinned to a single bound and
%% (ii) still TOUCH a bound. Condition (ii) is what keeps the step model
%% off smooth INTERIOR arcs: without it a segment running from 0.60*Tmax to
%% 0.64*Tmax -- entirely in the annulus interior, no switch anywhere near --
%% satisfies (i), and s lands inside [0,h], so it would be reconstructed as
%% a full Tmin->Tmax bang. Vacuous on this campaign's certified bang-bang
%% solution (G5: bound fraction 0.9917, one interior switch), but live the
%% moment drag or a pointing cone introduces a real interior arc.
isTrans = ~(allLo || allHi) && anyOn;
sSw = NaN;  A = Tmin;  B = Tmax;
if isTrans
    if m1 < m0, A = Tmax;  B = Tmin; end            % falling switch
    Isimp = (h/6)*(m0 + 4*mm + m1);                 % the NLP's own quadrature
    sSw   = (Isimp - h*B) / (A - B);                % s*A + (h-s)*B = Isimp
end
if isTrans && sSw >= 0 && sSw <= h
    if (ttc - (k-1)*h) <= sSw, Tmag = A; else, Tmag = B; end
else
    Tmag = L0*m0 + L1*mm + L2*m1;   % smooth/interior arc: quadratic as before
end
Tmag = min(max(Tmag, Tmin), Tmax);
%% NOTE (corrected, task-7b polish): an earlier version of this comment
%% claimed s falling outside [0,h] was what protected smooth interior arcs.
%% That was FALSE -- a reviewer demonstrated a 0.60->0.64*Tmax interior arc
%% whose s lands inside [0,h] and which the step model therefore rewrote as
%% a full Tmin->Tmax bang. The actual protection is the anyOn term in
%% isTrans above (the segment must touch a bound); the [0,h] range check is
%% a second, weaker guard, not the primary one.

Tv = dirT * Tmag;
end
