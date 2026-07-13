# Sundman-regularized multiple shooting (Z1-Sun) — build doc

Self-contained spec for the Sundman-regularized arcs, the structural fix for
Z1's conditioning noise floor (ZTL_RESULTS "Lessons learned" #9). Read that +
`Z0_BUILD.md` first.

## 0. Why (the diagnosis this addresses)

The physical-time MS trust region (geodesic) converged the 75 mN anchor to
1.4e-6 but floored there: cond(Jc)~1e9 -> linear-solve noise ~5-7e-7 in the
stiff perigee-continuity directions. The root is the perigee gravity-gradient
singularity: in PHYSICAL time the costate rate dlambda_r/dt = -G'lambda_v ~
(1/r1^3) lambda_v blows up at perigee (r1 = distance to Earth), making the STM
entries huge and the Jacobian ill-conditioned. Node placement can't fix it
(the amplification is un-splittable). SUNDMAN regularization changes the
INTEGRATION VARIABLE so the perigee passage is stretched and the rates are
tamed -- exactly why the DIRECT side (PSR, pSund=1.5) is well-conditioned.

## 1. The reparametrization

Independent variable t -> sigma in [0,1], via dt/dsigma = tauF * kappa(y),
kappa = r1^pSund, r1 = ||r - r_Earth|| = ||[r_x+muStar; r_y; r_z]||. tauF is
the total Sundman length (a scalar UNKNOWN, set by the fixed-tf constraint).
Under this map the augmented state gains the physical time t (needed for the
terminal t=tf condition); the costates are UNCHANGED (still the physical-time
PMP costates -- we only reparametrize the integration, so every RHS is just
scaled by tauF*kappa):

    Y = [r(3); v(3); m; lam_r(3); lam_v(3); lam_m; t]  in R^15
    dY/dsigma = tauF * [ kappa * f(y) ; kappa ]

where f(y) = the physical ztl_eom RHS (14). The costate rate near perigee
becomes dlambda_r/dsigma ~ tauF * r1^pSund * (1/r1^3) = tauF * r1^(pSund-3) --
milder for any pSund>0, and a knob (try 1.5 like PSR first; 2-3 tames the
1/r^3 costate singularity more, an indirect-specific option).

The switching function S(y) = 1 - ||lam_v||c/m - lam_m and the regime automaton
(on/medium/off at S=+/-eps) are UNCHANGED (S is a function of y only); events
fire on S in sigma exactly as in physical time.

## 2. Unknowns and residual (multiple shooting in sigma)

Nodes at FIXED sigma_k = (k-1)/M, k=1..M+1 (uniform in sigma == perigee-dense
in physical time, automatically). Node states carry Y (15). tauF is one global
unknown.

    z = [ lam0(7) ; Y_2(15) ; ... ; Y_M(15) ; tauF(1) ]     dim 15M-7
    node 1 (sigma=0): Y_1 = [rv0; 1; lam0; 0]  (r,v,m,t fixed; lam0 unknown)

Arc k integrates Y_k from sigma_k to sigma_{k+1} (needs tauF) -> F_k, with
per-arc STM Phi_k = dF_k/dY_k (15x15) AND w_k = dF_k/dtauF (15x1).

    R = [ F_k - Y_{k+1}          (k=1..M-1, continuity, 15 each) ;
          F_M(1:6) - rvf ;       (rendezvous, 6)
          F_M(14) ;              (lam_m(sigma=1)=0, 1)
          F_M(15) - tf ]         (physical time closes, 1)          dim 15M-7

Square (15(M-1)+8 = 15M-7). Jacobian (block-bidiagonal + a dense tauF column):
  continuity block k: d/dlam0 = Phi_1(:,8:14) [k=1] or d/dY_k = Phi_k [k>=2];
                      d/dY_{k+1} = -I(15);  d/dtauF = w_k.
  terminal block: d/dY_M = [Phi_M(1:6,:); Phi_M(14,:); Phi_M(15,:)] (8x15);
                  d/dtauF = [w_M(1:6); w_M(14); w_M(15)] (8x1).

## 3. Files

- `ztl_eom_sun.m` [dY,aux] = f(Y,P,regime): dY/dtau = [kappa*ztl_eom(y); kappa]
  (returns UN-normalized Sundman-time RHS; the flow applies tauF). P.pSund.
  CS-safe (reuse ztl_eom for yDot with nargout=1; kappa=r1^p is CS-safe).
- `ztl_A_sun.m`  A = 15x15 d(dY/dtau)/dY by complex step of ztl_eom_sun.
- `ztl_flow_sun.m` out = f(Y0,tauF,[s0 s1],P,wantSTM): integrate
  dY/dsigma = tauF*(ztl_eom_sun RHS), and if wantSTM the 15x15 Phi
  (dPhi/dsigma = tauF*A*Phi) AND w=dYf/dtauF (dw/dsigma = tauF*A*w + G,
  G=ztl_eom_sun RHS -- inhomogeneous). 15+225+15=255 ODEs. Events on S(y),
  saltation at eps=0 (as ztl_flow, but with the Sundman kappa factored into
  Sdot; the saltation matrix uses dS/dy unchanged since S(y)).
- `ztl_ms_residual_sun.m`, `ztl_ms_seed_sun.m` (chop the physical Z1 anchor or
  a lam0: integrate in sigma to get node states + tauF; continuity ~0 by
  construction).
- `ztl_ms_solve_tr.m` GENERICIZED to call prob.resFun (default ztl_ms_residual)
  so the Sundman residual drops in with the same trust-region + geodesic solver.
- `test_ztl_sun.m` gates: (S1) dY/dtau = kappa*[f;1] vs physical at sample pts;
  (S2) ztl_A_sun vs FD-of-field; (S3) Phi and w vs FD-of-flow on a short arc;
  (S4) seed continuity ~0 and t(sigma=1) ~ tf at the seed tauF; (S5) COND CHECK
  -- cond(J_sun) vs cond(J_phys) at the anchor (the whole point: expect a big
  drop).
- `z1_run_sun.m` driver.

## 4. Seed and gate

Seed from the physical-time Z1 anchor (`z1_anchor_75mN.mat`, lam0_BE): its
trajectory in sigma gives the node states; tauF seed = integral of 1/kappa...
actually tauF is found by integrating the physical trajectory and matching
t(sigma=1)=tf -- easiest: integrate dY/dtau in TAU until t reaches tf, record
tauF = that tau length, then rescale to sigma. Continuity ~0 by construction.

GATE Z1-Sun: ||R|| <= 1e-8. The load-bearing prediction (S5): cond(J_sun)
should be ORDERS below cond(J_phys)~1e9 -- if it drops to ~1e5-1e6 the
linear-solve noise floor drops to ~1e-10 and the trust region + geodesic
reaches 1e-8. If cond does NOT drop, Sundman does not address this floor and
the honest fallback is extended-precision linear algebra.

## 5. Pitfalls

- CS safety: same rules as Z0 (no abs/norm/real/min/max/sign in the dY/dtau
  path; kappa via sqrt(sum(dd.^2))^pSund, not norm()).
- tauF sensitivity w: the inhomogeneous term G (=dY/dtau) is easy to forget --
  without it the tauF column of J is wrong and the solve won't move tauF.
- The terminal now has 8 rows (added t=tf); the seed must satisfy t(1)=tf to
  ~integrator tol (tauF chosen for it), else the seed residual is dominated by
  the time row.
- pSund knob: start 1.5 (PSR); if cond(J_sun) still ~1e8, raise to 2-3 (the
  costate 1/r^3 wants more than the direct 1/r^2).
