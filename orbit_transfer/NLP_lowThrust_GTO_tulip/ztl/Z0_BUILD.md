# Z0 build document — ramp-family EOM + variational STM + event automaton

Self-contained build spec for increment Z0 of the ZTL campaign
(`../ifs/PLAN_PRONG_Z.md` §4). Written so a session with NO prior context can
execute it. Read `ZTL_RESULTS.md` "P0 CONCLUSION" first (2 minutes) — it
states why Z0 exists and what it must prove.

## 0. Context in three sentences

The min-fuel GTO→tulip campaign (CR3BP, 15 kg, 25 mN, ~40 revs) has a working
DIRECT pipeline (PSR) but every INDIRECT solver crawls or floors. P0
measured the cause: differencing-based Jacobians (complex-step AND finite-
difference) are structurally unusable for multi-rev CR3BP shooting
(derivative scales span 11 orders; CS corrupted at O(1) by the adaptive
integrator, FD takes secants across a curved valley). Z0 builds the cure —
the exact variational-STM machinery (Zhang 2015's ingredient (a)) for the
Bertrand–Épénoy ramp throttle family — and must pass a pre-registered
acceptance test: punch through the banked 75 mN floor iterate
(`results/p0i_fd_finish.mat`, ||R|| = 1.56e-2) to 1e-8 in O(10–30) Newton
iterations.

## 1. Problem constants and state layout

All problem constants come from `ztl_endpoints.m` (already built):
`[rv0, rvf, P] = ztl_endpoints()` gives the GTO departure state, the tulip
max-ydot rendezvous target (both ND rotating-frame), and
`P.muStar/.lStar/.tStar/.m0kg/.g0/.c/.Tmax25` (Tmax25 = ND thrust
ACCELERATION at m=1 for 25 mN; thrust level X mN → `(X/25)*P.Tmax25`).

Augmented state, ND rotating frame, mass as fraction (m0 = 1):

    y = [r(3); v(3); m; lam_r(3); lam_v(3); lam_m]   in R^14

This layout matches `../../lowThrust_GTO_tulip/lt_pmp_eom*.m` exactly (those
files also contain the analytic gravity-gradient G and Coriolis Hc blocks to
copy). Time domain (NOT Sundman) — deliberate, see PLAN_PRONG_Z §2.1.

## 2. The throttle family (BE ramp — do NOT use tanh anywhere)

Cost family: J_eps = (Tmax/c) * int [ u - eps*u*(1-u) ] dt, eps in [0,1].
eps=1 is min-energy (integrand = u^2), eps=0 is min-fuel (bang-bang).

Switching function (min-fuel convention, scale-fixed by the "1"):

    S = 1 - ||lam_v||*c/m - lam_m

Optimal throttle = clamped linear ramp with EXACT saturation boundaries:

    u*(S) = clamp( 1/2 - S/(2*eps), 0, 1 )        (eps > 0)
    u*(S) = 1 if S < 0, 0 if S > 0                (eps = 0)

Three regimes with exact boundaries at S = -eps and S = +eps:

    ON     : S <= -eps   -> u = 1
    MEDIUM : |S| < eps   -> u = 1/2 - S/(2*eps)   (u depends on y through S)
    OFF    : S >= +eps   -> u = 0

Thrust direction is always the primer: alpha = -lam_v/||lam_v||.

Dynamics (copy the gr/hv/G/Hc construction from `lt_pmp_eom_minfuel.m`):

    rDot   = v
    vDot   = gr(r) + hv(v) + u*(Tmax/m)*alpha
    mDot   = -u*Tmax/c
    lrDot  = -G(r)' * lam_v
    lvDot  = -lam_r - Hc' * lam_v
    lmDot  = -||lam_v||*u*Tmax/m^2

(The costate equations are the same in all regimes by the envelope theorem;
only u differs. In MEDIUM, u is a function of y, which matters for the
JACOBIAN A(y) but not for these expressions.)

Analytic pieces needed for events and saltation:

    dS/dy   : dS/dm = c*||lam_v||/m^2 ; dS/dlam_v = -(c/m)*lam_v'/||lam_v||;
              dS/dlam_m = -1 ; zeros elsewhere.            [1x14 row]
    Sdot    = -(c/m)*(lam_v'*lvDot)/||lam_v|| + (c*||lam_v||/m^2)*mDot - lmDot

Hamiltonian (for diagnostics; autonomous so it is conserved):

    Ht = (Tmax/c)*(u - eps*u*(1-u)) + lam_r'*rDot + lam_v'*vDot + lam_m*mDot

## 3. Why the derivatives here are trustworthy (the design's core)

Two rules make this machinery immune to the P0 failure mode:

1. **Complex step is applied ONLY to the field, never through an
   integrator.** `ztl_A(y,P,regime)` perturbs y componentwise by i*h
   (h=1e-100... use 1e-50; anything <=1e-20 is fine) and evaluates
   `ztl_eom` AT FIXED REGIME. Within a fixed regime the field is analytic
   (no abs/min/max/sign on the perturbed quantities — see §4 coding rules),
   so CS is exact to machine precision. There is no adaptive stepper in the
   loop to corrupt it.
2. **The STM solves the variational ODE** Phi_dot = A(y(t))*Phi alongside
   the trajectory (14+196 = 210 real ODEs). The derivative is the exact
   derivative of the continuous system, delivered at integrator tolerance.
   No differencing of the flow ever happens.

At eps=0 the field JUMPS at switches; there the STM is corrected by the
saltation matrix (§6). For eps>0 the field is continuous (u is C^0), so NO
saltation is applied at regime boundaries — but the integration must still
STOP and restart at each boundary (A jumps; smearing a step across it
degrades the STM's accuracy order).

## 4. File 1: `ztl_eom.m`

    function [yDot, aux] = ztl_eom(y, P, regime)
    % y      [14x1] (may be complex under CS probing)
    % P      struct: muStar, c, Tmax, eps
    % regime 'on' | 'medium' | 'off'
    % yDot   [14x1]
    % aux    struct (computed only if nargout>1): S, Sdot, u, Ht

CODING RULES (complex-step safety inside a fixed regime):
- Norms as `sqrt(sum(x.^2))`, NEVER `norm()` (norm uses abs -> kills CS).
- NO `real()`, `abs()`, `min/max`, `sign()`, or branching on y anywhere in
  the dynamics path. The regime argument replaces all branching:
  u = 1, or (1 - S/eps)/2... precisely `0.5 - S/(2*P.eps)`, or 0.
- `aux` MAY use real() (it is diagnostic, not differentiated).
- Guard: `assert(P.eps > 0 || ~strcmp(regime,'medium'))`.

Copy gr/hv/G/Hc literally from `lt_pmp_eom_minfuel.m` (they are already
CS-safe). Keep the MATLAB header convention (purpose/inputs/outputs/refs).

## 5. File 2: `ztl_A.m`

    function A = ztl_A(y, P, regime)
    % 14x14 exact Jacobian of ztl_eom's field at fixed regime, by complex
    % step: A(:,k) = imag(ztl_eom(y + 1i*h*e_k, P, regime))/h, h = 1e-50.

~15 lines. This is machine-precision because of §4's rules.

## 6. File 3: `ztl_flow.m` — the 3-regime event automaton

    function out = ztl_flow(y0, tspan, P, wantSTM)
    % Integrate y (+ Phi if wantSTM) from tspan(1) to tspan(2),
    % segment-by-segment between regime-boundary events.
    % out: yf [14x1], PHI [14x14] (I if ~wantSTM), events (struct array:
    %      t, S, Sdot, from, to, grazed), nSegs, tGrid/yGrid (coarse dense
    %      output for diagnostics), flag (0 ok, 1 graze, 2 maxSegs).

Algorithm:
1. Classify the current regime from S(y): regime boundaries decided by S vs
   +/-eps; if within `1e-12` of a boundary, use the sign of Sdot to pick the
   side being ENTERED.
2. Integrate `ode89` (RelTol `P.odeRelTol` [default 1e-13], AbsTol
   `P.odeAbsTol` [1e-15]) with DIRECTIONAL terminal events for the current
   regime:
   - ON: event value S+eps, direction +1 (leaving upward).
   - OFF: event value S-eps, direction -1 (leaving downward).
   - MEDIUM: two values [S-eps (dir +1); S+eps (dir -1)].
   - eps=0: ON: S, dir +1; OFF: S, dir -1. (No MEDIUM.)
   Directional events are what make the restart-at-boundary safe: at the
   start of the new segment the event value is ~0 but moving AWAY, so it
   cannot re-fire. Do not "nudge" time forward — that breaks the STM.
3. At each event: record (t, S, Sdot); **graze guard**: if
   |Sdot| < P.grazeFloor (default 1e-4) set flag=1 and record grazed=true.
   If eps == 0, apply saltation to the STM:
       Psi = eye(14) + ((fPlus - fMinus) * dSdy) / Sdot
   with fMinus = field in the OLD regime at the event state, fPlus = field
   in the NEW regime, dSdy the [1x14] row of §2, Sdot evaluated with the
   OLD regime's field. PHI <- Psi * PHI.  (For eps>0: NO saltation.)
4. Continue until t reaches tspan(2) or `P.maxSegs` (default 400) segments.

STM propagation: when wantSTM, integrate z = [y; Phi(:)] (210 states),
dz = [f; reshape(A*Phi, [], 1)]. A is evaluated by `ztl_A` at the CURRENT
y(t) inside the ode function — yes, that is a CS call per RHS evaluation;
it is 14 field evaluations, cheap (the field is ~50 flops).

NOTE: the ode89 RHS receives z REAL (the flow itself is never complex).

## 7. File 4: `test_ztl_z0.m` — unit gates (all must print PASS)

- **G1 (legacy equivalence at eps=1).** The legacy energy law
  (`lt_pmp_eom_energy`: u = sat(Tmax*(||lam_v||/m + lam_m/c))) is the BE
  eps=1 law under the EXACT costate map
      lam_BE = (2*Tmax/c) * lam_legacy
  (derivation: match u* of J1 = int u^2/2 vs J2 = (Tmax/c) int u^2; both
  clamp identically). Test: take the banked P0i iterate
  (`results/p0i_fd_finish.mat`, field anchor.lam0, 75 mN), map it, and
  verify (a) u agrees pointwise at t=0 to 1e-12, (b) terminal state from
  ztl_flow(eps=1) matches ode113(lt_pmp_eom_energy) at tf to <= 1e-8
  (both tight-tol). This validates EOM + automaton + the map in one shot.
- **G2 (A exactness).** ztl_A vs central-difference of the field (h=1e-7
  relative, field-only differencing IS valid) at 5 random states per
  regime: rel error <= 1e-6 (FD-limited; CS is the truth).
- **G3 (variational STM).** On a SHORT arc (0.05 ND, no events),
  Phi from ztl_flow vs central-difference of the FLOW (perturb y0
  componentwise, re-integrate): rel error <= 1e-6. Short arc = differencing
  is valid there; this is the last time the flow is ever differenced.
- **G4 (event integrity).** Construct an eps=0.5 arc that crosses regimes
  (e.g. the mapped P0i iterate at eps=0.5); verify (a) event times
  reproducible to <= 1e-10 under RelTol 1e-12 -> 1e-13, (b) S at each
  event equals +/-eps to <= 1e-11, (c) u is continuous across each event
  (|u_minus - u_plus| <= 1e-10).
- **G5 (saltation, eps=0).** On a 1-switch eps=0 arc (shrink the span so
  exactly one S crossing occurs): Phi WITH saltation vs central-difference
  of the flow across the switch: rel error <= 1e-5. Verify Psi != I
  actually applied (compare with saltation disabled: error should be much
  worse).

## 8. File 5: `z0_accept_75mN.m` — THE acceptance test (pre-registered)

Load `results/p0i_fd_finish.mat` (anchor: lam0 [7x1 legacy], Tmax_mN=75,
tf=7.2344, rv0, rvf, P). Then:

1. Map the seed: `lam = (2*Tmax/c) * anchor.lam0` (Tmax = 3*P.Tmax25, ND).
2. Residual (smooth energy problem, eps=1):
       R(lam) = [ yf(1:6) - rvf(:) ; yf(14) ]        [7x1]
   with yf from `ztl_flow([rv0; 1; lam], [0 tf], P, true)`; Jacobian
       J = [ PHI(1:6, 8:14) ; PHI(14, 8:14) ]        [7x7]
3. Sanity gate A: ||R(lam_mapped)|| must be ~1.5e-2 (same floor as P0i, up
   to integrator differences; if it is orders different, the map or the
   EOM is wrong — STOP and fix).
4. Newton with Armijo backtracking (reuse the p0i loop shape, J now from
   the STM; equilibrate columns; no truncation needed at cond ~1e4, keep
   the tsvd fallback anyway).
5. **GATE Z0 (the campaign fork):** ||R|| <= 1e-8 within 30 iterations.
   - PASS -> the 75 mN LADDER ANCHOR is converged; save
     `results/z0_anchor_75mN.mat` (lam0_BE, resNorm, throttle stats, dV,
     prop) and update ZTL_RESULTS.md + the ztl-p0-findings memory. Z3's
     energy ladder starts from this point (march thrust down 75 -> 25 mN,
     <= 10% steps, tf FIXED at 7.2344).
   - FAIL (floors again with an EXACT J) -> the wall is not derivative
     quality; record the trace and pivot to the direct-side dual seed
     (PLAN_PRONG_Z §6 R4/stopping-rule logic applies).

## 9. Budget and pitfalls

- Costs: plain flow ~1-3 s at 75 mN (13 revs); with STM (210 ODEs +
  14 CS field evals per RHS) expect ~10-60 s per Newton iteration.
  Acceptance test total: ~5-20 min. Run headless per the matlab-headless
  skill (R2025b ONLY: /Applications/MATLAB_R2025b.app/bin/matlab -batch).
- Pitfall 1: `norm()`/`abs()`/`real()` anywhere in the EOM path silently
  breaks ztl_A. Grep the finished ztl_eom for all three before running G2.
- Pitfall 2: non-directional events re-fire at segment starts. All events
  MUST set `direction`.
- Pitfall 3: do not reuse the legacy tanh law (`lt_pmp_eom_minfuel`) or the
  legacy energy law inside ztl_eom — the ramp family is the ZTL standard;
  legacy appears ONLY in G1 as an oracle.
- Pitfall 4: the banked lam0 is in the LEGACY costate convention. Everything
  downstream of Z0 uses the BE convention. The factor is (2*Tmax/c) with
  BOTH in ND units, and it is thrust-dependent — recompute per rung.
- Pitfall 5: ode89 needs the RHS as f(t,z); keep P bound via a closure, not
  globals. Local functions at the END of script files (MATLAB parser rule).
