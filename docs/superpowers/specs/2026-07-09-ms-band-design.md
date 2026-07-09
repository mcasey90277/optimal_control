# Design: `ms_band/` — indirect multiple shooting for the 1.01–1.11x transition band

**Date:** 2026-07-09
**Status:** approved (design review with user)
**Goal:** close the last open segment of the min-fuel dV–time front for the
GTO -> south-pole-tulip low-thrust transfer (Earth–Moon CR3BP, 15 kg / 25 mN /
Isp 2100 s) by producing **PMP-certified extremals at ~0.01x t_f spacing
through the transition band t_f = 1.01–1.11x min-time**, joining monotonically
to the min-time anchor (1.00x: 4.4665 km/s, 0 switches) and the certified
1.12x point (3.8278 km/s, 12 switches).

## 1. Background and why multiple shooting

The band is where the min-fuel control reorganizes from many-switch bang-bang
(~12 switches at 1.12x) down to always-burn min-time (0 switches at 1.00x).
It is bifurcation-rich — switches are born/die across it — and it defeated
every direct-collocation method tried (energy-backbone continuation times out
or hangs; from-scratch direct builds stall in restoration). Full account:
`NLP_lowThrust_GTO_tulip/LOW_THRUST_MINFUEL_CAMPAIGN.md`, sections
"Down-sweep CRACKED" and "transition band". The campaign's own preferred next
tool (option 1) and the verification plan's Layer 5
(`sundman_minfuel/OPTIMALITY_VERIFICATION_PLAN.md`) both point at indirect
multiple shooting, seeded from the KKT-dual costates the direct solver already
saves.

Reference methods:

- **Zhang, Topputo, Bernelli-Zazzera, Zhao (JGCD 38(8), 2015)** —
  `min_fuel_papers/2015-8-J-LowThrustMinimumFuelOptimizationRTBP.pdf`.
  Indirect min-fuel in the same CR3BP with Bertrand–Epenoy smoothing
  continuation; converged ~150-switch solutions. We adopt the smoothing
  philosophy: **the switch structure is implicit** (found a posteriori as the
  smoothed switching function sharpens), which is exactly what a
  bifurcation-rich band needs. We do NOT adopt their hard-switch event
  detection + STM jump matrices as the workhorse: the jump matrix carries a
  1/S-dot factor that blows up at grazing switches — precisely the band's
  generic event. (Kept as an optional finisher; see section 8.)
- **Leomanni, Bianchini, Garulli, Quartullo (2021)** —
  `min_fuel_papers/Leomanni2021_low_thrust_orbit_transfers_made_easy.pdf`.
  Regularized "ideal elements" direct method. Role here: fallback architecture
  if MS also stalls (slowly-varying elements tame the problem a different
  way); their Lyapunov initial-guess generator is not needed (we have better
  seeds: converged neighbors and KKT duals).

## 2. Chosen approach (decision record)

**Approach A — smoothed-control multiple shooting with eps-continuation,
two-sided t_f march.** Selected over (B) hard bang-bang MS with event
detection + STM jumps (fragile at grazing switches; heavy build) and (C)
switch-time parameterization (requires per-t_f structure guesses; structure
change IS the band). B survives as an optional per-point finisher; C's
bifurcation diagram falls out of A's S(t) profiles for free.

Success bar (user): dense certified front through the whole band.
Folder style (user): campaign-style production code + living campaign doc
(not a guided tutorial).

## 3. Formulation

Fixed t_f = factor x 6.290694 ND (tStar = 382981.289129055 s). Augmented
state y = [r(3); v(3); m; lambda_r(3); lambda_v(3); lambda_m] in R^14,
dynamics = `lowThrust_GTO_tulip/lt_pmp_eom_minfuel.m` (unchanged): CR3BP +
low thrust with primer direction alpha = -lambda_v/||lambda_v|| and
Bertrand–Epenoy tanh-smoothed throttle
u = (1 - tanh(S/(2 eps)))/2, S = 1 - ||lambda_v|| c/m - lambda_m.

Partition [0, t_f] into M arcs (default M = 24, a parameter) with joints
**uniform in Sundman tau** (dt/dtau = kappa = r1^1.5) computed from the seed
trajectory — joints cluster near perigee where sensitivity accrues. Time
remains the ODE independent variable (integration through perigee was never
the wall; the ~1e6 full-spiral sensitivity was, and partitioning reduces it
to ~10^(6/M) per arc). Escalation if per-arc conditioning still bites:
formally Sundman-transform the augmented system (documented, not built).

Unknowns (7 + 14(M-1)):
  Z = [lambda_0 (7); y_2 (14); ...; y_M (14)]   (interior joint values;
  node 1 state rv0, m0 = 1 is known)

Residuals (14(M-1) + 7), square system:
  - continuity: y(t_k^-; from arc k-1) - y_k = 0 at each interior joint
  - terminal: rv(t_f) - rv_f = 0 (6), lambda_m(t_f) = 0 (1, transversality)

## 4. Solver

- **lsqnonlin Levenberg–Marquardt** with SpecifyObjectiveGradient=true (LM
  measured 4x better than dogleg on this landscape in the single-shooting
  campaign). Fallback: fsolve trust-region with JacobPattern.
- **Jacobian:** per-arc complex-step STMs — 14 complex-step integrations per
  arc (tanh form is complex-step safe by design; established pattern in this
  repo), parfor over arcs, assembled into the sparse block-bidiagonal
  structure. Validation gate: matches finite differences to ~1e-6.
  Fallback if too slow: propagate STM by variational equations with
  complex-step D_yF only.
- **Integrator:** ode113, RelTol = AbsTol = 1e-12.
- **eps-march at fixed t_f** (`eps_march.m`): 1 -> 0.3 -> 0.1 -> 0.03 -> 0.01
  -> 3e-3 -> 1e-3 -> 3e-4 -> 1e-4, warm-started; bisect a failed eps step.
  Target floor 1e-4; accept 1e-3 with the bang fraction reported honestly.

## 5. Continuation strategy (two-sided, guarded)

- **Up-pass:** anchor at 1.01x seeded from the min-time indirect solution
  (machine-zero residual, matches pumpkyn to 8 sig figs) — solve at eps = 1
  first, then sharpen; march t_f up in 0.01 steps, halving to 0.0025 on
  failure. Note: the min-fuel problem is degenerate at exactly 1.00x
  (measure-zero coast), hence anchoring at 1.01x, not 1.00x.
- **Down-pass:** anchor at 1.11x seeded from the 1.12x direct solution's KKT
  duals, beta-scaled to physical costates (global sign -1; scale determined by
  matching the switching-law beta, then verified by propagating one arc —
  the OPTIMALITY_VERIFICATION_PLAN.md section F.6 route); march down.
- **Guard discipline** (campaign lesson): a step that fails its gate is
  discarded — it never poisons the warm start or overwrites the best result.
  Watchdog timeouts on every headless run (hangs, not MEX crashes, are the
  expected failure mode for MATLAB-native code).
- Where the passes meet: agreement is a free cross-check. A genuine fold (no
  extremal family threading some sub-interval) shows as a documented branch
  gap/overlap — reported as a structural finding, not forced.

## 6. Verification gate (acceptance = certificate, not convergence)

Per point, ALL of:
1. ||R|| <= 1e-9. For the indirect method the converged residual IS the
   first-order PMP certificate (primer direction, adjoint ODEs, transversality
   hold by construction).
2. S(t) regularity: min |S| on each interior arc bounded away from 0 (no
   singular arcs); finite-slope zero crossings (no Fuller/chattering).
3. dV monotone non-increasing vs accepted neighbors in t_f.
4. Band edges cross-validate against certified direct points:
   1.12x = 3.8278 km/s / 12 sw; 1.15x = 3.3696 km/s / 25 sw (reproduced by
   the MS machine itself in milestone M1).

Saved per accepted point (`msband_<factor>.mat`): Z, lambda_0, arc joints,
eps floor, ||R||, S(t) profile, dV, prop_kg, switch count/times, bang
fraction, diagnostics.

## 7. Folder layout

`NLP_lowThrust_GTO_tulip/ms_band/` (sibling of `sundman_minfuel/`):

| file | role |
|---|---|
| `setup_paths.m` | adds `../sundman_minfuel` (params, endpoints, duals) + `../../lowThrust_GTO_tulip` (EOM); no file copies — single source of truth |
| `arc_boundaries_tau.m` | tau-uniform arc joint times from a seed trajectory |
| `ms_residual.m` | core: MS residual + assembled sparse complex-step Jacobian |
| `ms_solve.m` | guarded single solve at fixed (t_f, eps); returns success flag; never mutates the input seed |
| `eps_march.m` | eps schedule controller with bisection + guard |
| `seed_from_mintime.m` | 1.01x up-anchor seed from the min-time indirect solution |
| `seed_from_duals.m` | down-anchor seed from beta-scaled KKT duals (verify by one-arc propagation) |
| `tf_march.m` | outer two-sided t_f driver: stepping, halving, watchdog, per-point save |
| `verify_band_point.m` | the section-6 gate as a function; per-point verdict |
| `plot_band_front.m` | extend the full-front figure with band points |
| `test_ms_jacobian.m` | regression: Jacobian vs FD |
| `test_ms_reproduce.m` | regression: reproduce min-time 1.00x and certified 1.15x |
| `MS_BAND_CAMPAIGN.md` | living campaign record (same style as existing docs) |

Deferred (build only if needed): `polish_hard_switch.m` — approach-B
finisher (freeze switch structure from the converged smoothed solution, one
hard-switch Newton polish to machine precision) if eps = 1e-4 points fail the
gate; full Sundman-transformed adjoint system if per-arc conditioning bites.

## 8. Milestones (hard gate before proceeding)

- **M1 — validate the machine on certified ground.**
  (a) `test_ms_jacobian` passes (~1e-6 vs FD).
  (b) `test_ms_reproduce`: MS converges min-time at 1.00x (dV 4.4665 km/s;
  always-burn so eps is irrelevant there) AND the certified 1.15x from a
  beta-scaled dual seed (||R|| < 1e-9, dV 3.3696 to 4–5 digits, 25 switches).
  Also resolves the two M1 unknowns: contents of
  `sundman_minfuel/ms_1.12.mat` (else regenerate duals via
  `solve_tf_minfuel.m`) and the practical complex-step-through-ode113 cost.
  No band work until M1 passes.
- **M2 — anchors.** 1.01x (up-seed) and 1.11x (down-seed) converge through
  the eps-march and pass the gate.
- **M3 — the march.** Two-sided fill at 0.01 spacing (refine 0.005/0.0025
  where the switch cascade is busy). Product: the certified band front + the
  switch-birth map (which apogee's coast appears at which t_f — read off the
  S(t) profiles).
- **M4 — assembly.** Extend `front_full_verified.png` with the band points;
  update `MS_BAND_CAMPAIGN.md` + the main campaign doc; optional back-feed of
  MS solutions as direct-solver seeds for an IPOPT cross-clean.

## 9. Lessons-learned constraints (binding, from the campaign do-NOT-redo list)

1. No single shooting over the full spiral (~1e6 sensitivity) — MS only.
2. No resampling of oscillatory trajectories for warm starts; seeds sample a
   solution at its own consistent times (MS joints are evaluated ON the seed
   trajectory).
3. Guarded continuation: failed steps are discarded, never propagated.
4. Acceptance is the PMP/monotonicity gate, not "residual got small".
5. LM over dogleg for this landscape; complex-step (not FD) derivatives where
   possible.
6. Costate scale from beta-matching, not first-principles de-scaling (the
   reviews' de-scaling dispute is moot under beta-matching).
7. Environment: MATLAB R2025b headless per the matlab-headless skill;
   multi-line -batch strings fail (write script files); IPOPT-style `fort.6`
   noise does not apply here but watchdog + per-run logs do.

## 10. Risks

- **Costate normalization near 1.00x.** The min-fuel costates lose their
  natural normalization as coasts vanish. Mitigation: anchor at 1.01x, eps = 1
  first; if LM wanders, add a one-shot normalization (fix ||lambda_v(0)|| from
  the min-time solution's scale) during the eps = 1 solve only.
- **Fold points in the band.** If extremal families genuinely fold, dense
  single-family coverage is impossible; the two-sided march exposes this as a
  branch gap/overlap. That outcome is acceptable: it is the structural
  answer, documented with the fold location.
- **Complex-step cost.** 14 x M integrations per Jacobian; if slow, variational
  -equation STM fallback (section 4).
- **eps floor vs certificate.** At eps = 1e-4 the throttle is near- but not
  exactly bang-bang. If gate item 2 marginally fails, the approach-B polish is
  the escalation, scoped in section 7.
