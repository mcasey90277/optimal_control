# Honest Evaluation — Min-Fuel ΔV vs t_f Front (GTO→Tulip, CR3BP)

**Date:** 2026-07-09 (written after the upper-band attack session)
**Audience:** any AI or reviewer picking up this campaign. Read
`LOW_THRUST_MINFUEL_CAMPAIGN.md` for the full history and
`direct/sundman_minfuel/OPTIMALITY_VERIFICATION_PLAN.md` for the verification plan.
This document is the candid self-assessment: what is trustworthy, what is not,
and what remains to build. Goals: (1) an **honest** ΔV vs t_f plot;
(2) **direct AND indirect** methods working at every t_f for min-fuel.

**Problem:** minimum-fuel low-thrust transfer, GTO → south-pole tulip orbit,
Earth–Moon CR3BP, 15 kg / 25 mN / Isp 2100 s, ~40-rev spiral. Min-time
baseline t_f = 6.290694 ND = 27.8845 d, ΔV 4.4665 km/s (never coasts).
Factors below are t_f / t_f^min. Code: `direct/sundman_minfuel/` (CasADi + IPOPT,
Sundman-regularized trapezoidal collocation, Bertrand–Épénoy energy→fuel
homotopy).

---

## Bottom line

The solver core is genuinely strong — the certification layer and the plot
are weaker than they look. Our "green" (PMP-certified) means **machine-tight
+ first-order-consistent with Pontryagin from the NLP's own KKT duals**. It
does NOT mean "on the true front." Proof from our own data: the grey 1.75×
point is a *feasible* trajectory (defect ~1e-14, endpoints enforced) at
**2.523 km/s** — *below* the certified-green 1.85× at **2.667 km/s**. A
certified extremal sitting above a feasible point at shorter t_f is (modulo
the phasing caveat below) a certified **local** minimum that is not the
global one. The current plot is honest about extremality, not optimality.

## What is trustworthy

1. **`casadi_minfuel_sundman.m` (core solver).** Sundman regularization
   dt/dτ = r₁^1.5 with τ_f held FIXED (a free τ_f makes a dense KKT column →
   MUMPS OOM — real, hard-won insight), cone-eliminated control
   (unit direction + throttle), exact AD Jacobian/Hessian, and two warm-start
   regimes: `warmTight=true` (mu_init 1e-4, bound_push 1e-9) for re-solving AT
   a bang-bang point, `warmTight=false` (mu_init 0.1, default push) for genuine
   continuation moves. Solutions converge to defect ~1e-14 with primer
   alignment ~0.06°. As *feasible extremal candidates*, the solutions are
   trustworthy.
2. **Energy-backbone continuation strategy.** The ε=1 (min-energy) problem is
   convex in the control and bifurcation-free: it threaded 1.15×→1.85× in 14
   steps (2026-07-09) and 1.15×→1.12× earlier, every step ~1e-14, zero drift.
   Bang-bang continuation across t_f fails (basin drift, MEX crashes); the
   smooth backbone + per-t_f independent sharpen is the validated recipe in
   BOTH directions. Backbone files banked: `energy_1.20.mat` … `energy_1.85.mat`.
3. **KKT-dual costate extraction + empirical-β switching-law check**
   (`verify_tf_front.m`, plan §D). Recovers the switching function
   S = 1 − β·W from the defect-constraint multipliers with one positive scale
   β pinned by least squares at switch intervals. Objectively flagged the
   scattered points. Transversality is gated RELATIVE
   (|λ_m(τ_f)|/max|λ_m| ≤ 1e-3) — an absolute gate is scale-dependent and
   wrongly failed larger-costate-scale solutions (fixed 2026-07-09).

## Honest problems

1. **Certification is self-consistency, not independence.** The PMP check
   uses the NLP's *own* duals. It catches convergence to non-extremals, but a
   discretization-level artifact fools the solution and its duals together.
   The independent adjoint check (Tier 1 of the verification plan) is OPEN.
2. **Local ≠ global; the plot doesn't yet distinguish.** The front is the
   lower envelope of ≥2 solution families:
   - many-switch family: optimal ~1.12–1.40×, ΔV min ≈ 2.961 at 1.40×, then
     FOLDS UP (1.45× rises to 3.098) — discovered 2026-07-09;
   - few-switch (~22-sw) family: holds 1.75–1.85× at ΔV 2.52–2.67.
   We only certify whichever family the seed finds. The branches cross near
   1.40–1.45×. Upper band 1.45–1.80× has NO certified envelope point yet.
3. **Discretization error unquantified.** Defect 1e-14 measures the DISCRETE
   trapezoid equations, not ODE accuracy. Trapezoid is O(h²) and smears
   switches. No mesh-refinement (N vs 2N) study exists for the min-fuel front;
   ΔV likely carries ~1e-3 km/s-level unquantified error, and switch counts
   (50/44/24…) partly reflect mesh chatter, not structure. Do not publish
   switch counts without a refinement check.
4. **Monotonicity caveat.** "ΔV non-increasing in t_f" strictly requires
   loitering at an endpoint. The tulip end is periodic in the rotating frame
   → argument holds modulo the tulip period; the GTO start precesses → small
   genuine phasing wiggles are possible between period multiples. The
   2.52→2.60→2.67 rise across 1.75/1.80/1.85× could be partly phasing, but
   0.14 km/s over 2.8 days is more likely local-minimum scatter. Unresolved.
5. **Code organization risk.** Eight overlapping drivers (`run_tf_sweep`,
   `run_tf_front`, `run_tf_2anchor`, `tf_step`, `energy_step`,
   `solve_tf_minfuel`, `direct_build_minfuel`, `build_energy_backbone`) with
   copy-pasted step logic and DIFFERENT homotopy schedules; historical scatter
   partly traces to this. Orchestration (watchdogs, retries) lives in
   throwaway shell scripts — a zsh `local`-line expansion bug
   (`local f=$1 out=...$f...` expands $f before assignment) clobbered all 13
   sharpen outputs on 2026-07-09 (~1.5 h compute; printed metrics survived).
   Smaller issues: `tfMin = Sm.tf/1.15` magic constant in three places;
   `%.2f` filenames collide on finer grids; `lamDef = lamAll(1:8N)` silently
   breaks if constraint order changes; `solve_tf_minfuel`'s schedule ends at
   ε=0.001 (never exactly 0 — tiny objective bias; reported ΔV is still real
   since computed from the flown mass).
6. **Sporadic uncatchable CasADi/IPOPT MEX fatal crashes** kill the whole
   MATLAB process (~1 in 10 solves). Mitigation that works: one solve per
   process + shell watchdog + one retry. Any orchestration must assume this.

## State of the data (2026-07-09, end of session)

| factor | ΔV km/s | sw | status |
|---|---|---|---|
| 1.00 (min-time) | 4.4665 | 0 | known anchor |
| 1.01–1.11 | — | — | HARD transition band, resists all methods (other terminal attacking) |
| 1.12 | 3.828 | 12 | certified green |
| 1.14 | 3.491 | 26 | certified green |
| 1.15 | 3.370 | 25 | certified green (the original certified solution) |
| 1.20 | 3.236 | 44 | certified green |
| 1.25 | 3.141 | 50 | certified green |
| 1.30 | 3.055 | 44 | solved 2026-07-09, .mat LOST to zsh bug — re-run cheap |
| 1.35 | 2.980 | 29 | same |
| 1.40 | 2.961 | 24 | same — minimum of the up-family |
| 1.45 | 3.098 | 24 | same — up-family folds; NOT on envelope |
| 1.50–1.70 | 3.2–5.6 | — | old scatter, non-extremal (grey) |
| 1.75 | 2.523 | 23 | FEASIBLE upper bound, fails switching-law check (grey) |
| 1.80 | 2.595 | 23 | same |
| 1.85 | 2.667 | 22 | certified green — but dominated by the 1.75× feasible point → local, not global |

Plot: `direct/sundman_minfuel/results/plots/front_full_verified.png`. Energy backbone
`energy_1.20…1.85.mat` all banked and machine-tight (the expensive part; any
single-t_f sharpen off it is ~10–30 min).

## Plan to the two goals

**Honest plot — three marker classes:** (a) feasible upper bound;
(b) direct-certified local extremal (KKT-dual PMP check); (c) direct+indirect
certified. Draw the envelope only through (b)/(c); a grey feasible point
below a green one is information, not noise.

**Direct v2 (consolidation):** one canonical `minfuel_at_tf(factor)`:
energy-backbone seed → tight re-clean → ε-sweep to 0. Branch enumeration per
t_f: sharpen from the energy seed AND from bang-bang continuation of both
neighbors; keep every certified extremal, record branch identity. N→2N mesh
refinement at 3–4 anchors to bound discretization error. Retire redundant
drivers; move watchdog/retry into a checked-in script (not scratchpad zsh).

**Indirect at each t_f (the missing independent certification):**
multiple-shooting PMP BVP seeded from the direct solution — the bridge that
makes indirect tractable where cold-start always failed here:
1. States seed from the direct trajectory; **costates seed from the de-scaled
   KKT duals** — the validated empirical-β fit IS the scale factor (this is
   the payoff of storing `lamDef` in every result file).
2. ~20–40 shooting segments kill the 1e6 single-shooting sensitivity
   (campaign lesson: single shooting fails on 40 revs for ANY objective).
3. Bang-bang control from the switching function via ODE event detection (no
   smoothing); fall back to ε-smoothed control-law homotopy if events chatter.
4. Fixed t_f, free final mass: transversality λ_m(t_f)=0; the cost multiplier
   normalization λ₀=1 fixes the costate scale (the "1" in S = 1 − (T/m)‖λ_v‖
   − (T/c)λ_m) — no normalization ambiguity.
5. Acceptance: indirect ΔV and switch times match direct to tolerance →
   independently certified point.

**Priority order on resume:** (1) re-run fixed up-sharpen → recover
1.30/1.35/1.40× with costates (~40 min, backbone banked); (2) few-switch
trace DOWN from 1.85× (and from the 1.75× feasible point) to settle
1.45–1.80×; (3) mesh-refinement spot checks; (4) build the indirect
multiple-shooting certifier; (5) merge with the lower-band result when the
other terminal's 1.01–1.11× attack lands.

**Session gotchas for AI operators:** run each solve in its own process with
a watchdog; never combine kill + launch in one shell command (SIGTERM races);
write scripts with the Write tool, not heredocs mixed with kills; in zsh
never `local f=$1 out=$DIR/x_$f.mat` on one line; the lower-band terminal
owns `energy_≤1.15` / `ms_*.mat` naming — the upper band uses `ms_up_*.mat`.

---

## Addendum (2026-07-09 late) — checker diagnosed; two-tier gate; front transformed

Written after the up-band campaign that followed the main text. Three
developments supersede parts of the above:

1. **The few-switch down-chain swept 1.45-1.80x** (neighbor continuation from
   the certified 1.85x): campaign-best **dV 2.434 km/s at 1.65x (46 d, 45.5%
   below min-time)**; switch count morphs 22->43 approaching the crossing at
   ~1.40-1.45x -- the "two families" are one connected structure over a fold.
2. **The mass certification failure above 1.25x was the CHECKER, not the
   solutions** (`diag_beta_checker.m`): the W^2-weighted LS beta estimator is
   fragile to ONE anomalous switch -- the FIRST -- whose implied scale decays
   smoothly with t_f (0.93 at 1.20x -> 0.12-0.27 in the dn family) while every
   other switch sits at a common scale to <0.8% MAD. Robust (median) beta:
   every envelope point 1.20-1.85x passes burn+coast at 99.8-100%.
3. **Two-tier gate implemented** in `verify_tf_front` (robust beta): tier 2
   FULL certified (first switch consistent within 10%) = {1.12, 1.14, 1.15,
   1.20, 1.85}; tier 1 INTERIOR certified (all gates pass, first switch
   flagged) = the entire rest of the envelope 1.25-1.80x incl. the 2.434
   minimum; tier 0 = 3 old scatter points. Honesty dividend: 1.25x DEMOTED
   from the old full-certified set (first-switch scale 0.84, previously
   smeared invisible by LS). The certified envelope now coincides with the
   best-known feasible envelope end-to-end (`results/plots/front_honest.png`).

**The open question is now precise:** the first-switch anomaly (departure-burn
cutoff), real and smoothly growing with t_f -- mesh under-resolution vs
genuine local non-extremality. That, plus the dominated 1.70-1.85x tail
(front(t_f>=1.65x) <= 2.434 by loiter-monotonicity), are the two sharpest
targets for the ms_band indirect arbiter.
