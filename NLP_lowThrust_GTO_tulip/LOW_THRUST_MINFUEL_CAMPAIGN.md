# Low-thrust min-fuel campaign — learnings & next steps

**Record of the July 7–8 2026 effort to solve the minimum-fuel GTO → tulip
transfer on the full ~40-revolution spiral, and what it would take to finish
the many-switch case.** This is the synthesis; per-piece detail is in
`MIN_ENERGY_NOTES.md`, `OVERNIGHT_STATUS.md`, and
`../lowThrust_GTO_tulip/gto_tulip_mintime_theory.pdf` §6.

---

## SOLVED (Jul 8 2026) — sharp many-switch bang-bang, machine-tight

The sharp global min-fuel solution on the full ~40-rev spiral is **solved and
certified**. Winning recipe = **Sundman regularization + energy→fuel homotopy +
no-resample seed**, all three together (each alone was insufficient):

| quantity | value |
|---|---|
| IPOPT status | `Solve_Succeeded` |
| max collocation defect | **2.4×10⁻¹⁴** (machine zero) |
| terminal rendezvous error | **0** (hits the tulip point exactly) |
| switches | **25** |
| bang-bang fraction (nodes at a throttle bound) | **99.4%** (59.4% full / 40.0% coast / 0.6% transition) |
| propellant (15 kg, 25 mN, Isp 2100 s) | **2.2640 kg** |
| ΔV | **3.3696 km/s** |

Cross-checks: the ΔV matches the earlier *uncertified* ~50-switch estimate
(3.37 km/s) to 4 figures; the propellant (2.264 kg) is **23% below** the
certified 3-switch burn-then-coast *local* optimum (2.950 kg), confirming the
many-switch burn-near-perigee / coast-near-apogee structure is the genuine
global optimum. Solution: `sundman_minfuel_certified.mat`.

**The recipe, and why each piece is load-bearing:**

1. **Sundman regularization** (`casadi_minfuel_sundman.m`). Change the
   independent variable time→τ with `dt/dτ = κ(r) = r₁^p` (p=1.5), carrying
   time as an 8th state. Every ODE is multiplied by κ, so the 1/r³ perigee
   terms that blew up IPOPT's exact Hessian become r₁^(p−3) — bounded. A
   *uniform* τ-mesh also concentrates nodes near perigee in physical time
   automatically. This alone lets the collocation defects reach 1e-10/1e-14
   (they floored at ~5e-3 on the raw problem). **Gotcha that cost hours:** make
   τ_f a *fixed constant*, not a decision variable — a free τ_f multiplies every
   defect → one dense KKT column → catastrophic MUMPS fill-in → OOM/SIGKILL at
   N≳400. Fixed transfer time is instead enforced by the t-state terminal
   condition t(τ_f)=t_f; the trajectory adjusts so ∫κ dτ = t_f.

2. **No-resample seed** (`run_sundman_from_seed.m` / `run_sundman_homotopy.m`).
   Seed from a *collocation-feasible* solution (the certified energy-seeded
   min-fuel `.mat`, which hits the target exactly), mapped to τ using **its own
   nodes** — σ = τ/τ_f, no pchip interpolation onto a uniform-σ mesh.
   Downsampling a 40-rev oscillatory trajectory left an irreducible ~1e-2 defect
   that pinned IPOPT in the restoration phase (false "locally infeasible" exit).
   Using the seed's own nodes makes the only initial infeasibility the small
   time-trap vs Sundman-trap mismatch, which IPOPT closes to 1e-14 in *normal*
   mode. This was the single change that flipped the energy solve from
   perpetual-restoration to `Solve_Succeeded`.

3. **Energy→fuel homotopy** (`run_sundman_homotopy.m`, `run_sundman_tail.m`).
   Objective `J(ε) = ∫s dt − ε∫s(1−s) dt`: ε=1 is ∫s² dt (energy, strictly
   convex → smooth ramp, no bang-bang, no restoration), ε=0 is ∫s dt (fuel,
   bang-bang). Solve ε=1 first (converges machine-tight from the no-resample
   seed), then step ε→0, warm-starting each solve from the last. The throttle
   sharpens continuously: edge 28%→75%→94%→99.6% as ε:1→0, switches 14→25,
   propellant 2.359→2.264 kg — **every point certified to ~3e-14**. Near ε=0 a
   coarse step (0.04→0.025) breaks the basin (inf_du→1e10), so the tail uses
   fine steps with a **guard** (a step that doesn't converge tight is discarded,
   never poisons the warm start or overwrites the best result).

This is the literal realization of the campaign's own two prescriptions —
"start smooth, then sharpen" (homotopy) and "regularize the dynamics" (Sundman)
— applied *together*. Neither alone sufficed: Sundman without the homotopy
thrashes into restoration on the bang-bang objective; the homotopy without
Sundman can't drive the perigee defects below ~5e-3. The two walls had to fall
at once.

## PMP verification via the costates (Jul 8 2026)

The solver now returns the **discrete costates** so the direct solution can be
checked against Pontryagin's principle. `casadi_minfuel_sundman.m` reads the
KKT multipliers of the dynamics-defect constraints (`out.lamDef` `[8xN]` =
`[λ_r;λ_v;λ_m;λ_t]` per interval; full stacked `out.lamAll`) — these are the
discrete adjoint, up to a positive mesh-weight scaling and a global sign. It
also computes two **scale-invariant** PMP checks in-solver: `out.primerAlignDeg`
and `out.lamMassEnd`. `run_tf_sweep.m` saves `lamDef` per t_f.

**Verified on the certified solution** (eps=0 re-solve, defect 2e-14):
- **Primer-vector condition** — the NLP thrust direction matches the costate
  primer `-λ_v/‖λ_v‖` to **0.058°** on every burn arc. A direct-method solution
  that never forms costates satisfies the PMP direction law to a hundredth of a
  degree — an independent direct-vs-indirect cross-check.
- **Transversality** — `λ_m(τ_f) = -1.7×10⁻⁷ ≈ 0` (free final mass).

Both are scale-invariant (direction cancels any positive weight; 0 stays 0), so
they need no calibration. **Still to close** for a full first-order certificate:
the **switching-function sign law** (`S = 1-‖λ_v‖c/m-λ_m < 0` on burns, `>0` on
coasts, zero-crossings at the 25 switches) — scale-DEPENDENT, needs the
duals de-scaled to node costates (undo the trapezoid weight, `τ_f`, `κ`), or,
scale-free, decode the throttle-bound multipliers in `lamAll` directly; and an
optional independent adjoint-ODE consistency check. See
`sundman_minfuel/TIER1_PMP_CERTIFICATION_SCOPE.md` for the full plan and the
(failed) continuous-costate route that motivated using the duals.

## ΔV–time front (tf-sweep) — a dense local-minimum landscape (Jul 9 2026)

`sundman_minfuel/run_tf_sweep.m` maps the min-fuel ΔV vs transfer-time trade by
t_f-continuation (fixed endpoints; only t_f varies, imposed through the carried
time state; the smooth energy solution is continued across t_f and re-sharpened
per t_f). Every point converges machine-tight and PMP-consistent (primer
alignment 0.06–0.17°). **But the ΔV values SCATTER — they do not trace a clean
monotone front:**

| t_f (×min / days) | ΔV (km/s) | switches |
|---|---|---|
| 1.10 / 30.7 | 3.785 | 26 |
| 1.15 / 32.1 | 3.457 | 38 |
| 1.23 / 34.3 | 3.074 | 39 |
| 1.33 / 37.1 | 3.009 | 24 |
| 1.45 / 40.4 | 3.195 | 20 |
| 1.60 / 44.6 | 4.185 | 28 |
| 1.80 / 50.2 | 2.708 | 24 |

(For reference the fine-schedule certified point at 1.15× is **3.370 / 25 sw** —
*better* than this sweep's 1.15× 3.457/38-sw, i.e. even the anchor fell into a
worse basin under the coarser sweep schedule.)

**Finding — the scatter is the physics, not a solver bug.** Each point is a
valid, distinct *local* minimum. The many-switch min-fuel problem has a **dense
set of local optima**: with ~40 apogees, there are combinatorially many
near-optimal choices of which apogees to coast through, each a separate
bang-bang basin. Single-thread t_f-continuation lands in whatever basin the
warm start drifts into at each t_f, so ΔV scatters instead of tracing one
family. This is the same local-minimum multiplicity that made the base problem
hard, now made visible. The **lower envelope** still shows the genuine trade —
min-time 4.47 (27.9 d) → ~3.0 (35 d) → **2.71 (50 d), ~39% below min-time** —
but the exact front curve cannot be pinned by single-thread continuation.

**To get an accurate front:** multi-start per t_f (K independent trials —
perturbed sharpening starts, varied schedules, and seeding from the *best-known*
certified solution rescaled — take the ΔV-minimum), then the best-of-K lower
envelope is the front. Finer t_f resolution helps only *after* each point is the
local best (more t_f alone just adds scatter). Per-t_f solutions (state,
control, costates) are saved in `tf_sweep_results.mat`; scatter plot
`attic/tf_dv_front.png`.

### Follow-up: certified-basin continuation (`run_tf_front.m`, Jul 9)

Continuing the CERTIFIED 25-switch solution in small (5%) t_f steps —
`run_tf_front.m` — produces a MUCH cleaner curve than the energy-continuation
scatter, but only over a limited band, and it exposes the basin structure
directly (data + plot `attic/tf_front.png`):
- **Good band, t_f = 1.15–1.45× (32–40 d):** a clean, smooth, decreasing front,
  ΔV 3.370 → **2.961 km/s** (min near 1.35×/37.6 d), then a slight uptick to
  3.063 at 1.45×. This is the trustworthy segment.
- **Drift, 1.50–1.70× (42–47 d):** ΔV climbs to an unphysical 5.58 km/s — the
  certified basin stops being near-optimal and the small-step continuation
  follows a suboptimal branch (edge% drops, throttle smears).
- **Basin jump, 1.75–1.85× (49–52 d):** the continuation falls into a
  *different, much better* family (22–23 switches, ΔV **2.52**→2.60→2.67) — the
  lowest ΔV found, ~44% below min-time.

**Reading:** the true front is monotone non-increasing in t_f (more time can
always coast more), so the drift hump and the 2.52-then-rising tail are both
continuation artifacts, not the front. The reliable takeaways: (i) the good band
gives a real 32–40 d segment; (ii) the 1.75×+ basin proves ΔV ≈ 2.5 km/s is
reachable at ~49 d. A single basin does NOT thread the whole t_f range — the
optimal family CHANGES with t_f (more coast ⇒ fewer, differently-placed
switches). The run CRASHED (CasADi/IPOPT MEX fatal error) on the down-pass, so
the 1.05–1.10× near-min-time points are missing; 15 up-pass points are saved in
`tf_front_results.mat`. Getting the true monotone front needs a scheme that
re-seeds the basin as t_f grows (or global multi-start per t_f).

### Per-t_f PMP verification of the front (`verify_tf_front.m`, Jul 9)

Each front point is checked against Pontryagin's first-order conditions using
its own KKT-dual costates (`out.lamDef`), via the validated empirical-β
switching-law route (`OPTIMALITY_VERIFICATION_PLAN.md` §D): recover
`S = 1 − β·W`, β pinned by LS at the switch intervals (absorbs the covector
scaling), then require burn-sign `S<0`≥99%, coast-sign `S>0`≥99%, β-spread ≤5%,
primer ≤0.2°, |λ_m(τ_f)|≤1e-3. This is an **objective per-t_f optimality test**
— it tells a genuine local extremal from a merely-converged-looking point.

Run on the 15-point front (`tf_front_verified.png`, green=certified):
**only 1.15–1.25× (32–35 d) certify** as first-order PMP extremals
(β-spread ≤2.3%, burn-sign ≥99.6%). The dV curve *looked* clean out to 1.45×,
but the switching law exposes that the light 4-step re-sharpen converged the
**primal** (dV, defect ~1e-14) without converging the **multipliers** — β-spread
climbs 0.6%→17% and burn-sign agreement falls 100%→56% past 1.25×. The primer
stays ~0.06° everywhere (scale-invariant), so **the switching law is the
discriminating test.** Lesson for the front: acceptance must be PMP-pass, not
just defect<1e-6, and re-sharpening must be fine enough to converge the costate
structure. `verify_tf_front.m` is the reusable checker (any results `.mat` →
per-point verdict table + colored plot); `run_tf_2anchor.m` continues the
certified AND the 1.75× low basin with a finer schedule to try to extend the
certified band and fill the 1.30–1.70× gap.

### Down-sweep CRACKED (Jul 9 2026) — energy backbone + tight re-clean

Extending the front DOWN in t_f (toward min-time) defeated every method until
this. The recipe below produced the first successful down-step in the campaign:
**1.13× = 3.5955 km/s, 24 switches, defect 4.5e-14, primer 0.097°** — converged,
PMP-consistent, monotone-correct (ΔV rises below 1.15×: 3.596 > 3.370, as
physics requires; ordering 3.596/1.13 > 3.370/1.15 > 3.236/1.20 > 3.141/1.25).

**Why the down-direction is uniquely hard, and what failed:**
- *Bang-bang continuation* (rescale a min-fuel solution to smaller t_f,
  re-sharpen): reducing t_f forces coast arcs to become burns (active-set
  change) — this **MEX-crashes** IPOPT (uncatchable, kills the process).
- *Fresh homotopy from a rescaled min-fuel seed*: rescaling the 1.15× solution
  to a shorter t_f makes it dynamically inconsistent, so the ENERGY solve either
  **plateaus at ~1e-3** or **blows up to inf_du ~1e10**.

**The working recipe** (two independent robust primitives, kept separate):
1. **Energy→energy continuation is crash-free DOWNWARD.** The ε=1 problem is
   convex in the control, so continuing the SMOOTH energy solution to a
   neighbouring (smaller) t_f converges to machine zero and never crashes —
   where every bang-bang method died. Build a *backbone* of energy solutions by
   stepping t_f down in small (~2%) steps from the cleanly-converged 1.15×
   energy, each warm-started (loose) from the previous. (`build_energy_backbone.m`)
2. **A tight re-clean is required before sharpening.** The loose continuation
   gives a defect-tight but *multiplier-inconsistent* energy solution; sharpening
   it directly blows up (inf_du ~1e10 on the first ε=0.6 step). Re-solving ε=1 AT
   the same t_f with the TIGHT warm start (no move ⇒ no wedge) cleans the duals,
   after which the fine energy→fuel sharpen converges every step to ~1e-14.
   (`solve_tf_minfuel.m`)

**Full down-recipe:** good energy → continue *energy* down in small steps
(smooth, crash-free) → tight multiplier re-clean at each t_f → fine sharpen to
bang-bang → PMP-verify + monotonicity. Because the backbone energy solutions are
built first (sequential, cheap), the per-t_f re-clean+sharpen steps are
INDEPENDENT and run in PARALLEL. Uncatchable MEX crashes are isolated by running
each t_f in its own process. This is the method for the whole down-band toward
Darin's min-time t_f (27.88 d, 4.4665 km/s, 0 switches — the known endpoint).

### Down-sweep status: a hard TRANSITION BAND at ~1.01–1.11× (Jul 9 2026)

Pushing the method above toward min-time hit a wall that BOTH continuation and
from-scratch solves share. Current map of the front (machine-tight solutions,
ΔV monotone-decreasing in t_f as physics requires):

| t_f× | days | ΔV (km/s) | switches | PMP |
|---|---|---|---|---|
| 1.12 | 31.2 | 3.8278 | 12 | **certified extremal** ✅ |
| 1.13 | 31.6 | 3.5955 | 24 | certified down-step (earlier build) |
| 1.14 | 32.5 | 3.4905 | 26 | **certified extremal** ✅ |
| 1.15 | 32.1 | 3.3696 | 25 | **certified extremal** ✅ |
| 1.20 | 33.5 | 3.2355 | 44 | **certified extremal** ✅ |
| 1.25 | 34.9 | 3.1409 | 50 | **certified extremal** ✅ |
| **1.01–1.11** | 27.9–30.8 | — | — | **transition band; resists all methods** ❌ |
| 1.00 (min-time) | 27.9 | 4.4665 | 0 | known endpoint ✅ |

Monotone ordering holds: 3.83(1.12) > 3.49(1.14) > 3.37(1.15) > 3.24(1.20) >
3.14(1.25). The **switch count drops toward the band** (25→12 from 1.15× to
1.12×) — direct evidence of the many-switch → min-time reorganization.
Full-front figure: `results/plots/front_full_verified.png` (green=PMP-certified, grey=not,
black=min-time), assembled by the `combine_front.m` pattern + `verify_tf_front`.

**Transversality nuance (RESOLVED):** the new down-points 1.12×/1.14× have
NEAR-PERFECT switching-law fits (100% burn-sign, 100% coast-sign, β-spread
0.5%) — as strong as the certified points. They initially colored grey because
`verify_tf_front` applied an ABSOLUTE transversality gate (|λ_m(τ_f)|≤1e-3), but
the costates are known only up to a positive scale, so an absolute gate is
scale-dependent. These fewer-switch down-band solutions carry a LARGER overall
costate scale (max|λ_m| ≈ 31 vs the certified points' smaller magnitudes), so
their raw λ_m(τ_f) ≈ −4×10⁻³ is a genuine zero relative to their own scale
(|λ_m(τ_f)|/max|λ_m| ≈ 1.3×10⁻⁴ at 1.12×, 2.0×10⁻⁴ at 1.14×). Fix: `verify_tf_front`
now gates on **relative transversality** |λ_m(τ_f)|/max|λ_m| ≤ 1e-3, which is
scale-invariant. Both points certify green; 1.85× (also fewer-switch) certifies
as well. Certified band is now **1.12×–1.25× + 1.85×**.

**What was tried on the band, and how each failed:**
- *Continuation (energy backbone, watchdog, 0.01 steps):* cleanly reached
  1.14/1.13/1.12×, then **1.11, 1.10, 1.09, 1.08, 1.07× all time out** (the
  watchdog auto-skips; each retry is a bigger jump from 1.12×, so it can't
  re-establish). 1.11× specifically *hangs* the solver (line-search stall).
- *Direct build from scratch (`direct_build_minfuel.m`):* a fresh burn+coast
  warm start built AT the target t_f → min-energy → homotopy. Tested at 1.07×,
  1.03×, 1.01×: the min-energy solve **stalls in restoration (inf_pr frozen
  ~0.3–0.6, inf_du → 1e8) or oscillates without converging.** The fresh warm
  start is not good enough — the same restoration wall the 1.15× seed took a
  hard multi-stage effort (+ no-resample) to beat. Counter-intuitively 1.01×
  (closest to min-time) is WORSE than 1.03×, likely near-degenerate (~1% coast
  ⇒ barely any structure to optimize).

**Why the band is hard (the physics):** 1.01–1.11× is the TRANSITION ZONE where
the min-fuel control reorganizes from the many-switch regime (~24 switches at
1.12×) down to the always-burn min-time limit (0 switches at 1.0×). Switches
disappear through this band, so it is **bifurcation-rich** — exactly the setting
where a direct-collocation basin fragments and continuation jumps or hangs.

**How to proceed (options, in preference order):**
1. **Indirect / multiple-shooting PMP, marching UP from min-time.** This is the
   tool matched to a bifurcation region: parameterize by the switch TIMES
   (which are few and shrinking toward min-time) and solve the TPBVP, seeded by
   the KKT-dual costates we already save (Layer 5 of
   `OPTIMALITY_VERIFICATION_PLAN.md`). Sundman-regularized MULTIPLE shooting
   (~20 arcs) keeps per-arc conditioning tame. Heavier, but structurally right.
2. **Fine-grained assault:** much smaller steps (0.005) from 1.12× with the
   watchdog + hang-timeout, accepting a low hit rate, to see whether ANY interior
   points are reachable and where exactly the switch count drops.
3. **Bank the certified band** (1.12–1.25×) + the min-time anchor as the result,
   documenting 1.01–1.11× as a genuine hard transition region — itself a
   publishable structural finding about the problem.

Infrastructure for the down-sweep: `build_energy_backbone.m` / `energy_step.m`
(chained energy continuation, process-isolated), `solve_tf_minfuel.m` (re-clean
+ sharpen), `direct_build_minfuel.m` (from-scratch burn+coast build), `tf_step.m`
(isolated single step). Watchdog + process isolation handle the uncatchable MEX
crashes and solver hangs. (Post-cleanup the canonical equivalents are
`minfuel_at_tf.m` + `orchestrate/*.sh`; the older drivers are atticked.)

## Up-band 1.30–1.45×: the fold region loses extremal support (Jul 9 pm 2026)

Branch enumeration above the certified 1.25× point, three independent seeds per
t_f — (a) the old up-pass bang-bang continuation, (b) the energy-backbone
sharpen (via `minfuel_at_tf` 'energy'), (c) a fresh neighbor-seed chain
continued from the certified 1.25× solution itself ('neighbor', light
schedule). Two decisive facts:

1. **The ΔV envelope is robust**: all three basins agree to ~0.01 km/s —
   1.30× ≈ 3.032, 1.35× ≈ 2.96–2.98, 1.40× ≈ 2.961 (the up-family minimum),
   1.45× ≈ 3.056. Solid feasible upper bounds, machine-tight, primer ~0.06°.
2. **The dual-inconsistency is a property of the REGION, not the seed**:
   burn-sign agreement matches across all three seeds to <0.2% at every
   factor and degrades smoothly toward the fold — 98.0% (1.30×), 92.6%
   (1.35×), 79.4% (1.40×), 71.6% (1.45×) — while β-spread grows 4.9→16%.
   Even the certified family's own gentle continuation loses switching-law
   support immediately past 1.25×. Not basin luck: the many-switch pattern
   itself stops being extremal-supported where the up-family folds
   (its ΔV turns back up at 1.45×).

Consequence: certified greens stay {1.12–1.25×, 1.85×}; 1.30–1.45× are
well-reproduced feasible envelope points, deliberately grey. The true optimum
in 1.30–1.80× most likely belongs to the FEW-SWITCH family that holds
1.75–1.85× (2.52–2.67, ~22 sw) — next probe: continue that family DOWN from
1.85×/1.75× into the gap. Independent arbiter remains the ms_band indirect
certifier. Diagnostic still open: where the burn-sign violations live
(near-switch vs mid-arc).

## The problem
GTO (350 × 35786 km, ω = −25°) → a south-pole tulip point in the Earth–Moon
CR3BP. 15 kg, 25 mN, Isp 2100 s (muStar = 0.012150585609624). Minimize
propellant at fixed transfer time. The full transfer is ~40 revolutions with
~40 perigee passes; the min-fuel optimum coasts near apogees and burns near
perigees, so a fully-resolved solution has **dozens of bang-bang switches**.

## Chronology — what we tried and what each taught

1. **Min-time (baseline, solved).** Indirect (complex-step single shooting)
   matches pumpkyn to 8 sig figs; direct NLP converges. Always-burn, no
   switches, so it is the easy case and the warm-start source for everything
   else.

2. **Min-fuel arrival leg (direct solved, indirect open).** Posed on the
   τ ≥ 4.0 tail with a phase-shifted target so a burn+coast start is feasible.
   Direct NLP converges machine-tight (single switch). Indirect TPBVP does
   **not** converge from the best seed — four seed strategies stalled at
   ‖R‖ = 1.55/0.83/0.33/0.14.

3. **Min-energy (the homotopy root).** Minimize ∫½u²dt → *continuous*
   saturated-ramp control, no bang-bang discontinuity. **Indirect converges
   to machine zero (‖R‖ = 6×10⁻¹²) on the leg** — where min-fuel's indirect
   never did. Full-spiral direct produces a coherent solution (defect 3×10⁻⁴).
   Full-spiral *indirect* does **not** converge (single shooting through ~40
   perigees = ~10⁶ sensitivity; ‖R‖ stalls at 0.236 even from a covector seed).

4. **Energy seeds fuel (the clean win).** Feeding the min-energy full-spiral
   solution into the min-fuel NLP as a warm start makes min-fuel **converge
   on the full spiral** — machine-tight (defect 2×10⁻¹⁵), 3 switches at
   tf = 1.15× min. This is the **certified full-spiral min-fuel result**
   (`minfuel_from_energy_seed.mat`). NB it is almost certainly a benign
   *local* optimum (burn-then-coast); the global optimum is many-switch.

5. **tf-continuation (grows switches).** Fixed tulip target, step tf up,
   warm-start each step from the last. Switch count climbs 0→2→6→…→~50.
   BUT the acceptance guard checked only the dynamics defects, not the
   throttle cone — so it accepted **cone-loose** solutions (cone off by
   ~2×10⁻²). Guard since fixed to check all constraints.

6. **Cone-elimination (fixes a real conditioning wedge).** The (w, s) control
   with cone w'w = s² is degenerate at coasts: s→0 forces w→0, where the
   thrust *direction* w/‖w‖ is undefined, which wedges the solver. Replacing
   it with a **unit direction α (‖α‖=1) + separate throttle s** decouples
   them and removes the wedge (seed feasibility 2.3×10⁻² → 1.7×10⁻³ by
   construction). Necessary, but **did not sharpen** the throttle by itself.

7. **Hermite-Simpson (4th-order, built & gradchecked).** Trapezoidal
   collocation penalizes sharp switches (a 0→1 jump between nodes creates a
   large defect, so the optimizer *smears* the throttle — the smeared
   throttle is the genuine trapezoidal optimum). HS should fix that. Built,
   Jacobian verified to 10⁻¹¹ — but with fmincon it made things *worse*
   (unit constraint drifted, feasibility 0.012), because…

8. **…the wall is the NLP solver, not the collocation order.** Every min-fuel
   collocation solve (trapezoidal, cone-eliminated, HS) plateaus at ~10⁻³
   with a smeared throttle under **fmincon interior-point + lbfgs** — the
   limited-memory Hessian can't resolve the coupling of a many-switch
   bang-bang problem. Min-*energy* and min-*time* (smooth control) converge
   fine; only min-*fuel* stalls. Textbook reason to move to an exact-Hessian
   solver.

9. **CasADi + IPOPT (installed, real bug fixed, but still stalls).**
   CasADi 3.7.0 at `~/casadi-3.7.0` (prebuilt arm64, loads despite the broken
   MEX *compiler*), bundles IPOPT and does automatic differentiation → exact
   sparse Jacobian AND Hessian. Found & fixed a genuine model bug (chained
   bound `lb<=X<=ub` is malformed in MATLAB → the box constraints were
   wrong). Even so IPOPT stalls on the many-rev min-fuel: defect floors at
   ~5×10⁻³ and the throttle thrashes.

## Why min-fuel is the hard one — the fundamental lesson

Min-fuel is dramatically harder than min-time and min-energy, and the root
cause is the **convexity of the Hamiltonian in the control**. How the control
enters the running cost decides everything:

- **Min-energy** — cost ∝ ∫u² dt is *strictly convex* in u. The optimal
  control is the unique interior minimizer: a *continuous* saturated ramp.
  Smooth control → smooth costate ODEs and a smooth shooting residual → large
  convergence basin. The genuinely easy problem.
- **Min-fuel** — cost ∝ ∫u dt is *linear* in u. Minimizing a linear function
  over a bounded set puts the optimum on a *bound*; the switching function's
  sign flips discontinuously → genuine **bang-bang**, non-smooth control.
- **Min-time** — also bang-bang in principle, but for this transfer maximum
  thrust is always optimal (coasting only wastes time), so it *degenerates*
  to always-on with **no switches**. Its ease is problem-specific (a
  switch-free optimum), not inherent smoothness.

So the operative distinction is not "min-fuel vs the rest" by objective name,
but that **min-fuel is the only one with genuine control discontinuities.**
That non-smoothness is what (a) shrinks the indirect shooting basin to
nothing and (b) makes low-order collocation *smear* the throttle — a jump
between two nodes creates a large defect, so the optimizer rounds the switch
off. (This is why the smeared throttle is the true trapezoidal optimum, not a
convergence failure.)

**But the campaign showed the difficulty is really TWO independent walls that
compound:**
1. *Objective-side* — the bang-bang non-smoothness above.
2. *Dynamics-side* — the 40-revolution CR3BP spiral: ~10⁶ shooting
   sensitivity from the perigee passes, and the 1/r³ near-perigee singularity
   that destabilizes even an exact-Hessian solver. This wall is independent
   of the objective; min-time and min-energy face it too.

Min-energy removes wall #1 — which is exactly why seeding min-fuel from
min-energy suddenly worked — but wall #2 remained and is what finally stopped
us (the regularization need). A many-revolution low-thrust transfer is
vicious precisely because it hands you *both* walls at once; min-time and
min-energy were "easier" only because they were fighting one wall, not two.

**The silver lining — and the practical prescription — follow directly:**
because min-energy is the *convex relaxation* of min-fuel, solve the smooth
(convex) problem first and deform toward the non-smooth one. That is the
whole basis of the energy→fuel homotopy (Bertrand–Épénoy), and it was by far
the most effective single move in this campaign. **Start smooth, then
sharpen.**

## The core learnings (reusable)

- **Smooth control is the homotopy root.** Min-energy (continuous control) is
  what makes min-fuel tractable. Solve the smooth problem, use it to seed the
  bang-bang problem. This unblocked the full-spiral min-fuel entirely.

- **Two independent difficulties, kept separate.**
  (i) *Objective-side* — bang-bang nonsmoothness; min-energy removes it.
  (ii) *Dynamics-side* — the ~40-perigee sensitivity / the CR3BP `1/r³`
  near-perigee singularity; independent of objective, and the deeper wall.

- **Single shooting is the wrong tool for 40 revs** (any objective): ~10⁶
  sensitivity, basin too small. Needs multiple shooting.

- **The exact Hessian can HURT near a singularity.** IPOPT's exact Hessian
  carries `1/r⁵` terms that blow up at each of the ~40 perigees (r ≈ 0.017),
  destabilizing its Newton steps — ironically fmincon's *crude* lbfgs Hessian
  was more stable because it's damped. The exact Hessian is only an asset
  once the problem is **regularized and scaled**.

- **The many-switch min-fuel optimum is genuinely global-hard.** The clean
  3-switch answer is a local optimum; the global optimum coasts at every
  apogee (~80 switches) and is what both solvers keep lunging at and missing.

- **Gotchas banked:** trapezoidal smears sharp switches (mesh, not solver);
  the throttle cone degenerates at coasts (use unit-direction); MATLAB does
  not chain `a<=x<=b` (write both sides); C++ MEX compilation is broken on
  this Mac but prebuilt binaries load fine; a continuation guard must check
  ALL constraints, not just the dynamics defects.

## Current state
- **SOLVED & certified (Jul 8):** sharp many-switch bang-bang full-spiral
  min-fuel — 25 switches, 99.4% bang-bang, propellant 2.2640 kg, ΔV 3.3696
  km/s, defect 2.4×10⁻¹⁴, terminal error 0. Via Sundman + energy→fuel homotopy
  + no-resample seed (see the SOLVED section up top).
  `sundman_minfuel_certified.mat`.
- **Certified (earlier):** min-energy → 3-switch machine-tight full-spiral
  min-fuel *local* optimum (defect 2×10⁻¹⁵, 2.950 kg — now superseded as a
  local minimum by the 2.264 kg global result). Min-energy leg indirect
  (6×10⁻¹²). Min-energy full-spiral direct (3×10⁻⁴).

## Next steps — polish & payoff (the hard problem is done)
The core numerical objective is met. Remaining work is packaging and payoff:

1. **Movie of the certified bang-bang solution** — the throttle/perigee-burn
   structure over the 40 revs (reuse the `movie/` animators; the solution is
   `sundman_minfuel_certified.mat`, states in Sundman τ with time as X(8,:)).
2. **tf sweep for the ΔV–time trade** — re-run the homotopy at several
   tf/tf_min values (larger tf ⇒ more coast ⇒ lower ΔV, more switches) to map
   the min-fuel Pareto front. The pipeline (`run_sundman_homotopy` +
   `run_sundman_tail`) now does one point robustly; loop it over tf.
3. **Independent verification (Tier-1 PMP certification)** — check the 25-switch
   structure against the PMP switching function S = 1 − ‖λ_v‖c/m − λ_m sign
   changes as a cross-solver certificate. **Scoped + attempted Jul 8; full record
   in `sundman_minfuel/TIER1_PMP_CERTIFICATION_SCOPE.md`.** FINDING: recovering
   costates from the PRIMAL (trajectory + primer directions) fails — the
   homogeneous costate map amplifies ~5×10¹¹ over the 40 revs and rendezvous
   leaves no BC on λ_r,λ_v, so the recovery is ill-posed (three methods tried,
   all fail; gravity-gradient model FD-validated, so it is conditioning not a
   bug). Corrected route: read the KKT DUALS from IPOPT (stable factorization),
   which needs `casadi_minfuel_sundman.m` to return `opti.dual` + a warm-started
   eps=0 re-solve to regenerate the `.mat`. `certify_minfuel_pmp.m` is in place as
   the scaffold (report/figure/switch-alignment work; recovery front-end pending
   the dual swap). Go/no-go on the pivot open.
4. **Optional indirect confirmation** — multiple shooting (partition the 40
   revs) could re-derive the same optimum from PMP for a belt-and-suspenders
   check, but is no longer *needed* — the direct result is machine-tight.

## File inventory (NLP_lowThrust_GTO_tulip/)
- Min-time direct: `lt_dynamics.m`, `nlp_constraints.m`, `solve_tfmin_nlp.m`,
  `NLP_lowThrust_GTO_Tulip.m`.
- Min-fuel (cone form): `lt_dynamics_throttle.m`, `nlp_constraints_minfuel.m`,
  `solve_minfuel_nlp.m`, `NLP_lowThrust_GTO_Tulip_minfuel.m`,
  `costate_seed_from_nlp.m`.
- Min-energy: `solve_energy_nlp.m`, `NLP_lowThrust_GTO_Tulip_energy.m`,
  `costate_seed_from_nlp_energy.m` (+ indirect in `../lowThrust_GTO_tulip/`).
- Cone-eliminated (unit-direction): `lt_dynamics_dirthrottle.m`,
  `nlp_constraints_minfuel_ue.m`, `solve_minfuel_nlp_ue.m`.
- Hermite-Simpson: `nlp_constraints_minfuel_hs.m`, `solve_minfuel_nlp_hs.m`.
- CasADi+IPOPT: `casadi_minfuel_trap.m` (needs `~/casadi-3.7.0`).
- **Sundman + homotopy (the solver that worked) — MODULARIZED into
  `sundman_minfuel/`** (self-contained library, code-reviewed Jul 8):
  `cr3bp_lt_params.m` (constants), `gto_tulip_endpoints.m` (BCs),
  `sundman_seed_map.m` (no-resample time→τ map), `casadi_minfuel_sundman.m`
  (core solver: Sundman-regularized, ε-parametrized energy→fuel objective,
  fixed τ_f), `sundman_homotopy.m` (guarded ε:1→0 sweep, folds the old
  from_seed/homotopy/tail trio), `run_certified_minfuel.m` (end-to-end entry),
  `README.md`. Certified result + seed `.mat` live in the folder.
  The original flat drivers (`run_sundman_from_seed.m`, `run_sundman_homotopy.m`,
  `run_sundman_tail.m`, `run_sundman_minfuel.m`,
  `test_sundman_eps1_noresample.m`) remain in the parent as the development
  record / diagnostics.
- tf-continuation: `tf_continuation_minfuel.m`, `tf_continuation_minfuel_fine.m`.
- Movies + solutions: `movie/` (min-fuel leg, coarse 6-switch, 53-switch,
  min-energy solo, three-way comparison; MP4 + GIF each).
