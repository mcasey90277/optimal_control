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
3. **Independent verification** — check the 25-switch structure against the PMP
   switching function S = 1 − ‖λ_v‖c/m − λ_m sign changes (recover costates from
   the converged NLP multipliers) as a cross-solver certificate.
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
