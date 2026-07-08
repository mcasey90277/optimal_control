# Low-thrust min-fuel campaign — learnings & next steps

**Record of the July 7–8 2026 effort to solve the minimum-fuel GTO → tulip
transfer on the full ~40-revolution spiral, and what it would take to finish
the many-switch case.** This is the synthesis; per-piece detail is in
`MIN_ENERGY_NOTES.md`, `OVERNIGHT_STATUS.md`, and
`../lowThrust_GTO_tulip/gto_tulip_mintime_theory.pdf` §6.

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
- **Certified:** min-energy → 3-switch machine-tight full-spiral min-fuel
  (defect 2×10⁻¹⁵). Min-energy leg indirect (6×10⁻¹²). Min-energy full-spiral
  direct (3×10⁻⁴).
- **Demonstrated but NOT certified:** the sharp ~50–80-switch global min-fuel.
  We reached ~50 switches (coherent to ~10⁻³, ΔV ≈ 3.37 km/s, ~25% below
  min-time — *directionally* right, numbers soft). No machine-tight, sharp
  bang-bang many-switch solution exists yet.

## Next steps — the regularization sub-project
Getting a machine-tight sharp many-switch solution is a dedicated
numerical-methods effort (standard practice for serious many-rev low-thrust),
not a solver flag:

1. **Remove the perigee singularity.** Sundman time transform (dt = r dτ) to
   slow the parameterization near perigee, or Kustaanheimo–Stiefel (KS)
   regularization to eliminate the `1/r` singularity from the dynamics. This
   is the single highest-leverage change — it is what makes the exact Hessian
   well-conditioned and lets IPOPT work.
2. **Explicit scaling** of variables and constraints for IPOPT (nlp scaling
   is not enough on the raw problem).
3. **Adaptive mesh refinement** — concentrate nodes where switches live
   (Betts Ch. 4), rather than borrowing the min-time integrator's density.
4. **Multiple shooting** for the indirect side — partition the 40 revs so
   each segment's sensitivity is ~10⁶/40; then the min-energy → min-fuel
   homotopy can certify to machine precision.
5. Then re-run the tf-continuation (with the fixed cone guard) inside
   IPOPT on the regularized, scaled problem for the sharp many-switch result.

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
- tf-continuation: `tf_continuation_minfuel.m`, `tf_continuation_minfuel_fine.m`.
- Movies + solutions: `movie/` (min-fuel leg, coarse 6-switch, 53-switch,
  min-energy solo, three-way comparison; MP4 + GIF each).
