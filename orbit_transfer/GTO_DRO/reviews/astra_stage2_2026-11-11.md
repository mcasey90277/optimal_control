# Astra review, GTO -> DRO stage 2 recipe, 2026-11-11 (GPT-6 Astra, raw API, 132 s, $0.27; verbatim)

## Verdict

**Do not launch this recipe unchanged as a catalog of minimum-time admissible transfers.** It currently produces a mixture of infeasible extremals, conjugate-refuted extremals, and candidates that pass some necessary tests. Those are different products and need separate acceptance states.

The pilot has already demonstrated two structural failures:

- Unconstrained refinement can discard the direct problem’s feasibility.
- First-success continuation can select branches roughly 12% slower than another available branch.

Neither is a rare numerical corner case. More compute with the same selection rules will reproduce them.

## Critical findings

### C1. The direct solve and indirect refinement solve different optimal-control problems

**Concrete failure:** A node-constrained direct trajectory seeds an unconstrained extremal that penetrates Earth. The refinement is not preserving admissibility, and its convergence cannot certify the original constrained problem.

#### Your interpretation is partly right, but too strong

If an Earth constraint genuinely influences the constrained optimum—with a nonzero state-constraint multiplier—then the usual unconstrained PMP equations generally do not describe that optimum. State-constraint multipliers, junction conditions, and possibly boundary-arc conditions are missing.

But the following conclusions **do not follow**:

- That there is no nearby admissible unconstrained extremal.
- That the refined trajectory is in the same homotopy class.
- That the refined trajectory is the unconstrained *minimum*, rather than merely an extremal.
- That the cell has no admissible transfer.
- That touching the floor alone establishes a nonzero multiplier.

**Demote the particular infeasible trajectory, not the cell.** Treating a demotion as the final verdict on a cell would discard potentially admissible alternative branches or constrained optima.

#### Cheap measurement before the fleet

Retain and inspect, for every relevant direct solve:

1. Earth and Moon constraint slacks at all nodes.
2. Their IPOPT multipliers, with objective, state, constraint, and time scaling documented.
3. NLP primal, dual, and complementarity residuals.
4. The transcription trajectory, controls, mesh, and defects.
5. Independently propagated clearance under the reconstructed direct control.

For a representative set of violated and nonviolated cells:

- Refine the mesh near first perigee.
- Repeat with slightly perturbed Earth floors.
- Check whether the active set, multipliers, trajectory, and objective sensitivity stabilize.
- Compare finite-difference change in optimal time with the appropriately scaled multiplier prediction.

A large multiplier from one coarse, poorly propagated NLP is **evidence about that discretized NLP**, not proof about the continuous optimum. An active state constraint can also have isolated contact rather than a finite boundary arc; do not assume the latter without inspecting the trajectory.

#### What to change

Choose explicitly between two architectures:

**A. Constrained direct optimization is authoritative.**  
Refine the mesh, validate continuous feasibility and endpoint closure, and retain the resulting constrained solution. Use unconstrained indirect refinement only for trajectories with convincingly inactive path constraints.

**B. Indirect solutions are authoritative.**  
Implement the correct constrained PMP formulation where necessary. A positional altitude constraint has higher-order boundary conditions because thrust acts through acceleration. This is not merely adding an altitude equation to the existing shooting residual.

A finite penalty is not an exact substitute: it can converge while violating the floor. A barrier or penalty continuation can be useful numerically, but requires a limiting study and an independent hard-feasibility test.

**Before the fleet:** implement A, or explicitly restrict the product to admissible unconstrained extremals and record the excluded cells as unresolved. Full constrained indirect shooting can come later. Do not spend the fleet assuming a single unconstrained `z8` represents every constrained optimum; constrained extremals may require additional arc, junction, and multiplier data.

---

### C2. Continuous path feasibility is not established, including for the final repaired entries

**Concrete failure:** Node constraints and minima over integrator output samples can both miss periapsis violations. The current sweep also appears to omit the **Moon floor entirely**.

For a rotating-frame primary at fixed position \(r_b\), define

\[
q_b(t)=\|r(t)-r_b\|^2,\qquad
\dot q_b(t)=2(r(t)-r_b)\cdot v(t).
\]

Locate interior minima using roots of \(\dot q_b=0\), and include endpoints. Use dense output and root refinement. Because an event detector can miss multiple roots inside a step, also control step size near either primary and check dense-output minima over relevant steps.

**Required measurements:**

- Minimum Earth and Moon altitudes and their times.
- Convergence of those minima under tighter ODE tolerances and smaller maximum steps.
- An altitude uncertainty estimate, or a conservative numerical margin.
- Endpoint and clearance convergence together, not clearance alone.

There is no universal “refine factor” that certifies this. Require a convergence result comfortably below the allowed clearance margin. A reported minimum of 240 km gives only 40 km of margin over the Earth floor; that margin has meaning only after numerical uncertainty is bounded.

#### Does the moving phase band identify the mechanism?

It is consistent with the constraint-removal explanation: thrust changes first-return timing, energy, angular momentum, and the departure phases that encounter dangerous perigees. It does **not** establish that explanation uniquely.

At 10 km/s, 292 s corresponds to roughly 2,900 km of travel. A coarse mesh near perigee is independently suspect, although travel distance is not an altitude-error estimate.

Your strongest evidence is that **45/46 refined violations are visible at the original node times**. Thus, the refined violations are not principally a sampling artifact in the refined trajectory. That does **not** tell you whether the original direct trajectory also violated the floor between nodes.

To separate the mechanisms, compare:

1. The transcription’s continuous state representation.
2. Independent propagation of its reconstructed control.
3. The refined indirect trajectory.

All three need continuous clearance measurements.

A perigee-radius condition is useful only if it refers to actual extrema of the controlled trajectory. An osculating Keplerian perigee proxy is not equivalent in a thrusting CR3BP trajectory.

**Before the fleet:** implement this verifier for both primaries and run it after **every replacement or repair**, not only before the repair chain.

---

### C3. The acceptance gates do not establish a rendezvous PMP solution

**Concrete failure:** A trajectory can pass a 100 km position gate and a solver-reproducibility test while having incorrect terminal velocity, incorrect free-time transversality, or poor continuous shooting closure.

A fixed-state rendezvous needs independently checked **position and velocity** residuals. State whether terminal mass is free or constrained and test the corresponding condition.

For the usual normal, autonomous, free-final-time formulation with a time-independent fixed target state:

- Check the Hamiltonian condition at final time.
- Check Hamiltonian constancy along smooth unconstrained arcs.
- If terminal mass is free and there is no terminal mass cost, check the associated mass-costate transversality condition.
- Check state, costate, stationarity, and multiple-shooting junction residuals.

The exact equations depend on your sign convention and normalization; they must be documented.

`tfMin(z8)` returning within \(10^{-6}\) of `z8` is not enough:

- What norm and scaling are used?
- Costates and time have different meanings and scales.
- Costates depend on PMP normalization.
- A solver can terminate with a small step and a non-small residual.
- A third-party implementation can share the same modeling or convention error.

**Fix:** accept on independently evaluated physical and PMP residuals, not primarily on \(|\Delta z|\). Keep the latter as a reproducibility diagnostic.

Set separate, purpose-based tolerances for terminal position, velocity, time, mass, and clearance. **100 km is a numerical screening tolerance, not self-evidently a rendezvous tolerance.**

Also verify whether “fixed arrival state” means a genuinely fixed point in the rotating frame or a target moving along the DRO with arrival time. Those are different endpoint problems and have different free-time conditions.

---

### C4. Conjugate-refuted trajectories cannot ship as minimum-time solutions

**Concrete failure:** The catalog labels as optimal trajectories that fail a necessary local-optimality test.

Assuming the conjugate test is correctly formulated and numerically resolved, an interior conjugate point refutes the relevant local-minimum claim. For a strictly path-interior trajectory, sufficiently small variations remain path-feasible, so the original altitude restrictions generally do not rescue that claim.

For a trajectory touching an active state constraint, however, an **unconstrained** Jacobi test is not automatically the correct constrained test. The admissible variations and junction conditions differ.

Passing a conjugate test proves neither:

- Global optimality.
- Best branch for the cell.
- Uniqueness.
- Feasibility.
- Robustness.
- Optimality for a different endpoint, free-time, mass, or state-constraint formulation.

Depending on the theorem and hypotheses actually implemented, absence of conjugate points can support a particular local sufficiency result. “Pass” without that theorem and its numerical validation should mean only **“not refuted by this test.”**

**Fix:**

- Exclude the 60 refuted entries from the minimum-time candidate product.
- Retain them separately as feasible extremals if they pass feasibility, with the refutation clearly exposed.
- Make refutation a repair queue at **every rung**, not just 15 N.
- Validate the conjugate detector under integration refinement and perturbation.

For selected refutations, an especially useful cross-check is to construct a variation suggested by the Jacobi calculation, restore endpoint feasibility, and demonstrate a shorter transfer.

The hypothesis gates also need clarification:

- \(\min\|\lambda_v\|>10^{-6}\) is scale-dependent and requires costate normalization and continuous minimization, not just samples.
- Being small is not itself proof of a singular arc; being nonzero only establishes regularity of the thrust-direction law under the chosen formulation.
- Explain exactly how “abnormal-lift dimension 1” relates to “no abnormal lift.” A one-dimensional multiplier space can describe uniqueness up to scaling of a normal lift; an actual abnormal lift is something else.

These should not be waivable blockers for an optimal-transfer catalog.

---

### C5. The “faster and loses no lower rung” repair rule can deliberately retain an inferior top rung

**Concrete failure:** A better admissible 15 N solution is rejected because it is a less successful numerical seed for subsequent solves. That confuses optimization quality with continuation convenience.

A repaired trajectory can also be faster but unacceptable because it violates a floor, fails endpoint closure, has inadequate mass, or fails the relevant optimality test. Speed alone is not sufficient.

**Fix:** separate the solution archive from the seed archive.

For each cell and rung, retain:

- The fastest independently validated candidate.
- Other branch-distinct candidates.
- Seeds useful for continuation, including slower or refuted extremals.

Do not overwrite already solved lower rungs merely because a new top rung was found. Keep the old branch as a seed while accepting the better top-rung result. A ladder is a computational device, not a requirement that all catalog entries belong to one branch.

Your observed 12% improvements are direct evidence that first-success selection is not a credible minimum-time search policy.

## Important findings

### I1. Removing the direct miss gate is defensible only if you stop treating refinement as preservation of the direct solution

**Concrete failure:** A large correction can change revolution count, encounter sequence, clearance, cost, and local branch. Neighboring cells can acquire unrelated solutions. The direct NLP’s multipliers and optimality properties then say little about the refined result.

A hard 1,000 km pre-gate is not justified by your evidence: it discards useful seeds and does not eliminate Earth penetration. But a 43,000 km mismatch cannot be ignored diagnostically.

Measure why it occurs:

- Transcription defects and convergence with mesh refinement.
- Control parameterization and reconstruction.
- Time rescaling and nondimensionalization.
- State/interpolation consistency.
- Control-frame conventions.
- Multiple-shooting versus continuous-propagation consistency.
- Accumulated sensitivity on multirevolution arcs.

**Fix:** use a large direct miss as a **seed-quality and branch-change flag**, not necessarily a rejection. Record pre/post time of flight, encounters, revolution counts, clearance, and trajectory differences. A newly validated extremal may still be useful, but it is not a faithful refinement merely because the solver converged.

Covector harvesting also needs a documented, transcription-specific mapping, including scaling and the treatment of path multipliers. Generic NLP duals are not automatically continuous costates.

---

### I2. The remaining basin search is too narrow, especially at low thrust

**Concrete failure:** Reverse warming can return the same refuted branch indefinitely, while lower rungs accumulate refutations without any search for alternatives.

For the five identical returns, changing the seed within the same basin is unlikely to help. Try, in roughly this order:

1. Seeds from other neighboring departure and arrival phases.
2. Both higher- and lower-thrust branches, with smaller intermediate thrust steps.
3. A controlled perturbation related to the detected negative/Jacobi direction, followed by constrained endpoint restoration.
4. Alternative time-of-flight and revolution-count seeds.
5. Pseudo-arclength continuation near folds.
6. A small, budgeted multistart set.

Reverse continuation from a lower rung is a cheap candidate repair, **not a universal remedy**. The 5 N rung has no lower fleet rung, although an auxiliary 4 N or 3 N continuation branch might supply seeds.

Do not infer too much from the increase from \(10/59\) to \(20/54\): the solved populations differ. Measure refutation changes on matched cells and, where possible, matched branches.

Extra-revolution seeds are worth trying, but a longer trajectory is not automatically safer; it may add dangerous perigee passages.

---

### I3. Easy phase-continuation repairs indicate initialization dependence, not automatically a bug

**Concrete failure:** Cold-start failure is being mistaken for cell difficulty or lack of feasibility.

Eight repairs succeeding is compatible with a poor cold-start basin and a useful neighboring solution. It becomes suspicious only if the repaired problem was not actually changed correctly.

**Cheap checks:**

- Recompute the requested departure state independently.
- Check its position, velocity, orientation, and mass against the flown initial state.
- Check arrival phase and endpoint identity.
- Check phase wrapping and indexing.
- Check that continuation changes the endpoint constraints, not just the stored label.
- Confirm the allowed control discontinuity at ignition and the intended always-on departure model.

A displacement of \(1/12\) of an orbital period need not be small near perigee. Record the actual state displacement. Use smaller phase steps if necessary.

For a few cells, approach the same phase from both directions. Different accepted times expose branch dependence that one-way continuation hides.

---

### I4. The deep-rung “wall” is a solver observation, not a reachability result

**Concrete failure:** Feasible low-thrust transfers could be rejected because the `z8` single-shooting representation is ill-conditioned. Conversely, a nonconverged multiple-shooting result could be excused as sensitivity.

To distinguish them:

1. Report every multiple-shooting junction residual in scaled and physical units.
2. Propagate each segment independently from its own junction state and costate.
3. Compare the stitched trajectory with progressively longer continuous propagations.
4. Integrate variational equations and measure relevant sensitivity growth.
5. Tighten tolerances, increase stored precision, and test higher-precision propagation on selected cases.
6. Compare with an independently mesh-refined direct solution.

If small residuals and tiny initial perturbations are amplified enough to explain the endpoint miss, conditioning is implicated. If segment closure or the final endpoint is poor already, the BVP has not been solved.

Also test smaller thrust steps, additional shooting segments, adaptive meshes, and longer flight-time seeds. The six selected cells and unswept solutions cannot establish a general wall at 2–3 N.

If a mathematically credible trajectory cannot be reproduced reliably from the distributed `z8`, **`z8` is an inadequate catalog representation**. Store a segmented trajectory, junction data, and propagation requirements.

---

### I5. The mass model and “always-on min-time” premise need an explicit audit

**Concrete failure:** The solver may produce transfers outside the available propellant budget or solve a more restrictive problem than the advertised minimum-time problem.

At 15 N,

\[
\dot m=-\frac{15}{1710g_0}\approx-8.95\times10^{-4}\ {\rm kg/s}.
\]

Exhausting the entire 150 kg would take about 46.6 h; any dry mass makes the available burn shorter. Report final mass, propellant use, dry-mass assumptions, and the nondimensional time conversion for every entry.

“Minimum time” does not, by itself, justify fixing thrust at maximum throughout for every mass-depleting, constrained rendezvous problem. If fixed magnitude is a prescribed model restriction, say **“minimum time within the always-on fixed-magnitude model.”** If the intended problem permits throttling or coasting, derive and test the throttle optimality conditions.

Likewise, do not use thrust-time monotonicity as a theorem for equality-fixed, always-on thrust. With throttleable thrust bounds and otherwise identical constraints, feasible-set inclusion supplies a monotonicity argument; your fixed-magnitude formulation does not automatically have that inclusion.

---

### I6. The final chain needs exhaustive invalidation rules, not a three-entry audit

**Concrete failure:** A repair made after the Earth sweep can reintroduce Earth penetration, Moon penetration, endpoint errors, or stale conjugate flags.

Every accepted replacement must invalidate and recompute all trajectory-dependent gates. Associate results with a hash of the actual trajectory/`z8`, model, normalization, and solver settings.

The reported final minimum Earth altitude must come from the **final versions of all entries**. A three-entry audit cannot establish that. Gate evaluation is cheap relative to optimization; run it exhaustively.

## Minor findings and reporting corrections

### M1. The `tf0` calibration statement is numerically wrong and statistically weak

The range \(0.2589\)–\(0.3238\) is not within 0.02 of 0.30. More importantly, accepted cold-start times are selected by the very initialization whose adequacy you are trying to assess. They cannot validate it independently.

**Fix:** compare a small spread of `tf0` values on matched cells, measuring accepted feasible time, branch identity, failure rate, and runtime. The faster repaired solutions are already evidence that the original distribution contains branch-selection bias.

### M2. “Solved” conflates incompatible categories

Report separately:

- Direct NLP converged.
- Indirect BVP converged.
- Endpoint validated.
- Both-primary clearance validated.
- Mass validated.
- Conjugate test passed, refuted, or inconclusive.
- Best known among searched branches.

On current evidence, **211 is at most the number not refuted by the stated conjugate test**, not a count of certified minimum-time transfers.

### M3. One sheet does not calibrate fleet coverage or cost

The pilot does not establish behavior for other orientations, DRO sizes, near-Moon encounters, different branch folds, or their repair burden. Running the 180° sheets first merely discovers a generalization failure earlier.

**Fix:** run a small stratified preflight across all orientation/size combinations, including perigee departures and difficult arrival phases. Measure end-to-end runtime including validation and repair, not just the fleet engine. The unused 300 s cap establishes neither reliability nor a useful runtime bound.

## What I would do before committing the fleet

1. **Install a single exhaustive final verifier:** endpoint position/velocity, mass, both-primary continuous clearance, PMP residuals, and numerical convergence.
2. **Retain direct trajectories and duals.** Diagnose a small representative set of the 46 demotions using perigee mesh refinement and floor perturbations.
3. **Stop discarding constrained direct solutions merely because unconstrained refinement cannot represent them.**
4. **Separate best solutions from continuation seeds.** Remove “loses no lower rung” as a reason to reject a better independently validated top rung.
5. **Make conjugate refutations a repair queue at all five rungs.** Keep unresolved refutations outside the minimum-time product.
6. **Run a small cross-sheet preflight and a matched-cell multiseed experiment.**
7. **Preserve enough data to reproduce and repair results without repeating the fleet.**

Full constrained indirect boundary-arc machinery and deeper-rung exploration can follow. Continuous feasibility, correct acceptance residuals, final-stage revalidation, and honest optimality labels cannot.