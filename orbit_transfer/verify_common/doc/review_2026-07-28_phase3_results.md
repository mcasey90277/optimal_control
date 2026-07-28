# Phase 3 results review (2026-07-28) — BOTH REVIEWERS: NOT DEFENSIBLE

Reviewed AFTER the measurements, with results and figures attached.
Both independently name the same weakest link: E = k*sum(h^3) extrapolated
from two points, contradicted by the study's own measurement that the worst
cell got 8x WORSE under redistribution (so the local coefficient is
phase-dependent and cannot be factored out).

---

## GPT-5.6-terra

- **[BLOCKING]** `“19 → 19 switches verifies one branch.”` -- No. The stated history already supplies a counterexample: 75 switches persisted across a real branch change. Switch count is only a coarse topology statistic; switch ordering, locations, arc identities, state/costate trajectory, and limiting continuous extremal may differ. Make this sound by continuation from one mesh to the other in small placement steps, with multiple seeds, matched switching sequence/times, and trajectory/costate proximity checks.

- **[BLOCKING]** `E = k Σh³; therefore E_production = 0.233 kg.` -- Two objective values identify only the *difference* under one redistribution, conditional on both solves representing the same continuous branch. They cannot establish a zero-intercept, linear, phase-independent global-error model. Plausible alternatives (`E∝(Σh³)^p`, a nonzero asymptote, or weighted `Σ C_k h_k³`) yield materially different extrapolations; even varying only \(p\) changes the inferred production error substantially. Falsify/validate with at least three same-branch meshes/redistributions spanning several values of `Σh³`, preferably including targeted phase redistributions, and fit competing models with uncertainty.

- **[BLOCKING]** `Constancy of median R_x/h³ supports the objective model E=kΣh³.` -- A state-residual scaling does not establish objective-error scaling. The missing quantity is the phase-dependent objective sensitivity and signed aggregation of defects. Moreover, the measured 8× worsening of the maximum residual directly rejects the constant-coefficient premise required to replace `Σ C_k h_k³` by `kΣh³`. The total mass may improve despite that, but it does not rescue the unweighted model.

- **[IMPORTANT]** `The residual is the trapezoid local truncation error.` -- It is consistent with a local discretization residual, but the evidence does not uniquely identify it as trapezoidal LTE. It can also include control reconstruction, switch placement, independent-variable conversion, and continuous evaluation/interpolation effects. A tiny discrete defect versus a larger continuous residual demonstrates a transcription-versus-continuous mismatch, not its exclusive cause.

- **[IMPORTANT]** `The three-row R_x/h³ constancy independently validates h³ scaling.` -- The rows deliberately preserve approximately 25 nodes/revolution, so they repeatedly sample closely related orbital phases with closely related step sizes. That is useful reproducibility evidence, but weak independence from orbital geometry. Test within smooth arcs using a log-log residual-versus-`h` slope; then repeat with materially different nodes/revolution, shifted longitude grids, and a controlled multiphase/fixed-switch formulation that preserves the branch.

- **[IMPORTANT]** `The residual is not driven by bang-bang structure.` -- Switch cells being 3.6× larger than interior cells means switching does contribute. A single largest interior residual does not show that bang-bang structure is unimportant to total objective error, particularly because defect multipliers may weight switch cells very differently. Evaluate signed multiplier-weighted contributions by cell class and hold the switch schedule fixed while refining smooth arcs.

- **[IMPORTANT]** `dJ = Σ λ_defᵀ d is a calibrated first-order mesh-error estimator.` -- The formula is conceptually valid only when `d` is precisely a perturbation of the NLP constraints under the same variables, scaling, quadrature, endpoints, and active set. Redistribution changes more than defect right-hand sides: it changes node locations, quadrature, control parametrization, and potentially the discrete extremal. A sign convention can explain a sign reversal, but cannot explain the 5.3× magnitude miss. The likely issue is model mismatch between the continuous residual and the actual perturbation induced by changing the transcription, compounded by finite/nonlinear reoptimization effects.

- **[IMPORTANT]** `The 5.3× miss means the sensitivity estimate is still useful for order of magnitude.` -- One failed controlled test establishes that it is uncalibrated for this perturbation; “order of magnitude” is not a validation criterion and can fail under other redistributions. Validate it with several small, signed perturbations along a continuation path and verify convergence of predicted/observed directional derivatives as perturbation size tends to zero.

- **[IMPORTANT]** `Voiding the 2.5 N and 1 N comparisons is correct.` -- Correct for claiming a pure mesh error from those two end points. Their raw mass changes combine mesh placement, basin/branch selection, and possibly switch-topology changes and should not be converted into error bars. A defensible route is continuation in placement with branch tracking, or a multiphase/fixed-switch-sequence formulation; otherwise report only the observed mixed effect.

- **[IMPORTANT]** `Basin fragility increases with switch count.` -- Not established. Two failed end-point comparisons at deeper rows are confounded by longer transfers, more opportunities for a switch displacement to cross a node, different basin geometry, and use of the same redistribution rather than a perturbation normalized by local mesh scale. Measure failure probability versus normalized node displacement across several rows and multiple randomized/graded redistributions before making this trend claim.

- **[BLOCKING]** `Basin selection dominates discretization by ~4.4× at 10 N.` -- The numerator is one observed improvement from a seed change and is at most evidence of a basin opportunity on that mesh; the denominator is a one-model, two-point extrapolation whose central assumption is contradicted by phase-dependent residual behavior. The ratio is therefore not a supported quantitative dominance claim. At most: the observed 1.015 kg basin improvement exceeds the *conditional* 0.233 kg estimate from the unweighted linear model.

- **[BLOCKING]** `The production mesh is good enough for mass and delta-V.` -- The evidence supports only that one 10 N redistribution produced a 0.146 kg mass difference. It does not establish a production-mesh error bound, cover the other thrust rows, or demonstrate delta-V robustness. “Good enough” needs a stated acceptance tolerance, repeated same-branch mesh perturbations, uncertainty bounds, and direct delta-V comparisons.

- **[IMPORTANT]** `The production mesh is not good enough for switch counts and times.` -- Topology changes under larger redistributions show sensitivity, but not yet an error magnitude in switch count or time. The 10 N case retained its count; the deeper cases cannot separate branch selection from mesh sensitivity. Report switch sensitivity as unresolved unless switch-time changes and topology reproducibility are measured under tracked continuation.

- **[IMPORTANT]** `Passing the PMP/self-gates and nearly constant λ_t shows a genuine continuous extremal.` -- These are necessary-condition consistency checks, not proof of a continuous extremal, local optimality, or branch identity. They are valuable verification of the checker and of residual evaluation, but the prose occasionally promotes them beyond what necessary conditions establish.

- **[MINOR]** `The 7×10¹¹ defect/residual ratio measures discretization error directly.` -- The two values are residuals of different equations/scalings and possibly different norms. The contrast is diagnostically useful, but its numerical ratio has no clear physical or convergence interpretation without normalization and a derivation linking them.

**VERDICT:** **NOT DEFENSIBLE.** The single weakest link is the conversion from one same-count redistribution difference to `E = k Σh³` and hence to **0.233 kg**: phase-dependent residual coefficients already contradict the unweighted model, while two points cannot distinguish it from plausible alternatives.

---

## Gemini 3.1 Pro

Here is the verification and validation review of the PMP residual results, focusing strictly on the numerical inferences and conceptual soundness of the claims.

- **[BLOCKING]** `Switch count 19 in both proves the two 10 N solutions are on the same branch.` -- An integer count does not uniquely define a discrete bang-bang topology (e.g., 9 thrust arcs + 10 coast arcs could shuffle their order or drop an arc while splitting another). The study's own history explicitly noted a case where 75 switches were held *across* a genuine branch change.
  **What would make it sound:** Verify the exact sequence of switch types (bang-to-coast, coast-to-bang) and confirm their orbital phases (true longitude) match monotonically between the two meshes.

- **[BLOCKING]** `Production-mesh discretization error is 0.233 kg, derived from E = k * sum(h^3).` -- This factors a single global constant $k$ out of the sum, implying the local truncation coefficient $C_i$ in $E \approx \sum C_i h_i^3$ is phase-independent. But in the same section, redistributing the nodes caused the maximum residual to get **8x worse**, directly proving that $C_i$ is highly phase-dependent (peaking at perigee/apogee geometry). Because moving the nodes changes the weighting of a highly variable $C_i$, the apparent $k$ is just an artifact of the two specific grids chosen and cannot be extrapolated to the continuous limit. 
  **What would make it sound:** Use at least a 3-point aligned ladder strictly verified to be on the same branch to fit an actual convergence order, or formulate an error model that integrates the known phase-dependent dynamics $C(L)$.

- **[BLOCKING]** `The first-order sensitivity estimate overpredicts by 5.3x, making it an order-of-magnitude bound.` -- In continuous optimal control, $\Delta J \approx \sum \lambda_{def}^T d_k$ is mathematically exact to first-order. A 5x discrepancy is not a loose bound; it is a fundamental violation of the estimator's premise. Because the optimal control is bang-bang and the grid is unaligned, the flow is discontinuous and the objective function is **non-smooth** with respect to node placement. Using smooth adjoint sensitivities across non-aligned switches calculates a mathematically fictitious gradient.
  **What would make it sound:** Abandon the first-order summation for estimating mass error in unaligned bang-bang collocation. It is theoretically invalid here.

- **[IMPORTANT]** `Constant median R_x/h^3 proves the residual is tracking local truncation, not switch structure.` -- The campaign strictly fixes the nodes-per-revolution ratio and uses true longitude as the independent variable. This guarantees the exact same orbital geometries are sampled at the exact same relative phases across all three rows. The constancy of the median is therefore a near-tautology measuring the smooth Keplerian dynamics of Earth orbit, not evaluating independence from thrust parameters. Furthermore, using a *median* dilutes the sparse but massive $O(1)$ jump errors at switch crossings (proven right in P3.1) into the far larger population of smooth interior flight steps.
  **What would make it sound:** Use the *mean* or integrated absolute residual—which captures the switch spikes—rather than the median, and test $R_x / h^3$ across runs with *differing* nodes-per-revolution.

- **[IMPORTANT]** `Basin selection dominates discretization error by ~4.4x at 10 N.` -- This compares a cleanly measured physical basin gain (+1.015 kg) against a theoretical mesh error (+0.233 kg) derived from an assumption ($k \sum h^3$) that the data itself falsified. Furthermore, it generalizes a sweeping campaign-level conclusion from $N=1$ row, side-stepping the fact that the meshes with higher switch counts failed the test entirely. 
  **What would make it sound:** Retract the specific 4.4x multiplier and state only that basin gains are the same order of magnitude as the *uncalibrated* mesh variations. 

- **[IMPORTANT]** `Basin fragility INCREASES with switch count (because deeper rows voided).` -- A 10 N trajectory has 19 switches and wide arcs; a 1 N trajectory has 171 switches and extremely narrow arcs. Applying the exact same global perturbation (a violent shift from uniform-longitude to uniform-time) is a profoundly larger topological shock to the 1 N mesh than the 10 N mesh. The voiding is an artifact of scaling. 
  **What would make it sound:** Normalize the grid perturbation magnitude against the duration of the shortest thrust arc in each row.

- **[MINOR]** `The production mesh is good enough for mass and delta-V but NOT for switch counts and times.` -- Mass error is literally the integrated absolute difference in control profiles. If the switch times are structurally wrong, the mass is necessarily suboptimal. Claiming mass is "good enough" while admitting switches are completely wrong requires proving that the objective function is exceptionally flat. Given that the first-order gradient check failed by 5x, objective flatness has not been rigorously established.
  **What would make it sound:** Reframe as: "The aggregate fuel use is relatively insensitive to exact switch timing, but the discrete times cannot be used for operational commanding."

### VERDICT
**NOT DEFENSIBLE.** 

**Single weakest link:** The assumption that $E = k \sum h^3$ permits a 2-point extrapolation to zero. The test affirmatively proved that the local error coefficient is wildly phase-dependent (spiking 8x upon redistribution), which prohibits factoring a global constant $k$ out of the discrete sum. The 0.233 kg error bar is a mathematical artifact of the two grids compared, not a true measurement of the absolute truncation error.
