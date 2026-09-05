## Findings

### PMP, smoothing, and propagation

- **[CORRECTNESS]** `cr3bp_minfuel_prop.m:65-86`, `cr3bp_minfuel_pmp.m:123-143` — **The Huber STM omits switch-time sensitivity.** For \(p<1\), the field jumps at \(h(y)=Q(y)-1=0\). Integrating \(\dot\Phi=A\Phi\) with the branchwise AD Jacobian does not capture the variation of the crossing time. At a transverse crossing, the required saltation update is
  \[
  \Phi^+=\left[I+\frac{(F^+-F^-)\nabla h^\top}
                              {\nabla h^\top F^-}\right]\Phi^- .
  \]
  Here \(F^-\) and \(F^+\) denote the fields immediately before and after the crossing. Step rejection does **not** supply this update. Consequently, the Huber shooting Jacobian is generally wrong on switching trajectories. This is an alternative explanation for stalled Newton residuals, directly undermining the claimed structural knockout.
  
  **Fix:** locate \(Q=1\) crossings with events, integrate smooth branches separately, and apply saltation updates. Detect grazing crossings separately when the denominator approaches zero. Validate the resulting endpoint STM against perturb-and-reintegrate sensitivities before rerunning the race. A discontinuous field does not, by itself, impose an unavoidable residual floor.

- **[ROBUSTNESS]** `cr3bp_minfuel_pmp.m:119-121`, `cr3bp_minfuel_pmp.m:141-143`, `cr3bp_minfuel_prop.m:73-74` — **The implemented epsilon family is not the stated tanh/logistic family, and it is not globally smooth.** With
  \[
  S=1-T_{\max}\left(\frac{\|\lambda_v\|}{m}+\frac{\lambda_m}{c}\right),
  \]
  the implemented law is
  \[
  \delta=\operatorname{clip}\!\left(\frac12-\frac{S}{2p},0,1\right),
  \qquad L=(1-p)\delta+p\delta^2 .
  \]
  This is a clipped-linear throttle, not
  \(\tfrac12[1+\tanh(-S/p)]\). The propagator uses this same field; there is no separate logistic implementation. The implemented running cost matches the clipped law, but the “exact Jacobian” claim needs qualification: a classical Jacobian generally does not exist at the clipping knees \(Q=1\pm p\). AD includes control derivatives away from those knees; that does not make the field globally smooth.
  
  **Fix:** identify the family explicitly as quadratic-cost/clipped-linear regularization in the experiment specification and certificate metadata. Add tests spanning both clipping knees, and qualify second-order results for this piecewise-smooth, bounded-control problem rather than treating it as a globally smooth logistic problem.

- **[CORRECTNESS]** `cr3bp_minfuel_pmp.m:58-62`, `cr3bp_minfuel_pmp.m:95-96`, `cr3bp_minfuel_pmp.m:126-128` — **Documented, accepted Huber parameters can violate the throttle bound.** The interface permits every \(\kappa>0\), but the low-\(Q\) branch clips to \(\kappa\), not to 1. For example, \(\kappa=2,\ Q=0.75\) produces \(\delta=1.5\). The correct minimizer in this example is \(\delta=1\).
  
  **Fix:** either restrict Huber parameters to \(0<\kappa\le1\), or support the documented range using
  ```matlab
  s = if_else(Q > 1, 1, fmin(fmax(ps*Q, 0), 1));
  ```
  Also require a finite scalar parameter and enforce the documented \(0<p\le1\) range for the epsilon family.

- **[CORRECTNESS]** `cr3bp_minfuel_prop.m:13-14`, `cr3bp_minfuel_prop.m:68-76` — **An incomplete integration can be returned as a completed segment.** MATLAB ODE solvers can issue a warning and return a partial trajectory when tolerances cannot be met. This function unconditionally takes the last returned sample as the state at `dt`; it never checks that `T(end)` reached `dt`. The header’s promise that integrator collapse throws is therefore unsafe. A shooting solver could enforce continuity against the wrong propagation duration.
  
  **Fix:** after either `ode113` call, verify completion to a roundoff-scale time tolerance and verify finite state/STM outputs. Throw a named error on incomplete integration so the shooting solver can reject the iterate.

- **[ROBUSTNESS]** `cr3bp_minfuel_pmp.m:114-116`, `cr3bp_minenergy_pmp.m:87-91` — **The primer regularization does not provide a physically admissible singular-arc treatment.** At \(\lambda_v=0\), `alpha` is the zero vector. If the mass-costate term requests positive throttle, the implementation consumes propellant while applying zero thrust acceleration. That is not a unit-direction control. Near this point, the direction derivative also becomes enormous; adding \(10^{-300}\) does not make the STM trustworthy.
  
  **Fix:** detect small primer norm on positive-throttle arcs and reject or classify the trajectory as requiring singular-control analysis. A genuine singular arc requires an admissible unit direction obtained from higher-order conditions, not a zero-vector fallback. Treat zero-throttle arcs separately, where direction is irrelevant.

### Conjugate-point test

- **[CORRECTNESS]** `ms_conjugate_test.m:116-136` — **The endpoint exclusion discards a genuinely interior crossing.** Samples stop at `tGrid(K)`, while the final time is `tGrid(K+1)`. Nevertheless, a flip ending at the last sample is marked `atFinal` and subtracted from `nCrossings`. Both ends of that bracket are strictly interior. For example, with three segments and opposite determinant signs at the two sampled junctions, the function returns `pass = true` despite detecting an interior crossing.
  
  The entire final segment is also unexamined.
  
  **Fix:** count every crossing between the currently sampled junctions as interior. Include propagation and sampling through the final segment, refine each crossing bracket, and classify a root as a final-endpoint case only when its estimated time coincides with the actual `tGrid(end)` within a specified tolerance.

- **[CORRECTNESS]** `ms_conjugate_test.m:128-136`, `ms_minfuel.m:67-72`, `cr3bp_minfuel_pmp.m:119-120` — **Structural rank deficiency on an initial coast is incorrectly counted as a conjugate-point refutation.** On a strict coast, \(\delta=0\) in a neighborhood of the extremal, so the state dynamics do not depend on the costates. If the trajectory initially coasts, then throughout that initial interval
  \[
  \Phi_{1:6,\,8:14}=0.
  \]
  The monitored seven-by-seven matrix therefore has rank at most one, with only the terminal-mass-costate row potentially nonzero. Its determinant is identically zero—not an isolated conjugate crossing proving loss of minimality. The code counts every zero sample as a focal point. Initial saturated arcs can likewise retain structural degeneracies.
  
  **Fix:** distinguish initial structural singularity from subsequent focal intersections; do not issue a refutation merely because the bounded-control endpoint map has not attained full rank. Use an active-set/switching-time second-variation treatment for production certification, and add analytic initial-coast and saturation tests. Return an indeterminate result when the assumptions needed for the regular Jacobi test are not established. An unconstrained LQ test does not exercise this failure mode.

- **[ROBUSTNESS]** `ms_conjugate_test.m:114-133` — **The determinant calculation is neither numerically stabilized nor an adequate absence-of-crossings test.** Taking the \(m\)-th root *after* `det(M)` only changes presentation; it does not prevent underflow, overflow, or sign corruption from ill-conditioned cumulative STM products. There are no finiteness checks. Junction-only sign comparisons also miss even-multiplicity zeros and crossing pairs between samples. Moreover, a zero sample between opposite signs can be counted both as a zero and as a sign flip.
  
  **Fix:** use stabilized variational propagation, dimensionally scaled matrices, and oriented signed-log-determinant calculations with singular-value diagnostics. Adaptively refine suspicious intervals and deduplicate root brackets. Reject nonfinite results and insufficient sampling as indeterminate, not PASS. Add tests for even-multiplicity crossings, closely spaced pairs, badly scaled blocks, and roots exactly at sampled junctions.

### Catalog integrity and certification

- **[CORRECTNESS]** `build_minfuel_catalog.m:31-43`, `build_minfuel_catalog.m:68-89` — **The packager does not enforce the advertised conjugate production gate.** If the rescue file is missing or its verdict is false, the original refuted record remains in `G`. Every record is then marked `has_solution = true` and appended to the catalog regardless of whether its verdict is PASS, FAIL, or missing (`-1`). The absent schema validator would have to provide the entire remaining defense; this function itself constructs those entries as solutions.
  
  **Fix:** resolve and validate the verdict before setting `has_solution` or appending costates/junctions. Require an explicit PASS plus a successful numerical solve. Abort or omit a failed/missing-verdict entry. For this particular advertised artifact, require the replacement and verify that exactly one original record was replaced.

- **[CORRECTNESS]** `build_minfuel_catalog.m:33-39`, `build_minfuel_catalog.m:80-88` — **Verdicts are not bound to the trajectory being packaged.** The lookup matches only objective kind, cell, and gamma. It does not match epsilon, junction states, time grid, dynamics constants, or solver version. Two different branches at the same cell and gamma can therefore receive the same verdict. The rescue path reduces the evidence to a Boolean override. Schema validation cannot detect a physically stale verdict from these keys.
  
  **Fix:** attach each verdict to a hash of the exact trajectory, time grid, smoothing parameters, boundary conditions, constants, and relevant code version. Match that identity before packaging, or recompute the test on the exact packaged trajectory. Store the actual diagnostics and crossing brackets, not just `conjPass`.

- **[ROBUSTNESS]** `build_minfuel_catalog.m:48`, `build_minfuel_catalog.m:74-78`, `ms_minfuel.m:19-20` — **Terminal reconstruction silently assumes a common uniform shooting mesh.** The solver accepts an explicit `seed.tGrid`, but the packager obtains `K` from the first record and propagates every final segment for `rec.tf/K`. For a valid nonuniform mesh, that is the wrong duration. Different mesh sizes also break the fixed array layout. The packaged junctions lack an accompanying per-entry time grid.
  
  **Fix:** retain the actual converged time grid for every entry and use `tGrid(end)-tGrid(end-1)` for terminal propagation. Store those grids with `Yj`. If the schema intentionally supports only common uniform meshes, assert that invariant for every record before packaging.

- **[ROBUSTNESS]** `build_minfuel_catalog.m:69-78`, `build_minfuel_catalog.m:110-113` — **The packager stores independently sourced physical quantities without checking their consistency.** It copies `g.mfFuel` while separately propagating a terminal state, copies `g.z` independently of `g.Y(:,1)`, and takes global thruster constants from the first pilot record even though propagation uses each record’s constants. Hard-coded dimensional units and Isp introduce another unchecked consistency requirement.
  
  **Fix:** assert that stored final mass agrees with `yT(7)`, initial costates agree with `g.Y(8:14,1)`, terminal position/velocity and \(\lambda_m\) meet tolerances, and all records share the catalog’s constants. In particular, enforce
  \[
  c_{\rm nd}
   = I_{\rm sp,s}\,g_{0,\rm km/s^2}\,
     \frac{TU_{\rm s}}{LU_{\rm km}}.
  \]
  Obtain dimensional constants from recorded provenance rather than independent literals. Recheck segment continuity before treating the packaged junctions as a certified record.

- **[ROBUSTNESS]** `ms_minfuel.m:75-82` — **The reported diagnostics undersample exactly the dynamics being used for certification.** `Hdrift` measures only differences between junction Hamiltonians; it can miss integration errors between nodes. `coastFrac` is an unweighted fraction of junction starts below a throttle threshold, not a time fraction. It can change solely through mesh redistribution and can badly alias short coast/burn arcs.
  
  **Fix:** evaluate Hamiltonian conservation on propagated interior samples, including the terminal point. Compute coast fraction by time integration,
  \[
  f_{\rm coast}
   =\frac1{t_f}\int_0^{t_f}
       \mathbf1_{\{\delta(t)<10^{-3}\}}\,dt,
  \]
  preferably using threshold-crossing times. Retain the present statistic only under an explicit name such as `coastJunctionFrac`.

- **[ROBUSTNESS]** `build_minfuel_catalog.m:95-107`, `ms_minfuel.m:70-72` — **A finite-epsilon verdict is not a certificate for the exact bang-bang fuel problem.** The packaged trajectories solve the finite-\(p\) quadratic regularization, and their STMs correspond to that same field. Even a repaired Jacobi test would apply to that finite-\(p\) extremal, not automatically to the \(p=0\) switching structure. The metadata acknowledges finite epsilon but then asserts final-mass convergence without recording an error estimate.
  
  **Fix:** describe entries as finite-epsilon approximations and store neighboring-epsilon mass changes, switching-time changes, and convergence diagnostics. Before claiming an exact min-fuel answer to a stated precision, perform a bang-bang switching-time refinement and appropriate sensitivity/second-order checks, or provide an independently justified error bound.

### Documentation

- **[STYLE]** `cr3bp_minfuel_pmp.m:2-79`, `cr3bp_minenergy_pmp.m:2-58`, `cr3bp_minfuel_prop.m:2-44`, `ms_conjugate_test.m:2-81` — **These headers lack traceable mathematical references.** Family names, project-local conventions, and an author/software mention do not identify the derivation or assumptions behind the regularization and mixed-boundary conjugate test.
  
  **Fix:** add explicit references with equation or section identifiers for the chosen regularization, CR3BP PMP conventions, discontinuous-flow sensitivities, and the fixed-time mixed terminal-condition test. Document where smoothness, controllability, and active-control assumptions are required.

## Claims not established

- **18 — “Epsilon ladder wins by knockout; Huber is structurally unsuited”: not established.** The Huber comparison uses an STM missing saltation sensitivity. A Newton residual floor with that Jacobian does not establish an unavoidable structural floor. The source also implements a clipped-linear, not logistic, epsilon law. Race data and driver logic are absent.

- **19 — “Seven certified min-fuel records”: not established as a certification claim.** The propagator can accept incomplete integrations; Hamiltonian and coast diagnostics are junction-only; finite-epsilon results do not themselves establish the exact fuel limit. The reported masses and continuation histories cannot be checked from this bundle.

- **20 — “13/14 conjugate verdicts, with one valid refutation”: not established.** The test explicitly discards one interior crossing bracket, omits the final segment, can mistake active-set structural degeneracy for focal points, and can miss crossings. These defects can produce both false passes and false refutations. The particular gamma-anomaly explanation therefore remains unverified.

- **21 — “Rescue validated and production-gated schema-v3 catalog shipped”: not established.** The rescue continuation is absent. The packager does not enforce PASS before inclusion, does not bind verdicts to trajectory identity, and assumes unverified mesh/physical-data consistency. Schema compatibility and rocket-equation derivation are outside the supplied source.

## QUESTIONS / COULD NOT VERIFY

- `catalog_schema.m` was not supplied. Its v1/v2 compatibility, required v3 fields, treatment of failed/unknown verdicts, and implementation of
  \[
  \Delta v_{\rm km/s}
    =I_{\rm sp,s}\,(0.00980665)\log(1/m_{f,\rm frac})
  \]
  cannot be verified.

- `ms_bvp.m`, the race/grid drivers, the rescue driver, and their result files were not supplied. Monotone epsilon stepping, retry seed reuse, final-depth acceptance, fixed-epsilon gamma continuation, node bookkeeping, and the reported numerical outcomes remain unverified.

- The claim that Huber \(\kappa=1\) reproduces the energy seed needs a normalization check in the missing driver. Since its objective is half the energy objective, the equivalent normal costates must also be halved:
  \[
  \lambda_{\rm Huber,\kappa=1}
   =\tfrac12\lambda_{\rm energy}.
  \]
  Reusing identical numerical costates does not reproduce the same throttle in general.

- The analytic conjugate-test tests and verdict artifacts were not supplied. In particular, it is unknown whether they cover the discarded last-interior bracket, final-segment roots, initial coasts, saturated arcs, and multiple crossings.