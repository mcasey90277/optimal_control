## Question 1 — Is it technically bulletproof?

**No.** Several central formulas are right, but this is a useful numerical audit of a candidate extremal—not a numerical certificate of strong local minimality. The largest problems are the unverified applicability of the reduced conjugate test, incomplete numerical coverage, and conflation of mathematical conditions with solver cross-checks.

### PMP formulation and computation

- **[CORRECTNESS]** `DRO_tulip/indirect/transfer_study.m:148–158` - **The Hamiltonian formula is correct under the stated minimum convention and normal normalization.** For the objective \(J=t_f=\int_0^{t_f}1\,dt\) in nondimensional time,
  \[
  H=1+\lambda_r^\top v+\lambda_v^\top
  \left(b(r)+Cv+\frac{T}{m}\alpha\right)-\lambda_m\frac{T}{c},
  \]
  with \(H\) minimized over the control. Thus `1 + Y(k,8:14)*F(1:7)` has the correct sign and state dimension. For this autonomous problem, free final time gives \(H(t_f)=0\), hence \(H\equiv0\) along an exact extremal. **Fix:** explicitly declare “normal multiplier \(p_0=+1\), minimum convention, objective = nondimensional elapsed time.” If the objective were dimensional seconds with the same multiplier normalization, its running cost would instead be `tStar`. These are equivalent optimization problems, but their normalized costates differ.

- **[MISSING-CONDITION]** `costate_common/ms_tfmin.m:61–79,100–105` - **The supplied source does not permit an end-to-end audit of the equations, shooting residual, or STMs.** `ms_bvp`, `mintime_prop_seg`, and the third-party `tfMinEoM`/`tfMinProp` implementations are unavailable. In particular, it is unknown whether their control set, Hamiltonian normalization, feedback derivatives, propagation tolerances, and stopping behavior match the claims. **Fix:** provide or pin those implementations and add contract tests: compare the library Hamiltonian with the explicit formula, its gradient with finite differences, and the segment STM with perturbed propagations. “All conventions match pumpkyn” is not verification that those conventions match the stated problem.

- **[MISSING-CONDITION]** `DRO_tulip/indirect/transfer_study.m:141–146` - **N1 does not independently check the adjoint equations, and a small BVP residual is not by itself “the statement that the first variation vanishes.”** It checks matching between integrations of the equations selected by the implementation. For this problem, with
  \[
  C=\begin{bmatrix}0&2&0\\-2&0&0\\0&0&0\end{bmatrix},
  \]
  the fixed-control adjoint equations are
  \[
  \dot\lambda_r=-D_rb^\top\lambda_v,\qquad
  \dot\lambda_v=-\lambda_r-C^\top\lambda_v,\qquad
  \dot\lambda_m=\frac{T}{m^2}\lambda_v^\top\alpha
              =-\frac{T}{m^2}\|\lambda_v\|.
  \]
  **Fix:** independently validate these equations and the state equations against the flown trajectory, preferably through an independent integration/defect calculation rather than noisy finite differences. Rename N1 “multiple-shooting matching and terminal residual.” Separately identify which PMP equations are enforced by construction and which are independently checked.

- **[MISSING-CONDITION]** `DRO_tulip/indirect/transfer_study.m:170–180` - **N5 is a normalization identity, not a test of the flown control or Hamiltonian minimality.** You construct `alpha = -lamv/|lamv|` and then check that its length is one. A propagator using the opposite direction could still pass this printed test. Moreover, N5 is omitted from `necessary`.

  For the sphere, the proposed law is indeed the **global** minimizer, not merely stationary. Writing \(\rho=\|\lambda_v\|>0\),
  \[
  H(\alpha)-H(\alpha_*)
    =\frac{T\rho}{2m}\|\alpha-\alpha_*\|^2,\qquad
  \alpha_*=-\lambda_v/\rho,\quad \|\alpha\|=1.
  \]
  **Fix:** retain this analytic explanation inline, but compare it with the control actually applied by the propagator. Log that control, or recover it from an independently validated trajectory derivative. Check sphere membership and the minimization gap
  \[
  \frac{T}{m}\left(\lambda_v^\top\alpha_{\rm flown}+\rho\right).
  \]
  Include the resulting Boolean in the PMP aggregate. Do not use `realmin` to conceal an undefined control at \(\rho=0\).

- **[MISSING-CONDITION]** `DRO_tulip/indirect/transfer_study.m:118–124,179–180` - **Basic admissibility and domain checks are missing.** Before using `log(1/mf)`, dividing by mass, or applying smooth optimal-control theory, establish finite real data, \(t_f>0\), successful propagation through the requested final time, \(m>0\), positive thrust/exhaust speed, and avoidance of the CR3BP singularities. For mandatory all-burn,
  \[
  m(t)=1-\frac{T}{c}t
  \]
  is an especially strong, cheap consistency check. **Fix:** add a feasibility block before the PMP block, including propagation coverage and mass-law residuals. If dry mass, collision radii, eclipses, or other mission constraints are intended, define and enforce them; otherwise explicitly state that they are excluded. Such constraints cannot be silently assumed satisfied.

- **[CORRECTNESS]** `DRO_tulip/indirect/transfer_study.m:79–90,148–168` - **The endpoint transversality conditions shown are correct for the actual fixed-phase problem, but “fixed endpoints” needs qualification.** Initial position, velocity, and mass are fixed; terminal position and velocity are fixed; terminal mass is free. Consequently \(\lambda_m(t_f)=0\) is correct for an interior, unconstrained terminal mass with no terminal mass cost. There are no additional zero-costate conditions for fixed position/velocity components. Also, `sD` and `sA` are parameters here, not optimization variables. **Fix:** state this endpoint manifold explicitly. If future scripts optimize either phase, add the corresponding endpoint tangent condition, such as
  \(\lambda_{rv}(t_f)^\top\,d\,stateA/ds=0\).
  If the target is time-dependent rather than a fixed rotating-frame state, replace the simple free-time condition by the appropriate moving-endpoint transversality condition.

- **[MISSING-CONDITION]** `DRO_tulip/indirect/transfer_study.m:141–173` - **The report never establishes whether the trajectory is one smooth all-burn arc, which determines whether corner conditions are relevant.** For a smooth sphere-control extremal with \(\rho>0\), there are no physical switching corners to check. Artificial multiple-shooting junctions are not physical corners: continuity of the full state and costate should already be enforced by the residual. **Fix:** declare and verify the smooth-arc assumption. If actual switching or multiphase junctions are allowed, add the applicable state/costate matching, Hamiltonian continuity at free interior switching times, and any constraint-induced jump conditions. Do not claim a generic “PMP checklist” while omitting that distinction.

- **[CORRECTNESS]** `DRO_tulip/indirect/transfer_study.m:118,143–144,207` - **The first- and second-order diagnostics are not all evaluated on the same numerical trajectory.** N1 and the conjugate test use the multiple-shooting junction solution and its STMs; N2–N5 and the lift-space test use single-shot propagation from `z8`. Small junction defects can amplify substantially in this problem—the solver’s own documentation emphasizes that. **Fix:** compare the flown augmented state against every shooting junction, with component scaling. Require that discrepancy, not merely terminal position error, to be within the error budget used for the second-order calculation. Recompute or validate the Jacobi propagation along the accepted common trajectory.

- **[CORRECTNESS]** `DRO_tulip/indirect/transfer_study.m:115–124,179–180` - **Solver convergence is printed but not enforced before expensive verification and optimality language.** A best iterate from an unsuccessful solve is still propagated, tested, and called “SOLVED.” The independently recomputed residual is useful, so `it.converged` need not be the sole authority—but disagreement needs handling. **Fix:** label this stage “solver result,” establish explicit numerical acceptance from finite residuals and trajectory consistency, and stop or mark subsequent tests `NOT EVALUATED` when acceptance fails. A convergence flag, a residual test, and a mathematical optimality condition are three different things.

- **[CORRECTNESS]** `DRO_tulip/indirect/transfer_study.m:74–89` - **`spline` does not construct a periodic spline.** It constructs a not-a-knot cubic spline; wrapping the argument with `mod` does not enforce periodic derivative matching and may also introduce a value jump if orbit closure is imperfect. **Fix:** use a genuinely periodic interpolation construction, or compute the selected endpoint by propagation to its phase. For the present fixed phases, the chosen interpolated vectors are still legitimate fixed endpoint data, so this does not automatically invalidate that BVP. It does invalidate the “PERIODIC SPLINE” description and leaves an unquantified difference between the specified endpoint and the actual periodic orbit.

### Legendre, switching, and normality

- **[CORRECTNESS]** `DRO_tulip/indirect/transfer_study.m:189–194` - **The strengthened Legendre formula is right, with important qualifications.** The ambient Euclidean Hessian in \(\alpha\) is zero because the Hamiltonian is affine in direction. The Hessian of the Hamiltonian **restricted to the sphere at the minimizing direction**, equivalently the constrained Lagrangian Hessian on \(T_\alpha S^2\), is
  \[
  \frac{T}{m}\rho\,I_2.
  \]
  Thus the stated reduction to \(\rho>0\) is valid only after establishing \(T>0\) and \(m>0\). **Fix:** print the actual coercivity measure `min((Tnd./mm).*nlv)`, its location, and an error-aware positive margin. A sampled floating-point value barely greater than zero is not numerical evidence of uniform strengthened Legendre.

- **[CORRECTNESS]** `DRO_tulip/indirect/transfer_study.m:196–199` - **The switching function has the correct sign for an optional-throttle problem, but it is not a control-selection condition for the stated mandatory all-burn sphere problem.** If throttle \(s\in[0,1]\) is admissible, minimizing over direction gives
  \[
  H=\text{drift terms}+1-TsQ,\qquad
  Q=\rho/m+\lambda_m/c.
  \]
  Then \(Q>0\) uniquely selects \(s=1\); \(Q=0\) permits any throttle, so even the helper’s strict “iff” wording needs qualification. If the only admissible control is a direction with \(s\equiv1\), there is no throttle decision.

  Moreover, for a consistent all-burn normal extremal,
  \[
  \lambda_m(t)=\int_t^{t_f}\frac{T\rho(\tau)}{m(\tau)^2}\,d\tau\ge0,
  \]
  so \(Q>0\) already follows from positive mass and \(\rho>0\). **Fix:** choose and declare the admissible control set. For the sphere model, label S2 a redundant consistency diagnostic. For a throttle-capable model, retain it as the active-set/regularity check and verify that the propagator actually follows that active set.

- **[CORRECTNESS]** `costate_common/mintime_hypothesis_gates.m:14–24,91–109` - **The lift-space construction is a useful strict-normality diagnostic, but its interpretation is overstated.** Given a genuine normal lift in \(S\) with \(\lambda^\top f=-1\), exact `dimS == 1` excludes a nonzero abnormal stationary lift satisfying the included endpoint conditions. That is a sound argument. However, \(S\) imposes collinearity, not the antiparallel sign required for Hamiltonian minimization. Consequently `dimS > 1` produces an abnormal **stationary** lift in the kernel of the Hamiltonian functional, not automatically an abnormal PMP-minimizing lift. **Fix:** describe this as a conservative linear-space test for absence of abnormal lifts, rather than an unconditional equivalence for full PMP lifts. Distinguish ordinary normality, already represented by \(p_0=1\), from strict normality/corank-one hypotheses.

- **[MISSING-CONDITION]** `DRO_tulip/indirect/transfer_study.m:201–204,211` - **The script ignores the self-consistency quantities on which the normality argument depends.** `gates.Hresid` and `gates.nullResid` are computed but neither printed nor gated. The independent fixed-control dynamics are useful here: if they disagree with the library propagation, the normality calculation may not concern the accepted lift at all. **Fix:** require an adequately resolved accepted member of the null space, a correctly normalized Hamiltonian, and a stable separation between the sixth and seventh singular values. Print the residuals, spectrum/gap, and actual rank threshold alongside `dimS`.

- **[CORRECTNESS]** `costate_common/lift_space_dim.m:10–18,46–49` - **The rank threshold is heuristic, and its stated safety guarantee is false.** The cap at `1e-3*sv(1)` does not make a noisy lift report zero nullity “never > 1”; several singular values can lie below that cap. Increasing the tolerance based on the residual can also turn model/integration error into an apparently acceptable null space. **Fix:** remove that guarantee. Establish an independent bound or convergence estimate for the constraint-matrix error, then demand that the numerical rank be unambiguous relative to it. Otherwise return `UNRESOLVED`, not a definitive integer classification.

- **[CORRECTNESS]** `costate_common/mintime_hypothesis_gates.m:73–109` - **The rank computation has substantial discretization error not controlled by its ODE tolerances.** It integrates an adjoint along separately PCHIP-interpolated state and costate data, then imposes alignment at only 200 times. Those interpolants are not necessarily an exact solution/control pair. Stacking more unweighted alignment rows also changes their weight relative to the single terminal condition. **Fix:** perform trajectory/interpolation/sample refinement studies, use appropriate quadrature and variable scaling for the constraint operator, and report rank stability across refinements. Tight `ode113` tolerances do not bound the error in the interpolated coefficients.

### BCT applicability and the conjugate calculation

- **[MISSING-CONDITION]** `DRO_tulip/indirect/transfer_study.m:183–211` - **“Bonnard-Caillau-Trélat” is not a sufficiently specified theorem invocation.** The report needs an exact theorem and a mapping from its hypotheses to this problem. Missing from that mapping are:
  - smoothness of the dynamics and minimized Hamiltonian in a neighborhood of the reference extremal;
  - exclusion of zero mass, CR3BP collisions, zero primer vector, and—if relevant—changes of throttle active set;
  - the correct endpoint manifold and associated regularity/corank conditions;
  - the precise normality requirement, including any requirement on restrictions to subarcs;
  - the appropriate Jacobi/focal problem and its short-time initialization;
  - absence of the relevant degeneracy throughout the required interval, including the terminal boundary;
  - the theorem’s topology of local optimality and any additional assumptions needed for its **strong**, rather than merely weak, conclusion.

  **Fix:** cite the theorem number and state these hypotheses explicitly. Some are structural facts of this model; others need tests. Four Boolean diagnostics are not a substitute for that correspondence.

- **[MISSING-CONDITION]** `costate_common/ms_tfmin.m:74–79` - **The equivalence of this reduced \(6\times6\) test to the required mixed-endpoint Jacobi problem is not established.** Excluding the mass-costate column and quotienting the position/velocity costate scale are plausible on a genuinely all-burn arc: neither affects the flown position/velocity trajectory. But that does not, by itself, prove that appending `flow6` gives the BCT test for fixed initial mass, free terminal mass, and free final time.

  Eliminating mass gives explicitly time-dependent six-state dynamics, since \(m(t)=1-Tt/c\). Also, terminal variations must respect the natural mass-costate condition as well as the free-time condition. **Fix:** derive the reduction from the full accessory problem. At the actual terminal point, verify its correspondence with
  \[
  \delta(rv)(t_f)+f_{rv}(t_f)\,\delta t_f=0,\qquad
  \delta\lambda_m(t_f)+\dot\lambda_m(t_f)\,\delta t_f=0,\qquad
  \delta H(t_f)=0,
  \]
  with fixed initial state. Do not “fix” this merely by adding a mass row: prove the appropriate reduction or use a validated mixed-boundary focal/second-variation method.

- **[CORRECTNESS]** `costate_common/ms_conjugate_test.m:168–247` - **Sign sampling cannot establish absence of conjugate times.** It misses two crossings within one segment and even-multiplicity zeros unless they happen to be sampled. With `zeroTol = 0`, almost any numerically represented touch will be missed. The helper honestly mentions some of this; the script discards that qualification. **Fix:** adaptively monitor singularity/conditioning, refine suspicious intervals, and locate candidate roots. For robust second-order work, use a suitable Jacobi-subspace/index method rather than determinant signs alone. Without validated interval coverage, report “no singularity detected,” not “no conjugate time.”

- **[CORRECTNESS]** `costate_common/ms_conjugate_test.m:196–204` - **Skipping everything before the first numerically full-rank sample can hide the very event being sought.** For this mandatory all-burn problem, there is no established initial coast. Small-time rank behavior can reflect controllability orders and conditioning, not absent control action; earlier singularities can also be skipped. **Fix:** justify a short-time interval free of relevant conjugate/focal points using the applicable theory or a controlled asymptotic/numerical initialization. Do not infer that justification from “a later sample became full rank.”

- **[CORRECTNESS]** `costate_common/ms_conjugate_test.m:160–165,201–247` - **The function can return `PASS` with inadequate coverage.** One live sample is enough; no bracket need be tested. A missing `Yend` silently omits the last segment while still allowing `PASS`. `tested` rejects `NaN` determinants but not all nonfinite values, and later near-rank-deficient samples do not automatically make the result unresolved. **Fix:** require explicit coverage through the requested final time, a justified starting interval, sufficient resolved samples, finite data, and a conditioning margin. Missing coverage must produce `UNDETERMINED`. The script must inspect coverage, not just `cj.pass`.

- **[CORRECTNESS]** `costate_common/ms_conjugate_test.m:166–193` - **Cumulative STM multiplication reintroduces the long-arc conditioning problem that motivated multiple shooting.** Equilibrating the final \(6\times6\) block does not recover information already lost while forming `PhiCum`. The positive row/column equilibration and LU determinant-sign calculation are sensible; they are not a cure for unstable variational propagation. **Fix:** use stabilized propagation of the relevant Jacobi subspace, such as continuous/repeated orthogonalization or an appropriate Riccati/symplectic method. Validate against tighter integrations and different segmentations, and check structural STM identities where applicable.

### Cross-check, thresholds, and verdict

- **[OVERCLAIM]** `DRO_tulip/indirect/transfer_study.m:21–23,175–180` - **The stated role of `tfMin` is correct; its placement in `necessary` is not.** Agreement with another shooting solver is an implementation cross-check, not a PMP condition. Conversely, disagreement does not establish that the candidate is not an extremal. These solvers also share propagation/dynamics and endpoint construction, so agreement cannot catch their common-mode errors. **Fix:** maintain separate aggregates such as `pmpSatisfied`, `crossCheckPassed`, and `secondOrderResolved`. Describe the comparison as independent **shooting-algorithm** agreement, not independent verification of the physical model.

- **[CORRECTNESS]** `DRO_tulip/indirect/verify_with_pumpkyn.m:59–71` - **`V.converged = true` means only “returned eight finite numbers.”** The wrapper obtains no convergence evidence from `tfMin`, and its flight check tests only position. Worse, `V.moved` is independent of the flight result, so the study’s N6 can pass even when the cross-check flight fails or never completes. **Fix:** separate `returnedUsableVector`, `solverConverged` if available, `agrees`, and `flightVerified`. Validate propagation completion, position and velocity matching, terminal mass transversality, and Hamiltonian residual. Use a full BVP residual when no reliable solver status is available.

- **[CORRECTNESS]** `DRO_tulip/indirect/verify_with_pumpkyn.m:59–68` - **The wrapper accepts an \(1\times8\) answer but subsequently concatenates it as though it were \(8\times1\).** `V.z8Pumpkyn = za(:)` normalizes the stored copy, not `za`; `[rv0(1:6); 1; za(1:7)]` can fail before `fenced` is called. **Fix:** assign `za = za(:)` immediately after validation and use that normalized vector throughout. Normalize endpoint states at the interface too.

- **[CORRECTNESS]** `DRO_tulip/indirect/verify_with_pumpkyn.m:55–60,98–110` - **Failure reasons and timeout guarantees are inaccurate.** Without a pool, the call is not capped at all. Any exception in the serial branch is swallowed and then described as exceeding the cap. **Fix:** return a structured status distinguishing timeout, exception, unusable output, and propagation failure; retain exception identifiers/messages. Either provide a real timeout mechanism or explicitly report that the serial call is uncapped.

- **[CORRECTNESS]** `DRO_tulip/indirect/transfer_study.m:145–173,193–211` - **The pass thresholds are mostly undocumented engineering choices, not defensible certification thresholds.** In particular:
  - `1e-8` in a position component is about **3.9 m** here.
  - N3 permits **100 km**, about \(2.6\times10^{-4}\) ND in position, and **10 m/s**, about \(9.8\times10^{-3}\) ND in velocity.
  - The mixed state/costate residual has no declared scaling.
  - Strict `> 0` checks provide no numerical safety margin.
  - `min|det|` is coordinate-dependent and not a conditioning certificate.
  - `|dz| < 1e-6` mixes final time and differently scaled costates.

  The loose single-shot arrival gate is particularly inconsistent with claiming a common accurate extremal for the second-order calculation. **Fix:** derive tolerances from physical endpoint requirements, propagation/interpolation error estimates, variable scaling, and sensitivity. Show refinement stability. Use positive lower margins for strict inequalities and uncertainty-aware rank/singularity decisions. A reasonable heuristic threshold can be useful, but call it that.

- **[OVERCLAIM]** `DRO_tulip/indirect/transfer_study.m:213–221` - **“Strict strong local minimizer … numerically certified at the sampled times” is too strong and mathematically ill-posed.** Strong local minimality is a property of the whole trajectory and its admissible neighborhood, not a property certified separately at sample times. “Same endpoints” also obscures that final mass is free. The failure branches overstate matters too: cross-check disagreement does not mean “NOT an extremal,” and failure of a sufficient hypothesis does not refute minimality.

  **Fix:** use a verdict such as:
  > Normal-PMP residuals satisfy the stated numerical tolerances. Selected regularity diagnostics passed; no candidate conjugate point was detected on the reported grid. Local optimality is not certified.

  Once the theorem mapping and numerical validation are genuinely established, state the precise conclusion: fixed initial state and terminal position/velocity, free terminal mass and time, specified admissible control set, and specified strong-local topology. No global optimality claim follows.

## Question 2 — Is it well written as a template?

**The numbered, visibly computed study format is worth keeping.** The main template problems are incomplete contracts, hidden numerical assumptions, repeated propagation, and a report whose certainty exceeds its calculations.

- **[STYLE]** `DRO_tulip/indirect/transfer_study.m:1–25,36–138` - **The overall narrative is good; retain it.** Orbit generation from parameters, explicit nondimensionalization, solving, implementation cross-checking, and mathematical diagnostics are understandable stages. The missing opening element is the precise optimization problem, and the missing intermediate element is numerical feasibility/acceptance. **Fix:** add a compact problem declaration before section 1 and a feasibility block immediately after solving. Keep the independent solver separate from the mathematical condition lists. This preserves the engineer’s style without hiding the machinery.

- **[STYLE]** `DRO_tulip/indirect/transfer_study.m:1–30` - **The main script violates the house header contract.** It has no INPUTS/OUTPUTS block, no array sizes, and no clear dependency or side-effect declaration. **Fix:** document edited parameter blocks as inputs; workspace results and the figure as outputs; `Y [N×14]`, `z8 [8×1]`, orbit tables `[N×6]`, and shooting starts `[14×K]`. List required MAT files, pumpkyn/pumpkynPie, CasADi, optional parallel facilities, and tested MATLAB/toolbox versions. Declare that the script clears the workspace, clears the console, and changes the path.

- **[READABILITY]** `DRO_tulip/indirect/transfer_study.m:41–90,130–134` - **Parameters are grouped, but the resulting objects are scattered.** Understanding the departure orbit requires tracking `dep`, `tD`, `rvD`, `infoD`, `ppD`, `stateD`, and `rv0`; the compatibility object `B` appears much later. **Fix:** retain parameter-driven construction but attach generated data to meaningful objects, for example `departure.params`, `.timeND`, `.stateND`, `.periodND`, `.stateAtPhase`, and `.endpointND`. Build the model/units object near the engine block. Use an adapter helper for the legacy `B` interface rather than making the teaching script revolve around that interface.

- **[CORRECTNESS]** `DRO_tulip/indirect/transfer_study.m:103–113` - **The advertised editable template has a hidden, incompletely checked seed operating point.** The guard checks departure period and petal count, but not branch, thrust, Isp, reference mass, mass ratio, or scale conventions. The library also discards the original time grid and the script invents a uniform one. An off-parameter trajectory can be a perfectly valid seed, but it is not a verified same-operating-point seed. **Fix:** store complete model metadata and the actual `tGrid` with each library entry. Either require matching metadata or explicitly label the entry as a continuation seed. Preserve the real endpoint instead of appending a duplicate start state unless the solver contract explicitly ignores that column.

- **[STYLE]** `DRO_tulip/indirect/transfer_study.m:150–173,233–234` - **Several short names impose unnecessary decoding, and two names are actively ambiguous.** `mm`, `nlv`, `z8`, and `it` force readers back to earlier assignments. Printed `lambda_0` is easily confused with the scalar PMP cost multiplier; `T.mfKg` actually stores propellant **used**, not final mass. **Fix:** prefer `massFrac`, `primerMag`, `shootingUnknowns`, `solveInfo`, and printed `lambda(0)`. Store both `finalMassKg` and `propellantUsedKg`. Define state/costate index groups once. The supplied code already respects the bans on bare `i`/`j` indices and `norm` shadowing; retain that discipline. `ii`/`jj` are not violations, though `idxDeparture`/`idxArrival` would teach better.

- **[STYLE]** `DRO_tulip/indirect/transfer_study.m:145–180,193–211` - **Repeated expressions and repeated thresholds make printed checks drift away from the aggregate verdict.** The omitted N5 Boolean is exactly the failure this style invites. **Fix:** calculate each metric and Boolean once, visibly inline, and delegate only presentation:
  ```matlab
  metric.maxAbsHamiltonian = max(abs(Hval));
  check.hamiltonian = all(isfinite(Hval)) && ...
      metric.maxAbsHamiltonian <= tol.hamiltonian;

  reportCheck("N2", "Hamiltonian", metric.maxAbsHamiltonian, ...
      tol.hamiltonian, check.hamiltonian);

  pmpSatisfied = all([check.boundary, check.dynamics, ...
      check.hamiltonian, check.controlMinimum, check.massTransversality]);
  ```
  Populate `tol` from a documented numerical policy. Keep engineering cross-checks out of this aggregate.

- **[READABILITY]** `DRO_tulip/indirect/transfer_study.m:201–228` - **“Every step is done HERE” is neither true nor a desirable literal template rule.** S3 delegates an adjoint integration and rank calculation; S4 consumes a precomputed conjugate verdict. Those algorithms should remain helpers. Inlining hundreds of lines of SVD, STM propagation, and root bookkeeping would obscure rather than teach the theorem. **Fix:** keep each condition’s mathematical definition, assumptions, extracted metric, tolerance, and acceptance decision inline. Put reusable numerical algorithms in documented, tested helpers returning evidence and status—not merely a Boolean. Explain this division in the header.

- **[READABILITY]** `DRO_tulip/indirect/transfer_study.m:139–228` - **The report is readable at first glance, but omits the information needed to interpret its PASS labels.** Thresholds, scaling, sample counts, coverage, and unresolved cases are largely absent. `cj.verdict` followed by another PASS/FAIL can also produce confusing output such as `ENDPOINT ... FAIL`. The inline/helper consistency assertion runs **after** the final verdict. **Fix:** print columns for ID, condition, value, tolerance/margin, units, and status. Use at least `PASS`, `FAIL`, `UNRESOLVED`, and `NOT EVALUATED`. Report conjugate coverage and minimum singular-value ratio; print normality residuals and rank gap. Run all consistency assertions before issuing a verdict. Required missing fields should be schema errors, not quietly replaced by `NaN` through `gv`.

- **[READABILITY]** `DRO_tulip/indirect/transfer_study.m:118,202,235` - **Repeated propagation makes the reader wonder which trajectory owns the reported numbers.** The script flies once, the hypothesis helper flies again, and the plot flies again; the shooting solution is a fourth representation. Repetition is not independence when it uses the same routine. **Fix:** create an explicit `flight` result containing time, augmented state, applied control if available, model metadata, and propagation settings. Let analysis and plotting accept that object. If a re-flight is intentionally a consistency check, label it as such and compare it numerically with the accepted flight.

- **[READABILITY]** `costate_common/ms_conjugate_test.m:10–51,155–215` - **The comments mix durable mathematics, implementation limitations, historical incidents, and review provenance.** The mathematics and limitations are valuable; dated references to individual reviews and measured catalog incidents make the algorithm harder to locate and the header longer to trust. Similar narrative clutter appears in `lift_space_dim` and `verify_with_pumpkyn`. **Fix:** keep equations, assumptions, failure modes, and exact contracts next to the code. Move review history and campaign anecdotes into version control, tests, or an audit document. A concise comment explaining *why* a step exists is more useful than recording who requested it.

- **[CORRECTNESS]** `costate_common/ms_tfmin.m:52–56` - **The self-resolving path is stale for this file’s supplied location.** Three `fileparts` calls on `.../costate_common/ms_tfmin.m` climb above the repository before appending `costate_common`. The script’s earlier `addpath` happens to mask this. **Fix:** remove helper-level path mutation and use one explicit project setup step, or resolve the actual sibling location correctly. Likewise, do not silently add `~/casadi-3.7.0` inside `fixedJacobian`; validate and report the configured dependency once.

- **[CORRECTNESS]** `DRO_tulip/indirect/verify_with_pumpkyn.m:43–45` - **The default-pool expression is evaluated even when the caller supplies `opts.pool`.** MATLAB evaluates function arguments eagerly, so `d('pool', gcp('nocreate'))` still invokes `gcp`, imposing a potentially unwanted Parallel Computing Toolbox dependency. **Fix:** branch explicitly:
  ```matlab
  pool = [];
  if isfield(opts, 'pool')
      pool = opts.pool;
  elseif exist('gcp', 'file') == 2
      pool = gcp('nocreate');
  end
  ```
  Document whether parallel execution is optional and how timeout behavior changes without it.

- **[STYLE]** `DRO_tulip/indirect/plot_transfer_3d.m:49–55` - **This loop is a good vectorization opportunity; the scalar library-call loops generally are not.** The script’s `ppval`-based phase functions already accept vector arguments. **Fix:** define that vectorized shape in the interface and write:
  ```matlab
  departureStates = B.stateD(ss);  % [6 x nSamples]
  arrivalStates   = B.stateA(ss);
  D = departureStates(1:3,:).';
  A = arrivalStates(1:3,:).';
  ```
  Keep ordinary `k` loops where the third-party API is scalar or where looping exposes the mathematics more clearly. Vectorization should improve the explanation, not merely reduce line count.

- **[OVERCLAIM]** `DRO_tulip/indirect/plot_transfer_3d.m:9–11,70,103–107` - **The plot’s claim to be the source of truth is false.** Geometry is repropagated, but time, delta-V, and propellant annotations are copied from `T`; they can disagree with that propagation. The plot also says “minimum-time transfer” regardless of the study verdict. **Fix:** obtain geometry and annotations from the same accepted `flight` object, or recompute all annotations from the re-flight and compare them. Label it “candidate minimum-time extremal” unless a stronger status has actually been established. For a reusable template, derive orbit names and legends from the orbit objects rather than hardcoding DRO/tulip.

- **[STYLE]** `DRO_tulip/indirect/dro_tulip_library.m:33–57` - **The library is a collection of file-specific assumptions rather than a validated reusable data contract.** It loads entire files, assumes nested field layouts, hardcodes phases and time scaling, and suppresses duplicates based only on phase. That is fragile template infrastructure. **Fix:** introduce a versioned entry schema containing model parameters, units, endpoint definitions, `z`, `Y`, `tGrid`, solver status, and provenance. Validate dimensions and finite data while loading; distinguish candidate seeds from independently audited results. Load only required variables where practical. The script should receive a well-defined seed object, not reconstruct one from historical MAT-file conventions.