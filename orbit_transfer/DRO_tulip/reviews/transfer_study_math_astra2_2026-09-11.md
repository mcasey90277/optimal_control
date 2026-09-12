## Executive assessment

The principal PMP corrections are mathematically sound. In particular, the full control-gap formula is correct **for a physically consistent control**, the mass-row throttle recovery has the correct normalization, the frozen-control adjoint comparison has the correct mass-costate sign, and **both kernel identities are valid for this 14-state all-burn system**.

However, three important gaps remain:

1. **The dense S4 gate does not require zero interior candidates**, contrary to your description. It can accept a true conjugate point classified as a “near-miss.”
2. **`certify_root` still enforces substantially less than `transfer_study`**: it does not enforce the new full-gap/throttle checks, X2, the strengthened S3 checks, or the dense conjugate scan.
3. The rank margin measures sensitivity to one numerical setting, **not the total error in the constraint matrix**.

None of these findings demonstrates that the reported anchor is nonoptimal. They demonstrate that the diagnostic machinery still has false-pass paths.

Line references below use the supplied filenames. `report_optimality.m`, `mintime_prop_seg.m`, and the underlying pumpkyn implementations were not supplied, so I cannot verify their internal changes.

---

## 1. Closure of findings 1–16 and 19–21

Here, **CLOSED** means the requested correction is present—not that a numerical diagnostic has become a mathematical certificate. Where the study and production certifier differ, I distinguish them explicitly.

- **[GAP]** `transfer_study.m:457–461`; `certify_root.m:269–280` — **1. CLOSED in the supplied callers.** Both now require `h6Ok` and a **strict** margin; `h6_margin.m:81–85` requires clearance above the supplied Hamiltonian residual. This is substantive, not cosmetic. The claimed change to `report_optimality` cannot be checked from this bundle.

- **[GAP]** `transfer_study.m:466–473` — **2. CLOSED as a wording correction.** V2 is now correctly described as a flight/exhaust-speed consistency check, not an independent thrust or mass-ratio check. “Shows up in N2” should preferably read “can show up in N2”: a scalar Hamiltonian projection is not a guaranteed detector of every field error.

- **[GAP]** `ms_conjugate_test.m:206–208,262–267,324–327`; `transfer_study.m:428–438` — **3. CLOSED.** Missing `Yend` prevents a free-time PASS, and `tFirstFullRank` is safely `NaN` when full rank is never attained. An already detected interior root can still produce FAIL despite incomplete coverage; that is logically appropriate.

- **[CORRECTNESS]** `pmp_pointwise_checks.m:117–125`; `transfer_study.m:342–345` — **4. CLOSED.** Absolute direction/full gaps and signed full-gap extrema are retained, and the study gates the absolute full gap. Negative full gaps are no longer clipped away.

- **[GAP]** `pmp_pointwise_checks.m:110–127`; `transfer_study.m:342`; `certify_root.m:182–187` — **5. PARTIALLY CLOSED.** The instrument and study now recover and gate both throttles and evaluate the direction-plus-throttle expression. But `certify_root` still gates only `PW.dirGap`, not `PW.fullGap` or either throttle error. Also, the reported “exact applied-field gap” remains conditional on agreement between acceleration-side and mass-side throttle; see §2.

- **[CORRECTNESS]** `ms_conjugate_test.m:283,291–298` — **6. CLOSED as a bracket-classification correction.** Opposite trustworthy nonzero signs, with the right endpoint at \(t_f\), imply an interior zero. The remaining issue is whether the numerical signs are trustworthy; `resolvedTol` alone does not establish that.

- **[GAP]** `ms_conjugate_test.m:223–231,340`; `transfer_study.m:435–438` — **7. CLOSED as requested diagnostics.** Both identities are measured and printed. Their mathematics is correct for the stated 14-state strict all-burn system. They are not acceptance gates, which is worth changing, but the adopted finding requested diagnostics.

- **[GAP]** `mintime_hypothesis_gates.m:111–122`; `transfer_study.m:366–373`; `certify_root.m:243–258` — **8. PARTIALLY CLOSED across the supplied pipeline.** X2 is substantive and correctly enforced by the study. The production certifier computes these fields but never gates them. Add `fieldErr` and `adjErrRef` to its required checks.

- **[GAP]** `transfer_study.m:410–418`; `certify_root.m:246–258` — **9. PARTIALLY CLOSED.** The study genuinely requires all four stated S3 checks. The production certifier still accepts `dimS == 1` without the lift residual, independent Hamiltonian residual, or `lift_margin`. The mathematical limitations of the error estimate also remain; they are not removed by requiring `LM.certified`.

- **[GAP]** `transfer_study.m:441–446`; `conj_spectrum.m:141–164,169–174` — **10. PARTIALLY CLOSED.** A dense scan is now enforced, but it does not close the sign test’s blind spots. The gate omits `nInteriorCand`, accepts “near-misses,” and ignores spectral candidates classified into entire start/end bands. A true conjugate point can pass; details below.

- **[ROBUSTNESS]** `validate_flight.m:73–92` — **11. CLOSED for the four requested additions.** Finite costates, sampled mass-law agreement, sampled positive mass, and strictly positive expected terminal mass are all checked. These remain sampled admissibility screens, not continuous trajectory validation.

- **[CORRECTNESS]** `transfer_study.m:199–208,543–555` — **12. PARTIALLY CLOSED.** The new propagation comparison is a real dynamics-based endpoint consistency test, unlike seam continuity. But it is an error estimate relative to a numerical reference trajectory, not a certified periodic-orbit error; moreover, multiplying a six-state position/velocity norm by `lStar` does not give the actual position error in kilometres. Report position and velocity discrepancies separately.

- **[GAP]** `transfer_study.m:281–285` — **13. CLOSED, wording only.** “PMP shooting equations” correctly replaces “first variation vanishes.” N6 supplies the control minimization condition.

- **[GAP]** `transfer_study.m:299–304,390–395` — **14. CLOSED, wording only.** The dependence of N2/N4 on the shooting equations and the implication of strict bang from the exact mass-costate equation, terminal transversality, and regularity are accurately explained.

- **[GAP]** `transfer_study.m:44–48,358–369` — **15. CLOSED in the study’s description.** X1 is now explicitly a second shooting implementation on shared physics; X2 is the independent field comparison. This does not turn a zero solver displacement into independent proof of convergence, but the shared-physics claim is corrected.

- **[GAP]** `transfer_study.m:491–504` — **16. CLOSED as the requested verdict rewrite.** Fixed endpoint position/velocity, initial mass one, free terminal mass/time, and fixed phases are explicit, and the text expressly disclaims certification. The theorem’s applicability still needs the qualifications discussed in §§3–4.

- **[GAP]** `transfer_study.m:435–445`; `conj_spectrum.m:98–104,141–144` — **19. PARTIALLY CLOSED.** The uncovered interval is honestly printed, and the dense scan adds samples inside the first segment. But its first sample is at \(t_f/192\), not arbitrarily near zero, and start-band spectral candidates are discarded rather than resolved. The short-time argument is still required.

- **[GAP]** `transfer_study.m:376–380,500–503` — **20. CLOSED as disclosure, not as a mathematical bound.** The absence of between-sample bounds is now explicit. No continuous positivity or root-exclusion proof has been added.

- **[GAP]** `transfer_study.m:96–98,307–311` — **21. CLOSED as wording.** N3 is correctly called a loose flown-arrival screen. N1 is a much tighter numerical boundary-value residual, but a small N1 residual still does not itself prove existence of an exact fixed-endpoint extremal.

---

## 2. Mathematical audit of the corrections and remaining defects

### 2.1 Full minimum-principle gap and recovered throttle

- **[GAP]** `pmp_pointwise_checks.m:110–122` — **The formula is correct, but its identification with the applied field’s Hamiltonian gap needs one additional consistency term.**

  For an admissible control \(b=u\alpha\), \(u=|b|\), the control-dependent Hamiltonian is
  \[
  H_c=\frac{T}{m}\lambda_v^\top b-\frac{T}{c}\lambda_m u,
  \qquad
  \min H_c=-T\max(Q,0).
  \]
  Therefore
  \[
  H_c-\min H_c
  =
  \frac{T}{m}\bigl(\lambda_v^\top b+|b|\rho\bigr)
  +T\bigl(\max(Q,0)-|b|Q\bigr),
  \quad \rho=|\lambda_v|.
  \]
  Your formula and its signs are correct.

  Likewise,
  \[
  u_{\rm mass}=-\frac{cF_m}{T}
  \]
  is correct for **mass fraction** with \(\dot m=-uT/c\); no additional mass factor belongs there.

  However, the actual recovered field has control contribution
  \[
  \frac{T}{m}\lambda_v^\top b-\frac{T}{c}\lambda_m u_{\rm mass}.
  \]
  Thus its gap is
  \[
  G_{\rm field}
  =G_{\rm code}
   +\frac{T}{c}\lambda_m\bigl(|b|-u_{\rm mass}\bigr).
  \]
  The new mass-throttle gate catches substantial inconsistencies, but it does not make the printed `fullGap` the exact field gap: the throttle tolerance is \(10^{-10}\), while the gap tolerance is \(10^{-12}\).

  **Concrete fix:** retain the decomposed control gap, but also compute the field gap using `uMass`, and gate acceleration/mass consistency explicitly. Label the two quantities accurately. Keep the signed extrema and absolute gates.

  One useful distinction: the direction part cannot be negative in exact arithmetic for \(T,m>0\). The negative gap associated with over-unit thrust arises in the **throttle part**, exactly as the corrected code now detects.

### 2.2 Coverage and the resolved final bracket

- **[GAP]** `ms_conjugate_test.m:262–265,283–295,328–331` — **The interior-root inference is correct; the numerical resolution policy is incomplete.**

  For a continuous determinant \(D\),
  \[
  D(t_a)D(t_f)<0,\qquad D(t_a)\ne0,\quad D(t_f)\ne0
  \]
  proves at least one root in \((t_a,t_f)\), not merely an endpoint possibility.

  But `sigR(end) > 1e-10` is a condition-number screen, not an error bound. Sign reliability requires the perturbation of the equilibrated matrix to be smaller than its smallest singular value. The test also does not require the other bracket endpoint to meet an analogous resolution criterion.

  More importantly, an **unresolved nonzero final determinant with the same sign as the preceding sample can still PASS**. `lastResolved` affects classification only when a sign change is found. The dense scan excludes \(t_f\), so it does not repair this endpoint gap.

  **Concrete fix:** make a numerically unresolved endpoint return UNDETERMINED regardless of whether a sign change is observed; resolve both bracket endpoints against a defensible matrix-error estimate. Apply similar handling to unresolved interior samples.

### 2.3 The right and left kernels with the mass costate present

- **[GAP]** `ms_conjugate_test.m:60–64,223–231` — **The implemented identities are correct for this system; the explanatory comment omits the mass-costate component.**

  Write the full canonical variables as \((q,m,p,\ell)\), with \(q=(r,v)\), \(p=(\lambda_r,\lambda_v)\), and \(\ell=\lambda_m\). On a strict all-burn arc,
  \[
  \dot m=-T/c,
  \]
  so fixed-time variations with respect to initial costates have \(\delta m=0\).

  Scaling \(p(0)\), while holding \(\ell(0)\) fixed, leaves the state and direction unchanged. Its propagated variation is
  \[
  W(t)=\bigl(0,0,p(t),\ell(t)-\ell(0)\bigr).
  \]
  Consequently,
  \[
  J(t)p(0)=0.
  \]
  The symplectic pairing of \(W\) with an initially vertical costate variation is initially zero and remains zero. At time \(t\), up to the chosen symplectic sign convention, that pairing is
  \[
  p(t)^\top\delta q+
  \bigl(\ell(t)-\ell(0)\bigr)\delta m
  =p(t)^\top\delta q.
  \]
  Hence
  \[
  p(t)^\top J(t)=0.
  \]

  **Your mass-costate argument is correct.** Holding \(\ell(0)\) fixed does not invalidate the left kernel.

  **Concrete fix:** amend the comment to show the extra mass-costate component and the use of \(\delta m=0\). Restrict this justification to strict all-burn variations; it does not automatically survive switches. Also make gross violations of these implementation identities block interpretation of S4 rather than merely print them.

### 2.4 Frozen-control adjoint comparison, including the mass row

- **[GAP]** `mintime_hypothesis_gates.m:35–39,119–121,189–194` — **The comparison is mathematically correct; the explanation should say “no derivatives of the optimizing control law,” not “does not depend on the control law.”**

  At the optimized control,
  \[
  H_x=f_x(x,u^*)^\top\lambda,
  \]
  where \(f_x\) is taken with control fixed. On this strict all-burn branch, \(u^*=1\) is locally constant and \(\alpha^*=-\lambda_v/\rho\) is independent of the state.

  In particular,
  \[
  \frac{\partial f_v}{\partial m}=-\frac{T}{m^2}\alpha,
  \]
  and therefore
  \[
  \dot\lambda_m
  =-\left(\frac{\partial f_v}{\partial m}\right)^\top\lambda_v
  =\frac{T}{m^2}\lambda_v^\top\alpha
  =-\frac{T\rho}{m^2}.
  \]
  The reference mass-costate row has the correct sign and magnitude.

  More generally, the envelope argument removes derivatives of the optimizing control from the derivative of the optimized Hamiltonian; it does **not** mean that \(H_x\) is independent of the **value** of the control.

  **Concrete fix:** correct that wording. No mathematical change to the reference adjoint calculation is needed. Also describe the reported error as a seven-row vector-norm comparison: it is not seven separately scaled relative row errors.

### 2.5 Endpoint propagation comparison

- **[CORRECTNESS]** `transfer_study.m:204–207,543–555`; `phase_state.m:63–64` — **The comparison is useful at the displayed phases, but its units, phase convention, and claim need correction.**

  For \(s\in[0,1)\), with the table beginning at zero and the same period used by both routines, propagating the first sample to \(sT_{\rm orbit}\) is the appropriate reference for the phase convention.

  What it measures is the discrepancy between:

  * the endpoint returned by the interpolant, and
  * a numerical CR3BP solution starting from the table’s first state.

  It does not independently establish that the first state is on an exact periodic orbit, that the period is exact, or that the reference propagation error is negligible. At \(s=0\), the comparison is identically zero and validates none of those things.

  There are also two concrete implementation defects:

  1. `e` is a six-state norm. `e*lStar` is **not the actual position error**. The reported 8.6 m is, at best, a position-error upper bound derived from that mixed norm.
  2. `phase_state` wraps with `mod(s,1)`, but the comparison helper does not. Negative phases are compared with the first sample; phases above one are propagated through additional periods.

  **Concrete fix:** use the same wrapped phase and verified period convention; report
  \[
  e_r=|\Delta r|\,l_*,\qquad
  e_v=|\Delta v|\,l_*/t_*,
  \]
  separately, validate propagation completion, and assess reference-propagation error. Call this an endpoint-consistency estimate unless the orbit and propagation are validated.

### 2.6 What Eckart–Young actually establishes

- **[CORRECTNESS]** `lift_margin.m:7–17,69–83` — **The header overstates the rank theorem.**

  For an \(m\times7\) matrix, \(\sigma_6\) is the distance to matrices of rank **at most five**. The distance to the nearest matrix lacking full column rank is \(\sigma_7\).

  If \(C_*\) is the exact constraint matrix and
  \[
  \|C_*-C_{\rm num}\|_2<\sigma_6(C_{\rm num}),
  \]
  then one obtains
  \[
  \operatorname{rank}C_*\ge6.
  \]
  One obtains rank exactly six only with the separate fact that an **exact nonzero lift** belongs to its kernel.

  A small residual for a supplied vector constructs a nearby matrix with a null vector; it does not construct an exact null vector of the original continuous problem.

  **Concrete fix:** state the lower-rank conclusion correctly, distinguish backward error from exact lift existence, and condition exact nullity one on an established exact normal extremal. The numerical calculation can remain a useful separation diagnostic.

- **[GAP]** `mintime_hypothesis_gates.m:93,125–134`; `transfer_study.m:105–110,410–412`; `lift_margin.m:60` — **The two-tolerance difference is an acceptable sensitivity estimate for one error component, not a total matrix-error estimate.**

  Changing `relTol` changes the frozen-control adjoint integration. It does **not** change the underlying `tfMinProp` flight or its pchip interpolants.

  Therefore the difference can miss common errors in:

  * the flown state and costate trajectory;
  * the interpolated direction and state;
  * endpoint and shooting data;
  * shared implementation or derivative errors;
  * accumulated errors in directions not exposed by the accepted lift residual;
  * matrix assembly and floating-point conditioning.

  Agreement of the lift residual at two tolerances shows that simply tightening that particular integration is not curing the residual. It does **not** prove that pchip is the source, nor bound its operator-norm effect on \(C\).

  Also, the supplied study actually compares **`1e-12` with `1e-9`**, at line 411—not `1e-10`. The separate observation about a floor at `1e-10` may be true, but it is not the setting used for the printed margin.

  **Concrete fix:** report the exact settings; refine the underlying flight/interpolation as well as the adjoint integration; use multiple resolutions or an independent construction to assess convergence. A proof needs a total error enclosure. A factor-of-ten policy applied to a two-run difference is not a validated safety factor.

- **[GAP]** `transfer_study.m:105,413`; `lift_margin.m:72`; `lift_space_dim.m:17–18` — **The \(10^{-4}\) lift threshold is defensible only as an empirical screen, and its meaning is inconsistent between callers.**

  On this anchor, \(2.2\times10^{-6}\) is about 45 times below the study threshold and very small compared with \(\sigma_6=0.678\). That supports using the threshold as a numerical rejection policy. It does not establish an exact lift or bound matrix error.

  The study checks
  \[
  \frac{\|C\lambda_0\|}{\|\lambda_0\|}<10^{-4},
  \]
  while `lift_margin` checks that same quantity against
  \[
  10^{-4}\sigma_1.
  \]
  These are not the same tolerance despite sharing one option value. The absolute stacked residual also depends on constraint weighting and sample count.

  Finally, the surviving claim in `lift_space_dim.m:17–18` that the cap makes a poor lift report zero nullity “never \(>1\)” is false.

  **Concrete fix:** choose and document one normalization for the lift residual, use it consistently, and remove the cap’s claimed guarantee. Calibrate acceptance against trajectory/interpolation refinement, not merely the observed floor.

### 2.7 A true conjugate point can be called “near-miss” and pass

- **[CORRECTNESS]** `transfer_study.m:442`; `conj_spectrum.m:155–174` — **The implemented dense gate is not the gate stated in your question.**

  The code requires:
  ```matlab
  CS.nInterior == 0 && CS.nZero == 0 && CS.multiplicity == 0
  ```
  Here `nInterior` means **coarse sign changes**, not interior candidates. It does **not** require:
  ```matlab
  CS.nInteriorCand == 0
  ```

  Furthermore, the plateau classifier is not a valid zero-exclusion test. Nested uniform grids need not produce a fixed reduction in the sampled minimum at every refinement level.

  A concrete failure mechanism is:

  * a genuine double conjugate point lies very close to a coarse-grid node;
  * locally two singular values behave like
    \[
    \sigma_5\sim a|t-t_c|,\qquad
    \sigma_6\sim b|t-t_c|,\qquad b/a<0.1;
    \]
  * the determinant has the same sign on both sides;
  * if the root is within \(h/32\) of that coarse node, the same node can remain nearest on the \(h\), \(h/4\), and \(h/16\) grids;
  * both refinement ratios can therefore be approximately one;
  * the sampled minimum can be below `tolCollapse` but above the `1e-8` override;
  * the candidate is called `near-miss`;
  * the multiplicity test also misses it because \(\sigma_6/\sigma_5<0.1\).

  All three terms in `s4Dense` can then pass despite the true rank loss.

  The singular-value ratio used for multiplicity measures relative sizes of two small singular values, **not whether both tend to zero**.

  **Concrete fix:** at minimum, make every unresolved interior candidate block PASS. Treat plateau as **UNRESOLVED**, not proof of a near-miss. To clear a candidate, establish a positive minimum against numerical error. Determine multiplicity from resolved rank loss, not a fixed ratio between the two smallest singular values.

- **[GAP]** `conj_spectrum.m:108–110,141–144,177–187` — **The scan has additional blind regions, and refinement discards useful sign information.**

  Entire candidate clusters touching the first segment are classified `start`; clusters reaching the last band are classified `endpoint`. Neither class is refined or counted in `nZero`. A genuine even-multiplicity conjugate point is not harmless merely because it lies in the first or last segment.

  The assertion that the determinant sign at \(t_f\) “carries no information” because of spectral grading is also wrong. Grading can make a sign numerically unreliable; it does not invalidate a resolved sign mathematically.

  Finally, `windowMin` returns only a singular-value minimum. Sign changes revealed by the finer integrations are discarded.

  **Concrete fix:** inspect candidates throughout \((0,t_f]\); exempt only an initial interval justified by a short-time argument. Carry determinant/sign information through refinement, include a resolved endpoint, and return UNRESOLVED when a candidate cannot be separated from zero. Replace “closes the blind spots” with “adds sensitivity to the blind spots.”

### 2.8 The production certifier still has the original enforcement problem

- **[GAP]** `certify_root.m:182–187,224–258,283–287` — **The stronger study verdict is not the verdict enforced for catalog certification.**

  The production path still omits:

  * `PW.fullGap`;
  * acceleration and mass throttle gates;
  * `g.fieldErr` and `g.adjErrRef`;
  * the accepted-lift residual;
  * the independent Hamiltonian residual as an S3 gate;
  * `lift_margin`;
  * `conj_spectrum`.

  Yet it returns `C.ok = true` and reason `"certified"`.

  This is not a new PMP error. It is an incompletely propagated correction—the same kind of enforcement gap found in the first review.

  **Concrete fix:** use a common decision routine for the study and production certifier, and add end-to-end mutation tests showing that failures of each newly required quantity block `C.ok`. Instrument-level mutation tests alone do not test caller enforcement.

### 2.9 The variational equations still lack an independent physics check

- **[GAP]** `mintime_rhs_point.m:22–23`; `mintime_hypothesis_gates.m:114–121`; `conj_spectrum.m:99–103` — **X2 validates the 14-vector field, not the STM generator used by S4.**

  Correct state and costate rows do not establish that the propagated 196 STM entries use the correct Jacobian. The frozen-control seven-state Jacobian is appropriate for the lift-space calculation, but it is not the full Jacobian of the optimized 14-state Hamiltonian field.

  In particular, the conjugate calculation must include
  \[
  \frac{\partial\alpha}{\partial\lambda_v}
  =-\frac{I-\alpha\alpha^\top}{\rho}.
  \]
  Both kernel identities can be small while other variational directions are wrong.

  **Concrete fix:** compare the pumpkyn STM generator against an independently differentiated optimized 14-state field, or independently integrate that variational system. Include full symplecticity and flow-derivative checks as supplementary diagnostics. This is an outstanding validation gap, not an assertion that the supplied STMs are actually wrong.

---

## 3. The verdict sentence

- **[GAP]** `transfer_study.m:491–504` — **The printed verdict is justified as conditional numerical evidence, not as an established application of the theorem.**

  The following statements are supported by the displayed run:

  * the diagnostics are consistent with a regular normal extremal;
  * no zero was detected by the particular scans;
  * the result is numerical evidence, not a certificate.

  The wording does not claim a global minimum or optimized phases, which is correct.

  The necessary qualification is that “if the stated hypotheses hold exactly” must include **the unresolved theorem-applicability and reduction arguments**, not merely replacing numerical tolerances by exact equalities. The current disclaimer emphasizes numerical limitations but does not explicitly name those structural theory gaps.

  **Concrete fix:** add:
  > Application of the theorem also remains conditional on the subarc-normality argument and on a valid free-mass, free-time second-variation reduction; these are not established by the diagnostics.

  I would also replace “the spectrum scan … refines candidates” with wording acknowledging that unresolved candidates can currently be classified away. For the displayed anchor there were no interior candidates, so that particular false-pass branch was not exercised.

---

## 4. Open theory items and what else is missing

The four recorded items are real and correctly motivated, but two need sharper formulations.

- **[GAP]** `mintime_hypothesis_gates.m:14–24`; `transfer_study.m:399–413` — **(a) Subarc normality: correctly open, but analytic continuation must address the boundary conditions as well as stationarity.** Establish that every relevant nontrivial subarc has the normality property required by the theorem. Analytic continuation can extend the frozen-adjoint stationarity identities along an exact analytic regular arc; it does not automatically transfer the full-arc condition \(\lambda_m(t_f)=0\) to a different subarc endpoint. The proof must account for that change and identify where H6 supplies the necessary nonvanishing normalization.

- **[GAP]** `ms_tfmin.m:77–82`; `ms_conjugate_test.m:29–32`; `transfer_study.m:420–425` — **(b) Six-state reduction: correctly open, and it must cover the competitor class, not just the reference dynamics.** Establish equivalence of the relevant second variations, endpoint conditions, critical directions, and conjugate tests after eliminating mass and its costate. Fixed-time all-burn costate variations have \(\delta m=0\), but free-time variations have a terminal mass contribution, and admissible competitors may reduce throttle. The argument must explain why strict bang allows those throttle directions to be handled without losing the claimed **strong**, rather than merely restricted or weak, local result.

- **[GAP]** `ms_conjugate_test.m:55–59`; `conj_spectrum.m:92–104` — **(c) Initial interval: correctly open.** Establish short-time full rank/no conjugacy and the appropriate initial sign or positivity on some \((0,\varepsilon]\), with a controlled remainder, then overlap that interval with the subsequently tested interval. Observing full rank at the first junction—or the first dense sample—does not establish what happened before it.

- **[GAP]** `transfer_study.m:376–380,500–503` — **(d) Between-sample bounds: correctly open, but include the conjugate matrix as well as scalar margins.** Establish continuous regularity, strict bang, clearance from singularities, and exclusion of determinant/rank loss on the whole required interval, including \(t_f\). On an exact all-burn arc, positive mass is especially simple to establish from the affine mass law; the other continuous statements still require bounds or analytical implications.

There are also additional steps between “every numerical line passes” and a theorem-backed conclusion:

- **[GAP]** `transfer_study.m:281–296`; `ms_tfmin.m:103–108` — **Existence of a nearby exact extremal is not established by a residual alone.** Establish an exact solution of the intended boundary-value problem near the numerical one, with enough control of state, costate, and final-time errors to transfer the strict margins and rank conclusions. A quantitative local existence argument or validated shooting solve is needed for certification; selecting the smaller residual from a more accurate integration mode is not a substitute.

- **[GAP]** `transfer_study.m:199–208`; `mintime_hypothesis_gates.m:125–147` — **Specify and control the exact boundary data and numerical objects to which the theorem applies.** If the endpoints are the literal supplied numerical vectors, say so. If they mean phases on exact DRO/tulip periodic orbits, establish the orbit/period/phase accuracy as well. Bound the resulting errors in the lift matrix and variational flow, not merely the residual along one supplied lift.

- **[GAP]** `transfer_study.m:376–378,497–498` — **Identify the precise theorem and verify its remaining structural assumptions for this formulation.** The audit must explicitly connect its regularity, normality, endpoint, control-set, and free-time hypotheses to this problem and establish the claimed strong-local topology. The matrix construction and H6 should appear in that connection, not be treated as self-validating because they produce sensible plots.

There is **no missing phase transversality condition** for the stated problem: the phases are fixed inputs. Nor should the review demand global optimality or exclusion of all competing transfer branches; the printed claim is local.

---

## 5. Per-condition verdicts

These grades assess the named numerical tests and their claimed role. **CORRECT does not mean interval-validated.** Changes refer to the first-review assessment summarized in the adjudication.

- **[GAP]** `pmp_pointwise_checks.m:75–81` — **N2: CORRECT.** \(H=1+\lambda^\top f=0\) is the correct autonomous free-final-time condition; the mathematics is unchanged, while X2 now supplies the missing independent field comparison in the study.

- **[GAP]** `pmp_pointwise_checks.m:83–84` — **N4: CORRECT.** Free terminal mass gives \(\lambda_m(t_f)=0\), provided the flight reaches the actual terminal time; unchanged mathematically.

- **[GAP]** `pmp_pointwise_checks.m:101–107`; `mintime_hypothesis_gates.m:119–121` — **N5: CORRECT.** The frozen-control state derivative gives the correct adjoint, including the negative mass-costate derivative; unchanged mathematically, with stronger independent checking.

- **[GAP]** `pmp_pointwise_checks.m:110–127`; `certify_root.m:187` — **N6: INCOMPLETE.** The admissible-control gap formula is correct and the study’s enforcement is substantially improved, but the reported exact field gap omits mass/acceleration inconsistency and the production certifier still ignores the new gates; the first-review enforcement gap is only partially closed.

- **[GAP]** `transfer_study.m:384–388` — **S1: CORRECT.** On a positive-mass all-burn arc the spherical Legendre form is \((T/m)|\lambda_v|I\), with continuous positivity still inferred only from samples; unchanged mathematically.

- **[GAP]** `transfer_study.m:390–397` — **S2: CORRECT.** \(Q>0\) is the correct strict-throttle condition and follows on an exact regular lift from the stated mass-costate equation and terminal condition; unchanged mathematically.

- **[GAP]** `transfer_study.m:410–418`; `lift_margin.m:60–94` — **S3: INCOMPLETE.** The study now enforces useful residual and separation diagnostics, but exact lift existence, total matrix error, and the subarc-normality connection remain open; enforcement improved without closing the theoretical assessment.

- **[GAP]** `transfer_study.m:420–446`; `conj_spectrum.m:141–174` — **S4: INCOMPLETE.** Coverage and final-bracket classification improved, but initial coverage, unresolved endpoint signs, ignored boundary-band candidates, plateau false negatives, and between-sample root exclusion remain; the original incompleteness persists.

- **[GAP]** `h6_margin.m:13–24,75–87`; `transfer_study.m:457–461` — **V1: CORRECT.** The reduced-Hamiltonian identity and strict clearance criterion are correct, and the previously missing enforcement is now present; mathematics unchanged, implementation finding closed.

- **[GAP]** `mintime_hypothesis_gates.m:111–122`; `transfer_study.m:370–373` — **X2: CORRECT.** It is a genuine independent seven-state/adjoint vector-field comparison with the correct frozen-control derivative, but does not validate the STM generator; this is a substantive new diagnostic addressing the first review.

**Bottom line:** the anchor still has persuasive numerical evidence of local optimality. The corrected study is appreciably more honest and stronger than the first version. It is not yet a certificate, and the most urgent remaining changes are to stop S4 from clearing unresolved spectral candidates and to make `certify_root` enforce the same checks as the study.