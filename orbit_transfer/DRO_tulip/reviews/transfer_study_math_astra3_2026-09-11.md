The widened gates are a substantial improvement, but the new candidate-resolution code still has false-clear paths. I am reviewing the supplied implementation and reported results—not reopening the four documented theory items.

## 1. `conj_spectrum`: candidate resolution

- **[CORRECTNESS]** `conj_spectrum.m:185–203, 250–274` -- **A bracket around the finest-grid argmin is not necessarily a unimodal bracket, and resolving one minimum does not resolve a cluster.** Golden section itself handles a V-shaped minimum perfectly well; differentiability is not required. The unsupported step is the assertion at line 191.

  At an isolated transverse, corank-one zero,
  \[
  \sigma_6(t)=|u^\top M'(t_*)v|\,|t-t_*|+O((t-t_*)^2),
  \]
  so a sufficiently small neighborhood is suitable. At a transverse corank-two zero, the two small singular values behave to first order like the singular values of
  \[
  (t-t_*)\,U_0^\top M'(t_*)V_0.
  \]
  If that projected derivative is nonsingular, both are locally V-shaped: a simultaneous collapse does **not**, by itself, invalidate golden section. But neither transversality nor a sufficiently small neighborhood is established here. Singular-value branch exchanges, closely spaced roots, higher-order contact, or numerical splitting of a multiple zero can make the bracket non-unimodal. A positive near-miss is smooth only locally when its smallest singular value is simple.

  **Concrete fix:** subdivide each cluster into search brackets for **all** detected local minima and sign-change brackets, not just the single finest-grid winner. Preserve ambiguous/multimodal cases as UNRESOLVED. Do not describe a local minimizer as resolution of the entire cluster without establishing the needed search assumption.

  The shifted grids defeat the **literal reuse of the same coarse node** in the old \(h/32\) example. They do not defeat sampling aliasing generally: fixed shifts can still sample a narrow zero less favorably than another positive well. For a single accurately evaluated, isolated V-shaped minimum, the nearest fine node and its neighboring fine steps do bracket the root; the minimization, rather than the shifts alone, supplies the improvement.

- **[CORRECTNESS]** `conj_spectrum.m:140–149, 175–203, 229–247` -- **The resolution discards evidence already computed, including refined sign changes.** `windowMin` obtains a spectrum but throws away the determinant. `locateMin` also throws away the determinant. Consequently, `cd.signChange` means “a coarse-inner-grid sign transition,” not “any sign change in the refinement window.”

  Furthermore, classification uses only the final `locateMin` result. A coarse or refined sample already below `zeroFloor` can be replaced by a larger located value and classified as a near-miss. This can happen when the minimizer chooses another well, or when segmented and single-window propagations disagree numerically. The shifted grids deliberately omit coarse nodes, making retention of previous evidence especially important.

  **Concrete fix:** carry a persistent record of:
  1. the smallest value and its spectrum across **every** evaluation;
  2. every numerically trustworthy opposite-sign bracket;
  3. disagreement between repeated evaluations of the same time using different propagation paths.

  Never upgrade an already observed floor-level candidate to CLEARED merely because a subsequent search returns a larger value. Disagreement beyond the measured numerical error should make the candidate UNRESOLVED.

- **[CORRECTNESS]** `conj_spectrum.m:120–149, 171–183, 235–247, 263–274` -- **The endpoint spectrum is stored, but it does not initiate a candidate, and the endpoint is not actually evaluated by either refinement grid or golden section.** Candidate construction uses only `1:N-1`. An endpoint window exists only if an inner cluster extends sufficiently far. Both shifted grids stop half a fine step short of the window endpoint. Golden section includes the endpoint in its interval but never evaluates it.

  Thus “the spectrum at \(t_f\) IS scanned” overstates the implementation. The separate junction test may block some such cases in `certify_root`, but that does not make `conj_spectrum.clear` an internally sound summary.

  **Concrete fix:** include the already computed \(t_f\) spectrum in candidate detection; explicitly evaluate both search endpoints; retain the \(t_f\) determinant with the same sign-resolution policy as other samples. A trustworthy opposite-sign bracket ending at \(t_f\) places a root strictly before \(t_f\); an unresolved or floor-level endpoint must not be discarded. Remove the incorrect numerical premise at lines 122–123 that its sign “carries no information.”

- **[GAP]** `conj_spectrum.m:157, 194–203, 266–274, 277–285` -- **Neither the zero floor nor the stopping tolerance is tied to the error in the matrix whose singular values are classified.** The relevant error is not a relative error in the full \(14\times14\) STM, nor the symplectic defect, nor the difference between two sampled minima.

  Let \(A(t)\) denote the matrix **after the chosen scaling**, and suppose its numerical error is estimated or bounded by \(\epsilon_A(t)\). Weyl’s inequality gives
  \[
  |\widehat{\sigma}_j(t)-\sigma_j(t)|\leq \epsilon_A(t)
  \]
  before accounting for the SVD’s own rounding error. A sensible relative floor therefore has the form
  \[
  \texttt{zeroFloor}(t)
  \sim
  \frac{s\,[\epsilon_A(t)+\epsilon_{\rm SVD}(t)
                   +L_A\,\delta t]}{\texttt{med}},
  \]
  where \(s\) is explicitly a safety factor and the last term accounts for minimum-location uncertainty when interpreting a located rank loss.

  **Concrete fix:** expose integration settings and compare the actual projected matrices at common times, especially candidate minima and \(t_f\), using tighter integrations and preferably another integration method. Start those comparisons from \(t=0\): both current refinements inherit the same stored prefix error, which their comparison cannot reveal. Measure errors in the STM block, the propagated state/flow column, matrix assembly, and scaling. Use a common quotient basis and controlled/common scalings for comparisons. Comparing singular values alone can conceal matrix errors.

  Return the final minimizer bracket and its convergence status. Line 270 implements an absolute nondimensional tolerance when `tB < 1`, **not** a tolerance relative to the initial bracket; translate its width into spectral uncertainty rather than assuming that \(10^{-9}\) in time is adequate.

  The anchor’s \(8.0\times10^{-5}\) is about **800 policy floors but only eight clearance thresholds**. It is encouraging, but it does not calibrate the floor. The appropriate measurement is whether the uncertainty in the relevant scaled matrix is comfortably below \(8.0\times10^{-5}\,\texttt{med}\).

- **[GAP]** `conj_spectrum.m:196–203, 219–226; certify_root.m:358–362` -- **“At the floor” establishes numerical rank ambiguity, not an exact zero or a mathematical refutation.** A sufficiently small positive near-miss also satisfies the ZERO rule. Conversely, an error perturbation can lift a true multiple zero into a positive computed minimum; if that error exceeds the assumed floor by enough, this code clears it.

  Also, `nSmall` measures **numerical corank at one approximate time**, not determinant-root multiplicity. A quadratic tangency can have determinant vanishing order two but only one vanishing singular value. At a corank-two zero, location error can leave \(\sigma_5\) above the floor while \(\sigma_6\) is below it.

  **Concrete fix:** distinguish “floor-level/possible rank loss” from a root established by trustworthy opposite signs. Both should block clearance, but only the latter supports the corresponding sign-crossing finding without further work. Call `nSmall` numerical corank, include location/error uncertainty when interpreting it, and do not equate the number of corank-two clusters with total conjugate multiplicity.

- **[CORRECTNESS]** `conj_spectrum.m:166–174, 204–226` -- **Contiguity does not establish that an entire cluster is structural startup.** A genuine even-corank zero, a tangency, or an undetected pair in the first segment can belong to the same low-spectrum cluster as the initial transient. It is then not refined and contributes nothing against `.clear`. A coarse sign change would still block via `nInterior`; the dangerous cases are precisely those without such a sign change.

  There is an especially simple diagnostic counterexample: if the block stays rank deficient throughout, every sample can belong to one `start` cluster, all determinants can be zero, and `.clear` becomes true. There is no “ever became testable” requirement.

  **Concrete fix:** require a finite, resolved departure from the structural startup region before reporting scan clearance. Refine suspicious structure away from \(t=0\) even when it remains contiguous with the startup cluster; otherwise retain it as UNRESOLVED. Report the actual uncovered time interval—`nStart = 1` does not report its extent. At minimum, separate “clear on the tested portion” from coverage status, and never clear an entirely rank-deficient scan.

- **[CORRECTNESS]** `conj_spectrum.m:277–285` -- **Column normalization preserves rank at each fixed time, but it need not preserve the limit \(\sigma_6(t)\to0\).** For finite positive column scalings, rank and determinant sign are unchanged. If every column norm remains bounded away from zero near a conjugate time, the interpretation is sound.

  A vanishing column is different. For example,
  \[
  M(t)=\operatorname{diag}\big((t-t_*)^2,1,1,1,1,1\big)
  \]
  loses rank at \(t_*\), but your normalized matrix is the identity at every other time. There is no spectral dip and no determinant sign change. The same issue can occur when a null direction happens to align with a chosen quotient column; it is therefore also basis-sensitive.

  **Concrete fix:** retain and inspect raw column norms against their **absolute numerical uncertainty**. A column approaching its uncertainty floor must itself generate an unresolved candidate. Freeze well-conditioned reference column scalings over each search window, or monitor an additional fixed-scaled spectrum. `realmin` prevents division by zero; it provides no numerical trust criterion.

- **[ROBUSTNESS]** `conj_spectrum.m:112–117, 237–245, 256–260; mintime_prop_seg.m:27–30` -- **The new scan can label an incomplete propagation as a completed time step.** `mintime_prop_seg` discards returned times and takes the final available row. The scanner then advances its nominal time unconditionally. An ODE call that returns early with finite output is not caught by a wall-time fence.

  **Concrete fix:** have the segment helper verify arrival at the requested duration, expected output dimensions, and finiteness of both state and STM. Propagation failure must abort the scan as unestablished, not supply a spectrum at the wrong time. Validate positive finite durations and integer sampling/refinement options before starting.

- **[GAP]** `tests/test_conj_spectrum.m:35–52, 71–86` -- **The existing controls do not exercise the newly claimed resolution mechanism.** The refuted 22.05-day fixture already has a coarse sign change; its ZERO classification is forced independently of golden section, the refinement ratios, and even-multiplicity detection.

  **Concrete fix:** extract the resolution logic so it can accept a synthetic matrix-valued function. Add mandatory tests for: a corank-two zero without a sign change; a quadratic touch; roots swept through all coarse/fine grid phases; multiple wells in one cluster; an endpoint-only candidate; a root merged with startup; a vanishing column; and controlled matrix noise around the floor. Assert `.clear`, not merely the existence of a ZERO candidate. The positive control should not silently disappear when an optional disk fixture is unavailable.

## 2. `certify_root`: enforcement, malformed data, ordering and caps

For **well-formed finite nonnegative scalar** inputs, `max(PW.fullGap, PW.fieldGap)` is the correct gated quantity: both fields are maxima of absolute gaps, so this enforces both requirements. Keep the two constituent values as provenance rather than retaining only their maximum.

The second gates build and spectrum scan are also **correctly fail-closed on timeout**: lines 325 and 370 return while `C.ok` is still false. Their exceptions are caught, and their failed outputs are not subsequently used. I found no new quantity being numerically gated before it is computed.

- **[CORRECTNESS]** `certify_root.m:218–228` -- **Gate 2b does not use the scalar validation applied to gate 5.** In MATLAB, `if ~(x <= tol)` is not a schema check:
  - `x = []` skips the refusal;
  - `x = [0 Inf]` also skips it, because the negated comparison is not true in every element;
  - scalar `-Inf` passes an upper-bound check;
  - scalar NaN and positive Inf are rejected by these particular comparisons, but that does not protect the aggregate inputs.

  `max(PW.fullGap, PW.fieldGap)` can also hide a NaN under MATLAB’s default missing-value behavior.

  **Concrete fix:** validate each required PW field as a real, finite, **nonnegative scalar before any aggregation or comparison**. Validate the individual full gap, field gap, acceleration throttle error, and mass throttle error; then form their maxima. Apply the same nonnegative-domain check to residual fields in gate 5—`scalar_verdict` currently validates representation, not the mathematical domain.

- **[CORRECTNESS]** `pmp_pointwise_checks.m:80, 101–138; mintime_hypothesis_gates.m:130–138` -- **NaNs and empty sampling can become perfectly finite zero residuals upstream, defeating even the new scalar guards.** The pattern
  ```matlab
  err = max(err, newErr);
  ```
  with an initial zero can discard a NaN. A NaN adjoint evaluation can therefore leave `adjErr` or `adjErrRef` equal to zero. Similarly, `pwOpts.nSample = 0` produces no N5/N6 evaluations and leaves their errors zero; `hRel = 0` can make the finite differences NaN and then have that failure hidden by the maximum.

  **Concrete fix:** reject nonfinite field and perturbed-field evaluations immediately; require positive finite FD steps and a positive integer sample count; require actual evaluation counts before allowing each gate to pass. Never initialize an “unperformed check” to a value that means success. Use explicit finite checks or NaN-propagating reductions, rather than relying on downstream scalar validation to reconstruct discarded failures.

- **[ROBUSTNESS]** `certify_root.m:326–332, 371–379; lift_margin.m:65–85` -- **Some new gate inputs remain unchecked or can terminate the sweep with an unnamed exception.** The two constraint matrices are checked only for equal shape, not for real finite numeric contents, seven columns, or sufficient rows. `lift_margin` calls `svd` before such validation. `LM.certified` is not read through a Boolean validator. The dense scan’s counts can be malformed or inconsistent with its scalar `.clear` flag without being rejected.

  **Concrete fix:** validate matrix schemas before SVD; validate `LM.certified` as a scalar Boolean; and catch/name failures from `lift_margin`. Validate dense counts as finite nonnegative integers and check consistency with `.clear`, together with an explicit completed/tested flag. Validate tolerance and mode options too: a NaN `marginMin`, for example, defeats `if M.margin < marginMin` in `lift_margin`.

- **[GAP]** `certify_root.m:363–384; tests/test_certify_enforcement.m:85–88` -- **The supposedly mandatory dense gate remains bypassable while returning “certified.”** `conjSpectrum = false` returns `C.ok = true` and the same certification wording. An empty diagnostic field is not sufficient protection for callers that consume only `.ok`; the test explicitly requires this bypass.

  **Concrete fix:** either forbid disabling the scan in certification mode, or return a distinct incomplete/diagnostic result with the full-study certification flag false. Change the test accordingly. A weaker certification mode must have a different explicit contract and must not be eligible for production writeback as a full pass.

- **[ROBUSTNESS]** `certify_root.m:290–318; tests/test_certify_enforcement.m:48–84` -- **The “first failed gate” promise and the mutation coverage are not as strong as stated.** Gate 5 validates all five summary fields before testing the X2 bounds. A malformed later `Hresid` can therefore mask an already out-of-tolerance `fieldErr`. More importantly, the tests do not independently establish enforcement of `adjErrRef`, the independent `Hresid`, both throttle constituents, malformed outputs, or second-build/scan timeout paths. Wrong-sign and half-mass-flow injections can be rejected by earlier Hamiltonian checks, so their rejection does not prove the later gate is wired.

  **Concrete fix:** validate and gate quantities in a documented order, or explicitly give malformed-result validation precedence over numerical gate ordering. Add a cheap gate-stack test seam with controlled returned summaries so each individual check can be made the first failing check. Include empty/NaN/Inf/vector mutations and named second-build/scan failures. The existing end-to-end physical mutations should remain as complementary tests.

  Separately, a conjugate verdict is currently refused before its hypothesis gates are evaluated. That is fail-closed for `.ok`, but a failure reason alone must not be interpreted as a physical conjugate refutation when the instrument’s prerequisites were never established.

- **[ROBUSTNESS]** `certify_root.m:145–147, 393–404; run_capped.m:44–50` -- **These are per-call fences, not a total certification wall-time limit, and `cancel(fut)` is not a general guarantee of killing and restarting a worker.** Without a pool, the named caps are entirely inactive. Native/uninterruptible work also requires more care than the header’s unconditional worker-kill claim suggests.

  **Concrete fix:** require a suitable pool for production sweeps, distinguish timeout from worker error in the returned status, and test the actual cancellation/recovery behavior of the deployed propagator. If a hard process limit is required, use process-level isolation/recovery. Add an overall budget if `wallSec` is intended to bound the whole certification. These limitations do **not** create a pass after a correctly reported second-build or spectrum timeout.

## 3. `ms_conjugate_test`: determinant-sign trust and kernel identities

- **[GAP]** `ms_conjugate_test.m:204–205, 250–263, 327–342` -- **\(10^{-10}\) is plausible protection against six-by-six LU roundoff, but it is not a derived trust threshold for the propagated matrix.** For pivoted LU of the equilibrated matrix \(A\), a standard backward-error model is
  \[
  |\Delta A_{\rm LU}|\lesssim \gamma_n\,|L|\,|U|.
  \]
  Thus the relevant relative LU error depends on the factorization growth, not just machine epsilon. At \(n=6\), \(10^{-10}\) will normally be generous for LU rounding alone. It says nothing about accumulated STM error, flow-column error, or errors amplified by equilibration.

  **Concrete fix:** measure the LU residual/growth and estimate the error in the **same equilibrated matrix**, including propagation and assembly. A defensible sign condition is
  \[
  \widehat{\sigma}_{\min}
  >
  \epsilon_{\rm matrix}+\epsilon_{\rm LU}+\epsilon_{\rm SVD},
  \]
  with the stated safety margin. This ensures that the perturbation cannot cross the singular set and change determinant sign. Translate that into a per-sample ratio threshold if the interface requires one. Keep determinant sign separately from its reconstructed magnitude; sign trust should not depend on successful exponentiation of a log-magnitude.

- **[CORRECTNESS]** `ms_conjugate_test.m:293–332` -- **The UNRESOLVED rule is applied after counting roots using the very signs declared untrustworthy.** For example, trustworthy positive samples surrounding one unresolved negative sample can generate interior crossings and force FAIL at line 330. That is not “an interior root was already found”; it is an interior root inferred from an untrusted sign. The final-bracket logic similarly tests resolution of the final sample but not necessarily of the sign against which it is compared.

  **Concrete fix:** classify samples as trusted-positive, trusted-negative, or unresolved **before** root counting. Opposite trusted signs, including across a run of unresolved samples, establish a crossing bracket. Equal trusted signs surrounding uncertain signs do not. A floating-point LU zero or a tolerance-defined zero alone should not become a proven interior touch. Otherwise return UNDETERMINED and refine. This preserves the intended priority of a genuinely established interior crossing over unrelated unresolved samples.

- **[GAP]** `ms_conjugate_test.m:233–241, 348–357; certify_root.m:266–273` -- **The quotient identities should be validity gates, not just printed diagnostics.** A large kernel residual means that the matrix being quotiented is inconsistent with the assumed invariance. Continuing to trust its determinant is then unjustified.

  **Concrete fix:** gate both identities before accepting a conjugate verdict. Calibrate tolerances against the propagated block/costate errors:
  \[
  r_R \lesssim \frac{\|\delta J\|}{\|J\|}
                  +\frac{\|\delta p_0\|}{\|p_0\|}
                  +O(u),
  \qquad
  r_L \lesssim \frac{\|\delta J\|}{\|J\|}
                  +\frac{\|\delta p_t\|}{\|p_t\|}
                  +O(u).
  \]
  A provisional \(10^{-10}\)–\(10^{-9}\) threshold has substantial headroom over the reported anchor residuals, but it must be checked against entry-specific numerical accuracy rather than treated as universal.

  Return explicit availability flags: skipped identities currently remain zero and can look like perfect checks. A failed or unavailable required identity should make the instrument UNDETERMINED, not label the trajectory physically refuted. These identities still do not discriminate all incorrect variational generators.

## 4. `test_stm_variational`: the symplecticity lesson and a cheap discriminator

- **[GAP]** `tests/test_stm_variational.m:60–81, 102–105` -- **The revised explanation is correct; the important qualification is that the frozen control is a prescribed time-dependent schedule.** At each time, the frozen generator is
  \[
  A_f(t)=J\,H_{zz}(z(t),\alpha(t)),
  \]
  with a symmetric Hessian, so \(A_f^\top J+JA_f=0\). Its fundamental matrix is symplectic even though it is not the derivative of the optimized flow. The control need not remain constant over the segment. There is therefore nothing surprising about the frozen generator having a symplectic defect as small as the true one.

  **Concrete fix / cheap discriminator:** directly test the generator block
  \[
  A(4\!:\!6,11\!:\!13)
  =-\frac{T}{m\rho}\left(I-\alpha\alpha^\top\right),
  \qquad \rho=\|\lambda_v\|.
  \]
  On this strict all-burn branch, this is precisely the block missing from the frozen-control generator. Its two tangential eigenvalues are \(-T/(m\rho)\), its radial eigenvalue is zero, and its trace is \(-2T/(m\rho)\). A tangential quadratic-form check or the full three-by-three identity is cheap and decisive; the radial-kernel check alone is not, because a zero block also passes it.

  Evaluate the actual variational RHS with an identity STM to expose the generator being integrated, rather than testing an independently constructed substitute. The reported FD agreement versus the \(1.4\times10^{-2}\) frozen discrepancy already provides strong evidence against this particular omitted-derivative mutation.

- **[GAP]** `tests/test_stm_variational.m:40–58` -- **One central-difference step on one anchor segment does not measure the STM accuracy needed by the spectrum floor.** The agreement is useful, but FD truncation and propagation error have not been separated; full-column normalization also does not directly assess errors in the small projected singular direction at a late candidate.

  **Concrete fix:** repeat with several step sizes and tighter propagation, observing the expected \(O(h^2)\) regime and eventual noise growth. Add generator-block checks along the trajectory and representative difficult entries. For floor calibration, compare full accumulated/projected matrices at candidate times, not just these first-segment columns. Retain the present test as a regression discriminator, not an accuracy certificate for the dense scan.

## 5. S3 lift-residual normalization

- **[GAP]** `mintime_hypothesis_gates.m:152–167; lift_margin.m:82–96` -- **The new quantity is a valid normwise backward error, but it is not sufficient evidence that these heterogeneous constraints are individually satisfied.** Indeed,
  \[
  \min_{(C+\Delta C)\lambda=0}
  \frac{\|\Delta C\|_2}{\|C\|_2}
  =
  \frac{\|C\lambda\|_2}{\|C\|_2\|\lambda\|_2}.
  \]
  That is mathematically correct for an **unstructured, global spectral-norm perturbation**. It allows errors in small rows to be judged against the largest rows anywhere in the stack. Repeating alignment rows also changes their weight relative to the single terminal-mass row.

  With \(\sigma_1=4.1\times10^3\), the \(10^{-6}\) gate allows
  \[
  \frac{\|C\lambda\|}{\|\lambda\|}\leq 4.1\times10^{-3}.
  \]
  That is not automatically a meaningful tolerance for either a local direction constraint or terminal transversality. The actual anchor value \(2.2\times10^{-6}\) is much better; the issue is what the acceptance rule permits.

  **Concrete fix:** retain the global backward error, but supplement it with:
  - maximum blockwise or componentwise scaled constraint residuals;
  - the alignment residual
    \[
    \frac{\|[\alpha_k]_\times\Psi_v(t_k)\lambda\|}
         {\|\Psi_v(t_k)\lambda\|}
    \]
    wherever the denominator is resolved and nonzero;
  - a separate terminal check on \(|e_m^\top\Psi(t_f)\lambda|\), using the normal-chart transversality tolerance/error scale.

  For the rank calculation, choose defensible block/row weights, including sampling weights if the stack is intended to represent an integral norm, and compare **the same weighted matrices** in `lift_margin`. Do not retain the old 2545-fold margin without recomputing it after changing weights.

  Also report the lift’s residual relative to \(\sigma_6\):
  \[
  \sin\angle(\lambda,v_7)
  \leq \frac{\|C\lambda\|}{\sigma_6\|\lambda\|}.
  \]
  On the anchor this is approximately \(2.2\times10^{-6}/0.678=3.2\times10^{-6}\), useful evidence that the supplied lift is close to the computed weakest right-singular direction. The \(\sigma_1\)-normalized number alone does not communicate that information.

**Bottom line:** I would not yet use this version to re-sweep the 115 entries and write production CLEARED/REFUTED classifications. The widened stack and its timeout returns are substantially better, but candidate resolution can still discard roots through single-well minimization, lost refinement/endpoint evidence, startup absorption, and time-dependent column normalization; malformed-data paths also remain. Fix those mechanisms, preserve ambiguity conservatively, calibrate the scaled-matrix error, and add adversarial synthetic resolution tests before the production sweep. The anchor’s located positive minimum is encouraging and the 22.05-day sign crossing remains useful evidence, but neither validates the new even-root clearance mechanism. The documented theory work can remain open as requested; these are implementation and numerical-classification changes that should precede the re-sweep.