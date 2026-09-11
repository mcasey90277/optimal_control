## A — Mathematics and consequences of the repairs

### A1 — N6

- **[REDUNDANT]** `transfer_study.m:250–262` - Your doubt is correct: N6 is another tautology, not merely a weak test. You construct the analytical minimiser rather than obtain the control actually applied by the propagator. Every sampled competitor must lose, and initialising `worst = 0` normally makes the reported result exactly zero. For unit directions,
  \[
  H(\alpha)-H(\alpha_*)=
  \frac{T}{m}\bigl(\lambda_v\!\cdot\!\alpha+\|\lambda_v\|\bigr)
  =\frac{\kappa}{2}\|\alpha-\alpha_*\|^2,\qquad
  \kappa=\frac{T}{m}\|\lambda_v\|.
  \]
  Here \(\kappa\) is a **curvature margin**, not a positive gap against every other direction. **Fix:** replace sampling with this exact optimality gap evaluated using the **applied** control—logged during propagation, or recovered from the implemented acceleration after subtracting an independently evaluated gravitational/Coriolis field—and check its admissibility separately. If you are not auditing the applied control, label directional minimisation **ANALYTIC**, remove the numerical gate, and promote S2’s throttle/bang minimisation check into N6. Keep S1: nonzero \(\lambda_v\) supplies regularity and strict curvature; it does not independently verify the implemented direction’s sign.

### A2 — N5

- **[FINE]** `transfer_study.m:226–246` - Yes: differentiating \(H\) in the state at fixed costate tests \(\dot\lambda=-H_x\); on the regular all-burn branch, fixed costate also fixes this direction, and your centred-difference formula is correct.

- **[CORRECTNESS]** `transfer_study.m:233–248` - The formula is sound, but one absolute step does not distinguish adjoint error from finite-difference error at the declared \(10^{-7}\) threshold. For reused cases, `m > 0` also does not ensure `m-hFD > 0`, and a perturbation can leave the current bang branch. **Fix:** explicitly freeze the applied control, use domain-safe component steps, and compare results at several steps, such as \(h/2,h,2h\). Report the differentiation uncertainty alongside the residual, or use an independently constructed analytic/AD state Jacobian. This remains an RHS-identity audit, not an independent test that the state dynamics model itself is correct.

### A3 — H6

- **[FINE]** `h6_margin.m:13–21` - For an exact normal all-burn extremal with \(T,c,m>0\), S1 and terminal transversality, \(\dot\lambda_m=-T\|\lambda_v\|/m^2<0\); therefore \(\lambda_m(0)\) controls the whole arc, and your strict inequality excludes zeros on **\([0,t_f]\)** without an all-arc search.

- **[FINE]** `h6_margin.m:21–24` - Two qualifications: on **\((0,t_f]\)** the equality case places the only zero at the excluded initial endpoint, so strict `<` is conservative there; the integral identity is correct because \((m\lambda_m)'=-TQ_{\rm mt}\), with \(m(0)=1\) and \(\lambda_m(t_f)=0\).

- **[CORRECTNESS]** `h6_margin.m:55–79` - The implementation needs a numerical clearance and validated options. With nonzero Hamiltonian residual,
  \[
  h(t)=H(t)-1+(T/c)\lambda_m(t),
  \]
  so an arbitrarily small positive \(\delta=1-(T/c)\lambda_m(0)\) is not enough to establish numerical separation from zero. Also, `marginMin < 1` can permit a genuine H6 violation, while failing a deliberately enlarged margin does **not** imply that \(h\) can vanish. **Fix:** require finite positive \(T,c\), require `marginMin >= 1`, compare \(\delta\) against the Hamiltonian/costate uncertainty allowance, and distinguish **mechanism excluded but requested headroom insufficient** from **mechanism not excluded**.

- **[CORRECTNESS]** `transfer_study.m:322–329`; `certify_root.m:220–234`; `report_optimality.m:86–109` - H6 is not consistently enforced. The study passes a missing `h6Ok`; `certify_root` ignores it entirely; the reusable report ignores it too. Consequently, a library entry can still be called certified with `g.h6Ok == false`. **Fix:** require a valid affirmative H6 result wherever the reduced conjugate instrument supports acceptance. Missing H6 should mean **UNRESOLVED/NOT CHECKED**, not pass. Keep its classification as an instrument-validity prerequisite, not a new theorem hypothesis.

### A4 — Repairs that remain incomplete or introduce conflicting behaviour

- **[FINE]** `transfer_study.m:264–272,331–338` - Excluding independent verification from `necessary`, and moving the consistency assertion before the verdict, are both correct; keep them.

- **[CORRECTNESS]** `report_optimality.m:71–84,102–118` - The separation repair did not reach the reusable report: independent agreement and witness flight still belong to `necessary`, so their failure still produces “NOT an extremal.” The report also has no records for the new N5/N6 diagnostics. **Fix:** give implementation cross-checks their own section/status, outside the PMP conjunction, and extend the diagnostic record to carry the new checks—or explicitly mark them not checked. Do not infer their completion from the older fields.

- **[OVERCLAIM]** `transfer_study.m:346–353`; `report_optimality.m:110–112` - The study’s conditional theorem statement and evidence qualification are substantially better, but “every hypothesis … was checked and holds” still asserts more than the sampled checks establish. The reusable report retains the old “numerically certified at the sampled times” claim unchanged. **Fix:** use the same wording in both: **“All required numerical checks passed. If the stated hypotheses hold exactly, this arc is a strict strong local minimizer …”**, followed by the sampling, rank and conjugate-instrument limitations.

- **[CORRECTNESS]** `transfer_study.m:303–329,355–359` - `ENDPOINT` becomes UNRESOLVED in the S4 row, then becomes “a sufficiency hypothesis fails” in the final verdict. An H6 failure is also incorrectly described there as a failed theorem hypothesis. **Fix:** preserve separate statuses through aggregation: **passed**, **failed**, and **unresolved/not checked**. An inconclusive conjugate result or unavailable reduction validity should end with **“sufficiency unresolved; no minimality claimed.”** `certify_root` may correctly refuse acceptance in these cases, but should preserve the reason rather than flatten everything into a numeric conjugate failure.

- **[CORRECTNESS]** `transfer_study.m:153–171` - “Admissible” is printed without enforcing either closest-approach quantity. The script also lacks the completed-flight and positive-final-time checks already present in `certify_root`; finite positive-mass samples can describe an integration that stopped early. **Fix:** share a flight validator that checks time coverage, positive \(t_f\), finite trajectory data and the all-burn mass law. Reject sampled body penetrations, and use collision events/dense-output closest-approach searches for between-sample clearance. Until then, label these values **sampled clearances**, not established admissibility.

- **[CORRECTNESS]** `transfer_study.m:95–113,372–385` - The periodicity repair is optional: missing or failing `csape` restores precisely the nonperiodic spline the comments reject. The seam diagnostic checks values, not the derivative discontinuity motivating the repair, and orbit closure is printed but not enforced. **Fix:** provide a toolbox-free periodic construction or fail explicitly when periodic interpolation is required. Validate closure before making a bounded endpoint reconciliation, and check seam values and derivatives. The fallback’s printed name is honest, but it does not fulfil the advertised interpolation contract.

- **[CORRECTNESS]** `transfer_study.m:146,180,187–194,271–272` - Thresholds are not yet a single source of truth. Verification runs before `tol` exists and uses its own defaults; the summary prints `tol.dz` but takes its pass/fail from `V.moved`. Editing only `tol.dz` can therefore produce a PASS beside a value exceeding the printed threshold. **Fix:** define tolerances before their first consumer, pass them into verification and validation, and derive displayed thresholds and statuses from the same configuration.

- **[CORRECTNESS]** `verify_with_pumpkyn.m:80–85,97–101`; `transfer_study.m:271–272` - Independent verification can report agreement as PASS even when the witness flight fails or is unavailable: `V.moved` is settled before that flight, and its failure does not affect the study’s summary. The witness flight also checks position only. **Fix:** retain `agreementOk` as one metric, but add a separate overall verification status requiring a usable, completed witness flight and position **and velocity** acceptance. Keep that status outside `necessary`; separation from PMP should not conceal an unsuccessful cross-check.

#### Supporting mathematical diagnostics

- **[CORRECTNESS]** `conj_spectrum.m:18–25,82–88` - A determinant near \(10^{-16}\) does **not** make its sign meaningless. With the stated singular values, \(\kappa_2(M)\approx1.7\times10^7\), so ordinary double-precision matrix evaluation can still resolve the sign comfortably; uncertainty in the integrated matrix is the relevant question. Excluding `t_f` also leaves the final sampling interval untested for sign changes. **Fix:** assess \(\sigma_{\min}(M)\) against matrix-error estimates/bounds, not determinant magnitude. If the endpoint is intentionally excluded, explicitly report that uncovered interval and do not treat this scan as covering \((0,t_f]\).

- **[OVERCLAIM]** `conj_spectrum.m:4–11,90–99` - A denser finite scan reduces the two blind spots; it does not close them. Moreover, a globally small \(\sigma_6\) plus \(\sigma_6/\sigma_5>0.1\) does not establish that two singular values locally collapsed. The option documentation gives the reciprocal inequality, which is impossible for descending singular values. **Fix:** call this a **candidate-detection scan**, correct the documented ratio, and refine/localise suspected events using each relevant singular value’s behaviour and numerical uncertainty. Do not equate absence of candidates with absence of conjugate times.

- **[OVERCLAIM]** `lift_margin.m:7–25,72–89` - The rigorous implication from \(\sigma_6>\|\Delta C\|\) is **rank at least six**, not rank exactly six. Exact rank six additionally needs an exact nonzero null vector; a small residual is not that witness. Your header correctly admits that the two-setting difference is only an estimate, but the returned “CERTIFIED … (Eckart–Young)” wording drops both qualifications. **Fix:** report **numerically supported nullity one / estimated rank margin**, and state exact rank conditionally on a genuine error bound and exact lift.

- **[CORRECTNESS]** `lift_margin.m:53–64` - The purported constructive lift check accepts `lam = zeros(7,1)` because its residual is zero after the `realmin` denominator guard. **Fix:** reject nonfinite or zero lift vectors before normalisation; a homogeneous null-space witness must be nonzero.

## B — Readability and reuse

- **[FINE]** `transfer_study.m:10–27` - The numbered, visibly computed study still has a good narrative; retain the departure → target → solve → audit → plot structure rather than replacing it with a front-door wrapper.

- **[READABILITY]** `transfer_study.m:121–152,226–232,267–270,320–345` - The mathematical narrative now competes with the history of the review. Dates, “YOUR objection”-style references, `FINDINGS 40`, and accounts of previous implementations make a reusable template read like a remediation log. **Fix:** retain the short mathematical reasons—why fixed-control differentiation, why separate first/second order, why conditional claims—and move incident history and review attribution into the audit document or regression tests.

- **[READABILITY]** `transfer_study.m:49–50,70–71,149–156,294,366–368` - The new objects do not yet establish ownership. `flight` is constructed but never consumed; the hypothesis routine and plotting routine still propagate again. Departure/arrival objects coexist with the original arrays, anonymous functions, `B`, and a later `T`. This is where a reader loses track of which representation is authoritative. **Fix:** make the grouped objects actual inputs to downstream consumers, including optional supplied-flight inputs to gates and plotting. Keep the independent witness flight distinct and identify the shooting mesh as the separate source of the conjugate diagnostic.

- **[READABILITY]** `transfer_study.m:281–325`; `report_optimality.m:89–96` - The reader must translate S1→H2, S3→H1, S4→H5 and S5→H6, while S2 is a PMP control requirement and S5 explicitly is not a theorem hypothesis. **Fix:** use stable diagnostic identifiers across the script and library, and distinguish **PMP**, **theorem hypotheses**, and **instrument validity** in the printed grouping. This can remain a short numbered study; it does not need another abstraction layer.

- **[READABILITY]** `transfer_study.m:314–318,349–354` - Coverage is more visible, but “first full-rank at 7,” for example, is an index, not a time span. The verdict also counts `size(it.Y,2)` rather than the instrument’s actual determinant evaluations. **Fix:** print the first and last evaluated times, the first resolved full-rank time as \(t/t_f\), the unresolved initial interval, and the actual determinant sample count. Report the rank threshold used beside the rank-related quantities.

- **[READABILITY]** `transfer_study.m:3,116–133` - “Edit the parameter blocks, press Run” promises broader reuse than the fixed operating-point seed gate permits; the reader discovers the restriction only after generating both orbits. **Fix:** state the seed prerequisite at the top and distinguish **running the shipped seeded example** from **studying a continuation-generated transfer**.

- **[READABILITY]** `h6_margin.m:44–46,58–59`; `mintime_hypothesis_gates.m:132` - `hMin` is the **maximum** reduced Hamiltonian, attained initially—the value closest to zero when H6 holds. **Fix:** rename it `hMax` and preferably expose the positive clearance `minusHMax`; that makes the inequality direction immediately readable.