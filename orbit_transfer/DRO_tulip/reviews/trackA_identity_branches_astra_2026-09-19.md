## Summary

The **phase-sensitivity signs, mass-costate equation, cancellation in \(\dot Q_{mt}\), and Hamiltonian gap identity are correct**. The principal problems are in what the numerical tests are claimed to establish:

- Endpoint samples of a slope bound do **not** produce a certified between-sample bound.
- A small flight-time residual does **not** identify a smooth costate branch.
- Arc attachment can silently revert to time-only matching, and its ambiguity handling is inadequate for suppressing new continuation runs.
- The audit is associated with the final catalog by location and count, not by content.
- The supplied evidence supports numerical consistency and improvement over the reference—not a completed mathematical sufficiency proof or global minimum-time optimality.

I reviewed the supplied source, not an executed checkout. The BCT/O–M papers and the referenced determinant-equivalence section were not supplied; I therefore distinguish checking the mathematics here from independently verifying a verbatim theorem quotation. Driver citations use the **new-side diff line numbers**. The unnumbered LaTeX is cited by section/paragraph.

---

## 1. Phase sensitivity and X3

### Derivation

- **[MATH: CORRECT]** `phase_sensitivity.m:9–22,73–74` — Both signs are correct.

  Use the augmented functional
  \[
  \mathcal A=t_f+\int_0^{t_f}\lambda^\top(f-\dot x)\,dt.
  \]
  On an extremal, its endpoint variation is
  \[
  \delta J
  =\lambda(0)^\top\delta x_0
   -\lambda(t_f)^\top\delta x_f
   +H(t_f)\,\delta t_f,
  \]
  where \(\delta x_f\) is the **total endpoint displacement**, including the effect of changing terminal time. Thus
  \[
  \boxed{\partial_{s_D}t_f=\lambda_{rv}(0)^\top x_D'(s_D)},\qquad
  \boxed{\partial_{s_A}t_f=-\lambda_{rv}(t_f)^\top x_A'(s_A)}.
  \]
  The mass contributions vanish because \(\partial_s m_0=0\) and \(\lambda_m(t_f)=0\). Free terminal time contributes no additional term because \(H(t_f)=0\). These conclusions use the prescribed phase endpoint maps, not a target with additional explicit dependence on \(t_f\).

- **[MATH: INCOMPLETE]** `phase_sensitivity.m:16–19,30–37` — “\(\lambda\) is the gradient of the cost-to-go” needs a differentiability qualification. It is correct at a differentiable value function with the appropriate optimal lift. More generally, these endpoint formulas give the derivative of the objective **along a differentiable, regular extremal branch**, even if that branch is not globally optimal. At a crossing of competing optimal branches, the minimum-time value can be nonsmooth.

  **Correction:** document the formula as a smooth-branch sensitivity; identify it with the minimum-time value gradient only when that value is differentiable and the branch realizes it.

- **[MATH: CORRECT]** `phase_sensitivity.m:38–41` — \(x'(s)=\tau f_{\rm ballistic}(x(s))\) is correct when phase is forward flow time divided by the orbit period. Using the derivative of the actual solver closure, as the code attempts, is important.

### What the finite-difference test establishes

- **[CORRECTNESS]** `phase_transversality_check.m:11–16,180–212`; `phase_sensitivity.m:24–28` — “No costate is read” and “validates \(\lambda(0)\) and \(\lambda(t_f)\) independently” overstate the test. `seed_from_z8` receives the stored costates, and the re-solve is an indirect PMP solve. The comparison observable uses only returned flight times, but the calculation itself uses costates and the same problem implementation. Moreover, each derivative checks **one scalar projection**, not the whole endpoint costate.

  **Concrete fix:** call this a *finite-difference endpoint-sensitivity consistency test*. Say: “The finite-difference observable does not use the returned costates.” Reserve “independent oracle” for the analytic single-integrator test or a genuinely independent formulation.

- **[ROBUSTNESS]** `phase_transversality_check.m:188–212` — Convergence of both displaced solves does not establish that they stayed on the original branch. A converged branch jump can invalidate the derivative comparison in either direction.

  **Concrete fix:** retain the perturbed solution vectors and residuals; check continuation continuity, preferably with a shooting-Jacobian tangent predictor. Return **UNRESOLVED—branch identity uncertain** rather than interpreting every converged discrepancy as a sensitivity failure.

- **[MATH: INCOMPLETE]** `phase_transversality_check.m:79,189,211` — \(h=2\times10^{-4}\) is a reasonable starting step, but the residual tolerance does not by itself justify its accuracy.

  For a smooth branch,
  \[
  \frac{T(s+h)-T(s-h)}{2h}
  =T'(s)+\frac{h^2}{6}T'''(s)+O(h^4).
  \]
  If each computed time has error at most \(\epsilon_T\), the additional derivative error is bounded by approximately \(\epsilon_T/h\). The relevant truncation quantity is **\(T'''\)**, not merely curvature \(T''\).

  If—only illustratively—\(\epsilon_T=3\times10^{-11}\), then
  \[
  \epsilon_T/h=1.5\times10^{-7}.
  \]
  But a shooting residual of \(3\times10^{-11}\) is not a time-error bound. Locally,
  \[
  |\delta t_f|
  \lesssim \|e_{t_f}^{\top}R_z^{-1}\|\,\|R\|,
  \]
  plus integration and endpoint-map errors.

  A seven-cycle sinusoidal dependence would have central-difference relative truncation around
  \[
  (2\pi\,7h)^2/6\approx1.3\times10^{-5},
  \]
  so the step is plausibly small relative to the tulip’s basic phase scale. That does not bound sharper branch structure or behavior near folds.

- **[ROBUSTNESS]** `phase_transversality_check.m:189–203` — The half-step is used only after failure; it supplies no truncation-error estimate when the first attempt converges. The report also always prints `X.delta`, even if a derivative actually used `delta/2`.

  **Concrete fix:** for validation cells, calculate both \(D_h\) and \(D_{h/2}\), record the actual step per derivative, and check second-order convergence. Under the smooth asymptotic model,
  \[
  |D_{h/2}-T'|\approx |D_h-D_{h/2}|/3.
  \]
  Tighten the solve tolerance separately to distinguish truncation from solve error.

- **[ROBUSTNESS]** `phase_transversality_check.m:272–275` — A purely relative gate with denominator floor \(10^{-12}\) is inappropriate near the phase-stationary points this tool is intended to investigate.

  **Concrete fix:** use
  \[
  |D_{\rm formula}-D_{\rm FD}|
  \le \text{absTol}+\text{relTol}\max(|D_{\rm formula}|,|D_{\rm FD}|),
  \]
  with `absTol` informed by finite-difference uncertainty.

- **[MATH: CORRECT]** `phase_transversality_check.m:221–234` — The three-valued aggregation is right: a resolved disagreement takes precedence over unresolved derivatives; unresolved or empty tests never pass. This assumes “resolved” also includes branch identity and adequate numerical accuracy.

  **For today’s supplied numbers:** six departure and five arrival derivatives resolved; all eleven agree within \(10^{-3}\); one arrival derivative is unresolved. The overall verdict is therefore **UNRESOLVED**, not PASS. Departure discrepancies are \(6.4\times10^{-8}\) to \(5.0\times10^{-7}\); the worst resolved arrival discrepancy is \(3.9\times10^{-4}\).

- **[ROBUSTNESS]** `phase_transversality_check.m:167–175` — The claimed \(10^{-10}\) accuracy of the departure tangent is unsupported by choosing a \(10^{-6}\) difference step.

  **Concrete fix:** differentiate the spline closure analytically, or verify the tangent with a step ladder. Also verify periodic wrapping at the closure boundary.

---

## 2. H2/H3 and between-sample bounds

### The differential identities

- **[MATH: CORRECT]** `mintime_hypothesis_gates.m:125–129` — In the stated rotating coordinates,
  \[
  C=\partial_v f_v=
  \begin{pmatrix}0&2&0\\-2&0&0\\0&0&0\end{pmatrix},
  \qquad C^\top=-C,\qquad \|C\|_2=2,
  \]
  and
  \[
  \dot\lambda_v=-\lambda_r-C^\top\lambda_v.
  \]
  Consequently the stated bound
  \[
  \left|\frac{d}{dt}|\lambda_v|\right|
  \le|\lambda_r|+2|\lambda_v|
  \]
  is valid.

  There is, however, a **strictly better bound at no additional computational cost**:
  \[
  \frac{d}{dt}|\lambda_v|
  =-\widehat{\lambda_v}^{\top}\lambda_r,
  \qquad
  \boxed{\left|\frac{d}{dt}|\lambda_v|\right|\le|\lambda_r|}.
  \]
  The Coriolis contribution vanishes because it is skew-symmetric. At a zero of \(\lambda_v\), use the corresponding almost-everywhere norm/Lipschitz statement rather than the divided formula.

- **[MATH: CORRECT]** `mintime_hypothesis_gates.m:43–45,126–127` — The mass-costate sign is correct:
  \[
  \dot\lambda_m=-H_m
  =\frac{T}{m^2}\lambda_v^\top\alpha^*
  =-\frac{T|\lambda_v|}{m^2}.
  \]
  Hence, with \(\nu=|\lambda_v|\),
  \[
  \dot Q_{mt}
  =\frac{\dot\nu}{m}
   -\frac{\nu\dot m}{m^2}
   +\frac{\dot\lambda_m}{c}
  =\boxed{\frac{\dot\nu}{m}}.
  \]
  Both mass terms cancel exactly.

- **[MATH: INCOMPLETE]** `mintime_hypothesis_gates.m:11–12` — “\(s^*=1\iff Q_{mt}>0\)” is correct only if it means **the unique minimizing throttle is 1**. At \(Q_{mt}=0\), every throttle minimizes the Hamiltonian, including 1.

  A useful additional consequence in this specific model is
  \[
  \lambda_m(t)=\int_t^{t_f}\frac{T\nu(\tau)}{m(\tau)^2}\,d\tau\ge0.
  \]
  Thus, for an exact all-burn extremal with positive mass and terminal transversality, H2 implies \(Q_{mt}>0\). H3 remains useful as a numerical consistency check, but the conditions are not independent here.

### The interval inequality and its implementation

- **[MATH: CORRECT]** `between_sample_bound.m:9–14` — If \(L\) genuinely bounds \(|v'|\) throughout the interval, then
  \[
  v(t)\ge\max\{v_k-L(t-t_k),\,v_{k+1}-L(t_{k+1}-t)\},
  \]
  giving
  \[
  \boxed{\min v\ge\frac{v_k+v_{k+1}-L\Delta t}{2}}.
  \]
  The endpoint compatibility condition \(|v_{k+1}-v_k|\le L\Delta t\) follows automatically from a valid Lipschitz bound. It is also a useful implementation check.

- **[CORRECTNESS]** `between_sample_bound.m:18–23,51`; `mintime_hypothesis_gates.m:89–92,123–131`; `certify_root.m:387–400` — The supplied \(L\) is not an interval bound. Taking the maximum of two endpoint bounds does not control an interior maximum. The helper acknowledges this, but its consumers promote the estimate back into a “LOWER BOUND over the whole arc.”

  For example, endpoint derivative samples cannot distinguish a constant function from
  \[
  v(t)=1-a\sin^2(\pi t/h),
  \]
  whose endpoint derivatives vanish but whose interior minimum is \(1-a\).

  **Concrete fix:** either:
  1. rename these outputs `minLamVEstimate`/`minQmtEstimate` and stop claiming continuous-time certification; or
  2. supply a verified interval slope bound and endpoint error enclosures.

  Seven thousand samples and a one-percent correction are reassuring numerical evidence, not the missing theorem.

- **[ROBUSTNESS]** `between_sample_bound.m:46–52` — The helper does not require finite real samples, finite real nonnegative slopes, or finite times. Missing/nonfinite interval values can contaminate or be omitted by reduction operations, depending on the operation’s behavior.

  **Concrete fix:** reject every nonfinite/nonreal input and every negative slope bound before forming intervals. If this becomes a certifying helper, accept **one bound per interval**, validate endpoint compatibility with error allowances, and subtract endpoint sample uncertainty.

### A cheap mathematically rigorous interval bound

- **[MATH: CORRECT]** `mintime_hypothesis_gates.m:128–131` — A comparison-system bound can replace endpoint maximization.

  Write
  \[
  K(t)=\partial_r f_v,\qquad
  a(t)=|\lambda_r(t)|,\qquad b(t)=|\lambda_v(t)|.
  \]
  If \(\|K(t)\|\le\kappa\) on an interval of length \(h\), then, in the appropriate norm-derivative sense,
  \[
  \dot a\le\kappa b,\qquad \dot b\le a.
  \]
  Starting with upper bounds \(a_0,b_0\),
  \[
  \begin{pmatrix}a(t)\\b(t)\end{pmatrix}
  \le
  \exp\!\left[
  \begin{pmatrix}0&\kappa\\1&0\end{pmatrix}(t-t_k)
  \right]
  \begin{pmatrix}a_0\\b_0\end{pmatrix}.
  \]
  In particular,
  \[
  a(t)\le
  a_0\cosh(\sqrt\kappa h)
  +\sqrt\kappa\,b_0\sinh(\sqrt\kappa h)
  =:A_k.
  \]
  Therefore valid choices are
  \[
  L_{\nu,k}=A_k,\qquad
  L_{Q,k}=A_k/m_{\min,k},
  \]
  with \(m_{\min,k}=m(t_{k+1})>0\) for the exact all-burn mass law.

  A conservative CR3BP Hessian bound is
  \[
  \kappa\le
  1+\frac{2(1-\mu)}{d_{1,\min}^3}
   +\frac{2\mu}{d_{2,\min}^3}.
  \]
  **The remaining requirement is real:** the primary-distance lower bounds and initial costate bounds must hold over the interval, not merely at its endpoints. Validated state enclosures, or independently established clearance bounds plus integration-error control, can supply them. Ordinary ODE tolerances alone do not make the result validated.

- **[MATH: INCOMPLETE]** `certify_root.m:387–391`; audit §gaps(b) — “Ten times the \(\lambda_m\) uncertainty” does not establish uncertainty margins for both tested quantities. With exact mass,
  \[
  \epsilon_Q\le\epsilon_{\lambda_v}/m+\epsilon_{\lambda_m}/c;
  \]
  uncertainty in mass adds another term. The error in \(|\lambda_v|\) is governed by the error in \(\lambda_v\), not by the error in \(\lambda_m\).

  **Correction:** propagate componentwise or normwise uncertainty into each gate. The observed margins \(0.2785\) and \(0.3126\) are large, but the stated justification of the \(10^{-5}\) floor is not sufficient.

---

## 3. Edge residuals and “branches”

- **[MATH: CORRECT]** `phase_edge_residuals.m:8–11,50` — On one \(C^3\) phase branch with exact \(G=T'\),
  \[
  r=-\frac{h^3}{12}T'''(\xi)
  \]
  for some interior \(\xi\), and the stated absolute bound follows. Numerical endpoint-time and gradient errors must be added to that bound.

- **[MATH: OVERSTATED]** `phase_edge_residuals.m:11–23`; `phase_transversality_check.m:17–21,121–128` — The converse inference is false. Small \(r\) does not establish a common branch, and a fixed 5- or 10-minute threshold does not establish departure from one.

  **Different branches can give small or zero residual because:**
  - distinct extremals can have identical time and phase sensitivities;
  - only a scalar costate projection is tested;
  - a time jump can cancel the trapezoidal sensitivity term;
  - branches can cross or approach within the threshold.

  An explicit cancellation example is
  \[
  T_1(s)=T_0,\qquad T_2(s)=T_0+b(s-h/2).
  \]
  Using branch 1 at \(s=0\) and branch 2 at \(s=h\) gives \(r=0\), although the two branch values at \(s=0\) differ.

  **One branch can give large residual because:** \(T'''\) is large, the grid is coarse, phase ceases to be a regular coordinate near a fold, endpoint derivatives are inaccurate, or the stored solutions/gradients are insufficiently polished.

- **[CORRECTNESS]** `phase_transversality_check.m:121–129,254–263`; `phase_branches.m:4–15,23–25` — The labels `SAFE`, `JUMP`, and `BRANCH` turn an inexpensive scalar consistency heuristic into a costate-interpolation guarantee.

  **Concrete fix:** rename them to, for example, `timeConsistentEdge`, `largeTimeResidual`, and `candidateComponent`. If downstream behavior actually requires common-branch evidence, add a full-solution continuation check.

  A useful check is to differentiate the shooting equations:
  \[
  z_s=-R_z^{-1}R_s,
  \]
  then compare the full \(z_8\) with a predictor/corrector or Hermite consistency residual. A midpoint solve and a conditioning check on \(R_z\) are much stronger evidence than \(T\) alone. Near folds, use pseudo-arclength rather than phase as the regular coordinate.

- **[MATH: OVERSTATED]** `phase_branches.m:4–8,47–62` — The graph algorithm correctly computes connected components of the supplied masks. But even genuinely local interpolation-compatible edges do not imply that arbitrary members of a component may be interpolated directly: the component may wind around the torus, follow a fold, or cover a nonconvex region.

  The measured result is therefore **525 and 9 threshold-accepted edges, forming 50 graph components**, not an established count of 50 mathematical solution branches.

- **[CORRECTNESS]** `phase_transversality_check.m:29–34,128–129` — Arrival edges are declared “not judged” because the axis is under-resolved, but nine such edges are nevertheless used to merge components. An accidental small residual in the under-resolved direction can join otherwise unrelated components.

  **Concrete fix:** do not use those arrival edges as positive branch evidence without refinement/full-solution checks. The large residual distribution is compatible with under-resolution, but does not establish that under-resolution is its sole cause.

- **[ROBUSTNESS]** `phase_edge_residuals.m:44–47` — The routine does not enforce finite, distinct, consistently circularly ordered phases. A singleton axis also gets spacing zero because `mod(1,1)=0`, making the self-edge trivially consistent.

  **Concrete fix:** validate the phase ordering and uniqueness; reject or mark singleton-axis residuals uninformative. Do not silently interpret arbitrary permutations as a nearest-neighbor phase grid.

---

## 4. Root identity and arc attachment

### What is already right

- **[MATH: CORRECT]** `family_map.m:137–142,276–278` — Comparing all seven initial costate components is appropriate for identifying the **full normal lift** at the same initial state and problem. Although \(\lambda_m\) does not affect the forced all-burn direction dynamics, it participates in normality and terminal transversality. Discarding it would lose information.

  A positive homogeneous scale does not affect direction. At a fixed initial state, a valid normal-chart direction also determines its scale through \(\lambda^\top f=-1\). The Euclidean angular metric is nevertheless coordinate-scaling dependent and can underweight small but sensitive components.

- **[MATH: CORRECT]** `family_map.m:271–274` — The final-point exact-hit case is guarded in the supplied version. If \(c=\text{numel}(q)\) is in `on`, `a=0` is selected without evaluating `s(c+1)`, and `c2` is clamped. There is no demonstrated end-index overrun here.

- **[MATH: CORRECT]** `family_map.m:280–282` — With finite nonnegative gaps and the default `tolDays=0.02`, a failing crossing cannot shadow a passing crossing: a passing score is at most \(0.02\), while a failing score is at least \(10^6\). This particular earlier failure mode is fixed.

### Remaining identity defects

- **[CORRECTNESS]** `family_map.m:263–282` — Missing, zero, or malformed costates can silently become a successful **time-only** match. In particular, a nonfinite costate can produce `cg=NaN`, and `isnan(cg)` explicitly permits the crossing. Interpolation cancellation to a zero vector also bypasses the direction comparison.

  **Concrete fix:** separate “not supplied,” “invalid,” and “valid comparison.” When attachment controls whether to suppress a new arc, require finite real nonzero costates on both sides. Missing/invalid data must give **UNKNOWN**, not “same root.” Keep any time-only fallback explicitly diagnostic.

- **[CORRECTNESS]** `family_map.m:264,280–282` — The comment says “costates first, then time,” but passing candidates are ranked solely by time gap. Exact score ties keep the first listed crossing. Two different families can both pass, yet the function reports one definitive family.

  **Concrete fix:** collect passing crossings, distinguish repeated hits of one family from multiple-family ambiguity, and return an ambiguity status. Use a documented normalized distance only to prioritize a confirming polish—not to make ambiguous identity disappear. Replace the \(10^6\) penalty with explicit pass/fail ordering.

- **[CORRECTNESS]** `family_map.m:127,137–142,276–278` — Taking `abs(rho)` discards information needed to interpret the homogeneous chart. If
  \[
  \lambda_h=\rho\lambda_n,\qquad \rho<0,
  \]
  then the same numerical ray becomes antiparallel and the direction gap is 2, not 0.

  More importantly, negative multiplication of PMP multipliers is **not** an admissible normal-chart gauge under the stated minimization convention: it reverses the minimizing thrust direction if the control is evaluated directly from \(\lambda_h\).

  **Concrete fix:** preserve signed \(\rho\); use \(|\rho|\) only for the displayed proximity-to-zero diagnostic. Require a documented positive normal chart, and exclude zero-crossing/near-zero chart segments from identity decisions. If the continuation formulation intentionally permits signed coordinate gauges, convert using the signed chart transformation and verify that its control law has the corresponding convention.

- **[MATH: INCOMPLETE]** `family_map.m:272–278` — Phase spacing \(2\times10^{-4}\) alone does not bound costate interpolation error. On a regular phase-parametrized arc,
  \[
  \|\lambda_{\rm interp}-\lambda(s)\|
  \le \frac{h^2}{8}\sup\|\lambda''(s)\|.
  \]
  Neither that curvature nor its amplification by normalization is bounded here. Near a fold, phase derivatives can become large. Varying homogeneous scale also changes the interpolation error.

  **Concrete fix:** use interpolation to produce a candidate, then polish that candidate at the requested phase. Convert to a common positive normal chart first where well-conditioned, or use the actual continuation chart with error control. Cache confirmations.

- **[MATH: INCOMPLETE]** `family_map.m:39–48,107–109`; `run_phase_torus.m:851–857` — The two identity tolerances differ enormously:
  \[
  1-\cos\theta\le10^{-4}
  \quad\Longleftrightarrow\quad
  \theta\lesssim0.01414\ {\rm rad}=0.81^\circ,
  \]
  corresponding to a unit-vector discrepancy of about **1.4%**. That is far looser than a \(10^{-6}\) relative normal-costate comparison.

  Also, \(0.02\) day is **28.8 minutes**, whereas the claimed \(10^{-5}\)-day arc interpolation error is \(0.864\) seconds. These may be useful candidate-search tolerances, but are not established root-identity tolerances.

  Distinct roots need not differ “at order one”; near a fold, they can be arbitrarily close. Conversely, the same root can move more under an ill-conditioned solve.

  **Concrete fix:** calibrate within-root discrepancies and nearest distinct-root separations after common polishing, including Jacobian conditioning. Use the loose test for candidate generation and a tighter, uncertainty-aware full-root test for suppressing continuation.

### Registry and `sameRoot`

- **[CORRECTNESS]** `run_phase_torus.m:842–843` — Registry de-duplication tests arrival phase but not the stored departure phase, although `registerRoot` receives and stores `sD0`.

  **Concrete fix:** include both phases and the problem identity. If this registry is intentionally spine-only, enforce that invariant on every loaded and inserted record.

- **[ROBUSTNESS]** `run_phase_torus.m:851–857` — `sameRoot` ignores flight time, lacks finite/nonzero input checks, and uses an asymmetric relative tolerance.

  There is an important mathematical qualification: for **exact**, same-initial-state, strict all-burn solutions, identical initial costates do determine the trajectory, and
  \[
  \dot\lambda_m=-T|\lambda_v|/m^2<0
  \]
  makes the zero \(\lambda_m(t_f)=0\) unique. Thus ignoring \(t_f\) is not an exact counterexample by itself in this model. But **approximate** costate equality does not supply a time-error bound.

  **Concrete fix:** require finite real normal-chart \(z_8\), compare \(t_f\) separately, and use symmetric absolute-plus-relative costate tolerances. Near the threshold, re-polish both candidates under the same settings rather than assuming identity.

The practical recommendation is a **two-stage identity rule**: cheap time/direction candidate search, followed by a fixed-phase confirming polish. That is much cheaper than unnecessary hours of continuation and much safer than suppressing a genuinely new arc.

---

## 5. Audit document and theorem scope

### BCT and O–M

- **[MATH: INCOMPLETE]** audit §hyp, BCT theorem quotation — I cannot independently authenticate the theorem number, exact quotation, or topology definitions without the paper. On the supplied statement, the document correctly distinguishes local \(C^0\) optimality from an unsupported claim of **strict** strong local optimality. It also correctly acknowledges that the stated strong-regularity hypothesis is stronger than a whole-arc numerical rank test.

  The claimed equivalence between Test 3 and the implemented determinant still depends on the unsupplied §equiv and on handling the nonautonomous mass reduction correctly.

- **[MATH: CORRECT]** audit §gaps(a) — \(S^2\) is open in itself as a control manifold. Thus the direction-only problem can fit the stated open-control setting, while the throttle extremal at \(s=1\) lies on the boundary of \([0,1]\). The original throttled problem is not directly covered merely by applying BCT to the direction-only problem.

- **[MATH: INCOMPLETE]** audit §gaps(a) — O–M’s mixed continuous/bang-bang framework is a natural relevant framework, and “cited, not instantiated” is an appropriate limitation. But a sphere-chart argument is not the only missing step. One must identify the applicable theorem and verify its endpoint qualifications, regularity, second-variation/coercivity requirements, free-time treatment, and the topology of the asserted minimum.

  With a uniform strict throttle margin and feasible first-order variations \(\delta s\le0\), the throttle contribution is
  \[
  \int_0^{t_f}-TQ_{mt}\,\delta s\,dt\ge0,
  \]
  and it vanishes only when \(\delta s=0\) almost everywhere. This explains why the critical directions have no throttle component. With no switches, there are no switching-time variables. The remaining accessory problem still includes state, direction, endpoint, and free-time variations—not merely the pointwise angular Hessian.

- **[MATH: CORRECT]** audit §gaps(a), Hamiltonian gap — For unit \(\alpha,\alpha^*\), \(m>0\), and \(\lambda_v\ne0\),
  \[
  \boxed{
  H(s,\alpha)-H(1,\alpha^*)
  =TQ_{mt}(1-s)
   +\frac{sT|\lambda_v|}{2m}\|\alpha-\alpha^*\|^2
  }.
  \]
  This follows from \(\alpha^*=-\widehat{\lambda_v}\) and
  \[
  1-\alpha^\top\alpha^*=\tfrac12\|\alpha-\alpha^*\|^2.
  \]

  The wording about “uniform lower bounds on its two coefficients” needs care: the second displayed coefficient includes \(s\) and vanishes at \(s=0\). H2 supplies the positive angular Hessian **at the reference throttle \(s=1\)**, not a uniform positive angular coefficient over all throttles.

### Remaining theorem overclaims

- **[MATH: OVERSTATED]** audit §gaps(c) — Prefixes and “a few interior windows” do not check a condition quantified over **every subinterval**.

  **Correction:** call those additional numerical probes. A proof exploiting analyticity of the dynamics, reference control, and stationary adjoints away from collisions and \(\lambda_v=0\) may be able to propagate an interval rank property, but that argument must actually be supplied.

- **[MATH: OVERSTATED]** `mintime_hypothesis_gates.m:14–21` — The linear space \(S\) imposes stationarity/parallelism, not the sign inequalities required for Hamiltonian minimization. Therefore `dim S = 1` is a sufficient exclusion of nonzero abnormal stationary lifts in the stated construction, but the claimed converse “no abnormal PMP lift iff \(\dim S=1\)” is not established: a vector in the linear kernel may violate direction or throttle minimization.

  **Correction:** distinguish *stationary lifts* from *PMP-minimizing lifts*. Keep the rank test as a sufficient exclusion test, with the separate numerical-rank qualifications.

- **[MATH: OVERSTATED]** audit §reductions, “The reduction is nevertheless sound, and what carries it is H3” — Strict first-order Hamiltonian minimization does not, by itself, instantiate a second-order sufficiency theorem for the full constrained-control problem. This older assertion conflicts with the new subsection’s honest admission that the reduction is not instantiated.

  **Correction:** say the reduction is justified at the level of critical directions under strict bang conditions; sufficiency for the original competitor class remains conditional on the applicable constrained-control theorem.

- **[MATH: WRONG]** audit §reductions, “the two Hamiltonians agree, \(\lambda_m\) playing the role of the time costate” — The equivalence requires a scaling and an added clock term. Put \(k=T/c\) and \(\tau=(m_0-m)/k\). The clock costate is
  \[
  p_\tau=-k\lambda_m,
  \]
  and
  \[
  H_7=H_6+p_\tau=H_6-k\lambda_m,
  \qquad H_6=1+\lambda_{rv}^{\top}f_{rv}.
  \]
  Thus \(H_7=0\) gives \(H_6=k\lambda_m\), not equality of the two Hamiltonians. The autonomous extended formulation is equivalent under this transformation; the unextended Hamiltonians are not identical.

- **[MATH: CORRECT]** audit §reductions, integral identity and invariances — For \(m_0=1\) and \(\lambda_m(t_f)=0\),
  \[
  \int_0^{t_f}TQ_{mt}\,dt=\lambda_m(0).
  \]
  More generally the right side is \(m_0\lambda_m(0)-m_f\lambda_m(t_f)\). Positive-scale trajectory invariance and independence from the initial mass costate are correct for the **forced all-burn trajectory map**. These identities are consistency checks, not guarantees that every underlying error will be detected.

### Fixed-phase scope and evidence

- **[MATH: CORRECT]** audit §gaps(d) — A fixed-phase second-order test does not certify a free-phase optimum. Vanishing of both phase sensitivities is the correct first-order condition on a regular branch.

- **[MATH: OVERSTATED]** audit §gaps(d) — “None of the 576 is [stationary] … so the free-phase minimum lies between grid points” goes too far. It establishes that none of the **stored roots** is a regular free-phase stationary candidate. It does not exclude an unrecorded better branch with a stationary root at the same grid phases, nor establish existence/location of the global orbit-to-orbit minimum.

- **[CORRECTNESS]** audit §gaps(d), X3 evidence paragraph — The claim that twelve derivatives agree is inconsistent with today’s supplied evidence: eleven resolved, one unresolved. The error range and “4–5 digits” description should also be updated.

**Most misleading sentence in the new subsection:** the §gaps(b) sentence asserting the interval inequality with \(L_k\) defined as the maximum of the two endpoint slope bounds. That inequality does not follow from that definition. The later disclaimer is valuable, but does not make the earlier mathematical assertion—or the consumer’s “whole-arc lower bound” language—true.

---

## 6. Catalog comparison and rebuild audit

### Catalog matching

- **[MATH: CORRECT]** `compare_phase_catalogs.m:81–87,169–180,125–128` — For these well-separated phase grids, circular nearest matching plus the bijection check correctly constructs the permutation. Reordering `entry_index` rather than the underlying `z8` columns is correct. Faster/slower signs are correct, with differences within `tolTfDays` deliberately treated as ties.

- **[CORRECTNESS]** `compare_phase_catalogs.m:72,153–165` — The problem-identity comparison is incomplete. Equal thrust, Isp, initial mass, orbit labels/parameters, and time scale do not establish equal dynamics and endpoint maps.

  **Concrete fix:** compare a canonical problem fingerprint including \(\mu\), relevant dimensional scales, orbit initial conditions/definitions, phase origins/orientations, and endpoint-closure identity. Explicitly select or require a single thrust rung.

- **[CORRECTNESS]** `compare_phase_catalogs.m:184–199` — Greedy largest-overlap matching does not minimize the number of family differences. For example,
  \[
  \begin{pmatrix}6&5\\5&0\end{pmatrix}
  \]
  gives 6 matched cells greedily but 10 under the optimal assignment. It does recognize an exactly identical partition under relabeling; the defect concerns the reported discrepancy count when partitions differ.

  **Concrete fix:** use maximum-weight bipartite assignment. Treat sentinel values such as “unattached” and “unidentified” separately: they are statuses, not ordinary family labels that can be permuted onto a known family.

- **[ROBUSTNESS]** `compare_phase_catalogs.m:97–108,145–148` — Nonfinite checking covers only common cells. `noWorse` can be true while a new-only entry is malformed. It also says nothing about certification of the compared roots.

  **Concrete fix:** validate every claimed occupied cell separately. Define `noWorse` explicitly as a comparison of valid recorded objective values within tolerance; require a content-bound audit before interpreting it as an accepted library improvement.

### Audit association

- **[CORRECTNESS]** `reproduce_library_70mN.m:204–215,245–272` — Exactly one audit in the last-round directory, with the right entry count, does **not** establish that it audited the final catalog. A stale audit of a different 576-entry catalog passes this association test.

  **Concrete fix:** store the audited catalog’s content hash and problem fingerprint in the audit, and verify them against the final catalog. Also validate that the audit covers the exact entry keys once each, rather than trusting only `nOk+nBad` and catalog metadata.

  Requiring exactly one pattern match is appropriately fail-closed against ambiguity; it is not a substitute for content identity.

- **[ROBUSTNESS]** `reproduce_library_70mN.m:43–44,162–168` — The advertised resumability conflicts with `.adoptArcs=false`: a second invocation rejects the arc files produced by the first invocation.

  **Concrete fix:** distinguish starting a fresh re-walk from resuming an initialized campaign with the same specification. Reject pre-existing untracked arcs at initialization, but permit verified campaign-owned arcs on resume.

### Interpretation of the full rebuild

The supplied outcome is strong practical evidence:

- complete coverage in both catalogs;
- 565 roots reproduced within the stated tolerances;
- 11 different roots, all substantially faster;
- no reported audited bad entries.

It is **not exact reproduction**: improvements of 1.95–3.93 days are not certification-noise ties. The 82 same-root family-label changes also show that provenance labels should be reported separately from mathematical solution agreement.

`R.ok=false` and `R.noWorse=true` are therefore sensible outcomes under the current definitions. I would expose separate results for **same solutions**, **same provenance partition**, and **audited no-regression/improvement**, rather than changing “reproduced” silently to mean any of those.

---

## Ranked top five

1. **Repair the continuous-time H2/H3 claim.** Endpoint slope samples produce an estimate, not a lower-bound certificate. Implement interval adjoint bounds with error enclosures, or downgrade the claim.
2. **Make arc identity explicit and fail closed.** No NaN/time-only fallback for suppressing discovery; preserve signed chart information; detect ambiguous matches; confirm candidates by polishing.
3. **Stop treating flight-time residual components as costate branches.** They are useful heuristics, not safe-interpolation or branch-identity certificates.
4. **Close or clearly retain the mathematical sufficiency gap.** The constrained-throttle reduction, every-subinterval regularity, and H5 remain conditional/assessed; the existing prose must not call them completed proofs.
5. **Bind the final audit to the actual final catalog.** Counts and filename patterns cannot establish that the accepted 576 entries are the ones audited.