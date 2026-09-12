## 1. Mathematics of the conditions

**Bottom line:** the Hamiltonian, mass-costate sign, switching-function sign, and spherical Legendre formula are correct for the stated problem. The quotient used by S4 is also defensible **for this all-burn reduction**, provided its additional assumptions hold. The main defects are the effective V1 gate, the incomplete coverage and numerical interpretation of S4, the calibration of S3, and overstated independence between checks.

### Reference calculation

Write
\[
q=(r,v),\qquad
C=\begin{pmatrix}0&2&0\\-2&0&0\\0&0&0\end{pmatrix},
\]
and let \(\Omega\) be the usual rotating CR3BP potential. The dynamics are
\[
\dot r=v,\quad
\dot v=\nabla\Omega(r)+Cv+\frac{uT}{m}\alpha,\quad
\dot m=-\frac{uT}{c}.
\]

Thus, with the normal multiplier normalized to one,
\[
H=1+\lambda_r^\top v+
\lambda_v^\top(\nabla\Omega+Cv)
+uT\left(\frac{\lambda_v^\top\alpha}{m}-\frac{\lambda_m}{c}\right).
\]

For \(\rho=\|\lambda_v\|>0\),
\[
\alpha_*=-\frac{\lambda_v}{\rho},\qquad
Q=\frac{\rho}{m}+\frac{\lambda_m}{c},
\]
so
\[
H(\alpha_*,u)=1+\lambda_r^\top v+
\lambda_v^\top(\nabla\Omega+Cv)-uTQ.
\]
Consequently, \(u=1\) minimizes \(H\) iff \(Q\ge0\); it is the unique minimizing throttle when \(Q>0\).

The adjoint equations, differentiating at fixed controls, are
\[
\begin{aligned}
\dot\lambda_r&=-\Omega_{rr}(r)\lambda_v,\\
\dot\lambda_v&=-\lambda_r-C^\top\lambda_v
             =-\lambda_r+C\lambda_v,\\
\dot\lambda_m&=\frac{uT}{m^2}\lambda_v^\top\alpha
             =-\frac{uT\rho}{m^2}\quad\text{at }\alpha_*.
\end{aligned}
\]

Because terminal mass is free, with no terminal mass cost or active terminal mass constraint,
\[
\lambda_m(t_f)=0.
\]
Its appearance in the dynamics does **not** change that boundary condition; it makes \(\lambda_m\) generally nonzero before arrival.

Because the problem is autonomous, the endpoint constraints have no explicit terminal-time dependence, and final time is free,
\[
H(t_f)=0,\qquad H(t)\equiv0,
\]
not merely “\(H\) is constant.”

For a unit tangent vector to \(S^2\) at \(\alpha_*\), the constrained second derivative is \(uT\rho/m\). At full thrust the spherical Hessian is therefore
\[
H_{\alpha\alpha}\big|_{T_{\alpha_*}S^2}
=\frac{T\rho}{m}I_2.
\]
This is the **constrained/Riemannian Hessian**, not the restriction of the zero ambient Hessian of a linear function.

The nondimensionalizations at `transfer_study.m:152–154` are consistent with these equations. Likewise, `fly_transfer.m:80` computes accumulated thrust delta-v,
\[
\int_0^{t_f}\frac{T}{m(t)}\,dt=c\log\frac1{m_f},
\]
not a difference of endpoint velocities.

### Findings on N2, N4, N5 and N6

- **[GAP]** `costate_common/pmp_pointwise_checks.m:60–62` — **N2 evaluates the right structural expression, including mass dynamics, but does not independently establish that the field is the stated physical field.** `1 + lambda'*F(1:7)` includes \(\lambda_m\dot m\), so there is no missing mass term here. However, both the Hamiltonian and the flight use `tfMinEoM`; a self-consistent wrong force or mass-flow model can satisfy this check. The actual `tfMinEoM` source is not supplied, so its state and adjoint formulas cannot be symbolically verified from this wrapper. **Fix:** compare its state rows against an explicit reference implementation of the equations above, including the Coriolis terms and \(\dot m=-uT/c\), and evaluate the corresponding reference Hamiltonian.

- **[GAP]** `costate_common/pmp_pointwise_checks.m:80–88,100–109` — **N5 is a field-consistency test, not a test of the derivative of the returned costate trajectory.** It compares the costate rows of `rhs` with a numerical derivative of the state rows of the *same* `rhs`. That is useful, but it neither independently validates the physical state equations nor bounds the ODE defect of `Y`.

  There is also an important qualification to the differentiation: the routine fixes the costate, but recomputes the feedback control as the state is perturbed. On the intended strict all-burn branch this is correct: \(\alpha_*=-\lambda_v/\rho\) is state-independent and \(u=1\) is locally constant. More generally, differentiation through an arbitrary erroneous feedback is not differentiation of \(H\) at fixed applied control; it can add \(H_a a_x\) terms.

  **Fix:** label this “adjoint-field consistency”; compare against the analytic adjoint equations above, with the applied control frozen; separately bound integration/continuous-output defects on the flown arc. Do not discard the latter requirement merely because naive differencing of sparse output is inaccurate.

- **[GAP]** `costate_common/pmp_pointwise_checks.m:90–95` — **The advertised direction-gap identity is conditional on actual unit throttle; the coupled acceleration/mass-flow control is not fully audited.** The recovered vector is
  \[
  b=\frac mT(F_v-F_{0,v})=u\alpha,
  \]
  not necessarily a unit direction. If \(u=1\), the computed expression is exactly
  \[
  H(1,\alpha)-\min_{\beta\in S^2}H(1,\beta)
  =\frac Tm(\lambda_v^\top\alpha+\rho).
  \]
  For arbitrary applied \(u\), the direction-only gap at that same throttle is instead
  \[
  \frac{uT}{m}(\lambda_v^\top\alpha+\rho).
  \]
  The gap against the **full** control minimum is
  \[
  H(u,\alpha)-\min_{\tilde u,\beta}H
  =\frac{uT}{m}(\lambda_v^\top\alpha+\rho)
   +T\{\max(Q,0)-uQ\}.
  \]
  When \(Q\ge0\), this becomes
  \[
  \frac{uT}{m}(\lambda_v^\top\alpha+\rho)+(1-u)TQ.
  \]

  The acceleration-norm check establishes only the acceleration-side throttle. It does not check that the field consumes mass at that same throttle. **Fix:** recover both
  \[
  u_{\rm acc}=\frac mT\|F_v-F_{0,v}\|,
  \qquad
  u_{\rm mass}=-\frac cT F_m,
  \]
  verify agreement and admissibility, and evaluate the complete gap above. With verified \(u=1\), the existing direction formula is correct.

- **[CORRECTNESS]** `DRO_tulip/indirect/transfer_study.m:250–251` — **“The first variation vanishes” is not the right description of this boundary-control problem.** At strict full thrust,
  \[
  H_u=-TQ<0.
  \]
  An inward throttle variation has \(\delta u\le0\), hence \(H_u\delta u\ge0\), generally not zero. The first-order control condition is a variational inequality/minimum principle. A shooting residual also does not establish that inequality by itself. **Fix:** say “the PMP shooting equations are satisfied to tolerance”; retain N6 as the separate control-minimization condition.

### S1 and S2

The quantities named by S1 and S2 are correct under \(T>0\), \(m>0\), and the applied minimizing direction. There is, however, a useful dependence between them and the other conditions:

\[
\lambda_m(t)
=\int_t^{t_f}\frac{u(s)T\rho(s)}{m(s)^2}\,ds\ge0.
\]
Therefore, if \(\rho(t)>0\), then \(Q(t)>0\). For this particular free-terminal-mass problem, strict bang follows from the exact adjoint equation, terminal mass transversality, admissibility, and S1. It is not an independent physical hypothesis once those facts hold.

This does **not** make its numerical evaluation worthless: it can expose violations of those facts.

### S3: what the lift-space rank establishes

- **[GAP]** `costate_common/mintime_hypothesis_gates.m:14–21,109–121` — **The constructed space is a space of stationary lifts; “no abnormal PMP lift iff \(\dim S=1\)” needs qualification.** The frozen-control Jacobian is the correct Jacobian for this construction. The constraints
  \[
  [\alpha(t)]_\times\lambda_v(t)=0,\qquad \lambda_m(t_f)=0
  \]
  correctly express angular stationarity and free-mass transversality.

  Also, for such a lift,
  \[
  \frac d{dt}(\lambda^\top f)
  =\frac Tm\lambda_v^\top\dot\alpha=0,
  \]
  because \(\lambda_v\parallel\alpha\) and \(\alpha^\top\dot\alpha=0\). Thus \(\lambda^\top f\) really is a conserved linear functional on this space.

  If an exact normal lift lies in \(S\), then \(\dim S=1\) does exclude any nonzero lift in the kernel of that functional. But if \(\dim S>1\), a kernel vector need not satisfy the PMP **minimization inequalities**: parallelism does not impose the required orientation of \(\lambda_v\), and the construction imposes no throttle inequality on alternative lifts. The linear-space statement is a corank-one/uniqueness condition, potentially stronger than absence of an abnormal *minimizing* lift.

  **Fix:** name the gate “corank-one stationary lift space,” state precisely which notion of abnormality the theorem uses, and retain \(\dim S=1\) as the sufficient exclusion condition rather than asserting the unqualified converse.

- **[GAP]** `costate_common/mintime_hypothesis_gates.m:97–121` — **A whole-arc rank calculation is not, without an argument, the strong-normality hypothesis on subarcs used in regular conjugate-point theorems.** In a general smooth problem, extra lifts can exist on an initial subarc and cease to satisfy stationarity later.

  There is a plausible mathematical repair here without testing infinitely many prefixes: away from collision, \(m=0\), and \(\rho=0\), the strict all-burn CR3BP extremal and its frozen-control lift constraints are analytic. Vanishing of an analytic stationarity constraint on an open subinterval extends along the connected arc. The mass costate is then determined by its terminal condition. This can establish the relevant subarc corank statement from the full-arc one.

  **Fix:** provide that analytic-continuation argument, including the required domain assumptions, or verify the precise subarc condition required by the cited theorem. Do not substitute a terminal rank label for it without explanation.

### S4 and V1: the quotient calculation

The quotient is **not automatically wrong because final time is free**. For this problem its algebra can be checked explicitly.

On a strict all-burn arc,
\[
m(t)=1-\frac Tc t.
\]
The six-state reduced Hamiltonian is
\[
h(t,q,p)
=p_r^\top v+p_v^\top(\nabla\Omega+Cv)
-\frac{T}{m(t)}\|p_v\|.
\]
It is time-dependent and homogeneous of degree one in \(p=(\lambda_r,\lambda_v)\).

Let
\[
J(t)=D_{p(0)}q(t).
\]
Homogeneity gives
\[
J(t)p(0)=0,\qquad p(t)^\top J(t)=0.
\]
Thus restricting initial costates to a complement \(P\) of \(p(0)\) correctly removes the scaling kernel.

The endpoint map with variable arrival time has derivative
\[
M(t)=[J(t)P,\ f_q(t)].
\]
The last column is the actual endpoint variation induced by varying time. It is **not** a zero-cost or endpoint-preserving gauge variation.

If \(E(t)\) is an orthonormal frame for \(p(t)^\perp\), then, up to orientation,
\[
\det M(t)
=\frac{p(t)^\top f_q(t)}{\|p(t)\|}
 \det\!\left(E(t)^\top J(t)P\right).
\]
Therefore, when \(h(t)=p(t)^\top f_q(t)\ne0\),
\[
\det M(t)=0
\quad\Longleftrightarrow\quad
\operatorname{rank}(J(t)P)<5.
\]

This establishes the algebraic role of H6. From the full normalized Hamiltonian,
\[
h(t)=-1+\frac Tc\lambda_m(t).
\]
Since \(\lambda_m\) decreases to zero, \(\lambda_m(0)<c/T\) ensures \(h<0\) throughout the arc.

The integral identity in `h6_margin.m:23` is also correct:
\[
\frac d{dt}(m\lambda_m)=-TQ,\qquad
\int_0^{t_f}TQ\,dt=\lambda_m(0),
\]
using \(m(0)=1\) and \(\lambda_m(t_f)=0\).

- **[GAP]** `costate_common/ms_conjugate_test.m:10–32,133–174` — **The determinant has the claimed reduced rank interpretation only under the reduction and regularity assumptions above; it is not a generic free-time Jacobi instrument.** Excluding \(\lambda_m(0)\) is legitimate locally on the strict all-burn branch because it has no effect on \(q(t)\). Removing the \(p(0)\) ray is legitimate before normalizing the lift. But the resulting six-state system is nonautonomous, and the endpoint mass is free; an autonomous fixed-point theorem cannot simply be imported by dimension counting.

  **Fix:** document the reduction of the free-mass second variation/critical cone to this six-state problem, and verify the two identities \(Jp(0)=0\) and \(p^\top J=0\) numerically as implementation diagnostics. State that the flow column represents variation of arrival time, not an endpoint-preserving “time reparameterization” degeneracy.

- **[CORRECTNESS]** `DRO_tulip/indirect/transfer_study.m:368–390` — **The effective V1 gate does not enforce the condition its own helper computes.** `h6_margin` requires both
  \[
  \text{margin}>1,\qquad \text{clearance}>H_{\rm resid}.
  \]
  The script instead accepts `h6Margin >= 1` and ignores `gates.h6Ok`, `h6Clearance`, and `h6Hmax`. Thus equality passes the script’s declared strict condition, and a negative reduced-Hamiltonian margin smaller than the measured Hamiltonian uncertainty also passes.

  **Fix:** use `gates.h6Ok` as the gate, retain strict comparisons, and require an explicit positive clearance exceeding the complete uncertainty budget. If a larger policy margin is requested, pass it into `h6_margin` rather than reimplementing its verdict.

- **[GAP]** `costate_common/h6_margin.m:13–21,59–85` — **H6 is correct only for the specified normal normalization and established monotonicity; it is not an invariant inequality on arbitrary costate representatives.** With objective multiplier \(p_0\),
  \[
  h=-p_0+\frac Tc\lambda_m,\qquad
  \lambda_m(0)<\frac cT p_0.
  \]
  The threshold \(c/T\) is correct because this solver fixes \(p_0=1\). Its use also relies on the correct minimizing direction, positive mass, and terminal transversality. A sampled Hamiltonian residual is not itself a rigorous error bound.

  **Fix:** state those dependencies and normalize before applying H6. Alternatively, formulate a correctly derived projected/mixed-boundary Jacobi test that does not multiply its intrinsic determinant by \(h(t)\); then rederive its validity conditions rather than carrying this particular spurious-zero mechanism into the instrument.

- **[GAP]** `DRO_tulip/indirect/transfer_study.m:316–390` — **The theorem application needs an explicit regular-extremal statement, not just the four gate names.** The required ingredients include a smooth collision-free positive-mass domain, the appropriate corank/normality condition, a correct free-mass/free-time second variation, an initial interval on which the second variation has the minimizing index, and absence of conjugacy through the terminal time.

  The strong-control part is available here: for every admissible control,
  \[
  H(u,\alpha)-H(1,\alpha_*)
  =\frac{uT\rho}{2m}\|\alpha-\alpha_*\|^2+(1-u)TQ
  \]
  when \(Q>0\). Uniform positive margins give a strengthened global Weierstrass inequality, not merely a local angular Hessian test.

  **Fix:** state the precise sufficiency theorem and show how each of these requirements follows from the exact problem facts or from a validated check. There is **no missing phase-transversality condition** for the fixed-phase problem actually solved.

## 2. The verdict sentence

The distinction is important:

- As a **conditional theorem statement about an exact regular extremal**, strict strong local minimality is a defensible conclusion after establishing the reduction and regularity conditions above.
- As a **verdict supported by these executed checks**, it is too strong: the program has not established S4, and its effective V1 gate is wrong.
- Merely replacing “strong” by “weak” would not repair that evidentiary gap.

- **[GAP]** `DRO_tulip/indirect/transfer_study.m:399–409` — **The numerical PASS branch is not a certificate of the hypotheses in its conditional sentence.** In particular, a sampled determinant with an omitted initial interval is not “no conjugate time in \((0,t_f]\).” The subsequent disclaimer is valuable, but does not make that numerical predicate equivalent to S4.

  **Fix:** either add validated root/trajectory enclosures, positive lower bounds for the strict margins, a certified rank separation, and a complete Jacobi-index/nonvanishing argument; or print something such as:

  > The numerical diagnostics are consistent with a regular normal extremal and with strict strong local minimality for fixed departure and arrival position/velocity, initial mass fraction one, and free terminal mass and time. No determinant zero was detected on the tested junction grid. This is numerical evidence, not a certificate.

  The exact conditional theorem can follow separately, with its assumptions explicitly listed.

- **[GAP]** `DRO_tulip/indirect/transfer_study.m:402–403` — **“Same endpoints” should not silently fix terminal mass, and “same endpoints and phases” is redundant here.** The actual problem fixes initial \(r,v,m\) and final \(r,v\), while leaving final mass free. It does not optimize departure or arrival phase.

  **Fix:** state those endpoint conditions explicitly. If the claim is later widened to free phases, add
  \[
  \lambda_{rv}(0)^\top \frac{dq_D}{ds_D}=0,\qquad
  \lambda_{rv}(t_f)^\top \frac{dq_A}{ds_A}=0,
  \]
  and the corresponding endpoint-manifold second-variation/focal conditions. The present point-to-point determinant cannot certify that widened problem unchanged.

## 3. Sampling and numerical failure modes

- **[GAP]** `costate_common/ms_conjugate_test.m:196–204` — **The beginning of the interval is discarded without proving it conjugate-free.** Choosing the first sample whose equilibrated singular-value ratio exceeds a threshold does not establish what happened earlier. The determinant starts singular at \(t=0\); it can become full rank, encounter a conjugate point, and become full rank again before the first accepted sample. Even when `firstFullRank == 1`, the interval from zero to that first junction remains untested.

  Moreover, the stated “initial coast” explanation does not apply to this certified strict all-burn branch. Poor small-time scaling is not evidence of a coast.

  **Fix:** establish a short-time positivity/nonconjugacy interval analytically or with validated bounds, and overlap that interval with the numerical continuation. Otherwise report the initial interval as uncovered and block S4.

- **[GAP]** `costate_common/ms_conjugate_test.m:164–165,244–255`; `DRO_tulip/indirect/transfer_study.m:389–390` — **Incomplete final-time coverage can return PASS.** Without `Yend`, the instrument stops at \(t_K<t_f\), but can still return `PASS`; the caller merely prints the uncovered interval and does not gate on it.

  **Fix:** require coverage through the intended \(t_f\), within a justified time tolerance, before S4 can pass. Return `UNDETERMINED` for incomplete coverage. For this caller, make the terminal state and all segment STMs mandatory.

- **[GAP]** `costate_common/ms_conjugate_test.m:183–247` — **Nonzero sampled determinant signs are not a nonconjugacy test, and near-singularity after the initial skip is not acted upon.** An even-multiplicity zero, two crossings within one segment, or an endpoint zero perturbed to a tiny same-sign number can pass. `rankTol` is used to select the beginning, not to flag unresolved near-rank-loss later. With `zeroTol=0`, an exactly conjugate endpoint generally need not produce an exactly zero floating-point determinant.

  **Fix:** treat any singular value comparable to its error bound as unresolved; refine using a continuous Jacobi representation. For certification, use interval nonvanishing/rank bounds or a validated Riccati/Maslov/Jacobi-index method. Finer sign sampling alone still does not exclude touches or close pairs.

- **[ROBUSTNESS]** `costate_common/ms_conjugate_test.m:225–237`; `DRO_tulip/indirect/transfer_study.m:353–362` — **Some status classifications are mathematically or operationally wrong.** Opposite nonzero signs at the last two sample times imply a zero **strictly inside** that bracket; it cannot be only at \(t_f\), whose sampled determinant is nonzero. Calling every last-bracket crossing “ENDPOINT” is unnecessarily inconclusive. Conversely, `UNDETERMINED` becomes the script’s default `FAIL`, and `firstFullRank = nS+1` can make the coverage print index past `cj.t`.

  **Fix:** distinguish a certified interior sign bracket from an unresolved terminal zero; propagate `UNDETERMINED` as unresolved; handle “never full rank” without indexing a nonexistent sample.

- **[GAP]** `costate_common/lift_space_dim.m:46–49`; `costate_common/mintime_hypothesis_gates.m:120–121` — **The rank threshold uses a one-vector residual as though it bounded the matrix error. It does not.** For the computed matrix,
  \[
  \sigma_{\min}(C)\le
  \frac{\|C\lambda_0\|}{\|\lambda_0\|}
  =\texttt{nullResid}.
  \]
  Therefore, absent the cap, `tol >= 10*nullResid` guarantees that at least one singular value is declared null. This incorporates the desired known-null direction into the threshold.

  More seriously, a small residual on \(\lambda_0\) gives no upper bound on errors in directions orthogonal to \(\lambda_0\). Such errors can destroy a second true null direction and create a false numerical `dimS == 1`.

  **Fix:** bound \(\|\widehat C-C\|\), including flight, interpolation, frozen-adjoint integration, and rounding errors. With an established exact normal lift, a bound
  \[
  \sigma_6(\widehat C)>\|\widehat C-C\|
  \]
  certifies rank at least six; the known null lift then makes the dimension exactly one. Refinement differences are useful evidence, but are not automatically rigorous error bounds.

  Importantly, **finite sampling is not itself a false-pass problem for exact lift-space rank**: a true continuous abnormal stationary lift would satisfy every sampled constraint. Exact sampled rank six, together with the known normal lift, already excludes it. The dangerous issue here is matrix/trajectory error and threshold selection.

- **[GAP]** `DRO_tulip/indirect/transfer_study.m:345–346,389–390` — **S3’s self-consistency diagnostics are printed but not required.** The script accepts `dimS == 1` even if the accepted lift has a bad `nullResid` or the independently assembled fixed field has a bad `Hresid`. The adaptive rank threshold can partially conceal the former.

  **Fix:** make finite, quantitatively acceptable lift and Hamiltonian residuals prerequisites for interpreting the rank. Require a resolved singular-value separation, not just a nullity produced by a selected tolerance.

- **[GAP]** `costate_common/ms_conjugate_test.m:168–173`; `DRO_tulip/indirect/transfer_study.m:246,326,352` — **S4 is not computed on the “one flight” used by N2–N6, and the variational equations are not independently audited.** It chains segment STMs about multiple-shooting junction states, whereas the pointwise checks use a single propagation and the gates make another propagation. A product of segment derivatives about slightly discontinuous junction states is not exactly the STM of the single flown trajectory.

  Also, N5 does not test the STM generator. S3 correctly needs the **frozen-control** Jacobian; S4 needs the derivative of the **optimized Hamiltonian flow**, including
  \[
  D_{\lambda_v}\alpha_*
  =-\frac1\rho\left(I-\frac{\lambda_v\lambda_v^\top}{\rho^2}\right).
  \]
  Confusing those two derivatives can invalidate S4 while leaving N5 satisfactory.

  **Fix:** verify junction-to-flight agreement and propagate/audit the variational system along the same reference arc. Test the STM against independent directional derivatives and Hamiltonian identities, including symplecticity in the full canonical formulation, the scaling kernel, and \(p^\top J=0\). Include chaining/conditioning error in the conjugate test.

- **[GAP]** `costate_common/pmp_pointwise_checks.m:58–75`; `costate_common/mintime_hypothesis_gates.m:79–86`; `costate_common/validate_flight.m:69–81` — **The continuous inequalities and equations are represented by samples with no between-sample bounds.** Specifically:
  - N2 uses all returned samples, not the continuum.
  - N5 and the applied-control part of N6 use at most 120 interior indices, excluding both endpoints.
  - S1 and S2 use sample minima.
  - Primary clearance uses sample minima.
  - N4 is an endpoint condition, so its uncertainty is endpoint propagation error rather than an unsampled interior extremum.

  A zero of \(\|\lambda_v\|\) is particularly easy to miss: its norm can be positive at both neighboring samples, while the normalized direction becomes singular between them. A close primary passage can likewise occur between safe samples.

  **Fix:** locate and bound extrema using dense output with error control, derivative bounds, or interval integration. Require margins exceeding those errors. For clearance, bracket closest approaches using
  \[
  (r-r_i)^\top v=0
  \]
  and evaluate their enclosed distances. Arbitrarily dense output is not a validated continuous bound.

- **[GAP]** `DRO_tulip/indirect/transfer_study.m:85–87,255–262,303–305` — **Neither the small internal residual nor the permitted flown miss establishes an exact nearby fixed-endpoint extremal.** The script allows a single-shot arrival error of 100 km and 10 m/s. Those are practical screening tolerances, not fixed-endpoint satisfaction for a local-minimum theorem. The observation that one integration mode returns a smaller residual does not establish that its residual is closer to the exact shooting residual.

  **Fix:** estimate propagation error and shooting conditioning; report both internal matching residuals and independently propagated boundary residuals. For a certificate, use a validated shooting/root enclosure or an appropriate Newton–Kantorovich/Krawczyk argument, then carry its uncertainty into the hypothesis checks. Tightening a residual alone does not address an ill-conditioned BVP.

- **[ROBUSTNESS]** `costate_common/pmp_pointwise_checks.m:77,93–95` — **Negative direction gaps are silently clipped to zero.** `dirGap` starts at zero and accumulates only a maximum. A slightly over-unit recovered thrust aligned with the minimizing direction produces a negative value; the looser `throttleErr < 1e-10` condition can permit it while the nominal gap tolerance is \(10^{-12}\).

  **Fix:** check both the minimum and maximum gap against a rounding/error budget, and make the throttle feasibility tolerance consistent with the induced Hamiltonian error. Once unit direction is independently established, the nonnegative form
  \[
  \frac{T\rho}{2m}\|\alpha-\alpha_*\|^2
  \]
  avoids cancellation near the minimizer.

- **[ROBUSTNESS]** `costate_common/validate_flight.m:60–81` — **The validator checks only state finiteness and the terminal mass law, not full augmented-flight admissibility.** It does not enforce real finite costates, positive mass at every sample, or the mass law throughout the trajectory. Nor does it explicitly require the expected all-burn terminal mass \(1-Tt_f/c\) to be positive; near exhaustion, its absolute mass tolerance can admit a positive reported mass with a slightly nonpositive expected mass.

  **Fix:** require real finite times and augmented states, consistent dimensions and time ordering, \(1-Tt_f/c\) bounded strictly above zero, and agreement with \(m(t)=1-Tt/c\) over the flight. If the clearance radii are actual path constraints rather than screening values, establish strict inactive clearance continuously or use the corresponding constrained PMP.

## 4. Tautologies, redundancies and common-mode checks

Not every redundant check should be removed. Many are valuable implementation or integration checks. They must not, however, be counted as independent mathematical evidence.

- **[GAP]** `DRO_tulip/indirect/transfer_study.m:265–277`; `costate_common/ms_tfmin.m:103–105` — **N4 is already a shooting equation, and N2 follows exactly from the terminal Hamiltonian equation plus autonomous canonical propagation.** On the exact same correctly propagated BVP solution, neither supplies a new independent mathematical restriction. Re-evaluating them on another propagation can expose numerical defects, which is worthwhile.

  **Fix:** describe them as independent evaluations of endpoint/transversality and Hamiltonian preservation/normalization, not as independent reasons that an extremal is minimizing.

- **[GAP]** `DRO_tulip/indirect/transfer_study.m:334–338` — **S2 is mathematically redundant for an exact admissible lift satisfying the minimizing direction, N4, N5 and S1.** The integral formula for \(\lambda_m(t)\) above forces \(\lambda_m\ge0\), hence \(Q>0\). It cannot independently fail while those exact premises hold.

  **Fix:** retain it as a valuable consistency and strict-margin diagnostic, but state the dependency. A failure means at least one of those premises is not established; it is not an additional independent physical mechanism.

- **[GAP]** `costate_common/lift_space_dim.m:47–48` — **The existence of at least one numerically declared null direction is often forced by the threshold construction.** The inequality \(\sigma_{\min}(C)\le\texttt{nullResid}\) makes this explicit. The complete `dimS == 1` test is not a tautology—the second-smallest singular value can also fall below the threshold—but the “we found the expected null direction” part is largely circular.

  **Fix:** report the known-lift residual separately from the independently resolved rank margin. Do not present threshold-induced nullity as confirmation of the assumed lift.

- **[GAP]** `DRO_tulip/indirect/verify_with_pumpkyn.m:78–105`; `costate_common/mintime_rhs_point.m:22`; `costate_common/ms_tfmin.m:103–105` — **X1 is a second shooting implementation, not an independent physics/PMP implementation.** Both solvers use `pumpkyn` propagation and equations. A common dynamics, Hamiltonian, mass-flow, or STM error can survive both. This is not the earlier copied-solver tautology, but the claim that it protects against a self-consistent solution of the wrong physical problem is overstated.

  Also, if the foreign solver returns its input on this particular call, its subsequent flight repeats the already successful flight. Perturbation tests showing that it moves other inputs do not establish convergence on this input.

  **Fix:** call it a cross-check of the nonlinear shooting implementation using shared dynamics. For independent physical verification, implement the state dynamics, adjoints and boundary residual independently.

- **[GAP]** `DRO_tulip/indirect/verify_with_pumpkyn.m:68–71,91–105` — **Reaching the target does not establish that the foreign solver solved the PMP BVP.** On a strict all-burn arc, changing \(\lambda_m(0)\) does not change \(r,v,m\); scaling the six position/velocity costates also leaves the direction and flight unchanged. Those changes can preserve perfect arrival while violating \(\lambda_m(t_f)=0\) and \(H(t_f)=0\).

  **Fix:** evaluate the foreign answer’s complete residual, including both terminal multiplier conditions and the control conditions, rather than inferring “solved” from arrival. Agreement with the original answer remains a useful but separate check.

- **[GAP]** `DRO_tulip/indirect/transfer_study.m:376–385`; `costate_common/pmp_pointwise_checks.m:70–72` — **V2 does not detect all of the wiring mistakes it claims to detect.** In particular, `PW.minQmt` contains no explicit `Tnd` or `mu` dependence: it is computed from the supplied `Y` and `cnd`. Passing the wrong thrust or mass ratio to `pmp_pointwise_checks` can change its Hamiltonian/adjoint results without changing this V2 comparison at all. Different flights can also share the same scalar minimum.

  **Fix:** compare actual provenance—initial state, costates, time interval, parameters, and trajectory agreement—not one scalar statistic. Retain the minimum comparison only as an additional consistency diagnostic.

- **[GAP]** `DRO_tulip/indirect/transfer_study.m:156–177`; `costate_common/phase_state.m:55–64` — **A periodic spline’s seam derivative is a construction property, not evidence that the endpoint lies on a CR3BP periodic orbit.** A correctly built periodic cubic enforces that seam behavior even for dynamically wrong orbit data. Closure and periodic interpolation do not check the interpolant’s interior accuracy, the relation \(\dot r=v\), or the coasting equations. For fixed phases, a globally periodic \(C^1\) interpolant is not itself an optimality-theorem requirement.

  **Fix:** validate the orbit propagation and endpoint interpolation error independently, preferably against dense physical propagation to the selected phase. Keep the seam check as an interpolation implementation check, not an astrodynamical certificate.

The revised applied-direction N6 check is **not** the old sphere-search tautology: injecting a wrong applied direction can genuinely open its gap. Its remaining limitations are the throttle/mass coupling and numerical issues identified above.

## PER-CONDITION VERDICT

- **N2 — CORRECT:** \(H\equiv0\) is the correct normal free-time condition, and `1 + lambda'*F(1:7)` includes mass dynamics, subject to the underlying field being the stated one.
- **N4 — CORRECT:** Free terminal mass gives \(\lambda_m(t_f)=0\), and the minimizing-direction mass-costate rate is \(-uT\|\lambda_v\|/m^2\).
- **N5 — INCOMPLETE:** The calculation is valid on the strict all-burn branch as an adjoint-field consistency test, but does not independently establish the physical equations or the flown trajectory’s adjoint defect.
- **N6 — INCOMPLETE:** The direction gap and switching-function sign are correct at verified full thrust, but the applied acceleration/mass-flow coupling and the complete control gap are not fully checked.
- **S1 — CORRECT:** For positive thrust and mass, \(\|\lambda_v\|>0\) is exactly positivity of the spherical angular Hessian at the minimizing direction.
- **S2 — CORRECT:** \(Q=\|\lambda_v\|/m+\lambda_m/c>0\) is the correct strict-full-throttle condition, although it follows from the other exact PMP conditions and S1 here.
- **S3 — INCOMPLETE:** The stationary lift-space construction is appropriate, but its interpretation and numerical rank need a resolved error margin and the required subarc-normality argument.
- **S4 — INCOMPLETE:** The quotient determinant has the correct reduced rank interpretation when \(h\ne0\), but sampled signs, skipped initial coverage, possible final coverage loss, and unaudited STM error do not establish absence of conjugacy.
- **V1 — INCORRECT:** The normalized H6 formula is correct, but the script’s effective gate accepts a non-strict ratio and ignores the helper’s required Hamiltonian-clearance verdict.