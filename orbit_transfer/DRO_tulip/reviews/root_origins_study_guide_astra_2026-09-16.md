The main physical argument is sound, but the document is not yet reliable as a rigorous guide. The clearest equation errors are the missing time scaling in Hermite–Simpson and the missing `unwrap` in the turns formula. The clearest data error is the handoff box mixing smoke-run measurements into the default run. Several passages also promise more than the code checks.

I checked the supplied script, log, and FINDINGS 75–77. The bodies of the collocation solver, covector mapper, and certifier are not included here, so their internal implementation cannot be independently verified from this bundle. I have distinguished mathematical verification from agreement with the campaign record.

## What I checked and found RIGHT

| Item | Verdict and reason |
|---|---|
| Hamiltonian sign convention | **RIGHT** for the normal minimum-principle convention \(H=L+\lambda^\top f\), \(L=1\), with controls minimizing \(H\). The minimizing direction is \(-\lambda_v/\|\lambda_v\|\) when the primer is nonzero and thrust is applied. |
| Free terminal mass | **RIGHT for this objective and endpoint set.** The objective has no terminal-mass term, and the fixed terminal position and velocity constraints do not depend on mass. Thus \(\lambda_m(t_f)=0\), provided terminal mass is genuinely free and no terminal-mass bound is active. Fixed terminal position and velocity do **not** require their costates to vanish. |
| Mass adjoint and throttle coefficient | **RIGHT:** \(\dot\lambda_m=-T_{\rm nd}q\|\lambda_v\|/m^2\), and after direction minimization, \(\partial H/\partial q=-T_{\rm nd}(\|\lambda_v\|/m+\lambda_m/c_{\rm nd})\). The qualification needed to conclude all-burn is addressed below. |
| Sundman \(\lambda_t=+1\) | **RIGHT under the normal minimizing Mayer convention**, with objective equal to the terminal clock state. Autonomy makes \(\lambda_t\) constant. The document should explicitly show that formulation; otherwise it invites double-counting the time objective. |
| Multiple-shooting counts | **RIGHT:** \(7+14(K-1)+1=14(K-1)+8=14K-6\). For \(K=24\), there are **330 unknowns and 330 equations**. Square does not mean nonsingular or well-conditioned. |
| ND thrust and exhaust speed | **RIGHT.** In particular, \(c_{\rm nd}=I_{sp}g_0t^\star/(1000l^\star)\), which is exactly the document’s expression after substituting \(g_{0,\rm nd}\). The factor 1000 is correctly placed. |
| Falling propellant and ideal \(\Delta V\) | **RIGHT at fixed exhaust speed and initial mass.** Let \(a=T_{\rm nd}t_f/c_{\rm nd}\). Then propellant fraction is \(a\), and \(\Delta V_{\rm nd}=-c_{\rm nd}\ln(1-a)\), with \(d\Delta V_{\rm nd}/da=c_{\rm nd}/(1-a)>0\). Therefore decreasing \(Tt_f\) decreases **both** propellant and ideal \(\Delta V\), while \(0\le a<1\). |
| Next-seed construction | **RIGHT.** The code preserves the banked position, velocity, and costate guesses, restores the fixed initial state, scales the time grid, and rebuilds the mass row for the new engine. It does not propagate that seed or rescale its costates. The resulting seed need not satisfy the new dynamics before shooting. |
| Polished grid | **RIGHT:** \(t_k^{\rm cert}=\sigma_k C.z(8)\), because the normalized breakpoints are preserved while polishing changes \(t_f\). Keeping the old absolute grid would be inconsistent. |
| Junction-array distinction | **RIGHT:** the harvested seed is reported as \(14\times25\), while the shooting output stores \(14\times24\) starts and 25 times. These are different objects, not contradictory counts. |
| Turns interpretation | **RIGHT in the prose and code:** sampled total absolute Moon-centered azimuth change in the rotating frame, with reversals adding—not a winding number. The displayed equation is missing an essential operation, as noted below. |

## A. Mathematical and physical findings

- **[MATH]** `root_origins_study_guide.tex:331–340` (The cold lottery — Hermite–Simpson) — **The defect and midpoint formulas omit \(t_f\) under the stated normalized-time definition.** You define \(f\) as the physical-ND-time dynamics through \(dx/ds=t_f f\), but then use \(f_k\) directly with an unspecified \(h\). The Simpson coefficients and midpoint signs are right; their time scaling is not complete.

  Replace the formulation with:
  ```latex
  Let $h=1/N$ and define the normalized-time right-hand side
  \[
  \bm F_j=T_j\,\bm f(\bm x_j,\bm u_j,q_j).
  \]
  Then
  \[
  \bm x_{k+1}-\bm x_k
  -\frac{h}{6}\left(\bm F_k+4\bm F_{k+\frac12}+\bm F_{k+1}\right)=0,
  \]
  \[
  \bm x_{k+\frac12}
  =\frac{\bm x_k+\bm x_{k+1}}{2}
  +\frac{h}{8}\left(\bm F_k-\bm F_{k+1}\right).
  \]
  The lifted time satisfies $T_{k+1}=T_k$, with
  $T_{k+\frac12}=(T_k+T_{k+1})/2$; at feasibility every $T_j=t_f$.
  ```
  Equivalently, apply Hermite–Simpson to the augmented state \((x,T)\) with derivative \((Tf,0)\). Make the midpoint variables part of the local stencil description rather than saying the entire defect touches only two variables/nodes.

- **[MATH]** `root_origins_study_guide.tex:136–140` (The question — all-burn argument) — **The signs and terminal-mass step are right, but the bold conclusion needs its assumptions and its zero-primer qualification.** “Away from a zero primer” establishes pointwise strictness where \(\|\lambda_v\|>0\); by itself it does not explain the almost-everywhere all-burn conclusion.

  Add:
  ```latex
  For the normal, unconstrained PMP considered here, with positive mass,
  no terminal-mass penalty or active terminal-mass constraint, and endpoint
  states independent of $t_f$,
  \[
  \lambda_m(t)=\int_t^{t_f}
  \frac{T_{\rm nd}q(\xi)\|\lambda_v(\xi)\|}{m(\xi)^2}\,d\xi\ge0.
  \]
  The throttle coefficient is strictly negative whenever
  $\|\lambda_v\|/m+\lambda_m/c_{\rm nd}>0$.
  Thus a nonzero primer is sufficient, but not necessary, for $q=1$.
  ```
  Also explain the zero set. On a regular arc with smooth CR3BP coefficients,
  \[
  \dot\lambda_r=-g_r^\top\lambda_v,\qquad
  \dot\lambda_v=-\lambda_r-g_v^\top\lambda_v.
  \]
  At a primer zero, \(\dot\lambda_v=-\lambda_r\). Simultaneous vanishing of \(\lambda_r,\lambda_v\) would force those costates to vanish identically; terminal \(\lambda_m=0\) would then give \(H=1\), contradicting free-time transversality. Consequently primer zeros are isolated under these assumptions, and \(q=1\) almost everywhere. Direction is undefined at a primer zero, but isolated values do not affect the trajectory.

  State explicitly that an active clearance constraint requires the appropriate constrained-PMP treatment; this is the unconstrained argument used by the shooter.

- **[MATH]** `root_origins_study_guide.tex:561–568` (The handoff — Sundman and \(\lambda_t\)) — **The answer \(+1\) is correct, but the derivation omits the change from the original running-cost formulation to a Mayer clock-state formulation.** Carrying the original running cost into the augmented Hamiltonian as well would count time twice.

  Replace with:
  ```latex
  On $\tau\in[0,1]$, write the equivalent Mayer problem
  \[
  \min t(1),\qquad
  x'=S\kappa(x)f(x,u,q),\qquad t'=S\kappa(x),\qquad t(0)=0.
  \]
  With the normal minimizing convention and no running cost,
  \[
  \mathcal H=S\kappa(x)\bigl(\lambda_x^\top f+\lambda_t\bigr).
  \]
  Since the dynamics and endpoint constraints have no explicit clock
  dependence, $\lambda_t'=-\partial_t\mathcal H=0$.
  The Mayer terminal condition is
  $\lambda_t(1)=\partial t(1)/\partial t(1)=1$.
  Hence $\lambda_t\equiv1$ and $\mathcal H=S\kappa H$.
  ```
  Use primes for \(\tau\)-derivatives here. Under the corresponding sign-reversed maximum-principle convention, the clock costate would be \(-1\); \(+1\) is not convention-independent.

- **[MATH]** `root_origins_study_guide.tex:578–580` (The handoff — multiple shooting) — **The counts are right, but the full Jacobian is not simply block-bidiagonal, and short segments do not guarantee good conditioning.** The continuity portion has that local structure, but the free-\(t_f\) unknown contributes a border column, and the endpoint conditions contribute a terminal block.

  Replace:
  > The continuity Jacobian has block-bidiagonal junction structure, bordered by the free-time column and terminal boundary rows. Its STMs cover individual segments rather than the whole arc. This reduces long-shot sensitivity amplification, but the full boundary-value Jacobian can still be ill-conditioned or singular.

  Likewise, at line 339 replace “the linear solver dies at \(N=400\)” with a specifically sourced implementation observation, or remove it. A scalar \(t_f\) produces an arrow/border structure; failure at a particular mesh size is not a mathematical consequence of that structure.

- **[MATH]** `root_origins_study_guide.tex:677–683` (The indirect ladder — turns) — **The equation does not implement the code because it differences raw `atan2` values without unwrapping.** Crossing the branch cut can add a spurious nearly-full revolution.

  Replace with:
  ```latex
  \[
  \theta_i=\operatorname{unwrap}\!\left[
    \operatorname{atan2}\bigl(y_i,x_i-(1-\mu^\star)\bigr)
  \right],
  \qquad
  \mathrm{turns}=\frac{1}{2\pi}\sum_i|\theta_{i+1}-\theta_i|.
  \]
  ```
  Call this a **sampled estimate** of angular total variation. Unwrapping assumes adequate angular sampling; it cannot recover unobserved reversals or arbitrarily large intersample angular changes. The azimuth also requires nonzero projected Moon-centered radius.

- **[MATH]** `root_origins_study_guide.tex:507–534` (The handoff — what multipliers mean) — **“Those numbers are the costates” is too exact, and “still converges to a slightly wrong answer” wrongly attributes a permanent shooting error to seed misstationing.** Defect multipliers become approximate continuous costates through a transcription-specific sign, scale, and station mapping. A half-step error can damage the seed or bias a direct-costate comparison; a correctly formulated shooting solve that subsequently converges must still satisfy its own BVP.

  Replace the ELI5 core with:
  > Those bookkeeping numbers contain approximations to the steering costates. After applying the right sign, scale, and time-station mapping, we can use them to fill in a shooting seed.

  Replace the station sentence with:
  > In this transcription’s covector mapping, defect multipliers are associated with segment midpoints. Associating them with nodes misplaces the harvested costate samples by half an interval and can degrade the seed or its comparison with a continuous solution.

  Retain the **RIGHT** observation that a constant \(\lambda_t\) row cannot detect a station shift. The mapper’s midpoint implementation itself is only supported here by the record, since its body is not supplied.

- **[MATH]** `root_origins_study_guide.tex:600–605` (The handoff — measured interpretation) — **A node-radius margin does not prove the entire path constraint is inactive or make harvested multipliers exact unconstrained costates.** Exact inactive-constraint multipliers vanish by complementarity, but this run prints a sampled radius diagnostic from a finite-tolerance NLP.

  Replace the parenthesis with:
  > The minimum reported node radius is 4514 km above the floor. This is strong sampled evidence that the clearance constraint is inactive, not a continuous-clearance proof. Inactivity of all relevant path-constraint rows is what permits an unconstrained covector handoff; the subsequent shooting solve independently checks the unconstrained BVP.

  Also replace “The seed was already a root to the solver’s precision” with:
  > The handoff reached a final residual below the requested \(10^{-11}\) target in one reported iteration.

  A final residual plus an iteration count does not, without the solver’s counting convention or initial residual, prove that no correction was made.

- **[MATH]** `root_origins_study_guide.tex:823–829` (The basin hunt — segment lengths) — **“No segment is longer than \(1/K\) of the arc” requires a uniform physical-time grid; the helper does not enforce that.** It uses the supplied increasing `tGrid`.

  Replace with:
  > Segment \(k\) has duration \(t_f(\sigma_{k+1}-\sigma_k)\). On a uniform normalized-time grid this is \(t_f/K\). Each segment is propagated independently from its banked start.

  Make the same conditional distinction for the \(k/K\) statement at line 624. Also call `seamKm` the worst **position** seam mismatch, not a check of all 14 state-and-costate components.

## B. Document–code correspondence

- **[MISMATCH]** `root_origins_study_guide.tex:147–187` (Flow diagram) — **The diagram is useful as dataflow, but misleading as actual control flow.** In particular:
  - The lottery executes before the direct ladder; it is not a parallel branch.
  - A direct-ladder refusal stops that ladder, but the handoff uses the **last accepted rung**. Only no accepted rung is fatal.
  - A refused/thrown Sundman re-solve falls back to the plain-time solution; it does not prevent harvesting.
  - The indirect ladder can be disabled, or stop short, and certification still follows.
  - Hunt warm-start failure skips the probes and yields `INCONCLUSIVE`.

  The certification-failure arrow **is RIGHT**: retain the pre-certification candidate and hunt from it as an uncertified warm start. Outcome B being decided after certification is also **RIGHT**.

  Redraw the principal control sequence as:
  > §1 → optional §2 / outcome A → §3 → §4 → optional §5 → §6 / outcome B → optional §7 / outcome C → §8.

  Add the fallback/refusal paths above, or explicitly label the existing diagram “dataflow, not execution order.” Change the unconditional “Sundman re-solve → harvest” label to “try Sundman; otherwise plain-time fallback → harvest.”

- **[MISMATCH]** `root_origins_study_guide.tex:171,892–908` (Flow and self-check) — **The fatal checks do not wait until §8, and B1/H1 are not the only ways execution can abort.** B1 is asserted immediately after §3; H1 immediately after §4. Orbit validation, usable-dual validation, and malformed passing-certificate output also have fatal assertions.

  Replace:
  > Among the named final-table gates, B1 and H1 are required to continue and are asserted when those stages finish. Other validation assertions and unhandled helper errors can also abort the run. §8 summarizes the gates that execution reached.

  Replace “nothing above can be believed” with “there is no usable root for the downstream walk.” A failed handoff does not invalidate the already-recorded cold lottery.

  Add `NOT RUN` to the outcome vocabulary at line 904; the code uses it when experiments are disabled.

- **[MISMATCH]** `root_origins_study_guide.tex:151–154` (Flow — costates and persistence) — **Costates do not come into existence only at §4, and the checkpoint does not contain “everything it had.”** The direct NLP already has multipliers; §4 maps them into a shooting representation. §7 harvests candidate costates again. Moreover, cold/direct records retain diagnostics, not their full trajectories, and per-exponent attempts are saved only when the enclosing rung is recorded.

  Replace with:
  > §4 is the first conversion from discrete multipliers to a state-and-costate shooting seed; §7 may perform further harvests. The script checkpoints completed meshes, rungs, handoff/certification stages, and probes. An interrupted current stage may be absent, and cold/direct checkpoints contain diagnostics rather than complete solution arrays.

  Do not imply atomic or guaranteed persistence: `saveq` warns on failure but does not supply that guarantee.

- **[MISMATCH]** `root_origins_study_guide.tex:202–207` (Inputs — reference guard) — **The guarded comparisons are right, but “every ‘vs reference’ line … is silent” is not what the script prints.** For a nonmatching problem, `relM` remains NaN while the cold table still prints its `vs ref` column and percentage field.

  Replace with:
  > Reference comparisons are computed only when the live problem matches. Otherwise the script explains that the reference belongs elsewhere; the cold table retains its reference column with NaN values, and the summary and certified-root comparison suppress the numerical comparison.

- **[MISMATCH]** `root_origins_study_guide.tex:283–287` (Orbits — seam gate) — **The prose says the seam’s value and derivative mismatches are both gated at \(10^{-6}\); the code gates only `seam.deriv`.** Orbit endpoint closure is separately gated at \(10^{-7}\).

  Replace with:
  > The script requires both sampled orbit closures to be below \(10^{-7}\), and the maximum interpolant seam-derivative mismatch to be below \(10^{-6}\). It does not separately gate a seam-value diagnostic here.

- **[MISMATCH]** `root_origins_study_guide.tex:425–428,451–455` (Direct ladder — acceptance predicate) — **The stated acceptance contract is stronger than `directOK`.** It gathers node throttles if `U` is available and midpoint throttles if `Um` is available; it does not require both histories to be present. It also supplies a common **core** predicate, not identical total acceptance: the ladder adds clearance and band screens, the hunt adds clearance, and the default lottery is unconstrained.

  Replace with:
  > All direct solves use the same core predicate. It checks both lower and upper throttle deviations over the available node and midpoint histories. In the normal solver output both are available, but this predicate does not independently require both arrays. Each stage may add further screens.

  The actual throttle inequalities are strict: \(1-\epsilon<q<1+\epsilon\). Accordingly, at lines 267 and 932–934 change “cannot see a sustained 0.1% under-throttle” to “accepts a sustained deficit just below 0.1%; unsampled deviations may also escape detection.” Exactly 0.1% is the mathematical rejection boundary.

- **[MISMATCH]** `root_origins_study_guide.tex:818–819` (Basin hunt — banking) — **Not every converged probe is banked as a harvested seed.** It must first pass `directOK` and clearance, and harvesting can still fail. The later three-count explanation is right but contradicts this sentence.

  Replace with:
  > Accepted probes retain their direct solution arrays. Harvesting is attempted when usable multipliers are available; successfully harvested seeds are banked for later shooting and certification.

  Add that every probe uses the **same held source root**, not the result of the preceding phase probe.

- **[PEDAGOGY]** `root_origins_study_guide.tex:685–689,772–775` (Budgets and poolless policy) — **The budget explanation omits a limitation that materially changes execution.** Without a pool, `fenced` invokes the function directly and its external cap is ignored. Also, 600 s is not a strict aggregate rung wall-clock limit: a shooting fence allows an additional 90 s, and elapsed time—including any overrun—is charged when the next remaining-budget calculation occurs. Direct NLP solves use CPU caps, not these hard wall fences.

  Add:
  > These are requested work budgets plus external fences where a pool exists, not a universal run-time guarantee. Without a pool, the external shooting and reconstruction fences disappear; internal solver limits remain where implemented. Certification is then skipped unless explicitly allowed unfenced.

  The distinction between “timeout or worker error” and a known CPU-cap exit is **RIGHT** and should remain.

## C. Numerical and historical findings

- **[MISMATCH]** `root_origins_study_guide.tex:599–600` (Handoff — Measured) — **This mixes smoke-run and full-default-run measurements.** The supplied default log reports a relative time difference of \(3.80\times10^{-8}\), with both node radii rounding to 6414 km. The approximately \(6.1\times10^{-7}\) and 4 km figures belong to the reduced smoke runs in FINDINGS 76–77.

  Replace with:
  > Default run: the Sundman re-solve differs in \(t_f\) by \(3.80\times10^{-8}\) relative; both minimum node radii round to 6414 km. The log does not resolve their exact radius difference.

  Although the default difference is technically “within \(6\times10^{-7}\),” that loose bound conceals its smoke-run provenance; the 4 km statement is directly wrong for the displayed default values.

- **[MISMATCH]** `root_origins_study_guide.tex:411–413,957` (Direct ladder — prediction accuracy) — **“To 1%” is not supported by the default table.** The final step has guess 0.7416 ND and answer 0.7575 ND:
  \[
  0.7575/0.7416-1\approx2.14\%.
  \]
  FINDINGS 75’s “every rung to 1%” sentence is itself inconsistent with this log.

  Replace with:
  > The empirical scaling is typically accurate to about 1% in this run; its largest printed inter-rung error is about 2.1%, on the final 0.75-to-0.5 N step.

- **[MISMATCH]** `root_origins_study_guide.tex:715–718` (Indirect ladder — historical routes) — **“Four routes were tried on this cell” includes a route that started from another cell.** FINDINGS 75 lists three anchor-cell routes plus the campaign route from the sheet’s fastest 0.5 N cell.

  Replace with:
  > Three reported routes from the anchor cell stalled at 0.12, 0.105, and 0.11 N. A fourth reported campaign route, starting from the sheet’s fastest 0.5 N cell rather than this anchor cell, reached 0.09 N before stalling at 0.067 N.

  Avoid “only from” as an impossibility claim about untried routes.

- **[MISMATCH]** `root_origins_study_guide.tex:734–748` (Certification — review attribution) — **Forwarding was a round-1 finding, not one of two round-2 discoveries.** FINDINGS 76 explicitly records the options-forwarding fix. FINDINGS 77 records the erroneous polished-grid adoption and its correction.

  Either change the introduction to “Two interface details matter here,” or attribute forwarding to round 1 and grid consistency to round 2. The technical explanations of both interfaces are otherwise **RIGHT**.

### Remaining number audit: RIGHT or source-qualified

| Location | Checked result |
|---|---|
| Lines 119–120 | CR3BP constants agree with the script at the stated precision. |
| Lines 210–255 | Operating point, reference metadata, and every displayed tolerance match the script. |
| Lines 258–268 | Historical slack range, \(N=800\) first-run rejection, \(1.3\times10^{-6}\) slack, and revised \(10^{-3}\) threshold agree with the record. The exact 0.1% boundary wording needs the correction above. |
| Lines 370–392 | **All six flight times, all six reference percentages, all six radii, and statuses match their supplied sources.** The historical altitude-to-radius additions are correctly rounded. The 1394 km node is inside the 1737.4 km Moon. The 7.8 m figure is correctly called a worst-interval error, not an endpoint miss. Reference \(4.0152\) ND / \(17.798\) d is correctly rounded. |
| Lines 390–392 | **37% spread, two clusters, none within 1%, outcome A supported:** correct for the two accepted current rows. The iteration-limit row is correctly excluded. |
| Lines 442–449 | **4.22 → 0.996 km/s and 22.3% → 5.8%:** correct default-run rounding. These are physical km/s after applying \(l^\star/t^\star\), whereas \(c_{\rm nd}\ln(1/m_f)\) alone is ND. Label that distinction explicitly. |
| Lines 491–499 | Historical coarse rung set, \(3.9\times10^{-14}\) defect, 54.7 d, 46.8 km/s, approximately 16×, and 94% agree with FINDINGS 75. Default eleven accepted rungs, worst slack, and 0.432 → 3.358 d agree with the log. “Safe periselene” should be weakened to the reported clearance/node diagnostic unless a continuous check is cited. |
| Lines 603–605 | **100% vote, \(\lambda_t=1.000000\), mapping flag true, \(4.59\times10^{-12}\), one iteration:** all match the default log. Only the inference about the *initial* seed is unsupported. |
| Lines 710–714 | **Seven thrust rungs, three Isp rungs, 11.119 → 10.606 d, five unconverged exponents on each failed rung, no fence event or policy skip, 0.81 → 1.31 turns, \(1.8\times10^{-4}\) km:** correct. “3 s” is the sum of the current log’s rounded one-second rows; FINDINGS 75’s 4 s belongs to an earlier run, not a contradiction. Prefer “about 3 s.” |
| Lines 829, 865–866 | **\(2.7\times10^{-8}\) km** correctly rounds the current \(2.68\times10^{-8}\). The earlier 61 s convergence and later 300 s CPU-cap event for \(s_A=0.7837\) are supported by FINDINGS 77 and correctly separated. |
| Lines 881–885 | **10.606 d source; 11.240 d / +6.0%; 9.081 d / −14.4%; third probe CPU-capped; two accepted, one faster, two banked:** all correct for the supplied full round-2 run. |
| Lines 804–805 | “Three of five families” agrees with the three named discoveries—fast2, direct18, direct11—in FINDINGS 75 and the run output, but conflicts with the script header’s “four of five.” Do **not** change it to four merely to match that header. Safest replacement: name the three documented discoveries and omit the disputed total. |

## D. Pedagogical findings

- **[PEDAGOGY]** `root_origins_study_guide.tex:302–327,399–404,956` (Cold lottery and direct ladder — ELI5/Intuition) — **Several simplifications become false guarantees.** “Every run that finishes says converged” contradicts the finished iteration-limit row. “The trouble is never whether it converges” contradicts both the lottery and the stall. “The optimiser never has to search” is not what warm starting guarantees. Mesh sensitivity motivates continuation; it does not prove continuation is uniquely necessary.

  Suggested replacements:
  > Two of the three cold solves passed the acceptance checks; the third stopped at its iteration limit. Even the accepted solves returned different flight times.

  > The direct solver needs no previous converged transfer, but it still constructs an initial numerical guess. Both convergence and which discrete solution it finds matter.

  > Warm starts often make the next solve easier, but they neither eliminate the nonlinear search nor guarantee that the same branch is followed.

  Replace “continuation is required equipment” with “this motivates informed seeding or continuation rather than trusting one cold solve.” Also remove “barely bends”: a short flight or sub-turn azimuth does not establish small trajectory curvature.

- **[PEDAGOGY]** `root_origins_study_guide.tex:725–730,777–783,925–946` (Certification and what is not established) — **The inspector ELI5 overstates the certificate, and the limitations list omits its most important scope restrictions.** “Checks every condition the theory demands” and “a true minimum” sound like unconditional mathematical certification, potentially global. The supplied material shows a numerical gate stack with hypothesis and validity checks.

  Replace with:
  > The inspector polishes the candidate, checks a separately flown arrival and PMP conditions, obtains a foreign-solver numerical witness, and applies the configured second-order and hypothesis tests. A pass supports local minimality under those tests’ assumptions; it is not a global-optimality proof.

  Add to the final warning list:
  - No globally fastest transfer is established, nor are all families enumerated.
  - Results concern the specified CR3BP, ideal engine/mass model, and fixed endpoint phases—not mission robustness or a higher-fidelity ephemeris model.
  - The \(t_f\)/radius agreement screens do not establish root or branch identity.
  - Accepted direct solutions and harvested hunt seeds are not thereby certified continuous extremals.
  - Numerical integration and sampled clearance diagnostics are not automatically rigorous error-enclosed continuous proofs.
  - The displayed Hamiltonian describes the normal formulation; its applicability and the certifier’s hypotheses must be maintained.

  The existing distinctions between a policy stall, an uncertified candidate, and failure to reach 70 mN are **RIGHT** and worth retaining.

- **[PEDAGOGY]** `root_origins_study_guide.tex:530–552` (Handoff — “same solution”) — **The surrounding claim still sounds stronger than its two scalar screens.** Agreement within 1% in time and minimum node radius is useful for rejecting a gross move, but cannot establish identical trajectories or branches. The code comment correctly admits there is no aligned trajectory comparison.

  Change the snippet title from “must be the same solution” to:
  > “Sundman re-solve: time/radius consistency screens”

  Add:
  > Passing these screens permits the handoff; it does not prove that the two discretizations represent the same continuous root.

- **[PEDAGOGY]** `root_origins_study_guide.tex:275–289,734–748` (Three-level structure) — **The orbit section has no Intuition box, while certification’s Intuition is largely change history rather than a first explanation of the pipeline.** This weakens the promised ELI5 → reason → rigor progression.

  Add an orbit Intuition explaining why closure alone is insufficient for a phase interpolant: phase-wrap smoothness matters when nearby endpoint phases seed related solves. For certification, lead with the interface reason:
  > A certificate applies to a particular polished trajectory and to the options actually consumed by the certifier. Therefore both the options and the returned state/time representation must remain consistent.

  Move reviewer attribution and “both mine” into a short provenance note. History helps when it explains a failure mode; repeated review chronology distracts when it replaces the explanation.

## E. Notation and presentation

- **[STYLE]** `root_origins_study_guide.tex:149,191,298,503` (Section numbering) — **There are two incompatible section-number namespaces without a consistent distinction.** The article’s automatically numbered sections differ from the script’s §0–§8, while bare `\S4` and `\S\ref{...}` refer to different systems. A reader following “at §4” can reach the wrong place.

  Use “script §4” for script labels and “Section~\ref{sec:handoff}” for article references, or align the article numbering with the script. Apply the convention to the diagram and captions as well.

- **[STYLE]** `root_origins_study_guide.tex:337,443–447,769` (Rigor notation) — **Several notation choices make correct ideas harder to follow.** \(T\) denotes thrust in one paragraph and lifted flight time in another; \(\Delta V\) is ND in its formula but dimensional in the immediately following measurement; `C.z_8` is neither quite mathematical indexing nor MATLAB syntax.

  Prefer a distinct lifted-time symbol such as \(\mathcal T_k\), write \(\Delta V_{\rm nd}\) before the physical conversion, and define \(t_f^{\rm cert}:=\texttt{C.z(8)}\), followed by \(t_k^{\rm cert}=\sigma_k t_f^{\rm cert}\).

I found **no undefined custom mathematical macros** in the supplied source. The snippets match the corresponding script slices; several deliberately end inside an `if` block, which is acceptable for excerpts, not evidence of erroneous code. I have not rendered the LaTeX, so I cannot certify float placement or figure fit.

## Ranked top five changes

1. **Repair the Hermite–Simpson equations:** explicitly include normalized-time scaling and the lifted-time state.
2. **Make the PMP and certificate claims properly conditional:** preserve the correct signs, finish the zero-primer argument, show the Sundman Mayer convention, and distinguish local numerical evidence from global proof.
3. **Correct the handoff evidence:** remove smoke/default run mixing and stop treating node clearance, one iteration, or two scalar agreement tests as stronger proofs.
4. **Redraw or relabel the flow diagram:** show actual ordering, plain-time fallback, last-accepted-rung continuation, early fatal checks, and optional/skipped stages.
5. **Finish the diagnostic rigor:** add `unwrap` to the turns equation, describe the bordered shooting Jacobian accurately, and replace the universal 1% scaling claim with the measured approximately 2.1% worst step.