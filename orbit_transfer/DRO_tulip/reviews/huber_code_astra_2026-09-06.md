## A. Conjugate-point test

- **[CORRECTNESS]** `ms_conjugate_test.m:164-185` -- “Not yet full rank” is not equivalent to “only structural zeros.” The skip can hide a genuine early conjugate point, including a rank loss during a saturated first arc, and returns `pass=true` even when **no sample was testable**. Fix: distinguish `structural`, `regular`, `rank-loss candidate`, and `unresolved`; never convert an entirely untested interval into PASS. Determine structural degeneracies from the active dynamics, not from a numerical rank threshold alone.

  Testing the `6x7` r,v block is a useful complement, **but not by demanding rank six**. On an uninterrupted saturated initial arc, `lambda_m(0)` has no state effect, and scaling `(lambda_r(0),lambda_v(0))` leaves the thrust direction unchanged. Thus the r,v block has rank **at most five**. Remove the persistent control-invisible directions and monitor loss of the resulting expected rank, with a nontrivial Jacobi-field check. An arbitrary square minor, or the smallest singular value of the unreduced `6x7` block, is not a valid replacement.

- **[CORRECTNESS]** `ms_conjugate_test.m:149-175` -- There is no sign baseline before the first retained sample. A conjugate point before that sample is invisible even without structural degeneracy. For example, the supplied scalar LQ problem with `T=4, K=1` returns PASS: its only sample is already beyond `pi`. Fix: establish the initial orientation from a justified small-time expansion or additional sufficiently early samples, and adaptively sample the first active interval. Until that interval is checked, report incomplete coverage.

- **[CORRECTNESS]** `ms_conjugate_test.m:156-161` -- The equilibrated determinant does **not** actually protect the sign test from raw-determinant underflow. If `det(M)` underflows to zero, line 161 multiplies the correct equilibrated sign by zero; the result is subsequently counted as an exact root. Overflow similarly corrupts reporting. Fix: store determinant sign independently of magnitude, preferably using pivoted LU with permutation parity. For reporting, compute
  \[
  \log|\det M|
  =\log|\det M_e|+\sum_i\log r_i+\sum_j\log c_j,
  \]
  and never feed the reconstructed, possibly underflowed magnitude back into crossing detection.

- **[ROBUSTNESS]** `ms_conjugate_test.m:134,158-169,193-196` -- `rankTol=1e-13` is not an error-calibrated structural-rank test. The STMs are propagated at `RelTol=1e-10`, and chaining/equilibration can amplify their errors. Worse, a structurally zero row contaminated by tiny independent errors is divided by its own tiny maximum and becomes order one. Conversely, a genuinely regular but poorly resolved block can remain skipped. Fix: retain physical/reference-scale absolute row norms, propagate or estimate STM error, and classify singular values relative to that error after scaling. Repeat borderline cases with tighter propagation and refined sampling; return UNKNOWN when classification changes. Merely changing `1e-13` to another universal constant is insufficient.

- **[CORRECTNESS]** `ms_conjugate_test.m:173-185` -- Exact zeros are neither deduplicated nor distinguished from sign crossings. For `[+,0,-]`, the code reports two crossings: one sign flip across the omitted zero and one zero sample. A run of zero samples is counted once per sample. An exact terminal zero also leaves `atFinal=false`. Fix: group contiguous near-zero samples into candidate root intervals, merge them with adjoining sign-change brackets, and refine each candidate. Report root count separately from multiplicity and numerical uncertainty.

- **[CORRECTNESS]** `ms_conjugate_test.m:176-185` -- A root **exactly at** `tf` does not, by itself, refute local minimality. The supplied LQ example at `T=pi` has
  \[
  \delta^2J[\eta]=\frac12\int_0^\pi(\dot\eta^2-\eta^2)\,dt\ge0
  \]
  for zero-endpoint variations, with equality for `eta=a*sin(t)`. It is minimizing, though not strictly. Fix: distinguish verified interior roots from endpoint degeneracy. Keep counting roots strictly inside the last bracket, but refine that bracket before assigning a refutation; a root localized only to the boundary should produce an endpoint-degenerate/undetermined verdict, not automatic FAIL.

- **[CORRECTNESS]** `ms_conjugate_test.m:173-181` -- The detector misses even-multiplicity roots without an exact sampled zero. This is separate from the documented possibility of two distinct roots in one segment. For two independent copies of the LQ example, `det(Phi_xl)=sin(t)^2`: the conjugate point at `pi` produces no sign change. Fix: monitor reduced-block singular-value minima or use a suitable Jacobi-curve/Maslov-index method, with adaptive localization. Label a sign-only result accordingly rather than treating it as an exhaustive conjugate-point test.

- **[CORRECTNESS]** `ms_conjugate_test.m:142-146,180-186` -- Missing `Yend` silently drops the final segment, yet the output can still PASS and can set `atFinal=true` for a bracket ending at `t_K`, not `tf`. Fix: require `Yend` for a complete free-time test, or explicitly return `complete=false` and `sampledThrough=t_K`. Compute `atFinal` against the actual final time, not merely the last retained sample.

## B. Huber propagation

- **[CORRECTNESS]** `cr3bp_minfuel_prop.m:82-90,115-127` -- The claimed “integrator collapse throws” contract is not enforced. MATLAB ODE solvers can warn and return a partial trajectory. Both paths then accept its last row as the requested endpoint; the event path also interprets an early return with no event as successful completion. Fix: validate finite output and attainment of the requested endpoint, or a verified terminal event, after **every** integration call. Explicitly throw on unexplained early termination.

- **[CORRECTNESS]** `cr3bp_minfuel_prop.m:106-123,132-145` -- Initialization at exactly `Q=1` always chooses `lo`, independent of the outgoing direction. It then relies on initial-zero event behavior to correct the branch. MATLAB initial-step event semantics are not a reliable substitute for branch initialization; restarting exactly on the same event surface also lacks a progress/duplicate-event guard. Fix: at a wall, choose the outgoing branch from `sign(dt)*dQdt`, explicitly handle tangency, and implement a tested restart policy that suppresses only the already-processed initial root—not a subsequent genuine crossing. Verify that the selected event is the terminal event corresponding to the returned endpoint.

  For these PMP equations, away from primer singularity,
  \[
  \dot Q=-\frac{T_{\max}}{m\,\rho}\lambda_v^\top\lambda_r,
  \qquad \rho=\|\lambda_v\|.
  \]
  This supplies a direct initialization test without estimating a derivative from output samples.

- **[CORRECTNESS]** `cr3bp_minfuel_prop.m:123-127` -- Backward propagation is accepted by the interface but broken by `te(end) < dt - 1e-14`. For `dt<0`, an interior event lies between `dt` and zero and fails that condition; propagation returns at the first event rather than at `dt`. Fix: either reject `dt<0` explicitly or implement direction-aware event ordering, branch transitions, and endpoint comparisons, with backward round-trip tests.

- **[ROBUSTNESS]** `cr3bp_minfuel_prop.m:139-143` -- The grazing scale is effectively an absolute `1e-12` test, not a meaningful relative transversality test. In this model the throttle-dependent terms cancel from `dot Q`, so normally
  \[
  n^\top F_+=n^\top F_-.
  \]
  Comparing one against the other supplies almost no independent scale. The guard also runs only with `needSTM=true`, and unnecessarily rejects tangencies when `p=1`, where there is no jump. Fix: route `p=1` through continuous propagation; for `p<1`, use an error-aware denominator test based on scaled terms such as `sum(abs(n.*Fm))`, event-state uncertainty, and the size of
  \[
  \frac{(F_+-F_-)n^\top}{n^\top F_-}.
  \]
  Separate “numerically unsafe saltation” from “physical grazing,” and report the latter consistently in state-only propagation too.

- **[ROBUSTNESS]** `cr3bp_minfuel_prop.m:123-127` -- A switching event exactly at a shooting boundary is silently assigned the incoming STM. At a discontinuous switch with fixed endpoint time, perturbations can put the switch just before or just after that endpoint, so an ordinary two-sided endpoint Jacobian generally does not exist there. Blindly adding saltation at every terminal event would not resolve this either. Fix: flag switch-on-boundary cases and move the shooting mesh away from them, or explicitly support one-sided/generalized derivatives with a consistent junction convention. Replace the fixed `1e-14` exclusion with a documented, time-scaled boundary-event policy.

## C. Continuous Huber law

- **[CORRECTNESS]** `cr3bp_minfuel_pmp.m:110-120,167-173` -- `delta` is not validated. Negative widths bypass the intended ramp and reintroduce a jump; zero widths introduce divisions by zero into the generated expressions. Nonfinite parameters are also admitted, and `eps` accepts `p>1` despite its advertised domain. Fix: validate finite scalar `p`, enforce the declared family-specific range, and require finite scalar `delta>0` for `huberc`. Handle the limiting `delta=0` as an explicit dispatch to `huber`, not through the ramp formula.

- **[CORRECTNESS]** `cr3bp_minfuel_pmp.m:166-169,201-203` -- `huberc.aux.sRaw` is not the claimed unclipped stationary throttle on the ramp. For example, at `p=0.3, delta=0.3, Q=1.15`, the stationary throttle is `0.65`, while `sRaw` reports `0.345`. Fix: report the piecewise inverse-law continuation
  \[
  s_{\rm raw}(Q)=
  \begin{cases}
  pQ,&Q<1,\\
  p+(1-p)(Q-1)/\delta,&Q\ge1,
  \end{cases}
  \]
  for `p<1`, followed by clipping for the actual throttle. Alternatively rename the existing field `sCoreRaw` and stop advertising it as the family's stationary throttle. Handle `p=1` separately.

## D. Continuation race

- **[CORRECTNESS]** `run_minfuel_race.m:133-139,183-188` -- Stage 2 forgets its already-accepted starting point: `xGood=[]` disables bisection on the first failed delta step. The inserted stage-2 pair itself is correct; the missing initialization prevents its use on the first gap. Fix: initialize `xGood=dLast` on entry to stage 2 and retain the corresponding seed as the accepted endpoint.

- **[CORRECTNESS]** `run_minfuel_race.m:201-203` -- Retirement is cleared merely because `deltaSched` is nonempty. An `eps`/`huber` arm, or a `huberc` arm with no eligible delta steps, can therefore finish reported as not retired despite having exhausted stage 1. Fix: preserve per-stage retirement state, and reset an active-stage flag only after confirming that stage 2 will actually run.

- **[CORRECTNESS]** `run_minfuel_race.m:183-196` -- `maxBisect` is documented as a per-gap insertion limit but implemented as one global budget of `maxBisect*numel(P.sched)`. One gap can consume many times its allowance, and stage 2 competes for a budget sized only by stage 1. Fix: attach a target-gap identifier and insertion counter to queued rungs, enforce `maxBisect` per gap, and maintain `A.nBisect` separately as an aggregate statistic. Count abandonment by target gap, not accidentally by its intermediate retries.

- **[CORRECTNESS]** `run_minfuel_race.m:225-237` -- Floor fallback does not clear or trim `A.delta`. After falling back to rung `kT`, line 237 combines the retained `p` with the **discarded deepest rung's delta**, potentially testing a different problem. `A.z` also remains the discarded solution when a tight fallback exists. Fix: trim every per-rung field through one shared routine, including `delta`, and restore `z` from the retained rung. Store `z` per rung rather than only as an overwrite-only aggregate.

- **[CORRECTNESS]** `run_minfuel_race.m:216-242` -- The acceptance verdict is not bound to the reported final solution. Successful refinement updates `A.Y` but not `seed`; fallback also leaves `seed` at the discarded loose rung. The acceptance call then re-solves from that stale seed, discards its solution, and attaches its acceptance result to the old `A.mf/A.Y/A.z`. It can converge to another branch. Fix: rebuild the acceptance seed from the retained/refined rung and either validate that exact solution without re-solving or replace the reported solution and all metrics with the accepted solve's outputs. Assert binding between the accepted and packaged solution.

- **[CORRECTNESS]** `run_minfuel_race.m:216-220` -- Successful floor refinement leaves `coastFrac`, `iters`, and `wall` describing the loose solve while `mf`, `Y`, `Hdrift`, and `condJ` describe the refined solve. This corrupts both the final report and per-rung comparisons. Fix: refresh solution-dependent metrics and record refinement effort explicitly, either as a separate attempt record or accumulated rung iterations/wall time.

- **[ROBUSTNESS]** `run_minfuel_race.m:157-166,211-213,217` -- Several potentially expensive/failing propagations run synchronously outside `run_capped`, without exception handling. One accepted loose predictor that collapses during its single-shooting flight can abort the entire race before the step is saved. These calls also defeat an arm-level interpretation of the hard timeout. Fix: run diagnostic/reconstruction propagation through bounded, failure-aware helpers; save the solve attempt before diagnostics and reject or mark failed diagnostics without losing the other arms.

## E. Switch diagnostics

- **[CORRECTNESS]** `huber_switch_diag.m:32-39` -- `unique(T)` removes exact duplicate times, but does not repair crossing classification. For sampled signs `[-,0,-]`, replacing zero by positive creates two fictitious crossings at a tangential touch. Conversely, a short excursion between output samples is missed. `tCross=T(kx)` is the left sample, not a located root. Fix: consume the propagator's event records, exposing them if necessary; distinguish crossing, touch, and boundary events using neighboring signs and analytic derivatives. For continuous families, locate roots through event detection/dense-output refinement. Near-duplicate times require tolerance-aware handling, not just exact `unique`.

- **[CORRECTNESS]** `huber_switch_diag.m:37-45` -- `gradient(Q,T)` at the **left bracket sample** is not a reliable saltation denominator at the crossing. A highly nonuniform grid can bias that secant estimate; nearly coincident times magnify cancellation. It can make a transversal crossing appear nearly tangent or miss a true small denominator. Fix: evaluate
  \[
  \dot Q(t_e)=\nabla Q(y_e)^\top F_-(y_e)
  =-\frac{T_{\max}\lambda_v^\top\lambda_r}{m\|\lambda_v\|}
  \]
  at the located event state, with the implementation's primer regularization treated consistently. Record event-location and derivative uncertainty.

- **[CORRECTNESS]** `huber_switch_diag.m:41-49` -- An extremum within `[0.9,1.1]` is not evidence that a grazing bifurcation is imminent. The sample extrema are unrefined, flat/noisy samples produce extra extrema, and the `1e-3*tf` exclusion can remove a genuine nearby extremum belonging to a short burn arc. Fix: locate extrema as roots of analytic `dot Q`, record `Q_ext-1`, curvature, and numerical uncertainty, and track them across the continuation parameter. A grazing claim needs evidence approaching
  \[
  Q(t_*,a_*)=1,\qquad Q_t(t_*,a_*)=0,
  \]
  with appropriate nondegeneracy checks such as `Q_tt != 0` and parameter unfolding. Keep “within 10%” only as a configurable screening window, not a bifurcation verdict.

## F. Verdict production and packaging

- **[CORRECTNESS]** `run_conj_fixedtf_sweep.m:109-115`; `build_minfuel_catalog.m:48-49,110-117` -- An unconverged re-solve can publish `conjPass=true`, and packaging/swap selection never checks `v.converged` or residual quality. A Jacobi calculation along inconsistent shooting segments is not a verdict on a feasible extremal. Fix: use a tri-state verdict; unconverged, incompletely tested, or numerically unresolved cases must be UNKNOWN. Require convergence, feasibility/continuity gates, and valid test coverage before either swapping or setting `conj_pass=1`.

- **[CORRECTNESS]** `build_minfuel_catalog.m:110-117` -- The binding assertion verifies only approximate initial costates and state rows. It does not verify `v.p == g.pDeepest`, the costate columns, free-time/quotient settings, final time, endpoints, constants, grid, or consistency of `g.Y(:,1)` with `g.z`. Close costates are not a conditioning-independent identity check: `delta y(tf) ≈ Phi(:,8:14)*delta lambda0` can be large, and near a conjugate point even a small perturbation can change the verdict. Fix: package the exact solution assessed by the sweep, or recompute the verdict on the packaged solution. Bind a full problem-and-solution fingerprint and retain the full test specification, including quotient, rank policy, coverage, and numerical settings. A bare `1e-8` relative costate tolerance cannot establish verdict equivalence.

- **[CORRECTNESS]** `run_conj_fixedtf_sweep.m:115` -- The recorded `.mf` is `it.Y(7,end)`, which under the supplied `info.Y` convention is the start of the last segment, not final mass. Fix: use `it.Yend(7)` or propagate the last segment to `tf`; name junction mass differently if that is what is intended.

## Regression coverage

- **[ROBUSTNESS]** `test_conj_fixedtf.m:75-77`; `test_huber_saltation.m:39-58` -- The coast fixture uses a singular matrix `Z=[0 0;0 1]`, which cannot be the fundamental matrix of a smooth Hamiltonian coast flow. The switch test exercises only one interior crossing and measures one globally normalized matrix error, which can hide bad smaller columns. Fix: use an invertible symplectic coast fixture—for example identity propagation for a deliberately frozen coast subsystem—and add early reduced-rank conjugacy, even-multiplicity roots, exact terminal roots, determinant underflow, and no-testable-sample cases. Add ascending/descending switches, initial/boundary switches, multiple switches, `p=1`, and near-grazing cases; compare scaled columns/directional derivatives over a finite-difference step-size sweep.

## QUESTIONS / COULD-NOT-VERIFY

- **[ROBUSTNESS]** `ms_conjugate_test.m:123,150-154` -- The free-time construction is algebraically consistent **if** the dynamics are autonomous, the supplied costate direction is an actual scaling invariance, and all STMs/flows belong to the same continuous extremal. In that case the column is exactly
  \[
  [\,\Phi_{x\lambda}(t,0)P,\ f_x(y(t))\,],
  \]
  and `P` does not act on the flow column. The min-time caller and `ms_bvp` are absent, so those premises cannot be verified. Add checks for the claimed kernel `Phi_xl*qDir≈0`, continuity defects, flow dimensions, and final-state provenance. Positive max-abs equilibration itself is sign-safe for finite matrices, including exact zero rows/columns; the defects are its use as a structural-rank oracle and the raw determinant retained afterward.

- **[ROBUSTNESS]** `cr3bp_minfuel_prop.m:114,135-143,177-179` -- No algebraic defect is evident in the incoming-branch selection or STM carry: `Fm` uses the branch actually integrated, the event STM is extracted from `ze`, saltation left-multiplies it, and the restart propagates `A*PHI`. Runtime correctness still depends on the initial-event behavior discussed above. Add explicit assertions that the event state corresponds to the stopped integration and that successive processed events make progress.

- **[ROBUSTNESS]** `cr3bp_minfuel_pmp.m:170-173` -- For `0<p<1`, `delta>0`, the supplied matching constant and ramp expressions check symbolically:
  \[
  L(s)=
  \begin{cases}
  s^2/(2p),&0\le s\le p,\\
  p/2+(s-p)+\dfrac{\delta(s-p)^2}{2(1-p)},&p\le s\le1,
  \end{cases}
  \]
  \[
  L'(s)=1+\frac{\delta(s-p)}{1-p},\qquad
  L(1)=1-\frac p2+\frac{\delta(1-p)}2.
  \]
  Thus `delta→0` recovers Huber, and `p→1` gives `L=s²/2` with throttle `clip(Q,0,1)`. No matching-constant error is established. Replace the `1e-300` denominator workaround with an explicit `p=1` branch and add symbolic/AD limit tests, particularly at the knees.

- **[ROBUSTNESS]** `run_minfuel_race.m:114`; `cr3bp_minfuel_pmp.m:144-158` -- Scaling costate rows in **all** junction columns is appropriate for a multiple-shooting predictor; scaling only the initial column would introduce unnecessary costate defects. For the exact energy-to-`p=1` Huber conversion, scaling all costates by `1/2` while preserving state rows is the correct transformation. What remains questionable is applying one scalar to every family: `eps,p=1` needs scale one. Make `seedScale` family-specific or derive it from the first rung.

- **[ROBUSTNESS]** `run_minfuel_race.m:223-224` -- `tern(okR,itR.normR,NaN)` is eagerly evaluated in MATLAB. If `run_capped` returns an empty/nonstruct `itR` on failure, this fallback logger itself crashes. Its implementation is absent, so the triggering contract cannot be verified. Replace the expression with an ordinary guarded assignment regardless.

- **[CORRECTNESS]** `build_minfuel_catalog.m:36-49` -- Higher final mass is a legitimate **selection preference** between feasible solutions of the same exact fuel problem; it is not evidence of local minimality, and this sampled necessary-condition PASS is not a local-minimum certificate. For finite `eps`, the solved objective is instead
  \[
  J_\epsilon=(1-p)\int s\,dt+p\int s^2\,dt,
  \]
  whereas final mass ranks only `int s dt`. Establish whether the swap ranks exact-fuel candidates or finite-smoothing objective values, require comparable problems and numerical uncertainty margins, and label the result “selected candidate; no detected conjugate point,” not “certified local minimum.” No defect is established in using higher mass merely as that explicitly limited preference.