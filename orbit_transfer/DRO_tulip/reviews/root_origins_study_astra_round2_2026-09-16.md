## 1. Verification of the round-1 fixes

The most consequential new defect is **adopting polished junction states without updating their time grid**. The disputed starts-only seed shape is **not** a numerical defect in the supplied consumer implementations.

- **[FIX OK]** `root_origins_study.m:440–445` -- The delta-V explanation now agrees with the measured default ladder. At fixed exhaust speed, the all-burn law gives \(1-m_f=Tt_f/c\), so the stated predictor \(t_f\propto T^{-0.6}\) gives decreasing propellant consumption as thrust decreases. The computed `dv` at line 394 also correctly uses the mass actually consumed. For editable ladders, make the explanatory exponent follow `dladder.tfExp` and print this fixed-Isp explanation only when Isp actually stayed fixed.

- **[FIX PARTIAL]** `root_origins_study.m:287,377,482,787,909–919` -- Every direct solve now uses `directOK`, including throttle, interpolation residual and lifted-time spread. However, missing `maxInterp` and `tfSpread` pass through defaults of zero, and the throttle test checks only undersaturation, not finite values or upper-bound violations. Complete the predicate as described in section 3.

- **[FIX OK]** `root_origins_study.m:255–265,279–282` -- The historical altitude/radius conversion is correct, both conventions are labelled, and the 7.8 m datum is correctly qualified as a worst-interval error. No further numerical correction is needed.

- **[FIX OK]** `root_origins_study.m:699–704` -- `m0kg`, `gateKm`, `gateVms` and `moonKmMin` are forwarded through the options that `certify_root` actually reads, and the printed units and meanings match. This fixes the configuration-forwarding defect.

- **[FIX PARTIAL]** `root_origins_study.m:746–755,789–815` -- The cross-phase “same basin / another basin” classification has correctly become a candidate report with an explicit same-phase-baseline requirement. But a candidate can still be counted without a banked seed, and infrastructure or feasibility failures can still be summarized as nonconvergence. Fix those separate paths below.

- **[FIX OK]** `root_origins_study.m:308–319` -- `uniquetol(..., tol.tfCluster*min(tf), 'DataScale', 1)` implements the intended absolute tolerance tied to the shortest accepted flight time. A long outlier no longer changes that tolerance. “Flight-time clusters,” rather than basins, is the right name. Minor wording correction: tolerance clustering is not an equivalence relation under arbitrary pairwise chaining, so describe this as MATLAB’s representative-based clustering with the stated absolute tolerance.

- **[FIX OK]** `root_origins_study.m:596–601` -- The depletion guard now correctly distinguishes the seed policy from exhaustion: `maxPropFrac = 0.6` retains 40% mass; exhaustion is at `cnd/Tnd`. Skips are counted.

- **[FIX PARTIAL]** `root_origins_study.m:609–622` -- NaN position misses no longer pass, velocity is gated, and the witness must finish at the requested time with positive final mass. Those specific defects are fixed. The routine still validates only the final flight row, not the complete returned arrays or their dimensions. This is a continuation screen—not full physical validation—and lines 659–661 now correctly say so. Add whole-array sanity checks without duplicating the entire certifier.

- **[FIX OK]** `root_origins_study.m:395–400` -- Direct refusal precedence is corrected: convergence/feasibility failure precedes clearance and the plausibility band. The implementation no longer labels an unconverged iterate a branch change. The contradictory comment at lines 378–385 still needs rewriting.

- **[FIX OK]** `root_origins_study.m:589–590` -- An Isp-only rung now gets one attempt. Since `Tcur/TN == 1`, changing the exponent would otherwise reproduce the same seed.

- **[FIX OK]** `root_origins_study.m:572–573,651–652` -- `turns` is correctly defined and computed as total variation of rotating-frame projected lunar azimuth divided by \(2\pi\). The printed qualification about reversals is correct. Older comments still call this revolutions or a winding number; those are not fixed.

- **[FIX PARTIAL]** `root_origins_study.m:584,768,940–949` -- Both continuation reconstruction and hunt preparation now use banked junction starts instead of a sensitive full-arc reconstruction. The basic segment propagation, forced departure state, output orientations and time offsets are correct. The helper needs input/output validation and an explicit seam-sample policy; details below. This is **not** an output-arity bug.

- **[FIX PARTIAL]** `root_origins_study.m:468,501,584,604–615,768` -- Pool creation now precedes the handoff shoot, and the expensive indirect propagations are fenced when a pool exists. However, the advertised per-rung wall budget is still exceeded by construction, and the poolless policy later collides with `certify_root`’s default refusal to run unfenced.

- **[FIX WRONG]** `root_origins_study.m:709–714` -- **The polished root is adopted incompletely: `C.Y` is paired with the pre-polish `itCur.tGrid`.** In `ms_bvp`, the normalized breakpoints stay fixed but their physical times scale with the newly solved `tf`. Therefore the certified junction states belong at
  \[
  t_k^{\rm cert}=\frac{t_k^{\rm input}}{t_f^{\rm input}}\,C.z(8).
  \]
  If polishing changes `tf`, the hunt propagates certified starts for the wrong segment durations, and `S.reached` stores an internally inconsistent root. The certificate itself remains valid; the downstream representation does not.

  **Concrete fix:** preferably return `C.tGrid = it.tGrid` from the certifier. With the current contract, reconstruct it using the preserved normalized grid:
  ```matlab
  sig = itCur.tGrid(:).' / itCur.tGrid(end);
  zCur = C.z(:);
  assert(numel(zCur) == 8 && size(C.Y,1) == 14 && ...
         size(C.Y,2) == numel(sig)-1);
  itCur.Y = C.Y;
  itCur.tGrid = sig*zCur(8);
  itCur.normR = C.normR;

  reached.z = zCur;
  reached.Y = itCur.Y;
  reached.tGrid = itCur.tGrid;
  reached.tfDays = C.tfDays;
  reached.polished = true;
  ```
  Do not silently retain old `Y` if a supposedly successful certificate lacks its promised junction states.

- **[FIX OK]** `root_origins_study.m:705–738` -- Certification prose now branches on `C.ok` before making a certified-root claim or reference comparison. The short-ladder branch no longer prints “certified” on certificate failure.

- **[FIX OK]** `root_origins_study.m:715–721` -- Reference agreement is now explicitly flight-time agreement, not root identity. The stated need for trajectory/costate comparison or a continuation connection is correct.

- **[FIX OK]** `root_origins_study.m:182,190,230,501–502,605–607` -- The two specifically identified hardcoded tolerances—seam derivative and shooting solve target—were moved into `tol` and are used. The distinction between solve targets and acceptance thresholds is sound. This does not establish the broader claim that every effective tolerance is exposed: H2/H3 and many certification settings still use independent constants/defaults.

- **[FIX OK]** `root_origins_study.m:421,435–438` -- Direct completion now tests acceptance of the final requested row, not just matching its thrust. The `B1b` label is consistent.

- **[FIX PARTIAL]** `root_origins_study.m:289–301,408–427,581–635,798–815` -- Returned/refused direct rungs and thrown hunt probes are now recorded. But indirect reconstruction failure exits before appending a rung, exponent attempts are retained only as aggregate counts plus the winner, and worker errors are counted as timeouts. “Every attempt is recorded with its reason” remains too strong.

- **[FIX PARTIAL]** `root_origins_study.m:702,842–848,952–970` -- The certificate, supplied certifier options, gates and outcomes are now saved; `saveq` no longer silently swallows save errors, and the final record line honors its return value. The archive still lacks the full run configuration/effective tolerances, and `certOpts` contains the live pool object. Save a value-only options snapshot and the full configuration; use a temporary file followed by replacement if checkpoint integrity matters.

- **[FIX PARTIAL]** `root_origins_study.m:20–24,57–62,127–141,155–159` -- Runtime prose was softened, but major introductory/settings comments were not. They still assert “Every mesh converges,” a “winding wall, not a sensitivity one,” and that a timed-out probe indicates “no nearby basin.” Replace these with the measured discrete-solve sensitivity and configured-search-stall descriptions already used later.

- **[FIX OK]** `root_origins_study.m:893–894` -- The direct wrapper now canonicalizes endpoints to rows, satisfying its advertised row-or-column wrapper contract and the direct solver’s row-input contract.

Two other applied improvements are still incomplete:

- **[FIX PARTIAL]** `root_origins_study.m:482–498,512–535` -- The Sundman re-solve now receives the shared feasibility screen, a flight-time comparison, explicit fallback reasons and a basic multiplier check. It still does not establish “the SAME solution”: no aligned trajectory comparison is made, the commented node-radius comparison is not enforced, and multiplier/control dimensions and finite values are not fully checked. Report time agreement as such; validate harvestability separately.

- **[FIX PARTIAL]** `root_origins_study.m:506–526` -- Sign-vote and time-costate diagnostics now have explicit displayed checks, and the half-step-association overclaim is corrected in the runtime output. But H3 passes missing/nonfinite values and ignores `lamTOK` in its verdict. That remains a fail-open evidence check.

## 2. Adjudication of the disputes

### (a) Junction-array shape

- **[PEDAGOGY]** `root_origins_study.m:694,931,940–945` -- **The triage is right about the supplied numerical path; my round-1 numerical-defect claim should be withdrawn.** The relevant operations are:

  - `ms_tfmin` overwrites only `seed.Y(1:7,1)`.
  - `ms_bvp` obtains `K` from `numel(seed.tGrid)-1`.
  - Its unknown packing reads column 1 and columns `2:K`.
  - `junctions(p)` reconstructs exactly K starts.
  - The last segment’s endpoint is propagated and passed to the terminal conditions; it is not read from `seed.Y`.
  - `flyFromJunctions` likewise obtains K from the grid and reads columns `1:K`.

  No shown consumer reads seed column K+1 or infers K from the number of junction-state columns. The `size(seedB.Y,2)` use at line 509 is merely display. Thus a **14×K starts-only array plus the K+1 time boundaries is lossless for these consumers**.

  **Concrete fix:** correct the library headers to accept/document K starts, optionally allowing a redundant endpoint column. Add shape checks against `numel(tGrid)-1`. Do not append an endpoint merely to satisfy the obsolete header.

  The actual representation defect in this revision is the **stale time grid after polishing**, not the missing endpoint column.

### (b) Free throttle and the \(10^{-3}\) gate

- **[PEDAGOGY]** `root_origins_study.m:176–188,352–357,917` -- **Leaving throttle free and checking the result is a defensible experiment; pinning it is not mandatory.** My round-1 recommendation was too categorical. A relaxed direct solve that demonstrably saturates can provide an appropriate all-burn seed.

  For the supplied regular, normal, positive-mass unconstrained PMP, there is also more analytical structure than “purely empirical saturation.” Writing throttle as \(q\),
  \[
  H=1+\lambda_r^\mathsf Tv+\lambda_v^\mathsf Tg
       -Tq\left(\frac{\|\lambda_v\|}{m}+\frac{\lambda_m}{c}\right).
  \]
  With free terminal mass,
  \[
  \dot\lambda_m=-\frac{Tq\|\lambda_v\|}{m^2},\qquad
  \lambda_m(t_f)=0,
  \]
  so \(\lambda_m(t)\ge0\). Away from zero-primer degeneracies, the throttle coefficient is strictly negative and minimization selects \(q=1\). This argument depends on those assumptions; it does not certify the returned discrete NLP vector.

  **Is \(10^{-3}\) defensible?** Yes, as an explicitly approximate **seed-compatibility screen**, supported by a later all-burn shoot and certification. No, as proof of exact saturation or of negligible trajectory error. The observed \(1.3\times10^{-6}\) slack makes the previous \(10^{-6}\) cutoff brittle, but does not uniquely justify \(10^{-3}\), nor prove that every future deviation below it is barrier slack.

  It fails to catch:

  - sustained 0.1% under-throttling;
  - throttle overshoots above one;
  - deviations between the sampled stations;
  - malformed throttle values hidden by summary reductions;
  - potentially amplified trajectory differences despite a small control discrepancy.

  The statement “a real switch takes the throttle to 0” describes nonsingular bang-bang solutions, not every discretization artifact or imperfect solution that should be rejected.

  **Better diagnostics:** retain the maximum sampled deficit, but supplement it with a **physical-time-weighted integrated deficit**
  \[
  D_q=\frac1{t_f}\int_0^{t_f}(1-q(t))_+\,dt,
  \]
  upper-bound violation, and physical-time occupancy below selected thresholds. Use Hermite–Simpson quadrature; for Sundman, include \(dt/d\tau\), rather than averaging stations uniformly.

  The corresponding missing propellant consumption is
  \[
  \Delta m=\frac{Tt_f}{c}D_q.
  \]
  Under a genuinely pointwise \(q\ge0.999\) bound, this is at most \(10^{-3}Tt_f/c\), but that is not an endpoint-error bound. An integral alone is also insufficient: it can hide a narrow deep dip. Report both integral and maximum statistics. To substantiate the barrier explanation, compare against tighter solves and, where available, bound-multiplier/complementarity information.

Also, **retain `returnModel=true`** in the present direct calls. Against the supplied producer, my earlier suggestion that it was unnecessary was not justified: constraint registration used for multiplier extraction is conditional on that option. Stripping the model afterward does not undo the need to request the registered solve.

## 3. New and remaining findings in the rewritten sections

### Flight reconstruction and resampling

- **[ROBUSTNESS]** `root_origins_study.m:940–949` -- `flyFromJunctions` trusts its grid and segment outputs without checking them. Validate a real finite, strictly increasing grid starting at zero; 14 state rows and K start columns; finite segment arrays; positive mass; matching time/state sample counts; and completion at each requested `dt`. This is boundary hardening, not evidence that the measured reconstruction failed.

  Several suspected issues are **not defects**:

  - Forcing `[rv0;1]` in column 1 matches `ms_tfmin`.
  - Restarting local time at zero is valid for this autonomous rotating-frame CR3BP.
  - Adding `it.tGrid(k)` produces the correct absolute times.
  - The output is samples-by-14, as `flight_to_junctions` requires.
  - The event function does not simply truncate the returned flight: `tfMinProp` restarts inside its `while tau(end) < tf` loop. A nonprogressing restart can hang, which is why the fence matters.

- **[ROBUSTNESS]** `root_origins_study.m:946,602–603` -- The seam rule drops the next segment’s initial sample, retaining the previous segment’s propagated endpoint. Thus resampling at a junction does not necessarily recover the exact banked start. For these accepted roots, the discrepancy should be at the shooting-residual level, so this is **not a demonstrated significant numerical error**.

  **Concrete fix:** measure the seam mismatch and document which side wins. To preserve banked starts exactly, omit each nonfinal segment’s endpoint and retain the following segment’s start. Keep the final endpoint. Do not silently conceal a large seam mismatch with interpolation.

- **[ROBUSTNESS]** `root_origins_study.m:770–777` -- `unique(...,'stable')` correctly preserves chronological order **when the input is already ascending**; it does not sort or validate a malformed flight. On a valid reconstruction, the deduplication, spline output transposes, primer columns 11:13 and control dimensions are correct. Complete the checks with finite times/states, strictly increasing deduplicated times, zero start, correct final time, and finite positive primer norms. `Inf > 1e-12` currently passes the primer test and can produce `Inf/Inf` NaNs.

- **[EFFICIENCY]** `root_origins_study.m:584,602–603` -- The ladder reconstructs and resamples a flight even though its segment count and normalized grid remain unchanged. It can preserve the banked starts directly, scale the grid to `tfGuess`, and replace the mass row using the new all-burn law. No endpoint column is needed by the next shoot. Keep segmentwise reconstruction for the hunt or an actual grid change; this removes unnecessary propagation and avoids perturbing the junctions through interpolation.

### Fences, witnesses and failure records

- **[ROBUSTNESS]** `root_origins_study.m:699–701` -- **The advertised poolless fallback cannot complete this script.** `fenced` runs directly with `pool=[]`, but `certify_root` asserts unless `allowUnfenced=true`; the script neither supplies that option nor catches the assertion. A long poolless ladder therefore aborts at certification instead of reporting an outcome.

  **Concrete fix:** expose an explicit poolless policy. Either refuse early, or record an unavailable certificate and continue reporting, or require deliberate user opt-in to unfenced certification. Do not silently enable unbounded calls.

- **[ROBUSTNESS]** `root_origins_study.m:584,592,604–615` -- The 600 s “per-rung budget” is still not a hard per-rung cap. Reconstruction gets a fixed 120 s; `capLeft` has a 30 s minimum; shooting gets another 90 s beyond that; and a successful shoot can start a further fixed 120 s witness after the deadline. Cancellation adds overhead too.

  **Concrete fix:** use one rung deadline. Recompute remaining time before reconstruction, shooting and validation; never start a stage without sufficient remaining budget. Separate any deliberate cancellation grace from the advertised work budget.

- **[ROBUSTNESS]** `root_origins_study.m:584–588,608,628–635,856–876` -- The taxonomy still conflates execution failure with timeout. `run_capped` returns `ok=false` for either timeout or worker error, so line 608 cannot label that Boolean “timeout.” Reconstruction failure is also described as exceeding 120 s even if the worker errored, and that rung is not recorded. `poolUnusable` can throw outside a per-rung catch.

  **Concrete fix:** extend the capped-call result with a reason code and exception details; record reconstruction, shoot and witness stages separately. Append the rung before stopping on reconstruction/infrastructure failure, and preserve one record per exponent, not just aggregate counters.

- **[ROBUSTNESS]** `root_origins_study.m:614–617,768,870–876` -- The `nout=2` wiring is correct: each call yields `{ok,time,state}`, including the anonymous forwarding handle around `tfMinProp`. No Hamiltonian output is required for the stated arrival witness, and `tfMinProp` computes its final Hamiltonian regardless of requested output count.

  The remaining issue is defensive validation: require a positive finite scalar `zt(8)`, an N×14 finite flight, N finite times, and the required endpoint before indexing it. The current `reached` test is a correct completion/final-state screen for normal outputs, not a complete flight validator.

### Acceptance and claims

- **[ROBUSTNESS]** `root_origins_study.m:909–919` -- `directOK` is fail-open for absent or empty `maxInterp` and `tfSpread`, because both default to zero. The supplied direct producer normally emits both, so this does **not** demonstrate a false pass in the measured run. It is nevertheless the wrong contract for a required diagnostic.

  `thrMin`’s default of zero is different: it correctly fails closed at the configured \(10^{-3}\) gate. But a present `thrMin=Inf`, or an all-above-one throttle history, passes its one-sided test.

  **Concrete fix:** require real finite scalar diagnostics in their mathematical domains; use NaN defaults or explicit required-field checks. Validate finite `X`, `U`, `Xm`, `Um`, and compute both throttle bounds. Save the individual diagnostic values and verdicts, not only the first refusal reason.

- **[PEDAGOGY]** `root_origins_study.m:507,519–526,833–845` -- H3 counts an unavailable or nonfinite time-costate check as a passed gate. On a Sundman mapping failure, it can even print “NOT EXPOSED without Sundman” despite having selected Sundman. A missing `lamTOK` also defaults to “the mapping’s own check agrees.”

  **Concrete fix:** use `PASS / FAIL / NOT APPLICABLE` explicitly. Plain time is not applicable and should be excluded from the checked-gate denominator; Sundman requires a finite `lamT` and the mapping’s promised validity flag. Never translate missing evidence into agreement.

- **[CORRECTNESS]** `root_origins_study.m:790–791,798–810` -- **An accepted hunt candidate need not actually be banked.** If multiplier extraction returns `lamDef=[]`—an explicitly supported producer failure mode—`okH` remains true, `seedH` remains empty, and the summary still says every candidate has a harvested seed.

  **Concrete fix:** track `directAccepted`, `harvestOK` and `seed` separately. Preserve the accepted numeric direct solution when harvesting fails, and print the actual number of banked seeds. Validate the complete harvest inputs before calling the helper; do not let a harvest exception erase an otherwise useful direct result.

- **[PEDAGOGY]** `root_origins_study.m:792,807–815` -- `nConv` actually counts accepted candidates, not solver convergence. It can be zero because all successful NLP solves failed feasibility, because every probe threw, or because no probes were configured. “No probe converged within its budget” does not distinguish these cases. A clearance-only refusal also produces the uninformative `refused: ` because `whyH` is empty.

  **Concrete fix:** report attempted / solver-converged / accepted / harvested counts. Use `INCONCLUSIVE` for infrastructure-only or zero-attempt runs, distinguish feasibility refusal from nonconvergence, and give clearance its own reason.

- **[PEDAGOGY]** `root_origins_study.m:94–100,203–204,288,325–328,715–719` -- The reference still has no immutable problem metadata. Editing orbit, phase, mass or target engine leaves “the campaign certifies … here” and the A/reference comparisons active.

  **Concrete fix:** store the reference together with its orbit definitions, phases, engine, mass and unit constants; enable comparisons only on a matching problem. Describe A as mesh sensitivity of accepted **discrete solutions**, not proven selection of distinct continuous extrema.

- **[PEDAGOGY]** `root_origins_study.m:674–677,705–738` -- B’s outcome is determined before certification. Its literal “reached the target engine” statement can remain true when certification fails, but it must not be interpreted as a certified end-to-end success. The hardcoded “cold 15 N start” also becomes false if the exposed ladder is edited.

  **Concrete fix:** print the actual initial engine and distinguish “target reached numerically” from “target root certified.” The default run’s `PARTIAL` outcome is appropriate.

- **[PEDAGOGY]** `root_origins_study.m:266–268,323–324,378–385,475–476` -- Several surviving sentences still overstate or contradict their own qualifications:
  - the historical node-safe, accurate discrete result is called a “minimum-time extremal” without the corresponding validation here;
  - accepted solver statuses and defect sizes are asserted to be “the same”;
  - leaving the time band is said necessarily to leave the branch, immediately before admitting a same-branch sharp turn can trip it;
  - a constant time-costate row is again said to check station association.

  **Concrete fix:** use the later, accurate formulations consistently throughout the file—not just in selected runtime paragraphs.

- **[ROBUSTNESS]** `root_origins_study.m:238–239,702,842–843` -- The final archive is still not a self-contained study record. It lacks the complete orbit/configuration objects, ladders, hunt settings and acceptance tolerances. Conversely, it includes the live `pool` through `certOpts`, which is runtime infrastructure rather than reproducible configuration.

  **Concrete fix:** save a value-only configuration snapshot, effective certifier settings, provenance, and explicit pre-/post-polish representations. Exclude the pool handle. Keep the accepted direct numeric data or harvested seeds needed to reproduce discoveries.

## 4. Audit of the measured default run

| Printed result | What the supplied evidence supports |
|---|---|
| **A SUPPORTED** | Yes, as discrete mesh sensitivity. Two accepted rows have 4.6809 and 6.4126 ND, giving approximately 37.0% spread and two clusters; neither is within 1% of the supplied default reference. The N=800 iterate is correctly excluded by both status and defect. This does not establish two validated physical basins. |
| **B PARTIAL** | Yes. Eleven direct rows are accepted; ten of twelve attempted indirect rungs are accepted; the held engine is 0.12 N / 900 s rather than 0.07 N / 900 s. |
| **The two final refusals were not timeouts** | Supported by the recorded paths: each reports five returned candidates failing the residual/finite-root acceptance test, with no `okRun=false` outcomes or seed-policy skips. This identifies the configured search’s numerical failure, not its mathematical cause or a physical wall. |
| **Certification at 0.12 N / 900 s** | Consistent with `C.ok` and the supplied certifier contract. The displayed residual/witness/conjugate/lift evidence is favorable, but the complete certificate—not these few printed fields—is the audit record for all gates. |
| **C CANDIDATES** | Yes. Three accepted discrete candidates at other phases, with two lower returned flight times, is exactly what the table shows. No basin identity is earned or claimed. Whether all three seeds were actually banked requires inspecting `S.hunt`; the output alone cannot establish that. |
| **6 of 7 machinery gates passed** | Consistent with the default code path: B1, B1b, H1, H2, H3 and C1 pass; B2 fails. The H3 unavailable-data defect does not appear in this run because a finite Sundman value is printed. |

Additional qualifications:

- **[PEDAGOGY]** `root_origins_study.m:404–406` -- The default output does **not display the first 15 N rung’s throttle slack**, precisely the rung that previously exposed the gate problem. Its acceptance implies it passed `directOK`, but the full table is not an independent printed audit of every throttle value. Print the statistic on the first rung too, or print the worst accepted-rung statistic afterward.

- **[PEDAGOGY]** `root_origins_study.m:653–654,706–708` -- Printed `0.0 km` and `0.000 km / 0.000 m/s` are rounded measurements, not exact zero. In particular, `0.0 km` only displays resolution of 0.1 km. Retain scientific-notation witness errors in the record or a detailed summary.

- **[PEDAGOGY]** `root_origins_study.m:20,57–59,136–141,158` -- The measured output directly contradicts the surviving “Every mesh converges” introduction and does not support the winding-wall or no-nearby-basin explanations. The rewritten runtime stall explanation is substantially more accurate than these comments.

I see **no contradiction between the numerical tables and the executable default branches**. The stale polished-grid defect could be numerically dormant here because the pre-certification residual is already below the certifier’s polish target; the transcript does not reveal how much `tf` changed. It still needs fixing.

The final `ORIGINS DONE` and `ORIGINS RUN EXITED` lines are not emitted anywhere in this file; they must come from the run wrapper or transcript processing. That is provenance, not a numerical inconsistency.

## Ranked top-5 changes next

1. **Synchronize the polished root’s grid and states** at lines 709–714; add a regression test in which polishing changes `tf` appreciably.
2. **Repair execution policy and records:** explicit poolless certification behavior, genuine remaining-time budgets, and separate timeout/worker-error/reconstruction-failure records.
3. **Make evidence checks fail closed:** required direct diagnostics, finite throttle histories and both bounds, plus H3’s explicit unavailable state.
4. **Guarantee actionable hunt records:** separate accepted NLP solutions from successful harvests, retain failed-harvest numeric solutions, and save the complete value-only configuration.
5. **Remove the remaining unsupported claims:** immutable reference metadata, numerical versus certified B outcomes, and consistent mesh-sensitivity/search-stall prose throughout the script.