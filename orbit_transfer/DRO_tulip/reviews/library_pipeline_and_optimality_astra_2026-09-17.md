## Scope and overall verdict

**The supplied pipeline supports a substantial numerical assessment of fixed-phase PMP extremals. It does not presently establish a rigorous certificate of strict strong local optimality.** There are also concrete defects that can lose candidates, reuse the wrong continuation lineage, mislabel resumed results, or prevent advertised configurations from running.

Two review limitations matter:

- `oc.ms_bvp` and `oc.ms_conjugate_test` are **not supplied**; the supplied files are delegates. Neither are `conj_resolve`, `mintime_prop_seg`, `validate_flight`, `h6_margin`, `phase_state`, or the final catalog builder. Their internal correctness cannot be inferred from their callers.
- This is a static review. No MATLAB executions or catalog measurements were available. Document references below use the supplied section/equation labels rather than invented line numbers.

Paths below are relative to `orbit_transfer`; unqualified application filenames are in `DRO_tulip/indirect/`.

# A. Q1 — Pipeline findings

### Entry points and configuration propagation

- **[CORRECTNESS]** `run_costate_library.m:177–187,626–638` — **The default 24-column arrival grid is rejected by the default front door.** Its last phase is
  \[
  0.0754+23/24=1.033733\ldots,
  \]
  but `checkPhaseList` requires every value to be below 1. Simply applying `mod` also violates the current strictly-increasing requirement. **Fix:** canonicalize generated grids consistently—e.g. `sort(mod(...,1))`—while keeping the continuation anchor phase separate from the first grid element; alternatively explicitly support cyclically ordered lists. **Test:** the documented `run_costate_library(struct('nD',24,'nA',24))` must get through grid validation without any solves.

- **[CORRECTNESS]** `run_costate_library.m:189,239–243,270–272` — **The requested departure origin is not passed to sheet construction.** The call omits `'sD',sD0`, so setup defaults to zero and a nonzero departure origin fails the subsequent identity check. **Fix:** pass the departure phase explicitly and provide a suitable anchor independently of the grid origin. **Test:** a campaign with `sD=[0.2,0.4,0.8]` must construct and certify its spine at 0.2, not zero.

- **[CORRECTNESS]** `run_costate_library.m:275–277,304–308` — **The advertised single-departure-phase case always returns `blocked`.** It deliberately empties `cols`, then the unconditional empty-column check returns before packaging. **Fix:** distinguish “no certified spine” from “no ribs required”; bypass the queue when `nPts==0` and package the sheet. **Test:** one certified `(sD,sA)` cell must produce a one-entry catalog without rib artifacts.

- **[CORRECTNESS]** `build_arrival_sheet.m:53`; `build_ribs.m:52–61`; `fill_holes_direct.m:96–100`; `arclength_arrival.m:112–155` — **Physics reconstruction unnecessarily depends on the default 70 mN anchor.** Both rib construction and hole filling reload and re-polish that anchor, even though they already have certified roots or need only endpoint closures. Another engine, departure phase, or orbit pair can fail before its own seed is used. The filler’s assertion that “only its layout is used” is false: setup re-solves it and enforces the time-change guard. **Fix:** use `physicsOnly=true` for ribs/filling and construct the minimal `anc.sD` from the sheet; separate residual/chart construction from anchor polishing for sheet rescans. **Test:** remove the default anchor file and build ribs from a self-contained non-default sheet.

- **[CORRECTNESS]** `build_70mN_library.m:73–78,105–114,154–155,199–204,219–220,251–252` — **The hand-driven chain still drops configuration at several interfaces.**
  - Anchor phases are calculated before grid overrides.
  - `setupOpts` omits `pmTulip`.
  - Sheet construction does not receive `arcDirUsed`, although plotting and the family map use it.
  - Rib construction receives `nPts` but neither `nD` nor explicit targets. A 24-row request therefore walks with the default `1/12` step and can revisit the same 12 phases.
  
  **Fix:** derive anchor specifications after overrides, forward the complete problem descriptor, pass the exact arc list, and pass `rib_targets(grid.sD,spinePhase,direction)`. **Test:** a nonuniform departure list, northern tulip, and non-default arc directory must survive all stages unchanged.

- **[ROBUSTNESS]** `run_costate_library.m:368–374,655–671,720–735` — **The generated finalizer does not reproduce the invocation.** In particular, `.onlyA` is omitted: a queue initialized for a subset can be reopened by the finalizer for all certified columns. The finalizer also forces stage switches and drops options such as `foreignPattern` and `pictures`. Both generated scripts hardcode one developer’s pumpkyn installation. **Fix:** serialize a validated campaign specification into a MAT file and have both jobs read it; make dependency startup configurable. **Test:** finalize a two-column subset campaign on a different installation path.

- **[ROBUSTNESS]** `run_costate_library.m:147–150,240–241,371–373,492–495`; `build_70mN_library.m:53` — **Only `outDir` is made absolute.** Relative arc, seed, and extra-rib paths can change meaning after generated jobs or the chain change directory. Missing extra ribs are additionally filtered out silently. **Fix:** canonicalize every artifact path at the front door and fail on explicitly supplied missing artifacts. **Test:** invoke from a directory different from both the application and output directories, using relative input paths.

### Continuation, folds, and crossings

- **[CORRECTNESS]** `build_arrival_sheet.m:70–78`; `crossings_from_arc.m:108–130`; `costate_common/arclength_ms.m:261–307` — **Default rescanning loses the walk’s fold-aware crossings.** The walker splits a step at a localized turning point and can find two crossings even when its stored step endpoints are on the same side of a level. The rescanner tests only adjacent stored `A.q` endpoints and replaces `A.crossings` wholesale. The localized turning-point state is not retained in `A.folds`. Thus a refined/rebuilt sheet can lose both candidates beside a fold, including candidates the original walk already found. **Fix:** persist turning-point `(q,p)` records and share one segment-splitting crossing extractor between online and offline operation. Until then, preserve and merge existing crossings rather than replacing them. **Test:** the unit-circle fixture crossing level 0.999 twice inside one step must produce the same two roots online and after rescanning.

- **[CORRECTNESS]** `crossings_from_arc.m:105–115` — **Exact terminal hits are direction-dependent, and the wrap search is artificially bounded.** A final segment with `s(i)>0, s(i+1)==0` is not selected; approach from below is selected. Exact interior hits can be emitted twice because deduplication includes the segment index. A one-point arc emits no anchor hit. Also, `nWrap=3` contradicts the arbitrary-period coverage claim. **Fix:** implement explicit endpoint ownership, handle isolated exact roots, and derive integer wrap bounds from `min(A.q)` and `max(A.q)`. **Tests:** increasing/decreasing terminal hits, a one-root arc, and an arc around unwrapped phase 5.

- **[CORRECTNESS]** `costate_common/arclength_ms.m:211–213,237–253` — **Fold resolution has both a stale-value bug and an ineffective tolerance construction.** The successful `abs(tq)<1e-10` exit occurs before assigning `tqLast`; first-iteration success leaves it `NaN`. Separately, using `10*sminX/sminAug` as the event tolerance confuses event residual with numerical uncertainty. For a scalar residual Jacobian `[a b]`,
  \[
  |\tau_q|=\frac{|a|}{\sqrt{a^2+b^2}}
           =\frac{\sigma_{\min}(R_x)}{\sigma_{\min}([R_x,R_q])},
  \]
  so that part of the “resolved” comparison cannot reject an inaccurate event location. **Fix:** assign the final tangent component before every exit; terminate localization using an independent arclength bracket/location error and an SVD backward-error estimate. **Test:** force localization to stop far from a scalar fold and require `resolved=false`.

- **[CORRECTNESS]** `build_arrival_sheet.m:53,70–78`; `certify_crossing.m:57–62`; `arclength_arrival.m:165–179` — **Every saved arc is reinterpreted using one freshly reconstructed chart, without checking its identity or mesh.** There is no comparison against `A.B.problem`, `A.anc.K`, or its departure phase. A different mesh causes a layout error; a same-size arc from another problem can be corrected under the new problem and subsequently attributed to the old family. **Fix:** validate complete arc identity before rescanning, reconstruct its chart from its stored mesh, and explicitly classify foreign-problem inputs as seeds rather than continuation evidence. **Test:** offer two same-size arcs differing only in departure phase or tulip branch; neither may silently become the other’s lineage.

### Sheet and rib interfaces

- **[CORRECTNESS]** `build_arrival_sheet.m:103–116`; `sheet_from_arcs.m:148–153`; `sheet_to_catalog_file.m:124–131` — **A seed can be certified at a nearby phase and exported as an exact-grid root.** Seeds are accepted within `1e-6` of a grid phase, but certification targets `lib(k).sA`, not the grid phase. The sheet then places the certificate in the nearby column. This is particularly consequential for the sensitive arrival orbit. **Fix:** use the nearby entry only as a seed; re-polish and certify against the exact selected grid phase and stamp that phase in the result. **Test:** perturb a seed phase by `5e-7`; the exported root must solve the exact grid endpoint.

- **[CORRECTNESS]** `rib_from_crossing.m:98–104,130`; `rib_targets.m:47–55` — **The rib walker discards the unwrapped-path contract.** `rib_targets` deliberately chooses a directed offset, but the walker uses shortest wrapped differences and resets `sD` modulo one after each target. For `sD=[0,0.1]`, a negative-direction target is `-0.9`; with `nD=2`, the walker changes its step to `+0.1`. It reaches the correct phase by the wrong continuation path, potentially another branch. **Fix:** retain an unwrapped current coordinate and use `target-cur`, with direction-preserving step clipping. **Test:** record every trial for sparse, nonuniform grids spanning more than half a turn.

- **[CORRECTNESS]** `sheet_from_arcs.m:116–118`; `build_ribs.m:81`; `sheet_to_catalog_file.m:126`; `rib_validate.m:59,75,92`; `package_phase_catalog.m:106–109` — **Root and identity matching are inconsistent across interfaces.** Examples include phase matching at `1e-6`, `1e-8`, and `1e-9`; problem matching at `1e-9` versus `1e-12`; and choosing a winning certificate by time alone within `1e-9` days although distinct equal-time roots are explicitly supported elsewhere. **Fix:** centralize separate contracts for phase identity, problem identity, and numerical root equivalence. Store a winner index/root ID in the sheet rather than recovering it from `TF`. **Tests:** two distinct equal-time roots; a rib whose identity passes the current validator but fails the current packager; a nearby-but-not-exact phase seed.

- **[CORRECTNESS]** `run_costate_library.m:238`; `rib_from_crossing.m:82–89`; `run_costate_library.m:321–332` — **Resume identity omits the spine root and certification policy.** An existing sheet is reused even when `.sheet` is on and new arcs/seeds have been supplied. Rib checkpoints and completed rib validation identify the problem and grid, but not `C0.z`, the chosen spine candidate, or gate policy/version. Replacing a spine with a faster family can therefore resume or reuse the old family’s rib. **Fix:** bind each artifact to input-content hashes, spine root ID, endpoint-rule version, and certification-policy hash. Invalidate descendants when any binding changes. **Test:** replace a sheet winner with a different certified root at the same phase and verify that old checkpoints cannot resume as the new rib.

- **[CORRECTNESS]** `build_arrival_sheet.m:90–92,115–126`; `build_ribs.m:88–91`; `direct_cell_solve.m:79–81` — **“One gate policy” is not propagated across the pipeline.** The sheet stores `S.policy`, but ribs pass only a new pool option; direct filling constructs another small default policy. A strict sheet policy therefore does not constrain the rest of its catalog. **Fix:** pass one immutable certification policy through all producers and store its ID in every certificate. **Test:** tighten a gate sufficiently to reject a fixture and verify identical rejection through crossing, seed, rib, and direct-harvest routes.

- **[CORRECTNESS]** `rib_validate.m:59–98` — **Publication validation does not establish the advertised certificate contract.** It checks that `.z` exists, not that it is a real finite eight-vector; it does not check each point’s `.sA` against the rib’s column, or require full-stack evidence. `NaN` problem values can evade the `abs(a-b)>tol` rejection because that comparison is false for `NaN`. Coverage is count-based: including the spine and omitting another target can still have the requested count. **Fix:** validate scalar types/domains first, validate complete point coordinates and root shape, and compare the exact set of target indices. **Tests:** `NaN` identity, `.ok=2`, malformed `z`, wrong point arrival phase, and a same-count/wrong-target rib.

- **[ROBUSTNESS]** `rib_validate.m:63–64`; `rib_from_crossing.m:121–124`; `run_costate_library.m:446–459` — **A genuine first-step stall cannot be published.** The walker returns a terminal empty rib, but validation rejects it. Repeated attempts can retire the column and prevent packaging the valid spine and other columns—the very holes intended for direct filling. **Fix:** distinguish a valid terminal “zero targets certified” result from malformed/infrastructure failure; report its zero coverage without treating it as missing. **Test:** a deterministic first-step certification refusal must yield a short published unit and a packageable catalog with holes.

### Packaging and durable provenance

- **[CORRECTNESS]** `package_phase_catalog.m:53–60,103–112`; `sheet_to_catalog_file.m:185–197` — **The packager accepts identity-free legacy ribs “on trust,” and accepts certification from a Boolean plus a finite vector.** Extra rib files bypass the work-unit validator. The compared identity also omits `muStar`, unit scales, endpoint construction, and gate policy. Consequently, a correctly shaped old artifact can enter a new catalog without satisfying its current certification contract. **Fix:** fail closed on missing identity; provide an explicit migration/re-certification tool. Require matching model, endpoint, and certificate-policy IDs. **Test:** offer a finite `.ok=true` legacy root with no identity/full-stack evidence; packaging must refuse it.

- **[CORRECTNESS]** `sheet_to_catalog_file.m:87–103,129–132,159–162,232–240` — **The phase-sheet export drops the junction states and most certificate evidence.** It exports `z8`, coarse conjugate verdict, three hypothesis numbers, and junction count, but not `Y`, H6, lift margin, pointwise checks, dense-scan coverage/status, or `fullStack`. Thus the shown interface cannot preserve the stated “z8 plus shooting junction states” deliverable or demonstrate which current gate stack each entry passed. **Fix:** export entry-indexed junction states and a versioned certificate record alongside each root. If a compact catalog intentionally excludes them, supply a content-addressed companion artifact and state that contract explicitly.

- **[CORRECTNESS]** `certify_root.m:127,188,232,252–254` — **Physical output units are not taken from the certified problem.** `m0kg` defaults to 150 independently of `B.problem.m0kg`; callers generally do not pass it. A non-150 kg problem can fly with the correct nondimensional propulsion yet report the wrong propellant and final masses. Length and time scales are likewise hardcoded rather than consumed from `B.problem`. **Fix:** derive all dimensional outputs and flight screens from the authoritative problem constants; treat options as assertions, not replacements. **Test:** a 300 kg fixture must report twice the dimensional mass for the same normalized mass fraction.

  The visible conversions themselves—`tf*tStar/86400`, velocity by `lStar/tStar*1000`, and `cnd*log(1/mf)*lStar/tStar`—are dimensionally correct. The direct solver also correctly converts the **1900 km lunar radius** screen to a **162.6 km altitude** constraint.

- **[CORRECTNESS]** `build_70mN_library.m:317–321,335–340,374–380`; `audit_phase_catalog.m:159–162` — **Saved audit/sweep evidence is not bound to catalog contents, and completion is conflated with success.** An audit-disabled chain accepts audit-file existence without reading its verdict or proving it covers the current roots. The sweep stage explicitly blocks incompleteness but does not explicitly block adverse `worstSpectrum`, H6, or lift results returned by a completed sweep. The unseen helpers may add safeguards, but this chain exposes none. **Fix:** require a catalog-content hash, complete entry census, certificate-policy ID, and aggregate pass verdict from both measurements. **Tests:** reuse a clean audit after changing one root; return a completed sweep with one adverse result; neither may leave the chain clean.

### Family attribution and direct filling

- **[CORRECTNESS]** `family_map.m:104–105,240–266,277–298`; `sheet_to_catalog_file.m:145–148` — **Family ownership is inferred from time, not root identity.** The implemented spine tolerance is **0.02 days**, not the methods document’s `1e-5` days; rib matching allows **0.03 days**. Distinct equal-time roots cannot be distinguished. Rib extrapolation uses coefficients valid only for the first three equally spaced offsets, while the pipeline now supports nonuniform grids. Direct-filled “ribs” have neither that spacing nor a spine lineage, yet receive the same inference. **Fix:** retain arc/root IDs and actual parent-root identity through continuation and ribs; mark direct roots unattached unless separately connected. For legacy reconstruction, compare full corrected roots and use actual phase offsets, returning “ambiguous” rather than nearest-time ownership. **Tests:** equal-time distinct roots, nonuniform ribs, and arbitrary direct-filled cells.

- **[CORRECTNESS]** `family_map.m:123–126,223–229` — **Family topology labels overstate the stored evidence.** Every fold candidate is treated as a fold regardless of `classified`/`resolved`; every ordinary terminal sample is labelled “walk budget … branch continues,” even when the stop was a stall. Taking `abs(rho)` also erases passage into the non-minimizing negative-multiplier chart. **Fix:** consume the actual fold verdicts, preserve signed `rho`, and classify stop reasons explicitly. **Test:** unresolved fold, Newton stall, and negative-`rho` fixtures must not become classified folds, budget continuations, or normal PMP branches.

- **[CORRECTNESS]** `fill_holes_direct.m:109–125,234–235,250–259` — **The direct filler can relabel old results as a new problem.** Resume loads `R` and `rec` without checking the saved `problem`, grid, root identity, or policy. Successful records mark cells done by old integer indices; old points are admitted without checking `.ok`; the next save stamps everything with the current `B.problem`. **Fix:** validate a full resume manifest before loading any result, map records by phase/root keys rather than indices, and never overwrite historical identity. Use an exclusive output lock and unique temporary names as well. **Test:** reuse the output filename for a changed grid, engine, or orbit branch; resume must refuse before admitting anything.

- **[CORRECTNESS]** `fill_holes_direct.m:118–119,132–155,220–224` — **A successful hole fill becomes permanently ineligible for later improvement in the same output file.** `done` is based on any prior success, not on current objective, seed availability, or current root. The improve pass queues such cells and then skips them. Already attempted cells also cannot be requeued when a better neighboring seed is discovered during the run. **Fix:** separate “filled” from “improvement attempt completed for these seed/root versions”; requeue when the relevant neighbor changes. **Test:** fill a cell, introduce a faster neighbor, rerun with `improveDays>0`, and verify that the cell is attempted again.

- **[ROBUSTNESS]** `arclength_arrival.m:147–148`; `crossings_from_arc.m:122`; `direct_cell_solve.m:61–68`; `costate_common/arclength_ms.m:151,327–374` — **Several expensive operations remain outside the hard-cap discipline.** Anchor polish, rescan correction, seed flight, and direct solve run synchronously. `maxCpuSec` is not a hard elapsed-time fence, and the arc deadline is checked only between outer continuation steps. The hand-driven long-arc calls also do not pass `partialFile` (`build_70mN_library.m:201–204`). **Fix:** fence these operations and publish bounded-stage heartbeats/checkpoints. **Test:** replace a propagator with a deliberately nonreturning evaluation and verify bounded failure with recoverable progress.

- **[ROBUSTNESS]** `build_arrival_sheet.m:98–101`; `family_map.m:101,304–306`; `audit_phase_catalog.m:59` — **Three MATLAB interface semantics are consequential.** A `for` loop iterates over the columns of its expression, so a column-cell `seedFiles` list causes this loop to process only its first file through `sf{1}`. Function arguments are evaluated before the call, so `pickTStar(S)` can throw even when `opts.tStarSec` is supplied, and `capped_pool()` is evaluated even when an audit pool is supplied. **Fix:** normalize seed lists with `(:).'`; use explicit `isfield` branches for potentially failing/expensive defaults. **Tests:** column-cell seed lists and an explicit time scale on a legacy sheet without stored `tStar`.

# B. Q2 — Sufficiency of the local-optimality tests

## (i) Hypothesis/gate table

Here **VERIFIED** means that the supplied caller enforces the stated **numerical check**. It does not mean an exact continuous-time theorem hypothesis has been proved.

| BCT hypothesis / sufficiency condition | Gate that checks it | Verdict | File:line |
|---|---|---|---|
| Normal PMP lift with normalized cost multiplier, fixed initial state, fixed final \(r,v\), free final mass and time | Shooting residual; pointwise \(H\), adjoint, and \(\lambda_m(t_f)\) checks | **VERIFIED**, numerically; exact root existence remains unvalidated | `ms_tfmin.m:104–122,138–148`; `certify_root.m:217–231,275–293` |
| Control globally minimizes the instantaneous Hamiltonian over \(S^2\times[0,1]\) | Analytic full gap; acceleration- and mass-side throttle checks | **VERIFIED** at evaluated samples | `pmp_pointwise_checks.m:139–154`; `certify_root.m:285–293` |
| Strengthened Legendre condition on \(T_\alpha S^2\) | Sampled \(\min|\lambda_v|>0\) | **VERIFIED** only as sampled positivity | `mintime_hypothesis_gates.m:113–117`; `certify_root.m:380` |
| Uniform strengthened Legendre margin over all of \([0,t_f]\) | No between-sample lower enclosure | **NOT CHECKED** | Same locations |
| Strict all-burn complementarity \(Q_{mt}>0\) | Sampled minimum | **VERIFIED** only at samples | `mintime_hypothesis_gates.m:115–118`; `certify_root.m:381` |
| Reduction from the original throttle-constrained problem to a smooth direction-control second variation | H3 supports the critical-cone argument, but no complete reduction theorem is instantiated | **ASSUMED** for the sufficiency claim | Audit `sec:reductions`, `sec:hyp`; `ms_tfmin.m:118–122` |
| No abnormal lift of the complete nominal trajectory | Sampled lift-constraint rank, accepted-lift residual, two-tolerance margin | **VERIFIED** as a numerical rank assessment, not an enclosed rank theorem | `mintime_hypothesis_gates.m:144–170`; `lift_margin.m:77–117`; `certify_root.m:382–410` |
| Any stronger normality/corank hypothesis on proper subintervals required by the chosen theorem | No prefix/subinterval test | **NOT CHECKED** | `mintime_hypothesis_gates.m:148–162` |
| Smooth flow, positive mass, avoidance of primary singularities throughout the interval | Delegated flight validation | **NOT CHECKED** as a continuous enclosure in the supplied code; helper unavailable | `certify_root.m:244–245,322–323` |
| Correct free-time, free-mass accessory problem and equivalence to the reduced determinant | Dense code explicitly forms \([\Phi_{rv}P,F_{rv}]\); theorem-equivalence argument remains incomplete | **ASSUMED** for sufficiency | `conj_spectrum.m:106–107,136–137`; audit `sec:equiv` and its repaired derivation |
| H6 excludes the reduced-Hamiltonian zero mechanism | H6 flag and strict margin enforced | **VERIFIED** as reported numerical inequalities | `mintime_hypothesis_gates.m:181–183`; `certify_root.m:418–433` |
| No conjugate/focal degeneracy anywhere in \((0,t_f]\), including startup and endpoint | Junction test plus candidate-based dense scan | **NOT CHECKED** continuously; startup coverage is not enforced | `certify_root.m:453–475`; `conj_spectrum.m:19–49` |
| Endpoint transversality and endpoint-curvature terms if either phase is an optimization variable | No phase stationarity equations or corresponding Jacobi boundary conditions | **NOT CHECKED**; unnecessary only for genuinely fixed-phase cells | `ms_tfmin.m:112,145–148`; `certify_root.m:195` |
| Independent reproduction of the extremal | Separate shooting solver, but shared pumpkyn dynamics and propagation | **VERIFIED** as solver agreement, **not** an optimality condition | `certify_root.m:298–331`; `ms_tfmin.m:14–20` |

## (ii) Findings

- **[CORRECTNESS]** `costate_common/ms_tfmin.m:112,138–148`; `arclength_arrival.m:165–166`; `certify_root.m:195` — **The solved and tested problem has fixed endpoint points, not freely selectable points on two orbits.** Continuation in phase changes the problem parameter; it does not put phase variations into a cell’s BVP. This is correct for a library entry claiming a minimum **conditional on its two specified phases**.

  If either phase is free in the claimed optimization problem, additional necessary conditions include
  \[
  \lambda_{rv}(0)^\top x_D'(s_D)=0,\qquad
  \lambda_{rv}(t_f)^\top x_A'(s_A)=0
  \]
  for the respective free phases. The accessory boundary conditions must admit the associated endpoint tangent variations and include endpoint-curvature terms. A point-target determinant does not supply that focal/manifold test. **Fix:** stamp the problem class as `fixed-phase-point-endpoints`; build a separate free-phase BVP and accessory problem when that is intended. **Test:** finite-difference the optimized time over neighboring phases and compare with endpoint costate sensitivities.

- **[CORRECTNESS]** `costate_common/conj_spectrum.m:106–107,136–137`; `costate_common/ms_tfmin.m:118–122`; audit `sec:equiv`, `prop:equiv`, repaired permanent-kernel derivation, and `sec:hyp` — **The rank identity is not yet a sufficient-optimality proof.** The document correctly withdraws its original non-immersion argument and identifies the permanent kernel. Quotienting that kernel establishes a square endpoint-map rank test. It does **not by itself** establish that this map’s rank losses are exactly the degeneracies of the constrained free-time/free-mass second variation, or that its index begins positive.

  The citation “BCT 2007 Thm. 3.x” is not an identifiable theorem statement. The **Bonnard–Caillau–Trélat smooth second-order theory** must be invoked with its actual normality/regularity, strengthened Legendre, endpoint, and no-conjugate-point hypotheses for the actual accessory problem. **Fix:** supply an exact theorem citation and a derivation matching this problem to it; independently construct the mixed endpoint accessory BVP and compare its nullity/index with the reduced determinant. Do not replace this missing link with agreement on one anchor.

- **[CORRECTNESS]** `costate_common/ms_conjugate_test.m:26–31`; `costate_common/ms_tfmin.m:121–122` — **The production junction test’s variational specification is not auditable from the supplied source.** The wrapper forwards to unseen `oc` code, and `ms_tfmin` supplies only a flow callback. The claimed state rows, costate columns, quotient direction, and free-time mode therefore depend on unseen defaults—especially significant after promotion to a problem-agnostic library. **Fix:** pass the full intended specification explicitly using the actual `oc` API, and save the resolved specification with the verdict. **Test:** inspect/assert the resolved configuration and run fixtures that distinguish fixed-time, free-time, and wrong-row tests.

- **[CORRECTNESS]** `costate_common/conj_spectrum.m:19–49`; `certify_root.m:453–475` — **A clear dense scan does not establish H5 on the whole interval.** It remains candidate-triggered finite sampling, with policy floors rather than bounds on the STM error. Narrow zero pairs and even-order zeros can evade candidate detection. More directly, `tUncovered` is copied into the certificate but never blocks acceptance. The original document’s `sec:sampling` caveat still applies; sampling at integrator steps would not eliminate it either. **Fix:** report an assessment rather than a proved certificate; for proof, cover the startup interval analytically or by validated asymptotics and enclose the relevant variational quantities on every remaining interval.

- **[CORRECTNESS]** `certify_root.m:453–469` — **The advertised multiplicity gate is not independently enforced.** `multiplicity` is validated as a count, but is absent from both the consistency expression and the acceptance condition. A result with `multiplicity>0`, zero other adverse counts, and `clear=true` passes this consumer. Whether the unseen resolver currently prevents that combination is a separate contract question. **Fix:** enforce the intended adverse condition explicitly, or rename/remove it if it is merely a diagnostic corank. **Test:** inject that exact combination through the existing `CS` override and require refusal.

- **[CORRECTNESS]** `costate_common/conj_spectrum.m:40–48`; `certify_root.m:436–440` — **A numerical floor hit or numerical corank is not a demonstrated conjugate zero or its multiplicity.** Refusing certification is appropriate. Calling every such result a refutation of local minimality is not. In particular, a degeneracy at \(t_f\) defeats this strict sufficient condition but does not, by itself, prove the trajectory is not a minimum. **Fix:** distinguish `numerically clear`, `unresolved`, `validated interior conjugate point`, and `endpoint degeneracy`; reserve refutation language for a justified necessary-condition violation.

- **[CORRECTNESS]** `costate_common/mintime_hypothesis_gates.m:110–118,145–154`; `costate_common/lift_margin.m:32–38,77–102` — **The lift margin estimates only one error component.** Both constraint builds use the same pumpkyn flight and the same interpolated state/control data; changing the adjoint integration tolerance does not bound their common errors. The helper’s own caveat is correct. The resulting finite matrix can support a useful rank assessment, but not a rigorous exclusion of continuous abnormal lifts. Also, the header’s “no abnormal lift iff dim S=1” at `mintime_hypothesis_gates.m:19` contradicts the corrected document: the converse fails without the PMP cone inequalities. **Fix:** use only the sufficient exclusion direction; independently refine the underlying flight/interpolation and ultimately bound total constraint-matrix error. Test proper subintervals if the selected theorem requires strong normality there.

- **[CORRECTNESS]** `certify_root.m:380–381`; `costate_common/mintime_hypothesis_gates.m:113–118`; audit `eq:Qmt`, `sec:reductions`, `sec:hyp` — **Strict sampled positivity is weaker than the uniform margins needed for the smooth/strict-bang argument.** Arbitrarily small positive sampled values pass, even below their numerical uncertainty.

  There is **not** a missing ordinary ambient \(H_{\alpha\alpha}\) test: \(H\) is linear in \(\alpha\), and the relevant constrained/Riemannian Hessian on the sphere is
  \[
  \left.H_{\alpha\alpha}\right|_{T S^2}
  =\frac{T|\lambda_v|}{m}I.
  \]
  Likewise, if uniform \(Q_{mt}>0\) is established, there are no throttle switches to regularize and no switching-time Hessian is automatically required. The exact gap is
  \[
  H(s,\alpha)-H(1,\alpha^*)
  =TQ_{mt}(1-s)
   +\frac{sT|\lambda_v|}{2m}\|\alpha-\alpha^*\|^2.
  \]
  This gives the needed quantitative strict-bang/Weierstrass information once its coefficients have uniform lower bounds. **Fix:** compute lower bounds above numerical uncertainty and explicitly connect the strict-bang critical cone, \(\delta s=0\), to the smooth direction-control second variation. If \(Q_{mt}\) reaches zero, refuse this certificate class and use the appropriate switching/singular-arc theory.

- **[CORRECTNESS]** `audit_phase_catalog.m:133–149` — **The audit fails open on the very gates it claims to repeat.** A failed or timed-out polish leaves `conjNow=NaN`; all the rejection comparisons then skip it. A failed hypothesis call similarly leaves `dimSNow=NaN`. Successful hypothesis output is not checked for positive `minLamV`, positive `minQmt`, H6, or lift margin. The dense scan and pointwise PMP stack are not repeated. Moreover, the polished root is discarded, so a conjugate result can refer to a moved root without checking that movement. **Fix:** require successful, finite, matching-root results from the same production certification policy, with no missing-gate pass path. **Tests:** force each external call to timeout, return malformed data, or return a different polished root; all must count as bad audit rows.

- **[CORRECTNESS]** `certify_root.m:126,244–271`; `audit_phase_catalog.m:47–58,127–131` — **The claimed flight accuracy does not match the acceptance implementation.** The production certifier gates the loose pumpkyn flight at 100 km and 10 m/s; the tight flight is used only for pointwise checks, not the arrival gate. It is not passed through `validate_flight`. The audit adds a 1 km screen but validates the witness flight only in position, not velocity, completion, mass law, or terminal PMP conditions. **Fix:** use one explicit flight-validation contract for the root and witness, including completed tight-flight endpoint reproduction, and store the actual numerical tolerances. A tighter flight should be validated before any terminal quantity is read from its last row.

- **[ROBUSTNESS]** `costate_common/pmp_pointwise_checks.m:125–137,156`; `certify_root.m:277–289`; `costate_common/conj_spectrum.m:146–151` — **Some independent numerical diagnostics are computed but not consumed.** `fdAgree` does not qualify the finite-difference adjoint check; `specConsistency` is not read by the certifier. Neither is itself an optimality theorem, and the independent adjoint field check provides additional protection, but their purpose is defeated if arbitrarily bad values are accepted. **Fix:** enforce sensible finite, uncertainty-scaled bounds or explicitly mark them advisory and remove them from any certificate claim.

- **[CORRECTNESS]** `certify_root.m:169–170,272,354,405,450,486` — **The test override can still produce a production-shippable certificate.** A warning is the only separation between injected instrument outputs and `C.ok=true`; the resulting certificate does not record that it was overridden. **Fix:** any nonempty override must force diagnostic-only status, or move the mutation seam into a test-only dependency injector whose results cannot be packaged.

### Can the stack pass a saddle, maximum, or non-strict minimum?

**Yes, the implemented stack cannot exclude those possibilities in general.** A PMP-minimizing control law and positive sphere curvature exclude a simple direction-wise Hamiltonian maximum; they do not establish positivity of the whole constrained trajectory second variation.

A saddle can survive if the relevant conjugacy is missed, startup coverage is omitted, or the wrong endpoint variation space is tested. A non-strict extremal can survive an unresolved exact degeneracy being treated numerically as nonsingular. Free-phase competitors are not tested at all.

Conversely, **under a correctly applicable sufficient theorem, uniform strengthened conditions, the correct endpoint accessory problem, and genuinely established absence of conjugate/focal degeneracy, such counterexamples would be excluded.** The problem is the gap between that conditional statement and these finite numerical measurements—not that strict local minima are unattainable in principle.

## (iii) Prioritized tests to add

1. **Fail-closed end-to-end gate mutation suite.**  
   Force every audit/certifier stage to fail, timeout, omit fields, return `NaN`, or return a different root. Include `multiplicity>0` with otherwise clear counts, nonempty override, and uncovered startup time.  
   **Failure means:** an uncertified or differently certified root can be advertised as passing.

2. **Independent accessory-problem and inertia test.**  
   Form the second variation with fixed initial state, terminal \(r,v\) fixed, terminal mass free, and free final time. Include the linearized terminal mass-costate and time conditions. Compute its reduced Hessian/inertia or an equivalent Riccati/Jacobi construction; compare degeneracies with \([\Phi_{rv}P,F_{rv}]\). Add phase tangent and curvature terms only for a separate free-phase problem.  
   **Failure means:** the reduced determinant is not testing the claimed variational problem, or an extremal has a negative/zero critical direction.

3. **Variational-flow validation independent of the shared STM implementation.**  
   Compare STM/Jacobian actions against independent AD or accurately scaled perturbations; check the mass-costate zero column, mass-row invariance, costate-scaling kernel, and Hamiltonian-flow identities. Use fixtures with an early conjugate point, a pair inside one coarse interval, an even-multiplicity degeneracy, and an endpoint degeneracy.  
   **Failure means:** second-order verdicts may be artifacts of the variational propagation.

4. **Continuous margin and coverage assessment, followed by validated bounds where proof is required.**  
   Locate minima of \(|\lambda_v|\), \(Q_{mt}\), mass, and primary distance using dense output and refinement; bound interpolation/integration uncertainty. Cover the initially singular determinant interval with a justified small-time analysis rather than a skip.  
   **Failure means:** the smooth strict-bang theorem’s hypotheses are unresolved or violated.

5. **Total-error lift-rank test.**  
   Refine the underlying flight, interpolation, adjoint integration, and sample set independently. For rigorous exclusion, bound the error of a suitable finite sampled constraint matrix and require its sixth singular value to exceed that bound, together with an exact/validated normal lift.  
   **Failure means:** abnormal-lift exclusion is unresolved—not necessarily that an abnormal minimizing lift exists.

6. **Mass-costate consistency identity.**  
   Test the document’s corrected identity, allowing terminal residuals:
   \[
   \int_0^{t_f}TQ_{mt}\,dt
   =m(0)\lambda_m(0)-m(t_f)\lambda_m(t_f).
   \]
   Also test \(\dot\lambda_m=-T|\lambda_v|/m^2\), monotonicity, and the reduced-Hamiltonian identity used by H6.  
   **Failure means:** a sign, mass normalization, transversality, or integration inconsistency.

7. **Validated root existence for a genuinely mathematical certificate.**  
   Apply interval Newton/Krawczyk or a comparable validated BVP method around the polished shooting solution, then validate the required second-order bounds.  
   **Failure means:** the numerical approximation has not been converted into a theorem about an exact admissible extremal. Failure to validate is not itself a refutation.

# C. Q3 — Generalization to any orbit pair

## DRO/tulip-specific inventory

| Location | Specific assumption or blocker |
|---|---|
| `build_70mN_library.m:60–78,141–155,249,360` | Hardcoded engine/tag, DRO/tulip descriptors, named anchors/families, filenames, and legacy period/petal-based picture path |
| `run_costate_library.m:173–201,195,520` | Only `tauDRO/NpTulip/pmTulip`; tulip-period formula; hardcoded `70mN` and `dro_tulip` artifact names |
| `arclength_arrival.m:63–86,101–106` | `ladder_endpoints` and a DRO/tulip-specific `B.problem`; periodic endpoint closures; tulip period locked by `Np` |
| `arclength_arrival.m:112–155,186` | Setup and arc operation coupled to a default campaign anchor |
| `build_arrival_sheet.m:94–96,124–127` | Default seeds from `dro_tulip_library`; saved endpoint construction not a portable, explicit orbit-pair definition |
| `build_ribs.m:55–59,66–73`; `rib_validate.m:89–94` | Reconstruction/identity fields enumerate only DRO and tulip parameters |
| `sheet_to_catalog_file.m:174–179` | Literal departure family `dro`, arrival family `tulip`, and legacy period fields |
| `package_phase_catalog.m:47,83,89–100` | Literal catalog name, `tau1_Np7` filename, description, and reconstruction instructions |
| `audit_phase_catalog.m:71–76` | Reconstructs only the legacy pair through `ladder_endpoints`; cannot independently consume an arbitrary pair descriptor |
| `family_map.m:118–123` | Family/anchor identity inferred from an `arrival_arc_*_(dn|up)` filename and assumed homogeneous-chart layout |
| `fill_holes_direct.m:73–74,97–100`; `direct_cell_solve.m:47–48,66` | Legacy direct-code location/name and legacy pair reconstruction; the numerical solve call itself takes arbitrary endpoint states |
| `certify_root.m:127,188` | Campaign-specific mass/unit defaults rather than authoritative model constants |
| `costate_common/get_family_orbit.m:76–129` | Implements DRO, Tulip, Halo, DPO, Lyapunov, GTO—not all nine documented families |
| `costate_common/survey_family_bounds.m:42–48,61–72` | Fixed Earth–Moon constants, lunar-vicinity criteria, and a periodic-orbit interpretation unsuitable as a universal admissibility policy |
| `costate_common/ms_tfmin_hom.m:10–15`; `ms_tfmin.m:92–96` | Legacy DRO names/comments, but the numerical problem definitions are endpoint-agnostic |
| `costate_common/harvest_ms_seed.m:62–72` | Not DRO-specific; it assumes the direct solver’s HS/Sundman output schema |

The names in generic numerical helpers are mostly migration debt, **not mathematical orbit restrictions**. Conversely, `ladder_endpoints`, legacy identity fields, default-anchor loading, and periodic phase handling are functional restrictions.

### Missing catalog families

`get_family_orbit.m:76–129` lacks **Pumpkin, LPO, Axial, and Cycler**, despite the problem-space document listing them. Add their getters, resolved parameter metadata, refinement, and validation. Pumpkin’s fixed \(2\pi\) period and Cycler’s resonance descriptor belong in their provider cases—not in continuation, certification, or packaging.

## What GTO changes—and what it does not

`get_family_orbit.m:89–127` constructs an **algebraic closed locus of rotating-frame states at fixed orientation**, using `fromPCI(0,...)` for every anomaly. It is explicitly not a propagated rotating-frame GTO trajectory.

That distinction yields two valid but different applications:

1. **Fixed-state/locus transfer library:** choose two endpoint states from such loci and solve the autonomous CR3BP transfer between them. This is compatible with the numerical shooting problem. A closed-locus parameter may legitimately wrap modulo one, but its scale must not be called a rotating-frame orbital period or an epoch offset. Differentiating the actual endpoint interpolant remains correct.

2. **Rendezvous with a physically evolving GTO:** the target must depend on epoch. In general,
   \[
   x_A(s+1)\ne x_A(s)
   \]
   when \(s\) advances one Kepler period along the physical rotating-frame trajectory. Periodic seam forcing, integer-shift crossing searches, toroidal neighbor selection, and mod-one identity matching are then wrong.

For moving-target rendezvous, the endpoint constraint is of the form
\[
x_{rv}(t_f)=x_A(t_0+t_f;\theta),
\]
not a fixed `rvf`. The free-time condition correspondingly includes target motion:
\[
H(t_f)-\lambda_{rv}(t_f)^\top \dot x_A(t_0+t_f;\theta)=0.
\]
The residual Jacobian and second variation must include this dependence.

**A two-closure periodic-orbit abstraction is sufficient for arbitrary fixed endpoint loci, but not for arbitrary physical moving-target transfers.** Do not conceal that distinction behind a fabricated GTO “period.”

Also, `survey_family_bounds`’ 100,000 km lunar-vicinity ceiling excludes Earth-side GTO members by policy. Closure alone would not prove that a closed algebraic locus follows the CR3BP flow; flow consistency must be a separate provider validation.

## Minimal abstraction

Use one pair object, with the existing `B` as its numerical-problem adapter:

```matlab
pair.dep = struct( ...
    'name',       'dro', ...
    'params',     depParams, ...
    'state',      stateD, ...          % 6x1 rotating-frame ND
    'dstate',     dstateD, ...         % derivative of that exact endpoint map
    'periodND',   periodD, ...
    'kind',       'periodic-flow', ... % or 'closed-locus', 'epoch-path'
    'wrap',       true);

pair.arr = ...;

pair.model = struct( ...
    'name', 'Earth-Moon CR3BP', ...
    'muStar', mu, 'lStar_km', lStar, 'tStar_s', tStar, ...
    'frame', 'synodic-barycentric');

pair.identity = ...;                  % resolved descriptors + endpoint-rule hash
```

For the minimal fixed-phase generator:

- `state(phase)`, period/parameter scale, and name are the essential endpoint data.
- `dstate` is needed by arrival continuation; compute it from the **same** stored representation as `state`.
- `kind` and `wrap` prevent a closed locus from being confused with a periodic flow.
- Add `d2state` when implementing free-phase endpoint second variations.
- Keep propulsion, phase choices, and certification policy in a transfer-problem object containing `pair`; they are not properties of an orbit.
- Save reconstruction descriptors and/or endpoint interpolation data, not only anonymous function handles. Bind certificates to that saved representation.

### Concrete signature and field changes

| File | Minimal change |
|---|---|
| `get_family_orbit` | Keep its three outputs; add optional `model` input and return resolved `kind`, `wrap`, actual period/parameter scale, and complete parameters in `info`. |
| `survey_family_bounds` | Add explicit `model` and `criteria` inputs; apply periodic-flow checks only to that endpoint kind. Return a result row for every requested member, including failures. |
| `run_costate_library` | `run_costate_library(problem, opts)`; grids remain in `opts`, while pair, propulsion, units, and policy come from `problem`. |
| `build_70mN_library` | Retain only as a thin preset constructing the 70 mN problem and calling the generic front door. Remove duplicated physics/grid plumbing. |
| `arclength_arrival` | Split setup into `arrival_setup(problem, root, sD, sA, opts)`; use `arclength_arrival(B, anc, opts)` for walking. No implicit anchor reload. |
| `arclength_ms` | No orbit-related signature change. Keep it parameter-agnostic. Persist localized event states for shared crossing extraction. |
| `ms_tfmin`, `ms_tfmin_hom` | No orbit-related signature change for fixed-phase points. Their current numeric endpoint interface is already appropriate. |
| `crossings_from_arc` | Existing signature can remain; obtain wrapping policy from `B.problem.pair.arr`, and require matching arc identity/chart. |
| `certify_crossing` | Existing signature can remain; use the supplied problem/chart, exact grid endpoint, and policy. |
| `sheet_from_arcs` | Existing numerical assembly can remain; require problem identity, exact phase keys, root IDs, and explicit winner indices. |
| `build_arrival_sheet` | `build_arrival_sheet(problem, opts)`; explicit arc/seed lists; no default family library or anchor dependency. |
| `rib_from_crossing` | Existing signature can remain through `B`; preserve unwrapped departure coordinates and add spine-root identity to checkpoints. |
| `build_ribs` | Reconstruct `problem` from the sheet descriptor, without anchor polishing; inherit the sheet’s policy. |
| `rib_targets` | Add endpoint-domain/wrapping policy; for nonperiodic domains use directed ordinary differences rather than mod-one differences. |
| `rib_validate` | Replace enumerated legacy physics fields with complete problem, grid, policy, and spine-root identities. |
| `certify_root` | Existing `B` argument can remain, but require `B.problem` for units, propulsion identity, endpoint class, and policy; remove dimensional defaults. |
| `sheet_to_catalog_file` | Emit generic departure/arrival descriptors, junction states, and versioned certificate records; legacy fields become optional compatibility aliases. |
| `package_phase_catalog` | Derive filenames/descriptions from pair identity; never reconstruct family identity from filenames. |
| `audit_phase_catalog` | Rebuild endpoints through the same versioned provider contract from catalog descriptors; validate the full stored certificate policy. |
| `family_map` | Use saved continuation/root IDs and lineage, not anchor-name regexes or time extrapolation. No physical orbit family needs special treatment. |
| `fill_holes_direct` | Reconstruct the generic problem from the catalog; topology-aware neighbors and identity-bound resume. |
| `direct_cell_solve` | Existing endpoint/`B` interface is nearly sufficient; rename/wrap the direct solver and forward the shared policy, units, and clearance requirements. |
| `harvest_ms_seed` | No orbit-pair change required; retain and expose dual-mapping diagnostics as seed diagnostics. They are not a replacement for subsequent certification. |

## Refactor order by effort

1. **Small:** introduce the pair/provider adapter, resolved identity, shared phase matching, and dimensional constants; remove family literals from metadata.
2. **Small–medium:** decouple physics/chart setup from anchor loading. This immediately enables single transfers, non-default phases, and alternate propulsion.
3. **Medium:** thread the shared problem/policy through sheets, ribs, certifier, direct filling, and audit; add the missing catalog families.
4. **Medium–large:** migrate artifact schemas and resume keys; preserve roots, junctions, certificate evidence, and lineage. This is essential for reliable reuse.
5. **Medium:** replace generated MATLAB source specifications with serialized campaign specifications and configurable dependency startup.
6. **Large, separate mathematical feature:** free-phase orbit-to-orbit optimality, including endpoint transversality and focal/accessory boundary conditions.
7. **Largest, separate problem class:** physical epoch-dependent GTO rendezvous and other nonperiodic/nonautonomous endpoints. This is not solved by renaming `stateD/stateA`.

# D. Ranked top five by consequence

1. **Strict-local-optimality overclaim:** continuous H5/startup coverage and the exact theorem/accessory-problem connection are not established. Existing passes are numerical assessments, not the claimed proof.
2. **Fail-open catalog audit:** failed second-order/hypothesis calls can produce clean audit rows, and the audit does not repeat the current full gate stack.
3. **Unbound reuse and provenance:** changed spines/policies can reuse old ribs, and direct-fill resume can stamp old results with a new problem identity.
4. **Loss of fold-adjacent candidates during rescanning:** rebuilding a sheet can discard roots already found by continuation and select a slower “fastest certified” candidate.
5. **Advertised generator configurations do not execute faithfully:** the default 24-column phase list fails validation, the single-departure case cannot package, and non-default problems remain coupled to the 70 mN anchor and legacy orbit plumbing.