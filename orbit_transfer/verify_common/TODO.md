# verify_common — TODO

Companion to `README.md`. Structure follows the campaign TODOs (Done / Open by
priority / Not-a-goal). Per-row results and the cross-campaign status live in
`../OPTIMALITY_CERTIFICATION.md`, not here — this file tracks work on the
**machinery**.

---

## Done

### 2026-07-25 — the layer, built and wired (plan `2026-07-25-foc-gate-layer.md`)
`foc_check` / `foc_manifest` / `foc_dual_to_costate` / `foc_ipopt_inertia` /
`foc_report`, plus wiring into all four campaigns and the explainer
`doc/first_order_checks.tex`. Coverage achieved on first run: Earth 2-body
9/9 rows (7/9 advisory pass), Earth CR3BP 4-row subset (3/4), tulip flagship,
ELFO three artifacts — **the ELFO campaign's first first-order gate of any
kind**. LEAD-0 (IPOPT `δ_w` inertia) ported campaign-wide in the same pass.
Report-only burn-in throughout.

### 2026-07-25 (Task 0, mesh-convergence-study) — mapped terminal covector, item 1.3 fully closed
Built the Hager (2000) mapped/transformed terminal covector (`rep.lamMassEndMapped`,
`foc_check.m` §(6)) and made it the gated transversality quantity, superseding
the I5 linear-extrapolation partial fix (kept as a reported-only companion
alongside the legacy raw one-sided dual). Re-ran earth 10/5/2.5 N and the
tulip flagship: all four land at KKT-stationarity floor (1e-17 .. 1e-25),
**retracting** the earth-2.5N transversality finding that the I5 fix could
not resolve — see `OPTIMALITY_CERTIFICATION.md` §A6 RETRACTION block. Test:
`tests/test_foc_terminal_covector.m`.

---

## Open — ranked by priority

### 1. Pre-promotion checklist (blocks turning burn-in into a hard gate)

These are the register's §A6 items, now **specified** by the 2026-07-25
three-way external review (GPT-5.6-terra + Gemini 3.1 Pro + host Claude;
verbatim reviews archived in `doc/review_2026-07-25_*.md`). Until they are
settled, a FAIL from this layer is a lead, not a verdict — and the layer
cannot be promoted.

**Outcome of the 2026-07-25 fix pass (items 1.1, 1.3, 1.4 — shipped, commit
`6b9c4a7`).** The fixes were built to *attribute* the outstanding advisory
FAILs, not merely move them, so each check now reports the corrected value
beside the legacy one:

| finding | legacy | corrected | attribution |
|---|---|---|---|
| tulip flagship Ṡ | 1.376e-07 | **2.701e+01** | **mesh artifact** — switches are transversal; flagship now advisory PASS |
| earth 5 N transversality | 3.189e-03 | **5.215e-04** | **endpoint bias** — now PASS |
| earth 2.5 N transversality | 1.265e-03 | **1.265e-03** | **REAL** — unchanged, so not the endpoint bias |
| ELFO energy seed (ε=1) | PASS incl. "sign law 100%" | PASS, 3 checks `--` | false PASS removed |

The 2.5 N result is the informative one: the endpoint bias only appears when
λ_m is still moving at t_f (a final *burn* arc), so a row whose final arc
coasts is untouched by the fix — which is exactly what happened. Its miss is
therefore a genuine property of that solution, pointing back at the
under-optimized-rung explanation (it is one of the rows warm re-solves
improve). `test_foc_mesh_invariance` locks 1.1 in both directions: across a 4×
refinement the normalized statistic holds (ratio 1.02) while the legacy one
falls by 4.1 — i.e. it was reporting `h`.

**What the review validated** (so it does not get re-litigated): the
switching-function recovery of Q1 is *exact* — the active bound multiplier
cancels `Sd` identically, so zeroing its row recovers magnitude and sign —
subject to the affine-throttle-cost conditions now in item 1.4; and the
dual→costate map is the correct *stationarity* combination (both reviewers
derived it, independently reproducing `DESIGN_dual_map`).

- [x] **1.1 Mesh-normalize the regularity statistic (register I1). DONE 2026-07-25.** Both
  reviewers independently confirmed the confound is real and that the current
  readings (tulip 1.4e-7 / 25 sw, ELFO front row 4.5e-5 / 50 sw) **cannot
  diagnose grazing**. Agreed fix — deweight, true divided difference,
  non-dimensionalize:
  ```matlab
  h = diff(sigma(:).');
  w = [h(1)/2, (h(1:end-1)+h(2:end))/2, h(end)/2];   % nodal trapezoid weights
  Sun = rep.Sd ./ w;                                  % deweighted switching fn
  swI = find(diff(burn) ~= 0);
  D   = abs(Sun(swI+1) - Sun(swI)) ./ h(swI);         % |dS/dsigma|, mesh-invariant
  Sref = median(abs(Sun(coastNodesAwayFromSwitches))); % exclude 1-node switch nbhd
  R    = (sigma(end)-sigma(1)) * D / Sref;            % dimensionless
  rep.sdotMinRel = min(R);
  ```
  Apply the same deweighting to the singular-arc test (item 1.3). **Threshold
  caveat (GPT):** no threshold on a *single* mesh is defensible — assert
  regularity only if `R` stays above the floor on **two** meshes. Recalibrate
  the 1e-3 gate on deweighted data before trusting it.
- [ ] **1.2 Make the direction check independent and sign-aware (register I2).**
  **Correction to the previously recorded fix:** excluding the cone multiplier
  before projecting does *not* restore independence — the projector annihilates
  the radial term anyway, `P_b(q + 2*mu*b) = P_b q`. Genuine independence needs
  a *different* gradient: build `q` from **defect duals only** (zero all other
  multiplier rows) or from a separately reconstructed Hamiltonian, then compare
  against the solved direction. For the min-vs-max branch, test the cone
  multiplier's sign: with `L = L0 + mu*(b'b - 1)` and `q + 2*mu*b = 0`, the
  minimizer has `b'q < 0` ⟺ `mu > 0` (reverse if the registered row is
  `1 - b'b`). Report normalized anti-alignment `-b'q/||q||`, which is more
  informative than tangency alone.
- [x] **1.3 Resolve the transversality endpoint bias (register I5). DONE 2026-07-25**
  (extrapolation, cheap fix) **+ SUPERSEDED 2026-07-25 (Task 0, mesh-convergence-
  study, principled fix).** Both reviewers confirmed the mechanism and first
  converged on the cheap linear-extrapolation fix:
  ```matlab
  lam(:,end) = LamDef(:,N)   + h(N)/(h(N-1)+h(N)) * (LamDef(:,N)-LamDef(:,N-1));
  lam(:,1)   = LamDef(:,1)   - h(1)/(h(1)+h(2))   * (LamDef(:,2)-LamDef(:,1));
  ```
  still live in `foc_dual_to_costate.m` as `lamX`/`rep.lamMassEndRel`, reported
  but no longer gated. The **more principled alternative** GPT also named —
  the exact discrete endpoint/mapped covector (Hager 2000, "Runge-Kutta
  Methods in Optimal Control and the Transformed Adjoint System," Numer. Math.
  87, 247-282), `lamHat_f = (h_N/2)*L_x(N+1) + (I - (h_N/2)*f_x(N+1))' *
  Lambda_N`, whose mass component is **zero by terminal-node stationarity**
  when the final mass is genuinely free — is now BUILT and is what
  `foc_check.m` gates on (`rep.lamMassEndMapped`, section (6); implemented by
  reading the already-assembled full-Lagrangian gradient `gL` at the
  terminal-node mass index, no new AD). It resolved what the extrapolation
  fix could not: earth 2.5 N's transversality "finding" (1.265e-3, unchanged
  by the I5 extrapolation because the final arc coasts) reads 2.268e-18 under
  the mapped covector — RETRACTED, see `OPTIMALITY_CERTIFICATION.md` §A6
  RETRACTION block. Regression test:
  `verify_common/tests/test_foc_terminal_covector.m`.
- [x] **1.4 Make the verdict ε-aware. DONE 2026-07-25.** `foc_check` folds `signPct`, the
  singular-arc count, and `sdotMinRel` into `rep.pass` **without knowing ε**.
  At ε>0 the throttle cost is quadratic, `Sd` is the derivative of the smoothed
  objective, interior `s` can legitimately give `Sd = 0`, and the bang-bang sign
  law does not apply. Add a manifest field (`throttleCostKind = 'affine' |
  'quadratic'`) or pass ε, and skip those three checks programmatically. The doc
  already claims these are advisory at ε>0 — the code does not honor it.
- [ ] **1.5 Relabel the δ_w verdict.** `certLocalMin` / "LOCAL MIN" overstates
  what zero inertia correction over a barrier-iterate tail proves: it is not a
  demonstration of reduced-Hessian positive definiteness on the critical cone,
  and does not address barrier/active-set limiting behavior or strictness.
  Relabel as an inertia-regularization *observation*, campaign-wide (the naming
  is inherited from `psr_ipopt_certify.m`).

### 2. Coverage holes

- [ ] **CR3BP deep rungs (0.2 / 0.1 N) — NOT RUN.** Deferred at build time;
  their warm re-solves are ~30k-node problems. 2.5 N also skipped. The CR3BP
  row of the register's coverage matrix is a 4-row subset, stated as such.
- [ ] **ELFO beyond the nominal 25 mN rung.** `run_foc_elfo` hardcodes
  nominal-rung `cBox` / `tfCapMult` / `maxIter`; on another rung the warm
  re-solve is a *different problem* and will trip the certified-quantity guard
  opaquely. Thread the artifact's own `fp` config before any ELFO ladder work.
- [ ] **Batch ELFO runs must isolate per process.** A combined single-process
  sweep of three artifacts died silently mid-solve (this codebase's documented
  MEX-fatal risk); re-run isolated, it converged cleanly. Mirror
  `elfo_batch.sh`.

### 3. Robustness and hygiene

- [ ] **`foc_report` thresholds are hardcoded to `foc_check`'s defaults.** No
  caller overrides them today, so labels and `rep.pass` agree — but a caller
  passing custom `opts` would get per-line PASS/FAIL text that silently
  disagrees with the verdict. Fix: record tolerances in `rep.tol` and have
  `foc_report` read them. *(Flagged independently by GPT-5.6-terra.)*
- [ ] **Harden the sign-resolution rule.** Choosing the smaller of
  `||grad f ± A'lam||` is safe *at a converged KKT point* (the wrong branch
  evaluates to `2*grad f`, which is large) — but GPT gives a valid
  counterexample at a **non-converged** point (`g_f=1, A'lam=0.2` → the wrong
  sign wins), and the sign is chosen *before* `kktStatInf` is gated. Add an
  ambiguity guard: require the chosen residual to be both small **and**
  decisively smaller than its opposite, else error "ambiguous dual convention"
  rather than silently signing `rep.lam`. (Gemini judged the rule "100% safe";
  that holds only post-convergence — hence the guard, not a rewrite.)
- [ ] **Assert defect-row ordering, not just variable layout.** The layout
  asserts confirm the X/U blocks numerically, but `reshape(lamAll(defRows), nx,
  Nseg)` silently corrupts costates if a registry ever groups defect rows
  non-interval-major. Store an explicit `creg.defectRowsByInterval` (nx×Nseg)
  map and index through it. *(GPT.)*
- [ ] **Input validation before differentiation:** `numel(sigma)==N1`, strictly
  positive mesh steps, manifest indices in range, finite and nonzero direction
  norms before `b/norm(b)`. *(GPT.)*
- [ ] **`earth_mee` in `mintime` mode is unhandled.** That solver registers
  `thrEq` (thr ≡ 1) with no `thrLo`/`thrHi`, while the manifest still sets
  `thrRow = 4`, so `Sd` would retain the equality multiplier and `signPct`
  would be meaningless. Unreachable today (`refresh_duals_mee` forces
  `fixedtf`); add a `thrEq`-present guard so it fails loudly if reached.
- [ ] **`dirTanMax` returns `[]` on an all-coast trajectory** with `dirRows`
  set. Unreachable with pinned endpoints, and it errors rather than passing
  silently, but it deserves an explicit guard if a coast-only manifest ever
  appears.
- [ ] **`man.autonomous` is declarative only** — nothing reads it. Either
  consume it (gate Hamiltonian constancy where it is true) or keep it as
  documentation, but do not let a future reader assume it routes behaviour.
- [ ] **`run_gergaud` live-solve FOC branch** and the `run_psr` /
  `run_elfo_minfuel` / `gen_elfo_mintime` wiring are static-checked only, not
  exercised end-to-end (each costs a full pipeline run). All are try/catch
  advisory, so failure degrades to a warning.

### 4. Strengthening the checks themselves

- [ ] **Independence (register G1).** Every gate here consumes the NLP's *own*
  multipliers — self-consistency, not an independent witness. The tulip's
  least-squares costate reconstruction is a second family and already
  disagrees on the flagship (raw-dual 0.058°/100% vs LS 2.0/40%, the LS side
  documented-unreliable over ~40 revs). A genuinely independent witness needs
  costates from an indirect solve.
- [ ] **Free-t_f Hamiltonian condition, value form (register G4, partial).**
  The dual form is confirmed — the ELFO min-time anchor reads
  `lamTimeEnd = -1.000`, exactly `lambda_t(t_f) = ±1`. But its
  `lamTimeCoV = 5.69e-2` is **not** small for an autonomous model; open
  candidates are the two-primary Sundman + `cScale` structure distributing the
  horizon condition differently than the plain reading assumes, or
  under-converged duals on that anchor. Derive the `freetf-cscale` horizon
  condition properly and replace the placeholder `horizonNote`.
- [ ] **Non-autonomous Hamiltonian test (register G6).** Under lunar gravity
  `H` is genuinely not constant, and the layer correctly refuses to gate
  constancy — but then checks nothing. Test `dH/dt = ∂H/∂t` instead; the
  adjoint recursion for it is already inside `kktStatInf`.
- [ ] **Continuous-residual check on the costates (the unmeasured half).**
  Every dual-side check we run is *discrete*: KKT stationarity IS the discrete
  adjoint recursion, so verifying it is near-tautological given that we DEFINE
  the costate to be the multiplier. It confirms the solver converged and that
  the multipliers belong to the model we think we solved (which is exactly what
  caught the `opti.dual` bug) — but it is silent on whether the discrete
  adjoint approximates the continuous one.
  **The asymmetry worth acting on:** on the PRIMAL side that discrete-vs-
  continuous gap has been measured — `../GTO_tulip/direct/PSR/psr_switch_hessian.m`
  forward-integrates the exact direct control from x0 over ~40 revs and diverges
  by ‖r‖~3, ‖v‖~5 while the defects read 1e-14. (That number is a statement
  about IVP CONDITIONING over 40 revs — the same wall that defeats indirect
  single shooting — NOT evidence the collocation solution is a poor
  approximation of the continuous BVP solution; cross-formulation agreement
  says otherwise: MEE 1377.10 vs Cartesian 1376.74 kg at 10 N, both inside
  HMG-2004's published band.) On the DUAL side the analogous gap has never been
  measured at all.
  **Cheap experiment:** interpolate λ(σ) from `rep.lam`, differentiate the
  interpolant, evaluate −∂H/∂x independently along the trajectory (AD on the
  same RHS the solver used), and report the residual — as a mesh-refinement
  study, not a single-mesh number. It can fail where KKT stationarity cannot,
  because it asks whether the transcription's adjoint approximates the
  continuous one rather than whether it is self-consistent. Pairs naturally
  with the switch-time mesh-band work, since both are "how much of the fine
  structure is real?" questions.
- [ ] **Mesh bands on the gate values.** Switch counts are known
  mesh-sensitive; the gate numbers are reported as single values, not bands.

---

## Not a goal — intentional scope boundaries

- **Replacing the physical verifiers.** `verify_pmp_mee` and
  `certify_minfuel_pmp` carry the interpretability (primer angle in degrees,
  the switching function in RTN, the Fig-16 plots) that a dimensionless
  residual cannot. This layer is a **gate over** them, and their disagreement
  is a useful signal — keep both.
- **Second-order certification.** Belongs to Part B of the register. The one
  second-order line printed here (`δ_w`) is a ported verdict, not new
  machinery; the live path to *strict* local minimality is the STM
  switching-time Hessian, which is not this module's job.
- **Forking `foc_check` per campaign.** The entire point is one instrument.
  New campaigns get a manifest.
