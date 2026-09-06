# orbit_transfer — TODO (program level)

Per-campaign items live in each subfolder's `TODO.md`; this file holds only
what spans campaigns. The organized long-form view (case table, verification
ladder, ordered roadmap) is `STATUS_AND_ROADMAP.md`. Task numbering follows the session task list where one
exists.

## Costate-catalog program (Darin roadmap)

- [ ] **Densify the red cells** of all three tulip catalogs (halo 46 pairs,
  DPO 65 pairs, DRO s_A=0.075 row + 6 stragglers — tasks #14/#15). The
  engines revisit any no-OK cell automatically; run `densify_ladder`-style
  passes with fresh attempt budgets.
- [x] **L1↔L2 halo-to-halo catalog** — DONE 2026-08-09, BOTH directions,
  shipped as deliverable 6 on the schema-v2 arrival-period axis (A: 1,952
  entries / 90%; B: 2,096 / 96%; directions measurably asymmetric). Still
  open from it: the lone 10 N tfMin-acceptance-NaN watch item, and
  densification of the remaining pairs.
- [x] **Thrust-depth extension to 0.5 N — DONE 2026-09-01** (roadmap §6
  step 4): all four catalogs extended onto rungs `[0.75 0.5]` N by the
  family-agnostic `extend_thrust_ladder` (warm continuation from each
  full-depth cell). 1,348 eligible → 1,314 at 0.75 N (97.5%) → 1,150 at
  0.5 N (85.3%); 2,464 new entries; median 1.2 s/cell-rung, one hard-cap
  worker kill in 2,662 solves. Catalogs repackaged (18,360 entries,
  schema v2 clean), FULL conjugate re-sweep 18,249/61/50 with verdicts
  written back, golden cells 20/20. **Finding: 60 of the 61 refutations
  sit on the two new rungs** (2.4% of deep entries develop interior
  conjugate points as t_f grows — consult `conj_pass`, not just
  `has_solution`, below 1 N); the 61st is the known DPO 15 N refutation,
  reproduced. Deliverable zips now predate both the extension and the
  verdicts — re-ship when Darin wants them.
- [x] **Deep probe (0.1 N + 25 mN, one cell) — DONE 2026-09-01** (roadmap
  §6 step 4, second half): `DRO_tulip/indirect/probe_deep_rungs` walked
  the fine sheet's fastest 0.5 N cell down to **0.09 N** (six rungs, all
  tfMin-accepted at |dz| = 0) and located the closure wall at
  **0.09 → 0.067 N**. It is a BASIN/topology wall, not sensitivity: only
  ~1.3 swept revs even at 90 mN, ms stalls at normR ~ O(1). 0.1 N target
  CLOSES; 25 mN needs winding-aware continuation (recorded route).
  Junction states banked for all closed rungs per the identifiability
  rule. Full record: `DRO_tulip/FINDINGS.md` §17.
- [ ] **Sub-0.5 N rungs** deferred by Darin "until customers need it" —
  the probe (above) now bounds what is reachable: catalog-recipe
  continuation closes to ~0.09 N on a well-behaved cell; 25 mN requires
  the winding-aware route first.
- [ ] Min-fuel catalogs: explicitly NOT yet (Darin, Aug 2026).
- [ ] **When min-fuel work restarts anywhere in this repo, check MfMax first.**
  It solves this class already and is built + validated here (runbook:
  `earth_elliptic_to_geo/indirect/mfmax/MFMAX_V1_RUNBOOK.md`). Four things
  worth borrowing: HOMPACK arclength continuation for energy→fuel instead of
  our hand-scheduled ε ladder; the start-at-the-target initial-condition
  homotopy (`lambda(1)`) as a seeding trick we do not use; freeing an unknown
  scalar by carrying it as a zero-derivative state with a normalized terminal
  condition (v1 does this for `t_f`); and its PMP-converged costates as clean
  shooting seeds where direct-KKT duals are unreliable.
- [~] **Min-energy → MIN-FUEL pipeline (roadmap §6 step 5)** — pilot DONE
  2026-08-14; **γ grid DONE 2026-09-01** (15 records attempted, band
  [1.1, 1.4] reliable; high-γ direct fails + a measured basin split at
  (6,8)); **energy→fuel RACE DONE 2026-09-02** (FINDINGS §18: ε ladder
  walks to ε ~ 0.002 — but its "Huber refuted with mechanism" verdict was
  **RETRACTED 2026-09-05, FINDINGS §22**: the Huber STM lacked the
  saltation matrix; corrected, Huber walks all 17 rungs with 0 fails to
  m_f within 4.5e-6 of ε — ε remains the shipped convention);
  **FIRST MIN-FUEL RECORD SET DONE 2026-09-02** (FINDINGS §19: 7/7
  backbone records to ε ~ 0.001–0.005, fuel gain +0.13…+0.95% of m₀
  growing with γ, coast up to 58%, junction states banked;
  `run_minfuel_grid` / `minfuel_grid.mat`). **Learned: ss acceptance is
  not a valid gate at deep ε** (shooting fragility over 20–30 d near-bang
  arcs) — certification = ms + absolute-H + endpoints + the conjugate
  verdict. **FIXED-t_f CONJUGATE TEST DONE 2026-09-02** (FINDINGS §20: the
  existing instrument with the fixed-tf spec, validated on the analytic LQ
  π-conjugate case; 7/7 energy + 6/7 fuel PASS at the time — **corrected
  2026-09-05, FINDINGS §22**: the instrument monitored the terminal block
  `[1:6 14]` instead of the full state `1:7`, and the one "refutation" was
  an initial-coast structural zero; fixed instrument re-sweep = 15/15
  PASS; the second-order gate stays production policy). **SCHEMA V3 + RESCUE DONE
  2026-09-02 (STEP 5 CLOSED)**: catalog_schema v3 (objective/γ axis, one
  catalog per objective, named axis3, stored mf_frac + deltav_from_mf,
  mandatory Yj junctions for minfuel; TDD 12/12, v1/v2 compat clean); the
  refuted (2,5)@1.4 REPLACED via γ-continuation at fixed deep ε
  (m_f = 0.949005, conj PASS, γ-monotonicity restored — recipe: descend ε
  once on a good branch, then move γ; family jumps measurably fail);
  first v3 artifact `costate_catalog_dro_tulip_minfuel.mat` (7 entries,
  FINDINGS §21). Remaining: high-γ band widening only.
  - [x] **Code-review P0s FIXED 2026-09-05** (three-way review, `DRO_tulip/
    reviews/minfuel_code_review_2026-09-05.md`; FINDINGS §22): Huber
    saltation (`test_huber_saltation`), conjugate block 1:7 + initial-coast
    skip + t_f sampled + last bracket counted (`test_conj_fixedtf` +3),
    catalog Isp label derived from c_nd + validator consistency check
    (`test_catalog_schema_v3` +1); catalog rebuilt (Isp 900 / 0.07 N,
    conj 7/7 bound to lam0); golden 20/20.
  - [ ] **Review P1/P2 still open** (same doc): `coastFrac` is a junction
    count, not a time fraction (§19's "58%" = 7/12 junctions) — integrate
    the indicator; persist `tGrid` per entry; race `maxBisect` is per-arm
    not per-gap; race resume never reloads its queue; `ode113` completion
    check in both propagators; `catalog_schema` `1./mf`, version whitelist,
    `conj_pass ∈ {0,1}` on solved entries; REFERENCES blocks on the
    `costate_common` math functions; store a neighbouring-ε m_f delta per
    entry.
  - [ ] **Re-sweep the min-time catalogs' conjugate verdicts** with the
    corrected instrument: the 18,249 verdicts were produced before the
    t_f sample existed (golden cells unchanged; the change can only ADD
    refutations). `conj_catalog_pass` unchanged; ~hours at K=24.
  - [x] **Huber over the whole grid — DONE 2026-09-05 (FINDINGS §23,
    `run_huber_race_grid`, `huber_grid.mat`):** 4 of 7 records reach
    p = 0.001 with 0–1 fails in 1–2 min (ε: 8–10 fails, 6–11 min), same
    solution to <7e-6, conj PASS; (6,8)@1.2 walls at p≈0.033; (1,2) fails
    its first rung from the energy seed at both γ — λ/2 seed fixes the
    first rung (κ=1 asymmetry confirmed) but (1,2) still walls at p≈0.77.
    Huber = fast-or-walls, ε = slow-but-arrives. ε stays shipped.
  - [x] **Huber wall mechanism — DONE 2026-09-06 (FINDINGS §24):** not a
    fold (cond(J) flat), not the gate (MfMax loose acceptance tested, net
    worse); **grazing bifurcations** of the switch structure — near-tangent
    crossing |dQ/dt|=0.045 on (6,8), Q_max=0.9907 about to cross on (1,2);
    clean cells show neither. Instrument: `indirect/huber_switch_diag`.
    λ/2 seed confirmed exact (2 iterations at p=1). MfMax review §2b.
  - [ ] **Huber follow-ups (from §23/§24):** (a) the cure is a FAMILY change
    at the corner: Huber→ε handoff when `huber_switch_diag` predicts a
    graze (min|dQ/dt| falling / Q extremum → 1), or a Huber-ε hybrid (jump
    replaced by a ramp of width δ — one field function), or step OVER the
    bifurcation (larger p step on failure, not bisection); (b) hedged race
    in the builder; (c) Huber on the cells ε cannot enter (high-γ band via
    MfMax idea 1.5 IC-homotopy seeds, deep-thrust wall, GTO_tulip
    1.01–1.11× band); (d) λ/2 as the standard Huber seed.
  - [ ] **Independent verifier for the fixed-t_f extremals (Mike,
    2026-09-06):** (i) MATLAB `bvp5c` on the same PMP BVP (collocation, not
    shooting — different algorithm, an afternoon; shares our field); (ii)
    MfMax **v0** (fixed-t_f, time domain) with a new CR3BP `user.f90`
    (Dhfun/Ufun/B2fun/Pfun hand-written in ND from the equations — a true
    foreign witness: Fortran, rkf45, hybrd, FD Jacobian), run in `ipar=3`
    single-shoot mode from our λ₀; agreement ceiling ~1e-6 |dz| (single
    shooting over 20–30 d), same as `ss_bvp_accept`, but independent.
    Closes survey gap 3.5 for case B.
  - [ ] **PLQ-penalty experiment for the energy→fuel homotopy** (Mike,
    2026-08-31 — do not forget): mlabTools `optim/` ships Rockafellar's
    piecewise linear-quadratic penalty class (`plq_make`/`plq_penalty`:
    Huber, quantile, Vapnik ε-insensitive, subgradients free from the QP
    dual; already on every session's path via startup). Try PLQ shapes as
    the smoothing family between the energy (L2) and fuel (L1) objectives —
    e.g. Huber's quadratic-core/linear-tails as a *parametrized* L2→L1
    bridge with κ as the continuation knob instead of the discrete ε
    ladder, and `vapnik2`'s flat dead-zone + quadratic tails for
    soft-constraint experiments. Compare against the hand-scheduled ε
    ladder AND MfMax's arclength continuation (the standing MfMax-first
    rule still applies).

## GTO→tulip costate catalog (Stage B — BUILT 2026-08-29)

- [x] **Stage B catalog — SHIPPED 2026-08-29**: `GTO_tulip/catalog/`,
  16 sheets (4 GTO orientations {0,90,180,270}° × 4 Np {5,7,9,12}), 12×6
  phasing grid, rungs `[15 12 10 7 5]` N at Darin-standard 150 kg / 1710 s.
  **2,625 entries, 840/1,152 pairs (73%), conjugate 2,625/0/0.** Orientation
  coverage 240/207/156/237 of 288 — a real π-dip at 180° (apogee toward the
  Moon). Diagnostic A traced the *initial* failure mode there to a
  cold-start defect (not a mesh problem); the residual gap that remains
  after the warm-recipe fix is *inferred*, not traced, to be genuine
  cold-basin difficulty of that geometry (see README Coverage section).
  Two drivers shipped: `gto_entry` (single cell) + `run_gto_catalog`
  (swath/regen); dedicated zero-safe `costate_catalog_pick.m` (the generic
  picker's log-distance metric breaks at orientDeg=0). README:
  `GTO_tulip/catalog/README.md`. Stage A (shared-engine min-time solve,
  25 mN) DONE 2026-08-25.
- [ ] **Deep-rung investigation (OPEN, split off 2026-08-27)**: the 3–1 N
  legs of the standard ladder produced **zero entries** on this fleet by
  both the warm and cold-mop-up recipes — a measured closure wall, not an
  attempt-budget artifact. Root-cause before attempting any deep-rung
  extension; `densify_ladder`-style grinding at deep rungs needs a
  `wallSec`/attempt cap first (measured: 13 h churn, 8.3 MB log, no wall).
- [ ] Add dual capture (`lamDef`) to `gen_tulip_mintime` so harvest seeding
  works (found missing 2026-08-25) — now scoped to a future **flagship
  (25 mN) regime** extension of this catalog only; the shipped v1 (1–15 N,
  few-rev) does not need it. **Binding rule for that future extension**
  (spec, verbatim): many-rev entries must ship full ms junction states, not
  bare z8 — bare z8 pins t_f only to ~1e-4 ND / arrival to ~100 km at ~40
  revs (measured on this campaign's own flagship two-root adjudication).
- [x] Library moves `seed_from_z8` + `ms_tfmin` → `costate_common` — DONE
  2026-08-26, golden_cells 20/20, seed bit-identical.
- [x] **Two-root question ADJUDICATED 2026-08-26 (same day): ONE extremal.**
  Head-to-head flight: trajectories coincide to 23.5 km / 28 d; cross-fly is
  symmetric (each z8 lands 4.3 km only in its home environment, ~125 km in
  the other). The 94 s difference is z8 IDENTIFIABILITY, not a faster basin
  — bare z8 pins t_f only to ~1e-4 ND at ~40 revs. Certified direct
  t_f = 6.290694 STANDS. **Stage B consequence: GTO-regime catalog entries
  must ship full ms junction states, not bare z8.** Record:
  `GTO_tulip/indirect/min_time/README.md`.
- [x] Migrations #2–#5 (duals_to_costates, ms_bvp + conjugate test +
  golden cells, flown-control verifier, foc covector unification) — done
  2026-08-08, each with an equivalence gate.
- [~] Promote to a MATLAB package — STARTED 2026-08-09 as the top-level
  cross-folder `../oclib/+oc` (moves 1–2 done: `oc.duals_to_costates`,
  `oc.fly_control`, both with bit-identical acceptance on orbit AND
  booster flagship data). Next: move 3 = `oc.ms_bvp` +
  `oc.ms_conjugate_test` promotion + cart-pole PMP-BVP integration demo.
  See `../oclib/README.md` and `../OCP_UNIFYING_MATH.md` §5.
- [x] Schema versioning — DONE 2026-08-09 (schema v2: `catalog_schema`
  validator + named derive registry + environment pinning; v1
  compatibility proven bitwise).
- [~] Acceptance-gate harness: `costate_common/ss_bvp_accept` (2026-08-14)
  is the generic single-shooting form, used by min-energy; still open:
  route the min-time tfMin gate through the same interface.
- [ ] Batched-driver template + monitor template into the library (the
  halo/DPO drivers differ only in names; monitors must watch process
  liveness, not just log errors).

## Library consolidation queue (recorded 2026-08-31, Mike-approved)

Promotion candidates under the reuse-twice rule, ordered by readiness.
Discipline per `consolidation-discipline`: measure before extracting, each
move carries an equivalence gate; refusals get documented too.

1. [ ] **Catalog consumer-helper suite** — `costate_catalog_pick`,
   `costate_lib_describe`, `costate_catalog_extremes`,
   `costate_catalog_extremes_movies` each exist twice in-repo
   (`DRO_tulip/indirect` master + `GTO_tulip/catalog` copy) plus vendored
   deliverable copies. The GTO copies are strict supersets (orientDeg-aware
   labels, zero-safe CIRCULAR picker metric — the period-keyed copies'
   log-metric is latent-unsafe at key 0/wrap — plus the ephemeris guard in
   the movies). Consolidate into `costate_common` with the GTO versions as
   the base; DRO copies become delegates. **Do AFTER the 0.5 N extension
   ships** (avoid churning the deliverable surface mid-campaign).
2. [ ] **oclib move 3**: promote `oc.ms_bvp` + `oc.ms_conjugate_test` +
   cart-pole PMP-BVP integration demo (already on the oclib roadmap above).
3. [ ] **`ladder_endpoints` ↔ `densify_ladder`** — the family-aware
   endpoint rebuild now exists as a function (`DRO_tulip/indirect/
   ladder_endpoints.m`, extend path) AND as densify_ladder's older inline
   block. Merge (reroute densify through the function) when densification
   work restarts, with its own bitwise gate.
4. [ ] **`run_capped`** (parfeval hard per-call timeout, born inside
   `extend_thrust_ladder` 2026-08-31) — the only fence that bounds a
   crawling integration (measured twice, >4 h each). Promote to
   `costate_common` on its second consumer (densify at deep rungs, the
   0.1 N / 25 mN probes, any grinding solver call).
5. [ ] Batched-driver + monitor templates (the standing item above —
   listed here so the queue is complete).
6. [ ] **`gto_entry` → generic `costate_common/catalog_entry`** — promotion
   path documented in `GTO_tulip/catalog/README.md`; waits for a second
   campaign to want single-entry solves.

- [ ] **HALO catalog driver path gap** (found by Stage B task 3, 2026-08-26):
  `run_halo_catalog`'s addpath set predates the oclib moves — a fresh run
  would fail resolving `oc.local_residual` via the dro_residual reroute.
  Add `oclib/` to its setup (one line) + verify-campaign fast tier before
  the next HALO catalog run. Same check for the DPO/HALO_HALO drivers.

## Certification (see OPTIMALITY_CERTIFICATION.md before planning any of this)

- [ ] STRICT second-order sufficiency is still nowhere; the live path is
  the STM switching-time Hessian (Part B). The conjugate-point test
  (necessary condition) now runs on catalog cells and is exposed to the
  Earth campaigns via `foc_check` opts.msInfo — wire an Earth-campaign
  multiple-shooting solution to actually exercise it.
- [ ] foc gate coverage is uneven across the four direct campaigns
  (see the Part A matrix); close the standing gaps.
- [~] **Continuous-residual (G1) gate everywhere** — engine promoted
  (`oc.local_residual`), DRO rerouted, earth-MEE gated (2026-08-25, first
  ever: switch-interval concentration finding). CR3BP-GEO wrapper DONE
  2026-08-26 (shared `verify_common/mee_residual`, switch-interval finding
  generalizes). Ladder sweep DONE 2026-08-26 (14/14,
  `verify_common/doc/g1_sweep_results.md`). Open: ELFO wrapper only; optional:
  the lunar-phase-feedback check (replace t-row by quadrature on one deep rung).
- [x] **Run `ms_conjugate_test` across the shipped catalogs and store the
  verdict** — DONE 2026-08-23 (`costate_common/conj_catalog_pass`, K=24,
  ~33 min): 15,896 entries, **15,895 pass / 1 FAIL / 0 unverified**; verdicts
  in the catalog .mats (`conj_pass` grids + `conj_test` provenance,
  schema-validated). The FAIL — DPO τ=2 → Np=7, (2/3,2/3), 15 N — is refuted
  as a local minimum. Residue: deliverable-2 fine-sheet library (1,105
  entries, old format) not swept; shipped zips predate the verdicts.

## Legacy campaign items (tracked in their own TODOs)

- GTO_tulip: 1.01–1.11× near-min-time band still open; indirect certification.
- GTO_ELFO: indirect Route C not started.
- earth_elliptic_to_geo: 0.5 N min-time wall; our own MATLAB indirect.
- earth_elliptic_to_geo_CR3BP: Phase-2 indirect not started.
