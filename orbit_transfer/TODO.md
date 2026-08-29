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
- [ ] **Sub-0.5 N rungs** deferred by Darin "until customers need it."
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
- [~] **Min-energy / fixed-t_f pipeline** — pilot DONE 2026-08-14 on three
  flagship DRO→tulip cells + a γ ladder (`DRO_tulip/run_minenergy_pilot`):
  fixed-tf `ms_bvp`, `ms_minenergy`, min-energy PMP field, generic
  single-shooting acceptance gate. Next: a t_f-multiplier axis in the
  catalog schema, a fixed-tf conjugate test, then energy→fuel homotopy on
  the same seeds (the min-fuel catalog route).

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
