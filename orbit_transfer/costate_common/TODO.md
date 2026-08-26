# costate_common — TODO

- [~] **Package promotion** — IN PROGRESS as the top-level cross-folder
  `../../oclib/+oc` (2026-08-09): `duals_to_costates` and the
  flown-control engine live there now (delegates here); `ms_bvp` +
  `ms_conjugate_test` are next (roadmap move 3, incl. a cart-pole PMP-BVP
  integration demo).
- [x] **Schema versioning** — DONE 2026-08-09: `catalog_schema.m`
  (versioned validator + named derive registry + environment pinning);
  v1 compatibility proven bitwise; first v2 catalogs = HALO_HALO
  deliverable 6.
- [~] **Acceptance-gate harness**: first concrete form DONE 2026-08-14 —
  `ss_bvp_accept` (single shooting on the same closures, K = 1) serves any
  cost with no pumpkyn twin; used by `ms_minenergy` (opts.accept). Still
  open: route the min-time tfMin gate through the same harness interface.
- [x] **ms_bvp fixed-tf variant** — DONE 2026-08-14 (`opts.fixedTf`; tests
  in `tests/`; golden cells bitwise unchanged). First consumer:
  `DRO_tulip/indirect/ms_minenergy` (min-energy pilot PASSED on flagship
  cells — see `DRO_tulip/FINDINGS.md`).
- [ ] **Fixed-tf conjugate-point test**: `ms_conjugate_test` quotients the
  free-time flow column; a fixed-tf variant (no flow column, λ-scaling
  invariance only partially present since L = s² fixes the scale) is
  needed before min-energy entries get a second-order verdict.
- [ ] **Conjugate test refinements**: junction-resolution disclaimer could
  be closed by dense STM sampling inside flagged intervals; an exact-zero
  det sample is currently counted as a crossing — refine to a bracketed
  root report.
- [ ] **Batched-driver + monitor templates**: halo/DPO shell drivers differ
  only in names; extract, and make monitors watch process liveness (a
  silent MATLAB death produced no log line and no alert, 2026-08-07).
- [ ] `survey_family_bounds` / older files still carry `%#ok` pragmas from
  before the no-pragma rule; strip on next touch.
- [ ] `golden_cells` engine cells converge in 1 iteration (wide Newton
  basin) — consider a rougher engineered seed so the iteration channel has
  more dynamic range.

- [x] Conjugate-point sweep at catalog scale — DONE 2026-08-23
  (`conj_catalog_pass`, 15,896 entries, 15,895/1/0; record in
  `DRO_tulip/process/COSTATE_LIBRARY_PIPELINE.md`). Standing rule: run it
  before packaging any new catalog/deliverable.
- [ ] Sweep the deliverable-2 fine-sheet library (1,105 entries, v2-library
  format) — needs a small adapter from the entry-array format to
  `conj_catalog_pass`'s sheet walk.
- [ ] Wire the conjugate sweep into the packaging path
  (`build_costate_catalog_family` or the deliverable checklist) so a catalog
  cannot ship without verdicts.

- [x] Library moves `seed_from_z8` + `ms_tfmin` -> here — DONE 2026-08-26
  (3 call sites rerouted; golden_cells 20/20; seed builder bit-identical;
  full test suite green).
