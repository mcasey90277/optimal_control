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
- [ ] **Acceptance-gate harness**: the tfMin acceptance gate is the one
  pipeline check not yet generic — generalize as a harness taking a
  per-family independent solver.
- [ ] **ms_bvp fixed-tf variant** (the header forbids faking it with tight
  guards); needed the day a fixed-time catalog appears.
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
