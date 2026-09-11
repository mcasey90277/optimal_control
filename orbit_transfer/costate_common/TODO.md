# costate_common — TODO

- [ ] **DEFERRED, with the measurement: a Hermite scheme for `periodic_pp`**
  (asked 2026-09-11 — is a choice of interpolant TYPE worth an option?).
  Measured first, on the tau = 1 DRO (105 samples) and the 7-petal tulip
  (1328): the periodic cubic's true error against propagation is **2.9 m
  median / 4.9 m max (DRO)** and **0.23 m / 1.07 m (tulip)** at interval
  midpoints, and `d/dt` of the interpolated position disagrees with the
  interpolated velocity by **1.5 mm/s** (DRO) / **0.23 mm/s** (tulip)
  against a ~500 m/s velocity scale. The certification gates are 100 km and
  10 m/s, and the worst audited miss is 0.29 km: interpolation is ~4 orders
  from binding, so NOT built.
  - The variant worth building when it IS needed is **Hermite**, not a
    higher-order spline for its own sake: the table already carries the
    derivative of its position rows (the velocity rows), and the
    acceleration is available from `cr3bp_field`. A spline discards both.
    Hermite would enforce them, cutting the error and making the
    interpolant self-consistent by construction. `opts.scheme` is already
    the hook, so deferring costs one switch case.
  - Ruled OUT on measurement, not taste: **trigonometric/FFT** (the natural
    periodic basis) needs uniform samples, and these tables come off the
    propagator with spacing ratios of **16384:1** (DRO) and 6221:1 (tulip);
    resampling would bake in the cubic's own error. **pchip / makima** are
    C1 only — worse for the seam derivative this function exists for.
  - TRIGGERS to revisit: an endpoint tolerance within ~1 km of the
    interpolation error; a continuation needing physically accurate
    derivatives better than ~1e-5 relative (today's 9e-13 is
    self-consistency against the same interpolant, not accuracy); or a much
    sparser orbit table (error scales as h^4).

- [ ] **Migrate the `interp1(..., 'spline')` endpoint sites onto `phase_state`**
  (2026-09-11, FINDINGS 44). About a dozen: `second_order_pass`,
  `conj_catalog_pass`, `gates_catalog_pass`, `DRO_tulip/direct/sweep_phasing_direct`,
  and the per-campaign `costate_catalog_extremes` / `*_example` helpers. One at a
  time, each with its own before/after check. NOT the shipped deliverable copies:
  those must stay self-contained, and the measured cost there is metres, only
  beside the seam.

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
