# costate_common — TODO

## Cleanup plan (2026-09-16), easiest to hardest

Measured state and the reasons for each step are in `README.md` → *State of
the folder*. Steps 1–9 change no numerical result and need no decision; steps
10–13 depend on step 0. Each step names the gate that says it is done.

- [x] **0. Two decisions (Mike) — DECIDED 2026-09-16.** (a) The admission
  rule becomes *generic by construction, with a test*: a file belongs here
  when nothing in it is specific to one campaign's orbits or bookkeeping and
  a test pins its contract, whether or not a second campaign calls it yet.
  The reason: DRO_tulip is the reference campaign new ones copy, so under
  *used by two campaigns* code born there could never be admitted. Clear
  stays: `ms_bvp`, `arclength_ms`, `conj_resolve`, `scalar_verdict`; clear
  goes: `rib_targets` (step 6). (b) Job control moves to a sibling folder,
  `orbit_transfer/campaign_common/` (step 11); the execution fences
  (`run_capped`, `capped_pool`, `current_pool`) stay here.

**A. Documentation and style — no behaviour change**

- [x] **1. README reflects the folder** — DONE 2026-09-16 (all 61 files in
  layers with measured consumers; tests classified; conventions with their
  current violations).
- [x] **2. Regenerate `../doc/library_catalog.md`** — DONE 2026-09-16
  (84 functions).
- [x] **3. Fix the junction contract text** — DONE 2026-09-16. On reading,
  the `info.Y` OUTPUT docs of `ms_bvp` and `ms_tfmin` were already right
  ([14 x K] starts); the stale text was the SEED docs, which did not say
  column K+1 is never read, and `DRO_tulip/indirect/probe_deep_rungs`,
  which described its stored `it.Y` as K+1. Seed-builder outputs
  (`seed_from_z8`, `seed_from_entry`, `harvest_ms_seed`,
  `flight_to_junctions`) and the min-fuel catalog's `.Yj` genuinely are
  K+1 (`build_minfuel_catalog` appends the endpoint) and were left alone.
  Gate met: Code Analyzer messages identical on all three files.
- [x] **4. House-style pass** — DONE 2026-09-16. 21 `%#ok` lines removed
  (suppression dropped, no code changed — no preallocation, so behaviour is
  untouched); `j` loop variables renamed; the 4 `% INPUTS:` headers now
  `%% Purpose` / `%% References` / `%% Inputs` / `%% Outputs` /
  `%% Revision History` (content kept; `gates_catalog_pass` now also
  documents `.lamVTol`, which the code already read). Gate met: Code Analyzer
  diff is exactly the 22 messages the 21 pragmas had suppressed (one
  `cr3bp_minfuel_prop` line grows two arrays), nothing else; touched tests
  pass (see README → Tests). The `AGROW` warnings are now visible, which is
  the point; preallocate them on next touch.

**B. Small removals and moves — one consumer each**

- [x] **5. Retire unused code** — DONE 2026-09-16. `conjugate_pole_predict`
  and its test DELETED (no caller; last present at commit f84e82b if the
  pole fit is ever wanted as a diagnostic). `cr3bp_field` became a local
  function of `tests/test_arclength_arrival`, its only consumer; the
  deferred Hermite item notes where to restore it from.
  `DRO_tulip/indirect/arclength_arrival`'s header still claimed its
  derivative came from `cr3bp_field`; corrected to the interpolant.
  Gate met: no code reference remains; `test_arclength_arrival` passes.

- [x] **6. Move `rib_targets` to `DRO_tulip/indirect`** — DONE 2026-09-16,
  with `test_phase_lists`. Gate met: the test passes from its new home.

- [x] **7. Move the campaign-only tests** — DONE 2026-09-16: 15 tests to
  `DRO_tulip/indirect/tests` (each bootstrap now climbs to `orbit_transfer`
  and names `costate_common` explicitly). `test_gto_family` was
  misclassified and STAYS: it tests this folder's `get_family_orbit`
  through a GTO_tulip fixture. Gate met: all 15 plus `test_phase_lists` and
  `test_arclength_arrival` pass from their homes, each run with
  `costate_common`, `DRO_tulip` and `DRO_tulip/indirect` removed from the
  path first, so every test had to find its own dependencies.

**C. Cut the library's dependency on DRO_tulip**

- [x] **8. Make `golden_cells` self-contained** — DONE 2026-09-16. The
  harvest cell's collocation inputs (`X(1:7,:)`, `lamDef`, `Um`, `tNodes`,
  endpoints, thruster, mass ratio) are copied into
  `golden_cells_data.mat` as `.harvestCell` (the file grows 1 KB -> 116 KB);
  the `DRO_tulip/indirect` addpath is gone (nothing it called lived there
  any more). Gate met: 20/20, and every cell's z, iters, normR, conjugate
  verdict, vote margin and lambda_t BITWISE equal to the pre-change run,
  executed with no DRO_tulip/HALO/DPO/GTO directory on the path and the
  gitignored `dsweep_12x12_cells.mat` renamed out of reach.
- [x] **9. Bring `ladder_endpoints` here** — DONE 2026-09-16, with
  `test_ladder_endpoints` (now in `tests/`, reading DRO_tulip and HALO_tulip
  results only as fixture data). `second_order_pass` no longer adds
  `DRO_tulip/indirect` to the path, and neither do its three tests, which
  had added it only to reach `ladder_endpoints` (and would have hidden a
  broken move). All eight DRO_tulip callers already put this folder on
  their path. Gate met: `test_ladder_endpoints`, `test_second_order_pass`,
  `test_second_order_sidecar_identity` and `test_second_order_parallel`
  pass with every DRO/HALO/DPO/GTO directory removed from the path first;
  `grep DRO_tulip costate_common/*.m` finds comments, the demo's and
  `golden_cells`' data-file references, and no addpath.
- [x] **10. Tests for the untested core** — DONE 2026-09-16, seven tests,
  each on a fixture with a known answer: `test_harvest_ms_seed` (real duals
  from `golden_cells_data`, bitwise against the pchip/extrap construction,
  sign flip, uniform mesh), `test_run_capped` (in-time, worker error,
  timeout within the cap, pool survives), `test_flown_control_error` (true
  flight of a quadratic throttle: zero error, injected misses read back,
  linear reconstruction caught; also covers `ctrl_quad` and
  `cr3bp_thrust_rhs`), `test_true_min_altitude` (a flyby whose periselene
  falls between nodes), `test_preflight_screen`, `test_newton_fixed_q`
  (analytic root, guards, cap), `test_survey_family_bounds` (real DROs,
  one admissible, one rejected, one unbuildable). Gate: all pass with no
  campaign directory on the path; one planted bug per function is caught
  by a failing ASSERTION (not a throw), files restored md5-identical. The
  one undetected mutant (dropping `newton_fixed_q`'s non-finite-RESIDUAL
  guard) is equivalent: a NaN residual always yields a NaN step, which the
  step guard catches on the same call.

**D. Structural moves — each needs a campaign-level reproduction, not only unit tests**

- [x] **11. Split out campaign orchestration** — DONE 2026-09-16. The nine
  files plus `run_campaign_workers.sh` and `campaign_supervisor.sh` and
  their two tests now live in `../campaign_common/` (README there). The
  fences stayed. Consumers repointed: `run_costate_library` (third code
  root threaded into the generated rib and finalize jobs; the supervisor is
  launched from the new folder), `build_ribs`, `build_arrival_sheet`,
  `fill_holes_direct`, and four archived `batch/torus/v*_job.m`.
  Gate met, in two parts: (a) `test_work_queue` and
  `test_campaign_processes` (real worker processes, kills, supervisor,
  finalizer) pass with `costate_common` and every campaign folder stripped
  from the path; (b) the torus3b 3 x 3 acceptance campaign was re-run
  end to end from the same spec (`results/torus3c_step11`, job script
  adding ONLY costate_common, so the driver had to find campaign_common
  itself): same two rounds, same sheet counts, the same new anchor
  (17.705 d at sA 0.7421), same holes census, fixed point after round 2,
  and its final catalog is BITWISE identical to the 2026-09-15 one in
  every field but the provenance note carrying the campaign tag
  (`test3b` -> `test3c`), which was the change.

- [x] **12. Bring the shared thrust ladder here** — DONE 2026-09-16.
  `thrust_ladder_library`, `casadi_mintime_dro`, `certify_dro_mintime` and
  `dro_residual` (with `test_dro_residual`) moved from `DRO_tulip` into this
  folder; their fixed-depth self-`addpath` guards are gone (the siblings are
  here now), and `dro_residual` keeps one for `oclib`. HALO, DPO, HALO_HALO,
  `probe_l1_l2_halo` and GTO's `setup_paths` no longer put ANY DRO_tulip
  folder on the path. Eight DRO scripts that had reached the engine through
  their own `direct/lib` and `direct/certify` now add `costate_common`.
  Placement note: `dro_residual`'s MEE sibling `mee_residual` lives in
  `../verify_common`, which argues for putting both there — refused for now
  because `thrust_ladder_library` → `certify_dro_mintime` → `dro_residual`
  would then make the two library folders mutually dependent. Unifying the
  two residuals is a separate item, below.
  Gate met: the per-campaign re-solve, run BEFORE and AFTER the move (the
  before-run proves the gate itself). One stored cell per campaign, its top
  three rungs re-solved from the sheet's own meta (each rung seeds from the
  one above, so a truncated ladder reproduces the stored run's first rungs),
  compared against the shipped sheet — HALO, DPO, HALO_HALO and GTO, plus
  `golden_cells` 20/20, `test_dro_residual` 3/3 and the CasADi engine's
  bitwise regression `test_minenergy_objective`.

**E. Rule-gated**

- [ ] **13. Promote the problem-agnostic engines to `../../oclib/+oc`:**
  `ms_bvp` and `ms_conjugate_test` first (named "next" below since
  2026-08-09), then `arclength_ms`, `newton_fixed_q`, `conj_resolve`,
  `lift_space_dim`. Blocked by `oclib`'s admission rule, not by effort: each
  needs a second TOP-LEVEL consumer (the planned cart-pole PMP-BVP demo, or
  booster_landing) and an equivalence gate, with a delegate left here as for
  `duals_to_costates`. Days, after a consumer exists.
  - [x] **`ms_bvp` DONE 2026-09-16** — the consumer arrived: the cart-pole
    PMP-BVP demo (`optimal_control_examples/ex3_cart_pole_pmp`) solves its
    Pontryagin BVP through the same engine with no orbit, no CR3BP quantity
    and no pumpkyn call, which is the second TOP-LEVEL consumer the rule
    asks for. Implementation now at `../../oclib/+oc/ms_bvp.m`; this folder
    keeps the delegate, so no caller changed. Gates, all run on the move:
    `golden_cells` 20/20 **bit-identical** to the pre-move capture (z, ‖R‖
    and iteration counts equal to the last bit, `max|Δz| = 0` on all four
    cells); the four-campaign ladder re-solve unchanged to the printed digit
    (HALO 4.551e-14, DPO 1.130e-14, HALO_HALO 3.478e-14, GTO 2.720e-14
    max |Δt_f|, OK flags matching); `test_ms_bvp_extra`,
    `test_ms_bvp_fixedtf`, `test_ms_tfmin_hom`, `test_arclength_ms`,
    `test_folder_rules` green; `test_cartpole_pmp` 9/9 through `@oc.ms_bvp`.
    No new folder-rules exemption was needed — the two `ms_bvp` tests call
    the delegate, so rule 6 is satisfied by the forwarding itself.
  - [x] **`ms_conjugate_test` DONE 2026-09-17** — the second consumer
    arrived: the cart-pole study scripts (`optimal_control_examples`,
    cart-pole sub-project A Task 4) call the Jacobi test directly on a
    fixed-t_f min-energy PMP extremal, which is the second TOP-LEVEL
    consumer the rule asks for. Implementation now at
    `../../oclib/+oc/ms_conjugate_test.m`; this folder keeps the delegate,
    so no caller changed (only two comment lines and one Revision History
    line differ from the pre-move file — no executable line moved). Gates,
    all run on the move: `golden_cells` 20/20 **bit-identical** to the
    pre-move capture (dro/halo/dpo/harv z, iters, normR all unchanged);
    the four-campaign ladder re-solve unchanged to the printed digit
    (HALO 4.551e-14, DPO 1.130e-14, HALO_HALO 3.478e-14, GTO 2.720e-14
    max |Δt_f|, dZ8 bit-identical too, OK flags matching); `test_conj_fixedtf`,
    `test_conj_coverage`, `test_conj_spectrum`, `test_conj_resolve`,
    `test_folder_rules` green, PASS/FAIL output byte-identical pre/post move.
    No new folder-rules exemption was needed — the five conjugate tests call
    the delegate, so rule 6 is satisfied by the forwarding itself, same as
    `ms_bvp`.
  - [ ] **The rest stay open, for the same reason they always were.**
    `arclength_ms`, `newton_fixed_q`, `conj_resolve` and `lift_space_dim`:
    no second top-level consumer yet.

---

- [ ] **`ms_bvp`'s `residual()` bare `catch` on `prob.prop` (found while
  building the cart-pole PMP-BVP demo, 2026-09-16 — not fixed here, `ms_bvp`
  itself is out of scope for that branch).** The catch has no identifier
  check and no logging: it treats a genuine coding bug inside a CONSUMER's
  propagator (undefined variable, missing struct field, dimension mismatch)
  exactly like a rejected iterate, returning the same `rejectR`/`eye(n)`
  pair either way. The shooting loop then simply fails to converge, with no
  trace of which case it was. This now reaches every consumer of the
  package, not only `orbit_transfer`: `optimal_control_examples/
  ex3_cart_pole_pmp/cartpole_pmp_prop.m`'s header documents the identical
  hazard from the cart-pole side and narrows its OWN try/catch (relabelling
  only `MATLAB:ode*` identifiers, rethrowing everything else unchanged) so
  it does not compound the problem, but that is a consumer-side mitigation,
  not a fix. Candidate fix: rethrow anything without the propagator's own
  documented collapse identifier, or at minimum log `err.identifier` /
  `err.message` before rejecting the iterate.

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
    acceleration is available from the ballistic CR3BP field (a local
    function of `tests/test_arclength_arrival` since 2026-09-16; promote it
    back into the library when this is built). A spline discards both.
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

- [x] **Migrate the `interp1(..., 'spline')` endpoint sites onto `phase_state`**
  — DONE 2026-09-11 (FINDINGS 44), 18 files in three tiers: the live
  instruments (`second_order_pass`, `conj_catalog_pass`, `gates_catalog_pass`,
  `audit_phase_catalog`), the ladder engines and phase sweeps (11 files incl.
  `thrust_ladder_library`, which the halo and DPO campaigns call unmodified),
  and the GTO campaign tools. Each site kept its previous ORIENTATION; the
  interpolant is now built once beside the orbit tables instead of per entry.
  Verified: gates/conjugate/audit recomputed equal on measured cells,
  golden_cells 20/20, Code Analyzer message-for-message against the committed
  versions.
- [x] **The ND propulsion conversion is one function** — DONE 2026-09-11
  (`nd_propulsion`): 28 files migrated, every value bitwise unchanged
  (the library holds the same expressions character for character).
  Verified: golden_cells 20/20, gates/conjugate/audit recomputed equal on
  measured cells, Code Analyzer message-for-message across all 28.
- [x] **The seed-from-flight cut is one function** — DONE 2026-09-11
  (`flight_to_junctions`): six engines plus `seed_from_z8`. The mass law
  (all-burn identity) is an option there, so the DERIVED-not-rescaled rule
  has one home. Gate: bitwise equal to the engines' inline form on a real
  flight; `golden_cells` covers `seed_from_z8`'s 3e-13 query-form shift.
- [ ] **`certify_dro_mintime` has no test of its own** (2026-09-16): it
  arrived with the ladder engine and is covered only through campaign
  re-solves, so it carries a `gap` exemption in `tests/test_folder_rules`.
  It composes four gates whose thresholds decide what a campaign certifies;
  a fixture test (a stored solution, one gate perturbed at a time) would
  pin them. Drop the exemption when it lands.

- [ ] **Reflow the four ladder-engine headers into the aligned columns**
  (2026-09-16): converted to `%% Purpose` on arrival with their documented
  content kept as prose; the Inputs/Outputs blocks are not in the house
  column format yet. Cosmetic, mechanical.

- [ ] **Unify `dro_residual` and `verify_common/mee_residual`** (2026-09-16,
  from step 12): both are thin layouts over `oc.local_residual` — the CR3BP
  Cartesian [r;v;m] split and the MEE one. They now live in different
  library folders. Either give `oc.local_residual` a layout descriptor and
  keep one caller-side splitter, or move both beside each other once the
  ladder's dependency direction allows it.

- [ ] **The shooting residual is PROPAGATION-MODE dependent, and every
  campaign quotes the tighter mode.** Measured 2026-09-11 on the 70 mN
  anchor: the SAME converged point gives `|R| = 7.84e-12` when the residual
  is asked for with its Jacobian (the integrator carries 210 states and
  takes finer steps -- the mode `ms_bvp` iterates in, and the mode every
  `tolR` gate is therefore evaluated in) and `1.41e-09` under plain 14-state
  propagation. A factor of 180, far above the 3e-11 the solver reports
  meeting. Nothing is wrong -- but "converged to 3e-11" means "in the
  Jacobian mode", and a re-evaluation of a stored solution will not
  reproduce it. Worth (a) saying so in the methodology doc, (b) deciding
  whether the catalogs should ALSO carry the state-only number, which is
  what a recipient re-evaluating an entry will see, and (c) checking whether
  the gap grows at deeper thrust, where arcs are longer.

- [ ] **Pointwise PMP checks for the FIXED-t_f objectives** (2026-09-11,
  Mike's question "can these be library generic functions?"). `pmp_pointwise_checks`
  already injects the vector field (`opts.rhs`), so it is generic across
  PROBLEMS -- it runs on every certified entry in DRO/HALO/DPO/HALO_HALO/GTO.
  It is NOT generic across OBJECTIVES: it hardcodes the min-time running cost
  (H = 1 + lambda.f), the free-time consequence H == 0, and the all-burn
  control recovery. So the min-energy and min-fuel catalog entries have NO
  pointwise check at all -- they carry the conjugate verdict and the
  hypothesis gates only. The fix is to inject the running cost and the
  optimal-control law the way the field already is (H constant rather than
  zero at fixed t_f; the smoothed argmin instead of bang). New capability,
  not a refactor; gate it on the min-fuel catalog's 18 entries.

- [ ] **Propulsion conversion: the sites NOT migrated, and why.**
  `cr3bp_common/cr3bp_lt_params` and the GTO `direct/` and `indirect/`
  campaign scripts still hold their own copies, because those modules do
  NOT have `costate_common` on their path (checked: only
  `GTO_tulip/catalog/setup_paths` adds it). Routing them here means a
  shared library depending on another shared library — an architecture
  call for Mike, not a refactor. Note those campaigns are already
  deduplicated LOCALLY through `cr3bp_lt_params`, so the duplication there
  is one file, not many. Also out: `GTO_tulip/attic/*` (archive) and the
  INVERSE conversions in `catalog_schema` / `build_minfuel_catalog`
  (c_nd → Isp), which are a different rule and would be their own addition.
- [ ] **Deliberately NOT migrated** (record, so it is not "found" again): the
  recipient-facing helpers that ship beside a catalog — `costate_lib_example`,
  `costate_catalog_example`, `costate_catalog_extremes(_movies)`,
  `costate_lib_extremes`, `costate_lib_v2_example`, the GTO
  `costate_catalog_gto_example` / `costate_catalog_extremes`, and the recipe
  STRINGS inside `build_costate_lib*` / `build_costate_catalog`. A shipped
  catalog's helpers must have no dependency on this library (GTO catalog
  README: "inlined, per deliverable-picker convention"), and the measured cost
  of the ordinary spline there is millimetres, only beside the seam. Revisit
  only if a deliverable is ever allowed to depend on `costate_common`.

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
- [x] `survey_family_bounds` / older files still carry `%#ok` pragmas —
  DONE 2026-09-16 with cleanup step 4.
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

## From the 2026-09-11 math review (FINDINGS 47) -- theory items, not code

- [ ] Audit document: prove the SUBARC corank-one condition from the
  whole-arc lift-space rank via analyticity of the strict all-burn CR3BP
  extremal (away from collision, m = 0, |lam_v| = 0): a stationarity
  constraint vanishing on an open subinterval extends along the arc.
- [ ] Audit document: write the reduction of the free-mass, free-time
  second variation / critical cone to the six-state NONAUTONOMOUS problem
  the quotient determinant tests; state that the flow column is the
  arrival-time variation, not an endpoint-preserving gauge.
- [ ] Conjugate instrument: close the interval before the first full-rank
  junction with a short-time sign expansion of det[J P, f] (currently
  reported as uncovered; the dense scan samples inside the first segment).
- [ ] Between-sample bounds for S1 (|lam_v| > 0), S2 (Q > 0) and the
  primary clearances: derivative bounds or interval integration; today they
  are sample minima and say so.
- [ ] `report_optimality`'s claimed-verdict sentence still reads "same
  endpoints"; align with the study script's endpoint statement.

## From the 2026-09-11 SECOND math review (Astra #2; FINDINGS 48)

- [ ] **RE-SWEEP the 70 mN library's second-order sheet** with the resolved
  `conj_spectrum` (located minima, shifted grids, UNRESOLVED class, endpoint
  clusters refined). The shipped sheet's 44 "near-miss" cells are OLD
  plateau verdicts; `conj_unresolved` is NaN for them. Clear the sidecar
  (or add a version key) so `second_order_pass` does not resume over stale
  records. Blocks the ship decision.
- [ ] Sharpened theory items (supersede the four above where they overlap):
  (a) subarc normality must carry the terminal condition lam_m(t_f) = 0 to
  a MOVED subarc endpoint, with H6 supplying the nonvanishing
  normalisation; (b) the six-state reduction must cover the COMPETITOR
  class, including reduced-throttle directions, to keep the claim strong
  rather than restricted/weak; (c) existence of an exact extremal near the
  numerical one (validated shooting or a Newton-Kantorovich/Krawczyk
  argument) -- a residual alone does not supply it; (d) between-sample
  bounds must include the conjugate MATRIX (rank/determinant), not only
  the scalar margins.
- [ ] `conj_spectrum` zeroFloor (1e-7 x median) and clearFactor (100) are
  POLICY values for the STM's accuracy; measure the STM error (two
  integration settings, as lift_margin does) and set the floor from it.
- [ ] The start transient (cluster touching sample 1) is reported as
  uncovered; the short-time sign/positivity expansion is still the only
  thing that would cover it.

## From the 2026-09-11 THIRD math review (Astra #3; FINDINGS 49) -- still open

- [ ] `conj_spectrum` floor from a MEASURED scaled-matrix error: safety
  factor x (matrix error + SVD error + slope x location uncertainty) / median,
  with the matrix error measured by re-propagating from t = 0 at a tighter
  setting or another integrator and comparing the PROJECTED matrices at the
  candidate times and t_f (both current refinements inherit the same stored
  prefix). Until then zeroFloor 1e-7 / clearFactor 100 are policy values.
- [ ] `ms_conjugate_test` resolvedTol from the same measurement: trust a
  sign only when sigma_min exceeds matrix + LU + SVD error with margin.
- [ ] STM discriminator: check the generator block A(4:6,11:13) =
  -(T/(m rho))(I - alpha alpha') along the trajectory (tangential eigenvalues
  -T/(m rho)); evaluate the actual variational RHS with an identity STM.
- [ ] S3 lift residual: add blockwise / componentwise residuals (alignment
  rows vs the terminal-mass row) and the angle bound |C lam|/(sigma_6 |lam|)
  beside the global backward error; recompute lift_margin under any change
  of row weights.
