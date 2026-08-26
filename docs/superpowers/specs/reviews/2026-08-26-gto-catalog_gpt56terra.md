## Findings

- **spec §3.3 / plan Task 6 — critical:** All orientation sheets share the same `tau_dep` (the GTO Kepler period) and `Np`; the existing nearest-`(tauDep,Np)` picker cannot select an orientation. **Fix:** add a GTO-specific orientation-aware picker that takes `orientDeg`, stores per-sheet `orientDeg`, and selects `(orientDeg,Np)` before grid lookup; do not overload `tau_dep`.

- **spec §2, §3.3 / plan Task 6 — major:** `dep_params` varies by sheet, but Task 6 describes it as a single catalog-level `{sma_km,ecc,orientDeg}` set. A packager default could reconstruct every sheet at `−25°` or another orientation. **Fix:** require and test sheet-local reconstruction recipes/parameters; package then replay one cell from at least two different orientations before any conjugate sweep.

- **spec §3.1 / plan Task 1 — major:** The flagship equality test is not independent: it compares the new implementation to `mintime_params`, using the same `orb2eci`/`fromPCI` convention. It only pins one point at epoch zero. **Fix:** add a geometric test that measures the Earth-relative perigee vector’s signed angle from the Earth→Moon vector in the rotating frame, plus an explicit assertion that epoch 0 aligns the inertial x-axis with Earth→Moon.

- **spec §3.1 / plan Task 1 — major:** The specification treats `orientDeg` as an epoch-zero rotating-frame quantity, but does not prohibit downstream code from interpreting `tau` as an epoch offset. That would rotate the Earth–Moon frame while retaining a fixed-orientation algebraic locus. **Fix:** audit/assert that `tau` is used only as a locus parameter for interpolation, never added to the CR3BP transfer start epoch; document this invariant in `get_family_orbit` metadata and tests.

- **spec §3.2 / plan Task 1 — critical:** A 25-point, highly eccentric GTO locus is too sparse for spline interpolation at the 12 departure fractions. This can silently create incorrect departure states between samples, especially around perigee. **Fix:** use a sufficiently dense mean-anomaly grid validated against direct algebraic states at every campaign grid point, or add an analytic GTO evaluator and have the engine call it instead of spline interpolation.

- **spec §3.2 / plan Tasks 1, 3 — major:** `mod(frac,1)` makes `sD_frac=1` identical to zero; the same may hold for arrival phase. If the grids are inclusive `linspace(0,1,n)`, one row/column is duplicate work. **Fix:** inspect the engine grid construction and assert unique wrapped fractions; use `[0:(n-1)]/n` or otherwise exclude the duplicate endpoint.

- **spec §2 / plan Task 3 — major:** The catalog omits the stated `−25°` flagship geometry and calls nearest-sheet selection adequate. Nearest `0°` is not an equivalent transfer problem and no orientation interpolation is established. **Fix:** add `−25°` as a shipped sheet, or replace one orientation only after explicitly changing the stated coverage and acceptance target.

- **spec §2 / plan Tasks 3–4 — major:** Four orientations at 90° spacing have no demonstrated continuation/interpolation validity. The pilot at `(0°,Np7)` cannot establish coverage for `90°/270°` or the longest target. **Fix:** run representative sentinel cells before the fleet for every orientation and for `Np=5,12`; add orientation resolution only where these reveal discontinuities or poor coverage.

- **context preflight fact / spec §3.2 / plan Task 2 — critical:** `preflight_screen` embeds cislunar admissibility assumptions, including a 500-km periselene floor, while this campaign departs a 350-km Earth-perigee GTO. It may reject valid cells, apply the wrong central-body safety test, or fail to protect Earth clearance. **Fix:** make screening body-aware for `dep_family='gto'`: explicitly enforce Earth-perigee/impact constraints and retain lunar-distance criteria only where physically intended; test representative perigee and apogee endpoints.

- **spec §3.2 / plan Task 2 — major:** Reusing the periodic-orbit engine is asserted rather than proven. Fixed-epoch GTO loci are closed only as a parameterized set, not as a rotating-frame flow orbit. **Fix:** add a compatibility test exercising interpolation, wrap, preflight, endpoint construction, and one warm-neighbor continuation across `sD≈0/1` before the pilot.

- **spec §4–5 / plan Task 3 — critical:** The pilot gate counts `any(Q.OK,3)`. A pair with only a 15-N solution passes, while the costly and scientifically important low-thrust rungs may all fail. **Fix:** gate on per-rung coverage, especially 1 N, ladder-complete-pair rate, and elapsed wall time; require thresholds before launching the remaining sheets.

- **spec §4 / plan Task 4 — major:** “Resume for free” is not established. Attempt counters are saved before solves; an interruption between solve completion and atomic result persistence can leave successful work unrecorded or mark cells exhausted. **Fix:** perform a deliberate kill-and-resume test during the pilot and verify no lost/duplicated entries, monotone attempt counters, and data-file—not log—recovery.

- **plan Task 4 — major:** The full campaign has no immutable run manifest. Pumpkyn/code/path changes during days of resume can mix solver behavior across sheets and invalidate comparisons. **Fix:** save MATLAB version, pumpkyn revision/path, git revision/diff status, engine options, and hashes of relevant functions into every sheet’s metadata before first solve.

- **plan Task 5 — minor:** The audit hard-codes three coordinates that may not have solutions, particularly the low-rung corner. **Fix:** select three available entries programmatically from required distinct sheet/rung classes, then record their actual indices and replay results.

## Top 3

1. Implement an orientation-aware sheet key/picker; `tau_dep` cannot encode orientation.  
2. Replace or validate the 25-sample spline locus before any solver run.  
3. Make GTO-specific preflight and pilot gates prove low-rung coverage, not merely one successful rung.
