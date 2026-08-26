**Finding 1: Orientation-Angle Unit Test False Positive**
* **Cite:** Spec 3.1 / Plan Task 1, Step 1
* **Issue:** The unit test verifies flagship equality and geometry exclusively at apsidal points (`mAnom` = 0 and π), where mean and true anomaly are numerically identical. If `orb2eci` expects true anomaly but the implementation directly feeds it the mean anomaly array, this test will falsely pass. This would silently generate a locus parameterized by true anomaly, entirely corrupting the time-fraction physics required by the engines. 
* **Severity:** Critical
* **Concrete Minimal Fix:** Add a non-apsidal check (e.g., at `mAnom = pi/2`) in `test_gto_family.m` against a known correct state to strictly pin the required Kepler mean-to-true anomaly solve.

**Finding 2: Anomaly Axis Interpolation Error**
* **Cite:** Spec 3.1 / Plan Task 1, Step 3
* **Issue:** The departure locus is generated using `M=25` uniformly spaced points in mean anomaly. For a highly eccentric GTO (e≈0.73), true anomaly changes extremely rapidly near perigee. Using a standard spline interpolation over just 25 time-spaced points will result in severe interpolation overshoot and velocity errors precisely where the physics is steepest, corrupting the transfer boundary conditions.
* **Severity:** Critical
* **Concrete Minimal Fix:** Increase the sample count `M` to ≥1000 in `get_family_orbit`, or explicitly build the locus using uniform spacing in true/eccentric anomaly (and map to a non-uniform `tau` array) so the spline adequately resolves perigee.

**Finding 3: Sheet Keying Collision**
* **Cite:** Spec 3.3 / Plan Task 6, Step 1
* **Issue:** The spec mandates placing the GTO Kepler period into the `tauDRO` (legacy `tau_dep`) slot. Because the Kepler period is identical across all four orientation sheets, the downstream picker (`costate_catalog_pick`, which selects by nearest `tau_dep` and `Np`) will be completely unable to distinguish the orientations, permanently rendering 75% of the computed catalog inaccessible to consumers.
* **Severity:** Critical
* **Concrete Minimal Fix:** Map `orientDeg` directly into the `tauDRO` / `tau_dep` slot in `build_gto_catalog.m` so the existing picker machinery can uniquely resolve the sheets, and document this override in the README.

**Finding 4: Grid Design Misalignment**
* **Cite:** Spec 2 / Plan Task 3, Step 1
* **Issue:** The planned 12-point anomaly grid is uniformly spaced in time fraction (`sD_frac`). For an eccentric GTO, uniform time spacing concentrates almost all geometric evaluation points near apogee. This directly contradicts the spec's intent to "double departure resolution where the physics is steepest" (perigee is left vastly under-resolved).
* **Severity:** Major
* **Concrete Minimal Fix:** Override the default `sD0` array in `run_gto_catalog.m` with a non-uniform grid that explicitly clusters points near perigee (e.g., using a cosine spacing or by mapping a uniform true-anomaly grid back to time fractions).

**Finding 5: Engine-Reuse Admissibility Screen Rejection**
* **Cite:** Plan Task 2, Step 4 / Context Line 12
* **Issue:** The `thrust_ladder_library` utilizes a `preflight_screen` originally tuned for cislunar space (≥500 km periselene floor). A standard GTO has a 350 km Earth perigee. If this screen generically evaluates nearest-body altitude, it will immediately reject valid perigee-departing transfers and fail the initial cold-start seed.
* **Severity:** Major
* **Concrete Minimal Fix:** Add logic inside `thrust_ladder_library`'s preflight screen to bypass or lower the Earth-altitude floor (e.g., to 100 km) specifically when `depFamily == 'gto'`.

***

**Top-3 List:**
1. **Sheet Keying Collision (Finding 3):** Ships a broken deliverable where downstream tools can only ever access the 0° orientation sheet.
2. **Anomaly Axis Interpolation Error (Finding 2):** Unresolved M=25 splines will permanently corrupt the boundary conditions of the solver at perigee.
3. **Orientation-Angle Unit Test False Positive (Finding 1):** Testing only apsidal points masks missing anomaly conversions, risking days of compute on a physically invalid parameterization.
