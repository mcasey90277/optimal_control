# costate_catalog_gto_deliverable7 — GTO → Tulip Minimum-Time Costate Catalog

Every entry in this catalog is a **converged minimum-time PMP solution**, not a
hint: `pumpkyn.cr3bp.tfMin` accepts each one **unchanged** (|Δz| < 1e-6 gate;
the worked example lands at |Δz| ≈ 8e-12). Give the picker five coordinates —
GTO departure **orientation**, tulip petal count, departure phase, arrival
phase, thrust — and you get the flight time and the full costate/final-time
vector `z8 = [λr(3); λv(3); λm; tf]` ready to fly or to seed `tfMin` at your
own endpoints.

This is the fifth catalog in the series (after DRO/HALO/DPO/L1↔L2-halo) and
the first with a **non-periodic departure orbit**. There is no departure
period to key sheets on, so the sheet axis is an orientation angle — read the
Orientation convention section before anything else.

## Contents

| File | Role |
|---|---|
| `costate_catalog_gto_tulip.mat` | The catalog: 16 sheets, 2,625 entries (schema v2, compact ND format, conjugate verdicts inside) |
| `costate_catalog_pick.m` | Five-coordinate lookup with the honesty contract (warns + reports what you actually got). **orientDeg-keyed — not interchangeable with the period-keyed catalogs' picker** |
| `costate_catalog_gto_example.m` | Worked example: pick → derive units → rebuild orbits → fly → tfMin → plots |
| `costate_lib_describe.m` | Prints the catalog's coverage, grids, orientation axis (orientDeg-aware labels) |
| `costate_catalog_extremes.m` | Min/max ΔV or flight time over any filtered slice, with plot |
| `costate_catalog_extremes_movies.m` | Renders the two extremes as side-by-side pumpkyn-style movies |
| `get_family_orbit.m` | Family name + params → orbit (the `'gto'` pseudo-family + tulip recipes) |
| `gto_extremes_movie.mp4/.gif` | Sample: cheapest (8.16 km/s) vs costliest (10.31 km/s) phasing on the orient=270°, Np=5 sheet at 5 N, shared clock |

All helper self-demos run on whichever `costate_catalog_*.mat` sits in the
current folder — this one included (type e.g. `costate_lib_describe` bare).

## Dependencies

pumpkyn + pumpkynPie on the MATLAB path (run their `startup.m`). The `.mat`
itself is dependency-free and self-describing.

## Quick start

```matlab
L = load('costate_catalog_gto_tulip.mat');
cat_ = L.costate_catalog_gto_tulip;
[tf_nd, z8, info] = costate_catalog_pick(cat_, 90, 7, 0.0, 0.0, 10.0);
% 90  -> departure orientation in DEGREES (sheets at 0/90/180/270)
% 7   -> tulip petal count
% then run costate_catalog_gto_example.m end-to-end
```

## Orientation convention (the one thing to internalize)

`orientDeg` is **the angle from the Earth→Moon line to the ellipse's perigee
direction, measured in the rotating frame's rotation sense.** It is the only
planar orientation angle a GTO needs in this frame — the rotating frame
already fixes the Earth–Moon line as the reference, so there is no separate
argument of perigee. The physical GTO itself never changes (350 × 35,786 km:
`sma_km = 24446`, `ecc = 0.724781150290436`, identical on every sheet); only
its orientation relative to the Moon varies, sampled at {0°, 90°, 180°, 270°}.

Departure **phase** is a separate axis: the **time fraction since perigee**
(mean-anomaly fraction) over one Kepler period (0.4402 d). The departure
state is algebraic — `get_family_orbit('gto', dep_params)` builds the locus
with no propagation.

**Schema consequence:** `sheets(k).tauDRO` / `.tau_dep` hold **orientDeg in
degrees, NOT a period** (the legacy field name is kept so the sheet layout
matches the sibling catalogs). Everything in this folder handles that
correctly; the *other* deliverables' pickers/describe tools do **not** — see
Gotchas.

## What it covers

- **Departure:** one fixed 350 × 35,786 km GTO at 4 orientations
  {0°, 90°, 180°, 270°}.
- **Arrival:** tulip orbits Np ∈ {5, 7, 9, 12}, pm = −1 branch (pm = +1 is the
  exact z-mirror: flip z components of states and costates).
- **Phasing:** 12 departure × 6 arrival phase fractions per sheet.
- **Thrust:** rungs **[15 12 10 7 5] N** at Isp 1710 s, m0 150 kg. The 3–1 N
  legs of the siblings' 9-rung ladder produced **zero entries** here by both
  warm and cold recipes — a measured closure wall, split off as future work,
  not silently absent.
- **Totals:** 16 sheets (4 orientations × 4 Np), **2,625 entries**, 840/1,152
  phase pairs solved (73%), 136 full 5-rung ladders.
- **Ranges** (measured over all entries by `costate_catalog_extremes`):
  t_f 1.039–3.172 d; ΔV 8.11–26.48 km/s; m_f 30.9–92.5 kg. The fastest entry
  is at 15 N and the slowest at 5 N, but ΔV orders the other way (lightest at
  5 N, heaviest at 15 N): for an all-burn min-time family t_f shrinks slower
  than 1/T, so the burned impulse T·t_f grows with thrust.
- **Coverage by orientation:** 0° → 83%, 90° → 72%, 180° → 54%, 270° → 82%.
  **The π-dip:** apogee-toward-the-Moon (180°) is the hard corner at every
  petal count — genuine basin difficulty of that geometry, not a solver
  artifact (the cold-start defect that originally masked this was found and
  fixed; the dip survived the fix).

## Entry-field reference

Top level: `name`, `description`, `created`, `provenance`, `schema` (v2),
`constants` (`muStar`, `lStar_km`, `tStar_s`), `thruster` (`isp_s`, `m0_kg`,
`c_nd`, thrust-ND formula), `rungs_N`, `sheets(16)`, `n_entries`, `env`
(solver/version pinning), `derive` (formula strings for every derived
quantity), `usage`, `conj_test` (conjugate-sweep provenance).

Per sheet `sheets(k)`:

| Field | Meaning |
|---|---|
| `tauDRO`, `tau_dep` | **orientDeg (degrees)** — the sheet key (legacy names) |
| `tau_arr`, `period_tulip_nd` | tulip period (ND) |
| `dep_family`, `dep_params` | `'gto'`, `{sma_km, ecc, orientDeg}` — reconstruction recipe |
| `arr_family`, `arr_params` | `'tulip'`, `{Np, pm}` |
| `Np`, `pm` | petal count, branch sign (−1 stored) |
| `sD_frac [12]`, `sA_frac [6]` | departure / arrival phase fractions |
| `has_solution [12×6×5]` | **the** coverage authority (dep × arr × rung) |
| `tf_nd [12×6×5]` | minimum flight time (ND) where solved |
| `entry_index [12×6×5]` | column into `z8` for that cell |
| `z8 [8×n]` | converged `[λr; λv; λm; tf]` per entry (tfMin-ready) |
| `conj_pass/_ncross/_atfinal` | per-entry conjugate-point verdicts |

## How it was made + verification

Per (sheet, phase-pair): direct collocation solve at the 15 N anchor, walked
down the 5-rung ladder with a thrust-locked warm recipe (`tf0 = 0.30`; a cold
mop-up pass covered cells the warm recipe missed), each rung's costates
refined by multiple shooting (`ms_tfmin`, pumpkyn-STM Jacobian) and passed
through three gates:

1. **ms residual** converged,
2. **flown arrival** — the PMP flight from z8 lands on the target (<100 km gate),
3. **tfMin acceptance** — pumpkyn's own solver accepts z8 unchanged
   (|Δz| < 1e-6; re-solve fidelity measured 9.6e-16 – 1.1e-8 across all
   2,625 entries in the conjugate sweep).

Beyond the three gates, this is the first catalog in the series with a **100%
conjugate-point census**: 2,625/2,625 pass (`ms_conjugate_test`, K=24,
free-time quotiented Jacobi field; zero interior sign changes anywhere).
Verdicts ride in the catalog itself. Independent 3-entry replay audit across
3 orientations: flown miss < 0.01 km, acceptance |Δz| < 4e-10, Earth
clearance clean. Minimum time = minimum ΔV at fixed thrust (continuous burn —
a theorem for this problem class, not an observation).

## Gotchas, honestly

- **Do not mix this folder's helpers with the other deliverables'.** The
  period-keyed copies of `costate_catalog_pick` / `costate_lib_describe` /
  `costate_catalog_extremes` shipped with DRO/HALO/DPO/L1-L2 will silently
  mis-serve this catalog: their log-distance sheet metric returns the wrong
  sheet at `orientDeg = 0` (log 0 = −Inf → NaN, and `min` skips NaN), and
  they read the degree-valued key as an ND period (bogus day conversions).
  This folder's copies are orientDeg-aware (circular ±180° wrap metric, true
  Kepler period from `dep_params.sma_km`) and remain correct on the
  period-keyed catalogs.
- **Coverage gaps are real** (27% of pairs unsolved; 180° is the weak
  quadrant): consult `has_solution`, and heed the picker's warnings — they
  tell you what you actually got.
- **Rebuild endpoints with `'spline'` interpolation** (the example shows it).
  The min-time flow is sensitive: a 0.03 km linear-interp error in the
  departure state grows to ~500 km at arrival; spline endpoints reproduce
  the audit's <0.01 km flown miss.
- Flight-time interpolation between rungs is linear and for `tf` only;
  **z8 is never interpolated** — you always get a stored converged vector at
  the nearest rung (`info.delivered.seedThrustN` vs `tfThrustN`).
- The movie renderer carries an **ephemeris guard**: if the Aerospace
  Toolbox's `planetEphemeris` is on your path without its data package
  (`aeroDataPackage`), the renderer temporarily hides it so pumpkyn's
  analytic Earth–Moon model takes over (path restored on exit). Without the
  guard that configuration errors inside `ShowEarth`.
- The spec's flagship geometry (orientDeg = −25°) is reachable via the
  picker's nearest-sheet honesty warning (it maps to the 0° sheet through
  the circular metric), not a stored sheet.

Casey / Koblick, 2026-08-31.
