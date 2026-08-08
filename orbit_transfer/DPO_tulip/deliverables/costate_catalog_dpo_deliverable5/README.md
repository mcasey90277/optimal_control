# costate_catalog_dpo_deliverable5 — DPO → Tulip Minimum-Time Costate Catalog

Every entry in this catalog is a **converged minimum-time PMP solution**, not a
hint: `pumpkyn.cr3bp.tfMin` accepts each one **unchanged** (|Δz| < 1e-6, observed
essentially 0). Give the picker five coordinates — DPO period, tulip petal
count, departure phase, arrival phase, thrust — and you get the flight time and
the full costate/final-time vector `z8 = [λr(3); λv(3); λm; tf]` ready to fly
or to seed `tfMin` at your own endpoints.

## Contents

| File | Role |
|---|---|
| `costate_catalog_dpo_tulip.mat` | The catalog: 16 sheets, 3,932 entries (compact ND format) |
| `costate_catalog_pick.m` | Five-coordinate lookup with the honesty contract (warns + reports what you actually got) |
| `costate_catalog_dpo_example.m` | Worked example: pick → derive units → rebuild orbits → fly → tfMin → plot |
| `costate_lib_describe.m` | Prints a catalog's coverage, grids, and field reference |
| `costate_catalog_extremes.m` | Min/max ΔV or flight time over any filtered slice, with plot |
| `costate_catalog_extremes_movies.m` | Renders the two extremes as side-by-side pumpkyn-style movies |
| `get_family_orbit.m` | Family name + params → propagated periodic orbit |
| `dpo_extremes_movie.mp4/.gif` | Sample: cheapest (0.76 km/s) vs costliest (7.10 km/s) transfer, shared clock |

All helper self-demos are catalog-agnostic: typed bare (e.g.
`costate_catalog_extremes`), they run on whichever `costate_catalog_*.mat`
sits in the current folder — this one included.

## Dependencies

pumpkyn + pumpkynPie on the MATLAB path (run their `startup.m`). The `.mat`
itself is dependency-free and self-describing.

## Quick start

```matlab
L = load('costate_catalog_dpo_tulip.mat');
cat_ = L.costate_catalog_dpo_tulip;
[tf_nd, z8, info] = costate_catalog_pick(cat_, 3.0, 9, 4.0, 10.0, 5.0);
% then run costate_catalog_dpo_example.m end-to-end
```

## What it covers

- **Departure:** Earth–Moon distant prograde orbits (planar family), periods
  τ ∈ {1, 2, 3, 4} ND = 4.4–17.7 days, from the measured admissible box
  (survey: τ 0.048–4.695 ND under the ≥500 km periselene / ≤100 Mm criteria;
  the upper edge is the 100 Mm criterion — DPOs grow with period).
- **Arrival:** tulip orbits Np ∈ {5, 7, 9, 12}, pm = −1 branch (pm = +1 is the
  exact z-mirror: flip z components of states and costates).
- **Phasing:** 6×6 grid per sheet (departure × arrival phase fractions).
- **Thrust:** rungs [15 12 10 7 5 3 2 1.5 1] N at Isp 1710 s, m0 150 kg.
- **Totals:** 16 sheets, 3,932 entries, 511/576 phase pairs solved (89%),
  238 full 9-rung ladders.
- **Ranges:** ΔV 0.76–7.10 km/s; flight time 0.24–3.58 days. Cheapest entry:
  τ=1 DPO → Np=12 tulip at 1 N, 0.7600 km/s, 6.65 kg propellant. The hard
  corner is the smallest DPO (τ=1) against the longest-period tulip (Np=12):
  27/36 pairs there vs 35/36 at the easy corners.

## Schema notes (read this)

- `sheets(k).tauDRO` is a **legacy field name** shared with the DRO catalog so
  the same pickers work on both: **here it holds the DPO period (ND)**, which
  is also the `getDPO` selector.
- `sheets(k).dep_family` / `dep_params` carry the reconstruction recipe; use
  `get_family_orbit('dpo', dep_params)`.
- Compact format: only canonical ND quantities are stored. Days, ΔV, masses
  are derived via the formula strings in `cat_.derive` (the example shows each).
- `sheets(k).has_solution` is **the** coverage authority — 11% of phase pairs
  have no solution at any rung; the picker warns and gives you the nearest
  solved pair.

## How it was made + verification

Per (sheet, phase-pair): direct collocation solve anchored cold at 15 N, then
warm-started down the rung ladder; each rung's costates refined by multiple
shooting (ms_tfmin, pumpkyn-STM Jacobian) and passed through three gates:

1. **ms residual** converged,
2. **flown arrival** — the PMP flight from z8 lands on the target (<100 km gate),
3. **tfMin acceptance** — pumpkyn's own single-shooting solver accepts z8
   unchanged (|Δz| < 1e-6; observed essentially 0).

Only three-gate entries are in the catalog. Minimum time = minimum ΔV at fixed
thrust here (continuous burn — a theorem for this problem class, not an
observation).

## Gotchas, honestly

- **Coverage gaps are real** (511/576 pairs): consult `has_solution` before
  trusting a request; the picker's warnings tell you what you actually got.
- Flight-time interpolation between thrust rungs is linear and for `tf` only;
  **z8 is never interpolated** — you always get a stored converged vector at
  the nearest rung (`info.delivered.seedThrustN` vs `tfThrustN`).
- Phases are days past each orbit's reference point (getDPO / getTulip
  outputs), modulo the period.
- Some low rungs are rejected pre-flight by the 500 km lunar-altitude floor
  (long spirals dip toward the Moon); the ladder stops a cell's descent there
  by design, which is one source of the partial ladders.

Casey / Koblick, 2026-08-08.
