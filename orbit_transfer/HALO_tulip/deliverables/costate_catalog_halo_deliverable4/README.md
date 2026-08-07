# costate_catalog_halo_deliverable4 — HALO → Tulip Minimum-Time Costate Catalog

Every entry in this catalog is a **converged minimum-time PMP solution**, not a
hint: `pumpkyn.cr3bp.tfMin` accepts each one **unchanged** (|Δz| < 1e-6, measured
max ~1e-9). Give the picker five coordinates — halo period, tulip petal count,
departure phase, arrival phase, thrust — and you get the flight time and the
full costate/final-time vector `z8 = [λr(3); λv(3); λm; tf]` ready to fly or to
seed `tfMin` at your own endpoints.

## Contents

| File | Role |
|---|---|
| `costate_catalog_halo_tulip.mat` | The catalog: 16 sheets, 3,980 entries (287 KB, compact ND format) |
| `costate_catalog_pick.m` | Five-coordinate lookup with the honesty contract (warns + reports what you actually got) |
| `costate_catalog_halo_example.m` | Worked example: pick → derive units → rebuild orbits → fly → tfMin → plot |
| `costate_lib_describe.m` | Prints a catalog's coverage, grids, and field reference |
| `costate_catalog_extremes.m` | Min/max ΔV or flight time over any filtered slice, with plot |
| `costate_catalog_extremes_movies.m` | Renders the two extremes as side-by-side pumpkyn-style movies |
| `get_family_orbit.m` | Family name + params → propagated periodic orbit (halo reconstruction needs it) |
| `halo_extremes_movie.mp4/.gif` | Sample: cheapest (0.65 km/s) vs costliest (7.28 km/s) transfer, shared clock |

## Dependencies

pumpkyn + pumpkynPie on the MATLAB path (run their `startup.m`). The `.mat`
itself is dependency-free and self-describing.

## Quick start

```matlab
L = load('costate_catalog_halo_tulip.mat');
cat_ = L.costate_catalog_halo_tulip;
[tf_nd, z8, info] = costate_catalog_pick(cat_, 2.8, 9, 4.0, 10.0, 5.0);
% then run costate_catalog_halo_example.m end-to-end
```

## What it covers

- **Departure:** Earth–Moon L2 SOUTHERN halos (pm = −1), periods τ ∈ {1.75,
  2.2, 2.8, 3.4} ND = 7.8–15.1 days, from the measured admissible box
  (≥500 km periselene, within 100 Mm of the Moon). The northern branch is the
  exact z-mirror: flip the z components of states and costates.
- **Arrival:** tulip orbits Np ∈ {5, 7, 9, 12}, pm = −1 branch (same mirror rule).
- **Phasing:** 6×6 grid per sheet (departure × arrival phase fractions).
- **Thrust:** rungs [15 12 10 7 5 3 2 1.5 1] N at Isp 1710 s, m0 150 kg.
- **Totals:** 16 sheets, 3,980 entries, 530/576 phase pairs solved (92%),
  242 full 9-rung ladders.
- **Ranges:** ΔV 0.65–7.28 km/s; flight time 0.68–20+ days. Cheapest entry:
  τ=1.75 halo → Np=12 tulip at 1.5 N, 0.6546 km/s, 5.74 kg propellant.
  Solvability improves with halo period (τ=3.4 sheets: 35–36/36 pairs).

## Schema notes (read this)

- `sheets(k).tauDRO` is a **legacy field name** shared with the DRO catalog so
  the same pickers work on both: **here it holds the HALO period (ND)**.
- `sheets(k).dep_family` / `dep_params` carry the full reconstruction recipe —
  a halo needs (tau, Lpt, pm), not just its period. Use `get_family_orbit`.
- Compact format: only canonical ND quantities are stored. Days, ΔV, masses
  are derived via the formula strings in `cat_.derive` (the example shows each).
- `sheets(k).has_solution` is **the** coverage authority — 8% of phase pairs
  have no solution at any rung; the picker warns and gives you the nearest
  solved pair.

## How it was made + verification

Per (sheet, phase-pair): direct collocation solve anchored cold at 15 N, then
warm-started down the rung ladder; each rung's costates refined by multiple
shooting (ms_tfmin, pumpkyn-STM Jacobian) and passed through three gates:

1. **ms residual** converged,
2. **flown arrival** — the PMP flight from z8 lands on the target (<100 km gate),
3. **tfMin acceptance** — pumpkyn's own single-shooting solver accepts z8
   unchanged (|Δz| < 1e-6; observed max ~1e-9).

Only three-gate entries are in the catalog. Minimum time = minimum ΔV at fixed
thrust here (continuous burn — a theorem for this problem class, not an
observation).

## Gotchas, honestly

- **Coverage gaps are real** (530/576 pairs): consult `has_solution` before
  trusting a request; the picker's warnings tell you what you actually got.
- Flight-time interpolation between thrust rungs is linear and for `tf` only;
  **z8 is never interpolated** — you always get a stored converged vector at
  the nearest rung (`info.delivered.seedThrustN` vs `tfThrustN`).
- Phases are days past each orbit's reference point (getHalo / getTulip
  outputs), modulo the period.
- The short flight times at high thrust (sub-day) mean departure phase matters
  enormously there: the ΔV spread across phasing at 15 N is several km/s, and
  it shrinks as thrust drops.

Casey / Koblick, 2026-08-07.
