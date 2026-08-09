# costate_catalog_halo_halo_deliverable6 — L1 ↔ L2 Halo-to-Halo Minimum-Time Costate Catalogs

Every entry is a **converged minimum-time PMP solution**, not a hint:
`pumpkyn.cr3bp.tfMin` accepts each one **unchanged** (|Δz| < 1e-6, observed
essentially 0). These are the first catalogs on the **schema-v2
arrival-period axis**: both ends are halos, so the picker's third argument is
the **arrival period (ND)**, not a petal count. Two catalogs ship — the
directions are independent problems (no symmetry maps one onto the other).

## Contents

| File | Role |
|---|---|
| `costate_catalog_halo_halo.mat` | **L1 → L2**: 8 sheets, 1,952 entries, 258/288 pairs (90%) |
| `costate_catalog_halo_halo_B.mat` | **L2 → L1**: 8 sheets, 2,096 entries, 277/288 pairs (96%) |
| `costate_catalog_pick.m` | Five-coordinate lookup; family-aware (arrival-period keying here); honesty-contract warnings |
| `costate_catalog_hh_example.m` | Worked example: validate → pick → derive via registry → rebuild both halos → fly → tfMin → reverse direction |
| `catalog_schema.m` | Schema validator + **named derive registry** (days/thrust/mass/ΔV formulas as code, one home) |
| `costate_lib_describe.m` | Coverage/grids/field reference for any catalog |
| `costate_catalog_extremes.m` | Min/max ΔV or flight time over filtered slices (use `.tauArr`, not `.Np`, here) |
| `costate_catalog_extremes_movies.m` | Side-by-side pumpkyn-style movies of the extremes |
| `get_family_orbit.m` | Family name + params → propagated periodic orbit (reconstruction for both ends) |
| `hh_extremes_movie.mp4/.gif` | Sample: cheapest vs costliest L1→L2 transfer, shared clock |

All helper self-demos are catalog-agnostic: typed bare, they run on whichever
`costate_catalog_*.mat` sits in the current folder.

## Dependencies

pumpkyn + pumpkynPie on the MATLAB path (run their `startup.m`). The `.mat`s
are dependency-free and self-describing (recipes + environment pinning inside).

## Quick start

```matlab
L = load('costate_catalog_halo_halo.mat');
cat_ = L.costate_catalog_halo_halo;
% third argument = ARRIVAL PERIOD (ND), not petal count:
[tf_nd, z8, info] = costate_catalog_pick(cat_, 2.7433, 2.8, 6.0, 3.0, 2.0);
% then run costate_catalog_hh_example.m end-to-end
```

## What it covers

- **Departures (A) / arrivals (B):** the complete admissible L1 southern set —
  only two members exist under the ≥500 km periselene / ≤100 Mm criteria:
  τ = 1.8037 ND (8.0 d, periselene 2,555 km) and τ = 2.7433 ND (12.2 d).
- **Arrivals (A) / departures (B):** L2 southern halos τ ∈ {1.75, 2.2, 2.8,
  3.4} ND (7.8–15.1 d). Northern pairs are exact z-mirrors: flip the z
  components of states and costates.
- **Phasing:** 6×6 grid per sheet. **Thrust:** [15 12 10 7 5 3 2 1.5 1] N at
  Isp 1710 s, m0 150 kg.
- **Ranges (L1→L2):** ΔV 1.01–8.40 km/s; flight 0.26–~3.1 d. Cheapest:
  τ=1.80 → τ=1.75 at 1 N, 1.0129 km/s, 8.8 kg. Fastest: same pair at 15 N,
  6.3 hours. (L2→L1: cheapest 1.0469 km/s, worst 10.16 km/s.)
- Solvability climbs with the far/circular L1 (τ=2.7433: 35–36/36 pairs
  everywhere); the Moon-diving τ=1.8037 is the hard end (26–33/36).

## Schema v2 notes (read this)

- `sheets(k).tau_dep` / `tau_arr` are the keys (**requested** family
  parameters); `tauDRO` is a legacy alias of `tau_dep`; **`Np` is NaN**
  (petal count is tulip-only) — filter extremes by `.tauArr`, never `.Np`.
- Full reconstruction recipes ride in `dep_family`/`dep_params` and
  `arr_family`/`arr_params` — use `get_family_orbit` for both ends.
- Derived values (days, ΔV, masses) come from the **named registry**:
  `catalog_schema('derive', cat, 'deltaV_kms', struct('tf_nd',...,'thrustN',...))`.
- `cat.env` records the MATLAB release and pumpkyn/pumpkynPie git revisions
  the sheets were built with.
- Validate any catalog with `catalog_schema('validate', cat)` → `{}` = clean.

## How it was made + verification

Per (sheet, phase-pair): direct collocation solve anchored cold at 15 N,
warm-started down the rung ladder; costates refined by multiple shooting
(pumpkyn-STM Jacobian) and passed through three gates: (1) ms residual
converged, (2) flown arrival <100 km (observed 0.00 km), (3) tfMin
acceptance unchanged (|Δz| < 1e-6, observed ~0). Only three-gate entries are
in the catalogs. Minimum time = minimum ΔV at fixed thrust (continuous burn —
a theorem for this problem class). A Jacobi (conjugate-point) necessary-
condition check passes on the golden benchmark cells of this pipeline.

## Gotchas, honestly

- **Coverage gaps are real** (A: 90%, B: 96% of pairs): `has_solution` is
  the authority; the picker warns and delivers the nearest solved pair.
- The directions are **separate catalogs** — a reverse trip needs the B
  catalog, and its tf genuinely differs from the forward trip's.
- Phasing dominates cost at high thrust: at 15 N the ΔV spread across the
  torus reaches several km/s; plan phasing before thrust.
- z8 is never interpolated: `info.delivered.seedThrustN` tells you the rung
  your seed actually came from; fly it there or let tfMin converge the
  correction at your thrust.

Casey / Koblick, 2026-08-09.
