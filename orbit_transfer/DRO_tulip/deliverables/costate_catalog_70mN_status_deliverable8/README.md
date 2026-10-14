# DRO -> Tulip 70 mN Minimum-Time Costate Library, with Optimality Status (deliverable 8)

Every entry here is a **converged minimum-time PMP extremal**, not a hint: hand its `z8` to `pumpkyn.cr3bp.tfMin` at the entry's own endpoints and it comes back unchanged. The library holds **1,152 primary transfers** (the fastest found in each cell of a 24 x 48 departure x arrival phase grid) and **656 alternatives** (every other transfer the builds found at a phase pair), each with an optimality **status** saying what was established about it, and the multiple-shooting **junction states** that fly a 16.14-44.05 day arc accurately.

Catalog `costate_catalog_dro_tulip_70mN` built 2026-10-02; this package generated 2026-10-13 by `build_dro_deliverable.m`. Every number below is computed from the shipped catalog and the audits bound to it by content key, not transcribed.

## Contents

| file | role |
|---|---|
| `costate_catalog_dro_tulip_70mN.mat` | the catalog (cleaned export; content keys `6a775f714b4cdf302a505aa56a8198ad` / `3a1c997b6def0ebabc1bfe044cb2ccd3`) |
| `get_family_orbit.m` | family name + parameters -> periodic orbit (the DRO and tulip recipes) |
| `phase_state.m` | state at a phase fraction of a periodic orbit |
| `periodic_pp.m` | the periodic spline phase_state is built on |
| `validate_flight.m` | the one admissibility check for a flown trajectory |
| `scalar_verdict.m` | reads a solver result as a gate verdict, or refuses |
| `catalog_schema.m` | the catalog schema: field list, validator, formula registry |
| `fly_transfer.m` | flight engine: one flight from z8, with mass and Delta-V |
| `run_capped.m` | hard wall-clock cap on a parfeval worker (used by tfmin_accepts) |
| `current_pool.m` | the open parallel pool or [], without throwing when the toolbox is absent |
| `costate_catalog_pick.m` | five-coordinate lookup (DRO period, petals, departure/arrival phase in DAYS, thrust) with the honesty contract; unchanged |
| `transfers_at.m` | every transfer recorded at a phase pair (fractions), primary and alternatives, ranked by t_f, with status |
| `status_key.m` | the status legend: codes, names, meanings, the order of the checks |
| `catalog_endpoints.m` | the certified endpoint rule, from the catalog's own orbit recipes |
| `fly_catalog_entry.m` | fly an entry: from its junction states when it has them, else from z8 |
| `tfmin_accepts.m` | the "pumpkyn tfMin accepts it unchanged" check (capped; status-1 rows skipped) |
| `entry_thrust_program.m` | flight engine: the transfer and thrust program, segment by segment from the junctions |
| `costate_lib_describe.m` | fact sheet: grids in ND and days, coverage, status counts |
| `costate_catalog_extremes.m` | min/max Delta-V or flight time over a slice (alternatives opt-in) |
| `plot_status_torus.m` | the status map of the phase torus |
| `costate_catalog_extremes_movies.m` | side-by-side movie of two transfers on one clock (pumpkynPie scene) |
| `print_transfer_summary.m` | prints a transfer's time of flight, Delta-V, propellant, flight samples |
| `plot_transfer_3d.m` | the 3D transfer figure the examples draw |
| `sample_status_pair.m` | picks the cell with the most status-2 alternatives and its pair |
| `costate_catalog_example.m` | worked example: pick, derive units, rebuild orbits, fly, tfMin, plot |
| `status_layer_example.m` | worked example: every transfer at a cell, statuses, fly two, the torus map |
| `status_pair_movie.mp4`, `status_pair_movie.gif` | the sample movie (see below) |
| `status_torus.png` | the status map of the phase torus (see below) |
| `README.md` | this file |

Every helper with a bare self-demo runs it on whichever `costate_catalog_*.mat` sits in the current folder (type e.g. `costate_lib_describe` with no arguments).

## Dependencies

- **pumpkyn + pumpkynPie** on the MATLAB path (run pumpkynPie's `startup.m`): the CR3BP dynamics, `pumpkyn.cr3bp.tfMin`, the orbit families, and the movie scene (`pumpkynPie.plot.SatelliteAnimator`).
- **Parallel Computing Toolbox: optional.** `tfmin_accepts` runs tfMin under a hard wall-clock cap (`run_capped`) when a pool is open; without the toolbox or a pool it calls tfMin directly.
- Nothing else: every other function the examples and helpers call is in this folder. The `.mat` itself loads with `load` alone.

## Quick start

Run the two worked examples from this folder (each about 10 seconds):

```matlab
costate_catalog_example   % pick by phases (days), fly, tfMin accepts it, plot
status_layer_example      % every transfer at one cell, statuses, fly two, the torus map
```

Or by hand:

```matlab
L = load('costate_catalog_dro_tulip_70mN.mat');
cat_ = L.costate_catalog_dro_tulip_70mN;
T = transfers_at(cat_, 0.0000, 0.0129);   % phases as fractions of a period
F = fly_catalog_entry(cat_, T(1));          % from its junctions
[ok, dz] = tfmin_accepts(cat_, T(1));       % pumpkyn tfMin accepts it unchanged
```

## What it covers

| axis | values |
|---|---|
| departure orbit | `dro`, recipe `get_family_orbit('dro', struct('tau', 1))`, period 1 ND |
| arrival orbit | `tulip`, recipe `get_family_orbit('tulip', struct('Np', 7, 'pm', -1, 'tau', 5.23598775598299))`, period 5.23599 ND |
| thrust | 0.07 N (one rung), Isp 900 s, m0 150 kg |
| phase grid | 24 departure x 48 arrival phases (fractions of each orbit's period) |
| primaries | 1,152 of 1,152 cells (100.0%) |
| alternatives | 656 (372 off the grid; 59 cells hold at least one, up to 11 per cell) |
| junction states | 1,152 of 1,152 primaries, 656 of 656 alternatives |
| t_f, primaries | 16.14 .. 23.15 days |
| t_f, alternatives | 16.23 .. 44.05 days |
| Delta-V, primaries | 0.6759 .. 0.9866 km/s (final mass 134.135 .. 138.941 kg) |
| Delta-V, alternatives | 0.6800 .. 1.9831 km/s (final mass 119.815 .. 138.877 kg) |

At fixed thrust the burn is continuous, so minimum time is minimum Delta-V: the ranges above order the same way. Delta-V and mass follow `cat.derive` (`mf = 1 - Tmax_nd*tf_nd/c_nd`, `dV = c_nd*log(1/mf)*lStar_km/tStar_s`).

## The statuses

A status records what was ESTABLISHED, not what might be true: a check that could not run gives the highest tier actually established, with the reason. An orbit transfer does not have to be optimal to be useful.

| code | name | primaries | alternatives | what it establishes |
|---|---|---|---|---|
| 4 | sufficient | 1,152 | 268 | first- and second-order sufficient conditions hold (the full certify_root stack passed) |
| 3 | necessary only | 0 | 110 | first-order necessary conditions hold (pointwise Pontryagin checks + tfMin witness); sufficiency not established |
| 2 | conjugate point found | 0 | 268 | an extremal with a conjugate point on the arc, under the hypotheses (gates + H6) that make it a refutation: NOT locally optimal |
| 1 | neither | 0 | 10 | flies to the target, but the necessary conditions are not established |

The checks, in order: certify_root: 1 polish, 2 flight (100 km / 10 m/s), 3 pointwise Pontryagin, 4 tfMin witness, 5 coarse conjugate test, 6 hypothesis gates (min|lam_v|, Q_mt, dim S, X2, lift margin), 7 H6, 8 dense conjugate scan.

A status-2 transfer (**conjugate point found**) is a converged extremal that tfMin accepts like any other; it is just not locally optimal. A faster status-2 alternative can still be the right choice for a mission, taken knowing exactly that.

### Borderline rows

5 alternatives are **borderline**: a borderline status (alternatives.borderline) is a lower bound: the status flipped between runs at a threshold, the claim both runs support is kept (status_meet on the lattice 1 < 3 < 2, 1 < 3 < 4), and the audit accepts the same root at any status whose meet with it is it.

| index | sD | sA | status | reason (first clause) |
|---|---|---|---|---|
| 21 | 0.0000 | 0.0233 | 3 | borderline: stored 2, re-audit 3 -> 3 (the common ground of the two runs) -- status not reproduced: stored 2, now 3 (conjugate test verdict 0 (FAIL)) |
| 219 | 0.0000 | 0.4087 | 1 | borderline: stored 3, re-audit 1 -> 1 (the common ground of the two runs) -- status not reproduced: stored 3, now 1 (tfMin witness exceeded its 300 s cap) |
| 419 | 0.0000 | 0.7525 | 3 | borderline: stored 3, re-audit 2 -> 3 (the common ground of the two runs) -- status not reproduced: stored 3, now 2 (conjugate test verdict 0 (FAIL)) |
| 491 | 0.0000 | 0.8879 | 3 | borderline: stored 3, re-audit 2 -> 3 (the common ground of the two runs) -- status not reproduced: stored 3, now 2 (conjugate test verdict 0 (FAIL)) |
| 544 | 0.0000 | 0.9504 | 3 | borderline: stored 3, re-audit 2 -> 3 (the common ground of the two runs) -- status not reproduced: stored 3, now 2 (conjugate test verdict 0 (FAIL)) |

## Phase conventions

- `transfers_at`, `catalog_endpoints`, `fly_catalog_entry`, `tfmin_accepts` and `plot_status_torus` take phases as **fractions of a period** in [0, 1), wrapping (s = 1 - 1e-9 finds s = 0), matched with a circular tolerance.
- `costate_catalog_pick(cat, tauDep, Np, depDays, arrDays, thrustN)` takes phases in **DAYS** (unchanged from the earlier deliverables): fraction x period_nd x tStar_s / 86400.
- Flight times are ND unless the field name says otherwise (`tfDays`, `dvKms`, `flyKm`, `flyVms`).

## Entry-field reference

Top level: `name`, `description`, `created`, `provenance`, `schema`, `constants`, `thruster`, `rungs_N`, `notes_key`, `families`, `sheets`, `n_entries`, `env`, `derive`, `conj_test`, `hyp_gates`, `usage`, `second_order`, `alternatives`, `status_key`, `build`.

Per sheet `sheets(1)`:

| field | meaning |
|---|---|
| `tau_dep` | departure-orbit period (ND) |
| `tau_arr` | arrival-orbit period (ND) |
| `tauDRO` | departure-orbit period (ND; legacy name, any family) |
| `dep_family` | departure family name, for `get_family_orbit` |
| `dep_params` | departure reconstruction parameters |
| `arr_family` | arrival family name |
| `arr_params` | arrival reconstruction parameters |
| `Np` | tulip petal count |
| `pm` | tulip branch (the pm = +1 branch is the z-mirror) |
| `period_tulip_nd` | tulip period (ND) |
| `sD_frac` | departure phase grid, fractions of the departure period |
| `sA_frac` | arrival phase grid, fractions of the arrival period |
| `has_solution` | **the** coverage authority (departure x arrival) |
| `tf_nd` | the primary's flight time (ND) |
| `entry_index` | column of `z8` / `junctions` for each cell |
| `z8` | converged `[lam_r; lam_v; lam_m; tf]` per entry, tfMin-ready |
| `conj_pass` | free-time conjugate test at certification (1 = none found) |
| `gate_min_lamv` | hypothesis gate: min magnitude of lam_v over the arc (> 0) |
| `gate_min_qmt` | hypothesis gate: min Q_mt (> 0, all-burn is the PMP control) |
| `gate_dimS` | hypothesis gate: lift-space dimension (1 = no abnormal lift) |
| `conj_K` | conjugate-test sample count |
| `family_index` | provenance: which continuation family found the entry (see `families.note`; read by no tool) |
| `entry_notes` | per-entry story, layout in `notes_key` |
| `conj_interior` | dense scan: interior sign changes (0 = none) |
| `conj_multiplicity` | dense scan: candidates with two or more collapsed singular values |
| `conj_interior_cand` | dense scan: located interior candidates |
| `conj_near_miss` | dense scan: near-misses (positive minimum located) |
| `conj_zero` | dense scan: zeros (at the floor or a sign change) |
| `conj_unresolved` | dense scan: unresolved candidates (inside the floor band) |
| `h6_margin` | H6 margin, (c/T)/lambda_m(0) (> 1 excludes the spurious-zero mechanism) |
| `lift_margin` | sigma_6 over the measured lift-matrix error (> 1 certifies dim S = 1) |
| `conj_min_rel` | closest near-miss, relative to the median sigma_6 |
| `status` | optimality status code (`status_key`) |
| `status_reason` | one-line reason for the status |
| `junctions` | multiple-shooting junction states per entry (state and costate rows x nodes) |

Per alternative `alternatives(k)`:

| field | meaning |
|---|---|
| `sD` | departure phase (fraction) |
| `sA` | arrival phase (fraction) |
| `iD` | departure grid index of its cell (0 = off the grid) |
| `iA` | arrival grid index of its cell (0 = off the grid) |
| `z8` | converged `[lam_r; lam_v; lam_m; tf]` |
| `junctions` | multiple-shooting junction states |
| `tf_nd` | flight time (ND) |
| `status` | optimality status code (`status_key`) |
| `status_reason` | one-line reason for the status |
| `inferred` | true when the status was classified from the certifier's record rather than read from its stamp |
| `conj` | conjugate test verdict (1 pass, 0 conjugate point, -1 undecided) |
| `conjVerdict` | conjugate test verdict, text |
| `minLamV` | hypothesis gate: min magnitude of lam_v over the arc |
| `minQmt` | hypothesis gate: min Q_mt |
| `dimS` | hypothesis gate: lift-space dimension |
| `h6Margin` | H6 margin |
| `liftMargin` | lift margin |
| `flyKm` | flown arrival miss, km (NaN = not measured) |
| `flyVms` | flown arrival miss, m/s (NaN = not measured) |
| `source` | how it was found (kind only, paths stripped) |
| `sheet` | sheet index (1) |
| `borderline` | true when the status is the meet of two runs (a lower bound) |

Legends that ride in the catalog: `status_key` (codes, names, meanings, check order), `notes_key` (the `entry_notes` layout), `second_order.meaning`, `conj_test.meaning`, `hyp_gates.meaning`, `families.note`, `derive` (every derived-quantity formula).

## How it was made, and how it was checked

Arrival axis by pseudo-arclength continuation (arclength_ms / arclength_arrival), departure axis by the bisecting walker (rib_from_crossing); assembled by sheet_from_arcs with the certified library as seeds; FINDINGS 37-38, 2026-09-09. Status layer (FINDINGS 98): candidates harvested from the builds, every inferred one re-certified, primaries matched to their junctions by root, assembled beside the unchanged primaries, audited fail-closed (CLEAN), borderline statuses relabelled to the meet of two runs. Exported for deliverable 8 on 2026-10-13.

Build record: 24 x 48 build 2026-09-20/21; refill 2026-10-02/03; merge with the 24 x 24 record 2026-10-03; status-layer backfill, re-certification and audit 2026-10-04/05. 2,944 candidates harvested, 653 re-certified.

Every transfer goes through the certifier's stack in order (`status_key.order`); its status records how far it got: certify_root: 1 polish, 2 flight (100 km / 10 m/s), 3 pointwise Pontryagin, 4 tfMin witness, 5 coarse conjugate test, 6 hypothesis gates (min|lam_v|, Q_mt, dim S, X2, lift margin), 7 H6, 8 dense conjugate scan.

Primaries, from the sheet grids:

- free-time conjugate test passed: 1,152 / 1,152 (`conj_pass = 1`)
- hypothesis gates: min magnitude of lam_v > 0 on 1,152, Q_mt > 0 on 1,152, dim S = 1 on 1,152
- dense conjugate scan: 0 with an interior sign change; smallest H6 margin 5.32, smallest lift margin 10.6

**Audits, bound to this file by content key** (fail-closed re-derivation from the catalog's own keys):

- primaries: 1,152 ok / 0 bad (audit of 03-Oct-2026 14:27:34)
- alternatives: 656 ok / 0 bad, over 8 audit chunks covering every row exactly once
- total: **1,808 ok / 0 bad**; catalog key `6a775f714b4cdf302a505aa56a8198ad`, alternatives key `3a1c997b6def0ebabc1bfe044cb2ccd3`

Worst flown miss recorded on the alternatives (`flyKm`/`flyVms`, 656 of 656 measured; the certifier's flight gate is 100 km / 10 m/s): **12.6 km / 0.84 m/s**. `fly_catalog_entry` flies from the stored junctions when an entry has them rather than one flight from `z8` over 16.14-44.05 days; the examples print the miss they get.

## The movie and the map

`status_pair_movie.mp4` / `status_pair_movie.gif`: one cell's primary beside its fastest status-2 alternative, same phases, one clock (1280 x 720), flown from their junctions.

`status_torus.png`: per cell, the primary's flight time, how many alternatives it holds, and which cells hold a status-2 transfer (45 cells); the 372 off-grid alternatives are drawn at their own phases.

## Gotchas, honestly

- **Interpolation between grid points is NOT included.** The junction-interpolation experiment is pending; until its hit rate is measured on this library, interpolate at your own risk and use the result only as a tfMin seed. Neighbouring cells can sit on different solution families (`family_index`).
- **Status-1 rows and tfMin.** Status 1 means the transfer flies to the target but the necessary conditions were not established, so tfMin is not expected to accept it: `tfmin_accepts` skips status-1 rows (ok = false, expected refusal) unless `opts.force` is set. There are 0 such primaries and 10 such alternatives.
- **372 of the 656 alternatives are OFF-GRID** (`iD = 0` or `iA = 0`: 63 off the departure grid, 309 off the arrival grid, 0 off both): found by ribs and sheet candidates at non-grid phases. `transfers_at` finds them only at their own phases; the torus map draws them as markers.
- **`source` wording is kind-only** (`sheet candidate, column n (re-certified)`, `rib ... point`, ...): machine paths were stripped at export; the content keys are unchanged because `source` is not hashed.
- **Phases: fractions vs DAYS** (see Phase conventions): the picker alone takes DAYS.
- **Borderline statuses are lower bounds** (the 5 rows above): the status flipped between two runs at a threshold and the claim both runs support is the one stored.

M. Casey / D. Koblick, Coorbital Inc.
