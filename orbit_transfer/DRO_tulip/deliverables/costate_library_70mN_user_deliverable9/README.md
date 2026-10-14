# DRO -> Tulip 70 mN Minimum-Time Costate Library

The library holds **1,152 primary transfers** (the fastest found in each cell of a 24 x 48 departure x arrival phase grid) and **656 alternatives** (every other transfer the builds found at a phase pair), each with an optimality **status** saying what was established about it, and the multiple-shooting **junction states**.

This is the **user kit**: the library, the 8 functions that describe it, select a transfer and show it, and one worked example. The movie is generated with pumpkyn's propagator from the stored junctions.

Catalog is `costate_catalog_dro_tulip_70mN`.

## Contents

| file | role |
|---|---|
| `costate_catalog_dro_tulip_70mN.mat` | the catalog (cleaned export; content keys `6a775f714b4cdf302a505aa56a8198ad` / `3a1c997b6def0ebabc1bfe044cb2ccd3`) |
| `get_family_orbit.m` | family name + parameters -> periodic orbit (the DRO and tulip recipes) |
| `costate_lib_describe.m` | fact sheet: grids in ND and days, coverage, status counts |
| `costate_catalog_pick.m` | five-coordinate lookup (DRO period, petals, departure/arrival phase in DAYS, thrust) with the honesty contract; unchanged |
| `transfers_at.m` | every transfer recorded at a phase pair (fractions), primary and alternatives, ranked by t_f, with status |
| `status_key.m` | the status legend: codes, names, meanings, the order of the checks |
| `costate_catalog_extremes.m` | min/max Delta-V or flight time over a slice (alternatives opt-in) |
| `plot_status_torus.m` | the status map of the phase torus |
| `costate_catalog_extremes_movies.m` | side-by-side movie of two transfers on one clock (pumpkynPie scene, flown with pumpkyn's own propagator from the stored junctions) |
| `costate_library_example.m` | worked example: describe, pick by days, tfMin accepts the pick, every transfer at the cell, the legend, extremes, the torus map, a movie |
| `status_pair_movie.mp4`, `status_pair_movie.gif` | the sample movie (see below) |
| `status_torus.png` | the status map of the phase torus (see below) |
| `README.md` | this file |

Every function with a bare self-demo runs it on whichever `costate_catalog_*.mat` sits in the current folder (type e.g. `costate_lib_describe` with no arguments).

## Dependencies

- **pumpkyn + pumpkynPie** on the MATLAB path (run pumpkynPie's `startup.m`): the CR3BP dynamics and propagator, `pumpkyn.cr3bp.tfMin`, the orbit families, and the movie scene (`pumpkynPie.plot.SatelliteAnimator`).
- The kit's own functions need base MATLAB only; every other function the example calls is in this folder. `pumpkyn.cr3bp.tfMin` (the example's section 3 and the recipe below) needs the **Optimization Toolbox** through pumpkyn's `fsolve`. The movie writes MPEG-4 with `VideoWriter`, which MATLAB does not provide on Linux (change the profile to 'Motion JPEG AVI' in `costate_catalog_extremes_movies.m` there).

## Quick start

Run the worked example from this folder (it prints a verdict line per section):

```matlab
costate_library_example   % describe, pick by days, tfMin accepts it, the cell's transfers, extremes, torus, movie
```

Or by hand:

```matlab
L = load('costate_catalog_dro_tulip_70mN.mat');
cat_ = L.costate_catalog_dro_tulip_70mN;
costate_lib_describe(cat_);                          % the fact sheet
[tf, z8, info] = costate_catalog_pick(cat_, 1, 7, 0.5, 3.0, 0.07);   % phases in DAYS; info.delivered says what you got
sh = cat_.sheets(1);
T = transfers_at(cat_, sh.sD_frac(1), sh.sA_frac(1));   % every transfer at a phase pair (FRACTIONS), with status
[eMin, eMax] = costate_catalog_extremes(cat_, struct('includeAlternatives', true), true, 'time');
plot_status_torus(cat_);                             % the status map
costate_catalog_extremes_movies(cat_, [], 'time', 'my_movie', struct('pair', {{T(1), T(end)}}));
```

## Handing an entry to pumpkyn tfMin

The endpoints are states on the two periodic orbits at the entry's phases (fractions of each period), rebuilt from the catalog's own recipe; `z8 = [lambda_r; lambda_v; lambda_m; t_f]` is the initial guess. For a certified entry tfMin's first residual is already at machine level and the vector comes back unchanged (`|dz| = 0` at the example's pick) or polished within a step or two: the plain `interp1` spline below differs from the certifier's periodic spline by at most 1.9e-10 ND (7.3e-05 km) at this library's grid phases, measured at build. The example prints `|dz|`.

```matlab
sh = cat_.sheets(1);
[tD, rvD] = get_family_orbit(sh.dep_family, sh.dep_params);   % one period of the DRO
[tT, rvT] = get_family_orbit(sh.arr_family, sh.arr_params);   % one period of the tulip
rv0 = interp1(tD/tD(end), rvD, sD, 'spline');                % sD, sA = phase fractions of the entry
rvf = interp1(tT/tT(end), rvT, sA, 'spline');
Tnd = (0.07/150)*cat_.constants.tStar_s^2/(cat_.constants.lStar_km*1000);   % thrust acceleration, ND
z = pumpkyn.cr3bp.tfMin(rv0, rvf, z8, Tnd, cat_.thruster.c_nd, cat_.constants.muStar);
```

A status-1 entry (see the statuses) flies to the target but its necessary conditions were not established, so tfMin is not expected to accept it and can grind for minutes: there are 0 such primaries and 10 such alternatives.

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

At fixed thrust the burn is continuous, so minimum time is minimum Delta-V: the ranges above order the same way. Delta-V and mass follow `cat_.derive` (`mf = 1 - Tmax_nd*tf_nd/c_nd`, `dV = c_nd*log(1/mf)*lStar_km/tStar_s`).

## The statuses

A status records what was ESTABLISHED, not what might be true: a check that could not run gives the highest tier actually established, with the reason. An orbit transfer does not have to be optimal to be useful.

| code | name | primaries | alternatives | what it establishes |
|---|---|---|---|---|
| 4 | sufficient | 1,152 | 268 | first- and second-order sufficient conditions hold (the full certify_root stack passed) |
| 3 | necessary only | 0 | 110 | first-order necessary conditions hold (pointwise Pontryagin checks + tfMin witness); sufficiency not established |
| 2 | conjugate point found | 0 | 268 | an extremal with a conjugate point on the arc, under the hypotheses (gates + H6) that make it a refutation: NOT locally optimal |
| 1 | neither | 0 | 10 | flies to the target, but the necessary conditions are not established |

The checks, in order: certify_root: 1 polish, 2 flight (100 km / 10 m/s), 3 pointwise Pontryagin, 4 tfMin witness, 5 coarse conjugate test, 6 hypothesis gates (min\|lam_v\|, Q_mt, dim S, X2, lift margin), 7 H6, 8 dense conjugate scan.

A status-2 transfer (**conjugate point found**) is a converged extremal that tfMin accepts like any other; it is just not locally optimal. A faster status-2 alternative can still be the right choice for a mission, taken knowing exactly that. `transfers_at` lists every transfer at a cell with its status; `costate_catalog_extremes` takes `filt.status` and `filt.includeAlternatives`.

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

- `transfers_at`, `plot_status_torus` and the movie take phases as **fractions of a period** in [0, 1), wrapping (s = 1 - 1e-9 finds s = 0), matched with a circular tolerance.
- `costate_catalog_pick(cat, tauDep, Np, depDays, arrDays, thrustN)` takes phases in **DAYS** (unchanged from the earlier deliverables): fraction x period_nd x tStar_s / 86400.
- Flight times are ND unless the field name says otherwise (`tfDays`, `dvKms`).

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

## The movie and the map

`status_pair_movie.mp4` / `status_pair_movie.gif`: one cell's primary beside its fastest status-2 alternative, same phases, one clock (1280 x 720), flown from their junctions with pumpkyn's propagator.

`status_torus.png`: per cell, the primary's flight time, how many alternatives it holds, and which cells hold a status-2 transfer (45 cells); the 372 off-grid alternatives are drawn at their own phases.
