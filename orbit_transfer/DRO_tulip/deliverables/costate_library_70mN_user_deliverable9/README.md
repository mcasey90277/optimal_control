# DRO -> Tulip 70 mN Minimum-Time Costate Library

The library holds **1,152 primary transfers** (the fastest found in each cell of a 24 x 48 departure x arrival phase grid) and **656 alternatives** (every other transfer found at a phase pair), each with an optimality **status** saying what was established about it, and the multiple-shooting **junction states**.

This is the **user kit**: the library, the 9 functions that describe it, select a transfer, plot it and film it, and one worked example. Flights use pumpkyn's propagator from the stored junctions.

Catalog is `costate_catalog_dro_tulip_70mN`.

## Contents

| file | role |
|---|---|
| `costate_catalog_dro_tulip_70mN.mat` | the catalog (content keys `6a775f714b4cdf302a505aa56a8198ad` / `3a1c997b6def0ebabc1bfe044cb2ccd3`) |
| `get_family_orbit.m` | family name + parameters -> periodic orbit (the DRO and tulip recipes) |
| `costate_lib_describe.m` | fact sheet: grids in ND and days, coverage, status counts |
| `costate_catalog_pick.m` | five-coordinate lookup (DRO period, petals, departure/arrival phase in DAYS, thrust); says exactly what it delivers |
| `transfers_at.m` | every transfer recorded at a phase pair (fractions), primary and alternatives, ranked by t_f, with status |
| `status_key.m` | the status legend: codes, names, meanings |
| `catalog_transfer.m` | one transfer's history (time, state, thrust direction, mass; ND and days / km / kg) and its 3D plot |
| `costate_catalog_extremes.m` | min/max Delta-V or flight time over a slice (alternatives opt-in) |
| `plot_status_torus.m` | the status map of the phase torus |
| `costate_catalog_extremes_movies.m` | side-by-side movie of two transfers on one clock (pumpkynPie scene) |
| `costate_library_example.m` | worked example: describe, pick by days, tfMin accepts the pick, the cell's transfers, one plotted, the legend, extremes, the torus map, a movie |
| `status_torus.png` | the status map of the phase torus |
| `README.md` | this file |

Every function with a bare self-demo runs it on whichever `costate_catalog_*.mat` sits in the current folder (type e.g. `costate_lib_describe` with no arguments).

## Dependencies

- **pumpkyn + pumpkynPie** on the MATLAB path (run pumpkynPie's `startup.m`).
- **Optimization Toolbox**, only for `pumpkyn.cr3bp.tfMin`.

## Quick start

```matlab
costate_library_example   % describe, pick by days, tfMin, the cell's transfers, plot one, extremes, torus, movie
```

Or by hand:

```matlab
L = load('costate_catalog_dro_tulip_70mN.mat');
cat_ = L.costate_catalog_dro_tulip_70mN;
costate_lib_describe(cat_);                          % the fact sheet
[tf, z8, info] = costate_catalog_pick(cat_, 1, 7, 0.5, 3.0, 0.07);   % phases in DAYS; info.delivered says what you got
sh = cat_.sheets(1);
T = transfers_at(cat_, sh.sD_frac(1), sh.sA_frac(1));   % every transfer at a phase pair (FRACTIONS), with status
F = catalog_transfer(cat_, T(1), struct('plot', true));   % its history (days, km, kg) and a 3D plot
[eMin, eMax] = costate_catalog_extremes(cat_, struct('includeAlternatives', true), true, 'time');
plot_status_torus(cat_);                             % the status map
costate_catalog_extremes_movies(cat_, [], 'time', 'my_movie', struct('pair', {{T(1), T(end)}}));
```

## Handing an entry to pumpkyn tfMin

`z8 = [lambda_r; lambda_v; lambda_m; t_f]` is the initial guess; the endpoints are the two periodic orbits at the entry's phase fractions `sD`, `sA`. For a status-4 entry tfMin's first residual is already at machine level: it returns z8 unchanged or polished within a step (`|dz| = 0` at the example's pick). A status-1 entry is not expected to be accepted and can make tfMin grind for minutes.

```matlab
sh = cat_.sheets(1);                                % z8, sD, sA: a transfers_at row's .z8 .sD .sA (days / period_days for the picker's output)
pD = sh.dep_params;  pD.muStar = cat_.constants.muStar;
pA = sh.arr_params;  pA.muStar = cat_.constants.muStar;
[tD, rvD] = get_family_orbit(sh.dep_family, pD);   % one period of the DRO
[tT, rvT] = get_family_orbit(sh.arr_family, pA);   % one period of the tulip
rv0 = interp1(tD/tD(end), rvD, sD, 'spline');
rvf = interp1(tT/tT(end), rvT, sA, 'spline');
Tnd = (0.07/150)*cat_.constants.tStar_s^2/(cat_.constants.lStar_km*1000);   % thrust acceleration, ND
z = pumpkyn.cr3bp.tfMin(rv0, rvf, z8, Tnd, cat_.thruster.c_nd, cat_.constants.muStar);
```

## What it covers

| axis | values |
|---|---|
| departure orbit | `dro`, `get_family_orbit('dro', struct('tau', 1))`, period 1 ND |
| arrival orbit | `tulip`, `get_family_orbit('tulip', struct('Np', 7, 'pm', -1, 'tau', 5.23598775598299))`, period 5.23599 ND |
| thrust | 0.07 N, Isp 900 s, m0 150 kg |
| phase grid | 24 departure x 48 arrival phases (fractions of each orbit's period) |
| primaries | 1,152 of 1,152 cells |
| alternatives | 656 (372 at off-grid phases; 59 cells hold at least one, up to 11) |
| junction states | 1,152 of 1,152 primaries, 656 of 656 alternatives |
| t_f | primaries 16.14 .. 23.15 days, alternatives 16.23 .. 44.05 days |
| Delta-V | primaries 0.6759 .. 0.9866 km/s, alternatives 0.6800 .. 1.9831 km/s |

At fixed thrust the burn is continuous, so minimum time is minimum Delta-V. `mf = 1 - Tmax_nd*tf_nd/c_nd`, `dV = c_nd*log(1/mf)*lStar_km/tStar_s` (`cat_.derive`).

## The statuses

A status records what was established about a transfer, not what might be true.

| code | name | primaries | alternatives | what it establishes |
|---|---|---|---|---|
| 4 | sufficient | 1,152 | 268 | locally optimal: first- and second-order sufficient conditions hold |
| 3 | necessary only | 0 | 110 | first-order necessary conditions hold; local optimality not established |
| 2 | conjugate point found | 0 | 268 | a converged extremal with a conjugate point on the arc: not locally optimal |
| 1 | neither | 0 | 10 | flies to the target; the necessary conditions are not established |

A status-2 transfer is a converged extremal that tfMin accepts; it is not locally optimal. A status-1 transfer flies to the target but tfMin is not expected to accept it and can grind for minutes. `transfers_at` lists every transfer at a cell with its status; `costate_catalog_extremes` takes `filt.status` and `filt.includeAlternatives`. Alternatives 21, 219, 419, 491, 544 carry a borderline status, a lower bound (`alternatives.borderline`).

## Phase conventions

- `transfers_at`, `catalog_transfer`, `plot_status_torus` and the movie take phases as **fractions of a period** in [0, 1), wrapping.
- `costate_catalog_pick(cat, tauDep, Np, depDays, arrDays, thrustN)` takes phases in **DAYS**: fraction x period_nd x tStar_s / 86400.
- Quantities are ND unless the field name carries a unit (`tfDays`, `dvKms`; `catalog_transfer`'s `tDays`, `rKm`, `vKms`, `massKg`).

## Entry fields

`sheets(1)`:

| field | meaning |
|---|---|
| `tau_dep` | departure-orbit period (ND) |
| `tau_arr` | arrival-orbit period (ND) |
| `dep_family` | departure family name, for `get_family_orbit` |
| `dep_params` | departure reconstruction parameters |
| `arr_family` | arrival family name |
| `arr_params` | arrival reconstruction parameters |
| `Np` | tulip petal count |
| `pm` | tulip branch (the pm = +1 branch is the z-mirror) |
| `sD_frac` | departure phase grid, fractions of the departure period |
| `sA_frac` | arrival phase grid, fractions of the arrival period |
| `has_solution` | **the** coverage authority (departure x arrival) |
| `tf_nd` | the primary's flight time (ND) |
| `entry_index` | column of `z8` / `junctions` for each cell |
| `z8` | converged `[lam_r; lam_v; lam_m; tf]` per entry, tfMin-ready |
| `junctions` | multiple-shooting junction states per entry (state and costate rows x nodes) |
| `status` | optimality status code (`status_key`) |
| `status_reason` | one-line reason for the status |

`alternatives(k)`:

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
| `borderline` | true when the status is a lower bound |

Every other field is explained by the legends inside the catalog: `status_key`, `notes_key`, `derive`, `second_order.meaning`, `conj_test.meaning`, `hyp_gates.meaning`, `families.note`.

## The map

`status_torus.png`: per cell, the primary's flight time, how many alternatives it holds, and which cells hold a status-2 transfer (45 cells); the 372 off-grid alternatives are drawn at their own phases. `costate_catalog_extremes_movies` films any two transfers.
