# DRO to ELFO minimum-time transfers

Low-thrust minimum-time transfers from a distant retrograde orbit (DRO) to an
elliptical lunar frozen orbit (ELFO) in the Earth-Moon CR3BP, built to answer
one question: how much harder or easier is DRO to ELFO than DRO to tulip?
Design: `docs/superpowers/specs/2026-11-19-dro-elfo-design.md`.

Regime: 150 kg, Isp 1710 s, thrust ladder 15 N down to 0.5 N (the catalog
regime, chosen because a matched DRO to tulip data set exists there; it is not
a vehicle design).

## The ELFO here

A Keplerian orbit about the Moon, default the IM shape from
`proj7/pipelines/gdop/im_elfo_optimum` (semi-major axis 12,000 km,
eccentricity 0.69, inclination 56.5 deg, argument of perilune 90 deg; period
0.308 ND = 1.365 d, radius 3,720 to 20,280 km). It is not a periodic orbit of
the CR3BP, so the arrival is a fixed curve of osculating states at one plane
orientation `raanDeg` (node angle in the rotating frame at the arrival
instant), parametrized by time since perilune. A catalog is a stack of phase
tori (DRO phase x ELFO phase), one sheet per (DRO size, `raanDeg`).

## Quick start

```matlab
cd orbit_transfer/DRO_ELFO
dro_elfo_transfer        % edit the knobs at the top: which DRO, which ELFO, phases, thrust
```

Default run (DRO period 1 ND to IM ELFO satellite 1, 5 N): 0.805 d, final mass
129.25 kg, delta-V 2.50 km/s, every gate and the conjugate test pass. The
script prints a verdict table, saves the result and the library-format entry
under `results/`, and draws the figure. One commented knob line switches the
target to the tulip.

## Files

| file | what it does |
|---|---|
| `dro_elfo_transfer.m` | user script: one transfer, verdicts, save, plot |
| `compare_dro_targets.m` | the matched comparison of tulip and ELFO sheets |
| `catalog/run_dro_elfo_torus.m`, `run_dro_elfo_batched.sh` | fleet front door and its resumable batched driver (phase-torus sheets) |
| `catalog/build_dro_elfo_catalog.m` | the chain: census, clearance, repairs, package, conjugate test, gates, audit, pictures |
| `catalog/dro_elfo_pick.m` | pick an entry by (DRO size, `raanDeg`, phases, thrust) |
| `lib/body_clearance_sweep.m`, `lib/repair_dro_elfo_cells.m` | located Moon and Earth clearance; cell repair |
| `tests/run_tests.m` | 8 suites |

Shared pieces written for this pair live in `costate_common/` and serve any
orbit pair: the `'elfo'` family in `get_family_orbit`, `pair_transfer`,
`plot_pair_transfer`, `transfer_difficulty`, `located_min_altitude`.

## Matched comparison (2026-11-19)

Same DRO (period 1 ND), same 6 x 6 phase grid, same cold recipe, same rungs,
fresh solves for all three sheets: tulip (control), ELFO at `raanDeg` 0 and 90.
36 cells and 396 possible entries per sheet. Reproduce the tables with
`compare_dro_targets` once the sheets are in `results/compare/`.

Share of cells accepted, and median transfer time (ND):

| thrust N | tulip | ELFO 0 | ELFO 90 | tf tulip | tf ELFO 0 | tf ELFO 90 |
|---|---|---|---|---|---|---|
| 15 | 1.00 | 0.97 | 0.97 | 0.098 | 0.091 | 0.091 |
| 10 | 0.97 | 0.94 | 0.97 | 0.121 | 0.112 | 0.112 |
| 5 | 0.97 | 0.94 | 0.97 | 0.175 | 0.168 | 0.165 |
| 3 | 0.97 | 0.94 | 0.97 | 0.237 | 0.221 | 0.215 |
| 2 | 0.97 | 0.92 | 0.92 | 0.305 | 0.285 | 0.283 |
| 1.5 | 0.94 | 0.86 | 0.92 | 0.365 | 0.346 | 0.339 |
| 1 | 0.92 | 0.75 | 0.83 | 0.468 | 0.491 | 0.449 |
| 0.75 | 0.89 | 0.67 | 0.72 | 0.595 | 0.595 | 0.611 |
| 0.5 | 0.72 | 0.39 | 0.50 | 0.816 | 0.787 | 0.845 |

| | tulip | ELFO 0 | ELFO 90 |
|---|---|---|---|
| accepted entries (of 396) | 372 | 334 | 350 |
| full ladders 15 to 0.5 N (of 36 cells) | 26 | 13 | 16 |
| lowest rung with at least 90% of cells | 1 N | 2 N | 1.5 N |
| refuted by the conjugate test (0 undecided) | 2 | 3 | 1 |
| below the 500 km Moon floor on the located minimum | 1 | 4 | 6 |
| lowest located Moon altitude (km) | 179 | -181 | -14 |
| accepted, clear of the floor and not refuted | 369 | 327 | 343 |
| full ladders after both checks | 23 | 12 | 15 |
| median revolutions about the Moon at 0.5 N | 0.34 | 0.69 | 0.75 |
| median delta-V at 5 N (km/s) | 2.40 | 2.30 | 2.25 |

Findings:

1. **Above about 1.5 N the two targets are equally easy.** 92 to 100 percent
   of cells solve on every sheet, and the ELFO transfer is 4 to 9 percent
   shorter with slightly less delta-V.
2. **Below 1.5 N the ELFO falls away faster.** At 0.5 N 72 percent of tulip
   cells are accepted against 39 and 50 percent for the ELFO; about half as
   many ELFO cells carry a full ladder.
3. **The ELFO transfers wind more at low thrust** (median 0.7 revolutions
   about the Moon at 0.5 N against 0.34), and transfer time varies more from
   cell to cell.
4. **Plane orientation matters:** `raanDeg` 90 is easier than 0 at every rung
   from 1.5 N down.
5. **Moon clearance is the ELFO's own hazard.** Ten accepted ELFO entries pass
   below the 500 km floor after refinement, two of them through the Moon,
   against one tulip entry. In each case the direct solution rode the floor
   and the refinement step, which has no floor, left it. All ten arrive at
   ELFO phase index 2 or 3 of 6. The chain demotes such entries.
6. **Second-order optimality does not separate the targets:** 6 of 1,056
   entries refuted, spread across the three sheets.
7. **Below 0.5 N neither target gets far with this recipe.** On four
   full-ladder cells per sheet: 0.3 N solved on 2 of 4 for all three; 0.2 N on
   0 (tulip), 1 and 1 (ELFO); nothing at 0.15 N or below. Twelve cells do not
   separate the targets there.

How ladders stopped: solver or acceptance failure, or the engine's
plausibility screen rejecting a direct solution. The screened solutions are 3
to 10 times longer than their neighbours' and repeat to four digits across
unrelated cells, so they are counted as real stops (the solver left the short
family), on both targets alike. The retry pass from a second transfer-time
seed applies only to cells with no accepted rung; one such cell remains (ELFO
0, cell (1,1)).

Limits of the measurement: one DRO size, a 6 x 6 grid, two plane
orientations, one recipe (top-rung seed 0.15 ND, 400 intervals, 500 km Moon
floor, no Earth floor), and cold basin selection is not reproducible to the
cell. The tulip control is a fresh solve, not the shipped catalog.

## Pilot sheet (2026-11-20)

One 12 x 12 sheet at `raanDeg` 0 (144 cells, 1,584 possible entries; 6.6 h on
one process), then the chain on a copy of it (14 min).

| thrust N | 15 | 12 | 10 | 7 | 5 | 3 | 2 | 1.5 | 1 | 0.75 | 0.5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| direct solution found (of 144) | 142 | 141 | 140 | 140 | 140 | 139 | 138 | 127 | 115 | 96 | 67 |
| accepted | 141 | 140 | 139 | 140 | 137 | 137 | 136 | 124 | 112 | 93 | 58 |
| median t_f (ND) | 0.091 | 0.104 | 0.114 | 0.138 | 0.166 | 0.221 | 0.284 | 0.344 | 0.474 | 0.593 | 0.829 |

- 1,357 entries accepted by the solve; 54 cells carry a full ladder. Most of
  the loss at low thrust is inherited: a ladder that stopped higher up never
  attempts the rungs below, so the 0.5 N row is a delivery count for this
  recipe, not a measure of what is solvable there.
- 13 entries in 9 cells are holes above an accepted rung (the ladder goes on
  past a rung whose refinement failed).
- Clearance: 27 accepted entries (2 percent, in 21 cells) pass below the
  500 km Moon floor after refinement, 11 of them through the Moon (lowest
  -1,035 km). 19 of the 27 arrive at ELFO phase 1/12 to 3/12 after perilune.
  No surviving entry lies between 500 and 550 km: entries either clear the
  floor or leave it outright. The chain demotes the 27.
- Repair solved both cells that had no accepted rung (20 entries).
- Packaged: 1,350 entries. Conjugate test 1,346 pass, 4 refuted, 0 undecided.
  Hypothesis gates 1,350 of 1,350. Audit through the picker clean.
- The independent re-flight misses the target by at most 9e-5 km (median
  4e-7 km), so the 100 km acceptance gate does not bind.

An outside review of the recipe (GPT-6.2 Astra, 2026-11-20) accepted it for a
seed library and asked that the comparison above be read as "this recipe
delivered fewer low-thrust ELFO entries", which is what it measures. Its
suggested follow-ups, not done: re-seed the first missing rung of a stopped
ladder from a phase neighbour at that rung; repeat the lowest-clearance cases
at 800 intervals; a half-step node angle on a paired subset.

## Library functions

Functions for reading the packaged catalog (`results/costate_catalog_dro_elfo.mat`,
one sheet per `raanDeg`), the ELFO counterpart of the DRO to tulip user kit
(deliverable 9). They read the catalog and fly entries with pumpkyn's
propagator; nothing is solved or re-certified. Sheets are matched exactly on
(DRO size, `raanDeg`), as by the picker, never to the nearest sheet. All in
`catalog/` unless noted.

| function | what it does |
|---|---|
| `dro_elfo_describe` | fact sheet: orbits, sheets, rungs, Isp, mass, phase grid; entries, full ladders and lowest-rung cells per sheet; t_f and delta-V per rung; conjugate-test counts |
| `dro_elfo_pick` | one entry by (DRO size, `raanDeg`, phases, thrust): the nearest cell holding an entry at that rung |
| `dro_elfo_transfers_at` | every rung recorded at the nearest phase cell of one sheet, high thrust first |
| `dro_elfo_fly` | flies one entry from its z8: history in ND and days / km / kg, arrival miss, located minimum Moon altitude, revolutions; optional 3D plot (rotating frame, Moon-centred) |
| `dro_elfo_extremes` | fastest and slowest (or least and most delta-V) entries under a filter on `raanDeg`, DRO size, thrust and conjugate verdict |
| `dro_elfo_torus_plot` | phase torus of one sheet at one rung (flight time and verdict map), or every sheet side by side on one colour scale (`'all'`) |
| `dro_elfo_movie` | movie (.mp4 and .gif, 1280 x 720) of one entry, or two side by side on one clock, in the pumpkynPie scene |
| `dro_elfo_entry`, `dro_elfo_find_sheet` | helpers: one catalog row with derived units; the exact sheet match |
| `dro_elfo_library_example.m` (top folder) | worked example of all of the above; outputs to `results/library_demo/` |
| `tests/test_dro_elfo_library.m` | test of the functions on the pilot catalog and a 2-sheet copy of it |

By hand:

```matlab
cd orbit_transfer/DRO_ELFO
setup_dro_elfo_paths
C = load('results/costate_catalog_dro_elfo.mat').costate_catalog_dro_elfo;
dro_elfo_describe(C);
e = dro_elfo_pick(C, 1, 90, 0.3, 0.6, 1);          % DRO 1 ND, raanDeg 90, phases, 1 N
F = dro_elfo_fly(C, e, struct('plot', true));       % F.tDays, F.rMoonKm, F.massKg, F.flyKm
T = dro_elfo_transfers_at(C, 1, 90, 0.3, 0.6);      % the thrust ladder at that cell
[eMin, eMax] = dro_elfo_extremes(C, struct('thrustN', 1));
dro_elfo_torus_plot(C, 'all', 1);
dro_elfo_movie(C, {eMin, eMax}, 'results/library_demo/fast_vs_slow');
```

The stored z8 belongs to its cell's phases (`e.sD`, `e.sA`); for other phases
it is a starting guess for `pumpkyn.cr3bp.tfMin`, as section 3 of the example
shows.

## Status

Built: user script, shared functions, front door and driver, comparison,
chain and picker; pilot sheet solved and run through the chain. Running: the
other 11 sheets (`raanDeg` 30 to 330), four processes in parallel. No
deliverable is packaged.
