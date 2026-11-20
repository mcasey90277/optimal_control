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

## Status

Built: user script, shared functions, front door and driver, comparison,
chain and picker. Not yet run: the 12 x 12 pilot sheet through the chain, and
the 12-sheet fleet (`raanDeg` every 30 deg). No deliverable is packaged.
