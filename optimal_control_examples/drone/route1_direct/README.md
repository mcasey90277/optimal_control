# route1_direct

Stage 1 of the drone program: the paper's multi-agent kinematic drone problem
(Zhang, Liu, Darbon & Karniadakis, arXiv:2610.30388, eq. 19) as one
CasADi/IPOPT trapezoidal NLP. Pairwise clearance `|x_i - x_j| >= 2r` is a hard
constraint at the nodes AND at the Hermite midpoints, or, in `barrier` mode,
the paper's own C1 log barrier in the cost (the smooth twin Route 2 shoots).
Five gates check each solve, and the folder reproduces the paper's Table 1.

## Run it

- One case with plots and a movie: edit section 1 of `run_drones.m` (layout,
  M, mode, N, T), then `run_drones`. It solves, runs the gates, prints the
  verdicts against the rigorous pair bound (and the paper's numbers for the
  `paper` layout), and writes a four-snapshot PNG and an MP4 to `results/`.
- The table: `drone_table1` or, detached, `batch/run_table1.sh` (about
  76 minutes; the M = 32 hard row alone is about 45: two ~11-minute
  uniform seeds, two ~2.5-minute alternating ones and an ~11-minute mirror). Writes `results/table1.csv` and the fixtures
  `data/drones_direct_M{4,8}.mat`; the log is `results/table1_batch.log`.
- Tests: `run_tests` from `drone/`.

## Table 1, our numbers

From `results/table1.csv` (batch of 2026-10-23, regenerated after the final
review). `pair bound` is the rigorous lower bound from the two-agent closed
form applied to the antipodal pairs.

| M | N | J hard | seed | gates | J barrier | J_LB | pair bound | OCFormer | SympOCNet | gap % | ours - OCFormer |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 400 | 0.2006416 | u+ | G1-G4 | 0.2015095 | 0.2000 | 0.2003202 | 0.2003 | 0.2001 | 0.32 | +3.4e-4 |
| 8 | 400 | 0.4038874 | a+ | G1-G4 | 0.4083980 | 0.4000 | 0.4006403 | 0.4036 | 0.4046 | 0.97 | +2.9e-4 |
| 16 | 400 | 0.8161853 | u- | G1-G4 | 0.8365288 | 0.8000 | 0.8012807 | 0.8307 | 0.8345 | 2.02 | -1.45e-2 |
| 32 | 100 | 1.2073306 | u+ | G1-G3 | 1.2969196 | 1.1556 | 1.1581178 | 1.2531 | 1.4413 | 4.48 | -4.58e-2 |

The gap is ours (hard) relative to J_LB. `seed` is the multistart seed the
row came from (`u` uniform bump, `a` alternating amplitudes 1.2 / 0.6 on
alternate agents of each ring; `+`/`-` the swirl sense); `gates` the gates
that ran (`gatesRun`), all of which pass. The paper's numbers are reported,
never used as gates.

**What "certified" means here.** The gates prove a GATED LOCAL MINIMUM: a
converged IPOPT solve with defects < 1e-9, clearance on a 10x finer
Hermite grid, J above J_LB, a mesh-doubling cost change < 1e-3 (M <= 16)
and a mirrored-swirl re-solve. That a row is the GLOBAL optimum rests on the
multistart evidence only: four cold seeds per hard row (costs per seed in
the batch log), the cheapest kept, plus any cheaper basin the mesh gate's
2N re-solve finds, which is adopted.

| M | u+ | u- | a+ | a- |
|---|---|---|---|---|
| 4 | 0.2006465 | 0.2006465 | 0.2006465 | 0.2006465 |
| 8 | 0.4044143 | 0.4044143 | 0.4039109 | 0.4039109 |
| 16 | 0.8162405 | 0.8162405 | 0.8162405 | 0.8162405 |
| 32 | 1.2073306 | 1.2073306 | 1.2272301 | 1.2272301 |

(seed costs at N = 100, cold; the rows then refine in N.)

## The finding

The two-agent closed form (`tests/test_direct_oracle.m`:
`ell = 2 sqrt(D^2 - R^2) + R (pi - 2 acos(R/D))`, `J* = ell^2/(4T)`) gives a
rigorous lower bound for every antipodal layout by relaxing it to its antipodal
pairs (the `pairBound` column). For M = 4 that bound is 0.2003202 and our
hard optimum 0.2006416 sits 3.2e-4 above it: that is what the cross-pair
avoidance costs (on the pair-optimal paths adjacent agents are 0.028 apart at
mid-time, less than 2r = 0.04), and every seed tried lands on the same
minimum.

- OCFormer's 0.2003 is a 4-decimal rounding and lies within 3e-5 of the
  bound 0.2003202: the rounding leaves at most 3e-5 for the cross-pair
  avoidance the gated optimum needs 3.2e-4 for.
- SympOCNet's 0.2001 is 2.2e-4 BELOW the rigorous bound. Two readings:
  either its trajectory violates the stated 2r clearance, or the paper's own
  trapezoid cost quadrature on its (unstated) time grid under-evaluates the
  swirl. The second is the likelier, given the paper's claim of positive
  clearance.
- At M = 8 the row is the Z4 "two squares" basin (alternate agents at radii
  0.056 and 0.028 at T/2), 0.4038874, which sits inside the paper's seed band
  0.4039 +/- 0.0001: at M = 8 the paper is consistent with the optimum. The
  regular octagon (0.4043878 at N = 400, the row before the final review) is
  a non-global local minimum 0.12 % higher; OCFormer's 0.4036 is 2.9e-4
  below our row.
- At M = 16 and 32 our hard optimum is 1.7 % and 3.7 % below OCFormer.

## Barrier rows

The barrier row is the hard basin's twin BY CONSTRUCTION: it is solved once,
at the hard row's final N, warm-started from the hard solution (IPOPT warm
start), and its mirror from the hard mirror. Its J INCLUDES the penalty
integral (trapezoid of the barrier at the nodes), measured in the batch log:

| M | J barrier | smooth part | penalty | % of J | min nodal h | branch (l = 1e-4) |
|---|---|---|---|---|---|---|
| 4 | 0.2015095 | 0.2006460 | 8.64e-4 | 0.43 | 8.0e-5 | quadratic |
| 8 | 0.4083980 | 0.4039461 | 4.45e-3 | 1.09 | 1.18e-4 | log |
| 16 | 0.8365288 | 0.8164603 | 2.01e-2 | 2.40 | 8.6e-5 | quadratic |
| 32 | 1.2969196 | 1.2087524 | 8.82e-2 | 6.80 | 6.9e-5 | quadratic |

So the barrier row sits 0.43 % (M = 4) to 7.4 % (M = 32) above the hard row,
and every row but M = 8 has some node on the barrier's quadratic branch
(`h <= l`), where the C1 extension, not `-eps log h`, is active.

## For Route 2

`out.lamDefect` (the defect multipliers, read from `opti.lam_g` by row range)
column k samples `+V` at the MIDPOINT of interval k: the adjacent-column
average `(mu_{k-1} + mu_k)/2` equals `V` at node k to 2.8e-17 on both barrier
fixtures (in hard mode to ~1e-7, except next to an active midpoint clearance,
whose multiplier enters V-stationarity: 1.8e-4 on the M = 8 hard fixture).
So the PMP costate is `lambda = -lamDefect`, stationed at the interval
midpoints, and the matching oclib map is `oc.duals_to_costates` with scheme
`'trapezoid-nodal'` (the adjacent-interval average at the N+1 nodes), NOT
`'trapezoid'` (left-node stationing, which injects a half-step phase error).

## Caveats

- The mesh gate's 2N re-solve is warm (IPOPT `warm_start_init_point`,
  `mu_init 1e-4`), and a cheaper 2N solution is ADOPTED as the row. Before
  the final review a cold restart at 2N repeatedly found the cheaper
  symmetry-broken Z4 basin at M = 8, and the loop read the cost drop as
  discretisation error and re-solved cold, landing back in the octagon. Now
  every hard row M <= 16 adopts its 2N solutions up to the N = 400 cap
  (the drops are 8.0e-6 to 5.4e-5 relative: ordinary refinement; the M = 8
  basin change now comes from the `a+` seed at N = 100); at the cap the
  2N = 800 change is 2.9e-6 to 3.6e-6 relative.
- M = 32 has no mesh gate (slow, per plan; `gatesRun` G1-G3). Its mirror
  re-solve lands 7.8e-5 above the row (`mirrorRel`, reported, not gated),
  although the multistart's own cold `u-` seed reached the row's cost
  1.2073306: on a ~240-iteration path the basin a cold solve picks is not
  reproducible to this level.
- G5 is a swirl-sense flip: it tests seed independence and the reflection
  symmetry, and cannot see a symmetry-broken basin.
- The cold barrier guess has coincident agents on ring layouts
  (`solve_drones_direct:guess` refuses it), another reason barrier rows are
  warm-started.
- M = 64 and 128 are stretch rows for the batch launcher and were not run.

## Files

| file | what |
|---|---|
| `solve_drones_direct.m` | one solve: modes `hard` / `barrier`, 2D (v_z pinned) or 3D with a ball obstacle, clearance at nodes and Hermite midpoints, straight-line-plus-swirl initial guess, defect multipliers by row range |
| `verify_drones_direct.m` | gates: G1 fine-grid clearance, G2 bound, G3 solver + defects, G4 mesh (2N), G5 mirror swirl; returns `.ok/.why` |
| `antipodal_pair_bound.m` | the rigorous lower bound from the two-agent closed form on a layout's antipodal pairs |
| `run_drones.m` | front door: section 1 settings, solve, gates with verdicts, panels + movie to `results/` |
| `drone_table1.m` | the Table 1 loop (4-seed multistart, warm mesh loop with adoption) -> `results/table1.csv`, fixtures `data/drones_direct_M{4,8}.mat` |
| `plot_drones.m`, `movie_drones.m` | the paper's four-snapshot figure; 1280x720 MP4 |
| `casadi_path.m` | puts CasADi on the path |
| `batch/run_table1.sh`, `batch/job_table1.m` | detached batch run with a log file (`JOB DONE` marker); the only practical way to run M = 64 / 128 |
| `tests/test_direct_oracle.m` | closed-form head-on case (G6), N = 100 and 400 |
| `tests/test_direct_small.m` | M = 4 circle: all gates, barrier pinned, failure path, sphere |
| `tests/test_table1_fixture.m` | committed fixtures, the CSV, the mesh-adoption decision |
| `tests/test_plot_drones.m` | headless plot and movie |
