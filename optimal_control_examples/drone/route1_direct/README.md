# route1_direct

Stage 1 of the drone program: the paper's multi-agent kinematic drone problem
(Zhang, Liu, Darbon & Karniadakis, arXiv:2610.30388, eq. 19) as one
CasADi/IPOPT trapezoidal NLP. Pairwise clearance `|x_i - x_j| >= 2r` is a hard
constraint at the nodes AND at the Hermite midpoints, or, in `barrier` mode,
the paper's own C1 log barrier in the cost (the smooth twin Route 2 shoots).
Five gates certify each solve, and the folder reproduces the paper's Table 1.

## Run it

- One case with plots and a movie: edit section 1 of `run_drones.m` (layout,
  M, mode, N, T), then `run_drones`. It solves, runs the gates, prints the
  verdicts against the rigorous pair bound (and the paper's numbers for the
  `paper` layout), and writes a four-snapshot PNG and an MP4 to `results/`.
- The table: `drone_table1` or, detached, `batch/run_table1.sh` (tens of
  minutes; the M = 32 hard solve alone is about 11 minutes). Writes
  `results/table1.csv` and the fixtures `data/drones_direct_M{4,8}.mat`.
- Tests: `run_tests` from `drone/`.

## Table 1, our certified numbers

From `results/table1.csv` (batch of 2026-10-23). `pair bound` is the rigorous
lower bound from the two-agent closed form applied to the antipodal pairs.

| M | N | J hard | J barrier | J_LB | pair bound | OCFormer | SympOCNet | gap % | ours - OCFormer |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 100 | 0.200646 | 0.201514 | 0.2000 | 0.200320 | 0.2003 | 0.2001 | 0.32 | +3.5e-4 |
| 8 | 400 | 0.404388 | 0.408398 | 0.4000 | 0.400640 | 0.4036 | 0.4046 | 1.10 | +7.9e-4 |
| 16 | 100 | 0.816240 | 0.836561 | 0.8000 | 0.801281 | 0.8307 | 0.8345 | 2.03 | -1.4e-2 |
| 32 | 100 | 1.207331 | 1.296920 | 1.1556 | 1.158118 | 1.2531 | 1.4413 | 4.48 | -4.6e-2 |

The gap is ours (hard) relative to J_LB. The paper's numbers are reported,
never used as gates.

## The finding

The two-agent closed form (`tests/test_direct_oracle.m`:
`ell = 2 sqrt(D^2 - R^2) + R (pi - 2 acos(R/D))`, `J* = ell^2/(4T)`) gives a
rigorous lower bound for every antipodal layout by relaxing it to its antipodal
pairs (the `pairBound` column). For M = 4 that bound is 0.200320, so the
paper's SympOCNet cost 0.2001 cannot satisfy the stated clearance, and
OCFormer's 0.2003 sits at the bound with no room for the cross-pair avoidance
that adjacent agents need (on the pair-optimal paths they are 0.028 apart at
mid-time). Our hard optimum 0.200646 is mesh-converged (N = 100 and 400 agree
to 5e-6) and swirl-independent.

At M = 16 and 32 our hard optimum is 1.7 % and 3.7 % below OCFormer. At M = 4
and 8 it is slightly above (+3.5e-4, +7.9e-4), consistent with the bound.

## Caveats

- M = 8 needed N = 400 to pass the mesh gate: changes of 1.3e-3 between
  doublings did not quarter, i.e. the doubled-mesh re-solve hops between
  near-equal swirl basins.
- M = 32 has no mesh gate (slow, per plan); its two swirl senses converge to
  local minima 7.8e-5 apart (`mirrorRel`; the lower one is the row).
- Barrier rows are warm-started from the hard solution, because the cold guess
  has coincident agents on ring layouts (`solve_drones_direct:guess` refuses it).
- Barrier at M = 32 is 7.4 % above hard.
- M = 64 and 128 are stretch rows for the batch launcher and were not run.

## Files

| file | what |
|---|---|
| `solve_drones_direct.m` | one solve: modes `hard` / `barrier`, 2D (v_z pinned) or 3D with a ball obstacle, clearance at nodes and Hermite midpoints, straight-line-plus-swirl initial guess, defect multipliers by row range |
| `verify_drones_direct.m` | gates: G1 fine-grid clearance, G2 bound, G3 solver + defects, G4 mesh (2N), G5 mirror swirl; returns `.ok/.why` |
| `antipodal_pair_bound.m` | the rigorous lower bound from the two-agent closed form on a layout's antipodal pairs |
| `run_drones.m` | front door: section 1 settings, solve, gates with verdicts, panels + movie to `results/` |
| `drone_table1.m` | the Table 1 loop -> `results/table1.csv`, fixtures `data/drones_direct_M{4,8}.mat` |
| `plot_drones.m`, `movie_drones.m` | the paper's four-snapshot figure; 1280x720 MP4 |
| `casadi_path.m` | puts CasADi on the path |
| `batch/run_table1.sh`, `batch/job_table1.m` | detached batch run with a log file (`JOB DONE` marker); the only practical way to run M = 64 / 128 |
| `tests/test_direct_oracle.m` | closed-form head-on case (G6), N = 100 and 400 |
| `tests/test_direct_small.m` | M = 4 circle: all gates, barrier pinned, failure path, sphere |
| `tests/test_table1_fixture.m` | committed fixtures and the CSV |
| `tests/test_plot_drones.m` | headless plot and movie |
