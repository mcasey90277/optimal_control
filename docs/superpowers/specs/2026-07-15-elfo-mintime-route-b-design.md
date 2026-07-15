# Design: GTO→ELFO min-time anchor (Route B — hard all-burn)

**Date:** 2026-07-15
**Status:** approved (brainstorming), pending implementation plan
**Scope:** Solve the GTO→ELFO **minimum-time** transfer to land a sharp
`tfMin_ELFO` anchor — the number that (a) fills the ELFO min-time goal-matrix
cell "by any means" and (b) relabels the mapped ELFO min-fuel ΔV–time front from
the borrowed `tf/tfMin_tulip` scale into ELFO's own `tf/tfMin_ELFO` units and
gives its 0-switch low-`tf` endpoint.

## Problem / goal

Both ELFO min-time cells are ⬜. The whole mapped ELFO min-fuel front is
currently plotted against `factor = tf/tfMin_tulip` (6.2907 ND) — a *borrowed*
scale, because ELFO's own min-time was never solved. The deliverable is one
trustworthy `tfMin_ELFO` (ND + days) with a machine-tight, verified always-burn
trajectory to the same terminal the front rendezvouses to.

Success criterion (set in brainstorming): **anchor by any means / fastest
trustworthy number** — we do not need both matrix cells, just a defensible
`tfMin_ELFO`.

## Context / decisions (from brainstorming)

- **Route A (energy→time continuation) was run first and walls.** Stepping the
  pinned `tf` down from the converged energy seed rides the throttle up cleanly
  to `tf = 6.2278 ND (27.61 d)` at edge **57.5%** (all rungs defect ~1.6e-15),
  then the next step (~6.04 ND) diverges (`inf_pr` frozen ~0.067, `inf_du` →
  1e10) — the same near-min-time conditioning wall tulip has below 1.12×. The
  smooth energy problem walls *before* it degenerates to all-burn (edge 57% ≠
  100%), so Route A yields only a **loose upper bound** `tfMin_ELFO ∈ (~5.5,
  6.23) ND`. (Calibration: tulip's energy continuation floored at 1.12× but its
  true min-time was 1.0× — the energy floor overestimated `tfMin` by ~12%; the
  same ~10–12% gap predicts `tfMin_ELFO ≈ 5.5 ND / 24–25 d`.) Route A's lowest
  converged rung `results/energy_elfo_f0990.mat` (6.2278 ND, edge 57.5%) is the
  **Route B warm start**.
- **Formulation = hard all-burn (`s≡1`), minimize `t(τ_f)`.** True min-time
  ansatz; tulip's min-time was 0-switch all-burn, so ELFO's should be too. If
  `s≡1` cannot reach `rvf` at any `t_f`, that *proves* ELFO min-time has coast
  arcs — reported as the finding, fall back to a throttle-free min-time (the
  pre-agreed escape hatch). No continuation: the objective **is** min-time, so
  IPOPT drives `t_f` down until the rendezvous BC can no longer be met all-burn.
- **Solver is a sibling, not a mode.** New `casadi_mintime_freetf.m` beside
  `casadi_energy_freetf.m`, following the campaign's sibling convention
  (`casadi_minfuel_sundman` → `casadi_energy_freetf` → this). Isolation +
  unit-testability outweigh the ~15-line dynamics duplication; the existing
  file's documented "min-time mode" (`tfTarget=[]` + `tfWeight>0`) is actually
  fuel+time, so it is *not* reused.
- **Anti-dive guard is diagnostic-only for now.** The GPT-flagged "dive at a
  primary to slow the clock" cheat is expected suppressed by `s≡1` +
  rendezvous BC. Gate on a post-solve diagnostic (`min(r1)` ≥ GTO perigee, `t(τ)`
  monotone, `cScale` in box); **only if it trips** add the hard path constraint
  `r1 ≥ r_perigee` and re-solve. (User decision: don't add the constraint
  pre-emptively.)
- **Same terminal as the front.** `rv0`, `rvf`, `τ_f`, and the two-primary clock
  params (`pSund`, `qSund`, `moonZone`) all come from the energy seed
  (`energy_elfo_f0990.mat`) so the anchor is consistent with the mapped front.
- **`factor` stays a labeling scale until the anchor exists**, then the front is
  re-plotted as `tf/tfMin_ELFO`. Out of scope to re-run the front here.

## Design

### (A) Solver — `elfo/casadi_mintime_freetf.m`

A sibling of `casadi_energy_freetf` with the identical two-primary Sundman clock,
`cScale` free-`t_f` slack (Betts banded free-time), and boundary conditions.
Exactly two departures:

- **Throttle hard-pinned `s≡1`.** Drop `s` from the control; control `u = α` (3),
  `‖α‖ = 1` the only control freedom. Thrust accel = `(Tmax/m)·α`,
  `mdot = −Tmax/cEx` (constant, max depletion). State stays
  `x = [r(3); v(3); m; t; cScale]` (9).
- **Objective = `X(8, end)`** (minimize physical `t(τ_f)`). No `∫s·dt`, no `ε`,
  no `tfTarget`. `cScale` floats (box `[0.10, 8]`), pinned from below by
  feasibility (must traverse GTO→ELFO within the fixed `τ_f` budget).

Unchanged from the energy sibling: `τ_f` fixed (= `tauf0` from the seed, the
~40-rev Sundman length); trapezoidal defects in σ (`dX/dσ = τ_f·cScale·κ·f`);
BCs `X(1:6,1)=rv0`, `X(7,1)=1`, `X(8,1)=0`, `X(1:6,end)=rvf`; two-primary
`κ = (r1^−q + (r2/D)^−q)^(−p/q)`; full gravity `muGain=1`; explicit two-sided
box bounds; unit-steering equality `sum(α²)−1 = 0`.

**Signature (mirrors the sibling):**
`out = casadi_mintime_freetf(sigma, rv0, rvf, Tmax, cEx, muStar, X0, U0, tauf0, opts)`
where `X0` is `[8×nN]` or `[9×nN]` (cScale row appended if absent), `U0` is
`[3×nN]` steering (or `[4×nN]` from an energy seed — throttle row dropped), and
`opts` carries `pSund, qSund, moonZone, cBox, tfCapMult, maxIter, warmTight`.

**Outputs:** `out.X [9×nN] .U [3×nN] .tf (=tfMin) .cScale .mf .maxDefect .maxUnit
.minR1 .tMonotone .success .ipoptStatus` plus KKT-dual costates `.lamDef` (for a
possible later PMP diagnostic — free to save, not required for the anchor).

### (B) Driver + verify — `elfo/gen_elfo_mintime.m`

Mirrors `gen_elfo_energy_gravhom` / `gen_elfo_minfuel`:

1. Load `results/energy_elfo_f0990.mat`; take `sigma, rv0, rvf, tauf0, pSund,
   qSund, moonZone, X` (states), and **override the throttle to `s≡1`** for the
   warm start (drop `U(4,:)`, keep steering `U(1:3,:)`; seed `t_f` from the
   seed's `X(8,end) = 6.2278`).
2. Call `casadi_mintime_freetf` → machine-tight all-burn min-time.
3. Save `results/mintime_elfo.mat` (`X, U, tf=tfMin, cScale, mf, maxDefect,
   minR1`, plus `rv0, rvf, tauf0, pSund, qSund, moonZone` for downstream reuse).
4. **Independent verify** (solver-free, the `verify_elfo_seed` pattern):
   re-integrate `dX/dτ` from `X(:,1)` under the same dynamics with `s≡1`,
   recompute the max defect and the rendezvous residual `‖X(1:6,end) − rvf‖`
   from scratch (no CasADi); confirm defect < 1e-8 and residual < 1e-8, and
   `mf = 1 − (Tmax/cEx)·tfMin` consistent.
5. Print `tfMin_ELFO` (ND + days) and the **relabel factor** `tf/tfMin_ELFO` for
   each existing front point.

### (C) Success criteria + guards

Accept the anchor only if ALL hold:
- IPOPT `Solve_Succeeded`; `out.maxDefect < 1e-8` (target 1e-10+);
  rendezvous residual < 1e-8; unit-steering `maxUnit < 1e-8`.
- **`s≡1` feasible** — an all-burn trajectory that hits `rvf`. If infeasible at
  every `t_f` (restoration failure / cannot close rendezvous): STOP, report
  "ELFO min-time is not all-burn," fall back to throttle-free min-time (separate
  effort, not built here).
- **Physicality / anti-dive:** `minR1 ≥ r_perigee_GTO` (no unphysical Earth
  dive), `t(τ)` monotone (automatic for `cScale>0`), `cScale` strictly interior
  to its box. If any trips → add `r1 ≥ r_perigee` path constraint, re-solve.
- **Sanity band:** `tfMin_ELFO < 6.23 ND` (strictly below the energy floor, as
  physics requires) and near the tulip-analogy estimate (~5.5 ND / 24–25 d). A
  wildly-off value flags a wrong `τ_f`/mesh, not a valid anchor.

**Cross-check (optional, cheap):** `pumpkyn.cr3bp.tfMin` (analytic-STM single
shooting) for a fully independent number. It floored ~1e-3 on 13-rev cases, so
treat agreement as confirmation and non-convergence as "single-shooting couldn't
close," NOT as invalidating the direct result.

### Data flow

```
energy_elfo_f0990.mat (6.2278 ND, edge 57.5%, s overridden to 1)
   → gen_elfo_mintime → casadi_mintime_freetf (min t(τ_f), s≡1)
   → mintime_elfo.mat  (tfMin_ELFO + all-burn trajectory)
   → independent verify (solver-free defect + rendezvous recompute)
   → report tfMin_ELFO (ND/days) + relabeled front factors
```

## Out of scope (untouched)

- Re-running / re-plotting the ELFO min-fuel front in the new units (follow-up).
- The indirect ELFO min-time cell (Route C); `casadi_energy_freetf`,
  `gen_elfo_minfuel`, PSR, and the front batch tooling.
- The throttle-free min-time fallback solver (built only if `s≡1` proves
  infeasible).

## Verification

1. **Solver unit** (`test_mintime_freetf.m`): construct on the seed; confirm the
   NLP builds with the `[9×nN]` state / `[3×nN]` control layout, `s≡1` wired
   (thrust term uses full Tmax), and one solve returns `Solve_Succeeded` with
   defect < 1e-8. A tiny finite-difference or CasADi-jacobian sanity on the
   dynamics block vs `casadi_energy_freetf` at `s=1` (they must match).
2. **Driver** (`gen_elfo_mintime`): produces `mintime_elfo.mat`; the independent
   solver-free verify recomputes defect < 1e-8 and rendezvous < 1e-8; guards
   evaluated and reported; `tfMin_ELFO` printed in the sanity band.
3. **Anchor sanity:** `tfMin_ELFO < 6.23 ND`, `mf` consistent with
   `1 − (Tmax/cEx)·tfMin`, `minR1 ≥ r_perigee`. (Requires CasADi / R2025b; if the
   run environment lacks them, say so — do not fake a pass.)

## Git

- Focused commits on `ifs-retarget` (solver, then driver, then the test), per the
  established flow.
- Result `.mat` and logs gitignored (`*.mat`, `results/logs/`), like the rest of
  `elfo/`.
