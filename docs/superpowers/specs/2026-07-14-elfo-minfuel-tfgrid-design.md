# Design: ELFO min-fuel tf-grid campaign

**Date:** 2026-07-14
**Status:** approved (brainstorming), pending implementation plan
**Scope:** Map the GTO→ELFO minimum-fuel ΔV–time front by sweeping the
energy→fuel homotopy across a grid of transfer times, mirroring the tulip PSR
batch infrastructure. Produces the ELFO min-fuel **tf-grid convergence map**
(which factors reach ε=0 bang-bang vs stall, switch counts, ΔV per factor).

## Problem / goal

The GTO→ELFO min-fuel-direct cell is 🟡: ε=0 solved at one factor (1.20×,
34-switch). The deliverable is a **band**, not a point — and the energy-seed
band is *wider* than the ε=0-convergent (fuel) band, so which factors actually
sharpen to bang-bang must be *mapped*, not assumed ([[minfuel-tf-grid-strategy]]).
The tulip target already has this via a batch trio (`psr_run_one` +
`psr_batch.sh` + `psr_collect_summary`); ELFO needs the parallel infrastructure.

## Context / decisions (from brainstorming)

- **Drive by `factor = tf / tfMin`**, full PSR parity. `tfMin` = the **tulip**
  min-time 6.2906939607 ND, reused as a shared reference scale — ELFO's own
  min-time is unsolved (both ELFO min-time cells are ⬜), so `factor` here is a
  labeling scale, not "×ELFO-min-time." `run_elfo_minfuel` already labels it
  this way (`t_f/tfMin_tulip`). **ELFO min-time is explicitly deferred** to a
  separate next objective (would fill 2 goal-matrix cells, give a real
  `tfMin_ELFO`, and anchor the 0-switch low-tf end of the front — but it is a
  known-hard open problem, not a prerequisite for this campaign).
- **Factor-key the energy seeds** (`energy_elfo_f<NNNN>.mat`, NNNN =
  round(1000·factor)), mirroring PSR's `energy_f<NNNN>.mat`. This reconciles the
  current ND-tf naming (`energy_elfo_tf<round(1000·tf)>.mat`) with PSR. No
  per-grid seeds exist on disk yet, so nothing to migrate; the base seed
  `energy_elfo_freetf.mat` keeps its special name as the lookup fallback.
- **`elfo_run_one` is a separate callable function**, mirroring how
  `psr_run_one.m` sits beside the interactive `run_psr.m`. The interactive
  `run_elfo_minfuel.m` script stays as-is; `elfo_run_one` calls the same
  `gen_elfo_minfuel` core, so duplication is minimal (seed-lookup + result-row).
- **Self-contained in `elfo/`**: result rows + summary land in `elfo/results/`
  (PSR uses a sibling `PSR_data/`; ELFO keeps everything under `elfo/`).
- **Process isolation is non-negotiable** — the uncatchable CasADi/IPOPT MEX
  FATAL crash (~1 in 10 solves) means each factor must run in its own
  `matlab -batch` process, exactly as `psr_batch.sh` does.

## Design

Two phases.

### Phase 1 — energy seed band (modify + run existing code)

**Modify `elfo/gen_elfo_energy_tfsweep.m`:**
- Re-parameterize opts from ND-tf to **factor**: `factorLo` / `factorHi` /
  `factorStep` (converting to `tfTarget = factor·tfMin` internally). Default band
  ≈ **1.11–2.00×, step ~0.08** (the current 7.0–12.5 ND / 0.5-ND tuning
  re-expressed); tunable after we see where energy stops converging. Keep the
  existing loose-continuation + step-halving mechanism unchanged.
- **Save seeds as `energy_elfo_f<round(1000·factor)>.mat`** (was
  `energy_elfo_tf<round(1000·tf)>.mat`); grid summary struct carries `factor`
  alongside `tf`.

**Modify `elfo/run_elfo_minfuel.m`:** update its per-factor seed lookup from
`energy_elfo_tf<round(1000·tf)>.mat` to `energy_elfo_f<round(1000·factor)>.mat`
(base-seed `energy_elfo_freetf.mat` fallback unchanged).

**Modify `elfo/gen_elfo_minfuel.m`:** docstring only — its two mentions of the
old `energy_elfo_tf####.mat` seed name → `energy_elfo_f####.mat`. (It takes the
seed via `opts.seedFile`, so no code path changes.)

**Run** `gen_elfo_energy_tfsweep` → banks `energy_elfo_f*.mat` per factor +
`energy_elfo_tfgrid.mat` (the energy band `[factorLo, factorHi]`).

### Phase 2 — fuel batch trio (build in `elfo/`)

**`elfo/elfo_run_one.m`** — `row = elfo_run_one(factor, opts)`:
- Resolve the seed `energy_elfo_f<round(1000·factor)>.mat` (error if missing).
- Call `gen_elfo_minfuel(struct('seedFile',seed,'target','ELFO',
  'epsMin',opts.epsMin, ...))` → the ε:1→0 homotopy solution file.
- Load the solution; compute the **result row**: `factor`, `tf` (ND), `tf_days`,
  `ok`, `epsReached` (did the homotopy hit `epsMin`?), `epsFloor` (the smallest
  ε converged if it stalled short), `dV`, `prop`, `switches`, `edge`, `defect`,
  `ipoptStatus`, `dataFile`, `err`.
- Save `row` to `elfo/results/elfo_result_f<NNNN>_minEps<eTag>.mat`. Resumable:
  skip the solve if that row already exists (unless `opts.rerun`).
- opts defaults mirror `psr_run_one` where applicable: `epsMin` [0],
  `movieMode` ['none'], `maxIter`, `resDir` [elfo/results].

**`elfo/elfo_batch.sh`** — near-clone of `psr_batch.sh`:
- Usage `elfo_batch.sh <epsMin> <factor1> [factor2 …]` or
  `elfo_batch.sh <epsMin> energy` (auto-discover: glob
  `elfo/results/energy_elfo_f*.mat`, parse the factor from each name).
- One `matlab -batch "cd('elfo'); elfo_run_one(<factor>, <opts>)"` **process per
  factor**; per-factor watchdog (`WATCHDOG_S`, default 1800 s → kill + continue);
  continue-on-crash (a MEX FATAL kills only that factor).
- MATLAB path `/Applications/MATLAB_R2025b.app/bin/matlab` (R2025a is broken —
  [[use-matlab-2025b]]); env overrides `MATLAB_BIN`, `MOVIE`, `WATCHDOG_S`.
- Logs to `elfo/results/logs/elfo_batch_<timestamp>.log`. Resumable (re-run the
  same command; finished factors skip instantly).
- After the loop, call `elfo_collect_summary`.

**`elfo/elfo_collect_summary.m`** — `res = elfo_collect_summary(epsMin)`:
- Scan `elfo/results/elfo_result_f*_minEps<eTag>.mat`, collect the `row`
  structs, sort by factor, print the tf-grid table (factor, tf/days, ok,
  ε-reached, switches, ΔV, prop, defect), save
  `elfo/results/elfo_batch_summary_minEps<eTag>.mat`. Clone of
  `psr_collect_summary.m`.

### Data flow

```
gen_elfo_energy_tfsweep  →  energy_elfo_f*.mat   (per-factor energy seeds)
        →  elfo_batch.sh  →  elfo_run_one(factor)/process
        →  gen_elfo_minfuel (ε:1→0)  →  minfuel_ELFO_tf*_minEps*.mat (solutions)
        →  elfo_result_f*_minEps*.mat (rows)
        →  elfo_collect_summary  →  the ELFO min-fuel tf-grid convergence map
```

The **ε-reached** column is the ELFO-specific payoff: the energy band is wider
than the fuel-convergent band, so the map records which factors sharpen all the
way to bang-bang (ε=0) and which stall (and at what ε-floor / switch count).

## Out of scope (untouched)

- **ELFO min-time** — the next distinct objective; not part of this campaign.
- `PSR/` and its trio (the template, unchanged).
- The `casadi_energy_freetf` / `gen_elfo_minfuel` solve cores (called, not
  modified — except the one docstring seed-name fix in `gen_elfo_minfuel`).

## Verification

1. **Phase 1:** run `gen_elfo_energy_tfsweep` on a short band (e.g. 1.15–1.30×);
   confirm it banks `energy_elfo_f<NNNN>.mat` files and an `energy_elfo_tfgrid.mat`
   with a sane band, each seed converged (maxDefect < 1e-6). Confirm
   `run_elfo_minfuel` finds a factor-keyed seed via the updated lookup.
2. **Phase 2 unit:** `elfo_run_one(1.20, struct('epsMin',0))` produces a solution
   + an `elfo_result_f1200_minEps0.mat` row with the expected fields; a second
   call with the row present skips the solve (resumability).
3. **Phase 2 batch:** `elfo_batch.sh 0 1.15 1.20` runs both factors in separate
   processes, survives a killed process (watchdog), and
   `elfo_collect_summary(0)` prints/saves the 2-row map. (Requires CasADi /
   R2025b; if unavailable in the run environment, say so — do not fake a pass.)

## Git

- One focused branch of commits on `ifs-retarget` (Phase-1 edits, then each
  Phase-2 file), main = latest per the established flow.
- Result `.mat` and logs are gitignored (`*.mat`, `results/logs/`), like PSR's.
