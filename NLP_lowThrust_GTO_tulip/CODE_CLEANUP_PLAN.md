# Code Cleanup & Reorganization Plan — NLP_lowThrust_GTO_tulip

**Date:** 2026-07-09. **Status:** PLAN — not yet executed.
**Companion docs:** `HONEST_EVALUATION_DV_TF_FRONT.md` (why: §5 code-organization
risk), `LOW_THRUST_MINFUEL_CAMPAIGN.md` (history).

## Diagnosis — how it got disorganized

1. **Three code generations flat in one folder.** The parent dir holds ~35
   loose `.m` files spanning: fmincon-era NLPs (`NLP_lowThrust_GTO_Tulip*.m`,
   `solve_*_nlp*.m`, `nlp_constraints*`, `lt_dynamics*`), Sundman prototypes
   (`run_sundman_*.m`, `test_sundman_*`), tf-continuation experiments
   (`tf_continuation_*`), plus logs, `.mat` snapshots, and `fort.6` junk.
2. **Stale duplicate core solver.** Parent-root `casadi_minfuel_sundman.m` is
   the OLD API (no `warmTight`; 7.9 KB vs the library's 11.3 KB). It already
   caused one shadowing bug ("Too many input arguments"); it will cause more.
3. **`sundman_minfuel/` mixes everything.** Library functions, 4 overlapping
   front drivers (`run_tf_sweep`, `run_tf_front`, `run_tf_2anchor`, `tf_step`)
   with *different* homotopy schedules, 15 `energy_*.mat`, result/plot files,
   and docs all in one flat dir.
4. **No result provenance.** `.mat` files carry no metadata (date, settings,
   seed lineage, git hash); `%.2f` filenames collide on finer grids; the
   canonical `tfMin` is recovered by the `Sm.tf/1.15` hack in 3 places.
5. **Orchestration is ephemeral.** Watchdog/retry loops live in scratchpad
   zsh scripts; a zsh `local`-expansion bug clobbered 13 result files
   (2026-07-09). Nothing checked in, nothing tested.
6. **Repo junk:** `fort.6` (IPOPT scratch), loose campaign logs, `.DS_Store`,
   LaTeX aux files, undecided untracked folders (`lieFiltering/`,
   `gauss_sum_curvature/`).

**Not a problem:** `ms_band/` (lower-band multiple-shooting attack, other
terminal) is well-organized — own campaign doc + unit tests. Model to follow.

## Constraints (hard)

- **Other terminal is live** in `ms_band/`, with `../sundman_minfuel` and
  `../../lowThrust_GTO_tulip` on its path, reading `cr3bp_lt_params`,
  `gto_tulip_endpoints`, and the dual `.mat`s. → NO renames/moves inside
  `sundman_minfuel/` or `lowThrust_GTO_tulip/`, and no touching `ms_band/`,
  until it is idle.
- **MATLAB hold** (user-ordered): no solver runs. File moves and new code are
  fine; smoke-testing waits.
- **Append-only artifacts:** certified `.mat`s are never modified, only moved
  with path fixes (note: `movie/animate_sundman_minfuel.m` loads the PARENT
  copy of `sundman_minfuel_certified.mat`; parent/lib copies verified
  identical, so dedupe must update that path in the same commit).

## Target layout

```
NLP_lowThrust_GTO_tulip/
├── README.md                      # rewritten: map + entry points
├── LOW_THRUST_MINFUEL_CAMPAIGN.md # history (stays)
├── HONEST_EVALUATION_DV_TF_FRONT.md
├── CODE_CLEANUP_PLAN.md           # this file
├── sundman_minfuel/               # CANONICAL LIBRARY (direct method)
│   ├── minfuel_config.m           # NEW: tfMin (explicit!), pSund, N, scheds
│   ├── minfuel_at_tf.m            # NEW: THE per-t_f driver (see below)
│   ├── aggregate_front.m          # NEW: combine+verify+plot, 3 marker classes
│   ├── casadi_minfuel_sundman.m   # core solver (unchanged)
│   ├── sundman_seed_map.m  sundman_homotopy.m  cr3bp_lt_params.m
│   ├── gto_tulip_endpoints.m  setup_paths.m
│   ├── verify_tf_front.m  certify_minfuel_pmp.m
│   ├── build_energy_backbone.m  energy_step.m  direct_build_minfuel.m
│   ├── run_certified_minfuel.m    # reproduces THE certified 1.15x result
│   ├── orchestrate/               # NEW: checked-in zsh (watchdog+retry)
│   │   ├── backbone_walk.sh  sharpen_batch.sh
│   ├── results/                   # NEW files land here (energy/ minfuel/
│   │   fronts/ plots/ logs/); existing .mats migrate in Phase 1
│   ├── test_minfuel_lib.m         # NEW: cheap non-solve checks
│   └── attic/                     # run_tf_sweep|front|2anchor.m, tf_step.m,
│       solve_tf_minfuel.m (superseded by minfuel_at_tf) + their .mats
├── attic/                         # parent legacy: fmincon-era solvers,
│   Sundman prototypes, tf_continuation_*, old logs/mats, OLD core-solver copy
├── movie/                         # unchanged (fix certified-.mat path)
└── ms_band/                       # other terminal's — untouched
```

## Code changes (the substance, beyond moving files)

1. **`minfuel_config.m`** — single source of truth. `tfMin = 6.2906939607`
   stored EXPLICITLY (kills `Sm.tf/1.15` in 3 places), `pSund`, canonical
   homotopy schedules (backbone / sharpen / re-clean), paths, IPOPT presets.
2. **`minfuel_at_tf(factor, opts)`** — consolidates `solve_tf_minfuel` +
   `tf_step` + the `step_tf`/`continue_anchor` copies. Options: seed source
   (`'energy'` backbone | `'neighbor'` bang-bang | file), schedule, branch
   tag, maxIter. Every output `.mat` carries a `meta` struct: date, git hash,
   solver settings, seed provenance, ipopt status, schedule actually used.
   Filenames `%.3f` (or milli-factor `f1200`) — no collisions.
3. **`lamDef` robustness** — `casadi_minfuel_sundman` records the defect
   constraint index range in `out` at construction instead of the implicit
   `lamAll(1:8N)` assumption.
4. **`aggregate_front.m`** — promotes scratchpad `combine_front.m` into the
   repo; produces the HONEST plot with 3 marker classes: feasible upper
   bound / direct-certified extremal / direct+indirect certified; envelope
   drawn only through certified points.
5. **`orchestrate/*.sh`** — the validated process-isolation + watchdog +
   one-retry pattern, checked in, factor list as arguments, `local`
   assignments on separate lines (the 2026-07-09 bug, documented in-file),
   logs to `results/logs/`.
6. **`test_minfuel_lib.m`** — runs WITHOUT solving: config consistency,
   seed-map roundtrip, filename encode/parse roundtrip, lamDef index range
   vs constraint construction order, schedule monotonicity. Cheap guardrail
   for every future refactor.

## Execution phases

**Phase 0 — safe NOW (additions + junk only; no renames in live folders):**
- Add `minfuel_config.m`, `minfuel_at_tf.m`, `aggregate_front.m`,
  `orchestrate/`, `test_minfuel_lib.m`, `results/` (empty, for new files).
  All NEW files — nothing the other terminal loads changes.
- Root junk: delete `fort.6` (×2), `.DS_Store`; clean stray LaTeX aux in
  `mpc_cart_pole/`; add `.gitignore` entries (`fort.6`, `.DS_Store`,
  `matlab_crash_dump.*`).
- Create parent `attic/` and move the PARENT-ROOT legacy files into it
  (fmincon-era solvers, prototypes, old logs/mats, and the stale
  `casadi_minfuel_sundman.m` — the shadowing hazard dies today). Safe: the
  parent root is NOT on the other terminal's path.
- Rewrite parent `README.md` as the folder map.

**Phase 1 — when the other terminal is idle (still no MATLAB needed):**
- Attic the superseded `sundman_minfuel` drivers; migrate `energy_*.mat`,
  `ms_*.mat`, front `.mat`s/plots into `results/` subdirs.
- Dedupe parent `sundman_minfuel_certified.mat` / `minfuel_from_energy_seed.mat`
  (identical copies verified) and fix the `movie/` load path in the same
  commit.

**Phase 2 — when MATLAB is unfrozen (validation):**
- `test_minfuel_lib` green; `minfuel_at_tf(1.20)` reproduces the banked
  backbone sharpen (dV/switches match logged 2026-07-09 values); then the
  physics campaign resumes on clean rails (recover 1.30/1.35/1.40× greens,
  few-switch down-trace from 1.85×, indirect certifier build).

## Decisions needed from Mike

1. `lieFiltering/` and `gauss_sum_curvature/` are untracked in the repo —
   commit them, or .gitignore?
2. Raw campaign logs (`sundman_n50.log`, `tf_continuation*.log`, …): attic or
   delete? (The MD docs carry all the learnings.)
3. Rename existing `energy_*.mat` to the new scheme in Phase 1, or leave
   legacy names and use the new scheme only for new files?
