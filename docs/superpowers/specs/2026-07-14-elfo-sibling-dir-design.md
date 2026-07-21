# Design: `elfo/` sibling deliverable directory

**Date:** 2026-07-14
**Status:** approved (brainstorming), pending implementation plan
**Scope:** Reorganize the GTO→ELFO direct min-fuel code into a self-contained
sibling directory that mirrors PSR's role for the tulip target.

## Problem

The GTO→ELFO direct-method code is scattered across three locations:

- **Working pipeline** — in `sundman_minfuel/`: `casadi_energy_freetf`,
  `gen_elfo_energy_gravhom`, `gen_elfo_energy_tfsweep`, `gen_elfo_minfuel`,
  `run_elfo_minfuel`, `elfo_export_data`, `verify_elfo_seed`,
  `smoke_energy_freetf`.
- **Endpoints + early/failed routes** — in `PSR/`: `gto_elfo_endpoints`,
  `probe_elfo_target`, `gen_elfo_energy_backbone`, `gen_elfo_energy_tangential`.
- **Min-time attempts** — in `min_time/` (out of scope; different method).

There is no single home for the ELFO deliverable, and the working pipeline
reaches into `PSR/` (for `psr_movie`) creating a cross-directory dependency.

## Context / constraints (established during brainstorming)

- **PSR is a self-contained, *vendored* tulip deliverable**, not "the
  direct-method directory." `PSR/setup_paths` adds only `PSR/` + `PSR/lib/`
  (19 frozen machinery copies snapshotted 2026-07-12) + pumpkyn. It was
  deliberately decoupled from the `sundman_minfuel/` dev library so ongoing dev
  edits can't break the deliverable.
- The name "PSR" = *PMP-Steered Refinement*, which is tulip-specific — ELFO
  transfers have **no** PMP-refinement step. So ELFO does not belong *inside*
  PSR under that name.
- Vendoring already carries a drift cost: `PSR/lib/minfuel_at_tf.m` is a Jul-12
  copy now stale vs the `sundman_minfuel/` master (certified-flag commit
  `a2364b2`), and `sundman_seed_map.m` was never vendored.
- **Dependency map of the ELFO pipeline** (verified 2026-07-14):
  - `casadi_energy_freetf` has **zero** shared-engine deps (takes everything as
    args) — moves cleanly.
  - The pipeline's only shared-engine touchpoints are `cr3bp_lt_params` +
    `minfuel_config` (both in `sundman_minfuel/`, already vendored into
    `PSR/lib`), plus `psr_movie`.
  - No non-ELFO file in `PSR/` references `gto_elfo_endpoints` or
    `probe_elfo_target` — safe to move them out of `PSR/`.

## Decisions (approved)

1. **Sibling `elfo/` dir**, not folding ELFO into PSR. One self-contained
   deliverable per target; the directory name matches the target rather than an
   absent refine step. PSR stays tulip-only.
2. **Shared-path machinery model, no re-vendoring.** `elfo/` references
   `cr3bp_lt_params` + `minfuel_config` from `sundman_minfuel/` on the path
   (single source of truth). No new drift surface. Accepted tradeoff: `elfo/`
   is coupled to the dev library, unlike PSR's frozen vendoring — the two
   deliverables use different disciplines, by design.
3. **Dedicated `elfo/elfo_movie.m`** (the ELFO-rendering path lifted out of
   `psr_movie`), so `elfo/` has no dependency on `PSR/`. Accepted tradeoff: a
   movie renderer exists in two forms; movie code is display-only and stable,
   so the mild duplication is acceptable in exchange for a clean split. (The
   rejected alternative was promoting the renderer to a shared
   `sundman_minfuel/lt_movie.m` and editing PSR to use it.)

## Target structure

New directory `GTO_tulip/elfo/`.

Files relocated with `git mv` (preserve history):

| File | From |
|---|---|
| `casadi_energy_freetf.m` | `sundman_minfuel/` |
| `gen_elfo_energy_gravhom.m` | `sundman_minfuel/` |
| `gen_elfo_energy_tfsweep.m` | `sundman_minfuel/` |
| `gen_elfo_minfuel.m` | `sundman_minfuel/` |
| `run_elfo_minfuel.m` | `sundman_minfuel/` |
| `elfo_export_data.m` | `sundman_minfuel/` |
| `verify_elfo_seed.m` | `sundman_minfuel/` |
| `smoke_energy_freetf.m` | `sundman_minfuel/` |
| `smoke_fixedtf.m` | `sundman_minfuel/` (free-tf solver smoke; sole external caller of `casadi_energy_freetf`) |
| `gto_elfo_endpoints.m` | `PSR/` |
| `probe_elfo_target.m` | `PSR/` |

New files created in `elfo/`:

- `setup_paths.m` — adds `elfo/` + `sundman_minfuel/` (for `cr3bp_lt_params` +
  `minfuel_config`) + pumpkyn (`~/Desktop/proj7/external/pumpkyn/src`).
- `elfo_movie.m` — the ELFO-rendering path from `psr_movie`, standalone.
- `README.md` — pipeline overview + the shared-path model + provenance.
- `results/` — ELFO outputs (moved from `sundman_minfuel/results/`).
- `attic/` — the two dead-end routes.

Moved to `elfo/attic/`:

- `gen_elfo_energy_backbone.m` (from `PSR/`)
- `gen_elfo_energy_tangential.m` (from `PSR/`)

Moved to `elfo/results/` (from `sundman_minfuel/results/`):

- `energy_elfo_freetf.mat`, `minfuel_elfo.mat`,
  `movie_ELFO_tf1p200_minEps0.gif`, `movie_ELFO_tf1p200_minEps0.mp4`

## Edits required inside moved files

`run_elfo_minfuel.m` is tied to `here = sundman_minfuel/` in three places that
must retarget once `here` becomes `elfo/`:

1. Line ~35 `cd(here); setup_paths(); addpath(fullfile(here,'..','PSR'));`
   → call `elfo/setup_paths` (which adds `sundman_minfuel/`); **drop** the
   `addpath(../PSR)` (no longer needed — `elfo_movie` is local).
2. Line ~35–36 `resDir = fullfile(here,'results')` → resolves to `elfo/results`
   automatically once the file lives in `elfo/` (seed reads
   `energy_elfo_freetf.mat` / `energy_elfo_tf%04d.mat` then come from
   `elfo/results`). Confirm the moved seed `.mat` files land there.
3. The `psr_movie(...)` call → `elfo_movie(...)`.

`dataDir = fullfile(here,'..','PSR_data')` (line ~110) stays as an in-place
data reference to the shared `PSR_data` store (same pattern PSR uses for
`sundman_minfuel/results`).

Any internal cross-calls among the moved files (e.g. `gen_elfo_minfuel` →
`casadi_energy_freetf`, `gen_elfo_energy_gravhom` → `gto_elfo_endpoints`)
resolve automatically since all callers and callees move together into `elfo/`
and onto the same path.

## Out of scope (untouched)

- `min_time/` ELFO code (indirect min-time — a separate method/exploration).
- `PSR/` itself and `PSR/lib/` vendoring (no de-vendoring, no renaming).
- The `cr3bp_lt_params` / `minfuel_config` masters (referenced, not copied).

## Verification

1. **Path resolution (clean):** from a fresh MATLAB path, run
   `elfo/setup_paths`, then `which -all` each pipeline function
   (`casadi_energy_freetf`, `gen_elfo_energy_gravhom`, `gen_elfo_energy_tfsweep`,
   `gen_elfo_minfuel`, `run_elfo_minfuel`, `elfo_export_data`,
   `verify_elfo_seed`, `smoke_energy_freetf`, `gto_elfo_endpoints`,
   `probe_elfo_target`, `elfo_movie`, `cr3bp_lt_params`, `minfuel_config`) and
   confirm each resolves to `elfo/` or `sundman_minfuel/` — no dangling refs,
   no `PSR/` resolution.
2. **Live smoke:** run `smoke_energy_freetf` headless (matlab-headless skill,
   requires CasADi at `~/casadi-3.7.0`) and confirm it reproduces its expected
   result (the free-tf form reproducing the f1.20 backbone at machine
   precision).
3. **Grep guard:** confirm no remaining reference to the moved files' old
   locations anywhere in the tree (`grep -rn "psr_movie" elfo/` returns
   nothing; no `sundman_minfuel/casadi_energy_freetf` style stale references).

## Git

- Use `git mv` for every relocation to preserve history.
- Single focused commit on the `ifs-retarget` branch (the reorg), message
  describing the sibling-dir + shared-path model.
- Update project docs (`CLAUDE.md` directory map, and the ELFO memory
  `elfo-retarget-open.md` file:line references) to point at `elfo/` in a
  follow-up if needed.
