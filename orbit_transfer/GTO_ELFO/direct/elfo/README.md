# elfo/ — GTO -> ELFO direct min-fuel pipeline

Self-contained deliverable for the minimum-fuel low-thrust GTO -> ELFO transfer
in the Earth-Moon CR3BP. Sibling of `PSR/` (the GTO -> tulip deliverable).

## Model: shared-path, not vendored

Unlike PSR (which vendors a frozen machinery snapshot into `PSR/lib`), this
directory references the two shared engine files -- `cr3bp_lt_params` and
`minfuel_config` -- from `../sundman_minfuel` on the path (`setup_paths.m`).
Single source of truth, no drift surface. Tradeoff: elfo/ tracks the dev
library, so a dev edit to those two files can change ELFO results.

## Pipeline

1. `gen_elfo_energy_gravhom.m` -- build the min-ENERGY seed via the two-primary
   gravity-homotopy ladder on the free-tf solver `casadi_energy_freetf.m`.
   (`gen_elfo_energy_tfsweep.m` builds the tf-grid of energy seeds.)
2. `gen_elfo_minfuel.m` -- sharpen energy -> fuel (epsilon 1 -> 0) to bang-bang.
3. `run_elfo_minfuel.m` -- end-to-end entry: solve -> export -> verify -> movie.
   - `elfo_export_data.m` -- costates from the two-primary KKT duals.
   - `verify_elfo_seed.m` -- solver-free seed verification.
   - `elfo_movie.m` -- control movie (copy of PSR's generalized psr_movie).
4. `gto_elfo_endpoints.m` / `probe_elfo_target.m` -- ELFO endpoints + geometry.

## Smoke tests

- `smoke_energy_freetf.m` -- free-tf form reproduces the f1.20 backbone.
- `smoke_fixedtf.m` -- pinned-tf leg-0 conversion is well-posed (no drift).

## Data

Results (`.mat`, movies) land in `elfo/results/`. The seed-data reads and the
`../PSR_data` reference are shared stores kept in place. Dead-end early routes
(`gen_elfo_energy_backbone`, `gen_elfo_energy_tangential`) are in `attic/`.

Full design rationale: `docs/superpowers/specs/2026-07-14-elfo-sibling-dir-design.md`.
