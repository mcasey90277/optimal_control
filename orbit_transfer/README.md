# orbit_transfer — optimal orbit-transfer campaigns

Container for all orbit-transfer work in this repo. Organization: **one folder
per transfer problem**, each split into `direct/` (collocation NLP → IPOPT) and
`indirect/` (PMP shooting) codebases, plus shared libraries, tutorials, and
reference material.

## Problem campaigns

| folder | problem | direct | indirect |
|---|---|---|---|
| `earth_elliptic_to_geo/` | elliptic (GTO-like) → GEO, Earth 2-body (HMG-2004 benchmark, 1500 kg, 10→0.1 N ladder) | **working, certified** — MEE/L-domain campaign, full Table-3 thrust ladder, `run_gergaud` front door | MfMax (Gergaud-group Fortran) built + validated as cross-check; our own MATLAB indirect is future work |
| `earth_elliptic_to_geo_CR3BP/` | same transfer **with lunar gravity** (Earth–Moon CR3BP) | not started — README/TODO scope the plan | not started |
| `GTO_tulip/` | GTO → south-pole tulip orbit, Earth–Moon CR3BP (15 kg, 25 mN) | **working, certified** — Sundman min-fuel engine, 25-switch flagship, ΔV–t_f front | built (ms_band, ifs, ztl, min_time) but not yet certified — the active goal |
| `GTO_ELFO/` | GTO → ELFO lunar frozen orbit, Earth–Moon CR3BP | **working** — front mapped, min-time anchor certified | not started (Route C open) |

Each campaign folder keeps `README.md` + `TODO.md` at its root, campaign
records in `process/`, technical notes in `doc/`.

## Shared library

| folder | what |
|---|---|
| `cr3bp_common/` | Single source for the CR3BP GTO problem definition — `cr3bp_lt_params`, `minfuel_config`, `gto_tulip_endpoints`, `gto_elfo_endpoints` — plus `setup_cr3bp_common()` (adds pumpkyn). Every GTO_tulip/GTO_ELFO module's `setup_paths` calls it. |

## Tutorials (guided build-it-yourself, with `mytry/` + verified checkpoints)

| folder | what |
|---|---|
| `min_energy_tutorial/` | Min-energy point-to-point transfer: indirect shooting + collocation + primer-vector verification (exercises PDF + reference solvers). |
| `lambert/` | Universal-variables Lambert solver incl. multi-revolution, validated vs pyKep. |

## Reference material

| folder | what |
|---|---|
| `min_fuel_papers/` | Source papers (HMG-2004 preprint, Caillau–Noailles 2001, Zhang 2015, MfMax manuals, …). |
| `min_fuel_paper/` | Our min-fuel paper outline (co-author Koblick). |

## Conventions

- MATLAB R2025b only; CasADi 3.7.0 at `~/casadi-3.7.0`; run modules from their
  own folder after calling that module's `setup_paths`.
- `.mat` results are gitignored campaign caches; committed figures live in each
  campaign's `results/` (under `direct/`).
- Cross-references between campaigns are deliberate and documented in each
  `setup_paths.m` header (e.g. GTO_ELFO reuses GTO_tulip's Sundman engine).
