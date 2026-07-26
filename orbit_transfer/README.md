# orbit_transfer — optimal orbit-transfer campaigns

Container for all orbit-transfer work in this repo. Organization: **one folder
per transfer problem**, each split into `direct/` (collocation NLP → IPOPT) and
`indirect/` (PMP shooting) codebases, plus shared libraries, tutorials, and
reference material.

## Problem campaigns

| folder | problem | direct | indirect |
|---|---|---|---|
| `earth_elliptic_to_geo/` | elliptic (GTO-like) → GEO, Earth 2-body (HMG-2004 benchmark, 1500 kg, 10→0.1 N ladder) | **working, certified** — MEE/L-domain campaign, full Table-3 thrust ladder, `run_gergaud` front door | MfMax (Gergaud-group Fortran) built + validated as cross-check; our own MATLAB indirect is future work |
| `earth_elliptic_to_geo_CR3BP/` | same transfer **with lunar gravity** (Earth–Moon CR3BP) | **working, certified** — full ladder 10→0.1 N, `run_cr3bp_geo` front door, reviewed technical note; Moon effect mapped vs thrust (~50→30 g) and vs lunar phase (sign-flipping quadrupole) | Phase 2, not started |
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
| `abstracts/` | Conference abstract drafts (e.g. cislunar-conference GTO→tulip transfer talk). |

## Conventions

- MATLAB R2025b only; CasADi 3.7.0 at `~/casadi-3.7.0`; run modules from their
  own folder after calling that module's `setup_paths`.
- `.mat` results are gitignored campaign caches; committed figures live in each
  campaign's `results/` (under `direct/`).
- Cross-references between campaigns are deliberate and documented in each
  `setup_paths.m` header (e.g. GTO_ELFO reuses GTO_tulip's Sundman engine).
- **Standard optimality report:** every production solve driver ends with
  `foc_report` — the fixed-format first-order block + `foc_<tag>.mat`
  sidecar; report-only burn-in, does not alter certified status. Core:
  `verify_common/`.
- **KKT multipliers: never `opti.dual()`.** It returns duals in CasADi's
  *canonicalized* constraint orientation, not the orientation of the `opti.g`
  row they pair with in `grad_f + A'*lam`, so for a row-wise-canonicalizing
  constraint (every collocation defect here) they come back sign-corrupted
  entry by entry. Record the group's `opti.g` row range at build time and
  index `opti.lam_g`. This cost two campaigns nine days of misdiagnosis in
  July 2026 — full post-mortem, the reusable
  "is-it-the-duals-or-my-derivation" diagnostic, and the weak-minimum
  warm-restart rules in
  **`earth_elliptic_to_geo/process/LESSONS_DUAL_EXTRACTION.md`**.

## Cross-campaign registers and lessons

Read before repeating an old mistake:

| file | what |
|---|---|
| `CODE_STRUCTURE.md` | **survey synthesis + extraction plan.** Consolidates the four campaign `doc/EXECUTION_PATHS.md` surveys: the measured comparison of the three code-sharing strategies this repo has tried (extend / vendor / fork), the standing rules derived from it, and a risk-ranked, gated extraction plan. Read before any refactor or library work |
| `OPTIMALITY_CERTIFICATION.md` | **live register** for goal-1 certification, both orders. **Part A (first order)** — instrument list + per-campaign **coverage matrix** (which gate runs on which transfer) + the five standing gaps. **Part B (second order)** — every experiment, the three blocking mechanisms, what is ruled out, open leads, and the decision on what to build next. Consult it *before* planning any optimality work: first-order coverage holes are invisible without the matrix (until 2026-07-25 ELFO had no gate at all — the FOC layer now covers all four campaigns, unevenly), and two campaigns have already built second-order instruments a third was about to re-plan |
| `earth_elliptic_to_geo/process/LESSONS_DUAL_EXTRACTION.md` | the `opti.dual` trap; how to tell a dual bug from a physics bug; why a physics-shaped correlation isn't evidence of physics; weak-minimum warm-restart guards |
| `earth_elliptic_to_geo/process/DEEP_THRUST_LESSONS.md` | making deep low-thrust rungs converge (maxIter, `liftDL`, why `scaleNLP` hurts) |
| `GTO_tulip/process/LADDER_PREP_PILOT_FINDINGS.md` | fixed-tau_f topology wall; why ELFO can ladder and the tulip cannot |
