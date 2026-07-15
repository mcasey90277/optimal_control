# Design: Explicit insertion points across the low-thrust pipelines

**Date:** 2026-07-15
**Status:** approved (brainstorming), pending implementation plan
**Scope:** Make the GTO departure (`rv0`) and the tulip/ELFO **insertion point**
(`rvf`) *explicit and changeable in the high-level drivers* of every low-thrust
pipeline (min-fuel, min-time, seed generation), instead of implicitly threading
them from whatever a seed `.mat` file happened to bake in. Option A: **keep the
current points** (no re-solve); just make them declared, guarded, and labeled.

## Problem / goal

An audit showed **no high-level driver declares its endpoints** — every one reads
`rvf`/`rv0` from a loaded seed (`E.rvf` / `S.rvf`). The endpoint helpers exist but
are bypassed by the production flow (`gto_tulip_endpoints`'s max-ẏ point is used
*only* to draw a movie backdrop). That disconnect is exactly how the "reference"
tulip min-time (max-ẏ, 6.2907 ND) drifted from the point the fronts actually use
(the slow campaign point, dMoon 28k, 5.8267 ND). Fix: declare the endpoints
explicitly in one place, guard against silent drift, and record them in outputs.

## Context / decisions (from brainstorming)

- **Option A — keep the current points.** They are legitimate points on their
  orbits (the tulip campaign point sits ON the tulip: 119 km / 0.1× the trace
  sample spacing). Not a re-solve; a declaration + guard + labeling task.
- **One hardcoded constant only.** Verified: all 20 tulip backbones share an
  identical `rvf` (spread 0.00), so the tulip `'campaign'` point is a single
  well-defined constant. The ELFO seed's `rvf` **exactly** equals
  `gto_elfo_endpoints('nearest', ref=tulip_campaign_rvf)` (diff 0.00), so ELFO is
  fully **reproducible from the helper** — no ELFO hardcode needed, and the
  drift-guard passes with zero re-solve.
- **Scope = all high-level drivers, interactive AND batch**, plus the min-time
  drivers and the seed generators (defense in depth: the guard then fires at
  every level).
- **`rv0` is already consistent** everywhere (matches to 0); the change centers
  on `rvf` (but the helper returns both for a single explicit block).

## The exact hardcoded constant (verified)

```
rv0 (GTO departure, all pipelines):
  [0.00349629072294633, -0.0072962582600817, 0, 4.19147893803368, 8.98865558978329, 0]
rvf tulip 'campaign' (slow south-pole insertion, dMoon 28.3k km, speed 0.31 ND):
  [1.00658107295709, 0.0425745746906059, -0.0557780910480905, -0.16004281347248, 0.0665702939657711, -0.260455693516549]
```

## Design

### 1. Single source of truth — `sundman_minfuel/insertion_states.m`

```
function [rv0, rvf, meta] = insertion_states(target, criterion)
```
- `target` : `'tulip'` | `'elfo'`
- `criterion` :
  - tulip: **`'campaign'`** (default; the hardcoded slow point above) · `'maxydot'`
    (`gto_tulip_endpoints` max-ẏ) · `'apoapsis'` (min-speed point on the tulip)
  - elfo: **`'nearest'`** (default; `gto_elfo_endpoints('nearest', ref=tulip
    campaign rvf)`) · `'apolune'` · `'perilune'` (both via `gto_elfo_endpoints`)
- Returns `rv0` (the shared GTO departure — hardcoded constant above, matching all
  pipelines), `rvf` (per target+criterion), and `meta` = struct with `.target
  .criterion .label` (e.g. `label='tulipCampaign'` / `'elfoNearest'`) for filenames
  and provenance.
- Lives in `sundman_minfuel/` (the shared engine dir, already on every pipeline's
  path). The ELFO branch reaches `gto_elfo_endpoints` by adding `../elfo` to the
  path internally (elfo callers already have it; tulip callers never hit that branch).

### 2. Explicit INSERTION block in every high-level driver

Each driver gets a visible, commented block near its parameters:
```matlab
% ---- INSERTION POINTS (edit here to retarget) --------------------------------
insertion = 'campaign';                       % tulip: 'campaign'|'maxydot'|'apoapsis'
[rv0, rvf, insMeta] = insertion_states('tulip', insertion);
```
(ELFO drivers use `('elfo','nearest')`.) Behaviour by driver role:

- **Consumers that load a seed** (`PSR/run_psr`, `PSR/psr_run_one`,
  `elfo/run_elfo_minfuel`, `elfo/elfo_run_one`, `sundman_minfuel/gen_tulip_mintime`,
  `elfo/gen_elfo_mintime`): declare via the helper, then **assert the loaded
  seed's `rvf` equals the declared `rvf`** (tolerance ~1e-12) — fail loud on drift.
  They still consume the seed; the declaration makes the target explicit + checked.
- **Seed generator that computes the target** (`elfo/gen_elfo_energy_gravhom`):
  replace its `gto_elfo_endpoints('nearest', ref=…)` call with
  `insertion_states('elfo','nearest')`; get the tulip ref via
  `insertion_states('tulip','campaign')`.
- **Seed generators that continue from a backbone** (`PSR/gen_energy_seed`,
  `sundman_minfuel/gen_tulip_energy_2p`, `elfo/gen_elfo_energy_tfsweep`): declare
  via the helper + assert the base backbone's `rvf` matches (same guard).

### 3. Endpoints + label into every saved output

Every result/data struct these drivers save gains `rv0`, `rvf`, and
`insertion` (the `meta.label` string). (Energy seeds already save `rv0`/`rvf`;
add the `insertion` label there too.)

### 4. Filename insertion tag

New outputs carry the label in the filename, e.g.
`minfuel_Tulip_tulipCampaign_tf…`, `mintime_tulip_tulipCampaign.mat`,
`psr_data_…_tulipCampaign_…`. Existing files keep their names (Option A: no
renames — the in-file `insertion` field is the provenance for those); the tag
applies going forward. Since `'campaign'`/`'nearest'` are the only points today,
the tag is future-proofing, but it makes any future retarget self-labeling.

## Files

- **Create:** `sundman_minfuel/insertion_states.m` (+ `test_insertion_states.m`).
- **Modify (explicit block + guard):** `PSR/run_psr.m`, `PSR/psr_run_one.m`,
  `elfo/run_elfo_minfuel.m`, `elfo/elfo_run_one.m`,
  `sundman_minfuel/gen_tulip_mintime.m`, `elfo/gen_elfo_mintime.m`.
- **Modify (declare + use/assert in seed gen):** `PSR/gen_energy_seed.m`,
  `sundman_minfuel/gen_tulip_energy_2p.m`, `elfo/gen_elfo_energy_gravhom.m`,
  `elfo/gen_elfo_energy_tfsweep.m`.
- **Modify (save endpoints+label + filename tag):** the export/save points these
  drivers reach (`elfo_export_data`, `psr_export_data`, the min-time savers, the
  batch summary savers) — add `rv0`/`rvf`/`insertion` to the saved struct and the
  tag to new filenames.

## Out of scope

- Re-solving any front or changing the actual insertion points (that would be
  Option B; a separate effort).
- The low-level solvers (`casadi_*`, `minfuel_at_tf`) — already parameterized by
  `rv0`/`rvf` args; they need no change (they receive the now-explicit values).
- Renaming existing seed/result files.

## Verification

1. **Helper unit test** (`test_insertion_states.m`): `insertion_states('tulip',
   'campaign')` returns exactly the hardcoded constant; `insertion_states('elfo',
   'nearest')` equals `gto_elfo_endpoints('nearest', ref=tulip campaign rvf)`;
   `rv0` matches the constant; `meta.label` strings correct. (No solve.)
2. **Drift-guard passes on existing seeds** (the load-and-assert path): running
   each consumer driver against its current seed does NOT trip the assert
   (verified pre-facts: all backbones + the ELFO seed match the helper to 0.00).
3. **No behavior change**: a min-time or min-fuel driver run produces the same
   `rvf` and result as before (the declaration equals what the seed held).
   (Requires CasADi/R2025b for the solve drivers; the helper + guard checks are
   CasADi-free and can be validated alone.)

## Git

- Focused commits on `main` (or a branch per the finishing flow): helper + test
  first (the guardrail), then drivers, then the save/tag changes.
- Result `.mat`/logs gitignored as usual.
