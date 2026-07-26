# GTO_ELFO — execution paths

Survey 2026-07-26, read-only. Tool: `docs/tools/callgraph.py`.

**Scale:** 24 `.m` + 3 shell scripts (+2 in attic). The tidiest campaign — and
the only one of the four with **zero orphan candidates**.

---

## 1. What it borrows, and what it forked

ELFO defines 24 files and borrows 14. But the split is informative:

| source | functions borrowed | top users |
|---|---|---|
| `cr3bp_common` | `cr3bp_lt_params` (18), `minfuel_config` (15), `cr3bp_fingerprint` (8), `thrust_tag` (5), `gto_elfo_endpoints` (2), `gto_tulip_endpoints`, `setup_cr3bp_common`, `check_cr3bp_fp` | everywhere |
| `verify_common` | `foc_check`, `foc_report`, `foc_manifest`, `foc_ipopt_inertia`, `setup_verify_common` | `run_foc_elfo` |
| **`GTO_tulip/direct/sundman_minfuel`** | **`insertion_states` — and nothing else** | 6 files |

`cr3bp_common` is doing its job well. The tulip dependency is the finding.

### The stale claim

`setup_paths.m`'s own header says this campaign *"reuses
`casadi_minfuel_sundman` / `insertion_states` / `minfuel_at_tf`, retargeted to
ELFO."* Measured: it **calls only `insertion_states`**. The other two appear
only in comments and provenance notes. ELFO did not retarget tulip's solver —
it **forked** it, into `casadi_energy_freetf` and `casadi_mintime_freetf`.

So `setup_paths` puts all 34 files of tulip's direct engine on the path to
reach one function, and the comment explaining why is no longer true.

### How big is the fork?

| | code lines (comments stripped) |
|---|---|
| `casadi_minfuel_sundman` (tulip) | 161 |
| `casadi_energy_freetf` (ELFO) | 188 |
| **identical lines in common** | **122 — 76% of the smaller** |

Longest shared run: **41 consecutive identical lines** (the IPOPT options
block). Then 16 (state box constraints), 12 (the output struct), 9 (warm-start
priming), 6 (the CasADi path bootstrap).

## 2. Execution paths

```
run_elfo_minfuel   [the min-fuel driver]
  elfo_find_energy_seed ─► picks a banked energy seed by t_f and fingerprint
  gen_elfo_minfuel      ─► casadi_energy_freetf   (eps 1 -> 0 sharpen)
  elfo_export_data      ─► data products
  run_foc_elfo          ─► foc_check / foc_report (verify_common)

gen_elfo_mintime   [the certified Route-B anchor]
  casadi_mintime_freetf   (s == 1 hard all-burn, t_f free via cScale)

gen_elfo_energy_gravhom  [the seed factory: 4-leg gravity-homotopy ladder]
  casadi_energy_freetf

Batch/ladder (shell, per-rung process isolation):
  elfo_energy_sweep.sh  ->  elfo_batch.sh 0 energy  ->  elfo_collect_summary
  elfo_movies.sh all    ->  elfo_render_movies -> elfo_movie   (no re-solve)
```

Third campaign with a working shell-script ladder pattern
(`run_cr3bp_ladder.sh`, `reproduce_table3.sh`, and these three) — the
requirement for a "clean script and shell script that reproduces thrust
ladders" already has three precedents to copy.

## 3. Findings worth acting on

1. **Narrow the tulip path dependency.** ELFO adds 34 files to the path for one
   function. `insertion_states` is endpoint/state machinery used by *both*
   campaigns, and `cr3bp_common` already holds `gto_tulip_endpoints` and
   `gto_elfo_endpoints`. **Move `insertion_states` into `cr3bp_common` and drop
   the tulip path entirely.** Small, well-justified, and it removes a shadowing
   surface: today an ELFO session that also touched PSR would have two
   definitions of `casadi_minfuel_sundman` on the path.
2. **Fix the stale `setup_paths` header** — it claims a reuse that does not
   happen and will mislead the next reader into thinking the solver is shared.
3. **Do NOT merge the two solvers.** The 24% that differs is the substantive
   part — the 9th state (`cScale`), free final time, the two-primary clock,
   the gravity homotopy. Merging would produce a heavily conditional solver and
   put two campaigns' certified results at risk for a cosmetic win.
4. **Do extract the incidental duplication.** The 41-line IPOPT options block
   is identical across both solvers and has no reason to be duplicated; the
   CasADi path bootstrap (6 lines) is in every solver in the repo. These are
   safe regardless of the fork decision, and gated trivially by byte-identity
   of results.
5. **18 of 23 non-test files unreached by any test** (~78%).
6. **Zero orphans** — the only campaign with none.

## 4. Library candidates

- **`insertion_states` → `cr3bp_common`** (finding 1). Highest value/risk ratio
  in this campaign: it removes a whole path dependency.
- **IPOPT options block** — 41 identical lines across `casadi_energy_freetf`,
  `casadi_mintime_freetf` and `casadi_minfuel_sundman`; almost certainly also
  in `casadi_lt_mee`. A single `ipopt_opts(profile)` helper.
- **CasADi path bootstrap** — 6 lines repeated in every solver.
- **`cr3bp_common` needs no work** — with 6 functions and up to 18 callers it
  is the proof that the shared-library approach works here when the shared
  thing is a *problem definition* rather than an algorithm.
