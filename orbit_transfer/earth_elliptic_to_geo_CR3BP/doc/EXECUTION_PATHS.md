# earth_elliptic_to_geo_CR3BP — execution paths

Survey 2026-07-26, read-only. Tool: `docs/tools/callgraph.py` (read its METHOD
AND ITS LIMITS section before trusting any single line).

**Scale:** 18 `.m` (15 non-test, 3 test) — the smallest campaign by an order of
magnitude, because it borrows almost everything.

---

## 1. The headline: this campaign is the *good* code-sharing pattern

It defines 18 files and **borrows 18 distinct functions** from
`earth_elliptic_to_geo` and `verify_common`. It duplicates **nothing** across
the campaign boundary. `setup_paths` delegates to the 2-body campaign's own
`setup_paths`, and the lunar physics enters the *shared* solver through an
opt-in `par.pert` branch in `lt_mee_rhs` rather than through a forked copy.

**The repo has now run this experiment twice, with opposite results, and that
is the most useful thing this survey has to say about the library question:**

| | strategy | outcome |
|---|---|---|
| **CR3BP** (here) | extend — path delegation + opt-in branches in the shared file | zero duplication; upstream fixes arrive automatically. The `opti.dual` fix reached this campaign for free, because it shares `casadi_lt_mee` verbatim |
| **PSR/lib** | vendor — frozen snapshot for self-containment | isolation works, but upstream fixes never arrive. Silently broke `run_psr` stage 5c, undetected for a day |

For the planned library work: **prefer opt-in extension of a shared file over
vendoring, and over inventing an abstraction.** The winning pattern here was
not "extract a common core" — it was "add a switch to the existing one."

### What it borrows

| from earth | used by | | from `verify_common` | used by |
|---|---|---|---|---|
| `kepler_lt_params` | 11 | | `foc_check` | 1 |
| `table3_certified` | 6 | | `foc_report` | 1 |
| `optdef` | 4 | | `foc_manifest` | 1 |
| `mee_res_to_cart_res` | 4 | | `foc_ipopt_inertia` | 1 |
| `casadi_lt_mee` | 3 | | `setup_verify_common` | 1 |
| `transfer_movie` | 3 | | | |
| `mee_seed`, `homotopy_mee`, `table3_recipes`, `verify_pmp_mee`, `switch_structure`, `hamiltonian_const_check`, `lt_mee_rhs` | 1–2 each | | | |

## 2. What it adds

| file | role |
|---|---|
| `lunar_params.m` | the new physics (μ_M, D_M, n_M, φ₀, gain) — hub, 5 callers |
| `run_cr3bp_geo.m` | **the front door**: parameter block, then stages A–E inline |
| `bridge_mu_continuation.m` | μ-continuation gain walk, 0 → 1, as a standalone stage |
| `solve_cr3bp_minfuel.m` | ε-sharpen stage; calls the bridge |
| `verify_cr3bp_pmp.m` | lunar-aware verification driver (physical + generic FOC) |
| `refresh_duals_cr3bp.m` | dual repair for banked artifacts |
| `compare_vs_2body.m` | the Moon-effect comparison harness |
| `sanity_bound.m` | tide/authority null model |
| `viz/` (5) | ladder, phase-sweep, quadrupole movie renderers |

## 3. Execution path

```
run_cr3bp_geo   [front door — the parameter block is section 1, "edit this section only"]
  stage A  two-pass seed .......... mee_seed            (earth)
  stage B  2-body energy, eps=1 ... casadi_lt_mee       (earth)
  stage C  mu-continuation 0->1 ... casadi_lt_mee       (earth, gain walk INLINE)
  stage D  eps-homotopy -> epsMin . homotopy_mee        (earth)
  stage E  data products .......... transfer_movie, mee_res_to_cart_res (earth)
     └─ verify_cr3bp_pmp
          ├─ refresh_duals_cr3bp ──► casadi_lt_mee      (earth)
          ├─ verify_pmp_mee                             (earth)
          ├─ foc_check / foc_report / foc_ipopt_inertia (verify_common)
          └─ hamiltonian_const_check, switch_structure  (earth)

solve_cr3bp_minfuel   [separate stage entry point]
  └─ bridge_mu_continuation ──► casadi_lt_mee, mee_seed (earth)
```

Ladder reproduction already exists: `run_cr3bp_ladder.sh` (per-rung process
isolation via the front door) — this is the pattern the tulip/ELFO campaigns
should copy for the ladder-script requirement.

## 4. Findings worth acting on

1. **The pipeline is implemented twice.** `run_cr3bp_geo` inlines stages A–D,
   while `bridge_mu_continuation` + `solve_cr3bp_minfuel` implement the same
   bridge-then-sharpen sequence as callable stages. The gain-walk logic appears
   in both (26 vs 49 lines mentioning `gsched`/`gain`). The README calls the
   latter "pipeline stages", but the front door does not call them. Decide
   which is authoritative and make the other delegate to it.
2. **The two-pass seed protocol exists in three places** —
   `run_transfer_mee` (earth), `bridge_mu_continuation`, and `run_cr3bp_geo`
   stage A. `bridge_mu_continuation`'s own header says it is *"mirrored
   VERBATIM from run_transfer_mee.m lines 132-161"*, which is an explicit
   copy-paste admission and therefore the **cleanest library-extraction target
   in the repo so far**: one named function, three call sites, a documented
   provenance trail, and no cross-campaign politics.
3. **13 of 15 non-test files are unreached by any test** (~87%) — the highest
   ratio of the campaigns surveyed. Mitigated by the front door being
   gate-guarded, but worth knowing.
4. **`setup_paths` uses a cwd trick** (both campaigns' files are named
   `setup_paths.m`; it `cd`s into the 2-body folder so cwd precedence resolves
   the other one). Documented as review amendment C, and it works — but it is
   the kind of cleverness that will confuse the next reader, and it makes the
   two campaigns' path setup non-composable in the general case.
5. Three viz/analysis scripts have no code callers but **are** referenced from
   docs and the TODO (`fig_ladder_dmf`, `fig_phi_sweep`,
   `phase_control_sensitivity`) — manual figure generators, not dead code.
   Third campaign in a row where folder-local orphan detection over-reported.

## 5. Library candidates (recorded, not acted on)

- **Two-pass seed protocol** — 3 call sites, explicit verbatim-copy comment.
  Highest-confidence extraction in the repo. Gate: byte-identical results on
  the earth 10 N row and the CR3BP 10 N row.
- **μ-continuation / gain-walk stage** — resolve the double implementation
  first (finding 1), then it is a single function.
- Everything already borrowed (§1) needs no work — it is the model.
