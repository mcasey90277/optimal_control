# earth_elliptic_to_geo — execution paths

Survey of what calls what, 2026-07-26. Produced by
`docs/tools/callgraph.py direct run_gergaud run_ladder reproduce_row`
(re-runnable; read that file's METHOD AND ITS LIMITS before trusting any single
line of it).

**Scale:** 136 `.m` under `direct/` — 96 non-test, 40 test.

**Headline finding: this campaign is already well-layered.** It is the most
mature of the four and should be the reference pattern the others are measured
against, not the one that gets restructured. What it needs is pruning at the
edges and a couple of misplaced files moved, not a reorganisation.

---

## 1. The layering

```
  frontdoor/  reproduce/          entry points (what a user runs)
        |
  drivers/                        solve modes: min-time anchor, fuel transfer
        |
  psr/                            switch-time refinement (optional stage)
        |
  core/                           transcription + dynamics + seed + homotopy
        |
  coords/  lib/                   pure utilities (frames, elements, optdef)

  verify/  verify/sosc/           optimality checks (1st + 2nd order)
  viz/                            rendering; consumes results, called by nothing
  cartesian_legacy/               retained only for the cross-formulation gate
  attic/                          off-path archive (see attic/README.md)
```

No layer reaches upward. `viz/` and `verify/` hang off the entry points rather
than sitting in the solve path — which is why a report failure can never fail a
solve.

## 2. The solve spine

Every path below bottoms out in the same four-function spine. If you understand
these, you understand the campaign:

| function | role |
|---|---|
| `mee_seed` | build the initial guess (calls `elements_to_cart`, `rtn_frame`, `lt_mee_rhs`) |
| `casadi_lt_mee` | **the transcription** — builds and solves the NLP; the single place the discrete problem is defined |
| `lt_mee_rhs` | the L-domain Gauss dynamics `dX/dL`, MX-safe; the only physics file |
| `homotopy_mee` | the ε: 1→0 energy→fuel sweep, repeatedly calling `casadi_lt_mee` |
| `interp_warmstart` | re-map a solution onto a new mesh for the next solve |

`casadi_lt_mee` is called from 13 places and is shared verbatim with the CR3BP
campaign (same file, `par.pert` switched on) — the single most load-bearing
file in the repo.

## 3. The three entry points

They differ by *use*, not by physics — all three drive the same spine.

**`run_gergaud`** (frontdoor) — one row, user-facing. Optionally anchor
(`run_mintime_mee`) → fuel (`run_transfer_mee`) → refine (`psr_mee_refine`) →
verify (`run_foc_mee`) → render (`transfer_movie`, `gergaud_plot`). This is the
template the other campaigns' entry scripts should follow.

**`run_ladder`** (drivers) — the thrust ladder. Loops
`run_mintime_mee` + `run_transfer_mee` per rung, warm-chaining downward.

**`reproduce_row`** (reproduce) — from-scratch reproduction with per-rung
recipes (`table3_recipes`), a keep-best-mass multi-start, and a one-sided
acceptance check (`verify_row` vs `table3_certified`).

```
run_gergaud
  ├─ run_mintime_mee ──┐
  ├─ run_transfer_mee ─┼─► casadi_lt_mee ──► lt_mee_rhs
  │    ├─ mee_seed ────┘        ▲
  │    ├─ homotopy_mee ─────────┘   (ε 1→0)
  │    └─ verify_sosc_mee ──► sosc_* (7 files)
  ├─ psr_mee_refine ──► run_transfer_mee, psr_switch_score_mee, psr_refine_sigma_mee
  ├─ run_foc_mee ──► refresh_duals_mee ──► casadi_lt_mee
  │                  verify_pmp_mee ──► mee_dual_to_costate, mee_primer_switch
  └─ transfer_movie / gergaud_plot ──► mee_res_to_cart_res ──► elements_to_cart
```

## 4. Hubs (change these carefully)

| function | callers | note |
|---|---|---|
| `kepler_lt_params` | 24 | campaign physics constants |
| `module_root` | 21 | path helper |
| `optdef` | 19 | option-with-default; pure utility |
| `casadi_lt_mee` | 13 | the transcription; **also used by the CR3BP campaign** |
| `interp_warmstart` | 8 | **cross-campaign library candidate** — see §6 |
| `elements_to_cart` | 7 | coords |

## 5. Findings worth acting on

1. **Ten `.m` scripts live in `direct/results/`** — `check_incl_M1` and nine
   `task7c_*` one-off drivers. Scripts in a *results* directory is a structural
   smell; `results/` should hold artifacts. They are all entry points (nothing
   calls them). Move to `attic/` or a `scratch/` folder.
2. **Five genuine orphan candidates**, confirmed repo-wide as having no callers
   anywhere: `check_incl_M1`, `diag_lamv_mag`, `diag_primer`, `diag_primer2`,
   `fig_switch_convergence`. The three `diag_*` are the superseded Campaign-B
   refutation experiments; `DESIGN_dual_map.md` cites them as the evidence
   base, so archive rather than delete.
3. **45 of 96 non-test files are not reached by any test** (~47%). Mostly
   `viz/`, `results/` scripts and one-offs, which is defensible — but worth
   knowing before relying on any of them.
4. **Do not trust a single-folder orphan list.** `hamiltonian_const_check` and
   `cart_to_elements` both looked orphaned in the folder-local graph and are
   called from the CR3BP campaign and the tests respectively.

## 6. Library candidates this survey exposes

For the later cross-campaign extraction phase — recorded here, not acted on:

- **`interp_warmstart` vs `PSR/lib/warmstart_on_mesh`** — both re-map a solution
  onto a refined mesh. Two implementations of one idea in two campaigns; the
  strongest single candidate.
- **`optdef`** — already a shared utility in `lib/`; the other campaigns each
  have their own `getdef`/`gd` local. Trivially unifiable.
- **The certified-quantity guard** (Solve_Succeeded + defect ≤ 1e-8 + one-sided
  final mass) now exists in four places across the repo.
- **`elements_to_cart` / `rtn_frame`** — pure coordinate utilities with no
  campaign dependence.

Every one of these needs a byte-identity or reproduce-the-certified-result gate
before consolidation, per the standing rule recorded for `PSR/lib`.
