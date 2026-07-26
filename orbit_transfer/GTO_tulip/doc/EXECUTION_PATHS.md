# GTO_tulip — execution paths

Survey 2026-07-26, read-only. Tool: `docs/tools/callgraph.py`.

**Scale: 266 `.m` — by far the largest campaign, and the indirect side (153) is
bigger than the direct side (81).** 37 more sit in two attics.

| area | files | status |
|---|---|---|
| `direct/sundman_minfuel` (+`refine`) | 23 + 11 | the certified min-fuel engine |
| `direct/PSR` (+`lib`) | 14 + 19 | switch refinement; **self-contained by vendoring** |
| `direct/movie` | 8 | renderers |
| `indirect/ms_band` | 57 | walled (near-min-time conditioning) |
| `indirect/ztl` | 47 | walled (min-time substrate dead) |
| `indirect/ifs` | 26 | walled (‖R‖≈0.023 terminal cluster) |
| `indirect/min_time` | 9 | validated (MS 4e-9) |
| `indirect/lowThrust_GTO_tulip` | 14 | base PMP shooting |

**130 of the 153 indirect files are the three walled attempts.** That is the
single largest block of code in the repo, and its recorded status is "did not
work." Not a cleanup target — it is the evidence base three campaign records
cite — but it should be understood as an archive with a live folder's name.

---

## 1. Headline: 21 duplicate basenames, and I undercounted them during the prune

**Correction.** In the 2026-07-26 prune I recorded: *"of the five names
duplicated between PSR/lib and sundman_minfuel/cr3bp_common, only
`gto_tulip_endpoints.m` is byte-identical; the other four DIFFER."* That was
true of the five names I compared and **materially understated the problem** —
I only diffed PSR/lib against `sundman_minfuel/` and `cr3bp_common/`, and
missed that PSR/lib also vendors from `sundman_minfuel/refine/` and from the
**indirect** `ms_band/`. The real count is 21 duplicate basenames, of which
**14 are byte-identical**.

The consequence matters for the library plan: I concluded then that the safe,
high-value consolidation targets did not exist here. They do — there are 14.

### The 14 identical copies (safe to consolidate)

| source | files vendored into `PSR/lib` |
|---|---|
| `sundman_minfuel/refine/` | `refine_loop`, `refine_sigma`, `prep_refine_seed`, `pmp_refine_indicator`, `warmstart_on_mesh` |
| **`indirect/ms_band/`** | `sms_eom`, `sms_jacobian_cs`, `sms_pack`, `sms_problem`, `sms_residual`, `sms_seed_duals`, `sms_unpack`, `verify_direct_pmp`, `beta_from_duals` |
| `cr3bp_common/` | `gto_tulip_endpoints` |

### The 6 that genuinely differ (need care)

| name | copies | note |
|---|---|---|
| `casadi_minfuel_sundman` | `sundman_minfuel` 16.6 kB vs `PSR/lib` 11.9 kB | 4.7 kB of divergence — the `returnModel`/`creg`, `vBox`/`rBox`, `boundSat` features never reached the vendored copy. **This is what silently broke `run_psr` stage 5c.** |
| `minfuel_at_tf` | `sundman_minfuel` vs `PSR/lib` | |
| `cr3bp_lt_params` | `PSR/lib` vs `cr3bp_common` | |
| `minfuel_config` | `PSR/lib` vs `cr3bp_common` | documented: `dirs` repointed |
| `pilot_rung_20mN` | tulip vs ELFO | same experiment name, different campaigns — legitimately distinct |
| `setup_paths` | 8 copies | correct by design; every module needs its own |

**A structural fact worth naming:** PSR — a *direct*-method pipeline — vendors
nine files from `ms_band`, an *indirect* campaign. The direct switch-refinement
machinery depends on the indirect multiple-shooting machinery. That is not
wrong (PSR is PMP-steered, so it wants the PMP/shooting tooling), but it is an
architectural dependency no folder name advertises.

## 2. Execution path

```
run_psr   [PSR — the direct min-fuel refinement pipeline]
  insertion_states ....... [sundman_minfuel]  -> cr3bp_lt_params, gto_tulip_endpoints
  minfuel_config ......... [PSR/lib]
  minfuel_at_tf .......... [PSR/lib]  -> casadi_minfuel_sundman, cr3bp_lt_params
  prep_refine_seed ....... [PSR/lib]  -> casadi_minfuel_sundman
  refine_loop ............ [PSR/lib]  (rounds: score -> refine sigma -> re-solve)
  psr_ipopt_certify ...... [PSR]      2nd-order delta_w observation
  [stage 5c FOC gate] .... prints a pointer; verify_common is outside PSR's
                           path boundary by design (fixed 2026-07-26)
  psr_export_data / psr_movie
```

Note `run_psr` crosses into `sundman_minfuel` for exactly one function
(`insertion_states`) — the self-containment claim in `PSR/lib/README.md` is
about `PSR/lib` covering the *machinery*, not about zero external references.

## 3. Findings worth acting on

1. **14 identical duplicates are the cleanest consolidation target in the
   repo** — no behaviour can change if the files are byte-identical. The
   blocker is not risk but path policy: PSR's self-containment is a deliberate,
   verified design (`requiredFilesAndProducts` under `restoredefaultpath`).
   Consolidating means deciding whether that guarantee is still wanted. My
   recommendation: keep the guarantee, but source the identical files from one
   place via an explicit `psr_vendor_check` test that asserts byte-identity
   against the originals and fails loudly when they drift — which is exactly
   the failure that broke stage 5c.
2. **`casadi_minfuel_sundman`'s 4.7 kB divergence is a live hazard**, not a
   historical one: the two copies will keep drifting, and the next upstream fix
   will silently miss PSR again.
3. **53 of 64 non-test direct files are unreached by any test** (~83%).
4. **7 orphan candidates, all in `movie/`** (`animate_*`, `gen_compare_data`,
   `psr_homotopy_movie`) — manual renderers. Fourth campaign running where the
   folder-local orphan list is all figure/movie scripts; the pattern is now
   reliable enough to state as a rule.
5. **Tool defect found and fixed here.** `callgraph.py` keyed by basename, so
   duplicate copies collapsed and their callers merged onto whichever the walk
   saw last — it reported `cr3bp_lt_params [PSR/lib]` as having 29 callers when
   that is the sum across two copies. It now detects and reports duplicates
   before the hub table, with a warning that the hub counts are untrustworthy
   for those names. **The earth and CR3BP surveys are unaffected** (no
   duplicate basenames in either tree).

## 4. Library candidates

- **The 14 identical files** — mechanical, gated by byte-identity. Highest
  confidence, blocked only on the path-policy decision in finding 1.
- **`warmstart_on_mesh` (here) vs `interp_warmstart` (earth)** — the
  cross-campaign candidate flagged in the earth survey; both re-map a solution
  onto a refined mesh. These are *not* byte-identical and live in different
  transcriptions, so this one needs real comparison, not a move.
- **`insertion_states` / `gto_tulip_endpoints`** — endpoint definitions already
  half-shared with `cr3bp_common`; finish the job.
