# PSR/lib — vendored machinery

These are **copies** of the machinery the PSR pipeline depends on, brought in on
**2026-07-12 (from git 5c0bdbc)** so that PSR would be self-contained: `run_psr`,
`psr_export_data`, and `psr_movie` reaching only `PSR/`, `PSR/lib/`, and the
external `pumpkyn` toolbox. That was verified at the time with
`matlab.codetools.requiredFilesAndProducts` under `restoredefaultpath`.

The originals were **kept in place** (copy, not move) so the IFS folder and the
`ms_band` / `sundman_minfuel` campaign scripts keep working unchanged.

## ⚠ Self-containment was lost on 2026-07-15 — read this before trusting the manifest

**It is no longer true that PSR reaches nothing in `sundman_minfuel/`.** Three
days after the vendoring, the insertion-points feature (commits `5089931`,
`a35a8ba`) added

```matlab
addpath(fullfile(here, '..', 'sundman_minfuel'));
```

to `run_psr`, `psr_run_one` and `gen_energy_seed`, to reach `insertion_states`.
`addpath` **prepends**, so that one line put the entire upstream folder *ahead*
of `PSR/lib` on the path. Measured 2026-07-26 on `run_psr`'s own path setup:

| name | resolves to | consequence |
|---|---|---|
| `casadi_minfuel_sundman` | **upstream** `sundman_minfuel/` | the vendored copy is **dead code** |
| `minfuel_at_tf` | **upstream** `sundman_minfuel/` | the vendored copy is **dead code** |
| `cr3bp_lt_params`, `minfuel_config`, `gto_tulip_endpoints`, `sms_*`, `refine_*` | `PSR/lib/` | as designed (upstream has no such file, or it is elsewhere) |

So PSR has been running the **upstream** solver since 2026-07-15, not the
frozen 2026-07-12 snapshot this directory implies. The vendoring did not fail
loudly; it was quietly overridden by a later, unrelated feature.

As of 2026-07-26 `insertion_states` moved to `cr3bp_common` and is **vendored
here** (file 20), so it is no longer the reason for that `addpath` — the
upstream *solver* now is, and the comment on each of the three call sites says
so. Removing the `addpath` would make the two dead files live again, i.e. it
would silently swap PSR's solver back to the frozen snapshot. **That is a
behaviour change on a certified pipeline and must not be done as cleanup.**
It is recorded as an open decision in `orbit_transfer/CODE_STRUCTURE.md`.

`PSR/tests/test_psr_vendor_drift.m` guards every file listed below.

## Manifest (20 files)

| file | origin | role |
|---|---|---|
| `casadi_minfuel_sundman.m` | sundman_minfuel | direct solver (CasADi+IPOPT, Sundman trapezoid); NLP dual extraction. **DEAD — shadowed by the upstream copy; see the warning above** |
| `cr3bp_lt_params.m` | sundman_minfuel | physics constants (muStar, lStar, tStar, Tmax, c, m0kg, Isp) |
| `gto_tulip_endpoints.m` | cr3bp_common | GTO start + south-pole tulip target (uses pumpkyn) |
| `insertion_states.m` | cr3bp_common | **added 2026-07-26.** Single source for the GTO departure + tulip/ELFO insertion states. Vendored when it moved out of `sundman_minfuel` into `cr3bp_common`, which PSR does not put on its path (doing so would shadow the PSR-owned `cr3bp_lt_params` / `minfuel_config` variants) |
| `minfuel_at_tf.m` | sundman_minfuel | canonical per-t_f driver (energy→fuel homotopy). **DEAD — shadowed by the upstream copy; see the warning above** |
| `minfuel_config.m` | sundman_minfuel | campaign constants + schedules. **EDITED**: `dirs` repointed to `../../sundman_minfuel/results` (energy backbones referenced in place) |
| `refine_loop.m` | sundman_minfuel/refine | PSR refinement loop. **PSR-owned** (carries the `outDir`/`solFile` additions) |
| `pmp_refine_indicator.m` | sundman_minfuel/refine | PMP switch-localization score (the refinement steer) |
| `refine_sigma.m` | sundman_minfuel/refine | mesh refiner |
| `warmstart_on_mesh.m` | sundman_minfuel/refine | no-resample warm start onto a refined mesh |
| `prep_refine_seed.m` | sundman_minfuel/refine | normalize a direct solution into refine-seed layout |
| `sms_seed_duals.m` | ms_band | KKT-dual → node-costate map (adjudicated mode 'd') |
| `beta_from_duals.m` | ms_band | costate scale (β) fit |
| `sms_eom.m` | ms_band | 16-dim Sundman PMP EOM (costate propagation) |
| `sms_problem.m` | ms_band | problem-struct factory |
| `sms_pack.m` / `sms_unpack.m` | ms_band | MS unknown ⇄ node layout |
| `sms_residual.m` | ms_band | MS residual (pulled in by the solver stack) |
| `sms_jacobian_cs.m` | ms_band | complex-step Jacobian (pulled in by the solver stack) |
| `verify_direct_pmp.m` | ms_band | first-order PMP verifier. **PSR-owned** (adjudication driven from `run_psr` `verifyOpts`) |

## Drift caveat

Two of these are **actively developed for PSR** and should be edited HERE, not
in the origin folders: `refine_loop.m` and `verify_direct_pmp.m` (and the
`minfuel_config.m` copy, which is deliberately different from its origin).
`casadi_minfuel_sundman.m` also carries a PSR-only addition (2026-07-12): it
returns `out.regHistory`, IPOPT's per-iteration Hessian regularization delta_w,
read by `psr_ipopt_certify.m` for the native-inertia local-min certificate. The
rest are **stable machinery** — if you ever need to sync a bug fix from the
origin (e.g. a `casadi_minfuel_sundman` fix), re-copy that one file and note it
here. The origins as of the copy were git 5c0bdbc.

## Not vendored (referenced in place)

- **Energy backbones / seed library**: `sundman_minfuel/results/energy/*.mat`
  and `.../minfuel/*.mat` (data, ~30 MB). The `minfuel_config` copy points
  there. If that tree is ever moved, update `minfuel_config.m` `dirs`.
- **pumpkyn toolbox**: `proj7/external/pumpkyn/src` (third-party, shared).
  Added to the path by `PSR/setup_paths`.

## Divergence from the originals (recorded 2026-07-26)

These are a **frozen snapshot**, and they have since diverged from the
`sundman_minfuel/` originals. That is the design working as intended, but it
has one consequence worth stating explicitly, because it has already caused a
silently-dead feature:

- `casadi_minfuel_sundman.m` here carries the `regHistory` capture but **not**
  the `returnModel` / `creg` constraint-registry hook added upstream on
  2026-07-25. That is why this copy could not serve the generic FOC gate — but
  it is *not* why the gate does not run in PSR. The gate does not run because
  `run_foc_tulip` and the whole `verify_common` layer sit outside PSR's path
  boundary, so the call threw "Undefined function". `run_psr` now prints a
  pointer instead of attempting it.

- **Which copy wins depends on the path state, and the two states disagree.**
  Measured 2026-07-26:

  | path state | `casadi_minfuel_sundman` resolves to |
  |---|---|
  | `PSR/setup_paths` alone | `PSR/lib/` ✅ as designed |
  | after the entry-point `addpath(../sundman_minfuel)` | **upstream** ❌ |

  An earlier note here claimed "PSR always runs these files, never the upstream
  ones." That was measured in the first state and is **wrong for real runs**:
  `run_psr`, `psr_run_one` and `gen_energy_seed` — every entry point — perform
  that `addpath` themselves, so the second state is the one that executes.
  `insertion_states` is unaffected: upstream no longer defines it, so the
  vendored copy here wins in both states.

Before folding these back into the shared sources, the tulip TODO's standing
rule applies: do it **with a reproduce-the-certified-result gate**, not as a
tidy-up.
