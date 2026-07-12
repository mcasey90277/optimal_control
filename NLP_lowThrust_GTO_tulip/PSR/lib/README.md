# PSR/lib — vendored machinery

These are **copies** of the machinery the PSR pipeline depends on, brought in on
**2026-07-12 (from git 5c0bdbc)** so that PSR is self-contained: `run_psr`,
`psr_export_data`, and `psr_movie` reach only `PSR/`, `PSR/lib/`, and the
external `pumpkyn` toolbox — nothing in `ms_band/` or `sundman_minfuel/`.
(Verified with `matlab.codetools.requiredFilesAndProducts` under
`restoredefaultpath`.)

The originals were **kept in place** (copy, not move) so the IFS folder and the
`ms_band` / `sundman_minfuel` campaign scripts keep working unchanged.

## Manifest (19 files)

| file | origin | role |
|---|---|---|
| `casadi_minfuel_sundman.m` | sundman_minfuel | direct solver (CasADi+IPOPT, Sundman trapezoid); NLP dual extraction |
| `cr3bp_lt_params.m` | sundman_minfuel | physics constants (muStar, lStar, tStar, Tmax, c, m0kg, Isp) |
| `gto_tulip_endpoints.m` | sundman_minfuel | GTO start + south-pole tulip target (uses pumpkyn) |
| `minfuel_at_tf.m` | sundman_minfuel | canonical per-t_f driver (energy→fuel homotopy) |
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
