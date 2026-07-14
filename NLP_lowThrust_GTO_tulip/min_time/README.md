# min_time — min-time (PMP, always-burn) CR3BP transfers

Self-contained module for **minimum-time** low-thrust CR3BP rendezvous, used as
the robust *root* for building min-energy / min-fuel solutions to new targets
(the tulip's own lineage: min-time → min-energy → min-fuel).

## Why min-time as the root

Min-time is **always-burn** (throttle ≡ 1, no switches on the interior for this
regime), so it avoids the two failure modes that stall the fixed-tf
energy-target homotopy when retargeting GTO→ELFO:
- **No throttle** → no saturation/edge sensitivity.
- **tf floats** → every intermediate target on a target-homotopy is reachable
  (just at a different tf), so it does not hit the fixed-tf "can't-reach-
  terminal" wall as the terminal moves into the Moon's gravity well.

## What pumpkynPie already provides

`pumpkyn.cr3bp.tfMin(rv0, rvf, [λ;tf], Tmax, c, muStar)` solves min-time PMP
rendezvous between **any** two rotating-frame states (analytic-STM single
shooting) — target-agnostic, so it handles an ELFO point directly. Its only
limitation is a hardcoded 100-evaluation budget that stalls hard multi-rev
cases. Its propagator/EoM (`tfMinProp`, `tfMinEoM`, `minDeltaV`) are reused
here.

## Files

- `setup_paths.m` — adds this dir + pumpkyn/src.
- `mintime_params.m` — CR3BP constants + GTO departure + tulip max-ẏ target
  (self-contained; mirrors `gto_tulip_endpoints`).
- `mintime_solve.m` — real-budget min-time TPBVP **single-shooting** solver:
  pumpkyn's exact residual/Jacobian (rendezvous, λ_m(tf)=0, H(tf)=0) via the
  propagated 14×14 STM, `lsqnonlin` TRR (vs pumpkyn's 100-eval cap).
- `mintime_ms_residual.m` / `mintime_ms_seed.m` — min-time **multiple-shooting**
  residual + block Jacobian (on pumpkyn per-arc STMs) and its arc-chopped seed.
- `mintime_ms_gate.m` — validates seed continuity + analytic block-J vs FD.
- `mintime_ms_solve.m` — converge the tulip min-time via MS (uses the genericized
  `ztl_ms_solve_tr` trust-region solver, `prob.resFun`).
- `mintime_ms_elfo.m` — homotope the MS target tulip→ELFO (predictor–corrector).
- `elfo_mintime.m` — single-shooting version of the tulip→ELFO homotopy.
- `direct_mintime_elfo.m` — DIRECT (fmincon `solve_tfmin_nlp`) attempt.

## Results (2026-07-13)

| method | tulip | ELFO retarget |
|---|---|---|
| single shooting (`mintime_solve`, `elfo_mintime`) | floors **~1e-3** (13-rev STM-product sensitivity) | stalls at s=0.05 |
| **multiple shooting** (`mintime_ms_*`) | **‖R‖=4e-9 ✓** (beats the wall, MS validated; J confirmed vs the solve) | homotopy fights min-time sensitivity even with predictor–corrector — impractically slow |
| direct fmincon (`direct_mintime_elfo`) | does not converge / scale (t_f plunges, infeasible at usable N) | — |

**Bottom line:** the min-time MS machinery WORKS (tulip 4e-9) and answers
"is indirect min-time viable? yes." But **retargeting to the ELFO via shooting
is impractical** (sensitivity); the ELFO seed is still open. See
`../elfo/ELFO_RETARGET.md`. Leading next candidate: a **direct min-time
collocation** (free t_f) built on `casadi_minfuel_sundman`.

## Run

```matlab
cd min_time
mintime_ms_gate     % validate machinery (seed + Jacobian)
mintime_ms_solve    % tulip min-time via MS -> results/mintime_tulip_ms.mat (4e-9)
```

## Caveat

The STMs are integrated across throttle-switch events without a saltation
correction, so J is exact only while the arc is all-burn (S<0 throughout);
`out.nSwitch` reports the converged arc's switch count. The tulip min-time is
all-burn (nSwitch 0), as expected.
