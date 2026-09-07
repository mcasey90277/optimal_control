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
- `direct_mintime_elfo.m` — DIRECT (fmincon `solve_tfmin_nlp`) attempt. **RETIRED to `../../attic/` 2026-07-26**: recorded here as non-converging, superseded by the certified `casadi_mintime_freetf` Route-B anchor, and in any case broken (its addpath pointed at a folder that does not exist, so none of its dependencies resolved).

## Results (2026-07-13)

| method | tulip | ELFO retarget |
|---|---|---|
| single shooting (`mintime_solve`, `elfo_mintime`) | floors **~1e-3** (13-rev STM-product sensitivity) | stalls at s=0.05 |
| **multiple shooting** (`mintime_ms_*`) | **‖R‖=4e-9 ✓** (beats the wall, MS validated; J confirmed vs the solve). **SUPERSEDED 2026-08-25 by the shared-engine root** (`mintime_ms_bvp_probe`, `costate_common/ms_bvp` via `ms_tfmin`, K=60): ‖R‖=2.0e-11 in 10.3 s, **t_f = 6.290694 ND — exactly the certified direct reference** (old root was 47 s off at 6.290815); flown arrival 4.3 km; tfMin accepts (|Δz|=5.1e-7, at the ss floor for a 6.3 ND arc); **conjugate test PASS, 0 crossings** — first second-order verdict in the 25 mN regime. Probe artifact: `results/mintime_ms_bvp_probe.mat` | homotopy fights min-time sensitivity even with predictor–corrector — impractically slow |
| direct fmincon (`direct_mintime_elfo`) | does not converge / scale (t_f plunges, infeasible at usable N) | — |

**Bottom line:** the min-time MS machinery WORKS (tulip 4e-9) and answers
"is indirect min-time viable? yes." But **retargeting to the ELFO via shooting
is impractical** (sensitivity); the ELFO seed is still open. See
`../../../GTO_ELFO/direct/elfo/ELFO_RETARGET.md`. That candidate was **built
2026-07-15**: `casadi_mintime_freetf` (Route B, free t_f via the `cScale` slack
state) anchors both targets — ELFO 6.0962 ND, and tulip **5.8267 ND to the
backbone rendezvous** (`../../direct/lib/gen_tulip_mintime`). Note that the two
tulip min-times are to *different* points: this module's `mintime_params`
targets the max-ẏ point (6.2907 ND, the value `minfuel_config.tfMin` still
carries), whereas the min-fuel front targets the backbone point (5.8267 ND) —
see the root README's *Factor scale* note.

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

## TWO-ROOT finding (2026-08-26, systematic-debugging record)

During the `seed_from_z8`/`ms_tfmin` library-move verification, the probe
converged to a DIFFERENT root than the day before — from bit-identical seed
and code (both verified). Root cause: the K=60 seed sits on a **basin
boundary between two genuine extremals**, and last-bit arithmetic
differences between MATLAB runtimes select the root (desktop session →
t_f = 6.290449990, 11 iters; `-batch` → t_f = 6.290693961, 16 iters; each
environment internally deterministic, ‖R‖ ≈ 2e-11 both).

- **Canonical root** (`results/mintime_ms_bvp_probe.mat`): t_f = 6.290694,
  = the certified direct reference.
- **Alternate root** (`results/mintime_ms_bvp_probe_altroot.mat`):
  t_f = 6.290450 — **94 s FASTER**, physically screened (periselene
  6,039 km, ΔV 4.4663 vs 4.4665 km/s), tfMin-accepted, conjugate test not
  yet run on it.

**ADJUDICATED 2026-08-26 (same day): there is ONE extremal, not two.**
Flying both roots head-to-head: trajectory separation max 23.5 km over the
28-day arc (median 1.3 km), periselene identical, λ_v0 directions identical
to 0.000°. The cross-fly is the proof — each root reaches the target only
in the environment that produced it (4.3 km at home; 123–132 km in the
other), symmetrically. The 94 s "faster root" is NOT a faster transfer; it
is the same extremal expressed in a different integrator's arithmetic.

**The real finding — z8 identifiability in the ~40-rev regime:** the
end-to-end flow amplifies last-bit arithmetic differences beyond the root
separation, so a bare z8 determines t_f only to ~1e-4 ND (~1 min) and the
arrival only to ~100 km. Each environment's Newton is internally tight
(‖R‖ ~ 2e-11) about ITS OWN rendering of the extremal. tfMin accepts both
because its own arithmetic re-polishes within the identifiability tube.

Consequences: (1) **the certified direct t_f = 6.290694 STANDS** — the
faster-basin scare is withdrawn; (2) **a GTO-regime costate catalog cannot
ship bare z8 entries** — it must ship the full multiple-shooting solution
(junction states), the refined-library lesson
(direct-duals-not-shooting-seeds) reappearing at the mN scale; (3) quote
flagship t_f from the ms solve + direct agreement, never from a single
end-to-end flight; (4) the conjugate test on the "alt root" is moot (same
extremal — the canonical PASS covers it).
