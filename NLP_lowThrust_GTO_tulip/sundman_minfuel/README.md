# sundman_minfuel — certified sharp bang-bang min-fuel GTO→tulip solver

Self-contained library for the **minimum-fuel low-thrust GTO → south-pole
tulip transfer** in the Earth–Moon CR3BP, on the full ~40-revolution spiral.
The optimal control is sharp bang-bang (~25 switches); this code solves it to
machine precision by combining three ingredients:

1. **Sundman regularization** — change the independent variable time → τ with
   `dt/dτ = κ = r₁^pSund` (r₁ = Earth distance), carrying time as an 8th state.
   Tames the near-perigee 1/r³ Hessian terms and auto-concentrates nodes at
   perigee. τ_f is held **fixed** (a free τ_f makes a dense KKT column → MUMPS
   OOM); fixed transfer time is enforced by the terminal condition t(τ_f)=t_f.
2. **Energy→fuel homotopy** (Bertrand–Épénoy) — objective
   `J(ε) = ∫s dt − ε∫s(1−s) dt`; ε=1 is energy (smooth), ε=0 is fuel
   (bang-bang). Sweep ε:1→0, warm-starting each solve.
3. **No-resample seed** — map a collocation-feasible time-mesh solution into τ
   using its **own** nodes (no interpolation), avoiding the ~1e-2 defect floor
   that pins IPOPT in restoration.

Solved via **CasADi + IPOPT** (exact sparse Jacobian and Hessian by AD).

## Certified result
`Solve_Succeeded`, max defect **2.4e-14**, terminal error **0**, **25
switches**, 99.6% bang-bang, propellant **2.2640 kg**, ΔV **3.3696 km/s**
(15 kg, 25 mN, Isp 2100 s, t_f = 1.15× min-time). For comparison the min-time
baseline (always-on) is 4.4665 km/s / 2.9247 kg — the min-fuel solution saves
~1.10 km/s (24.6%). Solution stored in `sundman_minfuel_certified.mat`.

## Files
| file | role |
|---|---|
| `cr3bp_lt_params.m` | Earth–Moon CR3BP + low-thrust ND constants (μ*, l*, t*, c, T_max) |
| `gto_tulip_endpoints.m` | departure GTO / arrival tulip states (needs pumpkyn) |
| `sundman_seed_map.m` | no-resample map of a time-mesh solution → Sundman coords |
| `casadi_minfuel_sundman.m` | **core solver**: Sundman collocation, ε-homotopy objective, IPOPT |
| `sundman_homotopy.m` | guarded ε:1→0 sweep (keeps best, discards loose steps) |
| `run_certified_minfuel.m` | **entry point** — reproduces the certified result end-to-end |
| `setup_paths.m` | adds pumpkyn (problem setup only) |
| `minfuel_from_energy_seed.mat` | collocation-feasible time-mesh seed |
| `sundman_minfuel_certified.mat` | the certified bang-bang solution |

## Run
```matlab
cd sundman_minfuel
best = run_certified_minfuel;      % ~10–15 min; writes sundman_minfuel_certified.mat
```
CasADi is loaded from `$CASADI_PATH` (default `~/casadi-3.7.0`). A movie of the
solution (with a running ΔV meter) is produced by
`../movie/animate_sundman_minfuel.m`.

## Costates & PMP verification
`casadi_minfuel_sundman` returns the **discrete costates** — the KKT multipliers
of the dynamics-defect constraints — so the direct solution can be checked
against Pontryagin's principle:
- `out.lamDef` `[8×N]` — `[λ_r; λ_v; λ_m; λ_t]` per interval (up to a positive
  mesh-weight scaling and a global sign); `out.lamAll` — the full stacked
  multiplier vector (includes the throttle-bound duals that encode the
  switching-function sign law).
- `out.primerAlignDeg` — mean angle between the NLP thrust direction and the
  costate primer `-λ_v/‖λ_v‖` on burn arcs (scale-invariant).
- `out.lamMassEnd` — terminal mass-costate proxy (transversality, ≈0).

**Verified** on the certified solution: primer alignment **0.058°**,
transversality **−1.7×10⁻⁷** — the direct solution meets the PMP direction and
transversality conditions to a hundredth of a degree. **Remaining:** the
scale-dependent switching-function sign law and an independent adjoint check.
Full plan and status in `TIER1_PMP_CERTIFICATION_SCOPE.md`; `run_tf_sweep.m`
saves `lamDef` per t_f.

## Notes
- The full method write-up and the "two walls" (dynamics vs objective) analysis
  are in `../LOW_THRUST_MINFUEL_CAMPAIGN.md`.
- **ΔV–time front.** `run_tf_front.m` continues the certified basin in small
  t_f steps (cleaner than `run_tf_sweep.m`'s energy-continuation, which
  scatters). The optimal *family changes* with t_f, so no single basin threads
  the whole range.
- **Down-sweep toward min-time** (the robust method — bang-bang continuation
  MEX-crashes going down): two phases. `build_energy_backbone.m` / `energy_step.m`
  chain the SMOOTH energy solution down in small t_f steps (convex ⇒ crash-free),
  saving `energy_<f>.mat` per t_f (process-isolated + watchdog to survive MEX
  crashes/hangs). Then `solve_tf_minfuel.m` (one per t_f, PARALLEL): tight
  multiplier re-clean of the backbone energy, then fine energy→fuel sharpen to
  certified min-fuel. First down-step: 1.13× = 3.5955 km/s / 24 sw / defect
  4.5e-14 (monotone-correct vs 1.15× = 3.370).
- **`direct_build_minfuel.m`** — solve a t_f from SCRATCH (fresh burn+coast →
  min-energy → homotopy, no continuation). Works only where the coast is small
  enough that the warm start is clean; hits the restoration wall in the hard band.
- **Hard transition band ~1.01–1.11×**: where the control reorganizes from
  many-switch to always-burn min-time. Both continuation and from-scratch solves
  fail there (bifurcation-rich). Certified front is 1.12×→1.25× + the min-time
  anchor; the band is an open item (indirect/multiple-shooting is the proposed
  next tool). Full account: campaign doc "Down-sweep CRACKED" + "transition band".
- **Per-t_f optimality verification.** `verify_tf_front.m` checks each solution
  against Pontryagin's first-order conditions from its own KKT-dual costates
  (empirical-β switching law) and colors the front by PMP-certified vs not.
  Use it to tell a genuine local extremal from a merely-converged point. Full
  five-layer plan in `OPTIMALITY_VERIFICATION_PLAN.md`.
