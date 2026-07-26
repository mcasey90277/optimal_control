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

## Results layout (2026-07-09 migration)

All result `.mat`s live under `results/` with collision-free milli-factor
names from `minfuel_config`: `results/energy/energy_f####.mat` (backbone),
`results/minfuel/minfuel_f####_<branch>.mat` (provenance-stamped, from
`minfuel_at_tf`), `results/minfuel/legacy_ms_f####.mat` (pre-migration
solutions, no meta), `results/fronts/` (front collections),
`results/plots/`, `results/logs/`. Only the two canonical artifacts stay in
the library root: `sundman_minfuel_certified.mat` and
`minfuel_from_energy_seed.mat`. **The canonical toolchain is
`minfuel_config` / `minfuel_at_tf` / `aggregate_front` +
`orchestrate/*.sh`** — the older drivers below (`solve_tf_minfuel`,
`run_tf_sweep/front/2anchor`, `build_energy_backbone`) still reference the
pre-migration filenames and are superseded (atticked in cleanup Phase 1; see
`../../process/CODE_CLEANUP_PLAN.md`).

## Notes
- The full method write-up and the "two walls" (dynamics vs objective) analysis
  are in `../../process/LOW_THRUST_MINFUEL_CAMPAIGN.md`.
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
  fail there (bifurcation-rich). Mapped front (all machine-tight, ΔV monotone):
  1.12×=3.83 (12 sw) / 1.14×=3.49 / 1.15×=3.37 / 1.20×=3.24 / 1.25×=3.14 km/s,
  + the min-time anchor (4.4665, 0 sw); switches drop toward the band. (1.12x
  switch count: see the Jul-10 adjudication — 10 PMP-certified switches + 1
  near-graze throttle dip miscounted by the s>0.5 threshold;
  `../ms_band/MS_BAND_CAMPAIGN.md`.) The band
  is an open item (indirect/multiple-shooting is the proposed next tool). Full
  color-coded figure `results/plots/front_full_verified.png`; account in the campaign doc
  "Down-sweep CRACKED" + "transition band".
- **Per-t_f optimality verification.** `verify_tf_front.m` checks each solution
  against Pontryagin's first-order conditions from its own KKT-dual costates
  (empirical-β switching law) and colors the front by PMP-certified vs not.
  Use it to tell a genuine local extremal from a merely-converged point. Full
  five-layer plan in `OPTIMALITY_VERIFICATION_PLAN.md`. **Transversality is
  checked RELATIVE** — |λ_m(τ_f)|/max|λ_m| ≤ 1e-3 — because the costates are
  known only up to a positive scale; an absolute gate on λ_m(τ_f) wrongly failed
  the fewer-switch down-band points (1.12×, 1.14×), whose overall costate scale
  is larger (their raw λ_m(τ_f) ≈ −4e-3 is a genuine zero against max|λ_m| ≈ 31).
  With the scale-invariant gate the full certified band is **1.12×–1.25× + 1.85×**.

## ELFO transfer (GTO → lunar elliptical frozen orbit)

The same campaign retargeted to a **lunar ELFO** (the south-pole nav orbit:
sma 12000 km, ecc 0.69, inc 56.5°, argp 90°) instead of the tulip. The lunar
capture leg needs a **two-primary** Sundman clock and a **free t_f**, so ELFO
runs on a redesigned solver rather than `casadi_minfuel_sundman`:

| file | role |
|---|---|
| `casadi_energy_freetf.m` | ELFO **core solver**: free-t_f (constant slack state `cScale`, banded KKT), two-primary clock κ=(r₁⁻�qᵁ+(r₂/D)⁻ᵠ)⁻ᵖ/ᵠ, Moon-gravity homotopy, Bertrand–Épénoy ε objective; optional pinned t_f (`tfTarget`) for well-posed energy |
| `gen_elfo_energy_gravhom.m` | make the GTO→ELFO **min-energy seed** via the gravity-homotopy ladder (gravity off → clock on → retarget → gravity on) |
| `gen_elfo_energy_tfsweep.m` | t_f-continuation mapping the ELFO energy band + banking a seed per grid t_f (built; not yet run) |
| `gen_elfo_minfuel.m` | ε:1→ε_min energy→fuel homotopy core (from an energy seed) |
| `run_elfo_minfuel.m` | **entry driver** (run_psr analog): solve → data export → verify → movie; target-tagged output names |
| `elfo_export_data.m` | data products → `../PSR_data/psr_data_<target>_...` (costates from the two-primary KKT duals) |
| `verify_elfo_seed.m` | independent endpoint + solver-free defect check (energy or fuel `.mat`) |

**Results (t_f = 1.20× = 33.46 d):** min-energy seed defect **1.8e-15**;
min-fuel **34-switch bang-bang**, 99.6% edge, propellant **~2.16 kg**, ΔV **~3.24
km/s**, defect **5.7e-14** — all independently verified. Output names carry the
target, e.g. `minfuel_ELFO_tf1p200_sw34_minEps0.mat`, so tulip and ELFO products
never collide. Full build record: `../elfo/ELFO_RETARGET.md`.

**Run:**
```matlab
cd sundman_minfuel
run_elfo_minfuel        % edit the PARAMETERS block: target / factor / epsMin / movieMode
```

**Deferred** (tulip single-primary machinery, not yet ported to the two-primary
free-t_f model): switch-localization mesh refinement and the full 16-dim per-arc
PMP-propagation certificate. The min-fuel **t_f-grid** convergence map (the energy
band is wider than the ε=0-convergent band) is the next campaign step:
`gen_elfo_energy_tfsweep` then `run_elfo_minfuel` per grid t_f.
