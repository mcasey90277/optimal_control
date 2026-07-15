# GTO low-thrust transfer — grand roadmap

**Goal:** BOTH direct and indirect methods for each of min-time, min-energy
(varying t_f), and min-fuel (varying t_f), to BOTH targets: the south-pole
**tulip** and the lunar **ELFO** (elliptical lunar frozen orbit, proj7
`im_elfo_optimum`: sma 12000 km, ecc 0.69, inc 56.5, argp 90).

## Status matrix (2026-07-15)

Legend: ✅ done/validated · 🟡 partial · ⬜ open

| problem | method | GTO→tulip | GTO→ELFO |
|---|---|---|---|
| **min-time** | direct | 🟡 `attic/solve_tfmin_nlp` (fmincon; converges at fine mesh, "easy case") | ✅ **hard all-burn `tfMin_ELFO` = 6.0962 ND = 27.02 d** (`elfo/casadi_mintime_freetf.m` + `gen_elfo_mintime.m`; s≡1, min t(τ_f); machine-tight def 1.7e-15, rendezvous 0, mass-identity 1e-16, independently verified; below the 6.23 ND energy floor). Anchors the ELFO factor scale + the front's 0-switch endpoint |
| | indirect | ✅ single-shoot = pumpkyn to 8 sig figs; **MS 4e-9** (`min_time/mintime_ms_*`) | ⬜ MS retarget fights shooting sensitivity (the direct anchor above is now a candidate seed) |
| **min-energy** (var t_f) | direct | ✅ energy backbones factor 1.12–1.95 (`sundman_minfuel/results/energy`) | ✅ **gravity-homotopy seed 1.8e-15** (`elfo/gen_elfo_energy_gravhom.m` → `elfo/results/energy_elfo_freetf.mat`, tf 33.5 d, 15.7% prop) |
| | indirect | ✅ Sundman-MS 75 mN anchor **4.8e-10** (`ztl/results/z1_sun_anchor_75mN.mat`); band via costates 🟡 | ⬜ (energy seed now exists; not yet run) |
| **min-fuel** (var t_f) | direct | ✅ PSR pipeline, 3- & 25-switch bang-bang certified, band [1.12,1.95] | ✅ **ΔV–time FRONT mapped** (2026-07-15): 11/14 factors ε=0 bang-bang, machine-tight (edge ~99.6%, def 1e-15..1e-12); **min 2.693 km/s @ 1.73×/48 d** (12.3% prop), monotone from 3.344 (1.11×/31 d) then flat; 3 gaps at the 1.65–2.0× fold (timed out even at 60 min = hard fold AT the optimum). Pipeline: `elfo/elfo_energy_sweep.sh`→`elfo_batch.sh 0 energy`→`elfo_collect_summary` (`results/elfo_batch_summary_minEps0.mat`) |
| | indirect | 🟡 IFS/ms_band: 1.12x = 10 switches certified; band = conditioning wall | ⬜ (energy seed exists; not yet run) |

## GTO→ELFO min-time ANCHORED (2026-07-15)

`tfMin_ELFO = 6.0962 ND (27.02 d)`, hard all-burn (s≡1), machine-tight and
independently verified. Built Route B after Route A (energy→time continuation)
walled: stepping the pinned t_f down from the energy seed rode the throttle up
cleanly to 6.2278 ND at edge 57.5%, then hit the near-min-time conditioning wall
(inf_du→1e10) — the same wall tulip has below 1.12× — giving only a loose upper
bound. Route B poses the true min-time directly: throttle pinned s≡1, t_f free via
the cScale slack, objective min t(τ_f); the energy 6.23 ND rung warm-starts it and
IPOPT drives t_f to the all-burn floor. Solver `elfo/casadi_mintime_freetf.m`
(sibling of `casadi_energy_freetf`), driver `elfo/gen_elfo_mintime.m`, result
`elfo/results/mintime_elfo.mat`. Loose and tight warm-start modes agree to 4 sig
figs (independent cross-check). Notable: ELFO's min-time is *slightly shorter*
than tulip's (6.0962 vs 6.2907 ND).

**Front relabel (the payoff):** the mapped ELFO min-fuel front was plotted vs
`factor = tf/tfMin_tulip`; it now relabels into ELFO's own units
`tf/tfMin_ELFO`. Endpoints: 1.11× tulip (6.98 ND) = **1.145× ELFO**; front
minimum 1.73× tulip (10.9 ND) = **1.79× ELFO**. Full record:
`elfo/ELFO_RETARGET.md` (Min-time anchor section). Still ⬜: the ELFO min-time
**indirect** cell (Route C).

## The ELFO-column blocker is CLEARED (2026-07-13)

The **GTO→ELFO min-ENERGY seed** — the one missing input the whole ELFO column
was blocked on — is MADE: `elfo/results/energy_elfo_freetf.mat`,
defect **1.8e-15**, independently verified. Built by the **gravity-homotopy
ladder** `elfo/gen_elfo_energy_gravhom.m` on the new free-t_f two-primary solver
`elfo/casadi_energy_freetf.m` (a GPT-5.6-terra + Gemini 3.1 Pro design
review killed the earlier direct-min-time-collocation plan as a detour and
prescribed this route instead). Full build record + the two extra fixes (pin
t_f; leg order clock-on-before-retarget with gravity off) in `elfo/ELFO_RETARGET.md`.

**Now open (unblocked):** min-fuel GTO→ELFO — re-run `casadi_energy_freetf` from
the energy seed with ε:1→0. Then the indirect ELFO cells.

## Key module map

- `sundman_minfuel/` — direct min-energy backbones + tulip min-fuel (Sundman
  collocation); the shared solver engine (`cr3bp_lt_params`, `minfuel_config`)
  that the per-target deliverables reference on the path.
- `elfo/` — GTO→ELFO direct deliverable (self-contained sibling of PSR;
  shared-path to the sundman_minfuel engine; reorg 2026-07-14):
  `casadi_energy_freetf.m` (free-t_f, two-primary clock, gravity homotopy);
  `gen_elfo_energy_gravhom.m` (4-leg ladder → energy seed); `gen_elfo_minfuel.m`
  (ε→0 fuel); `run_elfo_minfuel.m` (entry: solve→export→verify→movie, target-tagged);
  `elfo_export_data.m` (data products); `gen_elfo_energy_tfsweep.m` (tf-band map).
  Terminal-runnable, crash-robust + resumable: `elfo_energy_sweep.sh` (energy seed
  band) → `elfo_batch.sh 0 energy` (per-factor ε→0 bang-bang) → `elfo_collect_summary`;
  post-hoc movies `elfo_movies.sh all` (`elfo_render_movies`→`elfo_movie`, no re-solve).
  Build record `ELFO_RETARGET.md`.
- `PSR/` — PMP-Steered Refinement (direct GTO→tulip min-fuel pipeline).
- `ztl/` — indirect Sundman multiple-shooting (energy anchor 4.8e-10).
- `min_time/` — min-time (single + multiple shooting); tulip MS validated 4e-9.
- `ms_band/`, `ifs/` — indirect min-fuel band attempts (conditioning wall).
