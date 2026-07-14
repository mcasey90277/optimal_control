# GTO low-thrust transfer — grand roadmap

**Goal:** BOTH direct and indirect methods for each of min-time, min-energy
(varying t_f), and min-fuel (varying t_f), to BOTH targets: the south-pole
**tulip** and the lunar **ELFO** (elliptical lunar frozen orbit, proj7
`im_elfo_optimum`: sma 12000 km, ecc 0.69, inc 56.5, argp 90).

## Status matrix (2026-07-13)

Legend: ✅ done/validated · 🟡 partial · ⬜ open

| problem | method | GTO→tulip | GTO→ELFO |
|---|---|---|---|
| **min-time** | direct | 🟡 `attic/solve_tfmin_nlp` (fmincon; converges at fine mesh, "easy case") | ⬜ fmincon doesn't scale (t_f plunges) |
| | indirect | ✅ single-shoot = pumpkyn to 8 sig figs; **MS 4e-9** (`min_time/mintime_ms_*`) | ⬜ MS retarget fights shooting sensitivity |
| **min-energy** (var t_f) | direct | ✅ energy backbones factor 1.12–1.95 (`sundman_minfuel/results/energy`) | ✅ **gravity-homotopy seed 1.8e-15** (`sundman_minfuel/gen_elfo_energy_gravhom.m` → `results/energy_elfo_freetf.mat`, tf 33.5 d, 15.7% prop) |
| | indirect | ✅ Sundman-MS 75 mN anchor **4.8e-10** (`ztl/results/z1_sun_anchor_75mN.mat`); band via costates 🟡 | ⬜ (energy seed now exists; not yet run) |
| **min-fuel** (var t_f) | direct | ✅ PSR pipeline, 3- & 25-switch bang-bang certified, band [1.12,1.95] | 🟡 **ε=0 reached at tf=33.5 d (1.20×)**: 34-switch bang-bang, 14.5% prop, def 5.7e-15, verified (`gen_elfo_minfuel.m`→`minfuel_elfo.mat`). tf-GRID map pending (energy band ⊋ fuel band) |
| | indirect | 🟡 IFS/ms_band: 1.12x = 10 switches certified; band = conditioning wall | ⬜ (energy seed exists; not yet run) |

## The ELFO-column blocker is CLEARED (2026-07-13)

The **GTO→ELFO min-ENERGY seed** — the one missing input the whole ELFO column
was blocked on — is MADE: `sundman_minfuel/results/energy_elfo_freetf.mat`,
defect **1.8e-15**, independently verified. Built by the **gravity-homotopy
ladder** `gen_elfo_energy_gravhom.m` on the new free-t_f two-primary solver
`sundman_minfuel/casadi_energy_freetf.m` (a GPT-5.6-terra + Gemini 3.1 Pro design
review killed the earlier direct-min-time-collocation plan as a detour and
prescribed this route instead). Full build record + the two extra fixes (pin
t_f; leg order clock-on-before-retarget with gravity off) in `PSR/ELFO_RETARGET.md`.

**Now open (unblocked):** min-fuel GTO→ELFO — re-run `casadi_energy_freetf` from
the energy seed with ε:1→0. Then the indirect ELFO cells.

## Key module map

- `sundman_minfuel/` — direct min-energy backbones + min-fuel (Sundman collocation).
  ELFO: `casadi_energy_freetf.m` (free-t_f, two-primary clock, gravity homotopy) +
  `gen_elfo_energy_gravhom.m` (the 4-leg ladder) → `results/energy_elfo_freetf.mat`.
- `PSR/` — PMP-Steered Refinement (direct min-fuel pipeline) + ELFO retarget work.
- `ztl/` — indirect Sundman multiple-shooting (energy anchor 4.8e-10).
- `min_time/` — min-time (single + multiple shooting); tulip MS validated 4e-9.
- `ms_band/`, `ifs/` — indirect min-fuel band attempts (conditioning wall).
