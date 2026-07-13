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
| **min-energy** (var t_f) | direct | ✅ energy backbones factor 1.12–1.95 (`sundman_minfuel/results/energy`) | ⬜ target-homotopy stalls s=0.45 (Moon-ward) |
| | indirect | ✅ Sundman-MS 75 mN anchor **4.8e-10** (`ztl/results/z1_sun_anchor_75mN.mat`); band via costates 🟡 | ⬜ (blocked with direct) |
| **min-fuel** (var t_f) | direct | ✅ PSR pipeline, 3- & 25-switch bang-bang certified, band [1.12,1.95] | ⬜ **blocked on the ELFO energy seed** |
| | indirect | 🟡 IFS/ms_band: 1.12x = 10 switches certified; band = conditioning wall | ⬜ (blocked) |

## The one blocker for the whole ELFO column

Everything ELFO-min-fuel needs is a **GTO→ELFO min-ENERGY seed** (the homotopy
root the PSR fuel pipeline consumes; the pipeline is target-agnostic, so that
seed is the only missing input). Making that seed is the active problem — see
`PSR/ELFO_RETARGET.md` and `min_time/README.md` for the routes tried and their
walls. Leading candidate: a **direct min-time collocation** (free t_f, built on
`casadi_minfuel_sundman`) — the only approach with both direct-collocation
robustness and floating t_f.

## Key module map

- `sundman_minfuel/` — direct min-energy backbones + min-fuel (Sundman collocation).
- `PSR/` — PMP-Steered Refinement (direct min-fuel pipeline) + ELFO retarget work.
- `ztl/` — indirect Sundman multiple-shooting (energy anchor 4.8e-10).
- `min_time/` — min-time (single + multiple shooting); tulip MS validated 4e-9.
- `ms_band/`, `ifs/` — indirect min-fuel band attempts (conditioning wall).
