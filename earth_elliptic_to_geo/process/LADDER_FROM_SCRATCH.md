# Reproducing the min-fuel thrust ladder FROM SCRATCH

The one-page operational recipe: exact commands, per-rung strategy, expected
numbers, and the gotchas. For the full narrative see
`doc/campaign_reproduction_runbook.pdf`; for the deep-rung story see
`process/DEEP_THRUST_LESSONS.md`; the machine-readable per-rung parameters live
in `reproduce/table3_recipes.m` and the certified numbers in
`reproduce/table3_certified.m`.

## The problem
2-body (Earth-only `1/r²`, no Moon) **minimum-fuel low-thrust GTO→GEO** transfer,
reproducing Haberkorn–Martinon–Gergaud (JGCD 27(6), 2004). MEE / L-domain
collocation (true longitude `L` is the independent variable, time `t` a state,
total span `ΔL` a decision variable). Each rung is a **fixed-time** solve at
`t_f = c_tf · t_{f,min}` with `c_tf = 1.5`, run as a Bertrand–Épénoy
**energy→fuel homotopy** (`ε: 1 → 0`; `ε=1` smooth energy, `ε=0` bang-bang fuel).

## Prerequisites (every session)
```matlab
% MATLAB R2025b ONLY (R2025a license is broken)
cd ~/Desktop/optimal_control/earth_elliptic_to_geo
setup_paths          % adds module paths; the solver auto-adds CasADi from ~/casadi-3.7.0
```

## TL;DR — pick your goal
| Goal | Command |
|---|---|
| **Whole top ladder (10→0.5 N), crash-safe** | `cd reproduce && ./reproduce_table3.sh 10 5 2.5 1 0.5` |
| Whole top ladder, in-process (crash-free rungs only) | `reproduce_table3([10 5 2.5 1 0.5])` |
| **One rung, live solve** | `run_gergaud(struct('thrustN',5,'runMode','solve'))` |
| One rung, from-scratch + verified vs floor | `reproduce_row(5)` |
| **Deep rung 0.2 N** (warm-chain from 0.5 N) | `reproduce_deep_rung(0.2,'results/MEE_M2_0p5N.mat')` |
| **Deep rung 0.1 N** (warm-chain from 0.2 N) | `reproduce_deep_rung(0.1,'results/MEE_M2_0p2N.mat', struct('maxIter',5000))` |
| Campaign continuation builder (all rungs, warm-chained) | `run_ladder([10 5 2.5 1 0.5], struct())` |

## Expected results (the full certified ladder, `c_tf = 1.5`)
| T [N] | m_f [kg] | switches | revs | t_f [ND] | anchor strategy |
|---|---|---|---|---|---|
| 10  | 1377.10   | 19   | 7.33   | 33.33   | `coldB` (cold seed) |
| 5   | 1364.54   | 32   | 14.16  | 67.02   | `chain` from 10 N |
| 2.5 | 1369.79   | 76   | 27.84  | 133.88  | `chain` from 5 N |
| 1   | 1371.44   | 171  | 69.15  | 335.71  | `smallN_first` (low node density first) |
| 0.5 | 1375.28   | 362  | 138.60 | 669.42  | `R0law` (no anchor solve) |
| 0.2 | 1377.29   | 823  | 346.73 | 1673.55 | `R0law` + deep-rung recipe |
| 0.1 | 1377.29   | 1644 | 693.60 | 3347.10 | `R0law` + deep-rung recipe |

`m_f` is near-thrust-independent (~1377 kg) — the paper's Fig-23 result. The
min-time·thrust product is ~constant: `t_{f,min} ≈ 223.14 / T` ND (the "R0 law",
holds to <1% across the anchors), so `t_f` scales as `~1/T`.

## How a rung is built (two stages)
Both live in `drivers/`; `run_gergaud` / `reproduce_row` / `run_ladder` compose them.

**Stage A — min-time anchor `t_{f,min}(T)`** (`run_mintime_mee.m`), strategy per
rung from `table3_recipes.m`:
- `coldB` — cold tangential seed (10 N only).
- `chain` — warm-chain the previous rung's converged min-time trajectory, `ΔL`
  rescaled by the C-law `ΔL(T_new) = ΔL(T_prev)·(T_prev/T_new)` (5, 2.5 N).
- `smallN_first` — solve at low node density first, then refine (1 N; a raw cold
  seed stalls past ~3–4 revs).
- `R0law` — skip the anchor solve, use `t_{f,min} = 223.14/T` (0.5 N and deeper).

**Stage B — fixed-`t_f` fuel solve** (`run_transfer_mee.m`): seed (`mee_seed`) →
guarded `ε:1→0` homotopy (`homotopy_mee`) at `t_f = c_tf·t_{f,min}` → structure
report → certified-only save. Each rung k>1 warm-starts from the previous rung's
converged fuel trajectory (`ΔL` C-law rescaled).

## Deep rungs (0.2, 0.1 N) — the four levers
`reproduce_deep_rung.m` encodes the recipe that first cracked these (external
review 2026-07-19; full story in `process/DEEP_THRUST_LESSONS.md`). Each lever
is opt-in and **inert at feasible points** (10 N reproduces `m_f=1377.1012`, |Δ|=0):
1. **Rung-adaptive `dL` bound** — a fixed `dL≤2000` made them infeasible (0.2 N
   needs ΔL≈2168, 0.1 N ≈4335); now `max(2000, 5·dL0)`.
2. **`liftDL`** — lift scalar `ΔL` to a tied `N+1` sequence → block-banded KKT
   (the scalar was a dense Jacobian column → MUMPS/METIS crash at large N).
3. **Phase-correct β warm start** (`warmstart_phase_beta.m`) — recompute β from the
   tangential law at each node (a σ-interp aliases across rungs).
4. **ε-continuation with generous `maxIter`** — `maxIter=3000` (0.2 N) / `5000`
   (0.1 N); **under-iteration collapses the deep-ε tail** (this, not conditioning,
   was the final wall). `scaleNLP` is a KNOWN TRAP — leave it OFF.

## Gotchas (read before a long run)
- **Uncatchable MEX/CasADi crash at scale** — a `try/catch` inside MATLAB cannot
  catch it; it takes the whole process down. Use `reproduce_table3.sh` (one OS
  process per rung, relaunch-on-death) for the deep/hang-prone rungs. Per-step
  caches make a relaunch resume near where it died.
- **`-batch` block-buffers stdout** (~4 KB chunks) — the reliable live signal is
  the per-ε-step cache count (`ls <resDir>/*.mat | wc -l`) or the `[adap NN]` lines.
- **`c_tf` sits on a flat part of the fuel-vs-time curve**, so `m_f` is nearly
  `t_f`-insensitive — matching mass confirms the convention, not a 3-digit `c_tf`.
- **R2025b only.** `.mat` campaign caches are gitignored (never committed).

## Verify a solved rung
```matlab
verify_row(5)                                   % vs the table3_certified floor (one-sided)
hamiltonian_const_check('results/MEE_M2_5N.mat')% first-order PMP: H_t conserved (CoV ~1e-8)
verify_sosc_mee('results/MEE_M2_5N.mat')        % second-order (WEAK_MIN expected for bang-bang)
```

## Deeper references
- `doc/campaign_reproduction_runbook.pdf` — full whitepaper (anchors, recipes, numbers).
- `process/DEEP_THRUST_LESSONS.md` — deep-rung recipe + lessons + the arc.
- `process/DESIGN_thrust_ladder.md`, `reproduce/table3_recipes.m` — per-rung parameters.
- `reproduce/table3_certified.m` — the canonical certified numbers (source of truth).
