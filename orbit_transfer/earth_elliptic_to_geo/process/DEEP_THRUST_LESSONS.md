# Deep-thrust (0.2 N / 0.1 N) — lessons learned + reproduction recipe

**Result (2026-07-20): the FULL 10 → 0.1 N ladder is certified** — the recipe
below reached both deep rungs that were *"never attained / instant wall"* before
the 2026-07-19 external code review (GPT-5.6-terra + Gemini 3.1 Pro):
- **0.2 N** (`../direct/results/MEE_M2_0p2N.mat`): m_f=1377.29 kg, 823 sw, 346.7 rev,
  maxDefect **2.47e-13**, termErr 7.5e-36, incl 0°, `Solve_Succeeded`, edge 99.9%.
- **0.1 N** (`../direct/results/MEE_M2_0p1N.mat`): m_f=1377.29 kg, 1644 sw, 693.6 rev,
  maxDefect **5.04e-13**, termErr 0.00, incl 0°, `Solve_Succeeded`, edge 99.9%
  — reproduced by `reproduce_deep_rung(0.1,'../direct/results/MEE_M2_0p2N.mat')`
  (`maxIter=5000`, all 17 ε-steps ok=1), validating the driver + recipe at the
  last rung. m_f is near-thrust-independent across the whole ladder (10 N=1377.10
  … 0.1 N=1377.29), as the paper's Fig-23 predicts. (Caches gitignored, like all
  campaign `.mat`.)

**Switch-count mesh caveat (P0, 2026-07-21).** The switch counts above (823, 1644)
are 8-node/rev point estimates. The primal mesh-convergence study
`process/P0_SWITCH_MESH_CONVERGENCE.md` refined 0.2 N through 16/24/40 nodes/rev
(all ε=0-certified): **mass converges to 1375.8 kg** (the 8/rev 1377.29 is +0.11%
high) and **revs to 346.7** (mesh-invariant), but the **switch count converges to a
band ~866 ± 5, not 823** — the 8/rev value is a ~5% *undercount* (it does not diverge
with density: 24→40 rev goes down, staying in-band). Report deep-rung switch counts
as bands, not exact integers; 0.1 N (not yet refined) is expected to follow suit. The
external-review concern about 8/rev under-resolution is thereby resolved: the physics
(mass, revs) is trustworthy, only the exact switch integer was mesh-sensitive.

## Reproduce it
```matlab
setup_paths; addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));   % MATLAB R2025b
best = reproduce_deep_rung(0.2, '../direct/results/MEE_M2_0p5N.mat');       % warm from 0.5 N
% best.certified == 1, best.m_f_kg ~ 1377.29, best.maxDefect ~ 1e-13
```
`reproduce_deep_rung.m` is the committed driver; it encodes the four levers below.
0.1 N: `reproduce_deep_rung(0.1, '../direct/results/MEE_M2_0p2N.mat')` (warm-chain from the
0.2 N we just certified). It is a multi-hour run (0.2 N ≈ 34 min of solving over
17 ε-steps at N≈2773; 0.1 N is ~2× the revs). Runs are detached; observe via the
per-ε-step `[adap NN]` lines and the step-cache count (see "Observability" below).

## The four levers (each fixes ONE specific wall)

Every one is opt-in and **inert at feasible points** — the certified 10 N solve
reproduces byte-identically (`m_f=1377.1012`, |Δ|=0) with all of them on.

1. **Rung-adaptive `dL` bound** (`casadi_lt_mee.m`, commit `dbc17e3`).
   A fixed `dL ≤ 2000` (~318 rev) made 0.2 N (ΔL≈2168) and 0.1 N (≈4335)
   **structurally infeasible before any solve**. Now `dLub = max(2000, 5·dL0)`
   from the warm-start C-law estimate. *This was the single most concrete
   blocker* — found by GPT-5.6.

2. **`opts.liftDL`** (`casadi_lt_mee.m`, commit `8a3c78c`, plumbed `9933383`).
   The single scalar `ΔL` multiplies **every** defect ⇒ a DENSE Jacobian column
   ⇒ a bordered/"arrowhead" KKT hostile to sparse MUMPS factorization ⇒ the
   METIS crash at large N. Both reviewers pinned this independently. `liftDL`
   promotes `ΔL` to a tied `N+1` sequence (local bidiagonal chain) ⇒ block-banded
   KKT. Verified: max Jacobian column nnz **1546 → 19** at N=193 (~1700× at 0.2 N).

3. **Phase-correct β warm start** (`warmstart_phase_beta.m`, commit `d10f033`).
   The optimal β oscillates with orbital phase (L mod 2π); interpolating a source
   β sequence against the transfer fraction σ keeps the SOURCE rev-frequency, so
   onto a finer-rev rung it **phase-aliases** (0.5→0.2 N is 2.5× revs). Recompute
   β from the seed's tangential steering law (v̂ in RTN) fresh at each target
   node's own state+longitude — phase-correct, no aliasing. Found by Gemini.

4. **ε-continuation 1→0 with GENEROUS `maxIter`** (`homotopy_mee.m`, `adaptiveEps`,
   commit `ac0cade`). Do NOT attempt a direct ε=0 bang-bang solve — it thrashes
   at feasibility with no switching structure. Continue energy→fuel gently.
   **THE decisive detail: `maxIter` must be large enough that every ε-step
   converges to `ok=1` (defect ~1e-13).** At `maxIter=1500` the steps hit the cap
   under-converged (`ok=0`, best defect ~2.7e-8); that accumulated error
   collapsed the deep-ε tail into an IPOPT restoration spiral (`inf_du→1.8e8`)
   that OOM'd. At `maxIter=3000` every step converged and it marched cleanly to
   ε=0. **Under-iteration, not conditioning, was the final wall.**

## What did NOT work (and why) — do not repeat

- **`scaleNLP` (partial manual constraint scaling) is HARMFUL as implemented.**
  Scaling the O(tf) time-row defect by 1/tf while IPOPT's `gradient-based`
  auto-scaling is still on makes the two scalings FIGHT — the ε=1 energy solve
  went straight into a restoration failure (defect 82). The reviewers meant to
  *replace* auto-scaling with a *complete* user-scaling (all vars + all
  constraints + objective), not stack a partial manual scale on top. The
  `scaleNLP` code exists (`casadi_lt_mee`, opt-in, default off) but should NOT be
  used until a complete user-scaling is built. **OPEN ITEM.**
- **`adaptiveEps` bisection** was armed but **never triggered** at the winning
  config — with enough `maxIter`, no step failed. It is a correct safety net
  (bisect toward the last converged ε instead of propagating a bad iterate), kept
  for robustness, but the real fix was `maxIter`.
- **Direct ε=0 from an aliased warm start**: reaches feasibility (`inf_pr→5e-7`)
  but never resolves the ~350-switch bang-bang structure — churns and stalls.

## Observability lesson (for long detached runs)

MATLAB `-batch` **block-buffers stdout** when redirected to a file, so
per-iteration IPOPT output (and even the per-step summaries) **lag in ~4 KB
chunks** — not real-time. The RELIABLE progress signal is the per-ε-step cache
file (`<resDir>/<tag>_stepNN.mat`), flushed on save. Monitor via
`ls <resDir>/*.mat | wc -l` (converged-step count) and, lagged, the `[adap NN]`
summary lines. Diagnose a stop by: fresh `~/matlab_crash_dump*` (MEX/MUMPS
SIGBUS) vs. a restoration spiral in the log (`inf_du` blowup, `r` iters) vs. OOM
(RSS ballooning + no crash dump).

## The arc (for context)
Never-attained → (dL fix) feasible → (liftDL) 45k-var NLP iterates, no crash,
reaches feasibility → (phase-β + ε-cont, maxIter=1500) marches to ε=0.012 / 99.4%
bang-bang then the under-converged tail OOMs → (maxIter=3000) **every step ok=1 →
ε=0 certified.**
