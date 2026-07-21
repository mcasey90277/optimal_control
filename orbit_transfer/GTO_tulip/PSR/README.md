# PSR — PMP-Steered Refinement pipeline (entry point)

**This folder is the front door to the working direct-method pipeline** for
the CR3BP min-fuel GTO → south-pole-tulip transfer (15 kg, 25 mN, Isp 2100 s,
~40-rev spiral). One script runs the whole chain:

```
direct solve (energy→fuel homotopy)  →  PSR mesh refinement (sharp switch
times)  →  costate recovery + data export (../PSR_data/)  →  first-order
PMP certificate  →  control movie
```

## How to run

```matlab
cd GTO_tulip/PSR
edit run_psr          % set factor + seed + knobs in section 1
run_psr
```

Everything a user sets lives in **section 1** of `run_psr.m`: the transfer
time (`factor` = t_f / t_f_min), the seed policy, the homotopy schedule, the
refinement / verification knobs, and the movie mode. Outputs land in
`PSR/results/` under canonical names (`psr_direct_f1150.mat`,
`psr_refined_f1150.mat`, `refine_history_psr_f1150.mat`,
`verify_pmp_psr_refined_f1150.{mat,png}`, `psr_movie_f1150.*`). Stages skip
themselves when their output exists (resumable); `rerun*` flags force a redo.

## Seed policy (stage 2), spelled out

| seed | meaning | when |
|---|---|---|
| `'energy'` (default) | min-ENERGY backbone at this factor (`../sundman_minfuel/results/energy/energy_f####.mat`; exists for 1.12–1.15 by 0.01 and 1.20–1.85 by 0.05) — the canonical homotopy root | fresh solve at a new t_f |
| `'neighbor'` + `seedFactor` | an existing bang-bang solution at a nearby factor, time-rescaled + lightly re-sharpened | walking the dV-vs-t_f front |
| explicit path | any .mat with X/U (+ sigma, tauf0, rv0, rvf) | reruns, experiments |

Missing energy backbone → `../sundman_minfuel/orchestrate/backbone_walk.sh`.

## What each stage is (one line each)

1. **Direct solve** — `minfuel_at_tf`: CasADi+IPOPT trapezoid collocation on
   the Sundman mesh; Bertrand–Épénoy homotopy `J(ε) = ∫s dt − ε∫s(1−s) dt`
   marched ε: 0.6 → 0, deforming the smooth energy solution continuously into
   the bang-bang fuel solution.
2. **PSR refinement** — `refine_loop`: costates from the NLP's own KKT duals
   (adjudicated mode-'d' midpoint map) form the switching function S(τ);
   the mesh is refined where S localizes a switch worst; warm-started re-solve;
   repeat until switch times stabilize below a local mesh width. The indirect
   machinery *steers*, the direct solver *solves*.
3. **Costates + data export** — `psr_export_data`: PSR's NLP yields only raw
   interval duals; this stage runs the adjudicated mode-'d' dual→costate map
   (+β fit) and writes one self-contained product file to **`../PSR_data/`**
   (see below).
4. **Verification** — `verify_direct_pmp` (vendored in `lib/`): per-arc
   propagation of the full 16-dim state+costate system from the solution's own
   duals; primer alignment, switch-structure match, |H_t+λ_t| stationarity,
   λ_m(τ_f)=0 transversality. Appends its summary to the PSR_data file.
   Known-benign flags (terminal switch cluster, near-graze switches) are
   **adjudicated** via `verifyOpts.adjArcs`/`adjSwitches` in section 1 — the
   issued certificate lists them, so nothing is hidden. **First-order
   certificate = extremality only.** A second-order test (Jacobi /
   conjugate-point, or NLP reduced-Hessian SSC) is the planned upgrade to
   claim local minimality — currently NOT proven.
5. **Movie** — `psr_movie`: rotating-frame transfer colored burn/coast with
   primer thrust arrows + synced throttle strip + running ΔV meter
   (generalized from `../movie/animate_sundman_minfuel.m`).

## Data products (`../PSR_data/`)

Stage 4 writes `psr_data_tf<factor>_sw<k>.mat` (e.g. `psr_data_tf1p150_sw25.mat`;
factor with `.`→`p`, k = certified dual-S switch count). One file, two layers:

- **Seed-compatible top level** (`out, sigma, tauf0, rv0, rvf, factor`) — the
  exact layout `ifs_seed` / `verify_direct_pmp` / `sms_seed_duals` consume, so
  the file is a **ready-made IFS seed** with zero conversion.
- **Unpacked products** for independent analysis: `mesh` (σ, Sundman τ,
  physical t at nodes), `traj` (r, v, m, full X), `ctrl` (α, throttle s,
  switch times by certified dual-S crossings AND raw throttle crossings),
  `costate` (λ [8×nN], switching function S, β + spread quality diagnostic —
  note these are O(h) mesh-accuracy costates, see `../ifs/RESULTS_RUNG01_RUNG2.md`
  Rung A), `pmp` (transversality λ_m(τ_f), terminal rendezvous/time residuals,
  S-sign-law agreement, primer alignment), `scal` (dV, prop, m_f, defect),
  `const` (all physical constants — the file is self-contained), `provenance`
  (source, date, git hash, dual-map settings), and `verify` (certificate
  summary, appended by stage 5).

## Design note

**PSR is self-contained (as of 2026-07-12).** The entry drivers (`run_psr`,
`psr_export_data`, `psr_movie`) live at the top level; the 19 machinery files
they need are **vendored into `PSR/lib/`** (copies of `ms_band` +
`sundman_minfuel` machinery — see `lib/README.md` for the manifest, provenance,
and drift caveat). `setup_paths` adds only `PSR`, `PSR/lib`, and the external
`pumpkyn` toolbox. Verified with `requiredFilesAndProducts`: no PSR file reaches
`ms_band/` or `sundman_minfuel/`.

The originals were **copied, not moved**, so the IFS folder and the campaign
scripts keep working. The one input DATA dependency is **referenced in place**:
the min-energy backbones under `sundman_minfuel/results/energy` (the `lib`
copy of `minfuel_config` points there). Two vendored files are PSR-owned and
edited here going forward: `refine_loop` (carries the `outDir`/`solFile`
additions) and `verify_direct_pmp` (adjudication).

## Costs (this machine)

| stage | time |
|---|---|
| direct solve | ~30–90 min (13 IPOPT solves, N=4001) |
| refinement | ~20–60 min (≤4 re-solves + indicator) |
| verification | ~5–10 min |
| movie | seconds (`preview`) / ~15 min (`movie`) |

## Pointers

- Campaign record + open problems (1.01–1.11× band): `../LOW_THRUST_MINFUEL_CAMPAIGN.md`
- PSR design + headline results: `../sundman_minfuel/refine/{README,RESULTS}.md`
- Verifier provenance + dual-map adjudication: `../ms_band/MS_BAND_CAMPAIGN.md`
- Indirect finishing (IFS, research-grade, open): `../ifs/`
- Honest-evaluation notes on the dV-t_f front: `../HONEST_EVALUATION_DV_TF_FRONT.md`
