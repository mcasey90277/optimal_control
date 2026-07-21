# IFS — Indirect Finishing Solve

**Status:** machinery validated (unit tests all pass); full 1.12× convergence
still **OPEN**. Two build phases so far:
- **`RESULTS.md`** — the original single-shot attempt (cold direct seed, lsqnonlin):
  characterized the terminal-cluster / cold-seed shooting conditioning wall.
- **`RESULTS_RUNG01_RUNG2.md`** — the Rung 0+1+2 campaign (this repo's newer work):
  a scaled truncated-SVD solver (`ifs_solve2`) that *descends* the cold seed
  (1.96→~0.43) but doesn't converge; a working **min-time k=0 anchor**
  (`ifs_seed_mintime`, IFS's first end-to-end rendezvous convergence); and
  t_f-continuation (`ifs_tf_continuation` naive, `ifs_tf_arclength`
  pseudo-arclength) that reached a **min-time fold + costate-gauge degeneracy**.
  A GPT-5.6-sol consult (2026-07-12, `CONSULT_GPT56_{prompt,response}.md`)
  reframed the wall: the "gauge" is really a **two-BVP manifold separation**
  (gauge pinning retired), and the true bottleneck is **seed quality**. Next
  lever = a **backward adjoint sweep** seed (replace the KKT-dual map), then a
  min-energy anchor if a t_f gap remains.

## What this is

IFS is "point 3" of the direct↔indirect roadmap: take a good direct min-fuel
bang-bang solution (a PSR-refined direct-collocation result) and hand it to an
**actual indirect solver** that holds the switch *structure* fixed and places
every switch exactly at `S(τ)=0`, producing exact sub-mesh switch times, exact
costates, and a continuous-time first-order PMP certificate. It is the sibling
of the already-built **PSR = PMP-Steered Refinement**
(`../sundman_minfuel/refine/`, see that folder's `README.md`/`RESULTS.md`),
which upgrades the direct mesh but stays in decision-vector space; IFS goes
the rest of the way to an indirect (costate-carrying) solution.

## Method (one paragraph)

Hard-throttle Sundman-τ multiple shooting: each switch is an explicit shooting
node (not an event detected mid-arc), so there is **no saltation matrix** —
a switch's sensitivity is just the ordinary endpoint sensitivity of its two
neighboring arcs. Each arc runs at a *known constant* throttle (u=1 burn / u=0
coast, fixed by the direct solution's structure), so there is no smoothing
parameter ε and no 1/ε layer anywhere — this is what removes the crawl that
killed every prior smoothed indirect attempt (`ms_band`). Switch times are
parameterized by a **monotone stick-breaking** map (`τ_0 < τ_1 < … < τ_k <
τ_f` by construction) so an LM step can never reorder switches or push one
past the fixed terminal time. The residual is solved by **Levenberg-Marquardt**
(`lsqnonlin`) with a **per-arc complex-step sparse Jacobian** (each residual
block depends only on its own node/switch-time unknowns). Two terminal-BC
modes, both keeping the system square: `"rendezvous"` (the real problem — r,v
rendezvous + λ_m(τ_f)=0 transversality + fixed t_f) and `"fixedState"` (an
interior window with the end state pinned, used only as an early ground-truth
test — later found rank-deficient in the λ_m gauge and dropped as a gate, see
`RESULTS.md` §3a).

## File map

| file | role |
|---|---|
| `ifs_eom.m` | hard-throttle 16-dim Sundman PMP EOM (matches `sms_eom` with the entropy term dropped; throttle u is a fixed arc parameter) |
| `ifs_pack.m` / `ifs_unpack.m` | Z ⇄ (λ_0, {N_i node states}, {τ_i switch times}) |
| `ifs_residual.m` | R(Z) multiple-shooting residual (continuity + terminal + switch blocks) and its complex-step Jacobian; contains `ifs_int_arc`, the per-arc integrator with the zero/reversed-span guard |
| `ifs_taus.m` / `ifs_gseed.m` | switch-time parameterization (τ ⇄ unknowns) and its seed inverse; `mode='sigmoid'` (default, monotone stick-breaking) or `'direct'` (times are the unknowns, monotonicity by solver projection) |
| `ifs_seed.m` | builds Z from a direct/PSR `.mat`: switch times from `diag.tauCr`, node states by interpolation, costates via the mode-'d' dual→costate map (`sms_seed_duals`), arc throttles from the direct solution's sign pattern; `'full'`/`'window'` and `opts.tauParam` |
| `ifs_solve.m` | original LM driver (`lsqnonlin`, `Algorithm='levenberg-marquardt'`, `ScaleProblem='jacobian'`, complex-step Jacobian supplied) |
| `ifs_solve2.m` | **Rung-1 solver**: two-sided Jacobian equilibration + truncated/rank-revealing SVD Gauss-Newton step + Levenberg fallback + α-floor line search + adaptive truncation continuation. Same API as `ifs_solve`. Descends the cold seed where `ifs_solve` crawled/froze (see `RESULTS_RUNG01_RUNG2.md`). |
| `ifs_seed_mintime.m` | **k=0 all-burn min-time anchor** for continuation: min-time costates (`run_gto_tulip_indirect`) + `λ_T0=-Ht(0)`, `tauf` by integrating the hard burn to `t=tf`. IFS's first end-to-end rendezvous convergence. |
| `ifs_tf_continuation.m` | Rung-2 driver: naive t_f-march from the anchor tracking `max S`. Documents that naive stepping cannot cross the min-time fold. |
| `ifs_tf_arclength.m` | Rung-2b: pseudo-arclength (Keller) continuation across the fold — predictor/corrector on `[R; arclength]`. Mechanically correct; exposes the min-time costate-gauge degeneracy (see file header + results doc). |
| `run_combined_increment.m` | runs `ifs_solve2` on the gate(s) in both parameterizations vs the old lsqnonlin crawl. |
| `mint_easy_gate.m` | (finding) attempt to mint a 3-switch easy gate — the 3-switch local optimum is not a fixed point of the `eps=0` re-solve, so it can't be minted this way. |
| `ifs_certify.m` | post-hoc certificate: S=0 at each switch, per-arc sign-law sampling (re-integrates each arc — now guarded against zero-length spans), terminal residual, rendezvous transversality; reports (not acts on) vanishing-arc / missing-switch structure diagnostics |
| `run_ifs_1p12.m` | original 1.12× gate entry point (lsqnonlin), k≈10 switches, seeded from `../sundman_minfuel/results/minfuel/legacy_ms_f1120.mat`; saves the solve before certifying and tolerates a certify-time error |
| `test_ifs_eom.m`, `test_ifs_residual.m`, `test_ifs_jacobian.m`, `test_ifs_seed.m`, `test_ifs_solve.m` | unit tests (see `RESULTS.md` §2 for what each one certifies) |
| `setup_paths.m` | path setup (name-collides with `sundman_minfuel/setup_paths.m` — benign, callers `addpath(ifs)` first) |

## How to run

```matlab
cd NLP_lowThrust_GTO_tulip/ifs
setup_paths();
test_ifs_residual; test_ifs_jacobian; test_ifs_seed   % fast, seconds

% original gate (lsqnonlin) — characterizes the conditioning wall
res = run_ifs_1p12();                                  % ~30 min, 400 LM iters

% Rung-1 solver on the gate(s), both parameterizations (see RESULTS_RUNG01_RUNG2)
res = run_combined_increment();

% Rung-2 continuation from the min-time anchor
res = ifs_tf_continuation(1.00, 1.15, 0.01);           % naive: stalls at the fold
res = ifs_tf_arclength(1.00, 1.15, 0.10);              % pseudo-arclength: gauge-limited
```

`run_ifs_1p12` prints a `RUNG 2 SUMMARY` line and saves
`ifs_1p12_results.mat` (the solve `out`, plus `cert` if certify succeeded).

## Seed sources

- **Full 1.12× (Rung 2, the gate):** `../sundman_minfuel/results/minfuel/legacy_ms_f1120.mat`
  — k≈10 switches, already carries `out.lamDef` (zero prep needed).
- **25-switch 1.15× headline (Rung 3, not yet attempted):**
  `../sundman_minfuel/sundman_minfuel_certified.mat` via
  `prep_refine_seed.m` — k=25, 433 unknowns.
- **Interior window (dropped as a gate, kept as a rank-deficiency finding):**
  built in-line by `ifs_seed(matFile, struct('mode','window', ...))`.

## Status line

Machinery validated (EOM/residual/Jacobian/seed unit tests all pass). Full
1.12× convergence remains **OPEN**. Progress since the original attempt:
`ifs_solve2` (scaled truncated-SVD + adaptive truncation continuation)
*descends* the cold 1.12× seed (1.96 → ~0.43) but does not converge — the
cold-seed 40-rev shooting basin wall (the clean-band test showed this floor is
uniform across 1.12×/1.14×/1.25×, i.e. **seed quality, not the terminal
cluster**). Rung 2 built a working **min-time k=0 anchor** (IFS's first real
rendezvous convergence) and t_f-continuation; the min-time point is a **fold**
that naive stepping can't cross and pseudo-arclength wanders. A **GPT-5.6-sol
consult (2026-07-12)** reframed the "gauge": the near-null direction is a
**two-BVP manifold separation** (min-time H=0 vs min-fuel λ_m=0 at t_f,min),
*not* a scale gauge — so ‖λ0‖=1 pinning is retired. **Rung A (backward adjoint
sweep + adjoint smoother, `ifs_seed_adjoint.m`) was built and FALSIFIED
2026-07-12**: the adjoint flow amplifies 1.2e12 backward over 40 revs, and the
mesh-accuracy dual costates are dynamically inconsistent at that level — no
single 40-rev costate trajectory fits them (smoother-blended seed 2.13 vs dual
map 1.96; switching structure perfect at q=1.001±0.008 but arc-continuity rows
sit at data accuracy). Remaining levers: **Rung B — min-energy anchor +
s-homotopy with early structure lock** (compose `ms_band/sms_*` with IFS), and
bounded Rung-D substrate gains (PSR mesh / tighter direct tolerance). Drop the
min-time anchor. Full arc: `RESULTS_RUNG01_RUNG2.md`; plan:
`PLAN_OF_ATTACK_2.md`; consult: `CONSULT_GPT56_response.md`.

## Pointers

- **Current plan: `PLAN_OF_ATTACK_2.md`** (2026-07-12; supersedes
  `PLAN_OF_ATTACK.md`'s rung sequence; Rung A executed and falsified)
- **Rung B execution plan: `PLAN_RUNG_B.md`** (2026-07-12; self-contained
  handoff — smoothed-indirect anchor from energy duals, guarded ε-march to
  ε_lock, hard-throttle IFS finish; increments B0–B4 with gates)
- Design spec: `../../docs/superpowers/specs/2026-07-11-ifs-design.md`
- Full arc / diagnosed failure modes / gate numbers: `RESULTS.md` (this folder)
- Campaign context (direct-side record, two-walls analysis this build's walls
  echo): `../LOW_THRUST_MINFUEL_CAMPAIGN.md`
- PSR sibling (the direct-side mesh-refinement stage IFS follows):
  `../sundman_minfuel/refine/` (`README.md`, `RESULTS.md`)
- Diagnostics referenced in `RESULTS.md`: `../../.superpowers/sdd/ifs-diag-floor-report.md`,
  `ifs-diag-fullrank-report.md`, `ifs-diag-stall-report.md`
