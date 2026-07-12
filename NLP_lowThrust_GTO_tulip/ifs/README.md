# IFS — Indirect Finishing Solve

**Status:** machinery validated (unit tests all pass); full 1.12× convergence
**blocked** by terminal-cluster shooting conditioning; next lever = tighter
multiple shooting inside the sensitive arcs. See `RESULTS.md` for the full
arc.

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
| `ifs_taus.m` / `ifs_gseed.m` | monotone stick-breaking switch-time parameterization (τ ⇄ unconstrained g) and its seed inverse |
| `ifs_seed.m` | builds Z from a direct/PSR `.mat`: switch times from `diag.tauCr`, node states by interpolation, costates via the mode-'d' dual→costate map (`sms_seed_duals`), arc throttles from the direct solution's sign pattern; builds either the `'full'` problem or a `'window'` sub-arc |
| `ifs_solve.m` | LM driver (`lsqnonlin`, `Algorithm='levenberg-marquardt'`, `ScaleProblem='jacobian'`, complex-step Jacobian supplied) |
| `ifs_certify.m` | post-hoc certificate: S=0 at each switch, per-arc sign-law sampling (re-integrates each arc — now guarded against zero-length spans), terminal residual, rendezvous transversality; reports (not acts on) vanishing-arc / missing-switch structure diagnostics |
| `run_ifs_1p12.m` | Rung-2 gate entry point: full 1.12×, k≈10 switches, seeded from `../sundman_minfuel/results/minfuel/legacy_ms_f1120.mat`; saves the solve before certifying and tolerates a certify-time error |
| `test_ifs_eom.m`, `test_ifs_residual.m`, `test_ifs_jacobian.m`, `test_ifs_seed.m`, `test_ifs_solve.m` | unit tests (see `RESULTS.md` §2 for what each one certifies) |
| `setup_paths.m` | path setup (name-collides with `sundman_minfuel/setup_paths.m` — benign, callers `addpath(ifs)` first) |

## How to run

```matlab
cd NLP_lowThrust_GTO_tulip/ifs
setup_paths();
test_ifs_residual; test_ifs_jacobian; test_ifs_seed   % fast, seconds
res = run_ifs_1p12();                                  % the gate — ~30 minutes, 400 LM iters
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
1.12× convergence is **blocked** by terminal-cluster shooting conditioning
(cond(J)~5.9e9, GN-consistent, 6 orders short of the 1e-8 target after 400 LM
iterations) — a characterized negative result, not a code defect. Next =
tighter multiple shooting (interior shooting nodes inside the sensitive
arcs, especially the terminal switch cluster).

## Pointers

- Design spec: `../../docs/superpowers/specs/2026-07-11-ifs-design.md`
- Full arc / diagnosed failure modes / gate numbers: `RESULTS.md` (this folder)
- Campaign context (direct-side record, two-walls analysis this build's walls
  echo): `../LOW_THRUST_MINFUEL_CAMPAIGN.md`
- PSR sibling (the direct-side mesh-refinement stage IFS follows):
  `../sundman_minfuel/refine/` (`README.md`, `RESULTS.md`)
- Diagnostics referenced in `RESULTS.md`: `../../.superpowers/sdd/ifs-diag-floor-report.md`,
  `ifs-diag-fullrank-report.md`, `ifs-diag-stall-report.md`
