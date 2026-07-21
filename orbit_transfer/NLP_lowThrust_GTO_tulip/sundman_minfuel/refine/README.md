# refine/ — PMP-Steered Refinement (PSR) prototype

**Method name: PMP-Steered Refinement (PSR).** The direct collocation method
is the only solver; the indirect/PMP machinery never solves anything — it
*steers* where the mesh refines. This is "point 4" of the direct↔indirect
roadmap. Its sibling — **Indirect Finishing Solve (IFS)**, "point 3", where an
actual fixed-structure saltation-aware TPBVP takes over as the workhorse to
place S=0 exactly and diagnose switch-count changes — is a separate, unbuilt
method. Keep the distinction sharp: in PSR the indirect side *steers*; in IFS
it *solves*.

Point-4 prototype: adaptive mesh refinement for the Sundman min-fuel
bang-bang solver, driven by how well the collocation mesh **localizes PMP
switching-function zero crossings** (Option 1), with the Sundman-domain
Hamiltonian residual carried alongside as a **passive** diagnostic only
(Option 2 — a Hamiltonian-*driven* refiner — is not built here; this
prototype answers whether it is needed).

Each refinement round: measure the current solution's switch-localization
score (`pmp_refine_indicator`, reusing the KKT-dual costate recovery from
`ms_band/verify_direct_pmp.m`), bisect the worst-localized collocation
intervals (`refine_sigma`), build a no-resample warm start on the refined
mesh (`warmstart_on_mesh` — every original node's state/control is copied
verbatim; only inserted nodes are interpolated), and re-solve the direct
Sundman problem at `epsilon=0, warmTight=true`. The loop stops when switch
times stabilize (max move below the local mesh width, propellant drift below
`propTol`, switch count unchanged), at `maxRounds`, or if a re-solve fails to
converge tight. History is persisted to disk every round, so an IPOPT
MEX-crash (documented, uncatchable) never loses completed rounds.

Limitation: the stabilization check compares quantized bracket-midpoint
switch times (`tauSwitch`), not the sub-cell S=0 root (`diag.tauCr`), so the
acceptance test's resolution is only ~half a local cell width.

## File map

| file | role |
|---|---|
| `pmp_refine_indicator.m` | per-interval refinement score from the PMP switching function S(τ); also returns the passive Hamiltonian residual `Hres`/`HresMax` and switching-law violation count `nViol` |
| `refine_sigma.m` | bisects the top-K worst-scored collocation intervals, preserving every original node (no-resample discipline); guards `hFloor` / `maxAdd` |
| `warmstart_on_mesh.m` | builds the warm-start `(X0,U0)` on the refined mesh — originals verbatim, insertions pchip-interpolated (throttle step-held from the left) |
| `prep_refine_seed.m` | normalizes a certified solution into the seed layout `refine_loop` needs, regenerating KKT-dual costates (`out.lamDef`) via an `eps=0 warmTight` re-solve if the source `.mat` predates dual extraction |
| `refine_loop.m` | the driver: measure → refine → warm-start → re-solve, round by round; saves `refine_history_<tag>.mat` and `refine_<tag>.png` every round |
| `run_headline_1p15.m` | headline demonstration entry point on the certified 1.15× solution |
| `test_*.m` | unit/smoke tests for the above (Tasks 1–5) |

## How to run

Headline 1.15× demonstration (prepares `seed_1p15.mat` from
`../sundman_minfuel_certified.mat` if not already present, then up to 4
refinement rounds):

```bash
cd /Users/msc/Desktop/optimal_control/orbit_transfer/NLP_lowThrust_GTO_tulip/sundman_minfuel/refine
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('$(pwd)'); run_headline_1p15" 2>&1 | grep -v -i "home license\|personal use\|academic, research\|organizational use" | tee headline_run.log
```

Fast smoke test (2 rounds, the 10-switch 1.12× file, ~1 min/round — good for
iterating on the refinement logic itself without the 1.15× problem's longer
re-solve time):

```matlab
opts = struct('maxRounds', 2, 'tag', 'smoke_1p12', 'K', 6, 'maxAdd', 30);
history = refine_loop(fullfile('..', 'results', 'minfuel', 'legacy_ms_f1120.mat'), opts);
```

## Seed options

- **Default (headline): 1.15×** — `prep_refine_seed('../sundman_minfuel_certified.mat', 'seed_1p15.mat')`.
  The certified 1.15× solution was saved before dual extraction existed, so
  prep re-solves once (`eps=0 warmTight`, ~1 min) to regenerate `out.lamDef`
  before the seed can be used.
- **Fast alternative: 1.12×** — `results/minfuel/legacy_ms_f1120.mat`. Already
  carries duals; skips the prep re-solve. Used for the Task-5 smoke test and
  for quick iteration on `refine_loop` itself.

## Results and further reading

- `RESULTS.md` — the headline 1.15× run's actual per-round table, switch-time
  stabilization verdict, propellant drift, and the Option-2 escalation
  decision.
- Design spec: `docs/superpowers/specs/2026-07-11-pmp-mesh-refine-design.md`.
- Campaign record and next steps: `../../LOW_THRUST_MINFUEL_CAMPAIGN.md`.
