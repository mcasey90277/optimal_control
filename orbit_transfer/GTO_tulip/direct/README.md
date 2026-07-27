# GTO_tulip/direct — one campaign, one front door

Direct (collocation NLP) solution of the low-thrust **GTO → south-pole tulip**
transfer in the Earth–Moon CR3BP: 15 kg, 25 mN, Isp 2100 s, ~40-rev spiral,
Sundman-regularized mesh `dt/dτ = r₁^1.5`.

## Start here

```matlab
cd GTO_tulip/direct
setup_paths
edit run_gto_tulip      % set the parameters in section 1, then run
```

`run_gto_tulip` is **the** entry point. One parameter block, six stages, only
the first mandatory:

| stage | what | optional? |
|---|---|---|
| 1 | parameters | — |
| 2 | **direct solve** — energy→fuel homotopy (ε:1→0) via `minfuel_at_tf` | always runs |
| 3 | **PSR refine** — PMP-steered mesh refinement, sharpens switch times | yes |
| 4 | export data products (also a ready-made IFS seed) | yes |
| 5 | first-order PMP verification | yes |
| 6 | control movie | yes |

`epsMin` picks the objective along the one homotopy chain: `1` = min-energy
(smooth), `0` = min-fuel (bang-bang, the campaign objective), in between =
ε-optimal. The schedule is truncated to end exactly at `epsMin`, so a partial
run is a genuine prefix of the certified step sequence.

Sweep the ΔV/t_f front (this campaign's ladder — it sweeps **t_f, not thrust**,
because thrust is an open problem here; see below):

```bash
./run_tulip_front.sh              # default band, min-fuel
./run_tulip_front.sh 1.15 1.20    # just these factors
EPSMIN=1 ./run_tulip_front.sh     # min-energy sweep
```

## Verifying the campaign

Fast (seconds each, no solves) — run these after any structural change:

```matlab
cd GTO_tulip/direct; setup_paths
test_artifact_paths      % every literal .mat path in the campaign resolves
test_run_gto_tulip       % front-door guards + ONE copy of each shared name
test_minfuel_lib         % library guardrails
test_ladder_prep_tulip   % fingerprint / chain-helper rules
test_refine_sigma, test_warmstart_on_mesh, test_pmp_refine_indicator
```

Slow, but the one that actually proves the campaign still works:

```matlab
best = run_certified_minfuel(1500, '/tmp/repro_check.mat');   % ~30 min
% expect m_f ~0.8491, dV ~3.365 km/s, 2.26 kg, defect ~1e-14, certified=1
```

Pass an explicit `saveFile` as above: with no second argument it **overwrites
the published reference**. And judge it on mass/ΔV, not the switch count — a
re-solve is not bit-identical to the published artifact. Three measured
re-solves (2026-07-21, 2026-07-26 ×1) all land on **24** switches vs the
published 25, with mass agreeing to ~0.1%. The switch integer is basin-sensitive
even at fixed mesh; TODO C3 tracks turning it into a published band.

**Run the slow one after any move or path change.** The 2026-07-26 flatten
passed every fast test while `run_certified_minfuel` was broken — its seed load
still pointed at the pre-flatten location, and nothing else loads that file.
`test_artifact_paths` was written to close exactly that gap and now catches this
class in seconds, but it only checks *literal* paths; the reproduction is still
the ground truth.

## Layout

```
direct/
├── run_gto_tulip.m      THE front door (was PSR/run_psr.m)
├── run_one.m            function-call route (returns a struct)
├── run_batch.m          multi-factor batch driver
├── run_certified_minfuel.m   reproduce the flagship 1.15× result
├── run_tulip_front.sh, run_batch.sh, orchestrate/   sweep + batch shells
├── setup_paths.m        the ONE path setup for this campaign
├── lib/                 solver, homotopy, seeds, minfuel_at_tf, refine algorithm
├── certify/             PMP verification, second-order (SSOSC), front aggregation
├── viz/                 movies and frames
├── tests/               guardrail tests
├── results/             solver artifacts + caches (results/psr = pipeline intermediates)
└── data/                exported data products
```

## The mesh is inherited, not chosen — read before touching discretization

**There is no mesh parameter anywhere in this pipeline, and that is structural,
not an oversight.**

`N`, the normalized grid `sigma`, and the total regularized length `tauf0` are
all read off the energy backbone at the requested `t_f`
(`results/energy/energy_f####.mat`). Every downstream solve — the tight
re-clean, all 13 ε steps, any neighbour-seeded continuation — inherits them
unchanged. The flagship runs at **N = 4001 because that is what the backbone
has**, not because 4001 was chosen or justified.

Node placement follows from the Sundman law. With `dt/dτ = κ = r₁^1.5` and a
uniform σ grid,

```
Δt_k  ≈  (tauf0 / N) · r₁(x_k)^1.5
```

so density in time goes as `r₁^-1.5` — dense at perigee, sparse at apogee. The
mesh therefore adapts to the trajectory's **geometry** and to nothing else. It
does **not** adapt to the solution: the switching structure, where the control
is discontinuous and accuracy is hardest, has no influence on node placement.

**Why it is frozen — the no-resample discipline.** Interpolating a converged
solution onto a different mesh reintroduces a ~1e-2 defect floor that pins IPOPT
in restoration and never clears. The recipe is to map a collocation-feasible
solution into τ using its *own* nodes, no interpolation. That is what makes the
ε homotopy work, and its consequence is that `N` cannot be varied downstream:
changing `N` means resampling, and resampling is the thing that was outlawed.

**What this costs us.** *"Is this mesh adequate?"* **cannot be answered from
inside the pipeline** — there is nothing to vary. Answering it means generating
backbones at several densities with `gen_tulip_energy_2p` and comparing across
them (observed order of accuracy; do m_f and switch times converge). That study
is specified in `docs/superpowers/plans/2026-07-25-mesh-convergence-study.md`
and **has not been run**. The honest position today: N = 4001 is *sufficient in
practice* and *unjustified in theory*.

**The one mechanism that changes the mesh** is PSR refinement (stage 3), and it
is one-directional by design: it finds where the PMP switching function
localizes a switch worst, **adds** nodes there, and re-solves — always from the
backbone, never coarsening, never resampling existing nodes. It sharpens switch
times below the original mesh width. It does not establish that the underlying
mesh was adequate.

Full mathematical treatment: `../doc/gto_tulip_guide.pdf` §3.4.

## Why it looks like this (2026-07-26 flatten)

It used to be `sundman_minfuel/` and `PSR/` side by side, which read as two
competing pipelines. It was never two pipelines:

- `run_psr` **stage 2 has always called `minfuel_at_tf`** — the same driver the
  other entry points used. PSR was the front door all along; it was just named
  after its optional stage 3 and filed as if it were a peer of the library.
- `PSR/lib/` vendored **20 files** so PSR could be self-contained. That claim
  stopped being true on 2026-07-15, when an unrelated feature added an
  `addpath` that put the upstream folder *ahead* of the frozen copies —
  measured: PSR had been running the upstream solver for eleven days while its
  README said otherwise, and two vendored files were dead code. The vendoring
  is now dissolved; there is exactly one copy of everything, and
  `tests/test_run_gto_tulip.m` asserts it stays that way.
- The **PMP verifier** (`verify_direct_pmp`, `sms_*`) genuinely lives in
  `../indirect/ms_band/`. That is a real dependency — the direct pipeline
  checks itself with the independent indirect instrument — not duplication, so
  `setup_paths` names it instead of copying it.

## Two knobs that are open problems, not settings

| request | what happens | why |
|---|---|---|
| `factor` < ~1.12 | refused, lists the backbones on disk | the ε=1 energy backbone itself will not converge approaching min-time |
| off-nominal thrust | not exposed | fixed-τ_f topology wall; the 20 mN pilot (`lib/pilot_rung_20mN.m`) was an honest failure — `../process/LADDER_PREP_PILOT_FINDINGS.md` |

## Two certified optima at the flagship t_f

At factor 1.150 the energy-backbone route and the `run_certified_minfuel` chain
both converge machine-tight to 25 switches at the **same** t_f — and land on
**different** local optima:

| route | m_f | ΔV | propellant |
|---|---|---|---|
| `run_certified_minfuel` | 0.849066 | 3.3696 km/s | 2.2640 kg |
| energy backbone → `minfuel_at_tf` | 0.847086 | 3.4176 km/s | 2.2937 kg |

The seed route decides the basin. Stage 2 prints a basin check against the
certified reference whenever the t_f matches, so this cannot pass unnoticed.
Quote the `run_certified_minfuel` numbers as the campaign result. Open item in
`../TODO.md`.
