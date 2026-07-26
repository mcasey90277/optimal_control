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

## Layout

```
direct/
├── run_gto_tulip.m      THE front door (was run_psr.m)
├── run_one.m            function-call route (returns a struct)
├── run_batch.m          multi-factor batch driver
├── run_certified_minfuel.m   reproduce the flagship 1.15× result
├── run_tulip_front.sh, psr_batch.sh, orchestrate/   sweep + batch shells
├── setup_paths.m        the ONE path setup for this campaign
├── lib/                 solver, homotopy, seeds, minfuel_at_tf, refine algorithm
├── certify/             PMP verification, second-order (SSOSC), front aggregation
├── viz/                 movies and frames
├── tests/               guardrail tests
├── results/             solver artifacts + caches (results/psr = pipeline intermediates)
└── data/                exported data products
```

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
