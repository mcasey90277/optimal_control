# GTO_tulip — direct min-time/energy/fuel GTO→tulip solvers

Direct (collocation NLP) solvers for the low-thrust GTO → south-pole tulip
transfer in the Earth–Moon CR3BP (15 kg, 25 mN, Isp 2100 s, ~40-rev spiral).
Min-time reference: t_f = 6.290694 ND = 27.8845 d, ΔV 4.4665 km/s (to the
max-ẏ tulip point; the min-time to the front's own target is 5.8267 ND — see
**Factor scale** below). The flagship result is the **certified sharp bang-bang
min-fuel solution** at 1.15× (published row: 25 switches, defect 2.4e-14,
ΔV 3.3696 km/s; **best certified row, 2026-07-29: 24 switches, 2.2487 kg,
15.3 g better** — one of at least six optima at that t_f, see below) and the
ΔV-vs-t_f front built around it.

## Goals

1. **Perfect the direct code** — close the remaining direct open items (the
   1.01–1.11× near-min-time band).
2. **Get the indirect code working** — a certified indirect (PMP shooting)
   solve of the same problem; today the indirect campaigns are built but stall
   short of certification.

Concrete open items live in `TODO.md`.

## Three objectives, one pipeline

Min-time / min-energy / min-fuel are **not separate codebases** — they are one
homotopy chain through the same solver (which is why `direct/` is not split by
objective):

| objective | role in the chain | entry point |
|---|---|---|
| min-time | anchor: sets `t_f,min` (throttle ≡ 1 mode of the same core) | `direct/lib/gen_tulip_mintime` (direct); `indirect/min_time/` (PMP shooting) |
| min-energy | homotopy root (the SAME fuel solver at ε=1; smooth, big basin) | `direct/run_gto_tulip` with `epsMin=1`; `direct/lib/gen_tulip_energy_2p` |
| min-fuel | target (ε=0, bang-bang), reached by the ε:1→0 sweep | `direct/run_gto_tulip` (**the front door**, `epsMin=0`); `run_certified_minfuel` |

## Folder map

| where | what |
|---|---|
| `direct/` | **The whole direct campaign, one front door.** `run_gto_tulip` (solve → optional PSR refine → export → verify → movie), `lib/` (solver, homotopy, seeds, refine algorithm), `certify/` (PMP + second-order), `viz/`, `tests/`, `results/`, `data/`. Start at `direct/README.md`. Flattened 2026-07-26 from the old `sundman_minfuel/` + `PSR/` split — see that README for why. |
| `indirect/lowThrust_GTO_tulip/` | Base indirect campaign: PMP shooting w/ complex-step Jacobians, theory note + guided tutorial (`gto_tulip_mintime_theory.pdf`, `building_the_gto_tulip_solvers.pdf`). |
| `indirect/ms_band/` | Indirect multiple-shooting attack on the hard 1.01–1.11× transition band (own campaign doc + unit tests). |
| `indirect/ifs/` | IFS — Indirect Finishing Solve: direct-seeded indirect certification machinery. `IFS_data/` holds its gitignored caches. |
| `indirect/ztl/` | Zhang-thrust-ladder indirect probes (P0 findings recorded in its docs). |
| `indirect/min_time/` | PMP min-time root (always-burn shooting; seeds retargeting for tulip and ELFO). |
| `catalog/` | **Costate-catalog product** (Stage B, 2026-08-29): min-time PMP costates at the Darin-standard propulsion regime (150 kg / 1710 s, few-rev), NOT the 25 mN/~40-rev direct campaign above — a separate, non-periodic-departure (`'gto'` pseudo-family) member of the `costate_common` catalog line alongside DRO/HALO/DPO/L1↔L2-halo → tulip. 16 sheets (4 GTO orientations × 4 Np), 2,625 entries, 73% pair coverage, rungs [15 12 10 7 5] N, conjugate 2,625/0/0. Two drivers: `gto_entry` (single cell), `run_gto_catalog` (swath/regen). Start at `catalog/README.md`. |
| `process/` | Campaign narratives + plans (see Key documents below). |
| `doc/` | Technical notes, figures, briefing, and `doc/reviews/` (external code-review records). |
| `attic/` | Superseded code: fmincon-era min-time/min-fuel NLPs, Sundman prototypes, old continuation experiments. Do not use; see `attic/README.md`. The fmincon min-time formulation notes (density-matched mesh, throttle-on-bound gotcha, mesh-refinement table) are preserved in `attic/README_legacy_fmincon_era.md`. |

Shared problem definition (`cr3bp_lt_params`, `minfuel_config`,
`gto_tulip_endpoints`) lives in `../cr3bp_common/`; every module's
`setup_paths.m` pulls it in via `setup_cr3bp_common()`.

## Key documents

| doc | role |
|---|---|
| `process/LOW_THRUST_MINFUEL_CAMPAIGN.md` | Full campaign record: every method generation, what failed, why, and the winning recipe. Read first. |
| `process/HONEST_EVALUATION_DV_TF_FRONT.md` | Candid assessment: what the certification does/doesn't prove, branch structure of the front, open problems. |
| `process/CODE_CLEANUP_PLAN.md` | This reorganization (phased; Phase 0 done 2026-07-09). |
| `process/MIN_ENERGY_NOTES.md` | Min-energy (homotopy root) derivation notes. |
| **`doc/gto_tulip_guide.tex/.pdf`** | **Start here** — 9-page working guide: transfer problem, OCP, discretization (why Sundman, why fixed $\tau_f$, why the mesh is inherited), pipeline flow diagram, how to run, what is certified, what is not solved. |
| **`doc/run_gto_tulip_explained.tex/.pdf`** | Block-by-block walkthrough of the front-door *script*: each of its seven blocks, why it is written that way, and which lines encode a hard-won decision. Read after the guide. |
| `doc/sundman_minfuel_solution_note.tex/.pdf` | 7-page technical note: OCP, homotopy, Sundman regularization, IPOPT. |

## Entry points

**Start here — the front door** (house `run_gergaud` / `run_cr3bp_geo` pattern):

```matlab
cd direct
setup_paths
edit run_gto_tulip             % set factor / epsMin in section 1, then run
```

One editable parameter block, one solved transfer, artifacts and a figure under
a run name you choose. It drives `minfuel_at_tf` rather than re-implementing the
chain, so every run inherits the seed-fingerprint check, the schedule policy from
`minfuel_config`, the certification gate and the provenance stamp.

`epsMin` selects the objective along the one homotopy chain: `1` = min-energy,
`0` = min-fuel, in between = ε-optimal.

Two knobs are deliberately **refused** rather than silently failing — both are
open campaign problems, not settings:

| request | what happens | why |
|---|---|---|
| `thrustN` ≠ 25 mN | refused, redirected to `direct/run_tulip_ladder` | not a topology wall after all (2026-07-27): the ladder walks 25 → 20 mN in ~4% steps, every rung machine-tight, and sharpens 20 mN to ε=0 (11 sw, 2.104 kg, ΔV 3.113 km/s; certified-quality, full gate set pending). The fixed-τ_f solver stops at 21 mN and the free-time (`cScale`) solver at ~19.5 mN; that ceiling is **not** resolution (190 nodes/rev) and **not** t_f margin (tested) — cause open. `process/LADDER_PREP_PILOT_FINDINGS.md` records the original one-jump failure; `TODO.md` the ladder |
| `factor` < ~1.12 | refused, lists the backbones on disk | the ε=1 energy backbone itself will not converge approaching min-time (1.12× here ≈ 1.21× of the target-consistent min-time; see Factor scale) |

**At least six certified optima at the flagship t_f (2026-07-26 → 07-29).** At
factor 1.150 the energy-backbone route and the `run_certified_minfuel` chain
first showed two machine-tight optima 1.43% apart in ΔV; a 13-seed sweep on the
flagship's own mesh (`process/BASIN_1150_SWEEP.md`) then found **six**, with
24, 25 and 26 switches at the same t_f, and the best of them beats the
published flagship:

| route | m_f | sw | ΔV | propellant |
|---|---|---|---|---|
| **best certified — `direct/lib/sundman_minfuel_basin24_f1150.mat`** | **0.850087** | **24** | **3.3449 km/s** | **2.2487 kg** |
| `run_certified_minfuel` (published flagship) | 0.849066 | 25 | 3.3696 km/s | 2.2640 kg |
| energy backbone → `minfuel_at_tf` | 0.847086 | 25 | 3.4176 km/s | 2.2937 kg |

The 24-switch winner passed `run_foc_tulip` at parity with the flagship (KKT
5.5e-13, primer 4.8e-18, sign law 100%, Ṡ_min 28.3). Distinct seeds stay in
distinct basins; the **seed route** and even the **continuation path** decide
which you land on (third confirmation: the 20 mN rung reached by two paths,
4.4% apart in ΔV). The follow-up front sweep (`process/FRONT_SWEEP.md`) found
the same multiplicity concentrated at the short-t_f end (+53 g at 1.12×, +11 g
at 1.14×, +14 g at 1.30×, +9 g at 1.35×, nothing from 1.40× up; the four
improved rows certified 2026-07-30 as `direct/lib/minfuel_best_f####.mat`), so
the published front is a **lower bound in its left third and accurate in its
right two thirds**. `run_gto_tulip` cross-checks every certified run against
the reference when the t_f matches and says which basin it found.

**Publication rule (2026-09-06):** report the best certified row at each t_f
*together with* the multiplicity (six optima at 1.15×, switch count 24–26); the
multiplicity is itself a result. Do not quote the 25-switch 3.3696 km/s row as
"the" optimum.

**Factor scale (2026-07-15 finding, not yet rebased).** `cfg.tfMin = 6.2907 ND`
is the min-time to a *different* tulip point (max-ẏ,
`indirect/min_time/mintime_params`); the certified min-time to the rendezvous
the front actually targets is **5.8267 ND** (`direct/lib/gen_tulip_mintime`,
hard all-burn, mesh-invariant, 2026-07-15). Every "×" label in this campaign is
therefore against the wrong anchor by a factor 1.0796: the certified front
starts at **1.21×** the target-consistent min-time, the flagship 1.15× is
**1.24×**, and the unsolved "1.01–1.11×" band is **1.09–1.20×**. Physical t_f
and ΔV are unchanged. The config value is deliberately *not* changed here ---
file names, fingerprints and every stored artifact key on it --- so relabel at
presentation time until a coordinated rebase (`TODO.md`).

The lower-level entries remain available:

```matlab
run_certified_minfuel          % reproduce THE certified 1.15x result (~15 min)
minfuel_at_tf(1.30)            % solve one t_f from the energy backbone
aggregate_front                % collect + PMP-verify + honest 3-class front plot
test_minfuel_lib               % cheap no-solve guardrail checks
```

**Sweep the ΔV/t_f front** (the tulip analogue of the other campaigns' thrust
ladders — it sweeps `t_f`, *not* thrust, for the reason in the table above):

```bash
direct/run_tulip_front.sh              # default band, min-fuel
direct/run_tulip_front.sh 1.15 1.20    # just these factors
EPSMIN=1 direct/run_tulip_front.sh     # min-energy sweep
```

Process-isolated per factor, resumable, and a failed factor is recorded rather
than fatal — an uncertified factor is information about the front.

Batch orchestration (process isolation + watchdog + retry — required because
of sporadic uncatchable CasADi/IPOPT MEX crashes):

```bash
direct/orchestrate/backbone_walk.sh 1.15 1.20 1.25 1.30   # energy chain
direct/orchestrate/sharpen_batch.sh 2 1.30 1.35 1.40      # parallel sharpen
```

The folder root deliberately holds only `README.md` and `TODO.md`; campaign
records are in `process/`, technical notes in `doc/` (tidied 2026-07-21).

## Companions

- Indirect (PMP shooting) counterpart + theory note + guided tutorial:
  `indirect/lowThrust_GTO_tulip/`
- Problem source: pumpkynPie `Demos/LunaNet Analysis/lowThrust_GTO_Tulip.m`
