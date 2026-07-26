# GTO_tulip — direct min-time/energy/fuel GTO→tulip solvers

Direct (collocation NLP) solvers for the low-thrust GTO → south-pole tulip
transfer in the Earth–Moon CR3BP (15 kg, 25 mN, Isp 2100 s, ~40-rev spiral).
Min-time reference: t_f = 6.290694 ND = 27.8845 d, ΔV 4.4665 km/s. The
flagship result is the **certified sharp bang-bang min-fuel solution**
(1.15× min-time, 25 switches, defect 2.4e-14, ΔV 3.3696 km/s) and the
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
| `thrustN` ≠ 25 mN | refused, with the explanation | fixed-τ_f topology wall; the 20 mN pilot was an honest failure (`process/LADDER_PREP_PILOT_FINDINGS.md`) |
| `factor` < ~1.12 | refused, lists the backbones on disk | the ε=1 energy backbone itself will not converge approaching min-time |

**Two certified optima at the flagship t_f (measured 2026-07-26).** At
factor 1.150 the energy-backbone route the front door drives and the
`run_certified_minfuel` chain both converge machine-tight to 25 switches at the
*same* t_f — and land on **different** local optima:

| route | m_f | ΔV | propellant |
|---|---|---|---|
| `run_certified_minfuel` (flagship) | 0.849066 | 3.3696 km/s | 2.2640 kg |
| energy backbone → `minfuel_at_tf` | 0.847086 | 3.4176 km/s | 2.2937 kg |
| | **−0.00198** | **+1.43%** | **+0.0297 kg** |

Both are valid extremals; min-fuel here has closely spaced optima and the **seed
route** decides which you land on. `run_gto_tulip` cross-checks every certified
run against `sundman_minfuel_certified.mat` when the t_f matches and says plainly
which basin it found, so a front-door run can no longer quietly report numbers
1.4% off the headline result. **Quote the flagship numbers as the campaign
result.**

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
