# GTO_tulip — direct min-time/energy/fuel GTO→tulip solvers

Direct (collocation NLP) solvers for the low-thrust GTO → south-pole tulip
transfer in the Earth–Moon CR3BP (15 kg, 25 mN, Isp 2100 s, ~40-rev spiral).
Min-time reference: t_f = 6.290694 ND = 27.8845 d, ΔV 4.4665 km/s. The
flagship result is the **certified sharp bang-bang min-fuel solution**
(1.15× min-time, 25 switches, defect 2.4e-14, ΔV 3.3696 km/s) and the
ΔV-vs-t_f front built around it.

## Folder map

| where | what |
|---|---|
| `direct/sundman_minfuel/` | **THE canonical library** — Sundman-regularized CasADi+IPOPT solver, energy→fuel homotopy, energy-backbone continuation, PMP certification, front aggregation. Start at `direct/sundman_minfuel/README.md`. |
| `indirect/ms_band/` | Indirect multiple-shooting attack on the hard 1.01–1.11× transition band (own campaign doc + unit tests). |
| `direct/movie/` | Trajectory animations (certified solution with running ΔV meter). |
| `attic/` | Superseded code: fmincon-era min-time/min-fuel NLPs, Sundman prototypes, old continuation experiments. Do not use; see `attic/README.md`. The fmincon min-time formulation notes (density-matched mesh, throttle-on-bound gotcha, mesh-refinement table) are preserved in `attic/README_legacy_fmincon_era.md`. |
| `reviews/` | External code-review records. |

## Key documents

| doc | role |
|---|---|
| `LOW_THRUST_MINFUEL_CAMPAIGN.md` | Full campaign record: every method generation, what failed, why, and the winning recipe. Read first. |
| `HONEST_EVALUATION_DV_TF_FRONT.md` | Candid assessment: what the certification does/doesn't prove, branch structure of the front, open problems. |
| `CODE_CLEANUP_PLAN.md` | This reorganization (phased; Phase 0 done 2026-07-09). |
| `MIN_ENERGY_NOTES.md` | Min-energy (homotopy root) derivation notes. |
| `sundman_minfuel_solution_note.tex/.pdf` | 7-page technical note: OCP, homotopy, Sundman regularization, IPOPT. |

## Entry points

```matlab
cd sundman_minfuel
run_certified_minfuel          % reproduce THE certified 1.15x result (~15 min)
minfuel_at_tf(1.30)            % solve one t_f from the energy backbone
aggregate_front                % collect + PMP-verify + honest 3-class front plot
test_minfuel_lib               % cheap no-solve guardrail checks
```

Batch orchestration (process isolation + watchdog + retry — required because
of sporadic uncatchable CasADi/IPOPT MEX crashes):

```bash
direct/sundman_minfuel/orchestrate/backbone_walk.sh 1.15 1.20 1.25 1.30   # energy chain
direct/sundman_minfuel/orchestrate/sharpen_batch.sh 2 1.30 1.35 1.40      # parallel sharpen
```

Root-level `setup_paths.m` and the two root `.mat`s (`sundman_minfuel_certified`,
`minfuel_from_energy_seed`) are load-bearing for `direct/movie/` and legacy scripts —
dedupe scheduled for cleanup Phase 1.

## Companions

- Indirect (PMP shooting) counterpart + theory note + guided tutorial:
  `indirect/lowThrust_GTO_tulip/`
- Problem source: pumpkynPie `Demos/LunaNet Analysis/lowThrust_GTO_Tulip.m`
