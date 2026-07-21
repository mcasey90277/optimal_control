# GTO_ELFO — low-thrust GTO → ELFO transfers (Earth–Moon CR3BP)

Minimum-fuel / minimum-time low-thrust transfers from GTO to the ELFO
(elliptical lunar frozen orbit) benchmark target `im_elfo_optimum` (the proj7
lunar-nav constellation's IM-style orbit), in the Earth–Moon CR3BP. Same
spacecraft as the tulip campaign (15 kg, 25 mN, Isp 2100 s); same solver
architecture, retargeted.

## Status

- **Direct: working.** Full pipeline certified end-to-end (see
  `direct/elfo/README.md`):
  - min-ENERGY seeds via the two-primary gravity-homotopy ladder
    (`gen_elfo_energy_gravhom`, 2026-07-13), t_f-grid of seeds available;
  - energy → fuel sharpening to bang-bang (`gen_elfo_minfuel`);
  - **ΔV–t_f min-fuel front mapped** (2026-07-15);
  - **min-time anchor (Route B, all-burn direct): t_f,min = 6.0962 ND =
    27.02 d** (`casadi_mintime_freetf` + `gen_elfo_mintime`, machine-tight,
    verified) — the front is labeled against t_f/t_f,min.
- **Indirect: not started.** `indirect/` is a placeholder; the ELFO min-time
  indirect solve ("Route C") and an indirect min-fuel counterpart are open —
  see `TODO.md`.

## Three objectives, one pipeline

Min-time / min-energy / min-fuel share one solver core and form one homotopy
chain (so `direct/` is not split by objective):

| objective | role in the chain | entry point |
|---|---|---|
| min-time | anchor: `t_f,min` = 6.0962 ND (all-burn mode) | `direct/elfo/gen_elfo_mintime` |
| min-energy | homotopy root (same fuel solver at ε=1) + gravity-homotopy seeds | `direct/elfo/gen_elfo_energy_gravhom` |
| min-fuel | target (ε=0 bang-bang) via the ε:1→0 sweep | `direct/elfo/gen_elfo_minfuel`, `run_elfo_minfuel` |

## Folder map

| where | what |
|---|---|
| `direct/elfo/` | The working direct campaign: seed generation, homotopy sharpening, front, movies, smoke tests. Entry point: `run_elfo_minfuel`. |
| `indirect/` | Placeholder (README stub) for the indirect (PMP shooting) work. |

## Shared machinery (deliberate cross-references)

- `../cr3bp_common/` — shared CR3BP problem definition (`cr3bp_lt_params`,
  `minfuel_config`, `gto_elfo_endpoints`) + pumpkyn path.
- `../GTO_tulip/direct/sundman_minfuel/` — the Sundman min-fuel engine this
  campaign reuses (`casadi_minfuel_sundman`, `insertion_states`,
  `minfuel_at_tf`), retargeted to ELFO.
- `../GTO_tulip/indirect/min_time/` — the PMP min-time root used for
  retargeting experiments.

## Run

```matlab
cd direct/elfo
setup_paths
run_elfo_minfuel               % end-to-end: solve -> export -> verify -> movie
```
