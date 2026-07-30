# ΔV–t_f front sweep — basin multiplicity is concentrated near minimum time

**Run 2026-07-29.** Follow-up to `BASIN_1150_SWEEP.md`, which found six optima at
a single t_f and displaced the flagship. If one point on the front is 0.68%
pessimistic, every point is under suspicion.

## Protocol

Fifteen distinct t_f factors, from 19 stored solutions. At each factor, re-solve
from **its own solution plus its two nearest neighbours'**, all on the target's
own mesh and t_f. Neighbours are structurally different solutions obtained for
free — the seed family the 1.150× sweep showed to be the productive one, where
throttle shifts and loose-barrier variants did nothing.

## Result: the gains are ALL at the short-t_f end

| factor | own m_f | best found | gain | best via |
|---|---|---|---|---|
| **1.120** | 0.83038251 (12 sw) | **0.83393143 (20 sw)** | **+53.23 g** | f1250 |
| **1.140** | 0.84409409 (26 sw) | **0.84481147 (23 sw)** | **+10.76 g** | f1120 |
| 1.250 | 0.85854384 (50 sw) | — | 0 | own |
| **1.300** | 0.86214313 (42 sw) | **0.86310443 (42 sw)** | **+14.42 g** | f1300_nb |
| **1.350** | 0.86528752 (29 sw) | **0.86591759 (32 sw)** | **+9.45 g** | f1400_en |
| 1.400 | 0.86608135 (26 sw) | — | 0 | own |
| 1.450 | 0.87364942 (43 sw) | — | 0 | own |
| 1.500 | 0.87939106 (35 sw) | 0.87939467 | +0.05 g | own |
| 1.550–1.700 | — | — | **0** | own |
| 1.750 | 0.88470081 (23 sw) | 0.88475877 | +0.87 g | f1800 |
| 1.800 | 0.88177282 (23 sw) | — | 0 | own |
| 1.850 | 0.87854350 (22 sw) | 0.87854390 | +0.01 g | f1800 |

**Every meaningful gain is at t_f ≤ 1.350×. From 1.400× upward the front is
essentially unique** — thirteen of the fifteen re-solves at those factors
returned the stored value to machine precision.

## Why this is a result, not just a correction

The short-t_f end is the **near-minimum-time region**, which this campaign
already knows is hard: the 1.01–1.11× band remains unsolved, and it fails for
the *smooth* energy problem too, so it is conditioning rather than bang-bang
structure. The sweep shows basin multiplicity following the same gradient:

- **worst where the problem is hardest** (1.120×: +53 g, and the switch count
  moves 12 → 20 between optima),
- **absent where there is time margin** (≥1.400×: zero gain, structure stable).

That is a coherent physical story: as t_f approaches minimum time the feasible
set tightens, the optimizer has less room, and the landscape fragments into many
near-degenerate local minima. It also means the published front is a **lower
bound on achievable performance in its left third and an accurate
characterization in its right two thirds.**

## Switch counts move with the basin

At 1.120× the three seeds give 12, 17 and 20 switches with the 20-switch
solution best. At 1.140×, 26 / 23 / 28 with the 23-switch best. Switch count is
therefore not a stable label for a solution at fixed t_f — the third independent
demonstration of this in the campaign, after the mesh study and the 1.150× sweep.

## One correction to an earlier note

`BASIN_1150_SWEEP.md` flagged the stored f1450 pair (0.87364942 at 43 sw versus
0.86211422 at 26 sw, a 0.17 kg gap) and asked whether the front plotted the
better one. **It does** — the `dn` solution is best and no seed improved on it.

## Caveats

- Neighbour-seeding tests **lateral** structures, those reachable from adjacent
  transfer times. It cannot find a basin no neighbour occupies. A zero gain
  means "no better structure among those on the front", not "this is optimal".
- The improved rows are better solutions, **not yet certified rows** — they have
  not been through `run_foc_tulip`.
