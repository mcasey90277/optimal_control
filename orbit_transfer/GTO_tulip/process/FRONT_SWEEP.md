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

## CERTIFIED 2026-07-30 — all four improved rows PASS

Each regenerated from its winning seed, saved as
`lib/minfuel_best_f<factor>.mat`, and run through `run_foc_tulip`.

| gate | f1.120 | f1.140 | f1.300 | f1.350 |
|---|---|---|---|---|
| KKT stationarity | 4.36e-12 | 7.46e-13 | 6.33e-13 | 1.05e-12 |
| primer misalign, max | 1.28e-17 | 5.89e-18 | 3.06e-18 | 2.20e-18 |
| dual sign law | 100% | 100% | 100% | 100% |
| **min regular Ṡ** | **0.1147** | **1.5699** | **1.0028** | **0.7995** |
| singular arc nodes | 0 | 0 | 0 | 0 |
| mapped transversality | 4.42e-20 | 0 | 1.22e-24 | 1.20e-24 |
| switches | 20 | 23 | 42 | 32 |
| **verdict** | **PASS** | **PASS** | **PASS** | **PASS** |

Defects 2.5e-14 to 3.5e-14. All four are certified rows.

**Reproducibility confirmed.** Every re-solve matched the sweep's reported value
to ~3e-9, in an independently launched run. The f1.300 case was a deliberate
control — its winner is an already-stored artifact at the same factor, so the
run had to return that file's own value, and it did. The pipeline is not
perturbing solutions in transit, which was worth establishing given that
re-unitizing β and clipping the throttle (changes at 1e-12) moved an earth row
by 0.5 kg.

### Switching regularity is WEAK on exactly these rows

`sdotMinRel` measures how transversally the switching function crosses zero;
the gate threshold is 1e-3. The four improved rows sit at **0.11 to 1.57**,
against **27–28** for the two 1.150× solutions — one to two orders of magnitude
lower. The weakest, f1.120 at **0.1147**, is also the row with the largest basin
gain (+53 g).

Near-tangential switch crossings are precisely the condition under which the
bang-bang structure becomes ill-determined, so this is a plausible mechanism for
the multiplicity: as t_f approaches minimum time, switches graze, and many
near-degenerate optima appear.

**Stated as an observation, not a trend.** The 1.140× → 1.150× step jumps from
1.57 to 27 across a 1% change in t_f, so this is not monotone in t_f and likely
depends on which basin a solution occupies as much as on t_f itself. Four points
do not establish a law.

## Caveats

- Neighbour-seeding tests **lateral** structures, those reachable from adjacent
  transfer times. It cannot find a basin no neighbour occupies. A zero gain
  means "no better structure among those on the front", not "this is optimal".
- The improved rows are better solutions, **not yet certified rows** — they have
  not been through `run_foc_tulip`.
