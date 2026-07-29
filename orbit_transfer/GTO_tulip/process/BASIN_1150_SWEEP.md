# Basin sweep at t_f = 1.150x — the flagship is not the best optimum

**Run 2026-07-29.** Part A of `docs/superpowers/plans/2026-07-29-basin-and-multiphase.md`,
whose success criteria were fixed before the run.

## Why this was run

The GTO->tulip paper cannot report a final mass at 1.150x while this campaign's
own TODO records **two certified optima 1.43% apart in dV**. If a referee
reproduces and lands on the other, the result is worse than not publishing.

## Protocol

Thirteen seeds, every one solved on the **flagship's own mesh and t_f**
(N = 4001, t_f = 7.234298, tauf0 = 151.683747) so the comparison is
like-for-like. Seed families: four source solutions, each plus +/-4-node
throttle shifts, plus a loose-barrier variant of the flagship. Acceptance:
`Solve_Succeeded` and defect <= 1e-8. Comparison: higher final mass wins
(minimum fuel maximizes m_f), one-sided, per campaign convention.

## Result: there are at least SIX optima, not two

| seed | m_f | switches | defect |
|---|---|---|---|
| flagship (certified) | 0.84906583 | 25 | 1.9e-14 |
| flagship, shift -4 / +4 / loose | 0.84906583 | 25 | -- |
| psr_direct + both shifts | 0.84708554 | 25 | 2.9e-14 |
| psr_refined (projected) | 0.84805124 | 25 | 5.1e-14 |
| psr_refined, both shifts | 0.84797502 | 25 | 2.5e-14 |
| solve_minfuel, shift +4 | 0.84988608 | **26** | 4.1e-12 |
| **solve_minfuel (from the 14-sw artifact)** | **0.85008728** | **24** | 1.3e-13 |

The two-optima record was an undercount by a factor of three, and the switch
count is not even fixed across them: 24, 25 and 26 all appear at the same t_f.

## The certified flagship is beaten

|  | m_f (kg) | propellant (kg) | switches | defect |
|---|---|---|---|---|
| flagship | 12.735987 | 2.264013 | 25 | 2.04e-14 |
| **best found** | **12.751309** | **2.248691** | **24** | 1.30e-13 |

**15.322 g less propellant — 0.677% of the flagship's 2.264 kg.**

## Against the pre-registered expectations

1. *"The pair resolves in favour of 0.849066."* **Wrong in framing.** It does
   resolve that way *within the pair*, but the pair was the wrong question: a
   third structure beats both.
2. *"A seed from one optimum will not reach the other."* **Confirmed.** Flagship
   seeds all return 0.84906583 exactly; psr_direct seeds all return
   0.84708554. These are genuinely distinct basins.

## Caveats

- The `psr_refined` rows mix a mesh projection (4017 -> 4001) with a basin
  move and are **not clean data points**. The 24- and 26-switch results are
  clean, both native 4001.
- This is **best found under the seed families listed**, never "the optimum".
  A sweep raises the floor; it never proves a ceiling.
- The winner has **not** yet been through the campaign's certification gates
  (`run_foc_tulip`, `verify_pmp`). It is a better solution, not yet a certified
  row.

## Consequences

- **The flagship should not be published as-is.** The paper's honest position:
  at 1.150x the problem admits at least six local optima spanning 24-26
  switches; the best found under the documented search is reported; the search
  is reproducible. The multiplicity is itself a result worth reporting for a
  bang-bang problem of this size.
- **The dV-t_f front is likely an envelope over a much richer structure.** If
  one t_f carries six optima, the front as mapped is a lower bound on
  achievable performance, not a characterization of it.
