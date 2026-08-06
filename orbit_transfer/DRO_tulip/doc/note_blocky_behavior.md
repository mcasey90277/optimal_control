# Note for Darin: the "blocky behavior" is (mostly) not an artifact

*Casey, 2026-08-05. Companion figure: `family_walls_15N.png`.*

**Short version: block structure in a (phasing × t_f) map is the expected
signature of a multi-family solution landscape, and we can measure it in our
own verified library. Which block a given run lands in, however, depends on
the solver's seed — so independent solves produce patchwork maps even when
every individual solution is genuinely converged.**

## What we measured

In our DRO→tulip thrust-ladder library (every entry verified: multiple-
shooting residual ~1e-11, flown to the target, accepted unchanged by
`tfMin`), we computed the t_f change between every pair of adjacent phasing
cells at each thrust:

| Thrust | Median neighbor change | Neighbor jumps > 15% |
|---|---|---|
| 15 N | 8.5% | **51** (up to 75%) |
| 5 N | 7.4% | 58 |
| 1 N | 6.7% | 8 |
| 0.5 N | 5.6% | 5 |

Within one solution family, t_f varies smoothly (~5–8% per grid step). The
large jumps are **family walls**: adjacent phase pairs whose minimum-time
solutions are topologically different transfers (different near-side/far-side
routing, different revolution structure). The figure draws these walls in
black at 15 N — the map is visibly carved into blocks, including a 0.66-day
cell walled off beside 0.29-day neighbors.

Note the trend: **many walls at high thrust, few at low.** Short transfers
are acutely phasing-sensitive, so the fast landscape supports many distinct
locally-optimal families; a slow spiral averages the geometry and the
families merge.

## Why runs look blocky

Two mechanisms stack:

1. **The landscape genuinely has multiple families** (physics). We measured
   two feasible 15 N solutions at the *same* phase pair with t_f 0.40 vs
   0.60 days — both converged, both fly to the target.
2. **Which family a solve lands on depends on the seed** (method). Cold or
   independently-seeded solves scatter across families; neighboring grid
   points land on different ones; the map comes out blocky. This is not a
   convergence failure — every cell can be perfectly converged and the map
   still looks like patchwork.

## What to do about it

- **Continuation, not independent solves, for maps.** Walk seeds across the
  grid (each converged cell warm-starts its neighbors) and the map stays on
  one family per region; walls then appear only where the family genuinely
  ends. This is how our library was built.
- **Treat walls as information.** A cell that is much *slower* than its
  neighbors across a wall is usually sitting on a worse local optimum — its
  neighbor's solution, used as a seed, often finds it a faster family. (Our
  own 15 N map has a few such cells; re-seeding them is on our list.)
- **Diagnostic for any suspect run:** compute the neighbor-jump table above.
  Median ~5–10% with isolated >15% walls → multi-family landscape (expected).
  Jumps scattered randomly cell-to-cell with no block structure → an actual
  numerical artifact worth hunting.

If you send one of your blocky runs (any format with phasing + t_f, ideally
+ costates), we can overlay the wall diagnostic on it directly.
