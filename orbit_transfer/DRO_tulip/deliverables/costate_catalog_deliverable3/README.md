# DRO → Tulip Minimum-Time Costate CATALOG (deliverable 3)

The multi-orbit catalog: **3,936 converged minimum-time transfer solutions**
spanning **4 DRO periods × 4 tulip petal counts × a 6×6 phasing torus × 9
thrust levels**. Every entry is a root of the indirect (PMP) boundary-value
problem — hand its `z8` to `pumpkyn.cr3bp.tfMin` and it is accepted
unchanged in about a second. This supersedes nothing: deliverable 2 (the
τ=1/Np=7 pair at 12×12 phasing, 14 rungs) remains the fine-grained sheet for
that pair; this catalog is the coarse map across the orbit families.

## Contents

| File | Role |
|---|---|
| `costate_catalog_dro_tulip.mat` | the catalog — 16 sheets, 3,936 entries, COMPACT format |
| `costate_catalog_pick.m` | five-coordinate lookup (DRO period, petals, dep/arr phase, thrust) with honest warnings |
| `costate_catalog_example.m` | worked example: pick → derive → fly → tfMin → plot → period-axis view |
| `costate_lib_describe.m` | prints the catalog's facts: phase grids (ND + days), thrusts, periods, petals, coverage per sheet, best/worst sheet |
| `costate_catalog_extremes.m` | cheapest/dearest (or fastest/slowest) transfer under any filter (sheet, petals, thrust); side-by-side plot |
| `costate_catalog_extremes_movies.m` | renders the two extremes as side-by-side movies (.mp4 + .gif), same clock: the shorter finishes and freezes while the longer plays out |
| `sample_extremes_movie.mp4` | example output: min vs max ΔV phasing on the τ=2, Np=7 sheet at 5 N |

Requires **pumpkyn / pumpkynPie** on the MATLAB path. The data itself is
dependency-free.

## Quick start

```matlab
L = load('costate_catalog_dro_tulip.mat');
cat = L.costate_catalog_dro_tulip;
[tf_nd, z8, info] = costate_catalog_pick(cat, 2.0, 7, 3.0, 11.6, 5.0);
%                                        tau  Np  dep  arr  thrust(N)
% z8 -> pumpkyn.cr3bp.tfMin as-is;  info.delivered = what you actually got
```

Or run `costate_catalog_example` from this folder. For a fact sheet first:

```matlab
costate_lib_describe('costate_catalog_dro_tulip.mat')
```

## What it covers

- **DRO periods (τ = period in ND; getDRO selects by it):** 0.5, 1.0, 2.0,
  3.0 ND = 2.2, 4.4, 8.9, 13.3 days. The admissible family extends to
  τ ≈ 3.25 (the 100 Mm lunar-vicinity bound).
- **Tulip petal counts:** 5, 7, 9, 12 (periods locked by the
  (Np−2)/(Np−1)-month resonance: 20.9–25.3 days). **pm = −1 branch only**:
  pm = +1 orbits and their costates are exact z-mirrors (CR3BP symmetry) —
  flip all z-components.
- **Phasing:** 6×6 torus per sheet, grid origin (0,0), phases as fractions
  of each orbit's period.
- **Thrust:** 15, 12, 10, 7, 5, 3, 2, 1.5, 1 N at Isp 1710 s, m₀ 150 kg.
- **Coverage: 32–35 of 36 phase pairs on every sheet** (89–97%);
  `sheets(k).has_solution` is the per-(pair, rung) authority.

## The COMPACT format (data minimization)

Only canonical nondimensional quantities are stored — constants once at top
level, per-sheet phase fractions and rung availability/`tf_nd` lookup grids,
and the `z8` vectors themselves (which already contain t_f). Days, ΔV, and
masses are **derived**, and every formula rides along in `cat.derive`:

```
t_days  = t_nd * tStar_s/86400
Tmax_nd = (T_N/m0_kg)*tStar_s^2/(lStar_km*1000)
mf      = 1 - Tmax_nd*tf_nd/c_nd          (all-burn minimum time; exact)
dV_kms  = c_nd*ln(1/mf) * lStar_km/tStar_s
```

`cat.derive.orbit_reconstruction` holds the literal pumpkyn calls that
rebuild both orbits of any sheet.

## The honesty contract (the picker's warnings)

Whenever what is RETURNED differs from what you REQUESTED, the picker prints
a warning stating exactly what you are getting: the nearest **sheet** (if
your period/petal combination isn't stored), the nearest **grid pair** (if
your phases are off-grid), the nearest **solved pair** (if the requested
cell is unsolved), or a seed from a different **rung**. `info.delivered`
carries the same facts programmatically; `info.warned` flags any
substitution. Pass `warnFlag = false` to silence the console for batch use.
The right use of a substituted answer is as a **seed**: hand `z8` to `tfMin`
at your true endpoints/thrust and it converges the difference in seconds.

## How it was made, and how it was checked

Per sheet: a direct Hermite–Simpson collocation solve anchors each phase
pair cold at 15 N (nearly impulsive — converges in seconds), then thrust is
walked down with each rung warm-starting the next, so a pair stays on one
solution family. Every rung passes three gates before entering the catalog:

1. multiple-shooting refinement residual (`ms_tfmin`, analytic STM
   Jacobian): worst case ~1e-10
2. the PMP control law is FLOWN end-to-end and must arrive: worst case
   ~1e-4 km
3. `pumpkyn.cr3bp.tfMin` must return the entry unchanged: worst case
   ‖Δz‖ ~ 1e-9

Total build: ~26 h of unattended compute; ~25 s per verified solution.

## Gotchas, honestly

- **Coverage gaps are real**: 1–4 unsolved pairs per sheet, and lower rungs
  thin out (the 1 N rung holds fewer pairs than the 15 N rung on most
  sheets). `has_solution` is always the authority; the picker warns when it
  substitutes.
- **At a fixed thrust, minimum time = minimum ΔV** (continuous burn makes ΔV
  monotone in t_f). The metrics differ only across thrust levels.
- **Neighboring cells can sit on different solution families** — t_f can
  jump >15% across a "family wall". If an entry looks slow next to its
  neighbor, seeding tfMin from the neighbor may find the faster family.
- The coarse 6×6 grid is for interpolation + seeding, not final answers:
  grid spacing is hours-to-days of phase. tfMin closes the gap from any
  returned seed.

M. Casey / D. Koblick, Coorbital Inc. — August 2026
