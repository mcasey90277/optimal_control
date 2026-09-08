# DRO -> Tulip Minimum-Time Costate CATALOG (deliverable 3)

**4439 converged minimum-time transfer solutions** spanning **4 DRO periods x 4 tulip petal counts x a 6x6 phasing torus x 11 thrust levels**. Every entry is a root of the indirect (PMP) boundary-value problem: hand its `z8` to `pumpkyn.cr3bp.tfMin` and it is accepted unchanged in about a second.

Catalog built 2026-09-01; this package generated 2026-09-08 by `build_dro_deliverable.m`.
Every number below is computed from the catalog, not transcribed.

## What it covers

| axis | values |
|---|---|
| DRO period tau (ND) | 0.5, 1, 2, 3 |
| tulip petal count Np | 5, 7, 9, 12 |
| family branch pm | -1 (pm = +1 orbits and costates are exact z-mirrors) |
| phasing | 6x6 torus per sheet, phases as fractions of each orbit period |
| thrust (N) | 15, 12, 10, 7, 5, 3, 2, 1.5, 1, 0.75, 0.5 |
| propulsion | Isp 1710 s, m0 150 kg |
| time of flight | 0.049 .. 1.445 ND = 0.22 .. 6.41 days |

**Coverage: 4439 of 6336 grid slots (70.1%)**, 57%-80% per sheet.
`sheets(k).has_solution` is the per-(pair, rung) authority. Coverage
thins at low thrust, which is real and worth knowing before you plan
around it:

| thrust (N) | 15 | 12 | 10 | 7 | 5 | 3 | 2 | 1.5 | 1 | 0.75 | 0.5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| solved | 80% | 86% | 86% | 85% | 82% | 78% | 73% | 65% | 49% | 47% | 40% |

## How it was made, and how it was checked

Per sheet: a direct Hermite-Simpson collocation solve anchors each phase
pair cold at the highest thrust, then thrust is walked down with each rung
warm-starting the next, so a pair stays on one solution family. Every
entry passes three FIRST-ORDER gates:

1. multiple-shooting refinement residual (`ms_tfmin`, analytic STM Jacobian): ~1e-10
2. the PMP control law is FLOWN end to end and must arrive: ~1e-4 km
3. **`pumpkyn.cr3bp.tfMin` returns the entry unchanged** (|dz| ~ 1e-9) -- a
   foreign witness: different code, integrator and root finder

### Second-order verification (new since the previous deliverable)

**Conjugate-point test** (`conj_test`, 2026-09-08, K = 24): the free-time
quotiented Jacobi determinant, sampled at the multiple-shooting junctions.
`sheets(k).conj_pass`: 1 = no conjugate point found in (0, tf); 0 = a sign
change strictly inside, i.e. the entry is NOT a local minimum; -1 = not
decided. Result: **4427 pass, 11 fail, 1 undecided**.

**Sufficiency-hypothesis gates** (`hyp_gates`, 2026-09-07): the three hypotheses
the conjugate test assumes, checked per entry and stored as
`gate_min_lamv` (strong Legendre, min|lam_v| > 0), `gate_min_qmt`
(all-burn is the PMP control, min Q_mt > 0) and `gate_dimS` (no abnormal
lift, dim S = 1). **All entries pass all three.**

With `conj_pass = 1` these are the Bonnard-Caillau-Trelat sufficiency
hypotheses, verified numerically at the sampled times. What is NOT
claimed: the determinant is sampled at junctions, so a conjugate pair
inside one segment is invisible, and the second-order layer -- unlike
the first -- has no independent second implementation behind it.

## The COMPACT format

Only canonical nondimensional quantities are stored: constants once at top
level, per-sheet phase fractions and rung availability, `tf_nd` lookup
grids, and the `z8` vectors (which already contain t_f). Days, dV and
masses are DERIVED, and every formula rides along in `cat.derive`:

```
days_from_nd           t_days = t_nd * constants.tStar_s / 86400
dro_period_nd          sheets(k).tauDRO  (tau IS the period; getDRO selects by it)
phase_days             frac * period_nd * tStar_s/86400 (period_nd: tauDRO or period_tulip_nd)
thrust_nd              see thruster.thrust_nd_formula
m_final                mf = 1 - Tmax_nd*tf_nd/c_nd   (all-burn minimum time)
deltaV_kms             dV = c_nd*log(1/mf) * lStar_km/tStar_s
orbit_reconstruction   DRO: getDRO(tauDRO)->cont_np->prop;  tulip: getTulip(2*pi*(Np-2)/(Np-1), Np, pm)->cont_np->prop;  state at phase f: interp1(t, rv, mod(f,1)*t(end), 'spline')
pm_plus_one            z-mirror the pm=-1 orbit and costates (flip all z components)
```

## Quick start

```matlab
L = load('costate_catalog_dro_tulip.mat');
cat = L.costate_catalog_dro_tulip;
[tf_nd, z8, info] = costate_catalog_pick(cat, 2, 7, 3.0, 11.6, 5);
%                                        tau  Np  dep  arr  thrust(N)
% z8 -> pumpkyn.cr3bp.tfMin as-is;  info.delivered = what you actually got
```

Or run `costate_catalog_example`. For a fact sheet: `costate_lib_describe`.
Requires **pumpkyn / pumpkynPie** on the MATLAB path; the data itself is
dependency-free.

## The honesty contract (the picker's warnings)

Whenever what is RETURNED differs from what you REQUESTED, the picker
prints exactly what you are getting: the nearest sheet, the nearest grid
pair, the nearest SOLVED pair, or a seed from a different rung.
`info.delivered` carries the same facts programmatically and `info.warned`
flags any substitution. The right use of a substituted answer is as a
**seed**: hand `z8` to `tfMin` at your true endpoints and thrust.

## Gotchas, honestly

- **Coverage gaps are real** and concentrated at low thrust (see the table
  above). `has_solution` is always the authority.
- **Only the pm = -1 branch is stored.** The mirror branch follows from the
  CR3BP symmetry (flip all z components) but is not catalogued.
- **At fixed thrust, minimum time = minimum dV** (continuous burn makes dV
  monotone in t_f). The metrics differ only ACROSS thrust levels.
- **Neighbouring cells can sit on different solution families** -- t_f can
  jump >15% across a family wall. Seeding tfMin from a neighbour may find
  the faster family.
- The 6x6 grid is for interpolation and seeding, not final answers.

M. Casey / D. Koblick, Coorbital Inc.
