# The Costate-Library Pipeline

*Recorded 2026-08-04, after the first library (DRO→tulip, 12×12 phasing torus)
shipped. This is the production process for every future library.*

## The three steps

```
(1) DIRECT SOLVE          (2) REFINE: ms_tfmin        (3) TEST: pumpkyn tfMin
direct collocation   →    multiple shooting      →    feed every entry to the
sweep of the phasing      seeded by the FULL          single-shooting solver;
torus; save states,       state+costate trajectory;   PASS = returned unchanged
controls, costates,       machine-accurate            in seconds
tf per cell               [lambda0; tf] per cell
```

### Step 1 — Direct solve (the map-maker)

Hermite–Simpson collocation (N=800, Sundman mesh, lifted t_f, 500-km lunar
altitude floor) swept over the (departure phase × arrival phase) torus by the
wave engine (`direct/sweep_phasing_direct.m`): wave 0 knocks on cold cells,
wave 1 floods warm starts outward from the openers. A cell is **green only if
the flown control arrives** — the solved control history is integrated
end-to-end at high accuracy and must land within map tolerance (100 km /
10 m/s). Optimizer success certificates alone never color a cell.

**Data contract (non-negotiable, learned the hard way):** every completed cell
immediately saves its full record — states `X`, controls `U`/`Um`, time grid,
`tf`, sign-resolved costate trajectory `lamDef` (all nodes, from the Hager
covector mapping of the NLP defect multipliers), both phases, endpoint states,
and the orbit/thruster definition. Step 2 is only possible because the
*trajectories* of the costates are saved, not just their initial values.

### Step 2 — Refine with `ms_tfmin` (the precision-maker)

**Why single shooting cannot do this:** measured on all 132 green cells, the
PMP flight from the collocation-derived λ₀ misses the target by 12,700–558,000
km (median 47,500 km) — collocation duals are only ~1e-3 accurate and single
shooting amplifies that ~1e3× over a 4-ND spiral. `pumpkyn.cr3bp.tfMin` from
such seeds grinds for hours. Always run this cheap seed-residual probe (one
`tfMinProp` flight per cell) *before* spending solver time.

`indirect/ms_tfmin.m` splits the arc into K segments whose junction states are
unknowns (kills the amplification), seeds every junction from the cell's saved
state+costate trajectory, and solves the full system with trust-region-dogleg
and an **analytic block Jacobian assembled entirely from pumpkyn calls**
(`tfMinProp`'s 210-state mode returns segment STMs; `tfMinEoM` the Hamiltonian
row). Terminal conditions are tfMin's own: r(t_f)=r_f, v(t_f)=v_f,
λ_m(t_f)=0, H(t_f)=0. A **K-ladder (12→24→48)** escalates segmentation for
rough seeds; each rung costs seconds. Per-cell wall budget always; results
saved after every cell.

First-library outcome: 126/132 cells converged to residual ~1e-11, flight miss
~meters, 1–25 s per cell.

### Step 3 — Test against pumpkyn `tfMin` (the acceptance gate)

Every refined entry is fed to the real `pumpkyn.cr3bp.tfMin` as its seed.
**Pass = tfMin returns the seed essentially unchanged (|Δz| < 1e-6) in
seconds** — proving the entry is a root of the single-shooting equations, i.e.
directly usable in Darin's workflow with zero adaptation.

First-library outcome: **126/126 passed** (max |Δz| 2.8e-8, max wall 2.4 s).
Cross-method agreement: median |t_f(direct) − t_f(indirect)| = 2e-8 — an
independent certification of both methods at every phase pair.

## Packaging

`indirect/build_costate_lib.m` → a self-describing `.mat` structure
(`costate_lib_<pair>`): library name, both orbit definitions with
reconstruction recipes (DRO selected **by period** — `getDRO(tau)`'s τ *is*
the period in ND), CR3BP constants, thruster (dimensional + ND), per-entry
phases (fraction / ND / days), `lambda0`, `tf`, tfMin-ready `z8`, verification
metrics, and a `grid` availability view (`has_solution`, `tf_days`,
`entry_index`). Ships with `costate_lib_nearest.m` (torus-metric picker,
pumpkyn house style) and `costate_lib_example.m` (worked example).

## Axes of growth (the production roadmap)

- **Thrust axis** (in progress): continuation ladder per cell from the
  converged low-thrust entry toward Darin's 1–15 N / Isp 1710 s regime.
  Adaptive step ratio + ms_tfmin fallback seeded from the previous rung's
  propagated trajectory; basin jumps (revolution-count changes) are the
  expected failure mode of plain tfMin continuation.
- **Period axis**: repeat the whole pipeline per DRO period τ — each τ is one
  library sheet; together they form the ΔV–period–phasing map.
- **Other orbit pairs**: the pipeline is pair-agnostic; endpoints come from
  pumpkyn's nine catalogued families.

## Standing lessons baked into the pipeline

1. **Save every data product, per cell, at completion** — with orbit metadata.
2. **Green = continuously verified**, never an optimizer certificate.
3. **Prior runs are data**: rerun only what a reference run proved worthwhile.
4. **Measure seed quality before spending solver time** (the PMP-flight probe).
5. **Guard every solver against junk iterates** (sanity bounds + try/catch —
   an unguarded iterate once killed the MATLAB process via integrator
   collapse).
6. **Chain the stages automatically**; humans should read results, not
   babysit handoffs.


## Multi-family catalogs (added 2026-08-06)

The pipeline is now family-agnostic. Shared code lives in
`orbit_transfer/costate_common/` — the seed of an optimal-control library
alongside pumpkyn's astrodynamics:

- `get_family_orbit.m` — THE one place a family name + parameters becomes a
  propagated periodic orbit (dro / tulip / halo / dpo / lyapunov, all via
  pumpkyn getters + cont_np + prop). Engines build endpoints only through it.
- `survey_family_bounds.m` — the admissibility survey (>= 500 km periselene,
  <= 100 Mm from the Moon) for ANY family, with a periodicity guard: members
  whose seed did not truly converge are labeled NON-PERIODIC and excluded,
  never trusted (an interpolated L1 halo "member" showed 1e10-km excursions;
  and note the guard is written `~(err < tol)` so NaN closures fail it too).

**Measured Halo admissible box (2026-08-06,
`HALO_tulip/direct/results/halo_bounds.mat`):**

- **L2 southern is the workhorse**: continuously admissible for periods
  7.1–15.1 days (tau 1.61–3.42), periselene 2,800→49,900 km. Northern is the
  exact z-mirror.
- L1: mostly inadmissible (the family plunges below the lunar surface over
  most of its range); two clean members at tau = 1.80 and 2.74.
- L3: excluded wholesale — 0.39–1.05 Gm from the Moon, far outside the box.
- Short-period L2 (NRHO territory, < 7 d) falls below the 500 km floor at
  this grid; the true boundary deserves a finer survey before writing NRHOs
  off.

The Halo campaign (HALO_tulip/) proceeds on L2 southern, coarse periods
{1.75, 2.2, 2.8, 3.4} ND, reusing the ladder pipeline through
get_family_orbit.

## Fixed-time costs: the min-energy pilot (added 2026-08-14)

The three steps carry over unchanged in *shape* to a fixed-t_f cost; only
the bindings change. Proven on flagship 12×12 DRO→tulip cells with
MINIMUM ENERGY, J = ∫s² dt at t_f = γ·t_f^min (`DRO_tulip/run_minenergy_pilot`):

1. **Direct** — the same transcription (`casadi_mintime_dro`) with
   `objective='energy'` + `tfFix`, warm-started from the min-time cell
   stretched to γ t_f. Same Sundman HS mesh, floor, covector extraction.
2. **Refine** — `indirect/ms_minenergy`, a fixed-t_f binding of `ms_bvp`
   (`opts.fixedTf`): unknowns λ₀ (7), terminal r/v matched + λ_m(t_f) = 0,
   no H(t_f) = 0 (H is a first integral, reported and checked). Seeds from
   the harvest converge in 2–3 iterations at K = 12 — no K escalation.
3. **Accept** — pumpkyn has no min-energy solver, so the gate is the
   generic `costate_common/ss_bvp_accept`: single shooting (`ms_bvp` at
   K = 1) on the same closures, PASS = |Δz| < 1e-6 at the single-shooting
   residual floor (1e-6; with/without-STM full-arc flights differ ~6e-7).

Extra gates for a fixed-time cost: |ΔH| along the indirect flight
(absolute, < 1e-6), direct-vs-indirect J agreement (~1e-5, collocation
order), throttle interior somewhere (else it is min-time in disguise).
Result: 5/5 records pass every gate; full table in `FINDINGS.md`. Lesson to
carry: fixed-time problems at different γ are NOT nested (fixed rotating-
frame endpoints), so J and m_f can be non-monotone in γ — grid γ per cell
before declaring a catalog axis. Min-fuel entries are the energy→fuel
homotopy on these seeds.

## The conjugate-point sweep (added 2026-08-23)

The pipeline's first **second-order** pass over shipped product: every entry
of all five catalogs now carries a Jacobi (conjugate-point) verdict. Engine:
`costate_common/conj_catalog_pass`; instrument:
`costate_common/ms_conjugate_test` (free-time quotiented form, migration #3).

**Per entry:** rebuild both endpoint orbits from the sheet recipes (v1 sheets
fall back to dep = dro(tauDRO), arr = tulip(Np, pm)), fly the stored z8 with
`tfMinProp`, cut into K = 24 junctions, re-solve `ms_tfmin(conjTest)` seeded
AT the converged solution — 1 Newton iteration — then sweep the quotiented
determinant along the chained segment STMs.

**Result (33 min wall, 3 batched passes, ~0.14 s/entry):**

| catalog | entries | pass | fail | unverified |
|---|---|---|---|---|
| DRO → tulip | 3,936 | 3,936 | 0 | 0 |
| HALO → tulip | 3,980 | 3,980 | 0 | 0 |
| DPO → tulip | 3,932 | 3,931 | **1** | 0 |
| L1 → L2 | 1,952 | 1,952 | 0 | 0 |
| L2 → L1 | 2,096 | 2,096 | 0 | 0 |

**The one refuted entry:** DPO τ=2 → tulip Np=7, phases (2/3, 2/3), 15 N
(t_f = 0.115 ND, ΔV 5.12 km/s): one sign change strictly interior
(atFinal = 0) at re-solve fidelity 1.6e-12. A shipped, tfMin-accepted
extremal that is **not a local minimum** — its cell has a cheaper
neighboring transfer. Its lower-thrust rungs all pass, so it is an expensive
branch, not a corrupt solve. The single hit is also the sweep's best
validation: the test discriminates rather than rubber-stamps.

### Findings worth keeping

1. **The re-solve honesty gate is what makes 15,896 verdicts trustworthy.**
   A verdict is recorded only when the ms re-solve reproduces the stored
   entry (|z − z8| < 1e-6); measured 1e-12..1e-10 on every entry, zero
   exceptions. Seeded at its own solution, ms converges in 1 iteration —
   so the sweep certifies exactly the shipped solutions.
2. **Second order at catalog scale costs minutes, not campaigns.** 0.13–0.18 s
   per entry warm; the STMs ride free on the Newton Jacobian. There is no
   longer a cost argument for shipping entries without a Jacobi verdict.
3. **99.994% of accepted min-time entries are conjugate-point-free.** The
   three-gate acceptance pipeline almost always lands on locally minimizing
   branches — but not always, which is precisely why the verdict must be
   stored per entry rather than assumed.
4. **Storage:** per-sheet `conj_pass` / `conj_ncross` / `conj_atfinal` int8
   grids + a top-level `conj_test` provenance block, validated by
   `catalog_schema` (optional fields, v1/v2 compatible). Canonical .mats
   updated with `.bak_conj` backups; `*_conjprog.mat` sidecars kept.

### Standing rule going forward

**Run `conj_catalog_pass` before packaging any new catalog or deliverable.**
The sweep is cheap, resumable, and upgrades every entry from "certified
extremal" to "extremal with no conjugate point in (0, t_f)".

### Residue

- The deliverable zips (3–6) predate the verdicts; re-ship when Darin wants
  them (the DPO refuted entry is documented in `DPO_tulip/README.md`).
- The deliverable-2 fine-sheet library (1,105 entries, v2-library format)
  is not swept — needs a small format adapter.
- Free-time form only: min-energy (fixed-t_f) entries still have no
  applicable test (`costate_common/TODO.md`).
- Junction resolution: K = 24 sampling means a conjugate pair inside one
  segment can hide; "no sign change" is strong, not airtight.

Presentation deck (rigorous + intuitive walkthrough of the test):
`~/Documents/myLatex/conjugate_test_slides.html`.
