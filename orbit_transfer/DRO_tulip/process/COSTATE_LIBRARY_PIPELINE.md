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
