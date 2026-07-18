# earth_elliptic_to_geo — min-fuel low-thrust orbit-transfer pipeline

Direct-method (collocation NLP) reproduction of the **minimum-fuel** low-thrust
Earth orbit transfer of Haberkorn, Martinon & Gergaud, *"Low thrust minimum-fuel
orbital transfer: a homotopic approach,"* JGCD 27(6), 2004 — a 1500 kg satellite
from a low, elliptic, 7°-inclined orbit (P = 11625 km, e = 0.75) to equatorial
GEO. The pipeline produces one **row of the paper's Table 3** per thrust level
(m_f, propellant, ΔV, switch count, revolution count), plus a trajectory plot and
movie.

> **Two-body, Earth-centered — NOT CR3BP.** The only gravity is Earth's central
> field; there is no third body and **no Moon**. (The CR3BP + lunar-gravity
> low-thrust work is a *different* campaign, `../NLP_lowThrust_GTO_tulip/`.) Full
> methodology — coordinate systems, the optimal-control problem, the MEE
> formulation, discretization, and numerical lessons — is in
> **`doc/table3_method_note.tex`**.

## What the pipeline computes

The continuous optimal-control problem, per thrust level `T`:

```
minimize    fuel = m0 − m(t_f)              (equivalently  max m(t_f))
over        throttle δ(t) ∈ [0,1],  thrust direction β(t) ∈ S²  (unit vector)
subject to  MEE Gauss dynamics  dx/dt = a(x) + (T/m)·B(x)·(δβ)   [Earth 2-body + thrust]
            m(0)=m0,  x(0)=initial elliptic orbit,  x(t_f)=GEO (in elements)
            t_f = c_tf · t_{f,min}(T)        (fixed final time)
```

solved **directly** (transcribe → NLP → IPOPT) in Modified Equinoctial Elements
with true longitude `L` the independent variable and the total longitude span
`ΔL` a scalar decision variable. Bang-bang min-fuel is reached from a smooth
root by the Bertrand–Épénoy energy→fuel homotopy `J_ε = ∫δ − ε∫δ(1−δ)`, `ε:1→0`.

## How the pipeline fits together

```
run_gergaud.m                         ← FRONT DOOR (set thrust + endpoints, run)
     │
     ├── auto mode ─────────────────▶ load certified cache ─┐
     │                                                      │
     └── solve / probe mode:                                │
             run_mintime_mee.m   (free-L min-time anchor t_{f,min})
                    │  tfMinAnchor                           │
                    ▼                                        │
             run_transfer_mee.m  (fixed-t_f fuel solve)      │
                    │   seed → homotopy → report             │
                    │   ├── mee_seed.m        (warm start)   │
                    │   ├── homotopy_mee.m    (ε:1→0 sweep)  │
                    │   └── casadi_lt_mee.m   (NLP core) ◀── lt_mee_rhs.m (MEE dynamics)
                    ▼                                        │
             psr_mee_refine.m  (1 N / 0.5 N: switch-time mesh refinement)
                    │                                        │
                    ▼                                        ▼
             gergaud_row.m / gergaud_row_str.m  ◀───────────┘   (assemble + print Table-3 row)
                    │
                    ├── gergaud_plot.m            (static trajectory PNG)
                    └── mee_res_to_cart_res.m ─▶ transfer_movie.m   (mp4 + gif)
```

`run_ladder.m` is the batch analog of the solve path: it descends a whole
thrust list (10 → 5 → 2.5 → 1 → …) with warm-start continuation and prints the
full Table-3-style summary; `run_gergaud` is the single-row front door onto the
same machinery.

---

## Files

### `run_gergaud.m` — front door (start here)
PARAMETERS-block entry (edit the block and run, or call `row = run_gergaud(opts)`).
Set `thrustN ∈ {10,5,2.5,1,0.5,0.2,0.1}` N, the initial orbit (`P0_km,e0,i0_deg`,
default paper GTO) and final orbit (`Pf_km,ef,if_deg`, default GEO), `ctf`
(default 1.5). Three run modes:
- `auto` — reuse the certified cached row for default endpoints (instant for
  10/5/2.5/1/0.5 N), else solve.
- `solve` — run the live pipeline (required for custom endpoints).
- `probe` — force a live solve including the never-certified 0.2/0.1 N rungs,
  with an honest up-front conditioning-wall warning; reports the true
  `certified` flag rather than fabricating a row.

Emits the Table-3 row via `gergaud_row_str`, and (unless `returnOnly`) a plot
(`makePlot`) and movie (`makeMovie`) to `results/gergaud_<tag>.{png,mp4,gif}`.

### `run_mintime_mee.m` — free-longitude min-time anchor
`out = run_mintime_mee(thrustN, nodesPerRev, cfg)`. Two-basin (cold seed +
fuel-warm) keep-best min-time solve giving `t_{f,min}` for the fixed-time fuel
problem. `cfg` accepts `.xf` / `.initElems` (endpoint overrides), Stage A is
skipped for custom endpoints. Returns `.tfmin`, `.dL_mt`, `.revs`, `.solverOut`.

### `run_transfer_mee.m` — fixed-t_f fuel solve
`res = run_transfer_mee(cfg)`. Seed → guarded ε:1→0 homotopy → structure report
→ certified-only cache. `cfg`: `.thrustN .ctf .tfMinAnchor .xf .initElems .tag
.nodesPerRev .maxIter .warmStart`. Returns `res.report.{m_f_kg,switches,revs,
edge,incDeg,defect,certified}`, `res.fuel.{X,U,dL}`, `res.sigma`.

### `casadi_lt_mee.m` — solver core (NLP)
L-domain trapezoidal collocation with `ΔL` a scalar decision variable;
cone-eliminated control `[β(3);δ]`; modes `'mintime'` (δ≡1) and `'fixedtf'`
(ε energy→fuel objective). `opts.xf` [5×1] is the terminal target (default GEO
`[1;0;0;0;0]`); `opts.x0` [7×1] the initial MEE state. CasADi Opti + IPOPT
(`mumps_pivot_order=0`, tfTarget-relative `t` bound).

### `lt_mee_rhs.m` — MEE Gauss dynamics
`[dXdL, Ldot] = lt_mee_rhs(X, U, par)`. State `[P;ex;ey;hx;hy;m;t]`, control
`[β_RTN(3);δ]`; the paper's Gauss equations in RTN, converted `d/dt → d/dL` via
`Ldot`. Written MX-safe (inline floor-wrap of `L`, no `norm/abs/if` on state).
Carries the corrected `L̇` thrust term (paper's printed `1/m` is a typo → `Tmax/m`).

### `mee_seed.m` — warm-start seed
`[sigma,X0,U0,dL0,info] = mee_seed(par, opts)`. Constant-throttle ode113-in-`L`
propagation from the initial orbit at `L0=π`, sampled at uniform-σ nodes
(defect-free by construction). `opts.initElems` overrides the initial state
(default = paper literal); `opts.stopP` / `opts.nRev` set the span.

### `homotopy_mee.m` — energy→fuel sweep
`[best,tbl] = homotopy_mee(sigma,X0,U0,dL0,opts)`. Guarded ε:1→0 continuation
(loose first step, tight thereafter; never advance/cache on a failed step);
per-ε-step resume cache. Forwards `opts.xf` to every solve.

### `psr_mee_refine.m` (+ `psr_switch_score_mee.m`, `psr_refine_sigma_mee.m`, …)
PMP-steered switch-aware mesh refinement, ported from `../NLP_lowThrust_GTO_tulip/PSR/`.
Used at 1 N and 0.5 N to sharpen switch times below the base mesh width. **Does
not yet thread a custom terminal target** — `run_gergaud` skips PSR for custom
endpoints (see `TODO.md`).

### `run_ladder.m` — thrust-continuation batch orchestrator
`run_ladder(thrustList, cfg)`. Descends a strictly-decreasing thrust list, warm-
starting each rung from the one above, resume-safe per-rung cache, prints the
Table-3 summary + R0-law spread. The multi-row analog of `run_gergaud solve`.

### Row formatter — `gergaud_row.m` / `gergaud_row_str.m`
Pure functions: assemble the Table-3 row (prop = m0 − m_f, ΔV via the rocket
equation, `revs_paper` lookup) and format it as a fixed-width block. An
uncertified row gets an `UNCERTIFIED` banner; a certified-but-footnoted row
(0.5 N anchor-free, custom-endpoint) prints a `NOTE:` line.

### Visualization — `mee_res_to_cart_res.m`, `gergaud_plot.m`, `transfer_movie.m`
`mee_res_to_cart_res` reconstructs inertial `(r,v)` and rotates RTN `β`→inertial
from an MEE solution (the adapter that lets the Cartesian-format `transfer_movie`
render an MEE result). `gergaud_plot` writes a static trajectory PNG;
`transfer_movie` writes the throttle-colored mp4 + gif (accepts a struct or a
`.mat` path).

### `kepler_lt_params.m` — constants + canonical units
LU = GEO radius 42165 km, TU = √(LU³/μ) so μ=1, mass unit = m0. Isp default
2000 s. Single source of physical constants.

### `elements_to_cart.m` / `cart_to_elements.m`
Algebraic MEE ↔ inertial `(r,v)` maps (roundtrip-tested), used by the seed,
the reconstruction check, and the viz adapter.

### PMP verifier — `verify_pmp_mee.m` (+ `mee_dual_to_costate.m`, `mee_primer_switch.m`)
First-order Pontryagin certificate from the NLP's KKT duals (primer alignment,
switching-sign law). Diagnostic only; **primal certifications (defect/terminal)
never depend on it.** See `CAMPAIGN.md` footnote 5 for the open dual anomaly.

### Cartesian counterparts (original reproduction)
`casadi_lt_2body.m`, `run_transfer.m`, `run_mintime.m`, `homotopy_2body.m`,
`seed_2body.m`, `verify_pmp_2body.m` — the Sundman-regularized Cartesian solver
that reproduced the paper's 10 N case but stalled at 5 N; the MEE stack above
superseded it. Kept for the cross-formulation gate (see `CAMPAIGN.md`).

---

## Usage

```matlab
cd earth_elliptic_to_geo
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'))      % CasADi on path

% one Table-3 row, default paper endpoints, reuse the certified cache:
row = run_gergaud(struct('thrustN',10,'runMode','auto'));   % prints the row; row struct returned

% custom initial orbit, live solve + plot + movie:
run_gergaud(struct('thrustN',2.5,'runMode','solve','e0',0.6,'i0_deg',10, ...
                   'makePlot',true,'makeMovie',true));

% the whole ladder as a batch:
run_ladder([10 5 2.5 1])
```

Tests (fast, `matlab -batch`): `test_mee_xf`, `test_mee_seed_initelems`,
`test_mee_threading`, `test_mee_res_to_cart`, `test_gergaud_row`,
`test_run_gergaud_auto` (front-door suite); plus the campaign's no-solve guards
`test_params/elements/dynamics/terminal/seed/mee_rhs/mee_seed/...`.

## Related files (outside this directory)

- **Method note (full methodology + flow diagram):** `doc/table3_method_note.tex`
- **Campaign narrative + results tables + honesty footnotes:** `CAMPAIGN.md`
- **Outstanding work:** `TODO.md`
- **CasADi/IPOPT:** `~/casadi-3.7.0`
- **Parent campaign (CR3BP, shared solver architecture):** `../NLP_lowThrust_GTO_tulip/`
- **Problem source:** Haberkorn, Martinon & Gergaud, JGCD 27(6), 2004
  (`min_fuel_papers/Gergaud-Haberkorn-Martinon-JournalGuidance2004-preprint.pdf`)

## Results (headline)

Certified min-fuel ladder at c_tf = 1.5 (full table + honesty footnotes in
`CAMPAIGN.md`):

| T [N] | m_f [kg] | switches | revs (ours / paper) | movie |
|---|---|---|---|---|
| 10  | 1377.10 | 19  | 7.33 / 7.5  | `results/movie_MEE_10N.{mp4,gif}` |
| 5   | 1364.54 | 32  | 14.16 / 15  | `results/movie_MEE_5N.{mp4,gif}` |
| 2.5 | 1369.79 | 76  | 27.84 / 30  | `results/movie_MEE_2p5N.{mp4,gif}` |
| 1   | 1371.44 | 171 | 69.15 / 74.5 | `results/movie_MEE_1N.{mp4,gif}` |
| 0.5 | 1375.28 | 362 | 138.60 / 149 | (anchor-free, footnoted) |

R0 law `T·t_{f,min} ≈ 850 N·h` holds to 0.72% across the four certified anchors.
0.2 / 0.1 N were not attained (see `TODO.md`). The Cartesian formulation died
at 5 N; the MEE + ΔL formulation broke that wall.
