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
> low-thrust work is a *different* campaign, `../GTO_tulip/`.) Full
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
PMP-steered switch-aware mesh refinement, ported from `../GTO_tulip/PSR/`.
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
2000 s (the benchmark's exact value is **1994.8 s** — Caillau & Noailles 2001,
`orbit_transfer/min_fuel_papers/COCV_2001__6__239_0.pdf` p.255, δ=0.05112 km⁻¹s ⇒ Isp=1/(δg₀);
our default is 0.27% high). Single source of physical constants.

### `elements_to_cart.m` / `cart_to_elements.m`
Algebraic MEE ↔ inertial `(r,v)` maps (roundtrip-tested), used by the seed,
the reconstruction check, and the viz adapter.

### PMP verifier — `verify_pmp_mee.m` (+ `mee_dual_to_costate.m`, `mee_primer_switch.m`)
First-order Pontryagin certificate from the NLP's KKT duals (primer alignment,
switching-sign law). Diagnostic only; **primal certifications (defect/terminal)
never depend on it.** See `process/CAMPAIGN.md` footnote 5 for the open dual anomaly.

### SOSC certificate — `verify/sosc/` (NLP-level second-order local-min test)
A rigorous **second-order** certificate for a saved min-fuel row: it warm-re-solves
the NLP to recover a machine-tight KKT point, then classifies the reduced Hessian
on the critical cone via a direct null-space eig (`Z=null(A)`, `RH=Z'HZ`, `eig`) with
a `zt`-sensitivity gate. Verdicts: **PASS** (strict local min) / **WEAK_MIN** (PSD,
flat directions — no descent direction, not strict) / **FAIL** (proven saddle) /
**INCONCLUSIVE** (sign not resolvable, or too large to compute). Entry point
`verify_sosc_mee(<row.mat>)`; batch `recertify_table3([10 5 2.5 1 0.5])` writes
sidecar verdicts to `results/sosc/` (campaign caches untouched). **Key result:** the
10 N certified row is a **WEAK_MIN with 270 flat directions** — min-fuel bang-bang
extremals are *weak* (non-strict) minima, so strict SOSC is generically unreachable.
Gated into the driver/reproducer opt-in via `cfg.certifySosc` (default off); only a
FAIL verdict demotes a certified row. Design + method evolution + per-rung results:
`process/DESIGN_sosc.md` (esp. §12) and `process/PLAN_sosc.md`.

### Hamiltonian verification — `verify/hamiltonian_const_check.m`, `verify/hamiltonian_along_traj.m`
Independent **first-order (PMP)** checks from the recovered costates, complementary
to the second-order SOSC test. Because the L-domain carries **time `t` as a state**
(true longitude is the independent variable), the dynamics are autonomous in `t`, so
the time-costate `λ_t` is a first integral and `λ_t = −H_t` (energy–time conjugacy).
- **`hamiltonian_const_check(<row.mat>)`** — verifies `λ_t` (row-7 defect dual) is
  **constant** ⇔ the time-domain Hamiltonian `H_t` is conserved. KKT stationarity
  forces this exactly at a true KKT point; measured **CoV ~1e-9 (0.1 N) / 2e-8 (10 N)**,
  while a mass-costate control varies **56%** — so the test genuinely discriminates.
- **`hamiltonian_along_traj(<row.mat>)`** — reconstructs BOTH Hamiltonians node-by-node
  (`H_L = ℓ_L + λᵀ·dXdL`, `H_t = −λ_t`). Two subtleties, both handled: (1) the CasADi
  dual **sign** is pinned by transversality `H_t = dJ*/dt_f < 0` for min-fuel above
  min-time (the extremal identity does *not* discriminate sign); (2) the extremal
  identity `dH_L/dσ = ∂H_L/∂σ` (complex-step `∂/∂L`) is reported as a reconstruction
  check that holds only to **collocation order** (~1.6% at 26 nodes/rev), *not* a
  machine-precision claim. `H_L` breathes once per orbit (L explicit in the dynamics);
  `H_t` is flat.

### Movie / viz — `viz/hamiltonian_movie.m`, `viz/mee_res_to_cart_res.m`
`hamiltonian_movie(H, stem)` renders the two-panel synced animation (`H_L` breathing
vs `H_t` conserved, on one shared y-scale). The Cartesian adapter
`mee_res_to_cart_res(..., nDense)` gained an opt-in **render-densification** factor
(default `1` = byte-identical): the coarse mesh (~8 nodes/rev at deep rungs) draws each
orbit as a **polygon** (straight segments between nodes — a rendering artifact, not
physics); `nDense>1` pchip-interpolates the slow elements onto a fine σ grid and
re-maps to Cartesian at the fine longitudes, recovering the smooth spiral (throttle
held piecewise-constant so the switch count is preserved).

### Cartesian counterparts (original reproduction)
`casadi_lt_2body.m`, `run_transfer.m`, `run_mintime.m`, `homotopy_2body.m`,
`seed_2body.m`, `verify_pmp_2body.m` — the Sundman-regularized Cartesian solver
that reproduced the paper's 10 N case but stalled at 5 N; the MEE stack above
superseded it. Kept for the cross-formulation gate (see `process/CAMPAIGN.md`).

---

## Code layout

Code is organized into functional subfolders: `core/` `drivers/` `psr/`
`verify/` (with a nested `verify/sosc/` for the second-order certificate)
`frontdoor/` `reproduce/` `viz/` `coords/` `cartesian_legacy/` `lib/`
`tests/` `attic/`. At the module root sit only the two front-door docs
(`README.md`, `TODO.md`) and the two path helpers (`setup_paths.m`,
`module_root.m`). **Run `setup_paths` once per session** (after `cd`-ing into
this directory) before calling anything — it puts every subfolder on the MATLAB
path so functions resolve regardless of which subfolder they live in.

Documentation lives in two folders, split by maturity:
- **`doc/`** — polished deliverables: the method note and the campaign
  reproduction runbook (`.tex` + rendered `.pdf`) and external `reviews/`.
- **`process/`** — the working development record: the design specs and
  implementation plans (`DESIGN*.md`, `PLAN*.md`) and the campaign narrative
  with its honesty footnotes (`CAMPAIGN.md`).

## Usage

```matlab
cd earth_elliptic_to_geo
setup_paths                                           % put all subfolders on the path
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'))      % CasADi on path

% one Table-3 row, default paper endpoints, reuse the certified cache:
row = run_gergaud(struct('thrustN',10,'runMode','auto'));   % prints the row; row struct returned

% custom initial orbit, live solve + plot + movie:
run_gergaud(struct('thrustN',2.5,'runMode','solve','e0',0.6,'i0_deg',10, ...
                   'makePlot',true,'makeMovie',true));

% the whole ladder as a batch:
run_ladder([10 5 2.5 1])
```

Tests live in `tests/` (fast, `matlab -batch`): `test_mee_xf`,
`test_mee_seed_initelems`, `test_mee_threading`, `test_mee_res_to_cart`,
`test_gergaud_row`, `test_run_gergaud_auto` (front-door suite); plus the
campaign's no-solve guards `test_params/elements/dynamics/terminal/seed/
mee_rhs/mee_seed/...`. Each test locates the module root and calls
`setup_paths` itself, so it can be run standalone, e.g.
`matlab -batch "run('/abs/path/earth_elliptic_to_geo/tests/test_mee_xf.m')"`.

## Reproducing from scratch (best-found)

`reproduce_row.m` is a **keep-best-mass reproducer ENGINE**, not a bit-exact
replay: for one thrust rung `T` it re-solves the min-fuel transfer entirely
FROM SCRATCH (own `REPRO_`-prefixed tags, so it can never load or clobber the
campaign's production caches), then verifies the result against the campaign
floor. Because minimum-fuel means *maximize final mass*, its fuel stage runs a
multi-start (a seed set at the exact min-time-anchor `t_f`, falling back to a
tiny `t_f`-bracket) and keeps whichever certified candidate has the highest
final mass — the fuel bang-bang basin is razor-sensitive to `t_f` (a ~2e-5
change in `t_f` can flip the whole switch structure), so a single solve is not
trusted.

- **Entry points**, in order of preference for the deep/crash-prone rungs:
  - `reproduce/reproduce_table3.sh` — per-*process* watchdog: one MATLAB
    process per rung, relaunch-on-crash-or-hang (up to a per-rung attempt
    cap), because a CasADi/MUMPS MEX-level crash is not catchable from inside
    MATLAB. Usage (from `reproduce/`): `./reproduce_table3.sh` (default rungs
    `10 5 2.5 1 0.5`) or `./reproduce_table3.sh 10 5`; logs to
    `results/repro/reproduce_table3.log` (module root, resolved via
    `module_root()`).
  - `reproduce_table3.m(thrustList)` — thin **in-process** wrapper (`for T =
    thrustList, reproduce_row(T); end`, then prints the table). Convenient for
    the crash-free top rungs (10/5/2.5 N) or quick iteration; a fatal crash on
    a deep rung takes the whole process down with no auto-relaunch, so prefer
    the watchdog for 1 N and below.
  - `reproduce_row(T)` — the single-rung engine itself, callable directly for
    one row.
- **`results/repro/` namespace** — every rung writes
  `results/repro/REPRO_row_T<round(10*T)>.mat` (variables `row`, `anchor`,
  `sol`, `rep`; e.g. `REPRO_row_T100.mat` for 10 N, `REPRO_row_T5.mat` for
  0.5 N), kept separate from both the campaign's own `results/*.mat` caches
  and the driver-internal `REPRO_`-tagged per-stage cache files. `reproduce_row`
  also builds a `results/repro/` directory if absent, and `chain`-strategy
  rungs must be reproduced in order (`load_prev` errors loudly if a rung's
  predecessor hasn't been reproduced yet).
- **Recipe registry** — `table3_recipes.m` is a pure lookup: for each of the
  seven Table-3 rungs (10, 5, 2.5, 1, 0.5, 0.2, 0.1 N) it returns the exact
  proven anchor strategy (`coldB` cold min-time solve / `chain` warm-started
  from the previous rung / `smallN_first`, the 1 N low-node-density-then-
  mesh-refine anchor / `R0law`, the anchor-free `t_{f,min} ≈ 223.14/T`
  estimate), fuel-stage node density + seed throttle + warm-start source, and
  an optional PSR (post-solution mesh refinement) pass. 0.2 N and 0.1 N have
  registered recipes (`.seeded = true`) but have not been run to a certified
  row in this build.
- **One-sided verify, not a bit-exact match** — `verify_row.m` throws only if
  the reproduced final mass falls below `table3_certified(T).m_f_kg` minus a
  0.5 kg numerical-noise slack. A HIGHER mass always passes and is reported as
  an improvement; switch count and revolution count are reported for
  comparison but never gated, because the best min-fuel optimum can have a
  different (better) bang-bang structure than the campaign's row.
- **Table 3 gets updated with best-found numbers, not just replayed** — as
  each rung is re-run, `reproduce_table3_collect.m` prints the *updated*
  table (row + a BEAT/matched verdict against the campaign floor). This has
  already happened at 10 N: the engine independently found **18 switches /
  7.56 rev / 1378.46 kg**, beating the campaign's certified 1377.10 kg by
  1.36 kg (and closer to the paper's ~18-switch structure) — the campaign had
  under-optimized that rung, and the reproducer is expected to equal or beat
  it at every thrust level. The 5/2.5/1/0.5 N rungs have registered recipes
  and are expected to meet-or-beat their campaign floors the same way, but
  that is validated only when the full ladder is actually run (see `TODO.md`).

```matlab
reproduce_row(10)                    % one rung, live, verified against the floor
reproduce_table3([10 5 2.5])         % in-process, crash-free top rungs
```
```bash
cd reproduce && ./reproduce_table3.sh 10 5 2.5 1 0.5   # per-process watchdog, survives MEX crashes
```

> **Reproducing the ladder from scratch?** `process/LADDER_FROM_SCRATCH.md` is the
> one-page operational recipe — exact commands per rung, the full certified-number
> table (10→0.1 N), per-rung anchor strategy, the deep-rung four levers, and the
> gotchas. The canonical numbers themselves live in `reproduce/table3_certified.m`
> (all seven rungs).

## Related files (outside this directory)

- **Method note (full methodology + flow diagram):** `doc/table3_method_note.tex`
- **Campaign narrative + results tables + honesty footnotes:** `process/CAMPAIGN.md`
- **Outstanding work:** `TODO.md`
- **CasADi/IPOPT:** `~/casadi-3.7.0`
- **Parent campaign (CR3BP, shared solver architecture):** `../GTO_tulip/`
- **Problem source:** Haberkorn, Martinon & Gergaud, JGCD 27(6), 2004
  (`orbit_transfer/min_fuel_papers/Gergaud-Haberkorn-Martinon-JournalGuidance2004-preprint.pdf`)

## Results (headline)

Certified min-fuel ladder at c_tf = 1.5 (full table + honesty footnotes in
`process/CAMPAIGN.md`):

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
