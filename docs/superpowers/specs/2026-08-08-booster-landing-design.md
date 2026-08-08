# Booster Landing — 3-DOF Powered-Descent Guidance + Tracking (Design Spec)

**Date:** 2026-08-08
**Location:** `~/Desktop/optimal_control/booster_landing/`
**Status:** Approved design, pre-implementation

## Goal

Model a Falcon-9-class booster landing burn in MATLAB: solve the min-fuel
powered-descent guidance (PDG) problem two independent ways (direct collocation
and lossless convexification), certify they agree, then close the loop with a
time-varying LQR tracker and Monte-Carlo the landing under dispersions.
Campaign-style deliverable: one front door, `lib/ certify/ viz/ tests/ doc/
results/`.

## Scope decisions (locked)

| Decision | Choice |
|---|---|
| Model scope | 3-DOF point-mass landing burn (no attitude; 6-DOF and full return profile are future builds) |
| Architecture | Guidance (optimal trajectory) + tracking (TVLQR) + closed-loop Monte Carlo |
| Solvers | Both direct collocation (CasADi+IPOPT) and lossless-convexification, cross-validated |
| Fidelity | Phase 1 vacuum / flat Earth / constant gravity (convexification exactly valid); Phase 2 adds exponential atmosphere + quadratic drag to the collocation/tracking side |
| Format | Campaign-style code with front door `run_booster_landing.m` |
| Environment | MATLAB R2025b (`/Applications/MATLAB_R2025b.app/bin/matlab -batch`), CasADi+IPOPT already installed |

## Phase 1 problem statement (classic PDG, vacuum)

**State** (7): position r ∈ R³ (pad at origin, z up), velocity v ∈ R³, mass m.
**Control**: thrust vector T ∈ R³.

**Dynamics:**

```
ṙ = v
v̇ = g + T/m          g = [0; 0; -g0]
ṁ = -‖T‖/(Isp·g0)
```

**Constraints:**

- Thrust annulus (the famous nonconvex one): `0 < Tmin ≤ ‖T‖ ≤ Tmax`.
  Engine burns continuously — single burn, no on/off switching (classic PDG
  assumption; min-fuel throttle is then max–min–max bang-bang).
- Glideslope cone: `‖r_xy‖ ≤ z / tan(γ_gs)` with γ_gs the minimum elevation
  angle above horizon (default 30°, adjustable).
- Optional thrust-pointing cone: angle of T from vertical ≤ θ_max (default
  off; convex form `T_z ≥ cos(θ_max)·‖T‖`).
- Terminal: r(t_f) = 0, v(t_f) = 0. Mass floor m ≥ m_dry.
- Free final time t_f.

**Objective:** min fuel = max m(t_f) (equivalently min ∫‖T‖ dt).

**Nominal parameters** (public Falcon-9 estimates, all adjustable in
`booster_params.m`):

| Param | Value | Note |
|---|---|---|
| m_dry | 25 600 kg | F9 first-stage dry mass |
| m0 | 30 000 kg | dry + ~4.4 t landing propellant at burn start |
| Tmax | 845 kN | one Merlin 1D, sea level |
| Tmin | 0.40·Tmax = 338 kN | ~40% min throttle |
| Isp | 282 s | Merlin 1D sea level |
| g0 | 9.80665 m/s² | |
| r0 | [500; 100; 2000] m | post-entry-burn offset |
| v0 | [-30; 0; -180] m/s | descending ~180 m/s |
| γ_gs | 30° | glideslope elevation |

Physical payoff to verify: at landing mass, even min throttle gives T/W =
338k/(25.6 t·g0) ≈ 1.35 > 1 — the booster cannot hover, so the **hoverslam
emerges from the optimization** (terminal max-throttle arc arriving at zero
velocity exactly at the pad).

## Guidance solvers

### `solve_pdg_colloc` (house machinery)

- CasADi Opti, Hermite-Simpson collocation, time normalized to τ ∈ [0,1] with
  t_f a decision variable (house free-t_f pattern).
- Annulus handled directly: `Tmin² ≤ ‖T_k‖² ≤ Tmax²` (nonconvex NLP, IPOPT).
- Initial guess: straight-line r, linear v, m linear, T = mg (or warm-start
  from the convex solution — both paths supported; cold-start must work for
  the nominal case).
- N = 60 segments default.

### `solve_pdg_convex` (lossless convexification)

Açıkmeşe & Ploen 2007 / Blackmore–Açıkmeşe–Scharf 2010 change of variables:

```
u = T/m,  σ = Γ/m,  z = ln m
ṙ = v,  v̇ = g + u,  ż = -σ/(Isp·g0)
‖u‖ ≤ σ                          (relaxation; lossless at optimum)
ρ1·e^{-z̄(t)}[1 - (z-z̄) + (z-z̄)²/2] ≤ σ ≤ ρ2·e^{-z̄(t)}[1 - (z-z̄)]
```

with z̄(t) the max-thrust mass-depletion reference, ρ1 = Tmin, ρ2 = Tmax.
Objective: max z(t_f). All constraints convex (SOC + linear) for **fixed
t_f**; outer golden-section search on t_f.

- Discretization: same Hermite-Simpson grid, formulated in CasADi, solved with
  IPOPT. The problem is convex, so IPOPT's local optimum is global; a true
  conic solver (ECOS/CVX) is a documented drop-in, not part of this build.
- Losslessness check: ‖u‖ = σ on the optimal solution (within tolerance).

## Certification layer (`certify/`)

Gates, house style (report-only failures never silently pass):

- **G1 defects:** NLP defect norms at solution.
- **G2 continuous residual:** dense ode45 re-integration flying the
  interpolated control; position/velocity/mass error at t_f and along the
  path ("defect is not accuracy" lesson — measure the real residual).
- **G3 cross-method agreement:** |m_f(colloc) − m_f(convex)| < 0.1 kg;
  trajectory L∞ difference reported; t_f(free-NLP) vs t_f(golden-section)
  consistent.
- **G4 losslessness:** max |‖u‖ − σ| on the convex solution.
- **G5 PMP structure:** throttle profile is max–min–max (≤ 2 switches);
  primer vector (velocity costate from NLP duals — use `opti.lam_g` by row
  range, per the house sign-bug lesson) anti-parallel to thrust direction;
  switching function sign changes match throttle switches.

## Tracking layer

### `tvlqr_design`

- Linearize 7-state dynamics along (x*(t), u*(t)) → A(t), B(t) with control
  in thrust-vector components.
- Backward integration of the differential Riccati equation (ode45 on the
  vectorized P(t), terminal cost Q_f), gains K(t) stored on a dense time grid.
- Weights: position/velocity weighted for touchdown accuracy; mass state
  lightly weighted (observable but not directly controlled).

### `sim_closed_loop`

- ode45 truth sim: T_cmd(t) = T*(t) − K(t)·δx, magnitude clamped to
  [Tmin, Tmax] (saturation is part of the sim, not ignored).
- Terminates on z = 0 crossing (event function); reports touchdown position,
  velocity, propellant remaining.

### `run_monte_carlo`

Dispersions (defaults, adjustable): initial position ±100 m (1σ), velocity
±10 m/s, thrust-magnitude bias ±1.5%, Isp ±1%. N = 200 runs. Outputs:
landing-footprint scatter + 3σ ellipse, touchdown-speed histogram,
propellant-margin stats, success rate (inside pad radius 15 m, touchdown
speed < 2 m/s, propellant ≥ dry mass).

## Phase 2 — atmosphere (after Phase 1 is certified)

- Exponential atmosphere ρ(z) = ρ0·e^{−z/H} (ρ0 = 1.225 kg/m³, H = 8.5 km),
  quadratic drag D = −½ρ‖v‖v·Cd·A/m with F9 estimates (Cd ≈ 1.0 landing-leg
  config, A ≈ 10.75 m² for 3.7 m diameter).
- Collocation solver re-solves guidance with drag (opt-in flag in dynamics,
  house pattern from `lt_mee_rhs`); convex branch stays vacuum — the
  documented comparison quantifies what drag-free guidance misses.
- TVLQR re-linearizes with drag terms; Monte Carlo adds wind (constant +
  gust, ±10 m/s 1σ horizontal).
- Successive convexification (SCvx) is **future work**, out of scope.

## Layout

```
booster_landing/
├── run_booster_landing.m        % front door: params → solve both → certify →
│                                %   tvlqr → MC → plots/movie → results/
├── lib/
│   ├── booster_params.m         % all physical + solver params, one struct
│   ├── pdg_dynamics.m           % 3-DOF RHS, opt-in drag branch (Phase 2)
│   ├── solve_pdg_colloc.m
│   ├── solve_pdg_convex.m       % + golden-section t_f driver
│   ├── tvlqr_design.m
│   └── sim_closed_loop.m
├── certify/
│   └── certify_pdg.m            % gates G1–G5, single report struct
├── viz/
│   ├── plot_pdg_solution.m      % trajectory 3D + throttle + mass + primer
│   ├── plot_footprint.m         % MC scatter + ellipse
│   └── movie_landing.m          % house polished-graphics landing movie
├── tests/
│   └── test_booster_landing.m   % fast suite (below)
├── doc/
│   └── booster_landing_note.tex % theory note: PDG problem, convexification
│                                %   + losslessness proof sketch, TVLQR
└── results/                     % .mat solutions, PNGs, MP4, progress txt
```

## Conventions (Mike's standing rules apply)

- **Pumpkyn style throughout:** all MATLAB written per the house pumpkyn look
  (`%%`-delimited header blocks, column-aligned `=`, colon-terminated `%%`
  section comments, `if nargin==0` self-demos on library functions). Invoke
  the `matlab-pumpkyn-style` skill when writing each file.
- **Leverage pumpkyn where useful:** booster landing is not CR3BP, so the
  orbit-family machinery doesn't apply, but reuse pumpkyn/pumpkynPie
  utilities (plotting, integrators, param-struct patterns) wherever one fits
  rather than re-implementing. Document each reuse (or considered-and-refused
  reuse) in the note.
- **Run scripts Mike can run as easily as Claude:** `run_booster_landing`
  with no arguments does the complete nominal campaign from a fresh MATLAB —
  sets its own paths, prints progress, writes results/. Adjustable via an
  optional params/options argument, never by editing the file. Same contract
  as `run_cr3bp_geo` / `run_gto_tulip`.
- Standard MATLAB header block (purpose/inputs/outputs/references) on every
  function. No `i`/`j` loop variables.

## Test plan (`tests/`, fast — seconds to ~2 min)

1. Params sanity: T/W at min throttle > 1 at dry mass; m0 > m_dry.
2. Dynamics: complex-step vs analytic Jacobian of `pdg_dynamics`.
3. Coarse-grid smoke solve (N = 15) of both solvers converges.
4. Losslessness on the coarse convex solve.
5. Riccati: P(t) symmetric positive definite along the trajectory.
6. Closed loop, zero dispersion: tracks guidance to < 1 m / < 0.1 m/s at
   touchdown.
7. Glideslope satisfied on nominal solution (no violation > 1e-6).

The flagship reproduction (front door, full N, both solvers, gates green) is
the real acceptance test — re-run after any structural change (house lesson:
fast tests passing does not prove the campaign survived).

## Success criteria

- Both solvers converge on the nominal case; G1–G5 all green; final masses
  agree < 0.1 kg.
- Hoverslam structure visible: max–min–max throttle with terminal
  max-throttle arc to zero velocity at the pad.
- Monte Carlo ≥ 95% success at default dispersions.
- Landing movie renders clean (÷16 frame size, house H.264 lesson).
- Theory note compiles; aux files cleaned.

## Risks / mitigations

- **Nonconvex annulus needs a good initial guess** → straight-line init must
  work nominally; convex warm-start is the documented fallback.
- **Free-t_f scaling** → τ-normalized time, t_f bounded [5, 60] s.
- **Convexified mass bounds are Taylor approximations** → G3 agreement gate
  catches meaningful error; grid refinement documented if it bites.
- **TVLQR saturation near the hoverslam terminal arc** (T* at Tmax leaves no
  up-authority) → expected physics, shown in MC results, discussed in note;
  not "fixed" by weight tuning tricks.

## Future work (explicitly out of scope)

6-DOF rigid body (gimbal + grid fins), full return profile (boostback/entry),
successive convexification, MPC-in-the-loop guidance, real conic solver.
