# Drone dynamics library, multi-agent path planning (two routes) and a quadrotor — design (2026-10-23)

Reproduce the multi-agent "drone" trajectory problem of Zhang, Liu, Darbon
and Karniadakis, *A Transformer-based Approach for Real-time Multi-agent Path
Planning with Varying Initial/Final Conditions* (arXiv:2610.30388v1,
`papers/2610.30388v1.pdf`), with the classical machinery this repository
already holds, certify what the paper only approximates, and then put a real
vehicle under the same task. Four stages, in order, each a working tested
deliverable on its own:

0. **The 3D drone dynamics library** (Mike, 2026-10-23: "a helpful first
   step"). Three 3D models of increasing fidelity, one home for each one's
   constants, field, Jacobian and an independent physics oracle, nothing
   objective-specific in it. The paper's drone is the lowest rung; the
   quadrotor is the top.
1. **Route 1, direct.** The paper's problem (its eq. 19) on the kinematic
   model as a CasADi/IPOPT trapezoidal NLP with HARD pairwise clearance,
   plus the paper's own log-barrier form so Route 2 has a twin with the same
   objective. Reproduces their Table 1 (M = 4 .. 32 required, 64 and 128
   stretch).
2. **Route 2, indirect.** The barrier form's Pontryagin boundary-value
   problem (their eq. 7) shot with `oc.ms_bvp` at fixed final time, seeded
   from Route 1, and certified with `oc.ms_conjugate_test`.
3. **Quadrotor upgrade.** The same path-planning task on the rigid-body
   model: direct for M vehicles with clearance, indirect for one vehicle,
   both certified.

Each stage gets its own implementation plan. Plan 0 is written in full now;
plans 1-3 are written now too and re-checked when the stage before them
lands (a plan's seed route depends on what the previous stage taught).

Spec authority: this document. The paper is the source of the problem
statement and the numbers to compare against, not of the method.

## Stage 0 — The 3D drone dynamics library

`optimal_control_examples/drone/drone_dynamics/`, the `cartpole_common`
pattern generalised to a family. Every model is written in the control-affine
form `ẋ = F(x) + B(x) u` so one PMP field and one direct transcription
template serve all three, and every model carries its own oracle.

| model | state | control | what it is |
|---|---|---|---|
| `kinematic` | `p ∈ R³` | `v ∈ R³` velocity | the paper's drone: `ṗ = v`. The 2D problem is this model with `v_z` pinned to zero by bounds, not a separate model |
| `pointmass` | `[p; v] ∈ R⁶` | `a ∈ R³` thrust acceleration | the planner's quadrotor: `ṗ = v`, `v̇ = a − g e_3`; the model every differential-flatness planner uses |
| `rigid6dof` | `[p; v; Θ; ω] ∈ R¹²` | `f ∈ R⁴` rotor thrusts | the vehicle: 6-DOF rigid body with ZYX Euler attitude and a "+" rotor mixer |

**Interface, identical across models.**

```
p              = drone_params(model)                 constants, with .model and .nx .nu
[F, B]         = drone_field(x, p)                   F [nx x 1], B [nx x nu]; xdot = F + B u
A              = drone_state_jac(x, u, p)            d(F + B u)/dx [nx x nx], generated
[R, W, dRdTh]  = euler_zyx(Theta)                    rotation matrix, body-rate map, dR/dTheta {3}
```

`drone_field` is complex-step safe (no `abs`, `norm`, `max`, `min`, no
branch on a state value), so the generated Jacobians are tested against a
complex step and against central differences, as `cartpole_common` does.
`drone_state_jac` dispatches on `p.model` to one generated file per model
(`drone_state_jac_<model>_gen.m`, Symbolic Math Toolbox + `matlabFunction`),
and the generator `gen_drone_jacs.m` calls `drone_field` itself on symbolic
inputs so there is exactly one copy of the physics.

**The rigid-body model.** Crazyflie-class constants (Luis and Le Ny,
*Design of a Trajectory Tracking Controller for a Nanoquadcopter*, 2016,
Table 1; Förster, *System Identification of the Crazyflie 2.0 Nano
Quadrocopter*, ETH 2015):

```
m = 0.027 kg,  I = diag(1.395e-5, 1.436e-5, 2.173e-5) kg m²,
l = 0.046 m (arm),  kappa = 0.005964 m (rotor torque / thrust),
g = 9.81 m/s²,  fMax = 0.15 N per rotor   (hover: m g / 4 = 0.0662 N each)
```

State `x = [p; v; Θ; ω]`: inertial position and velocity, ZYX Euler angles
`Θ = [φ; θ; ψ]` (roll, pitch, yaw), body rates. Mixer (rotor 1 on +x_b,
2 on +y_b, 3 on −x_b, 4 on −y_b):

```
T   = f_1 + f_2 + f_3 + f_4
tau = [ l (f_2 − f_4);  l (f_3 − f_1);  kappa (f_1 − f_2 + f_3 − f_4) ]

ṗ = v
v̇ = −g e_3 + (T/m) R(Θ) e_3
Θ̇ = W(Θ) ω,   W = [1, sinφ tanθ, cosφ tanθ;  0, cosφ, −sinφ;  0, sinφ/cosθ, cosφ/cosθ]
ω̇ = I⁻¹ (tau − ω × I ω)
R = R_z(ψ) R_y(θ) R_x(φ)
```

Euler angles rather than a quaternion so the indirect BVP of Stage 3 stays
square (12 costates, no norm constraint, no gauge degeneracy for the
conjugate test). The singularity at `θ = ±90°` is a documented limit:
`drone_params('rigid6dof').thetaMax = 60°` is the box every direct solve
imposes and reports.

**The pointmass model.** `p = drone_params('pointmass')` carries `g` and
`aMax = fMax·4/m − g` (the rigid body's thrust-to-weight margin, so a
pointmass solution is one the vehicle could in principle fly).

**The physics oracle, `tests/test_drone_physics.m`.** Derived from mechanics
and geometry, never from the implemented rows (the house rule, paid for on
2026-09-17 when a cart-pole sign error survived nine self-comparing reviews):

1. `kinematic`: `F = 0`, `B = I`, exactly; the model is its own oracle and
   is tested only for shape and linearity.
2. `pointmass`: free fall with `a = 0` has `v̇ = −g e_3`; hover `a = g e_3`
   is an equilibrium of `v`; power balance `d/dt(½|v|² + g p_z) = a·v`.
3. `rigid6dof` hover: `f_k = m g/4`, `Θ = 0`, `ω = 0`, `v = 0` is an
   equilibrium (`F + B f = 0`).
4. `rigid6dof` free fall: `f = 0` gives `v̇ = −g e_3` at every attitude.
5. `rigid6dof` power balance, pointwise at random states and controls:
   `d/dt (½ m |v|² + m g p_z + ½ ω'Iω) = T (R e_3)·v + tau·ω`.
   The gyroscopic term does no work; a sign error in any row breaks this.
6. `rigid6dof` kinematic consistency, pointwise: `Ṙ = R [ω]_×` with
   `Ṙ = Σ_k ∂R/∂Θ_k Θ̇_k` (pins `W` independently of the rates row).
7. `rigid6dof` torque-free flight (`tau = 0`, any `T`) conserves the
   inertial angular momentum `R I ω` to integrator tolerance; a hover-trim
   flight stays put.
8. `euler_zyx`: `R'R = I`, `det R = 1`, `R(0) = I`, `R` matches the product
   of the three elementary rotations written out longhand, and the complex-
   step derivative of `R` matches `dRdTh`.

**Also in the library (geometry the problems share).**

```
[x0, xT, geom] = drone_layouts(which, M, r)        the paper's layouts + ours, as [3 x M] (z = 0 for 2D)
pairs          = pair_list(M)                      all pairs [2 x P]
[h, rho, n]    = clearance_terms(P, pairs, geom)   h = rho − 2r per pair, for one configuration P [3 x M]
[U, gU, HU]    = barrier_terms(P, pairs, geom)     the paper's C¹ log barrier summed over pairs (+ ball), grad, sparse Hessian
JLB            = lower_bound(x0, xT, T)            Σ |xT − x0|² / (2T)
```

`barrier_terms` is complex-step safe (branch on the real part of `h`), and
its analytic gradient and Hessian are tested against a complex step of `U`
and of `gU`. The barrier shape, from the paper's eq. 3b with `h = ρ − 2r`:

```
U(h) = −ε log(h)                                    if h > l
     = −ε log(l) + (ε/2) (((h − 2l)/l)² − 1)        if h ≤ l         (C¹ at h = l)
U'   = −ε/h                      |  (ε/l²)(h − 2l)
U''  =  ε/h²                     |   ε/l²                             (continuous at l too)
```

and for one pair with `n = (x_i − x_j)/ρ`:

```
∇_{x_i} U =  U'(h) n,          ∇_{x_j} U = −U'(h) n
∇²_{ii}  =  U''(h) n n' + U'(h) (I − n n')/ρ  = ∇²_{jj},   ∇²_{ij} = ∇²_{ji} = −∇²_{ii}
```

The ball obstacle `h = |x_i − c| − R_ball` enters the same sum with the
`i`-block only. Defaults `ε = 1e-5`, `l = 1e-4` (the paper's).

**Rule for the folder.** Objective-independent only, as `cartpole_common`:
no cost, horizon, boundary condition or solver anywhere in it. The barrier
is a penalty SHAPE on geometry; the objective that integrates it lives in
the routes.

## What the paper's problem is

M agents, single integrators in `d = 2` or `3` dimensions:

| | |
|---|---|
| state / control | `x_i ∈ R^d`, `ẋ_i = v_i` |
| cost | `J = ∫_0^T Σ_i ½|v_i|² dt`, **fixed** `T = 10` |
| boundary | `x_i(0) = x_i^0`, `x_i(T) = x_i^T` (points) |
| clearance | `|x_i − x_j| ≥ 2r` for all pairs, `r = 0.02` (`0.01` for M = 64, 128) |
| 2D layouts | M ≤ 16: one circle of radius 0.5, targets antipodal; M = 32: two circles (16 each, radii 0.33, 0.5); M = 64: two circles of 32 (0.33, 0.5); M = 128: four circles of 32 (0.35, 0.40, 0.45, 0.50) |
| 3D (ball) | obstacle `|x_i − c| ≥ R_ball`, `c = 0`, `R_ball = 0.1`; layouts NOT given in the paper |
| 3D (room) | box walls plus two box obstacles "as specified in the code release", NOT given |
| lower bound | `J_LB = Σ_i |x_i^T − x_i^0|² / (2T)` |

Their reported costs sit 0.2 % (M = 4) to 9.3 % (M = 128) above `J_LB`;
their "NLP ground truth" exists only for M ≤ 8. Their penalty is the barrier
above taken as the **max** over constraints at each instant (their eq. 3a).

**Scope decisions, from the paper's gaps.**
- The 2D Table 1 is the reproduction target. The 3D ball case is kept as a
  demonstration of obstacle handling with a layout WE choose (M = 16 on a
  sphere of radius 0.4, antipodal targets, documented as ours). The 3D room
  with box obstacles is dropped: its geometry is not in the paper.
- The barrier is a **sum** over pair terms, not the paper's max. The max
  makes the Hamiltonian field only piecewise smooth, which a PINN tolerates
  and a shooting method does not. The sum is a different penalised problem
  by at most `N_pairs · ε · |log h_typ|`, about 0.1 % of J at M = 16, and is
  reported beside the hard-constraint optimum, which is the reference in
  every table.

## Stage 1 — Route 1, direct (CasADi + IPOPT) on `kinematic`

**Transcription.** Uniform mesh of N intervals (default 100), `h = T/N`.
Decision variables `X ∈ R^{3M × (N+1)}` (agent-major rows) and `V` likewise.

```
defects      X(:,k+1) − X(:,k) − (h/2)(V(:,k) + V(:,k+1)) = 0,   k = 1..N
boundary     X(:,1) = x0,   X(:,N+1) = xT
2D           V(3i,:) = 0 for every agent (bounds), x0, xT have z = 0
cost         J = Σ_k (h/2)(L_k + L_{k+1}),   L_k = ½ Σ_i |V_i,k|²
clearance    |X_i,k − X_j,k|² − (2r)² ≥ 0   every pair, every node AND every
             interval midpoint of the cubic-Hermite interpolant,
             x_mid = (x_k + x_{k+1})/2 + (h/8)(v_k − v_{k+1})               (mode 'hard')
             (squared form: no square root, no singularity at coincidence)
ball         |X_i,k − c|² − R_ball² ≥ 0                                   (3D ball only)
```

A sampled path constraint holds between samples only to second order: on
a constraint arc of radius `2r` with sample spacing `s` (in space) the
interpolant cuts inside by up to `s²/(16r)` (the sagitta). That is `2e-4`
at N = 100 on the head-on case and cannot be driven to 1e-8 by any nodal
scheme, so the clearance gate REPORTS the fine-grid minimum against a
mesh-aware bar and the cost accuracy is carried by the closed-form oracle
and the mesh gate, where it belongs.

Mode `'barrier'` drops the inequalities and adds
`Σ_k (h/2)(U_k + U_{k+1})` with `U_k = barrier_terms(X(:,k))` to the cost.
`'hard'` is the reference; `'barrier'` seeds Route 2 and measures the bias.

**Initial guess and symmetry.** Straight lines plus a tangential bump
`A sin(π t/T) R_90 (x_i^T − x_i^0)/|x_i^T − x_i^0|`, `A = 0.05`, common
rotation sense (about `e_3`), so all agents swirl the same way. Every circle
layout is degenerate at the straight line (all agents meet at the centre at
`T/2`), and the mirror image of any solution is a solution of equal cost;
the common sense picks one. The mirror is solved too and must agree in cost
to 1e-6 relative (gate G5).

**Pair pruning (M ≥ 64, stretch only).** Constraint generation: start from
the pairs whose straight-line paths come within `4r`, solve, evaluate every
pair on the solution at a fine grid, add every pair closer than `2r + 1e-3`,
re-solve, at most five rounds. On the concentric-ring layouts all lines pass
the centre at `T/2`, so pruning is expected to keep nearly all pairs; it is
kept because it is the only route to M = 128 that does not start from
821,000 constraints. M = 64 and 128 run from a batch script under the
campaign discipline, never inside the test suite.

**Gates.**

| gate | what | bar |
|---|---|---|
| G1 clearance | min over pairs on a 10× finer cubic-Hermite reconstruction from (X, V) of `|x_i − x_j| − 2r` | ≥ −1e-2 · 2r ('hard', the sagitta bound at N ≥ 100; the number is printed in every table); reported ('barrier') |
| G2 bound | `J ≥ J_LB` | yes; `J − J_LB` reported |
| G3 solver | IPOPT status in {Solve_Succeeded, Solved_To_Acceptable_Level}; defects re-evaluated from the returned numbers | status ok, max defect < 1e-9 |
| G4 mesh | re-solve at 2N from the N solution | cost change < 1e-3 relative (M ≤ 16); reported above |
| G5 mirror | mirrored swirl's cost | 1e-6 relative |
| G6 oracle | two-agent head-on case vs the closed form below | 5e-4 relative at N = 100, 5e-5 at N = 400 (the discrete optimum cuts the arc's corner, so it sits BELOW J* by O(h²) and converges from below) |

**The closed-form oracle (G6).** `x_1: (−a,0,0) → (a,0,0)`,
`x_2: (a,0,0) → (−a,0,0)`, 2D. In centre-of-mass `c` and relative
`w = x_1 − x_2` coordinates the cost is `∫ |ċ|² + ¼|ẇ|² dt`; `c` is
optimally at rest; `w` travels from `(−2a,0)` to `(2a,0)` outside the disk
`|w| < 2r` in time T, so it follows the shortest such path at constant
speed. With `D = 2a`, `R = 2r`:

```
ℓ  = 2 sqrt(D² − R²) + R (π − 2 acos(R/D))
J* = ℓ² / (4T)
```

At `a = 0.4`, `r = 0.02`, `T = 10`: `ℓ = 1.602001`, `J* = 0.0641602`
(`J_LB = 0.064`). Physics the NLP does not know about: the one test that
cannot pass by the code agreeing with itself.

**Deliverables.** `solve_drones_direct` (one solve, both modes, 2D and 3D
ball), `verify_drones_direct` (G1-G5), `drone_table1` (the Table 1 loop
writing a tracked CSV beside the paper's numbers, committing the M = 4 and
8 fixtures), `plot_drones`, `movie_drones`, tests.

## Stage 2 — Route 2, indirect (`oc.ms_bvp`, fixed t_f) on `kinematic`

**PMP of the barrier form.** `H = ½|v|² + U(x) + λ'v` (minimising
convention, the cart-pole's):

```
∂H/∂v = v + λ = 0  →  v* = −λ       (∂²H/∂v² = I > 0: a minimum everywhere)
ẋ = −λ,    λ̇ = −∇U(x),    A = [0, −I; −∇²U(x), 0]
```

A `6M`-state autonomous Hamiltonian flow (2D: the z rows ride along at
zero). Unknowns `λ(0)`, conditions `x(T) = x^T`: square, smooth where
`h > l`, fixed time: `oc.ms_bvp` with `fixedTf`, `freeIdx0 = 3M+1 : 6M`.
`H` is constant along any trajectory of the flow (not zero: T is fixed), so
its spread measures integration quality, not optimality.

**Propagator.** `drone_pmp_prop(dt, y0, needSTM, geom)`: `ode113` at
1e-12/1e-14 with the variational equations alongside when `needSTM`
(`(6M)²` extra equations: 9,216 at M = 16), the cart-pole propagator's
collapse contract (throw on a truncated grid or non-finite result). The
Jacobian is the analytic `A` from `barrier_terms`, not a complex step.

**Seed.** The direct barrier solve at the same `(ε, l)`. Two seeds, compared
(a test of the covector rule on a new transcription, gate X1):
- `λ(t) = −V(t)` from the direct solution (exact at the optimum);
- `oc.duals_to_costates('trapezoid', mu = defect multipliers from
  opti.lam_g on the defect rows)`, sign resolved by correlation with `−V`.
The first seeds the engine.

**ε ladder (fallback only).** If the engine fails from the ε = 1e-5 seed,
walk ε from 1e-3 to 1e-5 in decades, seeding each rung from the previous
root, and report the rung that was needed.

**Certification.** `oc.ms_conjugate_test` with `stateRows = 1:3M`,
`costateCols = 3M+1:6M`, `quotientDir = []`, `freeTime = false`, `keepSTMs`
on the solve. The verdict is reported, never asserted PASS by a test unless
the measured root passes; the test asserts "not FAIL" on the M = 4 circle.

**Gates.**

| gate | what | bar |
|---|---|---|
| N1 converged | `info.converged`, residual | `|R| < 1e-10` |
| N2 terminal | re-flight from `λ(0)` at RelTol 2.5e-14 | miss < 1e-9 |
| N3 flown | `oc.fly_control` with the closed-form `v = −λ(t)` | miss < 1e-6 |
| N4 Hamiltonian | spread of H, relative | < 1e-9 |
| N5 adjoint | `λ̇` vs complex step of H in x | < 1e-9 relative |
| N6 clearance | min `h` along the arc vs `l` | reported; `> l` expected |
| X1 seed | the two seeds' agreement | < 1e-2 relative RMS |
| X2 direct | `J_indirect` vs `J_direct,barrier` | < 1e-3 relative |
| X3 segments | same `λ(0)` from K = 8 and K = 16 | < 1e-8 |
| S1 Legendre | `∂²H/∂v² = I` | identically, printed |
| S2 Jacobi | conjugate test verdict | not FAIL on M = 4 |

**Deliverables.** `drone_pmp_rhs`, `drone_pmp_prop`, `run_drones_pmp(M)`
(front door with `.ok/.why`), `drone_minenergy_study.m` (the numbered
study script in the cart-pole's shape), `results/route_compare.csv` with
columns `M, J_LB, J_hard, J_barrier_direct, J_indirect, OCFormer, SympOCNet,
conj`, tests.

## Stage 3 — Quadrotor on `rigid6dof`

- Direct, M vehicles: cost `J = ∫ Σ_i Σ_k f_{i,k}² dt`, fixed `T = 4 s`
  (the paper's 10 s is for unit-less kinematics; at 10 s a quadrotor barely
  leaves hover), hover-to-hover boundary (`v = 0`, `Θ = 0`, `ω = 0` both
  ends), clearance `|p_i − p_j| ≥ 2r` with `r = 0.06 m` (Crazyflie span plus
  margin), rotor bounds `[0, fMax]`, angle box `|φ|, |θ| ≤ thetaMax`.
  Layouts: M = 1 point-to-point `(−0.4, 0, 0.5) → (0.4, 0, 0.5)`; M = 4 on
  a circle of radius 0.4 at `z = 0.5`, antipodal. Trapezoid, CasADi/IPOPT,
  the Route 1 gates plus G1-local (`oc.local_residual` on the returned
  nodes with the trapezoid control reconstruction; bar 1e-6 per interval)
  and a bounds-margin report (rotor and angle boxes inactive on every
  shipped solution, or the solution is labelled bounded).
- Indirect, M = 1: `H = |f|² + λ'(F + B f)`, `f* = −B'λ/2` (unconstrained;
  fair only when the direct solution's bounds are inactive, which is
  reported), `λ̇ = −A(x, f*)'λ` with `A = drone_state_jac`. 24-state flow,
  `oc.ms_bvp` fixed-tf with `freeIdx0 = 13:24`, seed from the direct
  trapezoid duals via `oc.duals_to_costates('trapezoid')` with the sign
  resolved by correlation against the direct control (no thrust vote: the
  cart-pole's scalar-force pattern, applied to four rotors), conjugate test
  with `stateRows = 1:12`. The Route 2 gate table applies with `v` read as
  `f`.
- Also shipped: `pointmass` direct, M = 4, same layout, as the middle rung
  (no indirect, no new gates: one table row that shows what the attitude
  dynamics cost).
- Stretch, not planned: M = 4 indirect through the barrier.

**Deliverables.** `solve_quads_direct`, `run_quad_pmp`, `quad_pmp_rhs`,
`quad_pmp_prop`, `quad_minenergy_study.m`, `movie_quads`,
`results/quad_compare.csv`, tests.

## Layout

```
optimal_control_examples/drone/
  README.md                       what, how to run, the result tables, the caveats
  run_tests.m                     the folder's runner (registered in run_all_tests)
  drone_dynamics/                 STAGE 0: the library (objective-independent)
    drone_params.m  drone_field.m  drone_state_jac.m  gen_drone_jacs.m
    drone_state_jac_kinematic_gen.m  drone_state_jac_pointmass_gen.m  drone_state_jac_rigid6dof_gen.m
    euler_zyx.m
    drone_layouts.m  pair_list.m  clearance_terms.m  barrier_terms.m  lower_bound.m  hermite_resample.m
    README.md
    tests/  test_drone_physics (the oracle)  test_euler_zyx  test_drone_field  test_drone_state_jac
            test_drone_params  test_layouts  test_clearance_terms  test_barrier_terms
  route1_direct/                  STAGE 1
    solve_drones_direct.m  verify_drones_direct.m  drone_table1.m  plot_drones.m  movie_drones.m
    data/  results/  tests/
  route2_indirect/                STAGE 2
    drone_pmp_rhs.m  drone_pmp_prop.m  run_drones_pmp.m  drone_minenergy_study.m
    results/  tests/
  quadrotor/                      STAGE 3
    solve_quads_direct.m  run_quad_pmp.m  quad_pmp_rhs.m  quad_pmp_prop.m
    quad_minenergy_study.m  movie_quads.m
    data/  results/  tests/
```

Path policy: every entry point adds its own folder, `drone_dynamics/` and
`oclib/` relative to `mfilename('fullpath')`; CasADi is resolved the way
`cr3bp_minenergy_pmp` does it (from `~/casadi-3.7.0` if `casadi.Opti` is not
on the path). Nothing references `orbit_transfer`.

Fixtures: `.mat` files are git-ignored repo-wide; the ones the plans name
are force-added exactly as `ex3_cart_pole_pmp/data/*.mat` were, each with a
generator so it is reproducible.

## Testing strategy

- Physics first, mirror last: the mechanics oracle (Stage 0), the barrier
  derivatives against complex step (Stage 0), the closed-form head-on
  oracle (Stage 1) all run before any solve is trusted.
- Each solve ships its gates inside the front door (`.ok/.why`); the tests
  assert the front door's verdict AND re-derive the headline numbers.
- Every test that could pass by construction says so in its header (the
  cart-pole precedent) and is kept as a consistency check, not evidence.
- Slow things (M ≥ 32 direct, M = 64/128 stretch, the mesh gate at 2N for
  M = 16) are not in the test suite; they run from batch scripts and write
  the tracked CSVs.
- `run_all_tests.m` gains the `drone` suite after `conjugate_points`.

## Out of scope

The transformer or any neural surrogate, the paper's 3D room with box
obstacles, bounded-control (clipped) indirect solves (the saltation
machinery of the min-fuel line), rotor lag and aerodynamic drag, a
quaternion attitude variant, and a teaching write-up (a separate piece of
work on top of what lands here, as it was for the cart-pole).
