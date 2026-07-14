# CLAUDE.md

This file provides guidance to Claude Code when working with code in this repository.

## Project Overview

Educational implementations of **optimal control and estimation** methods. Core examples use MATLAB direct collocation (block move → cart-pole swing-up), based on Matthew Kelly's trajectory optimization tutorials. Also includes MPC and Lie group filtering (SO(3) error-state EKF) as LaTeX documents with planned MATLAB implementations.

## Directory Structure

```
optimal_control/
├── ex1_block_move/              # Minimum-energy point-to-point motion
│   ├── block_main.m             # Basic implementation (10 segments)
│   ├── block_move_2.m           # Cleaner version with plots (100 segments)
│   └── block_move_3.m           # Vectorized version (100 segments)
├── ex2_cart_pole_swing_up/      # Nonlinear underactuated system
│   ├── trap_collocation/
│   │   └── cart_pole_swing_up_trap_1.m   # Monolithic implementation
│   └── try2/                    # Modularized refactor
│       ├── cart_pole_swing_up_trap_2.m   # Main script using helpers
│       ├── nonlcon_cart_pole.m           # Equality constraints (BCs + defects)
│       ├── objective.m                   # Cost function (trap quadrature of u²)
│       ├── cart_accel.m                  # Cart acceleration from Lagrangian EOM
│       └── pendulum_accel.m             # Pendulum angular acceleration
├── lowThrust_GTO_tulip/         # CR3BP min-time GTO->tulip: theory note + guided
│                                #   tutorial + indirect (PMP shooting, complex-step)
│                                #   for min-time/energy/fuel + mytry/
├── NLP_lowThrust_GTO_tulip/     # Direct-NLP solvers: min-time/energy/fuel, cone-
│                                #   elimination, Hermite-Simpson, CasADi+IPOPT.
│                                #   Min-fuel CAMPAIGN record + next steps in
│                                #   LOW_THRUST_MINFUEL_CAMPAIGN.md (min-energy is
│                                #   the homotopy root; many-switch bang-bang open,
│                                #   needs regularized coords). movie/ has MP4+GIFs.
│                                #   Per-target min-fuel deliverables: PSR/ (GTO->
│                                #   tulip, w/ PMP refine) and elfo/ (GTO->ELFO);
│                                #   sundman_minfuel/ is the shared solver engine.
├── mpc_cart_pole/               # Model Predictive Control for cart-pole
│   ├── mpc_cart_pole_demo.m     # MPC controller (N=50 horizon, 20Hz control)
│   └── mpc_cart_pole_explained.tex  # Step-by-step code walkthrough
├── lieFiltering/                # Lie group attitude estimation on SO(3)
│   ├── SO3_Attitude_Filter.tex  # Error-state EKF on SO(3) manifold
│   └── simulating_IMU_data.tex  # Synthetic IMU data generation
├── papers/                      # Reference PDFs (Barfoot, Barrau & Bonnabel, etc.)
└── learning_docs/               # (empty, reserved)
```

## Method: Direct Trapezoidal Collocation

All examples use the same approach:

1. **Discretize time** into N segments with N+1 nodes
2. **Decision variables:** State and control at each node
3. **Dynamics constraints:** Trapezoidal defects between adjacent nodes
4. **Boundary constraints:** Initial and final conditions
5. **Objective:** Trapezoidal quadrature of running cost
6. **Solve:** Nonlinear program via `fmincon`

### Trapezoidal Defect Formula

For dynamics ẋ = f(x, u):
```
x_{k+1} - x_k - (h/2) * (f(x_{k+1}, u_{k+1}) + f(x_k, u_k)) = 0
```

### Trapezoidal Quadrature

For cost J = ∫ L(x,u) dt:
```
J ≈ Σ (h/2) * (L_k + L_{k+1})
```

## Example 1: Block Move

**Problem:** Move unit mass from x=0 to x=1 in 1 second, starting and ending at rest, minimizing ∫u²dt.

**Dynamics:**
```
ẋ = v       (velocity is derivative of position)
v̇ = u       (control force = acceleration, unit mass)
```

**Boundary Conditions:**
```
x(0) = 0,  v(0) = 0   (start at origin, at rest)
x(1) = 1,  v(1) = 0   (end at position 1, at rest)
```

**Decision Variables:** `X = [x; v; u]` or `X = [u; x; v]` depending on version

**Analytical Solution:** u(t) = 6 - 12t (linear control, crosses zero at t=0.5)

| File | Segments | Implementation Style |
|------|----------|---------------------|
| `block_main.m` | 10 | Loop-based, functions inside script |
| `block_move_2.m` | 100 | Loop-based, cleaner with plots |
| `block_move_3.m` | 100 | Vectorized (no loops), efficient |

### Running Block Move

```matlab
cd ex1_block_move
block_move_3  % recommended version
```

## Example 2: Cart-Pole Swing-Up

**Problem:** Swing inverted pendulum from hanging (θ=0) to upright (θ=π) by pushing the cart, minimizing ∫u²dt.

**Parameters:**
```matlab
m1 = 5;      % cart mass (kg)
m2 = 1;      % pendulum bob mass (kg)
L = 2;       % pendulum length (m)
g = 9.8;     % gravity (m/s²)
T = 5.0;     % time horizon (s)
N = 200;     % number of segments
```

**State Variables:**
- q1: cart position (m)
- q2: pendulum angle (rad, 0 = down, π = up)
- q1dot: cart velocity (m/s)
- q2dot: pendulum angular velocity (rad/s)

**Control:** u = horizontal force on cart (N)

**Decision Variables:** `X = [q1; q2; q1dot; q2dot; u]` (5 × num_nodes total)

**Dynamics:** Coupled nonlinear equations of motion:
```
q1ddot = (L*m2*sin(q2)*q2dot² + u + m2*g*cos(q2)*sin(q2)) / (m1 + m2*(1-cos²(q2)))

q2ddot = (L*m2*cos(q2)*sin(q2)*q2dot² + u*cos(q2) + (m1+m2)*g*sin(q2)) / (L*(m1+m2)*(1 - m2/(m1+m2)*cos²(q2)))
```

**Boundary Conditions:**
```
Initial: q1=0, q2=0, q1dot=0, q2dot=0  (cart at origin, pendulum down, at rest)
Final:   q1=0, q2=π, q1dot=0, q2dot=0  (cart at origin, pendulum up, at rest)
```

**Solution Strategy:**
1. Solve with relaxed force limits (±2000 N) to find feasible trajectory
2. Warm-start with tight limits (±40 N) for realistic solution

### Running Cart-Pole

```matlab
cd ex2_cart_pole_swing_up/trap_collocation
cart_pole_swing_up_trap_1  % takes ~1-2 minutes
```

## Solver Configuration

All examples use MATLAB's `fmincon` with SQP:

```matlab
options = optimoptions('fmincon', ...
    'Algorithm',             'sqp', ...
    'EnableFeasibilityMode', true, ...
    'SubproblemAlgorithm',   'cg', ...
    'Display',               'iter', ...
    'ConstraintTolerance',   1e-6, ...
    'MaxIterations',         1e7);
```

## Key Concepts

### Why Trapezoidal Collocation?

- Implicit method (more stable than explicit Euler)
- Second-order accurate
- Sparse constraint Jacobian structure
- Easy to implement and understand

### Decision Variable Ordering

Two common conventions:
1. **State-first:** `X = [x; v; u]` — groups by variable type
2. **Time-first:** `X = [x1,v1,u1, x2,v2,u2, ...]` — groups by time node

These examples use state-first ordering.

### Scaling Considerations

For better convergence:
- Normalize time to [0,1] when possible
- Scale state/control variables to similar magnitudes
- Use good initial guesses (critical for nonlinear problems)

## References

1. Kelly, M. "An Introduction to Trajectory Optimization: How to Do Your Own Direct Collocation." SIAM Review, 2017.
2. Betts, J. "Practical Methods for Optimal Control Using Nonlinear Programming." SIAM, 2010.
3. Tedrake, R. "Underactuated Robotics" (MIT Course 6.832)

## MPC Cart-Pole

Receding-horizon control for the cart-pole system. Predict-optimize-apply loop at 20 Hz.

- Prediction horizon: N=50 steps
- Weighted quadratic cost: Q on state error, R on control effort
- `mpc_cart_pole_explained.tex` walks through the code step by step

## Lie Filtering (lieFiltering/)

SO(3) attitude estimation using error-state EKF on the rotation manifold.

**Key mathematical objects:**
- SO(3) rotation group, so(3) Lie algebra
- Exponential map (so(3) → SO(3)) and logarithmic map (SO(3) → so(3))
- Hat/vee (wedge/vee) operators, skew-symmetric matrices
- Error-state formulation: small perturbations in the Lie algebra

**Documents:**
- `SO3_Attitude_Filter.tex` — Full EKF derivation; targeted at engineers who know Kalman filtering but not Lie theory
- `simulating_IMU_data.tex` — Generating synthetic accelerometer, gyroscope, magnetometer data for filter testing

## Future Extensions

Potential additions to this repository:
- Hermite-Simpson collocation (higher order)
- Multiple shooting method
- Pseudospectral methods (Legendre/Chebyshev)
- More examples: rocket landing, quadrotor, car parking
- Automatic differentiation for Jacobians

## Related Projects

- **Navigation** (`~/Desktop/navigation/`) — INS simulation, IMU modeling, SO(3) utilities in `lib/`
- **proj7** (`~/Desktop/proj7/`) — Tulip constellation GDOP/PDOP, cislunar SDA
- **Academic goals:** Mastering calculus of variations and optimal control (Goals_2026.tex)
