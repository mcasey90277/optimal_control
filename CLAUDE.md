# CLAUDE.md

This file provides guidance to Claude Code when working with code in this repository.

## Project Overview

MATLAB implementations of **direct collocation methods for optimal control**. Educational examples progressing from simple (block move) to complex (cart-pole swing-up). Based on Matthew Kelly's trajectory optimization tutorials.

## Directory Structure

```
optimal_control/
├── ex1_block_move/           # Minimum-energy point-to-point motion
│   ├── block_main.m          # Basic implementation (10 segments)
│   ├── block_move_2.m        # Cleaner version with plots (100 segments)
│   └── block_move_3.m        # Vectorized version (100 segments)
└── ex2_cart_pole_swing_up/   # Nonlinear underactuated system
    └── trap_collocation/
        └── cart_pole_swing_up_trap_1.m  # Full cart-pole swing-up
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

## Future Extensions

Potential additions to this repository:
- Hermite-Simpson collocation (higher order)
- Multiple shooting method
- Pseudospectral methods (Legendre/Chebyshev)
- More examples: rocket landing, quadrotor, car parking
- Automatic differentiation for Jacobians

## Related Projects

- **Lunar navigation (proj7):** `/Users/msc/Desktop/proj7` — uses similar optimization for trajectory design
- **Academic goals:** Mastering calculus of variations and optimal control (Goals_2026.tex)
