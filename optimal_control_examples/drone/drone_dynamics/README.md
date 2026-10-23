# drone_dynamics

The 3D drone plant family, in one place, shared by every drone example
(the paper's kinematic multi-agent problem, the point-mass planner model,
and the 6-DOF quadrotor), so each model's constants, field and Jacobian
have exactly one home. The `cartpole_common` rule applies: nothing here
knows a cost, a horizon, a boundary condition or a solver.

## Models (`drone_params(model)`, `drone_field(x, p)`, `drone_state_jac(x, u, p)`)

| model | state | control | dynamics |
|---|---|---|---|
| `kinematic` | `p` (3) | `v` (3) | `pdot = v` -- Zhang, Liu, Darbon & Karniadakis (arXiv:2610.30388) |
| `pointmass` | `[p; v]` (6) | `a` (3) | `vdot = a - g e3` |
| `rigid6dof` | `[p; v; Theta; omega]` (12) | `f` (4) rotor thrusts | rigid body, ZYX Euler, "+" mixer, Crazyflie-class constants (Luis & Le Ny 2016; Forster 2015) |

World frame: local-level, z up (ENU-compatible: x east, y north, z up); gravity is −g e3. NED callers rotate at their boundary.

All three are control-affine, `xdot = F(x) + B(x) u`, and complex-step
safe. `drone_state_jac` for `rigid6dof` is GENERATED (`gen_drone_jacs.m`,
Symbolic Math Toolbox) by differentiating `drone_field` itself; the
generated file is committed and `test_drone_state_jac` regenerates it and
compares byte for byte. The ZYX kinematics are singular at pitch 90
degrees; `p.thetaMax = 60 degrees` is the box the solvers impose.

## Geometry

| file | owns |
|---|---|
| `drone_layouts.m` | the paper's Table 1 layouts by M (inner ring radii inferred from its J_LB column), our `sphere` layout with a ball obstacle, the two-agent `headon` oracle layout |
| `pair_list.m` | all unordered pairs |
| `clearance_terms.m` | pair gaps `h = rho - 2r`, distances, unit vectors, ball gaps |
| `barrier_shape.m` / `barrier_terms.m` | the paper's C1 log barrier (eq. 3b) and its summed gradient and sparse Hessian |
| `lower_bound.m` | `J_LB = sum |xT - x0|^2 / (2T)` |
| `hermite_resample.m` | the trapezoid interpolant on a fine grid, for looking between nodes |

## The oracle

`tests/test_drone_physics.m` derives its expectations from mechanics and
geometry, never from the implemented rows: hover and free-fall limits, a
pointwise power balance `dE/dt = T (R e3).v + tau.omega` (the gyroscopic
term does no work), kinematic consistency `Rdot = R [omega]x`, and
conservation of the inertial angular momentum in torque-free flight. Every
other test compares the field against a differentiation of itself and says
so in its header.

## Consumers

`route1_direct/`, `route2_indirect/`, `quadrotor/` (Stages 1-3 of
`docs/superpowers/specs/2026-10-23-drone-multiagent-design.md`).
