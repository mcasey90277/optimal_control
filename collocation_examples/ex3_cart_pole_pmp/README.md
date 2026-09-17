# Cart-Pole PMP-BVP

The indirect (Pontryagin / costate-shooting) companion to
`ex2_cart_pole_swing_up`'s direct collocation solve of the same swing-up
problem. Same dynamics, same boundary conditions, same fixed horizon;
solved the other way, and cross-checked against the direct fixture.

## What this is

`ex2_cart_pole_swing_up` finds the swing-up trajectory by discretizing the
whole path and letting `fmincon` find state and control at every node
(direct method). This folder instead applies Pontryagin's Minimum
Principle: writes down the boundary-value problem the optimal trajectory
must satisfy, and finds the one unknown that trajectory needs — the
initial costate `lam(0)` — by multiple shooting. It reuses the physics
(`cartpole_field.m` reproduces `ex2`'s `cart_accel`/`pendulum_accel`
exactly, pinned by test) and the shared shooting engine `oc.ms_bvp`
(`oclib`), not the orbit-transfer campaigns' machinery.

## The problem

State `x = [q1; q2; q1dot; q2dot]` (cart position, pendulum angle from
hanging, their rates), dynamics affine in the control,
`xdot = F(x) + G(x) u`, cost `J = integral u^2 dt` over a FIXED horizon
`t_f = 5 s`, boundary conditions `x(0) = [0;0;0;0]` (hanging, at rest) and
`x(t_f) = [0; pi; 0; 0]` (upright, at rest). The control is unconstrained
(no bound enters the PMP formulation — the direct fixture only stays off
its own relaxed bound to approximate the same unconstrained optimum), so
stationarity of the Hamiltonian `H = u^2 + lam'(F + G u)` gives the control
in closed form,

    dH/du = 2u + lam'G = 0   =>   u* = -lam'G / 2,

with `d2H/du2 = 2 > 0` making it the minimizer everywhere — no switching
structure, which is what makes this problem a clean test of the shooting
engine rather than of switch-detection logic. Substituting `u*` into the
costate equation `lamdot = -dH/dx` gives the autonomous 8-state PMP flow
`y = [x; lam]` (`cartpole_pmp_rhs.m`). The boundary-value problem is then:
`x(0)` is already fixed (4 known values), so the 4 unknowns are `lam(0)`,
matched against the 4 terminal conditions `x(t_f) = x_f` at the fixed
`t_f` — square, four unknowns against four equations, with no free final
time and no extra parameters (`prob.nExtra = 0`, `opts.fixedTf = true` in
the `oc.ms_bvp` call).

## The six files

| file | owns |
|---|---|
| `cartpole_field.m` | The dynamics split `xdot = F(x) + G(x) u`; same equations and constants as `ex2_cart_pole_swing_up/try2`'s `cart_accel`/`pendulum_accel` (pinned by `test_cartpole_field`); complex-step safe. |
| `cartpole_state_jac.m` | The generated exact `4x4` state Jacobian `A = d(F + G u)/dx` at fixed `u`, derived symbolically and written out by hand as the house-styled wrapper around `gen_state_jac.m`'s generated body. Generated rather than complex-stepped because it is used *inside* a function that is itself complex-stepped (see `cartpole_pmp_prop.m`). |
| `gen_state_jac.m` | One-shot generator (Symbolic Math Toolbox + `matlabFunction`) of `cartpole_state_jac`'s raw body. Run once; the generated file is committed, not regenerated on every use. |
| `cartpole_pmp_rhs.m` | The 8-state PMP field `y = [x; lam]`: the closed-form control `u* = -lam'G/2`, the state rows `F + G u*`, and the costate rows `-A'lam`. This file, not `cartpole_field.m`, is where the optimality condition lives. |
| `cartpole_pmp_prop.m` | The propagator contract `oc.ms_bvp` expects (`prob.prop`): `dt, y0, needSTM -> yEnd, PHI`, integrating with `ode113` and taking `PHI` by COMPLEX STEP through `cartpole_pmp_rhs` (safe only because `cartpole_state_jac` is generated, not itself complex-stepped — a complex step inside a complex step would corrupt the inner derivative). Throws `cartpole_pmp_prop:collapse` on integrator failure, per the engine's contract; see its header for exactly which errors are relabelled and which are rethrown unchanged. |
| `gen_direct_ref.m` | One-shot generator of `data/cartpole_direct_ref.mat`: a direct trapezoidal solve of the SAME problem at a relaxed force bound (so it approximates the unconstrained optimum the PMP-BVP solves), harvesting the defect Lagrange multipliers the PMP seed is built from. Run once; the `.mat` is committed. |

`run_cartpole_pmp.m` is the front door and is not one of the six: it loads
the direct fixture, maps its defect multipliers to a costate seed via
`oc.duals_to_costates` (trapezoid station rule, sign resolved against the
direct solve's own control), shoots `oc.ms_bvp`, flies the recovered
control through the true dynamics (`oc.fly_control`) to check it arrives,
and reports `J`, the terminal and flown misses, and the stationarity
residual — optionally plotting the indirect solution against the direct
one.

## Running the solve

```matlab
cd collocation_examples/ex3_cart_pole_pmp
run_cartpole_pmp                                          % K=8, plots vs. the direct fixture
out = run_cartpole_pmp(struct('K', 16, 'plot', false));   % no plot, returns the struct
```

Headless, from the repository root:

```
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd collocation_examples/ex3_cart_pole_pmp; run_cartpole_pmp"
```

`run_cartpole_pmp.m` adds only `oclib` to the path beyond its own folder —
it needs nothing from `orbit_transfer` (see the header note on `.engine`
if you want to pass the `costate_common` delegate instead of `@oc.ms_bvp`;
that requires putting that folder on the path yourself).

## Running the tests

Each test is a plain function returning `ok`; run any of them directly
from `tests/` (each adds the paths it needs):

```matlab
cd collocation_examples/ex3_cart_pole_pmp/tests
test_cartpole_field        % cartpole_field vs. ex2's own helpers (~instant)
test_cartpole_state_jac    % cartpole_state_jac vs. complex-step and finite-difference (~instant)
test_cartpole_pmp_rhs      % the PMP field IS the PMP conditions (~instant)
test_cartpole_pmp_prop     % the propagator + STM contract (~instant)
test_direct_ref            % the committed direct fixture is a real, unconstrained solution (~instant)
test_cartpole_pmp          % the whole indirect solve end to end (~2 min: two full ms_bvp shoots, K=8 and K=16)
```

Or headless: `matlab -batch "cd tests; ok = test_cartpole_pmp"` (etc.).

## Fixture provenance

- `data/cartpole_direct_ref.mat` — generated by `gen_direct_ref.m` (run
  once; committed). Holds `tN, X, U, muDefect, J, p, tf, uMax, N`: the
  201-node direct trapezoidal solution, its control, and its trapezoidal
  defect multipliers, at `uMax = 2000 N` (the example's own bound is
  `40 N`; this relaxed bound is what makes the fixture approximate the
  UNCONSTRAINED optimum the PMP-BVP solves, per `test_direct_ref`'s check
  that the control stays off it).
- `data/cartpole_pmp_ref.mat` — a regression baseline only, holding
  `lam0, J`. No generator script: it was saved by hand from a converged
  `run_cartpole_pmp(struct('K', 8))` call. If the physics, the seed rule,
  or the engine changes on purpose and `lam0` is expected to move, re-save
  it the same way and note the change here — `test_cartpole_pmp`'s
  regression check has no other source of truth.

## Why this exists

This folder is the second TOP-LEVEL consumer of `oc.ms_bvp` (`oclib`),
alongside every `orbit_transfer` costate campaign — nothing here is an
orbit, a CR3BP quantity, or a pumpkyn call, which is exactly the property
`oclib`'s admission rule requires a second consumer to demonstrate before
a function may move there. That demonstration is the entire reason
`ms_bvp` was promoted from `orbit_transfer/costate_common/ms_bvp.m` (now a
thin delegate) to `oclib/+oc/ms_bvp.m`; `oclib/README.md`'s roadmap entry
and `orbit_transfer/costate_common/TODO.md`'s cleanup step 13 both name
this folder as the consumer that satisfied the rule. If this demo stops
running, the promotion loses its justification — which is also why it
stays a demo and not a campaign: no `results/`, no sweeps, nothing beyond
what it takes to keep exercising the engine honestly.
