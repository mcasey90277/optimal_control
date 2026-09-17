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

## The files

| file | owns |
|---|---|
| `cartpole_field.m` | The dynamics split `xdot = F(x) + G(x) u`; same equations and constants as `ex2_cart_pole_swing_up/try2`'s `cart_accel`/`pendulum_accel` (pinned by `test_cartpole_field`); complex-step safe. |
| `cartpole_state_jac.m` | The exact `4x4` state Jacobian `A = d(F + G u)/dx` at fixed `u`. A thin wrapper that unpacks `p` and calls the generated `cartpole_state_jac_gen.m`; there is no hand-copied body, so regenerating is one command. Generated rather than complex-stepped because it is used *inside* a function that is itself complex-stepped (see `cartpole_pmp_prop.m`). |
| `gen_state_jac.m` | Generator (Symbolic Math Toolbox + `matlabFunction`): calls `cartpole_field` ITSELF on symbolic inputs, so it holds no second copy of the physics, and writes `cartpole_state_jac_gen.m`. Re-run it whenever `cartpole_field` changes. |
| `cartpole_pmp_rhs.m` | The 8-state PMP field `y = [x; lam]`: the closed-form control `u* = -lam'G/2`, the state rows `F + G u*`, and the costate rows `-A'lam`. This file, not `cartpole_field.m`, is where the optimality condition lives. |
| `cartpole_pmp_prop.m` | The propagator contract `oc.ms_bvp` expects (`prob.prop`): `dt, y0, needSTM -> yEnd, PHI`, integrating with `ode113` and taking `PHI` by COMPLEX STEP through `cartpole_pmp_rhs` (safe only because `cartpole_state_jac` is generated, not itself complex-stepped — a complex step inside a complex step would corrupt the inner derivative). Throws `cartpole_pmp_prop:collapse` on integrator failure, per the engine's contract; see its header for exactly which errors are relabelled and which are rethrown unchanged. |
| `gen_direct_ref.m` | One-shot generator of `data/cartpole_direct_ref.mat`: a direct trapezoidal solve of the SAME problem at a relaxed force bound (so it approximates the unconstrained optimum the PMP-BVP solves), harvesting the defect Lagrange multipliers the PMP seed is built from. Run once; the `.mat` is committed. |

`run_cartpole_pmp.m` is the front door and is not one of these: it loads
the direct fixture, maps its defect multipliers to a costate seed via
`oc.duals_to_costates` (trapezoid station rule, sign resolved against the
direct solve's own control), shoots `oc.ms_bvp`, flies the recovered
control through the true dynamics (`oc.fly_control`) to check it arrives,
and returns its OWN verdict (`.ok`, `.why`) alongside `J`, the terminal
miss, the engine's own last-arc residual, the flown miss and the spread of
the Hamiltonian — optionally plotting the indirect solution against the
direct one. It does not report the residual `2u + lam'G`: `U` is built as
`-lam'G/2`, so that number is zero by construction and measures nothing.

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

Run them ALL through the runner, which fails the process if any check
fails — a bare `test_x` prints `FAIL` and still exits 0, which is fine
interactively and useless in automation:

```matlab
cd collocation_examples/ex3_cart_pole_pmp
run_tests                                      % ~3 min, errors on any failure
```

```
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('<this folder>'); exit(~run_tests())"
```

Individually, from `tests/` (each adds the paths it needs):

```matlab
test_cartpole_physics      % THE independent oracle: power balance, equilibria, energy (~instant)
test_cartpole_field        % cartpole_field vs. ex2's own helpers (~instant)
test_cartpole_state_jac    % the generated Jacobian vs. complex-step and finite-difference (~instant)
test_cartpole_pmp_rhs      % the PMP field IS the PMP conditions (~instant)
test_direct_ref            % the committed direct fixture: feasible, stationary, off its bound (~instant)
test_cartpole_pmp_prop     % the propagator, the STM, and the collapse contract by identifier (~1 min)
test_cartpole_pmp          % the whole indirect solve end to end (~2 min: two ms_bvp shoots, K=8 and K=16)
```

## Fixture provenance

- `data/cartpole_direct_ref.mat` — generated by `gen_direct_ref.m` (run
  once; committed). Holds `tN, X, U, muDefect, J, p, tf, uMax, N`: the
  201-node direct trapezoidal solution, its control, its trapezoidal
  defect multipliers and `solver` (fmincon's exit flag, first-order
  optimality and constraint violation, all asserted before the fixture is
  written), at `uMax = 2000 N` (the example's own bound is `40 N`; this
  relaxed bound is what makes the fixture approximate the UNCONSTRAINED
  optimum the PMP-BVP solves, per `test_direct_ref`'s check that the
  control stays off it — measured `max|u| = 68.5 N`).
- `data/cartpole_pmp_ref.mat` — a regression baseline only, holding
  `lam0, J`. No generator script: it was saved by hand from a converged
  `run_cartpole_pmp(struct('K', 8))` call. If the physics, the seed rule,
  or the engine changes on purpose and `lam0` is expected to move, re-save
  it the same way and note the change here — `test_cartpole_pmp`'s
  regression check has no other source of truth, so it pins a previous run
  of this same code and would preserve a mistake as faithfully as a fix.
  A wobble of order 1e-11 in `lam0` between runs is expected (a different
  solver stopping iterate on the same root) and is NOT a reason to
  regenerate: doing so on every wobble is how a tripwire stops being one.

## What it produces (measured 2026-09-17)

| quantity | value | what it means |
|---|---|---|
| `J` indirect / direct | 2779.381719 / 2781.109846 | 0.062% apart; the gap is the 201-node trapezoid's discretization error, and the indirect solve is the more accurate of the two |
| terminal miss | 3.67e-11 | re-flying `lam0` over the whole 5 s and comparing with `x_f`. Gate 1e-9 |
| engine residual | 7.27e-14 | the shooting's OWN last-arc terminal residual, independent of the reporting flight |
| flown miss | 1.67e-07 | flying the closed-form control `u = -lam'G/2` through the true dynamics. Gate 1e-6 |
| Hamiltonian spread | 5.77e-13 relative | `H` is constant (not zero: `t_f` is fixed) along an extremal of this autonomous flow |
| seed quality | correlation 0.9998, amplitude ratio 1.009 | the dual-derived costates agree with the direct solve's control in sign AND scale |
| `max\|u\|` | 68.5 N | well inside the fixture's relaxed 2000 N bound |

## A sign error lived here until 2026-09-17

The pendulum row of the dynamics carried the wrong sign — inherited from
`ex2_cart_pole_swing_up`, reproduced here on purpose (this folder's first
test pins `cartpole_field` to those helpers), and copied a third time into
the symbolic Jacobian's generator. With it, no sign convention of the
energy was conserved at `u = 0` and the *hanging* equilibrium was unstable.
GPT-6 Astra found it in review (`reviews/astra_review_2026-09-17.md`); it
survived everything else because the direct and indirect solves agreed to
0.075% while solving the same wrong problem, and every other check compared
the field against itself or against `ex2`.

What now prevents a repeat is `tests/test_cartpole_physics`: it derives its
reference from the GEOMETRY alone — bob at `(q1 + L sin q2, -L cos q2)` —
and demands the pointwise power balance `dE/dt = u q1dot`, the right
equilibrium characters, and energy conservation under no force. It fails
loudly on the old equations (337 W imbalance, 64 J drift over 3 s) and
passes on the corrected ones (5.7e-14 W, 1.6e-12 J). Any future edit to the
physics answers to it, not to another copy of itself.

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
