# cartpole_common

The cart-pole plant, in one place, shared by every example that solves it
(minimum energy, minimum time, minimum fuel) so the dynamics and constants
have exactly one home instead of one copy per example. A sign error in the
pendulum row survived until 2026-09-17 precisely because the equations
existed in several copies and every test compared one copy against another
rather than against an independent oracle -- see `test_cartpole_physics.m`.

## What lives here

| file | owns |
|---|---|
| `cartpole_params.m` | THE physical constants: `p = cartpole_params()` returns `struct('m1',5,'m2',1,'L',2,'g',9.8)` -- cart mass, bob mass, pendulum length, gravity. |
| `cartpole_field.m` | The dynamics split `xdot = F(x) + G(x) u`, so the control enters affinely. Complex-step safe (no `abs`/`norm`/`max`/`min`/branch on a state value). |
| `cartpole_state_jac.m` | The exact `4x4` state Jacobian `A = d(F + G u)/dx` at fixed `u`. A thin wrapper that unpacks `p` and calls the generated `cartpole_state_jac_gen.m` -- no hand-copied body, so regenerating is one command. |
| `cartpole_state_jac_gen.m` | The generated body (Symbolic Math Toolbox + `matlabFunction`); not hand-edited. |
| `gen_state_jac.m` | Generator: calls `cartpole_field` ITSELF on symbolic inputs, so it holds no second copy of the physics, and writes `cartpole_state_jac_gen.m`. Re-run whenever `cartpole_field` changes. |
| `tests/test_cartpole_physics.m` | THE independent physical oracle: power balance `dE/dt = u q1dot` derived from geometry alone, equilibrium character (hanging stable / upright unstable), energy conservation under no force. Every other test here compares the field against itself or against `ex2`'s helpers; this one doesn't, which is why it is what caught the 2026-09-17 sign error. |
| `tests/test_cartpole_field.m` | `cartpole_field` vs. `ex2_cart_pole_swing_up/try2`'s independently authored `cart_accel`/`pendulum_accel`. |
| `tests/test_cartpole_state_jac.m` | The generated Jacobian vs. a complex-step derivative and central finite differences of `cartpole_field`. |
| `tests/test_cartpole_params.m` | `cartpole_params` returns the specified constants, exactly the four fields the plant reads, and `cartpole_field` accepts what it returns. |

## The rule for adding to this folder

**Objective-independent only.** Anything that knows about a cost, a
boundary condition, a horizon, or a solve method belongs in an example
folder (`ex1_block_move`, `ex2_cart_pole_swing_up`, `ex3_cart_pole_pmp`,
`ex4_cart_pole_mintime`, and `ex5_cart_pole_minfuel` -- per
`docs/superpowers/specs/2026-09-17-cartpole-three-objectives-design.md`),
not here. This folder is the plant
and nothing else: the state, its derivative, and the derivative's Jacobian,
plus the physical oracle that pins them independently of any example's own
agreement with itself. If a function needs to know whether the problem is
minimum-energy, minimum-time, or minimum-fuel, it does not belong in
`cartpole_common`.

## Consumers

- `ex3_cart_pole_pmp/` (minimum-energy / PMP-BVP swing-up) -- the first
  consumer. Its files (`run_cartpole_pmp.m`, `run_tests.m`,
  `movie_cartpole.m`, `cartpole_pmp_rhs.m`, `cartpole_pmp_prop.m`) all add
  this folder to the path for `cartpole_field`/`cartpole_state_jac`
  (called directly by `cartpole_pmp_rhs.m`; the rest need it because they
  call things that call it). `cartpole_params()` itself has a single
  direct caller: the study script, `cartpole_minenergy_study.m` --
  `run_cartpole_pmp.m` instead reads `p` back out of the committed direct
  fixture.
- `ex4_cart_pole_mintime/` (minimum-time) and `ex5_cart_pole_minfuel/`
  (minimum-fuel) are expected to copy this folder's consumption pattern
  next, per
  `docs/superpowers/specs/2026-09-17-cartpole-three-objectives-design.md`
  (the reason this extraction happened before either was written, not
  after).

Nothing in this folder references `orbit_transfer`, a CR3BP quantity, or
pumpkyn in an executable line, and nothing here should start to.
