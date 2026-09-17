# Cart-pole PMP-BVP — design (2026-09-16)

An indirect solve of the cart-pole swing-up the `optimal_control_examples`
already solve directly: form the Pontryagin boundary-value problem and shoot
it with `ms_bvp`. Two purposes, in order:

1. a compact, tested demonstration that `ms_bvp` solves a problem with no
   orbit, no CR3BP and no pumpkyn in it;
2. the second TOP-LEVEL consumer that `oclib`'s admission rule requires
   before `ms_bvp` can move into `oclib/+oc` (costate_common TODO step 13).

The demo is not a tutorial. The teaching write-up (the house
`matlab-coding-tutorial` format) is a separate piece of work, later, on top
of what lands here.

## The problem

The same swing-up the direct example solves at its relaxed force bound:

| | |
|---|---|
| state | `x = [q1; q2; q1dot; q2dot]` — cart position (m), pendulum angle (rad, 0 = down), and their rates |
| control | `u` — horizontal force on the cart (N), **unconstrained** |
| cost | `J = ∫ u² dt` over a **fixed** `t_f = 5 s` |
| boundary | `x(0) = [0; 0; 0; 0]`, `x(t_f) = [0; π; 0; 0]` |
| constants | `m1 = 5 kg`, `m2 = 1 kg`, `L = 2 m`, `g = 9.8 m/s²` (the existing example's) |

The dynamics are affine in the control, `xdot = f(x) + g(x) u`, with

```
f = [q1dot;
     q2dot;
     (L m2 sin(q2) q2dot² + m2 g cos(q2) sin(q2)) / D1;
     (L m2 cos(q2) sin(q2) q2dot² + (m1+m2) g sin(q2)) / D2];
g = [0; 0; 1/D1; cos(q2)/D2];
D1 = m1 + m2 (1 - cos²(q2));
D2 = L (m1+m2) (1 - (m2/(m1+m2)) cos²(q2));
```

which is the existing `cart_accel` / `pendulum_accel` split at the `u` terms.

**PMP.** With `H = u² + λ'(f + g u)`:

- stationarity `∂H/∂u = 2u + λ'g = 0` gives the smooth minimiser
  `u* = -λ'g / 2` (the Legendre condition `∂²H/∂u² = 2 > 0` holds everywhere,
  so this is a minimum, and no switching structure exists to get wrong);
- the costate obeys `λdot = -∂H/∂x` evaluated at `u*`. The running cost has
  no state dependence, so this is `-(∂f/∂x + u* ∂g/∂x)' λ`.

Substituting `u*` closes an 8-state flow `y = [x; λ]`. Unknowns: `λ(0) ∈ R⁴`.
Conditions: `x(t_f) = [0; π; 0; 0]`, four of them. Square, smooth, fixed-time
— `ms_bvp`'s `opts.fixedTf` path.

**Why this formulation.** The bounded-control version (`|u| ≤ 40`, what the
direct example finally reports) makes `u*` a clipped function of `λ'g`, whose
switching structure needs the saltation machinery the min-fuel work uses.
That is a different demonstration. An unconstrained, fixed-time, smooth
extremal is the cleanest possible evidence that the ENGINE is generic, and
its cross-check against the direct solve is unambiguous.

## Files

New folder `optimal_control_examples/ex3_cart_pole_pmp/`:

| file | what |
|---|---|
| `cartpole_field.m` | `[F, G] = cartpole_field(x, p)` — the drift and the control column above, one home for the physics, complex-step safe (no `abs`, no `norm`, no branch on a state). |
| `cartpole_pmp_rhs.m` | `dy = cartpole_pmp_rhs(y, p)` — the 8-state PMP field with `u*` substituted; also returns `u*` on request. |
| `cartpole_pmp_prop.m` | `[yEnd, PHI] = cartpole_pmp_prop(dt, y0, needSTM, p)` — `ode45`/`ode113` at tight tolerances, integrating the 8×8 variational equations when `needSTM`, with `A(t) = ∂(rhs)/∂y` by complex step. Throws on integrator collapse (the engine's contract). |
| `run_cartpole_pmp.m` | The front door: build `prob` (`ny = 8`, `freeIdx0 = 5:8`, the three closures), seed, call `ms_bvp` with `fixedTf`, report `J`, the terminal miss and the flown miss, plot `x(t)` and `u(t)` against the direct solution. `nargin == 0` runs the whole thing. |
| `gen_direct_ref.m` | Produces the fixture below by solving the direct trapezoidal NLP at the relaxed bound (`fmincon`, the existing example's formulation) and saving nodes, control, defect multipliers and cost. Run once; committed so the fixture is reproducible rather than a mystery file. |
| `data/cartpole_direct_ref.mat` | That fixture, committed so nothing has to run `fmincon` to seed or to compare. Small (N = 200 nodes). |
| `tests/test_cartpole_pmp.m` | Below. |

Nothing moves or changes in `ex2_cart_pole_swing_up/`: `gen_direct_ref`
reuses its dynamics helpers and its formulation, and writes the fixture.

## Seeding

The campaigns' own pattern, on a different transcription:

1. read the fixture's defect multipliers;
2. map them to costates with `oc.duals_to_costates`, `scheme = 'trapezoid'`
   (the existing example is trapezoidal; the CR3BP campaigns are
   Hermite-Simpson, so this exercises a station rule the orbit work never
   uses in anger);
3. interpolate onto the shooting grid, take `λ(0)`, and shoot.

`K = 8` segments to start. Multiple shooting is not strictly needed on a 5 s
horizon, but the point is to drive the engine the way the campaigns do, and
K > 1 exercises the continuity blocks.

**Fallback if that seed does not converge:** a homotopy in the target angle,
solving `q2(t_f) = θ` for `θ` walked from a small angle up to `π`, each solve
seeding the next. Recorded here so the fallback is a decision already made,
not an improvisation. If the fallback is used, the spec's claim that duals
seed it is withdrawn in the README rather than quietly dropped.

## Verification

`tests/test_cartpole_pmp.m`, all with no campaign folder on the path:

1. **Physics agrees with the direct example.** `cartpole_field`'s
   `F + G u` reproduces `try2/cart_accel` and `try2/pendulum_accel` to
   1e-12 at a scatter of states and controls — the two demos provably solve
   the same problem.
2. **The Jacobian is right.** The complex-step `A` matches central finite
   differences to 1e-8 at several points along the extremal.
3. **The STM is right.** `PHI` matches a finite-difference sensitivity of
   the propagated map to 1e-6.
4. **The BVP converges**: terminal miss < 1e-9 in the four boundary
   components, and `ms_bvp` reports converged.
5. **The control is the PMP control**: `max |2u + λ'g|` along the solution
   < 1e-8 (stationarity actually holds, not just the endpoints matching).
6. **Flying it lands there**: `oc.fly_control` integrating `f + g u*(t)` from
   `x(0)` arrives within 1e-6 of the target — the same G1b idea the
   campaigns gate on.
7. **It agrees with the direct solution**: `J_indirect` within 1% of the
   fixture's `J`, and the control profiles within 2% RMS. The gap is the
   direct method's discretization error, not a tolerance to tune; the test
   states which side is expected to be more accurate (the indirect one).
8. **Regression**: stored `λ(0)` reproduced to 1e-8.

## The promotion (why this exists)

Once the demo passes:

1. `git mv orbit_transfer/costate_common/ms_bvp.m oclib/+oc/ms_bvp.m`;
2. leave `costate_common/ms_bvp.m` as a delegate, exactly as
   `duals_to_costates` is, so no orbit_transfer caller changes;
3. `run_cartpole_pmp` calls `oc.ms_bvp` directly (a top-level consumer that
   does not go through costate_common at all);
4. update `oclib/README.md`'s table with both consumers and the gate;
5. update `costate_common/README.md`, its TODO step 13, and
   `tests/test_folder_rules` (the delegate is exercised by the existing
   tests; no new exemption).

**Equivalence gates for the move**, all required:

- `golden_cells` 20/20 with unchanged numbers;
- the four-campaign ladder re-solve (`HALO`, `DPO`, `HALO_HALO`, `GTO`)
  matching the shipped sheets as it does today;
- `test_ms_bvp_extra`, `test_ms_bvp_fixedtf`, `test_ms_tfmin_hom`,
  `test_arclength_ms` green;
- the cart-pole test green through `oc.ms_bvp`.

## Out of scope

- The bounded-control (switching) version, and the minimum-time free-`t_f`
  version. Either would be a second demo later; neither is needed to earn
  the promotion.
- `ms_conjugate_test`'s promotion. The cart-pole extremal could drive its
  fixed-time form, but a second consumer for it needs a known-answer
  conjugate case, which is its own piece of work.
- The tutorial write-up.
- Any change to `ex2_cart_pole_swing_up/`.

## Risks

| risk | what happens | response |
|---|---|---|
| The dual-seeded `λ(0)` misses the basin | `ms_bvp` walks off or stalls | the angle homotopy above; recorded as a finding either way |
| Complex-step unsafe code creeps in (an `abs`, a `max`) | derivatives silently wrong, convergence mysterious | check 2 catches it; the field file says the rule in its header |
| The 1% cost agreement is optimistic | test 7 fails on a correct solve | the fixture records the direct solve's own defect and mesh; if the gap is discretization, refine the fixture rather than widen the tolerance, and say so |
| The demo drifts out of date once promoted | the second consumer stops being real | `oclib/README.md` names it as the gate; it is a test, so it runs |
