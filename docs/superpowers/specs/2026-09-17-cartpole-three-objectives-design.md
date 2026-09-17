# Cart-pole, three objectives — design (2026-09-17)

Minimum time and minimum fuel added beside the existing minimum-energy
swing-up, and all three given entry scripts in the `transfer_study.m` style:
numbered sections, the machinery visible rather than behind a front door, each
necessary and sufficient condition computed IN the script with its own
PASS/FAIL, and the inline numbers asserted against the shared instruments.

The teaching payload is the contrast. One plant, one set of boundary
conditions, three costs — and the costs are what decide whether the optimal
control is a formula, a switch, or a switch with a dead zone.

## The three problems

Shared: the cart-pole of `optimal_control_examples/`, `m1 = 5 kg`, `m2 = 1 kg`,
`L = 2 m`, `g = 9.8 m/s²`, `x(0) = [0;0;0;0]`, `x(t_f) = [0;π;0;0]`, dynamics
affine in the control, `xdot = f(x) + g(x) u`.

| | min-energy (exists) | min-time (new) | min-fuel (new) |
|---|---|---|---|
| cost `J` | `∫u² dt` | `t_f` | `∫\|u\| dt` |
| `t_f` | fixed, 5 s | **free** | fixed at `c·t_min` |
| bound | none (unconstrained) | `\|u\| ≤ 40 N` | `\|u\| ≤ 40 N` |
| Hamiltonian | `u² + λ'(f+gu)` | `1 + λ'(f+gu)` | `\|u\| + λ'(f+gu)` |
| optimal control | `u* = −½λ'g` | `u* = −u_max sign σ` | `u* = −u_max sign σ` if `\|σ\| > 1`, else `0` |
| along the arc | `H` constant | `H ≡ 0` | `H` constant |
| structure | none | bang-bang | bang-off-bang |

with the **switching function** `σ(t) = λ(t)' g(x(t))` throughout. The bound is
40 N, the ex2 tutorial's own; the min-energy solution peaks at 68.5 N, so the
bound is active in both new problems and the structure is genuinely exercised.

**Min-time.** `H = 1 + λ'(f+gu)` is linear in `u`, so the minimum over
`|u| ≤ u_max` is at a vertex: `u* = −u_max sign(σ)`. The problem is autonomous
with free final time, so `H ≡ 0` along the whole arc — that is one equation, and
it closes the system: unknowns are `λ(0) ∈ R⁴` and `t_f`, five, against four
terminal state conditions plus `H(t_f) = 0`.

**Min-fuel.** With `J = ∫|u| dt` the pointwise minimisation of
`|u| + σu` over `|u| ≤ u_max` gives three cases: `u = −u_max sign σ` when
`|σ| > 1`, `u = 0` when `|σ| < 1`, and an undetermined (singular) value when
`|σ| = 1` on an interval. `t_f` is fixed — with `t_f` free the problem is
ill-posed (more time always costs no more fuel), which is why it is set to
`c·t_min` with `c = 1.5` by default, from the min-time answer.

## Layout

The plant gains a third consumer, so it stops living inside one example:

```
optimal_control_examples/
  cartpole_common/           cartpole_field, cartpole_state_jac (+ its generator),
                             test_cartpole_physics (the power-balance oracle),
                             bangbang_prop (events + saltation; min-time and
                             min-fuel both need it)
  ex3_cart_pole_pmp/         min-energy machinery + cartpole_minenergy_study.m
  ex4_cart_pole_mintime/     min-time machinery  + cartpole_mintime_study.m
  ex5_cart_pole_minfuel/     min-fuel machinery  + cartpole_minfuel_study.m
  teaching_docs/             direct_and_indirect_cartpole.tex (updated at the end)
```

Moving the plant is the same rule the sign error taught: one home for the
physics, before the copies multiply. ex3's existing tests keep passing through
the moved files; that is the move's gate.

## Bang-bang: events and saltation

Decided (recommendation 1 of three, approved 2026-09-17): the switch is handled
**exactly**, by event detection, with the state-transition matrix corrected by a
saltation factor at each crossing. Not a `tanh` smoothing homotopy, which never
reaches the corner, and not switch-time extraction from the direct solve alone.

**The propagation.** Integrate with `odeset('Events', ...)` watching `σ(t) = 0`
(min-time) or `|σ(t)| = 1` (min-fuel). At each event, stop, flip the control
branch, restart. The event times are outputs, not unknowns.

**The saltation matrix.** Across a switch at `t_s` the flow is continuous but
its derivative is not, so the STM picks up a jump. With `F⁻`, `F⁺` the fields
just before and after and `n = ∇_y σ`,

```
Φ⁺ = [ I + (F⁺ − F⁻) nᵀ / (nᵀ F⁻) ] Φ⁻ .
```

Omitting it does not break the solve loudly — it degrades the Newton step, so
the symptom is slow or failed convergence, which is why this is the item most
worth testing directly. The orbit library's min-fuel work retracted a published
verdict over exactly this omission.

**Grazing.** If `σ` touches zero without crossing, `nᵀF⁻ → 0` and the saltation
factor blows up. The propagator must detect a near-zero denominator and refuse
by name (`:grazing`) rather than return a poisoned STM.

## The study scripts

One skeleton, three instances. Sections, following `transfer_study.m`:

```
0  tolerances          one struct, feeding every printed line AND every verdict
1  the plant           constants, and the physics oracle run here (not assumed)
2  the problem         cost, bound, boundary conditions, t_f rule
3  the seed            a direct collocation solve: costates AND, for the
                       bounded problems, the switch COUNT and rough times
4  solve               multiple shooting -> costates (+ t_f for min-time)
5  independent check   a second implementation must not move the answer
6  NECESSARY (N)       Pontryagin, one condition at a time
7  SUFFICIENCY (S)     Legendre/strict-bang, conjugate point
8  assert vs library   the inline numbers against the shared instruments
9  plot                states, control, and the switching function
```

Diagnostic IDs carry the same meanings as in `transfer_study`: **N** necessary,
**S** sufficiency hypotheses, **V** validity (can a verdict be believed),
**X** cross-check (a second implementation, or the direct solve).

Per-objective specifics of sections 6–7:

| ID | min-energy | min-time | min-fuel |
|---|---|---|---|
| N1 | BVP residual | BVP residual (incl. `H(t_f)=0`) | BVP residual |
| N2 | `H` constant along the arc | `max\|H\|` along the arc | `H` constant |
| N3 | terminal miss, re-flown | terminal miss, re-flown | terminal miss, re-flown |
| N4 | — (no transversality: `t_f` and `x(t_f)` both fixed) | `H(t_f) = 0` | — |
| N5 | adjoint equations, by independent differentiation | same | same |
| N6 | minimum-principle gap of the applied control | gap, **and** `u` on a vertex except at switches | gap, **and** `u ∈ {−u_max, 0, +u_max}` |
| N7 | — | switching function crosses zero exactly at the event times | `\|σ\|` crosses 1 at the event times; `\|σ\| ≤ 1` on coast arcs |
| S1 | `H_uu = 2 > 0` (strict Legendre) | strict bang: `min\|σ\|` away from switches > 0 | strict bang-off: `\|\|σ\|−1\|` away from switches > 0 |
| S2 | conjugate point (fixed-t_f Jacobi) | conjugate point (free-t_f form) | conjugate point (fixed-t_f) |
| V1 | the physics oracle passes | same | same |
| V2 | inline numbers vs the library instruments | same | same |
| X1 | direct vs indirect: cost and control | direct vs indirect: `t_f` and switch times | direct vs indirect: fuel and switch times |

**`strict`**, as in `transfer_study`: a failed N or X gate throws (the root or
the implementation is broken); a failed or unresolved S gate is a finding about
the trajectory, reported and not thrown, unless `strict` is set.

## The conjugate test, and a promotion it earns

Section 7 needs a Jacobi/conjugate-point test. `ms_conjugate_test` already
exists in `orbit_transfer/costate_common/`, supports both the fixed-time and
free-time conventions, and is the **next** promotion candidate in that folder's
cleanup step 13 — blocked only for want of a second TOP-LEVEL consumer.

These scripts are that consumer. **Decision: promote it**, exactly as `ms_bvp`
was promoted — move to `oclib/+oc/ms_conjugate_test.m`, leave a delegate in
`costate_common`, gate the move on `golden_cells` unchanged, the four-campaign
ladder re-solve unchanged, and the conjugate tests green. The study scripts then
call `oc.ms_conjugate_test` and keep `optimal_control_examples/` free of any
`orbit_transfer` dependency, which is the property that made the first promotion
legitimate.

The alternative — a local Jacobi test written into each script — was rejected:
three copies of delicate machinery, and the scripts' own rule is that anything
whose reimplementation would create a second unverified copy stays in the
library and is asserted against (section 8).

## Seeding

Each indirect solve is seeded by a direct collocation solve of the same problem,
as ex3 does — `oc.duals_to_costates` with the trapezoid rule, sign resolved by
the scalar correlation against the direct control.

The bounded problems need one thing more: the **switch structure**. The direct
solve's control history gives the number of sign changes and their approximate
times, which seed the event-driven propagation and, more importantly, tell the
study script how many switches to expect — a mismatch between predicted and
realised switch count is itself a reportable finding (N7).

Min-fuel additionally inherits `t_min` from the min-time study, so the
min-time sub-project must land first.

**Named fallback**, as in the ex3 spec, so it is a decision and not a scramble:
if a direct seed will not converge for the bounded problems, walk a homotopy in
the bound (`u_max` from large, where the answer resembles min-energy, down to
40 N), each solve seeding the next.

## Sub-projects, in order

**A — shared plant + min-energy study.** Extract `cartpole_common`; promote
`ms_conjugate_test`; write `cartpole_minenergy_study.m` against the already
working solve. Smallest, and it builds the skeleton the other two instantiate.
No new mathematics.

**B — min-time.** Direct solve with free `t_f` (time as a decision variable on a
normalised mesh); `bangbang_prop` with events and saltation; the 5-unknown BVP;
`cartpole_mintime_study.m`. Produces `t_min`, which C needs.

**C — min-fuel.** Direct solve at `t_f = 1.5·t_min`; the three-branch control;
singular-arc detection; `cartpole_minfuel_study.m`.

Each gets its own implementation plan. A is specified enough to plan now; B and
C are specified here at design level and get their plans when their predecessor
lands, so that B's plan can quote A's interfaces and C's can quote `t_min`.

## Verification

Every sub-project, before it is called done:
- the physics oracle passes (it is in `cartpole_common` and every script runs it);
- direct and indirect agree — on cost, and for the bounded problems on the
  switch count and switch times;
- the terminal miss and the flown-control miss are inside gates the script
  states;
- `run_tests` in each folder is green, and the ex3 suite still passes after the
  extraction;
- for the promotion: `golden_cells` 20/20 with unchanged numbers, and the
  four-campaign ladder re-solve at its recorded 1e-14 figures.

## Out of scope

- State constraints (a track limit on `q1`), which would add junction
  conditions and a different multiplier structure.
- Singular arcs as a *solved* case in min-fuel: they are **detected** and
  reported, not constructed. If the measured trajectory has one, that is a
  finding and a follow-on piece of work.
- A min-fuel front over many `t_f` values (the orbit campaigns' ΔV–t_f curve).
  One `t_f` here.
- Rewriting ex2. It stays as the direct-collocation reference.

## Risks

| risk | consequence | response |
|---|---|---|
| grazing switch (`σ` touches 0) | saltation denominator → 0, poisoned STM | detect and refuse by name; the study reports it as a structural finding |
| abnormal extremal in min-time (`λ₀ = 0`) | the PMP conditions hold vacuously | check the normalisation explicitly, as the orbit library's homogeneous chart does |
| chattering (many closely spaced switches) | event loop stalls | cap the event count and refuse; report the count the direct solve predicted |
| singular arc in min-fuel (`\|σ\| ≡ 1`) | control undetermined on an interval | detect via the dwell of `\|σ\|−1` near zero; report, do not construct |
| the direct seed misses the basin | no solve | the `u_max` homotopy above |
