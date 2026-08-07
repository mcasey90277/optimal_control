# Lessons Learned — Missile Trajectory Library

Running log. Newest entries at the top. Record what broke, what fixed it, and
what a future reader would otherwise rediscover the hard way.

## 2026-08-07 — Two files can each be correct and their composition still be wrong

State it in the general form:

> **When a quantity is represented twice — once in a state vector and once in a
> parameter struct — nothing in the arithmetic makes the two agree.** Each file
> can be correct under its own contract, and the composition still silently
> wrong. Correctness of parts is not correctness of the whole; the interface is
> where the defect lives, and no unit test can be positioned to see it.

The shape of it. Routine A integrates a quantity as a state. Routine B was
written before A existed, takes the same quantity from a parameter struct, and
never looks at the state. Compose them in one chain and the state component is
**inert** in B's phases: it rides along, it is plotted, it is reported, and it
drives nothing. Change it — at a staging event, a mass drop, a configuration
switch — and B keeps using the old parameter. The output does not diverge, does
not NaN, does not warn. It is a plausible trajectory of the wrong vehicle.

Three properties make this class of bug much worse than an ordinary one:

- **The evidence is present and misleading.** The state history shows the jump.
  A reader checking "did the mass drop at separation?" sees that it did, and
  concludes the dynamics used it.
- **Deleting the buggy code changes nothing.** Measured here: removing the
  staging link entirely left the flown trajectory **bit-identical** and changed
  only the recorded mass. A change with no effect is not a change a test can
  catch by comparing before and after.
- **Neither file is wrong.** There is no line to fix in A or in B. Reviewing
  either in isolation finds nothing, because in isolation there is nothing.

**The structural answer is a per-phase parameter set** — one carried on each
phase, so the question "which vehicle is this phase flying?" has exactly one
answer. Until that exists, the interim answer is a **guard at the adapter**: the
wrapper that lifts B into A's state space refuses to run unless the two
representations agree. Here `coorbital.eom.massConstant` raises
`coorbital:massConstant:massMismatch` when `x(7)` and `veh.mass` differ by more
than 1e-9 relative — loose enough to clear the residual an event solve leaves on
the state (measured 3.4e-10 kg on 900 kg), tight enough that a real staging jump
of hundreds of kilograms cannot pass.

**Two rules about that guard, both learned by getting them wrong first.**

- **Do not let the guard repair the disagreement.** Injecting the state value
  into a copy of the parameter struct is the obvious fix and it is worse than the
  bug, because it looks repaired. Mass was never the only thing that changed at
  separation: `Sref`, `CL` and `L/D` all still describe the jettisoned stack. A
  guard that silences the one symptom it can see, while leaving the rest wrong,
  has converted a detectable fault into an undetectable one. **Raise; do not
  reconcile.**
- **A tolerance that is only relative is not enough.** `abs(a-b) < tol*abs(b)`
  degenerates as `b` approaches zero. Floor it: never tighter than the same
  figure absolute.

And the general diagnostic, which is the transferable part: **for every quantity
in a composed system, ask how many places hold it and which one the arithmetic
actually reads.** If the answer is "two, and one of them," that is the bug,
whether or not it has bitten yet.

## 2026-08-07 — A closed-form reference cancels its own constants: the Hscale lesson, in a second place

The 2026-08-06 entry below records that an analytic reference computed from the
same constants as the model cannot detect a wrong constant. That was learned on
one constant, `Hscale`. It recurred immediately on a whole new parameter set, and
the recurrence is the point worth recording:

> **When you add a subsystem, its closed-form check is blind to that subsystem's
> constants for exactly the same reason — and you get no warning, because the
> check passes.** The lesson does not transfer itself. Re-run the diagnostic on
> every new parameter as it arrives.

The instance. The rocket equation, `deltaV = Isp*g0*ln(m0/mf)`, is the natural
validation for a powered phase and looks like a strong one. It is blind to
**every booster constant it contains**:

| Constant | Why it cancels |
|---|---|
| `Isp`, `g0` | Appear on both sides — inside `mdot = thrustVac/(Isp*g0)` on the simulation side, and inside `Isp*g0*ln(...)` on the reference side |
| `thrustVac` | In vacuum at zero incidence the whole speed equation is `dV/dt = T/m`, so the propagated result is `(thrustVac/mdot)*ln(m0/mf)`; substituting `mdot` removes it entirely. **Doubling the thrust halves the burn time and leaves the delta-V untouched.** |
| The mass ratio | The propagation starts at `m0` and terminates on an event set to `mf`, and the reference is built from the same mass fields. A wrong propellant load moves both sides identically. |

So the check validates the **shape** of the rocket equation — that the integral
of `T/m` really is an exhaust speed times a log mass ratio — and nothing
whatever about the numbers fed into it. It measures 2.4e-10 relative, which reads
as a very strong result and is a very weak one.

The defence is the same as for `Hscale` and there is no other: **pin each
constant exactly, at its source.** Six exact-equality pins carry this check, and
all six live together in `tests/test_constThrust.m:50-56` — `massDry`,
`massProp`, `thrustVac`, `Isp`, the payload `mass`, and `g0`. A comment in
`test_boost3DOF.m` names them and their line numbers, so the dependency is
written down rather than left to be rediscovered.

**One detail worth generalising: a range check is not a pin.**
`test_missileConst` bounds `g0` with `abs(c.g0 - 9.80665) < 1e-3`, which would
pass a `g0` wrong in the fourth decimal. Bracketing a constant to a plausible
band proves it is plausible; it does not prove it is the value the hand
arithmetic assumed. If a literal elsewhere in the suite was computed from a
constant, that constant needs equality, not an interval.

## 2026-08-07 — A reduction test validates only the terms that survive the reduction

State it in the general form:

> **A reduction test validates only the terms that SURVIVE the reduction.** Every
> term the reduction switches off is, by construction, invisible to it. And in
> its sharper form: **a BIT-EXACT reduction proves identity while providing zero
> evidence of correctness.** If the residual is exactly zero because the two
> routines evaluate character-identical expressions in the same order, then the
> check has confirmed that one file is a faithful copy of the other — and
> nothing whatever about whether either is right. Any error in the reference
> implementation is inherited and cancels invisibly.

This is the same disease as the shared-constant blindness recorded in the
2026-08-06 entry below, seen from the other side. There the shared object is a
*constant* that appears on both sides of a comparison and cancels; here it is a
whole *expression*. In both cases the tolerance is irrelevant: no tightening
recovers a sensitivity that is identically zero. Ask of any test that comes in
suspiciously clean: what would have to be wrong for this number to move?

Three entries in this log are now the same problem in different clothes — the
constant that cancels (2026-08-06), the expression that cancels (here), and the
whole *subsystem* whose constants cancel (the Tsiolkovsky entry above). Read them
together; the diagnostic is one question, not three.

Measured on `boost3DOF` reducing to `glide3DOF` with `env.prop` returning a dead
engine. Residual 0.000e+00 on all four states, `isequal` true, rotating and
non-rotating branches alike. That is a genuinely strong *structural* result —
transcription drift into the six shared equations is impossible — but the
mutation table shows exactly how much of the new physics it leaves uncovered:

| Mutation in `boost3DOF` | Reduction to glide | `test_boostEvents` | Force increment | Tsiolkovsky |
|---|---|---|---|---|
| thrust denominator frozen at liftoff mass | passes | passes | **fails, 35.2 % / 78.4 %** | **fails, 64.4 %** |
| `T*cos(alpha)` → `T` in `dV/dt` | passes | passes | **fails, 9.83e-3 rel** | passes (`alpha = 0` there) |
| `cos(sigma)` ↔ `sin(sigma)` on the thrust normal terms | passes | passes | **fails, 30.0 %** | passes |
| `dm/dt` sign flip | passes | **fails** | **fails, 2.00 rel** | — |

Three of the four mutations live entirely inside terms the reduction annihilates,
so the reduction cannot see any of them, and neither can an event test that only
watches the mass channel. The check that covers them is the **force increment**:
evaluate the same state twice with only the propulsion model changed, and require
the difference to equal the closed-form thrust projection. Because every shared
model is bit-identical across the two calls, the difference isolates the new
terms exactly, and it measures at 3.6e-16 relative in the clean case — the
sharpest instrument in the suite.

Two design rules follow, both learned by measurement:

- **Excite every coefficient you mean to test.** A zero angle of attack cannot
  distinguish `T cos(alpha)` from `T`; a zero bank angle cannot distinguish the
  flight-path channel from the heading channel. Choose states where no
  trigonometric factor is 1, 0, or equal to its partner, and *assert* that
  choice in the test so a later edit cannot quietly make the check vacuous.
- **Make the state disagree with the parameters.** Set the state mass to a value
  equal to no field, and no sum of fields, in either parameter struct. That is
  what turns "used the wrong mass" from a silent bug into a failing assertion.

A footnote on mutation design. The literal substitution of `veh.mass` for `x(7)`
is not runnable here: the boost phase is flown with `boosterDefaults`, which has
no `mass` field at all, so the mutant dies with `Unrecognized field name "mass"`
rather than producing a wrong number. That is a *loud* failure, and loud failures
are not what a mutation study is for. Freezing the denominator at a numeric
constant is the faithful mutant — it is what the bug would actually look like —
and it is the one that reveals the coverage gap. When a proposed mutation cannot
compile or cannot run, it has not proved the test bites; find the silent form.

## 2026-08-06 — An analytic test can never validate a constant it shares with the model

State it in the general form, because it is the most transferable thing this
milestone produced:

> **An analytic reference computed from the same constants as the model cannot
> detect a wrong constant.** The constant appears on both sides and cancels. No
> tolerance, however tight, recovers the sensitivity. The ONLY defence is to pin
> the constant at its source.

`test_missileConst` is therefore not a formality — it is the load-bearing test
for every physical constant in the library, and the analytic suite is built on
the assumption that it holds. Demonstrated three ways, all measured:

| Mutation | Analytic suite | Caught by |
|---|---|---|
| `Hscale` × 1.15 | both pass (Allen-Eggers agreement moves only +10.50 → +10.44 %) | nothing, originally — `test_missileConst` merely bracketed it to 6–9 km |
| `rho0` × 1.3 | both pass | `test_missileConst` pin; also `test_expAtmos` (sea-level pressure 87909.92 Pa) |
| `muE` × 1.02 | both pass | `test_missileConst` pin; also `test_sphereGrav` (surface gravity) |

Note *what* catches `rho0` and `muE` besides the constants pin: unit tests that
assert a hand-computed ABSOLUTE value which happens to depend on them.
`Hscale` had no such anchor anywhere in the suite — `test_expAtmos` checks that
density falls by one e-fold per scale height, which is true for *any* `Hscale`
and so is self-referential in exactly the way described above. That is the
precise reason `Hscale` was the one constant left exposed. It is now pinned to
7200 m within 1 m, with a comment saying why.

Generalising the diagnostic: for every constant, ask whether ANY test asserts a
number that would change if the constant were wrong. "Relative" checks — one
e-fold per scale height, energy conservation, a closed-form solution built from
the same constants — do not count.

**A wrong turn worth recording.** The first fix attempted was an anchor inside
`test_allenEggers` asserting the scale height against its hydrostatic value,
`H = R T / g`. That was wrong on two counts. First, Allen-Eggers needs only an
exponential *density* profile — temperature and pressure never enter the
solution, so hydrostatic consistency is not a precondition of the reference at
all, and the anchor would have false-alarmed on a perfectly legitimate change to
`T0`. Second, it did not assert much: `Rair*T0/g0` = 7317.81 m against
`Hscale` = 7200 m, a pre-existing 1.64 percent inconsistency, so a 3 percent
budget left only 1.8× headroom over a gap that was already there. The lesson
generalises: when a shared parameter turns out to be invisible, pin it where it
is *defined*; do not invent a physical-sounding cross-check in the test that
noticed the blindness.

**The same trap in its non-constant form**, caught earlier in the same task: the
equilibrium-glide reference must read `CL` from `veh.CL`, not from
`coorbital.aero.constLD`. Scaling `CL` by 1.2 *inside* `constLD` with the
reference reading it back out moves the mismatch from 14.2 percent DOWN to
9.9 percent — the mutant looks better than the original. Reading `CL` from the
vehicle instead, the same mutation drives the error from 1.51 to 11.52 percent.
Ask of every symbol in an analytic reference: where does the simulation get
this, and would a wrong value cancel?

## 2026-08-06 — Put the initial condition on the manifold, don't wait out the phugoid

The first draft of the equilibrium-glide test entered at 55 km, 5500 m/s,
gamma = -0.5 deg and then skipped the first 20 percent of samples as "the entry
transient." That entry state is 13 percent off the equilibrium speed at that
altitude (6320 m/s), so it excites a large phugoid, and at these altitudes the
oscillation damps over many hundreds of seconds. Measured decay of the max
relative error by decile: 0.165, 0.136, 0.142, 0.120, 0.102, 0.084, 0.075,
0.053, 0.049, 0.035. Skipping 20 percent leaves 14 percent error against a
5 percent tolerance; even skipping 80 percent leaves 4.9 percent. No honest
choice of skip fraction rescues it.

Placing the initial state ON the equilibrium manifold instead removes the
transient rather than waiting for it. Speed comes from the closed form, and the
flight path angle from the tangency condition `d/dt (V^2 - Veq^2(r)) = 0`, i.e.

    sin(gamma) = -(2 D/m) / (2 g + dVeq^2/dr)

with `dVeq^2/dr` a central difference on the closed form itself. That drops the
worst-case error over the WHOLE arc — no skip at all — to 1.5 percent, and it is
a stronger test: the arc now spans 45 to 20 km, a factor of 32 in density and
4.8 in speed, and the trajectory has to track a moving analytic curve the entire
way instead of merely ending up near it.

## 2026-08-06 — Allen-Eggers is a 10 percent check, not a 2 percent one

The design spec called for 2 percent agreement with the Allen-Eggers peak
deceleration. That is tighter than the approximation supports, and the reason is
quantifiable rather than hand-waved. At a 30 deg entry the closed form drops
gravity from the speed equation and freezes the flight path angle. Measured at
the peak: V = 3736 m/s against `Ve/sqrt(e) = 3639` (+2.7 percent, so +5.4 on
`V^2`), and gamma = -31.39 deg against the assumed -30 (+4.2 percent on
`|sin gamma|`). Product: +9.8 percent predicted, +10.50 percent measured. The
budget is set at 12 percent, and made one-sided BELOW as well — both neglected
effects can only add speed and steepen the path, so a simulated peak *under* the
analytic value means too much drag, not a better approximation. Note the brief's
10 percent would have failed by 0.5 points; the fix was to understand the excess,
not to nudge the number.

What keeps 12 percent from being toothless is beta-independence. The same entry
flown at two ballistic coefficients eight times apart (8571 and 1071 kg/m^2)
gives 508.1 and 502.2 m/s^2 — 1.2 percent apart — while the altitude of the peak
moves 15.0 km, from 4.86 to 19.86 km. A wrong drag area or vehicle mass cannot
move the prediction to match a wrong simulation.

## 2026-08-06 — Know what each analytic test does NOT cover

`Veq^2(r)` contains no drag term, so the equilibrium-glide test is blind to CD.
Confirmed by mutation: scaling the drag acceleration by 0.8 in `glide3DOF`
leaves `test_equilibriumGlide` passing (a wrong CD moves the vehicle along the
same curve at a different rate) while `test_allenEggers` fails by +38.3 percent.
Conversely a zero-lift Allen-Eggers entry says nothing about CL. The two tests
are complementary by construction, and neither should be described as validating
"the propagator" on its own.

## 2026-08-06 — Vacuum energy conservation is the test that matters

The glide EOM has six equations and many chances for a sign error. Propagating
a vacuum ballistic arc and checking that V^2/2 - mu/r holds constant to 1e-8
catches nearly all of them in one assertion, because any sign flip in gravity
or the centrifugal term breaks conservation immediately. Approximate checks
against analytic glide solutions do not — they tolerate a percent of error by
construction, which is enough to hide a wrong term.

## 2026-08-06 — `-batch` needs `run(...)`, not a bare relative path

Invoking the test harness as `matlab -batch "cd(...); tests/run_tests"` does
NOT work: MATLAB parses `tests/run_tests` as the expression `tests / run_tests`
(division of two undefined names) and errors with `Unrecognized function or
variable 'tests'` before `run_tests.m` ever executes. The working forms are
`matlab -batch "cd(...); run('tests/run_tests')"` or
`matlab -batch "cd(...); addpath('tests'); run_tests"`. Use `run('tests/run_tests')`
for all future headless invocations of this harness.

## 2026-08-06 — Constants live in one place

`missileConst` is the single source of truth for Earth and air constants.
pumpkyn's `getConst` was not used for these: it carries `g` and `R` but no
`muE`, no Earth radius, and two fields known to be wrong (`deg2ArcSec`, `c`).
Mixing the two would make it ambiguous which constant a routine actually used.
