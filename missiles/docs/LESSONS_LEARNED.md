# Lessons Learned — Missile Trajectory Library

Running log. Newest entries at the top. Record what broke, what fixed it, and
what a future reader would otherwise rediscover the hard way.

## 2026-08-06 — An analytic reference cannot see a parameter it shares with the code

The single biggest lesson from the analytic validation suite. A closed-form
solution is only a check on the parameters it does NOT share with the
simulation. Allen-Eggers gives peak deceleration as `Ve^2 |sin(gammaE)| / (2 e H)`.
Inflate `Hscale` by 15 percent and the prediction drops 13 percent — and so does
the simulation, by almost exactly the same amount, because Allen-Eggers *is* the
exact solution for that atmosphere. Measured: agreement moves from +10.50 to
+10.44 percent. The comparison is structurally blind to `H`, and nothing else in
the suite pinned it either (`test_missileConst` only brackets it to 6–9 km).

The fix is to anchor the shared parameter against something outside both sides.
Allen-Eggers assumes an isothermal atmosphere in *hydrostatic equilibrium*, and
for one of those the scale height is not free: `H = R T / g`. The library's
7200 m sits 1.64 percent from `Rair*T0/g0 = 7317.8 m`, so `test_allenEggers`
asserts that gap stays under 3 percent. That single line is what makes the
`Hscale` mutation fail; the deceleration comparison alone would have passed it.

The same trap, caught earlier in the same task: the equilibrium-glide reference
must read `CL` from `veh.CL`, not from `coorbital.aero.constLD`. Scaling `CL` by
1.2 *inside* `constLD` with the reference reading it back out moves the mismatch
from 14.2 percent DOWN to 9.9 percent — the mutant looks better than the
original. Reading `CL` from the vehicle instead, the same mutation drives the
error from 1.5 to 11.5 percent. Ask of every symbol in an analytic reference:
where does the simulation get this, and would a wrong value cancel?

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
