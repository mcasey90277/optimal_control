# Missile Trajectory Library — Code Organization

*Coorbital, Inc.*

MATLAB 3-DOF trajectory generation for hypersonic glide and ballistic vehicles
over a spherical, optionally **rotating** Earth. The library propagates
**prescribed-control** trajectories: give it a vehicle, a launch or entry state,
and a control schedule, and it integrates the flight and reports it. No
optimization — see [Out of scope](#out-of-scope).

Five milestones have shipped:

| Milestone | Delivers |
|---|---|
| 2026-08-06 | Validated glide propagator — `glide3DOF`, `phaseRun`, `HGV/run_glide` |
| 2026-08-07 | Powered boost, a descent phase, and full three-phase chains — `boost3DOF`, `constThrust`, `pitchProgram`, `massConstant`, `BM/run_ballistic` (boost → coast → impact), `HGV/run_boost_glide` (boost → glide → descent) |
| 2026-08-07 | **Point-to-point targeting and visualization** — `greatCircleBearing`, `rangeSolve`, the `+viz` package (`groundTrack`, `profilePlot`, `globe3D`, `globeMovie`), and `HGV/run_target`: give it a launch point and a destination and it solves the trajectory that connects them. See [Point-to-point targeting](#point-to-point-targeting) |
| 2026-08-07 | **Ballistic point-to-point targeting** — `BM/run_ballistic_target`, which ranges on the **loft angle** rather than on thrust termination and therefore has to deal with the two branches that come with it: a lofted arc and a depressed one for every range short of maximum. See [Ballistic point-to-point targeting](#ballistic-point-to-point-targeting-two-branches) |
| 2026-08-08 | **Two-axis targeting, and the closed-loop guidance seam** — `coorbital.util.aimSolve`, a damped Newton that solves the **launch azimuth beside a range control** against both components of the miss. Both targeting scripts now fly a **rotating Earth** and a **banked** track instead of refusing the first and warning about the second. `coorbital.guide.terminalConstraint` lands beside it as the interface PEG and VOA will both attach to; neither algorithm is implemented. See [Point-to-point targeting](#point-to-point-targeting) |

The design rationale is in `DESIGN.md`; the running log of what broke and why is
in `LESSONS_LEARNED.md`; `TODO.md` is the authority on what is *not* built. The
mathematics behind the targeting solve is `hgv_dynamics_note.tex` §8, and the
structural write-up of the two solvers is `software_design.tex`. This file
covers layout, use, and extension.

**Before you read a single trajectory number as a performance figure: every
vehicle and booster parameter in this library is a marked PLACEHOLDER.** See
[Placeholder parameters](#placeholder-parameters).

---

## Running

Everything below was executed against the code as committed. Every command shown
was run; every number shown came out of that run. The **targeting** figures, the
file counts and the suite result were re-measured on 2026-08-09 against two-axis
targeting; the `run_glide`, `run_ballistic` and `run_boost_glide` headlines are
from 2026-08-07, unmoved by that change and pinned to their last printed digit
by the suite; the movie figures are from 2026-08-07 and are marked where they
appear.

### The five entry scripts

| Script | Flight | Phases |
|---|---|---|
| `HGV/run_glide` | unpowered glide from an entry state | 1 (glide) |
| `BM/run_ballistic` | ballistic missile, pad to impact | 3 (boost, coast to apogee, descent) |
| `HGV/run_boost_glide` | boost-glide vehicle, pad to impact | 3 (boost, glide, terminal descent) |
| `HGV/run_target` | **launch point → destination point**, solved | 3 (boost to a solved cutoff, glide, descent) |
| `BM/run_ballistic_target` | **launch point → destination point**, solved, ballistic | 3 (boost to burnout or a solved cutoff, coast to apogee, descent) |

The first three fly *launch site + azimuth + control schedule → wherever the
physics puts it*. The last two invert that: they take two latitude/longitude
pairs and solve for the trajectory that connects them. Since 2026-08-08 both do
it the same way — the closed-form great-circle bearing is the **seed**, and
`coorbital.util.aimSolve` then solves the **launch azimuth beside a range
control** until both components of the miss are inside tolerance. They remain
separate scripts because the range control differs, and with it the whole
one-dimensional stage that produces the seed: `run_target` bisects the
thrust-termination time, which is monotonic in range, while
`run_ballistic_target` ranges on the **loft angle**, which is not, and so must
bracket and certify the max-range hump before it can bisect anything — which is
also why the ballistic one has two answers where the boost-glide one has one.

**What the second axis bought.** Both targeting scripts now fly a **rotating
Earth** and a **banked** track instead of refusing the first and warning about
the second: `earthSpin` true takes `run_target` from a 231551.628 m miss to
**4.361 m** and `run_ballistic_target` from 463211.19 m to **52.461 m**. Each
costs seven residual evaluations rather than the one a solve that converges at
its seed spends — six extra propagations — inside whole runs of 23 and 701. Both
figures, and what makes the second one large, are in
[Turning the Earth on](#turning-the-earth-on-what-it-costs-and-what-it-buys).

Each puts itself and the library root on the path, so it runs from anywhere
once MATLAB can see the file. Each is driven from a fenced
`%% USER PARAMETERS:` block at the top; nothing below that block should need
editing. Each takes an optional struct of overrides for named block entries, in
the block's own human units (degrees, kilometres), and a misspelt name raises
rather than silently leaving the shipped value in place:

```matlab
traj          = run_glide(struct('psiEntry',135,'latEntry',25,'showPlots',false));
[traj,info]   = run_ballistic(struct('separation',false,'showPlots',false));
[traj,info]   = run_boost_glide(struct('hHandoff',25,'showPlots',false));
[traj,info]   = run_target(struct('latTarget',30,'lonTarget',-140,'showPlots',false));
[traj,info]   = run_ballistic_target(struct('branch','depressed','showPlots',false));
```

The four chain scripts return a second output, `info`, carrying every number in
the printed summary at full precision. `run_boost_glide`'s and `run_target`'s
`info` additionally carry `phases`, `env`, `x0`, `boostVeh` and `glideVeh` —
everything an independent checker needs to re-integrate the same chain with a
different solver, which is exactly what `tests/test_fullChain.m` does.

Headless, which is how all five were verified:

```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles/HGV'); run_glide"

/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles/BM'); run_ballistic(struct('showPlots',false))"

/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles/HGV'); run_boost_glide(struct('showPlots',false))"

/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles/HGV'); run_target(struct('showPlots',false))"

/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles/BM'); run_ballistic_target(struct('showPlots',false))"
```

Either targeting script also takes `earthSpin` true, which since 2026-08-08 is
flown rather than refused:

```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles/BM'); run_ballistic_target(struct('showPlots',false,'earthSpin',true))"
```

### The glide entry script in detail

At the MATLAB prompt:

```matlab
cd ~/Desktop/optimal_control/missiles/HGV
run_glide
```

Called with an output — `traj = run_glide;` — it returns the trajectory struct;
called bare it prints the summary and returns nothing. All five scripts behave
this way.

The override struct exists so an automated test or a batch sweep can drive a
script at more than one operating point without editing it;
`tests/test_runGlide.m` is the reason it was added, and the throughput
requirement in `DESIGN.md` §1 is where it goes next. For a routine interactive
run, ignore it and edit the block.

### The test suite

```matlab
cd ~/Desktop/optimal_control/missiles
run('tests/run_tests')
```

Headless:

```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles'); run('tests/run_tests')"
```

**`run(...)` is required, not stylistic.** `matlab -batch "cd(...);
tests/run_tests"` does not work: MATLAB parses `tests/run_tests` as the
*expression* `tests / run_tests` and fails with `Unrecognized function or
variable 'tests'` before the harness ever loads. Verified — that is the exact
error it emits. The other working form is
`matlab -batch "cd(...); addpath('tests'); run_tests"`.

`run_tests` discovers every `tests/test_*.m`, runs it, prints `PASS`/`FAIL` per
file, and raises `missiles:testsFailed` at the end if any failed, so `-batch`
returns a non-zero exit code. Current state:

```
  PASS  test_aimSolve
  PASS  test_allenEggers
  PASS  test_boost3DOF
  PASS  test_boostEvents
  PASS  test_constLD
  PASS  test_constThrust
  PASS  test_equilibriumGlide
  PASS  test_expAtmos
  PASS  test_fullChain
  PASS  test_glide3DOF
  PASS  test_greatCircle
  PASS  test_greatCircleBearing
  PASS  test_missileConst
  PASS  test_phaseRun
  PASS  test_pitchProgram
  PASS  test_rangeSolve
  PASS  test_runBallistic
  PASS  test_runBallisticTarget
  PASS  test_runGlide
  PASS  test_runTarget
  PASS  test_saveFigure
  PASS  test_sphereGrav
  PASS  test_terminalConstraint
  PASS  test_viz

24 passed, 0 failed
```

305 s wall clock including MATLAB startup, measured 2026-08-09 on an Apple
silicon Mac, and zero warnings. The run was repeated and both were green.

**The suite got slower as it got sharper.** It was 42 s at 21 tests on
2026-08-07, and almost all of the increase is `test_runTarget` and
`test_runBallisticTarget`, which now fly whole two-axis targeting solves —
hundreds of trajectory propagations apiece — where they used to fly a handful.
That is the cost of covering the rotating and banked cases, and it is worth
paying: those two cases are the *only* coverage the two-axis path has from an
entry script.

Eighteen of the twenty-four are unit or analytic-validation tests — sixteen
exercise one library function apiece, and `test_allenEggers` and
`test_equilibriumGlide` check the propagator against a closed-form solution. The
other six test **compositions**, and they are the important ones:

| Test | What it composes |
|---|---|
| `test_runGlide` | `HGV/run_glide` end to end — unit conversion, the `greatCircle` call site, the peak search, the termination diagnosis, and the printed summary, which it parses rather than recomputes. Also guards the deliberate duplication between `vehicle_hgv` and `vehicleDefaults`. |
| `test_runBallistic` | `BM/run_ballistic` end to end — the boost → coast → descent chain, the staging link, and the Keplerian cross-check. |
| `test_runTarget` | `HGV/run_target` end to end — the **solved** azimuth, the bisection on cutoff time that seeds it, the separation link, the reachable-envelope refusal at BOTH ends of the band, and the printed summary. The miss it asserts is measured impact-to-target with `greatCircle` from the flown state, not read back out of the solver's own residual. It flies **three** configurations for that reason: the shipped zero-bank one, where the measured miss and the solver's residual agree to 9.3e-10 m and nothing can tell a measurement from a substitution; a **banked** case (part 6), 75° of terminal bank solved from a 21524.695 m seed miss to 54.795 m; and a **rotating** one (part 8b), flown rather than refused, 231551.628 m to 4.361 m. It also pins `aimIter == 0` on the shipped case — the correct assertion, and one that enshrines a coverage gap; see [What two-axis targeting left approximate](#what-two-axis-targeting-left-approximate). |
| `test_runBallisticTarget` | `BM/run_ballistic_target` end to end — the **solved** azimuth, the **certified bracketing of the max-range loft angle**, the two full-burn branch solves on the certified one-sided intervals either side of it, the **minimum-energy constrained minimisation**, and refusals on every path the script has: beyond maximum range, too close and an unreachable *branch* in part 8, the loft bracket in part 9, and the two-axis solve's own in part 15. Part 12 **flies a rotating Earth** and pins the deflection — where it used to assert a refusal. What it checks that no other test in this suite can: the branch is **re-derived from the flown loft angle against the certified maximiser interval** rather than read off the script's own label; the lofted and depressed arcs are asserted to differ by more than twice in apogee and 1.3 times in flight time (measured 2.55 and 1.48 at the vehicle's 6° clamp); the minimum-energy answer is asserted to be **below the burnout energy of both neighbouring feasible points and of both full-burn arcs**, by a margin far above the measured noise of the inner solve; and part 6b asserts the **vacuum equal-radius limit** — the same constrained minimisation applied to an impulsive burn at r = rE must return the classical γ\* and V\*, which it does to 1.6e−8 rad and 4.7e−16 relative. `assertBetweenArcs` is the assertion that a mode quietly falling back to *selecting* a branch cannot pass. Part 11 exercises the branch **measurement** on the one geometry where it can legitimately disagree with the label, and asserts `coalesced`. Part 13 flies **Plesetsk to New York**: refused on the shipped booster (7132.320 km against a 5211.525 km maximum) and, on a booster sized to the range, landing on the classical 1205.989 km / 26.012 min reference to +1.03 % and +4.42 %. It also pins the `alphaMax` story: the clamp is `BM/vehicle_bm`'s `alphaMaxDeg = 6`, a −40° `loftMin` must refuse with `coorbital:runBallisticTarget:maximumNotBracketed` for bracket width, and a 12° sensitivity study must cost 156.224 km of maximum range while buying reach down to 1684.117 km. |
| `test_fullChain` | The chain milestone's headline deliverable: boost, glide and descent on one seven-state vector. Re-integrates **across** each junction with `ode89` at 1e-12 — a different method at a hundred times the driver's tolerance — and asserts continuity in states 1–6 and the expected staging jump in state 7. |
| `test_boostEvents` | `eventBurnout` and `eventApogee` inside a live propagation, not just as scalar function calls. |

`test_runGlide` runs two propagations, the shipped due-east configuration and a
south-east one, because the shipped geometry is provably blind to a lat/lon
transposition at the `greatCircle` call site (with the entry point at the origin
the central angle is symmetric in the terminal latitude and longitude, so the
swap changes the summary by nothing at all). All four chain scripts ship an
off-axis geometry for the same reason, and both targeting tests prove their
geometry is discriminating — it asserts that every transposition of the `greatCircle`
arguments moves the range by more than fifty times the miss budget — rather
than assuming it.

`tests/test_aimSolve.m` is the solver's own coverage and it **integrates
nothing**: it drives `coorbital.util.aimSolve` against synthetic residuals with
closed-form roots, so the Jacobian, the SVD step, the conditioning ceiling, the
step cap, the line search and all seven refusals are exercised without a
propagation. It also pins the cost model — `nEval = 1 + 3*iterations + nShrink`
— and the stall exit, which is the exit guarding the whole "never a silent near
miss" contract and which had no coverage at all until 2026-08-08: deleting it,
gutting its message or ignoring `opts.maxHalve` entirely all left the suite
green. `tests/test_terminalConstraint.m` covers the five burnout constraints the
same way.

### Self-demos

Every library function except the three constant-return files runs a
demonstration when called with no arguments — `coorbital.atmos.expAtmos` plots the
density profile, `coorbital.aero.constLD` prints the coefficients,
`coorbital.prop.constThrust` prints thrust at sea level and in vacuum,
`coorbital.eom.massConstant` prints `dm/dt = 0.0 kg/s over a 7-state vector`.
Add the library root to the path first (`addpath ~/Desktop/optimal_control/missiles`).
This is the fastest way to see what a routine does without reading it.

---

## State and sign conventions

Read this before touching anything. These are the conventions a contributor is
most likely to get wrong.

```
x = [r, lon, lat, V, gamma, psi]      u = [alpha, sigma]       unpowered, 6-state
x = [r, lon, lat, V, gamma, psi, m]   u = [alpha, sigma]       powered,   7-state
```

| Symbol | Units | Meaning |
|---|---|---|
| `r` | m | **Geocentric radius**, not altitude. Altitude is `r - c.rE`. |
| `lon` | rad | Longitude, east positive |
| `lat` | rad | Geocentric latitude, north positive. **Singular at the poles** (`glide3DOF` errors rather than returning NaN). |
| `V` | m/s | Planet-relative speed. **Singular below 1 m/s** (errors). |
| `gamma` | rad | Flight path angle, **positive UP**. A descending entry has `gamma < 0`. Singular at ±90° (errors). |
| `psi` | rad | Heading, **clockwise from north**. 90° is due east. |
| `m` | kg | **Total mass currently carried** — payload plus whatever booster is still attached, propellant included. Present only in the seven-state vector. See [The seven-state convention](#the-seven-state-convention-and-the-inert-mass-trap), which you must read before building a chain. |
| `alpha` | rad | Angle of attack |
| `sigma` | rad | Bank angle. 0 puts all lift up; ±90° puts all lift into the turn. |

**The sign trap.** `coorbital.grav.sphereGrav` returns

```matlab
[gr,gLat] = sphereGrav(r,lat)
```

where **`gr` is a POSITIVE MAGNITUDE of inward attraction**, `mu/r^2` — it is
*not* a signed radial component. At the surface it returns `+9.7983`, not
`-9.7983`. The equations of motion are written to expect that: `glide3DOF` has
`Vdot = -aDrag - gr*sin(gamma) + ...`, so with `gamma < 0` on a descent the
gravity term *adds* speed, which is correct. Any replacement gravity model must
follow the same sign, and any new call site must subtract rather than add.
`gLat` is a signed **northward** acceleration.

Everything inside the library is **SI**: metres, m/s, radians, seconds, kg. The
only place human units appear is an entry script's `%% USER PARAMETERS:` block,
which converts to SI immediately afterwards in one place.

---

## Code organization

Everything reusable lives in one `+coorbital` package. Vehicle folders hold
**only** parameter files and entry scripts — no physics.

| Path | Contents |
|---|---|
| `+coorbital/+util/` | `missileConst.m`, `vehicleDefaults.m`, `boosterDefaults.m`, `greatCircle.m`, `greatCircleBearing.m`, `rangeSolve.m`, `aimSolve.m` |
| `+coorbital/+atmos/` | `expAtmos.m` |
| `+coorbital/+grav/` | `sphereGrav.m` |
| `+coorbital/+aero/` | `constLD.m` |
| `+coorbital/+eom/` | `glide3DOF.m`, `boost3DOF.m`, `massConstant.m` |
| `+coorbital/+guide/` | `prescribed.m`, `pitchProgram.m`, `terminalConstraint.m` |
| `+coorbital/+prop/` | `phaseRun.m`, `constThrust.m`, `eventAltitude.m`, `eventApogee.m`, `eventBurnout.m` |
| `+coorbital/+viz/` | `groundTrack.m`, `profilePlot.m`, `globe3D.m`, `globeMovie.m`, `saveFigure.m`, plus `private/` helpers |
| `HGV/` | `run_glide.m`, `run_boost_glide.m`, `run_target.m`, `vehicle_hgv.m` |
| `BM/` | `run_ballistic.m`, `run_ballistic_target.m`, `vehicle_bm.m` |
| `tests/` | `run_tests.m` plus twenty-four `test_*.m` |
| `docs/` | this file, `DESIGN.md`, `LESSONS_LEARNED.md`, `closed_loop_guidance.md`, the four plan briefs, the two LaTeX notes (`hgv_dynamics_note.tex`, `software_design.tex`), reviews |

Note the split: `+prop` holds both the **propagator** (`phaseRun`, the events)
and **propulsion** (`constThrust`). One package, two meanings of the word,
which is unfortunate but was not worth a rename; propulsion models are
identified by their signature, `[T,mdot] = f(t,P,veh)`.

### What each piece is for

**`util/missileConst`** returns a struct of nine constants — `muE`, `rE`,
`omegaE`, `rho0`, `Hscale`, `T0`, `Rair`, `gamAir`, `g0`. It is the **single
source of truth**; no routine anywhere hard-codes a physical constant. pumpkyn's
`getConst` was deliberately not used: it has no `muE`, no Earth radius, and two
fields known to be wrong.

**`util/vehicleDefaults`** returns mass, `Sref`, `CL`, `LD`, `noseRadius`. All
five are **PLACEHOLDER** open-literature values for a generic lifting body.
`noseRadius` is not read by anything: it is carried for the deferred
Sutton–Graves heating rate, and the file says so at the point of definition so
it is not mistaken for dead code.
`HGV/vehicle_hgv.m` starts from these and overrides — currently to the same
numbers, deliberately, so that real values have exactly one home when they
arrive. There is no `CD` field anywhere: `constLD` derives `CD = CL/LD`, so drag
cannot fall out of sync with lift.

**`util/greatCircle`** is the haversine central angle, chosen over the spherical
law of cosines because that form loses half its significant digits below about
1e-8 rad and terminal-phase ranges are short.

**`util/greatCircleBearing`** is the companion *initial bearing* of the same
arc, clockwise from north and wrapped into `[0,2π)` — which is exactly the
library's heading convention `psi`, so it drops straight into a launch state
with no sign flip. Since 2026-08-08 both targeting scripts use it as the **seed**
of the two-axis solve rather than as the answer. Bearing is **not** symmetric
under exchange of the two
points: LAX→JFK is 65.867° while JFK→LAX is 273.841°, not the 245.867° a naive
reversal would give. Its header states the three degenerate cases (coincident,
antipodal, polar) and what it returns at each.

**`util/rangeSolve`** bisects a scalar parameter until `fRange(x)` matches a
target value, treating `fRange` as an opaque, expensive black box — in the
intended use one call is a whole trajectory propagation. Bisection rather than
a secant or Newton method because there is no derivative, the function is only
piecewise smooth, and bisection cannot be thrown off a cliff by a local kink.
**An unreachable target is not an error**: it returns `converged = false` with
the achievable band in `info`, so the caller can phrase the refusal in its own
vocabulary. It never evaluates `fRange` twice at the same abscissa.

**`util/aimSolve`** is the two-dimensional sibling, and it is deliberately the
same contract: converge on the **achieved** residual rather than on the step,
never throw merely for failing to converge, and hand back an `info` struct rich
enough for the caller to write its own refusal. It drives a two-component
residual to zero in two controls — in the intended use, the launch heading `psi`
and a range control, against the down-range and cross-range miss in metres —
though nothing in the file knows that: `fResid` is any handle taking a `[2x1]`
and returning a `[2x1]`, treated as an opaque, expensive black box.

**The method had to be Newton, and that is a consequence and not a preference.**
`rangeSolve`'s entire safety argument is the bracket, and a bracket is a
one-dimensional object — the intermediate value theorem has no cheap planar
analogue. *Two coupled axes offer no bracket to preserve*, so every safeguard
bisection got for free is written out by hand: a condition-number ceiling read
off an SVD of the forward-difference Jacobian, tested **before** any step is
computed; a cap on the step applied as one scale factor on the whole step, never
component by component, which would rotate it off the Newton direction; and a
backtracking line search that treats a trial pair whose flight does not complete
as an **infeasible point** rather than as an error. The step comes out of that
same SVD — rotate the residual into the singular basis, divide by the singular
values, rotate back, `sDir = -(vJac*((uJac'*fCur)./sVal))` — rather than out of a
backslash: no second factorisation, no inverse formed anywhere, and no
near-singularity warning printed to the console from inside a loop a targeting
script runs thousands of times.

**The evaluation count is the cost, and it is part of the interface.** Each
residual evaluation is a full trajectory propagation. A clean iteration costs
**exactly three** — two difference probes and one accepted trial — because the
base point is never re-evaluated; the accepted trial value becomes the next base
value. So `nEval = 1 + 3*iterations + nShrink`, and `test_aimSolve` asserts it.
`dx0` is one step **per control** and a scalar is refused outright: a heading in
radians and a range control that may be an angle or a dimensionless fraction do
not share a scale, and silently broadcasting one step would bury a wrong
Jacobian column. `opts.maxStep` defaults to unbounded for the same reason in
reverse — a library function cannot know the caller's scales, and a wrong
default here would silently lengthen every solve. The callers set it;
`run_target` caps at 10° of heading and 5 % of the burn per iteration.

**Seven refusals, not one, because two of the remedies are opposites.**
`toleranceNotMet`, `singularJacobian`, `stepStalled`, `jacobianFailed`,
`stepUnderflow`, `jacobianNotFinite` and `stepNotFinite`. The last four all mean
loosely "the Jacobian could not be used", and collapsing them was considered and
rejected: `jacobianFailed` says **shrink** `dx0` and `stepUnderflow` says
**enlarge** it. A single message telling a caller to shrink when the fix is to
enlarge is worse than no message. The rule the library now applies generally:
*group failures by remedy, not by mechanism.* The one thing that does throw is
an unusable starting guess — a non-finite residual at `x0` raises
`coorbital:aimSolve:nonFiniteEval`, because there is no step to shrink there and
no earlier point to fall back on, which makes it an input error rather than a
failure to converge.

**`guide/terminalConstraint`** measures a burnout state against the five
terminal constraints **PEG and VOA both enforce** — radius, speed, flight path
angle, and the two components that fix the orbital plane — and returns a
**dimensionless** `[5x1]` residual that is zero exactly on the target manifold,
so a shooting or least-squares solver can drive it to zero without rescaling and
neither algorithm has to re-derive the five. It is the shared interface of the
closed-loop guidance seam; **neither algorithm exists**, and this file is the
only part of either that does. The mathematics is `hgv_dynamics_note.tex` §9,
"The five terminal constraints"; the design brief is
`docs/closed_loop_guidance.md`.

It takes **Cartesian position and velocity in the caller's own frame**, not a
library state, and that is deliberate rather than lazy. Converting the
seven-state spherical vector to an inertial Cartesian pair needs an inertial
frame, and this library has not defined one: with `env.omegaE` nonzero the
state's longitude is Earth-fixed and its speed is planet-relative, and there is
no epoch anywhere that says where the inertial x-axis points. Baking a guess
about that into the one interface both guidance laws share would put the guess
somewhere it could not be seen. The conversion belongs with the caller that owns
the epoch.

**`+viz`** is the plotting package: `groundTrack` (lat/lon track, one coloured
segment per phase, launch, target and impact marked), `profilePlot` (a chosen
subset of altitude, speed, Mach, dynamic pressure, load factor, mass and flight
path against time), `globe3D` (a still 3-D Earth with the trajectory arc) and
`globeMovie` (the same scene as an MP4 with the trajectory developing over
time). All four take the same `(traj,veh,env,opts)` argument list — several
of them do not read `veh` or `env` and say so in their headers — and every one
reads the trajectory and never writes it, so no figure can move a number in a
summary. **`globeMovie` is the one that styles itself dark**: `private/vizParent`
creates every figure on a **white** background, including `globe3D`'s still, and
`globeMovie` overrides it to black. It is also the only one that reaches for the
pumpkyn toolbox — the Earth texture and the starfield come from there when it is
on the path and degrade to a plain shaded sphere on black when it is not; both
paths are a fully working movie and need no network.

**`viz/saveFigure`** is the fifth public entry in the package and the odd one
out: it takes a figure handle and a path rather than a trajectory, and writes
the figure to an image file. It exists so the three decisions every caller needs
made — what resolution, what background, and what happens when the target
directory does not exist — are made **once** rather than in an inline
`exportgraphics` call per script with whatever arguments were to hand. 200 dpi
by default, an option because a slide and a paper want different answers; the
background defaults to **`'current'`, not white**, so a dark figure exports dark
and a light one exports light without the caller having to know which it has; a
missing directory is created rather than refused; the extension is optional and
defaults to `.png`, because callers hold a path *stem*; and the size on disk is
read back after the write, because `exportgraphics` returning quietly is not
evidence that pixels reached the file. It returns the absolute path actually
written. **It does not close the figure and modifies nothing about it** — that
is the same ownership contract the rest of the package keeps.

**`atmos` / `grav` / `aero` / `prop`** are the four swappable model families. See
[Adding a fidelity level](#adding-a-fidelity-level) — this split is the point of
the whole design.

**`eom/glide3DOF`** is the 6-state 3-DOF point-mass EOM (Vinh Eqs. 2.28–2.33).
It never names a model; it calls handles out of `env`. Rotating-Earth Coriolis
and centrifugal terms are written in from the start and gated on `env.omegaE`,
so setting that to zero recovers the non-rotating case *exactly* rather than
approximately. The three coordinate singularities (pole, vertical flight, zero
speed) raise identified errors rather than silently producing NaN.

**`eom/boost3DOF`** is the powered sibling: the same six equations plus a
seventh for mass. Thrust of magnitude `T` acts at angle of attack `alpha` from
the velocity vector, in the plane set by bank `sigma` — the same convention lift
already uses, so thrust and lift add directly. **The mass in every acceleration
denominator is the state `x(7)`**, so the vehicle feels itself getting lighter,
and `dm/dt = -mdot` closes the rocket equation. With `env.prop` returning zeros
the first six components reduce to `glide3DOF` bit-exactly.

**`eom/massConstant`** lifts a six-state EOM to seven states by appending
`dm/dt = 0`, and refuses to run if `x(7)` disagrees with `veh.mass`. Read
[the next section](#the-seven-state-convention-and-the-inert-mass-trap) before
using it.

**`prop/constThrust`** is a choked-nozzle rocket motor: constant mass flow with
an ambient back-pressure debit on thrust, `T = max(thrustVac - Aexit*P, 0)` and
`mdot = thrustVac/(Isp*g0)`. The flow is fixed because the throat is choked, so
back pressure changes the delivered thrust and not the flow — which is exactly
why deriving `mdot` from the *vacuum* thrust and the *vacuum* `Isp` is
self-consistent. When the clamp fires the model has left its domain of validity,
so `mdot` goes to zero **with** the thrust; a clamped engine that kept consuming
propellant would corrupt the mass history and fire the burnout event early.

**`util/boosterDefaults`** returns `massDry`, `massProp`, `thrustVac`, `Isp`,
`Aexit`, `Sref`, `CL`, `LD` for a single solid stage. All **PLACEHOLDER**. There
is deliberately **no `massWet` field**: it would have to mean either
booster-only or stack-including-payload, both readings are defensible, and the
ambiguity is precisely what the bookkeeping rule exists to remove. Add the two
masses that are there. The aerodynamic fields describe the **boosted stack**, a
slender body of revolution, not the waverider it carries.

**`guide/prescribed`** evaluates a control schedule `struct('tGrid',...,
'alpha',...,'sigma',...)`: linear inside the grid, clamped outside. Its
signature is `u = f(t,x,sched)` so a closed-loop guidance law drops in later
with the state already available.

**`guide/pitchProgram`** is the first guidance law that actually **reads the
state**. It interpolates a commanded pitch *attitude* `theta(t)` and returns
`alpha = theta - gamma`, taking `gamma` from `x(5)` — the classic prescribed-pitch
boost law, where the vehicle is flown to an attitude schedule and the angle of
attack falls out of the gap that opens between that attitude and the velocity
vector. Supply the optional `sched.alphaMax` whenever `theta` can run ahead of
`gamma`; a 90° commanded attitude at liftoff with `gamma` still near zero would
otherwise hand the equations of motion a physically absurd incidence. If the
field is absent, **no clamping is applied**.

**`prop/phaseRun`** is the driver. It integrates a `[1 x P]` array of phase
structs (`eom`, `guide`, `terminate`, `tspan`, and the optional `link`) with
`ode45`, hands the terminal state of each phase to the next, and records the
junction states so a Phase 2
optimizer can enforce them as linkage constraints. Default tolerances are
`RelTol = AbsTol = 1e-10`, overridable through `env.odeRelTol` / `env.odeAbsTol`.
Two behaviours worth knowing: each phase is referenced to its *own* `tspan(1)`,
so a phase given `tspan = [10 50]` contributes 40 s and opens no gap; and the
sample at a phase boundary is recorded once, carrying the *outgoing* phase's
control, so a control discontinuity does not appear in `traj.u`.

Two more were added in the boost-and-chain milestone. **Any state dimension:**
`nx` is taken
from `numel(x0)` and the control width `nu` is measured by calling phase 1's
guide once, before the loop, at `(tspan(1), x0)` — so a *stateful* guide sees one
evaluation more than the integrator makes and must tolerate it. And an
**optional per-phase `link`**, `xNext = link(xEnd)`, applied to a phase's
terminal state to form the next phase's initial state; that is how stage
separation drops the spent booster's mass. Absent or `[]` it is the identity
(`[]` is how a struct array mixes linked and unlinked phases, MATLAB requiring
every element to carry the same fields). The link is applied **before** the
junction is recorded, so `traj.junction(k).x` is the state on the **far** side of
a staging jump; the near side is the last `traj.x` row labelled phase `k`, at the
same instant on the cumulative clock.

**The three ODE events** are all terminal and all one-sided, so none can fire
from the wrong side:

| Event | Fires on | `direction` |
|---|---|---|
| `prop/eventAltitude(t,x,hStop)` | altitude descending through `hStop` | `-1`, so a lofted arc climbing back through it does not end the run |
| `prop/eventApogee(t,x)` | `gamma = x(5)` crossing zero descending | `-1`, so a still-climbing vehicle cannot trigger it. Reads only `x(5)`, so it works on six- and seven-state alike |
| `prop/eventBurnout(t,x,mBurnout)` | `x(7)` descending through `mBurnout` | `-1` |

`eventBurnout` takes the burnout mass as an **explicit scalar**, not reached for
inside the event from a booster struct. It must be `veh.mass + bst.massDry`:
comparing against `bst.massDry` alone would treat the payload as part of the
propellant budget and burn straight through it.

**`HGV/run_glide`** wires the single-phase glide together, propagates, and prints
a summary. It diagnoses *why* the propagation stopped, so a run truncated by the
time horizon can never be misread as a completed glide, and it flags when the
terminal Mach has dropped below the hypersonic regime in which constant `CL` and
`L/D` is defensible.

**`BM/run_ballistic`** and **`HGV/run_boost_glide`** are the chain scripts. Both
diagnose each phase's termination separately, both print a PLACEHOLDER caveat,
and both demonstrate the per-phase vehicle binding described next. `run_ballistic`
additionally carries a **Keplerian cross-check**: it re-propagates the burnout
state with aerodynamics off and compares against the closed-form free-flight
range of the same two-body arc, asserted at 1e-6 relative because both sides are
exact descriptions of one problem. It then re-propagates a *third* time with lift
alone suppressed, because switching `env.aero` off removes lift as well as drag
and the two have opposite signs — reporting their sum as "the drag effect" would
invite the wrong verdict whenever lift wins. Those differences are printed and
judged, never asserted.

---

## The seven-state convention and the inert-mass trap

**Read this before building a chain.** It is the hardest-won practical result of
the boost milestone.

### Why unpowered phases carry mass

`phaseRun` concatenates every phase into one `traj.x`, so **every phase in one
call shares one state dimension**. A powered phase needs seven states. Therefore
so does everything chained to it, including the unpowered coast, glide and
descent that carry no propellant and burn nothing.

The library solves this by **lifting the six-state EOM rather than mapping
between state definitions**: `coorbital.eom.massConstant(@coorbital.eom.glide3DOF)`
returns a seven-state EOM whose seventh derivative is identically zero. That is a
deliberate design change from `DESIGN.md` §2.2, which promised a
`+util/stateConvert.m` to translate at each junction. A uniform state vector
removed the need for it, so it was never written — recorded as a change, not left
as an unmet promise, in `DESIGN.md` §11.

### The trap

`glide3DOF` divides by **`veh.mass`**, not by the state mass `x(7)`. It is a
six-state routine and it is correct under its own contract: it never claimed to
read a mass state, because it has none.

Composed, though, there are **two sources of truth for one physical quantity and
nothing in the arithmetic makes them agree.** `x(7)` is carried, not read. A
chain that jettisons a stage through a `link` — dropping `x(7)` from the stack
mass to the payload mass — while continuing to hand the same `veh` struct to the
next phase will keep flying at the **old weight, silently**, with a mass jump
plainly visible in `traj.x` that the dynamics ignore. Nothing diverges, nothing
warns, and the trajectory looks entirely plausible.

Measured consequence in `run_ballistic`: deleting the separation link leaves the
flown trajectory **bit-identical** and changes only the recorded mass.

### The guard, and what it does not do

`massConstant` now refuses to run when `x(7)` and `veh.mass` disagree, to 1e-9
relative (never tighter than 1e-9 absolute — loose enough to clear the residual a
burnout-event solve leaves on the mass state, tight enough that a staging jump of
hundreds of kilograms cannot slip through). It raises
`coorbital:massConstant:massMismatch`. Verified — hand the wrapper a state mass
of 2400 kg with a 900 kg vehicle:

```matlab
f = coorbital.eom.massConstant(@coorbital.eom.glide3DOF);
f(0,[c.rE+40e3;0;0;4000;0;pi/2; veh.mass+1500],[0;0],veh,env)
```

raises `coorbital:massConstant:massMismatch` with the message:

```
The carried mass x(7) = 2400.000000000 kg disagrees with veh.mass =
900.000000000 kg by 1500 kg, over a 9e-07 kg budget. The wrapped equations of
motion divide by veh.mass and NEVER read x(7), so this flight would run
silently at the wrong weight. Rebuild the vehicle struct -- mass AND
aerodynamics -- after every staging event.
```

The guard deliberately does **not** repair the problem by injecting `x(7)` into a
copy of `veh`. That would leave `Sref`, `CL` and `LD` still describing the
jettisoned stack — worse than the bug it fixes, because it looks repaired.
Separation changes the whole airframe.

### What you must actually do

> **Any chain with a mass-changing link must bind a per-phase vehicle inside each
> EOM closure, rebuilding mass AND `Sref`, `CL` and `LD`.**

`phaseRun` carries one `veh` for the whole chain; that argument is therefore a
formality in a staged flight and is never read by the equations of motion. The
pattern, from `BM/run_ballistic.m`:

```matlab
eomFree = coorbital.eom.massConstant(@coorbital.eom.glide3DOF);
eomBst  = @(t,x,u,vehArg,envArg) coorbital.eom.boost3DOF(t,x,u,bst,envArg);
eomCst  = @(t,x,u,vehArg,envArg) eomFree(t,x,u,coastVeh,envArg);
```

`coastVeh` is built from the separation flag — `veh` when the booster is
jettisoned, `bst` when it is retained — and its `mass` field is then set to the
mass actually carried. **Both entry scripts demonstrate this**; `run_ballistic`
binds two vehicles (stack, re-entry body) and `run_boost_glide` binds two across
three phases (stack for boost; glide vehicle shared by glide and descent, since
nothing is jettisoned at the handoff — only the bank angle changes).

---

## Adding a fidelity level

**This is the architectural point of the library.** The equations of motion
never name a model. They call function handles carried in an environment struct:

```matlab
env.atmos  = @coorbital.atmos.expAtmos;
env.grav   = @coorbital.grav.sphereGrav;
env.aero   = @coorbital.aero.constLD;
env.prop   = @coorbital.prop.constThrust;   % powered phases only
env.omegaE = 0;
```

Raising fidelity is therefore a **one-line change in an entry script**, not an
edit to `glide3DOF`. Each family has a fixed signature, exactly as implemented:

| Family | Signature | Notes |
|---|---|---|
| atmosphere | `[rho,P,T,a] = f(h)` | `h` geometric altitude above the reference sphere (m); `rho` kg/m³, `P` Pa, `T` K, `a` m/s. Elementwise over `[N x 1]`. |
| gravity | `[gr,gLat] = f(r,lat)` | `r` geocentric radius (m), `lat` geocentric latitude (rad); `gr` positive-inward magnitude (m/s²), `gLat` signed northward (m/s²). Both `[N x 1]`. |
| aerodynamics | `[CL,CD] = f(alpha,mach,veh)` | `alpha` rad, `mach` dimensionless, `veh` the vehicle struct. |
| propulsion | `[T,mdot] = f(t,P,veh)` | `t` time since the start of the phase (s); **`P` ambient static PRESSURE (Pa), not altitude**; `T` delivered thrust (N), `mdot` positive propellant mass flow (kg/s) — the equations of motion apply the minus sign. Read only by `boost3DOF`. |

**Why the propulsion model takes pressure and not altitude.** Pressure is what
the physics needs: the back-pressure debit acts on the nozzle exit plane, and
`Aexit*(Pexit - P)` is a pressure term. Handing the model an altitude would force
it to hold its own atmosphere in order to get back to `P`, which would either
duplicate `env.atmos` or silently contradict it. The equations of motion already
have `P` in hand from `env.atmos`, so passing it costs nothing and keeps
propulsion **independent of which atmosphere is installed** — swap `expAtmos` for
US76 and every propulsion model follows automatically.

Earth rotation is not a handle — it is the scalar `env.omegaE`. Set it to
`c.omegaE` to turn on the Coriolis and centrifugal terms already present in
`glide3DOF` and `boost3DOF`, or leave it at 0. In every entry script the switch
is the boolean `earthSpin` in the user block, shipped `false`. In the two
**targeting** scripts, `HGV/run_target` and `BM/run_ballistic_target`, `true`
used to be **refused** outright, because their launch azimuth was a closed form
that is only the answer while the ground stands still. Since 2026-08-08 the
azimuth is solved and both scripts **fly** it. See
[Turning the Earth on](#turning-the-earth-on-what-it-costs-and-what-it-buys)
for what that costs and buys.

### Why `sphereGrav` returns a `gLat` that is always zero

Because **J2 does not**. An oblate gravity field has a genuine latitudinal
acceleration that a single "local up" scalar cannot represent. If the interface
had been a scalar `g = f(r)`, then adding `j2Grav` later would mean changing the
signature, and *every call site in `glide3DOF` would have to be rewritten* — in
three separate equations (`Vdot`, `gammadot`, `psidot`), each needing a new
projection of `gLat` onto the velocity frame, all under time pressure and with
no test that exercises the new term.

Instead the second output exists from the first commit, `sphereGrav` returns
`zeros(size(r))`, and the projection terms are already written into `glide3DOF`
and multiply out to zero under spherical gravity. Verified by injection: with an
otherwise identical state, replacing `env.grav` with a handle returning
`gLat = 0.01` and **no edit to `glide3DOF`** changes `psidot` by `-1.667e-6`
rad/s, exactly `-gLat*sin(psi)/(V*cos(gamma))`. The plumbing is live.

### Concretely: adding `j2Grav`

1. Write `+coorbital/+grav/j2Grav.m` with the signature `[gr,gLat] = j2Grav(r,lat)`,
   `gr` a **positive inward** magnitude and `gLat` signed northward, both
   elementwise over `[N x 1]`. Take `muE` and `rE` from
   `coorbital.util.missileConst` — do not hard-code them. J2 itself is a new
   constant, so add it to `missileConst` and **pin it** in `test_missileConst`
   (see the blind spot below).
2. Give it an `if nargin == 0` self-demo calling itself by full namespace.
3. Write `tests/test_j2Grav.m`. It must assert at least one **absolute** number
   that a wrong J2 would change — the equatorial-minus-polar surface gravity
   difference, say — not merely that `gLat` is nonzero.
3a. **Fly it on a NON-EASTWARD heading.** A due-east arc sets `cos(psi) = 0`,
   which zeroes the `gLat*cos(psi)` projections in `Vdot` and `gammadot` and
   leaves two of the three `gLat` terms in `glide3DOF` completely unexercised —
   the trajectory test would pass with either of them wrong, or missing. The
   shipped `run_glide` configuration is due east, so this cannot be inherited
   from it; override the heading. `tests/test_runGlide.m` flies a south-east
   case for exactly this reason and shows the pattern.
4. Change one line in the entry script:

   ```matlab
   gravFn = @coorbital.grav.j2Grav;
   ```

Nothing else moves. The EOM, the integrator, the guidance, the event functions,
and the plotting are untouched.

**One deliberate exception, which is not a counter-example.**
`BM/run_ballistic`'s Keplerian cross-check always re-propagates on
`coorbital.grav.sphereGrav` whatever `gravFn` is set to, because the closed form
it is compared against is built from `muE` alone and is two-body *by
construction*. The flown trajectory is oblate; only that one validation stays
Keplerian. Until 2026-08-07 the check inherited `env.grav`, so following the
recipe above aborted the script — measured at `2.9e-03` to `3.2e-03` relative,
on a 1e-6 assertion whose message insisted "a disagreement is a defect, not a
modelling difference". It was a modelling difference. Both entry scripts also
gate their printed validity paragraphs on the handles actually being the
defaults, so a swapped model no longer produces the line "the spherical
`j2Grav` carries no J2".

The same recipe with `[rho,P,T,a] = us76Atmos(h)`
adds a US76 atmosphere; the same recipe with `[CL,CD] = tableAero(alpha,mach,veh)`
adds a tabulated drag polar. Note that `alpha` currently has **no effect** on the
trajectory, because `constLD` ignores it by construction; it becomes live the
moment an alpha-dependent aero model is dropped in, and `run_glide` says so in
the user block.

---

## Point-to-point targeting

`HGV/run_target` takes a launch point and a destination point and solves for the
trajectory that connects them. Two unknowns have to be found — where to aim, and
when to stop thrusting — and there is exactly **one** geometry in which they
decouple. The library used to assume it. Since 2026-08-08 it solves them
together, in two stages.

**Stage one — the azimuth seed, and a bisection on thrust-termination time.**
`coorbital.util.greatCircleBearing` gives the initial bearing of the
launch-to-target great circle, clockwise from north, which is exactly the
heading state `psi`, so it drops into the launch state with no conversion. That
bearing is the whole answer under two conditions and no others: `env.omegaE = 0`
**and** every commanded bank zero. With both, the only horizontal force is a
lift vector held in the vertical plane, the ground track is the great circle the
vehicle left on, and matching the *distance* matches the *point*. Drop either
condition and the flown track is some other curve whose initial bearing is no
longer the aim point. It is still a good guess, and that is exactly how the
library now uses it: **as the seed**.

The range control bisected at that seed is the thrust-termination time. Of the
available control parameters, only this one is monotonic and single-valued:

| Parameter | Behaviour | Verdict |
|---|---|---|
| Thrust-termination time | less burn, less energy, shorter range | **chosen** |
| Pitch-program loft angle | two branches either side of a max-range hump | rejected — a root finder lands on whichever branch it started nearest |
| Glide handoff altitude | bifurcates on the glide phugoid's troughs; 30 km costs `run_boost_glide` 1881 km while every phase still reports nominal | rejected outright |

It also needs **no new machinery**: the boost phase already ends at `tspan(2)`
when the burnout event has not fired first, so a shortened `tspan` *is* a
commanded cutoff. `coorbital.util.rangeSolve` does the bisection; `run_target`
supplies the propagation.

**Early cutoff leaves unburned propellant, and it is thrown away.** At
separation the whole booster goes — dry structure and the propellant still in it
— so the post-separation vehicle is exactly the payload:

```matlab
        ph(1).link = @(x) [x(1:6); cfg.mGlide];
```

That is what lets one glide vehicle struct serve every cutoff, and it is what
`coorbital.eom.massConstant` checks on every derivative evaluation. Drop the
link and the carried mass is 17400 kg while the equations of motion divide by
900 kg; the guard raises `coorbital:massConstant:massMismatch` immediately
rather than flying the whole glide at nineteen times its weight.

**The reachable envelope is reported, not assumed.** The two bracket endpoints
are propagated before the solve begins, and their ranges *are* the envelope. A
target outside it is **refused in words**, with the band in kilometres and what
to change, at both ends:

```
    required range     18815.80 km    (great circle on the r = 6378.137 km sphere)
    reachable            202.86 to 7737.63 km
    shortfall          11078.17 km    (required minus the nearer edge of the band)
```

A refusal returns an empty trajectory and `info.refused = true`; it does **not**
throw, so a caller can read the band out of `info` without a `try`/`catch`. That
band is measured at the **seed** azimuth, so it is indicative for the azimuth
stage two settles on rather than exact — see
[What two-axis targeting left approximate](#what-two-axis-targeting-left-approximate).

### Stage two: two unknowns, two residuals, one damped Newton

**The residual is a miss vector resolved at the target**, not a range residual
beside a cross-track offset. Both components zero is exactly "the impact point
*is* the target", and the pair references nothing but two points — the target
and where the vehicle landed. The obvious alternative references a great circle
the flown track no longer follows, which is precisely what fails once the Earth
turns or a bank is commanded; it also couples both controls into both
components, where this form comes out **nearly decoupled** — cond(J) 8.07 to
8.57 on `run_target`'s two hard cases and 4.2 to 14.9 on the ballistic ones,
against a `condMax` ceiling of 1e8. The components are
built from `greatCircleBearing`'s own numerator and denominator rather than
from a bearing, because **no bearing exists at coincidence** —
`greatCircleBearing` rightly raises `degenerateArc` there — and a residual that
throws exactly where the solve converges is useless. Taken as components the
pair is smooth straight through coincidence and vanishes there, which is what a
residual must do at its own root. `hgv_dynamics_note.tex` §8 derives it.

`coorbital.util.aimSolve` then runs the damped Newton on the launch azimuth and
the range control together. `run_target` caps the step at 10° of heading and
5 % of the burn per iteration; the seed is already range-converged, so a step
anywhere near either cap means the Newton direction is wrong.

### Turning the Earth on: what it costs, and what it buys

Both targeting scripts fly `earthSpin` true. Until 2026-08-08 both **refused**
it, and before that `run_target` ran it anyway behind a caution and printed
three statements that contradicted each other around a number that was not a
miss. Measured 2026-08-09, all on the shipped geometries at a 1 km range
tolerance:

| Case | Miss at the seed | Miss achieved | Residual evaluations | Propagations, whole run |
|---|---|---|---|---|
| `run_target`, `earthSpin` true | 231551.628 m | **4.361 m** | 7 | 23 |
| `run_target`, 75° terminal bank | 21524.695 m | **54.795 m** | 4 | 18 |
| `run_target`, shipped | 511.243 m | 511.243 m, unchanged bit-for-bit | 1 | 14 |
| `run_ballistic_target`, `earthSpin` true | 463211.19 m | **52.461 m** | 7 | 701 |
| `run_ballistic_target`, shipped `'lofted'` | 779.491 m | **0.365 m** | 4 | 67 |
| `run_ballistic_target`, shipped `'depressed'` | 457.270 m | 457.270 m, unchanged bit-for-bit | 1 | 64 |
| `run_ballistic_target`, shipped `'minimum-energy'` | 39.009 m | 39.009 m, unchanged bit-for-bit | 1 | 733 |

**Read the cost column, because it is what decides whether rotation is
affordable.** A residual evaluation is one full trajectory propagation, and the
solve costs `1 + 3n` evaluations for `n` Newton iterations — so the hard rows
are two iterations and the 4-evaluation rows are one. The *rest* of each run is
stage one, and that is where `run_ballistic_target`'s cost lives: 62
propagations to bracket the hump and bisect both branches, against **693 to 731**
in `'minimum-energy'` mode. Two-axis targeting was cheap to add **here** because
it is eight propagations on top of seven hundred; that is not a claim that it is
cheap in general.

**Rotation deflects the vehicle. It does not move the target.** An earlier
edition of this file said the opposite — that the ground beneath the target
sweeps east under the vehicle, so the vehicle must be aimed where the target is
*going to be*, which would make the azimuth depend on the flight time and force
an outer iteration around the range solve. **That claim is withdrawn.** The
integrated state is planet-relative and the target is a ground point: both live
in the rotating frame and neither moves in it. What rotation does is bend the
vehicle off course, through the Coriolis and centrifugal terms that have been in
`glide3DOF` and `boost3DOF` from the first commit.

The distinction is not pedantic, because the two mechanisms predict different
**misses**: a target carried east under a fixed heading is largely a
*down-range* effect on an easterly course, while a deflected vehicle is
overwhelmingly a *cross-range* one. The decomposition settles it, on two
independent geometries:

| Rotating seed miss | Down-range | Cross-range | Ratio |
|---|---|---|---|
| `run_target`, 231551.628 m | −5858.645 m | +231477.499 m | 39.5 to 1 |
| `run_ballistic_target`, 463211.19 m | −17764.17 m | +462870.44 m | 26.1 to 1 |

Both scripts pin the ratio rather than the total, for exactly that reason —
`tests/test_runTarget.m` part 8b and `tests/test_runBallisticTarget.m` part 12.
The eastward-drift figure both scripts used to print (628 km on the rotating
`run_target` case) is **deleted**, computation and `info` field and all, because
it measures a mechanism the files no longer claim.

**Turning the Earth on also MOVES the reachable envelope**, and that is correct
physics rather than a solver artefact: an easterly launch is helped by the
planet's rotation. On the shipped ballistic geometry the maximum range rises
5211.5 → 5439.9 km and the depressed-branch floor 4708.5 → 5085.8 km, so the
shipped 4828.045 km target leaves the depressed band altogether and
`'depressed'` plus `earthSpin` is legitimately **refused**.

### Bank breaks the same assumption by the same route

A banked segment turns the heading, so the track leaves the arc it departed on
for a reason that has nothing to do with the planet. Rotation and bank are
therefore **one** failure of the closed form and not two, and one extra degree
of freedom closes both. Give `run_target`'s descent `run_boost_glide`'s 75°
terminal bank and the heading turns 166.28° between handoff and impact; the seed
misses by **21524.695 m**, of which 21515.222 m is cross-range, and the solve
brings it to **54.795 m** in four evaluations by aiming 0.3436° off the
great-circle bearing. `tests/test_runTarget.m` part 6 flies that exact case.

`run_target` still ships `descBank = 0` where `run_boost_glide` ships 75, but
**no longer for the reason it used to**. That zero was a targeting limitation
until 2026-08-08 — a banked descent was a case the script could only measure and
warn about. It is now simply the simplest default an example script should ship,
and the user block says so at the point of definition.

### What two-axis targeting left approximate

None of these is a defect. Each is a consequence of solving on a **surface**
with stage-one machinery whose guarantees are **one-dimensional**, and each is
stated in the code where it bites. Sizes measured 2026-08-09; `TODO.md` carries
the full list of seven.

**The envelope and the max-range bracket are computed at the SEED azimuth**, so
the too-far and too-close refusal gate is slightly approximate. Everything
`BM/run_ballistic_target`'s stage one certifies — the range hump, the certified
maximiser interval, the two monotone branch brackets, the reachable band — is
measured at one azimuth. Measured: **5.63° of aim moves the range hump 1.59°**,
about **0.28° of loft per degree of azimuth**, against a certified maximiser
interval **0.043° wide**; and turning rotation on moves the envelope by **30 to
42 km**. The alternative — re-bracketing the hump inside every Newton evaluation
— would multiply an already 700-propagation mode by the cost of a whole
bracketing and would still certify nothing, because the certification is a
one-dimensional argument.

**`'minimum-energy'` is minimised at the seed azimuth and only re-trimmed at the
solved one.** Stage two solves the azimuth beside the **cutoff fraction**,
holding the loft angle at the value the seed-azimuth minimisation settled on, so
the flown arc is not the minimum-energy arc at the azimuth it flies. Measured on
the rotating shipped run: minimised **−45.977595 MJ/kg** at cutoff 0.992294,
flown **−45.924475 MJ/kg** at 0.992621. The summary prints **two labelled
columns**, `MINIMISED (seed)` and `FLOWN (solved)`, differences the valley
evidence against the minimised one and says so, and reports the gap between
them as the size of the approximation. When the aim solve is a no-op — the
shipped non-rotating case — the two columns are exactly equal and the summary
says they are the same run rather than printing them twice in silence. Closing
this properly means minimising over the loft angle *at the solved azimuth*,
which is a nested solve and not a re-trim.

**The tolerance is sufficient but not necessary.** The convergence test is a
max-norm on the two components, so the per-component tolerance is the user's
range tolerance over √2 — 707.107 m for the shipped 1 km request. Components of
**900 and 100 m are a 906 m miss**, comfortably inside a kilometre, and they
**fail** the 707 m per-component test. The solver keeps working on a trajectory
the user would have accepted. That is the right way round — the alternative is
accepting one the user would not — but it should be read as a conservative test
and not as a measurement of the miss.

**And the shipped `run_target` case gives the two-axis path no coverage, and
structurally cannot.** On a non-rotating zero-bank run the cross-range component
is zero to rounding, so the miss *is* the down-range residual — which is the very
quantity the stage-one bisection has already driven inside the same tolerance.
The seed cannot be outside a tolerance the bisection was asked to meet, and
tightening it tightens the bisection in the same proportion: measured across a
**twentyfold** tightening, `aimIter` is 0 at every step. `tests/test_runTarget.m`
pins `aimIter == 0` there, which is the correct assertion and also enshrines the
gap. Coverage comes from the rotating and banked cases, and from
`tests/test_aimSolve.m`.

### The movie

`coorbital.viz.globeMovie` renders the trajectory developing over the Earth. It
is off by default because it is the expensive part of a run:

```matlab
run_target(struct('showPlots',false,'movieOn',true,'movieFrames',60, ...
                  'movieFile','/tmp/run_target.mp4'));
```

Run 2026-08-07, 60 frames at 1280×720 and 20 fps, twice. With the pumpkyn
toolbox on the path it took 13 s wall clock including MATLAB startup and wrote
6.4 MB over the `earth-clouds-4k` texture and the `starmap_4k` background;
without it, 2.1 MB over a plain shaded sphere on black, which is the documented
degradation and needs no network. Both are a working movie: the launch point,
the skipping glide and the terminal descent onto the target are visible in each,
and the three phases are coloured separately. Output goes to `tempdir` by
default — a movie is a build artefact and does not belong in the source tree.

#### Vertical exaggeration: true scale is the default

Both targeting scripts expose `altExag` in their USER PARAMETERS block, and both
ship it at **1 — true scale**. The globe used to exaggerate altitude because
otherwise the arc lay on the surface and showed nothing; the movie's **altitude
inset is true-scale already** and was added precisely so the globe need not lie
about it, and `globeMovie` drops the `(altitude exaggerated Nx)` caption clause
when the scale is 1, so a true-scale picture makes no claim it is not keeping.

The adaptive rule is still one word away: `altExag = 'auto'` selects `exagFor`,
which holds the apparent apogee at or under 0.3 rE, capped at 30 and floored at
2. Any positive number is used as given, and anything else raises — before the
solve, not after it.

Rendered and looked at 2026-08-07, 36 frames each at the new default:

| Run | apogee, as rendered | fraction of rE | true-scale verdict |
|---|---|---|---|
| `BM/run_ballistic_target`, lofted arc | 2118.45 km | 33.2 % | **reads well.** The arc stands clearly off the limb; the coast and descent legs are separately visible and the ballistic hump is obvious without any exaggeration. `'auto'` gives 2x, a barely perceptible difference |
| `HGV/run_target` | 118.09 km | 1.9 % | the ground track is crisp against the globe with both markers legible, but the **skip phugoid is invisible** on the sphere — it reads as a surface track. The inset carries it completely: six clear skips decaying from 121 km, in true scale, with axes. Ask for `'auto'` (16x) when the phugoid is the point of the picture |

Neither is unreadable, so neither default was reverted. What changed is where a
reader looks for altitude: the inset, which is honest, rather than the globe,
which was not.

**Read the first row's apogee as the arc that was rendered, not as the arc that
ships.** `BM/run_ballistic_target`'s demonstration target and angle-of-attack
clamp both moved on 2026-08-08 (see
[The clamp is a VEHICLE limit](#the-clamp-is-a-vehicle-limit-and-it-is-not-a-targeting-knob)),
and the shipped lofted arc now apogees at **1454.17 km, 22.8 % of rE** —
re-measured 2026-08-09. That is a third lower than the arc in the picture, and
the verdict above has **not** been re-checked at the new height. The `'auto'`
rule still returns its 2x floor there, since 22.8 % is already inside the 0.3 rE
target. `HGV/run_target`'s row is unaffected: its shipped run is bit-for-bit
what it was, apogee included.

### Saving the figures: `plotFile`

Until 2026-08-09 nothing in this library could write a figure to disk. Every
still that was in `results/` before that date was produced by an ad-hoc scratch
driver, not by an entry script — no entry script could have produced one.

`BM/run_ballistic_target` now exposes **`plotFile`** in its USER PARAMETERS
block, a path **stem without an extension**, shipped at
`fullfile(tempdir,'run_ballistic_target')` for the same reason `movieFile`
defaults there: a picture is a build artefact and does not belong in the source
tree. Each figure appends its own suffix and `.png`, and the three suffixes are
**fixed**, naming what each figure *is*:

| Figure | File |
|---|---|
| `coorbital.viz.profilePlot` | `<stem>_profile.png` |
| `coorbital.viz.groundTrack` | `<stem>_ground_track.png` |
| `coorbital.viz.globe3D` | `<stem>_globe.png` |

Fixed rather than derived from each figure's title, so a re-run **overwrites**
its own pictures instead of accumulating a pile of near-duplicates: a title
carries a miss distance, a minus sign and a degree symbol, and the sanitised
result is a different, truncated name on every run — which is exactly what the
`china_to_la_ground_track_10802_km_requi.png` left in `results/` by the old
scratch driver looks like.

```matlab
stem = '/Users/msc/Desktop/optimal_control/missiles/results/china_to_la';
[traj,info] = run_ballistic_target(struct('plotFile',stem));
info.plotFiles   % 1x3 cellstr: .../china_to_la_profile.png, _ground_track.png, _globe.png
```

Run 2026-08-09, that writes the three files and prints them, and the shipped
miss is unmoved at **39.0092687975735 m** — `saveFigure` reads a figure and
writes a file, and touches nothing the summary is computed from.

Saving happens whenever `showPlots` is true **and** the stem is non-empty; set
`plotFile` to `''` to draw the figures on screen without writing them. There is
deliberately **no second on/off flag** — the pattern is the one `movieFile` and
`movieOn` already establish in the same script, where the file name and the
switch are one decision rather than two that can fall out of step. `plotFile` is
in the override whitelist, so a mistyped option raises rather than silently
doing nothing. The paths are printed at the end of the run in the same style as
the existing `Movie:` line, and returned on `info.plotFiles` — present only when
the figures were drawn *and* saved.

**`HGV/run_target` has the same gap and is deliberately not wired.**
`coorbital.viz.saveFigure` was written for both callers and its header names
both; wiring the second is a scheduled follow-up in `TODO.md`, not an oversight.

---

## Ballistic point-to-point targeting: two branches

`BM/run_ballistic_target` takes a launch point and a destination and solves the
**ballistic** trajectory that connects them. Stage two is the same as
`HGV/run_target`'s in every respect — the same miss-vector residual, resolved in
the same frame, driven by the same `coorbital.util.aimSolve`. Stage one is not,
and the difference is the whole reason this is a separate script.

**The ranging control is the loft angle, and range is not monotonic in it.**
`run_target` bisects on thrust-termination time because less burn means less
energy means less range — one answer, and bisection is safe. A ballistic missile
is flown to a **loft angle** instead: the terminal attitude of the pitch
program. Range rises to a maximum at some max-range loft angle and falls away on
both sides, so **every range short of that maximum is reached by two
trajectories**:

| | depressed arc | max-range arc | lofted arc |
|---|---|---|---|
| loft angle | −126.73° | **−42.91°** | +19.43° |
| range | 4828.50 km | 5211.53 km | 4828.82 km |
| apogee | 569.18 km | 974.10 km | 1453.66 km |
| flight time | 1047.95 s | 1314.78 s | 1550.50 s |
| impact speed | 1781.32 m/s | — | 3012.98 m/s |
| impact angle | −22.97° | — | −43.40° |

Measured on the shipped geometry: 45°N 100°W to 62°N 28°W, 4828.045 km required
on a 40.555° **seed** azimuth, at the vehicle's 6° angle-of-attack clamp. **Both
arcs are always solved and always reported**, and only one is flown, because the
trade is the point — but only the flown one is then refined by the two-axis
solve, so the table above is the trade at the seed azimuth and the flown arc's
own solved numbers are reported separately. The depressed arc's commanded
−126.73° is not a pitch programme anyone would fly: it is where that branch goes
once the clamp saturates the achievable pitch-over, which is why `loftMin` ships
at −140°.

So the script **brackets the maximum first** — a 21-point coarse scan across the
loft bracket, adaptively refined until it *certifies* a single interior hump,
then golden-section refinement to 0.05°. What that returns is an **interval**,
`[aL,bL]`, and not a maximiser: a golden section can prove only that the
maximiser is inside its final bracket. So the loft axis is split on the
interval's **ends** — `[loftMin, aL]` and `[bL, loftMax]` — on which range
provably is monotonic (the precondition `coorbital.util.rangeSolve` documents and
does not check). Splitting at the interval's midpoint instead, which this script
did until 2026-08-08, leaves one bracket straddling the true maximum.

**All of stage one runs at the seed azimuth, and it has to.** Range against loft
is a *curve* at a fixed azimuth and a *surface* otherwise, and a golden section
on a surface certifies nothing — every guarantee the bracketing makes is a
one-dimensional argument. The consequence is item 1 of
[What two-axis targeting left approximate](#what-two-axis-targeting-left-approximate),
and it is why the stage-two Newton is **not** boxed to its branch's side of the
interval: the only interval available to draw the box on is the seed-azimuth
one, and at 0.28° of loft per degree of azimuth against a 0.043° interval a box
drawn there would exclude genuine roots. Three layers of *diagnosis* stand in
place of that prohibition — `opts.maxStep` caps the travel per iteration, the
branch is measured against the interval and the summary prints a caution naming
stage 2 when the two disagree, and a Newton that cannot converge refuses
outright rather than returning the arc it wandered onto.

**The branch is measured from the root's POSITION against that interval**: below
`aL` is depressed, above `bL` is lofted, inside is **coalesced** — at that range
the two arcs are the same arc and the distinction has no content. It used to be
measured from the flown apogee and flight time against the max-range arc's, and
neither of those is a branch invariant for a finite powered atmospheric arc:
drag, lift, burnout altitude and boost duration can make either non-monotone in
the loft angle, and near the maximum both differences vanish quadratically. They
are still printed, as the descriptive quantities they are.

**The steep arc arrives faster here, which is the opposite of the vacuum
intuition.** With the propellant load fixed both arcs leave burnout with almost
the same energy — −44.54 and −44.82 MJ/kg, 0.63 % apart — so in vacuum they
would arrive at almost the same speed. With an atmosphere the shallow −23.0°
depressed descent spends far longer in dense air than the steep −43.4° lofted
plunge and is braked to 59 % of its arrival speed. The script measures this and
prints it; it does not assert the textbook expectation.

### The three modes: two full-burn arcs, and a constrained minimisation

Two of the three fly the **full burn** and range on the loft angle alone —
`'lofted'` and `'depressed'`, one either side of the max-range loft angle. They
are what a fixed booster with no thrust termination can do.

The third states an optimisation problem and then solves it:

```
minimise    eps_BO = V_BO^2/2 - mu/r_BO
over        the loft angle and the cutoff fraction
subject to  R(loft,cutFrac) = the required range
```

Two parameters and one constraint leave a **one-dimensional feasible family**,
so the problem is well posed: parameterise the family by the loft angle, let the
cutoff be whatever makes the range, and minimise the burnout specific energy
along it.

**It used to be a gamma match, and that was the critical finding of the
2026-08-08 review.** The classical closed form

```
V*^2      = (mu/rE) * 2 sin(Lambda/2) / (1 + sin(Lambda/2))
gammaStar = 45 deg - Lambda/4
```

is derived for a free-flight arc whose two endpoints lie at the **same radius**,
with Λ the free-flight central angle between them. The script drove the burnout
flight-path angle to that γ\* with Λ taken as the **pad-to-target** angle, while
burnout sits downrange and 82 km above the impact sphere. The residual therefore
did not apply to the arc being flown, and the reported agreement with it verified
the wrong condition. γ\* is still computed and still printed — it is the right
yardstick for how far a finite boost is from the vacuum equal-radius idealisation
— but it is a **diagnostic beside** the achieved burnout gamma, not a residual.

**Neither full-burn arc is the answer.** At full burn the booster delivers a
fixed delta-V, so both arcs leave burnout with essentially the same energy —
0.63 % apart on the shipped case — and there is no minimisation to do. So
`'minimum-energy'` carries a **second control**: the burn is **cut short**, using
the same thrust-termination mechanism `HGV/run_target` ranges on.

| level | parameter | what it does | how it is bracketed |
|---|---|---|---|
| inner | cutoff fraction | enforces the **constraint** to `tolRangeMEKm` | the cutoff axis is *sampled* and a **sign change** of `R − R_req` is taken; monotonicity is then *verified* on the interval selected, not assumed across `[cutFracMin, 1]`. Where several roots exist the lowest-energy one wins |
| outer | loft angle | minimises **eps_BO** along the feasible family to `tolLoftMEDeg` | a coarse scan of `nScanME` points certifies a single valley — both ends of the family are full-burn arcs, so the energy is highest there — then a golden-section minimisation |

The **outer bracket is the two full-burn branch solutions**, already solved and
therefore free: between them and only between them the full burn reaches at
least the required range. Where a branch does not exist the user's own loft limit
serves at that end, which is admissible for the same reason — the full burn
*overshoots* there. An endpoint may itself be infeasible by a few hundred metres,
because `rangeSolve` places a branch within `tolRangeKm` **either side** of the
range; the scan tolerates that at the ends and refuses a hole in the middle.

**It is shown to be a minimum, and the noise is measured.** The summary prints
the burnout energy at the two neighbouring **feasible** points, and re-solves the
constraint ten times tighter at the settled loft angle to report the propagated
effect of the inner tolerance on the objective. On the shipped case:
−45.555275 MJ/kg against −45.486364 and −44.818885 at the neighbours, a valley
68 911 J/kg deep on its shallower side, against 101 J/kg of measured noise — a
ratio of 1.5e−3. Cost: 5 scan points and 16 golden-section steps, 669
propagations of a 733-propagation run, about 28 s. **All of that evidence is at
the seed azimuth**, differenced against the `MINIMISED (seed)` column and not
against the flown one — see
[What two-axis targeting left approximate](#what-two-axis-targeting-left-approximate).

**And the objective reduces to the classical arc in the vacuum equal-radius
limit**, which is the check that says it is the right objective. Take the boost
to be impulsive at r = rE in a vacuum; the feasible family is then the
one-parameter family of Keplerian arcs of central angle Λ from rE back to rE, and
minimising the same eps returns γ\* to **1.6e−8 rad** and V\* to **4.7e−16
relative**. `test_runBallisticTarget` part 6b writes that limit out in its own
hand and asserts it, on two different range angles.

**Against the classical arc on the shipped case**, as a diagnostic:

| | flown | classical | difference |
|---|---|---|---|
| apogee | 965.613 km | 952.696 km | **+1.36 %** |
| flight time | 21.2817 min | 20.1686 min | **+5.52 %** |
| burnout γ | 33.3290° | 34.1572° | **−0.828°** |
| burnout speed | 5682.3 m/s | 5807.2 m/s | −2.15 % |

**Those differences are physics, not solver error, and they are not residuals.**
The classical result assumes an *impulsive* burn *at* the impact radius in a
*vacuum*. This flight burns for 80.1 s and finishes 82.18 km up, and a Keplerian
arc from *that* burnout state apogees at 965.641 km — the flown figure to
0.027 km. Coast drag supplies the remainder.

**Minimum-energy requires `separation = true`.** A cut-short burn leaves
propellant in the booster, and the coast vehicle's mass is fixed before the
cutoff is known, so the *whole* booster goes overboard at cutoff — 150.7 kg of
unburned propellant with it on the shipped case, which sits at 93 % of maximum
range where there is little energy left to give back. The 3175 km case in part 7
cuts at 0.961 and throws away 1170.6 kg. `separation = false` cannot express that
and is refused for this mode.

**And it refuses what it cannot reach.** Plesetsk to New York is 7132.320 km
against this placeholder booster's 5211.525 km maximum, so it is declined
through the envelope path — and the refusal says *why* minimum-energy is the
harder ask, not the easier one: every trajectory it will consider carries less
energy than the full burn, and a burn cut shorter cannot fly further than the
full burn does. Give the same geometry a booster with the delta-V for the range
(52 t of propellant against 30 t, everything else unchanged) and the mode lands
on the reference: **1218.445 km of apogee against the classical 1205.989 km
(+1.03 %), 27.161 min against 26.012 min (+4.42 %), burnout γ 0.557° from γ\***.
That case is flown in `test_runBallisticTarget` part 13.

### The clamp is a VEHICLE limit, and it is not a targeting knob

`BM/vehicle_bm.m` carries `alphaMaxDeg = 6`, and **both** BM entry scripts read
it. Until 2026-08-08 `run_ballistic_target` carried 12° in its own user block
where `BM/run_ballistic` carried 6°, for nominally the same airframe — and the
12° had been chosen because it brought the demonstration target inside the
depressed branch. A control-authority limit chosen for reachability is not a
limit, and two limits for one vehicle made the two scripts' performance
non-comparable. Either script still accepts an explicit `alphaMax` override, for
a deliberate sensitivity study.

**6° is a PLACEHOLDER awaiting a qualification basis**, like every other number
in `vehicle_bm`. It is the value `run_ballistic` has always flown and the one
that was *not* chosen to make a feature work. What each value costs, measured:

| `alphaMax` | max-range loft angle | maximum range | depressed branch |
|---|---|---|---|
| 6° (the vehicle's) | **−42.907°** | **5211.525 km** | spans 4708.463–5211.525 km |
| 12° (sensitivity study) | +25.068° | 5055.302 km | spans 1684.117–5055.302 km |

Raising the clamp therefore **buys depressed-branch reach and pays about 156 km
of maximum range**. Neither figure is a reason to move it; the qualified number
is whatever the airframe is cleared for, and it does not exist yet.

Two consequences follow for the shipped configuration.

**`loftMin` ships at −140°, because that is what it takes to bracket the hump.**
The max-range angle at 6° sits 2.907° *below* the −40° this script used to ship,
so the coarse scan would find its largest range on an endpoint and `maxRangeLoft`
would refuse with `coorbital:runBallisticTarget:maximumNotBracketed`. The
refusal, and the bracket width that causes it, are both pinned in
`tests/test_runBallisticTarget.m` part 9.

**The demonstration target sits inside the 6° depressed band**, 4828.045 km
against a band of roughly 4708–5212 km, which is what lets all three modes fly on
the shipped configuration. It moved there from 3174.981 km when the clamp did:
at 6° the depressed branch cannot reach 3175 km at any loft angle, and that
geometry is now part 7's second minimum-energy case and part 8's
unreachable-branch refusal.

### Four refusal paths, none of which throws

`info.refusedWhy` names which gate stopped the run. Three of the four are
envelope refusals raised before anything is flown to the target:

| Case | What the summary says |
|---|---|
| Beyond maximum range | `BEYOND MAXIMUM RANGE by 7112.078 km`, with the 5211.525 km maximum, the −42.9074° loft angle achieving it, the 556.603–5211.525 km band, and each branch's own band and loft interval |
| Too close | `TOO CLOSE`, with the overflight in kilometres and what to change |
| **Branch unreachable** | the target is inside the envelope but only the *other* arc reaches it — `The depressed arc does not reach this target, though the other one does.` |

The third has no counterpart in `run_target`. Silently flying the branch that
does reach the target would be the worst outcome available: it converges, it
hits, and it is not the trajectory that was asked for. All three carry
`refusedWhy = 'envelope'` or `'minimumEnergy'`, return an empty trajectory with
`info.refused = true`, and none throws, so a caller reads the band out of `info`
without a `try`/`catch`.

The fourth and fifth arrived with two-axis targeting and with the certified
bracketing respectively. **`refusedWhy = 'aimSolve'`** fires when the stage-two
Newton cannot reach the tolerance: the run is refused with the best point it
found rather than returned as a near miss, which is the whole contract
`coorbital.util.aimSolve` was written to keep. **`refusedWhy = 'loftBracket'`**
fires when the coarse scan's very first propagation raises
`coorbital:pitchProgram:unreachableAttitude`.

**One caution for anyone writing a caller.** `'loftBracket'`'s `info` is *not* a
superset of the other three, and cannot be: it fires before `loftStarDeg`,
`rngMaxM`, `rngMinM`, `maxRange`, `classical` and the per-branch records exist.
It carries `bankAngleDeg`, `loftMinDeg`, `loftMaxDeg`, `alphaMaxDeg` and
`libraryErr` instead. This is deliberate and pinned by a test, but a caller
written against the other three as a template will break on it. The four are not
interchangeable.

`HGV/run_target` has two of the same paths, `'envelope'` and `'aimSolve'`, on
the identical contract.

### Per-phase vehicles, carried on `ph.veh`

This is the first entry script to use `coorbital.prop.phaseRun`'s optional
per-phase `ph.veh` rather than binding a vehicle inside an EOM closure and
ignoring the forwarded argument. The equations of motion are the library handles
themselves:

```matlab
         ph(1).eom = @coorbital.eom.boost3DOF;
         ph(1).veh = cfg.bst;
         ph(2).eom = cfg.eomCoast;          % massConstant(@glide3DOF)
         ph(2).veh = cfg.coastVeh;
```

`coorbital.eom.massConstant` still guards the mass half of the divergence on
every derivative evaluation, so a separation link that moved the state without
the vehicle following it still raises immediately.

---

## Validated results

All numbers below were produced by running the code, not copied from a report.

### The glide milestone

Re-run 2026-08-06.

| Check | Reference | Measured | Budget |
|---|---|---|---|
| Vacuum specific energy `V²/2 − μ/r` over a 300 s lofted arc | conserved exactly | `1.28e-13` relative drift | `< 1e-8` |
| Equilibrium glide `V(r)`, 45 → 20 km | closed form from `veh.mass`, `Sref`, `CL` | **1.51 %** worst case over the whole arc | `< 3 %` |
| Allen–Eggers peak deceleration, β = 8571 kg/m² | `459.85 m/s²` | `508.12 m/s²`, **+10.50 %** | `0 < rel < 12 %` |
| Allen–Eggers peak deceleration, β = 1071 kg/m² | `459.85 m/s²` | `502.16 m/s²`, **+9.20 %** | `0 < rel < 12 %` |
| β-independence of the peak (8× change in β) | must not move | `1.18 %` apart, while the peak *altitude* moves `15.01 km` (4.86 → 19.86 km) | `< 3 %` |
| Allen–Eggers peak *location*, in density | `rho_pk = beta*abs(sin(gammaE))/Hscale` | `+4.85 %` and `+4.36 %` | `< 6 %` |

The equilibrium-glide arc is not a soft one: density spans a factor of 32.2 and
speed a factor of 4.80 across it, and the trajectory has to track a moving
analytic curve the whole way.

The two analytic checks are **complementary, not redundant**, and neither
validates "the propagator" alone. `Veq²(r)` has no drag term, so equilibrium
glide is blind to `CD` — scaling drag by 0.8 leaves it passing while
Allen–Eggers fails by +38.3 %. A zero-lift Allen–Eggers entry says nothing about
`CL`. The Allen–Eggers tolerance is 12 % rather than the spec's 2 % because the
closed form drops gravity and freezes the flight path angle, which at a 30°
entry accounts for a predicted +9.8 % against the measured +10.50 %; it is also
one-sided **below**, since both neglected effects can only *raise* the peak, so a
simulated peak under the analytic value means too much drag, not a better
approximation.

### The boost and chain milestone

Re-measured 2026-08-07 by re-running the checks in `tests/test_boost3DOF.m` and
`tests/test_fullChain.m` outside the harness, so the figures below are this
document's own measurements and not the tests' comments.

| Check | Reference | Measured | Budget |
|---|---|---|---|
| `boost3DOF` reduces to `glide3DOF` with the engine dead, 4 states (equatorial/off-equator/rotating/steep-climb) | `glide3DOF` on the same state | **`0.000e+00`** — `isequal` true on all four | `< 1e-12` |
| **Thrust force increment**, 2 states, both angles of attack and both banks nonzero | `T cos α/m`, `T sin α cos σ/(mV)`, `T sin α sin σ/(mV cos γ)`, `−mdot` | `3.642e-16` worst relative | `< 1e-12` |
| Tsiolkovsky, full vacuum burn, no gravity, thrust along velocity | `Isp g0 ln(32400/2400)` = `6636.153368978 m/s` | `6636.153370556 m/s`, `2.377e-10` relative | `< 1e-8` |
| Mass linearity over the burn | `m = 32400 − 372.5886162804 t`, 80.518 s burn | `6.821e-15` worst relative | `< 1e-9` |
| Boost→glide junction, boost-glide chain, radius channel | `ode89` @ 1e-12 vs the driver's `ode45` @ 1e-10 | `1.340e-06 m` two seconds past the seam | `1e-3 m` |
| Glide→descent junction, radius channel | same | `3.120e-07 m` | `1e-3 m` |

Two things to understand about that table.

**The reduction is bit-exact, and that is worth exactly what it says.** Every
shared term in `boost3DOF` is a character-identical copy of the `glide3DOF`
expression evaluated in the same order, so with `T = 0` the residual is an IEEE
identity. It *proves* the six shared equations **are** `glide3DOF` — no
transcription drift, no dropped term, no altered sign — and it is simultaneously
**zero-bit evidence of correctness**, because any error in `glide3DOF` is
inherited and cancels invisibly. Mutations of the thrust terms survived the
entire suite before the force-increment check existed — **four of them measured
and tabulated** in the 2026-08-07 reduction-test entry of `LESSONS_LEARNED.md`,
eight in total according to the development log, which is not in this repository.
Take the four as the reproducible figure; the eight is provenance, not a
measurement, and unlike everything else in this section it cannot be re-run now
that the check exists. The force increment is
what covers them: evaluate the same state twice with **only** `env.prop` changed,
so every shared model is bit-identical across the two calls and the difference
isolates the new terms exactly. It is the sharpest instrument in the suite.

**Tsiolkovsky is blind to every booster constant.** `Isp`, `g0`, `thrustVac` and
the mass ratio all cancel between the propagation and the closed form — doubling
the thrust halves the burn time and leaves the delta-V untouched. It validates
the *shape* of the rocket equation and nothing about the numbers fed into it.
**Seven** exact-equality pins at `tests/test_constThrust.m:51-57` are the only
defence. Same disease as the `Hscale` blindness below; same cure.

**The junction reference is tolerance-ordered, not solver-named.** What makes
`ode89` @1e-12 a valid reference is that it is a hundred times tighter than the
driver, and `test_fullChain` asserts that ordering rather than leaving it in a
comment, so a later edit cannot downgrade the reference — or tighten the driver
through `env.odeRelTol` — and leave the check silently vacuous.

### What `run_glide` produces

With the shipped user block — 60 km entry, 6000 m/s, γ = −1°, due east from the
equator, zero bank, terminating at 5 km:

```
flight time         2073.77 s
ground range        6986.82 km   (great circle, r = 6378.1 km sphere)
central angle       62.7636 deg
terminal speed        321.01 m/s (Mach 1.01), gamma -12.60 deg
peak dynamic press.    34.59 kPa (t = 1178.6 s, h = 37.20 km)
peak deceleration       1.11 g   (10.87 m/s^2 sensed aero load)
samples                  3273    (ode45 adaptive)
```

Note that the run ends at Mach 1.01, well outside the hypersonic regime in which
holding `CL` and `L/D` constant is defensible — `run_glide` prints that warning
itself. And see [Placeholder parameters](#placeholder-parameters).

### What `run_ballistic` produces

Shipped block — 45°N, 100°W, azimuth 35°, jettison at burnout, impact at 0 km.
Run 2026-08-07:

```
boost           propellant exhausted at t =  80.518 s, m = 2400.0 kg (nominal)
coast           gamma reached zero at t = 838.208 s, apogee (nominal)
descent         reached 0.0 km at t = 1620.614 s (nominal)
liftoff mass       32400.0 kg   (payload 900 + dry 1500 + propellant 30000), T/W 2.990
burnout            t =  80.518 s, h = 102.000 km, V = 5768.20 m/s, gamma +47.600 deg
mass into coast      900.0 kg    (booster jettisoned)
apogee              1619.235 km  (exact event; 0.3569 of the range, so LOFTED)
flight time         1620.61 s
ground range        4536.36 km   (central angle 40.7508 deg)
impact speed        3156.49 m/s (Mach 9.96), gamma -47.725 deg
boost max q            84.13 kPa; boost sensed load 40.36 g at end of burn
re-entry max q       6281.46 kPa; re-entry decel 68.51 g at h = 1.66 km
Keplerian check     closed form 4467.5975 km vs propagated 4467.5975 km,
                    8.484e-13 relative (budget 1e-6, asserted)
drag effect           -0.0396 km (-0.0009 %); lift effect +0.1640 km (+0.0037 %)
impact speed, vacuum 5936.29 m/s vs flown 3156.49 m/s (-2779.80 m/s)
samples                 1053
```

The aerodynamic effect on *range* is tiny (0.003 %) because a steep ballistic
arc meets the atmosphere only on the terminal dive, where the force acts nearly
along the descent direction. It is not tiny in *speed*: drag removes 2780 m/s.

### What `run_boost_glide` produces

Shipped block — 20°N, 155°W, azimuth 60°, jettison at burnout, glide to a 15 km
handoff, 75° descent bank, impact at 0 km. Run 2026-08-07:

```
boost           propellant exhausted at t =  80.518 s (nominal)
glide           reached the 15.0 km handoff at t = 2120.526 s (nominal)
descent         reached 0.0 km at t = 2194.774 s (nominal)
burnout            t =  80.518 s, h = 46.998 km, V = 5840.75 m/s, gamma +9.758 deg
                   (a DEPRESSED burnout, which is what makes the glide possible)
glide              2040.01 s, 7527.52 km covered, peak altitude 168.28 km (first skip)
handoff            V = 694.89 m/s (Mach 2.19), gamma -4.213 deg
descent              74.25 s,   30.13 km covered
flight time         2194.77 s   (36.58 min)
ground range        7663.05 km  (central angle 68.8384 deg)
impact speed         401.57 m/s (Mach 1.27), gamma -32.228 deg
glide max q          260.19 kPa; peak aero load 8.33 g at h = 30.51 km
samples                 4469
```

The same handoff state flown on the *continued glide schedule* instead of the
descent bank takes 235.37 s and covers 96.44 km, against the descent's 74.25 s
and 30.13 km — the phase split is measured against that counterfactual, not
asserted.

**Two things about this run a reader will otherwise misread.**

**1. The reported aerodynamic peak is a LOAD FACTOR, not a deceleration.** The
8.33 g peak is `sqrt(aLift² + aDrag²)/g0` — what the structure feels. It breaks
down as 7.74 g of lift (92.8 % of the load) and 3.10 g of drag (37.1 %). **Only
the drag part brakes the vehicle.** An earlier version of this summary labelled
the 8.33 g figure a deceleration, which overstated the actual braking by 2.7×.
The per-phase mean drag decelerations — 0.266 g in the glide, 0.771 g in the
descent — are the deceleration numbers.

**2. `hHandoff` can bifurcate the trajectory.** The glide descends in a damped
skip phugoid, so a handoff altitude placed above a phugoid trough terminates the
glide a whole skip early. **Every phase still reports NOMINAL** — each event fired
exactly as asked — and the termination diagnosis cannot see it. Measured, all six
runs performed 2026-08-07:

| `hHandoff` | Ground range | Counterfactual rebound | Warning? |
|---|---|---|---|
| 15 km (shipped) | 7663.05 km | +0.00 km | no |
| 22 km | 7522.80 km | +0.60 km | **yes** |
| 25 km | 7490.28 km | +0.00 km | **no — blind spot** |
| 30 km | 5781.57 km | +31.33 km | **yes** |

At 30 km the run loses **1881 km of range, a quarter of the total**, and reports
nominal throughout. The script now detects this off a **counterfactual glide**:
it flies the handoff state on the *glide* schedule, and if the handoff caught a
trough the vehicle climbs back out of it, so a rebound above the handoff altitude
is the evidence. The propagation is needed anyway for the descent comparison, so
the check is free.

**The check is ONE-SIDED, and 25 km is the known blind spot.** It sees a trough
the continued glide climbs out of; it cannot see a handoff that truncated a
shallow skip without rebounding. The 25 km case loses 173 km against the shipped
run with no rebound and therefore no warning. Partial coverage was taken because
it costs nothing, not because it is complete.

### What `run_target` produces

Shipped block — launch 20°N 155°W, target 35°N 120°W, cutoff bracket 0.50 to
1.00 of full burn, 1 km range tolerance, zero bank throughout, Earth rotation
off. Re-run 2026-08-09, abridged from the printed summary:

```
boost            thrust terminated on command at t = 75.6820 s with 1801.8 kg
                 of propellant still aboard (nominal)
required range     3811.240 km     (great circle, launch to target)
achieved range     3811.751 km
range residual      +511.24 m      (achieved minus required, signed. NOT what the
                                    solve drives; see the two components below)
MISS DISTANCE        511.24 m      (impact to target, from the FLOWN state)
down-range miss     +511.24 m      (component along the intended course at the target)
cross-range miss      -0.00 m      (component across it. BOTH are what the solve drives)
tolerance            707.11 m      (on EACH component; 1000.00 m on the total miss)
seed azimuth      56.627204 deg    (the closed-form initial bearing, which is the GUESS)
launch azimuth    56.627204 deg    (SOLVED, +0.000000 deg from the seed)
seed miss            511.24 m      (down-range +511.24, cross-range -0.00)
flown azimuth     56.627204 deg    (-6.488e-12 deg from the commanded one)
solved cutoff       75.6820 s      (0.939941 of the 80.5178 s full burn)
seed cutoff         75.6820 s      (moved by +0.0000 s in stage 2)
propellant burned   28198.2 kg     of 30000.0; 1801.8 kg thrown away unburned
reachable           202.861 to 7737.630 km  (flown at the SEED azimuth)
iterations               10        (10 bisection steps seeding 0 Newton step(s))
propagations             14        (12 in the seed bisection, 1 in the Newton, plus
                                    one to re-fly the answer)
flight time         1626.09 s      (27.10 min)
impact            35.001317 N, 119.994629 W  (target 35.000000, -120.000000)
```

Three independent checks on that solve, all recomputed outside the script:

| Check | Reference | Measured |
|---|---|---|
| Impact-to-target distance, `greatCircle` on the flown terminal state | the reported 511.243460471 m miss | `511.243460471 m`, difference **`0.000e+00`** |
| Miss against the magnitude of the range residual — this is the zero-cross-range property, measured | equal | `9.286e-10 m` apart |
| `traj.x(1,6)` against the **solved** azimuth `info.psiLaunch` | equal | to 1e-13 rad, and the solved azimuth is asserted **different** from the seed by more than a degree on the rotating case, which is the one mutation a reinstated closed form would not survive |

The whole solve — fourteen propagations and the summary — takes about 2 s. A
propagation here is roughly 0.15 s, which is why the range tolerance can be
tightened almost for free: each halving of the bisection costs one more. Note
what the propagation line says about this case: **zero Newton steps.** The seed
miss is already inside the per-component tolerance, so the two-axis solve exits
at iteration zero having spent one evaluation, and the shipped run is
bit-for-bit what it was before `aimSolve` existed. That is a demonstration that
the second axis costs nothing when there is nothing to do — it is *not* evidence
that the solve works, and no tightening of the tolerance can make it so. The
rotating case (23 propagations, 231551.628 m → 4.361 m) and the banked one (18
propagations, 21524.695 m → 54.795 m) are where the evidence is.

### What `run_ballistic_target` produces

Shipped block — launch 45°N 100°W, target 62°N 28°W, loft bracket −140° to 85°,
`alphaMax` read from the vehicle (6°), 1 km range tolerance, 50 m constraint
tolerance, `branch = 'minimum-energy'`, zero bank throughout, Earth rotation
off. Re-run 2026-08-09, abridged from the printed summary:

```
required range      4828.045 km   (great circle, launch to target)
max-range loft      -42.9074 deg  (bracketed: 21 scan points + 13 golden steps)
CERTIFIED interval  -42.9290 to -42.8858 deg   ranges 5211.525 and 5211.525 km
MAXIMUM RANGE       5211.525 km   apogee 974.102 km, flight time 1314.777 s
reachable            556.603 to 5211.525 km   (flown at the SEED azimuth)

THE TWO ARCS THAT REACH THIS TARGET, both at the SEED azimuth
                          depressed           lofted
  loft angle (deg)        -126.7286          19.4336
  burnout gamma (deg)       21.8634          43.2151
  achieved range (km)     4828.5025        4828.8247
  miss (m)                 457.2699         779.4909
  APOGEE (km)              569.1767        1453.6613
  FLIGHT TIME (s)         1047.9491        1550.4998
  IMPACT SPEED (m/s)      1781.3246        3012.9766
  IMPACT ANGLE (deg)       -22.9679         -43.4028
  burnout energy (MJ/kg)   -44.5365         -44.8168
  loft from max (deg)       83.8212          62.3409
  bisection steps                 8               10
  propagations                   11               13

asked for   minimum-energy   flew minimum-energy (NOT ON A BRANCH)
                             MINIMISED (seed)     FLOWN (solved)
  loft angle (deg)                 -31.7656         held, same
  cutoff fraction (-)              0.994976           0.994976
  cutoff time (s)                   80.1132            80.1132
  BURNOUT ENERGY (MJ/kg)         -45.555275         -45.555275
  achieved range (km)             4828.0062          4828.0062
    THE TWO COLUMNS ARE THE SAME RUN. The seed azimuth was already inside the
    aim tolerance, so the two-axis solve converged at iteration 0.
  neighbouring feasible -45.486364 and -44.818885 MJ/kg; valley 68911.0 J/kg deep
                        on the near side and 736389.7 J/kg on the far side
  measured noise         1.010e+02 J/kg  = 1.47e-03 of the shallower side
  APOGEE                 965.613 km     classical  952.696 km   +1.36 %
  FLIGHT TIME            21.2817 min    classical 20.1686 min   +5.52 %
  burnout gamma           33.3290 deg   gamma* 34.1572 deg  (DIAGNOSTIC)
  propellant unburned     150.7 kg      thrown away with the booster
seed azimuth      40.555398 deg  (the closed-form initial bearing, which is the GUESS)
launch azimuth    40.555398 deg  (SOLVED, +0.000000 deg from the seed)
seed miss             39.01 m    (down-range -39.01, cross-range -0.00)
flown azimuth     40.555398 deg  (-1.134e-10 deg from the commanded one)
down-range miss      -39.01 m
cross-range miss      -0.00 m
MISS DISTANCE         39.01 m    (impact to target, from the FLOWN state)
impact            62.000071 N, 28.000731 W  (target 62.000000, -28.000000)
propagations             733     (37 bracketing the maximum, 1 re-flying it,
                                  11 depressed, 13 lofted, 669 minimum-energy,
                                  1 in the two-axis aim solve and 1 to re-fly
                                  the answer)
```

`info.nProp` used to return **57** where the summary printed 59, the difference
being the two propagations that re-create each branch's state history at its
converged loft angle. Both now read `mr.nEval + 1 + dep.nProp + lof.nProp`, and
`tests/test_runBallisticTarget.m` parses the printed count and compares it with
the returned one.

Four independent checks, all recomputed outside the script:

| Check | Reference | Measured |
|---|---|---|
| Impact-to-target distance, `greatCircle` on the flown terminal state | the reported 779.49087505 m lofted miss | `779.49087505 m`, difference **`0.000e+00`** |
| Same for the depressed arc | the reported 457.26985976 m | `457.26985976 m`, difference **`0.000e+00`** |
| Which branch was flown on a `branch = 'lofted'` run, re-derived from the flown loft angle against the certified maximiser interval −42.9290° to −42.8858° | `lofted` requested | +19.4831°, **above** `bL` — `lofted` |
| `traj.x(1,6)` against the **solved** azimuth `info.psiLaunch` | equal | to 1e-13 rad, with the solved azimuth asserted more than a degree off the seed on the rotating case |

Note the third row's loft angle. The **seed-azimuth** table above reports the
lofted arc at +19.4336°, because that is where stage one's bisection put it; the
flown arc is at **+19.4831°**, because stage two then moved it to close a
779.491 m seed miss to **0.365 m**. Those are two different numbers for two
different things and neither may stand in for the other.

Stage one costs **62 propagations** — both branches, the bracketing, the
certified interval — and stage two adds the aim solve and one re-fly on top:
**64** for a `'depressed'` run (3.3 s, no Newton step needed) and **67** for a
`'lofted'` one (one Newton iteration, four evaluations). A ballistic propagation
is roughly 0.05 s, three times cheaper than the boost-glide chain, which is what
makes a two-branch solve affordable at all. A `'minimum-energy'` run costs
**733** and about 28 s: every outer evaluation of the objective is a whole inner
feasibility solve.

### The branch is measured, and the measurement is now pinned

The script reports the branch it *measured* beside the branch it was *asked*
for, and that measurement had **no test behind it**: every case in the suite was
one where the two agree, so replacing the whole call with
`flownName = pickName; branchOK = true` left the suite green while
`info.branchMeasured`, `info.branchAgrees` and the summary's `MEASURED as` line
all quietly became restatements of the command.

Since 2026-08-08 the measurement is the flown loft angle's **position** against
the certified maximiser interval, which is the only branch invariant available.
Finding a case where that can *legitimately* differ from the label is
structural: each branch is solved on an interval lying entirely on one certified
side of the maximum, so a solution strictly inside its own bracket always
measures as the branch it was solved on. The exception is the **endpoint**: when
the required range falls within the range tolerance of the largest range the
search can certify, `coorbital.util.rangeSolve` short-circuits there, both branch
solves come back holding an arc INSIDE the unresolved interval, and the position
test can only answer **coalesced** — because at that range the two arcs are the
same arc.

Part 11 of `tests/test_runBallisticTarget.m` flies exactly that — 5212.791 km on
a 50 km tolerance against a 5211.525 km certified maximum — and asserts
`branchMeasured == 'coalesced'` against a `'lofted'` command, `branchAgrees`
false, `coalescedRq` true, and the printed caution. The caution diagnoses the
degeneracy ("landed INSIDE the certified maximiser interval … the bracketing is
NOT suspect") instead of blaming a bracketing step that did nothing wrong. With
the bypass in place, part 11 fails.

---

## Placeholder parameters

**Every vehicle and booster number in this library is a marked PLACEHOLDER
open-literature value.** They demonstrate that the machinery works and are of the
right order of magnitude. **No trajectory number in this document is a
performance prediction for any real vehicle.**

This table is the complete set — every field of all three parameter files.

| Struct | Values |
|---|---|
| `vehicleDefaults` / `vehicle_hgv` | `mass = 900 kg`, `Sref = 0.75 m²`, `CL = 0.35`, `L/D = 2.5`, `noseRadius = 0.05 m` |
| `vehicle_bm` | `mass = 900 kg`, `Sref = 0.385 m²`, `CL = 0.005`, `L/D = 0.02`, `noseRadius = 0.10 m` |
| `boosterDefaults`, mass and propulsion | `massDry = 1500 kg`, `massProp = 30000 kg`, `thrustVac = 950 kN`, `Isp = 260 s` (vacuum), `Aexit = 1.25 m²` |
| `boosterDefaults`, **stack aerodynamics** | `Sref = 1.77 m²`, `CL = 0.05`, `L/D = 0.25` |

**Do not read the booster aerodynamic triple as inert bookkeeping.** It sets the
boosted stack's drag through the whole powered ascent, so it moves the burnout
state every later phase inherits. And with `separation = false` it describes the
airframe for the *entire* unpowered flight as well: a slender body of revolution
at `L/D = 0.25` barely glides. Measured on `run_boost_glide` 2026-08-07 —
7663.05 km with separation against **2853.71 km without, a 62.8 % loss of
range**, driven by these three numbers.

`noseRadius` is the one field in the table no physics routine reads. It is
carried for the deferred Sutton–Graves heating rate and says so at its point of
definition; three tests assert only that it exists and passes through
`vehicle_hgv` unchanged.

Both chain scripts print `(PLACEHOLDER values)` on the vehicle and booster lines
of their own summaries, so the caveat travels with the output.

**What is actually pinned — corrected 2026-08-07.** An earlier version of this
section said `test_constThrust` "pins all of them exactly (lines 50–56)". It did
not: that span holds **seven assertions**, and they cover **six** of the eighteen
values in the table above — the five `boosterDefaults` mass and propulsion fields
plus the 900 kg payload mass. The seventh assertion pins `g0`, which is not a
parameter-file value at all; it comes from `missileConst`. Three gaps have since
been closed and the rest are stated honestly here:

| Values | Pinned where |
|---|---|
| `boosterDefaults` `massDry`, `massProp`, `thrustVac`, `Isp`, `Aexit`; payload `mass`; `g0` | `test_constThrust.m:51-57`, exact equality — seven pins |
| `boosterDefaults` `Sref`, `CL`, `LD`, plus the derived `CD = 0.20` | `test_constThrust.m`, own block immediately below — **added 2026-08-07.** Nothing in that file's hand arithmetic touches the triple, so every other check there stayed green at any positive values, while these three drive the 62.8 % range swing measured above |
| `vehicleDefaults` `Sref`, `CL`, `LD` | `test_runGlide.m`, exact equality — **added 2026-08-07.** `test_equilibriumGlide` builds its closed-form `Veq(r)` from these very fields, so they cancel between model and reference and it cannot see them. Same disease as the `Hscale` blindness below; same cure |
| `vehicle_bm` `mass`, `Sref`, `CL`, `LD` | `test_runBallistic.m`. `mass` is tied to the pinned 900 kg payload; the derived `CD = CL/LD = 0.25` is pinned at 1e-12; the ballistic coefficient `mass/(CD·Sref) = 9350.649350649351 kg/m²` at 1e-9; and **`CL` is pinned absolutely at 0.005**, which the other three do not cover — see below |
| `noseRadius`, both files | **Not pinned, deliberately** — no physics routine reads it. The tests assert only that it exists, is positive, and passes through `vehicle_hgv` unchanged |

Those pins are what stands between a mistyped constant and a suite that still
reports green.

**Why `vehicle_bm`'s `CL` needs a pin of its own, when `CD` and β are already
pinned.** Those two constrain only the *ratio*. Doubling `CL` to 0.010 and `LD`
to 0.040 together holds `CD = 0.25` and β = 9350.649350649351 kg/m² **exactly**,
so both joint pins stay green — and the trajectory still moves, because
`glide3DOF` reads `CL` **directly** for lift:

```matlab
             aLift = qbar.*veh.Sref.*CL./veh.mass;
             aDrag = qbar.*veh.Sref.*CD./veh.mass;
```

The split of a fixed `CD` between `CL` and `LD` is therefore **observable —
through lift, not through drag**. Measured on that mutation: impact speed moves
3156.49 → 3144.86 m/s and the pinned trajectory literal breaks at 3.683e-03
relative, with nothing in the failure saying why. An unexplained
trajectory-literal break is exactly what this section exists to prevent, so `CL`
is pinned at its own value and the cause fires first. `LD = CL/CD` is then fixed
by the two pins together.

An earlier version of this table asserted the opposite — that no trajectory
could see the split. That was wrong, and wrong in the instructive direction: the
reasoning stopped at the drag term.

Two artefacts of the placeholder motor are worth naming so they are not read as
physics. It is constant-thrust and does not throttle, so the sensed load climbs
to ~40 g at end of burn as the stack empties; a real stage would tail off or
stage before that, and both scripts say so. And the ballistic re-entry body's
`L/D = 0.02` gives a 68 g terminal deceleration at 1.66 km — plausible for a
high-β re-entry vehicle, but driven entirely by a placeholder.

---

## The validation suite's structural blind spot

The hardest-won lesson of the glide milestone, stated generally because it
transfers:

> **An analytic reference computed from the same constants as the model cannot
> detect a wrong constant.** The constant appears on both sides and cancels. No
> tolerance, however tight, recovers the sensitivity.

Proven by mutation, not argued:

| Mutation | Effect on `test_equilibriumGlide` + `test_allenEggers` | Contrast |
|---|---|---|
| `Hscale` × 1.15 | **both still PASS** — Allen–Eggers agreement moves only +10.50 → +10.44 % | |
| `rho0` × 1.3 | **both still PASS** | |
| `muE` × 1.02 | **both still PASS** | |
| drag acceleration × 0.8 in `glide3DOF` | equilibrium glide passes, Allen–Eggers **FAILS at +38.30 %** | a wrong *model* is caught; a wrong *constant* is not |

(All four re-measured 2026-08-06 against mutated copies of the library.)

Allen–Eggers *is* the exact solution for an exponential atmosphere, so `H`
cancels identically. What actually catches `rho0` and `muE` is not the analytic
suite at all but the unit tests that assert a hand-computed **absolute** value
which happens to depend on them (`test_expAtmos`'s 87909.92 Pa,
`test_sphereGrav`'s 9.7983 m/s²). `Hscale` had no such anchor anywhere, which is
precisely why it was the one constant left exposed.

**The only defence is pinning constants at their source.**
`test_missileConst` is therefore not a formality — it is the load-bearing test
for every physical constant in the library, and everything else is built on the
assumption that it holds. `Hscale` is now pinned to 7200 m within 1 m.

The working diagnostic, for anyone adding a constant: ask whether **any** test
asserts a number that would change if the constant were wrong. Relative checks —
one e-fold per scale height, energy conservation, a closed form built from the
same constants — do not count. And when a shared parameter turns out to be
invisible, pin it where it is *defined*; do not invent a physical-sounding
cross-check inside the test that noticed the blindness.

`LESSONS_LEARNED.md` has the full write-up, including the wrong turn taken first
and the non-constant form of the same trap (an analytic reference must read `CL`
from `veh.CL`, never back out of the aero model).

**The boost milestone found two more faces of the same problem**, both logged in
`LESSONS_LEARNED.md`:

- **A reduction test validates only the terms that survive the reduction.** The
  shared object that cancels is a whole *expression* rather than a constant. The
  `boost3DOF` → `glide3DOF` reduction is bit-exact and covers none of the thrust
  physics; four measured thrust mutations survived the whole suite before the
  force-increment check existed, tabulated in that entry.
- **Tsiolkovsky is blind to every booster constant.** `Isp`, `g0`, `thrustVac`
  and the mass ratio all cancel between propagation and closed form. **Seven**
  exact pins in `test_constThrust.m` are the entire defence — the `Hscale`
  lesson, applied to a second set of constants before it could bite. It was
  applied incompletely: the booster's *aerodynamic* triple went unpinned until
  2026-08-07, and it is not blind to Tsiolkovsky so much as absent from it. See
  [Placeholder parameters](#placeholder-parameters) for the corrected coverage.

The diagnostic is the same one in all four cases: **what would have to be wrong
for this number to move?**

---

## Conventions

MATLAB in **pumpkyn house style**, per `~/Desktop/proj7/doc/pumpkyn_style_guide.md`:

- `%%`-delimited header quartet — Purpose / Inputs / Outputs / (References) —
  with input and output columns aligned and every quantity carrying its size and
  units.
- `=` signs column-aligned, names right-justified against them.
- Colon-terminated `%%` step comments through the body.
- `if nargin == 0` self-demo in every library function, calling itself by full
  namespace.
- A `%% References:` block wherever the math has a source — Vinh for the EOM,
  Allen–Eggers for the entry solution, Sinnott for the haversine.

Beyond the style guide, four rules this library enforces:

- **SI throughout the library.** Metres, m/s, radians, seconds, kg. Human units
  (degrees, kilometres) appear *only* in an entry script's `%% USER PARAMETERS:`
  block and are converted immediately afterwards, in one place.
- **No `%#ok` pragmas**, of any flavour. Humans do not write them and they hide
  real problems. Unused output? Capture it with `~`. Growing array? Preallocate.
- **Never `i` or `j` as a loop variable** — they are the imaginary unit. Use `k`,
  `kb`, `kp`, `kt`, or a meaningful name.
- **Never `norm`** — use `vmag` with an explicit dimension. (Currently unused
  here; the rule stands for when vector work arrives.)

Re-verified 2026-08-09 across all **61** `.m` files under `missiles/` — 26
public library functions, 3 `+viz/private/` helpers, 5 entry scripts, 2 vehicle
parameter files and 25 under `tests/` (`run_tests.m` plus 24 `test_*.m`):

```bash
find /Users/msc/Desktop/optimal_control/missiles -name '*.m' | wc -l
```

Zero `%#ok`, zero `for i =` / `for j =`, and the single `norm(` hit is a comment
in `tests/test_viz.m` *explaining* why `norm()` is not used. A count refresh, not
a violation.

Three ruled exceptions to the self-demo rule: `missileConst`, `vehicleDefaults`
and `boosterDefaults` are pure constant returns whose entire content is visible
in the header, so a demo would print back what the reader is already looking at.
Every other library function has one.

Two standing cautions inherited from the design spec: **do not modify
`external/pumpkyn`**, and do not consume `getConst`'s known-bad `deg2ArcSec` or
`c` fields.

---

## Out of scope

Deliberately excluded. Each becomes its own plan.

| Excluded | Why / what it needs |
|---|---|
| **Closed-loop boost guidance** — PEG and VOA | Still prescribed pitch only, and **neither algorithm exists**: no `+guide` file for either, no test, no entry script that calls one. What landed on 2026-08-08 is the *seam*, not the laws — `coorbital.guide.terminalConstraint` evaluates the five terminal constraints both enforce, and a design spike measured the two questions gating them. The four gaps are named in `software_design.tex` §*The closed-loop guidance seam*; `docs/closed_loop_guidance.md` turns them into decisions with numbers attached; the mathematics is `hgv_dynamics_note.tex` §9. |
| **Multi-stage boosters** | One stage, one burn. `phaseRun`'s per-phase `link` is the mechanism a second stage would use; nothing else blocks it. |
| **Terminal guidance laws** | The descent phase flies a prescribed schedule, not a homing law. |
| **Aerothermal heating** | No Sutton–Graves anywhere. `noseRadius` remains a carried placeholder, flagged as such at its point of definition. |
| ~~**Rotating-Earth targeting**~~ | **Built, 2026-08-08.** Both targeting scripts fly `earthSpin` true instead of refusing it: 231551.628 m → 4.361 m on `run_target`, 463211.19 m → 52.461 m on `run_ballistic_target`. The *reason* this row gave was wrong as well as overtaken — there is no outer azimuth iteration, because the target does not move; the vehicle is deflected. See [Turning the Earth on](#turning-the-earth-on-what-it-costs-and-what-it-buys). |
| ~~**Cross-range steering with bank**~~ | **Superseded, 2026-08-08.** Cross-range is no longer measured and warned about — it is one of the two residual components `coorbital.util.aimSolve` drives to zero, and `info.crossWarn` is gone from both scripts. It is closed by the **launch azimuth**, not by bank: bank as a *control variable* is still not solved for, and remains out of scope. |
| ~~**Solving for loft angle rather than cutoff time**~~ | **Built.** `BM/run_ballistic_target` brackets the max-range hump and solves both full-burn branches — and, on `'minimum-energy'`, solves the loft angle and the *cutoff fraction* together for the textbook trajectory, which needs a burnout speed the full burn overshoots. See [Ballistic point-to-point targeting](#ballistic-point-to-point-targeting-two-branches). |
| ~~**`hgv_dynamics_note.tex`, `software_design.tex`**~~ | **Both written**, and both re-swept against two-axis targeting — the math note's targeting and guidance sections rewritten on 2026-08-08, `software_design.tex` now on its third edition, 2026-08-09. They were deliberately deferred until the interfaces had survived contact with working code. The `.tex` sources are tracked; the PDFs are gitignored, so a fresh clone builds them. |
| **Phase 2 optimization** against `orbit_transfer/verify_common` | Prescribed-control simulation first, deliberately: the throughput requirement (multiple trajectories per second) is what shapes the architecture. |
| Further fidelity increments — rotating Earth **on by default**, `j2Grav`, geodetic altitude, tabulated aero, US76 atmosphere | Rotation is now *supported* everywhere and still ships `false` in every user block, which is a default and not a limitation. The rest are each behind their own validation test. All are one-line swaps by design; see [Adding a fidelity level](#adding-a-fidelity-level). |

**Delivered since the first milestone**, and no longer out of scope: `boost3DOF`
with a prescribed pitch program, `constThrust`, the burnout and apogee events,
`BM/run_ballistic`, the descent phase, the full boost → glide → descent chain
in `HGV/run_boost_glide`, the **`+viz` package** — which was on this list for two
milestones and has now absorbed the plotting that was triplicated inline across
the entry scripts, and added `globeMovie` — **point-to-point targeting** in
`HGV/run_target` and `BM/run_ballistic_target`, **two-axis targeting** in
`coorbital.util.aimSolve`, which retired the rotating-Earth and cross-range rows
above, and — 2026-08-09 — **figure export**, `coorbital.viz.saveFigure` with
`BM/run_ballistic_target`'s `plotFile` on top of it, closing a gap nothing had
noticed: no entry script could write a picture to disk.
`util/stateConvert.m` was **not** built and will not be
— see [The seven-state convention](#the-seven-state-convention-and-the-inert-mass-trap)
and `DESIGN.md` §11.
