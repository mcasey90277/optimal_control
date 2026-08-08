# Missile Trajectory Library — Code Organization

*Coorbital, Inc.*

MATLAB 3-DOF trajectory generation for hypersonic glide and ballistic vehicles
over a spherical Earth. The library propagates **prescribed-control**
trajectories: give it a vehicle, a launch or entry state, and a control
schedule, and it integrates the flight and reports it. No optimization — see
[Out of scope](#out-of-scope).

Three milestones have shipped:

| Milestone | Delivers |
|---|---|
| 2026-08-06 | Validated glide propagator — `glide3DOF`, `phaseRun`, `HGV/run_glide` |
| 2026-08-07 | Powered boost, a descent phase, and full three-phase chains — `boost3DOF`, `constThrust`, `pitchProgram`, `massConstant`, `BM/run_ballistic` (boost → coast → impact), `HGV/run_boost_glide` (boost → glide → descent) |
| 2026-08-07 | **Point-to-point targeting and visualization** — `greatCircleBearing`, `rangeSolve`, the `+viz` package (`groundTrack`, `profilePlot`, `globe3D`, `globeMovie`), and `HGV/run_target`: give it a launch point and a destination and it solves the trajectory that connects them. See [Point-to-point targeting](#point-to-point-targeting) |
| 2026-08-07 | **Ballistic point-to-point targeting** — `BM/run_ballistic_target`, which ranges on the **loft angle** rather than on thrust termination and therefore has to deal with the two branches that come with it: a lofted arc and a depressed one for every range short of maximum. See [Ballistic point-to-point targeting](#ballistic-point-to-point-targeting-two-branches) |

The design rationale is in `DESIGN.md`; the running log of what broke and why is
in `LESSONS_LEARNED.md`. This file covers layout, use, and extension.

**Before you read a single trajectory number as a performance figure: every
vehicle and booster parameter in this library is a marked PLACEHOLDER.** See
[Placeholder parameters](#placeholder-parameters).

---

## Running

Everything below was executed on 2026-08-07 against the code as committed.
Every command shown was run; every number shown came out of that run.

### The five entry scripts

| Script | Flight | Phases |
|---|---|---|
| `HGV/run_glide` | unpowered glide from an entry state | 1 (glide) |
| `BM/run_ballistic` | ballistic missile, pad to impact | 3 (boost, coast to apogee, descent) |
| `HGV/run_boost_glide` | boost-glide vehicle, pad to impact | 3 (boost, glide, terminal descent) |
| `HGV/run_target` | **launch point → destination point**, solved | 3 (boost to a solved cutoff, glide, descent) |
| `BM/run_ballistic_target` | **launch point → destination point**, solved, ballistic | 3 (boost to burnout, coast to apogee, descent) |

The first three fly *launch site + azimuth + control schedule → wherever the
physics puts it*. The last two invert that: they take two latitude/longitude
pairs and solve for the trajectory that connects them. They are separate scripts
because they solve on **different controls** — `run_target` on thrust
termination, `run_ballistic_target` on the loft angle — and the ballistic one
therefore has two answers where the boost-glide one has one.

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

Headless, which is how all four were verified:

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

### The glide entry script in detail

At the MATLAB prompt:

```matlab
cd ~/Desktop/optimal_control/missiles/HGV
run_glide
```

Called with an output — `traj = run_glide;` — it returns the trajectory struct;
called bare it prints the summary and returns nothing. All three scripts behave
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
  PASS  test_allenEggers      PASS  test_greatCircle
  PASS  test_boost3DOF        PASS  test_greatCircleBearing
  PASS  test_boostEvents      PASS  test_missileConst
  PASS  test_constLD          PASS  test_phaseRun
  PASS  test_constThrust      PASS  test_pitchProgram
  PASS  test_equilibriumGlide PASS  test_rangeSolve
  PASS  test_expAtmos         PASS  test_runBallistic
  PASS  test_fullChain        PASS  test_runGlide
  PASS  test_glide3DOF        PASS  test_runBallisticTarget
                              PASS  test_runGlide
                              PASS  test_runTarget
                              PASS  test_sphereGrav
                              PASS  test_viz

21 passed, 0 failed
```

42 s wall clock including MATLAB startup, measured 2026-08-07 on an Apple
silicon Mac, and zero warnings. Fifteen of the twenty-one are unit or
analytic-validation tests — thirteen exercise one library function apiece, and
`test_allenEggers` and `test_equilibriumGlide` check the propagator against a
closed-form solution. The other six test **compositions**, and they are the
important ones:

| Test | What it composes |
|---|---|
| `test_runGlide` | `HGV/run_glide` end to end — unit conversion, the `greatCircle` call site, the peak search, the termination diagnosis, and the printed summary, which it parses rather than recomputes. Also guards the deliberate duplication between `vehicle_hgv` and `vehicleDefaults`. |
| `test_runBallistic` | `BM/run_ballistic` end to end — the boost → coast → descent chain, the staging link, and the Keplerian cross-check. |
| `test_runTarget` | `HGV/run_target` end to end — the closed-form azimuth, the bisection on cutoff time, the separation link, the reachable-envelope refusal at BOTH ends of the band, and the printed summary. The miss it asserts is measured impact-to-target with `greatCircle` from the flown state, not read back out of the solver's own residual — and it flies a **banked** case as well as the shipped zero-bank one, because at zero bank those two numbers agree to 9.3e-10 m and nothing can tell a measurement from a substitution. |
| `test_runBallisticTarget` | `BM/run_ballistic_target` end to end — the closed-form azimuth, the **bracketing of the max-range loft angle**, the two branch solves either side of it, the branch selector, and three separate refusals (beyond maximum range, too close, and an unreachable *branch*). What it checks that no other test in this suite can: the branch is **re-derived from the flown apogee and flight time** rather than read off the script's own label, the lofted and depressed arcs are asserted to differ by more than five times in apogee and two in flight time, and `minimum-energy` is flown at **two** targets because it takes the lofted arc at 3175 km and the depressed one at 4218 km — one case alone cannot tell a working rule from a hard-wired answer. It also pins the `alphaMax` guard: at `BM/run_ballistic`'s 6° clamp the depressed branch does not exist and the script must refuse with `coorbital:runBallisticTarget:maximumNotBracketed`. |
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
| `+coorbital/+util/` | `missileConst.m`, `vehicleDefaults.m`, `boosterDefaults.m`, `greatCircle.m`, `greatCircleBearing.m`, `rangeSolve.m` |
| `+coorbital/+atmos/` | `expAtmos.m` |
| `+coorbital/+grav/` | `sphereGrav.m` |
| `+coorbital/+aero/` | `constLD.m` |
| `+coorbital/+eom/` | `glide3DOF.m`, `boost3DOF.m`, `massConstant.m` |
| `+coorbital/+guide/` | `prescribed.m`, `pitchProgram.m` |
| `+coorbital/+prop/` | `phaseRun.m`, `constThrust.m`, `eventAltitude.m`, `eventApogee.m`, `eventBurnout.m` |
| `+coorbital/+viz/` | `groundTrack.m`, `profilePlot.m`, `globe3D.m`, `globeMovie.m`, plus `private/` helpers |
| `HGV/` | `run_glide.m`, `run_boost_glide.m`, `run_target.m`, `vehicle_hgv.m` |
| `BM/` | `run_ballistic.m`, `run_ballistic_target.m`, `vehicle_bm.m` |
| `tests/` | `run_tests.m` plus twenty `test_*.m` |
| `docs/` | this file, `DESIGN.md`, `LESSONS_LEARNED.md`, the three plans, reviews |

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
with no sign flip. Bearing is **not** symmetric under exchange of the two
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

**`+viz`** is the plotting package: `groundTrack` (lat/lon track, one coloured
segment per phase, launch, target and impact marked), `profilePlot` (a chosen
subset of altitude, speed, Mach, dynamic pressure, load factor, mass and flight
path against time), `globe3D` (a still 3-D Earth with the trajectory arc) and
`globeMovie` (the same scene as an MP4 with the trajectory developing over
time). Every one takes the same `(traj,veh,env,opts)` argument list — several
of them do not read `veh` or `env` and say so in their headers — and every one
reads the trajectory and never writes it, so no figure can move a number in a
summary. The Earth texture and starfield come from the pumpkyn toolbox when it
is on the path and degrade to a plain shaded sphere on a black background when
it is not; both paths are a fully working figure and need no network.

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

Two more were added this milestone. **Any state dimension:** `nx` is taken
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
`glide3DOF` and `boost3DOF`, or leave it at 0. In all three entry scripts the
switch is the boolean `earthSpin` in the user block, shipped `false`.

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
trajectory that connects them. That is **two separate solves**, and on a
non-rotating Earth they do not interact.

**Azimuth — closed form, no iteration.** `coorbital.util.greatCircleBearing`
gives the initial bearing of the launch-to-target great circle, clockwise from
north, which is exactly the heading state `psi`. It is written straight into the
launch state; `test_runTarget` asserts `traj.x(1,6)` equals it to machine
precision, and the shipped run reproduces it bit-identically.

**Range — bisection on thrust-termination time.** Of the available control
parameters, only this one is monotonic and single-valued:

| Parameter | Behaviour | Verdict |
|---|---|---|
| Thrust-termination time | less burn, less energy, shorter range | **chosen** |
| Pitch-program loft angle | two branches either side of a max-range hump | rejected — a root finder lands on whichever branch it started nearest |
| Glide handoff altitude | bifurcates on the glide phugoid's troughs; 30 km costs `run_boost_glide` 1882 km while every phase still reports nominal | rejected outright |

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
throw, so a caller can read the band out of `info` without a `try`/`catch`.

### Two limitations, printed in every summary

**The azimuth is exact only for a non-rotating Earth.** With `env.omegaE = 0`
the ground track of a zero-bank trajectory is a great circle, so the initial
bearing is the whole answer. Turn rotation on and it is not: over the shipped
1626 s flight the ground beneath a 35°N target sweeps **620 km** east, so the
vehicle would have to be aimed where the target is going to be. That makes the
azimuth depend on the flight time, the flight time on the cutoff, and the
cutoff on the azimuth — an **outer iteration** around the range solve, which
`run_target` does not have and which is out of scope for it. Setting
`earthSpin` true prints a caution saying the reported miss is against a target
that stood still.

**The miss is the residual of the range solve, along the great circle.**
Bisection matches a *distance*; nothing in it steers sideways. Cross-range comes
out at zero in the shipped configuration because every commanded bank angle is
zero — and that is a property of *that configuration*, not a general guarantee.
`run_target` **measures** the cross-track offset of the impact point from the
launch-to-target great circle every run and warns when it exceeds the range
tolerance. Give the descent `run_boost_glide`'s 75° terminal bank and it does:
the heading turns 166.42° between handoff and impact, the range solve still
converges to a 585.20 m residual, and the vehicle misses by **21524.70 m** — a
factor of 36.8 over the residual. The summary then prints a `*** WARNING ***`
and switches its limitations paragraph to the banked wording, because the
zero-bank explanation printed under a bank angle contradicts the warning
directly above it. `tests/test_runTarget.m` flies that exact case and pins all
of it; without it the measure-and-warn behaviour is unpinned, since at zero bank
the measured miss and the solver's own residual agree to 9.3e-10 m and no
assertion can tell them apart. That is why `run_target` ships
`descBank = 0` where `run_boost_glide` ships 75 — a targeting decision, not an
aerodynamic one, and the user block says so at the point of definition.

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

---

## Ballistic point-to-point targeting: two branches

`BM/run_ballistic_target` takes a launch point and a destination and solves the
**ballistic** trajectory that connects them. The azimuth solve is the same
closed-form `greatCircleBearing` call `HGV/run_target` makes. The range solve is
not, and the difference is the whole reason this is a separate script.

**The ranging control is the loft angle, and range is not monotonic in it.**
`run_target` bisects on thrust-termination time because less burn means less
energy means less range — one answer, and bisection is safe. A ballistic missile
is flown to a **loft angle** instead: the terminal attitude of the pitch
program. Range rises to a maximum at some max-range loft angle and falls away on
both sides, so **every range short of that maximum is reached by two
trajectories**:

| | depressed arc | max-range arc | lofted arc |
|---|---|---|---|
| loft angle | −14.36° | **25.07°** | 59.15° |
| range | 3174.31 km | 5055.30 km | 3174.87 km |
| apogee | 208.36 km | 1021.49 km | 2118.45 km |
| flight time | 659.01 s | 1327.81 s | 1816.01 s |
| impact speed | 924.20 m/s | — | 3502.79 m/s |
| impact angle | −18.23° | — | −62.74° |

Measured on the shipped geometry: 45°N 100°W to 62°N 60°W, 3174.98 km required
on a 39.198° azimuth. **Both arcs are always solved and always reported**, and
only one is flown, because the trade is the point.

So the script **brackets the maximum first** — a 13-point coarse scan across the
loft bracket, then golden-section refinement to 0.05° — which does two things:
it splits the loft axis into two intervals on each of which range *is* monotonic
(the precondition `coorbital.util.rangeSolve` documents and does not check), and
it produces the max-range arc's own apogee and flight time. Those two numbers
are the reference the flown branch is then **measured against**, because the
max-range arc sits between the two branches in both. A bracket lying on one side
of the maximum *should* keep the solution on that side; the script does not take
`should` for evidence, and prints the measured branch beside the requested one.

**The steep arc arrives faster here, which is the opposite of the vacuum
intuition.** With the propellant load fixed both arcs leave burnout with almost
the same energy — −45.15 and −44.91 MJ/kg, 0.55 % apart — so in vacuum they
would arrive at almost the same speed. With an atmosphere the shallow −18.2°
depressed descent spends far longer in dense air than the steep −62.7° lofted
plunge and is braked to 26 % of its arrival speed. The script measures this and
prints it; it does not assert the textbook expectation.

### What `minimum-energy` can and cannot mean, and why it says so at length

The textbook minimum-energy ballistic trajectory reaches a required range on the
**least burnout energy**, at a burnout flight-path angle of
γ\* = 45° − Λ/4. That presumes the energy is free to choose — you size the
booster to the range. This one is not: a fixed propellant load, burned to
exhaustion on both branches. **No member of the family is the minimum-energy
trajectory** unless the target sits at maximum range.

What is implemented is stated in the summary every run: the branch whose loft
angle lies **nearer the max-range loft angle**, on the grounds that the
max-range arc *is* the minimum-energy arc for its own range. Two other readings
are printed beside it and used by neither:

- **Least burnout energy read literally** would pick the *depressed* arc on the
  shipped case. That reading is rejected on purpose: with the propellant fixed, a
  lower burnout energy means a costlier ascent, not a cheaper trajectory.
- **Nearest γ\*** agrees with the implemented rule here, but it is measured on a
  burnout state the `alphaMax` clamp and not the script controls, so it is a
  diagnostic.

**And the rule is not always well conditioned.** The range-versus-loft hump is
roughly symmetric, so the two distances can fall nearly equal. Measured along the
shipped azimuth: 14 % apart at 3175 km (lofted wins), 4 % apart at 3789 km, and
by 4218 km the ranking has **flipped** and the depressed arc wins by 10 %. The
summary prints both distances and **cautions** inside a 5 % band, advising the
only sound course there — ask for `'lofted'` or `'depressed'` explicitly.
`tests/test_runBallisticTarget.m` flies all three of those targets.

### The clamp that decides whether the second branch exists at all

`run_ballistic_target` ships `alphaMax = 12°` where `BM/run_ballistic` ships 6°,
and that is a targeting decision, not an aerodynamic one. The loft angle is a
*commanded attitude*; what the vehicle achieves at burnout is limited by how fast
the angle-of-attack clamp lets the flight path be pushed over. Measured by
scanning the shipped bracket at 6, 8, 10 and 12°:

| `alphaMax` | burnout γ at the most depressed command | max-range angle inside the bracket? |
|---|---|---|
| 6° | 33.4° at a commanded −30° | **no** — range is monotone decreasing across the whole bracket |
| 8° | 23.6° at a commanded −30° | yes, near 5° |
| 10° | 15.0° at a commanded −30° | yes, near 25° |
| 12° | 7.6° at a commanded −30° | yes, at 25.07° |

At 6° the clamp holds burnout γ above the max-range value at *every* commanded
loft angle, **the depressed branch does not exist**, and the two-branch structure
this script is built around collapses. The bracketing step refuses with
`coorbital:runBallisticTarget:maximumNotBracketed` and names that cause first,
rather than reporting a one-branch problem as a two-branch one with one branch
quietly missing.

### Three refusals, none of which throws

| Case | What the summary says |
|---|---|
| Beyond maximum range | `BEYOND MAXIMUM RANGE by 7268.302 km`, with the 5055.302 km maximum, the 25.0684° loft angle achieving it, the 556.603–5055.302 km band, and each branch's own band and loft interval |
| Too close | `TOO CLOSE`, with the overflight in kilometres and what to change |
| **Branch unreachable** | the target is inside the envelope but only the *other* arc reaches it — `The depressed arc does not reach this target, though the other one does.` |

The third has no counterpart in `run_target`. Silently flying the branch that
does reach the target would be the worst outcome available: it converges, it
hits, and it is not the trajectory that was asked for. All three return an empty
trajectory with `info.refused = true`, and none throws, so a caller reads the
band out of `info` without a `try`/`catch`.

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
**Seven** exact-equality pins at `tests/test_constThrust.m:50-56` are the only
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
1.00 of full burn, 1 km range tolerance, zero bank throughout. Run 2026-08-07:

```
boost            thrust terminated on command at t = 75.6820 s with 1801.8 kg
                 of propellant still aboard (nominal)
required range      3811.240 km    (great circle, launch to target)
achieved range      3811.751 km
range residual        +511.24 m    (achieved minus required, signed)
MISS DISTANCE          511.24 m    (impact to target, from the FLOWN state)
launch azimuth       56.627204 deg (closed form, no iteration)
flown azimuth        56.627204 deg (-6.488e-12 deg from the commanded one)
cross-track              -0.00 m   (-4.064e-07 m)
solved cutoff         75.6820 s    (0.939941 of the 80.5178 s full burn)
propellant burned     28198.2 kg   of 30000.0 kg; 1801.8 kg thrown away unburned
reachable            202.861 to 7737.630 km
iterations                  10     (12 trajectory propagations, plus one to re-fly)
flight time          1626.09 s     (27.10 min)
impact               35.001317 N, 119.994629 W  (target 35.000000, -120.000000)
```

Three independent checks on that solve, all recomputed outside the script:

| Check | Reference | Measured |
|---|---|---|
| Impact-to-target distance, `greatCircle` on the flown terminal state | the reported 511.243460471 m miss | `511.243460471 m`, difference **`0.000e+00`** |
| Miss against the magnitude of the range residual — this is the zero-cross-range property, measured | equal | `9.286e-10 m` apart |
| `traj.x(1,6)` against `greatCircleBearing(launch,target)` | equal | **bit-identical**, `0` ulp |

The whole solve — twelve propagations, the re-fly, and the summary — takes about
2 s. A propagation here is roughly 0.15 s, which is why the range tolerance can
be tightened almost for free: each halving costs one more.

### What `run_ballistic_target` produces

Shipped block — launch 45°N 100°W, target 62°N 60°W, loft bracket −40° to 85°,
`alphaMax` 12°, 1 km range tolerance, `branch = 'minimum-energy'`, zero bank
throughout. Run 2026-08-07:

```
required range      3174.981 km   (great circle, launch to target)
max-range loft       25.0684 deg  (bracketed: 13 scan points + 13 golden steps)
MAXIMUM RANGE       5055.302 km   apogee 1021.486 km, flight time 1327.812 s
reachable            556.603 to 5055.302 km

                          depressed           lofted
  loft angle (deg)         -14.3602          59.1457
  burnout gamma (deg)       11.8153          62.8405
  achieved range (km)     3174.3149        3174.8651
  miss (m)                 666.0454         115.8585
  APOGEE (km)              208.3574        2118.4489
  FLIGHT TIME (s)          659.0085        1816.0100
  IMPACT SPEED (m/s)       924.2005        3502.7891
  IMPACT ANGLE (deg)       -18.2269         -62.7418
  burnout energy (MJ/kg)   -45.1546         -44.9073
  loft from max-range       39.4287          34.0773

asked for   minimum-energy   flew lofted   MEASURED as lofted
launch azimuth    39.197731 deg  (closed form, no iteration)
flown azimuth     39.197731 deg  (-2.200e-10 deg from the commanded one)
cross-track           -1.17e-05 m
MISS DISTANCE        115.86 m    (impact to target, from the FLOWN state)
impact            61.999681 N, 60.002110 W  (target 62.000000, -60.000000)
propagations              57     (29 bracketing the maximum, 1 re-flying it,
                                  14 depressed, 15 lofted)
```

Four independent checks, all recomputed outside the script:

| Check | Reference | Measured |
|---|---|---|
| Impact-to-target distance, `greatCircle` on the flown terminal state | the reported 115.85854573 m miss | `115.85854573 m`, difference **`0.000e+00`** |
| Same for the depressed arc | the reported 666.04542234 m | `666.04542234 m`, difference **`0.000e+00`** |
| Which branch was flown, re-derived from `max(traj.x(:,1))` and `traj.t(end)` against the max-range arc's 1021.486 km and 1327.812 s | `lofted` requested | apogee 2118.449 km and 1816.010 s, both **above** — `lofted` |
| `traj.x(1,6)` against `greatCircleBearing(launch,target)` | equal | **bit-identical**, `0` ulp |

The whole run — 57 propagations, both branches, the bracketing and the summary —
takes about 3 s. A ballistic propagation is roughly 0.05 s, three times cheaper
than the boost-glide chain, which is what makes a two-branch solve affordable at
all.

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
| `boosterDefaults` `massDry`, `massProp`, `thrustVac`, `Isp`, `Aexit`; payload `mass`; `g0` | `test_constThrust.m:50-56`, exact equality — seven pins |
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

The hardest-won lesson of this milestone, stated generally because it transfers:

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

Re-verified 2026-08-07 across all 39 `.m` files under `missiles/`: zero `%#ok`,
zero `norm(`, zero `for i =` / `for j =`.

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
| **Closed-loop boost guidance** — PEG and VOA | Prescribed pitch only. `hyperFLIGHT` uses both, and they are the natural next milestone. `pitchProgram` already reads the state, so the signature is ready. |
| **Multi-stage boosters** | One stage, one burn. `phaseRun`'s per-phase `link` is the mechanism a second stage would use; nothing else blocks it. |
| **Terminal guidance laws** | The descent phase flies a prescribed schedule, not a homing law. |
| **Aerothermal heating** | No Sutton–Graves anywhere. `noseRadius` remains a carried placeholder, flagged as such at its point of definition. |
| **Rotating-Earth targeting** | Both targeting scripts' closed-form azimuth is exact only at `omegaE = 0`. Rotation needs an outer azimuth iteration around the range solve — and in `run_ballistic_target` it would need one *per branch*, the two arcs flying for 659 s and 1816 s over the same geometry. See [Point-to-point targeting](#point-to-point-targeting). |
| **Cross-range steering with bank** | The range solve controls downrange only. `run_target` measures the cross-track offset and warns, but cannot close it. |
| ~~**Solving for loft angle rather than cutoff time**~~ | **Built.** `BM/run_ballistic_target` brackets the max-range hump and solves both branches. What remains out of scope is *choosing* between them on textbook minimum-energy grounds, which a fixed propellant load burned to exhaustion does not permit — see [Ballistic point-to-point targeting](#ballistic-point-to-point-targeting-two-branches). |
| **`hgv_dynamics_note.tex`, `software_design.tex`** | Written *after* the interfaces survive contact with working code. Two milestones in, they now have. |
| **Phase 2 optimization** against `orbit_transfer/verify_common` | Prescribed-control simulation first, deliberately: the throughput requirement (multiple trajectories per second) is what shapes the architecture. |
| Further fidelity increments — rotating Earth on by default, `j2Grav`, geodetic altitude, tabulated aero, US76 atmosphere | Each behind its own validation test. All are one-line swaps by design; see [Adding a fidelity level](#adding-a-fidelity-level). |

**Delivered since the first milestone**, and no longer out of scope: `boost3DOF`
with a prescribed pitch program, `constThrust`, the burnout and apogee events,
`BM/run_ballistic`, the descent phase, the full boost → glide → descent chain
in `HGV/run_boost_glide`, the **`+viz` package** — which was on this list for two
milestones and has now absorbed the plotting that was triplicated inline across
the entry scripts, and added `globeMovie` — and **point-to-point targeting** in
`HGV/run_target`. `util/stateConvert.m` was **not** built and will not be
— see [The seven-state convention](#the-seven-state-convention-and-the-inert-mass-trap)
and `DESIGN.md` §11.
