# Missile Trajectory Library — Code Organization

*Coorbital, Inc.*

MATLAB 3-DOF trajectory generation for hypersonic glide vehicles over a
spherical Earth. This milestone delivers a **validated prescribed-control glide
propagator**: give it a vehicle, an entry state, and a control schedule, and it
integrates the trajectory and reports it. No optimization, no boost, no
terminal phase — see [Out of scope](#out-of-scope).

The design rationale is in `DESIGN.md`; the running log of what broke and why is
in `LESSONS_LEARNED.md`. This file covers layout, use, and extension.

---

## Running

Everything below was executed on 2026-08-06 against the code as committed.

### The entry script

At the MATLAB prompt:

```matlab
cd ~/Desktop/optimal_control/missiles/HGV
run_glide
```

`run_glide` puts itself and the library root on the path, so it works from
anywhere once MATLAB can see the file. Everything a routine run needs is in the
fenced `%% USER PARAMETERS:` block at the top; nothing below that block should
need editing. Called with an output — `traj = run_glide;` — it returns the
trajectory struct; called bare it prints the summary and returns nothing.

It also takes one optional argument, a struct of overrides for named
USER PARAMETERS entries, in the block's own human units:

```matlab
traj = run_glide(struct('psiEntry',135,'latEntry',25,'showPlots',false));
```

Every entry in the block is overridable and nothing else is — a misspelt name
raises rather than silently leaving the shipped value in place. This exists so
an automated test or a batch sweep can drive the script at more than one
operating point without editing it; `tests/test_runGlide.m` is the reason it
was added, and the throughput requirement in `DESIGN.md` §1 is where it goes
next. For a routine interactive run, ignore it and edit the block.

Headless, which is how it was verified:

```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles/HGV'); run_glide"
```

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
  PASS  test_allenEggers      PASS  test_missileConst
  PASS  test_constLD          PASS  test_phaseRun
  PASS  test_equilibriumGlide PASS  test_runGlide
  PASS  test_expAtmos         PASS  test_sphereGrav
  PASS  test_glide3DOF
  PASS  test_greatCircle

10 passed, 0 failed
```

Nine of the ten test one library function apiece. `test_runGlide` is the odd
one out and the important one: it tests the **composition** — `HGV/run_glide.m`
end to end, including the unit conversion, the great-circle call site, the peak
search, the termination diagnosis, and the printed summary, which it parses
rather than recomputes. It also guards the deliberate duplication between
`vehicle_hgv` and `vehicleDefaults`. It runs two propagations: the shipped
due-east configuration, and a south-east one, because the shipped geometry is
provably blind to a lat/lon transposition at the `greatCircle` call site (with
the entry point at the origin the central angle is symmetric in the terminal
latitude and longitude, so the swap changes the summary by nothing at all).

### Self-demos

Every library function runs a demonstration when called with no arguments —
`coorbital.atmos.expAtmos` plots the density profile, `coorbital.aero.constLD`
prints the coefficients, `coorbital.eom.glide3DOF` prints one state derivative.
Add the library root to the path first (`addpath ~/Desktop/optimal_control/missiles`).
This is the fastest way to see what a routine does without reading it.

---

## State and sign conventions

Read this before touching anything. These are the conventions a contributor is
most likely to get wrong.

```
x = [r, lon, lat, V, gamma, psi]      u = [alpha, sigma]
```

| Symbol | Units | Meaning |
|---|---|---|
| `r` | m | **Geocentric radius**, not altitude. Altitude is `r - c.rE`. |
| `lon` | rad | Longitude, east positive |
| `lat` | rad | Geocentric latitude, north positive. **Singular at the poles** (`glide3DOF` errors rather than returning NaN). |
| `V` | m/s | Planet-relative speed. **Singular below 1 m/s** (errors). |
| `gamma` | rad | Flight path angle, **positive UP**. A descending entry has `gamma < 0`. Singular at ±90° (errors). |
| `psi` | rad | Heading, **clockwise from north**. 90° is due east. |
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
| `+coorbital/+util/` | `missileConst.m`, `vehicleDefaults.m`, `greatCircle.m` |
| `+coorbital/+atmos/` | `expAtmos.m` |
| `+coorbital/+grav/` | `sphereGrav.m` |
| `+coorbital/+aero/` | `constLD.m` |
| `+coorbital/+eom/` | `glide3DOF.m` |
| `+coorbital/+guide/` | `prescribed.m` |
| `+coorbital/+prop/` | `phaseRun.m`, `eventAltitude.m` |
| `HGV/` | `run_glide.m`, `vehicle_hgv.m` |
| `tests/` | `run_tests.m` plus ten `test_*.m` |
| `docs/` | this file, `DESIGN.md`, `LESSONS_LEARNED.md`, the plan, reviews |

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

**`atmos` / `grav` / `aero`** are the three swappable model families. See the
next section — this split is the point of the whole design.

**`eom/glide3DOF`** is the 6-state 3-DOF point-mass EOM (Vinh Eqs. 2.28–2.33).
It never names a model; it calls handles out of `env`. Rotating-Earth Coriolis
and centrifugal terms are written in from the start and gated on `env.omegaE`,
so setting that to zero recovers the non-rotating case *exactly* rather than
approximately. The three coordinate singularities (pole, vertical flight, zero
speed) raise identified errors rather than silently producing NaN.

**`guide/prescribed`** evaluates a control schedule `struct('tGrid',...,
'alpha',...,'sigma',...)`: linear inside the grid, clamped outside. Its
signature is `u = f(t,x,sched)` so a closed-loop guidance law drops in later
with the state already available.

**`prop/phaseRun`** is the driver. It integrates a `[1 x P]` array of phase
structs (`eom`, `guide`, `terminate`, `tspan`) with `ode45`, hands the terminal
state of each phase to the next, and records the junction states so a Phase 2
optimizer can enforce them as linkage constraints. Default tolerances are
`RelTol = AbsTol = 1e-10`, overridable through `env.odeRelTol` / `env.odeAbsTol`.
Two behaviours worth knowing: each phase is referenced to its *own* `tspan(1)`,
so a phase given `tspan = [10 50]` contributes 40 s and opens no gap; and the
sample at a phase boundary is recorded once, carrying the *outgoing* phase's
control, so a control discontinuity does not appear in `traj.u`.

**`prop/eventAltitude`** is a terminal, one-sided (`direction = -1`) ODE event,
so a lofted arc climbing back through the stop altitude does not end the run.

**`HGV/run_glide`** wires all of the above together, propagates, and prints a
summary. It diagnoses *why* the propagation stopped, so a run truncated by the
time horizon can never be misread as a completed glide, and it flags when the
terminal Mach has dropped below the hypersonic regime in which constant `CL` and
`L/D` is defensible.

---

## Adding a fidelity level

**This is the architectural point of the library.** The equations of motion
never name a model. They call function handles carried in an environment struct:

```matlab
env.atmos  = @coorbital.atmos.expAtmos;
env.grav   = @coorbital.grav.sphereGrav;
env.aero   = @coorbital.aero.constLD;
env.omegaE = 0;
```

Raising fidelity is therefore a **one-line change in an entry script**, not an
edit to `glide3DOF`. Each family has a fixed signature, exactly as implemented:

| Family | Signature | Notes |
|---|---|---|
| atmosphere | `[rho,P,T,a] = f(h)` | `h` geometric altitude above the reference sphere (m); `rho` kg/m³, `P` Pa, `T` K, `a` m/s. Elementwise over `[N x 1]`. |
| gravity | `[gr,gLat] = f(r,lat)` | `r` geocentric radius (m), `lat` geocentric latitude (rad); `gr` positive-inward magnitude (m/s²), `gLat` signed northward (m/s²). Both `[N x 1]`. |
| aerodynamics | `[CL,CD] = f(alpha,mach,veh)` | `alpha` rad, `mach` dimensionless, `veh` the vehicle struct. |

Earth rotation is not a handle — it is the scalar `env.omegaE`. Set it to
`c.omegaE` to turn on the Coriolis and centrifugal terms already present in
`glide3DOF`, or leave it at 0. In `run_glide` the switch is the boolean
`earthSpin` in the user block.

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
and the plotting are untouched. The same recipe with `[rho,P,T,a] = us76Atmos(h)`
adds a US76 atmosphere; the same recipe with `[CL,CD] = tableAero(alpha,mach,veh)`
adds a tabulated drag polar. Note that `alpha` currently has **no effect** on the
trajectory, because `constLD` ignores it by construction; it becomes live the
moment an alpha-dependent aero model is dropped in, and `run_glide` says so in
the user block.

---

## Validated results

Re-run 2026-08-06 against the committed code. All numbers below were produced by
running the code, not copied from a report.

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

**Caveat, and it matters.** Every vehicle parameter is a **PLACEHOLDER**
open-literature value (`mass = 900 kg`, `Sref = 0.75 m²`, `CL = 0.35`,
`L/D = 2.5`). These numbers demonstrate that the machinery works and are of the
right order; they are **not a performance prediction for any real vehicle**.
Note also that the run ends at Mach 1.01, well outside the hypersonic regime in
which holding `CL` and `L/D` constant is defensible — `run_glide` prints that
warning itself.

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

Verified: no `%#ok`, no `norm(`, and no `for i =` / `for j =` anywhere under
`missiles/`.

Two standing cautions inherited from the design spec: **do not modify
`external/pumpkyn`**, and do not consume `getConst`'s known-bad `deg2ArcSec` or
`c` fields.

---

## Out of scope

Deliberately excluded from this milestone. Each becomes its own plan.

| Excluded | Why / what it needs |
|---|---|
| **Boost phase** — `boost3DOF`, prescribed pitch program, `BM/run_ballistic.m` | Needs a Cartesian state with mass, hence `util/stateConvert.m`. `phaseRun` currently passes the state straight through because every Phase 1 phase shares the same 6-state glide vector. |
| **Terminal / descent phase** and the full boost → glide → descent chain | Waits on boost. |
| **`+viz` package** — ground tracks, 3-D trajectory over `pumpkyn.util.earth3D` | `run_glide` has inline plots for now. |
| **`hgv_dynamics_note.tex`, `software_design.tex`** | Written *after* the interfaces survive contact with working code, which is what this milestone establishes. |
| **Phase 2 optimization** against `orbit_transfer/verify_common` | Prescribed-control simulation first, deliberately: the throughput requirement (multiple trajectories per second) is what shapes the architecture. |
| Further fidelity increments — rotating-Earth validation, `j2Grav`, geodetic altitude, tabulated aero, US76 atmosphere | Each behind its own validation test. The plumbing is already in place; see [Adding a fidelity level](#adding-a-fidelity-level). |
