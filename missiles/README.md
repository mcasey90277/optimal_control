# Missile Trajectory Library

*Coorbital, Inc.*

A MATLAB **3-DOF point-mass trajectory generator** for hypersonic glide vehicles
and ballistic missiles over a spherical Earth. Give it a vehicle, a launch or
entry state and a prescribed control schedule and it integrates boost, glide and
descent on one state vector and reports the flight; give `run_target` two
lat/lon pairs and it solves for the trajectory that connects them. Every
propagator is checked against a closed-form solution. There is no optimization
here — this is the validated forward model that a Phase 2 optimizer will call.

> **Every vehicle and booster number in this library is a marked PLACEHOLDER
> taken from open literature.** No trajectory number anywhere — in this file, in
> a printed summary, in a plot — is a performance prediction for any real
> vehicle. See [Placeholders](#placeholders).

---

## Quickstart

Everything below was executed on 2026-08-07 against the code as committed, and
every number quoted in this file came out of one of those runs.

```bash
# glide only, from an entry state
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles/HGV'); run_glide"

# ballistic missile, pad to impact
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles/BM'); run_ballistic(struct('showPlots',false))"

# boost-glide vehicle, pad to impact
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles/HGV'); run_boost_glide(struct('showPlots',false))"

# launch point -> destination, solved
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles/HGV'); run_target(struct('showPlots',false))"
```

At the MATLAB prompt, `cd` to `HGV/` or `BM/` and type the script name. Each
script puts itself and the library root on the path. Each is driven from a
fenced `%% USER PARAMETERS:` block at the top; nothing below that block should
need editing. Each also accepts a struct of overrides for named block entries in
the block's own human units, and a misspelt name raises rather than silently
leaving the shipped value in place.

### The test suite

```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles'); run('tests/run_tests')"
```

Run 2026-08-07: **20 passed, 0 failed**, zero warnings, about 15 s wall clock
including MATLAB startup. `run(...)` is required, not stylistic — `-batch` parses
a bare `tests/run_tests` as the expression `tests / run_tests` and dies before
the harness loads.

### The movie

`run_target` can render the trajectory developing over a 3-D Earth. It is off by
default because it is the expensive part of a run:

```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles/HGV'); \
   run_target(struct('showPlots',false,'movieOn',true,'movieFrames',60, \
                     'movieFile','/Users/msc/Desktop/optimal_control/missiles/results/run_target.mp4'))"
```

Run 2026-08-07: 60 frames at 1280×720 and 20 fps, 2.4 MB, 10.6 s wall clock
including MATLAB startup. With the pumpkyn toolbox on the path the globe carries
a photographic Earth texture and a starfield; without it the run reports
`texture: 'plain', background: 'black'` and is still a fully working movie that
needs no network. `movieFile` defaults into `tempdir`.

**`missiles/results/` is where rendered output goes, and it is gitignored**
(`optimal_control/.gitignore:48`). Movies and frames are build artefacts — never
commit them.

---

## The four entry scripts

| Script | Flies | Phases | Headline, shipped user block |
|---|---|---|---|
| `HGV/run_glide` | unpowered glide from a 60 km / 6000 m/s entry state | 1 | **6986.82 km**, 2073.77 s, peak 1.11 g sensed aero load |
| `BM/run_ballistic` | ballistic missile, pad to impact | 3 (boost, coast to apogee, descent) | **4536.36 km**, 1620.61 s, apogee 1619.23 km |
| `HGV/run_boost_glide` | boost-glide vehicle, pad to impact | 3 (boost, glide, terminal descent) | **7663.05 km**, 2194.77 s, peak 8.33 g aero **load factor** |
| `HGV/run_target` | **launch point → destination**, solved | 3 (boost to a solved cutoff, glide, descent) | 20°N 155°W → 35°N 120°W: required **3811.240 km**, miss **511 m**, solved cutoff 75.6820 s |

The first three fly *launch site + azimuth + control schedule → wherever the
physics puts it*. `run_target` inverts that: closed-form great-circle bearing for
the azimuth, bisection on thrust-termination time for the range.

The 8.33 g on `run_boost_glide` is a **load factor**, what the structure feels
(7.74 g lift, 3.10 g drag) — not a deceleration. The script says so itself.

The three chain scripts return `[traj,info]`, with every printed number at full
precision in `info`.

---

## Layout

Everything reusable lives in one `+coorbital` package; the vehicle folders hold
only parameter files and entry scripts, no physics.

| Path | Contents |
|---|---|
| `+coorbital/+util/` | constants, vehicle and booster parameter defaults, great-circle range and bearing, the black-box range bisection |
| `+coorbital/+atmos/` | `expAtmos` — exponential isothermal atmosphere |
| `+coorbital/+grav/` | `sphereGrav` — spherical gravity, with the J2 latitudinal channel already plumbed and returning zero |
| `+coorbital/+aero/` | `constLD` — constant `CL` and `L/D`; `CD` is derived, never stored |
| `+coorbital/+eom/` | `glide3DOF` (6-state), `boost3DOF` (7-state powered), `massConstant` (lifts a 6-state EOM to 7) |
| `+coorbital/+guide/` | `prescribed` schedule interpolation, `pitchProgram` boost pitch attitude |
| `+coorbital/+prop/` | the multi-phase driver `phaseRun`, the `constThrust` motor, and the altitude / apogee / burnout ODE events |
| `+coorbital/+viz/` | `groundTrack`, `profilePlot`, `globe3D`, `globeMovie` (+ 3 `private/` helpers) |
| `HGV/`, `BM/` | 4 entry scripts and 2 vehicle parameter files |
| `tests/` | `run_tests.m` plus 20 `test_*.m` |
| `docs/` | design spec, code-organization README, lessons learned, three plan briefs, archived reviews |
| `results/` | rendered movies and frames — **gitignored** |

23 public library functions across the eight packages, 3 private helpers,
4 entry scripts, 2 vehicle files, 20 test files: 53 `.m` files in all.

The equations of motion never name a model — they call handles carried in an
`env` struct, so raising fidelity is a one-line change in an entry script rather
than an edit to `glide3DOF`. That is the architectural point of the library;
`docs/README.md` has the recipe.

---

## Validation

Reproduced 2026-08-07 by re-running each check. Budgets are the assertions
actually in the test files.

| Check | Budget | Measured |
|---|---|---|
| Vacuum specific energy `V²/2 − μ/r`, 300 s lofted arc | `< 1e-8` relative | **`1.281e-13`** |
| Equilibrium glide `V(r)`, 45 → 20 km, against the closed form | `< 3 %` | **1.51 %** worst over the whole arc |
| Allen–Eggers peak deceleration, two ballistic coefficients 8× apart | `0 < rel < 12 %`, one-sided below | **+10.50 %** and **+9.20 %** |
| `boost3DOF` → `glide3DOF` with the engine dead, 4 cardinal states | `< 1e-12` | **`0.000e+00`** — bit-exact |
| Thrust force increment, both `alpha` and `sigma` nonzero | `< 1e-12` | **`3.642e-16`** worst |
| Tsiolkovsky, full vacuum burn, thrust along velocity | `< 1e-8` relative | **`2.377e-10`** |
| Keplerian free-flight range, `run_ballistic`, asserted in the script | `< 1e-6` relative | **`8.484e-13`** |
| Cross-junction continuity, `ode89` @ 1e-12 vs the driver's `ode45` @ 1e-10 | `< 1e-3 m` on radius | **`1.340e-06 m`** (boost→glide), **`3.120e-07 m`** (glide→descent) |

Where those budgets live, if you want to read the assertion rather than trust
the table: `tests/test_glide3DOF.m:43`, `tests/test_equilibriumGlide.m:184`,
`tests/test_allenEggers.m:172` and `:177`, `tests/test_boost3DOF.m:180`, `:254`–`:275`
and `:333`, `tests/test_fullChain.m:238`, `BM/run_ballistic.m:496`.

**The two analytic entry checks are complementary, not redundant.** `Veq²(r)`
has no drag term, so equilibrium glide is blind to `CD`; a zero-lift
Allen–Eggers entry says nothing about `CL`. Allen–Eggers runs at 12 % rather
than a tight tolerance because the closed form drops gravity and freezes the
flight path angle, and it is one-sided **below** — both neglected effects can
only raise the peak, so a simulated peak *under* the analytic value means too
much drag, not a better approximation.

**And the structural caveat that shapes the whole suite:** an analytic reference
computed from the same constants as the model cannot detect a wrong constant —
it cancels. Physical constants are therefore pinned at their source in
`test_missileConst`, and the placeholder parameters are pinned in
`test_constThrust`, `test_runGlide` and `test_runBallistic`. `docs/README.md`
has the mutation table that proves it.

---

## Conventions

```
x = [r, lon, lat, V, gamma, psi]      u = [alpha, sigma]     unpowered, 6-state
x = [r, lon, lat, V, gamma, psi, m]   u = [alpha, sigma]     powered,   7-state
```

- `r` is **geocentric radius**, not altitude. Altitude is `r - c.rE`.
- `lon`, `lat` geocentric, east and north positive. Singular at the poles.
- `gamma` is the flight path angle, **positive up** — a descending entry has
  `gamma < 0`.
- `psi` is the heading, **clockwise from north**; 90° is due east.
- `m` is the **total mass currently carried**, payload plus whatever booster is
  still attached. Present only on powered chains, and every phase in one
  `phaseRun` call shares the state dimension, so the coast and glide phases
  carry it too.
- `sphereGrav` returns `gr` as a **positive inward magnitude**, not a signed
  radial component. The EOM are written to expect that.
- **SI throughout the library** — metres, m/s, radians, seconds, kg. Human units
  (degrees, kilometres) appear *only* in an entry script's `%% USER PARAMETERS:`
  block and are converted immediately afterwards in one place.

MATLAB in pumpkyn house style. Never `i` or `j` as a loop variable; no `%#ok`
pragmas.

**One trap worth knowing before you build a chain:** `glide3DOF` divides by
`veh.mass`, not by the state mass `x(7)`, so a staged chain must bind a
per-phase vehicle inside each EOM closure. `coorbital.eom.massConstant` guards
this and raises `coorbital:massConstant:massMismatch` rather than flying at the
wrong weight. The full write-up is in `docs/README.md`.

---

## Placeholders

Repeating it because it is the one thing that must not be lost: **every vehicle
and booster number is a marked PLACEHOLDER open-literature value.** They are the
right order of magnitude and they demonstrate that the machinery works. That is
all they are.

Both chain scripts print `(PLACEHOLDER values)` on the vehicle and booster lines
of their own summaries, so the caveat travels with the output. The complete
table of all eighteen values, and what each of them is actually pinned by, is in
`docs/README.md` under *Placeholder parameters*.

---

## Where to read next

| File | For |
|---|---|
| `docs/README.md` | **The detailed authority on delivered behaviour** — code organization, every interface, the seven-state convention, how to add a fidelity level, the full validated-results tables, out-of-scope list |
| `docs/DESIGN.md` | The pre-implementation design spec, deliberately left unedited, plus its two dated as-built sections (§10 the glide propagator, §11 boost and the full chain) recording where the spec and the code diverged |
| `docs/LESSONS_LEARNED.md` | The running log — what broke, what fixed it, and what would otherwise be rediscovered the hard way |
| `docs/plan_2026-08-06_glide_propagator.md`, `docs/plan_2026-08-07_boost_descent_chain.md`, `docs/plan_2026-08-07_targeting_and_viz.md` | The three milestone plan briefs, with the design decisions and out-of-scope lists behind each |
| `docs/reviews/` | Archived external reviews of the plan briefs |
| `TODO.md` | Open items — what is not built, what is unreviewed, and what is out of scope by design |

Two LaTeX notes are **forthcoming** and do not exist yet:
`docs/hgv_dynamics_note.tex` for the mathematics and `docs/software_design.tex`
for the software design. They were deferred by design until the interfaces had
survived contact with working code; they now have, and they are being written.
