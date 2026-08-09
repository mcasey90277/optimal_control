# Missile Trajectory Library

*Coorbital, Inc.*

A MATLAB **3-DOF point-mass trajectory generator** for hypersonic glide vehicles
and ballistic missiles over a spherical, optionally **rotating** Earth. Give it
a vehicle, a launch or entry state and a prescribed control schedule and it
integrates boost, glide and descent on one state vector and reports the flight;
give either targeting script two lat/lon pairs and it solves the **launch
azimuth and a range control together** for the trajectory that connects them.
Every propagator is checked against a closed-form solution. There is no
optimization here — this is the validated forward model that a Phase 2
optimizer will call.

> **Every vehicle and booster number in this library is a marked PLACEHOLDER
> taken from open literature.** No trajectory number anywhere — in this file, in
> a printed summary, in a plot — is a performance prediction for any real
> vehicle. See [Placeholders](#placeholders).

---

## Quickstart

Everything below was executed against the code as committed, and every number
quoted in this file came out of one of those runs. The **targeting** numbers,
the file counts and the suite result were re-measured on 2026-08-09; the
`run_glide`, `run_ballistic` and `run_boost_glide` headlines are from
2026-08-07 and are pinned to their last printed digit by the suite, which is
green as of 2026-08-09; the movie figures are from 2026-08-07 and are marked
where they appear.

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

# launch point -> destination, solved (boost-glide)
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles/HGV'); run_target(struct('showPlots',false))"

# launch point -> destination, solved (ballistic; branch is a user parameter)
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles/BM'); run_ballistic_target(struct('showPlots',false))"

# the same, over a rotating Earth
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles/BM'); run_ballistic_target(struct('showPlots',false,'earthSpin',true))"
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

Run 2026-08-09: **24 passed, 0 failed**, zero warnings, 305 s wall clock
including MATLAB startup. `run(...)` is required, not stylistic — `-batch` parses
a bare `tests/run_tests` as the expression `tests / run_tests` and dies before
the harness loads.

The suite got slower as it got sharper: it was about 15 s at 20 tests on
2026-08-07, and most of the increase is `test_runBallisticTarget` and
`test_runTarget`, which fly whole targeting solves — hundreds of trajectory
propagations each — rather than checking a function against a closed form.

### Saved figures

`BM/run_ballistic_target` writes its three figures to disk through
`coorbital.viz.saveFigure`. One user-block entry controls it, `plotFile`, a path
**stem without an extension**, defaulting to
`fullfile(tempdir,'run_ballistic_target')`:

```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles/BM'); \
   run_ballistic_target(struct('plotFile','/Users/msc/Desktop/optimal_control/missiles/results/china_to_la'))"
```

Each figure appends its own suffix and `.png` — `<stem>_profile.png`,
`<stem>_ground_track.png`, `<stem>_globe.png` — so the names are **fixed** and a
re-run overwrites its own pictures rather than piling up new ones. Saving
happens whenever `showPlots` is true and the stem is non-empty; setting
`plotFile` to `''` draws the figures without writing them. There is deliberately
no second on/off flag, exactly as `movieFile` has none beside it. The paths are
printed at the end of the run and handed back on `info.plotFiles`.

`HGV/run_target` has the same gap and is **not** wired to `saveFigure`; that is
a known follow-up and `TODO.md` carries it.

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
(`optimal_control/.gitignore:48`). Movies, frames and saved figures are build
artefacts — never commit them. The sanitised-title stills sitting there —
`china_to_la_ground_track_10802_km_requi.png` and the rest of the
`china_to_la_*` and `pacific_*` set — were written by an ad-hoc scratch driver
in an earlier session, because until 2026-08-09 no entry script could write a
figure at all. `plotFile` is how one does it now.

---

## The five entry scripts

| Script | Flies | Phases | Headline, shipped user block |
|---|---|---|---|
| `HGV/run_glide` | unpowered glide from a 60 km / 6000 m/s entry state | 1 | **6986.82 km**, 2073.77 s, peak 1.11 g sensed aero load |
| `BM/run_ballistic` | ballistic missile, pad to impact | 3 (boost, coast to apogee, descent) | **4536.36 km**, 1620.61 s, apogee 1619.23 km |
| `HGV/run_boost_glide` | boost-glide vehicle, pad to impact | 3 (boost, glide, terminal descent) | **7663.05 km**, 2194.77 s, peak 8.33 g aero **load factor** |
| `HGV/run_target` | **launch point → destination**, solved | 3 (boost to a solved cutoff, glide, descent) | 20°N 155°W → 35°N 120°W: required **3811.240 km**, miss **511.243 m**, solved cutoff 75.6820 s, 14 propagations |
| `BM/run_ballistic_target` | **launch point → destination**, solved, ballistic | 3 (boost to burnout or a solved cutoff, coast to apogee, descent) | 45°N 100°W → 62°N 28°W: required **4828.045 km**, miss **39.009 m** on the shipped `'minimum-energy'` branch, apogee 965.613 km, 733 propagations |

The first three fly *launch site + azimuth + control schedule → wherever the
physics puts it*. The last two invert that, and since 2026-08-08 both do it the
same way: the closed-form great-circle bearing is the **seed**, not the answer,
and `coorbital.util.aimSolve` runs a damped Newton on the launch azimuth *and*
a range control until both components of the miss — down-range and cross-range,
resolved at the target — are inside tolerance. The range control is the
thrust-termination time for `run_target`, and for `run_ballistic_target` it is
the loft angle on the two full-burn branches or the cutoff fraction in
`'minimum-energy'` mode.

**What two-axis targeting bought, measured 2026-08-09.** Both scripts now fly a
rotating Earth and a banked track instead of refusing the first and warning
about the second:

| Case | Miss with the seed bearing alone | Miss solved | Residual evaluations | Propagations, whole run |
|---|---|---|---|---|
| `run_target`, `earthSpin` true | 231 551.628 m | **4.361 m** | 7 | 23 |
| `run_target`, 75° terminal bank | 21 524.695 m | **54.795 m** | 4 | — |
| `run_target`, shipped | 511.243 m | 511.243 m, unchanged bit-for-bit | 1 | 14 |
| `run_ballistic_target`, `earthSpin` true | 463 211.19 m | **52.461 m** | 7 | 701 |
| `run_ballistic_target`, shipped `'lofted'` | 779.491 m | **0.365 m** | 4 | 67 |
| `run_ballistic_target`, shipped `'depressed'` | 457.270 m | 457.270 m, unchanged bit-for-bit | 1 | 64 |
| `run_ballistic_target`, shipped `'minimum-energy'` | 39.009 m | 39.009 m, unchanged bit-for-bit | 1 | 733 |

Read the cost column, because it decides whether turning rotation on is
affordable. A residual evaluation is one full trajectory propagation. The
solve costs `1 + 3n` evaluations for `n` Newton iterations, so the hard cases
above are two iterations and the 4-evaluation ones are one. The *rest* of each
run is the one-dimensional stage that produces the seed, and that stage is
where `run_ballistic_target`'s cost lives: 62 propagations for a full-burn
branch, but **693 to 731** for `'minimum-energy'`, because the loft angle is
not monotonic, the max-range hump has to be bracketed and certified before
either branch can be bisected, and minimum-energy then minimises the burnout
specific energy along the feasible family on top of that.

Two things a user turning rotation on should know. `run_ballistic_target`'s
reachable envelope is measured at the **seed** azimuth, so the too-far and
too-close refusals are slightly approximate — and rotation *moves* that
envelope: maximum range 5211.5 → 5439.9 km and the depressed-branch floor
4708.5 → 5085.8 km on the shipped geometry, which is why the shipped target
plus `earthSpin` is legitimately refused for `'depressed'`. That is correct
physics on an easterly launch, not a bug. And `run_target`'s shipped
non-rotating zero-bank case exits the two-axis solve at iteration zero, so it
demonstrates that the second axis costs nothing when there is nothing to do,
and nothing else. `TODO.md` carries both.

The 8.33 g on `run_boost_glide` is a **load factor**, what the structure feels
(7.74 g lift, 3.10 g drag) — not a deceleration. The script says so itself.

The four chain scripts return `[traj,info]`, with every printed number at full
precision in `info`. A refused run returns an empty `traj`, `info.refused =
true` and an `info.refusedWhy` naming which gate stopped it.

---

## Layout

Everything reusable lives in one `+coorbital` package; the vehicle folders hold
only parameter files and entry scripts, no physics.

| Path | Contents |
|---|---|
| `+coorbital/+util/` | constants, vehicle and booster parameter defaults, great-circle range and bearing, and the two black-box solvers — `rangeSolve`, the scalar range bisection, and **`aimSolve`**, the two-axis damped Newton that solves an azimuth beside a range control |
| `+coorbital/+atmos/` | `expAtmos` — exponential isothermal atmosphere |
| `+coorbital/+grav/` | `sphereGrav` — spherical gravity, with the J2 latitudinal channel already plumbed and returning zero |
| `+coorbital/+aero/` | `constLD` — constant `CL` and `L/D`; `CD` is derived, never stored |
| `+coorbital/+eom/` | `glide3DOF` (6-state), `boost3DOF` (7-state powered), `massConstant` (lifts a 6-state EOM to 7) |
| `+coorbital/+guide/` | `prescribed` schedule interpolation, `pitchProgram` boost pitch attitude, `terminalConstraint` the five burnout constraints PEG and VOA share (see `docs/closed_loop_guidance.md`) |
| `+coorbital/+prop/` | the multi-phase driver `phaseRun`, the `constThrust` motor, and the altitude / apogee / burnout ODE events |
| `+coorbital/+viz/` | `groundTrack`, `profilePlot`, `globe3D`, `globeMovie`, and `saveFigure` — the one export helper the whole package shares (+ 3 `private/` helpers) |
| `HGV/`, `BM/` | 5 entry scripts and 2 vehicle parameter files |
| `tests/` | `run_tests.m` plus 24 `test_*.m` |
| `docs/` | design spec, code-organization README, lessons learned, four plan briefs, the closed-loop guidance spike, two LaTeX notes (PDFs built locally, gitignored), archived reviews |
| `results/` | rendered movies, frames and saved figures — **gitignored** |

26 public library functions across the eight packages, 3 private helpers,
5 entry scripts, 2 vehicle files, 25 files under `tests/` (`run_tests.m` plus
24 `test_*.m`): **61 `.m` files in all**. Counted 2026-08-09 with

```bash
find /Users/msc/Desktop/optimal_control/missiles -name '*.m' -not -path '*/results/*' | wc -l
```

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
and `:333`, `tests/test_fullChain.m:238`, `BM/run_ballistic.m:545`.

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
| `docs/DESIGN.md` | The pre-implementation design spec, deliberately left unedited, plus its six dated as-built sections (§10 the glide propagator, §11 boost and the full chain, §12 targeting and visualization, §13 two-axis targeting, §14 the closed-loop guidance spike, §15 figure export) recording where the spec and the code diverged |
| `docs/LESSONS_LEARNED.md` | The running log — what broke, what fixed it, and what would otherwise be rediscovered the hard way |
| `docs/hgv_dynamics_note.tex` | **The mathematics** — frames, the 3-DOF equations term by term, the rotating-Earth projections, the seven closed-form references, the two-axis targeting solve, and a chapter separating targeting from the two specified-but-unimplemented closed-loop boost laws |
| `docs/software_design.tex` | **The structure** — package boundaries, data flow from the user block to a rendered movie, the two solvers, the entry-script shape, the four structural blindnesses, and the closed-loop guidance seam |
| `docs/plan_2026-08-06_glide_propagator.md`, `docs/plan_2026-08-07_boost_descent_chain.md`, `docs/plan_2026-08-07_targeting_and_viz.md`, `docs/plan_2026-08-08_aim_solve_and_closed_loop.md` | The four milestone plan briefs, with the design decisions and out-of-scope lists behind each |
| `docs/reviews/` | Archived external reviews of the plan briefs |
| `TODO.md` | Open items — what is not built, what is unreviewed, and what is out of scope by design |

Both LaTeX notes exist. The `.tex` sources are tracked; the PDFs are **not** —
`*.pdf` is gitignored at `optimal_control/.gitignore:8`, so a fresh clone has to
build them. Two passes, so the cross-references resolve, then clean the aux
files:

```bash
cd /Users/msc/Desktop/optimal_control/missiles/docs
/Library/TeX/texbin/pdflatex hgv_dynamics_note.tex && /Library/TeX/texbin/pdflatex hgv_dynamics_note.tex
rm -f *.aux *.log *.out *.toc *.lof *.lot *.fls *.fdb_latexmk *.synctex.gz
```
