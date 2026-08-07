# Point-to-Point Targeting and Visualization — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax.

**Goal:** An example script that takes a launch point and a destination point, solves for the trajectory that connects them, flies the full boost–glide–descent chain, plots the useful figures, and optionally renders a pumpkyn-style movie of the Earth with the trajectory developing over time.

**Architecture:** Two new library pieces make targeting possible — a great-circle initial bearing for the launch azimuth, and a one-dimensional range solver. A new `+viz` package, promised in the design spec since the first milestone and never built, absorbs the plotting currently triplicated inline across three entry scripts and adds the globe movie. The entry script is thin: it converts user inputs, calls the solver, calls the propagator, and calls `+viz`.

**Tech Stack:** MATLAB R2025b, the existing `missiles/+coorbital` library, and pumpkyn's `earth3D`, `sphere3D`, `stars3D`, `addFigureLogo` from `~/Desktop/proj7/external/pumpkyn/src/+pumpkyn/+util/`.

---

## Global Constraints

- pumpkyn house style: `%%` header quartet (`Purpose`/`Inputs`/`Outputs`/`Revision History`) with aligned three-column I/O blocks stating UNITS, closed by `%% ------------------------ Begin Code Sequence ---------------------------`; `=` signs column-aligned at column 20; colon-terminated `%%` step comments; `if nargin == 0` self-demo in every library function, calling itself by full namespace.
- `%% References:` whenever the math has a source.
- Never `i`/`j` as loop or index variables. Never `norm`.
- **No `%#ok` pragmas of any kind.** Unused outputs take `~`; growing arrays get preallocated.
- **No hard-coded physical constants** outside `coorbital.util.missileConst()`. Entry-script user blocks may hold user-chosen values.
- SI units in the library, radians for angles. Entry scripts may take degrees and kilometres in the user block and must convert immediately in one marked place.
- Author `Michael Casey`, date 08/07/2026, `Copyright 2026 Coorbital, Inc.`
- **Do not modify `+coorbital/+eom/glide3DOF.m` or `boost3DOF.m`.** Both are symbolically verified.
- **Never modify `~/Desktop/proj7/external/pumpkyn`** — it is a third-party dependency, read-only.
- Suite: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/missiles'); run('tests/run_tests')" 2>&1 | grep -vE "Home License|personal use|academic, research|organizational use"`
  A bare `tests/run_tests` does NOT work — `-batch` parses it as division. Allow 500 s.
- Baseline `16 passed, 0 failed`, zero warnings. The three existing entry scripts must keep their headline numbers unchanged: `run_glide` 6986.82 km / 2073.77 s / 1.11 g; `run_ballistic` 4536.36 km / 1620.61 s; `run_boost_glide` 7663.05 km / 2194.77 s / 8.33 g aero load.

---

## The targeting problem, and the design decisions taken

The library currently flies *launch site + azimuth + pitch program → wherever it lands*. The user wants *launch point + destination → trajectory*. That is two separate solves.

**Azimuth — closed form, no iteration.** The great-circle initial bearing from point 1 to point 2, measured clockwise from north, which is exactly the library's `psi` convention:

```
     psi0 = atan2( sin(dlon)*cos(lat2),
                   cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(dlon) )
```

Verified before this plan was written: LAX (33.94 N, 118.41 W) to JFK (40.64 N, 73.78 W) gives **65.867 deg** and **3978.8 km**, against the published great-circle values of about 66 deg and 3970 km.

**This is exact only because `env.omegaE = 0`.** With a non-rotating Earth the ground track of a zero-bank trajectory is a great circle, so the initial bearing is the whole answer. Under rotation the target moves during a 20-plus-minute flight and the azimuth would need an outer iteration. Rotation is out of scope here and the script must say so plainly in its summary, not bury it.

**Range — one-dimensional root solve on thrust-termination time.** Of the available control parameters:

| Parameter | Behaviour | Verdict |
|---|---|---|
| Thrust-termination time | Monotonic: less burn, less energy, shorter range | **Chosen** |
| Pitch-program loft angle | Two branches, lofted and depressed, either side of a max-range hump | Rejected as the default — a root finder lands on whichever branch it starts nearest, which is a bad default for an example script |
| Glide handoff altitude | Bifurcates on phugoid troughs (measured: 30 km loses 1881 km) | Rejected outright |

Thrust termination is also what real systems do for energy management, and it needs **no new library machinery**: the boost phase already ends at `tspan(2)` when the burnout event does not fire first, so a shortened `tspan` is a cutoff.

**Early cutoff leaves unburned propellant.** At separation, jettison the entire booster — dry structure and remaining propellant together — so the post-separation vehicle is exactly the payload:

```matlab
    phases(1).link = @(x) [x(1:6); veh.mass];
```

This satisfies the `massConstant` guard, which requires the carried mass to match the vehicle the next phase binds.

**The reachable envelope must be reported, not assumed.** Before solving, evaluate range at the minimum and maximum cutoff times to bracket. If the required range falls outside, fail with a message stating the achievable band in kilometres and what to change. A targeting script that silently returns the nearest miss is worse than one that refuses.

---

## File Structure

| File | Responsibility |
|---|---|
| `+coorbital/+util/greatCircleBearing.m` | Initial bearing, clockwise from north |
| `+coorbital/+util/rangeSolve.m` | Bisection on cutoff time to match a required range |
| `+coorbital/+viz/groundTrack.m` | Ground track on a lat/lon grid, with launch, target and impact marked |
| `+coorbital/+viz/profilePlot.m` | Altitude, speed, Mach, dynamic pressure, load factor vs time |
| `+coorbital/+viz/globe3D.m` | Static 3-D Earth with the trajectory arc |
| `+coorbital/+viz/globeMovie.m` | Animated globe, trajectory developing over time |
| `HGV/run_target.m` | The example entry script |
| `tests/test_greatCircleBearing.m`, `test_rangeSolve.m`, `test_viz.m`, `test_runTarget.m` | Tests |
| `docs/README.md` | Updated |

---

## Task 1: Great-circle bearing and the range solver

**Files:** Create `+coorbital/+util/greatCircleBearing.m`, `+coorbital/+util/rangeSolve.m`, `tests/test_greatCircleBearing.m`, `tests/test_rangeSolve.m`.

**Interfaces produced:**
- `psi0 = coorbital.util.greatCircleBearing(lat1,lon1,lat2,lon2)` — initial bearing (rad), clockwise from north, wrapped to `[0,2*pi)`. Elementwise, matching `greatCircle`'s broadcasting.
- `[tCut,rngAch,info] = coorbital.util.rangeSolve(rngReq,fRange,tLo,tHi,tolM)` — bisection on a scalar parameter. `fRange` is a handle taking the parameter and returning achieved surface range (m). Returns the parameter value, the achieved range, and an `info` struct with the bracket values, iteration count, and a `converged` flag.

- [ ] **Step 1: Write the failing tests**

`test_greatCircleBearing.m`, against hand-computed literals, all independent of the implementation:
- Due east on the equator, (0,0) to (0,30 deg): exactly `pi/2`.
- Due north, (0,0) to (10 deg, 0): exactly `0`.
- Due south, (10 deg,0) to (0,0): exactly `pi`.
- Due west, (0,0) to (0,-30 deg): exactly `3*pi/2` after wrapping.
- LAX to JFK: `65.867 deg` to within `0.01 deg`. Cite the source of that number in a comment.
- A point to itself must not return NaN.
- Wrapping: the result is always in `[0,2*pi)`, checked over a sweep.
- **Asymmetry:** bearing(A,B) is NOT generally bearing(B,A) reversed — unlike distance. Assert this on the LAX/JFK pair, because a reader who assumes symmetry will write a bug.
- **A CORRECTION to the line above, established during execution.** This plan originally called the asymmetry test "the one property that distinguishes a correct implementation from one that swapped its arguments." That is FALSE, and it was measured: swapping `lat1` and `lat2` in the formula moves both the forward and reverse bearings but leaves the offset-from-reciprocal at 27.974 deg in both cases, identical, so the asymmetry assertion cannot see the swap. A NEAR-MERIDIONAL pair is the real discriminator — (10N,0) to (80N,1E) gives 0.185 deg true against 178.952 deg swapped. Include such a case; that is the actual transposition detector.

`test_rangeSolve.m`, against a synthetic monotonic function with a known root so the solver is tested independently of any trajectory:
- Converges to the analytic root within tolerance.
- Reports `converged = false` rather than throwing when the requested value is outside the bracket, and the `info` struct carries the achievable band.
- Handles the exact-endpoint case.
- Respects the tolerance: a tighter `tolM` gives a tighter result.
- Errors with a named identifier if `tLo >= tHi` — a malformed bracket.
- **CORRECTION established during execution:** an earlier version of this line also demanded an error when both endpoints sit on the same side of the target. That is the SAME condition as the out-of-bracket case two lines above, which must NOT throw. Contradictory. The return-with-band behaviour governs; carry the named identifier and the achievable band inside `info` so the caller can print a useful refusal.

- [ ] **Step 2: Run and confirm they fail**
- [ ] **Step 3: Write both functions**
- [ ] **Step 4: Run the suite** — expect `18 passed, 0 failed`, zero warnings.
- [ ] **Step 5: Prove the tests bite**

Mutations, restoring byte-identically (md5) after each: swap `lat1` and `lat2` in the bearing formula — the asymmetry test must FAIL; change the bisection to return the midpoint without checking the tolerance — a test must FAIL; drop the wrap so a negative bearing can be returned — the wrapping test must FAIL.

- [ ] **Step 6: Commit**

---

## Task 2: The `+viz` package — static figures

**Files:** Create `+coorbital/+viz/groundTrack.m`, `profilePlot.m`, `globe3D.m`, `tests/test_viz.m`.

This package was in the design spec from the first milestone and was deferred twice. The three existing entry scripts each carry three inline `figure` blocks; that plotting is the thing to extract.

**Interfaces produced:** each takes `(traj,veh,env)` plus an options struct and returns a figure handle. They must not call `figure` implicitly when handed a target axes — accept an optional `'Parent'` so a caller can compose.

- [ ] **Step 1: Read the inline plotting in all three entry scripts** — `run_glide.m`, `run_ballistic.m`, `run_boost_glide.m`. Record in your report what is genuinely common and what is script-specific. Extract only the common part; do not force a shared abstraction over things that legitimately differ.
- [ ] **Step 2: Write the tests**

Plotting tests are mostly smoke tests, so make them earn their place: assert the figure has the expected number of axes, that every axis has non-empty `XLabel`, `YLabel` and `Title` strings, that units appear in the label text, that the data plotted matches `traj` (compare `XData`/`YData` of the primary line against the trajectory arrays), and that the ground track marks launch, target and impact. A plot test that only checks "no error was thrown" is worth little; these must check the figure says what it claims.

Must run headless under `-batch` without displaying. Use `'Visible','off'`.

- [ ] **Step 3: Write the three functions**
- [ ] **Step 4: Retrofit the three existing entry scripts to call `+viz`** — and confirm their headline numbers are UNCHANGED. This is the risky step: the numbers must not move, because plotting must not touch the trajectory. Re-run all three and diff the summaries against the values in Global Constraints.
- [ ] **Step 5: Run the suite** — expect `19 passed, 0 failed`, zero warnings.
- [ ] **Step 6: Prove the tests bite** — mutate a label to omit its unit, and plot the wrong array; both must FAIL. Restore md5-verified.
- [ ] **Step 7: Commit**

---

## Task 3: The globe movie

**Files:** Create `+coorbital/+viz/globeMovie.m`, extend `tests/test_viz.m`.

**Interface produced:** `mv = coorbital.viz.globeMovie(traj,opts)` writing an MP4 and returning the file path plus frame count.

- [ ] **Step 1: Read pumpkyn's renderers before writing anything**

`~/Desktop/proj7/external/pumpkyn/src/+pumpkyn/+util/earth3D.m`, `sphere3D.m`, `stars3D.m`, `addFigureLogo.m`. Read their actual headers for signatures and options — do not guess. `earth3D` takes an options struct with a texture type (`'clouds'`, `'night'`, `'day'`) and a `www` argument; `stars3D` takes figure and axes handles, a background mode, an image file and a Julian date.

- [ ] **Step 2: Decide and document the degradation path**

The textures are local files in the pumpkyn tree, but `earth3D` has a `www` argument and may attempt a fetch. The movie must work with no network and must work headless. Determine what actually happens in both cases and write the fallback: a plain shaded sphere if the texture is unavailable, `'black'` background if the starfield is not. State in the header what each fallback looks like so a user seeing a plain sphere knows why.

- [ ] **Step 3: Write the tests**

- A short movie (10-20 frames) renders headless and the file exists with non-zero size.
- The frame count matches what was requested.
- The trajectory data drawn in the final frame matches the full trajectory, and in an intermediate frame matches the trajectory truncated at that time — this is what makes it a movie of the trajectory *developing* rather than a static arc spun on a turntable.
- Cleans up its figure; no figures left open after the call.
- Deletes its output in the test, or writes to the scratch directory, so the repo stays clean.

- [ ] **Step 4: Write `globeMovie`**

Requirements:
- Convert `(r,lon,lat)` to Cartesian for plotting: `x = r cos(lat) cos(lon)`, `y = r cos(lat) sin(lon)`, `z = r sin(lat)`. With `omegaE = 0` these are planet-fixed coordinates; say so in the header, because a reader will assume ECI.
- The trajectory grows frame by frame, with the vehicle marked at the current position and a trailing path.
- Colour or mark the three phases distinctly, so boost, glide and descent are visible as such.
- A time and altitude annotation per frame.
- `drawnow` before every `getframe` — otherwise frames capture partially-rendered figures. This is a documented gotcha on this machine.
- Frame rate, size, and whether to spin the camera are options with sane defaults.

Budget roughly 2-3 s per frame at 1920x1080 on this machine; a 200-frame movie takes 8-10 minutes. Default to something short and say in the header how to make it longer.

- [ ] **Step 5: Render a real movie and inspect it**

Render at least 60 frames of the full chain. Extract a frame with `VideoReader` and `imwrite`, then LOOK at the PNG — do not just check the file size. Report what you see: is the Earth oriented sensibly, is the trajectory where it should be, are the phases distinguishable, is anything clipped or upside down.

- [ ] **Step 6: Run the suite** — expect `19 passed, 0 failed` (extending the existing viz test), zero warnings.
- [ ] **Step 7: Commit**

---

## Task 4: `HGV/run_target.m` — the example script

**Files:** Create `HGV/run_target.m`, `tests/test_runTarget.m`. Modify `docs/README.md`.

This is the deliverable the user asked for. `HGV/run_boost_glide.m` is the reference for entry-script quality; read it and match it.

- [ ] **Step 1: Write the script**

The user block must hold, in human units with a trailing comment giving units and a sane range for each: launch latitude and longitude, target latitude and longitude, the vehicle and booster selection, the handoff and stop altitudes, the glide and descent control schedules, the model selections, movie on/off with its frame count, and the range-solve tolerance.

Below the fence, in order: convert to SI; compute the required range and the launch azimuth from the two points; bracket the reachable envelope; solve for the cutoff time; fly the solved trajectory; report; plot; optionally render the movie.

The summary must report the targeting solve as its own block — required range, achieved range, miss distance, solved cutoff time as both an absolute time and a fraction of full burn, the reachable envelope, and the iteration count. Then the usual trajectory summary in the format the other scripts use, the model list, the PLACEHOLDER caveat, and the validity note.

**Two things the summary must state plainly**, because both are limitations a reader will otherwise not notice:
- The azimuth solution is exact only for a non-rotating Earth. Say it, and say what would change with rotation on.
- The miss distance is the residual of the range solve along the great circle. Cross-range miss is zero by construction here because bank is zero throughout and the ground track is a great circle — that is a property of this configuration, not a general guarantee.

- [ ] **Step 2: Run it and paste the complete console output into your report**

Then read that output as the end user. Is every number labelled and united? Is the targeting block clear about what was solved and how well? Fix anything you would not want to receive.

- [ ] **Step 3: Verify the targeting independently**

Take the achieved impact latitude and longitude from the run and compute the great-circle distance from there to the target with `coorbital.util.greatCircle`. Confirm it matches the reported miss distance. That is an independent check of the solve, computed from the flown trajectory rather than from the solver's own bookkeeping.

Also verify the azimuth: confirm the initial heading in `traj.x(1,6)` equals `greatCircleBearing(launch,target)` to machine precision.

- [ ] **Step 4: Write the integration test**

Pin the headline numbers to about 0.01 percent, read at full printed precision. Assert the miss distance is within the requested tolerance. Assert the reachable-envelope refusal path fires with a clear message for a target that is too far, and for one that is too close. Assert the impact point is within tolerance of the target, computed independently with `greatCircle`.

**Fly a test case that is not equatorial and not due east.** A due-east equatorial case is provably blind to a lat/lon transposition, because `cos(d) = cos(lat2)cos(lon2)` is symmetric at the origin. This project has paid for that lesson once already.

- [ ] **Step 5: Prove the test bites**

Mutations, restoring byte-identically after each: transpose the arguments to `greatCircleBearing` so the azimuth is wrong; return the bracket midpoint without iterating so the miss is large; skip the separation link. Each must FAIL.

- [ ] **Step 6: Update the README** — the new entry script, the two new util functions, the `+viz` package, and the targeting method with its non-rotating-Earth limitation. Every command printed must be run and shown to work.

- [ ] **Step 7: Run the full suite** — expect `20 passed, 0 failed`, zero warnings, and all four entry scripts producing their headline numbers.

- [ ] **Step 8: Commit**

---

## Task 5: `BM/run_ballistic_target.m` — the ballistic point-to-point script

**Files:** Create `BM/run_ballistic_target.m`, `tests/test_runBallisticTarget.m`. Modify `docs/README.md`.

The framework already exists. `BM/run_ballistic.m` flies boost -> coast -> impact and is validated against the Keplerian closed form to 8.5e-13, and Task 1 supplies the two pieces that turn "launch plus azimuth" into "launch plus destination". This task is the targeting wrapper applied to the ballistic chain.

**The one thing that genuinely differs from the HGV case.** A pure ballistic trajectory has TWO solutions for every range short of maximum — a lofted arc and a depressed arc, either side of the minimum-energy solution — and a ballistic-missile user expects to choose. So the ranging control here is the **loft angle**, the pitch program's terminal attitude, rather than thrust-termination, and the user block carries a branch selector:

```matlab
       branch = 'minimum-energy';   % 'minimum-energy' | 'lofted' | 'depressed'
```

- `minimum-energy` is the textbook default. Bracket the max-range loft angle first, then solve on whichever side of it the target falls.
- `lofted` and `depressed` solve the high and low branches. Each must REPORT which branch it actually landed on, measured from the flown apogee and flight time, rather than trusting that the bracket kept it there.
- Beyond maximum range, refuse with the maximum range stated, the loft angle achieving it, and the shortfall.

Report BOTH branch solutions when they exist, even though only one is flown, so the user can see the trade: give flight time, apogee, impact speed and impact flight-path angle for each. The lofted arc buys steeper, slower impact at longer flight time; the depressed arc is faster and flatter.

Everything else follows Task 4: same user-block discipline, same independent verification of the impact point with `greatCircle`, same off-axis test case, same mutation proofs.

---

## Out of scope

- Rotating-Earth targeting, which needs an outer azimuth iteration.
- Cross-range steering with bank, and any closed-loop terminal guidance.
- Solving for loft angle rather than cutoff time, and the two-branch structure that implies.
- Multi-stage boosters, throttling, and PEG or VOA boost guidance.
- Anything that modifies pumpkyn.
