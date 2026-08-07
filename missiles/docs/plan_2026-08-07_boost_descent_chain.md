# Boost, Descent, and the Full Three-Phase Chain — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the validated glide propagator to a full boost → glide → descent trajectory, with thrust and mass depletion during boost, a terminal descent phase, and entry scripts that fly the whole chain.

**Architecture:** Boost adds a seventh state, mass. Rather than build state-mapping machinery at phase junctions, every phase in a chain runs on the same 7-state vector: `boost3DOF` is natively 7-state, and the already-verified `glide3DOF` is reused untouched through a one-line `massConstant` adapter that appends `dm/dt = 0`. `phaseRun` generalises from a hard-coded 6 to `numel(x0)`. Propulsion joins atmosphere, gravity, and aerodynamics as an injected model handle, so a pressure-corrected or throttled engine drops in without editing the equations of motion.

**Tech Stack:** MATLAB R2025b (`/Applications/MATLAB_R2025b.app/bin/matlab`), `ode45`, the existing `missiles/+coorbital` library.

---

## Global Constraints

- Style: pumpkyn house style per `~/Desktop/proj7/doc/pumpkyn_style_guide.md` — `%%` header quartet (`Purpose`/`Inputs`/`Outputs`/`Revision History`) with aligned three-column input and output blocks stating units, closed by the `%% ------------------------ Begin Code Sequence ---------------------------` divider; `=` signs column-aligned with right-justified names at column 20; colon-terminated `%%` step comments; `if nargin == 0` self-demo in every library function calling itself by full namespace.
- Add a `%% References:` block whenever a function's math has a source.
- Never use `i` or `j` as loop or index variables — they are the imaginary unit. Use `k`, `ii`, `idx`, `kp`, `kt`, or a meaningful name.
- Norms via `sqrt(sum(x.^2,dim))`. **Never `norm`.**
- **No lint-suppression pragmas.** Never write `%#ok<NASGU>`, `%#ok<AGROW>`, or any other `%#ok` directive. If a variable is unused, delete it or capture it with `~`. If a loop grows an array, preallocate it.
- **No hard-coded physical constants** outside `coorbital.util.missileConst()`. Vehicle and booster properties belong in the vehicle structs. Analytic reference literals inside tests are legitimate and expected.
- SI units throughout the library: metres, seconds, kilograms, radians, kelvin, newtons. Entry scripts may accept degrees, kilometres, and tonnes in the user block and convert immediately.
- All angles in the state vector are radians.
- Every library function lives in `missiles/+coorbital/+<subpkg>/`, one function per file, file name equal to function name, called by full namespace.
- Author line in every header: `Michael Casey` with the date 08/07/2026, then `Copyright 2026 Coorbital, Inc.`
- **Do not modify `+coorbital/+eom/glide3DOF.m`.** Its six equations were re-derived symbolically term by term and every sign independently confirmed. It is reused, never edited. If you believe it needs a change, stop and report it.
- Run MATLAB headlessly as: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/missiles'); <command>" 2>&1 | grep -vE "Home License|personal use|academic, research|organizational use"`
- Run the suite as `run('tests/run_tests')`. A bare `tests/run_tests` does NOT work — `-batch` parses it as division of two undefined names. Allow 300 s; the suite integrates many trajectories.
- **The suite must stay green with zero warnings at every commit.** It currently reports `10 passed, 0 failed`. A MATLAB warning is a defect.

---

## Inherited interfaces — what already exists and is verified

| Function | Signature | Notes |
|---|---|---|
| `coorbital.util.missileConst()` | `c` | `muE`, `rE`, `omegaE`, `rho0`, `Hscale`, `T0`, `Rair`, `gamAir`, `g0`. Single source of truth. |
| `coorbital.util.vehicleDefaults()` | `veh` | `mass`, `Sref`, `CL`, `LD`, `noseRadius`. PLACEHOLDER values. |
| `coorbital.util.greatCircle(lat1,lon1,lat2,lon2)` | `d` | Haversine, stable at small separation. |
| `coorbital.atmos.expAtmos(h)` | `[rho,P,T,a]` | `h` geometric altitude (m). Isothermal, so `a` is a constant 316.9677 m/s. |
| `coorbital.grav.sphereGrav(r,lat)` | `[gr,gLat]` | **`gr` is a POSITIVE MAGNITUDE of inward attraction**, not a signed component. `gLat` identically zero. |
| `coorbital.aero.constLD(alpha,mach,veh)` | `[CL,CD]` | Scalar out. Ignores `mach` at this fidelity. `CD = CL/LD`. |
| `coorbital.eom.glide3DOF(t,x,u,veh,env)` | `xdot [6 x 1]` | `x = [r,lon,lat,V,gamma,psi]`, `u = [alpha,sigma]`. Verified symbolically. **Do not edit.** |
| `coorbital.guide.prescribed(t,x,sched)` | `u [2 x 1]` | Linear interpolation on `sched.tGrid`, clamped outside. Ignores `x`. |
| `coorbital.prop.eventAltitude(t,x,hStop)` | `[value,isterminal,direction]` | Terminal, `direction = -1` (descending crossings only). Uses `x(1)`. |
| `coorbital.prop.phaseRun(phases,x0,veh,env)` | `traj` | `t`, `x`, `u`, `phaseIdx`, `junction`. Each phase referenced to its own `tspan(1)`. |

State and sign conventions, restated because a boost EOM is where they get broken: `gamma` positive up, `psi` clockwise from north, `gr` a positive inward magnitude that the EOM subtracts.

---

## File Structure

| File | Responsibility |
|---|---|
| `+coorbital/+eom/massConstant.m` | Adapter: wraps a 6-state EOM as 7-state with `dm/dt = 0` |
| `+coorbital/+eom/boost3DOF.m` | 7-state powered-flight equations of motion |
| `+coorbital/+prop/constThrust.m` | Constant-mass-flow engine with ambient back-pressure correction |
| `+coorbital/+prop/eventBurnout.m` | ODE event: propellant exhausted |
| `+coorbital/+prop/eventApogee.m` | ODE event: flight-path angle crosses zero descending |
| `+coorbital/+guide/pitchProgram.m` | Prescribed pitch-attitude program; commands `alpha` from `theta - gamma` |
| `+coorbital/+util/boosterDefaults.m` | Booster parameters (wet/dry mass, thrust, Isp, exit area) |
| `+coorbital/+prop/phaseRun.m` | **Modified**: generalised to `numel(x0)` states |
| `BM/run_ballistic.m`, `BM/vehicle_bm.m` | Boost → ballistic coast → impact entry script |
| `HGV/run_boost_glide.m` | Full boost → glide → descent chain |
| `tests/test_boost3DOF.m` | Reduction-to-glide, Tsiolkovsky, mass linearity |
| `tests/test_boostEvents.m` | Burnout and apogee events |
| `tests/test_pitchProgram.m` | Pitch-program guidance |
| `tests/test_fullChain.m` | Three-phase integration test |
| `docs/README.md`, `docs/DESIGN.md`, `docs/LESSONS_LEARNED.md` | **Modified** |

---

## The boost equations of motion

State `x = [r, lon, lat, V, gamma, psi, m]`, control `u = [alpha, sigma]`.

Thrust magnitude `T` acts at angle of attack `alpha` from the velocity vector, in the plane set by bank angle `sigma` — the same convention lift already uses, so thrust and lift add directly:

```
dr/dt     =  V sin(gamma)
dlon/dt   =  V cos(gamma) sin(psi) / (r cos(lat))
dlat/dt   =  V cos(gamma) cos(psi) / r
dV/dt     =  (T cos(alpha) - D)/m - gr sin(gamma) + gLat cos(gamma) cos(psi)
dgamma/dt =  (T sin(alpha) + L) cos(sigma)/(m V) - (gr/V - V/r) cos(gamma)
             - gLat sin(gamma) cos(psi)/V
dpsi/dt   =  (T sin(alpha) + L) sin(sigma)/(m V cos(gamma))
             + V cos(gamma) sin(psi) tan(lat)/r - gLat sin(psi)/(V cos(gamma))
dm/dt     = -mdot
```

plus the same rotating-Earth Coriolis and centrifugal terms `glide3DOF` already carries, gated on `env.omegaE`.

**The mass in the denominators is the STATE `x(7)`, never `veh.mass`.** This is the single most likely bug in the file and it is silent: using `veh.mass` gives a plausible trajectory that simply ignores propellant burn. Task 6 mutates exactly this.

Setting `T = 0` and `mdot = 0` must reduce the first six components to `glide3DOF` exactly. That is the cheapest and strongest check available, because it cross-validates against a file whose every sign was confirmed symbolically.

---

## Task 1: Generalise phaseRun to any state dimension

**Files:**
- Modify: `missiles/+coorbital/+prop/phaseRun.m`
- Create: `missiles/+coorbital/+eom/massConstant.m`
- Test: `missiles/tests/test_phaseRun.m` (extend)

**Interfaces:**
- Consumes: `phaseRun(phases,x0,veh,env)` as it stands today.
- Produces: `phaseRun` accepting `x0` of any length `nx`, returning `traj.x` as `[N x nx]`. `coorbital.eom.massConstant(eomFn)` returning a function handle usable as a phase `eom`.

- [ ] **Step 1: Find every place phaseRun assumes six states**

Run: `grep -n "6" missiles/+coorbital/+prop/phaseRun.m`
Read each hit and classify it as a genuine state-dimension assumption or a coincidence. Record the list in your report before changing anything.

- [ ] **Step 2: Write the failing test**

Add to `missiles/tests/test_phaseRun.m`:

```matlab
%% A seven-state phase runs and returns a seven-column trajectory:
           eom7 = coorbital.eom.massConstant(@coorbital.eom.glide3DOF);
          ph7.eom = eom7;
        ph7.guide = @(t,x) coorbital.guide.prescribed(t,x,schedZero);
    ph7.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,20e3);
        ph7.tspan = [0 4000];
             x07 = [c.rE + 60e3; 0; 0; 6000; deg2rad(-1); deg2rad(90); veh.mass];
           traj7 = coorbital.prop.phaseRun(ph7,x07,veh,env);
    assert(size(traj7.x,2) == 7,'expected a 7-column trajectory, got %d',size(traj7.x,2));
    assert(all(abs(traj7.x(:,7) - veh.mass) < 1e-12), ...
        'mass must be exactly constant under massConstant');
```

- [ ] **Step 3: Run it and confirm it fails**

Expected: failure, because `massConstant` does not exist.

- [ ] **Step 4: Write massConstant**

```matlab
function eomFn = massConstant(baseEom)
%% Purpose:
%
%  Adapt a six-state atmospheric-flight EOM to the seven-state vector used by
%  a chain that includes a powered phase, by appending dm/dt = 0. Unpowered
%  phases carry mass so that every phase in a chain shares one state vector,
%  which removes the need for a state mapping at the junctions.
%
%% Inputs:
%
%  baseEom          Function handle             Six-state EOM with the
%                                               signature xdot = f(t,x,u,veh,env)
%
%% Outputs:
%
%  eomFn            Function handle             Seven-state EOM with the same
%                                               signature; component 7 is mass
%                                               (kg) and its derivative is zero
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
                 c = coorbital.util.missileConst();
               veh = coorbital.util.vehicleDefaults();
         env.atmos = @coorbital.atmos.expAtmos;
          env.grav = @coorbital.grav.sphereGrav;
          env.aero = @coorbital.aero.constLD;
        env.omegaE = 0;
              fD = coorbital.eom.massConstant(@coorbital.eom.glide3DOF);
              xD = [c.rE + 40e3; 0; 0; 4000; deg2rad(-2); deg2rad(90); veh.mass];
           xdotD = fD(0,xD,[0;0],veh,env);
    fprintf('dm/dt = %.1f kg/s over a %d-state vector\n',xdotD(7),numel(xdotD));
             eomFn = [];
    return;
end

             eomFn = @(t,x,u,veh,env) [baseEom(t,x(1:6),u,veh,env); 0];
end
```

- [ ] **Step 5: Generalise phaseRun**

Replace the hard-coded state width with `nx = numel(x0)`. The control width should likewise come from the first guide evaluation rather than being assumed to be 2. Update the header's `x0` and `traj.x` size annotations to `[nx x 1]` and `[N x nx]`, and add a `%% Notes:` sentence explaining that all phases in one call must share a state dimension, and that `coorbital.eom.massConstant` exists to lift a six-state EOM into a seven-state chain.

- [ ] **Step 6: Run the whole suite**

Expected: `10 passed, 0 failed`, zero warnings. **Every existing test must still pass unchanged** — this is a backward-compatible generalisation, and `run_glide` still runs six-state.

- [ ] **Step 7: Prove the generalisation is real**

Mutate `phaseRun` to hard-code 6 again in the `xSeg` preallocation. Confirm the new seven-state test FAILS. Restore, verify md5 byte-identical, rerun clean.

- [ ] **Step 8: Commit**

```bash
git add missiles/+coorbital/+prop/phaseRun.m missiles/+coorbital/+eom/massConstant.m missiles/tests/test_phaseRun.m
git commit -m "missiles: generalise phaseRun to any state dimension"
```

---

## Task 2: Propulsion model and booster parameters

**Files:**
- Create: `missiles/+coorbital/+prop/constThrust.m`, `missiles/+coorbital/+util/boosterDefaults.m`
- Test: `missiles/tests/test_constThrust.m`

**Interfaces:**
- Produces: `[T,mdot] = coorbital.prop.constThrust(t,P,veh)` and `bst = coorbital.util.boosterDefaults()`.

The propulsion family signature, joining the existing three:

| Family | Signature |
|---|---|
| propulsion | `[T,mdot] = f(t,P,veh)` — `t` time since phase start (s), `P` ambient pressure (Pa); returns thrust (N) and mass flow (kg/s, positive) |

Passing ambient pressure rather than altitude is deliberate: it is what the physics needs, the EOM already has `P` from `env.atmos`, and it keeps the propulsion model independent of the atmosphere model.

- [ ] **Step 1: Write the failing test**

`missiles/tests/test_constThrust.m` must check, against hand-computed literals:
- Vacuum thrust at `P = 0` equals `bst.thrustVac` exactly.
- Sea-level thrust equals `thrustVac - Aexit*101325`, computed by hand and written as a literal.
- `mdot` is `thrustVac/(Isp*g0)` and does NOT vary with `P` — a fixed propellant flow is the physical assumption, and back pressure changes thrust, not flow.
- `mdot > 0` (the sign convention is positive flow; the EOM applies the minus).
- Burn time `(massWet - massDry)/mdot` is positive and of a sane order; assert it against a hand-computed literal.

- [ ] **Step 2: Run it and confirm it fails**

- [ ] **Step 3: Write boosterDefaults**

A single-stage solid-equivalent booster sized to loft the existing 900 kg glide vehicle. Fields: `massWet`, `massDry`, `thrustVac`, `Isp`, `Aexit`, plus `Sref` and aero fields so the stack can fly with `constLD` during boost. **Every value is a PLACEHOLDER** and must be commented as such in the same style `vehicleDefaults` uses, with a header note saying they are open-literature order-of-magnitude figures and must be replaced before any result leaves the machine.

Choose values that give a physically sensible boost: a thrust-to-weight ratio at liftoff comfortably above 1 (aim for roughly 2 to 3), a burn time of order a minute, and a burnout speed high enough to set up the glide the existing entry script starts from. State the resulting T/W and burn time in the header as computed comments so a reader can sanity-check the set at a glance.

- [ ] **Step 4: Write constThrust**

`T = max(bst.thrustVac - bst.Aexit*P, 0)` — the `max` guards the unphysical negative-thrust case if someone supplies a huge exit area. `mdot = bst.thrustVac/(bst.Isp*c.g0)`, with `g0` from `missileConst`, never a literal 9.80665.

- [ ] **Step 5: Run the suite** — expect `11 passed, 0 failed`, zero warnings.

- [ ] **Step 6: Prove the test bites**

Mutate `Isp` by ×1.1 and confirm `test_constThrust` FAILS on `mdot`. Restore md5-verified.

- [ ] **Step 7: Commit**

```bash
git add missiles/+coorbital/+prop/constThrust.m missiles/+coorbital/+util/boosterDefaults.m missiles/tests/test_constThrust.m
git commit -m "missiles: constant-mass-flow propulsion model and booster parameters"
```

---

## Task 3: The boost equations of motion

**Files:**
- Create: `missiles/+coorbital/+eom/boost3DOF.m`
- Test: `missiles/tests/test_boost3DOF.m` (created in Task 6; write only a smoke assertion here)

**Interfaces:**
- Consumes: `env.atmos`, `env.grav`, `env.aero`, `env.omegaE`, and the new `env.prop`.
- Produces: `xdot = coorbital.eom.boost3DOF(t,x,u,veh,env)`, `x [7 x 1]`, `xdot [7 x 1]`.

- [ ] **Step 1: Read glide3DOF end to end**

You are writing its powered sibling. Match its structure, its variable names, its step-comment cadence, and its guard block exactly, so the two files read as a pair and a reviewer can diff them. Do not edit `glide3DOF.m`.

- [ ] **Step 2: Write boost3DOF**

Requirements, each of which a later task tests:
- The state is `[r, lon, lat, V, gamma, psi, m]`. The mass used in every acceleration denominator is `x(7)`. `veh.mass` must not appear anywhere in the file.
- Thrust and mass flow come from `env.prop(t,P,veh)`, with `P` from `env.atmos`.
- Aerodynamic forces use the same `env.aero` call as `glide3DOF`.
- Carry the same three singularity guards with the same error identifiers, adding a fourth: `coorbital:boost3DOF:massDepleted` when `x(7) <= 0`, because a negative mass silently inverts every acceleration.
- Carry the rotating-Earth terms gated on `env.omegaE`, identical in form to `glide3DOF`.
- `%% References:` citing Vinh for the flight mechanics and any standard text for the rocket equation.

- [ ] **Step 3: Smoke-check by hand**

Evaluate at a simple state — equatorial, due east, level, zero controls — and confirm `dm/dt` equals `-mdot` and `dV/dt` equals `(T - D)/m - gr sin(gamma)` computed separately. Paste both into your report.

- [ ] **Step 4: Run the suite** — must remain green; you have added no test yet.

- [ ] **Step 5: Commit**

```bash
git add missiles/+coorbital/+eom/boost3DOF.m
git commit -m "missiles: 3-DOF powered-flight equations of motion"
```

---

## Task 4: Burnout and apogee events

**Files:**
- Create: `missiles/+coorbital/+prop/eventBurnout.m`, `missiles/+coorbital/+prop/eventApogee.m`
- Test: `missiles/tests/test_boostEvents.m`

- [ ] **Step 1: Write the failing test**

`test_boostEvents.m` must assert:
- `eventBurnout` returns `[value,isterminal,direction]` with `value = x(7) - massDry`, `isterminal = 1`, `direction = -1` (mass only ever decreases, but the one-sided form documents intent and guards a future restart).
- `eventApogee` fires on `gamma` crossing zero **descending**, so `value = x(5)`, `direction = -1`. A vehicle still climbing must not trigger it.
- Behaviourally: propagate a boost phase and confirm it terminates with `x(7)` equal to `massDry` to within 1e-6 kg, and that the burn time matches `(massWet - massDry)/mdot` to 1e-6 s. That second assertion is the one that pins the event to the physics rather than to itself.
- Behaviourally: propagate a lofted unpowered arc and confirm `eventApogee` fires at the top, with `gamma` crossing zero and `dr/dt` changing sign there.

- [ ] **Step 2: Run it and confirm it fails**

- [ ] **Step 3: Write both events**

Follow `eventAltitude`'s structure exactly, including its `if nargin == 0` self-demo that prints the three returned values.

- [ ] **Step 4: Run the suite** — expect `12 passed, 0 failed`.

- [ ] **Step 5: Prove the tests bite**

Flip `direction` on `eventApogee` from `-1` to `+1` and confirm the test FAILS. Restore md5-verified. Then set `eventBurnout`'s `value` to `x(7)` alone (dropping `massDry`) and confirm the burn-time assertion FAILS. Restore md5-verified.

- [ ] **Step 6: Commit**

```bash
git add missiles/+coorbital/+prop/eventBurnout.m missiles/+coorbital/+prop/eventApogee.m missiles/tests/test_boostEvents.m
git commit -m "missiles: burnout and apogee ODE events"
```

---

## Task 5: Pitch-program guidance

**Files:**
- Create: `missiles/+coorbital/+guide/pitchProgram.m`
- Test: `missiles/tests/test_pitchProgram.m`

**Interfaces:**
- Produces: `u = coorbital.guide.pitchProgram(t,x,sched)` with `sched` carrying `tGrid`, `theta` (commanded pitch attitude, rad) and `sigma`.

Unlike `prescribed`, this guide **uses the state**: the commanded angle of attack is `alpha = theta(t) - gamma`, where `gamma = x(5)`. That is the classic prescribed-pitch boost law and it exercises the `guide(t,x)` signature that `prescribed` leaves unused.

- [ ] **Step 1: Write the failing test**

Assert against hand-computed values:
- At a grid point with `theta = 30 deg` and a state with `gamma = 10 deg`, the returned `alpha` is exactly `20 deg` in radians.
- Interpolation between grid points is linear in `theta`, verified at a midpoint.
- Outside the grid the schedule clamps, matching `prescribed`'s documented behaviour.
- `sigma` passes through unchanged.
- An `alpha` limit: the function must clamp `|alpha|` to `sched.alphaMax` if that field is present, and a test must confirm a commanded 90 deg attitude at `gamma = 0` clamps rather than returning a physically absurd angle of attack. State in the header what happens when the field is absent (no clamping).

- [ ] **Step 2: Run it and confirm it fails**

- [ ] **Step 3: Write pitchProgram**

- [ ] **Step 4: Run the suite** — expect `13 passed, 0 failed`.

- [ ] **Step 5: Prove the test bites**

Change `alpha = theta - gamma` to `alpha = theta + gamma` and confirm the test FAILS. Restore md5-verified.

- [ ] **Step 6: Commit**

```bash
git add missiles/+coorbital/+guide/pitchProgram.m missiles/tests/test_pitchProgram.m
git commit -m "missiles: prescribed pitch-attitude boost guidance"
```

---

## Task 6: Analytic validation of the boost EOM

**Files:**
- Create: `missiles/tests/test_boost3DOF.m`
- Modify: `missiles/docs/LESSONS_LEARNED.md`

This is the task that decides whether the boost dynamics are trustworthy. Three checks, in increasing strength.

- [ ] **Step 1: Reduction to glide**

With `env.prop` returning `T = 0, mdot = 0`, `boost3DOF`'s first six components must equal `glide3DOF` at the same state **to machine precision** (1e-12 relative or better), with `x(7)` set to `veh.mass`. Test at several states, including off-equator, non-eastward heading, nonzero bank, and with `env.omegaE` both zero and the real value, so the rotating branch is covered too.

This is the strongest structural check available: it validates every shared term against a file whose signs were confirmed symbolically, leaving only the thrust and mass terms to check separately.

- [ ] **Step 2: Tsiolkovsky**

In vacuum (`env.atmos` returning `rho = 0`, `P = 0`), with gravity zeroed (an `env.grav` returning `[0,0]`), and thrust along the velocity vector (`alpha = 0`), integrating from `massWet` to `massDry` must give

```
        deltaV = Isp * g0 * log(massWet/massDry)
```

to 1e-8 relative. Compute the reference from `boosterDefaults` and `missileConst`, never by rearranging the propagated result. Report the measured agreement as a percentage.

**Note the trap this test avoids and the one it does not.** It avoids self-reference by using a closed form. It does *not* escape the shared-constant blindness recorded in `LESSONS_LEARNED.md`: `Isp` and `g0` appear on both sides and cancel, so this test cannot detect a wrong `Isp`. Say so in the test header, and confirm `Isp` is pinned in `test_constThrust` — that pin is the only defence, exactly as `Hscale` was for Allen–Eggers.

- [ ] **Step 3: Mass linearity**

With constant mass flow, `m(t) = massWet - mdot*t` exactly. Assert against the analytic line to 1e-9 relative over the whole burn.

- [ ] **Step 4: Run the suite** — expect `14 passed, 0 failed`.

- [ ] **Step 5: Prove the tests bite — three mutations**

Restore byte-identically (md5) after each and put every output in your report:
1. In `boost3DOF`, replace `x(7)` with `veh.mass` in the acceleration denominators. **Tsiolkovsky must FAIL.** This is the silent bug the whole task exists to catch.
2. Change `dm/dt` from `-mdot` to `+mdot`. A test must FAIL.
3. Change the thrust term from `T*cos(alpha)` to `T` in `dV/dt`. The reduction-to-glide test will not catch this (thrust is zero there) — confirm which test does. If none does, that is a coverage gap: say so and add a nonzero-`alpha` thrust-projection assertion.

- [ ] **Step 6: Record the lesson**

Add to `LESSONS_LEARNED.md` whatever mutation 3 taught about which checks cover the thrust projection, in the general form: a reduction test validates only the terms that survive the reduction.

- [ ] **Step 7: Commit**

```bash
git add missiles/tests/test_boost3DOF.m missiles/docs/LESSONS_LEARNED.md
git commit -m "missiles: analytic validation of powered flight"
```

---

## Task 7: Ballistic-missile entry script

**Files:**
- Create: `missiles/BM/run_ballistic.m`, `missiles/BM/vehicle_bm.m`
- Test: `missiles/tests/test_runBallistic.m`

The simplest full chain, and the one with a closed-form cross-check: **boost → unpowered coast → impact.** Once the atmosphere is thin the coast is nearly Keplerian, so the range can be sanity-checked against a vacuum ballistic trajectory.

- [ ] **Step 1: Write the entry script**

Same standard as `run_glide.m`, which is the reference for what "clean entry script" means here — read it first and match it:
- A fenced `%% USER PARAMETERS:` block holding launch site, launch azimuth, pitch program, booster and vehicle selection, and stop altitude, in human units with a trailing comment giving units and a sane range for every entry.
- One clearly marked conversion block below the fence. Nothing below it needs editing for a routine run.
- Three phases: boost (terminated by `eventBurnout`), coast (terminated by `eventApogee`, recorded but not stopping the run), and descent to impact (terminated by `eventAltitude` at 0).
- A printed summary in the same format `run_glide` uses, reporting burnout conditions, apogee, range, flight time, peak deceleration, and terminal state, plus the model list, the PLACEHOLDER caveat, and a validity note.
- Report the termination reason for every phase; an unexpected termination must say so rather than quietly producing a short trajectory.

- [ ] **Step 2: Run it and read the output as the user**

Paste the complete console output into your report. Every number labelled, every unit stated. Fix anything you would not want to receive.

- [ ] **Step 3: Cross-check the range**

For a vacuum ballistic trajectory over a spherical Earth with burnout speed `V`, flight-path angle `gamma`, at radius `r`, the range angle has a closed form. Compute it independently and compare against the propagated ground range. Agreement will be imperfect because the real trajectory has drag and a finite boost arc — state the measured difference and judge whether it is the right size for those effects, exactly as the Sänger check was handled for the glide. Do not assert a tolerance you cannot justify from the physics.

- [ ] **Step 4: Write the integration test**

`test_runBallistic.m` pins the headline numbers to about 0.01 percent, asserts nominal termination, asserts mass is exactly constant after burnout, and asserts the three phases appear in `traj.phaseIdx` in order. Read the actual values from a run at full printed precision rather than copying rounded figures.

**Remember the lat/lon symmetry trap from the glide milestone**: a due-east equatorial trajectory is provably blind to a `greatCircle` argument transposition, because with the entry point at the origin `cos(d) = cos(lat2)cos(lon2)` is symmetric under exchange. Fly this test off-axis, or assert something that breaks the symmetry.

- [ ] **Step 5: Run the suite** — expect `15 passed, 0 failed`.

- [ ] **Step 6: Prove the test bites**

Transpose the `greatCircle` arguments and confirm `test_runBallistic` FAILS. Restore md5-verified.

- [ ] **Step 7: Commit**

```bash
git add missiles/BM/run_ballistic.m missiles/BM/vehicle_bm.m missiles/tests/test_runBallistic.m
git commit -m "missiles: ballistic entry script, boost to coast to impact"
```

---

## Task 8: Descent phase and the full boost-glide-descent chain

**Files:**
- Create: `missiles/HGV/run_boost_glide.m`
- Test: `missiles/tests/test_fullChain.m`

The descent phase is not new machinery — it is a third phase with its own control schedule and a lower termination altitude, which the phase driver already supports. What is new is flying all three together and checking the seams.

- [ ] **Step 1: Write the entry script**

Boost → glide → descent, on one 7-state vector, with `massConstant` lifting `glide3DOF` for the two unpowered phases. The user block must let the operator set each phase's schedule and handoff altitude independently. Follow `run_glide.m`'s layout exactly.

Descent is distinguished from glide by its control schedule — a steeper, higher-drag terminal segment — and by terminating at impact rather than at a handoff altitude. Document in the script what physically distinguishes the two phases, so a reader does not conclude the split is arbitrary.

- [ ] **Step 2: Run it and paste the full console output into your report**

- [ ] **Step 3: Write the integration test**

`test_fullChain.m` must assert:
- Time is strictly monotonic across both junctions, with no anomalous gap.
- The state is continuous across both junctions, compared against an **independent** reference. Use a different solver at a tighter tolerance — `ode89` at 1e-12 — because a same-solver re-integration from the same initial state reproduces bitwise output and silently restores circularity. This was established in the glide milestone; do not repeat that mistake.
- Mass decreases strictly during boost and is exactly constant afterwards.
- Mass at the first junction equals `massDry` to 1e-6 kg.
- Three phases appear in `traj.phaseIdx`, in order, each non-empty.
- The headline results pinned to about 0.01 percent.

- [ ] **Step 4: Run the suite** — expect `16 passed, 0 failed`, zero warnings.

- [ ] **Step 5: Prove the test bites**

Restore byte-identically after each: (i) drop the `massConstant` wrapper from the glide phase so mass is not carried — a test must FAIL; (ii) reorder the phases so descent precedes glide — a test must FAIL.

- [ ] **Step 6: Commit**

```bash
git add missiles/HGV/run_boost_glide.m missiles/tests/test_fullChain.m
git commit -m "missiles: full boost-glide-descent chain"
```

---

## Task 9: Documentation

**Files:**
- Modify: `missiles/docs/README.md`, `missiles/docs/DESIGN.md`, `missiles/docs/LESSONS_LEARNED.md`

- [ ] **Step 1: Update README**

Every command you print must be RUN and shown to work before you commit it. Every number must be produced by running the code. Cover: the two new entry scripts; the seven-state convention and why unpowered phases carry mass; the propulsion model signature added to the model table; the new events and guidance; and the updated validation results with their measured agreements.

- [ ] **Step 2: Update the DESIGN.md as-built section**

The existing dated as-built section records what the glide milestone shipped. Add a second dated entry for this milestone rather than editing the first — the record of what was true at each stage is the point. Note specifically that `stateConvert` remains unbuilt and explain why: a uniform seven-state vector removed the need for it, which is a design change from the original spec and should be recorded as one.

- [ ] **Step 3: Update LESSONS_LEARNED**

- [ ] **Step 4: Verify and commit**

Re-run every command in the README. Confirm the suite is green.

```bash
git add missiles/docs/README.md missiles/docs/DESIGN.md missiles/docs/LESSONS_LEARNED.md
git commit -m "missiles: document the boost, descent, and full-chain milestone"
```

---

## Out of scope for this plan

- PEG and VOA closed-loop boost guidance. Prescribed pitch only; `hyperFLIGHT` uses both and they are the natural next milestone.
- Multi-stage boosters. One stage, one burn.
- Terminal guidance laws. The descent phase flies a prescribed schedule, not a homing law.
- Aerothermal heating. `noseRadius` remains a carried placeholder.
- The `+viz` package, `hgv_dynamics_note.tex`, and `software_design.tex`.
- Fidelity increments: rotating Earth on by default, J2, geodetic altitude, tabulated aero, US76 atmosphere. All still one-line swaps by design.
- Phase 2 trajectory optimization.
