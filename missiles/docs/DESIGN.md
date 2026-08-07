# Missile Trajectory Library — Design Spec

*Coorbital, Inc. — 2026-08-06*

Design for a MATLAB library that generates boost / glide / descent trajectories for
hypersonic glide vehicles, with ballistic and ICBM variants sharing the same core.
Written to be read before the implementation plan; the math derivations live in
`hgv_dynamics_note.tex` and the interface detail in `software_design.tex`.

---

## 1. Purpose and scope

Generate full boost → glide → descent trajectories, in-house, as an alternative to
Darin's `hyperFLIGHT` tool.

**Phase 1 (this spec): prescribed-control simulation.** Given a vehicle, a launch
state, and a control schedule, propagate the trajectory and report it. No optimization.

**Phase 2 (planned, not specified here): optimal trajectories.** Solve for the control
history that maximizes range or satisfies terminal conditions, reusing the existing
`orbit_transfer/verify_common` collocation and PMP machinery.

### Capability target

`hyperFLIGHT` sets the bar (from *Coorbital HGV Capabilities*, slide 2):

| Phase | hyperFLIGHT | Phase 1 here | Later |
|---|---|---|---|
| Boost | VOA + PEG, 3-DOF | prescribed pitch program | PEG, then VOA |
| Glide | 6-state oblate-Earth, 3-DOF | 6-state spherical, 3-DOF | oblate + J2–J4 |
| Terminal | adaptive guidance on CL and bank | prescribed CL and bank | adaptive |

Slide 10 states a throughput requirement that shapes the architecture:
*"Computationally Efficient Open-loop Controllers Create Multiple Trajectories/Second"*
— 217 skip-glide trajectories in a single figure. **The propagator must be fast enough
to run in batch**, which is the main reason prescribed-control simulation comes first
and optimization second.

### Non-goals

- 6-DOF attitude dynamics. Everything here is 3-DOF point-mass.
- Aerothermal design. Heating rate is computed as a *constraint metric*, not a
  thermal-protection model.
- Any GUI. Entry scripts with a parameter block, not a front end.
- Re-implementing what pumpkyn already does correctly (see §3).

---

## 2. Architecture

```
missiles/
  +coorbital/
    +atmos/   expAtmos.m                    density, pressure, temperature, sound speed
    +aero/    constLD.m                     CL, CD from alpha and Mach
    +grav/    sphereGrav.m                  gravitational acceleration
    +eom/     glide3DOF.m boost3DOF.m ballistic3DOF.m
    +guide/   prescribed.m                  control schedule evaluation
    +prop/    phaseRun.m                    multi-phase driver, event handoff
              eventAltitude.m eventApogee.m eventBurnout.m
    +viz/     trajPlot.m groundTrack.m profilePlot.m
    +util/    stateConvert.m vehicleDefaults.m
  HGV/    run_glide.m      vehicle_hgv.m
  BM/     run_ballistic.m  vehicle_bm.m
  ICBM/   run_icbm.m       vehicle_icbm.m
  docs/   tests/
```

One shared library; the vehicle folders hold **only** parameter files and entry
scripts. Boost lives in the common library from the first commit because all three
vehicle classes use it.

### 2.1 Model injection — the extensibility mechanism

The equations of motion never name a model. They call handles carried in an
environment struct:

```matlab
env.atmos = @coorbital.atmos.expAtmos;    % → us76Atmos → nrlmsise
env.grav  = @coorbital.grav.sphereGrav;   % → j2Grav → j24Grav
env.aero  = @coorbital.aero.constLD;      % → table lookup
env.omegaE = 0;                           % → 7.292115e-5 rad/s to spin the Earth
```

Raising fidelity is a one-line change in an entry script. The EOM, integrator,
guidance, and plotting are untouched. Each model family shares a fixed signature so
the swap is total:

| Family | Signature |
|---|---|
| atmosphere | `[rho,P,T,a] = f(h)` — `h` geometric altitude (m), SI out |
| gravity | `[gr,gLat] = f(r,lat)` — radial and north components (m/s²), `[N x 1]` each |
| aerodynamics | `[CL,CD] = f(alpha,mach,vehicle)` |

The gravity signature returns **two** components deliberately. A scalar "local up"
value cannot represent J2, whose acceleration has a latitudinal term; committing to a
scalar now would force every EOM call site to change when `j2Grav` arrives.
`sphereGrav` simply returns `gLat = 0`.

The rotating-Earth Coriolis and centrifugal terms are **written into `glide3DOF` from
the start**, multiplied by `env.omegaE`. Setting it to zero recovers the non-rotating
case exactly; turning rotation on is a parameter change, not an edit. Same approach for
oblateness: `j2Grav` returns the same shape as `sphereGrav`.

### 2.2 Phase machinery

A phase is a struct:

```matlab
phase.eom       = @coorbital.eom.glide3DOF;
phase.guide     = @(t,x) coorbital.guide.prescribed(t,x,schedule);
phase.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,0);
phase.tspan     = [0 3000];
```

`coorbital.prop.phaseRun(phases,x0,vehicle,env)` integrates each in turn, using ODE
event functions to detect the handoff (burnout, apogee, altitude floor, impact) and
passing the terminal state of one phase as the initial state of the next. Phases may
use different EOM and different state definitions; `stateConvert` handles the mapping
at each junction, and the junction states are recorded so the optimizer in Phase 2 can
enforce them as linkage constraints.

---

## 3. What comes from pumpkyn

Reused directly, not rewritten:

| Need | pumpkyn routine |
|---|---|
| Geodetic ↔ ECEF ↔ ECI | `LLA2ECEF`, `ECEF2LLA`, `ECItoECEF`, `ECEFtoECI`, `LLA2ECI`, `ECI2LLA` |
| Time and sidereal angle | `juliandate`, `JD2GMST`, `JD2GAST` |
| Vector math | `vmag`, `bsxDot`, `bsxAng`, `multiplyDCM`, `Rx`, `Rz` |
| Physical constants | `getConst` |
| Earth and sky rendering | `earth3D`, `sphere3D`, `stars3D`, `addFigureLogo` |

**Not present in pumpkyn and therefore ours to write:** atmosphere, aerodynamics,
gravity harmonics (there is no J2 anywhere in the library — `getConst` has none),
atmospheric-flight equations of motion, a general integrator (`cr3bp.prop` is
CR3BP-specific), and all guidance.

Two standing cautions: **do not modify `external/pumpkyn`**, and do not consume
`getConst`'s known-bad `deg2ArcSec` or `c` fields.

---

## 4. Glide dynamics

The 6-state 3-DOF model, matching hyperFLIGHT's stated state count. Over a spherical
Earth with `omegaE = 0`:

```
state x = [r, lon, lat, V, gamma, psi]
control u = [alpha, sigma]        angle of attack, bank angle

dr/dt     =  V sin(gamma)
dlon/dt   =  V cos(gamma) sin(psi) / (r cos(lat))
dlat/dt   =  V cos(gamma) cos(psi) / r
dV/dt     = -D/m - g sin(gamma)
dgamma/dt =  (L cos(sigma))/(m V) - (g/V - V/r) cos(gamma)
dpsi/dt   =  (L sin(sigma))/(m V cos(gamma)) + (V/r) cos(gamma) sin(psi) tan(lat)

L = 0.5 rho V^2 S CL,    D = 0.5 rho V^2 S CD
```

`g` above is the radial component `gr`. The latitudinal component `gLat` is zero for
`sphereGrav` and enters `dV/dt`, `dgamma/dt`, and `dpsi/dt` through its projection onto
the velocity frame once `j2Grav` is in use — the terms are written now and multiply out
to zero under spherical gravity.

Rotating-Earth terms add to `dV/dt`, `dgamma/dt`, and `dpsi/dt` proportional to
`omegaE` and `omegaE^2`; they are present in the file, gated on `env.omegaE`.

Derived quantities reported but not integrated: Mach, dynamic pressure, load factor,
and stagnation-point heating rate (Sutton–Graves). These are the constraint metrics
Phase 2 will bound.

---

## 5. Validation

Analytic cases first, physical plausibility second. This is the reason for starting
simple: a spherical, non-rotating Earth with an exponential atmosphere has closed-form
solutions to check against, while an oblate rotating Earth with tabulated aero has
none.

| Check | Reference | Tolerance |
|---|---|---|
| Equilibrium glide: `V(h)` for constant L/D, small gamma | closed form | 1% over the glide |
| Allen–Eggers ballistic entry: peak deceleration and its altitude | closed form | 2% |
| Zero-lift vacuum: range and apogee of a ballistic arc | Keplerian | 1e-6 relative |
| Energy bookkeeping: work done by drag vs. mechanical energy lost | self-consistency | 1e-8 relative |
| Great-circle range, zero bank | spherical trigonometry | 1e-6 relative |

The vacuum and energy checks are exact and must hold to solver tolerance; they catch
sign errors and frame mistakes that the approximate checks would hide. Each becomes a
test in `tests/` with the expected value written into the file, following the verified
checkpoint pattern used in the tutorials.

---

## 6. Build order

1. `getConst` wrapper, `expAtmos`, `sphereGrav`, `constLD`, `vehicleDefaults`
2. `glide3DOF` — spherical, non-rotating; rotating terms present but gated off
3. `phaseRun` with altitude and apogee events
4. `HGV/run_glide.m` — prescribed constant CL and bank
5. Validation tests: vacuum, energy, equilibrium glide, Allen–Eggers
6. `+viz` — trajectory, ground track, altitude/velocity/Mach profiles
7. `boost3DOF` with a prescribed pitch program; `BM/run_ballistic.m`
8. Terminal/descent phase; full boost → glide → descent chain in one entry script
9. Fidelity increments, each behind its own validation: rotating Earth → J2 → oblate
   geodetic altitude → tabulated aero → US76 atmosphere

Items 1–5 are the first milestone: a validated glide propagator. Nothing after item 5
is started until those tests pass.

---

## 7. Conventions

MATLAB in pumpkyn house style, per `proj7/doc/pumpkyn_style_guide.md`:
`%%`-delimited header quartet with aligned input/output columns and units, `=` signs
column-aligned with right-justified names, colon-terminated `%%` step comments,
`if nargin == 0` self-demo in every library function calling itself by full namespace,
`[N x 3 x M]` N-last layout with explicit `dim`, norms via `vmag` and never `norm`,
and no `i`/`j` as loop variables.

Every function carries a `%% References:` block when its math has a source — Vinh for
the glide EOM, Allen–Eggers for the entry solution, Sutton–Graves for heating.

Entry scripts open with a fenced `%% USER PARAMETERS:` block holding launch site,
target, vehicle properties, control schedule, and model selections. Below that block,
nothing should need editing for a routine run.

---

## 8. Documentation

| File | Contents |
|---|---|
| `docs/DESIGN.md` | this document |
| `docs/README.md` | code organization, how to run, how to add a model |
| `docs/LESSONS_LEARNED.md` | running log — what broke, what fixed it |
| `docs/hgv_dynamics_note.tex` | frames, EOM derivation, atmosphere, aero, heating |
| `docs/software_design.tex` | architecture, interfaces, extension points, phase machinery |

The two LaTeX documents are written after the first milestone, when the interfaces have
survived contact with working code. `LESSONS_LEARNED.md` starts on day one.

---

## 9. Open items

- **Package name.** `+coorbital` matches `proj7/pumpkyn_style/+coorbital`. Alternative
  `+hyperfly` echoes hyperFLIGHT at the cost of consistency across repos.
- **Repository.** `missiles/` sits inside the `optimal_control` git repo, which has no
  remote configured. It is deliberately outside the `myLatex` tree, so none of this
  syncs to GitHub.
- **Vehicle parameters.** No representative HGV mass, reference area, or L/D has been
  chosen. Phase 1 uses plausible open-literature values, flagged as placeholders in
  `vehicle_hgv.m` until real numbers are supplied.
