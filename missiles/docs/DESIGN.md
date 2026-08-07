# Missile Trajectory Library — Design Spec

*Coorbital, Inc. — 2026-08-06*

> ## ⚠ READ THIS FIRST — this document is a pre-implementation spec
>
> Everything below §0 was written **before** the code existed and has **not**
> been revised to match it. Fifteen specific statements in it are now false:
> files that were never written, tolerances that were relaxed or replaced,
> validation checks that were delivered as something else. Taking any of it as
> a description of the shipped library will mislead you.
>
> **`README.md` is the authority on delivered behaviour.** This file is kept
> unedited as the record of what was *intended*, which is worth having — the
> gap between the two is the most interesting thing in the repo.
>
> Every known discrepancy is listed in **[§10 As built — 2026-08-06](#10-as-built--2026-08-06)**
> at the end. Read that section alongside any part of the spec you rely on.

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

---

## 10. As built — 2026-08-06

**Everything above this line is the pre-implementation spec, unedited.** This
section is the correction list, written after the first milestone shipped. It
does not revise the spec; it says where the spec and the code disagree and
which one won. `README.md` is the authority on delivered behaviour — where this
section and `README.md` differ, `README.md` is right and this section is stale.

Each entry was checked against the code on 2026-08-06, not carried over from a
review note.

### 10.1 Files the spec names that do not exist

| Spec says | As built |
|---|---|
| §0 preamble: "the math derivations live in `hgv_dynamics_note.tex` and the interface detail in `software_design.tex`" | **Neither file exists.** Both were deferred on purpose — see §8, which correctly says they are written *after* the interfaces survive contact with working code. The preamble states them as though they already do. The EOM sources are cited inline in `+coorbital/+eom/glide3DOF.m` instead. |
| §2 tree: `boost3DOF.m`, `ballistic3DOF.m` | Not delivered. Boost and ballistic are out of scope for this milestone. |
| §2 tree: `eventApogee.m`, `eventBurnout.m` | Not delivered. `eventAltitude.m` is the only event function. |
| §2 tree: `+viz/` with `trajPlot.m`, `groundTrack.m`, `profilePlot.m` | **No `+viz` package at all.** `run_glide` carries three inline `figure` blocks. |
| §2 tree: `+util/stateConvert.m` | Not delivered — see §10.3. |
| §2 tree: `BM/` (`run_ballistic.m`, `vehicle_bm.m`) and `ICBM/` (`run_icbm.m`, `vehicle_icbm.m`) | Neither folder exists. `HGV/` is the only vehicle folder. `tests/run_tests.m` already looks for `BM` and `ICBM` on the path and skips them silently, so the suite is ready for them. |

### 10.2 Files that exist and the spec's tree omits

The `+util/` line of the §2 tree reads `stateConvert.m vehicleDefaults.m`. As
built, `+util/` contains **`missileConst.m`**, **`vehicleDefaults.m`** and
**`greatCircle.m`** — and no `stateConvert.m`. The two additions are not
incidental:

- `missileConst.m` is the single source of truth for every physical constant,
  and `test_missileConst` is the load-bearing test of the whole suite. See the
  "structural blind spot" section of `README.md` for why: an analytic reference
  computed from the same constants as the model cannot detect a wrong constant.
- `greatCircle.m` is the haversine central angle, used for ground range.

### 10.3 Statements about behaviour that are false

| Spec says | As built |
|---|---|
| §2: "Boost lives in the common library from the first commit because all three vehicle classes use it." | **Boost is out of scope for this milestone**, so it lives nowhere. The reasoning still holds for when it arrives. |
| §2.2: "Phases may use different EOM and different state definitions; `stateConvert` handles the mapping at each junction" | `stateConvert` **does not exist**. `phaseRun` passes the 6-state vector straight through from one phase to the next, because every Phase 1 phase shares the same glide state. The junction states *are* recorded, as specified, so the Phase 2 linkage-constraint plan is unaffected. A boost phase carrying mass in a Cartesian state is what will force `stateConvert` to be written. |
| §4: "Derived quantities reported but not integrated: Mach, dynamic pressure, load factor, and stagnation-point heating rate (Sutton–Graves)." | Mach, dynamic pressure and load factor are computed and reported by `run_glide`. **Heating is not computed anywhere** — there is no Sutton–Graves code in the library. `vehicleDefaults` carries a `noseRadius` field for it, and now says at the point of definition that it is a placeholder for deferred work rather than dead code. |
| §7: "`if nargin == 0` self-demo in every library function" | Eight of the ten library functions have one. **`missileConst` and `vehicleDefaults` do not**, and deliberately: both are pure constant returns whose entire content is visible in the header, so a demo would print back what the reader is already looking at. This is a ruled exception, not an oversight. |
| §7: "Every function carries a `%% References:` block when its math has a source — Vinh for the glide EOM, Allen–Eggers for the entry solution, **Sutton–Graves for heating**." | The first two are delivered. There is no Sutton–Graves block because there is no heating code. |
| §6 build order, item 1: "`getConst` wrapper" | **`missileConst` is deliberately not a wrapper around pumpkyn's `getConst`.** `getConst` has no `muE`, no Earth radius, and two fields known to be wrong (`deg2ArcSec`, `c`). Wrapping it would have made the library's single source of truth depend on a struct that cannot supply its two most important members. `missileConst` is standalone and cites WGS-84 and US-76 directly. |

### 10.4 The validation table (§5) as delivered

Every row moved. The delivered budgets and the reasoning behind each are in the
"Validated results" table of `README.md`, which carries the measured numbers;
this is only the spec-to-code mapping.

| §5 row | Spec tolerance | As built |
|---|---|---|
| Equilibrium glide `V(h)` | 1 % over the glide | **3 %**, measured 1.51 %. The 1 % figure did not survive contact with a real arc: an open-loop constant-alpha glide phugoids about the equilibrium curve, and the quasi-steady assumption gives way at the low, slow end. 3 % is roughly twice the measured departure and still well below the 11.52 % a 20 % lift error produces. |
| Allen–Eggers ballistic entry: "peak deceleration **and its altitude**" | 2 % | Split into three delivered checks, none of them at 2 %. **(a)** Peak deceleration, budget **one-sided `0 < rel < 12 %`**, measured +10.50 % and +9.20 %. Not two-sided: the closed form drops gravity and freezes the flight path angle, and both neglected effects can only *raise* the peak, so a simulated peak *below* the analytic value means too much drag rather than a better approximation. **(b)** Peak *location*, budget **6 %** — and compared **in density**, not in altitude as the spec row implies, because the same discrepancy reads as 4.85 % and 4.36 % in density but 6.55 % and 1.52 % in altitude, so an altitude budget would mean different things at different ballistic coefficients. The spec gave the location no budget of its own. **(c)** A check the spec does not contain at all: **β-independence**, budget 3 %, an eightfold change in ballistic coefficient must move the peak by less than 3 % while moving its altitude by kilometres. That independence is Allen–Eggers' most distinctive claim and is now exercised directly. |
| "Zero-lift vacuum: **range and apogee** of a ballistic arc, Keplerian, 1e-6 relative" | 1e-6 | Delivered in `test_glide3DOF` as **conservation of specific energy `V²/2 − μ/r` over a 300 s lofted arc, at 1e-8** (measured 1.28e-13). **Range and apogee are never checked against a Keplerian reference.** Energy conservation is the sharper instrument for the sign and frame errors this row was meant to catch, and it needs no separate propagator to produce the reference; but it is a different claim, and anyone who wants the Keplerian range check should know it was not written. A rotating-frame Jacobi-integral analogue at the same 1e-8 was added alongside it, which the spec does not mention. |
| "Energy bookkeeping: work done by drag vs. mechanical energy lost, 1e-8 relative" | 1e-8, integrated over an arc | Delivered as an **instantaneous identity at 1e-12**: `dV/dt` from `glide3DOF` against `−D/m − g sin(gamma)` formed independently at a single state. Tighter, and it localises a fault to the equation rather than to somewhere along a trajectory, but it is a point check and not the integrated bookkeeping the row describes. |
| "Great-circle range, zero bank, 1e-6 relative" | 1e-6 | **Was** delivered only as a unit test of `greatCircle` itself; the *propagated* zero-bank arc was checked nowhere automated, and its range had been verified exactly once, by hand, by a reviewer. **As of 2026-08-06 this is closed**: `tests/test_runGlide.m` propagates the shipped zero-bank glide and asserts the reported range and central angle against pinned literals at 1e-4 relative, asserts the range against its own reported angle on the stated sphere, and asserts the invariance of both under moving the entry point and heading. Note the budget is 1e-4 against a *pinned measurement*, not 1e-6 against an independent closed form — spherical trigonometry gives the angle between two points, not the arc a vehicle flies, so there is no closed form to compare against here. |

### 10.5 Open items (§9), as resolved

| §9 item | Status |
|---|---|
| **Package name** — `+coorbital` or `+hyperfly` | **Resolved: `+coorbital`**, matching `proj7/pumpkyn_style/+coorbital`. |
| **Repository** — "`missiles/` sits inside the `optimal_control` git repo, **which has no remote configured**. It is deliberately outside the `myLatex` tree, so **none of this syncs to GitHub**." | **Resolved as to location** — `missiles/` is in-repo, in `optimal_control`, outside `myLatex`. **But the stated rationale is now false, and this one matters.** `optimal_control` *does* have a remote: `origin git@github.com:mcasey90277/optimal_control.git`. Being outside `myLatex` does not keep this off GitHub. Anyone pushing this branch is publishing the library. Nothing here is sensitive today — every vehicle number is a marked open-literature placeholder — but the assumption written into the spec is wrong and should not be relied on when real parameters arrive. |
| **Vehicle parameters** | **Still open, and the spec's text is still accurate.** Every field of `vehicleDefaults` and every override in `vehicle_hgv` is a marked PLACEHOLDER. The two files duplicate the same four values on purpose, so that real numbers have exactly one home when they arrive; `test_runGlide` now asserts they have not drifted apart. |

### 10.6 One thing built that the spec did not ask for

`HGV/run_glide.m` takes an optional struct of overrides for its
`%% USER PARAMETERS:` entries. It was added so `tests/test_runGlide.m` could
drive the entry script at a second operating point — necessary because the
shipped due-east equatorial geometry is provably blind to a lat/lon
transposition at the `greatCircle` call site. It also serves the batch
throughput requirement in §1. The fenced user block remains the interface for
a routine interactive run; see `README.md`.
