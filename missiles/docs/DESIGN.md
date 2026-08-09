# Missile Trajectory Library — Design Spec

*Coorbital, Inc. — 2026-08-06*

> ## ⚠ READ THIS FIRST — this document is a pre-implementation spec
>
> Everything below §0 was written **before** the code existed and has **not**
> been revised to match it. Twenty-five specific statements in it are now false
> — fifteen catalogued when §10 shipped, six more in §11.3, and four in §12.4:
> files that were never written, tolerances that were relaxed or replaced,
> validation checks that were delivered as something else, and a fidelity
> increment the spec called a one-line change that turned out to be one line in
> the propagator and a new solver in the targeting layer. Taking any of it as
> a description of the shipped library will mislead you.
>
> **`README.md` is the authority on delivered behaviour.** This file is kept
> unedited as the record of what was *intended*, which is worth having — the
> gap between the two is the most interesting thing in the repo.
>
> Every known discrepancy is listed in the dated as-built sections at the end,
> **one per milestone**, each written when that milestone shipped and never
> edited afterwards — the record of what was true at each stage is the point.
> Read them in order alongside any part of the spec you rely on:
>
> - **[§10 As built — 2026-08-06](#10-as-built--2026-08-06)** — the glide propagator
> - **[§11 As built — 2026-08-07](#11-as-built--2026-08-07)** — boost, descent, and the full chain
> - **[§12 As built — 2026-08-07](#12-as-built--2026-08-07)** — targeting and visualization
> - **[§13 As built — 2026-08-08](#13-as-built--2026-08-08)** — two-axis targeting
>
> **The later section is always the later record.** Several §10 entries
> ("boost is out of scope", "`stateConvert` will be written") were true when
> written and are now superseded; §11 says so at each point, and §12 and §13 do
> the same for their predecessors. **§12 and §13 were written together on
> 2026-08-09**, after their milestones rather than on the day, and they say so
> at their heads — they are reconstructions from the code and the commit record,
> not same-day records like §10 and §11.

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
| §2 tree: `+viz/` with `trajPlot.m`, `groundTrack.m`, `profilePlot.m` | **Delivered 2026-08-07**, after being skipped in two milestones. `groundTrack.m` and `profilePlot.m` are as specified; `trajPlot.m` was not written and `globe3D.m` was written instead, a static 3-D Earth with the arc over it being the picture a flat lat/lon grid cannot give at intercontinental range. **`globeMovie.m` was delivered alongside it on the same date** and this table omitted it until 2026-08-07: an MP4 of the same scene with the track GROWING frame by frame, the vehicle marked at its current position, a per-frame time/altitude readout and a phase legend, degrading to a plain grey sphere on a black background when pumpkyn's texture and starfield are not on the path. Two package-private helpers came with the package: `private/vizParent.m`, the only place in the package allowed to call `figure()`, and `private/earthSurface.m`, shared between `globe3D` and `globeMovie` so the still and the animation cannot end up rendering two different Earths. All four figure functions take the same `(traj,veh,env,opts)` arguments — `globe3D` and `globeMovie` read neither `veh` nor `env`, and say so in their headers. The three inline `figure` blocks in each entry script are gone. |
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

---

## 11. As built — 2026-08-07

**§10 above is the record of the glide milestone and is left unedited**, even
where this milestone has since overtaken it. This section is the correction list
for the second milestone: powered boost, a descent phase, and the full
boost → glide → descent chain. Same rules as §10 — it does not revise the spec,
it says where the spec and the code disagree and which one won, and `README.md`
remains the authority on delivered behaviour.

Every entry was checked against the code on 2026-08-07.

### 11.1 §10 entries this milestone supersedes

| §10 said, correctly at the time | As of 2026-08-07 |
|---|---|
| §10.1: "`boost3DOF.m`, `ballistic3DOF.m` — not delivered" | **`boost3DOF.m` is delivered.** `ballistic3DOF.m` was never written and will not be: a ballistic arc is `glide3DOF` at zero lift, so a third EOM file would have duplicated six equations to switch off two terms. `BM/vehicle_bm.m` carries `CL = 0.005`, `L/D = 0.02` instead, and `run_ballistic`'s cross-check forces `CL = 0` through a wrapper on the *aero handle* rather than through a separate EOM. |
| §10.1: "`eventApogee.m`, `eventBurnout.m` — not delivered" | **Both delivered**, both terminal and one-sided. See §11.3 for the one interface change. |
| §10.1: "`BM/` — neither folder exists. `run_tests.m` already looks for `BM` and `ICBM` on the path and skips them silently, so the suite is ready for them." | **`BM/` exists** with `run_ballistic.m` and `vehicle_bm.m`, and the suite picked it up exactly as predicted. `ICBM/` still does not exist. |
| §10.3: "Boost is out of scope for this milestone, so it lives nowhere. The reasoning still holds for when it arrives." | Boost now lives in `+coorbital/+eom/boost3DOF.m` in the common library, as §2 intended. Both `BM/` and `HGV/` use it. |
| §10.3: "`stateConvert` does not exist... A boost phase carrying mass in a Cartesian state is what will force `stateConvert` to be written." | **The prediction was wrong and `stateConvert` will not be written.** See §11.2 — this is the milestone's one genuine design change. |
| §10.2: "`+util/` contains `missileConst.m`, `vehicleDefaults.m` and `greatCircle.m`" | Add **`boosterDefaults.m`**. Still no `stateConvert.m`. |
| §10.3: "§7: '`if nargin == 0` self-demo in every library function' — **Eight of the ten** library functions have one. **`missileConst` and `vehicleDefaults` do not**… This is a ruled exception, not an oversight." | **Both counts have moved and there is now a THIRD exception.** Counted against the code on 2026-08-07: **17 library functions under `+coorbital/`, 14 of them with a self-demo, and three ruled exceptions — `missileConst`, `vehicleDefaults` and now `boosterDefaults`.** The reasoning is unchanged and now covers all three: each is a pure constant return whose entire content is visible in its header, so a demo would print back what the reader is already looking at. `README.md` ("Conventions") states this correctly; §10's "eight of ten, two exceptions" was the true count on 2026-08-06 and is left standing as that dated record, not corrected in place. |

Unchanged from §10 and still true: there is no
`hgv_dynamics_note.tex` or `software_design.tex`, no `ICBM/` folder, no
Sutton–Graves heating anywhere, and `missileConst` is still deliberately not a
wrapper around pumpkyn's `getConst`. **No longer true as of 2026-08-07: the
`+viz/` package now exists** — see the corrected row in §10.1.

### 11.2 The design change: `stateConvert` was replaced, not deferred

**Record this as a change to the design, not as an unbuilt promise.**

§2.2 specifies: *"Phases may use different EOM and different state definitions;
`stateConvert` handles the mapping at each junction."* That is not what was
built, and it is not what should have been built.

`phaseRun` concatenates all phases into a single `traj.x`, so every phase in one
call must share one state dimension. Rather than translating between a six-state
glide vector and a seven-state boost vector at each seam, **the whole chain runs
seven-state and the unpowered EOM is lifted to match**:

```matlab
eomFree = coorbital.eom.massConstant(@coorbital.eom.glide3DOF);   % appends dm/dt = 0
```

Why this is better than the spec's mapping layer:

- **A mapping at a junction is a place for state to be silently lost or
  reordered.** A uniform vector has no such place. There is nothing to get wrong
  at a seam because nothing happens at a seam.
- **`traj.x` stays a rectangular matrix** with one meaning per column across the
  whole flight. Plotting mass against time, or searching a peak across phases,
  needs no phase-aware indexing.
- **The optimizer story is unchanged.** §2.2's real requirement was that junction
  states be recorded for Phase 2 linkage constraints. They are, and they are now
  recorded in a single consistent coordinate system.
- **`stateConvert` would have had almost nothing to do.** The only genuine
  discontinuity in these chains is staging, which changes *one component* of an
  otherwise-continuous state. That is served by `phaseRun`'s optional per-phase
  `link` handle, `xNext = link(xEnd)`, in one line:
  `@(x) [x(1:6); x(7) - bst.massDry]`.

**The cost of the change, which is real and must be understood.** `glide3DOF`
divides by `veh.mass`, not by `x(7)`. So in an unpowered phase the carried mass
is **inert**: it rides along and drives nothing. Two sources of truth for one
physical quantity, with nothing in the arithmetic making them agree. A chain that
jettisons a stage while continuing to hand the same `veh` struct to the next
phase flies silently at the old weight, with a mass jump visible in `traj.x` that
the dynamics ignore. Measured: deleting the separation link in `run_ballistic`
leaves the flown trajectory bit-identical.

`massConstant` therefore raises `coorbital:massConstant:massMismatch` when `x(7)`
and `veh.mass` disagree. **The guard is not the fix; the per-phase vehicle is.**
Any chain with a mass-changing link must bind its own vehicle inside each EOM
closure, rebuilding mass **and** `Sref`, `CL`, `LD` — separation changes the whole
airframe, not just the weight. Both entry scripts do this. Full treatment in
`README.md`, "The seven-state convention and the inert-mass trap", and the
general lesson in `LESSONS_LEARNED.md`.

The spec's mapping layer remains the right answer if a future phase genuinely
needs a *different coordinate system* — a Cartesian ECI terminal phase, say.
Nothing here forecloses it; it simply was not needed for mass.

### 11.3 Other statements about behaviour that are now false

| Spec says | As built |
|---|---|
| §2 tree places `eventApogee.m` and `eventBurnout.m` alongside `eventAltitude.m` with no further detail | Delivered, but note `eventBurnout(t,x,mBurnout)` takes the burnout mass as an **explicit scalar argument**, not reached for inside the event from a booster struct. It must be `veh.mass + bst.massDry`; comparing `x(7)` against `bst.massDry` alone would treat the payload as part of the propellant budget and burn straight through it. |
| §2.1 lists three model families — atmosphere, gravity, aerodynamics | There are now **four**. Propulsion is `[T,mdot] = f(t,P,veh)`, taking ambient **pressure**, not altitude, because the back-pressure debit is a pressure term and passing `P` keeps propulsion independent of which atmosphere is installed. Read only by `boost3DOF`; `env.prop` is absent from unpowered runs and `glide3DOF` ignores it if present. |
| §2.2's phase struct has four fields: `eom`, `guide`, `terminate`, `tspan` | A fifth, optional, was added: **`link`**, a handle `xNext = link(xEnd)` applied to a phase's terminal state to form the next phase's initial state. Absent or `[]` is the identity (`[]` being how a MATLAB struct array mixes linked and unlinked phases). Applied **before** the junction is recorded, so `traj.junction(k).x` is the far side of a staging jump. |
| §4: control is `u = [alpha, sigma]`, and §7 implies guidance is schedule evaluation | Still `[alpha, sigma]`, but `guide/pitchProgram` is the first law that **reads the state**: it interpolates a commanded pitch *attitude* and returns `alpha = theta(t) - gamma`, taking `gamma` from `x(5)`. The `u = f(t,x,sched)` signature that §2.2's `prescribed` already had is what made this a drop-in. |
| §6 build order, item 7: "`boost3DOF` with a prescribed pitch program; `BM/run_ballistic.m`" and item 8: "Terminal/descent phase; full boost → glide → descent chain in one entry script" | **Both delivered**, in that order, and the chain script is `HGV/run_boost_glide.m`. Item 6 (`+viz`) was skipped at the time; it sat between them in the build order and was not a prerequisite for either. It was picked up on 2026-08-07 and the entry scripts were retrofitted onto it, their headline numbers unmoved. |
| §1 capability table, "Terminal: prescribed CL and bank" | Delivered as **prescribed bank only**. `constLD` ignores `alpha`, so the descent phase is flown by rolling to a large bank angle — 75° in the shipped block, leaving `cos(75°) = 0.2588` of the lift supporting the vehicle. Commanding `CL` directly awaits an alpha-dependent aero model, at which point the existing `descAlpha` entry becomes live with no code change. |

### 11.4 One thing built that the spec did not ask for

**`run_boost_glide` detects a bifurcation in its own user block.** The
glide descends in a damped skip phugoid, so the `hHandoff` altitude — advertised
as freely settable — terminates the glide a whole skip early if placed above a
phugoid trough. Measured: 7663.05 km of range at the shipped 15 km against
5781.57 km at 30 km, a loss of 1881 km, **with every phase reporting nominal**,
because every event fired exactly as asked. No termination diagnosis can see it.

The script now flies the handoff state on the *continued glide schedule* — a
propagation it needs anyway, to measure what the descent phase buys — and warns
when the vehicle would have climbed back above the handoff altitude. The check is
**one-sided and its blind spot is printed with it**: it cannot see a handoff that
truncated a shallow skip without rebounding, and 25 km is exactly that case.

This has no counterpart in the spec, which treats phase boundaries as given. It
is recorded here because it is the second time this library has found a
*silently wrong but internally consistent* result — the first being the
shared-constant blindness of §10 — and both were found by asking what the
existing checks could not see.

### 11.5 Open items (§9), as they now stand

| §9 item | Status |
|---|---|
| **Package name** | Resolved 2026-08-06: `+coorbital`. Unchanged. |
| **Repository** | §10.5 corrected the spec's claim that nothing syncs to GitHub: `optimal_control` **does** have a remote, `git@github.com:mcasey90277/optimal_control.git`. Still true, and now more consequential — this milestone adds booster performance parameters. They are all marked open-literature placeholders today, but the assumption written into the spec is wrong and must not be relied on when real numbers arrive. |
| **Vehicle parameters** | **Still open, and now larger.** `boosterDefaults` adds **eight** more PLACEHOLDER values — the five mass and propulsion fields (`massDry`, `massProp`, `thrustVac`, `Isp`, `Aexit`) **and the boosted-stack aerodynamic triple `Sref`, `CL`, `LD`**, which is not bookkeeping: it sets the stack's drag through the whole ascent and, with `separation = false`, the airframe for the entire unpowered flight. `vehicle_bm` adds a fourth vehicle parameter set. **Pinning, corrected 2026-08-07:** the five propulsion and mass fields are pinned exactly at `test_constThrust.m:50-56` (seven assertions, counting the 900 kg payload and `g0`); the aerodynamic triple was **not** pinned anywhere when this section was first written and now is, in its own block immediately below them, together with the derived `CD = 0.20`. `vehicleDefaults`' `Sref`, `CL` and `LD` are pinned in `test_runGlide`. Retuning any of these now fails loudly and says why rather than quietly changing every trajectory in the documentation. Both chain scripts print `(PLACEHOLDER values)` in their own summaries so the caveat travels with the output. |

---

## 12. As built — 2026-08-07

**§10 and §11 above are the records of the first two milestones and are left
unedited.** This section is the correction list for the third: point-to-point
**targeting** and the **visualization** package. Same rules as before — it does
not revise the spec, it says where the spec and the code disagree and which one
won, and `README.md` remains the authority on delivered behaviour.

**Written 2026-08-09, not on the day.** §10 and §11 were each written when their
milestone shipped. This one was not, and that is a difference worth declaring:
everything below was reconstructed on 2026-08-09 from the committed code and the
commit record, so it has the accuracy of a code reading rather than of a
same-day memory. Every entry was checked against the code on that date. Its
milestone spans `2026-08-07` through `2026-08-08` and it is dated by the
milestone, not by the writing.

### 12.1 §10 and §11 entries this milestone supersedes

| Said, correctly at the time | As of this milestone |
|---|---|
| §10.1: "`+viz/` with `trajPlot.m`, `groundTrack.m`, `profilePlot.m` — not delivered" | This one was already corrected **in place** in §10.1 on 2026-08-07, which is the single exception this document has made to its never-edit rule, and it is flagged there. `+viz` ships `groundTrack`, `profilePlot`, `globe3D` and `globeMovie`, plus three package-private helpers. `trajPlot.m` was never written. |
| §10.1: "`BM/` — `run_ballistic.m`, `vehicle_bm.m`" (§11.1: "`BM/` exists") | `BM/` gains a second entry script, **`run_ballistic_target.m`**. `ICBM/` still does not exist and the case for it is now weak: `run_ballistic_target` *is* the ICBM case, driven from a user block. |
| §11.1: "`+util/` — add `boosterDefaults.m`. Still no `stateConvert.m`" | Add **`greatCircleBearing.m`** and **`rangeSolve.m`**. Still no `stateConvert.m`, and §11.2's argument that it will never be written is unchanged. |
| §11.1: "17 library functions under `+coorbital/`, 14 of them with a self-demo, three ruled exceptions" | The count moved again with `+viz` and the two `+util` additions. The three ruled exceptions — `missileConst`, `vehicleDefaults`, `boosterDefaults` — are unchanged and the reasoning is unchanged. §11.1's numbers are left standing as that dated record. |

### 12.2 Files that exist and the spec's tree does not contain at all

The §2 tree has no targeting script of any kind. It lists `run_glide.m`,
`run_ballistic.m` and `run_icbm.m` — three scripts that fly *launch site +
azimuth + schedule → wherever the physics puts it*. Inverting that was not in
the spec, and two scripts now do it:

- **`HGV/run_target.m`** — launch point to destination for the boost-glide
  vehicle, ranging on **thrust-termination time**, chosen because range is
  monotonic and single-valued in it and bisection is therefore safe.
- **`BM/run_ballistic_target.m`** — the same problem for the ballistic vehicle,
  ranging on the **loft angle**, which is *not* monotonic. Range rises to a
  maximum and falls away on both sides, so every reachable range short of the
  maximum is flown by two arcs and the script has to bracket and **certify** the
  maximum before it can bisect either branch. It carries a
  `'minimum-energy' | 'lofted' | 'depressed'` selector, and it **measures**
  which branch it flew from the flown loft angle against the certified maximiser
  interval rather than trusting its own label.

Two library functions came with them, both in `+util`: `greatCircleBearing`,
the initial course of a great-circle arc, which raises `:degenerateArc` and
`:polarOrigin` rather than returning a plausible azimuth for a geometry that
has none; and `rangeSolve`, the scalar range bisection, which **refuses by
returning** — `converged = false`, the best-so-far point, and the bracket it
searched — rather than throwing or handing back a nearest miss.

### 12.3 The refusal contract, which the spec has no concept of

This is the milestone's one genuine addition to the architecture, and it is
worth recording as a design decision rather than as a feature.

A targeting script can be asked for something it cannot do: a target beyond the
reachable envelope, a target too close, a branch that does not reach. The spec
treats every entry script as a propagation that either runs or errors. Neither
is right here. Throwing loses the diagnostic — the caller wanted to know *how
far* out of reach the request was — and returning the nearest trajectory is
worse, because **a near miss looks exactly like a solution**.

What was built instead: an empty `traj`, `info.refused = true`, an
`info.refusedWhy` naming the gate, and the whole measurement record in `info`
regardless, so the caller can phrase its own refusal. The reachable envelope
comes free with the bracketing the solve does anyway, which is why it can be
printed on every run, converged or refused.

### 12.4 Statements about behaviour the milestone made false

| Spec says | As built |
|---|---|
| §5 validation table has no targeting row | Delivered checks the spec does not contain: the max-range loft angle is **certified** by adaptive refinement to a single interior hump and returns an *interval*, not a point; the flown branch is re-derived from the flown angle against that interval; and `'minimum-energy'` is checked in the **vacuum equal-radius limit**, where the constrained minimisation must reproduce the classical γ\* and V\* — measured 1.583e−8 rad and 4.698e−16 relative. |
| §1 capability table, "Terminal: prescribed CL and bank" (§11.3 already narrowed this to bank only) | `run_target` ships `descBank = 0` while `run_boost_glide` ships 75°. At the time this milestone shipped that was a **targeting limitation** — a banked descent turned the track off the departure arc and the script could only warn about the resulting cross-range. §13 removes the limitation; the shipped zero is now a default rather than a constraint. |
| §2.1: "`env.omegaE = 0` → 7.292115e-5 rad/s to spin the Earth", and §6 build order item 9 lists "rotating Earth" as a one-line fidelity increment | True of the **propagator** and false of the **targeting layer**, which is the distinction this milestone discovered and §13 closed. Turning rotation on is indeed one line in `glide3DOF` and `boost3DOF`. It is not one line in a script that takes its launch azimuth from a closed form derived for a non-rotating sphere: at the end of this milestone both targeting scripts **refused** `earthSpin = true` outright rather than returning a 231.6 km (HGV) or 463.2 km (BM) miss and calling it a solution. |
| §9 "Vehicle parameters" | `alphaMax` moved out of the entry scripts and became a **vehicle** property, `BM/vehicle_bm.m`'s `alphaMaxDeg = 6`, read by both `BM` entry scripts. It had been 12° in one and 6° in the other for one airframe, and the 12 had been chosen to bring a demonstration target inside the depressed branch — a placeholder tuned to make a demo work, which is the failure mode the whole PLACEHOLDER discipline exists to prevent. |

### 12.5 One thing built that the spec did not ask for

**`run_target` renders a movie.** `globeMovie` writes an MP4 of the trajectory
developing over a 3-D Earth, with the track growing frame by frame, a per-frame
time and altitude readout, and a phase legend — degrading to a plain grey sphere
on black when pumpkyn's texture and starfield are not on the path, so it needs
no network and no third-party toolbox to produce a working file. The spec's §2
tree has `trajPlot.m` and stops there.

It is off by default because it is the expensive part of a run, and its
altitude exaggeration ships at **1** — true scale — with `'auto'` available as
an opt-in adaptive rule. That default was itself a correction: exaggeration had
been on by default, which is a rendering choice quietly making a physical claim.

---

## 13. As built — 2026-08-08

**§10, §11 and §12 above are left unedited.** This section is the correction
list for the fourth milestone: **two-axis targeting**. Written 2026-08-09, one
day late, on the same terms as §12 — a reconstruction from the code and the
commit record, checked against the code on that date.

The closed-loop guidance seam specified in the same plan brief
(`docs/plan_2026-08-08_aim_solve_and_closed_loop.md`, Task 4) is a separate
piece of work with its own write-up at `docs/closed_loop_guidance.md`, and is
not covered here.

### 13.1 The milestone in one sentence

Both targeting scripts stopped **computing** the launch azimuth and started
**solving** it, beside the range control they already solved, against both
components of the miss instead of one.

`coorbital.util.aimSolve` is the new library function: a two-dimensional damped
Newton over a finite-difference Jacobian, with the step taken from an SVD of
that Jacobian rather than by backslash. It is the two-axis sibling of
`rangeSolve` and it honours the same contract for the same reason — converge on
the **achieved residual** and never on the step, never throw merely for failing
to converge, and hand back an `info` rich enough for the caller to write its own
refusal. A clean iteration costs exactly three evaluations, so the cost model is
`1 + 3n` and each evaluation is one full trajectory propagation.

### 13.2 §12 entries this milestone supersedes

| §12 said, correctly at the time | As of 2026-08-08 |
|---|---|
| §12.4: "at the end of this milestone both targeting scripts **refused** `earthSpin = true` outright" | **Both fly it.** `run_target`: 231 551.628 m of seed miss down to **4.361 m** in 7 residual evaluations. `run_ballistic_target`: 463 211.19 m down to **52.461 m**, also 7. The refusal, the `'earthSpin'` value of `info.refusedWhy` and the assertions gating the azimuth prose on `env.omegaE == 0` are all deleted rather than left as dead code. |
| §12.4: "a banked descent turned the track off the departure arc and the script could only warn about the resulting cross-range" | **Cross-range is a residual now, driven to zero.** At `run_boost_glide`'s 75° terminal bank `run_target` used to converge on range and land 21 524.695 m away; it now arrives at **54.795 m** in 4 evaluations. `info.crossWarn` is gone from both scripts, and so is the arcsin cross-track relation that produced it — nothing in the library computes it any more. |
| §12.2: `run_ballistic_target` "ranging on the **loft angle**" | Still true, and now only of **stage one**. Stage two solves the azimuth beside the mode's own range control: the loft angle on the two full-burn branches, the **cutoff fraction** in `'minimum-energy'` mode, whose loft is held at the value the minimisation settled on. |
| §11.1's and §12.1's library-function counts | Add `+coorbital/+util/aimSolve.m`, with `tests/test_aimSolve.m` beside it. |

### 13.3 The design decision: one solver closes two gaps

Recorded as a design change, because the two gaps did not look like one gap.

Earth rotation and a banked descent were separate limitations with separate
entries in `TODO.md`. They are the same limitation. The closed-form bearing is
the answer under exactly two conditions — `omegaE = 0` **and** `sigma ≡ 0` —
because only then is the flown ground track the great circle it departed on, so
matching the *distance* matches the *point*. Drop either condition and the aim
direction becomes an unknown. One extra degree of freedom, solved against a
residual that was already being measured, closes both.

**The mechanism that had been given for rotation was wrong, and its number
outlived its retraction.** Both scripts, and this project's documentation, used
to say that a turning Earth carries the target east under the vehicle so the
vehicle must be aimed where the target *will be*. It does not. The integrated
state is planet-relative and the target is a fixed ground point: neither moves
in the rotating frame. What rotation does is **deflect the vehicle**, through
the Coriolis and centrifugal terms that have been in the equations of motion
since the first commit.

The distinction is testable and was tested. A carried target is a down-range
effect; a deflected vehicle is a cross-range one. Measured seed misses:

| Script, rotating | Down-range | Cross-range | Ratio |
|---|---|---|---|
| `HGV/run_target` | −5 858.645 m | +231 477.499 m | 39.5 : 1 |
| `BM/run_ballistic_target` | −17 764.166 m | +462 870.439 m | 26.1 : 1 |

Both scripts deleted the eastward-ground-sweep figure they used to print
(628 km on the rotating `run_target` case) because it is the scale of the
mechanism the sentence above it had already retracted. **A retracted
mechanism's number left standing under a corrected sentence is the harder error
to notice, because the prose reads right.**

### 13.4 What did not move, and why that is the point

| Case | Before | After |
|---|---|---|
| `run_target`, shipped non-rotating zero bank | 511.243 m | **511.243 m, bit-for-bit**, 1 evaluation |
| `run_ballistic_target`, shipped `'depressed'` | 457.270 m | **unchanged bit-for-bit**, 1 evaluation |
| `run_ballistic_target`, shipped `'minimum-energy'` | 39.009 m | **unchanged bit-for-bit**, 1 evaluation |
| `run_ballistic_target`, shipped `'lofted'` | 779.491 m seed | **0.365 m**, 4 evaluations |
| `HGV/run_glide`, `HGV/run_boost_glide` | 6986.82 km, 7663.05 km | unmoved in every digit |

`aimSolve` tests its tolerance before doing any work, so a seed already inside
tolerance converges at iteration zero and hands the closed form back untouched.
Three of the four shipped targeting cases are that case. The fourth,
`'lofted'`, had a 779.491 m seed miss against a 707.107 m per-component
tolerance and therefore moved.

**That is also this milestone's coverage gap, and it is structural rather than
a tolerance that could be tightened.** On a non-rotating zero-bank run the miss
*is* the down-range residual the stage-one bisection has already driven inside
the same tolerance, so the seed can never be outside it: measured across a
twentyfold tightening of `run_target`'s tolerance, the solve still exits at
iteration zero every time. The shipped run a user gets by typing `run_target`
therefore exercises none of the Jacobian, the SVD, the step cap, the line
search or the seven refusals. Coverage comes from the rotating and banked cases
and from `tests/test_aimSolve`, which drives the solver against synthetic
residuals with closed-form roots and integrates no equation of motion at all.

### 13.5 What the milestone opened

Solving on a **surface** with machinery certified on a **curve** has
consequences, and they are open items rather than defects. All are in `TODO.md`
with measured sizes; the two that bear on the architecture:

- **The reachable envelope is measured at the seed azimuth.** Every guarantee
  `run_ballistic_target`'s stage one makes — the hump, the certified maximiser
  interval, the branch brackets, the reachable band and therefore the too-far
  and too-close refusals — is a one-dimensional argument at one azimuth, and is
  indicative rather than exact for the azimuth stage two settles on. Measured:
  the maximum moves about 0.28° of loft per degree of azimuth against a
  certified interval 0.043° wide, and turning rotation on moves the envelope by
  30 to 42 km. The visible consequence is that `'depressed'` plus `earthSpin`
  is refused for the shipped target — correct physics on an easterly launch,
  and pinned by a test so a reader does not blame the new solve.
- **The stage-two Newton is not boxed to the certified branch side**, and that
  is a decision. The only interval available to draw a box on is the
  seed-azimuth one, and at 0.28° of loft per degree of azimuth against a
  0.043° interval a box drawn there would exclude genuine roots. Three layers
  of **diagnosis** stand instead of a prohibition: a per-iteration step cap, a
  measurement of the flown angle against the interval with a printed caution
  when they disagree, and outright refusal when the Newton cannot converge.

### 13.6 Open items (§9), as they now stand

| §9 item | Status |
|---|---|
| **Package name** | Resolved 2026-08-06: `+coorbital`. Unchanged. |
| **Repository** | Unchanged from §11.5, and unchanged in consequence: `optimal_control` has a remote, so the spec's assumption that none of this syncs to GitHub is still false. This milestone adds no new parameters, only a solver. |
| **Vehicle parameters** | **Still open, and now with a sharper example of why it matters.** §12.4 records `alphaMax` being moved into `BM/vehicle_bm.m` after it was found set to 12° in one script and 6° in another for one airframe, the 12 having been chosen to make a demonstration target reachable. The shipped 6° is a placeholder awaiting a qualification basis, not a cleared limit, and it now says so in the vehicle file, in both `BM` user blocks and in the printed limitations. |
