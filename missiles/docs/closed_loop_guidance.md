# Closed-Loop Boost Guidance — Design Spike

**Status:** decisions, measured. No guidance law is implemented and none is designed here.
**Date:** 2026-08-08, revised 2026-08-09 against review. **Author:** Michael Casey.
**Scope:** settles the four questions that gate PEG and VOA, by experiment. It is the brief the implementation milestone is written from.

The mathematics is `docs/hgv_dynamics_note.tex` §*Guidance* (`sec:guidance`, `sec:termcon`, `sec:voa`, `sec:peg`) and is not re-derived here — it is cited. The four gaps are `docs/software_design.tex` §*The closed-loop guidance seam* (`sec:seam`). This document converts those gaps into decisions and attaches numbers to them.

**What is shipped alongside this brief:** `+coorbital/+guide/terminalConstraint.m` and `tests/test_terminalConstraint.m`. Nothing else. No PEG, no VOA, no driver, no dispersion machinery.

---

## Summary of the two decisions

**Decision 1 — guidance cadence.** A **driver above `phaseRun`, one `phaseRun` call per guidance cycle**, holding the solved parameter vector in the driver's own workspace. `phaseRun` is not modified, the state vector is not widened, and the purity contract is honoured literally rather than worked around: within a cycle the guide is a **constant**.

This **is option (a)'s honest version — "derivative zero, updated only at phase boundaries" — with the carrier relocated** from the integrated state vector to the driver's workspace. It is not a third mechanism; it is the same mechanism with the guidance parameters kept out of `traj.x`. What the measurement forces is the narrow conclusion that **no guidance update can fire inside a call to `guide(t,x)`**, which both (a)-as-specified and this satisfy. Two structural facts, not the measurement, decide between them; they are in *The two candidate designs* below.

**Decision 2 — the thrust-direction inverse.** Resolve the commanded direction on the velocity triad and use the **`atan2` form**, not the `arccos` form written in the math note:

```
    alpha = atan2(hypot(p.eGamma, p.ePsi), p.vHat)
    sigma = atan2(p.ePsi, p.eGamma)
```

Round-trip through the library's own boost equations over 720 cases: worst `|dalpha|` **2.22e-16 rad** for the `atan2` form against **1.21e-11 rad** for `arccos`. **The deciding argument is scale robustness, not conditioning** — `atan2` is homogeneous of degree zero and needs no unit vector, while `arccos` silently requires one and returns 6.522 deg for a true 6 deg at `|p| = 0.999`, an 8.7 % error. See the correction under *Why `atan2`, corrected*.

---

## The strategic frame, restated so the implementation does not relitigate it

**PEG targets the ballistic boost phase.** Its five terminal constraints *are* a free-flight-ellipse burnout state, which is the problem linear-tangent steering was invented for. The map from ground target to required burnout state already exists here and is validated: the Keplerian free-flight range from a burnout state (`docs/hgv_dynamics_note.tex` §`sec:kepler`, agreeing to 1.5e-15), plus `BM/run_ballistic_target.m`'s `minimum-energy` mode, which minimises burnout specific energy subject to a range constraint. That mode's output is exactly a `(Rd, Vd, gammaD, wHat)` tuple, which is exactly what `terminalConstraint` consumes.

**The HGV is out of scope.** When a glide phase decides the range there is no closed-form map from burnout state to impact: the burnout conditions producing a given impact point are a family, not a point, and choosing a member of that family is a modelling decision with no defensible default. `HGV/run_target.m` already solves the HGV aim problem the honest way — propagate the whole flight and measure where it lands. Nothing in this brief improves on that, and PEG would not.

---

## Question 1 — statefulness and cadence

### The measurement

Instrumented guide recording `(callIndex, t, x, u)` on every evaluation, wrapped around `coorbital.guide.pitchProgram`, flown on the shipped `BM` boost case: launch 45 N / 100 W, azimuth 35 deg, `omegaE = 0`, `boosterDefaults` stack, `vehicle_bm` payload, `RelTol = AbsTol = 1e-10`, boost phase alone, burnout on the mass event at **t = 80.517757895 s**. The accepted trajectory was recovered separately by running the same right-hand side through a bare `ode45` returning a solution structure, so `deval` could give the state actually flown at any recorded call time. Scratch only; `restoredefaultpath` before every run.

| Quantity | Measured |
|---|---|
| Returned samples (`traj.t`) | 493 |
| Total guide calls, one boost phase | **1323** = 1 sizing + 829 solver-stage + 493 reconstruction |
| Calls per second of flight | **16.43** |
| Accepted `ode45` steps | 123 |
| Stage calls per accepted step | 6.74 |
| Accepted step size, min / median / max | 0.00496 / 0.606 / **1.502** s |
| Rejected attempts (adjacent decreases in `t`) | **15** |
| Stage calls at a time below the running maximum | **172 of 829 (20.7 %)** |
| Largest single backward jump in `t` | **0.9218 s** |
| Distinct call times among 829 stage calls | **681** — 148 calls repeat a time already used, **with a different state** |
| Reconstruction tail equals `traj.t` | exactly, max difference **0 s** |
| Reconstruction times also seen as stage times | **123 of 493** |
| Worst radius deviation between a state shown to the guide and the state flown at that instant | **2.119 m** |
| Worst speed deviation, worst flight-path-angle deviation | **0.1279 m/s**, **1.02e-4 rad** |
| Median radius deviation, interior-stage calls | 0.0122 m |
| Worst radius deviation on rejected-attempt calls | 0.322 m |

### What that says, item by item

**Times are not monotonic, and the excursion is the size of a guidance cycle.** Fifteen rejected attempts put 20.7 % of the stage calls at times the integrator had already passed, the worst by 0.92 s. A guide that latched "solve when `t` crosses the next cycle boundary" would fire during a rejected attempt, on a trajectory that was never flown, and would then hold that solution across the retry. The failure is silent: the trajectory that comes out is plausible.

**The guide is called twice at the same instant with two different states.** 829 stage calls at only 681 distinct times. Dormand–Prince has two stages at `c = 1`, and their states differ — worst case **1.399 m of radius apart, at `t = 55.878095 s`**. **A latch keyed on time cannot work**, because time does not uniquely identify a call, and the two states it cannot distinguish between are metres apart, not centimetres.

**States shown to the guide are not states the vehicle occupied — and the worst offenders are the *accepted* steps.** The largest deviation, 2.12 m of radius and 0.128 m/s of speed, is an interior Runge–Kutta stage of an accepted step, not a rejected one, because rejected attempts are precisely the ones the integrator retries with a smaller step. This kills the intuition that "rejected steps are the problem". Every interior stage is an intermediate approximation, and there are 5.7 of them per accepted step.

**The reconstruction pass is worse than documented.** `phaseRun`'s header says it calls the guide once per returned sample after every stage call. It does — 493 calls, at times bit-exactly equal to `traj.t`. But only **123 of those 493 times were ever evaluated during integration**: `ode45`'s default `Refine = 4` puts 4·123 + 1 = 493 output points on the grid, and 370 of them are interpolation times the integrator never asked the guide about. A stateful guide is not merely advanced by the reconstruction pass; it is advanced at times that never existed.

### The two candidate designs, evaluated

**(a) Guidance state rides in the integrated state vector.** The measurement rules out the *continuous* version outright: a discretely updated quantity would have to be updated inside a `guide` or EOM call, i.e. at one of the 829 stage evaluations, 20.7 % of which are on retried arcs and 148 of which duplicate a time. It does **not** rule out (a) as the plan specified it — zero derivative, updates only at phase boundaries. That version and Decision 1 are the same mechanism: both need one integration segment per cycle, and both solve only on accepted terminal states. What separates them is two structural facts about `phaseRun`, neither of which is a measurement:

- **`phaseRun`'s link signature is `xNext = link(xEnd)` — there is no `t`.** A per-cycle solve placed in a link cannot know what time it is, so the clock would have to be smuggled into the state vector alongside the guidance parameters. That is a second non-state riding in `traj.x` to support the first.
- **The cycle count is not known in advance, and `phaseRun` cannot stop a chain early.** Burnout arrives when the mass event fires, and a pre-built array of sub-phases would keep integrating past it — every remaining phase would run, on an exhausted motor, with `boost3DOF` dividing by a mass state the burnout event was supposed to stop at.

Both are fatal to (a) in this library, and neither is a cost that could be paid down. Two further costs are worth recording accurately because they are easy to overstate: the event functions index positionally and are **unaffected** by appending `x(8:12)`, and the `massConstant` objection is right in its conclusion but wrong in its mechanism — that file's guard is the permissive `numel(x) < 7`, and what actually breaks is the hard-coded state width at `massConstant.m:112`.

**(b) `phaseRun` gains a guidance cadence.** Correct semantics, wrong location. It changes the contract of the one file the whole suite rests on, it introduces a third clock beside the two that already coexist, and — the fatal detail — a held command is a discontinuity `ode45` must be stopped and restarted across, which means either an event per cycle or **a sub-phase per cycle**. The sub-phase already exists. `phaseRun` need not know.

### Decision 1: a driver above `phaseRun`, one call per cycle

```
    z    = seed
    xNow = x0
    tNow = 0
    loop
        z     = solve(z, tNow, xNow)          % PEG's five parameters, or VOA's seven
        uHold = command(z)                    % p(0; z*) mapped to (alpha,sigma)
        ph.guide = @(t,x) uHold               % a CONSTANT: pure by construction
        ph.tspan = [tNow, min(tNow+dt, tMax)]
        seg   = coorbital.prop.phaseRun(ph, xNow, veh, env)
        append seg; xNow = seg.x(end,:)'; tNow = tNow + elapsed
        stop on the burnout event, or when t_go* <= dt
```

The command for cycle *k* is not known until cycle *k*−1 has flown, so the phases cannot be built up front and this must be a loop of `phaseRun` calls, not one call with a pre-built array. That is not a limitation; it is what "closed loop" means.

Everything the measurement objects to disappears: the solve happens **outside** any integration, exactly once per cycle, on the accepted terminal state of the previous cycle — a state the vehicle demonstrably occupied. Within a cycle the guide is a constant, so rejected steps, duplicate times and reconstruction calls are all harmless, and `traj.u` is the control that was flown. `phaseRun` is untouched, and the purity contract is satisfied in the strongest possible sense.

### The cost of Decision 1, measured

Same boost case; the guidance solve is *not* included, so this isolates the cost of the mechanism. Wall times are **medians of five repeats** — single-shot timings on this workload are noisy enough to invert the ordering of the two arms, and did in an earlier draft of this table.

| Δt | Cycles | Held-command Δr | ΔV | Δγ | Wall, held | Wall, live |
|---|---|---|---|---|---|---|
| — (one phase) | 1 | — | — | — | 0.0273 s (1.00×) | — |
| 2.0 s | 41 | **+670.279 m** | −2.629 m/s | +1.095e-2 rad | 0.0473 s (1.73×) | 0.0951 s (3.49×) |
| 1.0 s | 81 | **+325.553 m** | −1.282 m/s | +5.301e-3 rad | 0.0886 s (3.25×) | 0.1784 s (6.54×) |
| 0.5 s | 162 | **+162.537 m** | −0.640 m/s | +2.641e-3 rad | 0.1741 s (6.39×) | 0.3486 s (12.79×) |
| 0.25 s | 323 | **+81.437 m** | — | — | 0.3483 s (12.78×) | 0.7047 s (25.85×) |

Three things fall out, and they must not be confused with one another.

**The overhead is per cycle, not a fixed multiple.** It is linear in the cycle count at roughly **0.5 to 1.0 ms per cycle** (0.488 / 0.757 / 0.907 / 0.994 ms at the four cadences), so quoting a single "3.4×" is a statement about Δt = 1 s and nothing else. The live-schedule arm is uniformly *dearer* than the held arm, which is the sensible ordering — it evaluates the interpolating guide at every stage while the held arm returns a constant.

**The restart artifact is bounded by the integration tolerance and independent of the cycle count.** Cutting the boost into cycles while letting the guide keep reading the live schedule moves the burnout radius by **−2.79e-4 m at 21 restarts and −1.25e-4 m at 806 restarts** — flat, not accumulating. It does scale with the tolerance: at `RelTol = 1e-7` the same two cases give **−0.0635 m and −0.106 m**. So the honest statement is not "nil" but *bounded by, and of the order of, the integration tolerance already accepted elsewhere in the library*, and 806 restarts is no worse than 21. That is the whole risk of Decision 1 and it is measured to be absent.

**The zero-order hold is not nil, and it is first order in Δt.** Holding the command across the cycle costs **326 m of burnout radius at Δt = 1 s**. The four points give `Δr/Δt` of 335.1, 325.6, 325.1 and 325.8 m/s — flat, so the error is clean first order and the 2 s point is the only one carrying any second-order contamination. This is not a numerical artifact; it is what holding a command costs, and a real flight computer pays it too. Its consequence for validation is in the task list below: **the nominal closed-loop run will not reproduce the open-loop trajectory, and must not be validated against it.**

Burn time is bit-identical (80.517757895 s) in every row, because burnout is a mass event and `constThrust`'s `mdot` does not depend on steering.

---

## Question 2 — the thrust-direction to (alpha, sigma) inverse

### The forward map, and the inverse

`+coorbital/+eom/boost3DOF.m` applies thrust as (`docs/hgv_dynamics_note.tex` eq. `athr`)

```
    a_thr = (T/m) [ cos(alpha) vHat + sin(alpha)( cos(sigma) eGamma + sin(sigma) ePsi ) ]
```

so with `p` the commanded unit direction resolved on the velocity triad, `pv = p·vHat`, `pg = p·eGamma`, `pp = p·ePsi`:

```
    alpha = atan2(hypot(pg, pp), pv)          % NOT acos(pv)
    sigma = atan2(pp, pg)
```

`docs/hgv_dynamics_note.tex` eq. `dirtoalpha` gives `alpha = arccos(p·vHat)`. It is analytically identical and the wrong form to implement; the reason is under *Why `atan2`, corrected* below, and it is **not** the one this document gave in its first edition.

The frame chain a PEG command has to travel, and which the implementation owns:

```
    p_inertial --(rotate by thetaG about z)--> p_ecef --(lon,lat)--> [n,e,u] --(gamma,psi)--> [vHat,eGamma,ePsi]
```

It is a pure rotation at every stage — a *direction*, not a velocity, so no `omega x r` term enters and inertial-versus-planet-relative does not arise for the command itself. The triad is built from the *planet-relative* velocity, which is the right thing: `alpha` is then the aerodynamic incidence the equations and `env.aero` expect. Verified: triad Gram matrix minus identity worst element **4.44e-16**, `vHat × eGamma − ePsi` worst **1.11e-16**, and round trip through `thetaG ∈ {0, 0.4, 2.9}` rad recovers `alpha` to **1.25e-16 rad** and `sigma` to **2.0e-15 rad**.

There is genuinely no `thetaG` anywhere in the library — but that **does not block this leg**, and an earlier edition of this brief wrongly said it did. Gravity is `sphereGrav`, spherically symmetric with no J2, no longitude dependence and no third body, and the target is specified in the same rotating frame the state lives in. Nothing in the model can distinguish one choice of `thetaG(0)` from another, so **`thetaG(0) = 0` is a free convention with zero modelling content**. Adopt it, write it down once, and move on.

### Numerical verification

Forward truth was taken from `boost3DOF` itself, by differencing the derivative with the motor on against the derivative with the motor off. Aerodynamics, gravity, kinematics and the rotation terms cancel exactly in that difference, leaving the thrust acceleration resolved on the triad:

```
    a.vHat = d(Vdot),   a.eGamma = d(gammadot)*V,   a.ePsi = d(psidot)*V*cos(gamma)
```

720 cases: 4 states (2–120 km altitude, 200–6000 m/s, gamma from 0 to 80 deg, four headings, mid-latitude) × 18 values of alpha (0 to ±179 deg) × 10 values of sigma (0 to ±180 deg).

| Result | Value |
|---|---|
| Worst alpha error, **`atan2` form** | **2.22e-16 rad** (1.27e-14 deg) |
| Worst alpha error, `arccos` form | 1.207e-11 rad (6.92e-10 deg) |
| Worst sigma error | 8.97e-13 rad (5.14e-11 deg) |
| With `env.omegaE` on | `dalpha` **0**, `dsigma` 1.11e-16 rad |

### Why `atan2`, corrected

Conditioning near zero incidence, exact analytic forward, so only the inverse is being measured:

| alpha (rad) | `arccos` error | `atan2` error |
|---|---|---|
| 1e-10 | **1e-10** (100 % relative) | 0 |
| 1e-8 | **1e-8** (100 % relative) | 1.65e-24 |
| 1e-6 | 4.45e-11 | 0 |
| 1e-4 | 2.62e-13 | 0 |
| 1e-2 | 1.44e-15 | 0 |
| 1e-1 | 5.55e-16 | 1.39e-17 |

`arccos(1 − eps)` has unbounded derivative at the identity, so the `arccos` form loses its digits as `alpha → 0`. **That is true and it is not the reason to prefer `atan2` here, and the first edition of this brief said it was.** The correction, measured on the shipped ballistic ascent: `|alpha| < 1 deg` for only **7.01 s** of the 80.52 s burn, while **74.1 %** of the burn sits pinned at the 6 deg clamp, where the `arccos` error is **4.4e-16 rad**. On this vehicle, on this trajectory, `arccos` would have been fine.

**The deciding argument is scale robustness.** `atan2(hypot(pg,pp), pv)` is homogeneous of degree zero: it reads the *direction* of `p` and is blind to its length. `arccos(pv)` silently requires `|p| = 1` and degrades — or saturates — the moment it does not get one. Measured on a true 6 deg direction:

| Direction magnitude | `arccos` returns | `atan2` returns |
|---|---|---|
| 1e-9 | 90.000000 deg | 6.000000 deg |
| 0.999 | **6.522475 deg** (8.71 % error) | 6.000000 deg |
| 1 | 6.000000 deg | 6.000000 deg |
| 1.001 | 5.427660 deg | 6.000000 deg |
| 1e+9 | 0.000000 deg | 6.000000 deg |

The `atan2` spread across `|p|` from 1e-9 to 1e+9 is **exactly zero**. This matters because the direction arriving at the mapping is the output of a guidance solve and a chain of frame rotations, not a vector anyone has renormalised — and PEG's slerp itself needs a renormalised small-angle branch (`docs/hgv_dynamics_note.tex` eq. `slerp`) precisely because unit norm is not automatic. A form that is wrong by 8.7 % when handed a vector one part in a thousand off unit length is a trap; a form that cannot notice is not.

### Where it degenerates

1. **The forward map is two-to-one.** `(alpha, sigma)` and `(−alpha, sigma ± pi)` produce the same direction. Any inverse must pick a branch; the `atan2` form returns `alpha ∈ [0, pi]`. **This conflicts with the existing convention** — `pitchProgram` returns *signed* `alpha`, and `BM/run_ballistic.m` flies negative incidence on parts of its schedule. The implementation must state its branch in the header and the consumer must not assume it matches `pitchProgram`'s.
2. **`sigma` is undefined at `alpha = 0`, and MATLAB returns 0 silently.** `atan2(0,0)` is `0`, so a direction exactly along the velocity comes back as `sigma = 0` with no warning: measured **40 deg of error** against the bank actually commanded. There is no incidence to roll, so any `sigma` is correct for the *direction* — but it is not correct for the *lift vector*, which `sigma` also orients. Guard: below a threshold incidence, hold the previous cycle's `sigma` and say so, or refuse. Do not return the atan2 default.
3. **`sigma`'s sensitivity is `|dp| / sin(alpha)`.** Measured, one-ulp perturbation of `p`: `sigma` moves by 0.0179 deg at `alpha = 1e-12`, 1.79e-5 deg at 1e-9, 1.79e-8 deg at 1e-6, 1.79e-11 deg at 1e-3. With a PEG direction converged to a tolerance `d`, `sigma` is meaningless below `alpha ≈ d`. Set the threshold in item 2 from the guidance solve tolerance, not from `eps`.
4. **A reversed direction returns a plausible number, not a refusal.** `p → −p` on a true 6 deg command returns **174.0000 deg** — well formed, finite, and wrong. `alpha = pi` is thrust reversed along the velocity, physically absurd on a boost, and the `sigma` ambiguity returns there. **The guard must be a band, not a point**: refuse above some `alphaMax`-derived ceiling, not only at exactly `pi`.
5. **A zero direction vector returns `alpha = 0, sigma = 0` silently.** Measured: `p = [0;0;0]` gives `atan2(0,0) = 0` twice, so a failed or uninitialised guidance solve produces a perfectly plausible "thrust along the velocity, wings level" command with no NaN and no error. **A `|p|` guard is needed alongside the thrust guard below** — scale invariance is a virtue everywhere except at zero, where it becomes a way of not noticing.
6. **Vanishing thrust.** The direction is undefined when `T = 0`. Measured with the delivered thrust scaled down: at scale 1e-6 the inverse still recovers 6.0000 deg; at 1e-12 it has drifted to **5.99965 deg**; at exactly zero it returns **NaN**. Guard on `T/m` against a physical threshold, and raise rather than return NaN — `boost3DOF` refuses non-finite controls anyway, but it would refuse *downstream* of the mistake.
7. **`alpha` near ±90 deg is NOT degenerate and needs no guard.** For the `atan2` form it is an ordinary point: measured exactly 90.0000000000 deg with `sigma` recovered to 0 rad of error. Recorded because it is the natural thing to add a guard for, and adding one would refuse a perfectly well-posed command. (For the `arccos` form 90 deg is where the conditioning is *best*, which is the opposite of where it is needed — another way of seeing that the two forms fail in different places.)
8. **Feasibility is not the dynamics' job — and the limit is not on the vehicle being steered.** A 3-DOF point mass will fly any direction it is handed; there is no attitude state to violate, so **the mapping must carry the feasibility check the dynamics cannot**. But the limit is awkwardly placed: `coorbital.util.boosterDefaults()` has **no `alphaMaxDeg` field at all** — its fields are `massDry, massProp, thrustVac, Isp, Aexit, Sref, CL, LD` — while the 6 deg limit lives on `BM/vehicle_bm.m`, which is the *separated re-entry body*. The boost phase flies the stack. `BM/run_ballistic.m` reads `veh.alphaMaxDeg` from the payload and applies it to the boosted vehicle, which is defensible as a placeholder and is not what a header would lead a reader to expect. **Task 2 must say which struct it reads**, and if the answer is the booster then the field has to be added there first.

   And there is no margin to spend. On the shipped ballistic run the boost angle of attack spans **−6.0000 to +0.1214 deg**, is **negative for 84.2 %** of samples, and sits pinned at the −6 deg clamp for **59.68 s of the 80.52 s burn (74.1 %)**. The nominal open-loop trajectory is already saturated at its stated control-authority limit for three-quarters of the burn, so **PEG has essentially zero incidence margin on this vehicle**: almost any correction it wants to command is one the clamp will refuse. Whether that is a real vehicle limit or an over-tight placeholder is a question for the vehicle file, not for the guidance law — but the dispersion campaign must record clamp activity per run, because on this configuration the clamp, not the guidance, may well be what determines the result.

### What the 3-DOF point-mass model costs here

There are no attitude dynamics, no moments and no rate limits. Orientation is an instantaneous input, so **a commanded direction is achieved instantaneously and exactly**. Three consequences, stated plainly:

- **The cycle boundary is a step change.** PEG's slerp (`docs/hgv_dynamics_note.tex` eq. `slerp`) exists to turn the thrust direction smoothly over the arc, with a quadratic warp chosen so the command leaves `p0` at zero rate. This library commands `p(0; z*)` and holds it, then jumps at the next boundary. The *shape* PEG is built around is discarded; only the endpoint it produces is used. That is what the algorithm specifies (`pegcmd`) and it is still a loss of the fidelity the shape was designed for.
- **No pitch-rate limit can be violated.** The shipped attitude schedule averages 0.67 deg/s (89 deg to 34 deg over 82 s), comfortably inside any plausible booster limit — so the *nominal* is not the concern. A dispersed PEG correction demanding a large re-point is, and this model will execute it in zero time and report success.
- **`sigma` owns the lift vector as well as the thrust.** Quantified at 6 deg of incidence on the boosted stack: `a_lift` is 0.0716 m/s² at 2 km, **0.2115 m/s² at 20 km** (the worst case), 0.0070 at 60 km and 4.9e-6 at 120 km, against `a_thrust_normal = 4.327 m/s²`. So the lift the bank angle drags around is at most **4.9 %** of the normal thrust it is being commanded for. Note also that **`coorbital.aero.constLD` ignores `alpha` entirely** (`CL = veh.CL` always), so incidence-induced lift is not modelled at all — the coupling is only the rotation of a constant-magnitude lift vector, and a 60 deg incidence command would produce no aerodynamic penalty whatsoever. The feasibility check of item 8 is therefore the *only* thing standing between PEG and a physically absurd command.

---

## Question 3 — the dispersion design (specified, not built)

Flown against truth on the nominal plant, PEG converges to the trajectory the open-loop pitch program already flies, and reporting that agreement as validation would be the reduction-test mistake of `docs/software_design.tex` §`sec:blind`. The demonstration is a **dispersed Monte Carlo: same dispersion realisations, open-loop arm versus PEG arm, compared on burnout-state scatter.**

### The mechanism

The `env` struct already injects `atmos`, `grav`, `aero` and `prop` by handle, so a perturbed model is a wrapper handle with the same signature. Nothing in the library changes.

```
    envD.atmos = @(h)        dispAtmos(h, atmosFn, k.rhoMul, k.dHscale)
    envD.prop  = @(t,P,veh)  dispProp(t, P, veh, propFn, k.thrustMul, k.ispMul)
    envD.aero  = @(a,M,veh)  dispAero(a, M, veh, aeroFn, k.clMul, k.ldMul)
```

Initial-state error perturbs `x0`; a thrust misalignment perturbs the *achieved* control and so wraps the **guide** handle, not `env`.

### What to disperse

| Parameter | Distribution (1σ) | Wraps | Why it belongs |
|---|---|---|---|
| Thrust magnitude multiplier | N(1, 0.01), truncated ±3σ | `env.prop` | Motor-to-motor scatter; changes total impulse and burn time. The headline plant dispersion. |
| Isp multiplier | N(1, 0.005), ±3σ | `env.prop` | `mdot = T/(Isp g0)`, so the two together set the burn time, which sets the burnout event. |
| Density multiplier | log-normal, `ln k ~ N(0, ln 1.15)` | `env.atmos` | ~15 % density scatter. **Must perturb `P` by the same factor** — the isothermal atmosphere ties them, and `constThrust` reads `P` for back-pressure, so perturbing `rho` alone is internally inconsistent. |
| Scale-height perturbation | N(0, 2 % of `c.Hscale`) | `env.atmos` | Makes the density error altitude-dependent rather than a rigid bias, which a rigid multiplier cannot represent. |
| `CL` and `L/D` multipliers | N(1, 0.10) each | `env.aero` | Aerodynamic coefficient uncertainty. Small during boost by the numbers above; not small on the descent. |
| Dry mass, propellant load | N(1, 0.005) each | `veh` / `bst` structs | Note `massConstant`'s mass-agreement guard: the perturbed mass must be applied to **both** `x0(7)` and the vehicle struct, or the chain refuses — a feature here, not an obstacle. |
| Pad position | 100 m 1σ each in down/cross, 10 m in altitude | `x0(1:3)` | Ignition-state error. |
| Pad velocity | 0.5 m/s in speed, 0.05 deg in gamma, 0.1 deg in azimuth | `x0(4:6)` | Same. |
| Thrust misalignment | N(0, 0.2 deg) constant bias in `alpha`, N(0, 0.2 deg) in `sigma` | the **guide** handle | The one dispersion PEG can see and correct and the open-loop program structurally cannot. Expect the largest arm-to-arm separation from this. |

### Seeding

One substream per run: `RandStream('threefry','Seed',seedBase,'Substream',runIndex)`. Substreams are independent and order-free, so a `parfor` and a serial loop give identical results and any single run is reproducible in isolation from `(seedBase, runIndex)` alone.

**Draw every parameter from a fixed index in a fixed-length draw vector**, not sequentially as the code happens to need them. Otherwise adding a dispersion later shifts every subsequent draw and silently invalidates comparison with the earlier campaign. Reserve slack in the vector.

**Common random numbers across the two arms.** Both arms fly the *same* realisation of every dispersion, from the same substream. The comparison is then paired, which removes the sampling noise from the difference and is worth roughly an order of magnitude in the number of runs needed.

### What to record, per run

`runIndex`, `seedBase`, the full drawn parameter vector; burnout state (all 7 components) and burn time; **`terminalConstraint`'s five-vector `c` and the `ach` struct at burnout**; impact latitude/longitude, down-range and cross-range miss; peak dynamic pressure and peak axial load (feasibility); maximum commanded `|alpha|` and the fraction of the flight the `alphaMax` clamp was active; guidance cycles flown; per-cycle solver iteration counts and **a converged flag per cycle**; wall time. Failures are recorded, not dropped — a PEG arm that diverges on 3 % of realisations is the most important result the campaign can produce, and discarding those runs would hide it.

### The comparison

Three figures, in this order of importance:

1. **Burnout residual.** Empirical CDF of the Euclidean norm of `c` for the two arms, on log scale. This is the PEG-native metric: PEG targets `c`, not miss distance, and this is the plot that answers the question asked.
2. **Burnout-state scatter.** Two-panel scatter of `(dr, dV)` and `(dgamma, plane angle)` at burnout, open-loop in one colour and PEG in the other, with 1σ and 2σ covariance ellipses per population. Paired realisations may be joined by a line segment on a subsample, which shows the correction per realisation rather than only the aggregate.
3. **Impact scatter.** Down-range versus cross-range in km with CEP circles. Included because it is the quantity a reader will ask for, and captioned to say that it is downstream of a free-flight arc the guidance never saw.

### How many runs

The standard error of an estimated standard deviation is `sigma / sqrt(2(N−1))`: **N ≈ 51** for 10 % precision on the scatter, **N ≈ 201** for 5 %. For CEP the standard error is about `1.25 sigma / sqrt(N)`, so 200 runs give roughly 9 %. A 95th-percentile claim needs N ≥ 500.

**Recommendation: N = 200 for the headline comparison, N = 500 if any tail statistic is quoted.** Below about 50 the two populations' scatter statistics are not distinguishable and the plot will mislead. Paired sampling relaxes this for the *difference* but not for either arm's own scatter.

Cost, measured: `run_ballistic` warm is **0.192 s**, so 500 open-loop runs is about **1.6 minutes**. The PEG arm is a different animal: 81 cycles at Δt = 1 s, each a five-parameter least-squares solve whose every residual evaluation is a forward propagation. Cold-started at ~30 evaluations per cycle that is ~2400 propagations per flight, of order 100 s; warm-started from the previous cycle's `z*` — which is the entire reason PEG is affordable — perhaps 5 evaluations per cycle, of order 16 s per flight, so **500 flights is around 2.2 hours**. The cadence mechanism itself does not enter this budget: 0.09 s of overhead against a 16 s flight is four significant figures below the noise, and counting it here would be double-counting against a cost the same paragraph says is dominated by the solve. Budget it as a campaign: batch it, checkpoint it, make it resumable, and record partial results as they land.

### What this shows, and what it does not

**It shows** that PEG nulls a *plant* dispersion given perfect state knowledge. That is a real and useful experiment: the open-loop schedule was built for the nominal plant, and a 1 % thrust dispersion or a 0.2 deg thrust misalignment moves the burnout state by an amount the open-loop arm cannot recover and the PEG arm should largely remove.

**It does not show closed-loop performance.** There is no navigation model, no estimator, no measurement model. `guide(t,x)` receives the integrator's own truth, exact by construction. The result is therefore an **upper bound** on what a real system achieves, and the gap between it and reality is precisely the navigation error, which is not modelled anywhere in this library. Dispersing the plant while the guidance sees perfect state is a different experiment from closing the loop on an estimate, and the write-up must keep them separate in every sentence. `docs/software_design.tex` §`sec:seam` Gap 3 says this and it remains true after this spike.

**A second experiment is available, and it is cheap — but be precise about how big it is.** PEG's internal propagation is the two-body vacuum model of `docs/hgv_dynamics_note.tex` eq. `guidetwobody`, while the plant has an exponential atmosphere and drag. That is a genuine guidance-model/plant mismatch present *even on the nominal case with perfect state*, no dispersion required. Its size, measured over the **boost phase only** — which is the arc PEG steers:

- integrated drag deceleration **48.37 m/s**, peaking at **1.4601 m/s² at t = 34.02 s**;
- consistent with `run_ballistic`'s own attribution of −0.0396 km of range to drag.

**It is not 2780 m/s.** An earlier edition of this brief quoted that figure, which is `run_ballistic`'s vacuum-versus-flown *impact speed* — a whole-flight number, dominated by the terminal dive of the separated 900 kg re-entry body an hour after burnout, on a different vehicle from the one PEG steers. The boost-phase figure is **1.70 %** of it.

The rest of the mismatch is nil on the shipped case: PEG's `−mu r / r³` **is** `sphereGrav`, and `omegaE = 0`, so gravity and rotation contribute nothing. Task 6 therefore measures **drag plus the zero-order hold, and nothing else** — and the hold alone is 326 m of burnout radius at Δt = 1 s against a 48 m/s velocity perturbation from drag. **The hold may well dominate**, so task 6 must separate the two (run it at two cadences, or with drag switched off through `env`) or its result will be read as a drag floor when it is a hold floor. It still belongs before the dispersion campaign, and it still sets the floor; it is just a floor with two contributions and they must be reported apart.

---

## What is shipped with this brief

### `+coorbital/+guide/terminalConstraint.m`

The five terminal constraints as a dimensionless residual, zero on the target manifold, plus an `ach` struct of the achieved quantities. `docs/hgv_dynamics_note.tex` eq. `termcon`, evaluated and guarded, not re-derived. Both algorithms need it and neither should own it.

**Its interface was settled by this spike, with one deliberate omission.** It takes **Cartesian** `rVec`, `vVec` and a `des` struct, in whatever frame the caller's guidance problem is posed in — which for PEG and VOA is the inertial frame of the two-body equations they integrate. It does **not** take a library state, and it does not convert one. The conversion is a separate utility because its substance is the velocity term `v_inertial = v_rel + omega x r`, which belongs with whoever owns the chain's rotation rate rather than with the one function both guidance laws call in their inner loop. The axis convention is not the obstacle — `thetaG(0) = 0` is free, for the reasons given under Question 2 — but the conversion is real work with a real sign to get wrong, and it is task 1 below.

**And it carries one correction to the specification, which is the most important thing in this file.** Residuals four and five ask that the achieved angular momentum have no `q1` and no `q2` component. That is satisfied when `h_f` is **antiparallel** to `wHat` as well as parallel: a retrograde burnout at the right radius, the right speed, the right flight-path angle and in the right plane, flying round the wrong way, returns

```
    c = [0  2.2e-16  0  0  0]',   plane angle = 180.0000 deg
```

— five exact zeros. The null set of the specified five has **two sheets**, PEG's least-squares and VOA's shooting would both converge happily onto the wrong one, and §3's instruction to score dispersed runs on the norm of `c` would have recorded a wrong-way insertion as a perfect hit. The defect is inherited from the source deck. It is caught here, in the one place both algorithms pass through: `h_f · wHat <= 0` raises `coorbital:terminalConstraint:retrograde`, and a caller wanting the retrograde sheet negates `des.wHat`. A throw rather than a large residual because `coorbital.util.aimSolve` already establishes the contract that a residual which throws marks an infeasible point and shrinks the step — which is the behaviour wanted at the boundary between sheets. **Two further false zeros are guarded the same way**: a zero position or velocity vector is finite, real and three components long, and returns three of five residuals reading as satisfied with no NaN in `c` to give it away; and a radial achieved velocity gives `h = 0`, so both plane residuals are zero for the wrong reason.

`tests/test_terminalConstraint.m` asserts: zero residual on a state constructed to sit exactly on the manifold; each of the first three residuals measuring its own quantity with its stated scaling against hand-computed literals; **`c(3)`'s normalisation specifically, at a state with `r != Rd` AND `V != Vd`** — the only place the `Rd*Vd` and `rMag*vMag` normalisations differ, since everywhere else in the file they agree exactly and a mutation to the wrong one left the suite green; the plane residual magnitude equal to `sin(tilt)` for a rigid tilt of 0.1 to 20 deg out of the target plane; **invariance of `c(1:3)` and `hypot(c(4),c(5))` under a rigid rotation about `wHat`** — the continuous freedom whose transversality condition is VOA's seventh boundary condition, so the test exercises the derivation and not only the arithmetic; **that the specified five really are blind to the retrograde sheet, before asserting that the guard refuses it** — so the guard can never be mistaken for belt and braces — and that negating `wHat` makes the same state a perfect hit; cancellation of `bRef` from the constraint set while it *does* move `c(4)` and `c(5)` individually, and acceptance of a bRef 1e-9 long but perfectly transverse; internal normalisation of `wHat`; reshaping of row-vector inputs; and **fifteen refusals** across seven identifiers, including the radial-burnout case `gammaD = pi/2`, where the desired angular momentum is zero and the plane residuals have no scale.

---

## Ordered task list for the implementation milestone

1. **The state-to-Cartesian conversion**, `stateToCartesian` / `cartesianToState` and their tests. The substance is the velocity term `v_inertial = v_rel + omega x r`, which needs only `env.omegaE` — already present. The axis convention is *not* the work: `thetaG(0) = 0` is a free convention with zero modelling content, because `sphereGrav` is spherically symmetric with no J2, no longitude dependence and no third body, and the target is specified in the same rotating frame the state lives in. Write the convention down in one place and stop thinking about it. Everything downstream needs the conversion: `terminalConstraint` cannot be fed from a library state without it, and the Question 2 frame chain cannot be closed without it. It depends on nothing. **Do it first.**
2. **The direction-to-control mapping**, `+coorbital/+guide/thrustDirection.m` or similar. The `atan2` form; the branch convention stated in the header against `pitchProgram`'s signed `alpha`; the degeneracy guards of Question 2 — all of them, including the `|p|` guard and the reversed-direction band, and *not* a guard at `alpha = 90 deg`; and the feasibility check, whose first job is to **state which struct carries `alphaMaxDeg`**, since `boosterDefaults` does not have the field and the boost phase flies the booster. Record the clamp, because on this vehicle it binds for 74 % of the nominal burn. Testable entirely on its own against synthetic directions — reuse the motor-on/motor-off differencing harness from this spike, which is the strongest available forward truth.
3. **The target mapping**, ground target to `(Rd, Vd, gammaD, wHat)`, built on the validated Keplerian free-flight range and `BM/run_ballistic_target.m`'s `minimum-energy` mode. Ballistic only. Produces the `des` struct `terminalConstraint` already consumes.
4. **The closed-loop driver**, `+coorbital/+prop/closedLoopRun.m` or an entry script. Decision 1: a loop of `phaseRun` calls, one per cycle, holding `z*` in the driver's workspace, concatenating segments into one `traj`, and recording per-cycle solve diagnostics. `phaseRun` is not modified. Test it with a trivial "guidance law" that returns the pitch schedule's command — the live-schedule rows of the cadence table are the pinned expectation, and they say the burnout radius must be within 3e-4 m of the single-phase run at `RelTol = 1e-10`, independent of the cycle count. Pin it against the tolerance, not against a cycle count, because that is how it was measured to behave.
5. **PEG**, against 1–4. Nominal validation is that the Euclidean norm of `c` at burnout goes to zero across the cycles — **not** that the trajectory matches the open-loop one, which the zero-order-hold row of the cadence table says it will not (326 m at Δt = 1 s).
6. **The model-mismatch experiment.** Nominal case, perfect state, PEG's two-body vacuum internal model against the drag-carrying plant. Cheap, no Monte Carlo. It sets the floor for task 8, but the floor has **two** contributions — boost-phase drag (48.37 m/s of integrated deceleration) and the zero-order hold (326 m of burnout radius at Δt = 1 s) — and gravity and rotation contribute nothing on the shipped case. Separate the two, by running at two cadences or with drag switched off through `env`, or the result will be reported as a drag floor when it is a hold floor.
7. **VOA**, against 1–4, as the off-line benchmark PEG's slerp approximation is scored against. Twelve-state augmented system, seven shooting unknowns, the five constraints plus `H(t_f) = 0` plus the orthogonality condition of eq. `voaorth`. Not needed for the demonstration and correctly last among the algorithms.
8. **The dispersion campaign**, per Question 3. After task 6, so its result has a floor to be measured against.

Navigation is not on this list, per `docs/software_design.tex` §`sec:seam` Gap 3. It is larger than everything above put together, it has an obvious home elsewhere, and tasks 1–8 are worth having without it.

---

## Explicitly out of scope, and why

- **PEG or VOA for the HGV.** No closed-form map from burnout state to impact exists when a glide phase decides the range; the burnout states reaching a given target are a family and choosing among them is a modelling decision with no defensible default. `HGV/run_target.m` already solves that problem the only honest way available — propagate and measure.
- **Navigation, estimation, sensors.** No model, no filter, no measurement. Every result in this brief's Question 3 is therefore an upper bound, labelled as one.
- **Attitude dynamics, rate limits, autopilot.** 3-DOF point mass. Commanded directions are achieved instantaneously and exactly; PEG's slerp shape is discarded and only its initial direction is used, which is what the algorithm specifies and is still a fidelity loss.
- **Incidence-dependent aerodynamics.** `constLD` returns `veh.CL` regardless of `alpha`. A 60 deg incidence command costs nothing aerodynamically in this model. Feasibility must come from `veh.alphaMaxDeg`, and this is why task 2's check is not optional.
- **Throttling.** The only control is direction; `constThrust` does not throttle, and the 40.4 g end-of-burn load the shipped case already reports is an artifact of that. Both algorithms take `a_T(t)` as known and neither asks for the magnitude.
- **Multi-stage boost.** One powered phase. `phaseRun`'s `link` and per-phase `veh` would carry it, but nothing in this brief has been measured on a staged ascent.

---

## Reproduction

Every number above came from five scratch scripts run under `restoredefaultpath` with only `missiles/` and `missiles/BM` on the path — no copy of the tree anywhere reachable by `addpath(genpath(...))`. They are scratch and are not part of the library:

- `measure_calls.m` — the call-pattern census, Run 1 (`phaseRun`), Run 2 (bare `ode45` with `deval`), Run 3 (full chain).
- `measure_calls2.m` — adjacent-decrease counting, the raw stage-time sequence around the first rejection, and the deviation split by rejected versus accepted attempt.
- `measure_cadence.m` — the one-phase baseline against the per-cycle driver at Δt ∈ {2, 1, 0.5} s, with and without the command held.
- `verify_inverse.m` — the 720-case round trip through `boost3DOF`'s own thrust terms, the conditioning table, the bank ambiguity, vanishing thrust, and the lift-versus-thrust ratios.
- `verify_frame.m` — the inertial → Earth-fixed → local-horizon → velocity-triad chain, and the `run_ballistic` timing that sizes the campaign.
- `verify_review.m` — the second-edition numbers: the vehicle structs and the boost incidence statistics, the boost-phase drag budget, the scale-robustness table, the extra degeneracies, the worst duplicate-time spread, the restart artifact against cycle count and tolerance, the fourth zero-order-hold point, and the cadence timings as medians of five.

Suite state with `terminalConstraint` and its test added: **23 passed, 0 failed**, zero warnings, no pinned number moved.

## Revision note

This is the second edition. The first was reviewed and reproduced independently — the call census, the zero-order-hold table and the inverse verification all held — and seven things were corrected rather than added to: the retrograde false zero in `terminalConstraint` (Critical, and the only finding that would have produced a wrong answer downstream), two further false zeros on degenerate achieved states, a `c(3)` normalisation the test could not distinguish and a mutation proved it, the misattribution of a 2780 m/s whole-flight drag figure to the boost phase, the overreach of "neither documented option", the claim that a missing epoch blocked the frame conversion, and the conditioning argument for `atan2` — which is true and is not the reason. Each correction is marked where it lands.
