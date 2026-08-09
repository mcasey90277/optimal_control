# Closed-Loop Boost Guidance — Design Spike

**Status:** decisions, measured. No guidance law is implemented and none is designed here.
**Date:** 2026-08-09. **Author:** Michael Casey.
**Scope:** settles the four questions that gate PEG and VOA, by experiment. It is the brief the implementation milestone is written from.

The mathematics is `docs/hgv_dynamics_note.tex` §*Guidance* (`sec:guidance`, `sec:termcon`, `sec:voa`, `sec:peg`) and is not re-derived here — it is cited. The four gaps are `docs/software_design.tex` §*The closed-loop guidance seam* (`sec:seam`). This document converts those gaps into decisions and attaches numbers to them.

**What is shipped alongside this brief:** `+coorbital/+guide/terminalConstraint.m` and `tests/test_terminalConstraint.m`. Nothing else. No PEG, no VOA, no driver, no dispersion machinery.

---

## Summary of the two decisions

**Decision 1 — guidance cadence.** *Neither* documented option. A **driver above `phaseRun`, one `phaseRun` call per guidance cycle**, holding the solved parameter vector in the driver's own workspace. `phaseRun` is not modified, the state vector is not widened, and the purity contract is honoured literally rather than worked around: within a cycle the guide is a **constant**. The measurement below rules out every design in which a guidance *update* fires inside a call to `guide(t,x)`.

**Decision 2 — the thrust-direction inverse.** Resolve the commanded unit direction on the velocity triad and use the **`atan2` form**, not the `arccos` form written in the math note:

```
    alpha = atan2(hypot(p.eGamma, p.ePsi), p.vHat)
    sigma = atan2(p.ePsi, p.eGamma)
```

Round-trip through the library's own boost equations over 720 cases: worst `|dalpha|` **2.22e-16 rad** for the `atan2` form against **1.21e-11 rad** for `arccos`, and the `arccos` form's relative error reaches **100 %** at the small incidences boost actually flies.

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
| Worst `|dr|` between a state shown to the guide and the state flown at that instant | **2.119 m** |
| Worst `|dV|`, worst `|dgamma|` | **0.1279 m/s**, **1.02e-4 rad** |
| Median `|dr|`, interior-stage calls | 0.0122 m |
| Worst `|dr|` on rejected-attempt calls | 0.322 m |

### What that says, item by item

**Times are not monotonic, and the excursion is the size of a guidance cycle.** Fifteen rejected attempts put 20.7 % of the stage calls at times the integrator had already passed, the worst by 0.92 s. A guide that latched "solve when `t` crosses the next cycle boundary" would fire during a rejected attempt, on a trajectory that was never flown, and would then hold that solution across the retry. The failure is silent: the trajectory that comes out is plausible.

**The guide is called twice at the same instant with two different states.** 829 stage calls at only 681 distinct times. Dormand–Prince has two stages at `c = 1`, and their states differ: at `t = 5.979879583 s` the two calls saw radii 3.6 cm apart. **A latch keyed on time cannot work**, because time does not uniquely identify a call.

**States shown to the guide are not states the vehicle occupied — and the worst offenders are the *accepted* steps.** The largest deviation, 2.12 m of radius and 0.128 m/s of speed, is an interior Runge–Kutta stage of an accepted step, not a rejected one, because rejected attempts are precisely the ones the integrator retries with a smaller step. This kills the intuition that "rejected steps are the problem". Every interior stage is an intermediate approximation, and there are 5.7 of them per accepted step.

**The reconstruction pass is worse than documented.** `phaseRun`'s header says it calls the guide once per returned sample after every stage call. It does — 493 calls, at times bit-exactly equal to `traj.t`. But only **123 of those 493 times were ever evaluated during integration**: `ode45`'s default `Refine = 4` puts 4·123 + 1 = 493 output points on the grid, and 370 of them are interpolation times the integrator never asked the guide about. A stateful guide is not merely advanced by the reconstruction pass; it is advanced at times that never existed.

### The two candidate designs, evaluated

**(a) Guidance state rides in the integrated state vector.** The measurement rules out the *continuous* version outright: a discretely updated quantity would have to be updated inside a `guide` or EOM call, i.e. at one of the 829 stage evaluations, 20.7 % of which are on retried arcs and 148 of which duplicate a time. The only honest version is "derivative zero, updated only at phase boundaries" — which is a phase boundary mechanism wearing a state vector as a costume, and it drags along real costs: `nx` grows for every phase in the chain whether it guides or not, `massConstant`'s mass-agreement guard and the three event functions all index into a vector whose width has changed, and `traj.x` acquires columns that are not states.

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

Same boost case; the guidance solve is *not* included, so this isolates the cost of the mechanism.

| Δt | Wall time | vs. one phase | Cycles | Burnout Δr | ΔV | Δγ |
|---|---|---|---|---|---|---|
| — (one phase) | 0.0393 s | 1.00× | 1 | — | — | — |
| 2.0 s, held command | 0.1002 s | 2.55× | 41 | **+670.3 m** | −2.629 m/s | +1.095e-2 rad |
| 1.0 s, held command | 0.1330 s | 3.39× | 81 | **+325.6 m** | −1.282 m/s | +5.301e-3 rad |
| 0.5 s, held command | 0.1526 s | 3.88× | 162 | **+162.5 m** | −0.640 m/s | +2.641e-3 rad |
| 2.0 s, live schedule | 0.0708 s | 1.80× | 41 | +1.2e-6 m | −1.2e-6 m/s | −1.9e-11 rad |
| 1.0 s, live schedule | 0.1174 s | 2.99× | 81 | −1.1e-4 m | −9.6e-7 m/s | −1.8e-9 rad |
| 0.5 s, live schedule | 0.2233 s | 5.69× | 162 | −1.7e-4 m | −7.6e-7 m/s | −2.8e-9 rad |

Two things fall out, and they must not be confused with each other.

**The restart artifact is nil.** Cutting the boost into cycles while letting the guide keep reading the live schedule moves the burnout radius by at most **0.17 mm** and the speed by 1.2e-6 m/s. Restarting `ode45` 162 times does not damage the integration. That is the whole risk of Decision 1 and it is measured to be absent.

**The zero-order hold is not nil, and it is first order in Δt.** Holding the command across the cycle costs **326 m of burnout radius at Δt = 1 s**, halving cleanly to 163 m at 0.5 s and doubling to 670 m at 2 s. This is not a numerical artifact — it is what holding a command costs, and a real flight computer pays it too. Its consequence for validation is in the task list below: **the nominal closed-loop run will not reproduce the open-loop trajectory, and must not be validated against it.**

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

`docs/hgv_dynamics_note.tex` eq. `dirtoalpha` gives `alpha = arccos(p·vHat)`. It is analytically identical and numerically the wrong branch to implement; see the conditioning table.

The frame chain a PEG command has to travel, and which the implementation owns:

```
    p_inertial --(rotate by thetaG about z)--> p_ecef --(lon,lat)--> [n,e,u] --(gamma,psi)--> [vHat,eGamma,ePsi]
```

It is a pure rotation at every stage — a *direction*, not a velocity, so no `omega x r` term enters and inertial-versus-planet-relative does not arise for the command itself. The triad is built from the *planet-relative* velocity, which is the right thing: `alpha` is then the aerodynamic incidence the equations and `env.aero` expect. Verified: triad Gram matrix minus identity worst element **4.44e-16**, `vHat × eGamma − ePsi` worst **1.11e-16**, and round trip through `thetaG ∈ {0, 0.4, 2.9}` rad recovers `alpha` to **1.25e-16 rad** and `sigma` to **2.0e-15 rad**. **The library has no epoch, so `thetaG` is the caller's to supply** — the same gap flagged for `terminalConstraint` below.

### Numerical verification

Forward truth was taken from `boost3DOF` itself, by differencing the derivative with the motor on against the derivative with the motor off. Aerodynamics, gravity, kinematics and the rotation terms cancel exactly in that difference, leaving the thrust acceleration resolved on the triad:

```
    a.vHat = d(Vdot),   a.eGamma = d(gammadot)*V,   a.ePsi = d(psidot)*V*cos(gamma)
```

720 cases: 4 states (2–120 km altitude, 200–6000 m/s, gamma from 0 to 80 deg, four headings, mid-latitude) × 18 values of alpha (0 to ±179 deg) × 10 values of sigma (0 to ±180 deg).

| Result | Value |
|---|---|
| Worst `|dalpha|`, **`atan2` form** | **2.22e-16 rad** (1.27e-14 deg) |
| Worst `|dalpha|`, `arccos` form | 1.207e-11 rad (6.92e-10 deg) |
| Worst `|dsigma|` | 8.97e-13 rad (5.14e-11 deg) |
| With `env.omegaE` on | `dalpha` **0**, `dsigma` 1.11e-16 rad |

Conditioning near zero incidence, exact analytic forward, so only the inverse is being measured:

| alpha (rad) | `arccos` error | `atan2` error |
|---|---|---|
| 1e-10 | **1e-10** (100 % relative) | 0 |
| 1e-8 | **1e-8** (100 % relative) | 1.65e-24 |
| 1e-6 | 4.45e-11 | 0 |
| 1e-4 | 2.62e-13 | 0 |
| 1e-2 | 1.44e-15 | 0 |
| 1e-1 | 5.55e-16 | 1.39e-17 |

`arccos(1 − eps)` has unbounded derivative at the identity, so the `arccos` form returns garbage exactly where boost flies. The shipped ballistic ascent runs at incidences of a few degrees or less, clamped at `veh.alphaMaxDeg = 6 deg`; near thrust-along-velocity the `arccos` form has no digits left.

### Where it degenerates

1. **The forward map is two-to-one.** `(alpha, sigma)` and `(−alpha, sigma ± pi)` produce the same direction. Any inverse must pick a branch; the `atan2` form returns `alpha ∈ [0, pi]`. **This conflicts with the existing convention** — `pitchProgram` returns *signed* `alpha`, and `BM/run_ballistic.m` flies negative incidence on parts of its schedule. The implementation must state its branch in the header and the consumer must not assume it matches `pitchProgram`'s.
2. **`sigma` is undefined at `alpha = 0`, and MATLAB returns 0 silently.** `atan2(0,0)` is `0`, so a direction exactly along the velocity comes back as `sigma = 0` with no warning: measured **40 deg of error** against the bank actually commanded. There is no incidence to roll, so any `sigma` is correct for the *direction* — but it is not correct for the *lift vector*, which `sigma` also orients. Guard: below a threshold incidence, hold the previous cycle's `sigma` and say so, or refuse. Do not return the atan2 default.
3. **`sigma`'s sensitivity is `|dp| / sin(alpha)`.** Measured, one-ulp perturbation of `p`: `sigma` moves by 0.0179 deg at `alpha = 1e-12`, 1.79e-5 deg at 1e-9, 1.79e-8 deg at 1e-6, 1.79e-11 deg at 1e-3. With a PEG direction converged to a tolerance `d`, `sigma` is meaningless below `alpha ≈ d`. Set the threshold in item 2 from the guidance solve tolerance, not from `eps`.
4. **`alpha = pi`.** Thrust reversed along the velocity: same ambiguity, and physically absurd on a boost. Refuse.
5. **Vanishing thrust.** The direction is undefined when `T = 0`. Measured with the delivered thrust scaled down: at scale 1e-6 the inverse still recovers 6.0000 deg; at 1e-12 it has drifted to **5.99965 deg**; at exactly zero it returns **NaN**. Guard on `T/m` against a physical threshold, and raise rather than return NaN — `boost3DOF` refuses non-finite controls anyway, but it would refuse *downstream* of the mistake.
6. **Feasibility is not the dynamics' job.** `veh.alphaMaxDeg = 6 deg` in `BM/vehicle_bm.m` is the vehicle's control-authority limit. A 3-DOF point mass will fly any direction it is handed; there is no attitude state to violate. **The mapping must carry the feasibility check the dynamics cannot**, and must decide — explicitly — whether an infeasible command is clamped (with the clamp recorded, so the dispersed runs can report how often it bound) or refused.

### What the 3-DOF point-mass model costs here

There are no attitude dynamics, no moments and no rate limits. Orientation is an instantaneous input, so **a commanded direction is achieved instantaneously and exactly**. Three consequences, stated plainly:

- **The cycle boundary is a step change.** PEG's slerp (`docs/hgv_dynamics_note.tex` eq. `slerp`) exists to turn the thrust direction smoothly over the arc, with a quadratic warp chosen so the command leaves `p0` at zero rate. This library commands `p(0; z*)` and holds it, then jumps at the next boundary. The *shape* PEG is built around is discarded; only the endpoint it produces is used. That is what the algorithm specifies (`pegcmd`) and it is still a loss of the fidelity the shape was designed for.
- **No pitch-rate limit can be violated.** The shipped attitude schedule averages 0.67 deg/s (89 deg to 34 deg over 82 s), comfortably inside any plausible booster limit — so the *nominal* is not the concern. A dispersed PEG correction demanding a large re-point is, and this model will execute it in zero time and report success.
- **`sigma` owns the lift vector as well as the thrust.** Quantified at 6 deg of incidence on the boosted stack: `a_lift` is 0.0716 m/s² at 2 km, **0.2115 m/s² at 20 km** (the worst case), 0.0070 at 60 km and 4.9e-6 at 120 km, against `a_thrust_normal = 4.327 m/s²`. So the lift the bank angle drags around is at most **4.9 %** of the normal thrust it is being commanded for. Note also that **`coorbital.aero.constLD` ignores `alpha` entirely** (`CL = veh.CL` always), so incidence-induced lift is not modelled at all — the coupling is only the rotation of a constant-magnitude lift vector, and a 60 deg incidence command would produce no aerodynamic penalty whatsoever. The feasibility check of item 6 is therefore the *only* thing standing between PEG and a physically absurd command.

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

Cost, measured: `run_ballistic` warm is **0.192 s**, so 500 open-loop runs is about **1.6 minutes**. The PEG arm is a different animal: 81 cycles at Δt = 1 s, each a five-parameter least-squares solve whose every residual evaluation is a forward propagation. Cold-started at ~30 evaluations per cycle that is ~2400 propagations per flight, of order 100 s; warm-started from the previous cycle's `z*` — which is the entire reason PEG is affordable — perhaps 5 evaluations per cycle, of order 16 s per flight, so **500 flights is around 2.2 hours**. Plus the measured 3.39× cadence overhead. Budget it as a campaign: batch it, checkpoint it, make it resumable, and record partial results as they land.

### What this shows, and what it does not

**It shows** that PEG nulls a *plant* dispersion given perfect state knowledge. That is a real and useful experiment: the open-loop schedule was built for the nominal plant, and a 1 % thrust dispersion or a 0.2 deg thrust misalignment moves the burnout state by an amount the open-loop arm cannot recover and the PEG arm should largely remove.

**It does not show closed-loop performance.** There is no navigation model, no estimator, no measurement model. `guide(t,x)` receives the integrator's own truth, exact by construction. The result is therefore an **upper bound** on what a real system achieves, and the gap between it and reality is precisely the navigation error, which is not modelled anywhere in this library. Dispersing the plant while the guidance sees perfect state is a different experiment from closing the loop on an estimate, and the write-up must keep them separate in every sentence. `docs/software_design.tex` §`sec:seam` Gap 3 says this and it remains true after this spike.

**A second experiment is available and may be the more interesting one.** PEG's internal propagation is the two-body vacuum model of `docs/hgv_dynamics_note.tex` eq. `guidetwobody`, while the plant has an exponential atmosphere and drag. That is a genuine guidance-model/plant mismatch present *even on the nominal case with perfect state* — no dispersion required — and it will dominate the low-altitude part of the boost, where the measured drag on this vehicle is worth 2780 m/s of impact speed. Running the nominal case and reporting how much terminal residual PEG's own model error leaves is cheap, needs no Monte Carlo, and should be done **before** the dispersion campaign, because it sets the noise floor the campaign is measured against.

---

## What is shipped with this brief

### `+coorbital/+guide/terminalConstraint.m`

The five terminal constraints as a dimensionless residual, zero on the target manifold, plus an `ach` struct of the achieved quantities. `docs/hgv_dynamics_note.tex` eq. `termcon`, evaluated and guarded, not re-derived. Both algorithms need it and neither should own it.

**Its interface was settled by this spike, with one deliberate omission.** It takes **Cartesian** `rVec`, `vVec` and a `des` struct, in whatever frame the caller's guidance problem is posed in — which for PEG and VOA is the inertial frame of the two-body equations they integrate. It does **not** take a library state, and it does not convert one. Converting the seven-state spherical planet-relative vector to an inertial Cartesian pair needs an inertial frame, and **this library has not defined one**: with `env.omegaE` nonzero the longitude is Earth-fixed and the speed is planet-relative, and there is no epoch anywhere that says where the inertial x-axis points. That is the same missing `thetaG` as in the frame chain of Question 2. Baking a guess about it into the one interface both guidance laws share would put the guess where nobody would look for it. **The conversion belongs with the caller that owns the epoch, and defining the epoch is a task in the list below.**

`tests/test_terminalConstraint.m` asserts: zero residual on a state constructed to sit exactly on the manifold; each of the first three residuals measuring its own quantity with its stated scaling against hand-computed literals; the plane residual magnitude equal to `sin(tilt)` for a rigid tilt of 0.1 to 20 deg out of the target plane; **invariance of `c(1:3)` and `hypot(c(4),c(5))` under a rigid rotation about `wHat`** — the one-parameter freedom whose transversality condition is VOA's seventh boundary condition, so the test exercises the derivation and not only the arithmetic; cancellation of `bRef` from the constraint set while it *does* move `c(4)` and `c(5)` individually; internal normalisation of `wHat`; and seven degenerate inputs each raising its own identifier, including the radial-burnout case `gammaD = pi/2`, where the desired angular momentum is zero and the plane residuals have no scale.

---

## Ordered task list for the implementation milestone

1. **Define the inertial frame and the epoch.** One line in `missileConst` or a documented convention ("the inertial frame coincides with the Earth-fixed frame at t = 0 of the chain"), plus `stateToCartesian` / `cartesianToState` utilities and their tests. Everything downstream needs it: `terminalConstraint` cannot be fed from a library state without it, and the Question 2 frame chain cannot be closed without it. It depends on nothing. **Do it first.**
2. **The direction-to-control mapping**, `+coorbital/+guide/thrustDirection.m` or similar. The `atan2` form; the branch convention stated in the header against `pitchProgram`'s signed `alpha`; the five degeneracy guards of Question 2; the `veh.alphaMaxDeg` feasibility check with an explicit clamp-or-refuse decision and a recorded clamp flag. Testable entirely on its own against synthetic directions — reuse the motor-on/motor-off differencing harness from this spike, which is the strongest available forward truth.
3. **The target mapping**, ground target to `(Rd, Vd, gammaD, wHat)`, built on the validated Keplerian free-flight range and `BM/run_ballistic_target.m`'s `minimum-energy` mode. Ballistic only. Produces the `des` struct `terminalConstraint` already consumes.
4. **The closed-loop driver**, `+coorbital/+prop/closedLoopRun.m` or an entry script. Decision 1: a loop of `phaseRun` calls, one per cycle, holding `z*` in the driver's workspace, concatenating segments into one `traj`, and recording per-cycle solve diagnostics. `phaseRun` is not modified. Test it with a trivial "guidance law" that returns the pitch schedule's command — the live-schedule row of the cadence table is the pinned expectation, and it says the burnout state must be within 2e-4 m of the single-phase run.
5. **PEG**, against 1–4. Nominal validation is that the Euclidean norm of `c` at burnout goes to zero across the cycles — **not** that the trajectory matches the open-loop one, which the zero-order-hold row of the cadence table says it will not (326 m at Δt = 1 s).
6. **The model-mismatch experiment.** Nominal case, perfect state, PEG's two-body vacuum internal model against the drag-carrying plant. Cheap, no Monte Carlo, sets the noise floor for task 8.
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

Every number above came from four scratch scripts run under `restoredefaultpath` with only `missiles/` and `missiles/BM` on the path — no copy of the tree anywhere reachable by `addpath(genpath(...))`. They are scratch and are not part of the library:

- `measure_calls.m` — the call-pattern census, Run 1 (`phaseRun`), Run 2 (bare `ode45` with `deval`), Run 3 (full chain).
- `measure_calls2.m` — adjacent-decrease counting, the raw stage-time sequence around the first rejection, and the deviation split by rejected versus accepted attempt.
- `measure_cadence.m` — the one-phase baseline against the per-cycle driver at Δt ∈ {2, 1, 0.5} s, with and without the command held.
- `verify_inverse.m` — the 720-case round trip through `boost3DOF`'s own thrust terms, the conditioning table, the bank ambiguity, vanishing thrust, and the lift-versus-thrust ratios.
- `verify_frame.m` — the inertial → Earth-fixed → local-horizon → velocity-triad chain, and the `run_ballistic` timing that sizes the campaign.

Suite state with `terminalConstraint` and its test added: **23 passed, 0 failed**, zero warnings, no pinned number moved.
