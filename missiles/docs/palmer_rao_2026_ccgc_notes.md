# Palmer & Rao 2026: constrained computational guidance and control (CCG&C)

**Paper.** E. M. Palmer and A. V. Rao, "Method for Constrained Computational
Guidance and Control with Application to Hypersonic Entry," arXiv:2609.10813
[math.OC], posted 2026-09-09, 27 pp. Local copy (PDFs are gitignored):
`papers/2026_Palmer_Rao_Constrained_Computational_Guidance_Hypersonic_Entry_arXiv2609.10813.pdf`.

**Read by:** Claude, 2026-09-11. **Status:** notes only. Nothing is coded.

**Goal set by Mike 2026-09-11:** code the algorithm, reproduce the paper's
results, then upgrade the HGV code with it if that proves appropriate. Record
anything that transfers to the orbit-transfer work, DRO-to-tulip especially.

---

## 1. What the paper does

**The loop (their Ref. 20, Dennis, Hager and Rao, JGCD 2019).** Every guidance
cycle, `dT = 15 s`, re-solve the FULL remaining-horizon nonlinear optimal
control problem from the actual state. Transcription is LGR collocation
(GPOPS-II + IPOPT, full Newton). Fly the solved control, as a cubic spline, on
a *perturbed* plant (`ode113`, rel. tol 1e-8) until the cycle ends. Then take
the actual state as the next initial condition. The expired part of the mesh
is truncated and the rest remapped, so later cycles are smaller and faster.

**The problem it fixes.** If the reference solution rides a path constraint
(here the stagnation heating rate, active for 26.5 % of the flight,
`t in [168.69, 725.87] s`), the perturbed plant overshoots the limit. The next
initial condition is then outside the feasible set, and the NLP is
infeasible. The baseline holds the last feasible control until a solve
succeeds. With a +1 % density bias the baseline lost **30 guidance cycles**.
Over a 1000-run Monte Carlo it violated the heating limit on **every** run,
with a mean of 35 and a maximum of 44 infeasible cycles per run.

**The fix, CCG&C.** Keep every path constraint hard, and add a term to the
objective that buys MARGIN from the limit:

```
J_a = w1 * J  +  w2 * J_Pbeta,          w1 + w2 = 1
J_Pbeta = P_beta( c_i / eta_i )          (upper limit;  -c_i/eta_i for a lower one)
P_beta(y) = (1/beta) * log( integral_{t0}^{tf} exp(beta * y(t)) dt )    "log-integral-exponential"
```

`P_beta` is a smooth stand-in for `max_t y(t)`, which is nonsmooth. The
paper uses `beta = 5` and `eta_i = c_max,i`.

**Adaptive weights by threshold (Eq. 39).** At the start of each cycle the
weights are set from the ACTUAL constraint value at the end of the previous
cycle: `(w1, w2) = (1, 0)` if `c/c_max < xi`, otherwise
`(zeta1, zeta2) = (0.8, 0.2)`. Here `xi = 0.9`.

**The side effect, and a second fix (Eq. 40).** Switching the objective on and
off makes `alpha` and `sigma` jump between cycles, which excites phugoid
oscillations in `q` and sensed acceleration. So a third term is added whenever
the margin term is on:

```
w3 * integral ( C/(1 + exp(-k sin gamma)) - C/2 )^2 dt,   C = 2(1+e^-k)/(1-e^-k),   k = 3
```

It is a smooth, saturating penalty on flight-path angle, zero at `gamma = 0`.
They use `zeta3 = 5`.

**The demonstration problem.** Betts' Shuttle-class RLV maximum-crossrange
entry (Betts 2020, their Ref. 27). Point mass, spherical NON-rotating Earth,
exponential atmosphere. The state is `[r, theta(lon), phi(lat), v, gamma,
psi, alpha, sigma]`, and the controls are the RATES `u_alpha, u_sigma`, so
attitude stays continuous across cycles. There are three path constraints:
Sutton-Graves stagnation heating, dynamic pressure, and sensed acceleration.
Only heating is ever active. The objective is to maximize final latitude.

## 2. The numbers to reproduce

| quantity | CCG&C | CCG&C + phugoid penalty | CG&C baseline |
|---|---|---|---|
| reference `phi(t_f)` (offline, nominal) | 34.0 deg | same | same |
| MC mean `phi(t_f)` [deg] | 32.849 | 33.895 | 34.000 |
| MC max of max heating [MW/m^2] (limit 0.85) | 0.8304 | 0.8312 | **0.8618** |
| MC std of max heating [MW/m^2] | 2.17e-3 | 3.25e-4 | 1.95e-3 |
| infeasible guidance cycles per run | 0 | 0 | mean 35, max 44 |
| mean \|dh(t_f)\| / max [m] | 7.86 / 91.8 | 8.46 / 54.1 | 8.91 / 91.4 |
| mean solve time per cycle [s] | 0.213 | 0.436 | 0.090 (see note) |

The baseline's solve time is not comparable, because many of its early and
expensive cycles were infeasible and skipped. The single-case weight sweeps
are Table 4, over `(zeta1, zeta2)`, and Table 5, over `zeta3 = 1, 3, 5`.
Table 5's `zeta3 = 5` row is 33.898 deg with max heating 0.8299.

**The dispersion is a single constant density bias**,
`rho0_tilde ~ N(rho0, (0.01 rho0)^2)`, held for the whole run. There are no
navigation errors, no winds and no aero dispersions.

## 3. Reproduction spec, including the paper's defects

Constants are in Table 1 and limits in Table 2. `alpha` is in [-90, 90] deg,
`sigma` in [-90, 1] deg, `|u_alpha| <= 0.5 deg/s` and `|u_sigma| <= 5 deg/s`.
The limits are heating 0.85 MW/m^2, `q` 16.375 kPa and `A` 2.5 g. Nose radius
is 1 m. Items marked CHECKED were computed here on 2026-09-11.

1. **Eq. 28 swaps the coefficient labels.** Its printed right-hand sides give
   `C_D` the lift coefficients and `C_L` the drag ones. Use the numbers in
   Table 1: `C_L = -0.207 + 0.029245 a`, and
   `C_D = 0.0785 - 6.1593e-3 a + 6.2142e-4 a^2`, with `a` in degrees.
2. **Table 1's `K_q` unit is wrong.** It reads W/cm^2, but in SI (rho in
   kg/m^3, v in m/s, r_n in m) the formula gives **W/m^2**. CHECKED: 0.390
   MW/m^2 at the entry state (79.2 km, 7800 m/s), which matches the
   starting point of Fig. 4a.
3. **Table 3's radii are altitudes in disguise.** "645 km" and "639.6 km" are
   Betts' 260,000 ft and 80,000 ft entry and terminal altitudes, written as
   radii in units of 1e4 m. CHECKED: `r0 = 6450.248 km`, `rf = 6395.384 km`
   with `Re = 6371 km`. `v0 = 25,600 ft/s = 7802.88 m/s`,
   `vf = 2,500 ft/s = 762 m/s`, `m = 203,000 lb = 92,079 kg` and
   `S = 2690 ft^2 = 249.9 m^2` all match exactly.
4. **`alpha(t0)` and `sigma(t0)` are not given.** Fig. 6 shows about 20 deg
   and about -50 deg. Most likely they are free in the offline reference
   solve and then fixed to the actual state in each re-solve. This is a spec
   decision to record.
5. **`t_f` is free but never stated.** It is about 2100 s from the figures, as
   in Betts' free-final-time problem. Bound it in the NLP.
6. **`zeta3 = 5` breaks the stated rule** `w1 + w2 = 1` of Eq. 14. With
   `(0.8, 0.2, 5)` the weights are not normalized. Implement Eq. 42 as
   printed and note the inconsistency.
7. **`P_beta` is not close to the max at `beta = 5`.** Because
   `exp(beta*y) <= exp(beta*max y)`, the integral is at most
   `(t_f - t_0) exp(beta*max y)`, so `P_beta(y) <= max y + (1/beta)
   log(t_f - t_0)`. Over a ~2000 s horizon that offset is about 1.5, while
   `y = c/c_max` is at most 1. So as used, the term is a softmax-weighted
   time integral of the constraint. Its gradient weights are
   `exp(beta*y) / integral exp(beta*y)`, which puts only e^5 = 148 times more
   weight on the peak than on zero heating. The offset also changes as the
   horizon shrinks. This is not an error, since the gradient is what acts,
   but it changes how `beta` should be read. The normalized form
   `(1/beta) log( (1/(t_f-t_0)) integral exp(beta y) dt )` lies between the
   mean and the max. Test both.
8. **The mesh truncation and remap step is out of scope in the paper.** It is
   specified in Ref. 20 (Dennis, Hager, Rao 2019), and we do not have that
   paper yet. Substitute: re-solve on a normalized-time grid over the
   remaining horizon, warm-started by interpolating the previous solution.
9. **"Remaining horizon too small to solve" is not defined.** When that
   happens the previous control is held. Pick and document a rule.

**Tooling.** GPOPS-II is not installed here. CasADi 3.7.0 + IPOPT is, and
both booster_landing (Hermite-Simpson) and DRO_tulip direct already use it.
Plan: Hermite-Simpson first, reproducing the reference solution to 34.0 deg
and the active arc [168.69, 725.87] s. Build LGR collocation only if
transcription differences prevent a match. Implement `P_beta` as an extra
state `z' = exp(beta*y)` and `J_Pbeta = log(z(t_f))/beta`. Evaluate the
exponent stably (Blanchard, Higham and Higham 2021 is their Ref. 23). Overflow
is not a risk at `beta*y <= 5`.

**What "reproduced" should mean.** The qualitative claims must hold exactly:
zero infeasible cycles for CCG&C, and heating violated on every baseline run.
The table numbers should match to the extent a different transcription
allows. A 1000-run MC at about 0.2-0.4 s per solve and about 140 cycles per
run is roughly 10-15 CPU hours. It parallelizes with `parfor`.

## 4. What the paper does not show (baselines our reproduction should add)

- **Constraint tightening.** Re-solve with `c <= c_max - delta` (back-off),
  the standard robust-MPC remedy. It is one line and needs no new objective
  term. CCG&C's MC margin is about 2.2 % (0.8312 against 0.85). Does a 2 %
  tightening give the same crossrange with no phugoid? The paper never asks.
- **Continuous weights.** The threshold rule makes the objective jump between
  cycles, and the jump is what excites the phugoid. A smooth ramp in
  `c/c_max`, or a control-rate penalty (they suggest this in section 7),
  attacks the cause rather than the symptom.
- **Soft recovery.** An exact-penalty slack on the path constraint restores
  feasibility after a violation. This is better than holding a stale control.
- **Richer dispersions.** Density bias only is a lower bound on difficulty.
  Add navigation error, winds and aero coefficient dispersions.
- **Tail timing.** Only means and standard deviations are reported. The
  standard deviation (0.69-0.86 s) exceeds the mean, so report max and p99
  latency against the 15 s cycle.

## 5. Relation to our HGV code

**Where the library is (2026-08-09 state).**

- Guidance is open-loop. The `prescribed` and `pitchProgram` schedules are
  fixed before flight.
- The only closed loop is targeting (`aimSolve`), and it is closed *outside*
  the flight.
- `constLD` ignores `alpha`, so incidence costs nothing.
- Heating, `q` and g-load are not computed at all. `vehicleDefaults` reserves
  `noseRadius` "until Sutton-Graves heating arrives".
- `closed_loop_guidance.md` scoped the HGV OUT of PEG and VOA, because no
  closed-form burnout-to-impact map exists when a glide phase decides the
  range.

**Why CCG&C fits where PEG did not.** It is a numerical OCP re-solve, so it
needs no closed-form terminal map. It handles exactly the case the brief set
aside: guidance of the glide phase itself, with path constraints.

**Architecture fit is good.**

- **State.** The paper's `(r, theta, phi, v, gamma, psi)` is our
  `[r, lon, lat, V, gamma, psi]` exactly, `psi` clockwise from north in both.
  `glide3DOF` uses the same lift and bank form (`L cos sigma`,
  `L sin sigma/(V cos gamma)`) as Eq. 26 when `earthSpin = false`.
- **Cadence.** Decision 1 of `closed_loop_guidance.md` is a driver above
  `phaseRun`, one `phaseRun` per guidance cycle, with parameters held in the
  driver's workspace. That is the CCG&C loop's shape. The one difference is
  that CCG&C flies a time-varying schedule within a cycle. That is still pure,
  as `ph.guide = @(t,x) coorbital.guide.prescribed(t,x,schedFromSolve)`, with
  the schedule fixed for the cycle.

**What the upgrade needs, in order.**

1. **Heating, dynamic pressure and sensed acceleration as reported
   quantities.** Sutton-Graves, `q = rho v^2/2`, `A = sqrt(L^2+D^2)/g0`. This
   is useful even without guidance.
2. **An alpha-dependent aero model.** `coorbital.aero.polyAlpha` (a working
   name) with Betts' Shuttle polynomial, injected as a handle. It is the first
   non-placeholder airframe in the library.
3. **Paper constants without breaking "no hard-coded constants outside
   `missileConst`".** Ours differ from the paper's: `rE` 6378137 against
   6.371e6 m, `rho0` 1.225 against 1.2256, `Hscale` 7200 against 7254 m. And
   `expAtmos` reads `missileConst` directly. The reproduction needs a
   parameterized atmosphere handle or a constants override. Decide at spec
   time.
4. **An OCP solver layer.** CLAUDE.md says the library is "not a trajectory
   optimizer", so this is a real scope change and a new CasADi dependency.
   Recommendation: the transcription lives in `../oclib/+oc` (cross-folder,
   next to booster_landing's HS machinery), and missiles gets a thin guide
   wrapper.
5. **The guidance driver above `phaseRun`,** as specified in Decision 1.
   CCG&C becomes its first client, ahead of PEG.
6. **The MC dispersion harness** from `closed_loop_guidance.md` Question 3.
   The paper's density bias is its first dispersion.

**Is it appropriate for the HGV?** Probably yes for the glide phase, with
caveats. The HGV vehicle file is all placeholders, with no heating or
structural limits defined, so a constrained guidance law has nothing real to
respect until limits exist. Reproduce on the paper's Shuttle problem first:
it is fully specified and published, and it gives the library its first
external anchor for an optimized entry. Port afterward, with an HGV objective
such as terminal miss or energy at a target plus a heating margin, and judge
it against the simpler baselines in section 4 before adopting it.

## 6. Transferable ideas for orbit transfer (DRO-to-tulip)

1. **A riding constraint makes closed-loop re-solves infeasible, and a
   min-time arc rides one the whole way.** Min-time low thrust is all-burn:
   `u <= 1` is active for 100 % of the flight. A thrust or Isp deficit, or a
   navigation error, can leave the remaining min-time problem with no
   margin. Near `t_f` the target may even become unreachable. That is the
   exact analogue of the heating-riding arc. The fixed-`t_f` catalog entries
   at `t_f = gamma * t_f^min` already are the margin product: throttle
   authority held in reserve. The paper's framing, trading performance for
   margin so the re-solve stays feasible, is a reason to read `gamma` as a
   guidance margin.
2. **A smooth-min clearance objective for the lunar altitude floor.** The
   2026-08-03 floor experiments (DRO_tulip FINDINGS) found two problems. The
   floor-riding branch could not be resolved on a uniform mesh. And
   node-only enforcement let the true trajectory dip 120 km below a
   500 km floor. Adding `P_beta(-rho_Moon/eta)` as a weighted secondary
   objective, with the floor kept hard, pushes the arc off the floor. That
   gives a clearance-for-time trade curve, the analogue of the paper's
   crossrange-for-heating curve, and avoids the unresolvable riding arc.
3. **Shrinking-horizon re-solve is cheap with multiple shooting.** The
   analogue of mesh truncation and remap is to drop the expired shooting
   segments and re-seed from the remaining junction states. Seeded-at-a-root
   `ms_tfmin` solves converge in 1-2 Newton iterations (`seed_from_z8`). So
   closed-loop execution of a catalog transfer under dispersions is a
   natural experiment. `cond(J)` is known to lead the conjugate test (FINDINGS
   39), so it is a candidate real-time warning as the horizon shrinks.
4. **Smoothing literature.** Log-sum-exp is the standard smoothing of a max.
   It is the time-max counterpart of our eps, huber and huberc smoothings of
   `|u|`, and belongs in the smoothing-benchmark literature list. It is not a
   competitor for the L1 fuel term.
5. **Objective switching excites oscillation.** If a replanning loop ever
   switches cost terms on and off between cycles, penalize control rate or
   ramp the weights. The paper's phugoid is the entry-flight symptom of
   that.

## References worth fetching for the reproduction

- Dennis, Hager, Rao, "Computational Method for Optimal Guidance and Control
  Using Adaptive Gaussian Quadrature Collocation," JGCD 42(9), 2019 (the
  CG&C baseline and the mesh truncation and remap).
- Betts, *Practical Methods for Optimal Control and Estimation Using Nonlinear
  Programming*, 3rd ed., SIAM 2020 (the Shuttle max-crossrange problem).
- Palmer & Rao, "Adaptive Path Constrained Optimal Guidance with Application
  to Reusable Launch Vehicle Entry," AAS/AIAA SFM 2025 (preliminary version).
- Miller & Rao, JSR 59(3), 2022 (the phugoid penalty's source).
- Blanchard, Higham, Higham, IMA J. Numer. Anal. 41(4), 2021 (stable
  log-sum-exp).
