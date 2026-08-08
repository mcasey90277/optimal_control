# Two-Axis Targeting and the Closed-Loop Guidance Seam — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Targeting that works with Earth rotation on and with a banked trajectory, by solving two controls against two residuals instead of one against one; and the scaffolding that PEG and VOA closed-loop boost guidance will attach to.

**Architecture:** Today both targeting scripts solve ONE control (range) against ONE residual (down-range), and take the launch azimuth from a closed form that is exact only for a non-rotating, zero-bank great circle. Neither assumption survives contact with a rotating Earth or a banked descent, so both scripts currently refuse under rotation and merely warn under bank. The fix is a shared two-dimensional root solve on `(azimuth, range control)` driving `(down-range residual, cross-track)` to zero. That single utility closes both gaps, because they are the same gap: the flown ground track is not the great circle the closed form assumed.

**Tech Stack:** MATLAB R2025b, the existing `missiles/+coorbital` library.

---

## Global Constraints

- pumpkyn house style: `%%` header quartet (`Purpose`/`Inputs`/`Outputs`/`Revision History`) with aligned three-column I/O blocks stating UNITS, closed by `%% ------------------------ Begin Code Sequence ---------------------------`; `=` signs column-aligned at column 20; colon-terminated `%%` step comments; `if nargin == 0` self-demo calling by full namespace.
- `%% References:` whenever the math has a source.
- Never `i`/`j` as loop or index variables. Never `norm`.
- **No `%#ok` pragmas of any kind.**
- **No hard-coded physical constants** outside `coorbital.util.missileConst()`. Entry-script user blocks may hold user-chosen values.
- SI in the library, radians for angles; human units only in a user block, converted immediately in one marked place.
- Author `Michael Casey`, date 08/08/2026, `Copyright 2026 Coorbital, Inc.`
- **Never modify `~/Desktop/proj7/external/pumpkyn`.**
- Suite: `run('tests/run_tests')` under `-batch`; a bare `tests/run_tests` parses as division. Allow 900 s.
- Baseline `21 passed, 0 failed`, zero warnings. **`run_glide` 6986.82/2073.77/1.11 and `run_boost_glide` 7663.05/2194.77/8.33 must not move.** The two targeting scripts will move and their pins must be re-measured, not forced.

---

## Why one solver closes both gaps

The present scheme assumes the flown ground track is the great circle from launch to target, so the launch azimuth is that arc's initial bearing and only the range needs solving. Two things break it, and they break it the same way.

**Earth rotation.** With `omegaE` nonzero the Coriolis term deflects the trajectory. Measured: `run_target` misses by 231.6 km and `run_ballistic_target` by 315.7 km. Note the target does *not* move — the state is planet-relative and the target is a ground point, so both are fixed in the rotating frame. What changes is the path between them.

**Bank.** A banked segment turns the heading, so the track leaves the departure arc. Measured on the HGV geometry at a 75-degree terminal bank: 166 degrees of heading change and a 21.5 km miss on a 0.6 km range residual.

In both cases the down-range solve still converges — it is solving the wrong thing well. The missing degree of freedom is the aim direction. So:

```
    solve   (psiLaunch, rangeControl)
    such that   downRangeResidual = 0   and   crossTrack = 0
```

Cross-track is already measured by both scripts, so the second residual costs nothing new.

---

## Task 1: `coorbital.util.aimSolve` — the two-axis root solve

**Files:** Create `+coorbital/+util/aimSolve.m`, `tests/test_aimSolve.m`.

**Interface produced:**
`[xSol,fAch,info] = aimSolve(fResid,x0,dx0,tol,opts)` where `fResid(x)` takes a two-vector `[psi; control]` and returns a two-vector `[downRangeResidual; crossTrack]` in metres. Returns the solved controls, the achieved residuals, and an `info` struct carrying the iteration count, the evaluation count, the Jacobian history, and a `converged` flag.

Design points that matter, all learned the hard way on this project:

- **Each residual evaluation is a full trajectory propagation** costing order 0.15 s, so the evaluation count is the cost. A damped Newton with a finite-difference Jacobian costs three evaluations per iteration. Cache and report.
- **It must NOT throw when it fails to converge.** Return `converged = false` with the best-so-far point and enough in `info` for the caller to print a useful refusal — the same contract `coorbital.util.rangeSolve` already honours, and for the same reason: a targeting script that silently returns a near miss is worse than one that refuses.
- **Guard the Jacobian.** Near a range maximum `∂R/∂control` collapses and the system becomes ill-conditioned. Detect it, report the condition number, and refuse rather than taking a wild step.
- **Bound the step.** A full Newton step from a poor start can throw the propagation into a coordinate guard. Damp it, and shrink on failure rather than aborting.
- **A failed propagation is not a failed solve.** `fResid` may throw for a control pair that produces an incomplete flight. Catch it, treat that point as infeasible, and shrink the step.

- [ ] **Step 1: Write the failing tests, against synthetic residuals with known roots**

Test the solver independently of any trajectory:
- A well-conditioned linear system converges to the analytic root within tolerance, in the expected number of iterations.
- A mildly nonlinear system converges.
- A tighter tolerance gives a tighter result.
- A system with no root in reach returns `converged = false`, does NOT throw, and carries the best-so-far point.
- A singular Jacobian is detected and reported rather than producing an infinity.
- A residual that throws for some inputs is handled by step shrinkage, not by aborting.
- The evaluation count in `info` is truthful and endpoints are not re-evaluated.

- [ ] **Step 2: Run and confirm they fail**
- [ ] **Step 3: Write `aimSolve`**
- [ ] **Step 4: Run the suite** — expect `22 passed, 0 failed`, zero warnings.
- [ ] **Step 5: Prove the tests bite**

Mutations, restoring byte-identically with md5 verified: return the first Newton step without iterating; drop the damping so a poor start diverges; report `converged = true` unconditionally. Each must FAIL.

- [ ] **Step 6: Commit**

---

## Task 2: Two-axis targeting in `HGV/run_target.m`

**Files:** Modify `HGV/run_target.m`, `tests/test_runTarget.m`.

- [ ] **Step 1: Replace the closed-form azimuth with a solved one**

Keep the great-circle bearing as the INITIAL GUESS — it is the right starting point and it is exact in the non-rotating zero-bank case, which must still converge in very few iterations and land on the same answer it does today.

- [ ] **Step 2: Remove the rotation refusal**

The script currently refuses when `earthSpin` is true. That refusal was correct when it had no way to aim off; now it does. Replace it with a solved trajectory, and state in the summary that the azimuth is now solved rather than closed-form.

- [ ] **Step 3: Remove the bank warning, or demote it**

Cross-track is now driven to zero rather than measured and warned about. Keep the measurement in the summary — it is the evidence the solve worked — but it is no longer a caveat.

- [ ] **Step 4: Verify against the cases that previously failed**

The rotating case previously missed by 231.6 km and the 75-degree banked case by 21.5 km. Both must now hit within tolerance. Report both, and report how many propagations each took.

- [ ] **Step 5: Confirm the non-rotating zero-bank case is unchanged**

3811.240 km required, about 511 m miss. If the two-axis solve moves it, say by how much and why — a solved azimuth that lands on the closed-form value to machine precision is the expected result.

- [ ] **Step 6: Tests, then mutations**

Pin the rotating case and the banked case. Prove they bite: force the azimuth back to the closed form and confirm both new cases FAIL while the old case still passes.

- [ ] **Step 7: Commit**

---

## Task 3: Two-axis targeting in `BM/run_ballistic_target.m`

**Files:** Modify `BM/run_ballistic_target.m`, `tests/test_runBallisticTarget.m`.

Same as Task 2, with one complication worth naming up front: this script already runs a nested two-parameter solve for `minimum-energy` (loft and cutoff, against range and a burnout-energy objective). Adding azimuth makes it three controls. Think about the structure before writing — an outer azimuth loop around the existing solve is the obvious arrangement and probably the right one, but say why you chose what you chose, and report the propagation count, which will be the largest in the library.

All three branch modes must still work. The lofted and depressed modes are single-parameter in range, so they become straightforward two-axis solves.

- [ ] Steps mirror Task 2. Commit separately.

---

## Task 4: The closed-loop guidance seam

**Files:** Create `+coorbital/+guide/terminalConstraint.m`, `docs/closed_loop_guidance.md`. Modify `+coorbital/+prop/phaseRun.m` only if genuinely required — and if it is, say why before doing it.

This task builds the SCAFFOLDING for PEG and VOA. **It does not implement either algorithm**; that is the next milestone. The deliverable is an interface they can attach to without the library changing again, plus a design note that says what each will need.

Read `docs/Coorbital HGV Capabilities.pptx` slides 8 and 9 first — they are the specification. VOA is a 12-variable augmented state–costate system with 7 shooting unknowns minimising burnout time; PEG is a five-parameter vector linear-tangent law with a nonlinear correction and a closed-loop cycle of navigate, solve, command, execute. Both enforce the same five terminal constraints: radius, velocity, flight-path angle, and two orbital-plane components.

- [ ] **Step 1: Establish what the current seam already supports, by experiment**

A phase's guidance is a handle `guide(t,x)` called with the live state, and `pitchProgram` already reads `x(5)`. Write a scratch guide that solves a trivial terminal problem each call and confirm it runs inside `phaseRun` unmodified. Report what worked and what did not. **Do not change the library until you know what actually needs changing.**

- [ ] **Step 2: Write `terminalConstraint`**

A struct-returning function defining the five constraints in the deck's terms, with the residual function that measures a state against them. This is the shared piece both algorithms need and neither should re-derive.

- [ ] **Step 3: Establish whether guidance needs to carry state between calls**

PEG's cycle solves `z*` and commands `p0`, then executes for `Δt` before re-solving — that is per-cycle state, not per-call. `phaseRun` calls the guide once per `ode45` stage, many times per second of flight, and the reconstruction pass calls it again. A stateful guide would therefore see a different call pattern than it expects, and the library currently documents that a guide must be pure.

Determine what a correct implementation requires. The honest options are: hold guidance state in the integrated state vector (which the seven-state convention already anticipates); or give `phaseRun` an explicit guidance-update cadence distinct from the integration steps. **Recommend one, with reasoning, and write it up. Do not build it.**

- [ ] **Step 4: Write `docs/closed_loop_guidance.md`**

What PEG and VOA are, what they need from the library, what already exists, what is missing, and the recommended path for each. Name the specific gaps. This document is the brief for the next milestone, so it must be concrete enough to implement from.

- [ ] **Step 5: Commit**

---

## Task 5: Documentation

**Files:** Modify `README.md`, `TODO.md`, `docs/DESIGN.md`, `docs/README.md`.

- [ ] Update the entry-script table and the limitations: rotating-Earth targeting and cross-range steering now work; say what the residual accuracy is and what it costs in propagations.
- [ ] `DESIGN.md` is two as-built sections behind. Add them, dated, without editing the earlier ones.
- [ ] Refresh the stale counts flagged in `TODO.md`: test count, `.m` file count, and the claim that the two LaTeX notes do not exist.
- [ ] Strike what this milestone closes; add what it opens.
- [ ] Every command printed must be run and shown to work.

---

## Out of scope

- Implementing PEG or VOA. Task 4 builds the seam and the brief; the algorithms are the next milestone.
- Multi-stage boosters, throttling, terminal homing.
- Aerothermal heating.
- Fidelity increments: J2, oblate geodetic altitude, tabulated Mach-dependent aero, US76.
- Phase 2 trajectory optimization.
