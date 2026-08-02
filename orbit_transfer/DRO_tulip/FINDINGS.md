# DRO → tulip: first results (2026-07-31)

Goal cell: **direct / min-time** — the twin of `pumpkyn.cr3bp.tfMin`. Chosen
first because DRO→tulip is the only problem in the catalog where an
independently converged *indirect* solution already exists, so it is the one
place the direct method can be checked against a second opinion.

## The transfer, characterized from the reference solution

| | |
|---|---|
| t_f | 4.015242 ND = 17.798 days |
| revolutions about the Moon | **1.18** |
| Earth distance | 0.8535 – 1.1122 (never approaches Earth) |
| closest lunar approach | **0.0164 ND ≈ 4650 km altitude** |
| mass fraction | 1.0 → 0.918692 |

**No Sundman regularization is needed**, and that is a measured conclusion, not
an assumption: at 1.18 revolutions with no Earth passage, the stiffness that
forces regularization in GTO→tulip (a ~40-rev spiral through deep perigee)
simply is not present. Plain time-domain collocation is the right tool here.

## The direct formulation

Normalized time `s ∈ [0,1]`, `dX/ds = t_f·f(X,U)`, state `[r;v;m]`, control
`[alpha;thr]`, minimize `t_f`. Dynamics mirror `tfMinEoM` line for line.

**`t_f` had to be LIFTED.** A single scalar `t_f` couples to every defect,
giving one dense KKT column — and MUMPS does not survive it. The first version
killed MATLAB with a fatal error immediately after IPOPT printed the Jacobian
structure, at only N = 400. This is the same failure the tulip solver sidesteps
by fixing `tau_f` and the earth CR3BP campaign cured with `liftDL`. Replicating
`t_f` per node with local continuity (`T(k+1) = T(k)`) makes the arrowhead
banded; the solution is unchanged and the solve takes 6.8 s.

**Throttle left free, and it saturates.** `min(thr) = 1.000000` across the
trajectory, confirming the min-time bang-at-the-bound structure rather than
assuming it. Worth noting the reference *permits* a switch —
`tfMinEoM` sets `u = 0` when `S = -||lambda_v||c/m - lambda_m > 0` — but in the
converged solution `lambda_m` runs 5.50 → 0 while staying positive, so `S < 0`
throughout.

## THE RESULT: min-time here is not unique, and the reference is not the fastest

Solutions are machine-tight (defect 1e-14 to 1e-16, terminal error 0) but depend
strongly and non-monotonically on the mesh:

| N | t_f | vs reference | closest lunar approach |
|---|---|---|---|
| 400 | 4.143868 | +3.2% | 0.01682 (4819 km) |
| **800** | **3.881410** | **−3.3%** | **0.00559 (442 km)** |
| 1600 | 7.253098 | +80.6% | 0.01645 (4674 km) |
| *reference (indirect)* | *4.015243* | — | *0.0164 (4650 km)* |

**The N = 800 solution is faster than the reference minimum.** That is only
possible if the reference is a local, not global, minimum — and the mechanism is
visible in the last column: it buys speed with a far deeper lunar gravity
assist, 442 km altitude against 4650 km.

**So the problem as posed is ill-posed.** With no minimum-altitude path
constraint, solutions can approach the Moon arbitrarily closely, and "the
minimum time" is not well defined — it is a family parameterized by flyby
altitude, bounded only by the lunar surface. The indirect method never exposed
this because shooting cannot reach those solutions from a reasonable seed.

This is **basin multiplicity in a third independent context**, after min-fuel
GTO→tulip (six optima at one t_f) and min-fuel earth (every certified row
beaten). It is not a min-fuel phenomenon.

## RESOLVED (2026-08-02): the fast solutions are quadrature error, not better optima

The open question above — whether the 442 km solution is genuinely feasible or
is exploiting the discretization — has been answered. **It is exploiting the
discretization.**

### The measurement

For each interval, start at node `X_k`, apply the control linearly interpolated
between `U_k` and `U_{k+1}`, integrate the TRUE dynamics with `ode113` at
1e-11/1e-13, and compare against `X_{k+1}`. That is the true local error the
trapezoid rule commits — a different quantity from the NLP defect, which only
says how well the returned numbers satisfy the trapezoid *equations*.

On the N = 800, 442 km solution:

| quantity | value |
|---|---|
| NLP defect (trapezoid equations) | 1.4e-14 |
| continuous residual, median | 2.1e-07 |
| **continuous residual, max** | **1.5** |
| ratio median-continuous / discrete | 1.5e+07 |
| lunar altitude at the worst interval | 504 km |

A residual of O(1) in ND units is hundreds of thousands of km. **The 442 km
trajectory is not physical.**

### It is a clean power law, so the mechanism is understood

Residual vs lunar altitude over the whole trajectory fits

    R ~ altitude^(-4.24)

spanning eight orders of magnitude, 1e-8 at 1e5 km down to O(1) at periselene.
Trapezoidal LTE goes as `h^3` times the third derivative of position, and for
two-body-like motion that derivative scales as `mu^1.5 * r^-3.5`, predicting an
exponent near -3.5 against the measured -4.24. Observed, explained, not a
one-off.

### What it means

A minimum-time solver handed a mesh that under-resolves periselene will drive
the trajectory into exactly the region where its own error is largest, because
that is where the objective can be cheaply reduced. The solver did not fail; it
succeeded on the wrong problem.

**Transferable lesson, and it applies to every campaign in this repo: a
machine-tight defect is a statement about the discretization, not about the
trajectory.** Here the two differ by 1e7. Every campaign quotes 1e-14 defects as
evidence of a good solve. That evidence is necessary and nowhere near
sufficient. The continuous residual costs one `ode113` call per interval.

### Corollary 1: the "no Sundman needed" lesson is REFUTED, not merely qualified

Residual measured on all three converged solutions:

| N | closest approach | NLP defect | R median | **R max** | R max in km |
|---|---|---|---|---|---|
| 400 | 4819 km | 3.3e-14 | 1.9e-07 | **3.2e-02** | ~12,600 |
| 800 | 441 km | 1.4e-14 | 2.1e-07 | **1.5** | ~1e5 |
| 1600 | 4673 km | 3.9e-16 | 2.5e-08 | **2.9e-03** | ~1,100 |

The middle row is the pathological deep flyby. **The other two are the
well-behaved solutions at a sane periselene — and they are inaccurate too**,
by 12,600 km and 1,100 km at the worst interval.

So the claim is not "right about the reference, wrong about the deep flyby".
It is simply wrong: **a uniform-in-time trapezoidal mesh is inadequate for this
transfer at every density tried.** The original reasoning went from trajectory
SHAPE (few revs, no Earth passage) to mesh adequacy, and that inference does not
hold. We never measured it until now.

Density alone will not rescue a uniform mesh: panel (c) shows the grid is
uniform in TIME to a ratio of 1.000 but **353x non-uniform in lunar angle swept
per interval**. Adding nodes everywhere to fix error that lives in a few
intervals is the expensive route — which is the argument for regularization.

Open: which of Sundman / periselene-concentrated mesh / higher-order collocation
is right here. Untested.

### Corollary 2: the ill-posedness claim loses its evidence

The demonstration of ill-posedness was the 442 km solution being FASTER than the
reference. That solution is invalid, so it demonstrates nothing. The three rows
of the mesh table are three discretizations disagreeing, not three optima.

What survives is a physical argument, not a computed one: a deeper flyby really
does give a stronger assist, so it remains plausible that inf t_f is attained
only at grazing, and a floor is prudent regardless. But **we have NOT
demonstrated the ill-posedness.** Doing so requires a discretization that
resolves periselene, then showing t_f still falls monotonically as the floor is
lowered. Not run.

## The minimum-altitude path constraint

Built 2026-08-02, exposed as `opts.minAltKm` in `run_dro_tulip.m` (empty by
default, so the ill-posed baseline still reproduces). Imposed in squared form,
`||r_k - r_Moon||^2 >= rho_min^2`, one row per node — avoids the square root's
undefined gradient at rho = 0 and keeps an exact, well-scaled gradient
`2(r_k - r_Moon)`.

**Is a path constraint a direct-method advantage? Yes**, and the asymmetry is
larger than it looks. Direct: four lines of CasADi, N+1 extra Jacobian rows, the
same barrier machinery already handling `0 <= u <= 1`, nothing else changes.
Indirect: guess the constrained-arc structure a priori, derive the constraint
order (q = 2 here), append an `eta(t) >= 0` multiplier active only on the arc,
impose junction conditions where the **costates jump**, and solve a multipoint
BVP whose unknowns now include entry/exit times. The direct method turns a
change of mathematical object into a change of argument list.

**The caveat is real:** the constraint binds AT NODES ONLY. Between nodes the
trajectory may dip below the floor freely. That matters here more than usual,
because the constraint is active exactly at periselene, which is where the
collocation is least accurate (above). `min_k rho_k = rho_min` means "the nodes
clear the floor", not "the trajectory clears the floor".

### A pre-registered prediction that FAILED

Stated before running: with a floor, t_f should stop depending on N and should
RISE as the floor rises. Measured at 100 km:

| N | t_f | achieved alt km | defect | status |
|---|---|---|---|---|
| 400 | 4.357303 | 1048 | 1.4e-14 | Solve_Succeeded |
| 800 | 4.724866 | 4680 | 3.1e-16 | Solve_Succeeded |
| 1600 | 24.65 | 646 | 1.2 | Max_Iterations |

and at 500 km: N=400 gives 4.143868 (4818 km, converged, identical to the
unconstrained N=400); N=800 gives 9.574640 (1940 km) but Max_Iterations at
defect 4.3e-04, i.e. NOT converged.

Two observations. **The constraint was never active** — every converged run
clears the floor by 10x or more, so its multipliers are zero. **Yet the answers
differ from the unconstrained ones** (4.144 / 3.881 / 7.253 at the same meshes).
An inactive inequality cannot move a local minimum, but it can change which one
the interior-point method walks to, because the barrier terms depend on distance
to every inequality including the inactive ones.

So the N-dependence did NOT collapse. That separates two pathologies that were
being conflated: the floor removes the *unbounded* family of ever-deeper flybys
(a genuine repair), but it does not remove the *multiplicity* of local minima.
Different diseases, different cures.

## Diagnostics and movies (stage 4 of run_dro_tulip)

Off by default. `opts.plots` gives a six-panel figure: (a) Moon-centred
trajectory with lunar disc and floor to scale, (b) altitude vs time, (c) mesh
spacing in time AND in lunar angle swept, (d) continuous-time residual against
the NLP defect, (e) throttle, (f) residual vs altitude (the power law).
`opts.movies` gives `'traj'` (Moon-centred, base MATLAB only, always runs) and
`'scene'` (the pumpkynPie-lit starfield scene via
`pumpkynPie.plot.SatelliteAnimator`, driven by OUR direct solution; falls back
to `'traj'` if the class is absent). Nothing is written into pumpkyn/pumpkynPie.
Both render at 1280x720 (multiple of 16, per the H.264 shear lesson) and draw
the sub-node path by spline — a rendering choice, not a claim that the sub-node
path is dynamically correct.

## CERTIFIED (2026-08-02): Hermite-Simpson N=1600 reproduces the indirect answer

First certified direct min-time DRO->tulip solution.

```
===== CERTIFICATION: DRO->tulip min-time (hermite-simpson, N = 1600) =====
G1  continuous accuracy (worst interval)        3.2230e-01 km  (tol 1.0)   PASS
G2  agreement with the indirect reference t_f   1.8997e-06     (tol 1e-4)  PASS
G3  NLP defect                                  3.3307e-16     (tol 1e-9)  PASS
G4  control unit-norm error                     2.2204e-16                 PASS
G5  terminal boundary error                     0.0000e+00                 PASS
G6  throttle saturation (min u)                 1.0000e+00                 PASS
    t_f = 4.0152501 ND = 17.798 days   vs indirect 4.0152425
    min node altitude 4673.3 km, worst interval at 4673 km
```

**The direct method independently reproduces the indirect t_f to 6 significant
figures**, with worst-interval true error 0.32 km. That is the cross-validation
this campaign was built for, and it had never been possible before: the earlier
trapezoid answers were off by 3-80% and inaccurate by 1,100-12,600 km.

### The ladder is NOT monotone, and the middle rung went through the Moon

Hermite-Simpson, seeded from the indirect reference, unconstrained:

| N | t_f | rel err | worst true error | min node altitude |
|---|---|---|---|---|
| 400 | 4.6808938 | 1.66e-01 | 761 km | 4818 km |
| 800 | 4.1628670 | 3.68e-02 | **1,017,917 km** | **-719.6 km (INSIDE THE MOON)** |
| 1600 | **4.0152501** | **1.90e-06** | **0.32 km** | 4673 km |

The N=800 solve returned Solve_Succeeded with an HS defect at machine precision
and **a node 719.6 km beneath the lunar surface**.

### This DEMONSTRATES the ill-posedness — the earlier retraction is itself retracted

The ill-posedness claim was retracted earlier today because its only evidence
(the 442 km trapezoid solution) was an accuracy artifact. The N=800 HS result is
different in kind: a converged solve places a node INSIDE the Moon in order to
shorten the transfer. No quadrature argument is needed to reject it, and no
quadrature argument explains it away. **The unconstrained problem genuinely has
no minimum** — nothing in the formulation stops the trajectory from passing
through the primary.

So: ill-posedness DEMONSTRATED, by a stronger piece of evidence than the one
originally offered. The concern was right; the first proof was not.

### Consequences

1. **Order alone does not fix ill-posedness.** Fourth order made the N=1600 rung
   certifiable, and made the N=800 rung dive through the Moon more decisively.
   Accuracy and well-posedness are independent problems.
2. **The altitude floor and the accuracy study are not separate work items.**
   The unconstrained problem should not be used as the accuracy testbed at all;
   its minimizing sequence is not converging to anything physical.
3. **G2 needs rethinking for the CONSTRAINED problem.** It compares against
   pumpkyn's tfMin, which has no path constraint, so there is no indirect
   reference once a floor is imposed. Proposed split: G2a transcription fidelity
   (fix t_f = t_f_ref, solve for feasibility, require the trajectory to match the
   reference to G1 accuracy) and G2b optimizer quality (free t_f with a floor,
   require a well-resolved solution). Not yet built.

### Method notes worth keeping

- **Separated, not compressed, Hermite-Simpson.** The compressed form (midpoint
  interpolant substituted into f) solves at N=400 and kills MATLAB SILENTLY at
  N=800 — deep expression graph, MUMPS dies. Same failure and cure as the t_f
  lift and liftDL: lift Xm to a variable with an explicit interpolation
  constraint.
- **A warm start should be feasible for the constraints that define it.** Seeding
  Xm with the plain midpoint average violates the interpolation constraint. NOTE:
  fixing this did NOT change the N=400 answer (4.6808938 either way) — the
  hypothesis that it explained the basin difference was WRONG. It did change
  N=800. Recorded because the reasoning was wrong even though the change was
  right.
- The residual engine has its own test (certify/tests/test_dro_residual.m, 3/3):
  8.3e-16 on an exactly-integrated trajectory, injected 1e-5 error read at 1.0e-5.

## Next

1. **Redo the accuracy ladder WITH the altitude floor on.** The unconstrained
   problem is not a valid accuracy testbed — its minimizing sequence passes
   through the Moon. Accuracy and well-posedness turned out not to be separable
   work items.
2. **N = 3200 never ran** — the ladder was stopped one rung early, right after
   N = 1600 certified. Worth finishing to confirm the order estimate and to see
   whether the deep-flyby attractor reappears at higher resolution.
3. **Split G2 for the constrained problem.** pumpkyn's `tfMin` has no path
   constraint, so once a floor is imposed there is no indirect reference to
   compare against. Proposed: G2a transcription fidelity (fix t_f = t_f_ref,
   solve for feasibility, require a match to the reference at G1 accuracy) and
   G2b optimizer quality (free t_f with a floor, require a well-resolved
   solution).
4. **Mesh refinement and Sundman are now OPTIONAL, not required.** Fourth order
   alone reached 0.32 km. Either remains worth trying as a cheaper route to the
   same accuracy at lower N, but neither is on the critical path any more.
5. **The cell this campaign was built for is now unblocked:** compare the direct
   duals against the indirect costates via the Hager covector mapping. It needed
   a trustworthy direct solution sitting on the reference, and the certified
   N = 1600 Hermite-Simpson solve is exactly that. This has never been possible
   in this repo before.
