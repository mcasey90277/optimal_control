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
| **worst-interval POSITION error** | **1441 km** |
| worst-interval VELOCITY error | 1568 m/s |
| lunar altitude at the worst interval | 504 km |

**The 442 km trajectory is not physical.**

> CORRECTED 2026-08-02 after external review. This table originally reported a
> "continuous residual, max 1.5 ND = 1e5 km". That was the norm of the FULL
> seven-component state difference — position, velocity AND mass fraction —
> multiplied by lStar and labelled km. Dimensionally invalid. The real position
> error is 1441 km: still disqualifying, but not what was claimed.

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

### Corollary 2: the ill-posedness claim loses its evidence — and lost it TWICE

The demonstration of ill-posedness was the 442 km solution being FASTER than the
reference. That solution is invalid, so it demonstrates nothing. The three rows
of the mesh table are three discretizations disagreeing, not three optima.

What survives is a physical argument, not a computed one: a deeper flyby really
does give a stronger assist, so it remains plausible that inf t_f is attained
only at grazing, and a floor is prudent regardless. But **we have NOT
demonstrated the ill-posedness.** Doing so requires a discretization that
resolves periselene, then showing t_f still falls monotonically as the floor is
lowered. Not run.

**AND THEN WE MADE THE SAME MISTAKE AGAIN.** Later the same day the claim was
re-established on the Hermite-Simpson N=800 result — a converged solve with a
node 719.6 km inside the Moon. External review pointed out that this solution
has a **3,127 km position error**, so it is not a trajectory of the continuous
problem either, and a node inside the Moon in an infeasible DISCRETE solution
proves nothing about the continuous one. Withdrawn again.

Twice, one section apart, on the same class of evidence. The lesson is worth
more than the claim: **no conclusion about the continuous problem may be drawn
from a converged NLP whose continuous residual has not been measured** — not
about optimality, not about feasibility, not about well-posedness.

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
G1   local POSITION accuracy (worst interval)      3.2946e-03 km  (tol 1.0)   PASS
G1v  local VELOCITY accuracy (worst interval)      8.4152e-04 m/s (tol 1.0)   PASS
G1b  GLOBAL POSITION accuracy (flown end to end)   4.1773e-02 km  (tol 1.0)   PASS
G1bv GLOBAL VELOCITY accuracy (end to end)         1.8498e-03 m/s (tol 1.0)   PASS
G2   agreement with the indirect reference t_f     1.8997e-06     (tol 1e-4)  PASS
G3   NLP defect                                    2.9143e-16     (tol 1e-9)  PASS
G4   control unit-norm error                       2.2204e-16                 PASS
G5   terminal boundary error                       0.0000e+00                 PASS
G6*  throttle saturation, ADVISORY (min u)         1.0000e+00                 PASS
G8   lifted-time spread                            0.0000e+00                 PASS
G9   Hermite interpolation residual                4.4409e-16                 PASS
    t_f = 4.0152501 ND = 17.798 days   vs indirect 4.0152425
    local POSITION error: max 3.3 m, sum 26.8 m, global 41.8 m (mild accumulation)
    control reconstruction: min direction norm 1.0000, throttle overshoot 8.2e-09
```

**The direct method independently reproduces the indirect t_f to 5 significant
figures**, with 3.3 m worst-interval position error and 42 m end-to-end. That is the cross-validation
this campaign was built for, and it had never been possible before: the earlier
trapezoid answers were off by 3-80% and inaccurate by 1,100-12,600 km.

### The ladder is NOT monotone, and the middle rung went through the Moon

Hermite-Simpson, seeded from the indirect reference, unconstrained:

| N | t_f | rel err | worst POSITION error | min node altitude |
|---|---|---|---|---|
| 400 | 4.6808938 | 1.66e-01 | 7.38 km | 4818 km |
| 800 | 4.1628670 | 3.68e-02 | **3127 km** | **-719.6 km (INSIDE THE MOON)** |
| 1600 | **4.0152501** | **1.90e-06** | **0.0033 km (3.3 m)** | 4673 km |

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

## COSTATE COMPARISON (2026-08-02): the direct duals ARE the indirect costates

This is what the campaign was built for. Every other transfer here has one
method's answer only, so direct-derived costates can be checked for
self-consistency but never against an independent second opinion. DRO->tulip is
the exception: `tfMinProp` returns the whole indirect costate history.

Mapping: for defect constraints D_k = 0 with multipliers nu_k, stationarity of
the Lagrangian w.r.t. the interior states IS the discrete adjoint recursion, so
nu_k -> lambda(t_k) with no h-scaling (Hager 2000). Sign and scale were
**measured, not assumed** — sign resolved against the primal control (a
dual-convention-free quantity), scale by comparing magnitudes.

### Results, N=1600 Hermite-Simpson

| quantity | value |
|---|---|
| primer from duals vs the solution's own thrust direction | **max 1.2e-06 deg** |
| direct vs indirect lambda_v direction | median 0.0006 deg, max 0.033 deg |
| direct vs indirect lambda_r direction | median 0.0006 deg, max 0.086 deg |
| **scale factor lambda_v** | **0.999991** (CoV 4.1e-05) |
| **scale factor lambda_r** | **0.999995** |
| transversality lambda_m(t_f) | 8.2e-05 (indirect: 7.1e-10) |
| Hamiltonian, median abs(lambda'f + 1) | 9.8e-06 |

Three things are stronger than expected:

1. **The primer is exact to 1.2e-06 deg.** Not approximately — the KKT
   stationarity condition for the midpoint control IS the minimum condition
   alpha = -lambda_v/||lambda_v||, so the duals reproduce it identically.
2. **The scale factor is 1.000, not merely constant.** The two methods do not
   just agree up to normalization; they land in the *same* normalization. That
   was not designed in and it is not required by anything.
3. **lambda_r matches too** — median 0.0006 deg. lambda_r couples to the
   trajectory only indirectly and is the classic weak spot of shooting (it is
   where the IFS null direction sat, at 83%), so it is the sterner test.

### The one discrepancy, run to ground — and my first explanation was wrong

Raw Hamiltonian CoV is 5.6e-02, against a condition that wants `lambda'f == -1`
identically (H_PMP = 1 + lambda'f = 0 for free t_f). The deviation is confined to
the last ~8 intervals of 1600:

| interval set | median abs(H+1) | max abs(H+1) | CoV |
|---|---|---|---|
| all | 9.80e-06 | 1.024 | 5.6e-02 |
| trim 10 each end | 8.87e-06 | 0.173 | 6.2e-03 |
| trim 20 each end | 8.09e-06 | 0.023 | 2.5e-03 |

**I first attributed this to a terminal covector-mapping artifact** — the final
node carrying the endpoint equalities and contaminating the last multipliers —
and wrote that into the source and this file. **It is wrong.**

The test that refutes it: on this transfer the CLOSEST LUNAR APPROACH occurs at
**interval 1599 of 1600**, i.e. at the very end. Boundary and close approach are
confounded, and I attributed the deviation to the boundary without separating
them. Fitting |H+1| against altitude on intervals 1..1560 gives

    |H+1| ~ altitude^-4.73

and the blocks near the end come in at 0.5-1.1x that prediction:

| block | median actual | altitude-law prediction | excess |
|---|---|---|---|
| core 1..N-40 | 8.09e-06 | 1.02e-05 | 0.8x |
| N-39..N-20 | 1.59e-03 | 2.92e-03 | 0.5x |
| N-19..N-8 | 3.56e-02 | 5.22e-02 | 0.7x |
| **LAST 8** | **7.88e-01** | **7.40e-01** | **1.1x** |

No excess anywhere. **The end-interval deviation is ordinary truncation error at
periselene**, which on this trajectory happens to sit at t_f. Trimming the ends
does not remove an artifact — it removes the hardest part of the problem.

Third time this campaign has attributed something to the wrong cause by not
separating two confounded variables. Worth naming as a habit to break: when an
effect appears "at the end", check whether anything else is also at the end.


### STATE agreement, and why the raw number is misleading

The costate comparison above is the strong claim; the state comparison is the
one a reader asks for first, and it was NOT part of the original certification —
that rested on t_f plus endpoints, which the boundary conditions pin anyway.

| quantity | value |
|---|---|
| max position difference, matched ABSOLUTE time | 3.311 km |
| max position difference, matched FRACTIONAL time | **0.396 km** |
| median, fractional | 0.230 km |
| max velocity difference | 0.341 m/s |
| max mass-fraction difference | 3.2e-10 |
| thrust direction difference | median 5.3e-04 deg, max 0.101 deg |

**The raw 3.3 km is almost entirely a phasing offset, not a shape difference.**
The two solutions differ in t_f by 7.63e-06 ND = **2.92 seconds**. At the
trajectory's speed range (0.143–1.135 km/s) that alone predicts an along-track
displacement of 0.42–3.32 km, and the observed maximum is 3.311 km — the top of
the predicted range. Re-sampling at matched fractional time s = t/t_f removes it
and leaves 0.396 km, a factor of 8.4 smaller.

So: same trajectory to ~0.4 km peak / 0.23 km median over an 89,000 km,
17.8-day transfer, with the residual difference dominated by the last few
percent of the arc where the close approach sits.

### The Hager terminal covector (applied 2026-08-02)

lambda(t_f) is now read from the TERMINAL BOUNDARY multipliers rather than
extrapolated from the interval multipliers. Stationarity w.r.t. X(:,end) gives
nu_N'*dD_N/dX_end + nu_psi'*dpsi/dX_end = 0 with dD_N/dX_end = I + O(h), so
lambda(t_f) = -nu_psi to leading order. Requires `opts.returnModel = true`,
which is what registers the boundary constraint rows.

| check | value |
|---|---|
| lambda_v(t_f) angle to the indirect terminal costate | 1.57e-02 deg |
| relative error, full 6-vector | 4.96e-04 |

Note what it does NOT fix: the Hamiltonian deviation in the last intervals,
which is truncation error at periselene (see above), not a boundary artifact.
Applying the correction was still right — it is the correct reading of
lambda(t_f) — but it was proposed on a wrong diagnosis.

### Figures

- `direct/results/dvi_N1600_state.png` — trajectory overlay; position difference
  raw and phasing-removed; velocity; mass; thrust direction; difference vs
  altitude. Every panel that rises at late time also carries lunar altitude, so
  the two can be told apart by eye — the confound that produced a wrong
  diagnosis here is now visible in the plot itself.
- `direct/results/dvi_N1600_costate.png` — lambda_r, lambda_v and lambda_m
  overlaid (direct solid, indirect dashed; they are indistinguishable, including
  the large lambda_r excursions at 13.5 and 17.8 days); direction agreement;
  scale factor; Hamiltonian condition against altitude.
- Generator: `direct/viz/plot_direct_vs_indirect.m`.

### What this settles

The earlier claim that "the two methods agree" rested on t_f and the pinned
endpoints, and had to be walked back as not being trajectory-level agreement.
**It now is.** The two methods find the same extremal — same primer field, same
costate directions, same normalization — not merely the same final time.

Code: `direct/certify/costate_compare.m`. Figure:
`direct/results/costate_compare_N1600.png`. Data:
`direct/results/hs_N1600_duals.mat` (requires `opts.returnModel = true`, which
is what builds the constraint-row registry that locates the defect multipliers).

## COLD-START TEST (2026-08-03): accuracy survives, basin selection does not

Same solver, same Hermite-Simpson scheme, same meshes — but the crude internal
seed (states linear between endpoints, mass ramp 1→0.92, thrust along the chord,
tf0 = 4.0; zero content from the indirect solve, duals unseeded):

| N | t_f | vs ref | worst POS error | min node alt | G1 | G2 |
|---|---|---|---|---|---|---|
| 400 | 4.3806738 | +9.1% | 5.23 km | 4819 km | fail | fail |
| 800 | 4.4506628 | +10.8% | 6954 km | **-343 km (inside Moon)** | fail | fail |
| 1600 | 4.7832984 | **+19.1%** | **0.0078 km** | 4680 km | **PASS** | fail |

Reference (indirect, and the warm-started certified direct): t_f = 4.0152.

### What each row says

- **N=1600 is the decisive one.** The cold solve produced a genuinely accurate,
  physical extremal — 7.8 m worst-interval position error, G1 PASS, safe
  periselene — that is **19% slower than the reference**. The accuracy machinery
  works cold; **basin selection does not.** The direct method cold-finds *a*
  minimum-time extremal, not *the best known one*.
- **N=800 dove through the Moon from a cold seed too** (-343 km vs the
  warm-started run's -719.6 km). The N=800 mesh finding the unconstrained
  problem's hole is a property of the discretization + formulation, not of the
  seed. Third independent instance.
- **Three meshes, three different basins** (4.38 / 4.45-invalid / 4.78). Mesh
  density acts as a de facto random seed for basin selection.

### Consequences

1. **Seeding or continuation is required equipment on catalog pairs, not a
   convenience.** A cold direct solve certifies on accuracy while silently
   leaving 19% of the objective on the table — and nothing in the solver output
   distinguishes that from success. Only G2 (an independent reference) exposed
   it, and catalog pairs will not have one.
2. **This sharpens what the certified agreement means.** Warm-started from the
   indirect neighborhood, the direct method refines to the same extremal to
   0.4 km / 6e-4 deg. Cold, it does not find that neighborhood. The two-method
   agreement is a statement about refinement, not about global search.
3. **The indirect answer is the best known on this problem** — faster than
   everything the cold direct method found at any mesh. Worth remembering when
   tempted to treat the direct method as the global-search half of the
   partnership: on this problem neither method searches globally; the indirect
   one simply arrived with better converged costates (Darin's walk-down).
4. Caveat recorded: tf0 = 4.0 is a round number but was chosen knowing the
   answer is ~4. A tf0 sweep (2/3/5/6/8) would test how much that one scalar
   steers the basin; not yet run.

Data: `direct/results/cold_hs_N{400,800,1600}.mat`.

## THE PHASING MAP, COMPLETE (2026-08-04): 132/144 cells, 11.3 hours

The 12x12 torus sweep finished: **132 of 144 phasing pairs solved and
continuously verified** (flown-control gate <100 km / 10 m/s at N=800), 8 cold
openers + 124 warm-flood conversions, 262 total attempts, 676 min wall.

### The headline numbers

- **Global minimum of the map: t_f = 3.6566 ND = 16.21 days at
  (s_D, s_A) = (0.083, 0.409)**, dV 0.679 km/s, periselene 6,469 km, flown miss
  0.89 km. That is 8.9% faster than the indirect reference (4.0152) and 4.2%
  faster than the certified 3.8170 -- pending its own certification ladder.
- Two cost valleys: s_A ~ 0.41 and s_A ~ 0.66 (t_f 3.66-3.79 across many
  departure phases); ridges at s_A ~ 0.16 (5.7-5.9). Departure phase matters
  WEAKLY almost everywhere -- rows are near-uniform in s_D -- consistent with
  the departure axis being dynamically benign.
- **The 12 unsolved cells form exactly ONE ROW: s_A = 0.075** -- the demo's own
  arrival phase, the fastest point on the tulip (vArr = 1.119 km/s). Every cold
  seed and every warm start from adjacent rows failed all 4 tries there. The
  deepest crease of the map is precisely the slot the original heuristic chose.
  And we know solutions EXIST there: the two certified transfers (4.0152,
  3.8170) live at (0, 0.075) -- found earlier via the indirect seed and the
  floor-chain, routes the sweep does not use. Red means "unreachable by this
  sweep's seeds", not "no solution".
- 3 fast-arrival cells (vArr 0.74) DID open warm -- cold-seed difficulty was
  seed-conditioned, exactly as Mike suspected; only the vArr 1.12 column
  resisted everything.
- The costate catalog now holds **132 entries** (LAM0 in dsweep_12x12.mat):
  sign-resolved [lambda0; tf] per cell via the validated covector mapping.

### What the waves measured

- wave 0 (chord seed, full budget, 6x6 sublattice): 8/36 = 22% -- matches the
  transect's cold statistics.
- wave 1 (multi-source trajectory flood, half-size steps, <=4 tries/cell):
  124 conversions, most in 5-30 s; conversions kept landing deep into the
  retry tail (4th-neighbor attempts).
- Engineering ledger for production sweeps: FOUR budget layers were needed
  (solver iterations, solver CPU via IPOPT max_cpu_time, per-edge wall, and
  never-fly-garbage verification screening) -- each because a different stack
  layer could stall unboundedly.

### Follow-ups queued

1. Certification ladder on the map minimum (3.6566) and spot-winners.
2. The red row: solve its cells from the two existing certified solutions as
   seeds (the route that works there), completing the torus.
3. Re-solve the two suboptimal-basin outliers (7.51 at (0.417,0.409), 6.21)
   from their best neighbors -- best-per-cell keeps improvements.
4. Seed-sensitivity pass (dual-coast) now optional: only one row needs it.
5. The period axis (Darin's third dimension).

Figures: results/phase_torus_12x12_{torus.fig,torus.png,flat.png,tfmap.png}.
Data: results/dsweep_12x12.mat (map + catalog). Clean front door:
run_phase_sweep_ps.m; industrial: sweep_phasing_direct.m.

## CERTIFIED (2026-08-03): a flyable extremal 4.9% FASTER than the indirect reference

**t_f = 3.8169913 ND = 16.919 days, periselene 3,954 km, dV 0.7100 km/s** —
against the indirect reference's 4.0152425 / 4,673 km / 0.7485 km/s. All eleven
gates pass at Sundman N=6400: local error 4e-8 ND, **global single-shot
0.457 km / 0.045 m/s**, floor honoured between nodes, terminal error 2e-14.

**The indirect reference is no longer the best known flyable solution to this
problem.** Same endpoints, same physics, comfortable altitude — 4.9% faster and
5.1% cheaper. (Certification means transcription-accurate LOCAL extremal, as
always; global optimality is claimed by nobody.)

### Provenance — the discovery chain matters

1. The 500 km floor experiment's cold N=400 solve found a floor-riding basin at
   t_f≈3.70 — unresolved (983 km error), but pointing at real physics.
2. Sundman-regularized re-solve seeded from it converged to t_f=3.8169912 at
   N=800 — and drifted OFF the floor to 3,954 km. The deep flyby was never
   needed; the solver used the floor solution as a stepping stone to a basin
   the reference-seeded and cold solves never found.
3. Refinement N=800→1600→3200→6400: t_f stable to 8 digits throughout; global
   miss 29.4 → 7.3 → 1.8 → 0.457 km (×4 per doubling — discretization, not
   intrinsic sensitivity; the ~1500x error amplification along this trajectory
   is real but the input error shrinks faster).

**The seed lineage is entirely cold.** Step 1 was the crude straight-line seed
(linear states, chord thrust, tf0=4.0); every later step was seeded from its
predecessor. No indirect state, costate, or t_f entered the chain at any point
— the direct method beat the indirect method's converged answer with no help
from it.

So: the floor experiment discovered the basin, Sundman resolved it, the
refinement ladder certified it, and the G1b global gate — added only after
external review — is what forced the honesty at each step.

### Immediate consequences

- **For Darin's phasing map:** at the demo phasing itself, the map value is at
  most 3.8170, not 4.0152. Single indirect solves under-report the family even
  at their own anchor point.
- **Next: close the loop.** Map this solve's duals to costates, seed tfMin,
  and obtain the indirect twin of the NEW basin — the direct->indirect handoff
  in production for the first time on an answer the indirect method did not
  already have. That entry (pair, thrust, basin) goes into the catalog.
- The Sundman option (opts.sundman, kappa = rho_Moon^1.5, time as 8th state,
  auto re-sampling of time-uniform seeds) is now the recommended mode for
  anything that goes near the Moon: it reproduced the reference to 7 digits at
  HALF the mesh (N=800, 0.16 m local) before finding this.

Data: `direct/results/sundman_floor500_N{800,1600,3200,6400}.mat`.

## THE FLOOR-COST CURVE, COMPLETE (2026-08-03): cheap until the arrival geometry says no

Floor-ramp continuation from the 3.817 basin (Sundman N=800, every rung
constraint-ACTIVE, POS errors ~0.2 m):

| floor km | t_f | cost vs unconstrained |
|---|---|---|
| (natural 3954) | 3.8169913 | — |
| 4000 | 3.8169992 | +3 s |
| 4500 | 3.8180227 | +6.6 min |
| 4600 | 3.8184141 | +9.1 min |
| 4700+ | — | UNREACHABLE (3 solver configs failed) |

**Safety margin is nearly free** — until it abruptly isn't. The wall has two
layers, both set by the ARRIVAL GEOMETRY, not by the solver:

1. **Hard ceiling = the arrival endpoint's own altitude: 4,818.5 km.** The
   arrival state on the tulip is a fixed boundary condition; any floor above it
   makes the problem infeasible by construction. The 6000/8000/10000 km floor
   failures were CORRECT infeasibility detections — misdiagnosed at the time as
   restoration failures from bad seeds. (Second wrong diagnosis of the day
   corrected; the first was the "terminal covector artifact".)
2. **Practical ceiling ≈ 4,6xx km — the terminal-approach squeeze.** Every
   converged solution's closest approach occurs just before arrival: the
   natural approach dips ~140-220 km below the endpoint altitude before
   landing on it. At floor 4600 (218 km of dip allowance) continuation works in
   seconds; at 4700 (118 km) it fails under default barrier, reference-seeded,
   and mu_init=1e-6 warm-started configurations alike.

**The tie-in to the phasing sweep is the punchline:** both ceilings are
functions of the ARRIVAL PHASE. Pick a different arrival point on the tulip and
the endpoint altitude — hence the maximum enforceable floor — moves. "How much
lunar clearance can this transfer guarantee?" is a question about arrival
phasing, which is exactly the axis Darin's sweep varies. The sweep's periselene
channel should therefore also record the ARRIVAL-POINT altitude per phasing:
it is the feasibility ceiling for any altitude policy.

Data: floor_ramp_sundman_N800.mat, rampA_f4600.mat; failures logged in
scratch (floorhi/ramp2/probe4700).

## THE FLOOR EXPERIMENT (2026-08-03): 500 km, warm and cold, all meshes

Floor chosen far below the reference periselene (4,673 km) so it should be
INACTIVE in the reference basin — it exists to close the through-the-Moon hole.
G2 auto-disables under a floor (tfMin has no path-constraint capability).

| seed | N | t_f | POS err km | node alt | TRUE alt (between nodes) | status |
|---|---|---|---|---|---|---|
| warm | 400 | 4.6808938 | 7.4 | 4819 | 4674 | converged |
| warm | 800 | 4.6850952 | 0.23 | 4678 | 4673 | converged |
| warm | 1600 | 4.6853523 | **0.0071** | 4678 | 4673 | converged |
| cold | 400 | 3.6979385 | 983 | 577 | 506 | converged, unresolved |
| cold | 800 | 4.3884674 | 119 | **500** | **380 — VIOLATES floor** | converged, unresolved |
| cold | 1600 | 18.83 | 53,500 | 2220 | −1559 | **Max_Iterations** |

### Four findings

1. **The floor closes the hole.** Warm N=800, which unconstrained went 719.6 km
   inside the Moon, lands in the legitimate 4.685 basin.
2. **An INACTIVE floor still ejects the solver from the reference basin.** Warm
   N=1600 — same seed (the indirect trajectory itself), same mesh that
   certified unconstrained at t_f=4.0152501 — converged to 4.6853523 with the
   floor 4,178 km slack. Interior-point barrier terms involve distance to every
   inequality, active or not; even starting ON the reference trajectory did not
   hold the basin. Reproduces the 100 km trapezoid observation at the
   certifying mesh with the 4th-order scheme.
3. **The constrained optimum wants to ride the floor, and undercuts the
   reference.** Cold N=400 found t_f=3.698 (−8% vs the unconstrained reference)
   riding at 506 km — physics, not artifact: bounded family, minimizer on the
   boundary. But no mesh resolved it: 983 km error at N=400; at N=800 the nodes
   sit exactly at 500 km while the TRUE trajectory dips to 380 km — **the
   node-only enforcement caveat measured, a 120 km violation**; at N=1600 the
   solver could not converge at all in 8000 iterations.
4. **Bottom line: the constrained problem has NO certified solution.** Warm
   rows are well-resolved but in the wrong basin; cold rows chase the right
   (floor-riding) structure but cannot resolve it on a uniform mesh.

### What would fix it (untried)

- Reference basin under the floor: warm-start from the CONVERGED UNCONSTRAINED
  N=1600 solution (already feasible for the discrete dynamics) and/or IPOPT
  warm-start options (`warm_start_init_point yes`, small `mu_init`) so the
  barrier does not restart at 0.1 and shove the iterate.
- Floor-riding branch: this is where the periselene-concentrated mesh or
  Sundman regularization — optional for the unconstrained problem — becomes
  NECESSARY. The active-constraint arc sits exactly where a uniform mesh is
  weakest. Until then, "min-time with a 500 km floor" has no trustworthy value;
  we only know it is somewhere at or below ~3.7 ND if the floor-riding branch
  is real, and that the reference basin (4.685 under the floor barrier path) is
  an upper bound.

Data: `direct/results/floor500_{warm,cold}_N{400,800,1600}.mat`.

## The pumpkyn-style companion script (2026-08-03)

`direct/run_dro_tulip_ps.m` — one straight-line script in the style of
`demos/lowThrustDRO2Tulip.m`: same constants, same endpoint construction, same
figures. It runs the indirect solve for the reference, the direct HS N=1600
solve warm-started from it, prints the side-by-side comparison, and renders the
showMoon scene with both trajectories overlaid. Verified: reproduces the demo's
own propellant (12.196 kg) and dV (0.7485 km/s) exactly, and the comparison
numbers above (2.92 s, 0.396 km, 0.0056 m/s at matched s). The heavy machinery
(sweeps, gates, movies) stays in `run_dro_tulip.m`.

## Next

0. **Untried suggestions from the 2026-08-02 external review, recorded so they
   are not lost:** replace `opti.minimize(TF(1))` with `minimize(mean(TF))` so
   the cost gradient is spread across the lifted copies instead of loaded onto
   one (Gemini: better KKT conditioning; not tried, since changing the objective
   risks perturbing a certified result); and add a Betts-style polynomial
   residual alongside the re-integration, which is far cheaper and mirrors the
   mesh-refinement math, as a complement rather than a replacement.

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
   alone reached 3.3 m worst-interval position error. Either remains worth
   trying as a cheaper route to the same accuracy at lower N, but neither is on
   the critical path any more.
5. ~~The cell this campaign was built for~~ **DONE 2026-08-02 — see the costate
   comparison above.** Follow-on: apply the Hager TERMINAL correction so the last
   ~1% of intervals become informative, and port `costate_compare` to the other
   campaigns as a costate-catalog seeding check (it is the mechanism
   [[goal-costate-catalog]] needs: a direct solve now demonstrably produces
   costates good to 1e-3 deg, which is far better than any hand-built guess).

## MIN-ENERGY PILOT (2026-08-14): the fixed-t_f pipeline works end to end, 5/5

The first non-min-time run of the costate pipeline. Problem: fixed-t_f
minimum energy, J = ∫ s² dt (Bertrand–Epenoy ε = 1 — the same convention
GTO_tulip's energy→fuel homotopy starts from), throttle s ∈ [0,1], primer
direction, final mass free; t_f = γ · t_f^min(cell). PMP: same 14-state
field as min-time except the throttle law s* = clip((T/2)(‖λ_v‖/m + λ_m/c),
0, 1); terminal r, v matched + λ_m(t_f) = 0 (seven for seven); H is a
first integral, not zero.

Machinery added (all TDD, all tests green, golden cells bitwise unchanged):
`costate_common/ms_bvp` `opts.fixedTf` (drops the t_f unknown — a real
structural switch, not a guard) plus a Newton polish after an early fsolve
exit (its ‖JᵀR‖ < 1e-6 test fires with ‖R‖ ~ 4e-10 on short arcs; disabling
it instead made fsolve grind 1 → 7 iterations at the 1e-13 floor on the
golden DRO cell — the quality regression caught that); the min-energy field
`cr3bp_minenergy_pmp` (+ `_prop`) with exact CasADi-AD Jacobian, equal to
`pumpkyn.cr3bp.tfMinEoM` at saturation to 1.7e-13 (F) / 2.7e-13 rel (A);
`casadi_mintime_dro` `objective='energy'` + `tfFix` (default path bitwise
identical, trapezoid and HS); the binding `indirect/ms_minenergy`; and the
generic single-shooting acceptance gate `costate_common/ss_bvp_accept`
(K = 1 on the same closures — the role tfMin plays for min-time entries).

**Pilot** (`run_minenergy_pilot`, records in `direct/results/minenergy_pilot.mat`,
figures `minenergy_<iD>_<iA>_g<γ>.png`): three flagship 12×12 cells at
γ = 1.2 and a γ ladder on the golden cell. Warm start = the min-time cell's
trajectory stretched to γ t_f. Direct solves 6–150 s (N = 800 Sundman HS);
ms K = 12 converged every time in 2–3 iterations, 1–2 s.

| cell | γ | t_f [d] | J = ∫s² | m_f (energy) | m_f (min-time) | thr min | ‖R‖ ms | accept |Δz| | flown dir. [km] | indirect arrival [km] | J rel diff |
|---|---|---|---|---|---|---|---|---|---|---|---|
| (2,5) | 1.10 | 17.83 | 2.4724 | 0.9379 | 0.9260 | 0.243 | 2.6e-11 | 0.0 | 1.05 | 0.000 | 1.1e-5 |
| (2,5) | 1.20 | 19.45 | 1.9877 | 0.9421 | 0.9260 | 0.216 | 6.9e-12 | 0.0 | 1.10 | 0.000 | 2.7e-5 |
| (2,5) | 1.40 | 22.69 | 2.0452 | 0.9374 | 0.9260 | 0.192 | 3.6e-11 | 1.5e-10 | 3.77 | 0.000 | 7.5e-6 |
| (6,8) | 1.20 | 19.91 | 2.2828 | 0.9381 | 0.9242 | 0.378 | 6.5e-11 | 3.1e-10 | 12.18 | 0.031 | 2.8e-5 |
| (1,2) | 1.20 | 30.24 | 1.8160 | 0.9331 | 0.8849 | 0.118 | 4.4e-11 | 0.0 | 40.20 | 0.053 | 2.1e-5 |

All seven gates pass on all five: G1b flown direct control (< 100 km / 10
m/s), ms converged, single-shooting acceptance, |ΔH| along the indirect
flight 3e-9..4e-8 (integrator level), indirect flight lands, direct-vs-
indirect J agree to ~1e-5 (collocation-order), throttle interior. The direct
node throttle and the indirect flight's s(t) overlay to the eye
(`minenergy_2_5_g1.20.png`); saturation plateaus sit at periselene where
λ_r spikes. Final masses agree between routes to six digits.

### Two lessons that changed the gate definitions mid-pilot

1. **The single-shooting residual has a floor, and it is not ms_bvp's.** The
   same λ₀ propagated over the full arc WITH and WITHOUT the variational
   equations lands 6e-7 apart (ode113 RelTol 1e-10 × STM growth over 4–7
   ND). So "returned unchanged" must be judged at that floor: ss_bvp_accept
   now defaults tolR = 1e-6 (documented as the floor), verdict = |Δz| < 1e-6
   at that tolerance — which is also how tfMin's own gate reads. Before the
   change two cells "failed" acceptance with |Δz| = 3e-10 and 2e-11.
2. **Judge H conservation absolutely, not relatively.** Cell (1,2) has
   H = 0.093, so a 4e-8 absolute drift read as 4.6e-7 relative and "failed"
   the 1e-8 relative gate. The gate is now |ΔH| < 1e-6 absolute (measured
   3e-9..4e-8) — a PMP-shape sanity check, not a precision certificate; the
   precision certificates are the ms residual and the flown arrival.

### Observations worth carrying forward

- **J and m_f are not monotone in γ** on cell (2,5): m_f 0.9379 → 0.9421 →
  0.9374 and J 2.47 → 1.99 → 2.05 for γ = 1.1 → 1.2 → 1.4. Fixed rotating-
  frame endpoints mean the fixed-time problems at different γ are not nested
  (you cannot append a coast), and the γ = 1.4 solve, warm-started from the
  stretched min-time arc, may sit in a different family. A γ grid per cell —
  the "min-fuel = t_f-grid convergence map" lesson from GTO_tulip — will be
  needed before any min-energy catalog axis is declared; basin discipline
  applies here too.
- The min-energy m_f exceeds the min-time m_f by 1.2–4.8 % of m₀ at
  γ = 1.2 — a first quantitative "time-vs-fuel" number on this pair without
  yet solving min-fuel.
- Harvested collocation costates seed ms to 2–3 iterations at K = 12 for
  every cell; no K escalation was needed (the min-time torus needed 12→24→48
  on some cells). The smooth energy problem is the easy end, as predicted.

### What is NOT yet done (deliberately)

No second-order verdict for min-energy entries: `ms_conjugate_test` is the
free-time quotiented Jacobi test and does not apply verbatim (no flow
column; the λ scaling invariance is broken by L = s²). No catalog schema
axis for γ. No energy→fuel homotopy yet — that is the next step toward
min-fuel entries and reuses these seeds directly (`GTO_tulip`'s ε walk).

## 17. Deep-rung probe: the closure wall is a basin wall at 90 mN, not a sensitivity wall (2026-09-01)

Single-cell continuation walk (`probe_deep_rungs.m`) from the fine sheet's
fastest 0.5 N entry — cell (1,11), t_f = 3.18 d — down a 0.75-ratio rung
schedule toward 25 mN, same engine and gates as the catalogs, hard-capped
per call via `costate_common/run_capped`. The roadmap's "probe 0.1 N and
25 mN on one cell" item (§5A / §6 step 4).

**Result: six rungs closed, wall located between 0.09 and 0.067 N.**

| T (N) | t_f (d) | swept revs | perilune passes | ms normR | tfMin accept |
|---|---|---|---|---|---|
| 0.375 | 4.05 | 0.44 | 0 | 2e-13 | dz = 0 |
| 0.280 | 5.42 | 0.58 | 0 | 6e-14 | dz = 0 |
| 0.210 | 7.26 | 0.79 | 0 | 2e-14 | dz = 0 |
| 0.160 | 9.55 | 0.92 | 0 | 4e-13 | dz = 0 |
| 0.120 | 11.57 | 1.26 | 1 | 8e-14 | dz = 0 |
| 0.090 | 13.37 | 1.28 | 1 | 7e-14 | dz = 0 |
| 0.067 | — | — | — | stalls 0.8–2.4 | — |
| 0.050 | — | — | — | stalls 1.3–6.1 | — |

Three findings:

1. **The 0.1 N roadmap target CLOSES** (bracketed by clean 0.12 and 0.09 N
   entries); the 25 mN target does NOT close by single-step continuation.
2. **The wall is a BASIN wall, not a sensitivity wall.** Even at 90 mN the
   transfer sweeps only ~1.3 revolutions (petal-to-petal geometry, no
   spiral) — so pumpkyn tfMin single shooting accepted every closed rung at
   |dz| = 0 exactly, the regime where the GTO flagship's ~40-rev
   identifiability problem never appears. At 0.067 N every guess stalls at
   normR ~ O(1) (not divergence, not conditioning) — the stretched
   low-winding seed no longer matches a solution that presumably needs more
   winding. This is the free-time cousin of the fixed-τf tulip
   topology wall (ladder-prep P2, 2026-07): continuation cannot GROW
   winding from a topologically short seed.
3. **Junction states are banked** for all six closed rungs
   (`indirect/results/probe_deep_rungs.mat`, `R.Y` = full K+1 ms junctions
   per the identifiability rule) — ready seeds for any future deep-rung
   campaign.

**Follow-up route to 25 mN** (not attempted, recorded): a winding-aware
continuation — seed the sub-67 mN regime from a HIGHER-winding family
member (e.g. a longer-t_f branch at 0.09 N, or the multi-rev Lambert-style
initializers), or walk t_f upward at fixed thrust before descending
further. Finer thrust steps alone are unlikely to help across a topology
boundary.

## 18. The energy->fuel race: eps ladder wins by knockout; Huber-in-throttle is structurally unsuited (2026-09-02)

The pinned PLQ experiment (`run_minfuel_race`, cell (2,5) at gamma = 1.2,
both arms from the same min-energy seed, fixed-tf ms_minfuel + pilot gates,
pre-registered scoring).

**Arm A — Bertrand-Epenoy eps ladder: 1 -> 0.0017 in 15 clean steps** (2-5
Newton iterations, ~3 s each; bisection only below 0.002). m_f rose
MONOTONICALLY 0.942108 -> 0.947046 with zero basin flips, Hdrift ~1e-12
throughout, coast fraction growing 0.17 -> 0.33 as the bang structure
emerged, and the acceptance gate took the deepest solution at |dz| = 0.
m_f converged to six decimals by eps ~ 0.005: **the min-fuel answer is
m_f = 0.94705** (vs 0.94211 min-energy — the fuel saving is 0.49% of m0
at this cell, "free" for the same t_f).

**Arm B — PLQ Huber kappa walk: could not leave kappa ~ 1.** kappa = 1
reproduced the energy solution exactly (m_f 0.942108 — as it must: L =
s^2/2 is a positive scaling of the energy objective), and then EVERY step
below kappa ~ 0.91 failed — not by divergence but by a RESIDUAL FLOOR:
normR stalls at 4e-10..1e-6, just above tolR = 1e-10, exactly the
signature of the recorded structural defect (cr3bp_minfuel_pmp header):
for kappa < 1 the Huber throttle law JUMPS from kappa to 1 wherever an arc
crosses Q = 1, and the discontinuous field caps the achievable shooting
residual; the floor RISES as kappa falls (jump size 1 - kappa). Deepest
converged: kappa = 0.915. Loosening tolR would only postpone the wall.

**Verdict: the eps ladder is the production continuation family.** The
PLQ-in-throttle embedding is refuted with mechanism (switch-jump residual
floor + no-coast property), not just outscored — a clean negative result
for the experiment. PLQ penalties remain interesting for OTHER embeddings
(soft path constraints, direct-transcription objectives) where the control
box does not interact with an affine tail.

Race data (all steps, junction states, both arms):
`direct/results/minfuel_race.mat`.

## 19. The first min-fuel record set: 7/7 backbone records walked to eps ~ 0.001-0.005 (2026-09-02)

`run_minfuel_grid` (the race-winning eps ladder over every passing
min-energy backbone record; full ms junction states saved per the
identifiability rule; `direct/results/minfuel_grid.mat` + per-record mats).

| cell | gamma | eps reached | m_f fuel | m_f energy | gain [% m0] | coast |
|---|---|---|---|---|---|---|
| (2,5) | 1.10 | 0.0018 | 0.941107 | 0.937899 | +0.32 | 0.17 |
| (2,5) | 1.20 | 0.0017 | 0.947046 | 0.942108 | +0.49 | 0.33 |
| (2,5) | 1.40 | 0.0050 | 0.944025 | 0.937432 | +0.66 | 0.42 |
| (6,8) | 1.10 | 0.0012 | 0.928082 | 0.926771 | +0.13 | 0.08 |
| (6,8) | 1.20 | 0.0013 | 0.943936 | 0.938120 | +0.58 | 0.33 |
| (1,2) | 1.10 | 0.0034 | 0.901113 | 0.896934 | +0.42 | 0.25 |
| (1,2) | 1.20 | 0.0010 | 0.942523 | 0.933072 | +0.95 | 0.58 |

Findings:

1. **The fuel gain grows with gamma** (+0.13% of m0 at the tightest
   (6,8)/1.1 to +0.95% at (1,2)/1.2 with 58% coast): more time buys more
   coast buys more propellant, quantified per cell for the first time.
   Note m_f is NOT monotone in gamma within a cell ((2,5): 0.9411 ->
   0.9470 -> 0.9440) -- the gamma-basin structure survives into min-fuel.
2. **Single-shooting acceptance is NOT a valid gate at deep eps.** Only
   the 19.4-day (2,5)/1.2 passed (|dz| = 0); one more was a floor artifact
   ((2,5)/1.1 moved 4.4e-7 < tolDz but the ss residual could not reach
   1e-6), and the rest moved 2e-5..3.6 -- single shooting over 20-30 days
   of near-bang dynamics wanders, the min-time "shooting dies at depth"
   lesson reappearing in min-fuel. Entries are certified by ms convergence
   (1e-10), absolute H conservation, and the enforced endpoint match; the
   independent second-order verdict is Task 4's fixed-tf conjugate test.
3. The walk itself is robust: 2-10 endgame bisection failures per record,
   all in the eps < 0.005 sliver where m_f is already converged to ~1e-5.

## 20. Fixed-tf conjugate test: 13/14 verdicts pass; the one refutation explains the gamma anomaly (2026-09-02)

The fixed-final-time Jacobi test needed NO new instrument: the existing
ms_conjugate_test with the right spec -- freeTime = false (no flow
column), quotientDir = [] (the running cost breaks the min-time scaling
invariance), rows [1:6, 14] (the components vanishing under the terminal
conditions: r, v fixed, lam_m(tf) = 0), cols 8:14 (all seven initial
costates; lam_m is NOT degenerate here, the throttle law depends on it).
The monitored det equals the single-shooting BVP Jacobian at tf exactly.
Validated on the analytic LQ pi-conjugate case
(tests/test_conj_fixedtf, 5/5: crossings detected and bracketing pi at
junction resolution); wired into ms_minenergy and ms_minfuel as
opts.conjTest.

Verdicts (conj_fixedtf_verdicts.mat; every re-solve converged at
normR ~ 1e-11): all 7 MIN-ENERGY backbone records PASS; 6 of 7 MIN-FUEL
records PASS; **min-fuel (2,5)/gamma=1.4 (eps = 0.005) is REFUTED -- one
interior conjugate point.**

The refutation closes an open loop: that record was the grid's gamma
anomaly (m_f = 0.9440 at gamma = 1.4, WORSE than 0.9470 at gamma = 1.2
despite more time, and the shallowest eps its endgame reached). The
second-order test now identifies it as a non-minimizing extremal -- the
eps walk at gamma = 1.4 drifted into a weaker basin, consistent with the
gamma-basin structure measured on the energy side ((6,8) gate-1 splits).
Consequence for the future min-fuel catalog: the fixed-tf conjugate
verdict is a PRODUCTION gate (it caught exactly the entry a consumer
should not fly), and refuted cells should be re-walked from a different
gamma neighbor before packaging.

## 21. Gamma-continuation rescues the refuted record; schema v3 ships with the first min-fuel catalog (2026-09-02)

**The re-walk (production rule, first application).** A direct family jump
failed instructively (deep-fuel junctions stretched into an eps = 0.25
solve: too far in BOTH smoothing and time). The winning recipe is
GAMMA-CONTINUATION AT FIXED DEEP EPS: walk the healthy (2,5)@1.2 solution
1.2 -> 1.225 -> 1.25 -> 1.3 -> 1.35 -> 1.375 -> 1.4 (two bisections),
staying on the good fuel branch and moving only the time axis. Result at
gamma = 1.4: **m_f = 0.949005, conjugate test PASS** (0 crossings,
normR 1.3e-11) -- +0.50% of m0 over the refuted branch (0.944025), and
gamma-monotonicity RESTORED (0.9411 -> 0.9470 -> 0.9490): the anomaly was
the basin, not the physics. Continuation recipe note: eps-then-gamma
ordering matters -- descend eps once on a good branch, then move gamma at
fixed depth; do not re-descend eps from every gamma's own energy seed.

**Schema v3 + the first v3 catalog.** catalog_schema gains version 3
(objective/gamma axis): one catalog per objective; named .axis3 replaces
rungs_N; sheets carry tfmin_nd, p_floor, STORED mf_frac (coasts break the
all-burn identity; dV via the new deltav_from_mf derivation), lam0 [7 x n],
and -- REQUIRED for minfuel -- .Yj [14 x K+1 x n] junction states
(identifiability rule). Validator TDD 12/12 with all five shipped v1/v2
catalogs still validating clean (compat preserved). First v3 artifact:
`costate_catalog_dro_tulip_minfuel.mat` (7 entries, gamma {1.1,1.2,1.4},
fixed-tf conjugate verdicts inside incl. the rescued (2,5)@1.4;
build_minfuel_catalog.m is the packager). Step 5 is now CLOSED end to end.

## 22. Review fixes: the Huber "knockout" was our Jacobian; the conjugate refutation was a coast artifact; the catalog was mislabelled (2026-09-05)

A three-way external code review of the step-5 line (GPT-5.6-sol + GPT-6
Astra + host; `reviews/minfuel_code_review_2026-09-05.md`) found four P0
defects. All four are fixed, TDD (four new failing tests written first, all
green), the golden cells still 20/20, and the verdicts re-swept. Three of
the four change what sections 18, 20 and 21 claim.

**P0.2 -- section 18 RETRACTED in its mechanism.** The Huber throttle jumps
at Q = 1; `cr3bp_minfuel_prop` integrated Phi_dot = A Phi straight through
the jump, omitting the saltation update
Phi+ = [I + (F+ - F-) n'/(n'F-)] Phi-, n = grad Q. Measured against finite
differences on a one-switch arc: STM error 1.35e-3 (eps control case:
1.35e-7). With event-split propagation and the saltation matrix
(`tests/test_huber_saltation`: 2.1e-6, the FD floor) the SAME race cell
(2,5)@1.2 gives:

| arm | rungs | p deepest | m_f | fails | bisects | wall |
|---|---|---|---|---|---|---|
| eps (09-02) | 17 | 0.00168 | 0.947046 | 9 | 7 | 7.5 min |
| huber, no saltation (09-02) | 2 | 0.91 | 0.942176 | 8 | 6 | 3.7 min |
| **huber, saltation (09-05)** | **17** | **0.001** | **0.947041** | **0** | **0** | **1.2 min** |

Huber walks the full ladder with ZERO failures and lands within 4.5e-6 of
the eps mass. The "switch-jump residual floor" was Newton fighting a wrong
Jacobian; "Huber-in-throttle is structurally unsuited" is withdrawn. What
survives of section 18: the eps family has exact coast arcs at every
finite p and Huber does not (true, and visible: Huber's coast fraction is
0.00 until p = 0.001), and eps remains the shipped convention. Whether
Huber's cleaner walk (0 bisections vs 7) generalises is an open, cheap
experiment. Note also (Astra): both arms start at p = 1 from the energy
seed, but Huber kappa = 1 minimises at s* = Q, not Q/2, so its first rung
was a cold solve -- the race was never symmetric at the top.

**P0.3 + P0.4 -- section 20's block was wrong and its one refutation was
an artifact.** The fixed-tf test monitored `Phi([1:6 14], 8:14)` -- the
terminal shooting Jacobian, whose singularity means something only AT t_f.
An interior conjugate point is a Jacobi field vanishing in the FULL state
(mass included): only then does its zero-extension give an admissible
variation with zero second variation. Correct block: `Phi(1:7, 8:14)`
(sol; settled by Astra's admissible-variation argument after the host
defended the old block). Separately, on an initial COAST the state block
is structurally zero and the old instrument counted an exact-zero sample
as a focal point. The refuted record (2,5)@1.4 (p = 0.005) starts on a
coast: its first det sample is exactly 0 (sigma_min/sigma_max = 0) -- that
was the "conjugate point". `ms_conjugate_test` now: monitors the full
block, equilibrates before the sign test, skips samples until full rank is
first attained (`.firstFullRank`), samples THROUGH t_f (the old loop left
the final 1/K of every transfer unmonitored -- 8.3% at K = 12), counts the
last bracket (the old `atFinal` rule subtracted a strictly interior
crossing; it never fired in 14 + 4,405 verdicts), and echoes its spec.
`ms_bvp` returns `.Yend` so the free-time flow column exists at t_f.

Re-sweep (`run_conj_fixedtf_sweep`, 0.6 min; every verdict bound to the
lambda0 it was computed on): **15/15 PASS** -- 7 energy, 7 fuel grid
records, and the gamma-continuation rewalk. The section-21 narrative
"the refutation explains the gamma anomaly" is therefore wrong in
mechanism: the anomaly was a basin (the rewalk found the better optimum,
m_f 0.949005 vs 0.944025, and stays in the catalog on MASS, not on a
verdict), and the conjugate test never disagreed with the walk. The
second-order gate remains production policy; it simply had not yet caught
anything real.

**P0.1 -- the catalog mislabelled its propulsion.** `build_minfuel_catalog`
hardcoded `isp_s = 1710` (the min-time catalogs' value) beside
`c_nd = 8.673746`, which is Isp 900 s exactly -- the 12x12 torus substrate
(0.07 N / 900 s / 150 kg). A consumer recomputing c from `isp_s`, as every
other builder does, got delta-V 1.90x too large. Fixed: Isp and thrust are
now DERIVED from the record's `c_nd`/`Tmax_nd` (900 s, 0.070 N, asserted),
and `catalog_schema('validate')` rejects any catalog whose `c_nd` disagrees
with its `isp_s` (all five shipped min-time catalogs still validate
clean). The rebuilt `costate_catalog_dro_tulip_minfuel.mat`: 7 entries,
isp_s 900, thrustN 0.07, conj_pass 7/7 from the new sweep, verdicts bound
to lam0 and to rows 1:7 by assertion; reflight via the catalog's own
recipe: worst position miss < 0.2 km; delta-V 0.46-0.92 km/s.

**Not changed (P1/P2 of the review, recorded in TODO):** `coastFrac` is a
junction count, not a time fraction; the min-time catalogs' 18,249 verdicts
were produced by the pre-fix instrument (same PASS on the 20 golden cells;
the new t_f sample could only ADD refutations -- a re-sweep is cheap and
pending); the race's per-gap bisection cap; `tGrid` persistence.
