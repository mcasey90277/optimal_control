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

## 23. Huber over the whole grid: a second continuation family, cell-dependent, with a different failure character (2026-09-05)

Section 22 retracted the Huber "knockout" on one cell. This is the
robustness test: the saltation-correct Huber arm on all 7 min-fuel grid
records (same 17-rung schedule, gates and caps as the stored eps arms),
plus two lambda/2-seed controls, each deepest solution re-tested with the
corrected fixed-tf conjugate instrument (`run_huber_race_grid`,
`direct/results/huber_grid.mat`, 15.7 min total).

| cell | gamma | seed | rungs H/E | p_min H / E | m_f H / E | dm_f | fails H/E | bisects H/E | wall H/E [min] | coast H/E | conj |
|---|---|---|---|---|---|---|---|---|---|---|---|
| (2,5) | 1.10 | 1 | 17/19 | 0.001 / 0.0018 | 0.941104 / 0.941107 | -2.4e-6 | **0**/10 | 0/8 | **1.0**/8.8 | 0.17/0.17 | PASS |
| (2,5) | 1.20 | 1 | 17/17 | 0.001 / 0.0017 | 0.947041 / 0.947046 | -4.5e-6 | **0**/9 | 0/7 | **1.2**/7.5 | 0.33/0.33 | PASS |
| (2,5) | 1.40 | 1 | 17/14 | 0.001 / 0.005 | 0.944018 / 0.944025 | -6.2e-6 | **0**/8 | 0/6 | **1.9**/5.7 | 0.42/0.42 | PASS |
| (6,8) | 1.10 | 1 | 18/19 | 0.001 / 0.0012 | 0.928081 / 0.928082 | -8.4e-7 | 1/9 | 1/7 | **1.8**/10.5 | 0.08/0.08 | PASS |
| (6,8) | 1.20 | 1 | 13/19 | **0.033** / 0.0013 | 0.943751 / 0.943936 | -1.8e-4 | 9/9 | 7/7 | 6.6/6.2 | 0.00/0.33 | PASS (at 0.033) |
| (1,2) | 1.10 | 1 | **0**/17 | -- / 0.0034 | -- / 0.901113 | -- | 2/9 | 0/7 | 1.3/11.8 | --/0.25 | -- |
| (1,2) | 1.20 | 1 | **0**/19 | -- / 0.001 | -- / 0.942523 | -- | 2/2 | 0/2 | 1.6/2.5 | --/0.58 | -- |
| (2,5) | 1.20 | **0.5** | 17/17 | 0.001 / 0.0017 | 0.947041 / 0.947046 | -4.5e-6 | 0/9 | 0/7 | **0.8**/7.5 | 0.33/0.33 | PASS |
| (1,2) | 1.20 | **0.5** | 6/19 | **0.77** / 0.001 | 0.933445 / 0.942523 | -9.1e-3 | 11/2 | 9/2 | 8.3/2.5 | 0.00/0.58 | PASS (at 0.77) |

("rungs E" counts the eps arm's bisection inserts; the eps arms are the
09-02 records, not re-run.)

**Findings.**

1. **Where Huber works it is strictly better than eps on this protocol:**
   four of seven records reach p = 0.001 with 0-1 failures in 1.0-1.9 min
   against eps's 8-10 failures and 5.7-10.5 min, and land on the SAME
   solution (|dm_f| < 7e-6, coast fractions identical, conjugate PASS on the
   corrected instrument). The 09-02 "knockout" ran the wrong way.
2. **Where it fails it fails HARD, not softly.** eps never retires -- its
   failures are recoverable bisections. Huber hits Newton residual floors
   (normR stalls at 1e-3..3e-5 over 7 bisections on (6,8)@1.2 at p ~ 0.033;
   normR 1e-2..8e-3 on (1,2)@1.2 at p ~ 0.77) and retires. Those are
   walls with the Jacobian now correct -- a property of the field on those
   cells, not of the solver.
3. **The kappa = 1 seed asymmetry is real and measurable** (Astra's point):
   with the energy costates verbatim, (1,2) cannot solve its FIRST rung
   at either gamma (normR 1.5, Hdrift 0.25); with lambda/2 the p = 1 rung
   converges and the walk proceeds to p ~ 0.77. On the clean cell the
   lambda/2 seed gives the identical answer 33% faster (0.8 vs 1.2 min).
   Use lambda/2 as Huber's standard seed from the energy solution.
4. **The wall correlates with the coast structure the FUEL solution wants,
   not cleanly but suggestively:** the cell that never works, (1,2), is the
   longest transfer (30 d) whose eps solution coasts 25-58% of the arc;
   (6,8)@1.2 walls at 33% coast; the four clean records have eps coast
   0.08-0.42. Huber never coasts exactly (s = kappa Q > 0 wherever Q > 0),
   so on a cell whose extremal is mostly coast Huber must represent long
   near-zero-throttle arcs with a tiny positive s -- and its Q hovers near
   the jump at 1 across the whole arc. A grazing/near-grazing crossing
   (n'F- -> 0 in the saltation denominator) is the natural suspect for the
   residual floors. This is a HYPOTHESIS; the mechanism campaign
   (diagnostics along the walk: number and transversality of Q = 1
   crossings per segment, saltation condition numbers) would settle it.

**What this changes.** Huber is a legitimate second continuation family --
faster and failure-free on the cells it can enter, brittle on the
coast-dominated ones -- and the two families fail in complementary ways
(eps: slow, always recovers; Huber: fast, or walls). The obvious production
move is a **hedged race**: run both from the energy seed (Huber on lambda/2),
take whichever reaches the floor, and keep eps as the fallback that always
arrives. What it does NOT change: eps remains the shipped convention and
the catalog is unchanged. Recorded, not yet built: the mechanism study
(item 4), a Huber-then-eps handoff at the wall (Huber to 0.03, then switch
family on the same junctions -- the two limits are the same problem), and
Huber on the hard cells where eps itself fails (the high-gamma band, the
deep-thrust wall).

## 24. The Huber walls carry grazing-bifurcation signatures -- not folds, not the gate (2026-09-06; wording corrected in section 27)

MfMax review (`doc/mfmax_ideas_review.md`) offered two explanations for the
section-23 walls that its machinery would cure: a FOLD in p (arclength
continuation passes it) or an over-tight rung gate (its `homCI` accepts
continuation points at |S| < 1e-3). Both were tested on both walls
(`run_huber_wall_diag`, `huber_wall_diag.mat`; `ms_bvp` now reports
`cond(J)` at every final iterate, `run_minfuel_race` gained a loose-rung /
tight-floor gate with floor refinement). Both are ruled out.

**Not a fold.** cond(J) rises smoothly along the (6,8)@1.2 walk
(1.7e7 -> 1.3e8 over 13 rungs) and is FLAT across all nine failed iterates
(1.27-1.29e8). On (1,2)@1.2 (lambda/2 seed) it is 6.6-8.3e8 along the
accepted rungs and 8.6-9.1e8 at the failures near the branch; the 1e10-3e11
values occur only at iterates that had already diverged (normR 0.1-0.4).
A turning point would show dS/dz going singular on approach; nothing does.

**Not the gate.** With MfMax's own acceptance (normR < 1e-3, Hdrift < 1e-2)
the (6,8) arm accepted one extra rung (p = 0.03198 at normR 8.5e-4), which
then could not be tightened at the floor, so the arm fell back to its
deepest tight rung -- net one rung SHALLOWER than the tight protocol
(0.03409 vs 0.03275). On (1,2) every failure sat at normR >= 4.8e-3 and the
loose gate never fired. The walls are sharp in p (converges in 31
iterations at 0.03275, stalls at 0.03236) and the residual floor WORSENS
as p decreases (3e-5 at 0.0330, 2e-3 at 0.0300): the field, not the gate.

**What it is: a grazing bifurcation.** `huber_switch_diag` on the last
good rung of each wall, with three clean cells as controls:

| cell | gamma | p | Q = 1 crossings | min |dQ/dt| at a crossing | imminent graze |
|---|---|---|---|---|---|
| **wall (6,8)** | 1.2 | 0.0328 | 8 | **0.045** | -- |
| **wall (1,2)** | 1.2 | 0.774 | 2 | 2.51 | **Q_max = 0.9907 at t/t_f = 0.01** |
| clean (2,5) | 1.2 | 0.001 | 8 | 0.80 | none |
| clean (2,5) | 1.1 | 0.001 | 9 | 0.10 | none |
| clean (6,8) | 1.1 | 0.001 | 7 | 0.19 | none |

(6,8) has a crossing that is 2-18x less transversal than any on the clean
cells: the saltation denominator n'F- = dQ/dt is heading to zero -- two
crossings about to annihilate (a burn arc closing). (1,2) has a local
maximum of Q 0.9% below the jump at the very start of the transfer -- a
burn arc about to be BORN. In both cases, as p decreases a little further,
the switch structure changes at a tangency. At that point the trajectory
map is not differentiable in the unknowns (the saltation matrix
[I + (F+ - F-)n'/(n'F-)] blows up as n'F- -> 0); Newton stalls at a
residual floor with a perfectly finite cond(J), and no reparametrization
of p -- ladder, bisection, arclength -- passes it, because the solution
curve itself has a corner there. The eps ramp never meets this: its control
is continuous, so a new switch is born smoothly. The lambda/2 seed, by the
way, is confirmed exact: the p = 1 rung converges in 2 iterations from it.

**Consequences.**
1. Ideas 1.1 (loose gate) and 1.2 (arclength) of the MfMax review are
   closed for these walls; 1.3 (rescaling) is untested but cannot fix a
   corner; 1.5 (IC homotopy for the high-gamma band) is untouched by this.
2. The cure has to change the FAMILY at the corner, not the walk:
   (a) hand off to eps at the wall (continuous law, passes the structural
   change), optionally returning to Huber after -- the FINDINGS 23 handoff,
   now with a reason; (b) a Huber-eps hybrid: replace the jump by a ramp of
   tiny width delta (a PLQ with one more knee) so grazes become smooth --
   one field function, the same tests; (c) step OVER the bifurcation: when
   a rung fails, try a LARGER step past the graze and re-converge with the
   new structure -- our bisection does the opposite and walks straight
   into the singular point.
3. `huber_switch_diag` (indirect/) is the standing instrument: run it on
   every accepted rung during a walk and a wall becomes predictable
   (min |dQ/dt| falling, or an extremum of Q approaching 1) one or two
   rungs before it happens -- which is exactly when to hand off.

## 25. The Huber-eps hybrid passes both Huber walls -- and inherits eps's floor (2026-09-06)

Section 24 said the cure for Huber's walls had to be a FAMILY change at the
corner. `huberc` is that family: Huber's gentle core s = pQ below Q = 1,
then a CONTINUOUS ramp p -> 1 over Q in [1, 1+delta] instead of the jump,
s = 1 beyond; L(s) is the integral of the inverse law matched at s = p, so
H(s) is convex and the law is its exact argmin (`test_minfuel_pmp`: dH
5.6e-17; continuous at both knees; AD = FD; STM = FD through the generic
propagator with no saltation). Default delta = p. Same 17-rung schedule,
lambda/2 seed, tight protocol, on the two walls and the clean control:

| cell | gamma | family | rungs | p deepest | m_f | fails / bisects | wall |
|---|---|---|---|---|---|---|---|
| (6,8) | 1.2 | huber (saltation) | 13 | **0.0328 (wall)** | 0.943751 | 9 / 7 | 6.6 min |
| (6,8) | 1.2 | **huberc** | 18 | **0.0026** | 0.943921 | 9 / 7 | 6.1 min |
| (6,8) | 1.2 | eps (09-02) | 19 | 0.0013 | 0.943936 | 9 / 7 | 6.2 min |
| (1,2) | 1.2 | huber, lambda/2 | 6 | **0.774 (wall)** | 0.933445 | 11 / 9 | 8.3 min |
| (1,2) | 1.2 | **huberc** | 17 | **0.00138** | 0.942505 | 22 / 20 | 18.2 min |
| (1,2) | 1.2 | eps (09-02) | 19 | 0.0010 | 0.942523 | 2 / 2 | 2.5 min |
| (2,5) | 1.2 | huber (saltation) | 17 | 0.0010 | 0.947041 | **0 / 0** | 1.2 min |
| (2,5) | 1.2 | **huberc** | 15 | 0.0039 | 0.947028 | 8 / 6 | 5.8 min |
| (2,5) | 1.2 | eps (09-02) | 17 | 0.0017 | 0.947046 | 9 / 7 | 7.5 min |

**Both walls are gone.** On (6,8) huberc converges at p = 0.03 in 4
iterations where Huber stalled nine times, and continues 12x deeper; on
(1,2) the p ~ 0.77 wall is not there at all and the walk reaches the floor.
Final masses agree with eps to 1.5e-5 and 1.8e-5. This is the direct
confirmation of section 24's mechanism: remove the jump and the grazing
bifurcations stop being corners.

**What it costs.** With delta = p the ramp slope (1-p)/delta is ~1/p --
at the floor, 382 on (6,8) vs eps's 391 -- so huberc meets exactly the
steep-ramp difficulty eps has at small p, and it meets it slightly earlier
(its ramp is half eps's width at equal p). Visible in three places: (6,8)'s
new stall at p ~ 0.0025 is a Q MINIMUM at 1.0149 (a coast being born)
with a perfectly transversal crossing structure (min |dQ/dt| 1.19) -- a
steepness stall, not a graze; (1,2) needed 20 bisections in the
p = 0.5-0.6 and p < 0.002 bands; and on the CLEAN cell huberc stalls at
p = 0.0039 with 8 failures where plain Huber -- jump, saltation, no ramp
at all -- walked to 0.001 with none. Huber's cleanliness where it works
was BECAUSE it has no steep ramp; huberc gives that up to pass the walls.

**Where this leaves the three families.**

| | eps | huber (saltation) | huberc (delta = p) |
|---|---|---|---|
| walls (grazing corners) | none (continuous) | 3 of 7 cells | none measured |
| floor behaviour | steep ramp, 8-10 bisections | clean, 0-1 fails | steep ramp, 8-22 bisections |
| speed on clean cells | 6-11 min | **1-2 min** | 6 min |
| reaches the floor | 7/7 | 4/7 | 3/3 tested |
| m_f agreement | -- | < 7e-6 | < 2e-5 |

**Next knob, not yet turned.** delta need not shrink with p. A FIXED
delta (say 0.02-0.05) during the p-walk keeps the ramp gentle while the
core vanishes (at p -> 0 huberc becomes a one-sided eps of width delta);
a short second walk then takes delta down. Or delta = 2p to match eps's
steepness exactly. Or the hedge the numbers now argue for directly: run
HUBER, and switch to huberc the rung before `huber_switch_diag` predicts a
graze -- Huber's speed on the 60-70% of cells without grazes, huberc's
continuity where they occur, eps as the fallback that always arrives.
`run_minfuel_race` already carries `seedScale` and the loose/tight gates;
a `family schedule` (which family at which rung, with the diag as the
trigger) is the missing piece.

Records: `direct/results/minfuel_huberc_c{68,12,25}_g120.mat`.

## 26. The delta schedule: decouple the structural walk from the sharpening (2026-09-06; "family-independent floor" corrected in section 27)

Section 25 left huberc with eps's steep ramp at the floor because delta
shrank with p. `run_minfuel_race` now walks (p, delta) rungs: a delta RULE
for the p-walk (fixed, or a function of p) and an optional STAGE 2 that
holds p at its deepest value and walks delta down. Three rules on the two
former walls and the clean control (lambda/2 seed, tight protocol):

| cell | delta rule | p-walk fails | p reached | delta-walk fails | final delta | m_f | eps m_f | wall |
|---|---|---|---|---|---|---|---|---|
| (6,8) | delta = p (s.25) | 9 | 0.0026 (wall) | -- | -- | 0.943921 | 0.943936 | 6.1 min |
| (6,8) | delta = 2p | 9 | 0.0012 | -- | -- | 0.943929 | | 5.5 min |
| (6,8) | **delta = 0.03 fixed** | **0** | **0.001** | 8 | 0.0026 | 0.943930 | | 5.3 min |
| (1,2) | delta = p (s.25) | 22 | 0.0014 | -- | -- | 0.942505 | 0.942523 | 18.2 min |
| (1,2) | delta = 2p | 10 | 0.001 | -- | -- | 0.942510 | | 12.3 min |
| (1,2) | **delta = 0.03 fixed** | 15 | **0.001** | 12 | 0.0014 | 0.942510 | | 19.8 min |
| (2,5) | delta = p (s.25) | 8 | 0.0039 (stall) | -- | -- | 0.947028 | 0.947046 | 5.8 min |
| (2,5) | delta = 2p | 10 | 0.0019 | -- | -- | 0.947037 | | 7.2 min |
| (2,5) | **delta = 0.03 fixed** | **0** | **0.001** | 8 | 0.0038 | 0.947041 | | 8.0 min |

**1. With delta fixed, the p-walk is trivial.** On (6,8) and (2,5) it went
17/17 to p = 0.001 with ZERO failures, three Newton iterations per rung
at the end, and landed on eps's coast structure (0.33 on (6,8), exactly
eps's) -- the two cells where delta = p had walled or stalled. (1,2), the
30-day cell, still spent 15 failures in the p = 0.5-0.7 band, which every
rule and both other families also found hard (a structural change region
for that cell, not a delta effect); at least one of those was the
spurious tolR floor (normR 1.3e-10 vs tolR 1e-10, see 3).

**2. All the difficulty is in the sharpening, and it is the same for every
family.** Stage 2 costs 8-12 failures on every cell and stops at
delta ~ 0.0014-0.0038 -- the same place delta = p (0.0026), delta = 2p
(0.0022-0.0039) and eps's own ramp (p 0.001-0.0017) all stop. That is the
bang-bang-limit floor, not a property of Huber, huberc or eps: as the
ramp width goes below ~0.003 the switches become too sharp for Newton on
a 12-segment multiple-shooting mesh, whichever family produced them. And
by then m_f has CONVERGED: on (6,8) it moves 0.943928 -> 0.943930 from
delta 0.03 to 0.0026; on (2,5) 0.947041 at delta 0.0038 equals plain
Huber's floor value to 6 digits. Sharpening past delta ~ 0.003 changes
the answer at the 1e-6 level and costs most of the wall time.

**3. Two side observations.** (a) The fixed-delta floor solutions PASS the
single-shooting acceptance gate on (6,8) and (2,5) (`accept ok=1`,
|dz| = 0) -- most deep-eps entries do not (FINDINGS 19); a gentler ramp
at the same m_f is a friendlier solution for any single-shooting
consumer. (b) `ms_bvp`'s tolR = 1e-10 sits ON the residual floor of the
30-day arcs: three "failures" on (1,2) were normR 1.0e-10 .. 2.0e-10 and
one of them abandoned a gap. tolR 3e-10 (or floor-aware) is a one-line
fix with measurable effect on (1,2).

**What this changes.** The right decomposition of the energy -> fuel
homotopy is not "which family" but "structure, then sharpness":

  1. p-walk at FIXED delta ~ 0.03 (huberc): the switch structure forms with
     no drama -- 0 failures where the jump family walled;
  2. sharpen delta to ~ 0.003-0.005: m_f converged to 1e-6, ss-acceptable;
  3. STOP. Beyond that the problem is bang-bang with a fixed switch
     structure, which is a SWITCHING-TIME problem, not a homotopy -- the
     Maurer-Osmolovskii reduced problem the certification survey already
     wants for strictness (doc/extremal_and_local_min_survey.md, item 4).
     The homotopy's last, most expensive rungs are doing badly what a
     switching-time Newton would do exactly.

Plain Huber keeps one honour: on cells with no grazing bifurcation it
reaches p = 0.001 in 1-2 minutes with no ramp at all (s.23) -- the fixed-
delta recipe costs 5-8 minutes there. The family schedule of s.25 (Huber
first, huberc when `huber_switch_diag` predicts a graze) is still the
fastest route on the grid; the fixed-delta huberc walk is the ROBUST one.

Records: `direct/results/minfuel_huberc_c{68,12,25}_g120_{d2p,dfix}.mat`.

## 27. Corrections from Astra review #2 (2026-09-06): what sections 24-26 may and may not claim

A second GPT-6 Astra pass (code + theory; `reviews/huber_review_adjudicated_2026-09-06.md`)
and the resulting P0/P1 fixes change the wording, one number, and nothing
in the measured tables.

**Section 24 -- "grazing bifurcation" is a RISK SIGNATURE, not a proven
event.** (a) `Qdot` is throttle-independent and continuous across the
switch: `Qdot = -Tmax lam_v'lam_r/(m rho)` exactly (n'(F+-F-) = 9e-16
measured), so `huber_switch_diag` now uses it instead of a secant on the
ode113 grid. With the exact value the (6,8) wall's least-transversal
crossing is |Qdot| = 0.080 (not 0.045) against 0.108-0.80 on the clean
cells -- still the smallest, a 1.4-10x contrast, not 2-18x. The (1,2)
Q-maximum of 0.9907 at t/tf = 0.01 stands. (b) A grazing bifurcation
requires h = 0, hdot = 0, hddot != 0 on the branch plus a nonzero
unfolding in p; we observed the two indicators, not the event. (c) "No
parametrization of p passes a corner" is WRONG: the local behaviour at a
graze is typically square-root (switch-time sensitivity ~ mu^-1/2), which
w = sqrt(mu), explicit switching times, or a structure-change step can
traverse. The defensible statement: the smooth Newton-on-p walk stalls
there, and a continuous family passes it -- both measured. (d) Flat
cond(J) is weaker evidence against a fold than stated (finite-distance
samples at nonzero residual). The family change was A cure, not the only
one; it is the one that worked.

**Section 26 -- the sharpening floor is CONFIGURATION-dependent.** Three
cells on one solver setup (tolerances, scaling, K = 12, step policy)
establish an empirical floor for that setup, not a family-independent
one; ramp width in Q is not a temporal resolution (dt_ramp ~
delta/|Qdot|); segment count alone does not set switch resolution; and
delta -> 0 at fixed p = 0.001 recovers the Huber JUMP (s = pQ below Q = 1),
not exact bang-bang fuel. m_f convergence to 1e-6 is one scalar; switch
times may still differ. The switching-time endgame remains sensible; its
convergence certifies stationarity, not minimality, and the direction
control alpha(t) on burns must be in the accessory problem.

**The conjugate instrument** now reports `verdict` in {PASS, FAIL,
ENDPOINT, UNDETERMINED}: nothing-testable is no longer a PASS, a root on
the bracket ending at t_f is ENDPOINT (a weak, non-strict minimum -- to
examine, not a refutation), zero samples merge with adjacent sign changes,
and the determinant sign comes from a pivoted LU of the equilibrated block
with the magnitude kept as a log-determinant (an underflowing raw det no
longer reads as a root). Re-sweep on the fixed instrument: 15/15 PASS, all
tested, none endpoint-only; golden 20/20; the catalog rebuilt clean and its
verdicts are bound to converged re-solves at the entry's own p.

**Not changed by the review:** every measured table in sections 22-26
(the walls, huberc passing them, the delta schedule) and the shipped
catalog's entries.

## 28. The high-gamma band: the gamma-walk opens it, three families race in it, huberc alone reaches the floor at gamma = 2 (2026-09-07)

**The seed route.** The direct min-energy solve fails for gamma >= 1.7 on
(2,5) and gamma >= 1.4 on (6,8) and (1,2) (FINDINGS 19), so eps never had a
seed there. `run_gamma_walk` -- the MfMax homCI idea (a continuation in the
BOUNDARY DATA, here tf = gamma*tfMin) in the form FINDINGS 21 already
trusted -- walks the SMOOTH min-energy problem from each cell's deepest
PASS record, each step seeded by the previous junction states at the same
fractional times, geometric bisection on gamma, conj-tested at every step:

| cell | from | reached | walled at | steps | notes |
|---|---|---|---|---|---|
| (2,5) | 1.40 | **2.00** (1.47, 1.55, 1.70, 1.77, 1.81, 1.85, 2.00) | -- | 9, 2 bisections | cond(J) IMPROVES toward gamma 2 (3.4e6); 3 min total |
| (6,8) | 1.20 | 1.40 (1.30, 1.35, 1.40) | **1.43** | 7 | cond(J) 6e8 -> 5e9 -> 6.5e10 on approach |
| (1,2) | 1.20 | 1.247 (1.223, 1.247) | **1.27** | 6 | cond(J) 4e9 -> 1.2e10 on approach |

Twelve records, all conj PASS. Two walls are on the p = 1 smooth problem --
no switch structure exists there, so no throttle family can help -- and on
both, cond(J) grows by 1-2 orders of magnitude on approach, the signature
of a FOLD in gamma that the Huber walls did NOT show (FINDINGS 24). This
is where MfMax idea 1.2 (arclength) may finally apply; a target-orbit IC
homotopy and a larger K are the other candidates. Not attempted.

**The race** (`run_highgamma_race`: huberc fixed delta 0.03 + delta-walk
to 0.003, eps, huber with lambda/2; every deepest solution re-tested with
the corrected fixed-tf conjugate instrument; race tolR 3e-10 after the
1e-10 floor cost one arm -- see below):

| cell | gamma | huberc | eps | huber | m_f (best) | m_f energy | gain [% m0] |
|---|---|---|---|---|---|---|---|
| (2,5) | 1.47 | floor | 0.0089 | **floor, 0 fails, 1.1 min** | 0.947045 | 0.939708 | +0.73 |
| (2,5) | 1.55 | floor | 0.0017 | **floor, 0, 1.4 min** | 0.948688 | 0.941156 | +0.75 |
| (2,5) | 1.70 | floor | 0.0053 | **floor, 0, 1.7 min** | 0.949246 | 0.942154 | +0.71 |
| (2,5) | 1.77 | floor | 0.0045 | wall 0.051 | 0.949319 | 0.941559 | +0.78 |
| (2,5) | 1.81 | floor | 0.0040 | **floor, 0, 2.0 min** | 0.949306 | 0.940983 | +0.83 |
| (2,5) | 1.85 | floor | 0.0037 | **floor, 0, 2.1 min** | 0.949251 | 0.940300 | +0.90 |
| (2,5) | **2.00** | **floor (alone)** | wall 0.66 | wall 0.32 | 0.948585 | 0.937705 | +1.09 |
| (6,8) | 1.30 | floor | 0.0020 | **floor, 0, 1.1 min** | 0.950404 | 0.942911 | +0.75 |
| (6,8) | 1.35 | floor | 0.0080 | **floor, 0, 1.3 min** | 0.951979 | 0.944307 | +0.77 |
| (6,8) | 1.40 | floor | 0.0050 | **floor, 0, 1.4 min** | 0.952123 | 0.944933 | +0.72 |
| (1,2) | 1.22 | wall 0.81 | **0.0013** | wall 0.82 | 0.941482 | 0.931678 | +0.98 |
| (1,2) | 1.25 | wall 0.85 | 0.067 (re-run) | wall 0.85 | (0.939797) | 0.930044 | -- |

"floor" = p = 0.001 (huberc: then delta to 0.003). Every reached solution
is conj PASS on the corrected instrument. Where two or three families
arrive they agree to ~1e-5.

**Findings.**

1. **Ten new min-fuel records on cells the direct solver never touched**
   ((2,5) x 7 at gamma 1.47-2.00, (6,8) x 3 at 1.30-1.40), all three-family
   agreed where the families arrive, all conjugate-tested. These are the
   first min-fuel solutions above gamma = 1.4 in the program.
2. **At gamma = 2.0, huberc is the only family that reaches the floor.**
   eps stalls at p = 0.66 (m_f 0.9401 -- not close), huber at 0.32. The
   "cells eps cannot enter" test has an unambiguous answer on one cell.
3. **Huber's 0-fail / 1-2-minute pattern holds in the band** (8 of 10
   records where anything arrives), exactly as on the grid; it walled on
   (2,5)@1.77 and @2.00. The family schedule of FINDINGS 25 (Huber first,
   huberc at a predicted graze, eps fallback) is the reading of this table.
4. **The fuel gain saturates near gamma ~ 1.8 on (2,5)** (m_f 0.94931 at
   1.77, 0.94925 at 1.85, 0.94859 at 2.00 despite +1.09% over ITS energy
   value): more time stops buying propellant, and the gamma-basin
   structure of the energy problem survives into fuel. Physics, not
   numerics -- three families agree.
5. **(1,2) near its smooth wall is hard for every family.** At 1.22 only
   eps arrives (0.0013, 16 fails); at 1.247 nobody does (eps 0.067 with
   11 fails, the Huber pair 0.85). This is the 30-day cell one bisection
   from a fold-like wall; nothing to conclude about families from it.
6. **tolR = 1e-10 cost an entire arm**: the (1,2)@1.247 eps p = 1 rung
   "failed" at normR 1.0e-10 and the arm retired in 0.9 min. The race now
   runs at tolR 3e-10; the re-run is the row above. The gamma-walk records
   and all other arms were unaffected (their p = 1 rungs converged below
   1e-10).

**What this changes.** The high-gamma band is no longer "where the direct
solve fails"; it is where the boundary-data homotopy seeds and huberc
finishes. The ten records are catalog material (schema v3, one sheet,
axis3 = gamma); the builder currently reads only `minfuel_grid.mat` and
needs a second source. The two smooth walls at gamma 1.43 / 1.27 are the
next frontier, and for the first time the evidence points at a fold.

Records: `direct/results/minenergy_highgamma.mat` (12 records),
`highgamma_race.mat` (36 arms), `minfuel_hg_*.mat`.

## 29. Catalog v3.1: 18 min-fuel entries to gamma = 2, mixed continuation families, every verdict re-bound (2026-09-07)

The eleven high-gamma solutions of section 28 are packaged with the seven
grid entries: `costate_catalog_dro_tulip_minfuel.mat`, schema 3 minor 1.

**What changed in the schema (v3.1).** Entries from different continuation
families in one catalog: `sheets.family_code` (int8; `smoothing.codes` =
eps 1, huber 2, huberc 3) and `sheets.delta_floor` (finite for huberc
entries) beside `p_floor`; `smoothing.family = 'mixed'`. The validator
rejects a mixed catalog without `family_code`, an unknown code on a solved
entry, and a huberc entry without a positive `delta_floor`
(`test_catalog_schema_v3`, 17/17; the five shipped min-time catalogs still
validate clean). The reconstruction recipe is unchanged in form -- fly
`cr3bp_minfuel_prop` from `Yj(:,1,n)` over `tf_nd` -- with the entry's own
(family, p_floor, delta_floor); huber entries use the event-split
saltation propagator automatically.

**Which arm is packaged (`highgamma_select`).** The SHARPEST solution
among the arms that reached the floor (p <= 0.0015): ramp width of the
throttle law -- huber 0 (jump), huberc delta_floor, eps 2p; ties to eps.
Result: 8 huber (0 width), 2 huberc ((2,5)@1.773 delta 0.0093, @2.000
delta 0.0039 -- the two records where huber walled), 1 eps
((1,2)@1.223, the only arrival there). The other arriving arm of each
record is recorded (`altArms`) but not packaged; where two arrived they
agree to ~1e-5 in m_f.

**Every verdict re-bound.** `run_conj_fixedtf_sweep` now verdicts the
selected high-gamma solution in its OWN family (huber with saltation STMs,
huberc with the generic path): 26/26 PASS on the corrected instrument
(7 energy, 7 grid, 1 rewalk, 11 high-gamma), all tested, none endpoint-
only; the builder binds each verdict to the packaged entry by source,
lambda0, p, family and rows 1:7, and refuses unconverged or untested
verdicts.

**The catalog, entry by entry** (gamma axis now 13 values, 1.1 -> 2.0;
sheets sparse along gamma by design):

| cell | gamma | family | p_floor | delta | m_f | dV [km/s] |
|---|---|---|---|---|---|---|
| (1,2) | 1.10 | eps | 0.0034 | -- | 0.901113 | 0.919 |
| (1,2) | 1.20 | eps | 0.0010 | -- | 0.942523 | 0.522 |
| (1,2) | 1.223 | eps | 0.0013 | -- | 0.941482 | 0.532 |
| (2,5) | 1.10 | eps | 0.0018 | -- | 0.941107 | 0.536 |
| (2,5) | 1.20 | eps | 0.0017 | -- | 0.947046 | 0.480 |
| (2,5) | 1.40 | eps (rewalk) | 0.0017 | -- | 0.949005 | 0.462 |
| (2,5) | 1.473 | huber | 0.0010 | -- | 0.947039 | 0.480 |
| (2,5) | 1.55 | huber | 0.0010 | -- | 0.948679 | 0.465 |
| (2,5) | 1.70 | huber | 0.0010 | -- | 0.949240 | 0.460 |
| (2,5) | 1.773 | huberc | 0.0010 | 0.0093 | 0.949310 | 0.459 |
| (2,5) | 1.811 | huber | 0.0010 | -- | 0.949297 | 0.459 |
| (2,5) | 1.85 | huber | 0.0010 | -- | 0.949242 | 0.460 |
| (2,5) | 2.00 | huberc | 0.0010 | 0.0039 | 0.948585 | 0.466 |
| (6,8) | 1.10 | eps | 0.0012 | -- | 0.928082 | 0.659 |
| (6,8) | 1.20 | eps | 0.0013 | -- | 0.943936 | 0.509 |
| (6,8) | 1.296 | huber | 0.0010 | -- | 0.950396 | 0.449 |
| (6,8) | 1.347 | huber | 0.0010 | -- | 0.951971 | 0.434 |
| (6,8) | 1.40 | huber | 0.0010 | -- | 0.952114 | 0.433 |

(dV from `deltav_from_mf`, Isp 900 s.) Reflight of all 18 entries through
the catalog's own recipe (per-entry family, p_floor, delta_floor): worst
position miss 0.220 km ((6,8)@1.347, huber), median ~0.012 km.

**Reading the table.** On (2,5) the min-fuel m_f rises from 0.9411 (gamma
1.1) to a plateau of 0.9493 at gamma 1.7-1.85 and falls to 0.9486 at 2.0:
the propellant saved by more time saturates at about +0.8% of m0 over the
energy solutions and turns over -- a basin structure inherited from the
energy problem (FINDINGS 19, 28). (6,8) is still rising at its wall
(0.9521 at 1.40); (1,2) is barely started (0.9415 at 1.223). Where the
walls fall (FINDINGS 28: gamma 1.43 and 1.27, on the SMOOTH problem) is
now the limiting factor for this catalog, not the fuel homotopy.

## 30. The sufficiency-hypothesis gates at catalog scale: 18,360/18,360 pass -- after a rank-rule correction (2026-09-07)

`gates_catalog_pass` ran `mintime_hypothesis_gates` over every accepted entry
of the five min-time catalogs (DRO/HALO/DPO -> tulip, L1<->L2 halo):
one dense `tfMinProp` flight + one 7x7 adjoint integration per entry,
~0.75 s each, 4 h wall for 18,360 entries, sidecars `*_gatesprog.mat`.

**First census (fixed rank tolerance 1e-8): H2 0 fails, H3 0 fails, 512
"abnormal" (dim S = 0).** Every one of the 512 had `dimS = 0` -- not 2 --
with spectral gaps sv7/sv6 between 1.4e-7 and 2.5e-6 and lift residuals
|C lam0|/|lam0| between 1e-7 and 2.2e-6. dim S = 0 is impossible for a
solution: the accepted normal lift is a member of S by construction, and it
IS one numerically (nullResid ~ 1e-7). What the flag measured was the
rank tolerance: 1e-8 * sv1 is finer than the accuracy to which the known
member satisfies the constraints, so the one small singular value was
counted as "not small". A null space cannot be resolved finer than its
known member's residual.

**The rule (`lift_space_dim`, TDD 5/5):** dim S = #{sv < tol}, tol =
min(max(rankTol * sv1, 10 * nullResid), 1e-3 * sv1). The 1e-3 cap keeps a
noisy lift from ever inflating the count: a poor lift reads dim S = 0
(unresolved), never > 1. Golden cells unchanged (15/15; their residuals
are 4e-8..1.3e-7, below the old threshold anyway).

**Re-gated the 512 under the rule: 0 abnormal.** Final census, written
back into all five catalogs (`gate_min_lamv`, `gate_min_qmt`, `gate_dimS`,
`hyp_gates` provenance; backups `*.mat.bak_gates`):

| catalog | entries | min\|lam_v\| min / median | min Q_mt | max nullResid | max sv7/sv6 | dim S |
|---|---|---|---|---|---|---|
| DRO -> tulip | 4,439 | 9.3e-6 / 4.4e-3 | 1.5e-3 | 2.2e-6 | 2.5e-6 | 1, all |
| HALO -> tulip | 4,596 | 8.1e-6 / 2.8e-3 | 1.1e-3 | 1.4e-6 | 1.4e-6 | 1, all |
| DPO -> tulip | 4,457 | 2.4e-5 / 5.0e-3 | 9.7e-4 | 8.2e-7 | 1.1e-6 | 1, all |
| L1 -> L2 halo | 2,338 | 3.3e-6 / 2.4e-3 | 8.7e-4 | 7.9e-7 | 7.7e-7 | 1, all |
| L2 -> L1 halo | 2,530 | 2.8e-6 / 2.7e-3 | 1.8e-3 | 8.5e-7 | 8.4e-7 | 1, all |

**What this licenses (audit section 6).** For every entry with
`conj_pass = 1`: a normal extremal (H = 0 with lambda_0 = 1, checked to
~1e-10), no abnormal lift (dim S = 1 with a spectral gap of at least five
orders), strong Legendre and the all-burn control hold along the flown arc,
and the BCT determinant has no sign change at the sampled times. The
sampling caveat is the only one left on the min-time side.

**A graded quantity worth carrying: min|lam_v| falls with thrust.** The
primer never vanishes (H2 holds everywhere, floor 2.8e-6 against a
1e-6 gate), but it comes closest on the HIGH-thrust rungs: 20 entries
below 1e-5 and 593 below 1e-4, of which 248 sit at 15 N, 167 at 12 N, 122
at 10 N, and none below 2 N (medians 3e-4 at 15 N vs 1e-1 at 0.5 N; the
scale is pinned by H = 0, so these are comparable). A near-null of the
primer is where its direction turns fastest -- the fast reorientations of
short high-thrust transfers -- and it is where the Legendre form
(T/m)|lam_v| is weakest. Those are the entries a strict certificate would
find hardest, and the first to examine if the dense-det sampling upgrade
ever reports a near-root.

Reproduce: `gates_catalog_pass(catMat)` (resumes from the sidecar);
`test_lift_space_dim`, `test_mintime_gates`.

## 31. The REGIME MAP: 63 arms under one identical budget -- no family dominates, and every failure is a switch-structure failure (2026-09-07)

Goal (Mike, 2026-09-07): understand WHEN the eps, huber and huberc walks
work, and find transfers only one of them solves. Sections 23-28 could not
answer it: those arms came from three differently-tuned campaigns, so a
"wall" could always have been a budget or a seed. `run_regime_map` fixes
every knob -- same p schedule (17 rungs 1 -> 0.001), tolR 3e-10, wallSec
300, maxBisect 3, maxGaps 2, Hdrift gate -- and varies only what must vary
by construction: each family's correct warm seed (lambda/2 for huber and
huberc, whose p = 1 minimiser is s* = Q, not Q/2) and huberc's fixed-delta
p-walk plus delta sharpening (section 26). **63 arms = 21 (cell, gamma)
cases x 3 families**, cells (1,2) (6,8) (2,5), gamma 1.1 -> 2.0.

**Scoring: "reached p = 0.001" is NOT "found the optimum".** Measured here:
eps stops at p = 0.0016 on (2,5)@1.1 yet its m_f agrees with the two arms
that reached 0.001 to 2.4e-6 -- its remaining rungs were unnecessary, and
scoring that as a wall would invent a family difference. A family SOLVES a
case when its m_f agrees with the case's best arm to 1e-3 (0.15 kg of 150).
`regime_verdicts` also flags the trap in a relative criterion: the best arm
solves by construction, so a case whose winner never reached the limit
reads "furthest", never "only".

### The map

| | eps | huber | huberc |
|---|---|---|---|
| fails | **1 / 21** | 5 / 21 | 3 / 21 |
| uniquely solves | 3 (2 at the limit) | **0** | 1 |
| median wall where it solves | 8.0 min | **1.6 min** | 8.5 min |
| median failed rungs | 9 | **0** | 8 |

**Huber is the fast family, never the capable one: 0 of 21 cases are
huber-only, and where it works it walks with ZERO failed rungs in a fifth
of the time.** eps is the most robust; huberc sits between and is the only
family that reaches the bang-bang limit at gamma = 2 on (2,5).

**Genuine one-family-only cases (winner AT the limit):**

| case | verdict | losers, m_f deficit | switches: solved vs failed |
|---|---|---|---|
| (1,2) @ 1.223 | **eps ONLY** | huber 9.5e-3, huberc 9.3e-3 | 9 vs 2, 4 |
| (2,5) @ 2.00 | **huberc ONLY** | eps 8.5e-3, huber 6.3e-3 | 9 vs 5, 5 |

Two more cases have a single best arm that never reached the limit --
(1,2)@1.247 (eps furthest at p = 0.067) and (6,8)@2.0 (eps furthest at
p = 0.23): there NO family solved the problem. Plus one two-family case,
(1,2)@1.2 = eps + huberc, huber 9.1e-3 light. The remaining 16 cases are
solved by all three, agreeing to 8e-7 .. 2e-4.

**Different cells at the SAME gamma favour DIFFERENT families**: at
gamma = 2.0, (2,5) is huberc-only and (6,8) is eps-furthest with huberc
failing. The winner is not predictable from gamma, or from tf, or from
switch count -- which is the case for a family SCHEDULE rather than a
family choice.

### The mechanism: failures are switch-structure failures

Every failed arm stopped at a solution with FEWER switches than the case's
winner (2, 4 or 5 against 7-11), and **7 of the 8 failed arms stop at a
grazing-risk extremum of Q -- a local extremum sitting 0.3% to 2.6% short
of the threshold, i.e. a switch about to be born**:

| case | family | status | p_stop | switches | min\|dQ/dt\| | graze gap |
|---|---|---|---|---|---|---|
| (1,2)@1.2 | huber | failed | 0.774 | 2 | 2.48 | **0.0093** |
| (1,2)@1.223 | huber | failed | 0.818 | 2 | 2.21 | **0.0067** |
| (1,2)@1.223 | huberc | failed | 0.809 | 4 | 0.49 | none |
| (1,2)@1.247 | huber | failed | 0.846 | 2 | 1.88 | **0.0080** |
| (1,2)@1.247 | huberc | failed | 0.846 | 2 | 1.88 | **0.0081** |
| (2,5)@2 | eps | failed | 0.657 | 5 | 1.37 | **0.0260** |
| (2,5)@2 | huber | failed | 0.315 | 5 | 1.60 | **0.0073** |
| (6,8)@2 | huber, huberc | failed | 0.664 | 5 | 3.13 | **0.0030** |

The arms that SOLVED those same cases end with no grazing-risk extremum at
all, except eps on (1,2)@1.223 whose gap is 0.082 -- an order of magnitude
wider. So the wall is not a fold, not a tolerance and not the seed: **the
walk stops where the switch structure has to change and the family's law
cannot carry it through.** This generalises section 24, which had the
signature only for huber; eps and huberc wall the same way.

### Hypotheses, adjudicated

* **H1 (huber walls at grazing-risk structure changes): SUPPORTED, and it
  is not specific to huber.** 7 of 8 failures carry the signature. The one
  exception, huberc on (1,2)@1.223 (4 switches, no near-tangency,
  min|dQ/dt| = 0.49), is unexplained and is the open item.
* **H2 (eps fails at the floor on MANY-switch cells): REFUTED.** eps solved
  cases with 6 to 11 switches (median 8.5) and its single failure is a
  5-switch case; it solved at tf = 32.6 d and failed at 31.8 d. Neither
  switch count nor flight time explains it.
* **H3 (huberc's early walls unexplained): PARTLY EXPLAINED.** 2 of its 3
  failures carry the graze signature; (1,2)@1.223 does not.

### A correction to section 23

The 09-02 grid recorded huber reaching the floor on **4 of 7** records.
Under the identical budget with the lambda/2 warm seed (the Astra correction
of 09-05) and the mass criterion, huber solves **6 of 7** -- it fails only
(1,2)@1.2. The 4/7 figure measured the seed and the budget, not the family.

### Scope and what is NOT claimed

Three departure/arrival cells of one orbit pair (DRO -> tulip), one thrust
level, one K. "Fails" means under THIS budget: a bigger bisection budget or
a larger K may carry an arm further, and the two "furthest" cases show the
sample already contains problems no family solves. The huber-only hunt is
open -- 0 of 21 cases here, and it is the point of widening the cell sample
(regime map chunk 4).

Reproduce: `run_regime_map` (63 jobs, resumable, partitionable),
`regime_table`, `regime_verdicts`; tests `test_regime_features` (17),
`test_regime_verdicts` (12).

## 32. The 70 mN abstract case is SOLVED and certified (2026-09-08) -- but the "abstract is optimistic by 50%" verdict was WRONG: see the correction in section 35

> **CORRECTION (same day, section 35).** This section, and sections 33-34
> after it, compared the abstract against cell (1,11) of the 12x12 fine
> sheet -- arrival phase sA = 0.9087. **The abstract's case is a DIFFERENT
> arrival phase, sA = 0.0754**, and at that phase 70 mN reaches the tulip in
> **17.80 d for 0.7485 km/s and 12.20 kg (8.1% of wet mass)** -- the
> abstract's figures, exactly, and certified. Everything below about the
> abstract being optimistic is WITHDRAWN. The transfer physics in this
> section stands; it is a statement about sA = 0.9087 only.



The cislunar abstract ("MinTime Tulip Transfer", poster content approved for
public release) states a 150 kg spacecraft with a **70 mN** thruster at
**Isp 900 s** reaching a 7-petal tulip from a DRO in **~18 days for ~0.75
km/s and 12.2 kg**. That operating point was NOT in the catalog, which
floors at 0.5 N and runs Isp 1710 s: 70 mN is seven times below the lowest
catalogued rung. The abstract's figures follow from the all-burn identity
dV = (T/m0)*tf, not from a converged transfer.

### The certified solution

| quantity | value |
|---|---|
| thrust / Isp / m0 | 70 mN / 900 s / 150 kg |
| **t_f** | **5.96394 ND = 26.436 days** |
| **dV** | **1.1360 km/s** |
| **propellant** | **18.12 kg (12.1% of wet mass)** |
| ms residual | 2.45e-12 |
| flown arrival | 0.069 km |
| **pumpkyn tfMin** | **accepts, |dz| = 2.1e-10** |
| conjugate test | PASS |
| gates | min|lam_v| 0.599, min Q_mt 1.174, dim S = 1 |

Mesh independent: K = 24, 48 and 96 all converge to t_f = 5.96394 ND.
Cell (1,11) of the tau=1 / Np=7 fine sheet -- the abstract's own geometry.
Stored: `indirect/results/mintime_70mN_{certified,direct}.mat`.

**Against the abstract: t_f +47%, dV +51%, propellant +48%.** The abstract
understates the transfer in the OPTIMISTIC direction. The poster needs
26.4 days, 1.14 km/s and 18.1 kg, or a different operating point.

### How it was reached, and what failed first

Indirect thrust continuation from the banked 0.09 N solution walks to
**75.5 mN (15.15 d)** and then walls at 75.0 mN. The wall is real: it moved
82.5 -> 80 -> 75.5 mN as the machinery improved (t_f guess sweep, K 24->48)
and then K 48->96 did not move it at all. Approaching it,

| T (mN) | 90 | 85 | 80 | 78 | 76 | 75.5 |
|---|---|---|---|---|---|---|
| t_f (d) | 13.21 | 13.67 | 14.28 | 14.60 | 15.02 | 15.15 |
| \|lam_0\| | 5.4 | 5.7 | 12.9 | 21.3 | 38.6 | 45.9 |
| cond(J) | 3e6 | 7e6 | 2e11 | 1e8 | 3e12 | 2e12 |

-- costates diverging and the shooting Jacobian going singular, while the
conjugate test still PASSES and dim S stays 1 everywhere. So the short
family is not losing local optimality and is not acquiring an abnormal
lift; it is running into a singular Jacobian.

**What finally worked was DIRECT collocation at 70 mN seeded with the FULL
75.5 mN state/control history** (Hermite-Simpson, N = 800, consistent mass
profile, exact endpoints), then the standard covector harvest and
multiple-shooting polish. This is the catalog pipeline's own route; only
the continuation had been failing.

### Corrections to our own framing (GPT-6 Astra consultation, `reviews/mintime_wall_astra_2026-09-08.md`)

* **Loss of all-burn is EXCLUDED analytically, not just empirically.** With
  lam_m(t_f) = 0 and lam_m' = -(T/m^2)|lam_v|, lam_m(t) = int_t^tf
  (T/m^2)|lam_v| ds >= 0, so Q_mt = |lam_v|/m + lam_m/c >= 0 always. A
  coast arc cannot be the missing branch. Measured min Q_mt indeed GROWS
  (1.26 -> 4.22) into the wall.
* **"Stalls at ||R|| ~ O(1) but is not ill-conditioned" was
  self-contradictory** -- a well-conditioned square system has no interior
  nonzero-residual stationary point. Measuring cond(J) settled it: 1e12.
* **Our winding evidence was never valid.** Swept revolutions move only
  0.75 -> 0.77 across the whole family, and the 26.4 d solution has the
  SAME winding (0.76, net -0.22 rev) as the 15.15 d one. The 2026-09-01
  "continuation cannot grow winding" reading is not supported here; whatever
  ends the short family, it is not a winding change.

### Method lesson: fill deep rungs by DIRECT solve, not by continuation

A four-cell cold walk from 0.5 N down a 0.88-ratio ladder (17 rungs, t_f
guess sweep, K = 48) did far WORSE than the banked-seed chain on the very
same cell:

| cell | deepest closed, cold walk from 0.5 N | banked-seed chain |
|---|---|---|
| (1,11) | 143 mN | **75.5 mN** |
| (3,1) | 158 mN | -- |
| (1,9) | 480 mN | -- |
| (3,10) | 416 mN | -- |

So the wall a cell reports is a property of the SEED CHAIN as much as of the
cell: the same cell walls at 143 mN cold and 75.5 mN warm, and the direct
solve then reaches 70 mN outright. **Do not fill the catalog's low-thrust
coverage holes by extending the continuation ladder.** Solve the target rung
directly (Hermite-Simpson NLP seeded with the nearest converged full
trajectory), harvest, then polish -- which is what the pipeline was designed
to do and what closed 70 mN here.

### Open

**Whether 26.436 d is the MINIMUM time at 70 mN is not established.** It is
a certified extremal (PMP root, conjugate-passing, foreign-witness
accepted), but the jump from 15.15 d at 75.5 mN to 26.44 d at 70 mN is
large for a 7% thrust change, and the short family's fate at 75-70 mN is
unresolved. Fixed-t_f feasibility probes were unreliable (they failed even
1.7% below a converged solution, measuring the seed rather than
feasibility). The decisive experiment, recommended by Astra and not yet
run, is **scaled pseudo-arclength continuation in (z8, T)** on the
multiple-shooting unknowns: a fold shows tangent thrust-component -> 0 and
a sign change, with R_X losing one rank while [R_X R_T] stays full rank.

## 33. (arrival phase sA = 0.9087 ONLY -- see section 35) The 71.99 mN limit point -- NOT established as a fold, and the "minimum thrust" reading is REFUTED (2026-09-08; corrected same day after review)

> **CORRECTION.** This section first claimed a simple fold at 71.99 mN and
> concluded that the abstract's 70 mN is "below the family's limit". **Both
> claims were wrong.** The GPT-6 Astra code review
> (`reviews/arclength_code_review_astra_2026-09-08.md`, xhigh reasoning)
> refuted them, one on our own evidence:
>
> 1. **"Minimum thrust for this geometry is 72 mN" is contradicted by our own
>    certified 70 mN solution at the SAME endpoints (section 32).** A thrust
>    turning point on ONE extremal branch is not a feasibility threshold for
>    the problem. Conflating them was a straight logical error.
> 2. **A simple fold is not an endpoint of the solution curve** -- you round
>    it. We never rounded it; the arc approached asymptotically and stalled.
>    That is equally consistent with the normal-costate chart degenerating
>    (multipliers -> infinity as an abnormal configuration is approached), in
>    which there is NO finite corner in these coordinates.
> 3. **18.120 d is the last computed time, not a fold time.**
> 4. Our tangent (solve R_x v = -R_t, normalize [v;1]) is exactly the wrong
>    one near a fold: R_x is singular there. It should be the right null
>    vector of a FULL svd of [R_x R_T].
> 5. The R_T central difference uses h_theta ~ 1e-6 against residual noise
>    ~1e-12, giving derivative noise ~1e-6 in scaled coordinates -- LARGER
>    than the smallest augmented singular value (1.5e-7) we quoted as
>    evidence of regularity. That number may be noise.
> 6. Our step control inflates ds by 1.3x after every success regardless of
>    difficulty, and `A.ds` stores the NEXT proposal rather than the accepted
>    step -- so the recorded history cannot distinguish "arclength steps
>    collapsed" from "steps were fine but went entirely into costates".
>
> What SURVIVES: the 75.0 mN wall was a thrust-stepping artefact (arclength
> passed it with residuals at 1e-12), the roots found along the arc are
> accurate, and section 32's certified 70 mN solution is unaffected.
>
> The decisive test, from the review: along a regular normal branch with
> bounded multipliers, dt_f/ds = (integral of H_T) dT/ds with
> H_T = -(|lam_v|/m + lam_m/c) < 0, so **at a finite normal fold dt_f/ds must
> vanish together with the thrust tangent.** Our log's printed precision
> quantizes both increments to zero, so it cannot decide; a re-run logging
> |lam_0|, max|lam(t)|, the ACCEPTED step length, and the state-vs-costate
> split of each step will. If the arclength is being consumed by growing
> costates, the cure is not a better corrector but a homogeneous PMP
> normalization (rho^2 + |lam_0|^2 = 1, watching whether rho -> 0).
>
> The original text follows for the record.

### (superseded) The fast family has a MINIMUM-THRUST FOLD at 71.99 mN / 18.12 d

Section 32 left the decisive question open: does the fast DRO -> tulip family
fold near 75 mN (making the 26.44 d solution the answer at 70 mN), or does it
continue with a ~16 d transfer we failed to find? Those give OPPOSITE verdicts
on the cislunar abstract. **Pseudo-arclength continuation settles it.**

### The instrument

`arclength_thrust` (new) follows the solution CURVE in (X, T) instead of
stepping in thrust, so it can pass a turning point. Unknowns are the
multiple-shooting vector X = [lam_0; Y_2..Y_K; t_f], parameter is thrust;
everything runs in scaled coordinates (the blocks span costates O(10-50),
junction states O(1) and t_f O(3), so an unscaled arclength would be
meaningless and would change with K). R_T is a central difference at fixed
unknowns. The residual and its Jacobian are the PRODUCTION ones, reached
through a new `ms_bvp` **`assembleOnly`** mode -- continuation needs R and
R_X at points that are not roots, for a family of problems.

### The result: a fold, exactly as specified

Started from the converged 75.495 mN root and walked STRAIGHT THROUGH the
75.0 mN "wall" of section 32 -- residuals stayed at 1e-12. The wall was a
thrust-stepping artefact, not a boundary. The arc then converges on a limit
point:

| step | T (mN) | t_f (d) | sigma_min(R_X) | sigma_min([R_X R_T]) |
|---|---|---|---|---|
| 0 | 72.881 | 16.306 | 3.33e-06 | 7.02e-06 |
| 100 | 72.007 | 17.924 | 1.47e-08 | 7.96e-07 |
| 200 | 71.995 | 18.051 | 1.48e-09 | 4.04e-07 |
| 300 | 71.993 | 18.089 | 4.06e-10 | 2.64e-07 |
| 400 | 71.992 | 18.107 | 1.65e-10 | 1.95e-07 |
| 524 | 71.992 | 18.120 | 7.13e-11 | 1.48e-07 |

**R_X loses rank by five orders while the augmented [R_X R_T] falls by only
two and stays regular -- a rank-one deficiency with the augmented matrix
full rank. The tangent's thrust component goes to zero and thrust
asymptotes while t_f keeps growing.** That is the textbook simple fold, and
it is the criterion stated in advance (Astra, `reviews/mintime_wall_astra_2026-09-08.md`).

### What it means

**The fast family's minimum thrust is T* = 71.99 mN, at t_f* = 18.12 days.**
Below 72 mN this family DOES NOT EXIST. Everything else now follows:

* plain thrust stepping died at 75.0 mN -- 4% above the true limit, because
  a Newton corrector at fixed T cannot approach a fold;
* |lam_0| 5.4 -> 45.9 and cond(J) 3e6 -> 2e12 into the wall were the fold
  approaching, not a winding change (the winding never moved: 0.75 -> 0.77);
* the direct solve at 70 mN jumped to 26.44 d because at 70 mN the fast
  family is gone -- 26.44 d is on a DIFFERENT branch, and section 32's
  certified solution stands as the answer at that thrust.

### Consequence for the cislunar abstract (SUPERSEDED -- see the correction above)

The abstract asks for **70 mN in ~18 days**. Those are not compatible on
this family: 18.1 days is exactly the family's limit, but at **72 mN**, and
70 mN is 2.8% below the fold. The honest options are

1. **72 mN, 18.1 days** -- the family's minimum-thrust point, which very
   nearly matches the abstract's flight time and needs only the thruster
   figure changed; or
2. **70 mN, 26.44 days** -- the certified section-32 solution on the long
   branch, keeping the thruster and changing the schedule.

Option 1 preserves the abstract's headline and is the smaller edit. Neither
supports "70 mN in 18 days".

Reproduce: `arclength_thrust` from the 75.5 mN root
(`results/mintime_arclength_fold.mat`, logs `arclength{,2}.log`).

### Open

The arc CRAWLS at the fold (step length collapses as sigma_min(R_X) -> 0) and
never turned the corner, so the returning branch -- higher t_f at rising
thrust -- is unmapped. Turning it needs either a bordered/deflated corrector
at the limit point or a switch to t_f as the continuation parameter there.
The fold LOCATION is nonetheless established to 71.99 mN by five orders of
rank collapse.

## 34. (arrival phase sA = 0.9087 ONLY -- see section 35) The branch map: the fast family ends at 72.0 mN by LOSING NORMALITY, and 26.44 d is the only locally minimizing extremal found at 70 mN (2026-09-08)

Section 33's "fold at 71.99 mN" was retracted after review. What the
corrected instrument and the homogeneous chart then found is a fuller and
more useful picture than a single turning point. Figure:
`indirect/results/mintime_branch_map.png`; data `mintime_arclength_{diag,hom,70mN_dn,70mN_up}.mat`.

### Instrument (what changed since section 33)

* `arclength_thrust`: tangent from the FULL-SVD right null vector of
  [R_X R_T] (the old R_X \ R_T chart is singular exactly at a fold); honest
  logging of the ACCEPTED step, Newton count, |lam_0| and the costate
  fraction of each step in scaled coordinates; Newton-count step control;
  convergence tested after the last iterate; FD step for R_T moved onto its
  measured plateau (h_rel 1e-3: derivative stable to 1e-10 across three
  decades, sigma_min([R_X R_T]) identical to 4 digits -- the augmented
  regularity evidence was real, not noise); `opts.binding`, `opts.direction`.
* `ms_bvp` gains extra scalar unknowns (`prob.nExtra`, `prob.extraEq`);
  `ms_tfmin_hom` is the HOMOGENEOUS min-time binding: rho free,
  H(t_f) = rho + lam'f = 0, rho^2 + |lam_0|^2 = 1. Tests 7/7 and 6/6,
  fixed-tf 10/10 and golden 20/20 unchanged.

### 1. The discriminator settles the fast family: not a fold, a loss of normality

Re-run in the normal chart with the fixed instrument (401 roots, 75.5 ->
72.0 mN): the identity dt_f/dT = int_0^tf H_T dt holds to FOUR DECIMALS
along the whole arc (ratio 0.9999-1.0000 at every sampled point), so every
root is a regular normal extremal and the derivatives are exact. Along it
|lam_0| grows 46 -> 3409 (x74) with every accepted step at the maximum
length and 100% costate in scaled coordinates, and the slope dt_f/dT
diverges -25 -> -921 in lockstep with |lam_0| (H_T is homogeneous of degree
1 in lam), while T -> 72.00 mN and t_f -> 18.0 d stay finite. A vertical
tangent in (T, t_f) with UNBOUNDED multipliers: not a simple fold, which
has bounded multipliers. The branch reaches an abnormal configuration at
its turning point.

### 2. The homogeneous chart walks through it and maps the snake

On the sphere the multipliers cannot run away, and the arc from the 75.5 mN
root passes the point the normal chart could only approach:

| what | T (mN) | t_f (d) | rho |
|---|---|---|---|
| fast family, start | 75.495 | 15.15 | +0.0218 |
| **fast family ends: rho crosses 0** | **72.02** | **18.45** | 0 |
| abnormal connector (rho < 0, NOT min-time candidates) | 72.0 -> 74.1 | 18.5 -> 25.1 | < 0 |
| slow family, rho > 0 again, descends to its bottom | 70.68 | 43.0 | +0.0007 |
| further folds | 109.3, 83.0 | 26.8, 36.8 | +0.14, +0.011 |

The fast family's termination coincides EXACTLY with rho = 0: the turning
point is the abnormal point. Below 72.0 mN there is no fast solution.

### 3. The certified 26.44 d solution sits on its own branch, near its own fold

Continuing the certified 70 mN root in both directions:

| direction | folds (T mN / t_f d / rho) | lowest thrust |
|---|---|---|
| down | **69.57 / 26.63 / +0.0141** then up to 120 mN at 20.6 d | 69.57 mN |
| up | 71.49/25.8, then rho < 0 patch 70.4 -> 74.2, then 76.6, **59.36/44.4/+0.006**, 83.4, 67.6 | 59.36 mN |

So the 26.44 d point lies 0.5% above its branch's own thrust minimum -- a
GENUINE simple fold this time (rho bounded at 0.014) -- which is why a
direct solve seeded from the fast family landed on it: it was born there.
Upward, the branch snakes through eight folds and reaches 59.4 mN at 44 d.

### 4. Conjugate verdicts on the snake: the slow branches are NOT minima

| point | T (mN) | t_f (d) | rho | conjugate | interior crossings |
|---|---|---|---|---|---|
| certified 70 mN (section 32) | 70.00 | 26.44 | 0.0155 | **PASS** | 0 |
| its branch's fold | 69.57 | 26.63 | 0.0141 | FAIL | 1 |
| lowest thrust found | 59.36 | 44.43 | 0.0063 | FAIL | 2 |
| slow family bottom | 70.68 | 43.04 | 0.0007 | FAIL | 1 |
| 67.59 mN fold | 67.59 | 40.17 | 0.0165 | FAIL | 2 |

The fold point fails while the certified point 0.5% away passes: local
minimality is EXCHANGED across the fold, as it should be. Walked step by
step down the branch, the exchange is located: PASS at every root from
70.000 down to **69.580 mN**, FAIL from **69.569 mN** (the fold itself,
where the thrust tangent changes sign) onward. So the minimizing segment
of this branch is the piece ABOVE its fold, and the 70 mN certified point
lies 0.6% inside it. Every other extremal found at or below 70 mN is both
longer and a saddle.

### What this establishes, and what it does not

* **At 70 mN, 26.436 d is the only locally minimizing extremal in the
  mapped structure**, and the map says WHY nothing faster exists there: the
  fast family cannot reach 70 mN (it loses normality at 72.0), and the
  other branches that do reach it are saddles.
* **"70 mN in 18 days" is not available**: 18 d needs >= 72 mN.
* NOT a global certificate: a disconnected minimizing branch could exist.
  What the map adds over section 32 is that every branch REACHABLE by
  continuation from either known solution has been walked, and none beats
  26.44 d at 70 mN.
* rho on the slow branches is 0.0006-0.02: those extremals are nearly
  abnormal, the objective barely enters their Hamiltonian, which is
  consistent with their being saddles.

### For the poster

The defensible statement: at 70 mN the minimum-time transfer found is
26.4 days (1.14 km/s, 18.1 kg), certified as a local minimum by an
independent solver and the second-order test; the 18-day family requires
at least 72 mN and terminates there. Both numbers are on the map.


## 35. CORRECTION: the abstract's 18 days is RIGHT -- t_f at 70 mN varies ~50% with ARRIVAL PHASE (2026-09-08)

Sections 32-34 judged the cislunar abstract against cell (1,11) of the
12x12 fine sheet and concluded its 18 d / 0.75 km/s / 12.2 kg was
"optimistic by ~50%" and that "18 d needs >= 72 mN". **Both are withdrawn.**
They were statements about ONE arrival phase.

### What the abstract's case actually is

`sweep_phasing`'s defaults have always been the abstract's engine -- 70 mN,
Isp 900 s, 150 kg -- and its ANCHOR is Darin's demo phasing pair. Polished
through our own machinery:

| | anchor (sD 0, **sA 0.0754**) | cell (1,11) (sD 0, **sA 0.9087**) |
|---|---|---|
| t_f | **17.798 d** | 26.436 d |
| dV | **0.7485 km/s** | 1.1360 km/s |
| propellant | **12.20 kg (8.1%)** | 18.12 kg (12.1%) |
| ms residual | 5.2e-12 | 2.4e-12 |
| flown arrival | **0.0000 km** | 0.069 km |
| tfMin acceptance | **\|dz\| = 0 exactly** | 2.1e-10 |
| conjugate | PASS (0 interior) | PASS (0 interior) |
| gates | min\|lam_v\| 3.24, min Q_mt 3.83, dim S 1 | 0.599, 1.174, 1 |

The abstract quotes "approximately 18 days, about 0.75 km/s, 12.2 kg,
roughly eight percent of the wet mass". The anchor gives 17.80 d, 0.7485
km/s, 12.20 kg, 8.1%. **The abstract is accurate to every digit it states**,
and it is a certified local minimum. Stored:
`indirect/results/mintime_70mN_anchor.mat`.

### The real finding: arrival phase is the dominant variable

Same orbits, same engine, same departure phase; only the arrival phase
differs -- and t_f moves 17.80 -> 26.44 d (+48%), dV 0.75 -> 1.14 km/s
(+52%), propellant 12.2 -> 18.1 kg. **Arrival phase, not thrust, is what
sets the cost of this transfer at 70 mN.** For a constellation deployed by
staggered release into distinct phase slots, that spread IS the deployment
envelope, and it is the quantity a poster should show.

### How the error happened (and the warning that was already on file)

Cell (1,11) was chosen in section 32 because it was the FASTEST 0.5 N entry
of the fine sheet -- the best-conditioned start for a deep thrust walk. Its
orbits match the abstract (tau = 1 DRO, 7-petal tulip, 150 kg), which is
what was checked; its PHASES were never compared with the abstract's.
Astra's review had already said it plainly -- "changing arrival phase
changes the problem and can move, remove, or introduce branch folds" -- and
that sentence was recorded in section 34's own "open" list without being
applied to the abstract comparison.

**Standing rule: an operating point is (orbits, engine, DEPARTURE PHASE,
ARRIVAL PHASE). Matching the first two is not matching the case.**

### What sections 33-34 still establish

Everything measured there is correct AS A STATEMENT ABOUT sA = 0.9087: the
fast family at that phase ends at 72.0 mN by losing normality, the
homogeneous chart maps the snake, 26.44 d is the only locally minimizing
extremal at 70 mN there, and minimality is exchanged at the 69.57 mN fold.
That is a real and unusual piece of solution structure. It is not a
statement about the transfer in general, and the 6x6 sweep below shows why
that distinction matters.

## 36. The phase-sheet harness: six defects found by review, and the traversal method is wrong (2026-09-09)

`sweep_phase_mintime` (new 2026-09-08) maps the (departure phase x arrival
phase) sheet at 70 mN. It cost FIVE failed runs from our own bugs before a
GPT-6 Astra review at xhigh reasoning
(`reviews/sheet_harness_review_astra_2026-09-09.md`, 16k reasoning tokens)
found six more and judged the traversal method itself wrong.

### The defects (all fixed)

1. **The certification gates were computed and then IGNORED.** `tolDz` was
   never referenced; a `tfMin` exception left `dz = NaN` and the point was
   still accepted; a failed or missing conjugate verdict was accepted. The
   effective rule was "ms converged AND flown position < 100 km". Every
   stored point happens to carry dz = 0 and conj = 1, so **the table was
   true -- the harness's promise was not.** Now enforced.
2. **The anchor was recorded with FABRICATED diagnostics** (fly = 0, dz = 0,
   cj = 1) without checking that the stored certificate matched the
   requested operating point. Calling the sweep at another thrust would have
   labelled the old solution as an accepted anchor. Now the anchor is
   re-solved through the same gates and the operating point is asserted.
3. **The multiple-shooting warm start was destroyed at every sub-step**:
   `seed_from_z8` re-flew the whole trajectory from the predicted initial
   costates, reintroducing exactly the long-horizon sensitivity multiple
   shooting exists to remove. An MS seed does not need to satisfy the
   continuity defects -- removing them is the solver's job. Now the
   neighbour's junctions are kept.
4. **A failed retry replayed the entire edge** from the original point
   instead of resuming from the last accepted sub-step.
5. **Exact-budget rejection**: a point whose final sub-solve succeeded on
   call `maxSolve` was thrown away. The completion test now precedes the
   budget test.
6. **The reported sub-solve count was the final subdivision count**, not the
   total corrector calls -- so "368 sub-solves" understated the true cost.
   Both are now stored (`NSUB`, `NCALL`).

Also in the continuation tools: `ms_tfmin_hom` converted rho < 0 to the
normal chart via `abs(rho) > 0`. Dividing by a NEGATIVE rho is a negative
multiplier scaling and does not preserve the minimising direction; those are
not min-time candidates. Now gated on `rho > 1e-6` with `info.normalValid`.

### The strategic verdict

Astra's recommendation, which we accept: **replace the traversal, not the
engine.** Breadth-first / spine-and-ribs stepping cannot pass a branch
termination and cannot represent more than one candidate per grid point. The
right tool is **pseudo-arclength continuation in ARRIVAL PHASE over the full
multiple-shooting unknowns**, seeded from BOTH known families (the 17.8 d
fast anchor and the 26.4 d slow solution), tracing each curve in both
directions and recording every intersection with the requested grid levels,
with the cheap departure ribs hung off those intersections.

Arrival phase is an unusually cheap continuation parameter here: it enters
ONLY the terminal state-matching rows, so
`R_sA = [0; -x_A'(s_A); 0; 0]` with `x_A'(s_A) = P_A f(x_A(s_A))` for an
exact periodic orbit -- no finite differencing of the whole residual needed,
unlike the thrust driver.

### What is NOT settled

`cond(J)` reaching 5e10 while |lam_0| FALLS and the conjugate test passes is
not diagnosed. It is definitely NOT the normality loss found in thrust: on
the sphere that segment moves AWAY from rho = 0 (rho 0.0567 -> 0.0679).
Whether it is a fold in arrival phase, a mesh artefact at fixed K = 24, or
the shooting problem stiffening as the target moves around a large orbit
needs the separated conditioning diagnostics Astra lists.

### Standing caution

Missing sheet entries mean "this traversal did not reach them under its
rules" -- not infeasibility, not branch non-existence. Failure categories
(budget exhausted / Newton stalled / normality boundary / local test failed)
must be stored separately, and are not yet.

## 37. The arrival-phase sheet, rebuilt as pseudo-arclength continuation: the engine, the gate stack, and the first branch structure (2026-09-09)

Section 36 ended with the traversal method judged wrong. This is its
replacement, built to the Astra prescription and to TDD throughout: five
new units, five test suites, every one of them RED before it was GREEN.

### The pieces

| unit | what it owns | test |
|---|---|---|
| `costate_common/arclength_ms` | GENERIC pseudo-arclength continuation of `R(p,q) = 0` in any parameter, for any residual with a Jacobian | `test_arclength_ms` 9/9 (unit-circle fold fixture) + `test_arclength_ms_thrust` 5/5 (production regression) |
| `costate_common/cr3bp_field` | the ballistic CR3BP field, so the arrival-phase derivative is ANALYTIC | inside the arrival suite |
| `DRO_tulip/indirect/arclength_arrival` | the min-time DRO -> tulip problem bound to arrival phase in the homogeneous chart | `test_arclength_arrival` 8/8 |
| `DRO_tulip/indirect/certify_root` + `certify_crossing` | the gate stack for one candidate | `test_certify_crossing` 7/7 |
| `DRO_tulip/indirect/sheet_from_arcs` + `build_arrival_sheet` | crossings -> certified sheet | `test_sheet_from_arcs` 7/7 |
| `DRO_tulip/indirect/rib_from_crossing` | the departure-phase rib off a certified crossing | `test_rib_from_crossing` 7/7 |

### Four design decisions that were forced by measurement, not taste

**1. The tangent is a null vector, not a linear solve.** `v = -R_x \ R_q`
is exactly singular AT a fold, which is the one place the tangent matters.
The right null vector of the FULL SVD of `[R_x R_q]` is not. (Economy SVD
drops the extra right-null column of a wide matrix -- use the full form.)

**2. The arclength metric must not scale with the MESH.** With unit weights
the trajectory block's `K x 14 = 336` coordinates swamp the single phase
coordinate: an arclength step of 0.002 moved the arrival phase by 5e-5, and
60 steps advanced 0.0015 of a period. Weighting every junction coordinate
by `sqrt(K)` makes the block's contribution quadrature-like -- one
junction's worth rather than K junctions' worth -- and the same 60 steps
then cover 0.015.

**Measured, because the first version of this paragraph claimed more than
was checked.** The same 40-step arc at K = 24 and K = 48:

| | K = 24 | K = 48 |
|---|---|---|
| phase advance per unit arclength | 0.003413 | 0.002961 |
| arrival phase reached | 0.092202 | 0.088455 |
| `t_f` at the common phase 0.088455 | 17.455539 d | 17.455534 d |

So the SOLUTION CURVE is mesh-independent to 5e-6 days, as it must be --
both meshes discretize the same boundary-value problem. The
PARAMETERIZATION is not: the step advances 13% less phase on the finer
mesh, and step for step the difference reaches 2.96e-4 against a mean step
of 4.20e-4. What the weighting removes is the LEADING mesh dependence --
unweighted, doubling K would cut the phase advance by `sqrt(2)`, 41% --
leaving a 13% residual, plausibly from the junction states sampling the
same trajectory differently plus the fixed weights on `t_f` and `rho`.
"Comparable across meshes", not "identical".

**3. A fold is a rank statement, not a sign.** A sign change of the
tangent's q-component is only a CANDIDATE. It is called a fold when `R_x`
loses rank there (small `sigma_min`) while the augmented `[R_x R_q]` stays
regular -- and it is localized by an in-plane correction, not read off the
nearest accepted root.

**4. Every crossing of a grid level is kept.** After a fold an arc can
cross the same phase again, at a different `t_f`. Both are candidates; the
sheet stores all of them with their verdicts and takes the minimum over the
CERTIFIED ones. A grid point holding two local minima is a fact about the
problem, not a bug in the walker.

### The gate stack (certify_root)

A candidate is certified only if ALL of: the normal-chart polish converges
(tolR 3e-11); the flight from `z8` alone lands within 100 km AND 10 m/s of
the target -- **position and velocity**, the old harness gated position
only; pumpkyn `tfMin` returns within 1e-6 of the polished `z8` (a thrown
exception is a FAIL, not a pass); the free-time conjugate test returns PASS
(ENDPOINT = inconclusive = FAIL); and the min-time hypothesis gates hold
(`min|lam_v| > 0`, `min Q_mt > 0`, `dim S = 1`). Numbers are kept whether or
not a candidate passes, and the FIRST failed gate is named.

Regression: the anchor certifies at **17.7976 d / 0.7485 km/s / 12.20 kg**,
flown miss 0.000 km and 0.000 m/s, witness `|dz| = 0`, conj PASS, dim S = 1,
`min|lam_v| = 3.24`, `min Q_mt = 3.82`. Costates scaled by 3 are refused.

### The thrust regression: the generic engine IS the old instrument

Fed the production thrust residual with the scaling `arclength_thrust`
used, `arclength_ms` reproduces the archived 2026-09-08 diagnostic arc
**root for root**: max relative `t_f` error 0.0e+00 over the reference
roots in the first millinewton, `|lam0| = 73.8855` against the archived
73.8855. The generic engine is not a rewrite that happens to agree; on this
problem it is the same instrument.

### The rib is not symmetric in departure phase

`rib_from_crossing` walking the NEGATIVE departure direction out of the
anchor reproduces the 2026-09-09 certified ring exactly: **17.8775 d at
sD = 11/12** (ring 17.877) and **18.0688 d at 10/12** (ring 18.069). The
POSITIVE direction fails out of the anchor -- as it did in that sweep,
which recorded `(2,1) sD = 0.0833: FAILED from (1,1)`. Two lessons: a
reference number is only meaningful with its direction attached (the first
version of this test asserted the ring values in the +1 direction and the
walker was right to refuse), and a walker that only ever halves crawls at
the depth of its worst patch -- step RECOVERY after two clean sub-steps is
what makes the 1/1728-of-a-period patch affordable.

### First branch structure in arrival phase (arcs in flight)

Four arcs, two seeds x two directions, at 70 mN / Isp 900 s / sD = 0:

- from the **anchor** (sA 0.0754, 17.798 d): walking UP is smooth -- no fold
  in 0.0754 -> 0.24. Walking DOWN folds at **sA 0.0342**, then again at
  0.0706 and 0.0699 on the way back up, so the down direction returns
  through the anchor's own phase on a SECOND branch.
- from the **26.44 d solution** (sA 0.9087): folds at **0.9084 and 0.9103**,
  i.e. that solution sits on a narrow nose barely 0.002 of a period wide,
  with further folds near 0.9293. Both directions have since crossed the
  anchor's phase (0.0754 + 1), so the two seeds' branches will be compared
  at the same grid point.

`t_f` falls with increasing arrival phase out of the anchor: 17.7976 d at
0.0754 to 17.4018 d at 0.0907, agreeing with the independent fixed-step
walk (17.42 d at 0.0899). `rho` stays in 0.053 .. 0.068 throughout -- the
arrival direction does NOT lose normality, unlike the thrust direction of
section 33.

### The continuation's heuristics are guards and pace, not answers (added 2026-09-09)

`foldRatio`, `maxCorrFrac` and `newtonTarget` were CHOSEN, not derived. That
was the one risk of the five listed for outside review that no test covered,
so it was swept: 3 x 3 x 3 settings over two to three decades, on the cubic
whose fold positions are known exactly.

| knob | range swept | effect |
|---|---|---|
| `foldRatio` | 1e-3 .. 1e-1 (three decades) | **none at all** -- identical roots and identical verdicts in all nine (`maxCorrFrac`, `newtonTarget`) groups |
| `maxCorrFrac` | 1.5, 2, 4 | none measurable |
| `newtonTarget` | 2, 4, 8 | the only one that acts, and it acts as intended: it sets how large a step the corrector is asked to afford, hence how far a fixed step budget reaches (nt = 2 takes 61 small steps where nt = 4 takes 41, and stops on `nStep` rather than failing) |

Across all 27 settings a fold is reported exactly when the arc passed it, and
where both were found they sit at -2 and +2 to **7.4e-10**. So none of the
three moves a fold or changes which curve is traced. The reason `foldRatio`
is inert is that the rank collapse at a fold is orders of magnitude sharper
than anything in the range -- on the circle, sigma_min(R_x)/sigma_min([R_x R_q])
is 1.2e-3 at the fold against 8.7e-1 at a regular point.


## 38. The conjugate test is what separates the sheet: at every arrival phase the minimizer is the fastest extremal, and the slower ones are refuted (2026-09-09)

The first sheet assembled from four arrival-phase arcs produced 14
candidates at 70 mN and certified 2. That ratio is not a solver failure --
every one of the 14 converged, and 12 of them were thrown out by the
CONJUGATE TEST alone, after passing the multiple-shooting residual, the
flown arrival and the pumpkyn `tfMin` witness.

| arrival phase | certified t_f | refuted t_f (conjugate verdict 0) |
|---|---|---|
| 0.0754 | 17.7976 (the anchor, seeded) | 22.05, 28.81, 34.20 |
| 0.1587 | 16.2256 | 18.92, 30.24, 35.67 |
| 0.2421 | 16.8742 | 31.51, 36.93 |
| 0.3254 | -- | 38.24 |
| 0.9087 | 26.4361 (seeded) | **26.4537** |
| 0.9921 | -- | 27.47, 32.76 |

Three things this says.

**1. Every refuted candidate is SLOWER than the certified one at its phase.**
The pattern is a ladder: one fast branch and a stack of slow extremals
above it, exactly the structure the thrust-direction study found in section
34, where the slow branches were conjugate FAILs and the fast one was not.
The minimum-time minimizer is the fastest extremal, and the arcs walk
through the others on the way.

**2. The refutations are not sloppy solutions.** They fly to the target to
between 0.00 and 2.5 km with witness agreement `|dz|` of 1e-10 to 3.5e-8.
By every first-order measure they are extremals; the second-order test is
the only thing that separates them. Any pipeline that gates on "converged
and it flies there" -- which is what the old harness effectively did, see
section 36 -- would have shipped all 14.

**3. At the fold nose the test resolves 26 MINUTES.** At arrival phase
0.9087 the certified solution is 26.4361 d and a second root at 26.4537 d
is refuted: two extremals 0.0176 d apart, one minimizing and one not,
separated by the conjugate test alone. That is the nose of section 37 seen
from the other side -- the folds at 0.9084 and 0.9103 bracket the certified
solution, and the arc returning through the phase brings back its partner.

The practical consequence for the sheet: `S.TF` must be the minimum over
CERTIFIED candidates, never over converged ones, and every candidate is
kept with its verdict so this table can be read at all.

## 39. The departure-phase barrier is a CONJUGATE POINT, and the solver's "failure to converge" was that conjugate point announcing itself (2026-09-09)

The departure rib walks one way out of the anchor and not the other. On
2026-09-08 the positive sense was recorded as a solver failure -- the
2026-09-09 sweep log has `(2,1) sD = 0.0833: FAILED from (1,1)` -- and the
first rib test asserted the certified ring's numbers in that direction and
had to be corrected. Three probes turned that "failure" into a result.

### Probe 1: the stall is reproducible and step-independent

Walking the positive sense with a deep bisection budget stalls at
**sD = 0.0465** with `|R| = 5.5e-11` against a `tolR` of 3e-11; walking it in
steps of 1/96 instead of 1/12 stalls at **sD = 0.0467** with `|R| = 5.1e-11`.
Same place, same residual, whatever the step. That is not a step-size
problem: it is the solver's achievable floor at that phase sitting just above
a very tight number.

### Probe 2: with the plateau allowed through, the verdict changes

`certify_root` now sends a polish that plateaus below `tolRelax` on to the
gates (and says so in its reason -- FINDINGS: this also exposed `tolR` being
documented as an option and hardcoded in the call). The walk then gets past
the residual and stalls at **sD = 0.0466 with `conjugate test verdict 0`.**
The barrier is second-order, not numerical.

### Probe 3: the two are ONE phenomenon, measured

A conjugate point is a non-trivial solution of the linearised BVP, so the
multiple-shooting Jacobian must go singular as one is approached. It does:

| sD | t_f [d] | \|R\| | cond(J) | conj |
|---|---|---|---|---|
| 0.0000 | 17.7976 | 2.1e-11 | 4.40e9 | PASS |
| 0.0300 | 17.8963 | 6.4e-12 | 9.20e9 | PASS |
| 0.0380 | 17.9985 | 2.9e-11 | 2.79e10 | PASS |
| 0.0420 | 18.0992 | 7.1e-12 | 8.60e10 | PASS |
| 0.0440 | 18.1820 | 1.8e-11 | 2.35e11 | PASS |
| 0.0450 | 18.2416 | 2.6e-11 | 5.27e11 | PASS |
| 0.0460 | 18.3286 | 6.8e-12 | 2.14e12 | PASS |
| 0.0465 | 18.3973 | 2.3e-11 | 9.05e12 | PASS |
| 0.0466 | 18.4160 | 2.5e-11 | 1.46e13 | PASS |
| **0.0467** | 18.4379 | 1.4e-11 | **2.73e13** | **FAIL** |
| 0.0470 | -- | 1.2e+00 | 1.18e17 | (no convergence) |
| 0.0480 | -- | 2.9e+01 | 3.13e17 | (no convergence) |

`cond(J)` climbs **eight orders of magnitude over 0.047 of a departure
period**, the conjugate verdict flips between 0.0466 and 0.0467, and Newton
stops converging at all three thousandths of a period later. **The conjugate
point sits at sD = 0.04665 +- 0.00005.** The residual plateau of probe 1 was
this singularity, seen a few thousandths early and misread as a solver
limitation.

### What it means for the sheet

At 70 mN and arrival phase 0.0754, the branch continued from the anchor is
minimizing from the conjugate point at sD = 0.0467 backwards through 0 and
round through the negative sense.

> **Coverage figure corrected within the hour.** This section first said
> "roughly 63%", stopping at sD = 0.4167 because the 2026-09-08 sweep stopped
> there (`(7,1) sD = 0.5000: FAILED from (8,1)`). **That was the old walker's
> limit, not the problem's.** The rebuilt rib -- step recovery after two clean
> sub-steps, deeper bisection, and a polish plateau no longer read as failure
> -- walked straight through it and certified sD = 0.5000 (18.6369 d),
> 0.4167 (18.7501), 0.3333 (18.9106) and 0.2500 (19.0529). Coverage is at
> least **~80%** and the final figure follows the rib's own stop. Taking a
> previous run's stopping point as a property of the problem is exactly the
> mistake this section is about, and it was made twice in one page.

Whatever the final figure, the uncovered phases are not "unsolvable" -- they
are **not minimizing ON THIS BRANCH**, and whether a different certified
extremal covers them is exactly what the rest of the sheet is for.

Two lessons that generalise beyond this row.

1. **A residual that will not go below a floor, at a fixed place, independent
   of step size, is evidence about the PROBLEM, not the solver.** Loosening
   the tolerance was right, but for the opposite of the obvious reason: it did
   not rescue the walk, it let the walk report its real verdict.
2. **`cond(J)` is a leading indicator of the conjugate test.** It rises four
   orders of magnitude before the verdict flips. A continuation that watches
   it knows it is approaching a loss of minimality well before the test fires
   -- cheaper than the test, and available at every step.

## 40. Two outside reviews in one day: 84 findings, and in every confirmed case the instrument was right while the argument for it was wrong (2026-09-10)

The 70 mN catalog was built, packaged and shipped-ready before either review.
Both were commissioned deliberately: the first on the CODE CHAIN, the second
on the MATHEMATICS. Between them they found 84 things. **No catalog number
changed.** What changed is what we are entitled to say.

### Review 1 — the code chain (41 findings, 26 fixed)

Verdict: *"I would not release this with an unqualified CERTIFIED claim yet."*
It named three classes and all three were real: wrong-phase labelling, export
paths that trusted summaries instead of certificates, and an unenforced
foreign-solver convergence requirement. Adjudication:
`reviews/chain_review_adjudication_2026-09-10.md`.

The one I would keep in mind: **`dR/dsA` differentiated the CR3BP field while
the residual's target is a spline INTERPOLANT.** Different functions.
Correcting it moved the agreement from 7.4e-6 to 9.0e-13. My own test had
printed that 7.4e-6 and I recorded it without asking why it was not machine
precision.

### Review 2 — the mathematics (43 findings: 4 WRONG, 4 GAP, 3 FINE)

All four WRONG are now repaired, and two of the repairs left the document
**stronger** than before.

| finding | resolution |
|---|---|
| mass-adjoint sign | **DOCUMENT ONLY** -- code checked against the flight to ten digits, `max\|H\|` 3.3e-8. Corrected. |
| Proposition 1's "so" | **Rederived.** The construction has a PERMANENT kernel `w = lambda(0) + (1/b) e_m`, proved by our own two lemmas, so it is non-immersive at EVERY time. Rebuilt from the kernel, the same `[Phi_rv P, f_rv]` falls out. |
| scope: all-burn is not all competitors | **Re-founded on H3.** `H` is affine in throttle with slope `-T Q_mt`, so `Q_mt > 0` makes full throttle the unique minimiser with a margin, and a strictly bang component can be eliminated. H3 was carrying the reduction all along. |
| normality "iff" | **Downgraded to a sufficient exclusion**, which is the direction we use. The claim that abnormality needs `lam_v = 0` is simply false and is withdrawn. |

### H6, a new gate that came out of the fourth objection

In the reduction the reduced Hamiltonian is not conserved and `p'J = 0`, so
`det = 0` can mean `h(t) = 0` rather than a rank drop -- a spurious zero none
of H1-H4 excluded. In closed form
`lambda_rv . f_rv = -1 + (T/c) lambda_m`, zero exactly at `lambda_m = c/T`,
and `lambda_m` decreases monotonically to zero. So

> **H6: the spurious mechanism cannot fire iff `lambda_m(0) < c/T`.**

Measured over the shipped catalog: `c/T` = 49.38, `lambda_m(0)` in
[2.11, 10.88] across all 53 entries, **4.5x headroom, none can fire it.**

### Two exact identities

`d/dt(lam_m m) = -T Q_mt` with `lam_m(tf) = 0` and `m(0) = 1` gives

> **`integral of T Q_mt dt = lambda_m(0)`** -- verified 5.5004 vs 5.500404.

The total strict-bang margin IS the initial mass costate, so H6 also reads
"that margin must stay below `c/T`". One scalar ties H3, transversality and
the mass costate together.

### Two new instruments

- **`conj_spectrum`** -- dense singular-SPECTRUM scan on the instrument's own
  variational integration. Closes both blind spots of the sign test (two
  crossings in one segment; even multiplicity, where several values collapse
  and no sign changes). Its test runs a certified and a refuted entry
  together: 0 interior crossings against 1 at t/tf = 0.9219.
- **`lift_margin`** -- the dim S rank statement as an Eckart-Young MARGIN
  against a measured error, replacing a threshold whose claimed safety
  property it did not have.

### The lesson

When review 1 said the equivalence was "asserted, not established", I answered
that it was proved in a document the reviewer had not been given. **Both were
true.** The document existed and its central step did not follow, and I had
been citing it for weeks. *Having a proof written down is not the same as
having checked it.* The five rules in `doc/CERTIFICATION_DISCIPLINE.md` now
carry this one.

## 41. Second script review (GPT-6 Astra, xhigh): the sampled minimum-principle check was a second tautology, and a computed gate was not an enforced gate (2026-09-10)

The revised `transfer_study.m` and its ten callees went back to Astra with a
sharpened brief (my own doubt stated: N6 looked like a sampled restatement of
an analytic fact). Verdict: 11 CORRECTNESS, 6 READABILITY, 5 FINE, 3
OVERCLAIM, 1 REDUNDANT; $1.13. Full text in
`reviews/transfer_study_review2_astra_2026-09-10.md`. Everything actionable on
files the running sweep does not hold is fixed; the rest is in TODO with the
reason.

### The two findings that mattered

**N6 was a tautology, again.** The 400-point sphere sample compared the flown
direction `-lam_v/|lam_v|` against `H` evaluated on random unit vectors. `H`
is affine in the direction with coefficient `(T/m) lam_v`, so its minimiser
over the sphere is `-lam_v/|lam_v|` *by construction*; the sample could only
ever return zero. The replacement asks the question that CAN fail: is the
control the propagator **applied** the minimiser of OUR Hamiltonian? The
thrust acceleration is recovered as (powered field - coasting field), and the
exact gap `H(applied) - min H = (T/m)(lam_v . alpha_applied + |lam_v|)` is
evaluated directly: 7.5e-15 on the anchor, throttle 1 to 5.8e-15. The
mutation test injects a vector field that thrusts against the primer, and the
gap opens to 1e-2 while the shooting residual would happily converge to that
wrong problem. This is `pmp_pointwise_checks` (with N2, N4, N5 beside it);
the certifier now runs it on every flight and FAILS on it.

**H6 was computed and then ignored.** `mintime_hypothesis_gates` returned
`h6Ok`, `certify_root` stored the gates, and nothing read the flag: a library
entry could certify with `h6Ok == false`. Now gate 6 of the stack, with a
margin requirement (`h6MarginMin`, default 1x) so "excluded, but with no
headroom" is distinguishable, and a missing margin is a named failure. The
sweep so far reads 17-33x on every entry, so no shipped entry is affected --
which is luck, not design.

### Everything else that changed

- `validate_flight` (new, `costate_common`): ONE admissibility check --
  reached `t_f`, finite, all-burn mass law, clear of both primaries -- used
  by the certifier, the witness flight, `verify_with_pumpkyn` and the script.
  Before, the script and the certifier each had their own partial version.
- `verify_with_pumpkyn` has an OVERALL status: usable answer AND costates
  did not move AND its own flight is admissible and reaches the target in
  position AND velocity. Agreement alone was one metric of three.
- `report_optimality` has three groups (NECESSARY / SUFFICIENCY / CROSS-CHECKS)
  and four line states (PASS / FAIL / NOT CHECKED / UNRESOLVED). A failed
  cross-check no longer reads as "not an extremal"; an ENDPOINT conjugate
  verdict propagates as UNRESOLVED rather than as a refutation; the claim
  wording is one sentence shared with the script.
- `transfer_study`: thresholds in one block before anything runs; the
  periodic interpolant is REQUIRED (no silent `spline` fallback) and its seam
  is checked in value and derivative -- derivative 4.5e-16 against the 1e-2
  jump a not-a-knot spline has there; the value seam is the orbit's own
  closure (4.6e-8 on the tulip), which the interpolant cannot improve; orbit
  closure is asserted, not printed; N5 uses two coordinate-scaled step sizes
  (agreement 4.6e-9); diagnostic IDs are stable and grouped (N1-N6 Pontryagin,
  S1-S4 theorem, V1-V2 instrument validity, X1 cross-check); coverage is
  reported as t/t_f; the one `flight` object feeds every section and the
  figure (`plot_transfer_3d` accepts it and reports `flightSupplied`); review
  history is out of the comments.
- Tests: `test_validate_flight` (6), `test_pmp_pointwise_checks` (10, with
  the wrong-sign mutation), and extensions to `test_certify_crossing` (H6
  refusal), `test_report_optimality` (cross-check isolation, UNRESOLVED, H6
  NOT CHECKED), `test_verify_with_pumpkyn` (overall status),
  `test_plot_transfer_3d` (supplied flight). All green; the script re-run
  reaches the same anchor to the digit.

### Deferred, and why

`conj_spectrum` (drop the "det sign is meaningless" claim, report the
uncovered final interval, say "candidate-detection scan", fix the `multTol`
inequality), `lift_margin` ("rank >= 6" wording, reject zero or non-finite
`lam`), `h6_margin` (numerical clearance against the Hamiltonian residual) and
the `hMin` rename in `mintime_hypothesis_gates` all live in files the
second-order sweep is executing right now. Editing a function under a running
MATLAB job changes what its later calls do. They go in when the sweep flag is
set.

### One more front door

`build_70mN_library.m` runs the whole library chain -- anchors, four arcs,
sheet, ribs, package, audit, sweep, pictures, deliverable -- as a script in the
`transfer_study` style, each stage a switch so hours-long stages are reused
from their files. Dry-run on the existing results reproduces 115 entries /
115 of 144 cells / audit reused, end to end.

### Addendum: the sweep's multiplicity flags, read while it runs (88 of 115)

Five entries so far carry `multiplicity 1` from `conj_spectrum`: (8,1),
(9,1), (2,4), (3,8), (5,8), all with 0 interior sign changes. The sidecar
stores no location, so the flags were re-derived:

| entry | flagged sample t/t_f | sigma_6/sigma_5 there | reading |
|---|---|---|---|
| (8,1), (9,1), (5,8), (3,8) | 0.990-0.995 | 0.10-0.24 | the graded ENDPOINT collapse (FINDINGS 39-40): the last samples before t_f, where several singular values contract together |
| (2,4) | 0.573 | 0.18 | INTERIOR -- the case the instrument was built for |

The interior one was refined 16x (`nSub` 8 -> 32 -> 128). Its sampled
minimum went 1.23e-6 -> 6.68e-7 -> 6.68e-7 (relative 5.1e-4 of the median),
identical at the two finest resolutions, with no determinant sign change in
[0.50, 0.65] at any resolution. A true zero's sampled minimum keeps falling
with the spacing; this one PLATEAUS, so it is a finite near-miss of the
smallest singular value, not an even-order conjugate point. (3,8)'s
unflagged interior dip at 0.442 behaves the same way (8.97e-7 -> 8.59e-7).
The control (1,1) has no dip in that window at all (relative 3.3).

Two things follow. (2,4) sits at arrival phase 0.3254, the column whose
first candidate the sheet REFUTED by the conjugate test: the certified
entry is close, in the phase plane, to where the conjugate structure
changes, and its 5e-4 near-degeneracy says so. And `conj_spectrum` must
report the LOCATION and the refinement behaviour of every candidate, not a
count: as it stands the writeback will stamp `multiplicity 1` on five
entries whose flags mean two different things. That fix waits for the sweep
to release the file (TODO, top item); the count field should be read as
"candidates, see FINDINGS 41" until then.

## 42. The sweep read out, and the two instruments it corrected: the lift's "measured error" was the loose setting's own error, and a candidate scan must locate what it finds (2026-09-10)

The first second-order sweep completed 115/115 and wrote back: 0 interior
sign changes on every entry, worst H6 margin 4.5x, and **seven entries
uncertified on the lift margin** (4-9x against the 10x bar): (1,10) and six
of the seven entries in the 26-day column sA = 0.9087.

### The lift margin was measuring the wrong thing

`lift_margin` certifies dim S = 1 by Eckart-Young: sigma_6 must exceed the
MEASURED error in the constraint matrix, and that error was measured as the
difference between builds at integration tolerances 1e-10 and 1e-7. On
(1,10):

| setting pair | sigma_6 | sigma_7 | error estimate | margin |
|---|---|---|---|---|
| 1e-10 / 1e-7 | 0.988 | 1.6e-6 | 2.7e-1 | 3.7x, uncertified |
| 1e-12 / 1e-9 | 0.988 | 1.5e-6 | 2.3e-2 | 43x, certified |
| 1e-12 / 1e-10 | 0.988 | 1.5e-6 | 1.4e-3 | 718x, certified |

sigma_6 is 0.99 -- the rank statement could hardly be safer -- and the
"error" falls a decade per decade of the LOOSE setting. The two-build
difference is dominated by the looser build's own error, so with a 1e-7
second setting on a 26-day arc it was not an error estimate for the matrix
in use; it was a measurement of how bad 1e-7 is. The sweep now uses the
pair [1e-12, 1e-9] (`opts.relTolPair`), which is still conservative for the
1e-12 build it certifies. The seven "uncertified" verdicts were the
instrument, not the entries.

### The candidate scan now locates, classifies and refines

`conj_spectrum` returned a COUNT of "multiplicity" candidates; FINDINGS 41's
addendum found that count mixing the start-up transient, the graded endpoint
collapse and one interior near-miss. It now returns `candidates` -- each
with t/t_f, sigma_6/sigma_5 and a class (start / endpoint / interior) -- and
refines every interior one by re-integrating its window at 4x the sampling
from the stored state and STM. Discrimination, measured:

| entry | candidate at t/t_f | refined / coarse minimum | kind |
|---|---|---|---|
| refuted candidate at sA 0.0754 (22.05 d) | 0.92, with a sign change | 0.06 | ZERO |
| certified (2,4) | 0.573 | 1.00 (6.68e-7 at 32x and 128x alike) | near-miss |

A second gap surfaced while testing: the refuted entry's crossing was not a
candidate at all, because a simple zero straddled by two samples need not
dip below the depth threshold at either. A determinant sign change between
inner samples is now a candidate by definition. Multiplicity is counted on
INTERIOR candidates only, and the catalog writeback carries
`conj_interior_cand`, `conj_near_miss`, `conj_zero` beside it.

### The rest of the deferred review items

H6's clearance is judged against the arc's own Hamiltonian residual
(`h6_margin` `opts.Hresid`, wired from `mintime_hypothesis_gates`), `hMin`
is `hMax` (it is the largest value of the reduced Hamiltonian, attained at
t = 0), `lift_margin` refuses a zero or non-finite lift, and its certified
reason reads "nullity one numerically supported (rank >= 6)". Tests: five
instrument suites green plus a new one-entry `test_second_order_pass`.

The corrected re-sweep is running from a fresh sidecar; its census and the
redrawn torus follow.
