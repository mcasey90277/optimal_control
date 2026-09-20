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

### Addendum: the corrected re-sweep, and the near-miss RIDGE across the torus (2026-09-11)

Census 115/115, written back: 0 interior sign changes everywhere, every
lift certified with the tight pair (worst 28x -- the seven "uncertified"
verdicts are gone), worst H6 4.5x. **42 entries carry an interior
candidate**, and they are not scattered: they form a RIDGE. In arrival
columns 4-7 (sA 0.3254-0.5754) nearly every certified row shows one
near-degeneracy of the conjugate matrix, and its location walks smoothly
with the phases -- along a column it moves earlier as the departure phase
advances (column 4: t/t_f 0.573 at row 2 down to 0.510 at row 11; column
7: 0.458 down to 0.401), and across columns it moves earlier as the
arrival phase advances (row 3: 0.568 / 0.526 / 0.495 / 0.458 in columns
4-7). Column 9 carries a second, late one at t/t_f = 0.917 on rows 5-7.
Column 4 is the column whose first sheet candidate the conjugate test
REFUTED (FINDINGS 38), so the ridge is the certified sheet running close to
the surface in phase space where a conjugate point crosses the arc.

Every one of the 42 is a near-miss by hand refinement, and the four the
sweep labelled "zero" -- (6,5), (9,6), (10,7) at ratio 0.49 against a 0.5
threshold, and (3,7) at 0.29 -- all PLATEAU at the second level:

| entry | 8 / 32 / 128 samples per segment | second-level ratio |
|---|---|---|
| (6,5) | 4.71e-7 / 2.32e-7 / 2.32e-7 | 1.00 |
| (9,6) | 1.58e-6 / 7.69e-7 / 7.69e-7 | 1.00 |
| (3,7) | 6.33e-7 / 1.86e-7 / 1.49e-7 | 0.80 (a simple zero gives ~0.25) |
| refuted control | sign change | 0.06 at the first level |

So a single refinement level cannot separate a shallow near-miss from a
zero when the first step happens to land at 0.49. `conj_spectrum` now
refines TWICE (4x, then 16x) and calls a candidate a zero only when both
levels fall, or a sign change is present, or the minimum reaches 1e-8 of
the median. The 42 entries were re-measured under that rule; the
writeback below is the final one.

**Final writeback (2026-09-11):** the 42 entries re-measured under the
two-level rule -- **44 interior candidates, 44 near-miss, 0 zero**; worst
lift 28x, worst H6 4.54x, 0 interior sign changes. The catalog's
second-order fields (`conj_interior`, `conj_interior_cand`,
`conj_near_miss`, `conj_zero`, `h6_margin`, `lift_margin`) are the ones a
recipient should read; `phase_torus_findings.png` is redrawn from them.

## 43. The chain was not safe to rerun: three defects fixed before the live rerun (2026-09-11)

The TODO's top item is the live rerun of `build_70mN_library.m` (stages on,
not the reuse path). Reading the chain before running it found three ways it
would have damaged or misread the shipped 70 mN catalog. None changed a
catalog number; all three would have on the next run.

| defect | what would have happened | fix |
|---|---|---|
| packaging overwrites the catalog, and the shipped switches were package ON / sweep OFF | the second-order writeback (conj_interior, near-miss, zero, H6, lift) erased: catalogs are gitignored and `.bak_2nd` predates the writeback | `guard_catalog_overwrite` refuses by name unless the sweep stage is on, and backs up whatever it overwrites |
| the chain's sidecar pointer was `second_order_progress.mat` | that is the FIRST sweep's file: lift margins differ from the catalog's by up to 1.35e4, and it lacks the candidate fields. The v2 file matches the catalog on all six written-back fields, 115/115, to 0 | pointer is `second_order_progress_v2.mat` |
| the sidecar was POSITIONAL | a re-packaged catalog with a different entry set would inherit other entries' measurements on resume, silently | every record carries its cell key and z8; a resumed record must match (`second_order_pass:staleSidecar`); a pre-key sidecar is adopted only with `adoptLegacy` and only if every written-back value matches |

Also: `chainOverrides.outDir` / `.run` let a batch driver rebuild BESIDE the
shipped files, so a rerun is compared rather than trusted, and the catalog's
`second_order.meaning` now says two-level refinement (4x, then 16x), which is
what `conj_spectrum` does since FINDINGS 42.

Tests `test_second_order_sidecar_identity` (5) and
`test_guard_catalog_overwrite` (5), RED before GREEN; four mutations of the
new refusals all caught, files restored md5-identical. Commit `5df657c`.

The pattern is FINDINGS 40's and CERTIFICATION_DISCIPLINE rule 5 again, one
level up: the builder was fixed and audited, but the SCRIPT that drives the
builder had never been run on its live path, and its defaults were the
dangerous ones.

### The front door knew 10 of the 115 entries

`run_dro_tulip` seeded itself from `dro_tulip_library`, which gathers
solutions from the anchor and sweep FILES: 10 phase pairs. The other 105
certified, audited catalog entries were invisible to it, so asking for one
started a continuation walk -- minutes to hours -- to re-derive a transfer
already on disk. `dro_tulip_library` now lists the catalog on request
(`includeCatalog`; off by default, because `build_arrival_sheet` seeds from
the default list and seeding a sheet with its own output would be circular),
and the front door rebuilds the junction states it needs from z8 with
`seed_from_z8`. A catalog-only pair, (sD, sA) = (1/12, 0.2421), is now served
from the library and certifies in 32 s (pool start included) at the
catalog's own t_f, 18.603887 d, to the printed digit.
`test_run_dro_tulip_catalog` (6 checks) RED before GREEN; three mutations
(the listing branch, the front door's request, the z8 seed) all caught,
files restored md5-identical; `test_run_dro_tulip` still passes (173 s).
`transfer_study.m` followed the same day at Mike's request: section 4 asks
for the catalog and seeds from z8 when the entry has no junction states.
Run at (1/12, 0.2421) it reaches the full verdict (all required numerical
checks passed) at t_f 18.6039 d in 23 s; the anchor pair still takes the
file-backed route (17.7976 d, 10 s); an off-grid pair still refuses, now
under the identifier `transfer_study:noSeed` with its line breaks rendered
(the old `assert` printed a literal `\n`).

### The live rerun reproduces the shipped catalog exactly

`build_70mN_library.m` on its LIVE path -- sheet, package and audit on,
outputs to `results_rerun/` via `chainOverrides`, the shipped files untouched
-- ran 46 min in a clean `-batch` session (R2026a) under an OS watchdog and
ended `CHAIN DONE`. Against the shipped catalog:

| check | shipped | rebuilt |
|---|---|---|
| entries | 115 | 115 |
| cells present in only one | -- | 0 |
| max \|dz8\| over all 115 | -- | **0** (bitwise) |
| conjugate PASS | 115 | 115 |
| audit | 115 ok / 0 bad | 115 ok / 0 bad |
| certified arrival columns / candidates / certified | 11 / 32 / 11 | 11 / 32 / 11 |

**What this proves, and what it does not.** The 11 arrival-spine entries
were RE-CERTIFIED from the arcs by the gate stack and came back bitwise
identical, so the sheet stage is deterministic and the gate stack's current
code certifies exactly what shipped. The other 104 entries came from the
saved departure-rib files (`run.ribs` stays off -- the ribs take hours), so
for them the match proves the package and audit stages, not a re-walk of the
departure axis. Re-walking the ribs is the one piece of the live path still
unexercised.


## 44. The endpoint interpolant becomes a library function, and the two private copies disagreed (2026-09-11)

Mike's read of `transfer_study.m`: the private `periodicPP` should be its own
file, and probably a library function; and the four lines that build the two
phase-to-state closures should be a function too. Both were right, and the
evidence was stronger than the suggestion.

**There were already two copies, with different policies.** `transfer_study`'s
`periodicPP` REFUSED to run without the Curve Fitting Toolbox;
`arclength_arrival`'s `makePP` fell back to an ordinary spline IN SILENCE,
with a companion `ppDer`. Same rule, two files, opposite failure modes -- the
pattern `oclib/README.md` was written about.

**What the ordinary spline actually costs** (measured, tau = 1 DRO and
7-petal tulip, position):

| phase | DRO | tulip |
|---|---|---|
| 0 .. 0.9 | below a millimetre | below a millimetre |
| 0.999 | 3.2 m | 0.21 m |

So in VALUE the not-a-knot spline is harmless except beside the seam. The
difference that matters is the DERIVATIVE, which jumps by ~1e-2 there -- and
that is exactly what `dR/dsA` differentiates in the continuation.

**The two new units** (in `costate_common`, not `oclib`: both consumers are
inside `orbit_transfer`, and `oclib`'s admission rule wants a second
TOP-LEVEL consumer):

- `periodic_pp(t, y, opts)` -- the C1-periodic cubic, its derivative from its
  own coefficients, and two DIFFERENT seams: `seam.value` is the orbit
  table's own closure (no interpolant can improve it), `seam.deriv` the
  interpolant's derivative mismatch. The toolbox policy is now
  `opts.onMissing`, and the fallback warns.
- `phase_state(t, y, opts)` -- one orbit in, two closures out: the state at a
  phase FRACTION and `dx/ds`, which carries the period. This is the unit the
  four lines wanted to be, and the unit a dozen `interp1(..., 'spline')`
  sites across the catalogs should migrate to, one at a time.

**Equivalence, bitwise.** Four vectors captured from the pre-move closures --
arrival state at 0.0754 and at 0.9991 (next to the seam), the arrival phase
derivative at 0.3137, departure state at 0.9991 -- are reproduced with max
difference **0**, and are now the gate in `tests/test_phase_state`. The
production regression `test_arclength_ms_thrust` reproduces its archived arc
root for root (max rel err 0.0e+00, |lam0| 73.8855) and
`test_arclength_arrival` still measures `dstateA` against the interpolant
(2.6e-10) and the analytic `R_sA` against a finite difference (9.0e-13).
`transfer_study` reaches the same verdicts at the same t_f (18.6039 d
catalog pair, 17.7976 d anchor).

**A test threshold that was wrong, caught by its own failure.** The first
interpolation check asserted an error below 1e-7 at 60 intervals, where a
cubic's own error is ~3e-7. It now checks the ORDER instead -- halving h cuts
the error by 15.1, against the 16 a fourth-order method gives -- which is a
statement about the method rather than a magic number.

Three mutations (the periodic scheme, the derivative coefficient rule, the
per-phase scaling) each caught by the right test; files restored
md5-identical.

### Addendum: the backlog is closed, 18 files, three tiers (2026-09-11)

Mike: "can we fix that and have those sites now call our new library
functions." Done the same day, in tiers, each with its own check.

| tier | files | check |
|---|---|---|
| live instruments | `second_order_pass`, `conj_catalog_pass`, `gates_catalog_pass`, `audit_phase_catalog` | audit 4/4 with the shipped values; gates recomputed on 3 measured cells differ by **0** in min\|lam_v\|, min Q and dim S; conjugate verdicts on 2 cells unchanged; `test_second_order_pass` green |
| ladder engines + phase sweeps | `thrust_ladder_library` (the halo and DPO campaigns call it unmodified), `extend_thrust_ladder`, `densify_ladder`, `lowthrust_ladder`, `probe_deep_rungs`, `probe_abstract_case`, `sweep_phase_mintime`, `sweep_phasing`, `ms_refine_catalog`, both direct phase sweeps | Code Analyzer message-for-message against the committed versions (only difference: one alignment note whose line number shifted by the inserted lines); `golden_cells` 20/20 |
| GTO campaign tools | `audit_gto_entries`, `gto_entry`, `viz/gto_pilot_movie` | same message-for-message check, 0 new |

Three things worth keeping:

**The interpolant is now built once, beside the orbit tables, instead of per
entry inside the loop.** Several engines were re-fitting a spline through a
1328-point table for every cell.

**Orientation was preserved site by site, not normalized.** Some sites wanted
a row and some a column; a `.'` at the call keeps each exactly as it was.
Normalizing "while we are here" would have put a silent shape change into
engines that cannot be cheaply re-run.

**A first verification appeared to show verdicts moving** -- dim S off by 2,
a conjugate verdict off by -2. The check was wrong, not the code: the sweeps
walk the grid in their own order, so indexing by `find(has_solution)`
compared measured cells against cells still holding their -1 initializer.
Masking by the sweep's own attempt counter, every measured cell agrees
exactly. A difference of exactly 2 against a -1 initializer is worth
recognising on sight.

**What is deliberately NOT migrated**, recorded so it is not rediscovered as
a defect: the recipient-facing helpers that ship beside a catalog, and the
recipe strings inside the packagers. A shipped catalog's helpers carry no
dependency on `costate_common` (GTO catalog README: "inlined, per
deliverable-picker convention"), and the measured cost of the ordinary
spline there is millimetres, only beside the seam.

## 45. Two more duplications become library functions: the propulsion conversion and the seed cut (2026-09-11)

Mike, on the FINDINGS 44 leftovers: "if it makes sense to fix those two then
let's genericize and use library calls." Both made sense, and both were
bigger than the first estimate.

### The ND propulsion conversion (`nd_propulsion`, 28 files)

Three expressions -- ND gravity, ND exhaust speed, ND thrust acceleration at
unit mass fraction -- written out by hand in about twenty files, one of
which already commented that it used "the same ND thrust/exhaust conversion
`thrust_ladder_library` uses". Because every copy was textually IDENTICAL,
the library can hold the expressions character for character and every
migrated site keeps its value BITWISE. That is the whole equivalence
argument, and it is stronger than any tolerance.

Two closures came out of the call sites rather than from design:
`ndT` (thrust to T_nd) for the rung ladders, and `ndC` (Isp to c_nd) because
`probe_abstract_case` continues over Isp. Either input may be omitted:
`arclength_thrust` receives c as an argument and needs only the thrust
closure.

**Not migrated, and why:** `cr3bp_common/cr3bp_lt_params` and the GTO
`direct/` and `indirect/` campaigns do not have `costate_common` on their
path (only `GTO_tulip/catalog` does). Routing them here would make one
shared library depend on another -- an architecture call, not a refactor --
and those campaigns are already deduplicated LOCALLY through
`cr3bp_lt_params`. The inverse conversions (c_nd to Isp, in `catalog_schema`
and `build_minfuel_catalog`) are a different rule.

### The seed cut (`flight_to_junctions`, 7 files)

Cutting a flown trajectory into K+1 junction states, with the mass row
DERIVED from the all-burn identity when the seed is for a different thrust
or t_f. Six engines held it inline; `seed_from_z8` held a seventh copy.

**The one thing that was not bitwise, measured rather than assumed.** Six
sites queried the interpolant in NORMALIZED time; `seed_from_z8` queried in
ABSOLUTE time. On a real 70 mN flight the two differ by **3e-13 absolute /
5.5e-14 relative, all of it in the costate rows**. There is no formulation
that reproduces both, so the majority form was adopted and the difference
was carried into `golden_cells`, where it shows as the shooting residual
moving 5.6e-14 to 1.2e-13 against a 1e-10 bar with iteration counts
unchanged. A seed perturbation is not a solution change: the root is set by
the boundary-value problem, not by the guess.

### Two habits that earned their keep

**The orphan sweep.** After each migration, grep the edited files for the
variables whose definitions were removed. It found a packager still using
`g0` two lines below its deleted definition -- the kind of thing that runs
fine until the one branch that touches it.

**Code Analyzer message-for-message against the committed version**, rather
than "is it clean now". It caught the dead `unique` calls and grid vectors
the migration left behind in five engines (deleted), and it distinguished a
genuinely new message from the same message with a shifted line number.

Both of my first test thresholds were wrong, again, and their own failures
caught them: an absolute bar where rescaling costs one ulp, and a 1e-7
interpolation bar where a cubic on 60 intervals gives 3e-7.

## 46. The study script's seed lookup becomes two functions, and the third copy goes (2026-09-11)

Mike, reading `transfer_study` section 4: those lines should be a function,
"maybe just cosmetic, based on my human need for readability". It was not
cosmetic. That exact construction existed THREE times -- inline in the
script, as `run_dro_tulip`'s private `seed_of`, and inside
`build_arrival_sheet`'s seeding loop.

Two layers, because two different things were bundled in those 24 lines:

| unit | owns |
|---|---|
| `costate_common/seed_from_entry` | the CONSTRUCTION: a file-backed entry brings its own junction states (reused, final column appended, column 1 pinned to the actual departure state and the entry's costates); a catalog entry carries z8 only and is flown (`seed_from_z8`) |
| `DRO_tulip/indirect/dro_tulip_seed` | the LOOKUP: exact phase match, the operating point in all six fields, and the refusals |

**The refusal got better by being moved.** The inline version asserted "no
seed for this operating point and phase pair" for both failures at once, so
it could not say which had happened. There are now two identifiers --
`:noSeed` for the phase pair and `:operatingPoint` for the engine -- and the
test perturbs each of the six fields ALONE to prove each one refuses.

**Why the match stays exact, written into the function's header so it is not
"fixed" later.** A neighbouring phase would usually converge. At one arrival
phase this problem carries a ladder of extremals: a sheet column held 14
candidates, of which the conjugate test refuted 12 (FINDINGS 38). A
neighbour seed can land on a slower branch and pass every first-order check.
Walking there properly is `run_dro_tulip`'s job.

The script now reads `4. GET THE SEED` / `5. SOLVE`, and sections 5-8
renumbered to 6-9 (the N1-N6 / S1-S4 diagnostic IDs are independent of
section numbers and did not move). Section 4 went from 24 lines to 6, and
the operating-point check stays VISIBLE in the printed seed line rather than
disappearing into the function.

Verified: both entry routes bitwise against the inline constructions they
replace; `test_dro_tulip_seed` (6 checks) and `test_seed_from_entry` (6)
RED before GREEN; Code Analyzer message-for-message on all three edited
files, 0 new; the script reaches the same verdicts at the same t_f
(18.6039 d catalog pair, 17.7976 d anchor) and refuses (0, 0.1) with the new
message.

### Addendum: section 5 the same way (2026-09-11)

Mike, in a comment left in the file: "NEW SECTION - BUT MAKE IT FEWER LINES
OF CODE". Those 22 lines did four jobs -- gate the solve, fly once, validate,
derive and report -- and three of the four were duplicated elsewhere.

| unit | owns | consumers |
|---|---|---|
| `costate_common/fly_transfer` | fly the costates ONCE, validate through `validate_flight`, attach the flown miss, mass, propellant and Delta-V | the study script now; `certify_root`, `audit_phase_catalog` and the witness flight in `verify_with_pumpkyn` spell the same thing out |
| `DRO_tulip/indirect/print_transfer_summary` | the three numbers both consumers print | the study script AND `run_dro_tulip`, which printed them in a different format (Mike: "should also serve the front door") |

**Delta-V has two homes and they are now checked against each other.**
`catalog_schema`'s `deltav_from_mf` derivation is the catalog's rocket
equation; `fly_transfer` needs the same formula but takes physics, not a
catalog, so routing it through the schema would mean faking a catalog
struct. It writes the formula and the test asserts the two agree on the real
catalog (0.7484813765 vs 0.7484813765 km/s).

**What stayed in the script, deliberately:** both asserts. They are the
SCRIPT's policy -- nothing below is meaningful from a best iterate or an
inadmissible flight -- and a library that threw on admissibility could not
serve a certifier that must return a named reason instead.

Section 5: 36 lines to 21, of which the code is five statements; the rest is
the banner and the comments that say why the gates are there.

**A test that could not fail, caught by writing down what it was for.** The
first "inadmissible flight" check scaled the costates by 3 -- but the
min-time direction is invariant under a positive scaling, so the flight
stayed admissible and the check passed for the wrong reason. Flying past
mass depletion is the real refusal, and the validator names it.

### Addendum: section 7 runs the production instrument, and V2's "second opinion" was a copy (2026-09-11)

Mike asked for the four necessary-condition blocks to become named
functions, then asked the better question: "if we have
`pmp_pointwise_checks` in the library, why are we not using it?"

We were -- in `certify_root` on every certified entry, and in the script at
V2. The script ALSO kept its own copy of the same four computations, and V2
compared them. **V2 measured 0.0e+00.** Not "agrees to 1e-14": exactly zero,
because a copy of the same arithmetic on the same samples is not an
independent implementation. The duplication was buying readable math, not
cross-validation.

So section 7 now calls the instrument the catalogs are certified with, and
prints its numbers. Every value is unchanged (N1 1.41e-09, N2 3.32e-08,
N3 0.0000 km, N4 7.76e-10, N5 3.75e-10 over 120 samples, N6 7.47e-15), which
is the point: the script did not lose a check, it stopped keeping a mirror.

Three smaller things fell out:

- **N3 needed no function at all.** `fly_transfer` measured the flown
  arrival in section 5; the script now reads `flight.flyKm` instead of
  recomputing the norm. One fewer copy of that formula.
- **V2 is relabelled for what it is: a WIRING check.** Both instruments fly
  the same trajectory with the same settings, so it agrees to 0.0e+00 and
  cannot see an implementation divergence. What it can catch is the script
  handing one of them the wrong flight or the wrong thrust, which is the
  mistake that actually happens when a section is edited.
- **The generic question, answered precisely.** `pmp_pointwise_checks` is
  already generic across PROBLEMS (the field is injectable, and it runs on
  every campaign). It is not generic across OBJECTIVES: min-time's running
  cost of 1, the free-time H == 0, and the all-burn control recovery are
  baked in, so the min-energy and min-fuel entries have no pointwise check
  at all. Injecting the running cost and control law is recorded in
  `costate_common/TODO.md` as new capability.

The script is 63 lines shorter on this pass (429 total).

### Addendum: the residual depends on HOW the arc is propagated (2026-09-11)

Mike, on the study script: "why are we calling ms_tfmin again when we call
it on line 191? Can't we use the output?" We can, and chasing the difference
turned up something that qualifies every residual this repository quotes.

Section 7's N1 re-assembled the boundary-value residual and evaluated it at
the returned point. That reported **1.41e-09** while the solve itself
reported **7.84e-12** for the same solution. Not a different point, and not
the segment grid (checked: identical to 0). The difference is the
PROPAGATION MODE:

| the SAME converged point, residual evaluated | \|R\|_inf |
|---|---|
| with its Jacobian requested (210-state propagation, finer steps) | 7.844e-12 |
| with plain state propagation (14 states) | 1.414e-09 |

`ms_bvp` iterates with the Jacobian, so `info.normR` -- and therefore every
`tolR` gate in every campaign, including `certify_root`'s 3e-11 -- is the
first number. A recipient who re-evaluates a stored solution the cheap way
will see the second and find it 180x worse than advertised.

Nothing here is wrong, and no verdict moves: 1.41e-09 still passes the
script's 1e-8. But **"converged to 3e-11" means "in the Jacobian mode"**,
and that qualification was nowhere on file. `costate_common/TODO.md` carries
the follow-ups: say it in the methodology document, decide whether catalogs
should also carry the state-only number a recipient will reproduce, and
check whether the gap widens on the longer arcs at deep thrust.

N1 now reports `it.normR` -- the solve's own residual, at the point it
returned -- which is one line instead of three and is the number the gates
are actually set against.

### Addendum: section 8 reads the gates, and the ID letters say what they mean (2026-09-11)

Section 8 recomputed two of the four hypotheses it reports. `min|lam_v|`
(S1) and `min Q` (S2) were taken from the script's own flight while
`mintime_hypothesis_gates` -- the instrument that JUDGES them, and the one
`gates_catalog_pass` ran over all 18,360 catalog entries -- returns both,
with the times at which they occur. One call now serves S1, S2, S3 and V1.

Values unchanged (min|lam_v| 3.2359 at t/t_f 0.117, min Q 3.8248, dim S 1,
H6 margin 9.0x); S2 gained the location of its minimum, t/t_f 0.984, which
the gates were already computing and the script was throwing away.

**V2 is now a cross-check between two INSTRUMENTS rather than between a
script and a library**: `pmp_pointwise_checks` reads the flight this script
flew, `mintime_hypothesis_gates` flies its own from z8. It still agrees to
0.0e+00, because both propagate the same z8 deterministically, so it is
still labelled what it is -- a wiring check.

**And the ID letters are now spelled out in the header**, because Mike had
to ask what they stood for: N NECESSARY, S SUFFICIENCY (the theorem's
HYPOTHESES, not its conclusions), V VALIDITY (V1 the conjugate
instrument's own precondition, V2 the wiring), X CROSS-CHECK. A legend that
lists IDs without expanding them is not a legend.

## 47. Two outside math reviews of transfer_study: the conditions were right, the gates around them were not (2026-09-11)

GPT-6 Astra (xhigh) and Gemini 3.1 Pro reviewed the mathematics of the
study script's N/S/V/X conditions
(`reviews/transfer_study_math_{astra,gemini}_2026-09-11.md`; adjudication
in `reviews/transfer_study_math_adjudicated_2026-09-11.md`). Astra confirmed
every formula -- Hamiltonian, mass-costate sign, switching function,
spherical Legendre, the H6 identity and the free-time quotient's rank
interpretation -- and then found the script ENFORCING less than it printed.
Gemini graded everything CORRECT and its one code finding (an "undefined"
function that exists) was a bundle omission. Where they disagreed, Astra was
right each time it could be checked.

**Confirmed and fixed, in the shared instruments:**

- **V1 accepted equality and ignored the clearance.** The script gated on
  `h6Margin >= 1` and never read `h6Ok`; `certify_root` and
  `report_optimality` had the same `>=`. All three now take the helper's own
  verdict (strict margin AND clearance above the Hamiltonian residual). No
  entry moves (margins 9-33x), which is luck, not design -- the same
  "computed, not enforced" defect section 41 fixed in the certifier had
  been reintroduced in the script.
- **N6 checked the throttle on three of the four rows it enters.** The
  applied control is now recovered on the acceleration rows AND the mass
  row (`u_mass = -c F_m / T`), the gap is against the FULL control minimum
  (direction + throttle), and its signed minimum is kept -- the old
  accumulator started at zero and took a max, so an over-unit thrust along
  the minimiser produced a perfect zero. Two mutation tests: half mass
  flow (direction gap 6.7e-16 stays clean, mass-row throttle reads 0.500),
  1.001x thrust (signed gap -2.9e-3, was clipped).
- **S4 could PASS without covering the final segment**, and the script
  printed the uncovered interval without gating on it. A free-time test
  without y(t_f) is now UNDETERMINED with a reason; the script requires
  `.covered`. A sign change on the last bracket between RESOLVED nonzero
  samples is an interior root (FAIL), not ENDPOINT; ENDPOINT remains for an
  unresolvable final sign. The two exact identities J p(0) = 0 and
  p(t)'J(t) = 0 are measured: 3.8e-15 and 8.0e-14 on the anchor.
- **N2 and N5 tested pumpkyn's field against itself.** The gates now
  compare pumpkyn's state rows and adjoint rows, row by row along the arc,
  against the hand-written CR3BP + thrust field and its CasADi Jacobian:
  1.8e-16 and 9.3e-15. That is the independent physics check the X1
  cross-check was credited with and does not provide (both solvers
  propagate pumpkyn).
- **S3 printed its self-consistency residuals and required none of them.**
  The lift residual, the Hamiltonian residual and the Eckart-Young margin
  (`lift_margin`, 2545x on the anchor) are now all required. The rank
  threshold alone forces at least one small singular value
  (sigma_min <= nullResid by construction); the margin is what makes
  "exactly one" a measurement. The lift residual's floor is 2.2e-6 at both
  relTol 1e-12 and 1e-10 -- it is the pchip interpolation of the flown arc,
  not the integration tolerance -- so the tolerance is lift_margin's own
  1e-4, handed through as one value.
- **The dense spectrum scan is wired into S4** (192 samples, two-level
  refinement, 0 candidates on the anchor), the endpoint interpolant is
  compared against a propagation to the phase (2.2e-8 ND = 8.6 m on the
  tulip), and `validate_flight` checks costates, pointwise mass law,
  positive mass everywhere and a positive EXPECTED final mass.

**Wording corrected:** N1 is "the shooting equations are satisfied", not
"the first variation vanishes" (on a boundary control the first-order
control condition is the minimum principle, N6); N2/N4 are re-evaluations,
not independent evidence; S2 follows from N4 + N5 + S1 on an exact lift;
X1 is a second shooting implementation on shared physics; the verdict names
the endpoint conditions (fixed r, v at both ends, m(0) = 1, free terminal
mass and time, phases fixed inputs) and calls itself numerical evidence,
not a certificate.

**Reproduction:** the script runs end to end with every line PASS and the
same numbers as before the review (t_f 17.7976 d, 0.7485 km/s, min|lam_v|
3.2359, min Q 3.8248, H6 9.0x). New test `test_conj_coverage`; ten
existing tests re-run green.

**Open (theory, for `doc/mintime_second_order_audit.tex`):** the subarc
normality argument by analytic continuation of the strict all-burn
extremal; the explicit reduction of the free-mass second variation to the
six-state nonautonomous problem; a short-time sign expansion to close the
interval before the first junction; between-sample bounds for S1, S2 and
the clearances. Listed in `costate_common/TODO.md`.

## 48. Astra's second pass: the dense scan could clear a true conjugate point, and the certifier still enforced less than the study (2026-09-11)

The corrected script and its sixteen source files went back to GPT-6 Astra
at xhigh (461 s, $1.79; `reviews/transfer_study_math_astra2_2026-09-11.md`,
adjudication in `reviews/transfer_study_math_astra2_adjudicated_2026-09-11.md`).
It confirmed section 47's corrections -- the full-gap formula, the mass-row
throttle, the frozen-control adjoint comparison including the mass row, and
both kernel identities with the mass costate present -- and found three
things still wrong. All three verified against the code; all fixed in the
order Astra ranked them.

**1. The dense conjugate scan was not a zero-exclusion test.** The S4 gate
read coarse sign changes, zeros and multiplicity but never
`nInteriorCand`, so a candidate classified "near-miss" passed; and the
near-miss classification was a PLATEAU under nested 4x/16x refinement,
which a true double root within h/32 of a coarse node reproduces exactly
(the same node stays nearest at every level, both ratios read one).
Endpoint clusters were never refined at all. `conj_spectrum` now refines on
SHIFTED grids, brackets the finest-grid argmin one step either side and
LOCATES the minimum of sigma_6(t) by golden section; the located minimum is
judged against a numerical floor (1e-7 of the median, a policy value):
at the floor or with a sign change it is a ZERO, a hundred times above it a
NEAR-MISS with a measured positive minimum, and in between UNRESOLVED,
which blocks. Class is by contiguity, not by band: only a cluster touching
the first sample is the (uncovered) start transient; endpoint clusters are
refined through t_f. Multiplicity is the number of singular values at the
floor at a located zero. On the anchor the endpoint dip resolves to a
positive minimum of 8.0e-5 of the median at t/t_f = 0.9994 and clears; the
refuted sheet entry still reads ZERO. One second per scan.

**And the junction sign test could PASS on an untrusted sign.** The
resolved-last-bracket rule of section 47 acted only when a sign change was
seen; a same-sign tiny final determinant passed. Any live sample whose
sigma ratio is at or below `resolvedTol` now makes the verdict
UNDETERMINED unless an interior root was already found.

**2. `certify_root` gated `dirGap` and nothing else new.** The study
required the full gap, both throttles, X2, the lift residuals, the
Eckart-Young margin and the dense scan; the production certifier required
none of them and returned "certified" -- the computed-not-enforced defect
of sections 36, 41 and 47, one more time, in the file whose header says
that is the failure mode it exists to prevent. (The library chain does run
`second_order_pass` at the sweep stage, which Astra could not see, so the
shipped sheet did carry the dense scan and the margin; per-entry
certification did not.) Gates 2b, 5 and 7 now match the study, and
`test_certify_enforcement` proves the caller reads them: five tolerance
squeezes and four field injections each refuse by name (half mass flow
refuses with throttle error 0.500 stored; a wrong field handed to X2 alone
refuses on X2 at 3.7e-1; a clear factor of 1e12 turns the anchor's endpoint
near-miss UNRESOLVED and blocks), and the anchor certifies untouched with
every new number carried.

**3. One option value, two normalisations.** The script tested the lift
residual relative to |lam| and `lift_margin` tested it relative to
sigma_1 |lam|, both against "1e-4". The BACKWARD error
|C lam|/(|C||lam|) (`gates.nullResidRel`) is now the one normalisation in
the script, the certifier and `lift_margin` (liftTol 1e-6); on the anchor
it reads 5.3e-10 (sigma_1 = 4.1e3, |lam|-relative residual 2.2e-6). The endpoint
check multiplied a six-state norm by lStar and called it kilometres; it now
reports position and velocity separately -- and the number changed:
section 47's "2.2e-8 ND = 8.6 m" was the six-state norm, velocity-
dominated; the arrival POSITION error is 1.3e-10 ND (0.05 m) and the
velocity error 2.2e-5 m/s -- wraps the phase as `phase_state` does, and is
labelled an endpoint-consistency estimate against a numerical reference.
The `fieldGap` is evaluated too -- the field's own gap differs from the
reconstructed one by (T/c) lam_m (|b| - u_mass) -- and the gated number is
the worse of the two (half mass flow: 5.6e-2 against 1.9e-14).

**Corrected claims.** `lift_margin`'s header said sigma_6 > |dC| gives
rank six EXACTLY; Eckart-Young gives rank >= 6 (sigma_6 is the distance to
rank at most five), and exactly six needs the exhibited lift. The
two-tolerance difference is a sensitivity estimate of one error component
(the adjoint integration), not a total bound: flight, interpolants,
endpoints and assembly cancel in it. `lift_space_dim`'s header said its cap
made dim S "never > 1"; it does not. The envelope comment now says "no
derivatives of the control law", not "does not depend on it".

**The STM generator, checked independently** (`test_stm_variational`): all
14 columns of one segment's STM against central differences to 8.6e-8
(the lam_v columns, which carry d alpha*/d lam_v, to 5.2e-9), and
PHI' J PHI = J to 1.9e-10. A lesson came with it: symplecticity does NOT
discriminate the missing control-law derivative -- a frozen-control
generator is the linearisation of a fixed-alpha Hamiltonian and is
symplectic to the same 1e-11. The finite-difference columns are the
discriminating check; the frozen STM differs from the true one in the
lam_v columns by O(1).

**Verdict wording** now names the structural conditions the diagnostics do
not establish (subarc normality, the free-mass/free-time second-variation
reduction, existence of an exact extremal near the numerical one), in the
study and in `report_optimality`.

**Consequence for the 70 mN library:** its second-order sheet was swept
with the old plateau classification; the 44 near-miss cells are old-kind
evidence and `conj_unresolved` is NaN for them. A RE-SWEEP with the
resolved scan is required before the ship decision (TODO).

## 49. Astra's third pass, on the new code only: a located minimum is not a cluster, and the synthetic tests found the blind spot the reviewer described (2026-09-11)

Pass 3 (`reviews/transfer_study_math_astra3_2026-09-11.md`, adjudication
`..._astra3_adjudicated_...`) was scoped to what pass 2 had introduced, and
its bottom line was that the classifier was not yet fit to re-sweep the
library with. Items 1-7 of the adjudication were applied in its order;
item 8 -- a floor derived from a measured matrix error -- stays open.

**The resolution is now a pure function, and the tests drive it with
matrices whose rank loss is known.** `conj_resolve(Mfun, tGrid, opts)`
takes any matrix-valued function of time; `conj_spectrum` hands it the
propagated conjugate block. Everything it evaluates is KEPT: shifted grids
at 4x and 16x, then EVERY local minimum of the assembled record bracketed
and located by golden section with an a-posteriori unimodality check; a
bracket of opposite TRUSTED signs anywhere in the record is an established
root; the smallest value ever seen decides, and a later search returning a
larger value cannot upgrade it. The start cluster is refined past its
first local maximum -- the increasing prefix is the structural transient
and is the reported uncovered interval -- so a zero merged into the
transient shows as a dip after some growth. The last sample is a candidate
source and is evaluated. A vanishing raw column norm is a candidate in its
own right, because column normalisation hides that rank loss. A scan with
no trusted sample outside every cluster is not testable and not clear.

**The synthetic tests found a blind spot the review only gestured at.** A
V-shaped zero between two coarse samples reads far above the 1e-3 dip
threshold at both of them: with slope ~median/(0.1 t_f) the nearest sample
sits at 0.03 of the median. Without a sign change -- a corank-two zero, a
quadratic touch -- such a zero was INVISIBLE to candidate detection. The
fixture at 2.5 + h/3 failed; the on-node case had passed only because the
node happened to sit on the root. Local minima of the coarse spectrum
below 0.2 of the median are now candidates at any depth. On the anchor this
finds one more candidate, a shallow interior dip at t/t_f 0.776 (5e-3 of
the median), resolved as a near-miss and cleared; the scan takes 0.8 s for
511 evaluations. `test_conj_resolve`: 26 checks, every one asserting the
GATE -- transverse zeros at six grid phases including on a node and h/32
from one, corank-two zeros without a sign change on and off a node,
quadratic touches, two wells around a zero and the same wells cleared, a
zero at t_f and one fifth of a step before it, a zero merged into the
start (caught; uncovered stops at 0.148 before it at 0.2), a vanishing
column, a near-miss cleared at 1e-4, an unresolved minimum at 3e-6, and a
scan that never becomes full rank.

**The junction test classifies before it counts.** Every live sample is
trusted-positive, trusted-negative or unresolved first; a root is a
bracket of opposite trusted signs with unresolved samples between them
skipped; equal trusted signs around an unresolved run establish nothing
(the old "touch" is gone); an unresolved sample makes the verdict
UNDETERMINED unless a trusted bracket already refutes. ENDPOINT is no
longer a verdict. The two quotient identities are gates: a violation means
the block is not interpretable -- UNDETERMINED, never a refutation.

**The certifier validates the domain before it compares.** In MATLAB
`if ~(x <= tol)` is skipped by an empty x and by a vector with one passing
element, `-Inf` passes an upper bound, and `max` drops a NaN -- so a NaN
adjoint evaluation upstream read as a perfect zero error. Every gated
field is now a real finite NON-NEGATIVE scalar before any aggregation; the
pointwise checks and the gates poison a residual to NaN on any non-finite
evaluation and report the count; `lift_margin` validates its inputs before
the SVD; the dense counts must be finite non-negative integers and `.clear`
consistent with them; `mintime_prop_seg` asserts the propagation arrived.
`conjSpectrum = false` is a DIAGNOSTIC result (`ok` false, `okDiagnostic`
true), not a certificate. A TEST SEAM (`opts.override`, warns loudly) lets
each gate be the first failing one: 17 seam mutations -- empty, vector,
-Inf, NaN, negative, wrong type, inconsistent flags, and each new
tolerance alone -- refuse by name.

**Still open, recorded in TODO:** the floor as a measured matrix error
(the form: safety factor x error in the scaled matrix at the candidate
time, with the error measured by re-propagating from t = 0 at another
setting -- both current refinements inherit the same stored prefix); the
sign-trust threshold from the same measurement; the generator-block
discriminator A(4:6,11:13) = -(T/(m rho))(I - alpha alpha'); blockwise lift
residuals beside the global backward error, which at sigma_1 = 4e3 admits
too loose a residual on its own.

## 50. Seven hundred and fifty micrometres: a catalog is only as reproducible as its endpoint rule (2026-09-12)

The first live re-sweep of the 70 mN library with the resolved conjugate
classifier stopped at the audit, 56 minutes in, and never reached the sweep.
One entry of 115 came back BAD and `build_70mN_library` asserted on it. Two
separate problems, and both were worth the run.

### The gate did its job, and it cost a column

At arrival phase 0.8254 the sheet stage certified nothing. Its two candidates:

| candidate | t_f | verdict |
|---|---|---|
| crossing 9 | 25.503 d | dense scan not clear: **1 UNRESOLVED** (7 near-miss cleared) |
| crossing 10 | 25.507 d | flight inadmissible: Moon approach 1160 km < 1900 km |

The second is a pre-existing admissibility rule. The first is the new
resolver (FINDINGS 49) refusing to clear a candidate whose located minimum
falls inside the numerical floor band. So the rebuilt catalog carries 114
entries where the shipped one carries 115. **This is not a refutation**: an
UNRESOLVED verdict says the instrument cannot tell at the current floor, and
the floor is the policy value item 8 of the Astra-3 TODO is about.

### The audit failure was 0.75 mm

Entry (12,11) -- departure phase 11/12, arrival 0.9087, the longest arc in the
library at 26.4 days -- failed on witness disagreement, `|dz|` 7.4e-6 against
a 1e-6 tolerance. Its costates are **bitwise identical** to the shipped ones,
and the witness is deterministic (three runs, 7.398718e-06 every time).

The cause is the endpoint. `phase_state` (FINDINGS 44) replaced a privately
copied, silently-falling-back ordinary spline with the periodic cubic. The
two rules differ at exactly one grid phase:

| departure phase | periodic vs ordinary |
|---|---|
| 0 .. 0.8333 (columns 1-11) | 0, to machine zero |
| **0.9167 (column 12)** | **1.93e-12 ND = 0.75 mm** |

0.9167 is the phase nearest the seam on the coarse 105-sample DRO table, which
is where FINDINGS 44 predicted the ordinary spline would misbehave; the tulip
table has 1328 samples and never differs. Flying the stored costates proves
the mechanism exactly:

| departure endpoint used | flown miss |
|---|---|
| ordinary spline (what the rib was certified against) | **0.064220 km** = the shipped audit |
| periodic cubic (what the code gives today) | **2.124766 km** = the failed audit |

**0.75 mm becomes 2.1 km over 26.4 days, an amplification of 2.8 million**,
and the witness moves the costates by 7.4e-6 to absorb it. Four of the five
entries whose flown miss moved at all are in that one column. The re-certified
spine came back bitwise identical because sD = 0 is a KNOT of both
interpolants.

**The entries were never wrong.** Each still solves the problem it was
certified for. They solve a problem 0.75 mm away from the one the current code
poses. FINDINGS 44's bitwise-equivalence claim was true for the four phases it
captured -- none of which was a rib departure phase on the DRO table.

### What was done

**The data.** `repolish_endpoints` flies every stored rib point from the
CURRENT endpoints and re-polishes only those missing by more than a tolerance,
from the entry's own junction states, then re-certifies through the full gate
stack. Self-selecting, nothing hard-coded to a column: 104 points checked, 1
re-polished, |dz| 1.276e-08, miss 2.1248 -> 0.0644 km, 21 seconds, rib file
backed up before rewrite. A re-run on a clean library changes nothing.

**The chain.** A bad audit row no longer aborts it. The audit is diagnostic;
the SHIP gate is where badness belongs. Blockers are collected, the sweep and
pictures stages are individually fenced, the deliverable refuses by naming the
blockers, and the run ends `CHAIN CLEAN` or `CHAIN COMPLETE WITH n
BLOCKER(S)`. The batch verdict distinguishes a blocked library from a failed
chain. Losing an hour of sweep measurements to a diagnostic row is the defect
this fixes.

**Two guards.** `audit_phase_catalog` gains an ENDPOINT REPRODUCTION gate
(`tolFlyKm`, 1 km, against a worst legitimate 0.29 km): an entry that no
longer flies to its own endpoint is named as such, where it happens, instead
of surfacing as an unexplained witness disagreement. And `test_phase_state`
now pins the twelve departure and twelve arrival phases the library actually
uses, captured and verified bitwise, so a future change to the endpoint rule
fails in a test rather than in a catalog.

**The standing lesson.** A catalog entry is keyed by its PHASES, so it is only
as reproducible as the rule that turns a phase into a state. Any change to
that rule is a change to every stored entry, however small it looks: this one
was invisible in every direct measurement of the interpolant (the endpoint is
within 0.4 m of the orbit by propagation) and only appeared after 26 days of
amplification.

## 51. The re-sweep, read out: no zeros anywhere, and eight entries the instrument declines to clear AT t_f (2026-09-12)

The 70 mN library re-swept with the resolved classifier (FINDINGS 49), on the
repaired ribs (FINDINGS 50). `CHAIN CLEAN`, 114 entries audited 114 OK,
census complete, written back. Sweep wall time 3 h 40 m, about 2 minutes an
entry against the old classifier's seconds -- the cost of locating every
candidate instead of reading a plateau.

### The census

| measure | shipped (plateau rule) | rebuilt (resolved) |
|---|---|---|
| entries | 115 | 114 (column 10 lost both candidates, FINDINGS 50) |
| interior sign changes | 0 | **0** |
| located ZEROS | -- | **0** |
| multiplicity | 0 | **0** |
| cells with interior candidates | 42 (near-miss) | **109**, up to 7 each |
| cells with a cleared near-miss | 42 | **111** |
| cells UNRESOLVED | n/a | **8** |
| junction conjugate test | 115/115 PASS | 114/114 PASS |
| worst H6 / worst lift margin | 4.5x / 22x | 4.54x / 22x |

**The scan now finds far more candidates and clears nearly all of them.** 109
cells carry interior candidates against the old 42, because the coarse
LOCAL-MINIMUM rule added in FINDINGS 49 catches V-shaped dips that the depth
threshold alone missed -- the very blind spot the synthetic tests exposed.
Every one of them resolves to a located positive minimum. No entry in the
library produces a zero, an interior sign change, or a multiplicity event.

### The eight

All eight unresolved cells are the same phenomenon, and it is not scattered:

| cell | located minimum / median | t/t_f |
|---|---|---|
| (3,11) | 9.3e-07 | 1.0000 |
| (4,11) | 1.5e-06 | 1.0000 |
| (5,11) | 2.6e-06 | 1.0000 |
| (6,11) | 3.4e-06 | 1.0000 |
| (7,11) | 4.5e-06 | 1.0000 |
| (8,11) | 6.6e-06 | 1.0000 |
| (9,11) | 6.4e-06 | 1.0000 |
| (10,11) | 8.0e-06 | 1.0000 |

Every one is in arrival column 11 (sA 0.9087), the longest arcs in the library
at 26.4 to 27.7 days; every one is an ENDPOINT candidate whose minimum sits
exactly at t_f; and every one lands inside the floor band -- above the
zeroFloor of 1e-7, below the clear threshold of 100x that. The values rise
monotonically with departure index, which is the signature of a systematic
effect rather than eight independent near-conjugacies.

**This is the graded-endpoint collapse the instrument has always known about.**
`conj_spectrum`'s own header records the measurement: at t_f the hyperbolic
flow grades the whole spectrum down together, and a certified entry and a
refuted one measured sigma_min 1.37e-7 and 1.39e-7 there -- the endpoint value
does not discriminate. The old code therefore classified endpoint clusters and
never refined them. Astra's third review called that exclusion unsound (a real
even-order zero is not harmless for lying in the last segment), so they are
now refined -- and refinement lands them in the band where the floor cannot
separate grading from a zero. UNRESOLVED is the honest verdict: **not a
refutation, and not a pass.** The junction sign test still says PASS on all
114, and no interior structure appears anywhere.

### What this makes urgent

Open item 8 of the Astra-3 list -- a floor derived from a MEASURED error in
the scaled matrix rather than a policy value -- was a tidiness item this
morning. It is now the difference between a library that certifies 114 of 114
and one that certifies 106. The measurement to make is the error in the
projected matrix at t_f, from re-propagating from t = 0 at a tighter setting
or with a different integrator: if it is ~1e-5 of the median, these eight are
unresolvable in principle at this precision and must be reported so; if it is
~1e-9, a calibrated floor clears them and the library is uniform.

A second, cheaper discriminator is worth testing alongside: a located minimum
sitting exactly AT the interval endpoint with a monotone approach is what
grading looks like, whereas a conjugate point has its minimum strictly inside.
That distinction is suggestive, not decisive -- a conjugate point at t_f is
possible, and is what the junction test's retired ENDPOINT verdict used to
flag -- so it cannot be used to clear an entry on its own.

Nothing ships until one of those closes. The catalog carries the verdicts;
the deliverable stage was off.

## 52. The campaign chain reviewed: four defects that matched four paid-for incidents, and what a one-host queue actually needs (2026-09-13)

`run_costate_library.m` and the five campaign primitives under it
(`work_queue`, `campaign_worker`, `campaign_heartbeat`, `campaign_status`,
`run_campaign_workers.sh`) went to GPT-6 Astra at xhigh (59 KB bundle,
589 s, $1.47; `reviews/campaign_chain_astra_2026-09-13.md`). Verdict: "not
fit for an unattended multi-day campaign". Forty findings. Each was checked
against the code before anything was changed; the ones below were CONFIRMED
and fixed, the rest are adjudicated at the end.

**Two things fixed first, from Mike's reading of the entry script, before the
review landed.** (1) Where the phases are chosen was not visible: the arrival
grid came from `nA` in the sheet stage and the departure grid from `nD` deep
in the rib stage. Section 0 now takes `.sD`/`.sA` vectors (or derives them
from `.nD`/`.nA`/`.sA0`), refuses a non-uniform grid by name, and prints
both. (2) Which DRO and which tulip was not visible either: they came from
`arclength_arrival`'s own defaults. Section 0 now declares `tauDRO`,
`NpTulip`, `pmTulip` and the engine, prints them with the derived periods,
passes them to the sheet builder, and REFUSES to continue if the sheet on
disk was built for a different problem. Making the tulip knobs honest
exposed a latent defect: the setup accepted `NpTulip` but hardcoded the
period as `5*2*pi/6` and the branch as `-1`, so an 8-petal request would
have propagated an 8-petal orbit for the 7-petal period. The period is now
derived, `2*pi*(Np-2)/(Np-1)`, and the branch is read; the packager
(`sheet_to_catalog_file`) carried the same literal and now labels the
catalog with the sheet's own period. Bitwise identical at Np = 7 (the
rebuilt setup matches the stamped identity of the existing sheet on all 14
fields).

**The four confirmed headline defects.**

1. *The generated worker never beat.* `unitFcn = @(j, beat) build_ribs(...)`
   took `beat` and dropped it, so a claim was refreshed only between
   columns. Columns take 1.5 to 9 hours; the stale lease was 30 minutes.
   Any worker finishing early would have "reclaimed" a live column -- the
   exact mechanism that cost columns 13 and 16 nine hours each under the old
   launcher. Fix: `rib_from_crossing` takes `.progress` and calls it after
   EVERY solve (2-3 min); `build_ribs` forwards it; the job passes the
   worker's heartbeat there.
2. *Re-running destroyed live ownership and reset the retry budget.*
   `work_queue('init')` wiped every claim and every attempt file, and the
   script's own instruction was "re-run this function to package". So the
   documented workflow would have freed columns being walked and re-armed
   the livelock the attempt counter was built to stop. Fix: `init` opens an
   existing queue (claims and attempts kept, new units added) and creates
   only when none exists; `reset` is a separate explicit action that
   refuses while any claim is fresh.
3. *The watchdog measured lifetime, not inactivity.* A timer started at
   launch killed a healthy worker inside its fourth column regardless of
   progress. Fix: the launcher's watchdog reads the heartbeat AGE and kills
   only after WATCHDOG seconds of silence, which with per-solve beats is a
   hang and nothing else.
4. *Packaging followed launch with no barrier, and a file's existence was
   completion.* Fix: the entry script returns `out.state` in {pending,
   launched, blocked, packaged}; stage 4 runs only when every unit's
   artifact exists and none is claimed; rib files and the queue's own
   records are published atomically (write beside, move), so an interrupted
   save is not a finished unit.

**Also confirmed and fixed.** Ownership is now a TOKEN: beat and release
act only if the claim still carries the caller's token, so a worker
reclaimed while blocked in a solver cannot refresh, then delete, the new
owner's claim. Stale takeover is a rename (atomic; one reclaimer wins).
After winning a claim the done and attempt checks are made AGAIN. Attempt
records fail CLOSED (unparseable = blocked, unwritable = no claim). Status
gives ONE state per unit (a live final attempt is running, not retired) and
`nOpen`/`complete`/`finished` are stated rather than inferred from `nTodo`.
An empty or garbled heartbeat is UNKNOWN, never running; heartbeats are
written atomically. The worker's whole lifecycle is guarded (a fatal error
outside the unit try now writes a `fail` heartbeat), accounting happens once
at the work boundary and logging cannot turn a saved unit into a failed one,
and with `idleSec` the last live worker waits for held units instead of
leaving the tail unattended. The generated job derived `costate_common`
from the OUTPUT directory -- three `fileparts` of the sheet path gives
`DRO_tulip/costate_common`, which does not exist; it worked only because
`startup` had already put the real one on the path. Code roots now come from
the driver's own location; every path is absolute (MATLAB's `run` changes
directory to the job's folder); the shell command is quoted and MATLAB
literals escaped; `system`'s return code is checked; the launcher validates
its arguments, gives each launch an id so tags and logs are never reused,
and detects a worker that exits before its first beat. `test_work_queue` now
THROWS on failure and carries 28 checks, one per guarantee above (it also
caught my own bookkeeping error: a unit I thought untouched carried a prior
attempt, and the queue correctly retired it).

**Pushed back or deferred, with reasons.** A campaign manifest with fencing
generations and a coordinator-owned lease clock: this is one host, APFS,
five workers; the token plus the rename takeover plus the re-check close the
races that can actually occur here, and Astra's own text allows that "for
one host, kernel-held interprocess locks plus an exact-child supervisor may
be simpler". Separate-process barrier race tests: open; the claim test is
still sequential. Intra-column checkpoint/resume in the bisecting walker: a
real gap (a reclaimed column restarts from zero) and a bigger change; open.
Converting `build_70mN_library` from a script to a function: deferred; the
base-workspace handoff now saves and restores whatever was there. Tracking
calibration through the queue: deferred; its measurement is now persisted
so a later call cannot replace it with the 3600 s guess. `fmt_num` width
overflow: minor, open.

**State of the library.** 14 of 19 certified columns on disk; 13, 14, 16,
17, 18 in flight on the OLD explicit-range launcher (no queue), all five
alive at the time of writing. The new chain has been exercised end to end
with launch off against the live results directory and returns `pending`
with the launch command; it will be used for whatever the old workers leave
unfinished, and for every library after this one.

## 53. What the RUNNING 24x24 campaign needed to finish: a watchdog 27 minutes from firing, and a finalizer wired to the wrong library (2026-09-13)

Asked "are there fixes the currently running code needs to succeed", the
answer was yes, twice, and neither was in the rib workers' own code.

**1. The old launcher's lifetime watchdog.** `run_fine_ribs_range.sh` arms
`( sleep $SEC; kill $MPID )` per worker; these were launched with SEC =
21600 (6 h). Workers I and J started at 07:10:36, so both would have been
killed at 13:10 -- found at 12:43. I was half-way down column 16 after 5.5 h.
This is incident 2 of the discipline doc again. The four watchdog subshells
were killed with SIGKILL BEFORE their sleeps (killing the sleep first would
have let each subshell run straight on to its kill line); all four workers
confirmed alive. Hangs are covered by log age (the session monitor wakes on
20 min of silence).

**2. Static-range duplicates.** Worker H held [13 14]; worker K, launched
later, held [14]. H finished 13 at 12:25 and started 14, which K had walked
for 40 min: the skip-if-file-exists check runs only at column START. H was
stopped (column 13 on disk, nothing lost). Worker I holds [16 17] while L
walks 17, so a guard stops I as soon as it announces column 17, which it
prints only after column 16 is saved.

**3. The finalizer was the 12 x 12 library's.** `run_costate_library`
stage 4 handed `build_70mN_library` only outDir and switches, and the
chain's own parameter block pinned nA = nD = 12, the default sheet name
(not `arrival_sheet_70mN_nA24.mat`) and the rib glob `results/arrival_rib*`
(the 12 x 12 ribs). It would have asserted "sheet missing" hours after the
last column, or, given a matching name, packaged the wrong ribs at nD = 12.
Astra's pass-2 review flagged the same (K). Fixes: the chain accepts
`chainOverrides.grid`, `.sheetFile`, `.ribFiles` (byte-identical defaults;
the loaded sheet must have nA levels); the entry script passes the exact
grid, sheet and per-column rib files; it runs the chain SCRIPT inside a
local function, so its `clearvars` hits that workspace and not the shared
base (it wiped a test harness variable in the shared session today); the
working directory is restored by the caller; `packaged` requires a catalog
written by THIS call (mtime); the completion barrier (every column LOADS as
a rib, no live queue claim) now applies whether or not the rib stage runs;
and section 0 refuses phase vectors the builders cannot honour (anything
but the full 1/n lattice) and checks the sheet's arrival phases and
departure origin against the request.

**Verified on the finished columns, not assumed.** 15 of 19 columns
packaged through the entry script into a scratch folder: 329 entries (19
spine + 310 rib points) on the 576-cell grid, torus picture sensible, base
workspace and cwd intact. A detached sample job then ran the rest of the
chain on that catalog: audit 10/10 OK (including sD 0.9167 and 0.9583, the
phases nearest the seam where the 12 x 12 audit failed by 2 km: here 0.003
and 0.014 km), sweep 3 entries in 75 s on 2 workers, movie 8 frames in 15 s.
Scaled: audit ~1-2 h, sweep ~2-3 h, full movie ~20 min for ~450 entries.

**Armed.** `batch/fine_library_autochain.sh` waits for 19 columns on disk
AND no old rib worker alive (they hold no queue claims -- Astra pass 2, N),
then runs `batch/fine_library_finish_job.m`: package, audit, sweep, and the
full movie `results_fine/sweep_full_library_24x24.{mp4,gif}` at slow = 2,
with a verdict in `results_fine/FINISH_VERDICT.txt`. If the workers drain
with columns missing it writes BLOCKED and packages nothing.

Astra's pass 2 on the chain (`reviews/campaign_chain_astra_pass2_2026-09-13.md`,
206 KB, 393 s, $1.66) is not yet adjudicated beyond the items above; its
verdict is still "not fit for an unattended multi-day campaign", centred on
the queue's ownership protocol, which the running campaign does not use.

## 54. Astra's second pass on the chain: ownership by heartbeat age was an ABA race, so ownership is now a lock the process holds (2026-09-13)

`reviews/campaign_chain_astra_pass2_2026-09-13.md` (206 KB bundle with the
pass-1 review and my section-52 adjudication included, 393 s, $1.66).
Verdict again "not fit for an unattended multi-day campaign", and the
centre of it was right: the first queue transferred a claim whose beat was
30 minutes old, by renaming the stale directory. Two reclaimers could both
succeed (A renames B's claim and creates a fresh one; C, already past its
staleness check, renames A's fresh claim), and a slow-but-alive owner
would come back and refresh, then delete, its replacement's claim. The
tokens I had added did not close it: a check followed by an action is the
same stale observation. Astra's smallest closure was a kernel-held per-unit
lock retained through publication, with no age-based takeover at all; a
hung owner is killed by a supervisor and the kernel frees the lock.

**Primitives measured before anything was built on them.** (1) A java.nio
FileChannel lock from MATLAB is refused across processes while held,
acquired the moment the holder releases, and acquired after the holder is
killed with SIGKILL -- `unit_lock`. (2) MATLAB's `movefile` onto an
existing FILE does rename (the destination takes the source's inode; my
first reading said "copy" and was a harness bug: the captured text held
MATLAB's own status echo), but `movefile` onto an existing DIRECTORY moves
the source INSIDE it, which is how a takeover could have nested one claim
in another. java.nio `Files.move` with ATOMIC_MOVE is documented to rename
or throw, and throws "Is a directory" on that case -- `publish_atomic`.

**The model now.** A unit is owned by a lock its worker process holds
until it releases or dies; nobody can take it from a live owner. The
attempt is counted under the lock, before the work. The unit writes to an
attempt-specific temporary name the queue hands it; the WORKER validates
that file (`rib_validate` for ribs: loads, has certified points, carries a
problem identity; coverage reported) and moves it onto the output in one
rename while still holding the lock. "Returned normally" is not success. A
worker killed on its last attempt leaves the unit RETIRED, not the half-
open state Astra found (status said open, claim said skip). Heartbeat age
is an alarm; the launcher's per-child supervisor kills a worker silent for
`hangSec` = 2700 s (three times the solver's 900 s wall cap, the longest
silence a healthy walk can have), and has a startup deadline independent
of any heartbeat. The unowned calibration stage is gone. The queue's unit
set is fixed at creation and `open` is read-only. An immutable campaign
manifest (orbits, engine, phases, policy, code revision) is written on the
first call and checked on every later one. The entry script refuses to
launch while any foreign (old-launcher) worker is alive, refuses phase
vectors the builders cannot honour, checks the sheet's phases and origin
against the request, and names short columns as coverage blockers. Idle
workers wait for held units (polling every 15 s) so the tail is never
unattended. READY is a persistent file the worker writes after opening
the queue.

**Verified with processes, as Astra required.** `test_campaign_processes`
launches three real MATLAB workers through the real launcher against a
six-unit queue (5 s units; unit 3 takes 40 s and its owner is killed -9
mid-unit; unit 5 always throws). All eleven checks pass: every worker
READY; units 1, 2, 4, 6 computed EXACTLY ONCE; unit 3 taken over after the
kill and finished (two starts); unit 5 retired after exactly three
attempts with a `.failed` file and no stray `.part`; the queue finished
with no lock held; the killed worker reported running-stale, not done;
the survivors exited clean. Three runs failed first, each on the TEST or
the LAUNCHER, never on the queue: the cleanup deleted the evidence; the
waiter matched a job path that now travels in the environment, not argv
(rule 4 again); and READY was a heartbeat record the worker overwrote
within milliseconds, so a 2 s poll never saw it and reported workers that
had finished a whole campaign as "exited before ready". The 34-check
sequential suite passes; the entry script dry-runs both paths and refuses
a mismatched orbit, grid or manifest.

**Adjudicated and not done, with reasons.** Intra-column checkpoint/resume
(a reclaimed column restarts from zero; Astra: conditionally acceptable
once ownership is safe, and it is) -- open. A campaign-wide supervisor
that replaces lost capacity -- the idle-wait covers the tail with the
workers already launched; relaunching is a call to the entry script. The
hardcoded pumpkynPie bootstrap in the generated job is this machine's
environment. Separate-process race tests for simultaneous fresh claimers:
the lock makes that a kernel property rather than a protocol, and the
process test exercises it with three workers claiming from one queue.

**Also surfaced by the new barrier.** Six of the fifteen finished 24x24
columns are shorter than 23 points -- the walker stalled at a dense
conjugate scan (cols 6, 7), a polish that did not converge (12, 13, 15),
and, at column 21, a transversality margin of 1.45e-6 against a 1e-6
tolerance after only two points. The old barrier counted those files as
full columns. They are real, terminal, shorter columns; the ship decision
needs to see them.

## 55. Astra's third pass: the lock architecture is right, its implementation is not yet, and two blockers reproduce (2026-09-13)

`reviews/campaign_chain_astra_pass3_2026-09-13.md` (343 KB bundle incl.
certify_root and run_capped, 531 s, $2.42). Verdict: "not yet fit" -- but
"the rewrite does remove the pass-2 architecture's central ABA takeover
problem ... a one-host lifetime-lock model is appropriate. The remaining
work is no longer 'invent a sound ownership architecture'. It is
'implement the chosen architecture faithfully, correct the false timing
premise, and test its failure boundaries'."

Its four blockers, checked before being believed:

1. **A released claim keeps its authority -- REPRODUCED.** `release` closes
   the Java lock but cannot mutate the caller's struct, so `c.lock.held`
   stays true. In one session: A claims and releases; B claims the unit; A
   publishes with its old struct -> `ok = 1` and the output holds A's
   result; A's stale release deletes B's owner record. My sequential test
   passed only because it set `fake.lock.held = false` by hand.
2. **The registry can drop a live kernel lock -- REPRODUCED across
   processes.** Hold a lock; `clear unit_lock` (wipes the persistent map);
   probe the same file from the same process (opens and closes a second
   channel). A second MATLAB process then ACQUIRES the lock while the first
   still reports `held = 1`. POSIX fcntl semantics, exactly as predicted.
   Also: keys are raw path strings (aliases), and a repeated stale release
   can remove a newer registry entry.
3. **The 2700 s hang deadline is wrong -- CONFIRMED from the code.**
   `progress()` fires only after a whole `certify_root`, which runs seven
   separately capped stages: polish 900 + flight 300 + witness 300 +
   witness flight 300 + gates 900 + second gates 900 + dense conjugate scan
   900 = 4500 s, plus uncapped work, with no aggregate cap. "3 x the 900 s
   solve cap" multiplied the wrong number. Also `run_capped` requests
   `cancel(fut)` and does not verify the pool worker stopped.
4. **Reset is still check-then-act** -- by reading: it probes (acquires and
   RELEASES) then deletes attempt/owner records without holding the lock.

Also deterministic and small: the completion barrier's
`short(~good(ismember(cols, cols)))` indexes two arrays of different
domains and THROWS on the invalid-column path; `rib_validate` checks the
NAME of the identity variable, not its value or the column; the startup
deadline stops applying once any heartbeat exists; the finalizer's audit
stage is not fenced like the sweep; `package_phase_catalog` and friends
were not in the bundle.

**Effect on the live 24x24 run: none directly.** It runs on the OLD
launcher (no queue, no supervisor, watchdogs disarmed), and the finish job
packages through the entry script's barrier. The indexing bug would turn
an invalid final column into a thrown error instead of a named blocker --
nothing is packaged either way. An audit exception would stop the finish
job before the sweep and movie.

Not yet applied. Astra's five changes, in its order: lock ownership as a
lifecycle-controlled handle with a protected, identity-keyed registry and
exception-safe release, with reset holding the lock; stage-level progress
or a real aggregate deadline, plus verified pool cancellation; a campaign
controller lock serialising manifest, launch and finalization, with strict
unit-aware rib validation; an owning parent supervisor that reaps and
records exits; and a deterministic fault suite (stale handle, reset
boundary, commit/kill matrix, real supervisor kills, real pool lifecycle,
two controllers).

## 56. The pass-3 blockers closed: authority in a registry, beats at every stage, a parent supervisor, and the fault suite Astra asked for (2026-09-13)

Every blocker in section 55 is closed, and each closure is tested against
the sequence that reproduced it.

**Authority lives in the registry, not the struct.** `unit_lock` keeps a
process-wide, mlock'ed map keyed by file IDENTITY (device:inode) from the
lock file to the live Java objects and the current holder's token.
`holds`, `beat`, `publish` and `release` ask it; a released handle has no
authority whatever its copy of `held` says; release is idempotent and by
token, so a stale release cannot remove a newer holder's entry, and a lock
file deleted under a holder can still be released. The registry is
consulted BEFORE any channel is opened (POSIX drops every lock a process
holds on a file when any descriptor on it closes), and an overlapping-lock
report from Java parks the second channel rather than closing it and
errors by name. `reset` takes the unit lock and holds it while it deletes
records, and refuses ids not in the queue. Test: A claims and releases,
B claims, A's old struct cannot publish, beat or release (B's record
survives); `clear functions` leaves the registry intact; a path alias
maps to the same lock; and, across PROCESSES, a second MATLAB is still
refused after the owner clears its functions and probes its own lock --
the sequence that freed the lock this morning.

**Beats at every capped stage.** `certify_root` takes `.progress` and
ticks it after the polish and after each of the six fenced stages;
`rib_from_crossing` hands the worker's heartbeat down. The longest
silence between beats is now one stage cap (900 s) plus uncapped work,
and hangSec = 2700 s is three of those. `run_capped` VERIFIES a
cancellation: if the future has not finished 30 s after `cancel`, the pool
is deleted (its worker process with it) and the call errors as
infrastructure -- so a stuck native call fails the attempt for retry
instead of turning into short columns behind a poisoned one-worker pool.
`certify_root` refuses to run unfenced unless told to; the generated job
asserts a live pool before the worker reports READY.

**The launcher is the workers' parent.** Each worker runs under a
supervisor subshell that spawned it: it enforces the spawn-to-READY
deadline until the READY marker exists (whatever else is written), then
heartbeat inactivity; on either it sends TERM, waits up to 30 s, sends
KILL, kills the client's pool children, then `wait`s, so the exit code is
real and no recycled pid can be mistaken for the worker. Exit code and
reason go to `exit_<tag>`, the pid to `pid_<tag>`. READY is checked to
carry the spawned pid, and "attached then failed" is reported as such.

**One controller at a time.** The entry script holds `campaign.lock` for
the whole call, so manifest creation, launch and packaging are serialised.
Existing rib files are validated against their UNIT before the queue can
call them done (`rib_validate` now checks exactly one rib, the column, the
arrival phase, every point certified, points on the lattice and unique,
and the problem identity against the campaign's), and an invalid one is
quarantined by rename and the column re-queued. Foreign old-launcher
workers block packaging as well as launching. The barrier's index bug is
gone (one reason per column, parallel to the mask). 'packaged' is decided
by a RECEIPT the chain writes with this call's invocation id, not an
mtime. The chain takes the engine and orbits from the driver, gates
anchors/arcs/pool by the stages that need them, fences the audit like the
sweep, and closes only the figures it opened; the driver saves and
restores the figure default. The finish job appends verdicts and exits
non-zero unless packaged with a movie; the autochain checks the EXACT
expected rib files (`ribq/expected_units.txt`) and runs one instance at a
time. Audit- or sweep-only calls report state 'measured'.

**Verified.** `test_work_queue`: 49 checks. `test_campaign_processes`, four
phases with real MATLAB workers through the real launcher, 25/25: (0) the
cross-process registry check above; (1) three workers, unit 3's owner
killed -9 AFTER it wrote its temporary result -- the replacement published
its own result, not the dead owner's file; unit 5 wrote a partial result
and threw three times and left three `.failed` files; the only `.part`
debris belongs to the killed attempt; survivors exited rc 0, the killed
worker rc 137, all reaped by their supervisors; (2) maxAtt = 1, owner
killed mid-unit: RETIRED, finished, no restart; (3) a worker that beats
once and goes silent: the SUPERVISOR killed it after 69 s and recorded
"143 no heartbeat for 64s", the unit is abandoned and claimable, the
monitor alarms that no live worker remains, and no launched process
survived. The entry script dry-runs 'pending' against the live directory
(with the strict validator accepting all 16 real columns) and 'packaged'
with a receipt on the 15-column scratch copy, base workspace, cwd and
figure default intact.

**Still open (Astra's hardening list):** intra-column checkpoint/resume; a
campaign-wide supervisor that REPLACES lost capacity (idle workers cover
the tail, but a relaunch is a call to the entry script); status probes
are sampled, not a snapshot; the hardcoded pumpkynPie bootstrap; the
chain is still a script. The live 24x24 run is on the old launcher and
unaffected; its finish is armed on the new autochain.

## 57. The two gaps closed: a killed column resumes from its last point, and a supervisor keeps the campaign staffed and finalizes it (2026-09-13)

**Resume mid-walk.** `walk_checkpoint` (costate_common) saves the walker's
state after every ACCEPTED point -- the point index, the last certified
root and its multiple-shooting trajectory, the certified points so far,
the solve count -- atomically to `<output>.ckpt`. The checkpoint belongs to
the UNIT, not the attempt, and carries the walk's identity (arrival phase,
lattice, direction, targets, departure origin, problem). `rib_from_crossing`
resumes from it when the identity matches and ignores it by name when it
does not; `build_ribs` and the generated job pass it through; the worker
deletes it after the unit is published. Measured on a real rib: column 2
of the live 24x24 sheet, walked with a checkpoint, killed -9 after its
first accepted point; the next run logged "RESUMED at point 2 of 2",
finished, and its two points are BITWISE identical (max |dz| = 0, t_f to
all printed digits) to the same column walked unbroken by the live
campaign. The loss budget per kill drops from a whole column (up to nine
hours) to one point.

One bug caught on the way: `load` on a file named `.part` or `.ckpt`
without `'-mat'` reads it as ASCII and throws. `rib_validate` loads the
worker's temporary result, which is named `.part`, so every real unit
would have been refused publication. Both loaders now pass `'-mat'`, and a
test saves a valid rib under a `.part` name.

**Supervisor.** `campaign_supervisor.sh` is the long-lived controller: once
a minute it counts live workers from the launchers' pid files, reads the
queue's completion from files (every expected output exists, or the missing
ones have spent their attempts), relaunches the deficit through
`run_campaign_workers.sh` within a budget (3 x N launches; three
consecutive failed launches stop it), and when every output exists runs the
finalize job exactly once and exits with its status. All-retired remainder
is BLOCKED (exit 4, no finalizer). One instance per campaign directory.
`run_costate_library` with `.launch = true` now starts this supervisor, and
writes `finalize_job.m` -- itself, with the packaging stages on -- for it
to run. Measured (`test_campaign_processes` phase 4): four units, two
workers, one killed -9 mid-campaign; the supervisor launched a replacement
(8 launches in all across the test), every unit was published, the
finalizer ran exactly once, the verdict reads "FINISHED: all outputs
present; finalizer exited 0", and the lock was released.

**Tests.** `test_work_queue` 56 checks (checkpoint identity refusal, `.part`
validation). `test_campaign_processes` five phases, 26/26.

**Live 24x24 run:** unaffected; still on the old launcher, 16/19, finish
armed on the autochain. The supervised launch is for the next campaign.

## 58. The 24 x 24 library is built, audited, swept and filmed (2026-09-13)

The doubled-resolution 70 mN DRO -> 7-petal tulip library is complete:
`indirect/results_fine/costate_catalog_dro_tulip_70mN.mat`, and the full
phase-sweep movie `indirect/results_fine/sweep_full_library_24x24.{mp4,gif}`.
The finish job was started by `batch/fine_library_autochain.sh` at 16:49,
when the last rib column had landed and every old worker had exited, and
it finished at 21:45 with exit 0 and "CHAIN CLEAN: every stage that ran,
passed."

| | 12 x 12 (2026-09-11) | 24 x 24 (this) |
|---|---|---|
| certified entries | 115 | **406** of 576 cells |
| certified arrival columns | -- | 19 of 24 |
| t_f range | 16.23 .. 26.43 d | 16.23 .. 26.84 d |
| audit (rebuild from keys, fly, witness, conj) | clean | **406 OK / 0 bad** |
| second-order sweep | 0 sign changes | **0 interior crossings, 0 zeros, 0 unresolved**, worst H6 margin 5.32x, worst lift margin 23.7x |
| conj PASS | 115 | 406 |

Every column passed the unit-aware validator before packaging (nothing
quarantined), and the receipt ties the catalog to this invocation's 19 rib
files. The movie is 576 frames (406 certified, 170 gap frames for the
uncertified cells), 1280 x 720, 3 fps with slow = 2, 192 s; frames checked
at a certified phase pair and at a gap.

**Coverage blockers (data, not errors):**
- five arrival phases have no certified seed on the arrival sheet, so no
  rib: sA = 0.8254, 0.8671, 0.9504, 0.9921, 0.0337 -- the band past the
  fold where the fast family ends;
- ten columns stop short of 23 points where the walker could not certify
  the next point: polish non-convergence on most (cols 12-18), a
  transversality margin of 1.45e-6 against 1e-6 at col 21 (2 points), and
  at col 6 the dense conjugate scan found a ZERO -- a genuine end of the
  certified minimum along that rib, not a numerical stall.

The ship decision (deliverable zip) is Mike's; the audit and sweep that
gate it are clean.

**Cosmetic, open:** the torus picture rendered in the batch job has light
grey axis text on a white ground (the batch MATLAB appears to use the dark
theme's text colour); the same picture rendered in the interactive session
had black text.

**How long it took:** the arrival sheet 2026-09-12 afternoon; ribs from
18:27 on 09-12 to 16:44 on 09-13 (about 22 h wall, with the incidents and
hand fixes in sections 50-53); finish (package, audit 406, sweep 406,
pictures, movie) 4 h 56 min.

## 59. Why arrival phases 0.825-1.034 are missing: two ends of the fast family, one gate artifact, one unexplored interval (2026-09-14)

Asked why the 24x24 library has no certified transfers at sA = 0.825,
0.867, 0.950, 0.992, 1.034. The sheet's own records answer most of it;
GPT-6 Astra (`reviews/arrival_gap_astra_2026-09-14.md`, xhigh, 526 s,
$0.88) corrected two of my readings and asked for one experiment, which
settles the largest part.

**The fast family (A1, anchored at sA = 0.0754, 17.8 d) spans sA = 0.034
.. 0.841 and ends differently at each end.**
- High end: the +sA arc reached q = 0.8410 (t_f 25.72 d) and turned back
  with rho = 0.098 -- a FOLD of the endpoint projection, not a loss of
  normality. It crossed 0.8254 on the way out and back, which is where the
  two nearly equal extremals at that phase (25.503 and 25.507 d) come from
  (I had attributed them to the second family; wrong). One is refused for
  a 1160 km lunar approach (< 1900 km clearance), the other because the
  dense conjugate scan left one candidate unresolved -- consistent with a
  conjugate point approaching t_f near the fold.
- Low end: the -sA arc reached q = 0.0342 (t_f 19.32 d) with rho -> 0: a
  LOSS OF NORMALITY, the same way the family ends along the thrust axis at
  72 mN. The 0.0337 grid level lies 5e-4 in phase beyond it. So the fast
  family does not continue backward across the seam to 0.992/0.950
  (Astra's proposed search; the arcs had already made it, and it ends).
- 0.8671 lies strictly between the fast family's fold (0.841) and the
  second family's lowest reach (0.9084): no arc has entered that interval.
  Astra's intermediate-value argument does not apply because no arc joins
  0.909 to 0.825 -- but the interval is unexplored, not empty.

**The second family (A2, anchored at 0.9087, 26.43 d) is certifiable at
0.950, 0.992, 1.034 -- the refusal is an integration artifact, measured.**
The transversality gate reads lambda_m(t_f) off pumpkyn's tfMinProp flight
(ode45, RelTol 1e-10, AbsTol 1e-12) over 27-28 days. Re-integrating the
same initial state with ode113 at RelTol 1e-13 / AbsTol 1e-16:

| col | t_f [d] | gate value (ode45 1e-10) | tight (ode113 1e-13) | ms residual |
|---|---|---|---|---|
| 21 (certified) | 26.430 | 8.5e-7 | -3.7e-8 | 4.4e-13 |
| 22 | 26.702 | 1.437e-6 (refused) | -3.7e-8 | 1.2e-12 |
| 23 | 27.472 | 2.323e-6 (refused) | -6.4e-8 | 3.0e-12 |

The gate values reproduce the recorded refusals to three digits; at tight
tolerance the miss falls 40x, well inside 1e-6, and the certified 26.43 d
neighbour had passed by 15%. Astra's independent identity
p_m(0) = int_0^tf (T/m^2)|p_v| dt (exact when p_m(tf) = 0) holds on the
tight flights to five digits (7.0525 vs 7.0525; 3.5717 vs 3.5717). The
departure rib at 0.909 died on the same gate after two points. Their rho
is 0.0138 -> 0.0126, falling, so this family is nearly abnormal and may
end by normality loss somewhere past 1.03 -- to be found, not assumed.

**What is NOT the cause:** tulip geometry. The gap phases sit 25-38 Mm from
the Moon at 0.14-0.39 km/s; certified phases include a 6.6 Mm, 1.12 km/s
perilune pass.

**Astra's corrections, kept:** the value function is exactly periodic in
sA (the 17.8 d solution at 0.0754 is a solution at 1.0754), so a slower
family at a phase is never "the next turn" until the fast family's
continuation to that phase has been searched -- here it has, and it ends
by normality loss at 0.034. A fold is a singularity of the endpoint
projection with a terminal Jacobi degeneracy, not a reachable-set boundary.
"Certified" means the local second-order verdict; nothing here proves a
global minimum among disconnected extremals. Keep separate masks for
extremal found / numerically resolved / locally optimal / clearance-
admissible / best time found.

**Routes, in the order I would take them (not applied):**
1. Integrate the pointwise PMP checks at tight tolerance (ode113 1e-13)
   and report the loose-tight difference as the numerical uncertainty of
   each gate value, instead of trusting one ode45 1e-10 flight. Re-certify
   the three A2 candidates and re-walk the 0.909 rib. Hours.
2. Continue A2 downward from 0.909 into (0.841, 0.908): the interval no
   arc has entered; a crossing at 0.8671 either exists there or A2 ends.
3. Resolve the 25.503 d candidate at 0.8254 with a finer conjugate floor;
   if it is a conjugate point at t_f, that is the fold's signature and the
   phase's minimizer is on the other side (the Moon-grazing one, refused
   by a real constraint) -- then 0.8254 needs a clearance-constrained
   formulation, not a looser gate.
4. Direct multistart solves (collocation, free t_f) at 0.8671 and 0.8254
   with the clearance enforced, then polish and certify: the search that
   is blind to which branch the continuation happened to follow.
5. Map fold and normality-loss loci in (thrust, sA) before assuming a
   fixed-thrust gap is intrinsic.
6. Record the fast family's two ends (fold at 0.841 / 25.7 d, normality
   loss at 0.034 / 19.3 d) as library metadata: the t_f jump to the A2
   family at those phases is a property of the problem, not a hole.

## 60. Executing the gap plan: the corrected gate exposes saddles, the second family folds at 0.908, and the direct solver takes over (2026-09-14, in progress)

**Item 1, the certifier.** `certify_root` now runs the pointwise PMP checks
(H, adjoint, transversality, control law) on a tight flight -- ode113,
RelTol 1e-13 / AbsTol 1e-16 on the same field -- and records the loose
flight's value and the loose-tight gap as the check's numerical
uncertainty (`.lamMfLoose`, `.lamMfUnc`). The arrival gate keeps pumpkyn's
own flight. Re-certifying the four A2 candidates from their stored seeds:

| sA | t_f [d] | transversality now | verdict |
|---|---|---|---|
| 0.9087 | 26.430 | 3.7e-8 (unc 8.1e-7) | certified, as before |
| 0.9504 | 26.702 | 3.7e-8 (unc 1.4e-6) | **conjugate test verdict 0** |
| 0.9921 | 27.472 | 6.4e-8 (unc 2.3e-6) | **conjugate test verdict 0** |
| 1.0337 | 28.124 | 6.2e-8 (unc 2.2e-6) | **conjugate test verdict 0** |

So the gate artifact was real, and removing it recovers nothing: past
0.909 the A2 family is a saddle (a conjugate point before t_f). The
first-order gate had been hiding a second-order refutation.

**Item 2 is impossible.** The A2 continuation in -sA went 0.90873 ->
0.90837 and turned (tangent through zero at step 5): the family FOLDS at
0.9084. Between the fast family's fold at 0.841 (25.7 d) and this one at
0.908 (26.4 d) no known extremal family exists.

**Item 3.** The 25.503 d candidate at 0.8254 stays UNRESOLVED with the
conjugate floor lowered from 1e-7 to 1e-8: the near-zero is below 1e-8 --
a genuine conjugate point at (or just before) t_f, the terminal Jacobi
degeneracy Astra predicted next to the fold. Not a policy floor.

**Item 4, in progress.** The direct sweep harness (`sweep_phasing_direct`)
was the wrong tool: with a grid start different from its own anchor it
cold-seeds its first wave, and its first point at 0.0754 came back as a
55-day path through the Moon (five minutes per failure, twelve to go).
Stopped. Replaced by a focused probe: Hermite-Simpson + Sundman direct
solves at the five gap phases and two control phases, each warm-started
from the certified indirect trajectory on either side of the gap (0.7837
and 0.9087), with the 1900 km clearance enforced as a path constraint.

**Round 2 of the library, in progress.** `results_fine_v2`: the arrival
sheet rebuilt under the corrected certifier (every crossing re-certified),
the 18 valid ribs copied in, column 21 re-walked through the supervised
campaign (its rib had died on the same gate), finalized by the supervisor.

## 61. A faster family the continuation never found: 17-19 d transfers certified across the gap, seven days under the library's swept family (2026-09-14)

Direct Hermite-Simpson + Sundman solves (N = 800, IPOPT, the 1900 km lunar
clearance as a path constraint), warm-started from the certified 26.43 d
solution at sA = 0.9087, then harvested to multiple-shooting seeds
(`harvest_ms_seed`, sign vote -1 on every one) and put through the full
certifier (`certify_root`, corrected pointwise flight):

| sA | direct t_f | polished, CERTIFIED t_f | library's certified t_f |
|---|---|---|---|
| 0.7837 | 17.834 d | **17.8336 d** | 24.737 d |
| 0.8254 | 18.738 d | **18.7383 d** | -- (fold pair refused) |
| 0.8671 | 17.249 d | **17.2487 d** | -- (no extremal) |
| 0.9087 | 19.042 d | **19.0417 d** | 26.430 d |
| 1.0337 | 18.143 d | **18.1434 d** | -- (A2 saddle) |
| 0.9504 | 39.622 d (perilune on the floor) | polish did not converge | -- |
| 0.9921 | 72.068 d | polish capped | -- |

Every certified row passed the multiple-shooting residual, the flown
arrival, the tfMin witness, the pointwise PMP checks, the free-time
conjugate test, the sufficiency gates and the clearance; perilunes 6-8 Mm
except 0.8671 at 2.0 Mm (just clear). Also certified, from the fast
family's own seed at 0.7837: a 23.2086 d Moon-grazing transfer (perilune
1932 km), 1.5 d faster than the library's 24.737 d there.

**What this means.** The value function on the arrival axis is NOT the
family the pseudo-arclength continuation from the 0.0754 anchor traced.
That family is a local minimum everywhere it was certified, but from at
least sA = 0.78 to 1.03 a different family is 5-7 days faster and was never
reached, because the continuation follows its own sheet and both sheets
end in folds or normality loss before meeting. Astra's pass on the gap
said exactly this: "certified means local; a slower family at a phase is
never the answer until the fast continuation to that phase has been
searched" -- and a branch-blind direct solve is that search. The direct
sweep harness was the wrong tool (cold seeds); the right one was a direct
solve warm-started from the nearest certified trajectory of the OTHER
family.

**The 0.9504 / 0.9921 columns.** From the 0.9087 seed the solver lands on
a 39.6 d solution with the clearance ACTIVE (perilune on the 1900 km
floor) that the unconstrained PMP polish cannot reproduce -- an active
path constraint needs the constrained PMP, as Astra noted -- and on a 72 d
solution at 0.9921. Being re-probed from the new certified neighbours
(0.9087 at 19.04 d, 1.0337 at 18.14 d).

**Consequences for the library, and the plan.** The 406-entry 24x24
library and its round-2 rebuild are certified local minima of the SWEPT
family; over part of the range they are not the fastest certified
transfer. The torus is to be filled with the new family: a new anchor
(`results/mintime_70mN_anchor_fast2.mat`, sA 0.8671, 17.25 d),
pseudo-arclength arcs from it in both arrival directions, a sheet rebuilt
from ALL arcs so each column keeps its fastest certified solution, then
the departure ribs re-walked for every column whose spine changed, under
the supervised chain. Column 21's round-2 rib, under the corrected gate,
certified 5 points (was 2) and stalled at sD = 0.771 on an unresolved
conjugate candidate: that family's second-order margin ends there too.

## 62. Round 3 begins: the new family seeded into the sheet; 0.992 certified, 0.950 held by the lift margin (2026-09-14)

Third direct probe, seeded from the new family's certified neighbours:

| sA | seed | direct t_f | polished | verdict |
|---|---|---|---|---|
| 0.9504 | 0.9087 (19.04 d) | 19.988 d | 19.9881 d | REFUSED: lift margin -- dim S = 1 established by a factor 3.1, gate requires 10 |
| 0.9504 | 1.0337 (18.14 d) | 43.48 d | -- | slow branch, not polished |
| 0.9921 | 0.9087 | 64.82 d | -- | slow branch, not polished |
| 0.9921 | 1.0337 | 18.303 d | **18.3032 d** | **certified** |

So every one of the five gap phases now has a certified transfer of the
new family except 0.9504, whose 19.99 d candidate fails only the
normality-rank margin (a policy gate: the abnormal-lift dimension is 1 by
a margin of 3.1x, and the rule demands 10x; FINDINGS 30). It is kept as a
candidate, not a library entry.

**The new family's arcs.** From the 0.8671 anchor the +sA arc passed
1.034 with no fold by step 810 (still climbing, toward the phases the old
family holds at 16-24 d); the -sA arc folded five times in a tight
cluster at sA = 0.7995-0.8050 and is walking back up -- so the 17.83 d
solution at 0.7837 and the 18.74 d at 0.8254 are on sheets the direct
solver could jump to and this arc cannot reach. The family structure is
richer than one sheet.

**Round 3.** The seven direct-found certified solutions are stored in
`results/mintime_70mN_direct_certified.mat` and `dro_tulip_library` now
lists them (src 'direct_certified', z8 only; the seed is rebuilt by
seed_from_z8). The round-3 sheet (`results_fine_v3`) is being rebuilt from
all arcs plus these seeds; each column keeps its fastest certified
solution. Columns whose spine changes get fresh ribs under the supervised
chain; the old family's rib files are handed to the packager as extra rib
files, so per cell the fastest certified point of either family is kept
(`run_costate_library .extraRibFiles`).

## 63. Round 3 launched on the faster family; the seed phases had to be exact; two entries held by the lift margin (2026-09-14)

**A 3e-5 phase mismatch cost a sheet rebuild.** The direct probe's targets
were the grid phases rounded to four decimals (0.8671 for 0.867067, ...),
so the certified solutions sat 3.3e-5 off the lattice and the sheet
builder's seed filter (1e-6) took only the one that happened to be exact
(0.8254). Re-certified at the exact phases from the direct solutions'
harvested multiple-shooting seeds (a seed rebuilt from z8 alone did not
converge at 0.909 and 0.992):

| col | sA | t_f [d] | verdict |
|---|---|---|---|
| 18 | 0.783733 | 17.8342 | certified (also the 23.2092 d Moon-grazing one) |
| 19 | 0.825400 | 18.7383 | certified |
| 20 | 0.867067 | 17.2482 | REFUSED: lift margin 8.0x (rule: 10x); certified at 0.8671 |
| 21 | 0.908733 | 19.0174 | certified |
| 22 | 0.950400 | 19.9881 | REFUSED: lift margin 3.1x |
| 23 | 0.992067 | 18.3046 | certified |
| 24 | 0.033733 | 18.1431 | certified |

`results/mintime_70mN_direct_certified.mat` now holds the six certified
ones at exact phases. The 0.867067 refusal is a near miss of a policy
threshold on a noisy estimate (the abnormal-lift rank margin moved from
above 10 to 8.0 over 3e-5 in phase); it is recorded as a candidate, not
a library entry, and the estimator's jitter is an item for review.

**Round 3 (`results_fine_v3`).** The sheet was patched under the campaign
lock with the exact-phase certified solutions as candidates (each column
keeping its fastest): columns 18 (24.74 -> 17.83 d) and 24 (new, 18.14 d)
changed, 19 was already new (18.74 d); 21 and 23 missed in the patch job
(the z8-only seed) and join round 4. Ribs of the 18 unchanged columns
were reused from round 2; the supervised campaign is walking 18, 19 and
24 on four workers, with the round-2 rib files handed to the packager so
the earlier family's certified points survive per cell.

**Round 4.** When the new family's arcs finish: rebuild the sheet from all
arcs plus the six exact-phase seeds, reuse every rib whose spine is
unchanged, walk the rest (21 and 23 at least, more if the arcs lower the
0.075-0.74 spine), package with both families' ribs.

**Also seen along the arcs (in progress):** the +sA arc from the 0.8671
anchor passed 1.27 (= 0.27) with no fold and a smallest singular value
falling to 3e-9 -- an ill-conditioned stretch or a fold ahead; the -sA arc
folded five times in 0.7995-0.8050 and is walking an S-bend back up.

## 64. The faster family mapped by continuation, and its departure ribs (2026-09-14)

**The arcs from the 0.8671 anchor (17.25 d).** Both ran to their 4001-step
budget.
- +sA (`arrival_arc_fast2_up_long`): one smooth sheet from 0.8671 to a fold
  at q = 1.394 (= 0.394), t_f rising 17.25 -> 21.6 d and rho falling 0.137
  -> 0.05. Level crossings: 0.9087 at **17.92 d** (the library's old spine
  26.43; the direct-found 19.02), 0.9504 at **18.05 d** (no entry before),
  0.9921 at 18.36 d, 1.0337 at 18.66 d (the direct-found 18.14 is on
  another sheet), then 1.0754 at 18.86 d against the old family's 17.80 d
  -- so the crossover between the families lies near sA = 1.05-1.07, and
  from 1.075 to 1.37 (= 0.075-0.37) the old family is the faster one
  (16.2-18.3 d against 18.9-21.3 d).
- -sA (`arrival_arc_fast2_dn_long`): an S-bend. Five folds in 0.7995-0.805,
  then back up through the anchor's phase (19.3-19.9 d), a fold at 0.9635,
  and down again to end at 1.020 (25.3 d). Its crossings include **16.87 d
  at 0.8254** (the certified direct-found 18.74 there is on the other
  side of the bend), 20.4-20.5 d at 0.9504, and a slower 23-25 d branch.

**Round 3 ribs on the new spines** (four workers, supervised): column 18
(17.83 d) walked all 23 departure points -- COMPLETE; column 19 (18.74 d)
17 of 23, stalling at sD = 0.254 on polish non-convergence; column 24
(18.14 d) 1 of 23: it could not step 1.6e-4 in departure phase, so that
solution is isolated in sD. The sheet's fastest solution at a phase is not
necessarily the one with the widest certified neighbourhood.

**A finalizer bug, caught by the run.** The generated finalize job wrote
its extra rib files as a bare cell inside `struct(...)`, which in MATLAB
builds a struct ARRAY; the entry script's option reader failed on it in
six seconds. Fixed with double braces (`{{...}}`); round 3's finalizer was
patched by hand and rerun; the fix is in the generator for round 4.

**Round 4** is rebuilding the sheet from all six long arcs plus the exact-
phase seeds; every column keeps its fastest certified solution, changed
columns get fresh ribs, and the packager receives rounds 2 and 3's rib
files so per cell the fastest certified point of any family survives.

## 65. The arrival axis is complete: 24 of 24 phases certified (2026-09-14)

The round-4 sheet, rebuilt from all six long arcs (both families) plus
the exact-phase direct-found seeds, certifies a minimum-time transfer at
every arrival phase of the 24-grid at sD = 0. The five phases that were
empty or dominated yesterday now read:

| col | sA | before | now | source |
|---|---|---|---|---|
| 19 | 0.8254 | 18.738 (direct) | **16.870 d** | -sA arc crossing (the S-bend's lower sheet) |
| 20 | 0.8671 | none | **19.323 d** | -sA arc crossing (the 17.25 d anchor itself is refused by the lift margin at the exact phase) |
| 21 | 0.9087 | 26.430 | **17.922 d** | +sA arc crossing |
| 22 | 0.9504 | none | **18.054 d** | +sA arc crossing (the 19.99 d direct-found one fails the lift margin) |
| 23 | 0.9921 | none | **18.305 d** | +sA arc crossing / direct (equal) |
| 24 | 0.0337 | none | 18.143 d | direct-found (round 3) |
| 18 | 0.7837 | 24.737 | 17.834 d | direct-found (round 3) |

Columns 1-17 keep the original family (16.2-24.1 d): the new family's
+sA arc is slower there (18.9-21.3 d), so the crossover between the two
families sits between 0.784 and 1.075 on one side and near 0.78 on the
other. Everything the direct probe found, the continuation of the family
it exposed then found too, and in two columns found faster.

Round 4 walks the five changed columns' ribs on four supervised workers,
reuses the other 19 rib files (round 3's, which include the complete
23-point rib at 0.7837), and packages with rounds 2 and 3's ribs beside
its own so each cell keeps the fastest certified point of either family.
No blockers at launch.

## 66. The family map: which family each entry belongs to, and where each family ends (2026-09-14)

Item 5 of the arrival-gap plan. `family_map` reads the extremal families
off the stored arcs -- a family is one anchor, both walk directions --
and records per family the arrival-phase span its arcs reached, the
final-time span, every fold, the normality floor (min |rho| and where),
and the KIND of each end of the span: a fold (the branch turns back), the
walk budget (the branch continues, unmapped), a walk that ended with
|rho| -> 0, or the anchor itself (that direction was never walked). It
then attaches every certified root at every column to the arc that
passes through it (t_f interpolated along the arc; distinct roots at one
phase sit >= 0.05 d apart here, the arcs interpolate to 1e-5 d), and
every rib to the spine root it was walked from (t_f of its first three
points extrapolated quadratically back to the spine, matched to the
column's certified roots, refused when ambiguous). The catalog ships the
map at `cat_.families` and an int8 `family_index` per entry (1..n a
mapped family, -1 a certified root no mapped arc passes through, -2 a
rib whose spine root is unidentified, 0 no entry); the schema validates
the pair. The chain names its families in the anchors table, where the
fast2 anchor is now declared beside the other two.

Measured on the round-4 sheet (10 arcs, 24 columns):

| family | anchor | sA span | t_f span | ends |
|---|---|---|---|---|
| fast | 0.0754 (17.80 d) | 0.0342 .. 0.8410 | 16.22 .. 25.72 d | fold at 0.0342 (19.32 d); fold at 0.8410 (25.72 d) |
| A2 | 0.9087 (26.43 d) | 0.9084 .. 1.7109 | 25.98 .. 44.05 d | fold at 0.9084 (26.52 d); walk budget at 1.7109 (44.05 d) |
| fast2 | 0.8671 (17.25 d) | 0.7995 .. 1.3942 | 16.86 .. 25.25 d | fold at 0.7995 (20.36 d, the S-bend, min rho 3.8e-5 at 0.801); fold at 1.3942 (21.71 d) |

Columns 1-17 attach to `fast`, 19-22 to `fast2`, and three columns
attach to NOTHING: the direct-found roots at 0.7837 (17.83 d), 0.9921
(18.30 d) and 0.0337 (18.14 d). They are not on the fast2 arc -- fast2
passes 0.9921 at 18.36 d and 1.0337 at 18.66 d, both certified in the
same sheet as distinct roots (sheet_from_arcs merges only equal t_f AND
equal z8) -- and they are faster than it by 0.05-0.5 d. So the direct
solver found a FOURTH branch that no arc has walked, the one that won
the 0.7837 column outright. Two arcs from the 0.7837 root
(`arrival_arc_direct18_{up,dn}`) are walking now; the next sheet rebuild
takes them up automatically.

Where the families' ribs land in the packaged catalog (round-3 dry run,
432 entries): 382 on `fast`, 1 on `A2`, 26 on unattached roots (the 0.7837
column's whole rib and the 0.0337 and 0.8254 direct roots), 23 ribs
unidentified (-2: ribs from an earlier round whose spine root is not
among this sheet's certified roots -- expected, and now visible).

## 67. The fourth branch mapped: it owns 0.59-0.78 by 4-7 days (2026-09-14, evening)

Both arcs from the 0.7837 root (17.83 d) are in, 4001 roots each.

**up** (0.7837 -> fold at 1.2242 -> back to 0.917): t_f climbs from 17.8
to 24.8 d. Slower than the sheet's best at every grid level it crosses
except 0.8671 (19.01 vs 19.32 d, marginal). Its second fold at 0.8617 sits
where the chart is weakest (min |rho| 0.015 at 0.873).

**dn** (0.7837 -> fold cluster at 0.573-0.592 -> a slow sheet up to
1.028): on the way down the branch is the fastest thing yet seen on
0.62-0.74 --

| level | col | direct18 root | round-4 best | gain |
|---|---|---|---|---|
| 0.7421 | 17 | **17.17 d** | 24.07 d | 6.9 d |
| 0.7004 | 16 | **16.80 d** | 23.47 d | 6.7 d |
| 0.6587 | 15 | **17.28 d** (17.38, 17.39 on the return legs) | 22.82 d | 5.5 d |
| 0.6171 | 14 | **17.27 d** (17.86, 18.13) | 22.15 d | 4.9 d |

-- then it loses normality (min |rho| 3.5e-5 at 0.5785) in a cluster of
five folds between 0.5734 and 0.5917 and comes back up on a sheet at
32-38 d that is slower than everything (0.5754 at 32.6 d, 0.9921 at 37.6 d).
So the branch's fast sheet spans about 0.59 to 0.78, and the map now has a
4 d step between col 13 (0.5754, fast family 21.45 d) and col 14 (17.27 d):
the obvious place for a warm-started direct solve, from the 17.27 d root,
to look for a fifth root below 21 d.

The direct-found roots at 0.9921 (18.30 d) and 1.0337 (18.14 d) are on
NEITHER direct18 arc (up passes them at 20.9 and 21.6 d, dn at 37.6 d), so
they remain unattached -- a branch no arc has walked.

Housekeeping the day paid for: the first dn walk died in a MATLAB segfault
at step ~2160 with nothing on disk, so `arclength_ms` now writes the arc so
far atomically every 50 steps (`.partialFile/.saveEvery`); the relaunch
gave the table above from its partial file twenty minutes in. Round 2
closed clean (409 entries, audit 409/0, sweep 0 crossings, worst H6 5.3x,
lift 24x). Round 4's ribs finished (cols 19-23 walked to 23/16/22/21/20 of
23 points; 491 rib points + 24 spines = 515 entries packaged) and its
finalizer is auditing. Round 5 -- the sheet rebuilt with the direct18 arcs,
ribs re-walked wherever a spine changes, every earlier round's ribs offered
to the packager -- launched 18:25.

## 68. Round 5, and a fifth root at 0.4921 (2026-09-14, night)

Round 5's sheet (12 arcs, every earlier seed) certified 24/24 phases and
the fourth branch took five columns from the fast family: 0.6171 17.273,
0.6587 17.281, 0.7004 16.797, 0.7421 17.169 (were 22.2-24.1 d) and 0.8671
19.009 (was 19.323). Ribs for those five columns are walking; the other
19 ribs are reused from round 4. Round 3 closed clean (432 entries, audit
432/0, sweep 0 crossings); round 4's audit is at 500/515.

The direct probe at the branch's lower fold (0.5754, 0.5337, 0.4921,
seeded from the certified 17.27 d and 17.28 d roots):

| phase | seed | result |
|---|---|---|
| 0.5754 | 0.6171 | 20.00 d with the periselene AT the 1900 km floor (162 km altitude) -- a clearance-constrained arc, not an unconstrained extremal; the certifier would refuse it. The sheet's 21.45 d stands. |
| 0.5754 | 0.6587 | failed |
| 0.5337 | both | 76 d / 52 d, junk |
| 0.4921 | 0.6171 | **17.959 d, periselene 10,500 km** -- harvested, ms-polished, certified at the exact grid phase in 33 s: 17.9605 d, fly 0.000 km, conj PASS, lift 9.0e3x. The sheet has 20.220 d there. |

So the fast sheet is not one branch: below the fourth branch's fold at
0.59 there is another root 2.3 d under the fast family, and at 0.5754 the
fastest thing the direct solver finds grazes the Moon. The 0.4921 root is
in `mintime_70mN_direct_certified.mat` (the next sheet rebuild seeds from
it) and is the anchor of two more arcs (`arrival_arc_direct11_{up,dn}`,
launched 20:40, partial saves on) that will say whether it is the fourth
branch continuing past its fold on another sheet or a family of its own.

## 69. Round 5 packaged: 534 entries, four families in the map (2026-09-14, 23:00)

Round 5's ribs on the fourth branch's five columns: 0.6171, 0.7004, 0.7421
and 0.8671 walked COMPLETE (23/23 each -- the branch's ribs are long where
the fast family's stalled at 18-21 points); 0.6587 stalled after ONE point
("normal-chart polish did not converge" at sD = 0.9224 from the 17.281 d
spine). That column is the one where three roots sit within 0.11 d of each
other (17.281, 17.380, 17.390 -- the S-bend's three sheets crossing the
same level); the walker only follows the sheet's winner, so the packager
keeps round 4's 22.8 d rib there. Two remedies, in order of cost: walk the
rib from the column's OTHER certified roots when the winner's rib stalls
(a chain feature), or direct-solve the cells next to the spine seeded from
the neighbouring columns' rib points (17.3-17.4 d at the same sD).

Packaged: 534 entries (round 4: 515), family map with four families (fast,
A2, fast2, direct18; the direct18 span reads 0.5734..1.2242 with ten folds,
its lower end the 33.7 d slow sheet's fold). Audit running.

Also: `build_arrival_sheet` now ignores `<arc>.partial.mat` files -- a walk
in progress would otherwise enter the sheet as an arc (and twice, once it
finishes), and `family_map` would refuse the name.

## 70. The torus is full: 576 of 576 cells certified (2026-09-15, 09:10)

Round 6 closed clean (538 entries, audit 538/0, sweep 0 sign changes, worst
H6 5.32x, worst lift 11x); its ribs on the fifth branch's three columns
(0.4921 17.96, 0.5337 18.23, 0.5754 18.83 d) walked complete, 23/23 each.
That left 38 holes -- the tops of columns 6, 7, 15 (three cells each), the
whole of column 24 above its lone rib point (22 cells), and a few more --
which no continuation rib had reached.

`fill_holes_direct` (new) solved them cell by cell: a direct HS+Sundman
solve at the hole warm-started from up to three certified neighbours (the
same column's first -- they share the arrival geometry -- then the same
row's), harvested into a multiple-shooting seed and put through the full
gate stack; the fastest certified root kept and made a seed for the cells
still to come, so a column chains like a rib. **38 of 38 certified.**
Column 24 chained 22 cells at 17.7-18.5 d from its single rib point; the
(0.0417, 0.6587) cell certified at 17.25 d, the fourth-branch root the
walker could not step to. Two lessons cost an hour each: run scratch
MATLAB jobs under R2026a (R2025b has no Parallel Computing Toolbox, so
every certification died on `gcp`), and never seed from the neighbouring
COLUMN first -- a 16.8 d root at the next arrival phase sent IPOPT to
100-340 d junk and burned the 900 s cap, while the same column's root
converged in 1-3 min.

Round 7 (this morning) is round 6 re-packaged with the direct-hole rib
offered beside every earlier rib: **576 of 576 cells, 576 entries, t_f
16.23-24.91 d (mean 18.65, median 18.48)**, five families in the map (fast
44%, direct18 21%, fast2 13%, direct11 12%, A2 0; 22 cells on roots no arc
has walked, 40 rib points whose spine root the map could not identify --
mostly the direct-hole cells, which have no spine). Schema clean. Audit
and sweep running.

What "full" does not mean: the fastest known root everywhere. Column 15
(0.6587) reads 17.2 d at its spine and its first cell and 22-24.9 d
above that, because the fourth-branch rib stalled at once and the fast
family's rib filled the column; its neighbours are at 17-18 d. An
"improve" pass -- the same cell-by-cell solve at every cell slower than
its column neighbour by more than a threshold, chained upward from the
fast root -- is the next step, and the machinery is the hole filler with
one more selection rule.

## 71. Round 7 audited clean; the improve pass; the sidecar merge (2026-09-15, noon)

Round 7 (576 entries) audited **576 ok / 0 bad**. Its sweep refused its
own sidecar -- 538 records from round 6's sweep against 576 entries -- a
guard doing exactly what it was written for (FINDINGS 44: a positional
sidecar must never hand one entry's measurements to another). Records are
keyed by cell and z8, so the safe generalisation is a MERGE: reuse a
record only for the entry whose cell and z8 it measured, start every
other entry fresh, keep the old file, write the re-keyed one before
measuring. `second_order_pass` does that now (same-count mismatches are
still refused: they mean a different build); test case 6 covers it.

`fill_holes_direct .improveDays`: the same cell-by-cell solve at every
filled cell more than 2 d slower than a column neighbour, seeded from
faster neighbours only, a root kept only if faster than what the cell
holds, and the next slower neighbour queued when a cell improves -- so a
column is walked from one fast root. On column 15 (0.6587), where the
fourth-branch rib had stalled at once and the fast family's 22-25 d rib
filled the column, the chain climbed from the 17.25 d cell:

| sD | before | after |
|---|---|---|
| 0.1667 .. 0.5000 (9 cells) | 24.9 .. 24.0 d | **17.84 .. 17.06 d** |
| 0.9167 | 22.9 d | 18.63 d |
| 0.5417 .. 0.8750 (9 cells) | 23.9 .. 22.9 d | unchanged: the direct solve from the 17.06 d cell diverges (168-275 d) even at a 900 s cap, from below and from above |

plus (0.0417, 0.3671) 20.33 -> 18.31 d. Eleven cells improved by more than
two days; 49 direct points in the rib file. The nine stuck cells look
like a fold of the fourth branch's sheet in the departure phase between
sD 0.50 and 0.54 at this arrival phase -- the same wall the rib walker
hit -- so the fast family's 23 d roots may be the true minima there, or
the fast sheet may continue past a fold the direct solver cannot jump.
An arc in sD at fixed sA = 0.6587 would settle it; left open.

Round 8 = round 6 re-packaged with the improved rib file, then audit and
sweep with the merged sidecar; chained behind the last improve pass.

## 72. Round 8: the 24 x 24 library is finished (2026-09-15, 14:42)

Round 8 -- round 6 re-packaged with the direct-hole rib and the improved
column 15 -- closed clean: **576 entries, 576 of 576 cells; audit 576 ok /
0 bad; second-order sweep 576 done (sidecar merged: 527 records kept, 49
measured fresh), 0 interior conjugate crossings, worst H6 margin 4.44x,
worst lift margin 10.6x; schema clean.** t_f runs 16.23-23.86 d (mean
18.52, median 18.42). Families: fast 42%, direct18 21%, fast2 12.5%,
direct11 12%, A2 none; 22 cells on roots no arc has walked and 50
direct-solved cells with no spine to attribute. Every column is now within
2.5 d of its own minimum except column 15, whose nine stuck cells (sD
0.54-0.88) still carry the 23 d family (FINDINGS 71).

The library of record is `DRO_tulip/indirect/results/library_70mN_24x24_final/`
(catalog with the family map, receipt, keyed sidecar, the sD = 0 sheet
with every candidate, the direct-cell rib, the finalizer log, the torus
pictures); the rounds' rib files stay in `results_fine_v2..v6`.

What it took, from the 406-entry, 19-column library of 2026-09-13: five
families instead of one (three found by branch-blind direct solves
warm-started from the wrong family, then mapped by arcs), a family map
that attributes every entry and names each family's ends, the rib walker
resumable per point under a supervised queue, partial saves for the arcs,
a cell-by-cell direct filler with an improve pass, a sidecar that merges
across re-packagings, and the chain's git and pgrep calls made unable to
hang or lie. Open, none blocking: the nine column-15 cells (an sD-arc at
sA = 0.6587 would say whether 23 d is the minimum there), arcs from the
0.9921 / 0.0337 roots to attribute the last 22 unattached entries, the
deliverable zip and the ship decision.

## 73. run_phase_torus: the driver, its first acceptance run, and the Astra review (2026-09-15)

`run_phase_torus(spec)` is the one-call route to a phase torus over any
strictly increasing lists of departure and arrival phases (the lattice
assumption is out of the six files underneath; `rib_targets` owns the
rib-step rule). A 3 x 3 acceptance torus from one anchor with small budgets
ran the whole loop unassisted: arcs, sheet, ribs, finalizer, holes,
re-package, then DISCOVERY found the 17.705 d root at 0.7421 (the direct18
family) from the fast root, anchored it, and round 2 walked its arcs and
ribs and finished clean. Two defects surfaced on the way and were fixed:
the chain's grid-override block refused the new list fields, and a round
that packaged in-call was then judged by the previous attempt's FAILED
verdict file.

GPT-6 Astra (xhigh, 45 KB bundle, $0.96, 6 min) then reviewed the driver
and returned 21 findings (`reviews/run_phase_torus_astra_2026-09-15.md`).
Adjudication:

**Applied (17).** Campaign manifest in the state file, checked on resume
(a different grid, engine, orbit pair or tag is refused). Arcs live in the
campaign's own folder (`<outDir>/arcs`; `.arcDir` threaded through
run_costate_library, build_arrival_sheet, the chain, the packager and the
family map) so two campaigns can never consume each other's. Arc jobs
publish through a temp file and a rename and write a `.done` / `.fail`
verdict; the driver records their PIDs, waits on verdicts, fails fast on
a `.fail`, and kills what is left at the deadline. A live supervisor for
a round is adopted rather than relaunched. Anchor names are unique per
column (`d03_1`, `d03_2`, ...); the state is saved at every promotion.
"Spine unchanged" compares the root (t_f and z8), a lost spine is
reported, and a rib copied into a round for a spine that has since
changed is set aside. Every distinct certified probe root is registered
(`direct_certified.mat`, seeded into every later sheet); a spine root
found by the filler is registered too; promotion to anchor needs a gain
of `.acceptDays` (0.05 d, not 86 s) AND that no known family passes
through the root (`family_map` attachment) -- a root on an existing
family is a seed, not an anchor, so its arcs are not walked twice.
Discovery seeds come from every certified candidate (winning or not) of
columns within a circular phase radius, deduplicated, three per target;
`.probeAll` probes every column. The re-package stage asserts success.
One stopping rule: a round that registers no root and adds no anchor is
a fixed point ('done'); the last allowed round with pending work ends in
'budget' with the reason recorded; a finished campaign answers "nothing
to do" (plan mode too). The final folder is rebuilt fresh, required
products asserted, optional ones listed when absent. Shell paths quoted,
the startup folder and MATLAB binary are options, `.orbits`/`.engine`
default to the 70 mN campaign, `.rib.wallSec` reaches the rib jobs.
Discovery's direct solve is wrapped in try/catch. The phase lists must
be at least 1e-5 apart (the shared matcher resolution).

**Declined or accepted as limitations (4).** The clearance-floor
pre-check in `direct_cell_solve` stays a 5 km heuristic: the certifier
re-flies every candidate on ode113 and gates the continuous clearance at
1900 km itself, so the pre-check only spares a certification. The
direct-cell attempt has no outer wall deadline beyond the solver's CPU cap
and the certifier's wall cap; a hang in propagation or harvest would stall
the driver -- noted, not fenced (the fence's pool cannot be nested). The
arc jobs' level ladder now spans every whole-period copy the walk can
reach, but the sheet does not depend on it (its re-scan handles wrapping
itself). And the review's closing point stands as written: a certified
root faster than the spine is a minimum-eligible entry of the gate
stack, not a global-minimum certificate, and discovery is a heuristic
whose silence is evidence, not proof -- `.probeAll` is the wider net.

## 74. The study script reviewed: an enforcement gap in S4, a false branch, and a script that could not fail (2026-09-15)

GPT-6 Astra (xhigh, 53 KB bundle, $0.99, 7 min) reviewed
`transfer_study.m` -- the campaign's teaching artifact, whose contract is
that every necessary and sufficient condition is computed and gated IN the
script with its own PASS/FAIL. Fourteen findings
(`reviews/transfer_study_astra_2026-09-15.md`). Adjudication:

**Applied.**

*S4 could pass on a scan that was never testable.* The dense conjugate
scan's own verdict `CS.clear` was printed and never required; the
downgrade tested only the counters. `CS.clear` is
`testable AND nZero==0 AND nUnresolved==0 AND nInterior==0`, so a NOT
TESTABLE scan leaves every counter at zero and left S4 at PASS -- weaker
than `certify_root`, which refuses `~clear` outright. The script now
validates the scan's fields (scalar logicals, non-negative integer
counts), requires `.clear`, requires the flag to AGREE with its own
counters, and reports NOT TESTABLE as UNRESOLVED with the reason that no
zero was looked for, so none was excluded. This is exactly the failure
mode the study-script rule exists to prevent, in the script that states
the rule.

*A verdict branch claimed what it had not established.* Reaching the
`~crossCheck` branch establishes only `necessary`; the text said "every
PMP and sufficiency line passed, but the CROSS-CHECK failed" -- false
whenever a sufficiency gate had also failed. The branch now names the
actual sufficiency status (PASS / UNRESOLVED / FAIL) beside S4's and V1's.

*The script could not fail.* Every gate could print FAIL and the script
still exited 0. There is now a named gate table, a `studyOK` scalar, a
one-line self-check naming what did not pass, and an assertion -- with the
policy stated in section 0: a failed NECESSARY or CROSS-CHECK gate throws
(the root or the implementation is broken), while a failed or unresolved
SUFFICIENCY gate is a legitimate finding about this trajectory and is
reported. `selfCheck.strict = true` makes any outcome short of a claim an
error, for a regression harness. The plot is drawn before the assertion,
so a failing run still leaves the picture that explains it.

*The lunar clearance had no visible calculation.* It was enforced only by
the flight's admissibility flag. N7 now computes
`d_M(t) = lStar * |r(t) - (1-mu, 0, 0)|` in the script, reports its
minimum, where it occurs, the altitude and the margin, and gates it
against a threshold named in section 0 -- which also settles the
convention Astra asked about: 1900 km is a MOON-CENTRE distance, a
162.6 km altitude floor.

*And the new gate immediately caught its own author.* The first version of
N7 sliced `flight.Y(1:3, :)`, but `fly_transfer` returns the flight from
`pumpkyn.cr3bp.tfMinProp` as `[nTimes x 14]` -- times down the rows -- so it
measured the distance between the first three TIME SAMPLES and reported a
comfortable 332,999 km PASS at `t/t_f = 0.000`. The implausible location,
not the implausible number, is what gave it away: a lunar transfer's
closest approach is not at departure. Corrected to `flight.Y(:, 1:3)`, the
gate reads **6410.6 km from the Moon's centre at t/t_f = 0.999** (altitude
4673.2 km, margin +4510.6 km), which is the tulip end, where it belongs.
The lesson is the review's own: a gate that computes the wrong quantity
passes just as loudly as one that computes the right one, so every new
gate needs a number whose PLAUSIBILITY can be judged, not merely its
verdict. On the shipped anchor the whole script now reports 14 of 14 gates
passed, `studyOK = true`.

*Gate inputs are validated before they are compared.* A negative
"absolute" residual satisfies `< tol` and a NaN count makes `> 0` false;
either would let a broken instrument print PASS. Every scalar the
necessary block compares is now checked real, finite and non-negative
first, and the script fails closed on malformed data -- the policy
`certify_root` already followed.

**Declined, with reasons.**

*"Compute N2-N6 in the script instead of calling the library
instrument."* This reverses a documented decision recorded at the call
site: the script DID recompute all four and compare, and the comparison
measured 0.0e+00, because a copy of the same arithmetic on the same
samples is not an independent implementation. Running the instrument the
18,360 catalog entries were certified with is the stronger statement. The
pedagogical concern behind the finding is real, and the right answer is
not a tautological re-computation but showing the equations in the
narrative -- Astra's own derivation (including the `+C lambda_v` Coriolis
term in the velocity adjoint) is now quoted in
`doc/seeds_and_solution_route.tex`.

*"No conjugate time in (0, t_f]" overstates what a sampled scan
establishes.* Correct, and already said: the verdict block states in full
that this is numerical evidence and not a certificate, that positivity is
tested at sampled times with no between-sample bound, and that
application of the theorem remains conditional on the subarc-normality
argument and the free-mass/free-time reduction. The audit
`doc/mintime_second_order_audit.tex` is the home of that argument; the
script points at it. No wording change; a interval-arithmetic enclosure
is the open item, recorded here.

*Propagation-accuracy refinement studies, an FD step-size sweep, a
lift-margin error budget over trajectory and constraint refinement.*
All three are real gaps between "measured under this refinement" and
"bounded". They are studies, not fixes, and each is a half-day. Recorded
as open; the script's existing text already says the lift margin is a
measured sensitivity and not a proven rank separation.

## 75. Where the first root comes from: the cold lottery, two ladders, and the winding wall at the anchor cell (2026-09-15)

`indirect/root_origins_study.m` is the prequel to `anchor_study.m`. That script
starts from a certified root of a neighbouring transfer; this one starts from
nothing and measures the three mechanisms that make a first root: a cold solve
where the problem is easy, a ladder down the engine axis, and a branch-blind
basin hunt. Full default run: ~50 min under R2026a, log kept as
`origins_run6` in the session scratchpad, record
`indirect/results/root_origins_study.mat`.

**Provenance stated honestly at the top of the script.** The very first
DRO->tulip root on this machine was not built here: it arrived with
pumpkynPie as a converged indirect solution of exactly this cell (Darin's
walk-down), and it is the "reference" every early table is quoted against.
Everything in this script is the route for the cells that got no gift.

### The cold lottery, reproduced at the operating point

Same solver, same scheme, same t_f guess (4.0 ND); only the mesh changes.
70 mN / Isp 900 s / 150 kg, sD 0, sA 0.0754, unconstrained (as the 2026-08-03
record was). Reference: 4.0152 ND = 17.798 d.

| N | t_f (ND) | vs reference | periselene | status |
|---|---|---|---|---|
| 400 | 4.6809 | +16.6% | 6563 km | converged |
| 800 | 4.9909 | +24.3% | 2243 km | Maximum_Iterations_Exceeded |
| 1600 | 6.4126 | +59.7% | 1995 km | converged |

Spread 37.0%, two distinct basins at a 2% threshold, none of them the
campaign's root. This is the same phenomenon as the 2026-08-03 table
(4.38 / 4.45 / 4.78) with a wider spread; the conclusion is unchanged and now
reproducible from one script.

### The direct ladder, and a basin jump that looked perfect

15 -> 0.5 N at Isp 1710 s, 11 rungs, N = 800. t_f grew 0.432 -> 3.358 d and
dV fell 4.2204 -> 0.9960 km/s. The inter-rung guess `tf <- tf*(Tprev/T)^0.6`
predicted every rung to **1%** (printed rung by rung).

**A coarse rung set jumps basins, and nothing in the solver output says so.**
The first version of the ladder used [15 10 5 2 1 0.5] (ratios down to 0.4).
The 1 -> 0.5 N step returned a converged solution with defect 3.9e-14, safe
periselene and |u| = 1: **t_f 54.655 d, dV 46.8075 km/s -- 16x its own guess,
94% of the mass burned, a different winding number.** The ladder accepted it
and the handoff then ground on a 54-day many-revolution trajectory with ode45
tolerance failures.

The fix is the campaign's own `preflight_screen` rule, now in the script: a
rung is REFUSED if its t_f leaves [1/3, 3] x the guess that produced it. Two
details matter. The band applies to a STEP only -- the top rung's guess is a
cold round number with no branch behind it, and at 15 N the true answer is 41x
below it (the first version of the screen refused the 15 N rung and killed the
run). And a basin jump is a statement about the STEP, not the problem: the
remedy is a finer ratio, never a looser gate. The shipped rung set steps by
about 0.75 and does not jump.

### The seam: harvest with lambda_t exposed

The ladder runs in plain time, as the shipped catalog did. Under Sundman the
defect system carries the time state, so the multipliers carry one more row --
lambda_t, which PMP fixes at +1 when the objective is t_f. That row is a free
check on station association, sign and scale together, and it does not exist
without Sundman. The script re-solves the bottom rung ONCE in the Sundman
chart and harvests from that: **sign vote 100.0%, lambda_t = 1.000000**,
shooting |R| 4.59e-12 in ONE iteration.

### The indirect ladder, and the wall at THIS cell

Below ~0.5 N the walk continues by multiple shooting on banked junction
states. 0.375 -> 0.12 N at Isp 1710 (7 rungs), then the Isp stage 1710 ->
1400 -> 1150 -> 900 at fixed 0.12 N (3 rungs, **4 seconds total**, t_f falling
11.119 -> 10.606 d), then 0.11 and 0.10 N both REFUSE.

| route | deepest closed | wall |
|---|---|---|
| thrust only, Isp 1710 throughout, ratios ~0.5 | 0.16 N | 0.12 N |
| thrust only, Isp 1710 throughout, ratios ~0.85 | 0.12 N | 0.105 N |
| ratios ~0.85 + Isp stage at 0.12 N | 0.12 N (Isp 900) | 0.11 N |
| *campaign, from the FASTEST 0.5 N cell* | *0.09 N* | *0.067 N* |

Two findings. **A finer rung ratio buys real depth** (0.16 -> 0.12 N), because
a jump is a step property. **The Isp stage is free and helps but does not move
this wall**: it drops t_f 4.6% at fixed thrust in 4 s (faster depletion, a
lighter vehicle, more late-arc acceleration) and the next thrust rung still
refuses. Revolutions climb 0.81 -> 1.31 across the walk: this is a WINDING
wall, and **its depth is a property of the CELL**. The anchor phase walls near
0.11 N; the campaign reached 0.09 N because it started the deep walk from the
fastest 0.5 N entry of the sheet. Where you start the deep walk is a choice,
not a detail.

The root at 0.12 N / Isp 900 s **certifies**: t_f 10.6060 d, flown miss
0.000 km / 0.000 m/s, conjugate PASS, lift margin 32542x. B2 (target engine
reached) is reported FAIL and does NOT throw -- a wall is a finding, and
section 6 prints what was and was not demonstrated.

### The basin hunt, demonstrated live

From the 0.12 N root at sA 0.0754, warm-starting the DIRECT solver at three
other arrival phases:

| sA | t_f | vs source | reading |
|---|---|---|---|
| 0.3254 | 11.240 d | +6.0% | another basin, slower |
| 0.5754 | 9.081 d | **-14.4%** | **another basin, and FASTER** |
| 0.7837 | 10.411 d | -1.8% | same basin, continued |

Two of three landed in another basin and one came back 14% faster than the
family it was seeded from. That is the move that found fast2, direct18 and
direct11 (sections 61, 63, 68), reproduced from a cold start in one run.

### Two latent bugs found on the way

- **certify_root could not run without the Parallel Computing Toolbox.** It
  built its pool as `d('pool', gcp('nocreate'))`; MATLAB evaluates arguments
  eagerly, so gcp ran even when the caller had already passed a pool -- and
  gcp THROWS when PCT is absent or its licence is held by another MATLAB
  session on this machine. The desktop session is enough to take the seat,
  and then every `matlab -batch` job here loses certification at its first
  call. `verify_with_pumpkyn` had documented the eager half and still called
  gcp unguarded on the other branch. New `costate_common/current_pool` returns
  the open pool or [] without throwing; certify_root, verify_with_pumpkyn and
  second_order_pass route through it, certify_root honours an explicit
  `.pool = []` (which is what capped_pool returns and which means UNFENCED),
  and capped_pool no longer dies if a pool cannot be opened. Note for the
  record: **R2025b has no Parallel Computing licence here; R2026a does**, and
  the campaign launchers already point at R2026a.
- **fine_sheet_job's candidate counter shadowed its certified-column count**,
  so the verdict line reported the last column's candidate count. Fixed with
  the loop's own names (and the loop variable renamed off `j`).

## 76. root_origins_study reviewed by GPT-6 Astra: a wrong sentence about delta-V, an unchecked throttle, a table that mixed altitude with radius, and a "basin test" that was not one (2026-09-16)

`reviews/root_origins_study_astra_2026-09-16.md` (xhigh, 411 s, $1.18; the
prompt carried the script plus the contracts of every helper it calls, so the
reviewer could verify rather than speculate). Every finding was checked
against the code before anything changed. What follows is the triage.

### Verified wrong in the script, fixed

- **"dV rises slowly while t_f rises fast" was false, and the script's own
  table said so.** Ideal dV fell 4.2204 -> 0.9960 km/s down the direct
  ladder. At fixed exhaust speed the all-burn propellant fraction is
  `T*t_f/c`, and if `t_f ~ T^-0.6` then `T*t_f ~ T^0.4` FALLS as thrust
  falls. The script now prints both directions as measurements with that
  argument. (I wrote a slogan and did not read my own numbers.)
- **The throttle was never checked.** The direct solver leaves it FREE on
  purpose (its header: pinning it would be simpler; checking that it
  saturates is a genuine test of the formulation) and the script checked
  only `maxUnit`, the direction norm. Astra proposed `thrLock = true`; the
  campaign's own design says check, not pin. `directOK` -- one feasibility
  predicate for every direct solve in the script -- now gates on `thrMin`
  (measured `1 - u_min` 5.8e-7 and 3.6e-8 on the smoke rungs), and on
  `maxInterp` and `tfSpread`, which the solver computes separately from
  `maxDefect` and which no section had looked at.
- **The recorded lottery table mixed altitude with radius.** FINDINGS'
  2026-08-03 column is "min node alt"; I quoted -343 km as a "periselene"
  beside a live column that prints Moon-centre RADIUS. A radius cannot be
  negative. Script and document now convert (+1737.4 km) and label both.
- **Section 0 settings never reached certification.** `certify_root` reads
  `m0kg`, `gateKm`, `gateVms` and `moonKmMin` from ITS OWN options
  (defaults 150 / 100 / 10 / 1900) and not from the problem struct;
  editing `op.m0kg` or `num.clearKm` would silently have left the
  certificate on defaults. Forwarded, and printed on the certify line.
- **The hunt's "same basin / ANOTHER BASIN" labels were unsupported.** Each
  probe is a different boundary-value problem (the arrival state moved), so
  its t_f is not comparable to the source's as a basin test: one family's
  t_f varies with phase by far more than 2%, and two families can share a
  t_f. The hunt now reports a CANDIDATE ("X% faster/slower than the source
  at its own phase"), banks each converged probe as a harvested seed, and
  says what turns a candidate into a verdict: a same-phase baseline from
  `arclength_arrival` and a comparison of states and costates there --
  which is what sections 61, 63 and 68 actually did.
- **`uniquetol` did not implement the advertised relative test** (its
  default scales by the largest element, so a long outlier regroups the
  others). Now an explicit absolute tolerance, `tol.tfCluster` x the
  shortest converged t_f, and called a flight-time CLUSTER count, which is
  what it is.
- **The depletion guard's comment was wrong.** `m(t_f) = 1 - T t_f/c`, so
  exhaustion is at `t_f = c/T`; the `0.6 c/T` cut keeps 40% of the mass. It
  is a seed-propellant POLICY (`iladder.maxPropFrac`, recorded when it
  skips), not a physical bound.
- **The flown-miss gate was fail-open on NaN**, had no velocity gate, and
  did not check that the witness reached t_f with positive mass. Fixed;
  `gateVms` added.
- **Refusal reasons in the wrong order**: a band failure was reported before
  a convergence failure, so an unconverged iterate could be labelled a
  "basin jump". Precedence fixed; the band is now called what it is, a
  plausibility heuristic applied to steps only, not a branch detector.
- **An Isp-only rung that failed re-ran the identical solve five times**
  (the exponent has no effect on an Isp step). One attempt.
- **`revs` was total swept azimuth, not a winding number.** Renamed
  `turns`, defined on the printed header (reversals add).
- **The indirect ladder re-flew the whole arc from lambda(0)** to build the
  next seed, bringing back the amplification multiple shooting exists to
  avoid, and the banked junctions were never used. New `flyFromJunctions`
  rebuilds the trajectory segment by segment from `it.Y`; the hunt's warm
  start uses it too, with duplicate time samples removed and a zero-primer
  check, inside the try.
- **The pool was created AFTER the handoff shoot**, so that shoot and every
  full-arc propagation ran unfenced even with a pool. Pool first; the
  handoff, the reconstructions and the witness flights are fenced; the
  per-rung budget is passed into the fence caps.
- **The polished root was not adopted.** `certify_root` polishes; the
  certificate belongs to `C.z`/`C.Y`, but the hunt and the record used the
  pre-polish candidate. On `C.ok` the polished root is adopted and marked;
  the certificate, the options it ran with, the gate table and the three
  outcomes are saved. `saveq` reports failure instead of swallowing it.
- **"CERTIFIED root" was printed unconditionally** on the short-ladder
  branch. Section 6 now branches on `C.ok` first.
- **Two hardcoded tolerances** (seam 1e-6, solver target 1e-11) made "one
  source" false. Moved into `tol` with solve targets and acceptance
  thresholds kept apart and said to be different things.
- **"Direct ladder finished" tested thrust alone** and its label differed
  between the print and the gate table. Now the accepted-row index.
- Refused rungs and thrown probes were absent from the records; "of N
  probes" meant returned records. Every attempt is recorded with a reason
  taxonomy (timeout / unconverged / flown-miss / seed-policy skip / bad
  flight), and A, B and C each print an OUTCOME (SUPPORTED / PARTIAL /
  CANDIDATES / NOT OBSERVED / INCONCLUSIVE) separate from the machinery
  gates.
- Prose softened where the computation did not earn it: "winding wall, not
  a sensitivity one" -> "this search policy stalled on this cell"; "never
  loaded" -> the orbit getter refines a catalogued family seed, no TRANSFER
  is loaded; "near-impulsive, sub-revolution" -> the measured half-day.

### Verified NOT a numerical defect (the review's #2)

Astra ranked "the certification seed violates the junction-array contract"
second: `ms_tfmin` returns `info.Y` as 14 x K junction STARTS while the seed
contract documents 14 x (K+1). Measured: `reached.Y` is [14 24], `tGrid`
[1 25]. Then read the consumer: `ms_bvp` packs its unknowns from
`seed.Y(:,1)` and `seed.Y(:,2:K)` and never reads column K+1, and
`certify_root` touches only `seed.tGrid`. The K starts plus the fixed
departure state ARE the complete parameterisation; the endpoint column is
redundant. So feeding `it.Y` back is lossless and the certificates stand.
What IS wrong is the DOCUMENTED contract (`ms_bvp` header: ".Y [14 x K+1]")
against what `info.Y` returns, and my comment "full K+1 junction states".
Comment fixed; the library header mismatch is an open item, and anything
that indexes `.Y(:, end)` on an ms info output is reading the K-th START,
not the endpoint.

### Disagreed, with reasons

- `thrLock = true`: see above -- check, do not pin.
- "Certify periselene between nodes": the per-rung column is a node screen
  and now says so; the CERTIFIED root's clearance comes from
  `validate_flight` on the dense flight inside `certify_root`, which is the
  continuous check the reviewer asked for, just not per rung.

### Smoke after the changes (R2026a, reduced settings)

Lottery N = 200/400: 4.5510 / 4.6809 ND, 2 clusters, 0 within 1% of the
reference -> A SUPPORTED. Direct ladder 3/3 accepted, throttle saturated to
5.8e-7 / 3.6e-8. Sundman re-solve within 6.1e-7 of the plain-time rung's
t_f; floor slack 4514 km; H2 100%, H3 lambda_t 1.000000 with the mapping's
own check agreeing, H1 |R| 6.3e-13. Indirect 2/2. Certified at 0.375 N /
1400 s from the FORWARDED options. Hunt: 1 candidate, 13.7% slower, banked.
Full default run relaunched to refresh the record.

**Addendum, first full run after the review.** The new throttle gate at
`1 - u_min < 1e-6` refused the 15 N rung at N = 800 with `1 - u_min =
1.30e-6` (N = 400 had given 6.4e-7 on the same rung; 5.8e-7 and 3.6e-8 lower
down). That is interior-point bound slack at IPOPT's tol 1e-7, and it grows
with problem size and wherever the switching function is small -- not a
throttle dip, which would take u toward 0. The campaign's own
`certify_dro_mintime` carries the same check as G6*, ADVISORY at 1 - 1e-6 and
excluded from passAll, for the reason its comment gives: full throttle is an
empirical property of the extremals found here, not a theorem. Gate set to
1e-3 with that rationale in the tolerance block: above any barrier slack seen,
three orders below any dip that would change the problem the shooting solver
solves next. Rerun launched.

**Full default run of the reviewed script (2026-09-16, R2026a, ~55 min).**
Every number the pre-review run produced is reproduced: lottery 4.6809 /
4.9909 (iterate, excluded) / 6.4126 ND; direct ladder 11/11 accepted with
`1 - u_min` under the 1e-3 gate on every rung; handoff sign vote 100%,
lambda_t 1.000000 with the mapping's own check agreeing, |R| 4.59e-12 in one
iteration; indirect ladder 10 of 12 rungs accepted, 0.375 -> 0.12 N then the
Isp stage, then 0.11 and 0.10 N refused. The new per-rung taxonomy says WHAT
refused them: `unconverged 5` on each -- every exponent's shooting solve
returned without converging, none timed out, none was skipped by the seed
policy. So the stall is the shooting solver failing to converge from these
seeds, not a budget artefact. Certified at 0.12 N / 900 s on the forwarded
options: 10.6060 d, flown miss 0.000 km / 0.000 m/s, conjugate PASS. Hunt:
three converged candidates, two faster than the source (14.4% at sA 0.5754,
1.8% at 0.7837), banked as seeds; basin identity not claimed. Outcomes: A
SUPPORTED, B PARTIAL (0.12 N of 0.07 N), C CANDIDATES.

## 77. root_origins_study, Astra round 2: the fixes verified, one of them wrong, and the two disputes adjudicated (2026-09-16)

`reviews/root_origins_study_astra_round2_2026-09-16.md` (xhigh, 494 s,
$1.53). The bundle carried the current script, the round-1 review verbatim,
section 76's triage, the full default run's output and the complete bodies
of every helper the script calls, so the reviewer could verify the fixes
rather than re-review from scratch. Its verdict on the round-1 list: 14 FIX
OK, 9 FIX PARTIAL, **1 FIX WRONG**. Every item was checked against the code
before anything changed.

### The wrong fix, and it was mine

**Adopting the certifier's polished root with the pre-polish time grid.**
`certify_root` polishes with `ms_tfmin` and returns `C.z` and `C.Y`
(junction starts) but no grid; `ms_bvp` keeps the NORMALIZED breakpoints
fixed and re-solves t_f, so the junction TIMES belong to the polished t_f.
I paired `C.Y` with the old `itCur.tGrid`. The certificate was valid; the
record and the hunt's warm start were internally inconsistent by the
polish's change in t_f on every segment -- dormant in the default run only
because the pre-polish residual was already below the polish target. Fixed:
`itCur.tGrid = sig*C.z(8)` with the normalized grid preserved from the walk,
shapes asserted, and the relative change printed (it will read ~1e-12 when
the walk's own root was already polished, which is the honest number).

### The two disputes

- **Junction-array shape.** Astra withdrew its round-1 numerical claim
  after reading `ms_bvp`: no consumer reads column K+1 or infers K from
  `size(Y, 2)`; a 14 x K starts array plus K+1 times is lossless. It named
  the library headers as what to correct and said not to append an
  endpoint merely to satisfy an obsolete header. Agreed; still an open
  library item.
- **Throttle.** Astra agreed pinning is not mandatory and supplied the PMP
  structure that makes saturation expected here: with free terminal mass
  the mass costate satisfies lam_m(t_f) = 0 and lam_m' = -T q |lam_v|/m^2,
  so lam_m >= 0 and the throttle coefficient in H is strictly negative away
  from a zero primer -- q = 1. On the 1e-3 gate: defensible as an explicitly
  approximate SEED-COMPATIBILITY screen backed by the all-burn shoot and
  the certificate; not a proof of saturation. What it cannot catch: a
  sustained 0.1% under-throttle, overshoot above one, deviation between
  stations, malformed values hidden by a min. Applied: the throttle is now
  checked on BOTH sides (min below 1 - tol and max above 1 + tol refuse)
  over nodes and midpoints, and the history must be finite. Its suggested
  time-weighted integrated deficit is recorded as a better statistic for
  the certifier, not built here.

### The partial fixes, completed

- **Poolless policy.** The script's `fenced` degrades without a pool, but
  `certify_root` asserts unless `.allowUnfenced`; a poolless run would have
  aborted at section 6 after an hour. Now `num.allowUnfenced` (default
  false) is an explicit section-0 policy: without a pool the ladder and hunt
  run unfenced and say so, and certification is recorded as "not certified:
  no parallel pool" instead of throwing.
- **`run_capped` returns ok = false for a timeout OR a worker error.** I
  had labelled that "timeout". Now "fence: timed out or worker errored".
- **One deadline per rung.** Remaining budget is recomputed before the
  shoot and before the witness; no stage starts with under 10 s left; the
  +90 s cancellation grace is separate from the work budget.
- **No re-propagation on the ladder.** K and the normalized grid are the
  same rung to rung, so the banked starts ARE the next seed: grid scaled to
  the guess, mass row rebuilt. The 120 s fenced reconstruction per rung is
  gone, and so is the interpolation that perturbed converged junctions.
  `flyFromJunctions` stays for the hunt, which needs dense samples; its
  seam rule now keeps the BANKED start of the next segment (so a resample
  at a junction recovers it exactly), validates its inputs and segments,
  and returns the worst seam mismatch in km rather than hiding it.
- **`directOK` fails closed.** Required diagnostics default to NaN, not 0;
  the state array must be finite; the throttle history must be finite.
- **H3 has three states.** NOT APPLICABLE in plain time (kept out of the
  pass count), PASS/FAIL in the Sundman chart where a missing or non-finite
  lambda_t, or a mapping flag that is not true, FAILS. Missing evidence is
  never read as agreement.
- **The hunt tracks solver-converged, accepted and harvested apart.** An
  accepted direct solution whose multipliers did not come back keeps its
  numbers (X, U, t_f) so "anchor it" stays possible; the summary prints all
  four counts; a clearance-only refusal has its own reason.
- **Per-exponent attempt records** (guess, outcome, residual, misses,
  seconds) replace the aggregate counters; a rung that ran out of budget
  is recorded as that.
- **The reference carries its problem.** `ref` holds the orbits, cell,
  engine and mass it belongs to; every "vs reference" line is silent
  unless the live problem matches.
- **B's outcome is decided after certification**: SUPPORTED means a
  certified root at the target engine; "reached numerically" is its own
  state; the initial engine is printed from the table, not hardcoded.
- The Sundman re-solve now ENFORCES the node-radius agreement its comment
  claimed; the Sundman claim about lambda_t says sign and scale, not
  station association; the record is value-only (no pool handle) and
  carries the whole configuration and tolerances; the first rung's throttle
  slack and the worst accepted-rung slack are printed; misses are printed in
  scientific notation where "0.0 km" hid the resolution; the header no
  longer says every mesh converges, "revolutions" or "winding wall".

### Not applied

Astra's time-weighted throttle deficit and complementarity comparison
belong in the certifier, not a study script. Its per-exponent record is
in; a separate per-stage timing breakdown is not.

**Smoke after round 2 (R2026a, reduced settings):** lottery 2/2, 2 clusters,
A SUPPORTED; direct ladder 3/3 with the first rung's slack now printed
(6.4e-7) and the worst slack 6.4e-7 against the 1e-3 gate; Sundman re-solve
within 6.1e-7 in t_f and 4 km in node radius; H2/H3 PASS with the mapping's
flag true; H1 6.3e-13; indirect 2/2 with the worst accepted flown miss
8.4e-8 km printed in full; certified, polish moved t_f by 0 (the walk's root
was already at the polish target) and the grid rescaled with it; hunt: 1
attempted / 1 solver-converged / 1 accepted / 1 harvested, seam mismatch
5.2e-9 km on the re-flown source. Full default run relaunched.

**Full default run of the round-2 script (2026-09-16, R2026a, ~55 min).**
Every rung of the previous full run reproduces to the printed digits with
the no-propagation seeding (the banked starts, grid scaled, mass row
rebuilt): direct 11/11, worst accepted throttle slack 1.3e-6 against the
1e-3 gate (now printed, first rung included); handoff sign vote 100%,
lambda_t 1.000000 with the mapping's flag true, |R| 4.59e-12 in one
iteration; indirect 10 of 12 accepted; the two refusals at 0.11 and 0.10 N
now carry their per-exponent record -- all five exponents `unconverged` on
each, no fence event, no seed-policy skip. Worst accepted flown miss
1.78e-4 km / 1.76e-5 m/s, printed in full. Certified at 0.12 N / 900 s:
10.6060 d, flown miss 9.64e-5 km, conjugate PASS; polish moved t_f by 0 and
the grid was rescaled with it. Hunt: source re-flown from its junctions
with a worst seam mismatch of 2.7e-8 km; 3 attempted, 2 solver-converged,
2 accepted, 2 harvested, 1 faster (14.4% at sA 0.5754). The third probe
(sA 0.7837), which converged in 61 s to a -1.8% candidate on the previous
run, hit its 300 s CPU cap this time at a 55-day iterate; it is recorded as
"solver did not converge", which is the only thing a capped probe can be
said to show -- a lottery ticket lands on either side of its budget from
run to run, and the script no longer reads that as anything about basins.
Outcomes: A SUPPORTED, B PARTIAL (certified at 0.12 N), C CANDIDATES.

**The guide reviewed (2026-09-16).** `doc/root_origins_study_guide.tex` went
to GPT-6 Astra (xhigh, 578 s, $1.63) with the script, the full-run log and
these sections inlined; review at
`reviews/root_origins_study_guide_astra_2026-09-16.md`. It verified as RIGHT
the Hamiltonian signs, free terminal mass giving lam_m(t_f) = 0, the mass
adjoint and throttle coefficient, lam_t = +1 under the normal Mayer
convention, the 14K - 6 = 330 shooting count, the ND conversions, the falling
propellant AND falling ideal delta-V (d dV/da = c/(1-a) > 0), the next-seed
construction, the polished-grid formula, and every lottery, ladder, handoff
and hunt number bar two. Wrong, and fixed: the Hermite-Simpson equations
omitted the lifted-time scaling (F = T f); the turns formula omitted the
unwrap; **the handoff box quoted the SMOKE run (6e-7, 4 km) as the default
run (3.80e-8, 6414 -> 6414 km)**; "predicts to 1%" was 2.1% on the last
step (0.7416 -> 0.7575 ND), which this section's own "every rung to 1%"
also got wrong; "four routes on this cell" counted a route from another
cell; the forwarding fix was attributed to round 2 instead of round 1. The
diagram drew the lottery as a side branch and hid the plain-time fallback,
the last-accepted-rung continuation and the early fatal asserts; it is
redrawn in execution order. The all-burn argument now carries its
assumptions and the proof that primer zeros are isolated (vanishing
(lam_r, lam_v) forces H = 1). Provenance was moved out of the explanatory
boxes into its own section. **Family count reconciled:** the script header
and the seeds document said four of five were found by branch-blind direct
solves; this record names three (fast2 61, direct18, direct11 68), and all
three sources now say three.

## 78. The library generator reviewed as a whole: two Astra passes, two host forks, and what "certified" currently means (2026-09-17)

Mike asked for the DRO -> tulip costate-library generation code to go to
GPT-6 Astra at xhigh with an independent host review, and for thought on
three goals: sufficient tests for local optimality, a generator for any
orbit pair, and awareness of the cart-pole work in
`optimal_control_examples`. Two Astra passes ran through the raw API
(`run_phase_torus` with every callee interface, 228 KB, $1.86, 7 min; the
pipeline + gate stack + audit document, 439 KB, $2.88, 9 min) and two host
forks read the gate stack and the pipeline independently. Every adjudicated
claim was verified at its cited lines; the MATLAB semantics in dispute
(row orientation of the sheet fields, the launcher's `exec`, `pgrep`
self-match under `zsh -c`) were checked in the shared session.

Records: `reviews/run_phase_torus_review2_adjudicated_2026-09-17.md`
(19 driver findings, P0-P2) and
`reviews/library_pipeline_and_optimality_adjudicated_2026-09-17.md`
(21 pipeline findings, the Q2 gap table, the Q3 abstraction). Astra
transcripts beside them.

**A mis-bundled first launch.** The first driver bundle inlined
`DRO_tulip/run_costate_library.m`, the August thrust-ladder script, instead
of `DRO_tulip/indirect/run_costate_library.m`, the struct-driven front door
the driver calls. Two functions with one name, one folder apart; the
driver's `addpath` order makes MATLAB pick the right one, a reader will
not. Killed at 5 min, relaunched. Rename or delete the August script.

**What the passes converged on, independently.** (1) The sheet, ribs,
filler and discovery are all set up at the shipped operating point:
`run_costate_library` passes the sheet builder no `sD`, `anchorMat` or
`sA0`, so every campaign so far has worked because every campaign used
`sD0 = 0` and the 70 mN anchor. (2) Root identity is t_f alone in the
registry (1e-3 d), the family map (0.02 d) and the catalog winner
(1e-9 d); only the sheet compares z8. (3) The optimality certificate is
for the FIXED-PHASE point-to-point problem with SAMPLED H2/H3 positivity
and an ASSESSED (not proved) H5; a free-phase saddle passes by
construction, and the theorem the audit invokes is not cited. (4) The
catalog audit fails open: a polish that times out leaves `conjNow = NaN`
and every rejection is `isfinite(x) && ...`. The 576/0 figures of sections
70-72 stand as measurements of what the audit checked, not of what its
name says; re-run fail-closed before quoting them again.

**Driver-only defects, host-found and Astra-confirmed:** a failed arc can
never be retried (`.fail` tested before `.done`, never cleared); a
`budget` campaign cannot be resumed although its stop reason says to;
recorded PIDs are never read on resume (duplicate writers); discovery can
re-anchor a root the filler just anchored (`added` is required on the
spine path, not the discovery path); a resumed round drops the filler's
roots (`fh.nCert` counts this call only).

**Goal 2 (local-optimality tests), the answer.** Not sufficient for the
claim, sufficient for the assessment. Cheap and next: Lipschitz margins
for H2/H3 above the measured lambda_m uncertainty; the phase-transversality
cross-check `dt_f/ds_A = tau_A lambda_rv(t_f) . f_orb(x_f)` against the
sheet's own differences (zero cost, and the first check of lambda(t_f)
against anything); the fail-closed mutation suite through audit and
certifier; `multiplicity == 0` in the certifier's consistency test. Then
the short-time sign, an independent accessory-problem inertia test, and
one paragraph citing the mixed bang/smooth sufficiency theorem. Astra's
exact gap identity `H(s,alpha) - H(1,alpha*) = T Q_mt (1-s) + (s T
|lambda_v| / 2m)|alpha - alpha*|^2` is the quantitative form of H2/H3.

**Goal 3 (any orbit pair), the answer.** The thrust-ladder pipeline is
already pair-generic (`ladder_endpoints` / `get_family_orbit`); the
phase-torus pipeline is wired to `arclength_arrival`'s literal
`tauDRO/NpTulip/pmTulip` and every identity check downstream copies that
list. One `pair` struct (`dep`/`arr` with family, params, `state(s)`,
`dstate(s)`, period, `kind`, `wrap`; `model`; one identity string) as
`B.problem` is the refactor, and `physicsOnly` setup is its first step.
GTO is departure-only under this abstraction: as an arrival it is an
epoch-dependent rendezvous with a target-motion term in the free-time
condition, a different problem class. Four documented families are not
in `get_family_orbit` (Pumpkin, LPO, Axial, Cycler). The z8-only schema
needs junction states before any long-transfer pair.

**Goal 4 (the cart-pole work).** Its lessons are the same ones this review
found: an audit that re-runs the same instruments is a mirror, not an
oracle (the sign error that nine reviews missed); a gate that cannot fail
is carried as narrative, not as a gate (N6); every study script asserts its
inline numbers against the shared instrument. The phase-transversality
check is this pipeline's power-balance oracle: independent of the
certifier, derived from the problem's geometry, free.

**Goal 5 (replicate the study scripts under DPO_tulip).** Do the `pair`
refactor first, then instantiate DPO as its second consumer, the same
admission rule oclib uses; a copy made today carries the DRO literals and
the shipped-anchor dependency with it.

## 79. The P0 repairs, and one script that rebuilds the library (2026-09-18)

**The P0 items of section 78, fixed test-first.** Each was given a failing
check before the code changed; the front-door check reproduced the headline
defect live (a campaign at `sD(1) = 0.2` built its sheet at departure phase
0 and died on "the sheet was certified at departure phase 0.000000").

| # | defect | repair | where |
|---|---|---|---|
| 1 | the sheet, ribs, filler and discovery were set up at the SHIPPED operating point | the front door takes `.anchorMat`/`.anchorSA` and passes `sD(1)`, the anchor and its phase to the sheet; `commonOpts` carries them from the driver; ribs, filler and discovery ask for the closures only (`physicsOnly`) and no longer re-polish the 70 mN anchor | `run_costate_library` stage 1, `run_phase_torus>commonOpts, physicsOpts`, `build_ribs>setupFromSheet`, `fill_holes_direct>physicsFromCatalog` |
| 2 | a failed arc could never be retried; a MATLAB that died at launch left no verdict | `spawnArc` clears stale `.done/.fail/.pid`; each arc job runs under a small shell wrapper (`jobs/arc_*.sh`) that records MATLAB's pid and writes the `.fail` if MATLAB exits with no verdict | `run_phase_torus>spawnArc` |
| 3 | a `budget` campaign refused to resume | `campaignEnded`: `budget` ends the call only while `.maxRounds <= roundsDone` | `run_phase_torus>campaignEnded` |
| 4 | a resumed round dropped the filler's roots (`fh.nCert` counts this call only) | the round's own holes file is offered to the first package (`roundExtras`); spine roots are registered from the FILE | `run_phase_torus` round loop |
| 5 | recorded PIDs were never read on resume (two writers on one arc) | `adoptLiveJobs`: a live recorded job is waited on, not respawned; `killJobs` stops wrapper and MATLAB | `run_phase_torus>runArcs` |
| 6 | discovery could re-anchor a root the filler had just anchored | ONE rule for both paths, `promotionVerdict`: newly registered, beats the spine by `acceptDays`, on no known family; filler anchors are now counted | `run_phase_torus>promotionVerdict` |
| 7 | `packaged` stood in for audit and sweep success; the audit failed OPEN | `out.stages` (package/audit/sweep: true, false, NaN = not run), required by `assertPackaged` and by the generated finalizer; the audit treats a failed or timed-out re-polish or gates call, a moved root (`tolMove` 1e-6) and any non-pass of H2/H3/dim S/H6 as a BAD row | `run_costate_library>stageOutcomes`, `audit_phase_catalog` |

Tests (`indirect/tests/`): `test_run_phase_torus_p0` 20, `test_run_costate_library_seams` 10,
`test_fill_holes_physics` 5, `test_audit_fail_closed` 4, all green; the local
helpers are reached through a `localfunctions` seam (`run_phase_torus('localfunctions')`
and the same on the front door, the filler and the rib builder). Regression:
`test_phase_lists`, `test_rib_from_crossing`, `test_sheet_from_arcs`,
`test_family_map`, `test_crossings_from_arc` pass; an untampered audit of
entry 1 of the library of record is still a clean row (polished root
1.3e-10 from the stored one). Not yet exercised: a live round under the new
wiring (the 3 x 3 acceptance torus should be re-run), and the fail-closed
audit over all 576 entries.

A lesson from the patching itself: a substring replace of `capped(pool,`
also renamed the fence's own definition line and its `run_capped` call.
The test caught it; a diff against git, read line by line, confirmed the
repair. Bulk renames take a word boundary.

**One script rebuilds the library: `indirect/reproduce_library_70mN.m`.**
`run_phase_torus` is the engine but was not the reproduction: the record
was built by hand over eight rounds (section 10 of the runbook), the
driver's example spec has one anchor, and nothing compared a rebuild with
the record. The script is the build as one chain, in six numbered sections:
switches; the problem and the 24 x 24 grid; the FIVE families (anchor,
cell11, fast2, direct18, direct11 -- the chain script's table); the ten
walked arcs and the seven direct-found seed roots, adopted into the
campaign folder (`adopt_walked_arcs`; nothing is ever overwritten);
`run_phase_torus` with discovery off; and `compare_phase_catalogs`, cell by
cell against the record (coverage, t_f to 1e-6 d, z8 to 1e-6 relative,
families compared as a PARTITION so index renumbering is not a difference).
With no arguments it PLANS and creates nothing; `struct('go', true)` builds
(about 6-10 h with adopted arcs, about 30 h with `.adoptArcs = false`);
`struct('compareOnly', true)` re-runs the comparison on a finished folder.
`test_reproduce_library` (16 checks) pins the comparison on the record
against itself, on four planted differences, on relabelled families and on
a foreign grid, and the adoption and plan mode. The build itself has not
been run yet.

## 80. Track A: the costates tested against something that is not the costates; H2/H3 with a margin; two certifier holes closed; the theorem read (2026-09-18)

Section 78 asked whether the local-optimality tests are sufficient and
listed what would close the gap cheaply. Four items, each test-first.

### X3, the phase-transversality cross-check

In the library's convention (`H = 1 + lam.f`, normal chart `lam.f = -1`) the
costate is the gradient of the cost-to-go, so for endpoints sliding along
their orbits

    dT/ds_D = + lam_rv(0)   . x_D'(s_D)
    dT/ds_A = - lam_rv(t_f) . x_A'(s_A)

(`costate_common/phase_sensitivity`; the mass costate does not enter). The
sign and the transposes are pinned by an oracle with no orbit in it: the
minimum-time single integrator, `T = |x_f - x_0|`, where the formula matches
a numerical derivative of `T` to 1e-8 (`tests/test_phase_sensitivity`).

`indirect/phase_transversality_check` applies it to a catalog. On the 70 mN
library of record:

- **The exact test.** Re-solve the transfer at `s +/- 2e-4` in each phase
  and difference the two flight times; no costate is read. Twelve sampled
  derivatives (six cells, both phases): relative error 4e-8 to 9e-4, verdict
  PASS. This is the first check of `lam(t_f)` against anything but the solve
  that produced it. The verdict is three-valued -- a re-solve that does not
  converge is UNRESOLVED, not FAIL (two did not at the first attempt; a retry
  at half the offset with a longer budget resolved both).
- **The grid is the limit in arrival phase, not the costates.** Against the
  sheet's own central differences the arrival sensitivity looked 20% off;
  against the exact re-solve it is right to 5 digits. At 24 arrival phases
  the flight time is under-resolved along `s_A` (a 7-petal tulip): trapezoid
  edge residuals have a median of 311 min, against 0.5 min along `s_D`. The
  arrival edge map is returned and not judged.
- **Family labels hide jumps.** Along departure phase, 43 of 576 edges have
  a residual above 10 min -- the two cells are not neighbours on one smooth
  branch -- and **33 of the 43 sit inside one family label**. Column 15
  (sA 0.6587) jumps by 4 to 7 days between adjacent departure phases. This
  is the measured form of section 78's "identity by t_f alone": the family
  map attaches by flight time within 0.02 d and cannot see these.
  `costate_common/phase_edge_residuals` is the instrument; it costs nothing.
- **No entry is a free-phase optimum.** The fastest entry, (1,3), 16.226 d,
  has `(dT/ds_D, dT/ds_A) = (-0.92, -2.05)` d per unit phase; the entry
  nearest stationarity in both phases, (16,3), 16.600 d, still has
  (-0.60, -0.21). The orbit-to-orbit minimum lies between grid points, and
  the sensitivities say which way. A free-phase polish (two more unknowns,
  the two transversality equations) is the natural next step.

Result file: `indirect/results/phase_transversality_70mN_24x24.mat`.

### H2 and H3 with a margin, over the whole arc

The gate was `min_k |lam_v(t_k)| > 0` and `min_k Q_mt(t_k) > 0`: it passed
1e-300, passed values below the flight's own lambda_m error, and said
nothing between samples (FINDINGS 30 speaks of a 1e-6 gate the certifier
never had). One slope bound now covers both: in the CR3BP
`|d|lam_v|/dt| <= |lam_r| + 2|lam_v|`, and on an all-burn arc the two mass
terms of `dQ_mt/dt` cancel, leaving `dQ_mt/dt = (d|lam_v|/dt)/m` exactly.
`costate_common/between_sample_bound` turns samples plus a slope bound into
a lower bound over the interval (oracle: `1.2 + cos 3t`, true minimum 0.2
between samples; the coarse sampled minimum 0.2086 overstates it, the bound
does not). `mintime_hypothesis_gates` reports `.minLamVBound .minQmtBound
.dtMax`; `certify_root` and `audit_phase_catalog` gate the BOUNDS above
`hypFloor` = 1e-5. Measured on all 576 entries (about 7000 samples per
flight, steps <= 0.0035): the bound is within 1% of the sampled minimum,
smallest bounds 0.2785 and 0.3126, so no entry moves. It is an estimate on
a finely sampled flight, not a validated enclosure.

### Two holes in the certifier, found by mutation and shown live

Injected through the existing harness (`test_certify_enforcement`):
`multiplicity = 1` with every other dense-scan count zero returned
"certified"; and ANY test-seam override returned "certified", with nothing
in the certificate recording it. Now `multiplicity == 0` is part of the
scan's consistency test (it counts located zeros of corank >= 2, so it can
be nonzero only when `nZero` is; all 576 entries carry 0), and a non-empty
override makes the result DIAGNOSTIC ONLY.

### The theorem, read from the source

`doc/mintime_second_order_audit.tex` cited "BCT 2007 Thm. 3.x". The paper
was read (ESAIM COCV 13(2) 207-236). It is **Theorem 2.12** (with 1.13),
under (L) strong Legendre and (S) strong regularity = corank one on EVERY
subinterval; the free-final-time normal case is their **Test 3**, a zero of
`det(dx_1..dx_{n-1}, f)`, which IS the instrument's determinant. Three
corrections: the conclusion is "locally optimal in the C0 topology", not
"strict"; BCT assume an OPEN control set and say the bang-bang case is not
treated -- their orbit-transfer control lives on S^2, which is our problem
only after the throttle is fixed at 1; and (S) is checked here on the whole
arc only. The throttle case belongs to Osmolovskii-Maurer's theory of
controls with continuous and bang-bang components (2006, 2009, SIAM 2012;
references verified through Crossref), cited in the new subsection "What the
theorem does not cover, and what is now checked" as the framework, NOT instantiated (their continuous component is in R^m, ours on
S^2). The certificate sentence was rewritten to say what is measured.

### Still open on this track

Strong regularity on subintervals (checkable: lift-space rank on prefixes);
the short-time sign of the determinant and a bound on the uncovered prefix;
an independent accessory-problem inertia test; the free-phase second-order
condition; validated enclosures if a theorem is ever wanted. And one to
act on: the 43 departure-phase jumps should be reconciled with the family
map before the map is trusted for anchoring decisions.

## 81. The reproduction chain reviewed, rehearsed, and launched; the record passes the fail-closed audit (2026-09-18)

**Astra on the reproduction chain (xhigh, 154 KB, $1.46): NO-GO as
supplied.** Adjudication: `reviews/reproduce_library_adjudicated_2026-09-18.md`.
It confirmed the wiring (spec fields, arc adoption, the wrapper's shell
semantics, the family-partition comparison: zero differences iff the
partitions are identical) and found seven launch blockers, all fixed
test-first the same afternoon:

- **My own Track A defect.** `certify_root` added `.minLamVBound` /
  `.minQmtBound` only on the success path, so a failed and a passed
  certificate had different field sets; the sheet's seed loop and the rib
  walker append certificates raw, and MATLAB refuses that. With
  `librarySeeds = true` the full build would have thrown in its first
  sheet. Every field is now born in the initializer (`test_certify_schema`).
- The audit's new validators were weaker than the certifier's: `Inf` passed
  `> floor`, `logical(2)` is true, `logical(NaN)` THREW instead of making a
  BAD row. `realScalar` before every comparison.
- `compare_phase_catalogs` reported a MATCH for a catalog holding a NaN
  flight time (`NaN > tol` is false), reported "0 differences" when a family
  map was absent, compared libraries of different engines, and miscounted
  agreeing cells. Now `.nNonfinite`, `.familiesCompared` (required),
  problem identity asserted, `.nAgree`; costates compared without t_f in the
  norm.
- The script's verdict paired the final catalog with "the newest audit file
  anywhere"; `finalAudit` reads the last round of the driver's state and
  requires it to cover every entry.

Its findings on interrupted runs (non-transactional final publication, the
fixed-point decision using this invocation's counters, registration before
promotion, PID existence as ownership, mtime verdicts, arc identity on
adoption) are recorded as limitations, not fixed: a fresh uninterrupted
build does not exercise them.

**The 3 x 3 rehearsal (`results/torus3d_p0`, the 09-16 spec).** Fixed point
in 2 rounds, 2 anchors, 91 min; its catalog is identical to the 09-16 run's
in all 7 cells (zero deviation in t_f, costates and family). Seen working
live: the arc wrapper with a pid file per arc; the package/audit/sweep
assertion on every packaging; the round-2 sheet through the seed loop;
"1 new anchor" counted in round 1 (09-16 logged 0). It found one bug: the
packaging chain names its audit `audit_70mN.mat` whatever tag the driver
was given, so the script's lookup by the driver's tag only worked where the
two tags coincide -- as they would have in the full build, by accident. The
lookup is by pattern in the last round now. Final audit 7 ok / 0 bad.

**The library of record under the fail-closed audit: 576 ok / 0 bad.**
All 576 entries, 152 min, chunked and resumable
(`results/audit_failclosed_record_2026-09-18/`). The audit now treats a
failed or timed-out re-polish or gates call, a polished root more than 1e-6
from the stored one, a non-finite or malformed gate value, and any non-pass
of the H2/H3 whole-arc bounds (floor 1e-5), dim S = 1 or H6 as a BAD row.
Section 78's caveat -- that the record's "576 / 0" was measured by an audit
that failed open -- is closed: the count stands under the policy that cannot
skip a row.

**The full build was launched at 15:22** into
`results/reproduce_70mN_24x24` (five anchors, ten adopted arcs, seven seed
roots, discovery off). Its verdict belongs in the next section.

## 82. The library rebuilt by one script: every root of the record reproduced or beaten (2026-09-19)

`reproduce_library_70mN(struct('go', true))` ran unattended from
2026-09-18 15:22 to 09-19 15:34 (24 h 12 min) into
`indirect/results/reproduce_70mN_24x24`. One round, a fixed point
("24 certified spines, 24 walked, 0 new roots, 0 new anchors; status done").

| stage | result | time |
|---|---|---|
| arcs | 10 adopted, 0 walked; 7 seed roots in the registry | seconds |
| arrival sheet | 24 of 24 columns certified | 1 h 34 |
| ribs | 24 columns, 4 workers, none retired | ~9 h |
| package / audit / sweep | 519 of 576 cells; audit 519 / 0 (fail-closed); 0 interior crossings, worst H6 5.32x, lift 11x | ~4 h |
| holes + improve | 57 holes, 58 tried, **58 certified** | 7 h |
| re-package / audit / sweep | **576 of 576**; audit **576 ok / 0 bad** | 2 h 31 |

The 6-10 h estimate was wrong by a factor of three: the runbook's per-round
rib figures are for rounds that walk only the CHANGED columns, a rebuild
walks all 24; and the filler's 57 direct solves (about 9 min each) were not
budgeted. Quote a day, or raise `.nWorkers`.

**The comparison with the record** (`compare_phase_catalogs`,
`reproduce_result.mat`):

- coverage 576 / 576 in both; nothing non-finite;
- **565 cells hold the SAME root** (t_f within 1e-6 d, costates within 1e-6);
- **11 cells hold a different root, and in all 11 the REBUILD is faster**,
  by 1.95 to 3.93 d; the record is faster in none. Nine are rows 14-22 of
  the column at sA 0.6587 (record 22.95-23.86 d, rebuild 19.05-21.42 d):
  the "9 stuck col-15 cells" left open on 09-15, and the column where
  section 80's edge map measured jumps of 4 to 7 days between neighbours.
  The other two are (3, sA 0.2837) 19.862 -> 17.287 d and (4, sA 0.3254)
  19.758 -> 17.806 d. Fastest entry 16.226 d in both; mean t_f 18.5236 ->
  18.4678 d. The rebuild's filler, seeded from each hole's certified
  neighbours in one pass over a complete first-round catalog, found what
  eight hand-driven rounds had left;
- **82 cells differ in FAMILY LABEL only, with identical roots.** The record
  labels the ribs of the columns at sA 0.4921 and 0.7837 `direct11` and
  `direct18`; the rebuild labels them "rib unidentified" (-2); 20 cells go
  from "unattached" to "rib unidentified". That is `family_map`'s rib
  attachment by flight time within 0.03 d (sections 78, 80), now measured as
  NON-REPRODUCIBLE: the same roots, walked in one round instead of eight,
  get other labels. The family names differ too (`70mN_anchor` for `fast`):
  the driver does not pass the label table the chain script had.

So the script's strict verdict is FAIL ("the same library": 494 of 576
cells agree in everything) and the comparator's second line is the finding:
**NO WORSE than the record in any cell** -- covers it, slower nowhere,
faster in 11. For a minimum-time library that is a better library, and the
record should be superseded by it once the labels are dealt with.

Two defects found by the run itself, both in the comparator, both fixed
test-first: it demanded the two catalogs list their phases in the same
ORDER (the record lists arrival phases from the anchor's, 0.0754 ... 0.9921
0.0337; the driver sorts them), so the first verdict was "different grids"
-- cells are matched by phase now; and it could not say which catalog held
the better root (`.nNewFaster .nNewSlower .noWorse`).

Open: (1) decide whether `final/` of this run becomes the library of
record; (2) make family labels reproducible (root identity with costates,
lineage by root ID) or stop comparing them; (3) pass the family label table
through the driver; (4) the 11 improved cells deserve a look in the
second-order sweep written into the new catalog (it ran: 0 interior
crossings over all 576).

## 83. The rebuilt library adopted as the library of record (2026-09-19)

Mike's decision, 2026-09-19. `indirect/results/library_70mN_24x24_final/` now
holds the library rebuilt by `reproduce_library_70mN` (section 82): catalog,
receipt, fail-closed audit (576 / 0), sweep sidecar, sheet, direct-cell file,
pictures, the comparison with its predecessor, and a provenance README. The
hand-built library of 2026-09-15 is kept whole beside it as
`library_70mN_24x24_handbuilt_2026-09-15/`; nothing was deleted. The path is
unchanged, so the four test suites that read the record, the reproduction
script (a future rebuild is now compared with THIS library) and the runbook
need no edit; the four suites pass on the new record (54 checks), and the
installed catalog compares with the archive as in section 82 (565 same root,
11 faster, 0 slower, 0 missing, 82 label-only).

What changes for a reader of the catalog: the arrival phases are listed
SORTED (0.0337 first) where the hand-built file listed them from the anchor's
phase, so column k is not the column it was -- address cells by phase; and the
family names carry the campaign tag (`70mN_anchor` for `fast`).

`results/` is not tracked by git: the library of record exists on this disk
only. A deliverable zip, or tracking `library_70mN_24x24_final/` (3 MB without
the sheet), would fix that.

**The family label** (`family_index`). It says which continuation family -- one
connected branch of the minimum-time extremal set, traced by the arcs from one
anchor -- supplied an entry. The map (`family_map`) has two jobs: inside the
driver it decides whether a newly found root lies on a branch already walked
(no new arcs) or starts a new one (anchor it, walk it); in the catalog it is a
per-entry annotation. NOTHING reads the stored annotation: not the pickers,
not the audit, not the certifier. It matters for one future use, interpolating
a costate guess between neighbouring entries, which is only valid on one
smooth branch. For that use the label is the wrong instrument (it attaches by
flight time, and section 82 measured it non-reproducible, section 80 found 33
jumps INSIDE one label); the right one is already in hand: the trapezoid edge
residual of section 80, computed from the entries themselves.
