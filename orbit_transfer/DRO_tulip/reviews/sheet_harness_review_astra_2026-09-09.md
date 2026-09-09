## Bottom line

**Do not launch another production sweep unchanged.** The pasted harness has acceptance bugs, throws away the main advantage of multiple shooting, and does not implement the traversal described in its header.

My choice is **arrival-phase pseudo-arclength continuation in the full multiple-shooting variables, with multiple branch seeds and cheap departure ribs**. Replace `walkPhase`, not necessarily your dynamics or multiple-shooting engine.

Three important qualifications:

- **Your current code is already spine-and-ribs, not BFS.** Lines 108–156 do not walk arrival edges from every departure point. That criticism applies to an older implementation.
- **Pseudo-arclength passes regular folds, not the disappearance of an admissible branch.** You need the known slow-family seed as well as the fast-family seed.
- Your table establishes independently checked candidates. **It does not yet establish the global minimum-time envelope.** The observed 48% spread is real between those candidates; proving that it is the spread of the minimum-time value function requires more branch coverage.

Below I distinguish definite code defects from diagnoses that the supplied evidence cannot settle.

---

# A. Code review

## 1. The advertised certification gates are not enforced

This is the most serious correctness defect.

In `solvePoint`, lines 270–277:

- `tolDz` is **never used**.
- An exception from `tfMin` sets `dz = NaN`, but the point is still accepted.
- An arbitrarily large discrepancy between `za` and `zt` is accepted.
- A failed conjugate test, `cj = 0`, is accepted.
- A requested but missing conjugate result, `cj = -1`, is accepted.

The function finishes with:

```matlab
z = zt;  Y = it.Y;  ok = true;
```

Consequently, the actual acceptance rule is approximately:

> MS reports convergence, and the independently propagated final position is within 100 km.

That is not the acceptance rule stated in the header.

At minimum, enforce:

```matlab
witnessOK = isfinite(dz) && dz <= tolDz;
conjOK = ~doConj || ...
    (isfield(it,'conj') && isfield(it.conj,'pass') && it.conj.pass);
```

But **a small difference between two solver outputs is not itself proof that either output is a root**. An independent solver that fails and returns its initial guess gives `dz = 0`. Check its convergence status, if available, and independently evaluate the accepted solution’s terminal residual.

Also distinguish:

- `rootOK`: a usable continuation root;
- `physicalOK`: independently flown constraints pass;
- `localTestOK`: the requested conjugate test passes;
- `certified`: all required gates pass.

You may want to continue a mathematical branch after a local-optimality test fails. That does **not** justify recording it as an accepted minimum.

**This does not invalidate the independently certified points in your table. It invalidates the harness’s promise that every stored point received those certifications.**

### The flown gate is incomplete

Lines 267–269 test position only. Check, with explicit tolerances:

- final position **and velocity**;
- successful integration all the way to `tf`;
- `lam_m(tf)` and `H(tf)`;
- positive flight time and admissible mass;
- finite states and multipliers.

A 100 km position gate is a coarse sanity screen, not your “flown arrival approximately zero” certificate.

---

## 2. The anchor is silently accepted under potentially incompatible physics

Lines 93–104 always load `mintime_70mN_anchor.mat`, then record it with:

```matlab
fly = 0; dz = 0; cj = 1;
```

Yet callers can change:

- thrust;
- initial mass;
- Isp;
- orbit parameters;
- mesh specification.

There is no configuration match and no anchor re-solve.

For example, calling this with another thrust does not initiate thrust continuation. It simply labels the old solution as an accepted anchor under the new thrust.

**Fix:** validate the anchor’s exact endpoint states, full-precision phases, physical constants, normalization, and mesh contract against the requested run. Either reject a mismatch or explicitly re-solve/continue the anchor before starting.

Do not fabricate zero diagnostics. Store actual certification results or an explicit external-certificate reference.

---

## 3. You destroy the multiple-shooting warm start

This is the strongest actionable numerical problem in the harness.

You first predict only `z8` using a **single-shooting Jacobian**:

```matlab
[~, J, B] = sweep_phasing_shoot(...);
dz = J \ rhs;
```

Then `solvePoint`, lines 249–254, re-flies the entire trajectory from those initial costates:

```matlab
seed = seed_from_z8(zSeed, rv0(1:6), K_, ...);
```

That replaces the neighboring converged junctions with a fresh full-trajectory shot.

**You are repeatedly reintroducing the long-horizon sensitivity that multiple shooting was supposed to remove.**

After perturbing an ill-conditioned initial-costate estimate, a full re-flight can place distant junctions far from the neighboring BVP solution. It is not necessary for an MS initial guess to satisfy all continuity defects. Correcting those defects is its job.

### What to do instead

Predict the **entire packed MS unknown vector**:

- initial free costates;
- every junction state and costate;
- final time;
- homogeneous multiplier, when used;
- continuation parameter.

Use the uncondensed MS Jacobian for that prediction and correction. Keep the old junctions as the starting point if a predictor is unavailable.

The single-shooting predictor’s stated guarantee is also false:

> “a bad predictor must never be worse than no predictor.”

The check

```matlab
norm(dz) < 0.5*max(norm(z),1)
```

does not establish that. It uses an unscaled norm mixing time and costates, does not check the new residual, and does not enforce positive time.

I cannot verify the actual `J` and `B` implementation without `sweep_phasing_shoot`. The displayed sensitivity formula is consistent with the stated residual, assuming those derivatives are correct.

---

## 4. Failed refinement repeats already successful work

In `walkPhase`, every retry starts again from the original grid point:

```matlab
zc = z0; Yc = Y0;
```

If step 180 of 184 fails, the next attempt walks the whole edge again with 368 subdivisions.

That is not useful local adaptivity. It is whole-edge replay.

**Keep the last accepted continuation point. Shrink only the next attempted step.** Save those accepted intermediate points so a budget stop does not destroy hours of progress.

### The reported sub-solve count is wrong

On success:

```matlab
used = n;
```

But actual calls are counted by `nSolve`, including failed attempts and repeated prefixes.

Therefore, **if “368 sub-solves” came from this function’s log, it is not the total sub-solve count**. It is the subdivision count of the final successful traversal.

Return and store both:

- accepted continuation steps;
- total attempted corrector calls.

On failure, returning `used = NaN` discards the most important budget diagnostic.

### There is an exact-budget rejection bug

After each successful solve:

```matlab
if toc(tPt) > ptSec || nSolve >= maxSolve
    good = false;
    break
end
```

If the final required solve succeeds on call `maxSolve`, the point is rejected.

Check resource availability **before starting another solve**, and accept a completed target before deciding whether another call is permitted.

### Your defaults may exclude the known-good refinement again

For a 40,871 km chord:

- `stepKm = 2000` gives `n = 21`;
- with `maxBisect = 3`, the finest final pass is 168 subdivisions.

With `stepKm = 900`:

- `n = 46`;
- the finest final pass is 368 subdivisions.

Thus an 800-call budget does **not** ensure that the known-good subdivision level is reachable. The independent bisection cap can exclude it.

---

## 5. The clock backstop is not an end-to-end deadline

`wallSec` is passed to `ms_tfmin`. It does not visibly bound:

- the single-shooting predictor;
- `seed_from_z8`;
- the independent propagation;
- the independent `tfMin` solve;
- the extra propagation in `record`.

`ptSec` is checked only after those calls return.

**The code still cannot guarantee a wall-clock bound on one point.** Whether it can hang depends on the unshown routines.

Pass a shared deadline through every stage that can honor one. If a black-box routine cannot be interrupted, use process-level isolation/watchdog enforcement for a genuinely hard timeout.

There are also uncaught exceptions around the flown propagation and `record`. A successful MS solve followed by an integration exception can terminate the whole sweep rather than produce a classified failure.

Log stage, elapsed time, residual, current phase, step size, and failure reason during the walk—not just after an entire edge.

---

## 6. `K` is not actually controlling the mesh

`K` is passed into `solvePoint` and then ignored. Instead:

```matlab
K_ = size(YSeed,2);
```

This also exposes a contract ambiguity:

- Does `it.Y` contain `K` segment initial states?
- Or `K+1` nodes including the final endpoint?

The public seed contract elsewhere says `K+1`. The fallback here appends another column.

**I cannot call this a proven off-by-one without `ms_bvp`, but it must be asserted and tested.** Right now changing `opts.K` can change the metadata without changing the actual mesh.

Store the mesh explicitly, preserve it through continuation, and provide explicit pack/unpack/remesh functions.

---

## 7. Coverage is much more limited than “the sheet”

The arrival spine:

- goes in only one direction;
- stops at the first failed edge;
- does not use the independently known slow-family anchor;
- does not try an arrival continuation from another departure slice.

The departure ribs travel at most half a circle in each direction. If one direction fails early, the other direction is not allowed to take the longer route to the remaining points.

So missing entries mean:

> This particular traversal did not reach them within its rules.

They do not mean infeasibility, branch nonexistence, or even failure of continuation from every neighboring seed.

Store failure categories and branch identities. “Budget exhausted,” “Newton stalled,” “normality boundary,” and “local test failed” are different results.

---

## 8. Periodic interpolation needs validation

```matlab
interp1(..., mod(s,1)*period, 'spline')
```

does not by itself construct a periodic spline. Even if endpoint states agree, derivatives at the seam need not.

That matters especially when differentiating arrival phase.

Use a validated periodic orbit representation or accurate dense output with consistent derivatives. Check closure and derivative continuity, and check interpolation error against your terminal tolerances.

---

# Review of `arclength_thrust` and `ms_tfmin_hom`

The basic machinery is worth reusing:

- finite-differencing thrust at **fixed unknowns** is correct;
- the full augmented right-null-vector construction is correct;
- the homogeneous Hamiltonian and sphere normalization are algebraically consistent with positive homogeneity.

But several details need fixing before promoting this into the sheet engine.

## 9. Negative `rho` is incorrectly converted into a “normal” solution

In `ms_tfmin_hom`, line 79:

```matlab
if abs(info.rho) > 0
    info.zNormal = [p(1:7)/info.rho; ...];
end
```

The documentation correctly says conversion is valid for `rho > 0`. The code uses `abs(rho) > 0`.

**Dividing by negative `rho` is negative multiplier scaling. It does not preserve the minimizing thrust direction.**

Use a positive threshold and an explicit validity flag:

```matlab
normalValid = info.rho > rhoTol;
```

At very small positive `rho`, avoid automatically exporting enormous normal-chart multipliers as numerically usable data.

Also, `arclength_thrust` does not stop or classify the `rho = 0` event. Algebraic continuation through it can be useful diagnostically, but negative-`rho` roots must not enter the minimum-time candidate set.

Strictly, `rho = 0` is a possible abnormal PMP configuration—not proof that no physical transfer exists there.

---

## 10. Scaling is incomplete, and not mesh invariant

`Dx = max(abs(p0),1e-2)` provides column scaling. It does not provide:

- residual row scaling;
- a mesh-independent norm for distributed junction variables.

Duplicating trajectory nodes changes the Euclidean arclength even after componentwise amplitude scaling. Use a quadrature-weighted norm over normalized trajectory time, plus separate weights for scalar variables and the continuation parameter.

Your `fracCostate` diagnostic also excludes `tf`, thrust, and extras from its denominator. Therefore it is **not the fraction of total arclength consumed by costates**. It can read nearly one when most actual continuation motion is in time or thrust.

If that statistic is used to interpret branch geometry, fix its denominator to use the same complete metric as continuation.

---

## 11. A sign change is a fold candidate, not a certified fold

Lines 191–194 label every tangent-thrust sign change a fold.

Require:

- a localized parameter turning point;
- one-dimensional rank loss of the fixed-parameter Jacobian;
- full row rank of the augmented Jacobian;
- appropriate nondegeneracy/rank separation;
- protection against tiny-component sign noise and branch jumping.

A small tangent residual does not prove a unique one-dimensional nullspace. It can also occur when the augmented matrix has extra nullity.

---

## 12. The corrector needs globalization and branch-proximity control

The arclength corrector uses unrestricted Newton steps. It does not explicitly enforce positive time or mass, and it has no line search/trust region or accepted correction-distance limit.

That can produce failures—or convergence to a distant root that happens to satisfy the arclength plane.

Add:

- physical-domain checks;
- scaled merit-function globalization;
- a maximum correction relative to the predictor step;
- tangent-turning checks;
- classified exceptions and deadlines.

There is also a small step-controller edge case: a **successful** step can reduce the next `ds` below `dsMin`, after which the next correction loop never executes. Clamp successful proposals to the minimum allowed trial step; reserve “stalled below minimum” for a failed minimum-size attempt.

---

# B. What I would implement

## 1. A multi-branch arrival-phase PALC atlas

**I choose (a), with full MS variables and multiple seeds—not a single longer spine.**

Let \(p\) denote the packed multiple-shooting unknowns. At fixed thrust and departure phase, solve

\[
R(p,s_A)=0.
\]

There are \(n\) equations and \(n+1\) unknowns. Add the pseudo-arclength equation in properly scaled coordinates:

\[
\tau_k^\mathsf{T}(w-w_{\rm pred})=0.
\]

Use unwrapped arrival phase internally; apply modulo one only when evaluating the periodic endpoint.

### Arrival-phase differentiation is unusually cheap here

Because arrival is a **fixed state selected by phase**, changing \(s_A\) affects only the terminal state-matching conditions:

\[
R_{s_A} =
\begin{bmatrix}
0_{\text{continuity}}\\
-x_A'(s_A)\\
0_{\lambda_m}\\
0_H\\
0_{\text{normalization, if present}}
\end{bmatrix}.
\]

There is no need for the thrust driver’s two extra full residual evaluations to estimate this column.

For an exact periodic orbit,

\[
x_A'(s_A)=P_A f_{\rm CR3BP}(x_A(s_A)).
\]

If using an interpolant, differentiate that interpolant consistently.

### Concrete traversal

1. Seed the known fast and slow families. Use full-precision endpoint definitions, not rounded table phases.
2. Trace each arrival curve in both directions.
3. Locate **all** intersections with the requested arrival grid levels, including repeated intersections after folds.
4. Launch cheap departure ribs from those intersections.
5. Retain branch identities and all accepted candidates at a grid point.
6. Derive `S.TF` as the lowest time among discovered, accepted candidates.
7. When coverage is missing, start an additional slice/branch—not another blind replay of the failed spine.

I would use the homogeneous binding for branch tracing, after the fixes above, and convert to the normal chart for your existing certificates only while `rho` is safely positive.

Direct collocation would be my **targeted branch-discovery tool for holes and disconnected families**, not 144 unrelated primary solves. Neighbor-seeded collocation is still a local method and can follow the same wrong branch everywhere.

Following a two-parameter death locus is valuable later, especially if the existence boundary is itself part of the deliverable. It does not supply the interior time values, so I would not implement it first.

---

## 2. What does `cond(J) = 5e10` mean here?

**The supplied measurements do not identify the cause. Do not label it a fold, mesh failure, or ODE stiffness yet.**

They do strongly distinguish this segment from your previous normality-loss observation.

With your sphere normalization,

\[
\rho=\frac{1}{\sqrt{1+\|\lambda_0^{\rm normal}\|^2}}.
\]

The reported change \(17.6\to14.7\) corresponds approximately to

\[
\rho: 0.0567\to0.0679.
\]

On that measure, this segment is moving **away** from `rho = 0`.

### First identify which Jacobian is ill-conditioned

Log separately:

- the 8-by-8 single-shooting Jacobian;
- the uncondensed MS Jacobian;
- their scaled versions.

If the single-shooting Jacobian deteriorates while the scaled MS system remains healthy, that points toward long-flow sensitivity/condensation rather than a physical branch singularity.

Also inspect both singular-value extremes. A condition number can grow because the largest singular value grows, not because the smallest approaches zero.

### The decisive fold diagnostic

Form, with consistent scaling,

\[
A=W_R R_p D_p,\qquad b=W_R R_{s_A}s_*.
\]

Track:

- smallest singular values and rank separation of \(A\);
- smallest singular value of \([A\ b]\);
- arrival component of its oriented null tangent.

A regular arrival-phase fold should show:

- one rank lost by \(A\);
- augmented matrix \([A\ b]\) remaining full row rank;
- arrival tangent component reaching zero and reversing.

If the augmented matrix also loses rank, ordinary scalar PALC is no longer addressing a regular curve point.

### What else to measure

At a handful of accepted points:

- repeat on \(K=24,48\), using the same branch and consistent scaling;
- monitor segment STM growth;
- check Jacobian directional derivatives with a step-size study;
- monitor minimum \(\|\lambda_v(t)\|\), not just initial multiplier norm;
- monitor close approaches and integration tolerance sensitivity.

The direction law contains derivatives proportional to \(1/\|\lambda_v\|\). A small primer norm anywhere along the trajectory can cause difficulty without initial multiplier blow-up.

Also distinguish **unstable long-horizon sensitivity** from classical ODE stiffness. They are not synonymous.

A passing conjugate test does not bound numerical conditioning. It also does not exclude approaching a singularity that has not yet been reached; whether it detects the relevant degeneracy depends on its treatment of free final time.

**I would test the full-MS reformulation first. The current code contains enough numerical self-sabotage that failure to close an edge is not yet clean evidence about branch geometry.**

---

## 3. Can freeing arrival phase make it cheap?

**Freeing it as a continuation coordinate: yes, potentially much cheaper. Freeing it as an optimization variable: that changes your problem.**

If you truly optimize arrival phase, the terminal-manifold transversality condition includes

\[
\begin{bmatrix}\lambda_r(t_f)\\ \lambda_v(t_f)\end{bmatrix}^{\!T}
x_A'(s_A)=0.
\]

That selects stationary arrival phases. It does not produce the fixed-phase deployment envelope.

If you introduce \(s_A\) and then add \(s_A-s_{A,j}=0\), you have merely repackaged the original fixed-point problem. Its intrinsic singularity is still there.

PALC helps by allowing phase to move during tracing and avoiding a singular projection onto phase. You subsequently extract the prescribed-phase intersections.

### The largest immediate savings are elsewhere

Every internal step currently may involve:

1. a full single-shooting variational propagation for prediction;
2. a full re-flight to construct junctions;
3. another full independent flight;
4. an iterative independent single-shooting re-solve.

That is a lot of single shooting inside “multiple-shooting continuation.”

Use MS residual/physical checks for internal continuation points. Run the expensive independent certificates at grid intersections and selected diagnostic checkpoints. Do not certify every tiny sub-step and then ignore the certificate.

### A smaller exact simplification

Under your all-burn assumptions,

\[
m(t)=1-\frac{T_{\rm nd}}{c_{\rm nd}}t.
\]

Moreover, \(\lambda_m\) does not feed the direction-controlled \(r,v,\lambda_r,\lambda_v\) subsystem. It can be recovered by quadrature using \(\lambda_m(t_f)=0\).

Thus one can reduce to a 12-state propagation with seven shooting unknowns \([\lambda_{r0},\lambda_{v0},t_f]\), imposing the six endpoint conditions and terminal Hamiltonian condition.

That removes nuisance variables; it does **not** remove folds. I would make the existing MS phase continuation work before doing this additional refactor.

---

## 4. Replace the endpoint-chord subdivision rule

Your current rule does not actually target 900 km of endpoint motion per sub-step.

It uses the **chord between grid endpoints**, divided into uniform phase increments. On a large, curved, seven-petal orbit, chord length is not arc length, and phase speed is nonuniform. It also ignores terminal velocity motion.

Using your numbers:

- grid-endpoint chord: 40,871 km;
- `stepKm = 900` gives 46 subdivisions;
- phase increment: approximately \(0.00181\).

Where your measured \(0.0005\) phase increment produces about 1,100 km of motion, that subdivision corresponds locally to roughly **4,000 km**, not 900 km.

This alone can explain why the nominal geometric step setting repeatedly needs refinement. It does not explain every stall.

### The controller I would use

For each accepted full-MS point:

1. Compute a scaled tangent.
2. Predict by arclength \(h\).
3. Apply a globalized corrector.
4. Measure:
   - Newton effort;
   - predictor-to-corrector distance;
   - tangent turning angle.
5. Grow \(h\) conservatively after easy, low-curvature steps.
6. Shrink \(h\) locally after difficult or failed steps.
7. Never replay the accepted prefix.

Use a weighted six-state endpoint-motion estimate as an additional safety cap, not as the principal controller.

Keep separate:

- numerical step rejection;
- resource exhaustion;
- certification failure;
- physical/admissibility events.

A true conjugate-test failure or a localized `rho = 0` event is not something to cure by blindly doubling subdivisions.

---

# What I would do before another long run

1. **Fix acceptance and anchor validation.** Add tests where the witness throws, returns a large discrepancy, or the conjugate result fails/is missing.
2. **Test budgets with a mock corrector:** success exactly on the last permitted call, failure near the end of an edge, and resume from the last accepted intermediate point.
3. **Run one short arrival segment using the full MS tangent and stored junctions.** No single-shooting predictor and no full re-flight seed.
4. On that same segment, compare scaled MS conditioning at \(K=24\) and \(48\), and validate the analytic phase derivative.
5. Then trace both known families with phase PALC and extract grid intersections.

The nondimensional thrust conversion and the displayed \(\Delta V\) formula look consistent with normalized mass. Use the analytic mass law as a regression check. I cannot validate the underlying STM, conjugate implementation, or node-storage contract without the unshown engine routines.

**My recommendation is not “spend more compute on finer BFS.” It is: repair the certificates, stop collapsing MS back into single shooting, and build a branch-aware phase atlas.** That addresses the problem you actually have without pretending that continuation can manufacture a fast branch where none exists.