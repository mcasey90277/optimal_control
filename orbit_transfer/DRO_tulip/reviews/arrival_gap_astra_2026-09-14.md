## Bottom line

Your diagnosis is substantially right about **numerical conditioning and the incompleteness of the continuation search**, but I would challenge three central interpretations:

1. **The minimum-time problem is exactly periodic in \(s_A\).** An individual continued sheet of extremals need not be periodic, but the optimal value and the set of solutions are. In particular, the certified 16.4-day solution at \(s_A=0.0754\) is also a certified solution at \(s_A=1.0754\). A 28-day A2 solution near that phase cannot be interpreted as the only available “next turn” without searching the A1 branch backward across the phase seam.

2. **If a continuous A2 arc goes from \(0.9087\) to \(0.8254\), it necessarily passes through \(0.8671\).** Therefore “no extremal found at \(0.8671\)” is, on the stated ancestry, primarily a continuation/extraction question—not evidence of a dynamical hole. Exceptions would be an intervening inadmissible segment, a discontinuous branch jump, or a different phase lift.

3. **A phase fold is not, by itself, a reachable-set boundary or a winding-number transition.** For a normal branch it is usually a singularity of the endpoint/shooting projection, associated generically with a terminal Jacobi degeneracy and a change of local optimality.

I agree that the three A2 transversality refusals are **strong candidates for forward-integration evaluation artifacts**. But “the multiple-shooting residual is \(10^{-12}\)” is not sufficient to establish that: one needs a scaled error budget and a convergence experiment.

Finally, throughout this answer, “certified minimum” means the local optimality certified by your second-order theorem. None of the listed tests, including a local direct-solver witness, establishes global minimum time among disconnected extremals.

---

# 1. What actually determines the phase–time slope?

There is a useful exact identity here. It gives a much sharper interpretation than the helix picture.

Let
\[
y=(r,v),\qquad
f_0(y)=
\begin{pmatrix}
v\\ a_{\rm CR3BP}(r,v)
\end{pmatrix},
\]
and use the minimum-principle convention
\[
H=
1+p_r\cdot v+p_v\cdot a_{\rm CR3BP}
+\frac{T}{m}p_v\cdot u-\frac{T}{c}p_m.
\]
For a normal extremal,
\[
u=-\frac{p_v}{\|p_v\|},
\]
so
\[
H=1+p_y\cdot f_0(y)-\frac{T}{m}\|p_v\|
-\frac{T}{c}p_m=0.
\]

Here the objective multiplier is normalized to one. Other sign conventions give the same physical identities after translating the costates.

Because the arrival phase is the time phase of the unforced tulip,
\[
y_A'(s_A)=P_A f_0(y_A(s_A)),
\qquad P_A=\frac{5\pi}{3}.
\]

For any smooth branch of stationary solutions, with fixed departure state and free final mass, the endpoint sensitivity formula is
\[
\boxed{
\frac{dt_f}{ds_A}
=-p_y(t_f)\cdot y_A'(s_A).
}
\]
There is no \(p_m\,dm_f/ds_A\) contribution because \(p_m(t_f)=0\).

Using \(H(t_f)=0\),
\[
p_y(t_f)\cdot f_0(y_A)
=-1+\frac{T}{m_f}\|p_v(t_f)\|,
\]
and therefore
\[
\boxed{
\frac{1}{P_A}\frac{dt_f}{ds_A}
=
1-\frac{T}{m_f}\|p_v(t_f)\|.
}
\tag{1}
\]

This is valid on a smooth, unconstrained normal branch—not merely on the globally fastest branch.

## What your slope means

Your observed slope of approximately 15 days per phase unit gives
\[
\frac{15}{23.2}\simeq 0.65.
\]
Thus the terminal costate should satisfy approximately
\[
\boxed{
\frac{T}{m_f}\|p_v(t_f)\|\simeq 0.35.
}
\]

That is a directly testable prediction using your existing trajectories.

The slope is determined by the **terminal value sensitivity and the terminal thrust contribution to the Hamiltonian**. There is no universal \(0.65\), and no period-locking law selecting it.

Equation (1) also says:

- The slope need not be positive.
- On a strict-bang normal branch with nonzero terminal primer,
  \[
  \frac{dt_f}{ds_A}<P_A.
  \]
- At a stationary phase–time point,
  \[
  \frac{T}{m_f}\|p_v(t_f)\|=1.
  \]

Your minimum near \(s_A=0.159\) is consistent with that last condition. Strictly speaking, your A1 time history is not monotone over its whole reported interval: it first decreases from 16.4 to 16.23 days.

**Recommended diagnostic:** compare finite-difference or tangent-derived \(dt_f/ds_A\) against (1). It simultaneously checks the phase parametrization, multiplier normalization, continuation tangent, and terminal PMP data.

## Why autonomy does not make the time independent of phase

Autonomy removes absolute epoch from the problem. It does not make different points of a periodic orbit equivalent as fixed endpoint states.

At fixed departure state, changing \(s_A\) changes six terminal conditions. The associated cost is the projection
\[
-p_y(t_f)\cdot y_A'(s_A).
\]
That projection can have either sign.

Likewise, the reported arrival radii, speeds, and \(z\)-coordinates are too coarse to predict difficulty. What matters is the endpoint’s relationship to the thrust-driven transport geometry in **position–velocity space**, not simply whether it looks apolune-like.

---

# 2. Periodicity, the “helix,” and the most important missing search

Define the global value
\[
V(s_D,s_A)=\inf\{t_f:\text{admissible transfer to }y_A(s_A)\}.
\]
Since
\[
y_A(s_A+1)=y_A(s_A),
\]
you have exactly
\[
\boxed{
V(s_D,s_A+1)=V(s_D,s_A).
}
\]

The entire solution set is also periodic under relabeling \(s_A\mapsto s_A+1\).

An **individual lifted branch** may fail to close after one phase turn. Continuation around the phase circle can permute branches, or produce a sheet with additional dynamical excursions. A helical picture is therefore possible for a stationary solution sheet—but not for the single-valued global optimum.

Also, your quoted comparison is not at identical endpoints:
\[
1.034-1=0.034\ne 0.0754.
\]
The actual reported time difference is
\[
28.12-16.4=11.72\ {\rm d},
\]
over approximately \(0.959\) phase units. A full-turn increment near 12–13 days is an extrapolation, not an observed same-endpoint comparison.

## The omitted continuation direction

I would immediately continue A1 in **decreasing** arrival phase:

\[
0.0754\longrightarrow
0.0337,\ -0.0079,\ -0.0496,\ldots
\]

Modulo one, these are precisely the tail phases
\[
1.0337,\ 0.9921,\ 0.9504,\ldots
\]

This is potentially more consequential than making the 28-day solutions pass their numerical gate. It asks whether the tail columns admit a roughly 16–something-day continuation of A1 that your one-directional arrival-axis construction never explored.

There is no guarantee that A1 reaches all those phases: it might encounter a different fold, clearance boundary, or loss of optimality. But it is the natural search to perform before accepting the A2 times as the best available ones.

At the exact phase \(1.0754\), there is no uncertainty: A1 already supplies a 16.4-day normal certified solution. Moreover, if A1 is nondegenerate and has a positive clearance margin, the implicit-function theorem supplies normal nearby solutions on **both sides** of that phase.

Thus an A2 termination before \(1.0754\) cannot establish a genuine disappearance of normal sub-30-day transfers in a neighborhood of \(1.0754\).

---

# 3. What the fold at \(0.8202\) means

Write the fixed-phase shooting equations as
\[
F(z,s_A)=0,
\]
where \(z\) includes the shooting variables and \(t_f\).

At a generic simple fold,
\[
F_z \text{ has a one-dimensional kernel},
\]
while
\[
[\,F_z\;\;F_{s_A}\,]
\]
has full row rank. Along arclength \(\alpha\),
\[
\frac{ds_A}{d\alpha}=0
\]
at the fold, with a nonzero second derivative.

This means that the projection of the stationary-solution manifold onto the chosen phase coordinate ceases to be locally one-to-one. On one side there are typically two nearby stationary solutions; on the other, none from that local sheet.

## A normal fold has a distinctive time geometry

The endpoint sensitivity identity along arclength is
\[
\frac{dt_f}{d\alpha}
=
-p_y(t_f)\cdot y_A'(s_A)\frac{ds_A}{d\alpha}.
\]
If the normal costate remains finite, then at the phase fold,
\[
\boxed{
\frac{ds_A}{d\alpha}=0
\quad\Longrightarrow\quad
\frac{dt_f}{d\alpha}=0.
}
\]

So the full solution branch folds in shooting space, but its projection into \((s_A,t_f)\) is not generically an ordinary vertical turn with \(dt_f/d\alpha\ne0\).

Near a generic variational fold, the two time branches often share the same limiting slope, with their difference behaving like
\[
\Delta t_f=O(|s_A-s_A^\ast|^{3/2}).
\]
Nearly equal times are therefore compatible with a fold pair—but do not prove one.

Your two times at \(0.8254\) differ by about 5.8 minutes. To identify them as a fold pair, demonstrate that they coalesce in:

- full shooting variables;
- trajectory geometry;
- continuation arclength;
- a single kernel direction;
- second-variation index.

Notice also that \(0.8254>0.8202\). Those solutions are not simply the two local sides of the reported A1 maximum. They require another part of the global branch or another branch.

## Connection to conjugate points

Under the regularity assumptions behind your quotiented free-time test, singularity of the appropriate fixed-endpoint shooting map corresponds to a nontrivial Jacobi field satisfying the linearized boundary conditions.

Thus a generic fold is closely associated with a **terminal conjugate degeneracy**. A generic minimum–saddle pair changes Morse index by one across the fold.

That gives a concrete prediction:

> Track the smallest singular value of the correctly quotiented endpoint Jacobi map and the conjugate index toward the A1 fold. If it is an ordinary variational fold, the terminal degeneracy should appear there.

This prediction has qualifications: symmetries, chart degeneracies, unremoved gauge directions, or a higher-order singularity can contaminate the shooting rank test.

Your reported \(\rho=0.05\to0.10\) is good evidence against loss of normality there, assuming a consistent chart normalization. It is not, by itself, a classification of the singularity.

## “Augmented Jacobian” needs clarification

If the \(10^{-9}\) singular value belongs to \(F_z\), that is exactly where a fold signature should appear.

If it belongs to the fully bordered pseudo-arclength Newton matrix,
\[
\begin{pmatrix}
F_z&F_{s_A}\\
t_z^\top&t_s
\end{pmatrix},
\]
a simple fold should not make that matrix singular. Then I would investigate scaling, tangent accuracy, or an additional degeneracy. A numerical dip to \(10^{-9}\) is not itself proof of rank loss, but the distinction matters.

## Is it tangency to a minimum-time front?

Not ordinarily.

Where the value is smooth, tangency of the tulip to a level set of minimum time is
\[
\nabla_y V\cdot y_A'(s_A)=0,
\]
equivalently
\[
\frac{dV}{ds_A}=0.
\]
Your minimum near \(0.159\) is a candidate for that geometry.

The fold near \(0.8202\) instead concerns the singularity of the extremal endpoint projection—a caustic-type phenomenon. It need not lie on the global minimum-time front. A branch can cease to be globally minimizing at a Maxwell/cut point **before** it reaches a conjugate point or fold.

Nor does a fold imply a change of lunar winding.

---

# 4. Is A2 the next turn, and what does falling \(\rho\) mean?

## Same component versus different dynamical itinerary

The available data do not identify A2 as the next turn of A1.

Evidence for that claim would be an actual continuous connection between them in the full solution space, possibly through several folds. Similar time slopes are weak evidence: equation (1) allows unrelated branches to have similar slopes.

To characterize the dynamical difference, compare solutions at the **same phase**, using:

- Moon-centered angular accumulation in a specified projection;
- number and sequence of close lunar passages;
- crossings of selected Poincaré sections;
- residence near particular periodic orbits or invariant manifolds;
- thrust-work histories and Jacobi-constant histories;
- the continuation connection itself.

For example, with the conventional CR3BP Jacobi constant,
\[
\dot C=-2\frac{T}{m}v\cdot u.
\]
The distribution of this work over encounters can distinguish genuinely different transport mechanisms.

One caution: in this spatial problem, “one more revolution” is generally not a rigorous homotopy class. The exterior of a ball in three dimensions is simply connected. A projected winding count can change by moving out of the projection plane without hitting the Moon. It remains a useful itinerary label, but it is not usually a protected topological invariant.

## What small \(\rho\) establishes

If \(\rho\) is the objective multiplier in a fixed projective normalization, then \(\rho\to0\) means the multiplier approaches the abnormal set. In the normal chart \(p_0=1\), the costate norm may diverge.

That is a statement about **constraint qualification, endpoint-map degeneracy, and the relative weight of the time objective in the multiplier system**.

It is not directly a statement that:

- fuel is running out;
- the target is unusually far away;
- the target is approaching a special physical reachability boundary;
- all other normal branches disappear.

Every genuine minimum-time endpoint already lies on an appropriate minimum-time reachability frontier. Near-abnormality is a degeneracy of its supporting multiplier geometry, not merely “being near the boundary.”

You should also distinguish:

- an extremal admitting an abnormal lift;
- a strictly abnormal extremal admitting no normal lift;
- a limiting normal lift becoming abnormal.

These are not interchangeable. Your abnormal-lift dimension test is important precisely here.

## Your three \(\rho\) values do not predict an imminent zero

The trend
\[
0.0138,\quad0.0131,\quad0.0126
\]
over \(s_A=0.950,0.992,1.034\) is a decrease of about 9%, not an observed collapse.

A crude local linear extrapolation gives roughly
\[
\rho(1.075)\approx0.012,
\]
not zero. Such extrapolation is not a continuation guarantee, but the current numbers provide little affirmative evidence for an endpoint before \(1.075\).

The earlier 72 mN result is a reason to monitor the mechanism, not evidence that it must recur at these phases.

## A sharper abnormality diagnostic from the slope

If the homogeneous Hamiltonian uses objective multiplier \(\rho>0\), equation (1) becomes
\[
\frac{1}{P_A}\frac{dt_f}{ds_A}
=
1-\frac{T\|\lambda_v(t_f)\|}{m_f\rho}.
\tag{2}
\]

Consequently, if \(\rho\to0\) while the phase slope remains finite and near \(0.65P_A\), then
\[
\|\lambda_v(t_f)\|=O(\rho).
\]

If instead the homogeneous terminal velocity costate remains bounded away from zero, the slope must become large and negative as \(\rho\to0^+\).

This is a useful discriminator. Monitor together:

\[
\rho,\qquad
\|\lambda_v(t_f)\|,\qquad
\frac{\|\lambda_v(t_f)\|}{\rho},\qquad
\frac{dt_f}{ds_A},
\]
plus the endpoint-map singular values and abnormal-lift dimension.

A finite positive slope persisting into a genuine abnormal limit requires a special terminal multiplier degeneration. It should not be invisible in those diagnostics.

---

# 5. The \(0.8671\) “no extremal” result needs an immediate audit

This is the strongest logical issue in the reported continuation history.

If an A2 arc is continuous in unwrapped phase and joins
\[
s_A=0.9087
\quad\text{to}\quad
s_A=0.8254,
\]
the intermediate value theorem guarantees at least one point on that arc with
\[
s_A=0.8671.
\]

Folds do not remove this fact.

Therefore inspect the saved arc points and find every sign change of
\[
s_A(\alpha)-0.8671.
\]
For each bracket:

1. refine in arclength;
2. reproject onto \(F=0\);
3. locate the phase event;
4. perform fixed-phase correction if regular;
5. classify the resulting extremal.

Possible outcomes are informative:

- **A normal extremal appears:** grid extraction or correction failed.
- **The crossing is Moon-unsafe or conjugate:** an extremal exists, but not an accepted one.
- **The crossing has \(\rho<0\) or violates strict-bang assumptions:** the algebraic continuation passed through a PMP-inadmissible region.
- **There is no continuous crossing in the log:** the claimed branch ancestry contains a jump, restart, phase-wrap mismatch, or another change of problem parameters.

Until this is resolved, I would not call \(0.8671\) a physics gap.

Also, extending A1 immediately past its maximum at \(0.8202\) walks toward **smaller** phase. It cannot directly reach \(0.8671\). It may eventually find another turn, but that is a global continuation possibility, not the local implication of passing the fold.

---

# 6. The transversality refusal: likely numerical, but verify it properly

The pattern is persuasive:

- multiple-shooting residuals stay near \(10^{-12}\);
- the physical single-flight miss grows with duration;
- the same gate stops the \(0.909\) departure rib;
- the terminal mass-costate errors are only modestly above the threshold.

I would treat these as **numerical-validation failures pending diagnosis**, not evidence of a physical endpoint to the family.

But I would not yet label them harmless.

## Why a tiny multiple-shooting residual is insufficient

Let \(w=(x,p)\) denote the state–costate variables, and let
\[
d_k=w_{k+1}-\Phi_k(w_k)
\]
be the multiple-shooting junction defects.

To first order, the discrepancy between the globally flown trajectory and the shooting nodes satisfies
\[
e_{k+1}=D\Phi_k\,e_k-d_k+\eta_k,
\]
where \(\eta_k\) collects flow-evaluation errors.

Thus
\[
e_N\approx
-\sum_k\Phi_{N,k+1}d_k
+\text{propagated integration and initial-data errors}.
\tag{3}
\]

A \(10^{-12}\) defect can indeed produce a \(10^{-6}\) terminal error. But a small residual is only a backward-error statement; closeness to an exact extremal also depends on conditioning and integration accuracy.

Furthermore, a state arrival miss does not directly calibrate a mass-costate miss. They are different components of a highly anisotropic propagation map.

## The correct primary check

For a multiple-shooting extremal, evaluate the BVP and PMP on the **piecewise integrated multiple-shooting representation**, with explicit control of junction defects.

For terminal transversality, check the endpoint obtained by integrating the final segment from its initial node:
\[
p_m^{\rm last\ segment}(t_f).
\]

If you have an independently stored terminal costate constrained to zero, merely reading that stored component is insufficient. You must also verify the final segment’s matching defect.

The full validation should include:

- segment integration error estimates;
- scaled state and costate junction defects;
- terminal boundary conditions from the propagated last segment;
- Hamiltonian and control checks on each segment;
- consistency of costate normalization across all nodes;
- a global reconstruction/error assessment.

## A useful independent scalar identity

With the convention above,
\[
\dot p_m=-\frac{T}{m^2}\|p_v\|.
\]
Hence
\[
\boxed{
p_m(0)=\int_0^{t_f}\frac{T}{m(t)^2}\|p_v(t)\|\,dt
}
\tag{4}
\]
when \(p_m(t_f)=0\).

Evaluate this integral segment by segment, preferably independently of the mass-costate ODE integration. It will not magically eliminate conditioning, but it is an excellent cross-check for sign, scaling, quadrature, and cancellation problems.

## The experiment that would convince me

Take the certified \(0.909\) solution and the three rejected neighbors, and repeat:

1. single shooting with successively tighter tolerances;
2. higher-precision initial data and arithmetic;
3. an independent high-order integrator;
4. forward/backward integrations meeting at interior nodes;
5. multiple-shooting refinement with changed segmentation.

Then determine whether:

- the segmentwise endpoint transversality converges to zero;
- the global flown discrepancy decreases as predicted;
- equation (3), using measured variational propagation, explains its size and sign;
- the trajectories and times converge to a common limit.

If so, the gate-artifact diagnosis is compelling.

If the discrepancy plateaus despite higher precision, inspect chart conversion, stale or rounded initial costates, inconsistent dynamics, time scaling, and whether the displayed residual is scaled or omits the endpoint equation being checked.

## How should tolerance be set?

Not as an unexplained universal \(10^{-6}\), and not by blindly multiplying it by a large flow condition number.

Use three separate quantities:

1. **A declared stationarity tolerance** in a fixed physical/nondimensional normalization.
2. **An estimated numerical uncertainty** in the evaluated residual.
3. **A conditioning report** connecting the BVP representation to a flown reconstruction.

Costates scale with the multiplier convention. In a homogeneous chart, raw \(|\lambda_m(t_f)|\) is not comparable across different \(\rho\) unless normalization is fixed. Near abnormality, a residual small relative to \(\|\lambda\|\) can still be large in the normal chart.

I would report both a consistently scaled homogeneous residual and its normal-chart equivalent. For Hamiltonian bookkeeping, the mass-costate term has scale
\[
\frac{T}{c}\frac{\lambda_m}{\rho}.
\]
That is informative, although it is not a substitute for terminal stationarity itself.

**Amplification belongs in the numerical error budget, not in an unlimited acceptance allowance.** If you retain a single-flight \(10^{-6}\) requirement, increase precision until its numerical uncertainty is comfortably below \(10^{-6}\). Alternatively, certify the BVP representation with an explicit reconstruction/existence bound.

---

# 7. The \(0.8254\) candidates: distinguish three possibilities

The clear candidate with an unresolved conjugate test deserves immediate numerical refinement.

Use more than a raw determinant:

- a well-conditioned symplectic/Jacobi frame;
- QR or SVD stabilization;
- chart changes on the relevant Lagrangian Grassmannian;
- a Riccati or Maslov-index formulation where appropriate;
- crossing forms to classify candidate crossings;
- increased precision and interval/bracket refinement.

Check that the free-time/free-mass quotient removes the intended trivial directions and no more.

The outcomes have different meanings:

### A. No conjugate point, with a resolved positive margin

Then this is a genuine certification-policy hole. Certify it and seed ribs.

### B. A first conjugate point strictly before \(t_f\)

Then, under the theorem’s hypotheses, this candidate is not the required local minimum. The near-zero was real.

### C. Degeneracy exactly at \(t_f\)

Then strict second-order sufficiency fails, but nonoptimality does not automatically follow. A degenerate minimum may require higher-order analysis. Do not equate “not strictly certified” with “refuted.”

The Moon-unsafe candidate is a separate matter. If it passes the unconstrained second-order test, it may be the unconstrained local minimum whose safe alternative is not another smooth unconstrained extremal.

That leads to an important missing route.

## Clearance is a path constraint, not merely a filter

If \(1900\) km is part of the optimization problem, an optimal safe transfer may have active lunar-clearance contacts or boundary arcs.

For
\[
g(r)=R_{\min}^2-\|r-r_M\|^2\le0,
\]
the active constraint has higher relative degree with respect to thrust. The constrained PMP requires the appropriate path-constraint multipliers and contact/boundary conditions.

Such an optimum is generally **not covered by the smooth unconstrained PMP and conjugate test** you use for strictly clear trajectories.

Therefore a library formed by solving the unconstrained problem and rejecting unsafe solutions can have genuine coverage holes even when safe constrained minima exist.

At \(0.8254\), I would explicitly search for that possibility.

---

# 8. Routes to close the gap, ranked

I would use the following order. Some tasks can run in parallel.

## Priority 1: Exploit information already present

### 1A. Extract the \(0.8671\) crossings from the A2 arc logs

This has the highest information-to-computation ratio. On the stated continuous ancestry, the crossing must exist somewhere.

**Tells you:** whether the hole is extraction failure, inadmissibility, or incorrect branch bookkeeping.

### 1B. Continue A1 backward across the phase seam

Try the phases equivalent to \(1.034,0.992,0.950\), then farther backward.

**Tells you:** whether A2 is merely a slower local branch in phases where a much faster branch was never searched.

### 1C. Repair the transversality evaluation and rerun A2 ribs

Do the precision/error-budget experiment before globally changing tolerance.

**Tells you:** whether three arrival-column failures and much of the \(0.909\) rib failure are numerical.

### 1D. Resolve the \(0.8254\) conjugate candidate

This may close one column without finding any new trajectory.

**Tells you:** policy floor versus genuine loss of local optimality.

## Priority 2: Targeted continuation, not just more steps

### 2A. Continue A2 downward with exact phase-event extraction

Use small steps and conditioning-aware precision near turns. Preserve rejected extremals and continue them when mathematically legitimate; rejection from the certified library need not imply that a branch should be discarded as a search object.

### 2B. Continue A1 through its fold

Do this to characterize:

- the fold normal form;
- the Morse-index change;
- possible later turns;
- connection to A2.

I rank this below the seam search and the guaranteed A2 phase crossings for immediate gap closure. A1’s returning local branch initially moves away from \(0.867\).

### 2C. Build a two-parameter phase atlas

Do not make every departure rib depend on one arrival-axis seed.

Continue in \((s_D,s_A)\), use neighbors in both directions, and trace fold curves. A fold blocking one fixed-\(s_D\) slice may be bypassed in the two-parameter family.

**Tells you:** whether the library hole is an artifact of the chosen construction tree.

## Priority 3: Direct multistart solves at the difficult phases

Yes: collocation/IPOPT, followed by indirect polishing and your certification pipeline.

I would use:

- multiple initial times, roughly spanning the fast and A2 regimes;
- A1 seam seeds;
- A2 seeds;
- perturbed encounter sequences;
- the unsafe \(0.8254\) solution with clearance enforced;
- explicit \(t_f\le30\) day searches;
- mesh refinement and between-node clearance checks.

A direct solver is valuable because it need not remain on the branch your indirect continuation selected.

But failure of a collection of direct solves is not an infeasibility proof. Also, an active-clearance solution must be polished with the constrained PMP, not forced into the unconstrained shooting formulation.

## Priority 4: Thrust continuation, preferably in two parameters

Continuation from higher thrust is sensible, but I would map the \((T,s_A)\) geometry rather than use only isolated fixed-phase thrust homotopies.

Track:

- fold loci;
- normality-loss loci;
- clearance activation;
- terminal conjugate degeneracy;
- equal-time crossings between locally minimizing branches.

This distinguishes a fold barrier that can be bypassed from a genuine boundary of a particular normal sheet.

Higher thrust may improve the optimal value in the bounded-magnitude problem, but it does not guarantee that the specific branch you want survives continuously down to 70 mN.

## Priority 5: Homogeneous continuation to and through \(\rho=0\)

Useful for understanding branch organization and obtaining limiting abnormal trajectories.

Not a general gap-closing remedy.

Crossing to a negative objective multiplier in the same minimum-principle convention does not produce a valid minimum-time PMP extremal. A formal continuation through \(\rho=0\) can be mathematically informative while leaving the admissible multiplier cone.

A genuinely abnormal endpoint also needs an appropriate optimality analysis; your normal second-order certification cannot simply be carried through unchanged.

## Priority 6: Accept and certify a longer local minimum

This is appropriate for producing a feasible, locally optimal library entry after its numerical issues are resolved.

It is not yet appropriate to call it the minimum-time answer, or the necessary “next turn,” unless faster branches—especially the A1 seam continuation—have been addressed.

---

# 9. Is coverage of every phase a well-posed expectation?

There are three different expectations.

## A. A feasible transfer to every phase

This is a reachability question. Autonomy and periodicity do not guarantee it.

A phase can be unreachable within a prescribed horizon because of:

- insufficient time/thrust authority for the six-dimensional endpoint match;
- avoidance constraints;
- mass or other resource limits;
- dynamical transport barriers over that horizon;
- the restriction to full thrust if taken as an exact control-set restriction.

For your engine, fuel exhaustion is not close at 30 days. With \(m_0=1\),
\[
m(t)=1-\frac{T}{c}t,
\qquad \frac{T}{c}\simeq0.02024.
\]
At 28.12 days, \(t\simeq6.34\) ND, giving
\[
m_f\simeq0.872,
\]
or about 131 kg. The difficulty is transport and endpoint matching, not imminent depletion.

Because exact full-burn controls do not allow waiting, exact-time reachable sets need not be nested. For minimum-time questions, use the reachable-by-time set
\[
\mathcal R_{\le t}=\bigcup_{\tau\le t}\mathcal R_\tau.
\]

## B. An attained minimum-time transfer to every reachable phase

Existence requires hypotheses.

With a compact convex bounded-thrust control set, positive mass margin, trajectories confined away from singularities, and closed path constraints, standard existence arguments are available when the target is reachable within a finite bound.

With **exactly fixed thrust magnitude**, the control set is a sphere, not a ball. Relaxation/chattering issues must be considered; existence is not automatic from the usual convex-control theorem. Often one proves existence in the bounded-magnitude problem and then proves the optimum is full thrust almost everywhere.

Likewise, strict clearance \(d>1900\) defines an open constraint and can yield an unattained infimum. Closed clearance \(d\ge1900\) avoids that particular issue but may produce active-constraint minima.

## C. A strictly normal, unconstrained, second-order-certified minimum to every phase

This is a much stronger expectation and is not generally justified.

Even when feasible minima exist, some phases may have:

- multiple equal-time minima;
- a conjugate degeneracy at the endpoint;
- active lunar-clearance contacts;
- abnormal minimizers;
- loss of the regularity required by your particular certificate.

Abnormal-only optimality requires failure of an appropriate endpoint constraint qualification. A useful sufficient route to normality is regularity of the endpoint constraint map with respect to admissible control and final-time variations. A genuine strictly abnormal optimum lacks a normal multiplier representation.

But one branch approaching \(\rho=0\) is nowhere near a proof that a phase is reachable only abnormally—or that no normal sub-30-day solution exists. That would require a substantially more global reachability or branch-completeness argument.

---

## My revised diagnosis

I would replace your diagnosis with:

> The library has a branch-selection and numerical-validation gap, superimposed on genuine singularities of particular extremal families. A1 exhibits a normal projection fold, but neither that fold nor A2’s small objective multiplier establishes a gap in the global minimum-time solution. Three refusals are plausibly caused by unstable global evaluation of a well-resolved BVP; one candidate needs a sharper conjugate analysis; another violates a constraint whose active solutions may be absent from the formulation. The supposedly empty intermediate phase must first be reconciled with the continuous A2 arc ancestry. The unsearched backward continuation of A1 across the phase seam is a major alternative source of fast solutions.

One bookkeeping point also needs reconciliation before measuring progress: \(19\times24=456>406\), and the listed contiguous phases appear to comprise 18 columns plus the \(0.909\) column. I suspect “fully certified columns” means certified arrival-axis seeds rather than completed departure ribs. Keep separate masks for **extremal found**, **numerically resolved**, **locally optimal**, **clearance-admissible**, and **best time found**.

Those distinctions—and the two immediate continuation audits above—will tell you far more than increasing every arc’s step budget indiscriminately.