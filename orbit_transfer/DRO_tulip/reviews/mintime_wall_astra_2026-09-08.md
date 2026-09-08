**I would run two things next: a short, scaled pseudo-arclength continuation from your last converged solution, and a 70 mN direct solve initialized with the entire 75.5 mN trajectory.** The first diagnoses the branch; the second targets the answer. Stop spending runs on time-of-flight exponents.

Three important corrections to your framing:

1. **Your equations essentially rule out loss of all-burn structure on this normal extremal family.**
2. **A fold, a basin problem, and a change in winding are not mutually exclusive explanations.**
3. **A well-conditioned, square shooting system cannot have an interior, nonzero-residual stationary point of ordinary least squares.** If that is genuinely what you observe, something in “well-conditioned,” “stationary,” or “the Jacobian being used” needs revisiting.

## 1. What your existing evidence does—and does not—say

Your independent shooting verification is strong evidence that the converged trajectories are genuine PMP extremals. It does **not** establish their local or global time optimality.

### Loss of all-burn structure: effectively excluded here

With your sign convention and terminal condition,
\[
\dot\lambda_m=-\frac{T}{m^2}\|\lambda_v\|,\qquad
\lambda_m(t_f)=0,
\]
so
\[
\lambda_m(t)=\int_t^{t_f}
\frac{T}{m(s)^2}\|\lambda_v(s)\|\,ds\geq 0.
\]
Consequently,
\[
Q_{\rm mt}(t)=\frac{\|\lambda_v(t)\|}{m(t)}
+\frac{\lambda_m(t)}{c}\geq 0.
\]

For a nontrivial normal extremal, this is strictly positive in the interior. An extended interval with vanishing primer would force the spatial costates to vanish, inconsistent with the normal free-time conditions.

**Therefore, a coast arc is not the expected missing branch under your stated model.** If you measure negative \(Q_{\rm mt}\) on a converged solution, investigate a sign, terminal-condition, integration, or reconstruction error.

A different issue is **an isolated zero or near-zero of \(\|\lambda_v\|\)**. That makes the normalized direction and its derivatives problematic even though it does not imply a coast arc. Print the minimum primer norm and its time.

This conclusion depends on the stated problem: free terminal mass, no additional active constraints, and the normal \(L=1\) formulation.

### Winding: your current argument is not established

“Continuation cannot grow winding” is too strong.

For a **planar curve avoiding the Moon**, with fixed endpoints, the continuously unwrapped **net** angular displacement differs between homotopy classes by integer multiples of \(2\pi\). A continuous family cannot change that integer without encountering the excluded point or leaving the planar class.

But:

- Absolute swept angle is not that invariant.
- A three-dimensional Moon-centered path does not have the same planar topological classification.
- Projected winding can change through the projection axis without a lunar collision.
- A target period, especially for a seven-petal orbit, is not automatically one Moon-centered revolution.

Your \(0.75\to0.77\) variation suggests that the reported metric is not a strict fixed-endpoint planar winding invariant—or that some geometric/numerical detail matters.

**Measure these separately:**

1. Signed, unwrapped Moon-centered projected angle.
2. Total absolute angular sweep.
3. Minimum projected distance to the lunar axis and minimum lunar distance.
4. Ordered crossings of a chosen Moon-centered ray or Poincaré section.

Those distinguish “more looping” from a true planar homotopy change.

Also, **a different branch becoming faster does not make the old PMP root disappear**. A root solver does not reject an extremal because another extremal is better. An optimal-branch crossing alone cannot explain this wall.

## 2. The decisive branch experiment: scaled pseudo-arclength

Yes, I recommend it—but **as a diagnosis and branch-following tool, not a promise of reaching 70 mN**.

If your branch has a minimum-thrust fold at 75.3 mN and turns back toward larger thrust, pseudo-arclength will demonstrate that cleanly. It will not manufacture a root at 70 mN.

### Use your multiple-shooting unknowns, not a condensed single-shooting system

Let
\[
X=[\lambda_0;\,Y_2;\ldots;Y_K;\,t_f],\qquad R(X,T)=0.
\]
Keep \(K=48\) initially. Introduce fixed reference scales:
\[
X=X_{\rm ref}+D_Xu_X,\qquad
T=T_{\rm ref}+s_Tu_T.
\]
Also scale the residual:
\[
F(u)=D_R^{-1}R(X,T).
\]

The scales must make positions, velocities, mass, costates, time, and thrust numerically comparable. Do not use an unscaled Euclidean norm on this vector.

At a converged point \(u_n\), compute an oriented unit tangent
\[
A_n\tau_n=0,\qquad \|\tau_n\|_2=1,
\]
where
\[
A_n=
D_R^{-1}\begin{bmatrix}R_XD_X&R_Ts_T\end{bmatrix}.
\]

Predict:
\[
u_{\rm pred}=u_n+\Delta s\,\tau_n.
\]

Correct by solving
\[
\boxed{
F(u)=0,\qquad
\tau_n^\mathsf T(u-u_{\rm pred})=0.
}
\]

That is the arclength hyperplane you want. Orient the initial tangent toward decreasing thrust and subsequent tangents by positive dot product with the previous tangent.

If you compare different \(K\), give junction variables quadrature-like weights in the arclength metric so doubling the mesh does not double their influence. You can absorb those weights into the coordinates.

### Get \(R_T\) cheaply

For the first experiment, use a central finite difference **holding all shooting unknowns \(X\) fixed**. Check it at two perturbation sizes.

You are differentiating short, independent segment propagations—not a reoptimized solution. This is a reasonable place to use finite differences.

The production alternative is a segment sensitivity:
\[
\dot q=f_yq+f_T,\qquad q(t_{\rm segment,start})=0,
\]
including explicit thrust dependence in the terminal residual.

### What to record along the arc

At every **converged root**, record:

- \(T,t_f\);
- tangent thrust component \(\tau_T\);
- smallest singular values of the scaled \(R_X\);
- smallest singular value of the rectangular augmented matrix \(A\);
- primer minimum and winding diagnostics.

A simple fold has:

- \(\tau_T\to0\), then changes sign;
- \(R_X\) loses one rank;
- \([R_X\ R_T]\) retains full row rank.

With left and right null vectors \(w,v\), the usual additional nondegeneracy checks are
\[
w^\mathsf T R_T\ne0,\qquad
w^\mathsf T R_{XX}[v,v]\ne0.
\]
The second can be estimated directionally if needed.

**Crucially, conditioning at the stalled nonroot does not rule out a fold elsewhere on the root curve.**

## 3. What the different diagnoses would look like

| Hypothesis | Useful discriminating evidence |
|---|---|
| Simple thrust fold | Pseudo-arclength turns in \(T\); fixed-\(T\) Jacobian singular at the turning root; augmented Jacobian regular. |
| Ordinary basin/corrector failure | Arclength reaches a regular root at 75 or 70 mN, and fixed-\(T\) Newton succeeds when initialized there. |
| Different winding family needed | A separately found 70 mN solution has a different geometric/topological label. This does not by itself explain how the old family terminates. |
| Loss of all-burn | Excluded for regular normal extremals under your stated conditions; check primer regularity instead. |
| Loss of local optimality | Negative directions or zero crossings in the correctly constrained second variation—not merely a shooting residual or raw determinant. |

For conjugate-point work, I would **not begin with a raw determinant**. Its magnitude is a poor diagnostic, and the boundary conditions matter: fixed terminal \(r,v\), free terminal mass, and free terminal time.

Use an appropriately bordered Jacobi/second-variation construction, or exploit your direct machinery to examine the reduced Lagrangian Hessian on the linearized feasible subspace. Mesh convergence matters.

A conjugate point or loss of local optimality need not destroy the PMP branch. These are related but different questions, and a fold and a second-order degeneracy can overlap.

## 4. Direct collocation: yes, use it now

There is no compelling reason to insist on indirect continuation if your immediate deliverable is a 70 mN solution and the direct-to-indirect pipeline already works.

Your statement that direct collocation is robust “because it never integrates over a long horizon” is incomplete: your multiple shooting already avoids long individual propagations. The more relevant advantage is that the NLP can reshape the entire trajectory while managing feasibility and the objective, instead of correcting only a particular extremal family.

It remains a local nonlinear method. It can stall or stay in the same family too.

### First direct experiment

Initialize the 70 mN NLP with:

- the **full** converged 75.5 mN state/control history;
- a modest time stretch;
- an exactly consistent mass history for 70 mN;
- endpoints enforced exactly.

Do not reduce the seed to \(\lambda_0\) and a larger \(t_f\).

Run that alongside a geometrically different seed, described below.

### Costate harvesting checks

At this thrust, I would watch:

1. **Scaling and signs.**  
   Recover covectors from the actual discretized Lagrangian. Account for objective scaling, state scaling, defect scaling, and multiplier sign convention. Do not merely read raw IPOPT multipliers as continuous costates.

2. **Time normalization.**  
   With \(t=t_fs\), the defects contain \(t_f\). The free-time stationarity condition must be reconstructed consistently.

3. **Discrete versus continuous adjoints.**  
   Hermite–Simpson multipliers require the appropriate covector mapping. Interpolating raw midpoint or endpoint multipliers is not generally enough.

4. **Artificial active constraints.**  
   Path boxes, altitude limits, time bounds, or mass bounds used “for numerical safety” change the necessary conditions when active. Your current shooting equations then will not match the NLP.

5. **Primer resolution and control representation.**  
   Refine near small primer norm or rapid direction change. Check direction agreement,
   \[
   \alpha+\lambda_v/\|\lambda_v\|\approx0,
   \]
   away from primer zeros.

6. **Continuous validation.**  
   Check \(\lambda_m(t_f)\), \(H(t_f)\), Hamiltonian constancy, costate ODE defects, mass law, and independent propagated terminal errors.

A direct solution that is “feasible” only at collocation nodes is not ready for harvesting. Refine until the continuous defects and multipliers stabilize.

## 5. How to construct a genuinely higher-loop seed

**Increasing \(t_f\) alone does not insert a loop. Insert the geometry.**

Also, “one target period later is another admissible arrival” needs qualification. In autonomous CR3BP, a *ballistic* target-orbit segment returns to the same rotating-frame state after its period. Your spacecraft with nonzero thrust does not simply remain on that ballistic orbit. The appended trajectory is a useful initial guess, not an automatically feasible all-burn transfer.

### My preferred seed set

Build these full trajectories for direct collocation:

**A. Departure-loop seed**

1. Follow the ballistic DRO for one complete period.
2. Append the existing transfer.
3. Set the total time to approximately \(t_{f,\rm old}+\tau_{\rm DRO}\).
4. Rebuild mass using the target thrust over the whole interval.
5. Supply a smooth unit thrust-direction guess and allow the NLP to deform the ballistic loop.

Repeat with two departure periods if inexpensive.

**This is my first loop seed:** your departure period is much shorter than the target period, so it introduces a different itinerary without adding 23 days immediately.

**B. Arrival-loop seed**

1. Follow the existing transfer to the fixed target phase.
2. Append one full ballistic target period.
3. End at exactly the same target state.
4. Again rebuild mass and let the NLP deform the appended orbit under thrust.

**C. Forced itinerary seed**

If A and B immediately collapse back to the short solution, temporarily impose a small number of interior gates or section crossings that force the desired excursion. Solve that auxiliary problem, then continuously remove the gates.

The gate multipliers must vanish or the gates must be removed **before harvesting for your original shooting problem**.

These constructions change the guess’s itinerary. They do not guarantee the final solution preserves it. Track the crossing sequence and winding measure during optimization.

For a truly planar collision-free problem, a different homotopy class can provide a stronger separation. For unrestricted spatial motion, treat “higher winding” as a geometric itinerary, not an inviolable branch label.

## 6. What an \(O(1)\) stall actually tells you

By itself: very little. Trust-region methods are designed not to diverge.

But there is a sharp consistency test. For the **scaled residual actually minimized**,
\[
\phi(u)=\frac12\|F(u)\|^2,\qquad
\nabla\phi=J_F^\mathsf TF.
\]
For a square, nonsingular \(J_F\),
\[
J_F^\mathsf TF=0\quad\Longrightarrow\quad F=0.
\]

So a genuine interior nonroot least-squares stationary point requires singularity. A nonroot stall with a well-conditioned Jacobian instead suggests termination/scaling issues, constraints, inaccurate derivatives, integration noise, or failure to take a useful step.

### Print this in one failed run

- Raw and scaled residual norms, separated into:
  - each segment’s continuity block;
  - terminal position;
  - terminal velocity;
  - \(\lambda_m(t_f)\);
  - \(H(t_f)\).
- \(\|J_F^\mathsf TF\|\), with its own meaningful scaling.
- Largest and smallest singular values of **the actual scaled Jacobian**.
- Trust radius, accepted step size, termination reason, and actual/predicted reduction ratio.
- Newton correction \(\delta u\), plus:
  \[
  \|J_F\delta u+F\|,\qquad \|\delta u\|.
  \]
- \(\phi(u+\eta\delta u)\) for several decreasing \(\eta\).
- Directional Jacobian checks using finite differences at two step sizes—especially along the proposed Newton step and smallest right singular vector.
- Minimum mass, minimum primer norm, minimum lunar distance, and segment integration failures.
- Projection \(w_{\min}^{\mathsf T}F\), indicating how much residual lies in the hardest-to-correct left-singular direction.

For a correct smooth Jacobian, a Newton direction satisfies
\[
\nabla\phi^\mathsf T\delta u=-\|F\|^2.
\]
Thus sufficiently small steps should decrease the merit function. Failure of that test is a strong implementation/numerical clue, not evidence of winding.

## Ranked by payoff per effort

1. **One instrumented stall run:** scaled gradient/SVD, directional derivative checks, Newton merit-function probe. Cheap and potentially decisive.
2. **Short pseudo-arclength run from 75.5 mN:** determine whether the branch turns, continues, or encounters primer/numerical pathology.
3. **70 mN direct solve from the complete 75.5 mN trajectory:** run in parallel with item 2.
4. **Direct solve with one inserted DRO period**, then an arrival-period seed if necessary.
5. **Second-order optimality analysis** once you have the competing branches or a clearly identified degeneracy.
6. **More segments and more scalar \(t_f\) guesses:** low expected payoff now.

**My working diagnosis is “an unresolved branch-following or corrector failure,” not “proven winding wall.”** Your next experiment should measure the solution curve itself. Pseudo-arclength answers that question; direct collocation attacks the 70 mN deliverable.