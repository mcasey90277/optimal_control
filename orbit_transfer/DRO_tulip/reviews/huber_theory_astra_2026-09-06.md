## Bottom line

1. **Claim 1 is not yet a sufficient-certificate claim.** The endpoint manifold, free-time reduction, endpoint-map regularity, and meaning of normality must be audited. Integrator-resolution determinant sampling is not a proof of absence of conjugate/focal points.
2. **Claim 2 mischaracterizes continuous clipping.** At regular bang–interior junctions, branchwise differentiation can give the correct first variational flow, without a nontrivial saltation matrix or an additional independent switching-time Hessian. Sufficiency must address the full constrained second variation.
3. **Claim 3 is a plausible diagnosis promoted to a false theorem.** Neither grazing nor exclusion of a fold is established. The claim that no arclength formulation can pass is wrong.
4. **Claim 4 provides useful numerical evidence, not a family-independent floor or established bang-bang structure.** Switching-time refinement is sensible, but its convergence does not itself certify optimality.
5. **Claim 5 is correct for \(0<p<1,\ \delta>0,\ s\in[0,1]\).** Its limits require explicit treatment.

---

## 1. Normal min-time, free terminal mass, no conjugate time

### Verdict: the underlying sufficient-condition strategy is legitimate; the stated hypotheses and proposed numerical gates are insufficient.

### A. Free terminal mass changes the boundary-value second variation

Your terminal constraint is the manifold
\[
M_f=\{(r_f,v_f,m):m>0\},
\]
not a fixed point in seven-dimensional state space. Transversality is
\[
\lambda(t_f)\in N^*M_f,
\qquad\text{hence}\qquad \lambda_m(t_f)=0.
\]

The Jacobi/accessory problem must use this endpoint manifold and allow variation of terminal time. In a physical-time linearization, the terminal state variation includes
\[
\delta x(t_f)+f(t_f)\,\delta t_f.
\]
Its \(r,v\) components must vanish; its mass component is not prescribed. Linearized transversality must also be included.

**A point-to-point conjugate test with all seven terminal state variations fixed is not the appropriate sufficient test.** It tests a smaller class of competitors. In geometric language, one needs the conjugate/focal or Lagrangian-intersection construction appropriate to the actual endpoint conditions. Terminology varies; the boundary conditions do not.

For a curved endpoint manifold there is also an endpoint-curvature contribution to the second variation. Here \(M_f\) is affine, so that particular contribution vanishes. The free mass variation and mixed boundary conditions remain.

There is an important simplification specific to case A:
\[
m(t)=m_0-\frac{T}{c}t.
\]
Thus mass can be eliminated, yielding a six-state, time-dependent problem. At a varied terminal time,
\[
\delta m_f=-\frac{T}{c}\delta t_f.
\]
Free terminal mass therefore does **not** introduce an independently controllable terminal coordinate in the all-burn model. But retaining mass can introduce structural rank deficiencies unless the Jacobi construction accounts for this. A generic seven-state determinant is particularly suspect.

**Required audit:** write down the exact reduced exponential map or accessory boundary-value problem and show that its kernel is precisely the admissible null space of the second variation. “Free-time quotient, as in HamPath” is not that derivation.

### B. Strong Legendre plus normality is not the whole theorem

Your direction Legendre calculation is correct:
\[
H_{\alpha\alpha}\big|_{T_\alpha S^2}
=\frac{T}{m}\|\lambda_v\|\,I
\]
at the minimizing direction, assuming \(T>0,\ m>0\). On a compact interval, continuity and absence of primer zeros provide a uniform positive bound.

However, the geometric conjugate-time sufficiency results also require the regularity assumptions supporting their construction: smooth dynamics away from collisions and zero mass, a smooth regular minimizing feedback, and the appropriate rank/corank or strong-normality assumptions on the endpoint map. A Jacobi determinant is meaningful only after those issues and the structural null directions have been resolved.

The safe theorem-level statement is:

> Under the hypotheses of the applicable free-time, endpoint-manifold second-order theorem, positivity/coercivity of the relevant second variation—equivalently, under its Jacobi/Morse-index hypotheses, absence of the relevant conjugate/focal degeneracy—together with the required Hamiltonian minimization conditions yields a strict strong local minimum.

In this problem the unique global direction minimizer helps with the strong/Weierstrass part. But **Legendre and a determinant check do not replace the other hypotheses**.

Likewise, “an interior conjugate time refutes minimality” needs the relevant index theorem and regularity hypotheses. A zero of an incorrectly reduced determinant, a redundant mass direction, or a numerically near-zero determinant does not refute anything.

### C. Normality: distinguish two meanings

If an **exact** feasible extremal admits a PMP multiplier with \(\lambda_0=1\), satisfying the adjoint equation, minimum condition and transversality, then it has a normal lift. That is sufficient to call it normal in the usual nonexclusive sense.

But:

- Finite shooting variables alone establish nothing.
- Numerical acceptance establishes an approximate normal witness, not an exact one.
- Existence of a normal lift does **not** exclude a different abnormal lift for the same state-control trajectory.
- Testing the same costate direction after dropping the “1” excludes only that particular proposed abnormal multiplier.
- If the conjugate theorem requires strict normality, strong normality, or a corank-one endpoint map, a normalized normal BVP does not automatically establish it.

Also audit “costate scaling quotient.” Once \(\lambda_0=1\), scaling \(\lambda\) alone is not a PMP symmetry:
\[
H=1+h,\quad H=0
\quad\not\Rightarrow\quad
1+a h=0.
\]
Homogeneous geometric constructions can use projectivization or an energy normalization, but the equivalence to your normalized free-time BVP must be demonstrated.

### D. Integrator-resolution sampling is not a sufficient certificate

It can miss:

- two crossings between samples;
- an even-multiplicity zero with no sign change;
- nearly singular blocks obscured by determinant scaling;
- a terminal focal degeneracy;
- integration error in the Jacobi flow.

An adaptive ODE mesh controls the estimated local error in the integrated variables. It does not automatically isolate all zeros of a determinant.

For a **numerical sufficiency assessment**, use a geometrically appropriate rank/index method, singular values or Riccati charts, crossing localization, and demonstrated integration/mesh convergence. For a **rigorous certificate**, add validated enclosures and an exact-solution existence argument, or another genuinely rigorous error-control framework.

A sampled \(\min\|\lambda_v\|\) gate has the same logical issue: it needs a between-sample bound to prove nonvanishing.

### E. Is \(C^0\)-local minimality what the catalog needs?

It is valuable: strong local optimality allows competitors with nearby state trajectories but controls that need not be uniformly close. For variable final time, specify the topology using a common time parametrization and proximity of final times.

It does **not** imply:

- global optimality;
- optimality against a distant transfer branch;
- a useful-sized basin of local optimality;
- robust optimality under model or endpoint perturbations.

Also, the certificate concerns the **all-burn problem**. If the catalog advertises minimum time with throttle \(s\in[0,1]\), you must establish the all-burn reduction or test optimality against throttle variations too. With mass depletion, that equivalence should not simply be assumed.

**References:** Bonnard–Caillau–Trélat’s 2007 second-order/conjugate-point treatment; Agrachev–Sachkov, *Control Theory from the Geometric Viewpoint* (2004), for regular extremals, exponential maps and conjugate points; Vinter, *Optimal Control* (2000), for endpoint constraints and strong optimality. Cite the actual theorem and verify its hypotheses against your boundary conditions—not merely the software implementing one specialization.

---

## 2. Smoothed fuel with clipped throttle and branchwise AD

### Verdict: partly right about the critical cone; wrong to identify every bound–interior junction with a bang-bang switching-time problem.

For the continuous clipped laws—your eps family and huberc with \(\delta>0\)—the minimized canonical field \(F(y)\) is generally **continuous and piecewise smooth**.

At a transversal branch boundary, the saltation formula is
\[
\Xi
=I+\frac{(F^+-F^-)\,n^\top}{n^\top F^-}.
\]
Because \(F^+=F^-\) there,
\[
\Xi=I.
\]

Consequently, integrating
\[
\dot\Phi=F_y(y(t))\Phi
\]
using the derivative of the branch actually occupied by the nominal trajectory, and matching \(\Phi\) continuously, can give the **correct derivative of the full flow**. Switching times do move under perturbation; their first-order effects cancel at a continuous junction.

This is not merely differentiation of an artificially restricted problem with junction times fixed.

### What “frozen active set” can legitimately mean

On an open bound-active arc with strict Hamiltonian complementarity, a critical first-order variation has
\[
\delta s=0
\quad\text{a.e. on that arc}.
\]
That is the proper active-set restriction in the critical cone. Isolated regular junctions have measure zero and do not invalidate it.

But strict complementarity does not hold *at* the junction. The useful hypotheses are:

- strict sign in the interiors of active arcs;
- isolated transversal threshold crossings;
- positive interior throttle curvature;
- regularity of the other control components and endpoint constraints;
- treatment of degeneracies such as the irrelevant direction \(\alpha\) on a coast.

Distinguish:

1. **Freezing the branch for differentiation along the nominal flow:** often correct.
2. **Forbidding perturbed trajectories from changing branch or moving junctions:** not the full perturbation problem.

### Is a separate switching-time Hessian necessary?

**Not as a universal statement.** For continuous saturation of a strictly convex control cost, coercivity of the full constrained second variation can be analyzed with critical-cone/accessory-system or Riccati methods. A separate bang-bang switching-time Hessian is not automatically an additional missing term.

Such finite-dimensional switching conditions are central for genuine discontinuous bang-bang controls, and can appear in mixed-arc formulations. But a theorem must match the actual junction type.

Your original Huber family is different: despite \(p>0\), it is linear above \(s=p\), the minimizer is nonunique at \(Q=1\), and the selected throttle jumps. It cannot be grouped with continuous eps clipping merely because its parameter is positive.

Finally, an STM is not itself a second-order test. You still need the correct endpoint boundary conditions, critical directions and Jacobi/index interpretation. A fixed-seven-state endpoint test again omits the actual free-mass terminal condition.

**Correct replacement:** branchwise AD can correctly generate the variational flow for continuous clipped feedback under regular crossings. Whether a resulting Jacobi test is necessary or sufficient depends on the accessory boundary-value problem and the constrained second-order theorem being implemented.

**References:** Vinter (2000); Bonnans–Shapiro, *Perturbation Analysis of Optimization Problems* (2000), for critical cones and complementarity; Osmolovskii–Maurer, *Applications of Second Order Optimality Conditions in Optimal Control Engineering* (2012), for the specific bang-bang and mixed-structure second-order frameworks. Do not use “Maurer–Osmolovskii–Büskens conditions” as an undifferentiated theorem covering all clipped laws.

---

## 3. Huber walls as grazing bifurcations that no continuation can pass

### Verdict: grazing is plausible but unproved; several claimed consequences are false.

### A. What would establish a generic grazing event?

With
\[
h(y,p)=Q(y,p)-1,
\]
you need a solution on the branch satisfying
\[
h=0,\qquad \dot h=0,
\]
and a suitable nondegeneracy condition such as
\[
\ddot h\ne0,
\]
plus nonzero unfolding in the continuation direction. For discontinuous fields, specify the relevant one-sided derivatives and solution convention.

Neither observation establishes this:

- \(|\dot Q|=0.045\) is small relative to selected controls, but nonzero.
- \(Q_{\max}=0.9907\) is below threshold. Without its derivative along the solution branch, you do not know whether it will cross.

These are **grazing risk indicators**, not detection of a bifurcation. The magnitude of \(\dot Q\) also depends on time normalization.

### B. Actual local behavior: often square-root, not universally a corner

A generic local contact model is
\[
h(t,\mu)=a\mu-b(t-t_*)^2+\text{higher-order terms}.
\]
On the side where two roots exist,
\[
t_\pm(\mu)=t_*\pm\sqrt{a\mu/b}+\cdots.
\]
Thus switching-time sensitivity can diverge as \(\mu^{-1/2}\).

For a discontinuous vector field, the resulting endpoint correction can generically scale like
\[
C\sqrt{\mu_+}.
\]
That is continuous but not locally Lipschitz, rather than a finite-slope \(|\mu|\)-type corner. However:

- leading coefficients can vanish;
- composition of the entry and exit maps matters;
- endpoint projections can cancel leading terms;
- sliding, nonuniqueness or nonshrinking excursions change the analysis.

There is **no universal regularity conclusion from the saltation denominator alone**. The numerator and geometry also matter.

Your PMP switch is not an arbitrary relay. After minimizing direction, its canonical Hamiltonian has the form \(H_0+\phi(Q)\). At \(Q=1\), the vector-field jump is proportional to the symplectic gradient of \(Q\), hence
\[
\nabla Q^\top(F^+-F^-)=0.
\]
The normal velocities agree on the two sides. This structure is relevant to the grazing analysis; generic Filippov slogans are not a substitute for it.

For continuous clipping, a grazing activation often gives a correction of order \(\mu_+^{3/2}\): the excursion duration is \(O(\sqrt\mu)\), while the vector-field difference is \(O(\mu)\). The endpoint map may therefore be \(C^1\) but not \(C^2\). “Continuous control means everything is smooth” is also too strong.

### C. “No arclength continuation can pass” is wrong

Classical smooth pseudo-arclength continuation may lose its usual guarantees at a genuinely nonsmooth point. That does not mean all continuation formulations fail.

Possible formulations include:

- explicit switching times and arc durations;
- a variable \(w=\sqrt{\mu}\);
- active-structure changes at zero arc duration;
- nonsmooth or complementarity formulations;
- one-sided branch tracking and restart beyond the event.

Even a square-root graph can become smooth under a different parameter. A corner does not imply that the solution set cannot be traversed.

The valid statement is:

> The current smooth Newton/p-continuation formulation may be unsuitable near a change in switch structure.

The conclusion that the family **must** change is unsupported.

### D. Flat \(\operatorname{cond}(J)\) does not exclude a fold

At a simple smooth fold, the exact \(R_z\) becomes singular along exact solutions approaching the fold. That prediction is useful—but the reported experiments do not meet those conditions convincingly.

You have:

- finite-distance samples;
- nonzero residuals at failed iterates;
- potentially inaccurate event sensitivities;
- scaling-dependent condition numbers already around \(10^8\).

A nearly constant condition number at stalled iterates is not a reliable exclusion test. At most it says you did not observe the expected divergence in the matrix you computed.

Use scaled singular values, verify derivatives independently, track the augmented continuation tangent, and test whether its \(p\)-component approaches zero.

### E. Other events?

Two switches colliding and a burn arc shrinking to zero are often **the same generic grazing event**, not competing explanations. Other possibilities include:

- a switch entering or leaving through \(t=0\) or \(t_f\);
- a different solution-branch fold;
- higher-order contact or an incipient singular arc;
- event-order changes;
- integration/event-location or Jacobian errors.

An extremum near the start of the trajectory deserves an explicit endpoint-event check.

**References:** Nordmark, “Non-periodic motion caused by grazing incidence in an impact oscillator” (1991), for the square-root grazing paradigm; di Bernardo–Budd–Champneys–Kowalczyk, *Piecewise-smooth Dynamical Systems* (2008), for discontinuity maps and grazing; Kuznetsov, *Elements of Applied Bifurcation Theory*, for folds and continuation. Impact results illustrate a mechanism; they are not automatically the theorem for your canonical PMP switch.

---

## 4. Continuous ramp, universal sharpening floor, switching-time endgame

### Verdict: promising engineering, overgeneralized interpretation.

### A. The source material already contradicts two headline claims

The fixed-\(\delta\) p-walk did **not** have zero failures on every cell: \((1,2)\) had 15.

Nor are the reported terminal masses uniformly shown to have converged to \(10^{-6}\). Some inter-family differences remain around \(10^{-5}\). Small movement along the last few rungs is evidence of stabilization, not an error bound relative to the bang-bang optimum.

Moreover, endpoint mass is one scalar. Switch times and trajectories can differ materially while final mass is nearly unchanged.

### B. No family-independent floor has been established

Three cells using essentially the same solver, integration settings, tolerances, scaling and step policy establish a **configuration-dependent empirical floor**.

A ramp width in \(Q\)-space is not a temporal resolution:
\[
\Delta t_{\rm ramp}\approx
\frac{\delta}{|\dot Q(t_s)|}
\]
for a transversal crossing. For huberc,
\[
\frac{ds}{dQ}=\frac{1-p}{\delta},
\]
but the effective variational coefficients also involve \(Q_y\), the vector fields and the trajectory.

The floor can depend on:

- ramp slope and crossing speed;
- crossing curvature near grazing;
- state and costate scaling;
- integrator tolerance and maximum step;
- accuracy of STM integration through derivative discontinuities;
- shooting-node placement and segment length;
- unstable propagation and nonnormal amplification;
- Newton globalization, continuation step size and residual tolerances;
- genuine degeneracy of the limiting BVP.

A 12-segment multiple-shooting formulation is **not inherently unable to resolve switches narrower than some universal \(\delta\)**. Each segment contains an adaptive continuous integration. More segments can improve propagation conditioning, but segment count alone does not set switch resolution.

Test the hypothesis by varying tolerances, maximum step, scaling, segment count and event-aligned nodes separately. Compare actual \(Q\)-band widths and slopes, not raw family parameters.

Also, at fixed \(p=0.001\), taking \(\delta\to0\) gives the original Huber-type jump, **not exact bang-bang fuel**: the low branch still permits \(s=pQ\).

### C. Switching-time Newton is a sensible candidate, not a proved endgame

Before fixing a switching structure, establish evidence for:

- finitely many separated limiting switches;
- transversal switching;
- strict switching-function sign between switches;
- no singular arc or unresolved small arc;
- correct coast/burn ordering and endpoint conditions;
- regularity of the reduced endpoint constraints.

Then formulate the exact limiting bang-bang BVP using arc durations or switch times, possibly with multiple shooting inside arcs.

But direction \(\alpha(t)\) remains a continuously varying control on burns. A throttle switching-time Hessian alone does not automatically account for all direction variations. Those must be included in the accessory problem or rigorously eliminated, including their coupling to switch-time variations.

### D. Newton convergence is not a strict second-order certificate

Convergence proves, numerically, stationarity/root finding. Certification additionally requires:

- the full PMP conditions and regular switching;
- endpoint constraint qualification;
- positive definiteness/coercivity of the appropriate **Lagrangian** second variation on admissible critical directions;
- the hypotheses transferring reduced-problem positivity to the asserted local minimum of the original control problem.

In particular, test the reduced Hessian on the kernel of the endpoint-constraint Jacobian, not the raw objective Hessian.

The survey’s broader assertion that fuel linearity makes every NLP Hessian intrinsically incapable of detecting strictness is also false. The Lagrangian includes nonlinear dynamics and endpoint constraints; its reduced Hessian need not vanish merely because \(L(s)=s\). Discretization artifacts and inadequately represented switches are real issues, but that proposed impossibility argument is invalid.

**References:** Maurer–Büskens–Kim–Kaya, “Optimization methods for the verification of second order sufficient conditions for bang-bang controls,” *Optimal Control Applications and Methods* 26 (2005), 129–156; Osmolovskii–Maurer (2012). Match the theorem to the mixed continuously controlled direction/bang-bang throttle problem.

---

## 5. Exact argmin for huberc

### Verdict: correct for \(0<p<1,\ \delta>0\), with the control domain explicitly restricted to \([0,1]\).

Define
\[
L(s)=
\begin{cases}
\dfrac{s^2}{2p},&0\le s\le p,\\[2mm]
\dfrac p2+(s-p)+
\dfrac{\delta(s-p)^2}{2(1-p)},&p<s\le1.
\end{cases}
\]
Then
\[
L'(s)=
\begin{cases}
s/p,&0<s<p,\\[1mm]
1+\dfrac{\delta(s-p)}{1-p},&p<s<1.
\end{cases}
\]
Value and first derivative agree at \(p\), and
\[
L''=
\begin{cases}
1/p,\\
\delta/(1-p),
\end{cases}
\]
are positive. Thus \(L\) is \(C^1\), piecewise quadratic and strictly convex—indeed strongly convex on \([0,1]\).

The necessary and sufficient optimality condition is
\[
0\in L'(s)-Q+N_{[0,1]}(s).
\]
Solving gives
\[
s(Q)=
\begin{cases}
0,&Q\le0,\\
pQ,&0<Q<1,\\
p+\dfrac{(1-p)(Q-1)}{\delta},&1\le Q<1+\delta,\\
1,&Q\ge1+\delta.
\end{cases}
\]
The minimizer is unique, including at both knees.

### Qualifications

- **“Saturated \(L=L(1)\)” must not mean an unconstrained constant extension for \(s>1\).** That extension is not convex at \(1\). The convex extended-value formulation is \(L+I_{[0,1]}\).
- **At \(p=1\):** the ramp interval disappears. Define \(L=s^2/2\) and \(s=\operatorname{clip}(Q,0,1)\) explicitly. Do not rely on evaluating a formula containing \(1/(1-p)\).
- **At \(p=0\):** the core formula is singular, but the limit is
  \[
  L(s)=s+\frac{\delta}{2}s^2,\qquad
  s(Q)=\operatorname{clip}\!\left(\frac{Q-1}{\delta},0,1\right).
  \]
- **At \(\delta=0,\ 0<p<1\):** strict convexity is lost on \([p,1]\). At \(Q=1\),
  \[
  \arg\min_{[0,1]}(L-sQ)=[p,1].
  \]
  Your original jump law chooses one minimizer there. For fixed \(Q=1\), the positive-\(\delta\) family selects \(p\); along \(Q=1+\theta\delta\), its limiting selection can be any point of \([p,1]\).

### Is it a known regularization?

It is naturally a **convex piecewise-linear-quadratic modification of the Huber penalty**, made strictly convex by adding quadratic curvature to its linear tail on a bounded interval. Equivalently,
\[
s(Q)=\nabla\big(L+I_{[0,1]}\big)^*(Q).
\]
This is standard convex-conjugate construction of a continuous monotone feedback from a strongly convex penalty.

I would not assign the exact two-parameter formula a standard proper name without a specific source. Huber’s original loss supplies the core/linear-tail construction; Rockafellar, *Convex Analysis* (1970), supplies the conjugate/subgradient argument; Rockafellar–Wets, *Variational Analysis* (1998), treats the PLQ framework.

**Build on the exact convexity result. Treat the claimed universal floor, diagnosed grazing bifurcations, and catalog-wide strict-minimum certificates as hypotheses still requiring validation.**