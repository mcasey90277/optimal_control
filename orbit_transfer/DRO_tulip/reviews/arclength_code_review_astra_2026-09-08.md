## Bottom line

**I would not yet call the 71.99 mN feature a simple fold.** Your results establish that thrust stepping failed and that continuation found additional accurate roots. They do **not** establish a finite turning point at \(t_f=18.120\) d.

My main conclusions are:

1. **Your pseudo-arclength corrector is mathematically correct. Your tangent calculation is the wrong one to rely on near a fold.** Replace it with the null vector of the full augmented Jacobian.
2. **An ordinary, finite, well-resolved simple fold should not cause pseudo-arclength to crawl indefinitely.** Investigate derivative accuracy, actual accepted step lengths, and—especially—costate growth.
3. **A normal-costate chart approaching infinity is a serious alternative explanation.** Your earlier costate divergence makes this more than a generic caveat.
4. **“The minimum thrust for this transfer geometry is 72 mN” is already contradicted by your 70 mN solution**, assuming it uses exactly the same endpoints and admissible-control model.
5. The 26.436 d solution is a well-validated candidate local minimum. **Nothing supplied certifies global minimum time.** An overnight search can substantially strengthen the evidence, but normally cannot establish globality.

This is a source review and mathematical assessment; I have not executed the code or audited the underlying `pumpkyn` dynamics or conjugate-test implementation.

# A. Implementation review

## 1. Tangent: replace lines 123–124

You currently compute
\[
v=-R_x^{-1}R_\theta,\qquad
\tau=\frac{(v,1)}{\|(v,1)\|},
\]
where \(\theta=T/s_T\).

Away from singular \(R_x\), this is mathematically a valid tangent. At a simple fold, however, it is an **invalid coordinate chart for calculating that tangent**: \(dX/dT\) diverges while the arclength tangent remains finite.

### Minimal replacement

You already perform an SVD of the augmented matrix. Use its right null vector:

```matlab
B = [Jx Rt];                         % n-by-(n+1)

[~, SA, VA] = svd(B);                % FULL SVD, not 'econ'
sminAug = SA(n,n);
tauNew = VA(:,end);                  % the structural right null vector

if isempty(tau)
    if tauNew(end) > 0
        tauNew = -tauNew;
    end
elseif tauNew' * tau < 0
    tauNew = -tauNew;
end

tau = tauNew;
```

**The full-SVD detail matters:** an economy SVD of this wide matrix can omit the extra right-nullspace vector. At your \(n=330\), a full SVD is entirely reasonable, particularly compared with the propagation cost.

Alternatively, after initialization, a bordered tangent solve is standard:
\[
\begin{bmatrix}
B\\
\tau_{\rm old}^{T}
\end{bmatrix}q
=
\begin{bmatrix}0\\1\end{bmatrix},
\qquad \tau=q/\|q\|.
\]
It remains nonsingular when \(B\) has full row rank and the previous tangent is not orthogonal to the current one.

### Does your current formula force the thrust component to keep its sign?

**No.** The initialization chooses decreasing thrust, but your subsequent dot-product orientation permits \(\tau_T\) to change sign. The formula can, in exact arithmetic, produce correctly oriented tangents on opposite sides of a fold.

Thus:

- It is numerically fragile and undefined at the exact fold.
- It **does not algebraically prohibit turning the corner**.
- Normalization can sometimes preserve an accurate direction despite a very large \(v\). You cannot infer tangent corruption solely from a small \(\sigma_{\min}(R_x)\).

Log both:
\[
\frac{\|B\tau\|}{\|B\|\|\tau\|},
\qquad
\cos^{-1}(\tau_k^T\tau_{k-1}).
\]
The first tests the computed null vector against the computed Jacobian—not against errors in that Jacobian.

## 2. Corrector and Jacobian: the formulation is correct

Your equations are
\[
F(w)=
\begin{bmatrix}
R(w)\\
\tau_k^T(w-w_{\rm pred})
\end{bmatrix}=0,
\]
with Jacobian
\[
F_w=
\begin{bmatrix}
R_x&R_\theta\\
\tau_k^T
\end{bmatrix}.
\]

Lines 160–163 implement exactly that.

At a simple fold, \(R_x\) is singular but this bordered Jacobian is nonsingular, provided the tangent is appropriate. **You already have the essential Keller-style pseudo-arclength corrector.** You do not need to add deflation or a special fold-crossing equation to make an ordinary fold traversable.

Keep the hyperplane tangent fixed during each corrector. Do not differentiate it with respect to the Newton iterate.

### What I would change

**Add globalization.** Plain full-step Newton can leave the local convergence basin. Backtracking on a sensibly scaled augmented residual is a straightforward improvement. Treat invalid propagations as failed trial steps, not as equations with a meaningful placeholder Jacobian.

**Check the last Newton update.** Your loop does not test convergence after the update on iteration `nMax`. A point that converges on that update is reported as a failure and causes unnecessary halving.

**Require coordinate accuracy as well as residual accuracy.** A small residual alone is not enough when the bordered system is poorly conditioned. Monitor the final scaled Newton correction or estimated correction, with separate tolerances for the BVP and arclength equation.

For perspective, a smallest singular value of \(1.5\times10^{-7}\) can amplify a residual component of \(10^{-9}\) into a correction of order \(7\times10^{-3}\) in the corresponding coordinates. This is only a linearized illustration, but it is much larger than `dsMin = 1e-4`.

Also reject nonfinite `Jxi` and `Rti`, not just nonfinite `Ri`.

## 3. Step control: crude, but not intrinsically fold-blocking

Halving after failure is standard. A tangent predictor followed by your corrector has the correct near-fold behavior; it does **not** need a special thrust-reversal predictor.

The problem is that every success increases the step by 30%, whether convergence took one iteration or twelve. That tends to produce repeated overshoot–halve cycles.

Use actual difficulty and geometry, for example:
\[
ds_{\rm next}
=
\operatorname{clip}\!\left[
ds_{\rm used}
\sqrt{\frac{N_{\rm target}}{\max(N_{\rm Newton},1)}}
\right],
\]
with bounded change factors, plus a tangent-angle or predictor-correction-size limit. A target of roughly four Newton iterations is a reasonable starting choice, not a universal optimum.

Record:

- `dsUsed`: the predictor length that actually succeeded;
- `dsNext`: the next proposed length;
- retry count and Newton iterations;
- \(\|w_{k+1}-w_k\|\);
- predictor correction \(\|w_{k+1}-w_{\rm pred}\|\);
- tangent angle;
- conditioning of the bordered corrector.

**Your current `A.ds` is not the accepted step length.** After success you enlarge `ds`, and the next stored row contains that enlarged proposal. This obscures whether the arc is actually collapsing.

Two very different situations can look like “thrust has stopped moving”:

- **Actual arclength steps collapse:** investigate numerical accuracy, curvature, globalization, or loss of augmented rank.
- **Arclength steps remain substantial but mostly change costates:** the curve may be escaping to infinity in your coordinates.

The supplied table does not distinguish them.

## 4. Scaling: usable initially, but neither mesh-independent nor robust to costate growth

The transformation
```matlab
Dx = max(abs(p0), 1e-2);
Jx = Jp .* Dx(:)';
```
is internally consistent. Scaling thrust by `TN0` is also mathematically valid.

However:

### The metric depends on the mesh

Repeating junction variables in a Euclidean norm means their aggregate contribution grows roughly as \(\sqrt K\) for similar trajectory variations. Per-component normalization does **not** remove that dependence, despite the comment in the header.

A better trajectory metric is of the form
\[
\|\delta w\|^2 =
a_\lambda\|S_\lambda^{-1}\delta\lambda_0\|^2
+\sum_{k=2}^{K}\omega_k
  \|S_Y^{-1}\delta Y_k\|^2
+a_f(\delta t_f/s_f)^2
+a_T(\delta T/s_T)^2,
\]
where \(\omega_k\) are normalized mesh/quadrature weights and the block weights are explicit choices.

### Initial component magnitudes are not always good characteristic scales

A component initially near zero gets a scale of \(10^{-2}\), even if it naturally becomes order one later. Different components of the same physical block can receive very different weights for accidental reasons.

Use characteristic block scales, or representative trajectory magnitudes with physically motivated floors.

### Fixed initial scales cannot handle unlimited multiplier growth gracefully

If costates grow by orders of magnitude, your arclength can become overwhelmingly costate arclength. Small physical changes may then require huge travel in \(w\).

**A fixed nonsingular scaling cannot create or destroy a finite mathematical fold.** It can dramatically affect practical progress, conditioning, and reported singular values.

Do not silently update scales inside Newton. If you rescale between continuation windows, transform/recompute the tangent and restart the local metric consistently.

## 5. Finite-difference \(R_T\): conceptually correct, accuracy not demonstrated

You correctly differentiate at fixed unknowns:
\[
R_\theta=s_T\,\frac{\partial R}{\partial T_N}.
\]
The thrust nondimensionalization is automatically included by `mkRes`, and the finite difference also includes explicit thrust dependence in the terminal Hamiltonian.

The concern is numerical resolution. Near 72 mN:
\[
h_T\approx 7.2\times10^{-8}\ {\rm N},
\qquad h_\theta\approx10^{-6}.
\]
Residual-evaluation noise of order \(10^{-12}\) can therefore produce derivative noise around \(10^{-6}\) in scaled coordinates. Correlation between evaluations may improve this, but **you have not established it**, and your smallest augmented singular value is already below that illustrative noise level.

Run a step-size study over several scaled perturbations—say \(10^{-3}\) through \(10^{-6}\)—with tightened propagation tolerances. Look for a stable derivative plateau and stable tangent/fold diagnostics.

There is an additional subtlety: your finite differences request only the residual, whereas Newton requests residual plus STM. Those propagation modes can use different adaptive integration histories. Check that their residuals agree to the accuracy required here.

### Preferred production fix

Integrate thrust sensitivities on each segment:
\[
\dot S=f_Y S+f_T,\qquad S(0)=0.
\]
For the final segment:
\[
R_T^{\rm terminal}=g_Y S+g_T.
\]
Include the explicit \(H_T\) contribution in \(g_T\). This removes subtraction noise and makes the augmented Jacobian much more credible.

Also make derivative evaluations explicitly reject failed propagations. If both perturbed calls return the same rejection residual, your current code can report a spurious zero derivative.

## 6. Engine and binding: the main assembly looks right

For the stated autonomous dynamics, fixed endpoints, and fixed normalized segment fractions:

- Packing and continuity blocks are consistent.
- `J(rows,n) += Fh*dsg(k)` is the correct segment-duration derivative.
- `dgdy*Fh*dsg(k)` is the correct terminal time column.
- The eight terminal conditions match seven initial costates plus free final time.
- Holding junction unknowns fixed when differentiating thrust is correct.

I do not see a missing time-column term for the problem as stated.

Important caveats:

### Expensive residual-factory side effect

Every `mkRes(T)` invokes `assembleOnly`, which first propagates the **original seed with STMs**, just to return a residual handle. You then propagate the requested iterate separately.

That is a substantial avoidable cost, especially inside finite differences and Newton. Expose a factory that returns the closures without evaluating the original seed, or accept thrust as an argument to the residual.

### Hidden guards must be logged

The residual contains:

- time limits relative to the original seed;
- `pMax = 1e7`;
- propagation-collapse rejection.

The time upper bound seems unlikely to explain 18 d, but a costate bound or propagation issue could. Return a failure reason. Do not treat all failures as evidence of difficult geometry.

### Do not use `fixedTf=true` in this binding to switch continuation coordinates

`ms_tfmin` still supplies eight terminal conditions, including \(H=0\), and its output packing assumes trailing \(t_f\). Its binding has not been adapted to the engine’s fixed-time option.

To continue using \(t_f\) as the coordinate, retain the original extremal equations, make thrust unknown, and add
\[
t_f-t_{\rm target}=0.
\]
That is different from formulating a genuinely fixed-time optimal-control problem.

# The critical alternative: a fold at finite coordinates, or a normality limit at infinity?

This is the first scientific diagnosis I would make after replacing the tangent.

Your Hamiltonian uses the normal convention
\[
H=1+\lambda^T f.
\]
That fixes the objective multiplier to one. It is a valid chart for normal extremals, but it can become singular as an extremal approaches an **abnormal multiplier configuration**.

Then normal costates can diverge while:

- physical trajectories remain finite;
- thrust approaches a finite value;
- time approaches a finite value or changes slowly;
- physical motion per unit costate-based arclength vanishes;
- finite computed augmented Jacobians remain full rank.

That can look exactly like an endlessly approached “fold,” but there is no finite corner to round in the current coordinates.

**Your table does not prove this alternative either.** But the reported history of costate growth means you should actively test it.

Log
\[
\|\lambda_0\|,\quad
\max_t\|\lambda(t)\|,\quad
\min_t\|\lambda_v(t)\|,
\]
along with state/time versus costate contributions to each accepted arclength step.

## A useful model-specific sensitivity check

For a smooth family of normal stationary solutions with fixed endpoints,
\[
\frac{dt_f}{ds}
=
\left(\int_0^{t_f} H_{T_N}\,dt\right)
\frac{dT_N}{ds},
\]
using consistent nondimensional units and thrust conversion.

For your all-burn minimum-direction law,
\[
H=1+\lambda_r^Tv+\lambda_v^Tg
-\frac{T_{\rm nd}}m\|\lambda_v\|
-\frac{T_{\rm nd}}c\lambda_m,
\]
so
\[
H_{T_{\rm nd}}
=-\frac{\|\lambda_v\|}{m}-\frac{\lambda_m}{c}.
\]
With free final mass,
\[
\lambda_m(t_f)=0,\qquad
\dot\lambda_m=-\frac{T_{\rm nd}\|\lambda_v\|}{m^2},
\]
hence \(\lambda_m(t)\ge0\) along this regular all-burn model.

Consequently, along a regular normal branch, \(dt_f/dT<0\), with a finite slope when the relevant multipliers remain bounded.

**At a finite normal thrust fold, \(dt_f/ds\) must therefore vanish too.** A picture in which thrust turns while time passes through with nonzero arclength derivative is not consistent with bounded normal multipliers in this formulation.

Check this identity numerically. It is a strong discriminator among:

- a genuine finite normal fold;
- multiplier blow-up/normality loss;
- inaccurate tangents or derivatives.

## If multipliers are diverging

Use a homogeneous PMP representation:

- introduce objective multiplier \(\rho\);
- use \(H=\rho+\lambda^Tf\);
- impose a normalization such as
  \[
  \rho^2+\|\lambda_0\|^2=1;
  \]
- retain the nonnegative-objective-multiplier condition for admissible minimum-time PMP candidates.

Initialize by setting
\[
a=(1+\|\lambda_0^{\rm old}\|^2)^{-1/2},\qquad
\rho=a,\qquad
\lambda^{\rm new}(t)=a\,\lambda^{\rm old}(t)
\]
throughout the trajectory.

For a positively homogeneous all-burn costate system, this preserves the state trajectory and steering. It replaces an infinite-multiplier chart with bounded initial multipliers and lets you inspect whether \(\rho\to0\).

Continuing algebraically into \(\rho<0\) does not make those points valid minimum-time PMP candidates.

# Fold detection and localization

Your sign-change test is a reasonable **crossing detector**, once tangents are reliable. It is not an adequate early-warning or localization procedure.

Also, **18.120 d is the last reported time, not an established fold time.**

For a simple fold at a finite point, verify:
\[
R=0,\qquad
\operatorname{rank}R_x=n-1,
\]
and, with right/left null vectors \(\phi,\psi\),
\[
\psi^T R_\theta\ne0,\qquad
\psi^T R_{xx}[\phi,\phi]\ne0.
\]

The first nonzero condition makes the full solution set a regular curve; the second makes the thrust turning point quadratic rather than degenerate.

### In flight

Monitor:

- signed \(\tau_\theta\);
- the two smallest singular values of \(R_x\);
- relative conditioning and numerical accuracy of \([R_x\ R_\theta]\);
- estimated thrust curvature \(d\tau_\theta/ds\);
- boundedness of the unknowns.

A local prediction
\[
s_{\rm fold}-s_k\approx
-\frac{\tau_\theta}{d\tau_\theta/ds}
\]
is useful when curvature is resolved and nonzero. Small \(\tau_\theta\) alone is not a fold certificate.

### Localize a candidate

Solve the extended system
\[
R(x,\theta)=0,\qquad
R_x(x,\theta)\phi=0,\qquad
c^T\phi=1,
\]
where \(c\) is a fixed normalization vector near the candidate right null vector.

That finds an isolated fold in a one-parameter family. It requires reliable Jacobian derivatives/Hessian-vector products. Use it after establishing that the candidate is finite and numerically resolved—not to force an asymptote into a fold interpretation.

# B. How I would proceed

## 1. Turning and mapping: implement robust standard pseudo-arclength

**My choice is your existing bordered pseudo-arclength formulation, with full-SVD tangents, accurate thrust derivatives, and a globalized corrector.**

Why:

- It is exactly the tool needed for an ordinary fold.
- You have nearly all of it.
- It avoids selecting another coordinate that may itself turn.
- At this problem size, robust dense linear algebra is affordable.

I would not introduce deflation to cross a simple fold. Deflation is potentially useful later to discover additional roots.

I would not make switching to \(t_f\) the main method. It works only where \(\tau_{t_f}\ne0\), and the sensitivity identity above shows why it is not guaranteed to regularize a finite normal thrust fold.

**If the corrected method still makes substantial arclength progress only through growing costates, stop shrinking steps and “trying harder to turn.” Reformulate the multiplier normalization.**

## 2. Is 70 mN reachable, and are there intermediate branches?

**Numerically, yes: your 70 mN trajectory is already the evidence.** Independent propagation and indirect acceptance are much stronger evidence than continuation failure is evidence of infeasibility.

Whether it is the shortest solution is unresolved.

A verified local minimum of thrust on one extremal branch means only that the branch initially returns toward **higher** thrust. Turning that fold would not immediately deliver a lower-thrust solution. There could be:

- further folds reconnecting to lower thrust;
- another connected component;
- additional local minima at intermediate times;
- nonminimizing extremals connecting pieces of the diagram.

There is no basis to assert that intermediate branches must exist—or that they do not.

The winding observations are useful morphology, **not branch identification**. Distinct stationary solutions can share winding, and fractional winding of open trajectories is not a complete topological classification.

### Structured overnight search

I would run three parallel tracks.

**Track 1 — Continue both validated families.**

- Restart the fast family before the suspected limit.
- Continue the 70 mN family toward higher thrust as well as lower thrust.
- Map \((T,t_f)\), multiplier norms/normality, and local-optimality diagnostics.
- Retain nonminimizing extremals in the branch map; label them rather than discarding them.

Starting from both ends is much more informative than another 500 one-sided steps.

**Track 2 — Direct multistart at 70 mN.**

Use full-trajectory seeds from both families, not just initial costates. Include:

- several initial flight times below 26.436 d;
- smooth perturbations of thrust-direction histories;
- retimed trajectories;
- different close-approach/waypoint patterns;
- other winding patterns, not just the observed one.

Polish every distinct promising solution with multiple shooting, independent propagation, and your second-order tests.

**Track 3 — Time-cap feasibility searches.**

Solve direct problems with
\[
t_f\le t_{\rm cap}
\]
for several caps below the incumbent, using multistart or continuation of the cap. A feasibility-restoration objective can help generate seeds.

Do not infer infeasibility from local solver failure.

Also, do not assume feasibility is monotone in **exact fixed flight time**: with forced full burn, “arrive earlier and wait” is not automatically admissible. Feasibility under a time **upper bound** is monotone by set inclusion; feasibility at a prescribed exact time need not be.

For additional organization, arrival-phase homotopy and deflated stationary-root searches can help find other branches. Any solution used in the final comparison must return to the originally specified phase.

## 3. Is the fold physical or geometric? Poster wording

These are not cleanly separate categories. A fold of an extremal branch is a geometric property of the specified dynamical optimization problem. It may depend strongly on endpoint phase.

**Changing arrival phase changes the problem and can move, remove, or introduce branch folds.** With arrival phase as a second parameter, you can continue a locus of folds rather than one isolated turning point.

But the proposed poster claim is not defensible:

> “The minimum thrust for this transfer geometry is 72 mN.”

Your same-endpoint 70 mN solution contradicts it.

Before localization, say something like:

> “The continued short-time extremal branch approaches 71.99 mN, where the present formulation becomes strongly ill-conditioned; a distinct 26.436 d solution is obtained at 70 mN.”

After verifying a finite fold:

> “A local thrust turning point is identified on the continued short-time extremal branch for the specified departure and arrival phases. This is not a global feasibility threshold.”

And avoid calling an entire continued stationary branch “minimum-time” without checking which portions remain locally minimizing.

## 4. Can you certify 26.436 d as global minimum overnight?

**Not from continuation, multistart, and conjugate tests alone.**

Taking your tests at face value, the supplied evidence supports accurate feasibility and local optimality under the tested model. It does not exclude a disconnected faster basin.

A genuine global certificate needs:

1. a validated feasible upper bound; and
2. a global lower bound excluding all shorter admissible transfers, to a stated tolerance.

Possible routes include:

- verified reachability exclusion;
- a verified HJB subsolution/value-function lower bound;
- rigorous global optimization with validated dynamics and a justified representation of all admissible controls;
- sufficiently tight certified relaxations.

For nonlinear CR3BP dynamics with mass and free thrust direction, those are substantial research/computational tasks. Certifying only a finite collocation NLP does not automatically certify the continuous-time control problem.

**Realistic overnight deliverable:** a well-documented best-known solution, supported by structured multistart, two-sided continuation, mesh/tolerance refinement, independent propagation, and local sufficiency tests.

Report the search envelope and number of distinct branches/candidates. Say “lowest time found” or “best-known locally minimizing solution,” not “global minimum.”

# What I think is wrong in the present reasoning

The important corrections are:

- **Small \(R_x\), apparently regular augmented Jacobian, and small thrust tangent do not yet establish a simple fold.**
- **A simple fold is not an endpoint of the solution curve.**
- **The last computed \(t_f\) is not the fold time.**
- **The existing tangent formula is fragile, but does not itself force one-sided continuation.**
- **A passing conjugate test does not establish globality, and singularity should not automatically be dismissed as unrelated to local-optimality degeneracy.**
- **Same winding does not mean same extremal branch.**
- **A branch-local thrust minimum is not a feasibility threshold. Here the 70 mN solution explicitly disproves a 72 mN threshold for the same full transfer problem.**

**First action:** make the full-SVD tangent change, log actual accepted steps and multiplier norms, and check the time–thrust sensitivity identity. Those inexpensive changes should tell you whether you need ordinary continuation improvements—or whether you have been trying to turn a corner that lies at infinity in the normal-costate chart.