Here is your domain-aware review.

### Q1 — The Dual $\mapsto$ Costate Map

Your diagnosis of the KKT architecture is correct: mathematically, the primer vector should align perfectly, and the anomaly is purely an extraction artifact introduced by IPOPT's row-scaling of the highly coupled `cScale` constant-state equation.

**The Algebra:**
Let $\mathcal{L}$ be the NLP Lagrangian. Ignoring path bounds on the states, we isolate $\alpha_k$, which appears *only* in the velocity block of the defect equations $C_k^{(v)}$ and the norm constraint:
$$C_k^{(v)} = v_{k+1} - v_k - \frac{1}{2}d\tau \left[ cS_k \kappa_k f_v(x_k,u_k) + cS_{k+1} \kappa_{k+1} f_v(x_{k+1},u_{k+1}) \right] = 0$$
Where $f_v(x_k,u_k) = g(r_k) + \frac{T_{max}}{m_k}s_k \alpha_k$.

The stationarity condition for $\alpha_k$ evaluated on an internal node ($1 < k < N+1$) is:
$$\frac{\partial \mathcal{L}}{\partial \alpha_k} = \underbrace{-\lambda_{v,k-1}^T \frac{d\tau}{2} \left[ cS_k \kappa_k \frac{T_{max}}{m_k} s_k \right]}_{\text{from } C_{k-1}^{(v)}} \underbrace{-\lambda_{v,k}^T \frac{d\tau}{2} \left[ cS_k \kappa_k \frac{T_{max}}{m_k} s_k \right]}_{\text{from } C_k^{(v)}} + 2\mu_{\alpha,k} \alpha_k = 0$$

Because $cS_k$, $\kappa_k$, $T_{max}$, $m_k$, and $s_k$ are all purely given as evaluated *scalars* at node $k$, they completely factor out of the vector sum. Rearranging:
$$\alpha_k \propto \left( \lambda_{v,k-1} + \lambda_{v,k} \right)$$
*Notice that neither $\kappa$, $cScale$, nor $d\tau$ (since it is uniform) shift the direction*. The mathematically required NLP costate for the primer *is* the strict node-centered average of the adjacent interval velocity duals.

**Why it reads 17–24 degrees off:**
When `cScale` is packaged as $X(9,:)$ with equality constraints $cS_{k+1} - cS_k = 0$, IPOPT scales rows via `nlp_scaling_method`. Since $cS_k$ multiplies *every single dynamics row* in the Jacobian, IPOPT heavily rescales the $cScale$ constraints against the physical constraints to keep the KKT matrix well-conditioned. CasADi's `opti.dual()` un-scales the multipliers, but because of numerical mixing at 1e-9 tolerance between the time-domain costate $\lambda_t$ and the "slack clock" costate, the retrieved velocity block `lamDef(4:6,:)` evaluates to a rotated projection of the true costate space. 

**The Decisive Test:**
Your candidate fix is perfectly sound. By mathematically collapsing to the 8-state model (absorbing $c^*$ into $d\tau$), you eliminate row 9 and its scaling cross-talk. 
- **Decisive Test:** Run the `cScale`-eliminated 8-state NLP warm-started exactly at the previous solution, force 0-1 iterations, and extract `lamV`. If the primer alignment immediately plummets to ~0.06 deg (matching the CR3BP sibling), the anomaly is confirmed as an IPOPT scaling artifact of the slack state. 

---

### Q2 — The Low-Thrust Min-Time Stall

Your diagnosis is exactly right. Interior-point methods use `warm_start_init_point = 'yes'` and steep bound/multiplier pushes to safely initialize the barrier $\mu$ close to zero. If the iterate has a defect of 5e-3, it is severely primal-infeasible. Applying tight NLP parameters here forces the line-search to fail step-acceptance criteria, immediately throwing IPOPT into the restoration phase which gets trapped in a local minimum of constraint violation.

**Ranked Candidate Causes for the 5 N failure:**

1. **`warmTight=true` locking out the Infeasible Basin**
   *Why it happens at 5 N but not 10 N:* At 10 N, 4.5 revolutions means the initial Cartesian phase error (winding difference) from the cold seed is easily pliable. At 5 N, 9 revs means the coordinate wrap spirals map tightly over one another; a 5e-3 defect here crosses highly disjoint topology, so the tight barrier immediately traps.
   *Decisive Test:* Run your proposed fix (`warmTight=false` until defect < 1e-6). 
2. **Spurious Free-L Manifold Traps (Phase Ambiguity)**
   *Why:* The terminal manifold (5 constraints) does not restrict angular momentum direction ($r \times v$) nor integer phase wraps. At high revs, the manifold gradients can attract the solver toward a localized (e.g. out-of-plane or inverted) topology.
   *Decisive Test:* Set `term=fixed`, pinning the exact final $(r,v)$ of the Route-B energy seed, and try the min-time solve. If `fixed` converges cleanly but `manifold` stalls, the manifold branch traps the optimizer.
3. **Conditioning Growth in Cartesian Coordinates**
   *Why:* Trapezoidal collocation condition numbers scale roughly as $O(N^2 \cdot Revs^2)$ in Cartesian coordinates.
   *Decisive Test:* Examine the linear solver output in the iteration log; if `alpha_primal` and `alpha_dual` are frequently dropping below 0.1, the matrix is locally singular.

---

### Q3 — Baseline Integrity

The Cartesian baseline is highly solid. The constraints, dynamics, and structural implementation accurately reproduce the mathematical target.
* **Quadrature Int[s]dt:** By the chain rule, $dt = cS \cdot \kappa \cdot \text{tauf0} \cdot d\tau_{unit}$. The discrete evaluation $\sum \frac{d\tau}{2} (w_k + w_{k+1})$ where $w_k = cS_k \kappa_k \Phi(s_k)$ perfectly aligns with exact trapezoidal integration over $\tau$.
* **Units:** Your canonical nondimensionalization ($AU \approx 0.224 \text{ m/s}^2$ and $c \approx 6.38$) corresponds precisely to the specific equations and problem boundaries given.
* **[ROBUSTNESS]** `earth_elliptic_to_geo/lt2b_rhs_time.m:21` — `rn2 = r(1)^2... + 1e-12`. While AD-safe, ensure that upstream functions consuming terminal state evaluation don't unintentionally fail tight terminal gap tolerances ($\text{termErr } < \text{1e-8}$) due to gradient jitter on the fractional offset.

---

### Q4 — Vet the Go-Forward Plan (MEE + True Longitude)

The plan to adopt Modified Equinoctial Elements (MEE) with True Longitude ($L$) as the independent variable is the mathematically pristine standard for long-duration low-thrust, **but requires three specific structural revisions:**

**(a) $L$ as the Independent Variable**
$L$ is unequivocally monotone for this regime. $dL/dt = h/r^2 > 0$ strictly, as typical thruster accelerations ($10^{-4}$ g) cannot stall orbit momentum.

**(b) The Clock / Sundman Role**
You do *not* need a Sundman transform. $dt = \frac{1}{\dot{L}} dL$ acts as a built-in Sundman
Here are the three specific structural revisions required to implement the MEE + $L$ formulation, followed by the completion of Q4.

**The Three Structural Revisions for the MEE+$L$ Transcription:**
1.  **Domain Mapping:** Map the free terminal longitude $L_f$ to a fixed unit grid $\sigma \in [0,1]$ using $L(\sigma) = L_0 + \sigma \Delta L$. The total angle $\Delta L$ becomes a scalar decision variable.
2.  **State Augmentation:** Physical time $t$ becomes a strict state equation: $\frac{dt}{d\sigma} = \frac{\Delta L}{\dot{L}(x,u)}$.
3.  **Objective Measure:** The fuel/energy objective integrand must be rescaled by the spatial measure: $J(\epsilon) = \int_0^1 \Phi_{homotopy}(s) \frac{\Delta L}{\dot{L}(x,u)} d\sigma$.

---

### (b) The Clock / Sundman Role
$L$ entirely subsumes the Sundman transform. You do not need a separate $\tau$ or `cScale`. By dividing the Gauss Variational Equations by $\dot{L} = \frac{h}{r^2} + \text{thrust}_W$, the system naturally slows its spatial steps at apogee (where $r$ is large, so $dL/dt$ is small), automatically concentrating mesh points exactly where the apogee-centered burns occur. 
**If $\dot{L} \to 0$:** This would mean angular momentum is stalling to zero or reversing. For low thrust ($10^{-4}$ g) in Earth orbit (1 g), the thrust normal component $w$ is orders of magnitude too weak to overcome $\frac{h}{r^2}$. $\dot{L}$ remains strictly positive.

### (c) Singularities and Conditioning Traps
*   **Targeting $e \to 0, i \to 0$:** MEEs are specifically designed to be perfectly smooth and non-singular here (the coordinates evaluate to exactly $0$).
*   **True Singularity:** MEEs are singular at retrograde equatorial orbits ($i = 180^\circ$, where $hx, hy \to \infty$). You are moving $7^\circ \to 0^\circ$, so this is safely bounded.
*   **[ROBUSTNESS] Conditioning Trap:** $\Delta L$ will be roughly $754 \times 2\pi \approx 4700$ radians. If you use raw $L$ inside trig functions (e.g., $\cos(L)$ for the $r$ expansion), finite differencing or large state-scale disparity will damage the KKT matrix. Wrap the evaluation to $mod(L, 2\pi)$ inside your dynamics, or ensure CasADi is explicitly outputting analytical derivatives for the trig evaluations so large operands do not cause catastrophic cancellation.

### (d) Node-per-rev FLOOR and Collocation Methods
*   **The Floor:** To capture a minimum-fuel sequence of coast-burn-coast-burn, you need at bare minimum ~3 nodes per segment to avoid immediate infeasibility bounds. This dictates an absolute floor of **$12$ to $16$ nodes per revolution**. At $754$ revs, your mesh sits at $9,000 - 12,000$ points.
*   **Does Hermite-Simpson (HS) or Pseudospectral (PS) help?** **NO.** 
    *   PS assumes infinite smoothness; applying it to bang-bang control guarantees severe Gibbs ringing around the switch points unless knot sequences are placed exactly on the implicit switches (intractable for 1,500+ switches).
    *   HS improves state integration but *still smears* discontinuous jumps linearly across the interval. 
*   **The "Trapezoidal Smears" Lesson:** You cannot fix smearing on a rigid grid. The required countermeasure is to **not solve it as pure bang-bang**. You must rely heavily on the $\epsilon$-homotopy. Stop the continuation at a small positive $\epsilon$ (e.g., $10^{-4}$) rather than pure $0$. This slightly rounds the corners of the burn, absorbing the smearing into physical fuel efficiency rather than KKT constraint violations.

### (e) What else will bite us / Defeating the Paper's Assertion
*   **The Paper is Outdated:** The 2004 assertion that direct methods are "predictably unsuitable" was true for dense-Jacobian sequential quadratic programming (SNOPT) and finite differencing. A 100,000-variable NLP with structured sparsity, fed with CasADi's exact algorithmic differentiation (AD) exact-Hessians into MUMPS/IPOPT, solves in minutes on modern hardware. Your MEE + $L$ plan fundamentally defeats their argument.
*   **[SEVERITY] ROBUSTNESS -- Trap:** The multiple-revolution basin. Even in MEE, an initial guess with 700 revolutions provided to IPOPT will experience extreme local minima trapping against $701$ or $699$ revolutions. A cold continuous solve at $0.1$ N will fail. 
    *   **Fix:** You must execute a geometric continuation strategy on `Tmax`. Solve at $10$ N $\to$ warm start $5$ N $\to$ warm start $2.5$ N $\dots \to 0.1$ N. Let the optimizer naturally stretch $\Delta L$ step-by-step to reach the $754$ rev basin safely.
