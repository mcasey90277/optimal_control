Here are the explicit answers to your highest-priority questions, followed by a prioritized list of code findings.

### Q1. HERMITE-SIMPSON CORRECTNESS AND ORDER
**(a) Correctness of equations:** Yes. The equations, signs, and the factors of `1/8` and `1/6` exactly match the standard Separated Hermite-Simpson formulation. 
**(b) Genuinely fourth-order:** Yes. This formulation perfectly preserves the $O(h^4)$ global convergence properties of the Hermite-Simpson method.
**(c) The empirical puzzle ($O(h^5)$ measured):** You are experiencing **(iii)**. Your measurement script `dro_residual.m` re-integrates *each individual interval* from $t_k$ to $t_{k+1}$ and compares the result to $X_{k+1}$. This isolates and measures the **Local Truncation Error (LTE)** of the interval. For a scheme with a global error of $O(h^p)$, the local step error is mathematically $O(h^{p+1})$. Since Hermite-Simpson is globally 4th order, its LTE is strictly $O(h^5)$. A 2x refinement yielding a ~32x ($2^5$) reduction confirms your implementation is executing flawlessly. Your documentation is mathematically standard (schemes are classified by their global order), but you should clarify that the *gate* measures the $h^5$ local error.
**(d) Free $U_m$ vs averaged:** Your reasoning is entirely correct. Making $U_m$ a free variable is strictly necessary. Averaging $U_m = (U_k + U_{k+1})/2$ degrades the integration rule, dropping the scheme's global accuracy to 2nd order. Furthermore, the arithmetic average of two unit vectors lies strictly inside the unit sphere, which would either force a dynamic re-normalization (altering the order properties) or violate the physical bounds of the control problem. 

### Q2. CONTROL RECONSTRUCTION IN THE ERROR MEASUREMENT
**(a) Implied control representation:** Yes. The Hermite-Simpson state quadrature inherently assumes that the system vector field (and thus the continuous control) can be represented by a quadratic interpolating polynomial over the subinterval.
**(b) Basis correctness:** Yes. Your Lagrange polynomial `u(w) = (2w-1)(w-1)*ua + 4w(1-w)*um + w(2w-1)*ub` is algebraically perfect. It evaluates exactly to $u_a$, $u_m$, and $u_b$ at $w = 0, 1/2, 1$.
**(c) Trapezoidal control:** Yes. Trapezoidal collocation strictly implies a piecewise linear control between nodes.
**(d) Re-normalization inconsistency:** **Yes, an inconsistency is introduced.** The NLP enforced `norm(u) == 1` strictly at the nodes (and midpoint). However, a polynomial traversing through points on a sphere will chord *inside* the sphere between those points. By re-normalizing dynamically inside `local_rhs`, you are integrating a slightly different vector field than the polynomial approximation the NLP implicitly solved. However, this inconsistency is an unavoidable artifact of direct collocation (polynomials cannot naturally parameterize spherical manifolds). Doing it this way is standard practice; it just contributes an $O(h^3)$ perturbation to the defect.
**(e) Metric choice (Re-integration vs. Betts residual):** Re-integration measures the physical local drift in physical units (km, km/s). Betts evaluates the continuous algebraic defect of the interpolating polynomial inserted into the vector field ($f(x(t), u(t)) - \dot{x}(t)$). 
*Tradeoffs:* The Betts residual is vastly faster to compute (no ODE solver needed) and directly mirrors the mesh-refinement math. However, its units are abstract (e.g., $km/s^2$ acceleration defect). Your re-integration approach is computationally heavy but gives absolute physical intuition. For a final pass/fail gate, re-integration is often preferred by practitioners for its unquestionable physical meaning.

### Q3. GATE DESIGN
**(a) G1 max-of-local vs. global error:** Your suspicion is correct: **max-of-local drastically understates the true global error.** Local errors accumulate, and in the CR3BP regime, they can amplify exponentially via Lyapunov divergence. A local error of 1 km per step over 1000 steps could easily manifest as thousands of kilometers of actual terminal miss distance. However, *true global error* (integrating from $t_0$ to $t_f$) is often fundamentally uncomputable for highly unstable trajectories because the open-loop integration will blow up regardless of discretization. Max-local is a highly standard, pragmatic gate for chaotic regimes, but you must acknowledge it is only a local confidence check. 
**(b) G2 structural flaw:** If a path constraint is provided, the unconstrained indirect solution is invalid. You should simply automate the bypass. Inside `certify_dro_mintime.m`, before processing gates, add: `if isfield(o,'minAltKm') && ~isnan(o.minAltKm), tfRef = []; end`. This properly self-disables G2 when you run the constrained problem.
**(c) Missing/vacuous gates:** 
*   **Vacuous:** G5 (Terminal Boundary Error) checks `norm(X(:,end) - rvf)`. Since `X(:,end) == rvf` is an explicit hard equality constraint in the NLP, G5 is basically just re-checking that IPOPT met its feasibility tolerance (which G3 already monitors). It is harmless but redundant.
*   **Missing:** You verify `thrMin >= 1 - 1e-6` at the *nodes* (G6), but since Hermite-Simpson implies a quadratic control interpolant, it is mathematically possible for a quadratic passing through $1, 1, 1$ to bow downwards between nodes if the derivatives dictate it (though unlikely for bang-bang). It is safer to pull the continuous min-throttle directly out of the ODE sampling pass.

---

### FINDINGS (Code Review)

- **[CORRECTNESS]** `certify_dro_mintime.m:150` — `ode113` restricted to exactly 9 output points. By passing `linspace(a,b,9)` as the `tspan`, you force `ode113` to return ONLY those 9 specific points, silently dropping the dense integration steps it calculates internally. This defeats the purpose of looking for a continuous path-constraint violation between nodes.
  *Fix recommendation:* Call `ode113` with `[a b]` so it returns its natural dense grid, or use `deval(sol, linspace(a,b,100))` to sample the interpolant much more finely.

- **[CORRECTNESS]** `casadi_mintime_dro.m:306` — `out.Xm` is not preserved. Because you evaluate `out.Xm = []` and never capture it from `sol.value(Xm)`, the defect check at line 313 computes `Xmn` using the exact interpolation algebraic rule. This masks any solver tolerance violations specific to `Xm`.
  *Fix recommendation:* Save it via `out.Xm = full(sol.value(Xm));`. Then use `out.Xm` directly when re-computing `Fmn` at line 315 rather than rebuilding `Xmn` algebraically.

- **[EFFICIENCY]** `casadi_mintime_dro.m:242` — Objective function puts the entire gradient on a single lifted node variable `TF(1)`. This forces the solver to propagate the cost gradient through $N$ equality constraints, degrading the conditioning of the KKT matrix.
  *Fix recommendation:* Use `opti.minimize(mean(TF));` or `opti.minimize(sum(TF)/nN);`. It yields the exact same mathematical problem but gives the NLP structurally symmetric dual paths.

- **[EFFICIENCY]** `casadi_mintime_dro.m:269` — Warm-start function evaluation is placed in a native MATLAB `for` loop. CasADi functions invoked iteratively in MATLAB code are notoriously slow.
  *Fix recommendation:* Vectorize the initial state mapping using CasADi's `.map()` feature. Replace the loop with: `F0 = full(fdyn.map(nN)(X0, U0));`.

- **[ROBUSTNESS]** `casadi_mintime_dro.m:208` — Strict inequality used for constraints (`> 0.05`). NLP solvers (including IPOPT) fundamentally operate on closed sets and treat `>` identically to `>=`. Writing strict inequalities is misleading to readers who might assume the solver prevents the variable from resting exactly on the boundary.
  *Fix recommendation:* Change to `>= 0.05`. (Also applies to `TF(:) >= 0` at line 207). 

- **[STYLE]** `certify_dro_mintime.m:73` — `altKm(min(R.kWorst+1, numel(altKm)))`. The worst error interval occurs between nodes, but this logs the altitude *at the right-hand node* of the worst interval. The variable is confusingly named `worstAltKm` implying it was the lowest altitude seen, when it is just a point-reference for where the worst local error occurred.
  *Fix recommendation:* Rename to `altAtWorstErrorKm` or similar to clarify intent.
