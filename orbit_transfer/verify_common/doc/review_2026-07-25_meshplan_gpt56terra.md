## Q1 — primal order

| Quantity | Arbitrary uniform mesh, unaligned simple switches | Switches as true mesh breakpoints / phase boundaries |
|---|---:|---:|
| State, global norm | **generically \(O(h)\)** | **\(O(h^2)\)** for trapezoidal collocation |
| Free final mass / \(\Delta V\) | **generically \(O(h)\)** | **\(O(h^2)\)** |
| Terminal state error | **\(O(h)\)** if measured by continuous propagation | **\(O(h^2)\)** |

A switch-crossing interval has an \(O(h)\) integral error: its wrong thrust duty-cycle occupies \(O(h)\) time with an \(O(1)\) RHS jump. One such interval dominates the \(O(h^2)\) smooth-arc errors. Thus rows 32–33 should **not** predict \(p\approx2\) for a uniform nodewise bang-bang transcription.

If terminal equalities are imposed in the NLP, “terminal state error” is identically solver tolerance and is **not a convergence observable**. Replace it with terminal error after independently integrating the reconstructed control under the continuous dynamics, or a dense-trajectory error against a reference.

With an actual breakpoint at every switch, each arc is smooth and trapezoidal’s second order is available again. Merely placing *many* nodes near an estimated switch is not the same guarantee: the switch must remain aligned to a mesh breakpoint to \(o(h)\), ideally be a decision variable / phase boundary.

Relevant theory: Dontchev, Hager & Veliov, *Second-Order Runge-Kutta Approximations in Control Constrained Optimal Control*, SIAM J. Numer. Anal. 38 (2000), 202–226, doi:10.1137/S0036142999351765. Its second-order result requires its regularity/coercivity hypotheses; it is not a blanket result for an unaligned nodewise bang-bang solution.

## Q2 — costates and defect multipliers

For a smooth regular problem, a **properly scaled and covector-mapped** Runge–Kutta multiplier/costate has the order of the RK transcription: \(O(h^2)\) at mesh points for trapezoidal. Raw defect multipliers are not automatically nodal samples of the continuous costate: their sign, scaling, and endpoint interpretation depend on the defect convention.

For unaligned bang-bang switches, expect the same generic event-location limitation: **\(O(h)\) globally**, absent additional structure. Do not claim a universal multiplier order without first verifying that your multiplier reconstruction is the discrete adjoint consistent with your precise CasADi defects and quadrature.

There is no universal “trapezoidal boundary order is one” theorem. What often looks like boundary order reduction is instead an **unmapped raw final defect multiplier**. The transformed/covector-mapped terminal adjoint satisfies the discrete terminal transversality condition exactly, to NLP KKT tolerance.

Reference: Hager, *Runge-Kutta Methods in Optimal Control and the Transformed Adjoint System*, Numerische Mathematik 87 (2000), 247–282, doi:10.1007/s002110000178.

## Q3 — \(\lambda_m(t_f)=0\)

Your conclusion is correct **only for the raw last defect multiplier**. With
\[
D_N=x_{N+1}-x_N-\frac{h_N}{2}(f_N+f_{N+1}),
\]
stationarity in \(x_{N+1}\) gives, up to your sign convention,
\[
0=\phi_{x_f}+\frac{h_N}{2}L_{x_{N+1}}
+\left(I-\frac{h_N}{2}f_{x_{N+1}}\right)^T\Lambda_N.
\]
For free final mass and \(\phi_m=0\), the raw \(\Lambda_{N,m}\) is indeed generally \(O(h)\). That is a **known endpoint representation offset**, not evidence that continuous \(\lambda_m(t_f)\neq0\).

Define instead the mapped terminal covector
\[
\widehat\lambda_f =
\frac{h_N}{2}L_{x_{N+1}}+
\left(I-\frac{h_N}{2}f_{x_{N+1}}\right)^T\Lambda_N .
\]
Then \(\widehat\lambda_{m,f}=0\) by the discrete KKT system, apart from IPOPT tolerance. Compare a continuous-costate reconstruction with this mapped quantity, not with raw `Lambda_{N,m}`.

Therefore: if `lamMassEndRel` is raw \(\Lambda_{N,m}\), an observed \(p\approx1\) is expected and should **not** be called a real first-order anomaly. It is not an estimator of transversality error. If it is already the mapped covector, \(O(h)\) is not the default prediction; investigate the mapping, endpoint constraint multipliers, and normalization.

## Q4 — Richardson

Richardson is defensible only after establishing a fixed smooth asymptotic regime:

- **Defensible:** smooth-arc quantities on switch-aligned meshes; state/functionals after phase-wise reparameterization; mapped costates away from switches; a quantity empirically showing stable topology and a consistent signed expansion.
- **Not defensible as a default:** switch count, unmatched switches, raw switch times on a grid, raw endpoint defect multipliers, sign percentages, KKT residuals, and values across a branch/topology change.
- **Uniform unaligned bang-bang ladder:** it may empirically fit an \(O(h)\) expansion, but moving switch locations make the leading coefficient oscillatory with mesh phase. A three-level Richardson value can be badly misleading.

Use a continuous-dynamics propagation defect estimator as the primary estimator: reconstruct control, locate its switches, integrate each arc accurately, and compare propagated states against the collocation polynomial/nodes. Report local errors, especially in switch cells. Adjoint-weighted residual estimates are also appropriate for target functionals, provided jumps are included explicitly.

## Q5 — warm-chain bias

Yes. Warm chaining can lock a discrete branch and its inherited switch topology. Passing `Solve_Succeeded`, defects, and NLP FOCs only establishes that each level is a stationary point of its own NLP; it does **not** establish that all levels represent the same continuous branch or the best discrete approximation.

Add at least one independent branch check for every headline row:

1. Fine solve from the warm chain **and** from an independently constructed fine-grid seed.
2. Perturb final time / switch locations / throttle seed in both directions.
3. Run fine-to-coarse continuation as a reverse check.
4. Compare objective, switch topology, and dense propagated trajectory.
5. Where ambiguity persists, solve competing fixed-topology / multiphase formulations.

A differing switch count with comparable objective is a branch finding, not a failed solver gate.

## Q6 — switch structure

Your treatment of count is right: require stabilization, but also report unmatched switches and their physical significance.

For an unaligned nodewise bang-bang transcription, the apparent switch location is grid-quantized, so its generic accuracy ceiling is **\(O(h)\)** for a transversal/simple switch. No \(O(h^2)\) switch-time claim is justified merely because the state scheme is trapezoidal.

For a simple switch treated as an optimized phase boundary/event time, with smooth arcs and a nonzero derivative of the switching function, **\(O(h^2)\)** switch-time convergence is the natural trapezoidal expectation. Degenerate/grazing switches, near-singular arcs, clustered switches, and topology changes invalidate that prediction.

Agamawi, Hager & Rao, *Mesh Refinement Method for Solving Bang-Bang Optimal Control Problems Using Direct Collocation* (2019), arXiv:1905.11895, is directly relevant: it detects switches from the switching function and uses them to create mesh structure. It supports treating switch localization as a dedicated event-detection problem, not as an ordinary scalar Richardson observable.

## Design corrections

- Replace the table’s blanket \(p\approx2\) for `mf`, `dV_kms`, and terminal error with **\(p\approx1\) unaligned; \(p\approx2\) only with verified switch-aligned phase boundaries**.
- `kktStatInf`, `sdotMinRel`, `signPct`, and raw dual diagnostics are solver/discretization diagnostics, not convergent physical quantities; do not put them through the same `CONVERGED` verdict.
- Use physical local step sizes, especially near each switch, not node multipliers alone. Longitude and Sundman meshes make “\(h=1/N\)” physically nonuniform.
- Three levels yield one fragile slope. Use at least four, preferably \(1,2,4,8\), and report sliding three-level slopes plus signed differences.
- For PSR, state the falsifiable claim correctly: it may improve constants and recover second order **only if it makes switches true retained breakpoints**, not simply dense neighborhoods.

