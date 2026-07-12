## Diagnosis and Recommendations

---

### Q1 — The fold + gauge diagnosis is partially wrong, and the error matters

**What you're diagnosing correctly:** There is a genuine near-null direction in the Jacobian at the min-time anchor, and the pseudo-arclength corrector is drifting along it.

**Where the diagnosis goes wrong:** The near-null direction is *not* a scale gauge of the min-fuel problem. It's the manifold-separation direction between two different BVPs.

At t_f = t_f,min the min-time and min-fuel problems share the same optimal *trajectory* (all-burn is the only feasible control), but their transversality conditions are structurally different:

| Problem | Terminal transversality |
|---|---|
| Min-time (free t_f) | H(τ_f) = 0; all costates free at boundary |
| Min-fuel (fixed t_f) | λ_m(τ_f) = 0; Hamiltonian *not* pinned to zero |

When you evaluate R_min_fuel at the min-time costates and get ‖R‖ ~ 2e-7, you're not measuring distance along the min-fuel branch to a fold — you're measuring the transverse separation between two closely tangent but distinct costate manifolds. The pseudo-arclength corrector, starting *off* the min-fuel manifold, must simultaneously (a) correct onto it and (b) advance along it. With cond ~2e9 and a loose stopping tolerance forced by the ~1e-6 anchor accuracy, those motions are not separated, and λ0 drift occurs at fixed t_f.

**The "fold":** A genuine switch-birth bifurcation IS expected near t_f,min on the min-fuel branch (the branch is near-vertical in the (t_f,λ0) plane), but it lives on the min-fuel manifold. You can't see it cleanly because you never fully landed on that manifold.

**Pinning ‖λ0‖ = 1 will not fix this.** That addresses the scale gauge of the *min-time* problem, which is a different null direction. The min-fuel problem has no scale gauge; the cost breaks that symmetry. The corrector will still drift along the manifold-separation direction. This is the concrete thing your Q1 diagnosis gets wrong, and it explains why standard gauge fixes would be ineffective.

**The seed-accuracy and geometric problems are coupled**, not independent: you seeded from a close-but-wrong BVP solution, and the loose tolerance means the corrector never exits the near-tangent manifold-crossing region.

---

### Q2 — Min-energy anchor is strictly superior; here is the concrete path

Yes. The min-energy (quadratic cost) anchor avoids every pathology at once:
- No scale gauge (the quadratic cost makes the Hamiltonian definitively nonzero and uniquely normalized)
- S(τ) is smooth and strictly bounded away from 0 at the energy optimum; no switch-birth singularity
- dλ0/dt_f is bounded and nondegenerate along the energy-optimal arc family
- Your direct solver already uses Bertrand-Épénoy, so the min-energy KKT duals are better-conditioned than at min-fuel (the multipliers against a quadratic cost have smoother mesh dependence)

**Continuation parameter:** Use the Bertrand-Épénoy cost homotopy s ∈ [0,1]:

```
J_s = (1 − s) · ½∫(T/T_max)² dt  +  s · ∫(T/T_max) dt
```

At s = 0: smooth throttle with unique finite costates, no switches. At s = 1: bang-bang. The stationarity condition gives a smooth throttle T*(τ) = clamp(0, 1, −m·S / (2(1−s)·c)) for s < 1.

**Switch-structure locking to prevent basin collapse as s → 1:**

1. During s-continuation, monitor T*(τ) on each arc between existing switch nodes.
2. When T* > 0.95 uniformly on an arc → promote to burn arc; when T* < 0.05 uniformly → promote to coast arc.
3. At s = s_lock (when *all* arcs satisfy T* ∈ [0, 0.05] ∪ [0.95, 1.0] uniformly), freeze the switch structure, extract the costates, and hand off to your existing fixed-structure IFS.
4. s_lock is operationally determined by watching max_τ(|S|) on each arc type staying > 0.1 for at least two consecutive continuation steps.

You've built most of this on the direct side. The missing piece is extracting the IFS seed at s_lock from the continuation rather than from KKT duals.

---

### Q3 — Gauge fix for the min-time anchor (if you insist)

The scale-gauge fix (‖λ0‖ = 1) addresses the wrong direction as argued above. If forced to use this anchor, the minimal correct fix is:

**Add the min-fuel Hamiltonian value as an explicit equation.** Compute H0_ref from your direct solution's initial conditions (first arc), then augment the residual:

```
[R_min_fuel(λ0, t_f);  H(τ=0) − H0_ref]
```

This explicitly separates the manifold-landing correction (the H-equation) from the arclength advance, so the pseudo-arclength tangent's t_f-component is recoverable. The extended system is now (n+2) equations in (n+1) unknowns + arclength, yielding a clean 1D null space.

However: this still requires ‖R_min_fuel‖ < ~1e-10 at the anchor for the corrector to distinguish manifold-landing from fold-navigation. With your current ~1e-6 accuracy at the min-time anchor point, expect ‖R‖ to floor at ~1e-6 again. The fundamental issue is that you need the min-fuel costates at t_f,min — and those are not the min-time costates.

---

### Q4 — This is the highest-value thing to build next

Yes, KKT-dual-to-costate failure for a 40-rev indirect problem is expected. The KKT multipliers from an N-point collocation scheme are O(h) approximations to the true adjoints. Round-trip shooting sensitivity over 40 revs amplifies a relative error ε in λ0 by a factor of order exp(N_rev · σ_Lyapunov · t_rev). For a CR3BP spiral near the primary, the Lyapunov exponent in the state × costate space is large enough that even a 0.01% error in λ0 produces O(1) endpoint residuals. The KKT duals are typically accurate to 0.1%–1% depending on mesh density — squarely in the danger zone.

**Best first build: backward adjoint sweep on the direct state trajectory**

The procedure replaces the KKT dual extraction with a backward integration that is accurate to integration tolerance rather than mesh-spacing tolerance:

1. Accept the converged direct solution state history x_dir(τ) on its mesh (dense enough for accurate interpolation).
2. At τ = τ_f, initialize costates from the min-fuel transversality conditions: λ_r and λ_v from the terminal constraint multipliers (dual of the rendezvous equations in the direct NLP), λ_m(τ_f) = 0 exactly, λ_t(τ_f) = 0.
3. Integrate the *adjoint equations only* backward from τ_f to 0, using x_dir(τ) as the driven state. Do not re-integrate the state equations. Drive the adjoint ODEs with the direct solution's state via dense output interpolation. Use a tight tolerance (rtol ~ 1e-11).
4. Record the adjoint trajectory at every shooting node. These costates, by construction, satisfy the adjoint DEs and terminal transversality on the direct solution's state trajectory.
5. The IFS sees a residual of order O(direct_defect_size² + integration_error), which for a good collocation mesh is ~1e-8 to 1e-10 — well inside the basin.

**This is categorically different from KKT dual extraction.** You are not sampling the discrete Lagrange multipliers; you are re-integrating the continuous adjoint equations using the direct solution as the forcings. The only error sources are (a) direct solution defect quality and (b) adjoint ODE integration tolerance, both of which you control.

The other options (rev-by-rev continuation, a_T stepping) are more expensive and don't exploit the existing direct solution. Build the adjoint sweep first.

---

### Q5 — Hard throttle is correct given a sufficient seed; ε is for cold starts

ε-smoothing is necessary for cold seeds over 40 revs. It is not necessary when the warm-start lands inside the basin. With the backward adjoint sweep giving ~1e-9 residual seeds at each shooting node, your hard-throttle IFS with complex-step Jacobians and the existing Gauss-Newton/LM fallback should converge in 5–20 iterations. The "crawl" in your prior ε attempt was the symptom of cold-seed + tight basin, not a property of the bang-bang problem per se.

If the adjoint sweep alone isn't sufficient (residual stalls again but at a much lower level, say ~1e-4), the remaining gap is almost certainly the direct-solution defect at the shooting-node boundaries. In that case, first tighten the direct collocation tolerance before extracting the adjoint sweep, rather than introducing ε.

---

### Overall verdict on the IFS framing

The framing is sound. The failure is not architectural — it's that the costate seed quality bottleneck was misidentified as a geometric fold/gauge issue rather than a BVP-manifold separation issue. Fix sequence:

1. **Build the backward adjoint sweep** to replace KKT dual extraction → expect cold IFS convergence on existing direct solutions.
2. **If continuation across t_f is still needed** (gap in the direct solution band), switch to the min-energy anchor with s ∈ [0,1] homotopy → avoid the min-time anchor entirely.
3. **Drop the min-time anchor.** It is degenerate for two compounding reasons (wrong BVP and near-tangent manifolds at t_f,min), and no amount of gauge pinning recovers a clean continuation from it.
