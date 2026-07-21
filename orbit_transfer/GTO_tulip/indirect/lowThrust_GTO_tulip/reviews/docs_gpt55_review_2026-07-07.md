The core PMP Hamiltonian, switching sign convention, costate ODEs, and `lambda_m(tf)=0`, `H(tf)=0` conditions are mostly correct. Main issues are consistency and pedagogy.

1. **`gto_tulip_mintime_theory.tex` — Nondimensional units, lines 77–79 — severity: major**  
   **Issue:** The note says `\tilde c = I_{sp} g_0` “in ND”, but that is dimensional exhaust velocity.  
   **Fix:** Write  
   ```tex
   \tilde c = I_{sp} g_0\,t^\star/(1000\,l^\star) \approx 20.24 .
   ```  
   In the tutorial checkpoint, replace the algebraically equivalent but opaque expression with the same direct formula.

2. **`gto_tulip_mintime_theory.tex` — Direct route, lines 288–339; `building_the_gto_tulip_solvers.tex` — Phases D–E, lines 386–446 — severity: critical**  
   **Issue:** The theory note describes a 4-control NLP with variables `(w,s)`, throttle cone `w^T w=s^2`, `nZ=11(N+1)+1`, and mass flow `-Tmax s/c`. The tutorial actually builds a reduced full-thrust NLP with only a 3-vector unit direction, `w^T w=1`, `nZ=10(N+1)+1`, and constant mass flow `-Tmax/c`. These are different transcriptions.  
   **Fix:** Either make the theory note match the tutorial’s full-thrust reduced NLP, explicitly justified by `S<0` everywhere, or teach the full 4-control throttle NLP in the tutorial. Do not mix counts, constraints, and Jacobian blocks from both.

3. **`gto_tulip_mintime_theory.tex` — Direct route, lines 288–337; `building_the_gto_tulip_solvers.tex` — Phase E, lines 435–444 — severity: major**  
   **Issue:** The theory note derives defects for equal segments `h=tf/N`, including  
   ```tex
   \partial d_k/\partial t_f = -(f_k+f_{k+1})/(2N).
   ```  
   The tutorial uses nonuniform fixed mesh fractions `sigma_k`, so the actual formula is different.  
   **Fix:** Generalize the theory note to  
   ```tex
   h_k=t_f(\sigma_{k+1}-\sigma_k),\qquad
   \partial d_k/\partial t_f
   =-\frac{\Delta\sigma_k}{2}(f_k+f_{k+1}).
   ```  
   Mention that the equal-spacing formula is the special case `Delta sigma = 1/N`.

4. **`gto_tulip_mintime_theory.tex` — Geometry / terminal manifold, lines 111–117 and 158–169 — severity: major**  
   **Issue:** The documents repeatedly say “rendezvous on a tulip orbit,” but the OCP fixes one preselected tulip state. A true rendezvous with an unspecified tulip phase would have an additional phase/time manifold condition.  
   **Fix:** Say “rendezvous with a fixed sampled state on the tulip orbit.” If the intended target is any point on the tulip, add a phase variable and the corresponding transversality condition along the target-flow tangent.

5. **`gto_tulip_mintime_theory.tex` — Terminal-manifold remark, lines 158–169 — severity: minor**  
   **Issue:** “Transversality will pin the costate components normal to `M`” is misleading. The normal components `lambda_r, lambda_v` remain free multipliers; only the tangent mass component must vanish.  
   **Fix:** Replace with: “The costate may have arbitrary normal components associated with the six fixed rendezvous constraints; orthogonality to the free mass tangent gives `lambda_m(tf)=0`.”

6. **`building_the_gto_tulip_solvers.tex` — Introduction/results, lines 91–95 and 521–531 — severity: major**  
   **Issue:** “The two methods … check each other to five significant figures” conflicts with the reported direct result: `6.288574` vs `6.290694`, a relative difference of about `3.4e-4`, not five significant figures.  
   **Fix:** Change to “agree to mesh accuracy” or “agree to about `3e-4` relative at `N=12000`, improving under refinement.”

7. **`building_the_gto_tulip_solvers.tex` — Checkpoint E, lines 478–482 — severity: minor**  
   **Issue:** The text says the defect table converges at trapezoidal `O(h^2)` “per-segment” rate, but the table roughly drops by a factor of 8 per mesh doubling, consistent with local trapezoidal defect `O(h^3)`. Global trajectory accuracy is `O(h^2)`.  
   **Fix:** Say: “The interpolated warm-start defects decrease with the trapezoidal local `O(h^3)` defect scaling; the trajectory converges globally as `O(h^2)`.”

8. **`gto_tulip_mintime_theory.tex` — Ballast-exploit pitfall, lines 315–323 — severity: minor**  
   **Issue:** The pitfall is directionally right, but it should make clear that the exploit comes from decoupling thrust-vector magnitude from mass flow.  
   **Fix:** Add: “A physical throttle-vector model must either enforce `||w||=s` or use mass flow proportional to the actual thrust magnitude; the exploit appears when `s` can exceed `||w||`.”

9. **`building_the_gto_tulip_solvers.tex` — Phase B hint, lines 299–301 — severity: major**  
   **Issue:** “MATLAB’s relational `if S > 0` compares real parts of complex numbers” is brittle and version-dependent; it may fail or obscure the complex-step intent.  
   **Fix:** Tell learners to write  
   ```matlab
   if real(S) > 0
       u = 0;
   else
       u = 1;
   end
   ```  
   This is explicit and complex-step safe away from switches.

10. **`building_the_gto_tulip_solvers.tex` — Checkpoint E, lines 458–465 — severity: minor**  
   **Issue:** The sparse-gradient finite-difference checkpoint gives a verified number but no random seed, variable ranges, or mass/primary-avoidance safeguards, so learners cannot reproduce it exactly.  
   **Fix:** Provide a complete small `N=5` test script with fixed `rng`, positive masses, states away from primaries, random unit controls, and the exact reported tolerance.

11. **`gto_tulip_mintime_theory.tex` — KKT/costate remark, lines 360–366 — severity: minor**  
   **Issue:** “KKT multipliers … are the costates in disguise” is true only after the appropriate discrete covector scaling/sign convention, especially with trapezoidal defects and nonuniform `sigma`.  
   **Fix:** Add: “after the usual discrete covector mapping, including quadrature weights and sign conventions.”

12. **`gto_tulip_mintime_theory.tex` — PMP setup, lines 176–182 — severity: minor**  
   **Issue:** The note implicitly assumes the normal PMP multiplier convention `lambda_0=1`.  
   **Fix:** Add one sentence: “We use the normal convention `H=L+\lambda^T f`; abnormal cases are not considered for this transfer.”

## 3 highest-priority edits

1. Align the direct NLP formulation across both documents: 4-control throttle NLP vs reduced full-thrust 3-control NLP.  
2. Fix the nonuniform-mesh defect/Jacobian formulas in the theory note to match the tutorial.  
3. Correct the nondimensional exhaust-velocity formula and the overclaim that the direct/indirect answers agree to five significant figures.
