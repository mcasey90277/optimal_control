## Mathematical correctness

### 1. **P0 — The mass-costate equation has the wrong sign**
**Eq. `costate`, lines 153–156; closed-form \(\dot Q\), lines 636–640.**

From eq. `H`,
\[
\dot\lambda_m=-H_m=\frac{sT}{m^2}\lambda_v\cdot\alpha
=-\frac{sT}{m^2}|\lambda_v|
\]
after minimizing over direction. Replace the last equation in `costate` with this **negative** expression.

This also resolves an internal contradiction: the stated
\[
\dot Q=-\frac{T\lambda_v^\top\lambda_r}{m|\lambda_v|}
\]
is correct with the corrected costate equation, because the mass-depletion and mass-costate terms cancel and the Coriolis matrix is skew-symmetric. With the printed positive sign, an additional \(2sT^2|\lambda_v|/(cm^2)\) remains.

### 2. **P0 — The min-time Hamiltonian normalization changes without being declared**
**OCP, lines 125–126; transversality, lines 172–174; terminal residual, line 331; lift argument, line 452.**

For the stated Mayer formulation \(\phi=t_f,\ L=0\), free-time transversality is \(H(t_f)+1=0\), not \(H(t_f)=0\). The latter, and the later statement \(\lambda\cdot f=-1\), use the equivalent running-cost formulation \(L=1,\ \phi=0\).

**Fix:** Declare immediately before eq. `H`:
> “From this point onward, min-time uses \(L\equiv1,\ \phi=0\); consequently \(H(t_f)=0\) and \(\lambda\cdot f=-1\).”

Also state the admissible physical domain \(m(t)>0\), excluding primary collisions. If a dry-mass bound is imposed, include it in the OCP and qualify \(\lambda_m(t_f)=0\) by its inactivity.

### 3. **P0 — The smoothing formulas need domain corrections, and huberc does not approach bang-bang as \(p\to0\) with fixed \(\delta\)**
**Smoothing families, lines 539–553; Fig. `families` caption, line 534.**

The **printed huberc cost is a consistent Legendre pair on its stated positive-\(Q\) branches**: at \(s=p\), both pieces have value \(p/2\) and derivative \(1\); their curvatures are \(1/p\) and \(\delta/(1-p)\). Thus it is \(C^1\) and strictly convex for \(0<p<1,\ \delta>0\). That part should be retained.

However:

- Both Huber laws print \(s=pQ\) for all \(Q<1\), giving inadmissible negative throttle when \(Q<0\). Such values can occur during shooting iterations even if absent on accepted extremals.
- For huber at \(Q=1\), the argmin is the interval \([p,1]\), not a unique throttle.
- At fixed \(\delta>0\), huberc has
  \[
  L(s)\longrightarrow s+\frac{\delta}{2}s^2,\qquad
  s^*(Q)\longrightarrow\operatorname{clip}\!\left(\frac{Q-1}{\delta},0,1\right),
  \]
  not the bang-bang law.

**Fix:** Add the branch \(s=0\) for \(Q\le0\), specify huber’s set-valued minimizer at \(Q=1\), and replace “all three tend to bang-bang as \(p\to0\)” with:
> “eps and huber approach bang-bang as \(p\to0\); huberc requires both \(p\to0\) and \(\delta\to0\), away from the indifference surface \(Q=1\).”

Likewise, qualify the ELI5 “never in between” claim at lines 95–98: it excludes singular intervals with \(Q\equiv1\).

### 4. **P0 — Multiple-shooting breakpoints are used as segment durations**
**Eq. `msres`, lines 333–338; Algorithm `ms_bvp`, lines 355–358.**

The algorithm defines \(\sigma=\mathrm{tGrid}/t_f\), so \(\sigma_k\) is a normalized **breakpoint**, but propagates segment \(k\) for \(\sigma_k t_f\). With a grid starting at zero, the first propagation has zero duration.

**Fix:**
\[
\tau_k=\frac{t_k}{t_f},\qquad
\beta_k=\tau_{k+1}-\tau_k,\qquad
\varphi_k(Y_k)=\varphi(\beta_k t_f,Y_k).
\]
Use \(\beta_k\), not \(\sigma_k\), in propagation and time derivatives. With \(F\) explicitly defined as the **14-state PMP vector field**, the complete blocks are
\[
R_{k,Y_k}=\Phi_k,\quad R_{k,Y_{k+1}}=-I,\quad
R_{k,t_f}=\beta_kF(y^-_{k+1}),
\]
and
\[
R_{K,Y_K}=Dg\,\Phi_K,\qquad
R_{K,t_f}=Dg\,\beta_KF(y^-_{K+1}).
\]
For the first junction, select only the seven free costate columns. The continuity residual itself is otherwise correctly structured.

### 5. **P0 — A root bracket ending at \(t_f\) is not an endpoint conjugate point, and an endpoint conjugate point does not establish a weak minimum**
**Conjugate instrument, lines 441–444; ELI5, lines 405–412.**

A sign change between the last two samples generally brackets a root **strictly before** \(t_f\). Calling this `ENDPOINT` and describing it as “a weak, non-strict minimum” is unjustified; even a root exactly at \(t_f\) does not by itself establish minimality.

**Fix:**
> “A sign change in any interval, including the final interval, brackets a candidate conjugate time. Refine the bracket before classifying it as interior or terminal. A terminal conjugate point makes the strict sufficiency test inconclusive; it does not establish a weak minimum.”

The ELI5 should describe an **infinitesimal focusing variation**, not necessarily a second finite trajectory that intersects the first. Its great-circle analogy is useful if it preserves the distinction between reaching the first conjugate point and passing it.

### 6. **P1 — The lift-space linear algebra is correct, but the identification with PMP abnormal lifts needs an additional argument**
**Abnormal-lift probe, lines 449–457.**

For the linear space \(S\) as defined and the nonzero functional \(\ell(\lambda)=\lambda\cdot f\), the normal member implies
\[
\dim\ker(\ell|_S)=\dim S-1.
\]
But \(S\) imposes only \(\lambda_v\parallel\alpha\), whereas a minimizing PMP lift also satisfies directional and throttle inequalities. A nonzero element of \(\ker\ell\) is therefore an abnormal **stationary** lift; the text has not shown that it, or its negative, is an admissible minimizing abnormal lift.

**Fix:**
> “\(\dim S=1\) excludes every nonzero stationary abnormal lift and hence every PMP abnormal lift. The converse requires proving that a nonzero kernel element satisfies the PMP minimum condition, or explicitly adopting the stationary-lift definition used by the cited theorem.”

Give that theorem’s precise definition. Also show how \(C\) is assembled—for example, by propagating the fixed-control adjoint fundamental matrix and stacking direction-alignment constraints and the terminal mass-costate row.

### 7. **P1 — The conjugate instruments need their reductions and applicability hypotheses stated locally**
**Free-time construction, lines 415–447; fixed-time construction, lines 697–706.**

The free-time displayed matrix needs \(P\in\mathbb R^{6\times5}\) spanning \(\lambda_{rv}(0)^\perp\), not the complement of the seven-component \(\lambda(0)\) described immediately above. Explain that the mass-costate column is discarded and mass has no fixed-time costate sensitivity on the strict all-burn branch. Also, the determinant is a scalar; it is the matrix inside it that belongs to \(\mathbb R^{6\times6}\).

**Fix:** Define those dimensions explicitly and summarize the image-space lemma and its hypotheses rather than relying entirely on the external audit.

The fixed-\(t_f\) choice \(\Phi_{1:7,8:14}\) is appropriate for an **interior fixed-full-state necessary-condition probe**; it should not be replaced by the terminal mixed block \([1{:}6,\ 14]\). However, obtaining the correct flow derivative does not alone establish a Jacobi theorem for constrained, piecewise-smooth or discontinuous feedback. State the required regularity, treatment of structurally singular intervals, and the switching/Jacobi result invoked for huber. Skipping initial singular samples must remain a diagnostic convention unless its validity is established for that arc structure.

The saltation formula itself has the correct sign for an autonomous surface, identity state reset and a transverse crossing. State \(n^\top F^-\ne0\) explicitly alongside it.

### 8. **P1 — The rank tolerance is a heuristic, and its cap contradicts its stated rationale in some cases**
**Rank rule, lines 455–457; spectral-gap claim, lines 460–464.**

The outer cap \(10^{-3}\sigma_1\) can force the tolerance below the residual-based term \(10\|C\lambda_0\|/\|\lambda_0\|\). Moreover, that residual is an empirical accuracy indicator, not a rigorous matrix-error bound; the reported gap alone does not establish the exact rank.

**Fix:**
> “This is a residual-informed numerical rank heuristic. If the residual-based scale exceeds the cap, or the singular-value separation is not robust to propagation accuracy and constraint scaling, report normality as UNDETERMINED.”

Document normalization of \(C\), and report \(\sigma_6,\sigma_7\), the selected tolerance and its sensitivity—not just the ratio.

## Logical consistency and data flow

### 9. **P1 — “Identifiability” confuses numerical conditioning with mathematical uniqueness**
**Data model, line 202; seed reconstruction, lines 576–578; reflight, lines 693–695; Fig. `files` caption, line 741.**

For a well-posed smooth flow—or a specified transverse hybrid flow—exact initial data determine the trajectory. Junctions do not create mathematical identifiability; they preserve a numerically usable representation of an unstable boundary-value solution. The document then says the deliverable is reconstructed by flying from its first junction alone, making the claimed role of all the junctions unclear.

**Fix:**
> “Store junctions to preserve numerical fidelity and support stable multiple-shooting reconstruction; initial-costate-only reconstruction can amplify rounding and integration errors.”

Specify separately the actual delivery/reconstruction procedure, the full-horizon reflight check, and the junction-based re-solve. Explain why energy records can safely discard junctions while fuel records cannot.

### 10. **P1 — Verdict binding omits data that change the trajectory and its STM**
**Binding, lines 707–710; verdict-file fields, line 736.**

The listed binding keys include \(\lambda_0,p,\) family and rows, but omit \(\delta\), which changes the huberc field. They also do not explain how the re-solved trajectory is reconciled with the packaged junctions; matching an initial costate alone is especially weak after the document’s conditioning discussion.

**Fix:** Bind the verdict to a complete solution/configuration identifier covering endpoints, \(T,c,t_f\), family, \(p,\delta\), junction data/grid and instrument version. State whether the re-solved junctions replace the catalog candidate or must pass an explicit agreement test.

## Honesty about numerical evidence

### 11. **P1 — The “proven versus assessed” distinction is not maintained in the strongest claims**
**“Certify,” lines 471 and 693–695; PASS licenses, line 764; acceptance, lines 370–375.**

Finite residuals and unvalidated propagation provide evidence for an approximate extremal, not a proof of an exact one. Dense samples do not establish inequalities “along the arc,” a numerical rank does not prove absence of abnormal lifts, and a sampled determinant does not establish absence of conjugate points—even “at junction resolution” can conceal the actual sign-based criterion.

**Fix the PASS row to say, approximately:**
> **Min-time:** “A numerically resolved normal PMP candidate, independently re-solved; sampled regularity and all-burn gates pass, numerical lift-space dimension is one, and the stated determinant detector finds no conjugate-point evidence.”  
> **Min-fuel:** “A numerically resolved extremal candidate of the specified smoothed problem, meeting the stated reflight tolerance, with no detected interior conjugate point under the sampled instrument.”

Also require successful foreign-solver termination and an independently evaluated terminal residual. Merely returning \(|\Delta z|<10^{-6}\) could mean that the foreign solver did not move. Define the norm and scaling of \(\Delta z\).

### 12. **P1 — Several conclusions are stronger than the measured experiments support**
**Multiple-shooting explanation, lines 307 and 317–321; huberc discussion, lines 568–571; future-work claims, lines 774–783.**

Multiple shooting avoids explicitly forming a long-horizon sensitivity product in each residual evaluation; it does not make the global problem necessarily well-conditioned or eliminate instability. Likewise, continuous huberc feedback still has derivative corners, observed success on two walls is not a general cure, and a switching-time Newton solve is not itself a second-order certificate.

**Concrete replacements:**
- “Multiple shooting reduces long-arc amplification within each residual evaluation; the assembled boundary-value Jacobian may remain ill-conditioned.”
- “huberc removes the throttle jump and passed the two tested Huber-wall cases.”
- “Switching-time refinement would provide a candidate on which the required reduced second variation and switching regularity conditions could be tested.”

Treat \(m_f\) agreement at \(10^{-6}\) as observed stabilization across accepted rungs, not an error bound relative to the exact fuel optimum. The catalog-wide conjugate counts at lines 458–460 should also be explicitly marked **legacy-instrument results**, not current corrected-instrument coverage.

## Reimplementation completeness

### 13. **P1 — The direct seed generation and endpoint/mesh construction remain black boxes**
**Dynamics, lines 102–114; direct collocation, lines 258–283; fuel pipeline, lines 484–487; seed and glossary, lines 576–578 and 794.**

A reader cannot reconstruct either direct NLP or the covector harvest from “Hermite–Simpson” and \(\lambda\approx\pm\mu/w\). The document omits the actual defect convention, midpoint control representation, objective quadrature, constraint scaling and scheme-specific dual map; the direct min-energy formulation is barely described at all.

Add one compact transcription subsection giving those equations and the exact mapping used. Also specify:

- How orbit phases map to fixed rotating-frame endpoint states, including phase origin and grid indexing; distinguish these from free-phase endpoints.
- How \(K\) is selected or increased, and whether breakpoints remain uniform.
- The nondimensional conversions for thrust and exhaust speed, and either explicit \(g(r)\) or a complete definition including the mass parameter and primary locations.
- Whether \(t_f^{\min}\) means a proven minimum or the accepted, numerically assessed min-time reference value.

These are essential algorithm inputs, not implementation trivia.

## Figures

### 14. **P1 — Several diagrams omit required inputs or imply stronger processing than the text describes**
**Figs. `mt-pipeline`, `mt-verify`, and `mf-pipeline`.**

- **Min-time pipeline, lines 229–241:** `tfMin` has no incoming solution/configuration arrow. Add \(z_8\), endpoints and engine parameters as its inputs, and show its convergence/residual output. Change “Every catalog entry has passed through … direct” to explain that later thrust rungs may use continuation from a directly seeded ancestor.
- **Min-time verification, lines 400 and 458–460:** “What a min-time entry has passed” suggests every entry passed the lower row, despite reported refuted/unverified entries. Label it “additional assessment and stored verdict,” with non-PASS outcomes shown.
- **Fuel pipeline, lines 486–487:** `ms_minenergy` does not have only seven unknowns once junctions are included. Label “seven initial costates plus internal junctions,” and show the missing direct-solution/dual-harvest-to-MS-seed edge. Its output at finite \(p,\delta\) should be labeled a **smoothed fuel candidate**.

## Overall verdict

The document has a strong organizing idea and contains valuable algorithmic details, especially the event-split saltation treatment, continuation strategy and explicit recognition that fixed-time interior conjugate testing is necessary-only. It does not yet achieve its stated purpose as a mathematically reliable, reimplementable companion: the mass-costate sign, Hamiltonian convention, segment-time definition and huberc limit need correction, while the abnormal-lift and conjugate-point conclusions need tighter hypotheses. After those repairs, a compact direct-transcription/seed specification and consistently numerical wording for verification would make this a substantially stronger senior-engineer review document.