**Bottom line:** the all-burn mass elimination and the column-space lemma are sound. The identification with the BCT conjugate-time condition is not. The central problem is the passage from an autonomous problem with a terminal manifold to a **projected, non-autonomous endpoint map**. There is also a definite sign error in the mass adjoint.

- **[FINE] Mass reduction —** “For every admissible all-burn trajectory, \(m(t)=m_0-Tt/c\)” is correct; with \(t_f<m_0c/T\), eliminating mass is a smooth equivalence, and correctly pulled-back second variations—not merely feasible sets—agree.

- **[GAP] Mass reduction — the correct reduction does not establish the claimed theorem transfer.**  
  The unsupported passage is:
  > “Consequently the point-target, free-time second-order theory applies … Nothing below depends on which form is used.”

  Write \(k=T/c\), \(y=(r,v)\), and introduce the clock
  \[
  \tau=\frac{m_0-m}{k},\qquad p_\tau=-k\lambda_m.
  \]
  The reduced and extended Hamiltonians are
  \[
  K(\tau,y,p)=1+p\cdot F(\tau,y,u^*),\qquad
  H_{\rm ext}=K+p_\tau.
  \]
  Thus they do **not** simply “agree”: one includes a clock costate. The original terminal conditions become
  \[
  p_\tau(t_f)=0,\qquad K(t_f)=0.
  \]
  There is generally no condition \(K(0)=0\).

  The linearized terminal condition includes
  \[
  \delta p_\tau(t_f)+\dot p_\tau(t_f)\,\delta t_f=0,
  \]
  where the first term is the fixed-time variation. This condition is absent from the tangent-space construction used in Proposition 1.

  Importantly, **the flow column itself is not the error**:
  \[
  \partial_t y(t;p_0)=F(t,y(t),u^*(t))
  \]
  remains correct for a non-autonomous system. One should not append an ad hoc \(\partial_tF\) correction to that column. The issue is the parameter domain, terminal transversality, and the second-order theorem attached to the resulting map.

  **What would settle it:** derive the reduced free-time accessory problem, including terminal-time terms, and identify its degeneracies with the proposed determinant. Alternatively, retain the clock and use the appropriate terminal-manifold construction. Mass elimination can avoid an explicit focal construction, but it does not avoid its mathematical content.

- **[WRONG] Mass reduction/BCT identification — the projected determinant can acquire zeros that are not an interior loss of optimality for the final-target problem.**  
  The passage being attacked is:
  > “Nothing below depends on which form is used.”

  There is a concrete obstruction. For the correctly reduced Hamiltonian, put
  \[
  h(t,y,p)=p\cdot F(t,y,u^*).
  \]
  It is homogeneous of degree one in \(p\). Consequently, if
  \(J(t)=D_{p_0}y(t)\), then
  \[
  p(t)^\top J(t)=0.
  \]
  This follows directly by differentiating \(p(t)\cdot\delta y(t)\) along the Hamiltonian variational equations; it is initially zero and remains zero.

  Therefore, whenever \(\operatorname{rank}J=5\),
  \[
  \det[J P,F]=0
  \quad\Longleftrightarrow\quad
  p(t)\cdot F(t)=h(t)=0.
  \]
  In the ordinary autonomous normal point-target problem, this extra mechanism is excluded by \(h\equiv-1\). In your reduced problem it is not:
  \[
  \frac{d h}{dt}
  =-\frac{Tk}{m(t)^2}|\lambda_v(t)|,
  \qquad h(t_f)=-1.
  \]
  Thus \(h\) can cross zero before \(t_f\). None of H1–H4 excludes that.

  Here is an explicit counterexample to the proposed general transfer of theory. Consider
  \[
  \dot z=(2,1,0)+\frac1m u,\qquad
  \dot m=-1,\qquad u\in S^2,\quad m(0)=1,\quad z(0)=0.
  \]
  Set
  \[
  a(t)=\frac1{1-t},\qquad A(t)=-\log(1-t).
  \]
  The constant control \(u=-e_1\) gives
  \[
  z(t)=(2t-A(t),t,0).
  \]
  Choose the target \(z_f=z(0.9)\). Every all-burn competitor satisfies
  \[
  z_1(t)\ge 2t-A(t)>2(0.9)-A(0.9)
  \qquad(0\le t<0.9).
  \]
  Hence this trajectory is globally minimum-time. Equality at \(0.9\) forces \(u=-e_1\) almost everywhere, so it is unique.

  Nevertheless, the direction-parametrized endpoint map is
  \[
  z(t;q)=(2t,t,0)-A(t)q,\qquad q\in S^2.
  \]
  Its proposed free-time determinant is proportional to
  \[
  A(t)^2\bigl(2-a(t)\bigr),
  \]
  and has a sign-changing zero at \(t=1/2\), strictly before the globally optimal terminal time.

  This example has a normal lift \(p=e_1/8\), positive angular Legendre Hessian, and—with the correct mass-adjoint equation—
  \[
  \lambda_m(t)=\frac{10-a(t)}8,\qquad Q=\frac54>0.
  \]
  The zero comes from \(h(1/2)=0\), not from loss of rank of the fixed-time wavefront.

  This is **not a CR3BP counterexample**. It refutes the general inference that deterministic mass reduction makes this projected determinant the autonomous BCT test. The same possible extra factor \(h(t)\) exists in the CR3BP formulas.

  **What would settle it:** a problem-specific theorem distinguishing genuine accessory-problem degeneracies from these projection zeros. Better determinant sampling cannot repair this issue.

- **[FINE] Lemma and the displayed rank identity —** The column-space proof is correct, and using \(\operatorname{Im}A=\operatorname{Im}B\) to conclude \(\operatorname{rank}[A,f]=\operatorname{rank}[B,f]\) is entirely legitimate here.

- **[WRONG] Proposition 1 — the “so” after the rank identity does not follow from the stated BCT construction.**  
  The disputed passage is:
  > “so \(t\) is a BCT conjugate time iff the instrument’s \(6\times6\) determinant … vanishes.”

  Your stated autonomous BCT construction has a **seven-dimensional domain**:
  \[
  \mathbb R\times\mathcal L,\qquad \dim\mathcal L=6.
  \]
  It is structurally non-immersive at every time, not merely at isolated conjugate times.

  Indeed, let \(b=\dot m=-T/c\), and let \(e_m\) denote the mass-costate coordinate vector. The nonzero vector
  \[
  w=\lambda(0)+\frac1b e_m
  \]
  satisfies
  \[
  dH_0(w)=-1+\frac1b b=0,
  \qquad
  \Phi_{x\lambda}(t,0)w=0.
  \]
  Thus \(w\in T\mathcal L\) is a permanent kernel direction. Consequently,
  \[
  \operatorname{rank}
  [\Phi_{x\lambda}|_{T\mathcal L},f]\le6<7
  \]
  for every \(t\).

  Dropping the mass output and testing rank six does not preserve the original immersion condition. It creates a different map. Quotienting the permanent kernel can produce a meaningful six-dimensional parameter space, but you must then prove which variational problem that quotient describes.

  Also, the explanation
  > “the mass row of \(\Phi_{x\lambda}\) is identically zero”
  
  is not, by itself, the explanation for singularity of the **free-time** matrix: its mass entry in the flow column is \(b\ne0\). The permanent tangent kernel above is the relevant explanation.

  At the actual terminal time, imposing \(\lambda_m(t_f)=0\) can remove this extra parameter freedom. That may support a terminal boundary-map singularity test. It does not establish the claimed first-conjugate-time test along the entire arc.

  **What would settle it:** explicitly construct the correct six-dimensional domain, including terminal transversality, and prove that its singularities are those of the relevant second variation. The rank identity alone cannot do that.

- **[FINE] Quotient remark, narrowly construed —** The \(H=0\) slice **has been handled explicitly and correctly for the image identity**; scaling need not be a symmetry of the BVP, and no additional \(H=0\) argument is needed for that lemma.

- **[IMPRECISE] Quotient remark — which determinant is “the same”?**  
  The passage is:
  > “The determinant is therefore the same function of \(t\) up to a nonzero factor.”

  This can be made correct for a **defined quotient map**. Choose a fixed complement to the permanent kernel in \(T\mathcal L\), for example the vectors
  \[
  \left(Pw,-\frac{Pw\cdot f_{0,rv}}b\right).
  \]
  In this basis, the projected quotient differential is exactly the instrument’s matrix. Other fixed quotient bases change its determinant by a nonzero constant.

  But it is not a nonzero multiple of the original seven-state BCT determinant: that determinant is identically zero. The algebra establishes coordinate equivalence for the projected quotient map, not equivalence of that map with the required BCT variational construction.

  **What would settle it:** name the quotient and its endpoint map precisely, and separate this coordinate-invariance statement from the unproved theorem identification.

- **[WRONG] PMP equations — the mass-adjoint sign is wrong.**  
  The document states
  > “\(\dot\lambda_m=T|\lambda_v|/m^2\).”

  From the displayed minimized Hamiltonian,
  \[
  H=1+\lambda_r\cdot v+\lambda_v\cdot(g+h)
       -\frac{T|\lambda_v|}{m}-\frac Tc\lambda_m,
  \]
  the canonical equation is
  \[
  \boxed{\dot\lambda_m=-\frac{T|\lambda_v|}{m^2}}.
  \]
  With the printed positive sign, the asserted Hamiltonian conservation fails:
  \[
  \frac{dH}{dt}=-\frac{2T^2}{c\,m^2}|\lambda_v|
  \]
  along the stated all-burn equations.

  With the corrected sign and \(\lambda_m(t_f)=0\),
  \[
  \lambda_m(t)=
  \int_t^{t_f}\frac{T|\lambda_v(s)|}{m(s)^2}\,ds\ge0.
  \]
  This also supplies the natural all-burn PMP argument: \(Q_{mt}>0\) whenever \(|\lambda_v|>0\).

  **What would settle it:** correct the equation and check every assertion involving Hamiltonian conservation, clock-costate identification, and switching-function values. The \(rv\) trajectory invariances survive because \(\lambda_m\) is decoupled, but the claimed Hamiltonian lift does not survive the printed sign.

- **[GAP] Scope of the mass reduction — all-burn competitors are not all competitors of the stated problem.**  
  The initial problem allows \(s\in[0,1]\), whereas the reduction concerns:
  > “every admissible all-burn trajectory”.

  For the original admissible class, a first-order variation gives
  \[
  \delta m_f
  =-k\,\delta t_f-k\int_0^{t_f}\delta s(t)\,dt,
  \]
  with one-sided throttle variations at \(s=1\). Thus the extra mass variation really exists outside the all-burn restriction.

  Strict \(Q>0\) along the reference does not make these competitors disappear. It can support a strict-bang sufficiency argument. Alternatively, an all-burn theorem for relevant minimizers, together with the needed existence argument, can justify restriction. Neither is equality of feasible sets.

  **What would settle it:** state whether the certificate concerns the restricted all-burn problem or the original throttle problem, and prove the transfer between them.

- **[GAP] Sufficiency theorem and hypothesis dictionary — not a complete, usable statement of BCT 2007.**  
  The passage is:
  > “BCT 2007 Thm. 3.x … reads: if the strong Legendre condition holds … and there is no conjugate time … then … strict strong … local minimizer.”

  **Strong \(C^0\)-local sufficiency before the first conjugate time is the appropriate kind of conclusion in the regular autonomous setting.** “Strong \(C^0\)” is not, by itself, the objection. But the document does not supply a faithful complete theorem statement: “Thm. 3.x” is a placeholder, and the regularity/corank framework and endpoint-map assumptions are material.

  In particular, these are different assertions:
  1. a normal multiplier exists;
  2. no abnormal minimizing multiplier exists;
  3. the relevant endpoint-map singularity has corank one, with the interval-wise regularity required by the conjugate-point theory.

  The table treats them too loosely. A full-horizon probe with terminal mass transversality is not automatically the corank condition for an autonomous seven-state **point-target** endpoint map—indeed that map has the permanent degeneracy exhibited above. Nor is full-horizon information, without an argument, interval-wise information.

  Another concrete mismatch is that autonomous normality supplies \(p\cdot f=-1\) for the state space to which the theorem is applied. Your six-state projection does not retain that property.

  **What would settle it:** give the exact BCT theorem number and statement, including the control class, corank assumptions, endpoint conditions, and topology for variable-duration trajectories. Then provide a hypothesis-by-hypothesis mapping to the **correct reduced construction**. As written, I would not accept the advertised sentence as a complete statement and application of BCT 2007.

- **[WRONG] Normality subsection — abnormality does not require \(\lambda_v\equiv0\), and the kernel argument omits PMP inequalities.**  
  The first false passage is:
  > “an abnormal lift … would require \(\lambda_v\equiv0\)”.

  Abnormality means the **cost multiplier** is zero. It does not imply that the control costate vanishes. An abnormal extremal can have \(|\lambda_v|>0\) and obey the same unique direction law.

  The second overclaim is:
  > “the kernel is exactly the abnormal lifts … no abnormal lift exists iff \(\dim S=1\).”

  Your \(S\) imposes linear adjoint and stationarity conditions, but only
  \(\lambda_v\parallel\alpha\). A minimizing lift also needs
  \[
  \lambda_v=-\rho\alpha,\quad \rho\ge0,
  \qquad Q_{mt}\ge0,
  \]
  together with the other PMP requirements. A nonzero vector in the kernel of \(\lambda\mapsto\lambda\cdot f\) need not satisfy those inequalities, even after one global sign change.

  Therefore \(\dim S=1\), with a normal element having nonzero functional value, is a **sufficient algebraic exclusion** of abnormal lifts. The claimed converse is not established.

  **What would settle it:** distinguish stationary multipliers from minimizing PMP multipliers, and test the kernel against the appropriate cone rather than only its dimension.

- **[GAP] Meaning of PASS — an initial interval of genuine rank deficiency cannot simply be skipped.**  
  The passage is:
  > “an initial arc where the state block is structurally rank-deficient having been skipped”.

  Rank deficiency at \(t=0\) is expected. Genuine rank deficiency on a positive-time interval would, under the document’s own definition, mean singularity throughout that interval. Skipping it is incompatible with claiming that the only remaining issue is sparse sampling.

  **What would settle it:** distinguish small-time numerical ill-conditioning from exact deficiency and provide a short-time rank/asymptotic argument covering the omitted interval.

The defensible algebraic conclusion is that the instrument tests singularities of a particular **projected, projectivized extremal endpoint map**. The document proves a useful representation of that map’s differential. It does **not** prove that its zeros are exactly the BCT variational conjugate times for the stated transfer problem.