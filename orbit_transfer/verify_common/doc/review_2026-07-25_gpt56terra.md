## Q1 — Switching-function recovery

- **[CORRECTNESS]** `foc_check.m:102-108` — Zeroing *only* the throttle-box rows is the right discrete recovery **provided the fuel objective is affine in throttle, the clock factor is positive, and those rows are the complete set of constraints acting directly on that throttle variable**.  With \(\mathcal L=f+s_{\!L}\lambda_g^Tg\), write
  \[
  0=\partial_{s_k}\mathcal L
   =S_k^d+\partial_{s_k}\mathcal L_{\rm box},\qquad
  S_k^d:=\partial_{s_k}\mathcal L_{\setminus{\rm box}}.
  \]
  Thus the active bound multiplier exactly cancels \(S_k^d\); removing its row recovers both the magnitude and sign of what it absorbed. For \(0\le s\le1\), minimization gives \(S^d<0\Rightarrow s=1\), \(S^d>0\Rightarrow s=0\).
  
  It is not valid as a bang-bang check for a weakly active bound (\(S^d\simeq0\)), fixed throttle represented by coincident lower/upper bounds (dual split is nonunique), unregistered duplicate/direct throttle constraints, or \(\varepsilon>0\) quadratic homotopy legs. There, this quantity is the **derivative of the smoothed discrete objective**, not an affine switching function: interior \(s\) can correctly have \(S^d=0\), and its sign does not imply the stated bound law. Add a manifest/runtime flag such as `throttleCostKind='affine'`; skip `signPct`, singular, and switch checks unless it is affine.

## Q2 — Defect duals and endpoint costates

- **[CORRECTNESS]** `foc_dual_to_costate.m:15-17` — The stated weighted average is not an exact discrete-adjoint identity, and, if \(\Lambda_k\) are regarded as interval-midpoint samples, its nonuniform weights are reversed. Linear interpolation from interval midpoints gives
  \[
  \widehat\lambda_k=
  \frac{h_k\Lambda_{k-1}+h_{k-1}\Lambda_k}{h_{k-1}+h_k},
  \]
  not the implemented \((h_{k-1}\Lambda_{k-1}+h_k\Lambda_k)/(h_{k-1}+h_k)\).

  More fundamentally, no exact nodal costate map can be obtained from only \(\{\Lambda_k,h_k\}\). For \(d_k=x_{k+1}-x_k-\frac{h_k}{2}(F_k+F_{k+1})\), with \(D_k=F_{x,k}\), state stationarity at an interior node is
  \[
  0=r_k+
  (I-\tfrac{h_{k-1}}2D_k)^T\Lambda_{k-1}
  -(I+\tfrac{h_k}2D_k)^T\Lambda_k ,
  \]
  where \(r_k\) includes nodal running-cost and other constraint contributions. The proper discrete adjoint boundary/costate quantities therefore require the dynamics Jacobian and those contributions; averaging is only a visualization/interpolation convention.

- **[CORRECTNESS]** `foc_check.m:146-150` — KNOWN-3 is real. `lam(:,N+1)=Lambda_N` is interval-centred to leading order, hence samples approximately \(t_f-h_N/2\), not \(t_f\). The cheapest useful continuous estimate is last-midpoint extrapolation:
  \[
  \widehat\lambda(t_f)=\Lambda_N+
  \frac{h_N}{h_{N-1}+h_N}(\Lambda_N-\Lambda_{N-1}),
  \]
  subject to smoothness and with mesh-convergence required. The cheapest **correct discrete endpoint covector** is instead
  \[
  p_{N+1}^{d}=(I-\tfrac{h_N}{2}D_{N+1})^T\Lambda_N
  \]
  plus terminal-objective and terminal-constraint terms under the selected sign convention. For a genuinely free final mass with no terminal mass cost/constraint, its mass component is zero by terminal-node KKT stationarity. Report this separately as a discrete consistency identity; use extrapolation only as a continuous-transversality estimator. Do not keep treating the one-sided interval dual as \(\lambda_m(t_f)\).

## Q3 — Mesh-normalized switch regularity

- **[CORRECTNESS]** `foc_check.m:158-170` — Replace the current statistic with a deweighted, divided-difference slope. Define nodal trapezoid weights
  \[
  w_1=h_1/2,\quad w_k=(h_{k-1}+h_k)/2,\quad
  w_{N+1}=h_N/2,
  \]
  and \(\widehat S_k=S_k^d/w_k\). If the clock factor varies by node, use \(\widehat S_k=S_k^d/(w_k\Xi_k)\); if it is a positive constant, omitting \(\Xi\) only changes the common scale. For a burn/coast transition between \(k_i,k_i+1\),
  \[
  D_i=\left|\frac{\widehat S_{k_i+1}-\widehat S_{k_i}}
                       {\sigma_{k_i+1}-\sigma_{k_i}}\right|,\qquad
  R_i=\frac{(\sigma_{N+1}-\sigma_1)D_i}{S_{\rm ref}},
  \]
  where \(S_{\rm ref}=\operatorname{median}\{|\widehat S_j|:j\text{ is outside a one-node switch neighbourhood and }|\widehat S_j|>\text{noise floor}\}\). Report \(\min_iR_i\). This is dimensionless and invariant to local mesh refinement for a smooth transversal crossing.

  No universal numerical threshold is physically defensible from one mesh. A provisional \(R_i>10^{-3}\) may be retained only after recalibration on deweighted data; it must be paired with refinement stability, e.g. require \(R_i\) to remain above a chosen problem-calibrated floor on two meshes. The present \(1.4\times10^{-7}\) and \(4.5\times10^{-5}\) cannot diagnose grazing.

## Q4 — Direction minimum condition

- **[CORRECTNESS]** `foc_check.m:114-123` — Removing the cone multiplier before projection does **not** make the tangent residual independent: for \(c(\beta)=\beta^T\beta-1\),
  \[
  P_\beta(q+2\mu\beta)=P_\beta q,\quad P_\beta=I-\beta\beta^T.
  \]
  It remains bounded by the full KKT residual when calculated from the same AD gradient and same dual vector. To make it operationally independent, construct \(q_k=\partial_{\beta_k}\mathcal L_{\setminus{\rm cone}}\) from an independent defect-dual/costate reconstruction or a separately reconstructed Hamiltonian, then compare it with the solved direction.

  It must also test the radial branch. With the explicit convention
  \[
  \mathcal L=\mathcal L_0+\mu(\beta^T\beta-1),\qquad q+2\mu\beta=0,
  \]
  and a positive thrust coefficient, the minimizer has \(\beta^Tq<0\), equivalently \(\mu>0\); the maximizer has \(\beta^Tq>0\), equivalently \(\mu<0\). Here \(\mu=s_{\!L}\lambda_{\rm cone}\) after resolving the global sign. If the registered row is \(1-\beta^T\beta=0\), reverse the multiplier-sign test. Apply it only on burns with \(\|q\|\) safely above noise, and report normalized anti-alignment \(-\beta^Tq/\|q\|\), not merely tangency.

## Q5 — Lagrangian-sign resolution

- **[ROBUSTNESS]** `foc_check.m:91-93` — Selecting the smaller of the two residuals is a useful diagnostic but not a safe convention resolver. If the true-plus residual is \(r=g_f+A^T\lambda\), the wrong-minus residual is \(2g_f-r\). For scalar \(g_f=1\), \(A^T\lambda=0.2\), the true residual is \(1.2\) while the wrong residual is \(0.8\), so the code chooses the wrong sign. At an exact solution the choice is unique only when \(g_f\) and \(A^T\lambda\) are materially nonzero; zero objective gradient/zero multiplier cases are intrinsically ambiguous.

  Determine and store the raw IPOPT/CasADi convention once, using a tiny CI-tested calibration NLP with a known nonzero active multiplier and recorded raw row orientation. At runtime, use that convention and fail/report “ambiguous dual convention” if the selected residual is not small **and** decisively smaller than its opposite.

## Q6 — Document/code drift

- **[CORRECTNESS]** `doc/first_order_checks.tex:190-198` — “State-column stationarity **is** the discrete adjoint recursion” is overstated. It is the full NLP state stationarity equation, which includes endpoint/path/boundary multipliers and only induces a discretized adjoint after defining a particular discrete costate. It is not automatically the trapezoidal continuous adjoint stated there.

- **[CORRECTNESS]** `doc/first_order_checks.tex:194-198` — Free-\(\Delta L\)/`cScale` stationarity is a discrete scalar NLP condition, not automatically the physical \(H(t_f)=0\) transversality condition. This is especially important for `cScale`; the code itself admits the value-form derivation is pending at `foc_check.m:178-181`.

- **[CORRECTNESS]** `doc/first_order_checks.tex:330-350` — The displayed dual map reproduces the incorrect nonuniform weighting in code and calls it a “standard” discrete-adjoint result. Replace it with the qualification and formulas in Q2.

- **[ROBUSTNESS]** `doc/first_order_checks.tex:300-304` — “Equals continuous \(S\) up to a positive factor” needs conditions: affine throttle cost, positive clock/span factor, complete removal of every direct throttle-bound row, and no weak/fixed bound degeneracy. It is false for the stated \(\varepsilon>0\) quadratic legs.

- **[CORRECTNESS]** `doc/first_order_checks.tex:404-423` — \(-\lambda_t=H\) and \(\lambda_t(t_f)=\pm1\) are not established by this implementation and are not generally valid for a Sundman-reparameterized system with `cScale`. They require a defined extended Hamiltonian, a specified costate convention, and a derivation of the time-state/cScale transformation. The code computes only a mapped defect-dual row and its CoV; it never computes \(H\), verifies this identity, or gates the endpoint value.

- **[ROBUSTNESS]** `doc/first_order_checks.tex:261-263` — The \(\sqrt3\,\texttt{kktStatInf}\) bound is correct only for exactly three direction rows and a unit direction. The generic code does not assert either condition before normalizing; document it as \(\sqrt{n_\beta}\,\texttt{kktStatInf}\) under \(\|\beta\|=1\).

- **[CORRECTNESS]** `doc/first_order_checks.tex:321-324` — The text says \(\varepsilon>0\) sign-law output is “advisory-only,” but `foc_check.m:187-192` folds `signPct`, singular arcs, and `sdotMinRel` into `rep.pass` without knowing \(\varepsilon\). Either skip these checks programmatically or remove that claim.

## Additional material defects

- **[CORRECTNESS]** `foc_check.m:158-165` — Singular-arc detection uses weighted \(S^d\), so its near-zero criterion changes with local mesh density even away from switching slopes. Run detection on \(\widehat S=S^d/(w\Xi)\), use a physical-\(\sigma\)-length threshold in addition to node count, and exclude a switch neighbourhood before classifying a singular run.

- **[ROBUSTNESS]** `foc_check.m:68-76,128-133` — The layout assertions establish only that the first presumed blocks numerically equal `X(:)` and `U(:)`; they do not validate defect-row ordering. `reshape(lamAll(defRows),nx,Nseg)` silently corrupts costates if registry rows are grouped/reordered rather than interval-major. Store `creg.defectRowsByInterval` as an `nx`-by-`Nseg` map and index through it.

- **[ROBUSTNESS]** `foc_check.m:69,75,129` — No assertion verifies `numel(sigma)==N1`, strictly positive mesh steps, `numel(xall)` capacity, valid manifest indices, or finite/nonzero direction norms. Add these before differentiation and before `b/norm(b)`.

- **[ROBUSTNESS]** `foc_check.m:122,191` — An all-coast throttle solution makes `max(tanAbs(burn))` empty, which can break the scalar advisory expression. Return `NaN`/skip the direction branch when there are no burn nodes.

- **[CORRECTNESS]** `foc_report.m:49-53` — The printer hard-codes defaults even when `foc_check(...,opts)` used different tolerances, so it can print PASS while `rep.pass` used a different gate. Save effective tolerances in `rep` and consume those.

- **[CORRECTNESS]** `foc_ipopt_inertia.m:6-15,74-81` — Zero IPOPT `delta_w` over a tail is not a local-minimum/SOSC certificate. It only says IPOPT’s linear-system inertia correction was not applied at those barrier iterates; it does not prove positive definiteness of the final reduced Hessian, account for barrier/active-set limiting behavior, or establish strictness. Relabel as an inertia-regularization observation, not `certLocalMin` / `LOCAL MIN`.
