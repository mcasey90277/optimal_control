No CRITICAL findings.

- **[MAJOR]** `gto_tulip_mintime_theory.tex:459-470` (§6, “the ‘+1’ anchors...”) - “In min-time, \((\lambda,H)\) can be rescaled” conflicts with the earlier normal convention \(H=1+\lambda^T f\), \(\lambda_0=1\). Fix by scoping the statement: the min-time **control law/switching sign** is costate-scale invariant, while the normalized full TPBVP scale is fixed by \(H(t_f)=0\); min-fuel makes even the throttle law scale-sensitive via the bare \(+1\).

- **[MAJOR]** `gto_tulip_mintime_theory.tex:432-435` (§6, “we measured both methods failing...”) - The later “honest accounting” documents indirect seed stalls, but not a direct-method failure on the full min-fuel spiral. Fix by either removing “both methods” or adding the direct full-spiral evidence: mesh, warm start, exit flag, feasibility/objective behavior.

- **[MAJOR]** `building_the_gto_tulip_solvers.tex:603-614` (Phase G checkpoint) - The phase is introduced as the arrival-leg min-fuel problem, but the checkpoint suddenly asks for full-transfer min-fuel shooting without a driver/signature or setup path. Fix by adding a `solve_minfuel_indirect.m` exercise/signature and explicit full-transfer inputs, or mark this as an optional diagnostic after the leg machinery exists.

- **[MAJOR]** `building_the_gto_tulip_solvers.tex:667-682` (Checkpoint H2) - `costate_seed_from_nlp` is expected but never assigned as an exercise before it is used. Fix by adding a build step with signature and required rows: linear costate STM along NLP arc, primer-direction cross-product equations, \(\lambda_m(t_f)=0\), and switch-based sign/scale anchoring.

- **[MAJOR]** `building_the_gto_tulip_solvers.tex:695-700` (“What comes next”, item 1) - “Either should convert... into a converged indirect solution” overpromises an open problem. Fix to “natural escalation paths to try/expected to improve the basin”; do not imply verified convergence.

- **[MINOR]** `building_the_gto_tulip_solvers.tex:147-154` vs `623` - Min-fuel filenames are inconsistent/vague (`solve_minfuel_nlp.m + driver` vs `NLP_lowThrust_GTO_Tulip_minfuel.m`), and several named direct min-fuel reference files are not present in this checkout. Fix names and say explicitly which are learner-created vs provided references.

- **[MINOR]** `building_the_gto_tulip_solvers.tex:625` - “control \(u=[\mathbf w;s]\)” collides with scalar throttle \(u\) used throughout PMP. Fix by naming the 4-vector \(q\), \(U_k\), or `ctrl`, reserving \(u\) or \(s\) for throttle.

- **[MINOR]** `gto_tulip_mintime_theory.tex:496` - “lighter is better here” is imprecise for a min-fuel objective, which rewards larger final mass. Fix: “lighter is dynamically advantageous, even though the objective penalizes propellant loss.”

- **[MINOR]** `building_the_gto_tulip_solvers.tex:623` and `652` - Very long `\texttt{...}` strings in prose/tcolorbox titles are overfull-prone. Fix with displayed `verbatim`, shorter aliases, or line-broken filename lists.

**Verified correct**

- Min-fuel cost \(J=\int (T_{\max}/c)u\,dt=m(0)-m(t_f)\).
- Hamiltonian coefficient and switching function \(S=1-\|\lambda_v\|c/m-\lambda_m\); burn for \(S<0\), coast for \(S>0\).
- Primer direction \(\alpha=-\lambda_v/\|\lambda_v\|\).
- Costate ODE signs, including \(\dot\lambda_m=-\|\lambda_v\|uT_{\max}/m^2\).
- Fixed-\(t_f\) transversality: six rendezvous residuals plus \(\lambda_m(t_f)=0\), no \(H(t_f)=0\).
- Bertrand–Epenoy smoothing \(u_\varepsilon=(1-\tanh(S/(2\varepsilon)))/2\).
- Direct min-fuel transcription: fixed \(t_f\), controls \([\mathbf w;s]\), cone equality \(\mathbf w^T\mathbf w=s^2\), objective \(-m_{N+1}\), no \(t_f\) Jacobian column.
- Quoted leg results match the verified values: \(N=3000\), max defect \(\sim2\times10^{-15}\), propellant \(1.0622\) kg vs \(1.0650\), burn fraction \(76.9\%\), switch at \(2.2907\) ND.
