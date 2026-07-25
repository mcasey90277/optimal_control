Here is the adjudication of the known issues and your priority questions, conforming strictly to the required review format.

- **[CORRECTNESS]** `foc_check.m:103` — **Q1. SWITCHING-FUNCTION RECOVERY**. Yes, it is mathematically CORRECT. 
  - **Derivation**: For a simple bound $s\in [0,1]$, KKT states $\nabla_s \mathcal{L}_{full} = \nabla_s \mathcal{L}_{sans} + s_{Lag} \lambda_{box} \nabla_s g_{box} = 0$. Since the constraint gradient is $\pm 1$, the gradient missing those rows precisely evaluates the marginal Hamiltonian physical cost $S_d$ which the bound multipliers absorbed, including the sign ($\pm S_d = -s_{Lag} \lambda_{box}$).
  - **Edge-failures**: The only condition where this silently fails is if *other* arbitrary constraints act directly on the throttle (e.g., system-wide rate path limits $\int \dot{s}^2$, inter-engine mixture constraints). For strictly independent lower/upper bounds, it is exact.
- **[CORRECTNESS]** `foc_dual_to_costate.m:45` — **Q2. DUAL-TO-COSTATE MAP**. The weighted formula mathematically perfectly matches the discrete adjoint of the trapezoidal rule precisely because of how KKT acts.
  - **Derivation**: Taking $\frac{\partial \mathcal{L}}{\partial x_k} = \Lambda_{k-1} - \Lambda_k - A_k^T(\frac{h_{k-1}}{2}\Lambda_{k-1} + \frac{h_k}{2}\Lambda_k) = 0$ and dividing by the span $w_k = \frac{h_{k-1}+h_k}{2}$ identically isolates the continuous $\dot{\lambda}_k$ via the relationship: $\frac{\Lambda_{k-1}-\Lambda_k}{w_k} = A_k^T \left( \frac{h_{k-1}\Lambda_{k-1} + h_k\Lambda_k}{h_{k-1}+h_k} \right)$.
  - **Endpoint / KNOWN-3 Bias**: Yes, relying on one-sided $\lambda_{N+1} = \Lambda_N$ causes the explicit failure because $\Lambda_N$ essentially prices at $t_f - h/2$, missing the endpoint by $O(h) \dot{\lambda}$.
  - **Fix**: Linearly extrapolate the interval duals based on their discrete geometric midpoints:
    `lam(:,N+1) = LamDef(:,N) + h(N)/(h(N-1)+h(N)) * (LamDef(:,N) - LamDef(:,N-1));`
    `lam(:,1) = LamDef(:,1) - h(1)/(h(1)+h(2)) * (LamDef(:,2) - LamDef(:,1));`
- **[ROBUSTNESS]** `foc_check.m:168` — **Q3. KNOWN-1 CONCRETELY**. $S^d_k$ is weighted by the local trapezoidal interval size $w_k$, scaling inconsistently along non-uniform/PSR densified meshes. The dimensionless regularity constraint should measure rate of change evaluated identically matching the mesh physics.
  - **Fix**: Non-dimensionalize `Sd` back to the continuous curve via node step-weights before applying finite differences:
  ```matlab
  h = diff(sigma(:).');
  w = [h(1)/2, (h(1:end-1) + h(2:end))/2, h(end)/2];
  S_unw = rep.Sd ./ w; % Unweighted representation of S
  coastS = abs(S_unw(~burn));
  scaleS = median(coastS(coastS>0)); 
  if isempty(scaleS) || isnan(scaleS), scaleS = max(abs(S_unw)); end
  
  swI = find(diff(burn) ~= 0);
  if ~isempty(swI)
      sdot = abs(S_unw(swI+1) - S_unw(swI)) ./ h(swI) / max(scaleS,1e-30); 
      rep.sdotMinRel = min(sdot); % Rate independent of mesh spacing
  end
  ```
- **[CORRECTNESS]** `foc_check.m:116` — **Q4. KNOWN-2 CONCRETELY**. Evaluating $v_{full}$ limits your tangential residual directly to `kktStatInf` by identity. To make it genuinely independent and accurately locate the minimum branch, you must compute $v_{unc}$ cleanly avoiding unrecorded constraint multipliers. Because the Manifest lacks a concrete index for the cone constraint, isolate it directly retaining purely dynamic defects:
  - **Fix**: 
  ```matlab
  lamDefOnly = zeros(size(lamAll)); lamDefOnly(defRows) = lamAll(defRows);
  g_def_only = gf + s*(A.'*lamDefOnly);
  for k = 1:N1
      b = U(man.dirRows,k);  b = b/norm(b);
      v_unc = g_def_only(uix(man.dirRows,k));
      tanAbs(k) = norm(v_unc - (v_unc.'*b)*b);
      % Physical branch check against primer vector phase
      if (v_unc.' * b > tolStat && burn(k)), warning('Maximizing direction branch identified'); end
  end
  ```
- **[ROBUSTNESS]** `foc_check.m:92` — **Q5. SIGN RESOLUTION**. Autoresolution is entirely robust. IPOPT barrier solves $f + s^\ast \lambda^T g=0$ to machine $O(10^{-14})$. The reversed incorrect condition necessarily evaluates identically to $2 \nabla f$ which is universally measurable and $\gg 0$ for non-trivial continuous boundary models (like minimum-time endpoints or bounded free-mass limits). It is $100\%$ safe. For a static alternative, `CasADi` natively formulates using $s = +1$; explicitly binding this negates the empirical check.
- **[STYLE]** `verify_common/doc/first_order_checks.tex` — **Q6. DOC-vs-CODE DRIFT**. 
  - _KKT discrete adjoint (stationarity is transversality)_: TRUE.
  - _Tangential bounded by $\sqrt{3} kktStatInf$_: TRUE. Maximum geometric extension of a 2-Norm compared to the $\infty$-Norm evaluation bounds of three dimensions is definitively identically $\sqrt{3}$.
  - _Min-time transversality limit_: TRUE.
  - _$S^d$ equals continuous function up to positive factor_: OVERSTATED/FALSE. While functionally strictly positive ($w_k > 0$), it obfuscates that the factor varies unpredictably node-over-node driving KNOWN-1 issues—masking geometric truth within finite differences.

**Loose Findings**
- **[ROBUSTNESS]** `foc_check.m:160` — Singular arc detection mechanism will routinely exhibit false-positive FAIL triggers in mesh densified layers resulting from uncorrected $S^d_k$ massing. Refined blocks mathematically drop $w_k$ bringing perfectly compliant regions actively beneath $10^{-6}\cdot \text{scaleS}$. Fix using identical unweight computation variables isolated explicitly in Q3 (`abs(S_unw) < max(1e-6*scaleS, 1e-14)`).
- **[STYLE]** `foc_check.m:185` — Typo in array definition assignment string logic; using `checks=[checks, {'singularArc','sdotRegular'}];` overrides column consistency if `checks` happens to establish as vertical over diverse iteration cycles. Recommend forcing row layout tracking continuously throughout.
