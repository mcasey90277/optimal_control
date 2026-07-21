## BOTTOM LINE

Proceed with the KKT-dual pivot, but not with the diagnosis as written. The forward-shooting objection is right; the “primal recovery is structurally ill-posed because endpoint rendezvous gives no λ boundary condition” is overstated. Method (2) is not genuinely doomed; it is a badly posed **penalized min-norm LS** for an almost homogeneous system. A hard-constrained/stabilized formulation can recover a costate direction if the trajectory/control data are PMP-consistent. The KKT route is sound and cleaner for certification, but the current code/comments do **not** give the correct node covector mapping; multipliers are interval KKT variables, not automatically continuous node costates.

- **[CORRECTNESS]** `GTO_tulip/sundman_minfuel/TIER1_PMP_CERTIFICATION_SCOPE.md:141-148` -- Method (2) is diagnosed too strongly. The failure of `A\b` only proves the **soft LS/min-norm formulation** is wrong/ill-conditioned, not that primal costate recovery is impossible. Concrete fix: enforce recursion as hard constraints and fit primer rows only in the nullspace. Form `R y = 0`, compute a stabilized basis `y = Z z` for the recursion manifold by multiple-shooting/QR/SVD or block Riccati, then solve
  ```text
  min_z || P_burn Z z ||^2
  s.t.  n' Z z = 1,  -alpha_k' lam_v,k >= 0 on burn nodes
  ```
  or solve a constrained generalized eigen/SVD problem for the smallest primer residual subject to normalization. This avoids the min-norm collapse.

- **[CORRECTNESS]** `GTO_tulip/sundman_minfuel/TIER1_PMP_CERTIFICATION_SCOPE.md:143-147` -- The “no boundary condition on λr, λv” statement is only half right. For fixed initial/final position/velocity, those endpoint state variations are fixed, so PMP imposes no transversality values on λr, λv at either end. But that does **not** make the inverse problem structurally ill-posed; the interior primer-direction constraints plus adjoint dynamics can define a well-posed overdetermined homogeneous BVP up to scale/sign. Concrete fix: rephrase as “endpoint transversality does not pin the costate scale/direction, so the chosen interior-data inverse problem must be solved with hard dynamics and explicit normalization/sign constraints.”

- **[CORRECTNESS]** `GTO_tulip/sundman_minfuel/certify_minfuel_pmp.m:141-148` -- `A\b` minimizes all recursion and primer rows simultaneously after row scaling. That allows violation of the adjoint recursion and lets the solver trade physical propagation for norm reduction. Concrete fix: do not penalize recursion. Either eliminate it (`Lam_k = Phi_k Lam_1` with stabilized products / periodic rescaling / multiple shooting) or solve equality-constrained LS:
  ```text
  min_y ||D y||^2
  s.t. R y = 0,
       q' y = 1,
       sign constraints if needed.
  ```

- **[CORRECTNESS]** `GTO_tulip/sundman_minfuel/certify_minfuel_pmp.m:127-138` -- Primer direction rows use `(I - alpha alpha') lam_v = 0` and one normalization `-alpha' lam_v = 1`, but do not enforce the sign at other burn nodes. Thus the wrong-sign solution can satisfy all projector rows and show primer error `2.0`. Concrete fix: add sign inequalities `-alpha_k' lam_v,k >= gamma_k > 0` on selected burn nodes or use signed residuals after hard recursion, e.g. `lam_v,k + rho_k alpha_k = 0`, `rho_k >= 0`.

- **[CORRECTNESS]** `GTO_tulip/sundman_minfuel/certify_minfuel_pmp.m:97-104` -- The hand-built adjoint map is not the discrete adjoint of the actual NLP. The NLP dynamics include thrust, mass, carried time, guarded radii, and kappa dependence (`casadi_minfuel_sundman.m:80-87`), while the certifier uses only the continuous 6x6 CR3BP λr/λv block. That is acceptable for a **continuous-PMP check along a frozen state**, but not for claiming dual consistency with the transcription. Concrete fix: separate claims: continuous adjoint verification uses `lt_pmp_eom_minfuel.m:73-75`; discrete KKT certification must use CasADi Jacobians of the exact defect `D`.

- **[CORRECTNESS]** `GTO_tulip/sundman_minfuel/TIER1_PMP_CERTIFICATION_SCOPE.md:169-170` -- “De-scale by trapezoid weight / τf / κ” is not precise and is partly misleading. For the defect
  ```text
  D_k = X_{k+1} - X_k - h_k/2 (F_k + F_{k+1}),  h_k = tauf Δσ_k,
  F = κ f,
  ```
  with Lagrange multipliers `ν_k` on `D_k = 0`, the discrete stationarity gives node covectors by endpoint sums, not by simply dividing each `ν_k` by `h_k κ_k`. Interior node `i` satisfies approximately
  ```text
  0 = w_i ∂L_i/∂X_i
      + ν_{i-1}
      - (h_{i-1}/2) F_X(i)' ν_{i-1}
      - ν_i
      - (h_i/2) F_X(i)' ν_i
      + boundary/path/bound terms.
  ```
  A consistent node covector proxy is therefore closer to a left/right or averaged defect multiplier with sign convention, e.g. `λ_i^- ≈ -ν_{i-1}`, `λ_i^+ ≈ -ν_i` depending on whether the augmented Lagrangian is `J + ν'D` or `J - ν'D`; the `hκ` factors already live inside the defect Jacobian. Validate sign by primer alignment.

- **[CORRECTNESS]** `GTO_tulip/sundman_minfuel/casadi_minfuel_sundman.m:191-210` -- Current dual interpretation is already overclaimed. `lamDef = reshape(lamAll(1:8*N),8,N)` returns interval defect multipliers, not node costates, and `lamMassEnd = lamDef(7,end)` is not the final mass transversality. Final mass transversality comes from stationarity of `X(7,end)` including the last defect multiplier, objective endpoint weight, and active bounds if any. Concrete fix: expose separate handles for defect, unit sphere, bounds, and boundary constraints; reconstruct stationarity/costates from the full KKT system.

- **[ROBUSTNESS]** `GTO_tulip/sundman_minfuel/casadi_minfuel_sundman.m:174-175` -- Using `opti.lam_g` as a flat vector relies on constraint ordering. It is fragile once constraints are inserted or changed. Concrete fix: store constraint expressions:
  ```matlab
  defectCon = D(:) == 0;
  unitCon = ...
  opti.subject_to(defectCon);
  ...
  out.nuDef = reshape(full(sol.value(opti.dual(defectCon))),8,N);
  ```
  Do not slice `lamAll(1:8*N)` unless also asserting the ordering.

- **[CORRECTNESS]** `GTO_tulip/sundman_minfuel/casadi_minfuel_sundman.m:110-111` -- The unit-sphere control constraint matters for KKT stationarity in α. The primer condition is not just “defect multipliers imply `alpha = -lam_v/||lam_v||`”; stationarity includes the sphere multiplier:
  ```text
  0 = objective_α + defect_α terms + 2η_k α_k
  ```
  Since objective has no α, projection orthogonal to α should recover primer alignment. Concrete fix: check the projected α-stationarity using defect multipliers and `η`, not only angle of raw `ν(4:6)`.

- **[CORRECTNESS]** `GTO_tulip/sundman_minfuel/certify_minfuel_pmp.m:159-171` -- The λm reconstruction assumes `λ_m(tf)=0` and then scales using switch nodes, but the NLP has fixed terminal position/velocity/time and free final mass with possible active mass bounds (`casadi_minfuel_sundman.m:114-123`). If `m(tf)` is near/at a bound or IPOPT bound multipliers are non-negligible, `λ_m(tf)=0` is false in the NLP KKT sense. Concrete fix: check final mass lower/upper bound multipliers and stationarity before using free-mass transversality.

- **[ROBUSTNESS]** `GTO_tulip/sundman_minfuel/certify_minfuel_pmp.m:80-81` -- Switch indices from `diff(s > 0.5)` identify interval boundaries, but scaling `W(swIdx)` uses the left node only. For bang-bang trapezoid controls, the switching condition should be checked by interpolation or bracketing across `[k,k+1]`. Concrete fix: estimate zero/switch locations on intervals and fit scale using both endpoint values or interval midpoint controls.

- **[CORRECTNESS]** `GTO_tulip/sundman_minfuel/TIER1_PMP_CERTIFICATION_SCOPE.md:157-160` -- IPOPT duals are not “immune” to dynamic range. Sparse KKT factorization is much better than forward shooting, but multipliers can still be affected by NLP scaling, active bounds, degeneracy at bang-bang controls, and nonunique multipliers. Concrete fix: state the KKT route is the most direct discrete certificate, then verify with stationarity residuals, primer projection residuals, switching signs, and mesh refinement/dual consistency.

- **[ROBUSTNESS]** `GTO_tulip/sundman_minfuel/casadi_minfuel_sundman.m:143` -- IPOPT gradient-based scaling can change reported multiplier scaling conventions unless handled carefully. CasADi usually returns multipliers in the original NLP convention, but this should be numerically verified. Concrete fix: finite-difference or AD-check stationarity:
  ```text
  ∇J(z) + J_g(z)' λ + bound multipliers = 0
  ```
  before interpreting λ physically.

- **[CORRECTNESS]** `lowThrust_GTO_tulip/lt_pmp_eom_minfuel.m:66-75` -- The continuous PMP equations are in physical time τ as used by that file, while the NLP uses `sigma` with `dX/dsigma = tauf κ f` (`casadi_minfuel_sundman.m:106-107`). Covectors themselves do not require division by `κ`; the independent-variable change multiplies the Hamiltonian, and the adjoint w.r.t. σ is `dλ/dσ = -tauf ∂(κ f)'/∂x λ - tauf ∂(κ L)/∂x`. Concrete fix: when comparing duals to continuous PMP, map by sign/node convention first; use σ-adjoint equations including κ derivatives for discrete residuals, or physical-time equations only after converting derivatives, not by arbitrary `κ` de-scaling.

- **[ROBUSTNESS]** `GTO_tulip/sundman_minfuel/TIER1_PMP_CERTIFICATION_SCOPE.md:162-173` -- The dual plan omits control-bound multipliers for `s ∈ [0,1]` (`casadi_minfuel_sundman.m:117-119`). At bang-bang nodes, throttle stationarity is an inequality-complementarity condition, not `S=0` except at switches/singular arcs. Concrete fix: use defect duals to compute the switching function and verify complementarity:
  ```text
  S_k > 0 => s_k = 0,
  S_k < 0 => s_k = 1,
  S_k ≈ 0 only near switch intervals,
  ```
  while also checking IPOPT bound multipliers have the corresponding sign.

- **[STYLE]** `GTO_tulip/sundman_minfuel/casadi_minfuel_sundman.m:49-51` -- The docstring says `.lamDef` are “discrete costates ... up to positive mesh-weight scaling and global sign.” That is too vague and currently wrong enough to mislead future certification. Concrete fix: rename to `.nuDef` until a verified covector reconstruction is implemented.
