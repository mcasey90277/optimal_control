## Errors

### E1. Fixing both \(\tau_f\) and physical \(t_f\) changes the optimization problem

**File:** `direct/lib/casadi_minfuel_sundman.m`, lines 13–18, 133–138  
**Quote:** “The total regularized length tau_f is held FIXED”; “the trajectory adjusts so that int(kappa dtau) = tf.”

The physical-time endpoint is enforced correctly, but this is **not generally equivalent to the stated fixed-physical-time minimum-fuel problem**.

Once \(dt/d\tau=\kappa(x)\) is prescribed, every physical trajectory determines its own regularized duration:
\[
\tau_f=\int_0^{t_f}\frac{dt}{\kappa(x(t))}.
\]
Holding that duration fixed adds the isoperimetric constraint
\[
\int_0^{t_f}\kappa(x(t))^{-1}\,dt=\tau_{f,0}.
\]
Different physically admissible trajectories need not have the same integral. Enforcing terminal time does not remove this extra restriction.

Consequently, machine-tight defects and excellent discrete KKT residuals can certify the **extra-constrained NLP**, without certifying the advertised physical problem.

**Fix:** Use the sparse constant-scale-state construction already implemented in the min-time solver, but retain the fixed terminal physical time:
```matlab
% x = [r; v; m; t; cScale], with cScale free and positive
cs = x(9);
fdyn = Function('f', {x,u}, ...
    {[cs*kappa*[v; accel; mdot; 1]; 0]});

gint = Function('g', {x,u}, ...
    {[cs*kappa*s; cs*kappa*s*(1-s)]});

% Keep:
opti.subject_to(X(8,end) == tf);
% Do not fix X(9,1) or X(9,end).
```
Keep `tauf0` fixed as a coordinate-domain length. Update state sizes, bounds, defect extraction and manifests to nine states. Verify that the artificial `cScale` bounds are inactive.

Alternatively, free scalar `tauf`; sparsity concerns justify a different implementation, not fixing a physically trajectory-dependent quantity.

---

### E2. The FOC gate does not verify the signed primer minimum condition or complete KKT conditions

**Files:**

- `verify_common/foc_check.m`, lines 146–155, 300–309
- `direct/lib/casadi_minfuel_sundman.m`, lines 237–254
- `direct/certify/run_foc_tulip.m`, lines 100–107

**Quotes:**

- “minimum condition, direction part: tangential dL/dbeta”
- `v = gL(uix(man.dirRows,k));`
- `if mean(cang) < 0, cang = -cang; end`
- “mirrors cert.passed's own primer/sign-law thresholds”

There are three distinct problems.

1. **Tangential stationarity is not Hamiltonian minimization.** Both parallel and antiparallel steering give zero tangential residual. In fact, here it is a projection of the already assembled full stationarity residual. A value such as \(5\times10^{-18}\) is not a normalized primer-direction error and cannot be compared directly with the LS direction error of 2.

2. **The engine chooses the costate sign from the observed steering.** For
   \[
   \mathcal L=J+\sum_k\mu_k^\top D_k,\qquad
   D_k=X_{k+1}-X_k-\cdots,
   \]
   the conventional minimizing-Hamiltonian costate has the sign opposite to the effective defect multiplier. That sign should follow from the Lagrangian convention, not be flipped to improve the primer result. The current procedure can conceal steering toward a Hamiltonian maximum.

3. **`rep.pass` does not require primal feasibility, dual feasibility or complementarity.** Nor does it itself require successful optimization. The wrapper's guard partly mitigates this for this particular call chain, but the generic gate can award an advisory PASS without full KKT conditions.

The cross-check thresholds also are not equivalent: mean angle below \(1^\circ\) is substantially weaker than maximum unit-vector error below \(10^{-2}\), approximately \(0.573^\circ\). `cert.passed` additionally includes switch matching and scale spread.

**Fix:**

- Retain full-Lagrangian stationarity as one check.
- Add feasibility, signed inequality-multiplier and complementarity checks using `opti.g`, `opti.lbg`, `opti.ubg` and sign-resolved multipliers.
- Compute the direction coefficient with the direction-normalization and direction-box contributions removed. For a nonzero coefficient \(q_k\), check
  \[
  \alpha_k=-q_k/\|q_k\|,
  \]
  not merely its tangential projection.
- Resolve the multiplier convention from the NLP Lagrangian, once.
- Compare maximum signed, normalized direction errors on matching node masks. Retain midpoint angles only as diagnostics.

The full AD assembly `gf + s*A.'*lamAll` does include all constraint rows; the incomplete registry does **not** truncate that calculation. That part is correct.

---

### E3. The reported mapped terminal covector includes terms that can cancel a nonzero physical transversality residual

**File:** `verify_common/foc_check.m`, lines 181–194, 213–217  
**Quote:** “gL already sums any OTHER constraint touching X(massRow,N+1) (an active box bound, say)”

That is precisely why the complete `gL` entry is **not**, in general, the mapped terminal covector.

Let \(h=tauf\,\Delta\sigma_N\), \(F=dX/d\tau\), \(\ell\) be the running integrand in \(\tau\), and \(\mu_N=s\,\Lambda_N\). With the defect convention used here, the conventional terminal costate is
\[
p_f
=-
\left(I-\frac h2F_x(X_f,U_f)\right)^\top\mu_N
-\frac h2\ell_x(X_f,U_f).
\]
In the absence of other terminal dependence, stationarity gives \(p_f=\Phi_x\). For free final mass and no terminal mass cost, its mass component is zero.

But if a mass box is active, full stationarity can be zero because the box multiplier cancels a **nonzero** mapped mass costate. Reporting that cancellation as free-mass transversality earns an inappropriate PASS.

Also, the comment “CANNOT fail unless kktStat already has” is false for the **relative** test: division by a small costate scale can make it fail even when the absolute stationarity test passes.

**Fix:** Assemble the terminal mass contribution from the running objective and defects only, and separately verify that final mass really is free:
```matlab
mu = s * lamAll;
gObjDef = gf + A(defRows,:).' * mu(defRows);
imf = (N1-1)*nx + man.massRow;

rep.lamMassEndMapped = abs(gObjDef(imf)) / max(scale,1e-30);
```
This expression is appropriate for these two solvers because neither objective has a terminal **mass** cost. A generic endpoint covector implementation must explicitly separate Mayer and running costs.

Report terminal box slack and multiplier contributions. Do not label a bound-constrained final mass “free”.

---

### E4. The LS certifier's disagreement is not established to be merely accumulated trapezoidal error

**File:** `direct/certify/certify_minfuel_pmp.m`, lines 7–12, 95–115, 173–206  
**Quotes:**

- “only a 2nd-order … discretization of the CONTINUOUS adjoint”
- `Mst(:,:,k) = [zeros(3), -G; -I3, -Hct];`
- “scale normalization rho(kb)=1”
- `qint = kappa .* (nlamv .* s .* Tmax ./ m.^2);`

The CR3BP six-costate matrix has the correct signs for the **ordinary physical-time adjoint**, and the stated trapezoidal recursion is a correct discretization of
\[
\frac{d\lambda_{rv}}{d\tau}=\kappa M\lambda_{rv}.
\]
However, that is not generally the adjoint system of the actual fixed-\(\tau_f\) NLP.

For the extended Sundman Hamiltonian
\[
K=\kappa\bigl(\ell_{\rm phys}+p_y^\top f+p_t\bigr),
\]
the position adjoint contains
\[
p_r'
=-\kappa G^\top p_v
-\nabla_r\kappa\,
  \bigl(\ell_{\rm phys}+p_y^\top f+p_t\bigr).
\]
The second term vanishes for the correctly free-regularized-horizon representation of the original problem, where \(K=0\). It need not vanish in the current fixed-\(\tau_f\) formulation. The LS fit omits it.

There are additional issues:

- **`rho(kb)=1` does not hard-normalize the costate.** Its relation to \(\lambda_v\) is only penalized. A zero costate remains feasible for the equality and bound constraints, with a nonzero fitting residual.
- **The mass reconstruction assumes the primer fit succeeded.** Generally,
  \[
  \dot p_m=\frac{sT_{\max}}{m^2}p_v^\top\alpha.
  \]
  Replacing \(p_v^\top\alpha\) by \(-\|p_v\|\) is invalid on arcs where the recovered direction disagrees with the control.
- **`W(swIdx)` samples the left node of a switching interval**, not a continuous switch at which \(S=0\).
- `lsqlin`'s exit status is ignored; `cert.passed` does not gate the recursion residual.
- Solving a sparse global BVP avoids explicit products of unstable transition matrices, but does not make the problem “immune” to conditioning.

**Fix:**

1. First decide which problem is being certified. For the current NLP, reconstruct the full extended adjoint, including \(\nabla\kappa\) terms. For the intended physical problem, fix E1 and then validate a physical-time reconstruction under mesh refinement.
2. Normalize the actual covector, for example:
   ```matlab
   % Add a hard equality:
   -alpha(:,kb).' * Lam(4:6,kb) == 1
   ```
   implemented as a row of `Aeq`.
3. During diagnostic reconstruction, use the actual dot product:
   ```matlab
   pmDot = s .* Tmax ./ m.^2 .* sum(lamv .* alpha,1);
   % Integrate pmDot backward with pm(tf)=0.
   ```
4. Locate/interpolate switching times, or use interior-throttle stationarity nodes, and quantify the mesh error in scale recovery.
5. Gate solver status, finite positive scale, recursion residual and nonzero primer norm.

The normalized switching formula itself,
\[
S=1-\frac{c_{\rm Ex}\|p_v\|}{m}-p_m,
\]
is sound when costates are normalized to running cost \((T_{\max}/c_{\rm Ex})s\), and primer alignment holds. The global switch-based scale can implement that normalization.

**Bottom line:** The LS route is fundamentally reasonable as a multiple-shooting-like consistency diagnostic. This implementation's failure cannot be attributed confidently to “40 revolutions plus second-order truncation” without fixing the model mismatch and performing a refinement study.

---

### E5. The time-costate interpretation is wrong or incomplete; `cScale` does not explain a 5% variation

**Files:**

- `verify_common/foc_check.m`, lines 170–175, 283–290
- `GTO_ELFO/direct/elfo/casadi_mintime_freetf.m`, lines 104, 138–143

**Quotes:**

- “fixed t_f: H=const generally nonzero; constancy via lamTimeCoV”
- “value-form H(tf) check reported via lamTimeEnd … derivation pending”

For the min-time solver, `t` does not occur in the dynamics. With no active time bounds, discrete stationarity therefore gives, **exactly**,
\[
\mu_{t,k-1}-\mu_{t,k}=0,\qquad
1+\mu_{t,N}=0.
\]
Thus the sign-resolved **defect multipliers** satisfy \(\mu_t=-1\). In the conventional minimizing-Hamiltonian costate convention \(p=-\mu\), the Mayer endpoint condition is instead
\[
p_t(t_f)=+1.
\]

Neither statement alone is the free-physical-time condition \(H_{\rm phys}(t_f)=0\).

For this formulation,
\[
K=c_s\kappa(p_y^\top f+p_t),\qquad
p_{c_s}'=-\kappa(p_y^\top f+p_t).
\]
With `cScale` free at both endpoints and its bounds inactive, \(p_{c_s}(0)=p_{c_s}(\tau_f)=0\). Together with autonomous-Hamiltonian constancy, this yields \(K=0\), hence
\[
1+p_y^\top f=0,
\]
the conventional physical min-time Hamiltonian condition.

**Fix:**

- Report `lamTimeDefectEnd = -1` explicitly in defect-multiplier convention, not ambiguously as the physical costate.
- Check raw interval time multipliers for constancy.
- Check the `cScale` endpoint/adjoint conditions and inactive scale bounds.
- Treat a direct \(H_{\rm phys}=0\) reconstruction as an additional mesh-dependent check.

A 5% variation can result from an incorrect node mapping, active/artificial time-bound contributions, accumulated stationarity error or numerical conditioning. It is **not an intrinsic effect of introducing the constant scale state**. A correctly normalized weighted average of constant interval multipliers also remains constant on a nonuniform mesh.

`foc_dual_to_costate.m` was not supplied, so its mapping cannot be verified here. That missing dependency is particularly important for the observed CoV.

For the existing fixed-\(\tau_f\) min-fuel problem, constant \(p_t\) also does not establish constant ordinary physical Hamiltonian: generally \(K/\kappa=H_{\rm phys}+p_t\), with \(K\) potentially nonzero. This is another consequence of E1.

---

### E6. “Apoapsis” is implemented as minimum rotating-frame speed

**File:** `cr3bp_common/insertion_states.m`, lines 67–71  
**Quote:** `[~, idx] = min(vecnorm(tr(:,4:6), 2, 2)); % slowest point`

Minimum rotating-frame speed is not the definition of lunar apoapsis and need not coincide with it in CR3BP.

**Fix:** If lunar apoapsis/apolune is intended:
```matlab
rMoon = [1-p.muStar, 0, 0];
rhoMoon = vecnorm(tr(:,1:3) - rMoon, 2, 2);
[~, idx] = max(rhoMoon);
rvf = tr(idx,1:6);
```
Refine the sampled maximum if an accurately located event is required.

Otherwise rename the criterion and label to `minRotatingSpeed`; do not call it apoapsis.

---

### E7. The seed mapper can return a zero “unit direction”, giving a singular normalization constraint at initialization

**File:** `direct/lib/sundman_seed_map.m`, lines 16–17, 28, 34–35  
**Quote:** `alpha = w ./ max(sqrt(sum(w.^2,1)), 1e-9);`

For an admitted cone seed, coasting nodes can have \(w=0\). The result is then \(\alpha=0\), not a unit vector. At that initial value, the gradient of
\[
\alpha^\top\alpha-1=0
\]
is zero. This is an initialization defect, unlike the perfectly regular unit-sphere constraint at a feasible point.

**Fix:**
```matlab
nw = vecnorm(w,2,1);
good = nw > 1e-12;

alpha = zeros(size(w));
alpha(:,good) = w(:,good) ./ nw(good);
alpha(:,~good) = repmat([1;0;0], 1, nnz(~good));
```
Using neighboring valid directions is a better warm start, but any genuinely unit fallback satisfies the representation. Validate consistency between positive throttle and nonzero cone direction.

---

### E8. The min-time driver can publish a failed iterate and claims more than an all-burn local solve establishes

**Files:**

- `direct/lib/gen_tulip_mintime.m`, lines 70–89, 102–104
- `GTO_ELFO/direct/elfo/casadi_mintime_freetf.m`, lines 8–9, 38

**Quotes:**

- “MATCH: direct solve CERTIFIES the indirect reference.”
- “true direct min-time”
- `save(outFile,...)`
- `.tf(=tfMin)`

The driver prints acceptance diagnostics but does not enforce them before declaring a match or saving. The solver explicitly returns debug iterates after failure.

Furthermore, a local optimum of a **hard all-burn** transcription is not by itself proof of the unrestricted minimum-time problem. One must justify that full thrust satisfies the throttle minimum condition, or compare with an unrestricted-throttle solve. Even then, first-order conditions do not prove global minimum time.

**Fix:** Before declaring a match or publishing an anchor, require successful convergence, full feasibility, unit control, endpoint accuracy, positive monotone time and inactive artificial boxes. For example:
```matlab
assert(strcmp(out.ipoptStatus,'Solve_Succeeded') && ...
       out.maxDefect < 1e-8 && out.maxUnit < 1e-8 && ...
       rferr < 1e-8 && out.tMonotone && ~out.boundSat.hit, ...
       'gen_tulip_mintime:uncertified', ...
       'Refusing to publish an uncertified min-time candidate.');
```
Then add the appropriate first-order and mesh checks.

Use “converged all-burn local min-time candidate” and “agrees with the reference” unless the stronger claims have separate supporting evidence.

---

### E9. The moved min-time driver contains a wrong solver path and a conflicting results root

**File:** `direct/lib/gen_tulip_mintime.m`, lines 33–35  
**Quote:**
```matlab
addpath(fullfile(here, '..', 'elfo'));
resDir = fullfile(here,'results');
```

With this file in `GTO_tulip/direct/lib`, the added directory is `GTO_tulip/direct/elfo`, not the supplied solver location `GTO_ELFO/direct/elfo`.

The result directory is `GTO_tulip/direct/lib/results`, whereas `minfuel_config` establishes `GTO_tulip/direct/results` as the campaign tree. A preconfigured MATLAB path can hide the former error.

**Fix:**
```matlab
repoRoot = fileparts(fileparts(fileparts(here)));
addpath(fullfile(repoRoot,'GTO_ELFO','direct','elfo'));

cfg = minfuel_config();
resDir = cfg.dirs.root;
```
Resolve setup paths explicitly rather than relying on `cd(here)`. If `lib/results` is intentionally a separate tree, document it and align the generating/consuming drivers accordingly.

## Imprecisions

### I1. The \(10^{-12}\) guards soften gravity; they are not exact CR3BP

**Files:**

- `direct/lib/casadi_minfuel_sundman.m`, lines 110–111
- `GTO_ELFO/direct/elfo/casadi_mintime_freetf.m`, lines 92–94
- `direct/certify/certify_minfuel_pmp.m`, lines 102–106
- `direct/lib/sundman_seed_map.m`, lines 37–46

**Quote:** “Earth distance (guarded)”

Adding \(10^{-12}\) to squared distance gives a softening length \(10^{-6}l_\star\), approximately **0.390 km**. It changes both gravity and the Sundman clock.

At GTO perigee its relative gravity effect is only approximately \(5\times10^{-9}\), so this is unlikely to explain the reported LS failure. Nevertheless, defects of \(10^{-14}\) are defects of the **softened model**, not errors relative to exact CR3BP.

The LS gravity gradient and seed clock are unsoftened. Thus the seed mapper's “exact discrete inverse” is not exactly the inverse of the engine's clock.

**Fix:** Centralize the softening parameter and use the same model everywhere. For a softened gradient:
\[
G=\operatorname{diag}(1,1,0)
-(1-\mu)\left[\frac I{d^3}-\frac{3dd^\top}{d^5}\right]
-\mu\left[\frac I{r^3}-\frac{3rr^\top}{r^5}\right],
\]
where \(d^2=dd^\top+\delta^2\), \(r^2=rr^\top+\delta^2\).

For the seed clock:
```matlab
r1 = sqrt((Xseed(1,:)+muStar).^2 + ...
          Xseed(2,:).^2 + Xseed(3,:).^2 + soft2).';
```
Document the model perturbation, or remove it consistently if exact CR3BP is required.

---

### I2. Interval duals are not arbitrary mesh-scaled costates, and the primer/mass quantities are proxies

**Files:**

- `direct/lib/casadi_minfuel_sundman.m`, lines 61–66, 238–256
- `GTO_ELFO/direct/elfo/casadi_mintime_freetf.m`, lines 187–198

**Quotes:**

- “up to a positive mesh-weight scaling and a global sign”
- `lamMassEnd = lamDef(7,end)`
- “independent optimality certificate”

The extraction order is correct: `D(:)==0` is the first constraint call, and MATLAB/CasADi column-major stacking gives consecutive interval blocks.

For these **unscaled state-increment defects**, however, the interval dual already has costate scaling. It should not generally be divided by interval length. The node control-stationarity combination is, at an interior node,
\[
\bar\mu_k
=\frac{h_{k-1}\mu_{k-1}+h_k\mu_k}{h_{k-1}+h_k},
\qquad p_k=-\bar\mu_k.
\]
Endpoint combinations are one-sided.

`lamDef(7,end)` belongs to the **last interval**, \([\sigma_N,\sigma_{N+1}]\). It is not the exact terminal covector; an \(O(h)\) endpoint offset is unsurprising. Comparing interval duals to averaged endpoint controls is likewise only a midpoint approximation.

**Fix:** State the actual multiplier convention and node map. Rename outputs, for example, `lamMassLastInterval` and `primerMidMeanDeg`. Use the signed nodal minimum-condition check from E2 and the terminal map from E3 for certification. Do not call a diagnostic derived from the same NLP duals independent evidence.

---

### I3. Neighbor time rescaling is a legitimate infeasible guess, not a Sundman transformation

**File:** `direct/lib/minfuel_at_tf.m`, lines 77–80  
**Quote:** “rescale time state to new t_f”

With positions and `tauf0` unchanged, the right-hand side of the time defect is unchanged. Multiplying the time state by \(q=t_f/t_{f,\rm old}\) therefore creates time defects proportional to \(q-1\).

That is acceptable for a continuation **initial guess**, but it does not preserve Sundman feasibility or represent a dynamically consistent time dilation. Rescaling `tauf0` alone would fix the time defects while disturbing the other state defects; it is not a general cure.

**Fix:** Document explicitly:
```matlab
% Endpoint-consistent, dynamically infeasible continuation guess.
% Only the carried-time state is rescaled; the next solve repairs defects.
```
For the formulation corrected under E1, initialize and continue the scale state deliberately.

The output formulas at lines 131–132 are correct:
- `dV` is ideal accumulated thrust \(\Delta V\) in **km/s**, not the norm of net inertial or rotating-frame velocity change.
- `prop_kg` is consumed propellant for initial normalized mass one.

---

### I4. The switching regularity statistic is mesh-deweighted, but it is not physical \(\dot S\)

**File:** `verify_common/foc_check.m`, lines 223–245, 249–278  
**Quote:** `Shat = rep.Sd ./ w;`

Removing the trapezoidal **sigma** weights is a real improvement. For the single-primary solver, however,
\[
\frac{Sd_k}{w_k}\approx tauf\,\kappa_k S_{\rm phys,k},
\]
not \(S_{\rm phys,k}\). With a scale state, include `cScale` as well.

A positive smooth factor preserves zero crossings and simple-versus-multiple zeros in the continuum. But finite-mesh slope magnitudes and coast-normalized thresholds can be strongly influenced by the clock, particularly around perigee. The near-zero three-node test also does not represent a mesh-independent physical arc duration.

**Fix:** Either rename these outputs as **sigma-Hamiltonian regularity diagnostics**, or use
```matlab
Sphys = Sd ./ (w .* tauf .* kappa);  % also .* cScale where applicable
Sdot  = diff(Sphys) ./ diff(X(timeRow,:));
```
Normalize with physical transfer duration and a documented physical-switching-function scale. Require a refinement comparison before claiming regularity.

Also replace line 140's exact-sign/burn-threshold classification with tolerance-aware tests: negative/positive switching coefficient on upper/lower bounds, and near-zero coefficient for interior throttle. Interior throttle at an isolated discrete switching node is not automatically a sign-law failure.

---

### I5. The clock/Hessian and strict-convexity claims overstate what the formulas establish

**File:** `direct/lib/casadi_minfuel_sundman.m`, lines 8–10, 22, 185, 209–214  
**Quotes:**

- “gravity Hessian terms (~1/r^3)”
- “energy, strictly convex”
- “At a genuine local min … reduced Hessian is PD”

The \(r^{-3}\) scaling describes a gravity **Jacobian**; second derivatives of acceleration scale as \(r^{-4}\). Multiplying by a state-dependent clock also generates derivative-of-clock terms. The stated Hessian exponent is therefore not a correct general Hessian estimate.

At \(\epsilon=1\), the integrand is strictly convex in scalar throttle for fixed positive clock. The complete NLP is not convex, and smooth ramps or absence of restoration are not guaranteed.

Finally, zero IPOPT Hessian regularization is not a sufficient second-order certificate. Genuine non-strict minima can have semidefinite reduced Hessians; coast-node direction variables supply obvious flat directions here.

**Fix:** Replace with qualified statements about improved conditioning and scalar-throttle smoothing. Describe `regHistory` as an inertia/solver diagnostic, not a native proof of SOSC. An SOSC claim requires analysis on the relevant critical cone/reduced space, including flat control directions.

---

### I6. The seed time parameter and endpoint selection need narrower descriptions

**Files:**

- `direct/lib/sundman_seed_map.m`, lines 19, 31–32, 40–53
- `cr3bp_common/gto_tulip_endpoints.m`, lines 6, 25–28
- `direct/run_certified_minfuel.m`, lines 54–62

**Quotes:**

- “seed node parameter … any monotone range”
- `tSeed = sgNorm * tf`
- “arrival at the maximum-ydot point”
- “endpoint check vs pumpkyn”

`sgNorm` must be affine in physical seed time. An arbitrary monotone trajectory parameter is insufficient. The reciprocal-of-mean clock formula is otherwise the correct discrete inverse.

Endpoint pinning occurs **after** computing the clock. Changing an endpoint position can invalidate the initially exact time trapezoid even after softening is made consistent.

The tulip maximum is a maximum over returned propagation samples, not an event-refined continuous maximum. Moreover, `insertion_states('tulip','campaign')` and `'maxydot'` are deliberately distinct endpoint choices.

The flagship endpoint “check” only prints differences and catches all exceptions; it is not a drift gate.

**Fix:** Specify physical-time mesh coordinates, pin positions before computing the clock, and say “sampled maximum-\(\dot y\) point”. If endpoint identity is required, assert against the declared `campaign` endpoint; keep optional toolbox regeneration as a separate diagnostic.

The external `orb2eci`, `fromPCI`, `getTulip` and `prop` implementations were not supplied. Their argument conventions, frame conversion and orbit-family parameter consistency cannot be independently verified from these files alone.

---

### I7. Some output labels imply stronger certification or freer bounds than the code provides

**Files:**

- `direct/lib/minfuel_at_tf.m`, lines 39–42, 117–119
- `direct/lib/casadi_minfuel_sundman.m`, lines 259–272
- `GTO_ELFO/direct/elfo/casadi_mintime_freetf.m`, lines 41–45, 206–225

**Quotes:**

- “certified … at least one schedule step converged tight”
- “BCs pin the endpoints by construction”
- “widen via opts”

The driver correctly requires convergence at the requested schedule endpoint, not merely any step. A custom schedule ending above zero can nevertheless produce `certified=true`; that is schedule completion, not pure-fuel certification.

The endpoint-bound exclusion rationale is too broad: final mass is free, and in the min-time problem final time and scale are free. The diagnostics omit terminal mass and the time cap.

The min-fuel solver exposes only position/velocity box options; the min-time solver exposes neither position nor velocity nor mass box options. “Widen via opts” is therefore not generally actionable.

**Fix:** Separate `scheduleConverged`, `epsReached` and `fuelFirstOrderVerified`. Check artificial bounds on every **free** component, including terminal components. Expose actual bound options or give the correct source-edit instruction.

## Readability

### R1. The min-time driver's header describes a different seed and clock

**File:** `direct/lib/gen_tulip_mintime.m`, lines 5–17, 23 versus 44–66  
**Quote:** “SINGLE-primary clock (moonZone=0)”; “energy_f1120.mat”; `moonZone(=0)`

The body defaults to a tagged **two-primary** energy seed and `moonZone=0.15`.

**Fix:** Rewrite the purpose, warm-start, options and saved-field descriptions to match the body. Document `.moonZone` as an input option. Replace the unsupported “28k km … overflows … Hessian” causal assertion at lines 60–62 with a measured solver-conditioning/factorization observation unless that mechanism has actually been isolated.

---

### R2. Several cross-references and default locations are stale

**Files, lines and fixes:**

- **`cr3bp_common/cr3bp_ipopt_opts.m`, line 6:**  
  Quote: `GTO_tulip/direct/sundman_minfuel/casadi_minfuel_sundman.m`  
  Replace with `GTO_tulip/direct/lib/casadi_minfuel_sundman.m`.

- **`direct/lib/casadi_minfuel_sundman.m`, lines 42–44:**  
  Quote: “adaptive barrier”  
  Replace with “monotone barrier, larger `mu_init`, default warm-start pushes”, matching `cr3bp_ipopt_opts`.

- **Same file, line 25:**  
  Quote: `run_sundman_*`  
  Reference the actual homotopy/per-\(t_f\) drivers.

- **Same file, lines 85–86:**  
  Quote: this file's own pathname as the pattern it “mirrors”.  
  Replace with the intended `PSR/lib/...` precedent, consistent with lines 212–213, or remove the self-reference.

- **`direct/certify/certify_minfuel_pmp.m`, line 46; `run_foc_tulip.m`, lines 20–21:**  
  Quote: “beside this file/one”  
  Replace with `../lib/sundman_minfuel_certified.mat`.

- **`direct/certify/run_foc_tulip.m`, line 53:**  
  Quote: “this folder's own sundman_homotopy.m”  
  Correct to the library location.

The certifier's `../../min_energy_tutorial/primer_check.m` reference also needs validation against the current tree; the referenced file was not supplied.

---

### R3. Several functions do not meet the required header contract

**Files and lines:**

- `direct/lib/sundman_seed_map.m`, lines 2–29: no **REFERENCES** section.
- `direct/run_certified_minfuel.m`, lines 2–34: no **REFERENCES** section; `best` has no output dimensions or field-size reference.
- `cr3bp_common/minfuel_config.m`, lines 86–89: `local_fparse` has only a one-line purpose.
- `direct/lib/minfuel_at_tf.m`, lines 158, 192, 204, 228: all four local helpers lack full INPUTS/OUTPUTS/REFERENCES headers.
- `direct/certify/certify_minfuel_pmp.m`, lines 273–275: `ternary` has no header.
- `direct/certify/run_foc_tulip.m`, lines 136–138: incomplete helper header.
- `verify_common/foc_check.m`, lines 327–328: no helper header.
- `direct/lib/gen_tulip_mintime.m`, lines 110–111, and `casadi_mintime_freetf.m`, lines 244–245: no helper header.

**Fix:** Add the house-style sections to every function, including local helpers. For trivial utilities, `REFERENCES: None` is sufficient. Mark parameter/report structs as scalar structs and specify sizes for their important array fields.

---

### R4. Prohibited loop indices and inconsistent assignment alignment remain

**File:** `direct/certify/certify_minfuel_pmp.m`, lines 215–216, 263–264  
**Quote:** `for i = 1:numel(swIdx)`

**Fix:**
```matlab
for kSwitch = 1:numel(swIdx)
    ...
end
```
Update both indexed uses. These are the only `i`/`j` loop-variable violations in the supplied code; `iiR`, `jjR`, `gi` and `gj` are not violations.

**File:** `cr3bp_common/cr3bp_lt_params.m`, lines 30–39  
**Quote:** `p.muStar = ...` versus `p.thrustN = ...`

Align the parameter assignments consistently. The dense multi-assignment lines in `minfuel_at_tf.m`, particularly lines 71–80, would also benefit from one aligned assignment per line.

---

### R5. The LS header still describes superseded algorithms and an unconditional REVIEW policy

**File:** `direct/certify/certify_minfuel_pmp.m`, lines 4–5, 118–130 versus 149–183, 241–252  
**Quotes:**

- “returns a REVIEW verdict, not PASS”
- “Solve min_y ||P y||²”
- “immune to that dynamic range”

The body uses signed, bounded `lsqlin` fitting with auxiliary magnitudes, not the projector-only saddle-point objective described immediately above it. It can print PASS.

**Fix:** Describe the actual signed constrained-LS problem and distinguish historical results from enforced policy. If this is intentionally a PARTIAL diagnostic, make that status explicit in the returned struct and printed verdict, rather than relying on the flagship currently failing its numerical thresholds.

## Checklist verdicts

- **A — I1:** Centrifugal term, Coriolis term, Earth/Moon vectors and gravity signs **PASS**. The guards define softened CR3BP and are inconsistent with the unsoftened LS/seed formulas.
- **B — E1, I1, I6:** State and objective clock factors, time ODE, mass ODE and sigma trapezoid **PASS**; fixed `tauf0` plus fixed physical `tf` is not the advertised equivalent formulation.
- **C — PASS:** `IntF - epsilon*IntS` has the correct sign and scale; at one it is \(\int s^2dt\), at zero it is proportional to propellant. `Tmax/m`, `Tmax/c`, `c` and the unit conversions are correct.
- **D — PASS:** At \(\|\alpha\|=1\), the equality gradient is nonzero and the \(\pm1.1\) box is inactive. Coast directions create flat directions, not this alleged sphere/box constraint-qualification defect. See E7 for zero-direction initialization.
- **E — E2, E3, I2:** First-`8N` column-major extraction **PASS**; signed primer certification and terminal-covector interpretation require correction.
- **F — I3:** `dV` and `prop_kg` **PASS**; neighbor time scaling is only an intentionally infeasible initial guess.
- **G — E4, R5:** The physical adjoint matrix/recursion is recognizable and correctly signed, but the fitted problem, normalization and diagnostic claims need correction.
- **H — E2, E3, E5, I4:** Full AD stationarity assembly is correct for the supplied models; the PASS criteria and physical interpretations are insufficient. The unprovided manifest and dual-mapping helper remain unverified.
- **I — E5, E8:** The `cScale` dynamics and Mayer objective **PASS**. Raw time-dual \(-1\) is correct under the stated defect convention, but is not alone \(H(t_f)=0\); a 5% variation is not intrinsic to this formulation.
- **J — R1–R5:** House-style, stale-header and cross-reference defects remain.

**Clean portions/files:** `cr3bp_lt_params.m` is mathematically and dimensionally clean; only alignment/header wording needs attention. `minfuel_config.m` is clean within its documented nominal-anchor and post-default-override contract, apart from the local-helper header. `cr3bp_ipopt_opts.m` is functionally clean; its stale solver pathname needs updating. No rotating-frame force-sign error was found in either engine.

## Verdict

**The chain is not yet mathematically certified for the problem as stated.**

The principal issue is not a force sign, thrust unit or quadrature factor: it is the extra trajectory constraint introduced by fixing the regularized duration. The existing excellent residuals can be entirely genuine while certifying that restricted problem.

The next priorities are to distinguish signed Hamiltonian minimization from tangential KKT stationarity, compute terminal covectors without bound cancellation, and resolve the LS model mismatch before explaining its failure as accumulated discretization error.

The sparse `cScale` machinery already present provides a practical route to correcting the main formulation issue without reintroducing a dense free-duration column.
---

## Run metadata
GPT-6 Astra via the raw OpenRouter API (`astra_raw.sh`; crush cannot carry Astra on bundles this size): HTTP 200, wall 264 s, 42,335 prompt / 11,658 completion (3,106 reasoning), finish_reason=stop, cost $1.112. Bundle 127 KB = prompt (`direct_core_chain_astra_prompt_2026-09-06.txt`) + 14 files, line-numbered.

## Host adjudication (2026-09-06, Claude; nothing applied yet)

| item | host verdict | evidence |
|---|---|---|
| **E1** fixed τ_f + fixed t_f adds the isoperimetric constraint ∫dt/κ = τ_f0 | **CORRECT in principle, and the campaign already has evidence it binds.** TODO P0 (2026-07-26): freeing τ_f via `cScale` at the flagship t_f, ε=1, converged with **cScale = 1.0051** — if the constraint were slack, cScale would have stayed at 1. Magnitude of the effect on m_f at ε=0: UNMEASURED. First-order checks cannot see this (the control minimisation is unaffected; only the costate ODE gains a −∇κ·(K/κ) term), which is also the most economical explanation of the recorded raw-dual PASS vs LS-reconstruction FAIL on the flagship (E4). | needs the decisive experiment below |
| E2 tangential direction check is sign-blind; engine picks the costate sign from observed steering | correct as coded (`v − (vᵀb)b` in foc_check; `if mean(cang)<0, cang=-cang` in the engine). The register's own I2 already flagged sign-blindness. Fix: resolve the multiplier sign once from the Lagrangian convention; report a signed normalised primer error. | foc_check.m lines 146–155; engine 244 |
| E3 mapped covector reads the full gL entry incl. any active terminal box multiplier | correct in principle; no practical effect on tulip rows (mass box 0.3 ≤ m ≤ 1 is inactive at m_f ≈ 0.85), but the logic is wrong in general and the "cannot fail unless kktStat has" comment is false for the relative test. | foc_check.m 181–217 |
| E4 LS certifier fits the *physical* adjoint, not the fixed-τ_f NLP's | follows from E1; the disagreement on the flagship is then EXPECTED, not "40 revs of truncation". Also real: rho=1 soft normalisation, mass reconstruction assumes the primer fit, W(swIdx) at the left node. | certify_minfuel_pmp.m |
| E5 λ_t = −1 is the defect-multiplier condition, not H_phys(t_f)=0; the 5% lamTimeCoV is not intrinsic to cScale | plausible; the mapping helper (`foc_dual_to_costate`) was not in the bundle. OPEN. | LEAD-4 in the register |
| E6 'apoapsis' = min rotating speed | correct as coded; only the `'campaign'` criterion is used in production. Rename or recode. | insertion_states.m 67–71 |
| E7 zero unit direction at coast nodes in the seed map | correct as coded (`w./max(‖w‖,1e-9)` → 0). Initialisation wart; fix is a two-line fallback. | sundman_seed_map.m 28 |
| E8 min-time driver prints but does not enforce acceptance | correct; add the assert. "true direct min-time" → "converged all-burn candidate". | gen_tulip_mintime.m |
| **E9** `addpath(here/../elfo)` from `direct/lib` → `GTO_tulip/direct/elfo`, which does not exist; `resDir = lib/results` also nonexistent | **CONFIRMED** (`ls direct/elfo` fails). Masked by setup_paths. Second instance of the 2026-07-26 flatten breaking a literal path that `test_artifact_paths` doesn't cover (it checks .mat paths, not addpath targets). | verified 2026-09-06 |
| I1 1e-12 softening (≈0.39 km) differs between engine, LS certifier and seed clock | correct; effect ~5e-9 at perigee, but the inconsistency should be centralised. | |
| I2 interval duals already carry costate scaling; lamMassEnd is the last-interval dual | correct; matches the register's I5 (superseded by the mapped covector). | |
| I3 neighbor time rescale is an infeasible continuation guess | correct; document as such. | |
| I4 Ṡ statistic is σ-Hamiltonian, not physical Ṡ (misses τ_f·κ) | correct; the register's I1 already suspected a discretisation confound. | |
| I5–I7, R1–R5 | accept; R2's stale paths are real (`sundman_minfuel/` references survive the flatten). | |

**The decisive experiment for E1** (not run): solve the flagship (25 mN, factor 1.150, ε: 1 → 0, seeded from `sundman_minfuel_certified.mat`) with the free-time solver `casadi_energy_freetf` (moonZone ≤ 0, i.e. the same single-primary clock) so that τ_f is free through `cScale` while t(τ_f) = t_f stays pinned, and compare m_f with the fixed-τ_f result in the same basin (0.849066). `run_tulip_ladder` with `thrustStop = thrustStart = 0.025`, `epsMin = 0` is essentially this run. If m_f moves by more than the basin-sweep resolution (~1e-8 reproducibility) the certified rows are extremals of the restricted problem and every published m_f needs re-solving free-τ_f; if it does not, E1 is a formal defect with no practical consequence at this t_f — either way the docs must stop calling the fixed-τ_f trick "equivalent".
