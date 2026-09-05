## A. PMP / variational equations

- **[CORRECTNESS]** `costate_common/cr3bp_minfuel_pmp.m:10-14,118-121` — The implemented “eps” family is **not** the stated tanh/logistic Bertrand–Epenoy smoothing. It uses the polynomial homotopy
  \[
  L_\epsilon(s)=(1-\epsilon)s+\epsilon s^2,\qquad
  s^*=\operatorname{clip}\!\left(\frac{Q-(1-\epsilon)}{2\epsilon},0,1\right).
  \]
  It is internally consistent with that running cost, but it is a different regularization with exact coast arcs and hard clipping. Either rename the experiment and claims, or implement the advertised law
  \[
  s^*=\tfrac12\!\left[1+\tanh\!\left(\frac{Q-1}{\epsilon}\right)\right]
  \]
  together with the matching entropy term
  \[
  L_\epsilon(s)=s+\frac{\epsilon}{2}\big[s\log s+(1-s)\log(1-s)\big].
  \]

- **[ROBUSTNESS]** `costate_common/cr3bp_minfuel_pmp.m:114-116,133-141` and `costate_common/cr3bp_minenergy_pmp.m:87-100` — The primer regularization is not Hamiltonian-consistent. With \(n=\sqrt{\lambda_v^T\lambda_v+10^{-300}}\) and \(\alpha=-\lambda_v/n\), the implemented Hamiltonian contains
  \[
  \lambda_v^T\frac{sT\alpha}{m}
  =-\frac{sT}{m}\frac{\lambda_v^T\lambda_v}{n},
  \]
  but `Q` assumes the coefficient is \(-sTn/m\). At \(\lambda_v=0\), the code can consume propellant while producing zero thrust acceleration. Reject `norm(lv)<primerTol` as singular, or consistently use `primer=dot(lv,lv)/nlv` in `Q` and the throttle law.

- **[CORRECTNESS]** `costate_common/cr3bp_minfuel_pmp.m:123-128,143` and `costate_common/cr3bp_minfuel_prop.m:15-17,84-86` — The Huber vector field jumps at the state-dependent surface \(Q=1\), but the STM integrates only \(\dot\Phi=A\Phi\). It omits the saltation matrix, so the shooting Jacobian is wrong whenever a segment crosses the switch:
  \[
  \Phi^+=\left[I+
    \frac{(F^+-F^-)\nabla(Q-1)^T}
         {\nabla(Q-1)^TF^-}
  \right]\Phi^-.
  \]
  Add event localization, split propagation at each crossing, and apply this update. This omission can itself produce the reported Newton residual floor, so the Huber “structural knockout” is not established.

- **[ROBUSTNESS]** `costate_common/cr3bp_minfuel_pmp.m:95-96,126-128` — Huber accepts any `p>0`; for `p>1` and \(Q<1\), `s=min(max(pQ,0),p)` can exceed one. Enforce `0 < p <= 1`, or cap the result at one.

- **[ROBUSTNESS]** `costate_common/cr3bp_minfuel_pmp.m:119-120,143` and `costate_common/cr3bp_minenergy_pmp.m:89-102` — At `sRaw=0` or `sRaw=1`, hard clipping makes the field only piecewise differentiable. CasADi supplies a branch derivative, not the claimed “exact” classical Jacobian. Detect clipping-surface crossings and use piecewise/saltation treatment, or replace clipping with a differentiable saturation.

- **[ROBUSTNESS]** `costate_common/cr3bp_minfuel_prop.m:65-75` and `costate_common/cr3bp_minenergy_prop.m:59-69` — Neither propagator verifies that `ode113` reached `dt`; MATLAB can return a partial trajectory after a warning, which is then accepted as the endpoint. Add `T(end)==dt` within tolerance, finite-output checks, and terminal events for \(m\le m_{\min}\) and primary collision.

The unregularized expressions away from switches are otherwise consistent:
\[
Q=T\left(\frac{\|\lambda_v\|}{m}+\frac{\lambda_m}{c}\right),\qquad
\dot\lambda_m=-\frac{sT\|\lambda_v\|}{m^2}.
\]

## B. Conjugate-point test

- **[CORRECTNESS]** `DRO_tulip/indirect/ms_minenergy.m:87-95` and `DRO_tulip/indirect/ms_minfuel.m:67-72` — The fixed-time test monitors rows `[1:6,14]`, i.e. the **terminal shooting residual Jacobian** \(\partial[r_f,v_f,\lambda_m(t_f)]/\partial\lambda_0\). That is not the interior fixed-initial-state Jacobi projection. Interior conjugacy requires rank loss of
  \[
  \frac{\partial x(t)}{\partial\lambda(0)}
  =\Phi(1\!:\!7,8\!:\!14).
  \]
  Use `stateRows=1:7`, `costateCols=8:14`, `quotientDir=[]`, `freeTime=false`. Test `[1:6,14]` separately only as terminal BVP-Jacobian conditioning.

- **[CORRECTNESS]** `costate_common/ms_conjugate_test.m:113-124` — The final segment STM is never chained: the loop ends at `K-1`. Therefore every conjugate point in the last segment is invisible. Chain all `K` STMs and sample `info.tGrid(2:K+1)`; separately classify an actual rank loss exactly at \(t_f\).

- **[CORRECTNESS]** `costate_common/ms_conjugate_test.m:129-136` — `atFinal` is assigned when a sign-changing interval ends at the **last sampled interior junction**, not at \(t_f\). That crossing is strictly interior, yet it is subtracted from `nCrossings`, potentially turning a refutation into `pass=true`. Count every bracket wholly inside \((0,t_f)\); set `atFinal` only after endpoint root localization or an endpoint rank test.

- **[CORRECTNESS]** `costate_common/ms_conjugate_test.m:122-133` — Determinant sign changes miss even-multiplicity rank losses and any pair of crossings inside one segment. Exact `dets==0` is also not a floating-point rank criterion. Equilibrate the block and monitor \(\sigma_{\min}(M)/\sigma_{\max}(M)\) with an explicit tolerance and adaptive interval refinement.

- **[CORRECTNESS]** `costate_common/ms_conjugate_test.m:129-133` — A node-aligned root can be double-counted. For signs `(+ ,0, -)`, the compressed nonzero samples contribute one sign flip and `nnz(dets==0)` contributes another. Merge contiguous zero samples and adjacent sign-change brackets into one localized root.

- **[ROBUSTNESS]** `DRO_tulip/indirect/ms_minenergy.m:84-95` and `DRO_tulip/indirect/ms_minfuel.m:64-72` — The conjugate test runs even when `info.converged` is false. Refuse to issue a verdict unless the BVP residual and propagation checks pass.

- **[CORRECTNESS]** `costate_common/tests/test_conj_fixedtf.m:35-56` — The scalar LQ test cannot distinguish the wrong terminal block from the correct state-projection block, and it tests neither endpoint roots, node-aligned roots, even-multiplicity losses, nor two crossings within one segment. Add multidimensional fixtures that distinguish `Phi(1:7,8:14)` from `[Phi(1:6,8:14);Phi(14,8:14)]`, plus those missing cases.

## C. Continuation and acceptance

- **[CORRECTNESS]** `DRO_tulip/run_minfuel_grid.m:47-58` — Any arm with at least one converged rung is packaged, even if it retired far above the requested final epsilon or failed `acceptOk`. Require an explicit final verdict such as `finalReached && acceptOk && msResidualPassed && HdriftPassed` before appending to `G`.

- **[CORRECTNESS]** `DRO_tulip/run_minfuel_race.m:144-162` — Failure of the deepest single-shooting gate does not invalidate the arm; the solution remains exported with its mass and junction states. If that gate is genuinely required, set `A.valid=false` and make grid/catalog packaging reject it. Otherwise stop calling it an acceptance gate.

- **[ROBUSTNESS]** `DRO_tulip/run_minfuel_race.m:129-139` — `maxBisect` is documented as a per-gap limit, but the implementation allows `maxBisect*numel(P.sched)` total insertions—51 under defaults. Track a per-gap counter and reset it only after success or abandonment.

- **[CORRECTNESS]** `DRO_tulip/run_minfuel_grid.m:34-35` — A MATLAB string verdict such as `"FAIL"` enters the non-`char` branch and can be accepted if stale `gates` are all true. Normalize with `strcmpi(string(r.verdict),"PASS")`; legacy gates should be migrated explicitly rather than overriding a recorded verdict.

- **[ROBUSTNESS]** `DRO_tulip/run_minfuel_grid.m:55-57` and `DRO_tulip/build_minfuel_catalog.m:74-78` — `G.tGrid` is always discarded and the catalog reconstructs the last segment with `rec.tf/K`. This is valid only for a uniform grid. Store the actual final `tGrid` and use `g.tGrid(end)-g.tGrid(end-1)`.

- **[ROBUSTNESS]** `DRO_tulip/run_minfuel_race.m:86-97,142` — “Save after every step” is not restartability: the saved queue, current seed and counters are never loaded. An interruption reruns the arm from the energy seed and may enter a different basin. Persist and validate the continuation state, then resume it.

## D. Catalog integrity

- **[CORRECTNESS]** `DRO_tulip/build_minfuel_catalog.m:29-42` — The gamma rescue is an unauthenticated MAT-file substitution. The builder checks only `conjPass`; it does not verify final epsilon, BVP residual, endpoints, acceptance, dimensions, or that `kBad` exists uniquely. Validate all of these and bind the artifact to the solution using a fingerprint of `z`, `Y`, `tf`, smoothing family and `pDeepest`.

- **[CORRECTNESS]** `DRO_tulip/build_minfuel_catalog.m:59,80-89` — Missing conjugate verdicts leave `conj_pass=-1` on solved entries, yet packaging succeeds. Lookup also selects the first match by cell/gamma without proving it belongs to the packaged branch or epsilon. Require exactly one branch-bound verdict in `{0,1}` for every solved entry and abort otherwise.

- **[CORRECTNESS]** `costate_common/catalog_schema.m:142-154` — The validator describes `conj_pass` as an `int8` grid but checks only shape and provenance. It accepts doubles, NaNs, arbitrary values, and unknown `-1` on solved entries. Enforce the type, allowed values, and `conj_pass(has_solution) ~= -1`.

- **[CORRECTNESS]** `costate_common/catalog_schema.m:118-140` — Core v3 invariants are unvalidated: `lam0` need not have seven rows; `tf_nd`, `tfmin_nd`, and `entry_index` shapes are unchecked; entry indices can duplicate or omit columns; `axis3.name` need not be `"gamma"`; and `tf_nd=gamma*tfmin_nd` is not enforced. Add these checks, including that solved indices are exactly a permutation of `1:nSolved` and unsolved indices are zero.

- **[CORRECTNESS]** `costate_common/catalog_schema.m:184-192` — Both delta-v formulas use matrix division (`1/mf`) rather than elementwise division. Grid inputs can error or produce matrix-algebra results. Use
  ```matlab
  cnd .* log(1 ./ mf) .* lStar ./ tStar
  ```
  and validate `0 < mf <= 1`. The scalar dimensional conversion to km/s is otherwise correct.

- **[ROBUSTNESS]** `costate_common/catalog_schema.m:73-80` — Unknown future versions such as schema 4 or 99 are silently validated as v3. Require an integer version in `{1,2,3}` and reject unsupported versions.

- **[ROBUSTNESS]** `costate_common/catalog_schema.m:101-104,146-148,160-163` — Validation can throw after detecting missing sheet fields because later checks still dereference `arr_family` or `has_solution`. Guard dependent checks so malformed catalogs return diagnostics rather than crashing.

- **[ROBUSTNESS]** `DRO_tulip/build_minfuel_catalog.m:48,60,78` — The junction count is taken from `G(1)` and imposed on every entry. Mixed segment counts will break packaging. Either declare and validate a common grid or store per-entry junction arrays in cells.

## E. Code quality and tests

- **[EFFICIENCY]** `DRO_tulip/run_minfuel_race.m:110-113,132` and `DRO_tulip/build_minfuel_catalog.m:73` — Result vectors, queues, and `lam0` grow dynamically. Preallocate to the schedule/bisection upper bound and trim, or use indexed record arrays.

- **[STYLE]** `costate_common/cr3bp_minfuel_pmp.m:2-79`, `costate_common/cr3bp_minenergy_pmp.m:2-58`, `costate_common/cr3bp_minfuel_prop.m:2-44`, `costate_common/cr3bp_minenergy_prop.m:2-39`, `costate_common/ms_bvp.m:2-94`, `costate_common/ms_conjugate_test.m:2-81`, and `costate_common/catalog_schema.m:2-62` — These mathematical top-level functions omit the required `REFERENCES` section. Add authoritative CR3BP/PMP, regularization, multiple-shooting, conjugate-point, and rocket-equation references.

## Claims not established

- **Claim 18 is not established.** `cr3bp_minfuel_pmp.m:123-143` and `cr3bp_minfuel_prop.m:84-86` omit saltation across the Huber switch. The observed Newton residual floor can therefore be an incorrect-Jacobian artifact, not proof that the family is structurally unusable. The “eps” arm is also not the advertised tanh/logistic regularization.

- **Claim 19 is not established as a certified record set.** `run_minfuel_grid.m:47-58` publishes incomplete or acceptance-failed walks, and `run_minfuel_race.m:144-162` does not make final acceptance binding. The reported numerical entries may be real, but the code does not enforce the stated certification.

- **Claim 20 is not established.** `ms_minenergy.m:87-95`, `ms_minfuel.m:67-72`, and `ms_conjugate_test.m:113-136` use the wrong interior block, omit the final segment, suppress a last-interior crossing, and rely on determinant sign changes. Consequently neither the 13 passes nor the single refutation is reliable.

- **Claim 21 is not established.** `build_minfuel_catalog.m:29-42,80-89` trusts an external rescue artifact and can ship missing or stale conjugate verdicts; `catalog_schema.m:118-154` validates neither branch identity nor mandatory solved-entry verdicts. No source-controlled gamma-continuation driver was present.

## Questions / could not verify

- The referenced MAT artifacts and direct-collocation/IPOPT code were not included, so the reported masses, residuals, monotonicity, dual-sign convention, and endpoint data could not be independently reproduced.
- No gamma-continuation rewalk source file was present; only conditional consumption of its MAT artifact appears in `DRO_tulip/build_minfuel_catalog.m:29-42`.
