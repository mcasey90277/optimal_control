# GPT-6 Astra review of the cart-pole PMP-BVP demo (2026-09-17)

Raw OpenRouter API (crush wedges on Astra above ~15 KB), 70 KB bundle of all 14
files, reasoning effort high: HTTP 200, 398 s, 23,494 prompt / 13,847 completion
tokens (7,379 reasoning), $0.99.

## Host verification of the headline finding (done before anything else)

Astra's first finding -- the pendulum row has the wrong sign -- is CONFIRMED,
and it is upstream of this demo: `ex2_cart_pole_swing_up/try2/pendulum_accel.m`
lacks the leading minus of Kelly (SIAM Review 2017). Measured with u = 0 from
x0 = [0; 0.7; 0.3; -0.2] over 3 s, using the ORIGINAL ex2 helpers:

| energy convention (coupling, potential) | drift |
|---|---|
| +, + | 6.2 J |
| +, - (hanging-zero) | 64.0 J |
| -, + | 31.9 J |
| -, - | 101.1 J |
| pendulum row NEGATED, hanging-zero energy | **1.6e-12 J** |

No sign convention of the energy is conserved by the committed dynamics, so
this is a genuine error, not a relabelling of which equilibrium is "up".
Released at rest at q2 = 0.05 the committed field moves AWAY from zero
(0.05 -> 0.091 in 0.5 s): "hanging" is unstable. The direct and indirect
solves agree to 0.075% because they solve the SAME wrong problem; the
cross-check validated consistency, never physics.

Not affected: the `ms_bvp` promotion and its gates (the engine is
problem-agnostic), the PMP derivation, the STM assembly and the
nested-complex-step arrangement, all of which Astra judged sound.

---

Static review only: the shooting engine, direct-example helpers, and `.mat` fixtures were not supplied, so I cannot reproduce the reported numbers. There is nevertheless a decisive physics error visible in the supplied equations.

## Principal findings

- **[CORRECTNESS]** `cartpole_field.m:58–60` — **These are not the dynamics of the specified hanging-zero cart-pole.** At zero force and near the alleged hanging equilibrium, the code gives
  \[
  \ddot q_2 \simeq \frac{(m_1+m_2)g}{Lm_1}q_2=+5.88q_2,
  \]
  making hanging unstable. Near \(q_2=\pi\), it makes the angular mode stable—the opposite of the specified physics.

  More specifically, take the conventional bob coordinates
  \[
  r_b=(q_1+L\sin q_2,\,-L\cos q_2).
  \]
  The mechanical equations are
  \[
  (m_1+m_2)\ddot q_1+m_2L\cos q_2\,\ddot q_2
       -m_2L\sin q_2\,\dot q_2^2=u,
  \]
  \[
  L\ddot q_2+\cos q_2\,\ddot q_1+g\sin q_2=0.
  \]
  With \(D=m_1+m_2\sin^2q_2\), the supplied cart acceleration is correct for this convention, but the angular acceleration must be
  \[
  \ddot q_2=
  -\frac{L m_2\cos q_2\sin q_2\,\dot q_2^2
       +u\cos q_2+(m_1+m_2)g\sin q_2}{LD}.
  \]
  **Both `F(4)` and `G(4)` need their signs reversed.** Merely relabelling zero as upright does not repair all the supplied centrifugal and coupling signs.

  The same error is reproduced in `gen_direct_ref.m:101–102` and `gen_state_jac.m:33–34`. Matching the inherited helpers therefore establishes compatibility, not physical correctness. Fix the equations consistently, regenerate the Jacobian and both fixtures, and rerun both methods. Add an independent Lagrange-equation or mechanical-power test, not another copy of the acceleration formulas.

- **[CORRECTNESS]** `README.md:89–99`, `tests/test_cartpole_field.m:84–85` — **The documented test commands do not turn failed checks into failed MATLAB processes.** Every test’s `chk` follows the same pattern: print `FAIL`, return `false`, and continue. Both `test_cartpole_field` and `ok = test_cartpole_pmp` can therefore finish normally after failed checks; `matlab -batch` need not return failure.

  Returning a Boolean is a valid local interface, but an asserting runner is missing from the supplied module. Add a runner that executes all six tests and finally calls `assert(all(results))`; document that runner for automation. Alternatively, migrate to `matlab.unittest`. This is distinct from whether the individual conditions are mathematically meaningful.

- **[STYLE]** `cartpole_pmp_rhs.m:12–25`, `README.md:34–43` — **The PMP substitution argument is correct; the trajectory-optimality language needs qualification.** For the normal Hamiltonian,
  \[
  \bar H(x,\lambda)=H(x,\lambda,u^*)=
  \lambda^\mathsf TF-\tfrac14(\lambda^\mathsf TG)^2,
  \]
  and
  \[
  \bar H_x=H_x+H_u u_x^*=H_x
  \]
  because \(H_u=0\). Thus `-A.'*lam`, with `A` differentiated at fixed control and then evaluated at \(u^*\), is exactly right. The outer derivative used for the STM must, and does, differentiate the dependence on \(u^*\).

  There is no endpoint exception: the state endpoints and final time are fixed, so the endpoint variations vanish. There is no condition \(\lambda(t_f)=0\), and **no condition \(H(t_f)=0\)**. The autonomous Hamiltonian is constant, not necessarily zero.

  However, \(H_{uu}=2>0\) establishes the unique **pointwise Hamiltonian minimizer**, not a minimizing trajectory. The nonlinear BVP can have multiple normal extremals, including saddles. Say “normal PMP extremal” rather than implying sufficiency. Also distinguish the four-dimensional initial-costate parametrization from the additional interface-state unknowns introduced by multiple shooting. Being square and smooth does not guarantee a nonsingular shooting Jacobian.

- **[ROBUSTNESS]** `run_cartpole_pmp.m:68–82` — **The sign choice is sensible in principle, but its implementation assumes vector orientation and accepts essentially no confidence criterion.** For real sampled controls, the sign of their dot product selects whichever global multiplier sign minimizes the unweighted squared control mismatch. That is the right scalar replacement for the primer vote.

  However, `uImplied` is explicitly a row while `uDirect` follows `tS`’s interpolation shape. A column `tS` can produce an outer product through implicit expansion at line 77, instead of the intended samplewise product. Moreover, `sum(...) ~= 0` admits arbitrarily weak correlation and can admit non-finite correlation; it does not establish meaningful sign information.

  Force both vectors to columns, check finiteness and nonzero norms, and use a normalized dot product with a documented confidence threshold. For nonuniform stations, use quadrature weights. Check the amplitude mismatch separately: a sign vote cannot detect an incorrect costate scale or station mapping. Not recomputing `uImplied` after the flip does **not** invalidate the existing nonzero check, but storing the corrected implied control would make diagnostics clearer.

  PCHIP interpolation and short endpoint extrapolation are defensible **for a seed**, not as a covector-mapping theorem. Validate station ordering and coverage, report the extrapolation distance, and compare endpoint extrapolation schemes if the seed becomes fragile.

- **[ROBUSTNESS]** `cartpole_pmp_prop.m:95–99` — **Checking that `ode113` reached the requested endpoint is necessary, but the endpoint tolerance is too permissive for its purpose.** A truncated return less than \(10^{-9}\) seconds short is accepted for ordinary segment lengths. For a sufficiently short nonzero `dt`, even negligible progress can satisfy the test. Near a collapse, a tiny time deficit need not imply a tiny state error.

  Use the solver’s successful-endpoint convention—normally the requested endpoint is returned exactly—or a few-ULP comparison, rather than a \(10^{-9}\)-second acceptance band. Validate real, finite, correctly sized inputs before the `dt == 0` shortcut, which currently returns malformed or non-finite inputs without checking them. Test the completion gate independently with finite truncated returns, including one just inside the present tolerance.

- **[ROBUSTNESS]** `cartpole_pmp_prop.m:28–36,87–92` — **`MATLAB:ode*` is not the right semantic boundary between numerical collapse and programming error.** Bad options and derivative-size mismatches can carry ODE-prefixed identifiers, yet backtracking cannot repair them. Conversely, an identifier alone is not a guarantee of where a failure originated.

  More importantly, the comment says the engine uses a **bare catch**. If that is accurate, rethrowing the original identifier here still does not prevent the engine from swallowing the programming error.

  Establish a specific numerical-failure contract: convert only recognized recoverable integration failures to `cartpole_pmp_prop:collapse`, preserve their causes, and make the engine catch only that identifier or an explicitly documented family. Rethrow configuration and programming errors through the engine. Add integration tests proving both behaviors. The engine body is unavailable, so its actual discrimination cannot be verified here.

- **[ROBUSTNESS]** `run_cartpole_pmp.m:100,113–132` — **The reporting integration bypasses the collapse protection and the front door does not enforce its claimed success gates.** It consumes `info.Y` without first checking convergence, calls raw `ode113`, does not verify that the returned time reaches `tf`, and then reports `X(:,end)` as the terminal state. The interpolation clamp can subsequently hold the last available control past a truncated arc.

  Check engine success before extracting the answer, reuse a checked integration/completion helper for reporting, and explicitly verify `t(end) == tf` before quadrature or flying the control. Either enforce the terminal/flown gates here or return an explicit unsuccessful status. The tests enforce some gates; `run_cartpole_pmp` itself currently only measures and reports them.

- **[ROBUSTNESS]** `cartpole_pmp_prop.m:77–81,110–119` — **The STM assembly is correct, but its accuracy is not guaranteed by the tolerance settings alone.** The implementation uses 8 state/costate variables plus 64 STM entries, seeds `eye(8)`, consistently uses MATLAB column-major `reshape`, and applies the correct left multiplication \(D f(y)\Phi\).

  The complex-step path through `cartpole_pmp_rhs`, `cartpole_field`, and `cartpole_state_jac` contains no evident non-holomorphic operation. In particular, both transposes in `cartpole_pmp_rhs.m:56–58` are correctly nonconjugating. The argument-count branches are harmless. `abs`, `max`, PCHIP, and the control-query clamp elsewhere are **outside** this differentiated path.

  Ordinary nested `imag(f(x+ih))/h` differentiation would indeed corrupt the shared imaginary channel. The implementation avoids that nesting. This does not forbid multicomplex or suitable automatic differentiation.

  One missing contract restriction matters: the outer complex-step formula assumes a **real base state**. The propagator is not itself safely complex-steppable when `needSTM=true`. Validate/document real inputs. Keep the generated inner derivative; there is no mathematical reason to replace it.

- **[ROBUSTNESS]** `run_cartpole_pmp.m:87–100,113–115`, `tests/test_cartpole_pmp.m:104–106` — **This is a genuine non-orbital engine consumer, but the tests validate much less of multiple shooting than the promotion argument suggests.** The supplied problem has a real 8-dimensional propagator, four free initial costates, four terminal constraints, a full node seed, a requested segment grid, and a fixed-final-time option. That is substantive problem adaptation, not merely a renamed orbit wrapper.

  However, all reported trajectories are reconstructed from **only the first costate column**. An engine regression leaving bad interface states, ignoring the requested segment grid, or returning stale noninitial nodes could remain invisible if it still supplies the correct `lam0`. The K=16 check neither verifies convergence nor checks its terminal miss.

  Independently recompute every segment continuity residual from the engine’s returned nodes; check actual grid size, fixed horizon, fixed initial-state components, and the terminal residual for both K values. Add a perturbed-seed case, a nonuniform grid case, and an engine-integration test for rejected propagations. A spy can confirm the requested segment structure and STM requests. These tests cost extra propagations, but not necessarily extra full solves. This demo does not establish coverage of free time, extra parameters, switching, or other engine paths.

- **[STYLE]** `run_cartpole_pmp.m:119–121,135–138`, `README.md:60–62` — **The reported stationarity residual is an unlabelled construction identity.** `U` is assigned `-Lam.'*G/2` immediately before computing `abs(2*U + Lam.'*G)` from the identical operands. Barring exceptional arithmetic, that residual is zero even for completely wrong costate evolution.

  Label it “control-reconstruction arithmetic consistency,” or omit it from the headline numerical evidence. The end-to-end centred-FD test is explicitly and honestly labelled this way; the headline `statMax` is not. In the supplied source I see only one explicitly labelled arithmetic-consistency test, not two. The direct fixture’s pinned-boundary caveats are also honest, but those are separate contract checks.

- **[STYLE]** `gen_state_jac.m:4–6,19,37–40` — **The generator does not write the file its header says it writes.** It writes `cartpole_state_jac_raw.m`; a human must then transfer/wrap the body. That manual step is also where stale generated derivatives can enter.

  Correct the documentation, or commit the generated raw function and have the public wrapper call it. Better still, construct the symbolic dynamics by calling `cartpole_field` with symbolic inputs rather than maintaining another physics copy. That change costs little and makes regeneration reproducible; retain a separate mechanics-based oracle in the tests.

- **[EFFICIENCY]** `run_cartpole_pmp.m:125–127` — **The flown-control callback repeatedly constructs a PCHIP interpolant.** Precompute `pp = pchip(t,U)` once and evaluate it with `ppval` in the callback. This avoids repeatedly processing 2001 samples during adaptive integration without changing the intended interpolant. The clamp is not a complex-step hazard, but should only operate after successful endpoint validation.

## Check-by-check test audit

The distinctions below matter: a check can be redundant without being impossible to fail, and a test against a production function’s output is not automatically a tautology merely because it restates that function’s specification.

- **[ROBUSTNESS]** `tests/test_cartpole_field.m:59–71` — **Useful compatibility tests, but no independent physical oracle and an overstated complex-safety check.**
  - **Line 59, dynamics versus helpers:** catches altered signs, denominators, coefficients, or control insertion relative to the inherited helpers. The \(10^{-12}\) tolerance is meaningful for these samples. It cannot catch the shared mechanical error identified above.
  - **Line 60, `G` versus helper control differences:** catches an incorrect affine-control column, including unwanted control in the kinematic rows. This is genuinely independent of the module’s `G` calculation.
  - **Line 61, helper second differences:** catches nonlinear control dependence introduced into the **helpers**. It does not directly test a cart-pole-module implementation change; it validates an oracle assumption. It is not a tautology because those helpers are called independently.
  - **Lines 65–66, kinematics:** catches swapped rates, incorrect kinematic rows, and control leakage. A valid exact contract check.
  - **Lines 70–71, imaginary propagation:** catches complete removal of the angular imaginary perturbation from either output, or non-finite evaluation at this point. It does **not** establish complex-step safety: conjugation, partial `real`/`abs` use, and many wrong derivatives still propagate some imaginary part.
  
  Add the independent mechanical-power/linearization checks, and replace the safety claim with derivative comparisons against real finite differences over multiple components and states.

- **[ROBUSTNESS]** `tests/test_cartpole_state_jac.m:59–66` — **The real-input derivative checks are sound; the “no nested-step damage” check does not test its claim.**
  - **Line 59, Jacobian versus complex step:** catches incorrect generated coefficients, fixed-control differentiation errors, and reshape/indexing errors relative to `cartpole_field`. The \(10^{-12}\) threshold is strong at the sampled scales.
  - **Line 60, Jacobian versus central differences:** catches discrepancies independently of complex-step evaluation, including some non-holomorphic field changes. The \(10^{-7}\) absolute threshold is reasonable here. Contrary to the header at lines 8–10, it cannot expose a physical error shared by the field and symbolic derivation: it differentiates that same field.
  - **Line 64, real/finite:** catches complex contamination and non-finite entries at one real sample; even a zero matrix passes.
  - **Line 66, finite complex evaluation:** catches exceptions or non-finite output only. A nested-complex-step implementation can return a finite but entirely wrong matrix. So can `real(A)`. This is an unlabelled weak smoke test with an incorrect explanatory label.

  Compare `imag(A(x+ih e_j,u))/h` with real finite differences of `A`, and include perturbations of **`u`**, which becomes complex on the outer PMP differentiation path. Rename the current finite-value test if retained.

- **[ROBUSTNESS]** `tests/test_cartpole_pmp_rhs.m:61–69` — **These mostly test real implementation outputs effectively, but do not establish trajectory minimality.**
  - **Line 61, stationarity:** catches a wrong sign, missing factor \(1/2\), wrong costate indices, or an incorrect returned control. Unlike `out.statMax`, this checks `u` returned by the function under test; it can fail under a production mutation.
  - **Line 62, positive Hamiltonian gaps:** can catch sufficiently wrong returned controls. But once exact stationarity passes, each gap is algebraically \(d^2\); this adds no independent evidence of optimality. Small control errors can pass the gap test, though line 61 catches them. Label it as a pointwise Hamiltonian-minimization check, not a trajectory-minimization test.
  - **Line 63, state rows:** catches wrong state dynamics or use of a different control in `dy` than the returned `u`.
  - **Line 64, costate rows:** catches the costate sign, missing terms, wrong generated Jacobian use, and erroneous differentiation through only the dynamics after control substitution. This is a good independent check against fixed-control differentiation of the scalar Hamiltonian. \(10^{-9}\) is not worryingly loose at these sample scales.
  - **Lines 68–69, zero costates:** catches constant offsets in the control/costate equations and a completely frozen state. `any(dy0(1:4) ~= 0)` does **not** verify the claimed uncontrolled drift; arbitrary nonzero state rows satisfy it.

  Compare the zero-costate state rows explicitly with `F`. Add some solution-scale costates; current random costates are \(O(1)\), versus the reported \(O(10^3)\) solution.

- **[ROBUSTNESS]** `tests/test_cartpole_pmp_prop.m:36–65` — **Good local STM consistency coverage, but weak failure discrimination and no independent propagation-time oracle.**
  - **Line 36, shapes:** catches most dimension errors. `numel(yE)==8` accepts a row vector despite the stated column-vector contract. Use `isequal(size(yE),[8 1])`.
  - **Line 46, STM versus finite differences:** catches wrong STM multiplication order, flattening, seeding, omitted derivative terms, and complex-step corruption when they affect this trajectory. This is a real derivative test, not a construction identity.
    
    However, the relative metric normalizes every error by the **largest** STM entry, potentially hiding errors in small entries. Also, central differences at \(h=10^{-6}\) amplify propagation error by roughly \(10^6\). `RelTol=1e-12` is not a guarantee of \(10^{-12}\) global state error, so the \(10^{-6}\) comparison is plausible but not certified. Use componentwise absolute-plus-relative tolerances, a step-size sweep, and representative points on the converged arc.
  - **Line 49, zero-time identity:** catches a broken zero-time branch. It does not verify identity seeding in the positive-time integration; that is a separate code path.
  - **Line 53, two half steps:** catches restart inconsistencies and excessive integration error. A consistently time-scaled or time-reversed flow still satisfies this property; even the identity map does.
  - **Line 54, STM composition:** catches many composition/variational inconsistencies. It is a meaningful semigroup test, not a tautology, but a coherently wrong flow and its derivative can pass.
  - **Line 57, STM/no-STM agreement:** catches differing state evolution between the two branches. It does not independently establish that either evolves the requested field. It also never checks that the no-STM second output is `[]`.
  - **Line 65, blow-up throws:** catches returning normally for this particular extreme input. **Any exception passes**, including an indexing bug or deliberate `error('bug')`. It does not test the advertised identifier discrimination, and covers only `needSTM=false`.

  Add a short-time check against `cartpole_pmp_rhs`, or a separate propagation oracle, plus deterministic completion-gate tests. Require the exact collapse identifier and test programming-error propagation in both branches. An additional symplectic-STM check is useful structural coverage, but is not a replacement for the derivative comparison.

- **[ROBUSTNESS]** `tests/test_cartpole_pmp.m:52–106` — **The K=8 endpoint gate is strong; several other gates are diagnostics or branch baselines, not optimality or engine-structure tests.**
  - **Line 52, `info.converged`:** catches a false success flag from an otherwise running solve. It trusts the engine’s own definition; a wrongly hard-coded `true` passes.
  - **Line 53, terminal miss:** is the strongest end-to-end check here. It catches an incorrect initial costate, wrong endpoint formulation, many fixed-time regressions, and false convergence that fails fresh propagation. Add an explicit final-time check so “last returned point” cannot substitute for \(t_f\).
  - **Lines 66–67, Hamiltonian constancy:** genuinely catches many costate/control inconsistencies and integration errors. It is not true merely because `U` was reconstructed. However, every accurately integrated canonical arc conserves it, including nonminimizing arcs and arcs of the wrong physical model. Some coordinated errors, such as scaling the entire canonical vector field, also preserve it. The \(10^{-9}\) relative gate is meaningful against the reported \(2.4\times10^{-11}\); explicitly reject non-finite values and handle a nearly zero Hamiltonian scale.
  - **Lines 86–88, centred-FD stationarity:** honestly labelled arithmetic consistency. With the current reconstruction it is zero algebraically, up to cancellation. It can catch a stale/misindexed stored `U` versus `(X,Lam)`, but cannot validate the PMP evolution. “Cannot fail by construction” is slightly too absolute: storage mutations and floating-point cancellation can fail it. The \(10^{-6}\) threshold reflects cancellation, not optimality accuracy.
  - **Line 90, flown terminal miss:** catches control sampling/interpolation errors, stale control/state pairing, and inconsistencies between state evolution and recovered open-loop control. This is useful and not tautological. It shares `cartpole_field`, so it cannot validate the physical equations. The \(10^{-6}\) gate is sensible for this interpolated-control check, but cannot substantiate the \(10^{-9}\) shooting claim on its own.
  - **Lines 93–94, direct cost within 1%:** catches large cost-scale errors and substantially different-cost branches. It admits discrepancies about thirteen times the reported 0.075% gap. That is a coarse cross-method gate, not a precision or minimality check. Establish a mesh-refinement error budget before tightening it.
  - **Line 98, control RMS within 2%:** catches wrong sign, appreciable amplitude errors, timing shifts, and a different control branch. It permits localized discrepancies and offers no optimality proof. Its tolerance is acceptable as a coarse shape test.
  - **Lines 101–102, stored `lam0`:** catches branch changes and numerical/model changes that move the root beyond \(10^{-8}\). It is a strong regression pin, but its source of truth is a previous run of the same implementation. It preserves existing mistakes just as effectively as correct results.
  - **Lines 105–106, K=16 agreement:** catches segment-count sensitivity in the returned initial costate. It passes if K is ignored, and can pass despite K=16 reporting failure or returning bad interface states. Check convergence, actual segment count, continuity, fixed time, and terminal miss independently for K=16.

  None of these checks distinguishes a saddle from a local minimum if it shares the stored branch. A projected second-variation/Hessian test or a conjugate-point calculation would address local minimality; additional direct mesh refinement and multistart would strengthen, but not prove, global optimality.

- **[ROBUSTNESS]** `tests/test_direct_ref.m:36–69` — **This validates a feasible discretized fixture against inherited equations, not an unconstrained optimum or trustworthy multiplier seed.**
  - **Lines 36–38, shapes:** catches missing/wrong-sized arrays. Arbitrary or non-finite multiplier values pass; only `muDefect`’s dimensions are inspected.
  - **Lines 39–40, horizon:** catches changed endpoint times and `tf`. It does not check interior times for finiteness, ordering, duplicates, or the documented uniform mesh.
  - **Lines 41–42, initial state:** catches incorrectly pinned initial boundary data. The caveat that this is not propagation accuracy is honest.
  - **Lines 43–46, terminal state:** catches a grossly wrong pinned terminal target. Again, not tautological under fixture mutation. Its \(10^{-6}\) tolerance admits a boundary error \(10^4\) times the generator’s requested constraint tolerance; that is unnecessarily loose for this committed fixture.
  - **Lines 47–49, inactive bound:** catches contact with the declared force bound, provided the metadata are valid. It does **not** establish optimality: any feasible suboptimal control comfortably inside the bound passes. Non-finite or incorrectly enlarged `uMax` also undermines the check.
  - **Line 68, trapezoidal defects:** catches corrupted nodes, controls, time increments, and disagreements with the helpers. The \(10^{-12}\) gate is strong for the committed data. It proves discrete feasibility—not continuous integration accuracy or freedom from a shared physics error.
  - **Line 69, positive finite cost:** catches zero, negative, NaN, and infinite `J`. Any unrelated positive number passes. This is a very weak smoke test, not cost validation.

  Pin the specified parameters and metadata, check finiteness of every fixture array, and recompute `J = trapz(R.tN,R.U.^2)`. Validate multiplier stationarity, not just multiplier shape. For example, at an interior node of the uniform mesh, the stored defect convention implies
  \[
  2hU_j-\frac h2 G_j^\mathsf T(\mu_{j-1}+\mu_j)\simeq0.
  \]
  Interior adjoint KKT residuals should also be checked.

  Finally, `gen_direct_ref.m:49–50` accepts any positive `fmincon` exit flag and discards its diagnostics. Save and check feasibility and first-order optimality before committing a fixture. An inactive bound makes a properly stationary discrete solution locally unconstrained; it does not transform mere feasibility into an optimum.

- **[ROBUSTNESS]** `tests/test_cartpole_field.m:41–57`, `tests/test_cartpole_state_jac.m:46–57`, `tests/test_cartpole_pmp.m:65,84` — **Several maximum-error accumulators lack explicit finite-value protection.** MATLAB reductions can omit missing/NaN values, so a mixture of finite and non-finite components can yield an apparently acceptable error statistic. The isolated finite-value checks elsewhere do not cover every random sample or every reported trajectory.

  Assert `all(isfinite(...))` on each tested array before reducing errors, and use explicit missing-value behavior where appropriate. A numerical test should never rely on how `max` happens to treat NaNs.

## Plausibility of the reported numbers

- **[STYLE]** `README.md:24–28`, `tests/test_cartpole_pmp.m:22–24` — **The reported numbers are numerically plausible, but they do not validate the specified mechanics or the claimed discretization explanation.** \(J\approx2202.87\ \mathrm{N^2\,s}\) over five seconds implies an RMS force of about \(21.0\ \mathrm N\), compatible with a \(56\ \mathrm N\) peak. A 0.075% direct/indirect gap is plausible for a 201-node trapezoidal transcription, but attributing the entire gap to discretization requires mesh refinement and quadrature-error checks.

  A \(2.5\times10^{-11}\) endpoint miss and \(2.4\times10^{-11}\) relative Hamiltonian variation are also credible for a smooth, well-resolved PMP arc with these tolerances. They can occur for the supplied **wrong dynamics**.

  For the actual hanging-to-upright mechanics, the bob gains \(2m_2gL=39.2\ \mathrm J\). That is not comparable directly with \(\int u^2dt\), which is not energy. The relevant independent physical check is
  \[
  E(t_f)-E(0)=\int_0^{t_f}u\,\dot q_1\,dt.
  \]
  Add this power-balance check using energy derived from the geometry. It would expose the shared sign error that the present cross-checks preserve.

## Verdict

**As written, this does not solve the classic cart-pole problem as specified: the angular dynamics have the wrong sign.** Conditional on accepting those inherited equations, the PMP derivation and STM assembly are sound, and the reported results are consistent with accurately finding a normal extremal—not with proving a minimum.

The module genuinely exercises a non-orbital, fixed-time shooting problem. Its tests would catch many algebraic, STM, and K=8 endpoint regressions. They are **not yet strong enough to support the broader engine-regression claim**: documented test execution need not fail automation, multiple-shooting continuity is unchecked, K=16 success is unchecked, error discrimination is untested, and the physical oracle shares the central mistake.