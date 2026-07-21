## Findings

- **[CORRECTNESS]** `lowThrust_GTO_tulip/solve_minfuel_indirect.m:30` - The optional `epsSchedule` default is guarded with `nargin < 8`, but the function has 8 required inputs plus the optional 9th. Calling the documented 8-argument form leaves `epsSchedule` undefined at line 48. Fix to `if nargin < 9 || isempty(epsSchedule)`.

- **[CORRECTNESS]** `lowThrust_GTO_tulip/run_gto_tulip_minfuel.m:64` - `solve_minfuel_indirect` is called with arguments shifted: `rvf` is passed as `m0`, `tf` as `rvf`, `lamSeed` as `tf`, etc. This will fail at runtime and also poses the wrong shooting problem. Fix to pass the mass fraction explicitly, e.g. `solve_minfuel_indirect(rv0, 1, rvf, tf, lamSeed, Tmax, c, muStar, schedule)`. Also rename line 54’s `m0` to `m0kg` to avoid confusing kg mass with nondimensional mass fraction.

- **[ROBUSTNESS]** `lowThrust_GTO_tulip/lt_pmp_eom_minfuel.m:64-65` - `alpha = -lambda_v / sqrt(sum(lambda_v.^2))` has no guard for near-zero primer magnitude. Bad shooting iterates or singular/coast-like segments can produce `NaN/Inf` and kill `ode113`. Add a complex-step-safe guard using `real(lamvMag)` in the branch; for very small `lamvMag`, return a bounded fallback direction or stop with a clear diagnostic.

- **[ROBUSTNESS]** `GTO_tulip/costate_seed_from_nlp.m:75` - If no NLP node satisfies `s > 0.9`, `burnIdx` is empty and `burnIdx(1)` errors. This can happen for low-throttle smoothed/interior-point solutions or failed NLP outputs. Add an explicit check and either lower/adapt the burn threshold or throw a descriptive error before using `burnIdx(1)`.

- **[ROBUSTNESS]** `GTO_tulip/costate_seed_from_nlp.m:84-90` - The scale gauge claims to enforce `S = 0` at the first throttle switch, but `swIdx` is only the node before a `s > 0.5` logical change, not the switch time. If there is no switch, `swIdx` is empty and `lamSeed` can become empty. Fix by bracketing the crossing and interpolating the actual `s = 0.5` time; add a no-switch fallback/error.

- **[ROBUSTNESS]** `GTO_tulip/NLP_lowThrust_GTO_Tulip_minfuel.m:129-150` - The warm start is described as “exactly feasible,” but the controls are clipped to `s = 0.98` on the burn and `s = 0.02` on the coast while the states come from full-burn/coast propagation. Therefore the trapezoidal dynamics are not exactly consistent with the supplied controls. Either use exact controls `s = 1/0` for a truly feasible burn+coast seed, or re-propagate the warm-start states with the clipped interior controls and update the target accordingly.

- **[CORRECTNESS]** `lowThrust_GTO_tulip/run_gto_tulip_minfuel.m:79` and `GTO_tulip/NLP_lowThrust_GTO_Tulip_minfuel.m:202` - `nCoasts = sum(abs(diff(sign(S)))/2 == 1)/2` can return fractional coast counts, e.g. one burn-to-coast switch gives `0.5` coast arcs. Count contiguous intervals with `S > 0` or `u < 0.5` instead.

## Verified correct

- `lowThrust_GTO_tulip/lt_pmp_eom_minfuel.m:63-75` - Min-fuel switching function, smoothed throttle law, thrust direction, and costate ODE signs are PMP-consistent for fixed-`tf`, free-final-mass min fuel.
- `lowThrust_GTO_tulip/lt_pmp_eom_minfuel.m:50-64,73-80` - Complex-step-safe patterns are used: `sqrt(sum(x.^2))`, nonconjugating `.'`, and no complex-valued branch conditions.
- `lowThrust_GTO_tulip/shoot_residual_minfuel.m:34-46,55-56` - The 7 residuals and complex-step Jacobian column pattern are correct.
- `GTO_tulip/lt_dynamics_throttle.m:45-66` - Dynamics and Jacobians are correct for the scaled-control convention `w = s*alpha`, cone `w.'*w = s^2`, acceleration `Tmax*w/m`, and mass flow `-Tmax*s/c`.
- `GTO_tulip/nlp_constraints_minfuel.m:48-56,65-104` - Fixed-`tf` trapezoidal defects, cone equality, sparse triplet counts, and gradient transpose convention for `fmincon` are correct.
- `GTO_tulip/solve_minfuel_nlp.m:55-84` - Fixed-time final-mass objective, endpoint bounds, and return of `lambdaS.eqnonlin` are correct.
- `GTO_tulip/costate_seed_from_nlp.m:36-72` - Costate transition ODE, burn-node primer cross-product rows, terminal `lambda_m(tf)=0` row, and smallest-singular-vector extraction are structurally correct.
