# oclib — the cross-folder optimal-control library (`+oc`)

The shared home for optimal-control machinery used by MORE THAN ONE
top-level campaign folder (`orbit_transfer`, `booster_landing`, future
`missiles` Phase 2). Created 2026-08-09 per the plan in
`../OCP_UNIFYING_MATH.md` §5. MATLAB package form: `addpath('.../oclib')`
once, then call `oc.<function>`.

**Admission rule (measure before extracting):** a function moves here only
when a second top-level consumer exists and the move carries an equivalence
gate (bit-identical outputs on both consumers' stored flagship data).
Campaign folders keep thin delegates at the old paths so no caller breaks.

## Contents

| function | what | consumers | equivalence gate |
|---|---|---|---|
| `oc.duals_to_costates` | The covector mapping: defect KKT multipliers → costate samples. Owns ALL scheme-specific station-association rules (Hermite–Simpson **midpoints**, trapezoid left-nodes, trapezoid-nodal weighted average), the empirical sign vote, the λ_t = +1 check. | `orbit_transfer` (catalog harvest, foc gates — via the `costate_common` delegate) + `booster_landing` (G5 primer gate) | orbit: golden cells + bitwise harvest A/B; booster: G5 primer angle reproduced at the acos machine floor |
| `oc.local_residual` | The local-residual engine (G1): per-interval re-integration RESTARTING from the transcription's own left node; returns the raw per-component miss at the right node (dimensional splitting stays with the caller — CR3BP vs MEE layouts differ). Sibling of `fly_control`, which CARRIES the state. | `orbit_transfer` DRO gate (`dro_residual`, rerouted) + `orbit_transfer` earth-MEE gate (`mee_residual`, new). **Admission-rule note:** both consumers are inside `orbit_transfer` — the second *top-level* consumer (booster R1's tracking-error gate) is anticipated, not yet real. Recorded rather than hidden. | DRO: old-vs-new max diff 1.7e-18 on both stored reference solutions + 3/3 instrument harness; unit test `tests/test_local_residual` (exactness, injected-error locality, shape) |
| `oc.fly_control` | The flown-control engine: integrate an RHS closure (dynamics + the caller's control reconstruction) over a node grid, per-interval-restart or single-span, configurable integrator/tolerances. Control reconstruction stays with the consumer — it is domain policy (annulus splits, throttle clamps). | `orbit_transfer` G1b (`flown_control_error`, perInterval/ode113) + `booster_landing` G2 (span/ode45) | orbit: globKm 0.891913 / 8.443927 at d = 0.000e+00; booster: G2 residuals 0.0088486 m / 0.000299075 m/s / 1.02479e-05 kg identical |
| `oc.ms_bvp` | The generic multiple-shooting two-point BVP engine: the problem arrives as three closures (`prop`/`rhs`/`terminal`), the engine owns only structure — the unknown vector `[y1(freeIdx0); y2..yK; tf]`, the block-bidiagonal-bordered Jacobian assembled from segment STMs, the trust-region-dogleg solve, the iterate guards, and `opts.fixedTf` / `prob.nExtra` as real structural switches. Problem physics never enters: the CR3BP min-time binding (`ms_tfmin`), the fixed-t_f min-energy/min-fuel bindings and the cart-pole PMP both hand it closures and nothing else. | `orbit_transfer` (every costate campaign, via the `costate_common` delegate) + `optimal_control_examples/ex3_cart_pole_pmp` (directly, as `@oc.ms_bvp`) | orbit: `golden_cells` 20/20 **bit-identical** pre- vs post-move (z, ‖R‖, iteration counts to the last bit) PLUS the four-campaign ladder re-solve of one stored cell each — HALO 4.551e-14, DPO 1.130e-14, HALO_HALO 3.478e-14, GTO 2.720e-14 max \|Δt_f\|, OK flags matching; cart-pole: `test_cartpole_pmp` 9/9 (terminal miss 2.47e-11, Hamiltonian constant to 2.4e-11, λ(0) regression 8.2e-12) |
| `oc.ms_conjugate_test` | The Jacobi (conjugate-point) test on a converged multiple-shooting extremal, both time conventions: free-final-time (`det([Φ_xλ·P, f(t)])` on the costate-scaling quotient) and fixed-final-time (same instrument, no flow column, `quotientDir = []`, rows 1:7 = the full state — a Jacobi field must vanish in r, v AND m). A necessary-condition check (Jacobi), not a sufficiency proof; segment STMs come free from `ms_bvp`'s `opts.keepSTMs`. | `orbit_transfer` (every costate campaign, via the `costate_common` delegate). **Admission-rule note:** the second *top-level* consumer — the cart-pole min-energy study script landing later in this same sub-project (Task 4, `optimal_control_examples`), calling `oc.ms_conjugate_test` directly on a fixed-t_f PMP extremal with no orbit, no CR3BP quantity and no pumpkyn call in it — is anticipated, not yet real at this commit. Recorded rather than hidden. | orbit: `golden_cells` 20/20 **bit-identical** pre- vs post-move (dro \|z-z8\| 3.7e-11/iters 1/normR 1.2e-13, halo 4.4e-11/1/6.0e-14, dpo 1.2e-10/1/2.9e-14, harv \|z-zref\| 4.5e-12/iters 2 — every number unchanged) PLUS the four-campaign ladder re-solve of one stored cell each, unchanged to the printed digit — HALO 4.551e-14, DPO 1.130e-14, HALO_HALO 3.478e-14, GTO 2.720e-14 max \|Δt_f\| (dZ8 1.689e-12/1.466e-12/2.354e-13/2.895e-13, OK flags matching); `test_conj_fixedtf`, `test_conj_coverage`, `test_conj_spectrum`, `test_conj_resolve`, `test_folder_rules` all green, PASS/FAIL output byte-identical pre- vs post-move |

## Tests

`tests/` holds one unit test per function, each on a problem with an exact
answer, so they run from a fresh clone with no campaign data:
`test_duals_to_costates` (station rules for the three schemes, h-scaling,
the sign vote in both orientations, the lambda_t check, refusals),
`test_fly_control` (exact flow in both modes, control looked up from global
time, solver option, refusal), `test_local_residual` (exact nodes, injected
error stays local). `ms_bvp` and `ms_conjugate_test` are the exception to
the one-test-per-function shape: their tests stayed in
`orbit_transfer/costate_common/tests` (`test_ms_bvp_extra`,
`test_ms_bvp_fixedtf`; `test_conj_fixedtf`, `test_conj_coverage`,
`test_conj_spectrum`, `test_conj_resolve`) and call the delegate, so they
exercise the forwarding and the engine at once — moving them would have cost
that coverage for nothing. Each returns `ok`. Every function except
`ms_conjugate_test` also carries a `nargin == 0` self-demo; `duals_to_costates`'s
demo reads a gitignored DRO_tulip results file, so it runs only on the
machine that holds it. `ms_conjugate_test`'s only `nargin` check is a
default-argument test (`nargin < 2`), not a demo — it needs a converged
multiple-shooting `info` struct with segment STMs to run against, which its
unit tests in `costate_common/tests` supply as fixtures instead.

## Why this exists (the one-sentence version)

Booster G5 compared a Hermite–Simpson segment dual against the *node*
control instead of the *midpoint* control it is the multiplier for — the
exact bug `orbit_transfer` had already found by review and fixed. Two
folders, one rule, no shared home → the bug was independently re-invented.
This package is the home.

## House style

Pumpkyn conventions: `%% Purpose / Inputs / Outputs / Revision History`
headers, `nargin==0` self-demos, no Code Analyzer pragmas. Physics never
lives here — forward models stay per-domain; this package owns
transcription-side and PMP-side *structure* only.

## Roadmap (from OCP_UNIFYING_MATH.md §5, in order)

1. ~~`oc.duals_to_costates`~~ (this move)
2. ~~`oc.fly_control`~~ (done)
3. ~~(part) local-residual engine~~ — `oc.local_residual` landed 2026-08-25
   (see table). ~~`oc.ms_bvp`~~ landed 2026-09-16: the cart-pole PMP-BVP demo
   (`optimal_control_examples/ex3_cart_pole_pmp`) was built as the cross-folder
   integration test, and being a second TOP-LEVEL consumer with no orbit,
   no CR3BP quantity and no pumpkyn call in it, it is what admitted the
   engine (`costate_common/TODO.md`, cleanup step 13). ~~`oc.ms_conjugate_test`~~
   landed 2026-09-17 (see table for the admission-rule note and the
   equivalence gate) — the delegate stays at the old path. `ms_tfmin` and the
   seed builder `seed_from_z8` moved into `costate_common` on 2026-08-26.
   **Still open:** `arclength_ms`, `newton_fixed_q`, `conj_resolve` and
   `lift_space_dim` — no second top-level consumer yet, so the admission
   rule still refuses them.
4. (deferred) transcription defect builders — build for the next new
   campaign, adopt backward if the diff supports it
