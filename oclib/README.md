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
| `oc.fly_control` | The flown-control engine: integrate an RHS closure (dynamics + the caller's control reconstruction) over a node grid, per-interval-restart or single-span, configurable integrator/tolerances. Control reconstruction stays with the consumer — it is domain policy (annulus splits, throttle clamps). | `orbit_transfer` G1b (`flown_control_error`, perInterval/ode113) + `booster_landing` G2 (span/ode45) | orbit: globKm 0.891913 / 8.443927 at d = 0.000e+00; booster: G2 residuals 0.0088486 m / 0.000299075 m/s / 1.02479e-05 kg identical |

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
3. `oc.ms_bvp` + `oc.ms_conjugate_test` — promotion unchanged, plus a
   cart-pole PMP-BVP demo as the cross-folder integration test
4. (deferred) transcription defect builders — build for the next new
   campaign, adopt backward if the diff supports it
