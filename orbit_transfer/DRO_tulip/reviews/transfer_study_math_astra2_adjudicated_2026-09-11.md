# transfer_study math review, second pass: adjudication and corrections (2026-09-11)

Review adjudicated: `transfer_study_math_astra2_2026-09-11.md` (GPT-6 Astra,
xhigh, 461 s, $1.79, on the corrected sources + the first adjudication + the
run output). Host: Claude (Fable 5.1). Every finding was checked against the
code before a verdict; the fixes were applied in the order Astra ranked them.

## Summary

Astra confirms the first round's corrections: the full-gap formula, the
mass-row throttle, the frozen-control adjoint comparison including the mass
row, and both kernel identities for the 14-state system are correct.
Per-condition: N2, N4, N5, S1, S2, V1, X2 CORRECT; N6, S3, S4 INCOMPLETE.
Three real gaps remained, all verified against the code and all closed here:

1. The dense S4 gate never read `nInteriorCand`, so a "near-miss" passed,
   and the plateau classifier was not a zero-exclusion test.
2. `certify_root` enforced less than the study (no full gap, throttles, X2,
   lift residuals, `lift_margin`, dense scan).
3. The lift residual carried two normalisations under one option value;
   the endpoint check multiplied a mixed norm by a length.

## Findings: verdict and action

| # | Finding | Verdict | Action |
|---|---|---|---|
| 2.7 | Dense gate omits `nInteriorCand`; plateau is not exclusion; nested grids can keep the same nearest node (root within h/32) so a true double root reads "near-miss" | CONFIRMED | `conj_spectrum`: refinement on SHIFTED grids, then golden-section location of the sigma_6 minimum; kind = zero (at the floor or a sign change), near-miss (located minimum >= 100x floor), UNRESOLVED (between). `.clear` is the gate; unresolved blocks. Multiplicity = number of singular values at the floor at a located zero |
| 2.7b | Start/endpoint bands discarded unrefined; "sign at t_f carries no information" is wrong; refinement discards sign information | CONFIRMED | class by CONTIGUITY: only a cluster touching sample 1 is the (uncovered) start transient; endpoint clusters are refined through t_f. Header corrected: the scan adds sensitivity, does not close the blind spots |
| 2.2 | Unresolved same-sign final determinant can PASS; `resolvedTol` only acts on a sign change | CONFIRMED | `ms_conjugate_test`: any live sample with sigma ratio <= resolvedTol makes the verdict UNDETERMINED unless an interior root was already found |
| 2.8 | Production certifier gates `dirGap` only; no X2, lift residuals, `lift_margin`, dense scan | CONFIRMED (mitigated at sweep level by `second_order_pass`, which Astra did not see) | `certify_root`: gates 2b/5/7 widened to the study's set; `test_certify_enforcement` (NEW): 5 tolerance squeezes + 4 field injections each refuse by name; the anchor certifies untouched |
| 2.1 | `fullGap` is not the field's gap when the two throttles disagree: differs by (T/c) lambda_m (|b| - u_mass) | CONFIRMED | `.fieldGap` evaluated; the gated number is max(fullGap, fieldGap); mutation test: half mass flow opens fieldGap 1e-3 while fullGap stays 1e-16 |
| 2.6a | `lift_margin` header overstates Eckart-Young: sigma_6 gives rank >= 6, exactly 6 needs the exhibited lift | CONFIRMED | header rewritten; the reason string says so |
| 2.6b | Two-tolerance difference is a sensitivity estimate of ONE component; the printed pair is 1e-12/1e-9, not 1e-10 | CONFIRMED | stated in the header and printed with the settings; the study labels it a sensitivity estimate |
| 2.6c | Lift tolerance: study tests residual/|lam| < 1e-4, `lift_margin` tests residual < 1e-4 sigma_1 | CONFIRMED | ONE normalisation: the backward error |C lam|/(sigma_1 |lam|) (`gates.nullResidRel`), tol 1e-6 in the study, the certifier and `lift_margin` (liftTol default 1e-6) |
| 2.6d | `lift_space_dim` header claims the cap makes dim S "never > 1" | CONFIRMED | header corrected; states that dim S >= 1 is threshold-induced |
| 2.5 | Endpoint error: six-state norm times lStar is not a distance; phase not wrapped; identically zero at s = 0 | CONFIRMED | position and velocity reported separately (km, m/s), phase wrapped with mod as `phase_state` does, propagation completion asserted, labelled a consistency estimate |
| 2.3 | Kernel identities correct; comment should show the mass-costate component | CONFIRMED (correct as implemented) | comment amended |
| 2.4 | Envelope wording: "no derivatives of the control law", not "does not depend on it"; the error is a 7-row vector norm | CONFIRMED | comment rewritten, mass row derivation spelled out |
| 2.9 | X2 validates the field, not the STM generator behind S4 | CONFIRMED | `test_stm_variational` (NEW): all 14 STM columns vs central differences; symplecticity PHI' J PHI = J; and a frozen-control generator (no d alpha*/d lam_v) shown to BREAK symplecticity, so the check discriminates |
| 3 | Verdict must name the structural theory gaps, not only numerical tolerances | CONFIRMED | verdict adds: subarc normality, free-mass/free-time second-variation reduction, existence of an exact nearby extremal; `report_optimality` aligned (its "same endpoints" wording replaced) |
| 4 | Open theory items correctly stated; two sharpened (subarc normality must carry lambda_m(t_f) = 0 to a moved endpoint with H6 as the normalisation; the reduction must cover reduced-throttle competitors to keep the claim STRONG); plus: existence of an exact extremal near the numerical one | ACCEPTED | `costate_common/TODO.md` updated |
| 1.2 | "shows up in N2" should read "can show up in N2" | CONFIRMED | wording |

## Consequence for the shipped 70 mN library

The library's second-order sheet (`second_order_pass`, FINDINGS 42) was
swept with the OLD plateau classification: 44 cells carried a "near-miss"
verdict that was a plateau, not a located positive minimum, and endpoint
clusters were not refined. `second_order_pass` now records
`conj_unresolved` (NaN for entries swept before today) and the writeback
`meaning` says so. A RE-SWEEP is required before the ship decision; until
then the library's near-miss column is evidence of the old kind.

## Verification

- `transfer_study.m` end to end after all fixes: every line PASS, verdict
  claimed with the new wording. Same solution as before both reviews
  (t_f 17.7976 d, 0.7485 km/s, min|lam_v| 3.2359, min Q 3.8248, H6 9.0x).
  New readings: lift backward error 5.3e-10 (sigma_1 4.1e3), Eckart-Young
  margin 2545x at relTol 1e-12 vs 1e-9, dense scan 192 samples with 1
  endpoint candidate resolved to a located minimum of 8.0e-5 x median at
  t/t_f 0.9994 (cleared), 1 start transient (uncovered), 0 unresolved;
  endpoint interpolant 0.05 m / 2.2e-5 m/s at the arrival phase.
- `costate_common/tests`, all PASS: `test_pmp_pointwise_checks` (18, incl.
  fieldGap mutation), `test_validate_flight` (11), `test_lift_margin`,
  `test_lift_space_dim`, `test_gates_h6_wiring`, `test_mintime_gates`,
  `test_stm_variational` (NEW: 14 FD columns 8.6e-8, lam_v columns 5.2e-9,
  symplecticity 1.9e-10, frozen generator differs O(1) in lam_v columns),
  `test_conj_fixedtf` (17), `test_conj_coverage` (12), `test_conj_spectrum`
  (refuted entry still ZERO), `test_report_optimality`, `test_fly_transfer`,
  `test_h6_margin`, `test_certify_enforcement` (NEW: anchor certifies with
  every new number; 5 squeezes + 4 injections refuse by name),
  `test_certify_crossing` (anchor certifies through the widened stack).
- Lesson recorded: symplecticity of the STM does NOT discriminate a
  generator missing d alpha*/d lam_v (a frozen-control generator is a
  fixed-alpha Hamiltonian linearisation and is symplectic to 1e-11); the
  finite-difference columns do.
- NOT re-run here: the 70 mN library's second-order sweep. Its near-miss
  column is old-classification evidence until re-swept (TODO).
