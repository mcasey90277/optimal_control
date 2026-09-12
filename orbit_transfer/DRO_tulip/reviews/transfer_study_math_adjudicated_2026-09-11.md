# transfer_study math review: adjudication and corrections (2026-09-11)

Reviews adjudicated: `transfer_study_math_astra_2026-09-11.md` (GPT-6 Astra)
and `transfer_study_math_gemini_2026-09-11.md` (Gemini 3.1 Pro). Host: Claude
(Fable 5.1). Every finding was checked against the code before a verdict.

## Summary

Astra's review confirms the mathematics of every condition (Hamiltonian,
mass-costate sign, switching function, spherical Legendre, the H6 identity)
and finds real defects in how the script ENFORCES them. Gemini grades all
nine conditions CORRECT, including two the script itself labels as sampled
evidence, and its one code-level finding is wrong. Where they disagree,
Astra is right in every case that could be checked.

All code-level findings are applied below. The theory items are recorded as
open work for `doc/mintime_second_order_audit.tex` (see TODO).

## Astra findings: verdict and action

| # | Finding | Verdict | Action |
|---|---|---|---|
| 1 | V1 gate `h6Margin >= 1` ignores `h6Ok`, `h6Clearance`; equality passes | CONFIRMED | script, `certify_root`, `report_optimality`: gate on `h6Ok` AND strict margin; clearance printed against the Hamiltonian residual |
| 2 | V2 claims to catch a wrong thrust; `minQmt` carries no `Tnd`/`mu` | CONFIRMED | comment corrected: V2 catches flight/exhaust-speed wiring; thrust and mu errors surface in N2 |
| 3 | S4 PASS possible with incomplete final coverage; caller does not gate; `firstFullRank = nS+1` indexes past `cj.t` | CONFIRMED | `ms_conjugate_test`: missing `Yend` -> UNDETERMINED with `.covered`, `.reason`, `.tFirstFullRank`; script gates on coverage |
| 4 | Negative direction gaps clipped (`dirGap` starts at 0, max only) | CONFIRMED | signed `.gapMin`/`.gapMax` kept; gate is `abs` |
| 5 | N6 recovers the throttle on the acceleration rows only; full control gap not evaluated | CONFIRMED | `u_mass = -c F_m/T` recovered and gated; full gap (direction + throttle) evaluated; mutation tests for both |
| 6 | Last-bracket sign change between resolved nonzero samples is an interior zero, not ENDPOINT | CONFIRMED | FAIL when the final determinant is resolved (`resolvedTol`); ENDPOINT otherwise; `nEndResolved` reported; fixed-tf test expectation updated |
| 7 | Verify `J p(0) = 0` and `p(t)' J(t) = 0` as diagnostics | ADOPTED | measured in `ms_conjugate_test` (`kernelRight`, `kernelLeft`): 3.8e-15 and 8.0e-14 on the anchor |
| 8 | N2/N5 test pumpkyn's field against itself; compare against an explicit reference implementation | CONFIRMED | `mintime_hypothesis_gates` compares pumpkyn's state rows and adjoint rows against the hand-written field and its CasADi Jacobian: 1.8e-16 / 9.3e-15 on the anchor; script line X2 |
| 9 | S3 self-consistency residuals printed but not required; threshold forces at least one small singular value | CONFIRMED | script requires lift residual, Hamiltonian residual AND the Eckart-Young `lift_margin` (2545x on the anchor) |
| 10 | Sign sampling at junctions cannot see touches or close pairs | CONFIRMED (already disclosed) | `conj_spectrum` dense scan (192 samples, two-level refinement) wired into S4 and gated |
| 11 | `validate_flight` checks only endpoint mass, not costates, pointwise mass, positive expected mass | CONFIRMED | all four added, with tests |
| 12 | Periodic-spline seam is a construction property, not an orbit check | CONFIRMED | endpoint interpolant compared against a CR3BP propagation to the phase: 2.2e-8 ND (8.6 m) on the tulip |
| 13 | "First variation vanishes" is the wrong description of N1 | CONFIRMED | wording: the PMP shooting equations are satisfied; the control condition is N6 |
| 14 | N2/N4 not independent evidence; S2 follows from N4+N5+S1 | CONFIRMED | wording in the script |
| 15 | X1 is a second shooting implementation on shared physics | CONFIRMED | wording; X2 is the physics check |
| 16 | "Same endpoints and phases" redundant; final mass free not stated; numerical PASS is not a certificate | CONFIRMED | verdict rewritten per Astra's proposed sentence |
| 17 | Subarc normality needs the analytic-continuation argument | OPEN (theory) | TODO for the audit document |
| 18 | Six-state nonautonomous reduction with free mass needs an explicit second-variation reduction | OPEN (theory) | TODO for the audit document |
| 19 | Initial interval before the first junction is uncovered | CONFIRMED | reported as uncovered; the dense scan samples inside the first segment; short-time sign expansion is a TODO |
| 20 | No between-sample bounds on S1, S2, clearance | CONFIRMED (disclosed) | stated in the section header and the verdict; validated enclosure is future work |
| 21 | N3's 100 km / 10 m/s do not establish a fixed-endpoint extremal | CONFIRMED | wording: a screen on the single-shot flight; N1 is the statement |

## Gemini findings

| Finding | Verdict |
|---|---|
| X1 physics tautology | CONFIRMED (same as Astra 15) |
| `Yend` missing truncates the final segment silently | CONFIRMED in the instrument; cannot occur on this call path (`ms_bvp` always supplies `Yend` under `keepSTMs`); fixed regardless (Astra 3) |
| Undefined function `scalar_verdict` | WRONG. `costate_common/scalar_verdict.m` exists; the reviewer's bundle omitted it |
| S3 sampling can only raise the rank, so it fails conservatively | INCOMPLETE. True of sampling, not of matrix error; Astra's Eckart-Young objection is the right one, now enforced |
| S4, V1, verdict sentence all CORRECT / justified | OVERRULED by Astra 3, 1, 16 |

## Verification

- `costate_common/tests`: `test_pmp_pointwise_checks` (17/17, two new
  mutations), `test_validate_flight` (11/11), `test_conj_coverage` (NEW,
  11/11), `test_conj_fixedtf` (16/16, one expectation updated),
  `test_gates_h6_wiring`, `test_mintime_gates`, `test_fly_transfer`,
  `test_certify_crossing`, `test_report_optimality`, `test_h6_margin`: all
  PASS.
- `transfer_study.m` end to end: every line PASS, verdict claimed; t_f
  17.7976 d, Delta-V 0.7485 km/s, min|lam_v| 3.2359, min Q 3.8248, dim S 1
  (margin 2545x), 0 crossings on 24 junctions and 192 dense samples, H6
  9.0x with clearance 0.889 against |H| 3.3e-8 -- identical to the
  pre-review numbers. No shipped catalog verdict moves: the 70 mN library
  has 0 sign changes, so the last-bracket rule change touches no entry.
- Code Analyzer clean on the edited files (two pre-existing "might be
  unused" notes in `ms_conjugate_test` untouched).
