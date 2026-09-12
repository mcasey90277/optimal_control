# transfer_study math review, third pass: adjudication and corrections (2026-09-11)

Review adjudicated: `transfer_study_math_astra3_2026-09-11.md` (GPT-6 Astra,
xhigh, 550 s, $1.73), scoped to the NEW code of pass 2: the candidate
resolution in `conj_spectrum`, the widened `certify_root`, the unresolved
rule in `ms_conjugate_test`, the symplecticity lesson, the lift
normalisation. Host: Claude (Fable 5.1). Bottom line as received: not fit to
re-sweep the library yet. Items 1-7 were applied in Astra's order; item 8
(floor from a measured matrix error) stays open.

## Findings: verdict and action

| # | Finding | Verdict | Action |
|---|---|---|---|
| 1.1 | Golden section on the finest-grid argmin resolves one minimum, not a cluster; bracket need not be unimodal | CONFIRMED | `conj_resolve` (NEW, pure): EVERY local minimum of the assembled record is bracketed by its neighbours and located; unimodality of everything evaluated in the bracket is checked afterwards; a non-unimodal bracket is UNRESOLVED |
| 1.2 | Refinement discards evidence (determinants, earlier minima); a floor-level sample can be upgraded by a later search | CONFIRMED | every evaluation is kept; trusted opposite-sign brackets are counted across the whole record; the smallest value ever seen decides, never a later larger one |
| 1.3 | t_f never evaluated by the refinement; endpoint spectrum does not initiate a candidate | CONFIRMED | the last sample is a candidate source; endpoint windows run to t_f; bracket endpoints are in the record |
| 1.4 | Floor and stopping tolerance not tied to the matrix error | ACCEPTED, OPEN | still a policy value; the form Astra gives (safety factor x measured scaled-matrix error at the candidate time) is item 8 in TODO |
| 1.5 | "At the floor" is numerical rank ambiguity, not an exact zero; `nSmall` is numerical corank | CONFIRMED | zero split into `.nZeroSign` (trusted bracket, established) and `.nZeroFloor` (possible rank loss); both block; `.nSmall` documented as numerical corank |
| 1.6 | Contiguity can absorb a real root into the start transient; a rank-deficient scan reads clear; extent not reported | CONFIRMED | the start cluster is refined past its first local maximum (the increasing prefix is structural and is the reported `.tUncovered`); a scan with no trusted sample outside every cluster is `.testable = false` and not clear. Synthetic case: a zero at t = 0.2 inside the transient is caught, uncovered stops at 0.148 |
| 1.7 | Column normalisation hides a vanishing column | CONFIRMED | raw column-norm ratio recorded per evaluation; below `colTol` it is a candidate and UNRESOLVED. Synthetic case: diag((t-t*)^2, 1, ...) is caught |
| 1.8 | `mintime_prop_seg` takes the last row of a possibly incomplete propagation | CONFIRMED | asserts arrival at the requested duration, dimensions and finiteness; errors otherwise |
| 1.9 | Tests do not exercise the resolution mechanism | CONFIRMED | `test_conj_resolve` (NEW, 26 checks) on synthetic matrices: transverse zero at six grid phases incl. on-node and h/32, corank-two zero without sign change on- and off-node, quadratic touch, two wells around a zero and the same wells cleared, zero at t_f and h/5 before it, zero merged with the start, vanishing column, near-miss cleared at 1e-4, unresolved at 3e-6, never-full-rank; every case asserts `.clear` |
| + | (found by the synthetic tests) a V-shaped zero between two coarse samples can read far above the dip threshold at both and, without a sign change, was invisible | CONFIRMED | local minima of the coarse spectrum below `tolLocal` (0.2) are candidates at any depth |
| 2.1 | Gate 2b uses plain comparisons: `[]`, `[0 Inf]`, `-Inf` slip through; `max` hides NaN | CONFIRMED | every gated field validated as a real finite NON-NEGATIVE scalar before aggregation or comparison (`nonneg_verdict`); documented that validation precedes numerical ordering |
| 2.2 | `max(err, new)` drops a NaN; `nSample = 0` / `hRel = 0` leave zero residuals | CONFIRMED | `pmp_pointwise_checks` validates the sample count and step, poisons residuals to NaN on any non-finite evaluation and reports `.nonfinite`; the gates do the same for X2; the certifier refuses on either |
| 2.3 | `lift_margin` svd before validation; `LM.certified` not validated; dense counts unchecked | CONFIRMED | `lift_margin` validates shape, realness, finiteness, options; certifier catches its refusal by name, requires `certified` to be a scalar logical, validates every count as a finite non-negative integer and `.clear` as consistent with them |
| 2.4 | `conjSpectrum = false` returns "certified" | CONFIRMED | DIAGNOSTIC contract: `C.ok = false`, `C.okDiagnostic = true`, `C.fullStack = false`, reason says so; cannot be packaged |
| 2.5 | Mutation coverage: field injections are rejected by N2 first, so later gates unproven | CONFIRMED | TEST SEAM `opts.override` (warns loudly) lets each gate be the first failing one: 17 seam mutations incl. empty, vector, -Inf, NaN, negative, wrong type, inconsistent flags, each refused by name |
| 2.6 | Per-call fences are not a total wall-time limit; unfenced without a pool | ACCEPTED, unchanged | documented; the batch launcher's OS watchdog is the backstop |
| 3.1 | `resolvedTol` is an LU-roundoff screen, not a derived trust threshold | ACCEPTED, OPEN | still 1e-10; the derivation (sigma_min above matrix + LU + SVD error) is item 8's sibling |
| 3.2 | UNRESOLVED rule applied AFTER counting with untrusted signs | CONFIRMED | trusted-positive / trusted-negative / unresolved classified FIRST; roots are brackets of opposite trusted signs; equal trusted signs around unresolved samples establish nothing; ENDPOINT verdict retired |
| 3.3 | Kernel identities should gate | CONFIRMED | `kernelTol` 1e-8; a violation is UNDETERMINED (not interpretable), never a refutation; `.kernelChecked` / `.kernelOk` reported |
| 4.1 | Symplecticity lesson correct; frozen control is a prescribed schedule | CONFIRMED | test comment amended |
| 4.2 | Cheap discriminator: the generator block A(4:6,11:13) = -(T/(m rho))(I - alpha alpha') | ACCEPTED, OPEN | TODO (the FD columns already discriminate) |
| 5 | Backward error admits |C lam|/|lam| up to 4e-3 at sigma_1 = 4e3; add blockwise residuals and the angle to v_7 | ACCEPTED, OPEN | TODO; the anchor's 5.3e-10 / 2.2e-6 are far inside either rule |

## Verification

- `costate_common/tests`, 15/15 PASS in one batch: `test_conj_resolve` (NEW,
  26 synthetic checks), `test_conj_fixedtf`, `test_conj_coverage`,
  `test_conj_spectrum` (the anchor CLEAR with one interior near-miss at
  t/t_f 0.776 and one endpoint near-miss; the refuted 22.05 d entry refuted
  by a trusted sign bracket), `test_pmp_pointwise_checks`, `test_lift_margin`,
  `test_lift_space_dim`, `test_mintime_gates`, `test_gates_h6_wiring`,
  `test_stm_variational`, `test_report_optimality`, `test_fly_transfer`,
  `test_validate_flight`, `test_certify_enforcement` (5 squeezes + 4 field
  injections + 17 seam mutations, each refused by name; the anchor
  certifies untouched with fullStack true; conjSpectrum off is diagnostic
  only), `test_certify_crossing`.
- `transfer_study.m` end to end: every line PASS, same solution (t_f 17.7976
  d, 0.7485 km/s, min|lam_v| 3.2359, min Q 3.8248, H6 9.0x, lift backward
  error 5.3e-10, Eckart-Young 2545x). S4 junction test: PASS, "no sign
  change among trusted samples, none unresolved", identities 3.8e-15 /
  8.0e-14 now gated. Dense scan: 192 samples, 1 interior + 1 endpoint
  candidate, both near-miss with a located positive minimum, 0 unresolved,
  start transient uncovered.
- Anchor scan cost: 0.8 s, 511 evaluations.
- NOT done: the library re-sweep (the launch kit is ready; the classifier is
  now the reviewed one) and open item 8 (floor from a measured matrix error).
