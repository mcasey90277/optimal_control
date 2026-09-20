# Track A mathematics, root identity and the branch map: GPT-6 Astra (xhigh), adjudicated (2026-09-19)

**Bundle.** `phase_sensitivity`, `phase_edge_residuals`, `phase_branches`, `between_sample_bound`, `mintime_hypothesis_gates`, `phase_transversality_check`, `family_map`, `compare_phase_catalogs`, `reproduce_library_70mN`, the day's driver diff, the certifier's H2/H3 hunk, and the audit document's theorem statement, new subsection and reductions section. Raw API, xhigh, 132 KB, 440 s, $1.52. Transcript: `trackA_identity_branches_astra_2026-09-19.md`. None of this code had had an outside review.

## Confirmed CORRECT by the review
The phase-sensitivity signs and derivation (mass costate does not enter; free final time adds no term); `x'(s) = tau f_ballistic`; the CR3BP adjoint equation and `||d f_v / d v|| = 2`; the mass-costate sign `lam_m' = -T|lam_v|/m^2` and the EXACT cancellation `dQ_mt/dt = (d|lam_v|/dt)/m`; the interval inequality given a true slope bound; the third-order edge residual on a smooth branch; comparing all seven costate components; the exact-hit index guard and the pass/fail scoring in `attachPoint`; the Hamiltonian gap identity; `S^2` as an open control manifold and "cited, not instantiated" for Osmolovskii-Maurer; the integral identity `int T Q_mt dt = lam_m(0)`; the fixed-phase scope; the comparator's phase permutation and faster/slower signs. It also gave a strictly better slope bound for free: the Coriolis matrix is skew, so `|d|lam_v|/dt| <= |lam_r|`.

## Accepted and fixed today (test-first; all suites green in a clean batch session)
| finding | host verdict | what changed |
|---|---|---|
| A small flight-time residual does not identify a costate branch; "SAFE / BRANCH" turned a scalar heuristic into an interpolation guarantee | **CONFIRMED BY MEASUREMENT**: across time-consistent s_D edges of the record the costates change by 32% (median), up to 88% | `phase_branches` -> `phase_components`; outputs renamed `.timeConsistentD/.A .component .nComponent`; `.costateJumpD` reported; nothing is called "safe" |
| The under-resolved arrival axis was "not judged" yet its 9 accepted edges merged components | CONFIRMED | components are joined through s_D edges only (each lies in one arrival column) |
| `attachPoint`: NaN or missing costates silently became a flight-time match; two families passing reported one; `abs(rho)` discards the chart's sign (rho < 0 makes the stored rows antiparallel) | CONFIRMED | four-valued status `attached / ambiguous / none / unknown`; non-finite caller costates or a costate-free arc -> `unknown`, fam 0 (the root is anchored, not suppressed); the signed chart is used; the nearest crossing's costate gap is reported when nothing passes |
| `sameRoot` asymmetric, no finite checks, ignores t_f; the registry dedups without the departure phase | CONFIRMED | symmetric tolerance, finite and non-zero required, t_f compared; both phases in the dedup; a non-finite root is refused |
| X3: "no costate is read / independent oracle" overstates; purely relative gate fails near stationary phases; the half-step retry is not recorded | AGREED | renamed a finite-difference consistency test; gate `abs + rel`; `.stepD/.stepA/.absErr*` recorded |
| Greedy family pairing overstates differences (`[6 5; 5 0]`: 10 vs the optimal 6); status codes are not families | CONFIRMED | `matchpairs` assignment; codes < 1 compared as statuses |
| `.adoptArcs = false` refused its own arcs on resume | CONFIRMED | allowed when the campaign's state file exists |
| Audit document: "the two Hamiltonians agree" is WRONG (`H_7 = H_6 - k lam_m`, `p_tau = -k lam_m`; host re-derived it); the L_k inequality does not follow from end samples; "reduction is sound" overstated; X3 numbers stale; "minimum lies between grid points" overstated; (S) "checkable" overstated; gap coefficient wording | ALL ACCEPTED | `doc/mintime_second_order_audit.tex` corrected; H2 => H3 on an exact all-burn extremal noted |

## Accepted, deferred while the 96-phase sheet job ran; APPLIED 2026-09-20 (test-first, 115 checks green in a batch process)
- `mintime_hypothesis_gates`: use `|lam_r|` as the slope bound; rename `.minLamVBound/.minQmtBound` to `...Estimate`; header "no abnormal lift iff dim S = 1" -> a sufficient exclusion of stationary lifts.
- `certify_root` / `audit_phase_catalog`: the refusal text must say "lower-bound ESTIMATE"; the floor's justification.
- `between_sample_bound`: validate finite inputs and non-negative slopes.
- Bind the audit to the catalog by a content key (hash of has_solution, tf, z8), not by location and count.

## Recorded, not done
A rigorous interval bound (comparison system `a' <= kappa b, b' <= a`, needs primary-distance enclosures over the interval); a confirming polish before suppressing a new arc (two-stage identity); Richardson `D_h` vs `D_h/2` on validation cells; branch-identity check of the displaced re-solves; an analytic departure tangent; a canonical problem fingerprint in the comparator; calibrating `tolCos = 1e-4` (0.81 degrees, ~1.4% of a unit vector) against measured within-root and nearest-distinct-root separations.

## The experiment the review provoked
Whether an interpolated costate is USABLE is an experiment, not a residual: interpolate z8 at an s_D midpoint and polish. The first attempt ran unfenced in the shared interactive session and hung it for over 45 minutes (the campaign rule it broke: solver calls are fenced, in a batch process). To be rerun properly.
