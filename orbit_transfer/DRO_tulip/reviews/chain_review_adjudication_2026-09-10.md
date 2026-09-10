# Adjudication — Astra chain review, 2026-09-10

Review: `chain_code_review_astra_2026-09-10.md` (GPT-6 Astra, xhigh, 14,371
reasoning tokens, 505 s, $1.49). 41 findings on the 15-file chain that builds
the 70 mN DRO → tulip costate catalog.

Astra's verdict: *"I would not release this with an unqualified CERTIFIED
claim yet."* That verdict was correct. Its three named classes — wrong-phase
labelling, export paths that trust summaries instead of certificates, and an
unenforced foreign-solver convergence requirement — were all real.

| severity | found | fixed | recommendation | refuted |
|---|---|---|---|---|
| CORRECTNESS | 8 | 5 | 3 | 0 |
| SILENT-ACCEPTANCE | 11 | 11 | 0 | 0 |
| ROBUSTNESS | 16 | 8 | 8 | 0 |
| STYLE / EFFICIENCY | 6 | 2 | 4 | 0 |

Nothing was refuted outright. One finding was **empirically bounded** rather
than accepted at face value — see the witness note below.

## Fixed (commits b92db9b, aad6019, 212a15f, 14f9c7e)

**Silent acceptance — all eleven.**

1. The front door could certify one phase and report another (`round(1/|delta|)`
   does not reproduce an arbitrary phase: 5/12 becomes 2 and walks a half
   period). The walker now takes EXPLICIT targets, and the front door asserts
   the certified phase matches the request instead of overwriting it.
2. Export from a located CERTIFICATE, not a finite summary entry; `t_f` read
   from the certificate's own `z(8)` so the two representations cannot
   disagree; `ok = true` with `z(8) = NaN` no longer passes (and no longer
   displaces a good entry, since `existingTF <= NaN` is false).
3. The spine went to the caller's `sD0`, which built the grid, so the lookup
   selected row 1 by construction and any sheet landed at departure zero. The
   grid ORIGIN and the CERTIFIED phase are now separate things.
4. Packaging took engine, orbits and phases from fresh option defaults. A
   versioned **problem identity** is now stamped at setup, stored on the
   sheet, and used at packaging; options may ASSERT a value, and disagreement
   is an error.
5. Ribs could be generated at a different operating point from the sheet they
   joined. `build_ribs` reconstructs from the sheet's identity and asserts it;
   every saved rib carries it; the packager checks it.
6. Malformed verdicts passed plain `if` tests. `scalar_verdict` (18 checks): a
   verdict is a real finite SCALAR or it is not a verdict. `[]` made `v ~= 1`
   false; `Inf > 0` read as a satisfied gate; NaN tripped no bound.
7. The witness need not have converged — see below.
8. A returned flight array was assumed to be a completed flight. Now: positive
   `t_f`, `t_end` actually reaching it, finite final state, and the all-burn
   mass law with positive final mass.
9. Certification policy applied to crossings but not to library seeds. One
   policy now serves both.
10. Seeds were accepted on `C.ok` alone; and a REFUSED candidate could dedupe
    a later SUCCESSFUL one, dropping the certificate that would have filled
    the cell. Merges now require matching `ok`.
11. A persistent staging directory plus a wildcard could ship stale sheets.
    Staging is now a fresh directory.

**Correctness — five.**

12. `dR/dsA` differentiated the CR3BP field while the residual's target is an
    INTERPOLANT. Different functions. Now the piecewise polynomial's own
    coefficients are differentiated, periodic where the toolbox allows.
    Measured: 7.4e-6 → 9.0e-13 (Richardson, so the check tests the derivative
    and not the finite difference's truncation). **My own test had reported
    7.4e-6 and I recorded it without asking why it was not machine precision.**
13. A level sitting exactly on a turning point was recorded ZERO times. Endpoint
    ownership is now explicit.
14. "Fold localized" meant only that some corrected point lay on the curve.
    Now the EVENT must be resolved, with a scale-aware tolerance — forced by
    measurement: at one cubic fold `|tau_q|` bottoms out at 4.1e-6, exactly
    `sigma_min(R_x)`, while the fold POSITION is exact to 1e-9.
15. A fold now requires the augmented matrix to be numerically REGULAR with a
    real gap, not merely a larger `sigma_min` than `R_x`.
16. The fixed-q Newton had no locality guard, so both halves of a split could
    return the same root and lose the second crossing.

**Robustness — eight.** Reverse continuation defaulted to `qStop = [inf inf]`
and stopped at step 0, so a reverse arc silently did nothing. Walking from a
sweep-library entry hit an unsupported file layout — an explicit root is now
accepted. The last converged crossing was taken even if it failed
certification while an earlier one would have passed; all are now tried,
fastest first. Chart layout, `Dx` positivity and step ordering are checked
rather than assumed.

## The witness finding, bounded by measurement

Astra: agreement between the foreign solver's output and our input does not
prove an independent solve; a solver returning its initial guess on stagnation
gives `dz = 0` and passes. **Correct in general.** Measured for this solver on
2026-09-10: perturbing `lam0` by 1.5x, 3x and 10x moved its answer by 9.3, 37
and 4612, and an impossible target moved it by 5.6. It does not return its
input. So our `dz = 0` results are genuine agreement, and the shipped entries
stand. The hole is closed anyway: the WITNESS's own solution is now flown and
must itself reach the target, which is a stronger statement than two matching
vectors.

## Recommendations, not taken (with reasons)

- **Continuous-time sufficiency.** Sampled minima do not prove positivity
  between samples, and a thresholded rank does not prove `dim S = 1`. Correct.
  The right response is language, not code: the catalog says *numerically
  certified under a stated policy*, and `doc/mintime_second_order_audit.tex`
  already states what a rigorous version would need. Interval arithmetic on
  this problem is a research project, not a patch.
- **"Minimum time" must not imply a global optimum.** Agreed and now stated
  wherever the sheet is described: these are fixed-phase local candidates and
  "fastest found". Phase transversality is NOT enforced — continuation in a
  phase does not impose it — and that is now written down.
- **Crossing enumeration is resolution-dependent.** Two turns inside one step
  can leave endpoint tangent signs unchanged. Documented as such rather than
  claimed exhaustive; curvature-based subdivision is the fix if we ever need
  the guarantee.
- **Residual row scaling.** Only the unknowns are scaled, not the residual
  rows. A real gap; deferred because changing the merit function invalidates
  every step-control constant currently validated by test.
- **End-to-end wall budgets** and fencing the continuation's own residual
  calls. The per-call fence is in; a budget spanning a whole arc is not. This
  is the next real robustness item.
- **`parfeval` cancellation as a demonstrated hard kill.** Trusted on the
  strength of `run_capped`'s own history, not independently proven here.
