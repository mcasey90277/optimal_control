# Research memo: closing the three open gaps in the second-order certificate

**Status: findings and proposals, nothing implemented yet.** Written after the
2026-09-10 Astra review of `transfer_study` named three gaps: the conjugate
test's sampling, the rank threshold behind `dim S`, and whether the reduced
6x6 determinant is provably the right Jacobi problem.

---

## Gap 3 first, because it is mostly already closed

Astra: *"The equivalence of this reduced 6x6 test to the required
mixed-endpoint Jacobi problem is not established."*

**It is established, in `doc/mintime_second_order_audit.tex`, which was not in
the review bundle.** I sent code only. Two results there answer it:

- **The endpoint is not mixed.** For an all-burn arc `m(t) = m0 - Tt/c`
  whatever the control does, so the free mass direction is never exercised by
  an admissible variation: `dm(tf) = -(T/c) dtf` is *determined* by `dtf`. The
  competitor set is identical to that of a 6-state **point-target** problem, so
  the point-target free-time theory applies and no focal construction is needed.
  (This is exactly the opposite of the fixed-time min-fuel line, where the
  throttle varies and the mass really is free.)
- **Proposition 1** proves the instrument's determinant equals the BCT
  immersion test: the column spaces coincide, so the ranks do, so the sign
  changes are the same events.

**Action: send the document, not the code.** The proof deserves an adversarial
read by someone who has it in front of them. One caveat I want checked: the
reduced 6-state system is non-autonomous through `m(t)`, and the free-time
construction's flow column is built from the 6-state field. The audit asserts
this is immaterial; that assertion is the thing to attack.

---

## Gap 1: what the sampling misses — now measured, not speculated

I rebuilt the state-transition block **by finite differences of the flow in
`lambda_0`** (an independent construction from the variational equations the
instrument integrates) and evaluated the quotiented 6x6 at 400 dense samples
along the certified anchor, recording the whole singular-value spectrum rather
than only the determinant.

**Result 1 — the early-arc skip is vindicated, quantitatively.** Sixteen of
eighteen determinant sign changes lie in the first 6.5% of the arc, where the
block is structurally near-zero and `|det|` is 1e-8..1e-6 while `sigma_min`
oscillates over three decades. These are noise, and they are exactly what the
instrument's `firstFullRank` skip exists to discard. Dense sampling *without*
that skip would report sixteen spurious conjugate points.

**Result 2 — and something the junction sampling cannot see.** In the last
0.5% of the arc `sigma_min` collapses from a median of 1.7e-3 to **6.0e-8**,
with two determinant sign changes at `t/tf` = 0.9950 and 0.9975.

| region | median sigma_min | min sigma_min | det sign changes |
|---|---|---|---|
| whole arc | 1.7e-3 | 6.0e-8 | 18 |
| `t/tf` > 0.05 | 1.7e-3 | 6.0e-8 | 6 |
| last 0.5% | -- | 6.0e-8 | 2 |

The 24 junctions are spaced about 4% apart, so this structure falls between the
last two samples. **The instrument reports PASS with `min|det|` = 2.4e-4; the
dense scan says the matrix comes within 6e-8 of singular before `tf`.**

> ## ADJUDICATED, same day: it is NOT a conjugate time, and no shipped claim
> ## is affected.
>
> The discriminator was a POSITIVE CONTROL -- run the same endpoint spectrum on
> entries the conjugate test REFUTED and see whether it separates them from a
> certified one. It does not:
>
> | entry | verdict | sigma_min at t_f | last gap / median gap |
> |---|---|---|---|
> | anchor, 17.80 d | **certified** | 1.37e-07 | 3.76 |
> | 22.05 d, sA 0.0754 | refuted | 1.70e-08 | 34.9 |
> | 18.92 d, sA 0.1587 | refuted | **2.17e-04** | 0.91 |
> | 31.51 d, sA 0.2421 | refuted | 1.39e-07 | 0.88 |
>
> Two facts settle it. A refuted entry (31.51 d) has the SAME endpoint
> sigma_min as the certified one, 1.39e-7 against 1.37e-7 -- so the endpoint
> value carries no information about the verdict. And a refuted entry
> (18.92 d) has a perfectly healthy endpoint spectrum, sigma_min 2.2e-4 with
> no anomalous gap -- so its conjugate point is INTERIOR, where the instrument
> does look and did find it.
>
> The endpoint smallness is therefore a property of the transition matrix at
> t_f shared by certified and refuted arcs alike: the geometric grading of a
> hyperbolic flow's STM (sigma = 2.45, 5.6e-2, 1.3e-2, 3.9e-4, 1.8e-5, 1.4e-7
> -- each roughly 30x the next, no single collapsed direction). It is the same
> phenomenon as the 4.4e9 shooting conditioning, and it is exactly why the
> instrument treats an endpoint zero as INCONCLUSIVE rather than a refutation.
>
> Two method notes worth keeping. The finite-difference construction that
> raised the alarm is accurate to only ~2e-6 relative at h = 1e-6 (verified by
> convergence against the integrated STM at four step sizes, clean h^2), which
> on some entries is ABOVE the sigma_min being measured -- the original scan
> was reading its own noise on those. And the determinant near t_f is the
> product of six graded singular values, about 1e-16 here, so its SIGN is
> meaningless there whatever the geometry does.

**The original two candidate explanations, retained for the record:**

1. **A conjugate time at or just before `tf`**, which would violate the BCT
   hypothesis and invalidate the minimality claim for this entry.
2. **The known endpoint conditioning.** The shooting Jacobian at this operating
   point has `cond(J)` = 4.4e9, i.e. a smallest-to-largest singular value ratio
   near 2e-10; the dense scan's worst ratio is 2.5e-8. Same story, and a
   nearly-singular endpoint constraint is expected to show here.

There is a decisive discriminator: **`tf` is a conjugate time if and only if
the shooting Jacobian is singular.** Ours converged with a finite condition
number, which argues for explanation 2 -- but that argument must be made
carefully, since 4.4e9 is not small.

### Proposals for gap 1, in order of cost

**A. Monitor the singular-value SPECTRUM densely, not the determinant.**
The determinant conflates magnitude with conjugacy: it is 7e-19 where the
matrix is merely small, and its sign is meaningless in the early arc. The
spectrum separates them. It also sees **multiplicity**: a conjugate time of
multiplicity k sends k singular values to zero together, which is precisely
the even-order case a sign test is blind to. Here `sigma_5` never drops below
1.7e-5, so no multiplicity-2 crossing occurred on this arc -- a statement the
determinant cannot make at all. Cheap: the variational integration is already
dense, only the evaluation is discarded.

**B. Cross-check with an independent instrument: the Morse index.** Discretize
the second variation on the existing mesh and count its negative eigenvalues.
By the Morse index theorem the count equals the number of conjugate points with
multiplicity, so zero negative eigenvalues means no conjugate point. This is a
sparse symmetric inertia computation, shares no code path with the
determinant, and sees multiplicity by construction. Two instruments agreeing is
worth more than one instrument sampled more finely.

**C. Rigour, if we ever need it: validated integration.** Enclose the
variational equations with interval or Taylor-model arithmetic and obtain a
sign-definite enclosure of the determinant, or a positive lower bound on
`sigma_min`. Only this can *prove* absence. Tools: INTLAB (MATLAB), CAPD
(C++), Julia's Taylor models. This is a project, not a patch, and I would not
start it before A and B.

---

## Gap 2: the rank threshold behind `dim S`

Astra: *"The rank threshold is heuristic, and its stated safety guarantee is
false."* Correct on both counts. The current rule caps the tolerance at
`1e-3*sigma_1` and the header claims this prevents a noisy lift from reporting
a spurious nullity; it prevents no such thing.

**Proposal D: replace the verdict with a certified margin.** The mathematics is
already available and is a theorem rather than a heuristic:

- `dim S >= 1` is **constructive** -- we exhibit the lift, so nothing is
  assumed.
- `dim S <= 1` follows from Eckart-Young: the distance from the constraint
  matrix to the nearest rank-deficient one is `sigma_6`. So if
  `sigma_6 > ||dC||`, where `||dC||` bounds the numerical error in the matrix,
  then `rank = 6` **exactly** and `dim S = 1` is certified.

What is missing is `||dC||`, and it is obtainable: integrate the adjoint at two
tolerances and take the difference as a measured error estimate, or bound it
rigorously with the validated integration of proposal C. Then report the
**margin** `sigma_6 / ||dC||` per entry instead of a pass or fail against an
arbitrary number. On the anchor the reported values are `nullResid` 2.1e-6 and
a singular-value gap of 3.4e-7, so the margin is currently unknown rather than
comfortable -- which is the point.

---

## What I would do next, in order

1. ~~**Adjudicate the endpoint collapse.**~~ **DONE, same day** -- see the box
   above. It is conditioning, not conjugacy, established by a positive control
   rather than by the argument originally proposed: a REFUTED entry has the
   same endpoint sigma_min as the certified one (1.39e-7 vs 1.37e-7), and
   another refuted entry has a healthy endpoint spectrum with an interior
   conjugate point. No shipped claim is affected.
2. **Implement A**, dense spectrum monitoring, and re-sweep. Cheap, and it
   converts the sampling caveat from "unknown" to "measured".
3. **Implement B**, the Morse-index cross-check, on the golden cells.
4. **Send the audit document** for the adversarial read that gap 3 deserves.
5. Leave C and the rigorous half of D scoped but unstarted.
