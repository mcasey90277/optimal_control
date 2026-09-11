# Certification discipline — how not to repeat 2026-09-09/10

Three external reviews in two days found the same class of defect three
times. Not three different mistakes: **one habit, three instances.**

| date | the defect | the habit |
|---|---|---|
| 09-09 | the phase sweep COMPUTED the conjugate verdict and the witness agreement, then ignored both | computed ≠ enforced |
| 09-10 | `certify_root` passed `wallSec` and I called it a fence; it is an in-process check that cannot bound one grinding evaluation | assumed ≠ checked |
| 09-10 | `dR/dsA` differentiated the CR3BP field while the residual used a spline; my own test reported 7.4e-6 agreement and I wrote the number down | agreed-ish ≠ correct |

Each cost real time: twelve hours of hung compute, a lost rib, and a catalog
that had to be audited after shipping rather than before.

## The five rules

**1. A gate is enforced where it is computed, and a test proves it refuses.**
Every gate in `certify_root` now has a test that feeds it something bad and
requires a NAMED refusal. A gate with no refusal test is decoration: the
2026-09-09 sweep's gates all "passed" because nothing ever consulted them.

**2. When a numerical check agrees to less than machine precision, ask why.**
7.4e-6 between an analytic derivative and a finite difference is not
round-off; it is a different function. The number was on screen for a day
before an outside reviewer asked the question I should have asked. If a
check cannot reach machine precision, either explain the floor (truncation,
conditioning, rank resolution) or find the bug.

**3. Identity travels with the data; it is never re-derived at a boundary.**
Engine, orbits and phases now ride from setup to sheet to catalog as a
versioned `problem` struct. Options may ASSERT a value; disagreement is an
error. Re-deriving at the boundary is how a real trajectory acquires the
wrong Isp, and it was harmless here only because the defaults happened to
match — luck, not design.

**4. Every external call is fenced by something that can kill it.**
An in-process wall check fires between iterations and cannot bound a single
evaluation. Only `run_capped`'s cancellation can. Exceeding a cap is a NAMED
failure, never a silent pass — a candidate whose witness could not be run is
not a candidate that passed.

**5. Audit the ARTIFACT, not only the builder.**
Fixing a builder says nothing about a file already shipped.
`audit_phase_catalog` re-derives every entry from the catalog's own keys and
flies it, exactly as a recipient would. **A deliverable does not ship until
its audit is clean**, and `build_dro_deliverable` enforces that rather than
trusting anyone to remember.

## Two habits of language

- Say **"numerically certified under a stated policy"**, never "certified"
  unqualified. Sampled minima do not prove positivity between samples and a
  thresholded rank does not prove a dimension.
- Say **"fastest found"**, never "minimum time", for a sheet value. These are
  fixed-phase local candidates; phase transversality is not imposed, and
  continuation in a phase does not impose it.

## When to get an outside review

**Before shipping, not after.** The chain review cost $1.49 and 8 minutes and
found 41 findings including three release blockers, on code that already had
nine green test suites. Green tests prove the code does what I thought;
they cannot tell me what I failed to think of.

## Two rules added 2026-09-10 (FINDINGS 41)

6. **A gate that is computed is a gate that is enforced.** `h6Ok` was
   returned by the gates, stored on every entry, and read by nothing; an
   entry could certify with it false. If a quantity is worth computing in
   the gate stack it is worth failing on, with a named reason and a margin --
   and its absence is a failure, never a pass.

7. **A check that cannot fail is not a check.** Twice now a numerical
   "verification" evaluated a formula against the formula's own minimiser
   (N5 with `|alpha| = 1`, N6 with the sphere sample). The test of a check is
   its mutation: name the bug class it exists for, inject that bug, and watch
   it fail. `pmp_pointwise_checks` ships with a wrong-sign vector field in its
   test for exactly this reason.
