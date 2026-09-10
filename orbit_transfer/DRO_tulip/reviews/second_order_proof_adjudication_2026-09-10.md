# Adjudication — Astra review of the second-order proof, 2026-09-10

Review: `second_order_audit_proof_review_astra_2026-09-10.md` (GPT-6 Astra,
xhigh, 18,128 reasoning tokens, $1.19). This one reviewed the MATHEMATICS in
`doc/mintime_second_order_audit.tex`, not code — the earlier chain review had
called the equivalence "asserted, not established" only because I sent it code
without the document.

**Headline: the document is weaker than I believed, the code is not affected,
and one of our gates is safe in the direction we actually use it.**

| verdict | count | what it touches |
|---|---|---|
| WRONG | 4 | one is a document typo; three are real derivation errors |
| GAP | 4 | scope and completeness of the theorem transfer |
| IMPRECISE | 1 | wording |
| FINE | 3 | the Lemma, the rank inference, the H = 0 handling |

## What is NOT affected

**No catalog number changes.** The audit of the shipped 53 entries verified
each against its own claims and passed 53/53; nothing here touches that.

**The mass-adjoint sign error is in the DOCUMENT ONLY** — and it was checked,
not assumed. The document printed `lam_m_dot = +T|lam_v|/m^2`; the canonical
equation gives the negative sign. The code carries the negative form, matching
to ten digits, and `max|H|` along the arc is 3.3e-8, which it could not be if
the sign were wrong. Corrected in the document. Worth recording anyway: this
is the document I reason from, and the earlier algorithm document had the same
error fixed weeks ago without this one being checked.

## What must be repaired, and how bad each is

**1. Proposition 1's conclusion does not follow from its own premise.
(Serious, but the instrument probably survives.)**

Astra exhibits a PERMANENT kernel direction of the construction as the
document states it: with `b = -T/c`,

    w = lambda(0) + (1/b) e_m,     dH_0(w) = 0,     Phi_xlambda(t,0) w = 0.

**Our own two lemmas prove this**, which is why I accept it without a
numerical check: `Phi_xlambda * lambda(0) = 0` is the scaling lemma, the
`lam_m` column is identically zero, so any multiple of `e_m` added to
`lambda(0)` stays in the kernel; and `dH_0(w) = f.lambda(0) + (1/b) f_m =
-1 + 1 = 0`. So the 7-dimensional map is non-immersive at EVERY time, not
merely at conjugate times, and "rank drop = conjugate time" cannot be read off
it.

The instrument quotients exactly this away, so it is almost certainly
computing the right object — but the DERIVATION has to be rebuilt from the
effective 6-dimensional domain rather than asserted on a 7-dimensional one.
Encouragingly, Astra rates the Lemma and the rank inference **FINE**, so the
hard half is already sound; what is missing is the correct identification of
the domain and its terminal transversality (`lam_m(t_f) = 0`, which is exactly
what removes the extra freedom).

**2. Scope: all-burn competitors are not all competitors. (Real, and it
narrows the claim.)** The stated problem allows throttle in [0,1]; the mass
reduction covers all-burn trajectories only, where
`delta m_f = -k delta t_f - k \int delta s`. Strict `Q > 0` does not make the
throttled competitors vanish — it can *support* a strict-bang argument, which
the document does not make. So what we certify is minimality among **all-burn**
trajectories with the same endpoints. Extending to all admissible controls
needs the bang argument written down.

**3. Normality: the justification is wrong, the GATE IS SAFE.** The document
says an abnormal lift "would require lambda_v == 0". False — abnormal means
the *cost* multiplier vanishes, and an abnormal extremal can have
`|lambda_v| > 0`. The claimed "no abnormal lift **iff** dim S = 1" is not
established either, because S imposes only `lambda_v || alpha` and omits the
PMP inequalities. **But the direction we use is the safe one**: Astra states
that `dim S = 1` with a normal element of nonzero functional value is a
*sufficient algebraic exclusion* of abnormal lifts. Our gate excludes; it does
not claim the converse. Fix the wording, keep the gate.

**4. The projected determinant can acquire structural zeros.** Homogeneity of
degree one in the multiplier gives `p(t)' J(t) = 0`, so the projected object
has degeneracies of its own that are not losses of optimality. This needs
working through against our actual quotient before I can say whether it bites.

## Suggestions

1. **Rewrite the derivation** starting from the permanent kernel: state it,
   quotient it explicitly, define the effective 6-dimensional domain with
   `lam_m(t_f) = 0` included, then apply the (already-endorsed) image lemma.
   My expectation is that the same 6x6 determinant falls out, which would
   vindicate the instrument and fix the proof at once. That expectation is
   what to test, not to assume.
2. **Restate the claim's scope**: minimality among ALL-BURN competitors, until
   the strict-bang argument is written.
3. **Reword the normality subsection** as a sufficient exclusion, and drop the
   `lambda_v == 0` sentence entirely.
4. **Work finding 4 through** against our quotient.
5. Only then send it back for a second read.

## The lesson worth keeping

The earlier review said this equivalence was asserted rather than established.
I answered that it was proved in a document the reviewer had not been given.
**Both were true**: the document exists, and its central step does not follow.
Having a proof written down is not the same as having checked it, and I had
been citing this one for weeks.
