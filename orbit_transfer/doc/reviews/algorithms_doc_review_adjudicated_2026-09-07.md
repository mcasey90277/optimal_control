# Algorithm document — GPT-6 Astra review, adjudicated (2026-09-07)

Reviewer: GPT-6 Astra via the raw OpenRouter API (185 s, $0.60; raw output in
`algorithms_doc_astra_raw_2026-09-07.md`). Document: `doc/algorithms_orbit_transfer.tex`
(review draft of 2026-09-07, 16 pp). Fourteen findings; the host adjudicated each
against the code before editing. Verdict: 12 accepted (4 P0 + 8 P1), 1 accepted
with a correction to the reviewer's premise (#5), 1 partly declined (#10).

| # | Astra | verdict | what changed |
|---|---|---|---|
| 1 | P0 mass-costate sign wrong | **ACCEPTED — a real error in the document** (the code is AD from H and was never wrong). lam_m' = −sT\|lam_v\|/m² ≤ 0, so lam_m ≥ 0 ahead of its terminal zero. With the wrong sign the closed-form Q̇ would not have cancelled. | eq. `costate` corrected, remark added |
| 2 | P0 min-time normalization undeclared (Mayer φ=t_f gives H(t_f)+1=0, the text uses H=0) | ACCEPTED | declared before eq. `H`: running-cost form L≡1 from there on; admissible domain m>0, no dry-mass bound stated |
| 3 | P0 s=pQ negative for Q<0; huber set-valued at Q=1; huberc does NOT tend to bang-bang at fixed δ; ELI5 "never in between" ignores singular arcs | ACCEPTED on all four (the code clips at 0; the document omitted it) | max(pQ,0); [p,1] at Q=1; "eps and huber → bang-bang as p→0, huberc needs δ→0 too" in text and caption; ELI5 caveat |
| 4 | P0 breakpoints σ_k used as segment durations | ACCEPTED | β_k = τ_{k+1} − τ_k throughout eq. `msres`, the Jacobian blocks (∂R_k/∂t_f = β_k F, last row Dg Φ_K and Dg β_K F) and Algorithm 1 |
| 5 | P0 a sign change in the last bracket is not an endpoint conjugate point; "weak minimum" unjustified | ACCEPTED with a correction to the premise: the instrument never PASSES an ENDPOINT verdict (pass = strcmp(verdict,'PASS')), so no catalog verdict was inflated; but the "weak, non-strict minimum" reading in the document AND in the instrument's header was unjustified. 2 of 18,360 catalog entries carry the flag. | document: ENDPOINT = inconclusive, refine the bracket; `ms_conjugate_test.m` header comment corrected (no behaviour change); ELI5 rewritten around an infinitesimal focusing variation |
| 6 | P1 lift-space "iff" too strong — S is the space of STATIONARY lifts | ACCEPTED | dim S = 1 ⇒ no abnormal lift (sufficient, conservative); C assembly written out |
| 7 | P1 P ∈ R^{6×5} spanning lam_rv(0)^⊥, det vs matrix, hypotheses local, n^T F^- ≠ 0, fixed-t_f regularity | ACCEPTED | all stated; coast skip labelled a diagnostic convention |
| 8 | P1 rank rule is a heuristic; cap can undercut the residual term | ACCEPTED (the code already flags nullResid > 1e-4 as unresolved, never counts it as dim 1) | stated as a residual-informed heuristic; unresolved case reported |
| 9 | P1 "identifiability" conflates conditioning with uniqueness; reflight from the first junction contradicts the junction story | ACCEPTED | reworded everywhere: junctions preserve a numerically usable representation; the 18 fuel entries are short enough (≤30 d) for lam_0-only reflight, the 40-rev flagship is not |
| 10 | P1 binding omits δ; reconciliation of re-solve vs packaged junctions unstated | PARTLY DECLINED: the builder looks up by source+cell+γ and then requires the re-solved initial costate to agree with the packaged one to 1e-8 relative — that agreement pins (p, δ, family) far more tightly than a key would (different (p,δ) move lam_0 by ≫1e-8). ACCEPTED that the document must say so. | binding paragraph rewritten: lookup keys, 1e-8 agreement, rows 1:7, converged re-solve, packaged junctions kept |
| 11 | P1 "certify"/PASS wording slides from assessed to proven; tfMin |Δz| alone could mean "did not run" | ACCEPTED | PASS-row rewritten in numerical language; acceptance = tfMin convergence AND |Δz|<1e-6 |
| 12 | P1 overclaims: MS "well-conditioned", huberc "cure", switching-time = certificate, m_f 1e-6 as error bound, legacy counts | ACCEPTED on all five | softened as proposed; conjugate counts labelled legacy-instrument |
| 13 | P1 direct transcription / dual map / endpoints / K / units / g(r) / t_f^min are black boxes | ACCEPTED | new §2.1 "Conventions a reimplementer needs" (units, g(r), phase→endpoint, HS separated form with lifted t_f copies, dual map, K ladder, residual and reflight definitions) |
| 14 | P1 figures: tfMin lacks inputs; "every entry passed direct" wrong for continuation rungs; verify figure implies universal level-2 pass; ms_minenergy has junction unknowns; fuel output is a smoothed candidate | ACCEPTED | fig 1 input edge + caption; fig 3 caption (61 refuted / 50 unverified); fig 4 labels |

What Astra got wrong: nothing material. Its one premise error (#5, that the
instrument called ENDPOINT a minimum) was the document's fault — the header
comment said exactly that.

Lessons: (i) a sign error survived the document's own derivation because the
closed-form Q̇ was quoted from the code, not re-derived from the printed
equations — re-derive every quoted identity from the printed ones; (ii) the
smoothing families' clip at 0 and huberc's δ→0 requirement were in the code
and in FINDINGS §26 but not in the document — an exposition that omits a
branch of the law is wrong, not merely incomplete.
