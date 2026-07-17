# DESIGN — Campaign B: Resolve the PMP dual/primer anomaly

**Goal.** Turn the open finding — a reproducible 10–24° primer misalignment in the
KKT duals of `casadi_lt_2body` solutions — into a *resolved* mechanism with a
verified dual-to-costate map, restoring campaign-grade first-order PMP
certification (primer <1°, switching-sign ≥99%, relative transversality ≤1e-3)
and unlocking the paper's Fig 16 analog (switching function ψ and ‖Bᵀp‖ plots,
their assumptions H1/H2).

Date: 2026-07-17. Status: approved design, pre-implementation.
Prerequisite for: Campaign A (thrust ladder) inherits whatever verification
front-end this campaign produces. Recommended to run FIRST (cheaper, no long
solves, hardens the tooling A will use).

---

## 1. Evidence base (all preserved; do not re-run to believe)

From Task 12 + five controller refutation experiments (scripts in
`results/dual_anomaly/`, artifacts `results/M1_dualpolish.mat`,
`results/M1_wideal_probe.mat`):

| experiment | result |
|---|---|
| raw primer from `lamDef(4:6,:)` | M0 13.3°, M1 20.5°, M2 18.8°, M2-N1200 23.5° (bulk, not outliers — M1 median 17.3°) |
| tight ε=0 dual-polish re-solve at the solution | primer IDENTICAL to 4 decimals (20.4645°) — duals converged & deterministic |
| node-centered averaging of adjacent interval duals | no effect (17.26°→17.22°) |
| α box bounds widened ±1.01→±2 | identical (20.4645°) |
| ‖λ_v‖ magnitudes on burns | healthy (30–58, 83% of total dual magnitude); misalignment uniform across magnitude quartiles |
| **discriminating clue** | M2-N1200 **transversality passes at 9.5e-11** while its primer reads 23.5°; the tulip solver (8 states, no cScale, same naive map) reads **0.06°** |

Reading: the duals are the true KKT multipliers of the NLP; the naive map
(interval dual rows 4:6 ≡ continuous λ_v at the node) is wrong for THIS
transcription. Prime structural difference from the tulip solver: the **cScale
slack state** (row 9) multiplying every dynamics row, plus the ∂κ/∂r chain terms
it drags into the adjoint.

## 2. Hypotheses, ranked

**Superseded 2026-07-17 by the three-way core review** (`doc/reviews/2026-07-17_triage.md`,
GPT-5.6-terra + Gemini 3.1 Pro + host). Both reviewers derived, independently and before
any experiment was run, that the original H1 below is **dead as algebra**:
`cScale`'s defect row `c_{k+1} − c_k` has derivative **zero** w.r.t. α and v, and the
factor `q_k = c_k·κ_k` multiplying every dynamics row is a **common node-local scalar
that cancels out of the direction** (Gemini: *"neither κ, cScale, nor dτ shift the
direction"*). cScale cannot rotate the primer through the stationarity algebra — full
stop. The anomaly therefore lives in the **dual extraction**, not the discrete adjoint
map. Hypotheses below are the review's replacement, not a repair of the original.

- ~~**H1 (primary): transcription dual-map error.**~~ **REFUTED (algebra).** The claim
  was that the discrete adjoint of the 9-state cScale-augmented trapezoid differs from
  the plain 8-state one via a κ/cScale/mesh-weighted combination — the ms_band
  campaign's situation. It does not: cScale's row and the shared `q_k` scalar are both
  provably inert to the α/v stationarity equations. Retained here struck-through as a
  record of what was tested and killed, per house honesty rules (§6).

- **H1 (new) — scaling cross-talk (Gemini's mechanism).** IPOPT's `nlp_scaling_method`
  row-scales the cScale-coupled KKT system; if `opti.dual` un-scales imperfectly, the
  retrieved λ_v is a "rotated projection" of the true costate. cScale survives as a
  suspect here — **not as an algebra term in the direction, but as the structural
  source of the row-scaling asymmetry** (its defect row has a different natural scale
  than rows 1–8). Decisive/cheapest test: **T0** below.

- **H2 (new) — extraction convention/order (terra's mechanism).** What `opti.dual`
  returns is not the multiplier convention or ordering the diagnostic assumes — or, more
  seriously, the returned duals are not a usable KKT certificate at all for this
  transcription. Decisive test: **T1** below (the tangential Lagrangian residual —
  discriminates "wrong convention" from "not a certificate").

Host adjudication (triage): H1(new) and H2(new) are compatible, not competing — Gemini
names a specific corruption channel, terra names the test that tells corruption apart
from a genuinely broken map. Both survive pending T0/T1.

### Evidence table addendum — algebraic refutation of the old H1

| experiment | result |
|---|---|
| **cScale stationarity derivation** (terra + Gemini, independent) | `∂(c_{k+1}−c_k)/∂α_k = ∂(c_{k+1}−c_k)/∂v_k = 0`; `q_k = c_k·κ_k` is a common node-local scalar factor in the dynamics rows and cancels from the α/v-stationarity direction. **cScale cannot rotate the primer as algebra.** Confirms independently that experiment (b) (node-centered averaging, 17.26°→17.22°) already showed the naive interval→node indexing fix is *necessary but not sufficient* — consistent with the corruption being extraction-side, not map-side. |

The §1 evidence table itself is unchanged and now better contextualized: the five
refutation experiments there still stand as the empirical record; this addendum is the
algebraic explanation for why none of them (including the averaging fix) closed the gap.

## 3. Approach — T0 → T1 → T2, cheapest-decisive-first (per triage)

Superseded 2026-07-17: the original Step 1/Step 2 order is replaced by the
triage's ranked test order. T0 is new (falls out of Gemini's scaling mechanism
and nobody had named it); T1 is terra's decisive discriminator, promoted ahead
of the 8-state test because it is cheaper and settles H1(new) vs H2(new)
directly; T2 is the old "Step 2" re-scoped — it is now a *formulation-
equivalence / scaling-isolation* test, not a coupling test (the coupling
hypothesis it was built to test is dead).

**T0 — scaling off (new, cheapest; directly tests H1-new).** Set
`opti.solver('ipopt', ..., struct('nlp_scaling_method','none'))` (currently
unset ⇒ IPOPT default `'gradient-based'`, `casadi_lt_2body.m:132`), re-solve
warm at an already-converged solution, re-check the primer with the existing
naive map. One option line. If Gemini's row-scaling mechanism is right, this
moves the primer immediately; if it doesn't move, H1(new) is weakened (not
killed — un-scaling could still be imperfect in a way this toggle doesn't
touch) and T1 becomes the load-bearing test.

**T1 — terra's tangential Lagrangian residual (decisive; discriminates H2-new
from "duals are not a KKT certificate").** Assemble the FULL NLP Lagrangian
gradient from the returned duals — defect duals (`conDef`), cone multipliers
(‖α‖=1), terminal/equality multipliers, and bound multipliers — and evaluate
the **tangential** component `(I − α_kα_k')·∂L/∂α_k` at interior burn nodes.
Tangential, not raw, because the cone constraint's own multiplier absorbs the
radial component of ∂L/∂α_k by construction; a raw residual would be
nonzero even for a correct certificate.
- Residual machine-small AND primer still rotated ⇒ **H2(new) confirmed**: the
  duals ARE a valid KKT certificate, so the corruption is in *extraction*
  convention/order (which rows/order `opti.dual` hands back) — chase the
  convention against low-level `nlpsol` `lam_g`.
- Residual NOT small ⇒ the returned `opti.dual` set is not assembling into a
  valid Lagrangian stationarity point at all for this transcription — escalate
  per §6 (time-boxed, becomes its own derivation campaign).

**T2 — 8-state elimination (RECAST: equivalence test, not a coupling test).**
For a CONVERGED fixed-t_f solution, cScale is a known constant c*. Re-pose
with cScale eliminated (absorb c* into the clock: τ_f₀' = c*·τ_f₀; 8 states,
tulip-style transcription), warm-start AT the solution, solve (≈0 iterations),
extract duals from the reduced problem. Under the old H1 this was billed as
"cScale row confirmed as the culprit." That framing is retired: cScale cannot
couple through the algebra (§2), so a clean primer out of T2 does NOT indict
cScale-coupling — it isolates whether the 9-state transcription's scaling/
extraction path is the difference-maker versus the algebraically-equivalent
8-state one. **Terra's caveat carried over verbatim in substance: T2 "does not
make a single-interval map correct"** — even if T2 comes back clean, the
step-weighted averaging fix (below) is still required on the 9-state solver,
because T2 sidesteps the interval→node indexing question rather than
resolving it.

**[CORRECTNESS] — mandatory regardless of T0/T1/T2 outcome.** The nodal primer
must use the **step-weighted adjacent-interval average**, one-sided at the
endpoints:

    lam_v,k = (h_{k-1}*Lam_v,k-1 + h_k*Lam_v,k) / (h_{k-1} + h_k),   1 < k < N+1
    lam_v,1 = Lam_v,1;   lam_v,N+1 = Lam_v,N

not interval k's dual assigned directly to node k (the current
`verify_pmp_2body.m:34-43` / `casadi_lt_2body.m:162-166` behavior). On our
**uniform** τ-mesh this collapses to the simple adjacent average, which is
exactly refutation experiment (b) in §1 — and that experiment already showed
the fix is **necessary but NOT sufficient** (17.26°→17.22°, nowhere near the
10–24° gap). Ship the weighted form anyway (it is provably correct for
non-uniform meshes and costs nothing here); do not expect it to close the
anomaly on its own.

**Step 3 — Institutionalize.** Wire whichever of T0/T1/T2 resolves the
mechanism into `verify_pmp_2body` as its front-end (plus the mandatory
weighted-average map above regardless). Re-verify M0 / M1 / M2 / M2-N1200 from
their banked .mats.

**Step 4 — Payoff: Fig 16 analog.** With trustworthy costates, plot ψ(t) =
1 − (T/(Isp g₀))p_m − (T/m)‖Bᵀp‖-equivalent and ‖λ_v‖(t) against the throttle
(our Cartesian switching function S = 1 − β·W with the now-correct costates),
verifying pinpoint zeros at switches and no singular arcs — the paper's H1/H2.

## 4. Module changes

| module | change |
|---|---|
| `casadi_lt_2body.m` | T0: `nlp_scaling_method` option threaded to `opti.solver` call (`:132`) — one line. + option `pinCScale` (fix cScale ≡ given constant, drop row-9 freedom; used only by T2's extraction re-solve) — OR a sibling `extract_duals_8state.m` if touching the core is riskier |
| `results/dual_anomaly/` | + `diag_kkt_audit.m` (T1 reproducer: full Lagrangian assembly from all dual groups + tangential residual dump) |
| `verify_pmp_2body.m` | + the mandatory step-weighted adjacent-average primer map (replaces `:34-43`'s direct interval→node assignment) + dual-extraction front-end (mode option: 'raw' legacy / 'reduced' new default once T0/T1/T2 pick a mechanism); gates unchanged, NEVER weakened |
| `fig_switching.m` | NEW: the Fig-16 analog figure from a verified-dual solution |
| `README.md` / `DESIGN.md` | finding status update once resolved |

## 5. Milestones & gates

- **D0** — T0 verdict: primer re-checked with `nlp_scaling_method='none'` at a
  converged solution. *(Gate: moves ⇒ H1-new strongly supported, go straight to
  institutionalizing a re-solve-with-scaling-off front-end; does not move ⇒ H1-new
  weakened, proceed to D1.)*
- **D1** — T1 verdict: tangential Lagrangian residual `(I−αα')∂L/∂α_k` quantified
  on burn nodes across the full dual set (defect + cone + terminal + bound
  multipliers). *(Gate: an unambiguous H2-new-confirmed vs not-a-certificate
  branch decision, with numbers.)*
- **D2** — T2 verdict: primer from the 8-state-elimination re-solve's duals on
  M1, read as an equivalence/isolation result (NOT a coupling confirmation —
  see §3). *(Gate: clean primer narrows the corruption to the 9-state
  scaling/extraction path; not clean ⇒ the corruption is present even in the
  algebraically-simpler transcription, escalate per §6.)*
- **D3** — Verifier green across milestones: M0/M1/M2/M2-N1200 with primer <1°,
  burn & coast sign ≥99%, relative transversality ≤1e-3, β-spread reported.
  *(Gate: ver.pass=1 on at least M1 and M2; any residual failure re-characterized
  honestly, not gate-weakened.)* Includes the mandatory weighted-average map
  regardless of which of D0/D1/D2 resolved the mechanism.
- **D4** — Fig-16 analog produced (`results/fig_switching.png`): ψ zero-crossings
  align with throttle switches (nodeTol ≤1), no finite singular intervals.

## 6. Verification & honesty rules

- Gates are the tulip campaign's; they may not be loosened to declare success.
- If T0/T1/T2 ALL fail to localize the corruption (scaling toggle doesn't move
  it, tangential residual is small yet primer stays rotated in a way T2 can't
  explain, 8-state re-solve reproduces the same rotation): STOP, document (the
  finding graduates from "suspected extraction artifact" to "genuine open
  problem"), and escalate to a full symbolic discrete-adjoint derivation as its
  own future campaign (ms_band precedent). Time-box: do not exceed ~2 days
  total without a controller decision.
- Every experiment writes a script + artifact under `results/dual_anomaly/`.
- Per triage action 9: do NOT gate the Campaign A MEE validation (`DESIGN_thrust_ladder.md`
  Gate P2a) on the raw-dual primer metric until this campaign closes — the m_f/
  structure gates stand on their own for Campaign A.

## 7. Estimate

T0: minutes (one option line, re-solve at a converged point). T1: ~half a day
(reproducer + full-Lagrangian assembly is new work, not previously scoped).
T2: ~half a day (as before, now framed as isolation not confirmation). Step 3:
~half a day. Step 4: ~2 h. Total ≈ 1.5–2 days, low compute — the reordering
changes emphasis, not the budget.

## 8. References

- Task 12 report + Fix sections (`.superpowers/sdd/task-12-report.md`).
- `results/dual_anomaly/diag_*.m` (the five refutation experiments).
- ms_band campaign: `NLP_lowThrust_GTO_tulip/ms_band/MS_BAND_CAMPAIGN.md`
  (adjudicated midpoint dual map — the method precedent).
- Paper Fig 16 (H1/H2 verification target).
- **`doc/reviews/2026-07-17_triage.md`** (this file's §2/§3 rewrite) + the two
  full reviews alongside it, `2026-07-17_core_review_gpt56terra.md` and
  `2026-07-17_core_review_gemini.md` — the source of the H1 refutation and the
  T0/T1/T2 test order.
