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

- **H1 (primary): transcription dual-map error.** The discrete adjoint of the
  9-state cScale-augmented trapezoid differs from the plain 8-state one; the
  correct nodal costate is some κ/cScale/mesh-weighted combination of interval
  duals (the ms_band campaign's situation exactly — resolved there by deriving
  the "midpoint-principled map, mode 'd'").
- **H2 (secondary): extraction/scaling artifact.** CasADi `opti.dual`
  conventions or IPOPT's internal constraint scaling leaking into the returned
  multipliers (would show up as a nonzero explicit stationarity residual).

## 3. Approach — four steps, two of them decisive experiments

**Step 1 — Minimal reproducer + explicit KKT audit.** Tiny instance of the SAME
transcription (N≈40, ~1 rev, fixed-t_f, eps=0, fixed terminal). Dump the full
constraint Jacobian, objective gradient, duals, and bound multipliers; evaluate
the α-stationarity residual ∂L/∂α_k numerically at the converged point.
- Residual ≈ 0 with primer ≠ 0 ⇒ H1: read the correct dual→costate map off the
  stationarity equations symbolically (write out ∂L/∂v_k, ∂L/∂α_k including the
  cScale row and ∂κ/∂r terms).
- Residual ≢ 0 ⇒ H2: compare `opti.dual` vs low-level `nlpsol` `lam_g`; retest
  with `nlp_scaling_method='none'`; chase the convention.

**Step 2 — The 8-state A/B (decisive, parallel to Step 1).** For a CONVERGED
fixed-t_f solution, cScale is a known constant c*. Re-pose with cScale
eliminated (absorb c* into the clock: τ_f₀' = c*·τ_f₀; 8 states, tulip-style
transcription), warm-start AT the solution, solve (≈0 iterations), extract duals
from the reduced problem. Primer collapses to ~0.1° ⇒ cScale row confirmed as
the culprit AND this re-solve IS the practical fix (a dual-extraction front-end,
the exact analog of the tulip campaign's "ε=0 re-solve to regenerate duals").

**Step 3 — Institutionalize.** Wire the verified map (Step 1 route) or the
8-state extraction re-solve (Step 2 route) into `verify_pmp_2body` as its
front-end. Re-verify M0 / M1 / M2 / M2-N1200 from their banked .mats.

**Step 4 — Payoff: Fig 16 analog.** With trustworthy costates, plot ψ(t) =
1 − (T/(Isp g₀))p_m − (T/m)‖Bᵀp‖-equivalent and ‖λ_v‖(t) against the throttle
(our Cartesian switching function S = 1 − β·W with the now-correct costates),
verifying pinpoint zeros at switches and no singular arcs — the paper's H1/H2.

## 4. Module changes

| module | change |
|---|---|
| `results/dual_anomaly/` | + `diag_kkt_audit.m` (Step 1 reproducer + residual dump) |
| `casadi_lt_2body.m` | + option `pinCScale` (fix cScale ≡ given constant, drop row-9 freedom; used only by the extraction re-solve) — OR a sibling `extract_duals_8state.m` if touching the core is riskier |
| `verify_pmp_2body.m` | + dual-extraction front-end (mode option: 'raw' legacy / 'reduced' new default); gates unchanged, NEVER weakened |
| `fig_switching.m` | NEW: the Fig-16 analog figure from a verified-dual solution |
| `README.md` / `DESIGN.md` | finding status update once resolved |

## 5. Milestones & gates

- **D1** — KKT audit: stationarity residual quantified on the reproducer;
  H1-vs-H2 verdict with numbers. *(Gate: an unambiguous branch decision.)*
- **D2** — 8-state A/B verdict: primer from reduced-transcription duals on M1.
  *(Gate: <1° ⇒ H1 confirmed + fix in hand; else Step 1's derived map must do it.)*
- **D3** — Verifier green across milestones: M0/M1/M2/M2-N1200 with primer <1°,
  burn & coast sign ≥99%, relative transversality ≤1e-3, β-spread reported.
  *(Gate: ver.pass=1 on at least M1 and M2; any residual failure re-characterized
  honestly, not gate-weakened.)*
- **D4** — Fig-16 analog produced (`results/fig_switching.png`): ψ zero-crossings
  align with throttle switches (nodeTol ≤1), no finite singular intervals.

## 6. Verification & honesty rules

- Gates are the tulip campaign's; they may not be loosened to declare success.
- If BOTH hypotheses fail (Step 1 residual ≈ 0, Step 2 primer still off, no
  derivable map): STOP, document (the finding graduates from "suspected map
  error" to "genuine open problem"), and escalate to a full symbolic discrete-
  adjoint derivation as its own future campaign (ms_band precedent). Time-box:
  do not exceed ~2 days total without a controller decision.
- Every experiment writes a script + artifact under `results/dual_anomaly/`.

## 7. Estimate

Steps 1+2: ~half a day (no long solves — reproducer is tiny, extraction
re-solves are ~0-iteration). Step 3: ~half a day. Step 4: ~2 h.
Total ≈ 1.5–2 days, low compute.

## 8. References

- Task 12 report + Fix sections (`.superpowers/sdd/task-12-report.md`).
- `results/dual_anomaly/diag_*.m` (the five refutation experiments).
- ms_band campaign: `NLP_lowThrust_GTO_tulip/ms_band/MS_BAND_CAMPAIGN.md`
  (adjudicated midpoint dual map — the method precedent).
- Paper Fig 16 (H1/H2 verification target).
