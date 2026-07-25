# Lessons — the dual-extraction anomaly (Campaign B)

**Open:** ~2026-07-16 (earlier in Cartesian form). **Resolved:** 2026-07-25.
**Blast radius:** first-order PMP certification in `earth_elliptic_to_geo` AND
`earth_elliptic_to_geo_CR3BP`; PSR blocked in the CR3BP campaign; a published
technical note carrying the failure as a standing limitation.

Companion to `DEEP_THRUST_LESSONS.md`. That one is about making hard solves
converge; this one is about **not spending nine days on the wrong hypothesis**.

---

## 1. What it actually was

`opti.dual(con)` returns the multiplier for CasADi's **canonicalized**
constraint orientation, not the orientation of the `opti.g` row it pairs with
in `∇f + Aᵀλ = 0`. Minimal demonstration — the same constraint, three writings:

| written as | `opti.dual` stationarity residual | `opti.lam_g` |
|---|---|---|
| `x - c == 0` | 0.0 | 0.0 |
| `c - x == 0` | **10.0** | 0.0 |
| `x == c` | 0.0 | 0.0 |

The collocation defect `X(:,k+1) - X(:,k) - (h/2)·ΔL·(F_k + F_{k+1}) == 0`
canonicalizes **row by row**, so the corruption was entry-wise: magnitudes
exact to 1e-16, ~44–60% of signs inverted.

**Fix:** record the constraint group's `opti.g` row range at build time and
take duals from `opti.lam_g(defRows)`. Correct by construction — those rows
*are* the orientation that pairs with `lam_g` in stationarity.

Result: `verify_pmp_mee` went from FAIL to PASS on **9/9** certified rows —
the whole 10 / 5 / 2.5 / 1 / 0.5 / 0.2 / 0.1 N ladder plus both PSR rows —
and on the CR3BP 10 N and 5 N rows. Every one reads primer median
**0.000°**; sign agreement 100% everywhere except the 1 N PSR row at 99.94%.
Headline deltas: 10 N 32.370° / 78.35% → 0.000° / 100.00% with
switch-alignment error 21.2 → 0.15; 1 N PSR 59.967° → 0.000°; CR3BP 10 N
32.4° / 78.4% → 0.000° / 100.00%. Sweep table:
`results/verify_pmp_all.mat` via `verify/run_verify_pmp_all.m`.

---

## 2. Why it took nine days — the expensive lessons

### L1. The working implementation was already in the repo. Diff siblings before theorizing.

Every CR3BP-family solver — `casadi_minfuel_sundman` (tulip), both ELFO
`casadi_*_freetf` solvers, `psr_second_order` — already extracted duals via
`opti.lam_g` with a row-range reshape. **Only the two earth solvers used
`opti.dual`.** The tulip solver's primer read 0.06° the whole time.

Worse: that 0.06° was *recorded in the design doc as the "discriminating
clue"* and then explained with a structural theory (tulip has 8 states and no
`cScale`; the MEE/Cartesian solvers have an extra slack row). The real
explanation was six lines of extraction code in a sibling folder.

> **Rule:** when one implementation of a shared concept works and another
> doesn't, diff the two implementations **before** building a theory about why
> their *problems* differ. A structural difference between two models is a
> seductive explanation precisely because there is always one available.

### L2. A physics-shaped correlation is not evidence of physics.

The misalignment correlated with orbital eccentricity — genuinely, and
strongly: measured `corr = +0.93` (this was quoted qualitatively for over a
week and only *measured* on the day it was solved). That correlation drove
every hypothesis toward eccentricity-dependent transcription effects.

But a per-row sign flip on a defect whose rows have eccentricity-dependent
magnitudes will *always* produce an eccentricity-correlated error. The
correlation was a downstream shadow of the corruption, not a fingerprint of
its cause.

> **Rule:** a correlation with a physical quantity constrains *where the error
> shows up*, not *what causes it*. Rank hypotheses by mechanism, and check the
> plumbing before the physics.

### L3. The decisive test was specified, ranked first, estimated at half a day — and not built.

The 2026-07-17 three-way review (GPT-5.6-terra + Gemini 3.1 Pro + host) got
this **right**. It killed the original `cScale` hypothesis as algebra, named
H2(new) — "the returned duals are not the convention the diagnostic assumes,
chase it against low-level `nlpsol` `lam_g`" — and specified test T1 (assemble
the full Lagrangian from all dual groups, check the tangential residual) as
the decisive discriminator.

T1 was never built. The campaign ran five refutation experiments around the
symptom instead, and the anomaly stayed open across two campaigns and one
published note. When T1 was finally built it took under an hour and localized
the bug on the first run.

> **Rule:** when a triage names a decisive test, build the decisive test. A
> ranked hypothesis list with the top test unexecuted is not a diagnosis.

### L4. A false structural claim in a header was load-bearing.

`verify_pmp_mee.m`'s header asserted this transcription "has no cone
constraint on beta … `||beta||=1` is enforced by construction, not as an NLP
constraint with its own multiplier." It is not: `casadi_lt_mee` imposes
`beta(1)^2+beta(2)^2+beta(3)^2 == 1` as a real constraint with a real
multiplier. Harmless numerically (that multiplier is purely radial and drops
out of the tangential test), but it made the multiplier inventory look shorter
than it was and reinforced the "nothing left but the transcription" framing.

> **Rule:** header claims about *structure* ("there is no X here") should be
> asserted in code or verified when written. Prose facts rot silently and get
> reused as premises.

### L5. "Certified" was doing double duty.

Certification rested on four primal NLP metrics; the PMP gate failed. Both
statements were true and both were reported honestly — but "certified" in the
campaign's shorthand quietly came to mean "certified except for the thing that
never passes," which lowered the cost of leaving it broken.

> **Rule:** a check that is always red stops being information. Either fix it,
> or move it out of the gate set and say why in one line — don't carry it as
> permanent yellow.

---

## 3. The reusable diagnostic

When a first-order optimality check fails on a converged NLP, this separates
"my duals are wrong" from "my analytic reconstruction is wrong" in one run:

1. Get the raw multipliers: `lam = full(sol.value(opti.lam_g))`.
2. Assemble the full Lagrangian gradient by AD:
   `gL = grad_f + s·Aᵀ·lam`, choosing `s = ±1` by whichever gives the smaller
   `‖gL‖∞` (CasADi's overall Lagrangian sign is a convention, not a bug).
3. Read off the control block of `gL` and project out the directions absorbed
   by norm/cone constraints — for a unit-vector control, the **tangential**
   part `(I - ββᵀ)·∂L/∂β`.
4. Interpret:
   - tangential ≈ 0 **and** your analytic condition fails ⇒ the duals are a
     valid certificate; **your reconstruction or its inputs are wrong**.
   - tangential ≉ 0 ⇒ the multipliers are not a stationarity point; suspect
     extraction, convergence, or the model you rebuilt.

On the 10 N row this read `‖∇ₓL‖∞ = 1.5e-14` and tangential `8e-17` — which
immediately cleared `mee_primer_switch`'s derivation (correct all along) and
pointed at the inputs. It also eliminated the last structural suspect for
free: all 194 `Ldot`-guard multipliers were zero, so the one control-dependent
constraint that *could* have rotated the primer wasn't active.

Implementation: `results/dual_anomaly/diag_t1_beta.m`.

---

## 4. Standing rules for this codebase

- **Never use `opti.dual()` for a multiplier you will differentiate against.**
  Record the `opti.g` row range at build time; index `opti.lam_g`. This is
  already house practice in the CR3BP campaigns — it is now uniform.
- **Fixing a solver does not fix banked results.** Every `.mat` cache carries
  the derived quantities as of solve time. Dual repair needs an explicit
  refresh path: `verify/refresh_duals_mee.m`, `refresh_duals_cr3bp.m`.
- **Guard warm-restart refreshes on certified quantities, not node drift.**
  These extremals are weak minima (the SOSC study found 270 flat directions at
  10 N), so the solver legitimately slides within a flat optimal set — node
  drift reached 0.32 on rows whose mass was unchanged. Gate on status,
  machine-tight defect, and a **one-sided** final-mass check (minimum fuel
  means maximize m_f, so higher mass is never a failure — same convention as
  `reproduce/verify_row.m`).
- **A warm re-solve is a free under-optimization probe.** Refreshing duals
  found better optima on **six of nine** rows: 5 N +0.235 kg, 2.5 N +0.51 kg,
  1 N +0.44 kg, 1 N PSR +0.47 kg, 0.5 N +0.02 kg, 0.5 N PSR +0.02 kg. The
  three that did **not** move are exactly 10 N, 0.2 N and 0.1 N — the rungs
  built or rebuilt most recently and most carefully (drift ≤1e-2, Δm_f ≤1e-10).
  That pattern is independent confirmation of the reproducer engine's "the
  campaign under-optimized the warm-chained rungs" finding, and a concrete
  reason to revisit the mid-ladder Table-3 entries.

---

## 5. Artifacts

| file | what |
|---|---|
| `results/dual_anomaly/diag_optidual_minimal.m` | minimal reproduction of the CasADi behavior (3-variable NLP, three constraint writings) |
| `results/dual_anomaly/diag_t1_beta.m` | the T1 test: full-Lagrangian tangential β-stationarity; also measures the eccentricity correlation |
| `results/dual_anomaly/diag_rawdual.m` | structural fingerprint (magnitudes equal, signs differ) + the fix test |
| `tests/test_dual_extraction.m` | regression guard: locks the CasADi behavior **and** asserts β-stationarity of delivered duals at a converged point |
| `verify/refresh_duals_mee.m`, `../../earth_elliptic_to_geo_CR3BP/direct/refresh_duals_cr3bp.m` | repair banked caches |
| `verify/run_verify_pmp_all.m` | ladder-wide PMP sweep |
| `DESIGN_dual_map.md` | status banner with the full evidence chain; body preserved as the pre-resolution record |
