# Booster Landing — Final Whole-Branch Review: Fix-Wave Report

Date: 2026-08-09
Scope: apply all findings from the final whole-branch review in one pass, one
commit. Doc/gate-reporting edits only; the one code guard is drag/cone-path
only. Vacuum-certified numbers unchanged.

## I1 — G5 primer sub-check invalid under a finite pointing cone

**Files:** `certify/certify_pdg.m`, `certify/print_certify_report.m`,
`doc/booster_landing_note.tex`

- `certify_pdg.m`: added a header note ("G5 primer INVALID under a finite
  pointing cone") explaining why the primer-vector necessary condition
  (thrust parallel to the velocity costate) does not hold once
  `P.theta_max_deg` is finite — the optimal control then sits on the cone
  boundary (KKT projection of the unconstrained direction), so `T` is not
  parallel to `lambda_v` by design, not by bug.
- Code change: the G5 block now branches on `isfinite(P.theta_max_deg)`.
  When the cone is active, `rep.G5_primer_deg` is still computed and
  reported (info only), `rep.G5_primer_mode = 'skipped-cone'` is set
  (mirrors the G3/G4 `'skipped'` pattern used when there's no convex twin),
  and `rep.G5_pass = rep.G5_structOk` (primer no longer gates). When the
  cone is `Inf` (every certified result in this campaign), behavior is
  byte-identical to before, including the drag-loosened 10 deg bound —
  that bound now only ever applies in the `Inf`-cone branch, so a cone run
  cannot also be scored against it.
- `print_certify_report.m`: added a `'skipped-cone'` branch that prints the
  primer row as an info row (value shown, threshold `(skipped-cone)`, no
  verdict), same shape as the G3/G4 `skipped` rows above it.
- `doc/booster_landing_note.tex`: extended the pointing-cone future-work
  bullet (§Future work) with a paragraph on the projection/skip semantics,
  so a future reader who enables the cone knows G5's primer check is
  cone-aware, not silently wrong.
- **Verification (measured, not assumed):** ran both solvers at
  `P.theta_max_deg = 15` (nominal grid, `N=60`/`Nconv=120`). Both riders
  bind the cone at the boundary: collocation max angle 14.9998 deg, convex
  14.9963 deg (agreement ~3.5e-3 deg, matching the review's ~4e-3 deg
  claim). The OLD naive primer check on that same solution reads 7.32 deg
  at the nominal grid (6.95 deg at a coarser N=40 check) — comfortably over
  the vacuum <1 deg gate on a fully-converged optimum, confirming this was
  a real false-fail risk, not a hypothetical one. `P.theta_max_deg`
  defaults to `Inf` everywhere in this campaign's shipped results, so no
  certified number changes; `test_certify_nominal` re-run green (below).

## I2 — Spec self-contradiction (README-designated tiebreaker)

**File:** `docs/superpowers/specs/2026-08-08-booster-landing-design.md`

- Success-criteria bullet: `"final masses agree < 0.1 kg"` →
  `"final masses agree < 1.0 kg (adjudicated: the residual is the measured
  Taylor-bound model error, ~0.43 kg)"` — matches the README's flagship
  `G3 |dmf| = 0.43 kg (gate: < 1.0 kg)` and `certify_pdg.m`'s adjudicated
  1.0 kg gate.
- Same bullet: `"terminal max-throttle arc to zero velocity at the pad"` →
  `"...to the adjudicated 1.5 m/s terminal descent rate at the pad"` —
  matches `P.vf = [0;0;-1.5]` and its own adjudication note two lines
  below in the same bullet, which the old text directly contradicted.
- Conventions section (M2): qualified the `if nargin==0` self-demo
  convention to `"where meaningful (several library functions are nullary
  or campaign-internal; none shipped a self-demo — convention relaxed at
  final review)"` — the spec no longer states an unqualified rule that
  every shipped file in fact deviates from.

## M1 — README task pointer

**File:** `README.md`

`task-{1..10}-{brief,report}.md` → `task-{1..12}-{brief,report}.md`
(verified: 12 task briefs/reports actually exist in
`.superpowers/sdd/2026-08-08-booster-landing/`).

## M3 — `solve_pdg_convex` single-`tf` path

**File:** `lib/solve_pdg_convex.m`

The `opts.tf` fast-path built `tf_curve` from `sol.mf` unconditionally,
contradicting the function's own docstring ("`v` (column 2) is `mf` for
BOTH code 2 and code 3, `-Inf` otherwise"). Fixed to capture `v0` from
`mf_or_neginf(sol)` and store that in column 2, matching the golden-search
path just below it. `test_convex_lossless` (which drives this file) re-run
green.

## M4 — `certify_pdg` header: "G2ff = second half of G2"

**File:** `certify/certify_pdg.m`

Added a sentence to the G2 entry in the gate enumeration at the top of the
file: G2ff is framed explicitly as the second half of G2 (same
reconstructed control, feedforward-feasibility question instead of
continuous-residual question), not an unrelated bolt-on gate.

## M5 — README: `.mc` presence parity with `.mcD`

**File:** `README.md`

The `R` returns... line listed `.mc` unconditionally alongside
`.P .solC .solV .rep .ctrl .out0 .when`. Corrected to state `.mc` is present
only when `cfg.doMC` is true — the same gating already documented (and
correctly worded) for Phase 2's `.mcD`.

## M8 — Front-door Phase-2 stage banner

**File:** `run_booster_landing.m`

`=== [P2] Drag-on re-solve (warm-started) ===` → `=== [P2 +1] ...` — makes
explicit that Phase 2 is one stage appended after the 5-stage Phase-1
sequence, rather than an unnumbered stage that looked orphaned next to
`[1/5]`..`[5/5]`.

## Bonus — `tvlqr_design` zetaA floor assert

**File:** `lib/tvlqr_design.m`

The derived-`QA` path already documented (in a comment) that `zetaA >=
1/sqrt(2)` is required for `qVel >= 0` (LQR's Butterworth floor) but never
checked it. Added
`assert(zA >= 1/sqrt(2) - eps, 'tvlqr_design:zetaAFloor', ...)` right where
`zA` is read, before it's used to derive `qV`, so a below-floor
`opts.zetaA` throws a clear error instead of silently producing a negative
weight.

## Test results (MATLAB R2025b, `-batch`, this run)

| Test | Result |
|---|---|
| `test_certify_nominal` | PASS (both blocks: coarse `all_pass`=PASS, nominal `all_pass`=PASS; G3 `|dmf|`=0.4337 kg coarse / 0.4279 kg nominal, both < 1.0 kg; G5 primer 0.853 deg coarse / 0.574 deg nominal, both < 1 deg — vacuum path, `theta_max_deg=Inf`, byte-for-byte unaffected by the I1 guard) |
| `test_convex_lossless` | PASS — `gap=4.81e-05`, `mf=26435.3 kg` |
| `test_run_front_door` | PASS — `tf=16.596 s`, `fuel=3535.4 kg`, gates ALL PASS, MC 6/6 landed (100%) |

No regressions; all vacuum-certified numbers (tf, mf, gate values) match
the README's flagship table.

## Documentation compile

- `doc/booster_landing_note.tex` recompiled 3x with `pdflatex`
  (2x required by the task; a 3rd pass was run because the 2nd pass still
  emitted "Label(s) may have changed" — 3rd pass came back clean, 23 pages,
  no errors, no undefined references). Aux files cleaned (`.aux .log .out
  .toc .fls .fdb_latexmk .synctex.gz`).
- `verify-paper` (`~/ai_council/venv/bin/python
  ~/Documents/myLatex/tools/paper_verify.py booster_landing_note.tex`):
  citations 10 verified / 1 to-review (`betts2010`, pre-existing WARN,
  untouched by this fix wave) / 0 failed / 1 skip (non-DOI software cite,
  expected); figures 4/4 OK. No FAIL.

## Files touched

- `booster_landing/certify/certify_pdg.m`
- `booster_landing/certify/print_certify_report.m`
- `booster_landing/doc/booster_landing_note.tex`
- `booster_landing/lib/solve_pdg_convex.m`
- `booster_landing/lib/tvlqr_design.m`
- `booster_landing/run_booster_landing.m`
- `booster_landing/README.md`
- `docs/superpowers/specs/2026-08-08-booster-landing-design.md`
