# As-built package review (2026-07-22) — triage

Reviewers: GPT-5.6-terra + Gemini 3.1 Pro (opencode-headless; Gemini needed
3 attempts — large-bundle + `-c` triggers a provider 400, inline-only works),
host Claude adjudicating. Scope: the as-built Phase-1 package — front door
`run_cr3bp_geo.m`, campaign drivers, the `lt_mee_rhs` pert branch, tests,
and `homotopy_mee`'s cache machinery — steered at the two defect classes the
front-door validation had just caught live (cache-trust semantics S1, silent
basin selectors S2) plus user-facing robustness, as-built-vs-oracle drift,
and data-product honesty.

## Verdicts

GPT: **Hold** (16 findings). Gemini: **Ship with fixes** (11 findings,
heavily convergent with GPT). Host disposition: **all substantive findings
fixed in one wave (A1–A13 + B1–B4), validated, shipped** — the cold
front-door run reproduces the certified 1377.1545 kg / 19 sw / 4.2e-15
after the wave.

## Applied (A-wave, from GPT; B-wave, Gemini-new)

- **Fail-closed caching:** `homotopy_mee` gains opt-in `fpStrict` (missing/
  partial fingerprints error; cached FAILED steps re-solve) — default off,
  2-body campaign byte-identical; CR3BP callers set it. Front-door
  `check_fp_local` hard-errors on missing fields (B1). fp completeness:
  `x0Elems, maxIter, N, gainSched` added (A1). Mesh-consistency asserts on
  checkpoint load against partial deletes (B3 — same class as the v2
  poisoned-cache crash).
- **Basin guards:** seed knobs recipe-exact incl. no `max(60,·)` floor
  on-table (B2); rev-count warn vs `table3_certified.revs` (A4);
  certified-reference warning when a gate-passing solve lands >0.1 kg below
  the known result (A8 — the 1.84 kg silent-basin class).
- **Robustness/interface:** parameter validation block (A7); narrow catch
  on `table3_recipes` ids (A6); gainSched validation + recipe-sourced bridge
  seeds (A10); stronger final gates incl. maxUnit/termErr/boundSat (A5);
  ctfEff metadata under explicit tfTargetTU (A9); lunar-separation
  post-check (A11); oracle uses the production-guarded d³ + an MX-graph
  evaluation test of the pert branch (A12); `lamDef` relabeled honestly as
  defect-constraint duals — GPT's proposed converter `mee_dual_to_costate`
  does not exist (A13, modified).

## Rejected (do not act)

- **(Gemini top-3) remove the `+1e-12` d³ guard** — backwards: the guard
  protects unconstrained IPOPT trial iterates (the documented
  Invalid_Number_Detected class behind `LdotFloor` in the same file);
  physical bias at ≥8 LU separation ~1.6e-14; the test now uses the guarded
  expression.
- (Gemini) catch-filter keyed on `MATLAB:UndefinedFunction` — would mask a
  missing function on the path; the narrow known-id catch (A6) is correct.
- (Gemini) hard [6.5, 9] rev-window assert — 10 N-specific; the scaled
  vs-certified-revs warning generalizes across rungs.

Raw reviews: `2026-07-22_pkg_review_{gpt56terra,gemini}.md`. Fix commit:
see git log (`fix(cr3bp-geo): package-review fixes A1-A13 ...`).
