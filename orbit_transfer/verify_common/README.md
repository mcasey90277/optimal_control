# verify_common — shared first-order optimality (FOC) gate layer

One generic, AD-based first-order optimality checker that runs identically on
**all four** transfer campaigns, plus the fixed-format report every production
solve driver now emits. Built 2026-07-25 (plan
`docs/superpowers/plans/2026-07-25-foc-gate-layer.md`).

**Why it exists.** Before this layer, each campaign hand-derived its own PMP
verifier in its own coordinates (Gauss matrices, RTN components, the `K_L`
clock-coupling term, Sundman `r^1.5`). That is why there were four of them,
why they disagreed in coverage, and why GTO→ELFO had **no first-order gate at
all**. The fix was not a fifth verifier but a change of method: verify in the
**transcription's own variables**, by differentiating with CasADi the exact
`opti.f` / `opti.g` that were solved. A checker that differentiates the solved
model cannot drift out of sync with the physics — and needs to know nothing
about MEE, Sundman, or the Moon.

## Files

| file | what |
|---|---|
| `foc_check.m` | **The core.** Assembles the full Lagrangian gradient by AD from `opti.f`/`opti.g`/`opti.lam_g`, resolves the ± sign convention, and evaluates every first-order condition: KKT stationarity, the minimum condition (direction + throttle sign law), nodal costates, mass transversality, time-costate behaviour, singular arcs, regular switching. Returns one `rep` struct. |
| `foc_manifest.m` | Per-campaign bookkeeping — `nx`/`nu`, which state row is mass/time, which control rows are direction/throttle, autonomy, horizon kind. Empty `[]` fields mean **skip that check**, never "pretend it passed". Names: `earth_mee`, `earth_cr3bp`, `tulip`, `elfo_fuel`, `elfo_mintime`, `toy`. |
| `foc_dual_to_costate.m` | Interval defect duals → nodal costates, step-weighted (one-sided at the endpoints). Generic in `nx` (7/8/9 in use). |
| `foc_ipopt_inertia.m` | Interprets IPOPT's `δ_w` Hessian-regularization history: converge with `δ_w=0` and the reduced Hessian is PD without correction. LEAD-0 port of the tulip campaign's `psr_ipopt_certify`, which was the one second-order instrument that actually delivered. |
| `foc_report.m` | **The standard report.** Fixed-format block printed by every driver + a `foc_<tag>.mat` sidecar. Identical layout across all four campaigns so reports are comparable at a glance. |
| `setup_verify_common.m` | Path setup. Self-contained — adds only this folder, no campaign paths, no CasADi. |
| `doc/first_order_checks.tex` | **Read this first.** Every check explained twice: formal condition + intuition, then code + how to read the output. Annotated real report block, thresholds table, and the open caveats. Compile: `pdflatex first_order_checks.tex` (×2). |
| `tests/` | `test_foc_check_toy` (self-contained bang-bang OCP, no campaign deps), `test_foc_check_10N` (integration vs the certified Earth row), plus unit tests for the map, manifests, inertia interpreter, and report. Cache-dependent tests SKIP when the `.mat` is absent. |

## Running it

From each campaign's own folder, after its `setup_paths` and
`addpath(fullfile(getenv('HOME'),'casadi-3.7.0'))`:

```matlab
% Earth 2-body            (earth_elliptic_to_geo/direct/)
[rep, ver] = run_foc_mee('results/MEE_M2_10N.mat');
T = run_verify_pmp_all();          % whole ladder, generic + physical side by side

% Earth CR3BP             (earth_elliptic_to_geo_CR3BP/direct/)
ver = verify_cr3bp_pmp(struct('thrustN',10,'phi0',0));   % result in ver.foc

% GTO->tulip              (GTO_tulip/direct/)
rep = run_foc_tulip();             % + LS-vs-raw-dual costate cross-check

% GTO->ELFO               (GTO_ELFO/direct/elfo/)
rep = run_foc_elfo('mintime', fullfile('results','mintime_elfo.mat'));
```

Each entry point warm re-solves at the artifact's saved primal (to attach the
live model and re-derive multipliers), **guards that re-solve on certified
quantities** — solver success, defect ≤ 1e-8, one-sided final mass, or
non-degraded `t_f` for min-time — then checks and reports. Production drivers
(`run_gergaud`, `run_cr3bp_geo`, `run_psr`, `run_elfo_minfuel`,
`gen_elfo_mintime`) emit the report automatically after a live certifying
solve, try/catch-guarded so a report failure can never fail a solve.

## Policy: report-only burn-in

`rep.pass` is **advisory**. Every printed block says so. Nothing here can
demote a certified row or loosen an existing gate — the physical verifiers
(`verify_pmp_mee`, `certify_minfuel_pmp`) and the campaigns' primal
certification are untouched. Promotion to a hard gate is a recorded future
decision, blocked on the pre-promotion checklist in
`../OPTIMALITY_CERTIFICATION.md` §A6.

Two consequences worth internalizing:

- A FAIL is a **finding to investigate**, not a verdict on the solution — and
  the burn-in's own record argues for suspicion: of the five advisory failures
  this layer has produced, **all five turned out to be its own artifacts**
  (two mesh-normalization, three endpoint-representation). The one standing
  disagreement is the tulip LS-vs-raw-dual cross-check, where the LS side is
  independently documented as unreliable over ~40 revs.
- Conversely, a PASS here is first-order only. It certifies **extremal**, not
  **minimal** — see Part B of the register.

## Rules for extending this

- **Never `opti.dual()`.** Record the constraint group's `opti.g` row range at
  build time and index `opti.lam_g`. See `../README.md` Conventions and
  `../earth_elliptic_to_geo/process/LESSONS_DUAL_EXTRACTION.md` for the
  nine-day misdiagnosis this rule exists to prevent.
- **Adding a campaign = adding a manifest**, plus a `returnModel`/`creg` hook
  in its solver (additive, default-off, byte-identity tested). Do not fork
  `foc_check`.
- **Empty manifest fields skip checks.** If a condition does not apply, say so
  with `[]` and let the report print `--`. Never make a check pass vacuously.
- `foc_check`'s layout and manifest-semantic asserts (X-block-then-U-block;
  mass non-increasing; time non-decreasing) exist so a wrong manifest fails
  **loudly** rather than producing plausible nonsense. Keep them.

## Which file says what

Three documents, deliberately non-overlapping — check the right one before
adding to any of them:

| file | owns | does NOT own |
|---|---|---|
| `doc/first_order_checks.tex` | the **mathematics**: what each condition is, why, and how to read its report line. Present tense; corrections live in its provenance appendix | per-row results, open work |
| `../OPTIMALITY_CERTIFICATION.md` | the **cross-campaign status**: coverage matrix, per-row numbers, blocking mechanisms, open leads, decision record. Both orders | how an instrument works (link the header), module-local chores |
| `TODO.md` | **work on this module**: pre-promotion checklist, coverage holes, robustness items | results (they belong in the register) |

Function headers remain authoritative for how a given function works. None of
the three restates a derivation that lives in a header.

## Related

- `../OPTIMALITY_CERTIFICATION.md` — the cross-campaign register: per-row
  results, coverage matrix, both orders, open leads, and the §A6 pre-promotion
  checklist. **The register is the live status; this folder is the machinery.**
- `../earth_elliptic_to_geo/process/LESSONS_DUAL_EXTRACTION.md` — why the
  multiplier rule above exists, and the reusable
  "is-it-the-duals-or-my-derivation" diagnostic.
- `TODO.md` — open work on this module.
