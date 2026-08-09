# costate_common — TODO

- [ ] **Promote to a proper MATLAB package** (`+oc` or similar) once the
  interfaces settle; the SDD is the architecture document. Until then this
  flat folder + `addpath` is the contract.
- [ ] **Schema versioning** (accepted debt, SDD): version the sheet-.mat
  schema and the compact catalog format before a second data consumer
  exists; replace `cat.derive` executable strings with a named formula
  registry; pin environment (MATLAB/CasADi/IPOPT/pumpkyn revisions) in
  campaign metadata enforcement, not just recording.
- [ ] **Acceptance-gate harness**: the tfMin acceptance gate is the one
  pipeline check not yet generic — generalize as a harness taking a
  per-family independent solver.
- [ ] **ms_bvp fixed-tf variant** (the header forbids faking it with tight
  guards); needed the day a fixed-time catalog appears.
- [ ] **Conjugate test refinements**: junction-resolution disclaimer could
  be closed by dense STM sampling inside flagged intervals; an exact-zero
  det sample is currently counted as a crossing — refine to a bracketed
  root report.
- [ ] **Batched-driver + monitor templates**: halo/DPO shell drivers differ
  only in names; extract, and make monitors watch process liveness (a
  silent MATLAB death produced no log line and no alert, 2026-08-07).
- [ ] `survey_family_bounds` / older files still carry `%#ok` pragmas from
  before the no-pragma rule; strip on next touch.
- [ ] `golden_cells` engine cells converge in 1 iteration (wide Newton
  basin) — consider a rougher engineered seed so the iteration channel has
  more dynamic range.
