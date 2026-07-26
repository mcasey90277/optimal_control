# GTO_ELFO — TODO

Goal: **both methods working.** Direct is in good shape; the open work is
mostly on the indirect side.

## Direct — polish

- [ ] **Narrow the tulip path dependency (survey 2026-07-26,
  `doc/EXECUTION_PATHS.md`).** `setup_paths` adds all 34 files of
  `GTO_tulip/direct/sundman_minfuel` to reach exactly **one** function,
  `insertion_states`. Move that function into `cr3bp_common` (where
  `gto_tulip_endpoints` and `gto_elfo_endpoints` already live) and drop the
  tulip path. Removes a real shadowing surface: today an ELFO session that has
  also touched PSR carries two definitions of `casadi_minfuel_sundman`.
- [ ] **Fix the stale `setup_paths` header.** It claims the campaign "reuses
  casadi_minfuel_sundman / insertion_states / minfuel_at_tf, retargeted to
  ELFO". Measured: only `insertion_states` is called. The solver was FORKED
  (`casadi_energy_freetf`, `casadi_mintime_freetf`), not retargeted. The
  comment will mislead the next reader into thinking the solver is shared.
- [ ] **Extract the incidental duplication, not the solver.** `casadi_energy_freetf`
  shares **122 identical code lines (76%)** with tulip's
  `casadi_minfuel_sundman`, including a **41-line identical IPOPT options
  block**. Do NOT merge the solvers — the 24% that differs is the substantive
  part (9th `cScale` state, free t_f, two-primary clock, gravity homotopy) and
  merging risks two campaigns' certified results. DO extract the options block
  and the 6-line CasADi path bootstrap, which are incidental and appear in
  every solver in the repo.


- [x] **First-ever FOC (first-order-condition) gate landed (2026-07-25, FOC
  gate layer Task 9).** `run_foc_elfo.m` (new, `verify_common`-backed) runs
  KKT stationarity / min-condition / sign-law / transversality / regular-
  switching gates on three artifacts at the nominal 25 mN rung: energy seed
  (ADVISORY PASS, KKT 5.89e-13, eps=1 sign-law caveat), min-time anchor
  (ADVISORY PASS, KKT 7.25e-14 — `lamTimeEnd=-1.000` exactly confirms the
  free-t_f dual-form transversality condition λ_t(t_f)=±1; `lamTimeCoV=5.69e-2`
  left as an OPEN lead, see `../OPTIMALITY_CERTIFICATION.md` LEAD-4), and the
  1.33×/50-switch front row (ADVISORY FAIL — sdotMinRel 4.485e-5 near-graze
  switches + δ_w 1.98e-6 borderline, same weak-minimum pattern as the tulip
  flagship, not a new defect). `foc_ipopt_inertia` (LEAD-0 port) also now
  gives this campaign its first second-order verdict: energy seed + min-time
  anchor LOCAL MIN, front row NOT CERTIFIED. Report-only burn-in — advisory,
  does not gate certified status. Not yet covered: other thrust rungs, the
  PSR-equivalent refinement path, and a full ladder sweep. Details:
  `.superpowers/sdd/2026-07-25-foc-gate-layer/task-9-report.md`.
- [ ] **Near-min-time end of the ΔV–t_f front.** The front is mapped and
  labeled against t_f,min = 6.0962 ND (27.02 d), but the transition region
  just above min-time deserves the same scrutiny the tulip front got (the
  tulip campaign's 1.01–1.11× band was hard even for the smooth energy
  problem — check whether ELFO shows the same wall).
- [ ] **Switch-count bands.** Report ELFO switch counts as mesh-convergence
  bands, not integers (lesson from earth_elliptic_to_geo's P0 study).
- [ ] Keep `elfo_export_data` / movies current as new front points land.
- [ ] **Thrust ladder.** Same goal as the tulip and earth-GEO campaigns: sweep
  T_max around the nominal 25 mN with per-rung min-time anchors +
  fixed-c_tf fuel solves (thrust-continuation warm-chaining, certified-only
  caching), and check the T·t_f,min ≈ const law analog for the ELFO target.
  Port the `../earth_elliptic_to_geo/` ladder recipe.

## 2026-07-21 review follow-ups (doc/reviews/2026-07-21_triage.md)

Fixed same day: factor semantics rebased on `tfMin_ELFO` with tf-nearest seed
selection (C1); resume path certification-gated + single-trajectory saves
(C2); acceptance gates require `Solve_Succeeded` (C3).

- [x] **tf-sweep tight re-clean (C4).** DONE (ladder-prep T5): `solve_tf` now
  tight-re-cleans a clean loose probe before banking.
- [ ] **Mesh-band repeat + filename hygiene (C5).** Repeat ≥1 front point at a
  refined mesh and report switch counts as bands; parameterize N with t_f;
  drop the switch integer from identity-bearing filenames at the next format
  change. (Still open — not part of ladder-prep.)
- [x] **Ladder-prep trio (C6, feeds the thrust-ladder goal above):** DONE
  (ladder-prep package, Tasks 1–6): per-rung thrust via `minfuel_config(over)`
  + fingerprints on seeds/checkpoints/rows + `elfo_find_energy_seed` fp filter;
  scaled `cBox`; `boundSat` on both freetf solvers; `gen_elfo_mintime` anchor
  save gated on certification; `chain_rung_seed_elfo`. **Pilot gate: ELFO
  20 mN PASSED** — certified=1, defect 1.5e-15, 11 switches, first certified
  off-nominal rung in either campaign; certified at ~0.99× naive 1/T min-time
  (first `T·t_f,min` deviation signal). Writeup:
  `../GTO_tulip/process/LADDER_PREP_PILOT_FINDINGS.md`. ELFO can ladder the
  near-nominal band now (the `cScale` slack avoids the winding wall the tulip
  engine hit).

## Indirect — get it working

- [ ] **Route C: ELFO min-time indirect.** The direct min-time anchor (Route B)
  is certified; the indirect counterpart was scoped but never built. Start
  from `../GTO_tulip/indirect/min_time/` (PMP always-burn shooting root,
  pumpkyn `tfMin` machinery) retargeted to `gto_elfo_endpoints`.
- [ ] **Indirect min-fuel.** After Route C: PMP bang-bang shooting seeded from
  the certified direct solution (the same direct-seeded strategy as
  `../GTO_tulip/indirect/ifs/` — reuse its lessons, including the
  terminal-cluster conditioning failure mode).
- [ ] Replace the `indirect/` README stub with real structure once the first
  solver lands.
