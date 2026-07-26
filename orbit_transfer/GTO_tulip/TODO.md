# GTO_tulip — TODO

Two standing goals (2026-07-21): **(a) keep perfecting the direct code,
(b) get the indirect code working.** Full history behind every item:
`process/LOW_THRUST_MINFUEL_CAMPAIGN.md`.

## (a) Direct — perfect what works

- [x] **Front door.** DONE 2026-07-26: `direct/sundman_minfuel/run_gto_tulip.m`
  (house `run_gergaud` / `run_cr3bp_geo` pattern) + `run_tulip_front.sh`,
  guarded by `test_run_gto_tulip.m` (no-solve). This was the last campaign
  without one. It drives `minfuel_at_tf` rather than re-implementing the chain,
  and **refuses** the two knobs that name open problems below — off-nominal
  thrust (topology wall) and factors under ~1.12 (no energy backbone) — with
  messages that route to the relevant script and findings doc. The sweep script
  walks **t_f, not thrust**, for the same reason.

- [ ] **Two certified optima at t_f = 1.150×, 1.43% apart in ΔV.** Found
  2026-07-26 while smoke-testing the new front door. The energy-backbone route
  (`minfuel_at_tf`) and the `run_certified_minfuel` chain both reach 25
  switches, machine-tight, at the SAME t_f, but different local optima:
  m_f 0.847086 vs 0.849066 (ΔV 3.4176 vs 3.3696 km/s, +0.0297 kg). The seed
  route decides the basin.

  **MEASURED 2026-07-26 across every banked multi-route factor.** The spread is
  real, is NOT confined to 1.15, and is far larger than first estimated:

  | factor | routes (switches) | ΔV spread |
  |---|---|---|
  | 1.150 | certified-chain (25) vs energy (25) | **1.43%** |
  | 1.300 | en (42) vs nb (42) | **0.76%** |
  | 1.350 | en (29) vs nb (28) | 0.06% |
  | 1.400 | en (26) vs nb (26) | 0.02% |
  | 1.450 | dn (43) vs nb (26) | **9.84%** |

  All are ε=0, machine-tight (defect ≤5e-14), ~99.6% bang-bang, at identical
  t_f and endpoints — genuinely distinct extremals, not convergence noise. The
  1.450 pair differs by 17 switches.

  **The published front is NOT biased by this** — the earlier worry that it
  "may sit systematically ~1% high" is REFUTED: `aggregate_front` builds its
  envelope as `min(ΔV)` per factor over *all* points, so it already takes the
  best known. Open questions that remain: is any of these the global optimum;
  is there a cheap route-diversity policy (`minfuel_at_tf` keep-best-mass
  across routes, cf. the earth campaign's razor-thin-basin lesson); and does
  the spread correlate with switch count rather than factor (1.450's 43-vs-26
  suggests topology, not tolerance).
  `run_gto_tulip` now cross-checks and reports the basin, so this cannot pass
  unnoticed, but the underlying question is unresolved. **Related evidence:**
  the C3 mesh-band item below already recorded a re-solve landing 24 switches /
  ΔV 3.3660 ("basin-sensitive even at fixed mesh") — that one differed by
  ~0.1%, this one by 1.43%, so the spread is wider than a switch-integer
  wobble and the two should be investigated together.

- [ ] **Close the 1.01–1.11× near-min-time band.** The certified front starts
  at 1.12×; the energy backbone itself is only generatable for
  t_f ≈ 1.12×–1.95× (the band is hard even for the SMOOTH energy problem —
  a near-min-time conditioning wall, not just bang-bang structure). Candidate
  attacks: direct continuation from the min-time anchor downward, or the
  indirect band campaign below.
- [x] **PSR/lib de-dup — DONE 2026-07-26, by dissolving it.** The 2026-07-26
  flatten merged `sundman_minfuel/` and `PSR/` into one `direct/` campaign with
  a single front door, and deleted the 20 vendored files outright. There is now
  exactly one copy of every shared name, asserted by
  `direct/tests/test_run_gto_tulip.m` check 3.

  The de-dup was safe precisely because the vendoring had already stopped
  working: an `addpath` added on 2026-07-15 put the upstream folder ahead of
  `PSR/lib`, so PSR had been running the UPSTREAM solver for eleven days while
  its README claimed self-containment, and `casadi_minfuel_sundman` /
  `minfuel_at_tf` in `PSR/lib` were dead code. Removing dead code is not a
  behaviour change — which is why no re-certification was required.

  The PMP verifier (`verify_direct_pmp`, `sms_*`) was NOT re-vendored: it lives
  in `indirect/ms_band/` and `direct/setup_paths.m` names it. The direct
  pipeline checking itself with the independent indirect instrument is a real
  dependency, not duplication.

- [ ] **Front hygiene.** Keep `aggregate_front`'s honest 3-class front current
  as new t_f points land; switch counts reported as bands (mesh-sensitivity
  lesson from the earth_elliptic P0 study applies here too).
- [ ] **Free-span reformulation — the unblock for the thrust ladder.** The tulip
  solver holds `τ_f` FIXED (a free scalar would make a dense KKT column and kill
  the MUMPS factorization at N≈4000). The cost is that a chained seed's winding
  number cannot grow, so a lower-thrust rung — which needs MORE revolutions —
  cannot be reached. The 20 mN pilot failed there, at the ε=1 energy step. This
  is a topology wall, not tuning.

  **What is reusable, measured 2026-07-26.**

  *Not* the earth CR3BP campaign, despite it having walked a full 10 N → 0.1 N
  ladder on a free `ΔL` span. Its parameterization is **Earth-centred MEE**, and
  the tulip target sits **28,303 km from the Moon** — inside the lunar SOI
  (~66,000 km). Earth-centred elements are the wrong coordinates for the arrival
  leg. Only the *principle* transfers: a free scalar coupled to every node must
  be replicated per node and tied by local continuity constraints (`liftDL`).

  **ELFO is the real reuse, and it is close to a drop-in.**
  `casadi_energy_freetf` already has BOTH things the tulip solver lacks:
  1. **free time via the `cScale` slack state** (Betts) — banded KKT, and the
     revolution count can grow. This is what let ELFO's 20 mN rung CERTIFY where
     tulip's failed, in the same dynamical setting.
  2. **the two-primary Sundman clock** `κ = (r1^-q + (r2/D)^-q)^(-p/q)` — and the
     tulip campaign needs this independently of the ladder: its single-primary
     `κ = r1^1.5` gives a COARSE mesh at the arrival leg (r1 ≈ 398,000 km) which
     is exactly where the Moon dominates and the dynamics are fastest.

  Critically, `casadi_energy_freetf(sigma, rv0, rvf, ...)` takes the endpoints as
  **arguments**, and `insertion_states('tulip','campaign')` supplies tulip's. So
  the solver can be pointed at this transfer without modification.

  **Phased plan.**
  - **P0 (cheap, decides everything).** Run the ELFO free-t_f solver on TULIP
    endpoints at nominal 25 mN and the certified 1.15× horizon. Does it close at
    ε=1? If yes, tulip gains a free-span solver essentially for free and the
    ladder becomes an ELFO-style chain. If no, learn where it breaks before
    committing to anything larger.
  - **P1 Seeding.** P0 needs an ε=1 root in the 9-state free-t_f, two-primary
    representation. The existing tulip energy backbone is 8-state, single-primary
    — its node placement does NOT transfer (the two-primary clock redistributes
    the mesh, so the no-resample discipline is violated). Either re-map with an
    explicit resample and accept the defect hit, or build a tulip analogue of
    `gen_elfo_energy_gravhom`.
  - **P2 Ladder.** With P0/P1 done, chain rungs by HOLDING THE FACTOR, not t_f
    (t_f,min scales ~1/T, so the same t_f at lower thrust can sit below the new
    rung's minimum time — the first 20 mN pilot attempt was genuinely infeasible
    for exactly this reason). `chain_rung_seed_tulip` already does the
    time-rescale; it would need the 9-state layout.
  - **P3 Gate.** Certified rung = ε=0 reached, machine-tight defect, clean
    `boundSat`, fingerprint recorded — the ladder-prep gates already exist.

  **Risks to name up front.** ELFO's success does not guarantee tulip's: the
  targets differ in geometry (the tulip insertion is a high-inclination
  south-pole orbit). Adding `cScale` changes the state dimension, so every
  downstream consumer of the tulip 8-state layout (refine, PMP verifier, export,
  movies) needs auditing. And a second tulip solver would be a FORK — if it is
  built, `test_solver_fork_parity` should be extended to cover it rather than
  letting a third copy drift.

- [ ] **Thrust ladder (Table-3 analog for the tulip).** Port the
  `../earth_elliptic_to_geo/` ladder recipe — per-rung min-time anchor,
  thrust-continuation warm-chaining, certified-only caching, R0-law check —
  to the tulip problem: sweep T_max around the nominal 25 mN, map
  t_f,min(T) and m_f(T), and test whether a T·t_f,min ≈ const law holds in
  the CR3BP. Prior art: the indirect ztl campaign was a thrust-ladder
  attempt (its P0 findings — ~75 mN sweet spot, fixed-t_f ladder argument —
  constrain the design); the earth-GEO MEE ladder is the proven direct
  machinery. FIRST TASK — pumpkyn stress probe (ztl-style): run pumpkyn
  `tfMin` at 22.5 and 20 mN, cold AND warm-chained from the 25 mN solution,
  to find where the indirect single-shooting min-time substrate dies. Do not
  block on it: per-rung anchors can be DIRECT all-burn solves (ELFO Route-B
  precedent) or R0-law estimates (earth 0.5 N precedent) — min-time + high
  revs is the regime where direct shines (no switch structure) and single
  shooting is weakest (STM sensitivity grows with revs).

## 2026-07-26 prune findings

- [x] **`run_psr` stage 5c never worked — FIXED.** The FOC/KKT gate added
  2026-07-25 called `run_foc_tulip(finalSol)` from inside PSR. PSR is
  self-contained by design, `run_foc_tulip` and `verify_common` are outside its
  path boundary, so the call threw `Undefined function` and the surrounding
  catch turned it into an advisory warning on **every** run — the gate never
  executed once. (Its reviewer had flagged the wiring as "static-checked only,
  not exercised live"; that risk was real.) Now prints a pointer telling the
  user to run the report from `sundman_minfuel`, which preserves the isolation
  `PSR/lib/README.md` promises. Second, independent reason it could not have
  worked: the vendored solver lacks the `returnModel`/`creg` hook `foc_check`
  needs.
- [x] **`indirect/min_time/direct_mintime_elfo.m` retired to `attic/`.** Its
  own README already logged it as non-converging and it was superseded by the
  certified `casadi_mintime_freetf` Route-B anchor — and it was broken besides:
  its `addpath` pointed at `indirect/attic`, which does not exist, so none of
  its four dependencies resolved.
- [ ] **Attic is 1.6 MB / 30 files and off-path** (verified: nothing on any
  campaign path resolves into it). Left as the archival record it is meant to
  be; `attic/README.md` documents it. Prune further only if it starts costing
  something.

## 2026-07-21 review follow-ups (doc/reviews/2026-07-21_triage.md)

Fixed same day: acceptance gates now require `Solve_Succeeded` (C1) and
`certified` requires the requested homotopy endpoint (C2).

- [ ] **Mesh-band study for the flagship (C3).** Re-solve the certified 1.15×
  solution at ≥2 finer meshes; publish the switch count as a band (the
  earth-campaign P0 protocol). LIVE EVIDENCE (2026-07-21 gate-hardening
  regression): a fresh re-solve of the certified recipe landed 24 switches /
  ΔV 3.3660 vs the published 25 / 3.3696 — same mass to ~0.1%, switch integer
  basin-sensitive even at fixed mesh. **CONFIRMED again 2026-07-26** (post-flatten
  reproduction): 24 switches, m_f 0.849279, ΔV 3.3644 km/s, 2.2608 kg, defect
  1.3e-14, certified. Three data points now say the same thing — mass and ΔV
  reproduce to ~0.1%, the switch INTEGER does not. Report it as a band. Note the
  2026-07-26 re-solve was *marginally better* than the published artifact
  (m_f 0.849279 vs 0.849066), so the published 25-switch row is not the best
  point in its own neighbourhood either. Also relax `certify_minfuel_pmp`'s strict
  integer PMP-crossing match (node-grazing switches fail it spuriously).
- [x] **Bound-saturation diagnostic + box widening (C4).** DONE (ladder-prep
  T3): `casadi_minfuel_sundman` now emits `out.boundSat` + a saturation
  warning; velocity/position boxes are opt-in-widenable via the trailing
  `opts` arg (default byte-identical). The 20 mN tulip pilot's
  `boundSatWorst=vBox` is the first lead this diagnostic produced.
- [x] **Ladder-prep trio (C5/C6, feeds the thrust-ladder goal above):** DONE
  (ladder-prep package, plan `docs/superpowers/plans/2026-07-21-ladder-prep.md`,
  Tasks 1–6): per-rung thrust via `minfuel_config(over)` + `cr3bp_fingerprint`/
  `check_cr3bp_fp`/`thrust_tag` on every cache boundary; `chain_rung_seed_tulip`
  (time-rescale + no-resample re-map, fresh `tauf0`); `rF` fallback added to
  `gen_tulip_energy_2p` step_solve. **Pilot gate:** ELFO 20 mN PASSED
  (certified, defect 1.5e-15); tulip 20 mN was an HONEST FAILURE — fixed-τf
  topology wall (chained 25 mN winding can't grow to the 20 mN/1.15× rev count;
  the ε=1 energy step won't close). Full writeup:
  `process/LADDER_PREP_PILOT_FINDINGS.md`. NEXT (ladder campaign): the cheap
  widened-`vBox` disambiguation probe, then the θ-domain/Δθ-free spiral
  reformulation if the wall is confirmed (tulip has no CR3BP MEE+ΔL analog).

## (b) Indirect — get it working

Today: machinery built and validated, no certified indirect solve yet.

- [ ] **IFS to certification.** `indirect/ifs/` (direct-seeded indirect
  finishing solve) stalls at ‖R‖ ≈ 0.023 on terminal-cluster conditioning at
  the full 1.12× problem. Next move (per campaign record): retarget to a
  clean-band t_f whose switches are non-grazing, then walk back to 1.12×.
- [ ] **ms_band.** Multiple-shooting attack on the 1.01–1.11× band — blocked by
  the same near-min-time wall; revisit after IFS certifies anywhere.
- [ ] **Use the ztl P0 findings** (recorded in `indirect/ztl/`): the min-time
  substrate is dead as a ladder start; fixed-t_f ladder argument; the
  cold-landscape asymmetry ("arrives-warm easy"); ~75 mN sweet-spot signal.
  These constrain which indirect strategies are worth another attempt.
- [ ] **Success bar:** an indirect (PMP shooting) solution of the 1.15×
  min-fuel problem certified against the direct result (25 switches,
  ΔV 3.3696 km/s) — or an honest documented refutation of why single/multiple
  shooting cannot close it at this scale.

## Housekeeping

- [ ] `indirect/min_time/` serves both tulip and ELFO retargeting — if it grows,
  consider promoting it to a shared home (noted in the restructure spec).
