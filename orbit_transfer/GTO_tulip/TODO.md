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

- [x] **RESOLVED 2026-07-29 — there are at least SIX optima, not two, and the
  flagship is beaten.** A 13-seed sweep on the flagship's own mesh found
  converged optima at m_f 0.84708554, 0.84797502, 0.84805124, 0.84906583,
  0.84988608 and 0.85008728, with switch counts of 24, 25 AND 26 at the same
  t_f. The best (24 switches, from the 14-switch artifact) beats the certified
  flagship by **15.3 g of propellant, 0.677%**. Distinct seeds stay in distinct
  basins, so these are real and separate. Record:
  `process/BASIN_1150_SWEEP.md`. **Open follow-up:** certify the winner through
  `run_foc_tulip`/`verify_pmp` before it replaces the flagship, and re-run the
  sweep across the ΔV–t_f front, which is probably an envelope over a much
  richer structure.
- [ ] ~~**Two certified optima at t_f = 1.150×, 1.43% apart in ΔV.**~~ Found
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
  - **P0 — RUN 2026-07-26. Result: `cScale` transplants cleanly; the two-primary
    clock does not.** Pointed `casadi_energy_freetf` at tulip endpoints, nominal
    25 mN, the certified 1.15× horizon, seeded with the existing 8-state
    single-primary energy backbone (the solver appends `cScale` itself, so no
    seed conversion was needed — the state-dimension obstacle P1 assumed turned
    out not to exist for P0):

    | variant | iters | status | maxDefect |
    |---|---|---|---|
    | single-primary (`moonZone ≤ 0`) + `cScale` | 199 | **Optimal Solution Found** | **6.8e-14** |
    | two-primary (`moonZone = 0.15`) | 1403 | **Restoration Failed** | 3.5e-01 |

    The control is machine-tight and hits t_f exactly with `cScale = 1.0051`, so
    adding the free-time slack state to this transfer costs nothing. The
    two-primary failure is not marginal — defect 0.35 means genuinely infeasible,
    confirming the predicted node-placement incompatibility: a mesh laid out by
    `κ = r₁^1.5` cannot be warm-started into a clock that wants nodes
    redistributed toward the Moon.

    The two-primary clock was assumed necessary for tulip (coarse mesh at the
    arrival leg). P0 shows it is not necessary at nominal thrust — so the ladder
    was attempted next with single-primary + `cScale`.

  - **P0b — RUN 2026-07-26. `cScale` alone does NOT clear the wall.** Chained the
    same backbone to 20 mN holding the FACTOR (t_f 7.234 → 9.043 ND), same mesh
    (4001 nodes, same σ, same τ_f), single-primary, ε=1:

    | iter | inf_pr | inf_du | lg(rg) |
    |---|---|---|---|
    | 0 | 7.4e-03 | 8.2e-02 | — |
    | 154 | 1.3e-01 | 1.1e+03 | — |
    | 454 | 1.1e-01 | 5.4e+12 | 8.9 |
    | 848 | 6.2e-02 | 3.6e+14 | 11.2 |

    Primal infeasibility never came below ~2e-2, dual infeasibility diverged
    through fourteen orders of magnitude, and IPOPT piled on Hessian
    regularization (lg(rg) → 11) before MATLAB died with a **fatal MEX error** at
    iteration 848 — the campaign's known sporadic CasADi/IPOPT crash, and the
    reason every batch script runs one process per item. The 15 mN case never
    started. **The rung did not close.**

    **Two competing explanations, not yet separated:**
    1. **Resolution, not topology.** `cScale` scales `dx/dτ = c·κ·f`, so it *does*
       buy geometric arc per unit τ — the winding number can grow. But the node
       budget did not: ~25% more revolutions were asked of the same 4001 nodes,
       so each revolution is under-resolved and the defects cannot be met. If
       this is the cause, the fix is **scale N with the revolution count**, not
       reformulate. Cheap to test.
    2. **Genuine span limitation.** Fixed τ_f caps the arc regardless of `c`, and
       free span (free τ_f, or a Δθ-domain formulation) is genuinely required.

    Explanation 1 is the cheaper hypothesis and should be tested first: re-run
    20 mN at N scaled by t_f (≈5000 nodes) before concluding anything about
    topology. Note the earlier 20 mN pilot failure used the fixed-t_f solver at
    fixed N too, so it does not discriminate between these either.

  - **A1 — RUN 2026-07-26. THE 20 mN RUNG CLOSES. Neither explanation was
    needed: it was the CONTINUATION STEP SIZE.** Both the original pilot and my
    P0b jumped 25 → 20 mN in one 20% step. Walking the same path in ~4% steps,
    chaining each converged rung into the next (single-primary + `cScale`,
    ε=1, N=4001 throughout, one process per rung):

    | rung | status | maxDefect | cScale |
    |---|---|---|---|
    | 24 mN | Solve_Succeeded | 1.44e-14 | 1.0145 |
    | 23 mN | Solve_Succeeded | 2.98e-14 | 1.0199 |
    | 22 mN | Solve_Succeeded | 3.00e-14 | 1.0219 |
    | 21 mN | Solve_Succeeded | 3.46e-14 | 1.0211 |
    | **20 mN** | **Solve_Succeeded** | **2.33e-14** | 1.0222 |

    Every rung machine-tight, at the SAME node count that failed in one jump.
    So there is no resolution wall at 20 mN and no evidence of a topology wall
    either — the recorded "fixed-τ_f topology wall" was, at least at this rung,
    an artifact of an oversized continuation step.

    **What is NOT yet established, and matters:**
    1. **This is ε=1 (energy) only.** A certified min-fuel rung needs ε→0, and
       the energy band has always been wider than the ε=0-convergent band in
       this campaign. The pilot's recorded failure was also *at the ε=1 step*,
       so the comparison is apples-to-apples — but a certified 20 mN min-fuel
       rung does not exist yet.
    2. **Two variables changed at once.** A1 used small steps AND `cScale`; the
       pilot used a big step AND fixed τ_f. Which mattered is not separated.
       `cScale` stayed within 2% of 1.0 on every rung, which weakly suggests it
       is NOT doing the work — the discriminating run is small-step continuation
       with the ORIGINAL fixed-τ_f `casadi_minfuel_sundman`. If that also reaches
       20 mN, no reformulation is needed at all and the free-span item can be
       closed.

  - **B1/B2 — RUN 2026-07-26. BOTH variables mattered, and the 20 mN rung is
    now a min-fuel solution.**

    **B1, the discriminator** — same ~4% continuation with the ORIGINAL
    fixed-τ_f `casadi_minfuel_sundman` (no `cScale`):
    24 mN ok (2.9e-14), 23 mN ok (3.8e-14), 22 mN ok (5.6e-14), **21 mN FAILS**.

    So neither single explanation was right. Step size mattered — fixed-τ_f
    reaches 22 mN in small steps where a single 25→20 jump failed outright. And
    `cScale` mattered — fixed-τ_f stops at 21 mN while the free-time formulation
    continued to 20 mN. The fixed-τ_f limitation is real, just **further down
    than the pilot suggested and not a wall at 20 mN**.

    **B2, ε→0 at 20 mN** from the A1 energy solution, campaign schedule, 13 steps,
    every one converged:

    | | m_f | ΔV | propellant | switches | defect |
    |---|---|---|---|---|---|
    | 25 mN (published) | 0.849066 | 3.3696 km/s | 2.2640 kg | 25 | 2.0e-14 |
    | **20 mN (new)** | **0.859722** | **3.1127 km/s** | **2.1042 kg** | 11 | 4.5e-14 |

    Physically consistent: lower thrust → longer t_f (9.043 vs 7.234 ND) → less
    propellant and lower ΔV, moving toward the low-thrust limit.

    **What this is and is NOT.** It reached ε=0 machine-tight, which meets the
    `minfuel_at_tf` definition of certified. It has NOT been through the
    campaign's full gate set — no `boundSat` check, no fingerprint recorded, no
    `certify_minfuel_pmp` PMP propagation, no FOC gate. And it was produced by
    the ELFO solver with `moonZone=-1`, not by tulip's own machinery. Call it a
    certified-quality solve pending gating, not a gated campaign rung.

    The 11-switch count is worth checking rather than assuming: the 25 mN
    reference has 25. Fewer, longer arcs at lower thrust is plausible but is a
    topology change, and this campaign has already been shown to have multiple
    certified optima at one t_f.

  - **Ladder front door built + run 2026-07-27 (`run_tulip_ladder.m`).** Six
    rungs 25 → 20 mN in 4% geometric steps, every one machine-tight, final rung
    sharpened to ε=0.

    **THIRD confirmation of basin multiplicity, and the sharpest.** At IDENTICAL
    t_f (9.0429), thrust and endpoints, the ladder's 20 mN result differs from
    the 1 mN-step chain run earlier the same day:

    | continuation path | m_f | ΔV | prop | sw |
    |---|---|---|---|---|
    | 1 mN steps (A1/B2) | 0.859722 | 3.1127 km/s | 2.1042 kg | 11 |
    | 4% geometric (ladder) | **0.865477** | **2.9753 km/s** | **2.0178 kg** | 10 |

    4.4% in ΔV, 86 g of propellant, and the ladder's is BETTER. The first two
    confirmations showed the SEED ROUTE selects the basin; this one shows the
    **continuation path** does too — same endpoints, same horizon, only the
    intermediate rung spacing differs. Consequence: a single-path ladder is an
    upper bound on ΔV, not the front. Build ladders as envelopes over paths, the
    way `aggregate_front` already does over routes.

  - **Descent — RUN 2026-07-27. Ceiling found, and there are TWO stacked
    limits, not one.**

    | | lowest reached | evidence |
    |---|---|---|
    | ε=1 **energy** | **19.5 mN** | defect 5.1e-14, revs 21.1, 190 nodes/rev |
    | ε=0 **min-fuel** | **20 mN** | defect 4.5e-14, 11 sw, 2.1042 kg, ΔV 3.1127 km/s |

    5% steps stopped at 19 mN; halving to 2.5% bought one more rung (19.5 mN)
    and then stopped at 19.01 mN. But 19.5 mN **will not sharpen** — the ε
    homotopy fails at its FIRST step (ε=0.6, defect 2.2e-04, edge 1.6%). So the
    fuel ceiling (20 mN) sits above the energy ceiling (19.5 mN), which is this
    campaign's long-recorded pattern: the energy band is wider than the
    ε=0-convergent band.

    **The energy ceiling is NOT resolution.** 190 nodes/rev at 19.5 mN is
    generous, not starved — so "add nodes" is not the fix, and node count is not
    the argument for MEE phasing.

    **Likely cause of the energy ceiling: the t_f rule, not the formulation.**
    Rungs hold the factor by scaling t_f as `T_src/T_new`, i.e. assuming
    `t_f,min ∝ 1/T` exactly. That is an approximation, and the error compounds
    on the way down — by ~19 mN the requested t_f may sit BELOW the true
    `t_f,min(19 mN)`, which is genuinely infeasible. Four things support it:
    (a) the pilot findings already recorded this exact mode ("the first pilot
    run requested 0.92× tfMin(20 mN) and was genuinely infeasible");
    (b) the failure signature is infeasibility, not stiffness — `inf_pr` never
    descends while `inf_du` diverges to 7e11; (c) it is not resolution; and
    (d) **`cScale` climbs monotonically across the whole descent** — 1.0145
    (24 mN) → 1.0222 (20) → 1.0258 (19.5) — i.e. the free time scale is
    progressively stretching, exactly what it does when the requested horizon
    gets tight. That last signal is only visible BECAUSE the free-time state
    exists; the fixed-τ_f solver had no way to express "I need more time", which
    may be why it quit earlier still, at 21 mN.

    **TESTED 2026-07-27 — HYPOTHESIS REFUTED.** Ran it as two clean steps:
    (A) stretch 20 mN from factor 1.15 → 1.25 at fixed thrust — **converged**,
    defect 2.6e-14; (B) step 20 → 19 mN at the wider horizon — **failed
    identically** to the 1.15 case (`inf_pr` starts 4.8e-02 and never descends,
    `inf_du` → 1.4e11, MEX crash, no artifact). Giving the rung ~9% more
    transfer time changes nothing.

    So the ~19 mN energy ceiling is **not** the scaling law. Also withdrawn: the
    `cScale`-climbs-when-tight argument. Step A gave the rung MORE time and
    `cScale` still rose (1.0222 → 1.0297) — since `c = t_f / ∫κ dτ`, it rises
    whenever t_f outpaces the τ-integral, which happens both when the horizon is
    tight AND when a longer one is simply requested. It is not a tightness gauge
    and should not be read as one.

    **What is left.** Not resolution (190 nodes/rev), not t_f margin (tested),
    not step size at this scale (2.5% steps also stop). Both solvers wall in the
    same neighbourhood — fixed-τ_f at 21 mN, free-time at ~19.5 mN — which
    points at something about the PROBLEM near 19–21 mN rather than either
    formulation: plausibly a topology change (revolution count or lunar-encounter
    geometry) that smooth continuation cannot cross. Untested. The next probes
    would be (i) inspect how the trajectory geometry changes across 22 → 20 mN
    in the solutions we have, and (ii) try a discrete jump with a re-generated
    backbone at 19 mN rather than continuation.

    Artifacts: `results/ladder_2026_07_27/` (rung solutions + summaries);
    movie `results/minfuel_tulip_20mN.{mp4,gif}`, stills `_early/_mid/_late.png`.

    **Consequence for the plan.** The thrust ladder is NOT blocked — 25 → 20 mN
    is a 20% reduction on an item recorded as walled since 2026-07-21. The
    free-span reformulation is not justified; do not build it. Next: (i) the
    t_f-margin test above, (ii) gate the 20 mN rung properly (boundSat,
    fingerprint, PMP, FOC), (iii) port `cScale` into `casadi_minfuel_sundman` as
    an OPT-IN branch before any rung is banked as a campaign artifact. Next: (i) the fixed-τ_f
    small-step control run, (ii) ε→0 sharpening at 20 mN, (iii) extend the
    continuation downward to find where it actually breaks — with the
    expectation that the real ceiling is set by nodes/rev in Cartesian
    (~100/rev vs MEE's ~25), putting the practical floor near 5–10 mN rather
    than earth's 0.1 N.

    **Also measured:** the two-primary clock cannot be reached by warm-starting
    from a single-primary mesh (defect 0.35, Restoration Failed). Any move to
    the two-primary clock needs its own seed — a tulip analogue of
    `gen_elfo_energy_gravhom` — and is a separate question from the ladder.
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
