# IFS Plan of Attack (2026-07-11)

Synthesis of the IFS docs (`README.md`, `RESULTS.md` incl. the post-merge
scaled-SVD investigation), the GPT-5.6-sol methodology review
(`reviews/gpt56sol_2026-07-11.md`), the campaign record
(`../LOW_THRUST_MINFUEL_CAMPAIGN.md`), the PSR sibling
(`../sundman_minfuel/refine/`), and the two reference papers in
`../../min_fuel_papers/`. Written so the next session can pick this up cold.

## 1. Where we stand

- **PSR (direct workhorse) works and is merged.** Certified many-switch
  bang-bang solutions exist (1.12×–1.25× + 1.85×), with KKT-dual costates.
- **IFS machinery is validated** (EOM/residual/Jacobian/seed unit tests all
  pass; four failure modes diagnosed and fixed) but the full 1.12× (k=10)
  gate **does not converge**: LM crawl 1.96 → 0.023 over 400 iterations,
  cond(J) ≈ 5.9e9, GN-consistent, 6 orders short of the 1e-8 target.
- **The conditioning is diffuse and shifting, not surgical** (post-merge
  scaled-SVD finding): near-null direction is 76% on switch #4's node (the
  shallowest S-crossing, |dS/dτ|=0.11) at the seed, then 83% on the initial
  position costate λ_r0 at mid-crawl. The near-double-root terminal pair was
  **refuted** as the culprit; τ_f over-constraint was **ruled out**
  empirically. This is the textbook weak spot of indirect shooting: λ_r is
  only indirectly coupled to the trajectory and weakly determined by a short
  first arc.

## 2. What the two reference papers say about exactly this

### Zhang, Topputo, Bernelli-Zazzera, Zhao, JGCD 38(8), 2015 (the existence proof)

PDF: [`../../min_fuel_papers/2015-8-J-LowThrustMinimumFuelOptimizationRTBP.pdf`](../../min_fuel_papers/2015-8-J-LowThrustMinimumFuelOptimizationRTBP.pdf)
("Low-Thrust Minimum-Fuel Optimization in the Circular Restricted Three-Body
Problem," DOI 10.2514/1.G001080)

Same problem class (Earth–Moon CR3BP, indirect min-fuel, GTO → L1 halo).
Their hardest case — 0.6 N, 140.3 days, **~150 bang-bang switches over ~150
revolutions** — converges. How:

- **Single shooting**, unknowns = λ_i ∈ R⁷ only (no node states, no switch
  times as unknowns). Residual: r,v match (6) + λ_m(t_f)=0 (1). t_f fixed
  = c_tf × t_fmin.
- **Exact analytic STM** via variational equations (14+196 ODEs), with a
  **saltation/jump matrix** Ψ = I + (ẏ⁺−ẏ⁻)(∂S/∂y)/Ṡ at every event-detected
  switch, composed into Φ(t_f,t_i). (Ψ's 1/Ṡ blows up at grazing switches —
  their known weakness, consistent with our wall; they never hit it in their
  cases.)
- **Machine-precision hybrid switch detection**: Newton on S(t_sw)∓ε with
  analytic Ṡ (4–5 iterations to 1e-15), bisection fallback; integration run
  as a 3-regime automaton (On/Medium/Off) so no step ever smears across a
  switching layer. RK7(8), RelTol=AbsTol=1e-14.
- **Never solve anything cold.** Double continuation: thrust ladder
  10 N → 0.3 N for the min-time and min-energy families (each step seeded by
  the last), then Bertrand–Épénoy ε-march 1 → 0 at the target thrust with a
  quadratically clustered schedule ε_j=(j²−1)/(N²−1), N=10.
- **Their stated knob when convergence fails: move t_f** (c_tf 1.1 → 1.6 as
  thrust drops). "Varying c_tf is a practical method to overcome convergence
  problems."
- No costate normalization, no Sundman, no scaling tricks — robustness comes
  entirely from exact derivatives + exact switch location + tight tolerances
  + continuation.
- Published converged λ_i values for all cases (Table 5) — useful sanity
  scale: λ_r,i = O(10), λ_v,i = O(0.03), λ_m,i ≈ 0.12.

**Key reconciliation with ms_band:** their smoothed-ε integration never hits
the 1/ε wall that killed the ms_band eps-march because the Medium band is
integrated as its own regime with exact band-boundary detection — the
integrator and STM are always smooth-per-arc. The ms_band failure integrated
the smoothed EOM blindly. IFS's hard-throttle explicit-node scheme is the
ε=0 endpoint of their automaton with nodes instead of events.

### Leomanni, Bianchini, Garulli, Quartullo 2021 (arXiv:2101.08160, JSR)

PDF: [`../../min_fuel_papers/Leomanni2021_low_thrust_orbit_transfers_made_easy.pdf`](../../min_fuel_papers/Leomanni2021_low_thrust_orbit_transfers_made_easy.pdf)
("Optimal Low-Thrust Orbit Transfers Made Easy: A Direct Approach")

"Made easy" = **made direct**: GPOPS-II hp-adaptive collocation + IPOPT +
ADiGator, no costates anywhere, no assumed switch structure (throttle is a
free bounded control; the hp mesh sharpens the switches). Plus:

- **Ideal elements + ideal-anomaly Sundman-type regularization** (Hansen):
  slowly-varying nonsingular states, independent variable whose rate is
  perturbation-independent, sparse input matrix. Two-body+J2 only — the
  element set presupposes an osculating conic, so it dies near the Moon /
  tulip endgame. Not our coordinate path.
- 3-parameter Lyapunov (Q-law-like) initial guess; 405-rev, 869-switch
  min-fuel solved on a laptop in ~3–5 h.
- Transferable hygiene: keep stiff smoothed terms out of the cost integrand;
  q = w/‖w‖ with convex ‖w‖²≤1 instead of ‖q‖=1; event-aware initial mesh
  densification (literally the PSR idea).
- **Verdict for us: it's a vote for PSR, not a fix for IFS.** Costate
  conditioning, shallow crossings, λ_r observability — those pathologies
  don't exist in their formulation because there are no costates.

### The synthesis

IFS's formulation (hard throttle, explicit switch nodes, no saltation
needed) is sound — externally confirmed, unit-proven. What IFS is missing
relative to Zhang 2015 is **everything around the solve**: continuation into
the answer instead of one cold 400-iteration LM run, and a solver that steps
*through* weak directions instead of crawling against them. Our own
diagnosis agrees: residual GN-consistent, null direction diffuse and
shifting → no surgical fix exists; the cure is seed quality + step strategy.

## 3. The rungs

### Rung 0 — a missing easy gate + a better seed (cheap, do first)

**FILE-FACT CORRECTION (verified 2026-07-11 by inspecting the `.mat`s):** the
seed files named in the first draft do not fit `ifs_seed`'s required layout
(`out.X[8xN]`, `out.U[4xN]`, `out.lamDef[8x(N-1)]`, top-level `factor`;
consumed via `sms_seed_duals`). Ground truth:
- `../sundman_minfuel/minfuel_from_energy_seed.mat` — this is the genuine
  **3-switch** local optimum, BUT it is the OLD 7-state fmincon format
  (`nlp.X[7x4001]`, `nlp.U`, `nlp.eqnonlin`; no `out.*`, no `lamDef`, no
  `factor`). `ifs_seed` **cannot read it** as-is.
- `../sundman_minfuel/results/solve_minfuel_f1150.mat` — correct `out.X/lamDef`
  layout, BUT it is **14 switches** (not 3) and has **no top-level `factor`**
  (`sms_seed_duals` calls `sms_problem(S.factor,…)` → would error).
- The only turnkey `ifs_seed`-compatible bang-bang rendezvous file is
  `../sundman_minfuel/results/minfuel/legacy_ms_f1120.mat` (has `factor=1.12`,
  12 switches) — i.e. the *existing* 1.12× gate. There is **no ready
  few-switch compatible file**.
- `../sundman_minfuel/refine/refine_history_smoke_1p12.mat` is 607 B — it
  carries only per-round metadata + `tauSwitch[1x12]`, **not** a refined
  trajectory (no X/U/lamDef). It cannot seed IFS directly.

**Revised, actually-buildable Rung 0:**
- **Mint an easy few-switch rendezvous gate** by a one-shot prep re-solve
  (the `prep_refine_seed` pattern; CasADi/IPOPT confirmed loadable): re-solve
  a low-switch 1.15× case at `eps=0, warmTight` to regenerate `out.X(8)/
  lamDef/factor`, giving a k≈3 full-rendezvous seed IFS can read. Converges →
  machinery validated end-to-end on a real rendezvous problem. **Caveat (be
  honest):** fewer switches removes the terminal-cluster issue but NOT the
  long weakly-observable first arc — the scaled-SVD investigation found λ_r0
  was the *dominant* null direction mid-crawl, and that arc exists at k=3 too.
  So a 3-switch gate is a real machinery test, not necessarily an "easy" solve.
- **Better 1.12× seed = inject PSR's refined `tauSwitch`** into the
  `legacy_ms_f1120` seed's switch times (small `ifs_seed` change to accept an
  override), OR re-run PSR persisting the full refined solution. The 607 B
  history file alone is not enough.
- **Reality check on Rung 0 in isolation:** a better seed does not fix the
  conditioning (this plan's own thesis). Rung 0 delivers a *validated easier
  gate + cleaner seed*; the actual convergence test needs Rung 1's solver
  change. Run them together, not Rung 0 alone expecting convergence.

### Rung 1 — the conditioning package (one unit of work)

1. **Canonical scaling**: physical state scales, reciprocal-costate column
   scales, per-block row scales (continuity vs terminal vs switch). All
   conditioning judgments on the scaled J from here on. (RESULTS.md lever 1,
   GPT rec 4.)
2. **Custom damped-Newton/LM loop with truncated-SVD (rank-revealing) step**
   on the scaled Jacobian + line search, replacing `lsqnonlin`. J is only
   178×178 — an SVD per iteration is free. Direct answer to "GN-consistent
   but crawling": step freely through the near-null λ_r direction instead of
   letting LM damping fight it.
3. **Retire the stick-breaking sigmoid** (`ifs_taus.m`): dτ/dg → 0 as a gap
   closes — the parameterization degenerates exactly where the problem is
   hardest (GPT red flag). Use switch times directly (the zero-span guard in
   `ifs_int_arc` already exists) or gap variables with a floor.
4. Optional, high-value: swap per-arc complex-step (~25 ode113 calls/arc)
   for **per-arc variational-equation STMs** à la Zhang. ~10–20× cheaper
   Jacobians turn the 30-min/400-iteration budget into thousands of
   iterations.

### Rung 2 — continuation in t_f (Zhang's medicine, and the real prize)

Once *any* full-rendezvous IFS point converges (3-switch 1.15×, or 1.12×
post-Rung-1): continue in t_f in small steps, re-seeding from the converged
(λ_0, nodes, τ_i); detect folds (S = Ṡ = 0) to add/delete switch pairs;
pseudo-arclength near folds. This is simultaneously RESULTS.md lever 2 and
the tool the campaign already identified for the **open 1.01–1.11×
transition band**: min-time (1.00×, k=0) already converges indirectly, so
IFS + t_f-continuation marching **up** from min-time is the
bifurcation-tracking instrument the band needs. IFS stops being just a
finisher and becomes the band-mapper.

Fallback within this rung: ε-continuation done *right* (Zhang's regime-aware
3-arc automaton with exact band-crossing detection — not the blind smoothed
integration that killed ms_band). Heavier build; only if t_f-continuation
stalls.

### Rung 3 — structural (only if Rungs 0–2 stall)

- **Interior non-switch shooting nodes** inside multi-rev arcs (at k=10 the
  arcs span several revolutions each — per-arc sensitivity is genuinely
  large). RESULTS.md lever 3. GPT caution: helps long-arc STM conditioning,
  does not cure vanishing-arc rank loss (which was refuted anyway).
- Limit case: **indirect collocation** — solve the 16-dim PMP system + S=0
  switch conditions as collocation defects on the PSR mesh, seeded by direct
  states + dual-mapped costates. The direct solver's machine-tight
  convergence on identical dynamics at 4001 nodes is proof the collocated
  Jacobian is tractable where the shooting one isn't.

### Cheap formulation audits along the way (both from the GPT review)

- **Free τ_f** via σ ∈ [0,1] with an unknown Sundman-length scale, so the
  certificate applies to the true fixed-t_f problem, not the
  fixed-Sundman-length one. (The empirical residual-split test ruled out
  over-constraint as the *stall* cause; the *equivalence* question for the
  certificate remains.)
- Confirm the "1" in S and λ_m(τ_f)=0 come from one consistent fuel-objective
  normalization.

## 3b. Chosen first increment: Rung 0 + Rung 1 together (2026-07-11)

Decision: build and run Rung 0's easy gate and Rung 1's solver as **one
increment**, because Rung 0 alone cannot converge (a better seed does not fix
conditioning) and the actual open question is "does IFS converge a real
rendezvous case once the step strategy is fixed?" Concrete pieces:

1. **`ifs_solve2.m`** — a custom damped Gauss–Newton solver replacing
   `lsqnonlin`, keeping the same `(Z0, prob, opts)` signature and the existing
   `ifs_residual` analytic Jacobian. Per iterate: column scaling
   (reciprocal-magnitude of the current iterate, with a floor) + Jacobian-row
   equilibration; **SVD of the scaled Jacobian**; a **truncated / rank-revealing
   Gauss–Newton step** (singular values below `relTol·σ_max` are dropped, so the
   step moves fully along well-determined directions and does NOT crawl against
   the near-null λ_r0 direction); Levenberg damping as a fallback when the
   truncated step fails to reduce ‖R‖; backtracking line search. Records the
   singular-value spectrum / numerical rank per iterate for diagnosis.
2. **Direct-τ parameterization** — `prob.tauParam = 'direct'` makes the switch
   times themselves the unknowns (monotonicity via a min-gap projection in the
   solver), retiring the stick-breaking sigmoid whose dτ/dg→0 compounds the
   crawl. `'sigmoid'` stays the default so the existing unit tests are
   unchanged; `'direct'` is toggled on for the gate runs.
3. **Minted 3-switch easy gate** — the genuine 3-switch local optimum lives
   only in the old 7-state `minfuel_from_energy_seed.mat`, so mint an
   `ifs_seed`-compatible seed by mapping it to Sundman 8-state
   (`sundman_seed_map`) and an `eps=0, warmTight` re-solve
   (`casadi_minfuel_sundman`, confirmed loadable) to regenerate
   `out.X(8)/lamDef/factor`. Saved as `ifs/seed_3sw_1p15.mat`.
4. **Run on both gates** — `ifs_solve2` on the minted 3-switch gate (machinery
   checkout on a real rendezvous problem) and on the existing 1.12× gate
   (`legacy_ms_f1120.mat`), comparing the ‖R‖ trajectory against the old
   lsqnonlin crawl (1.96 → 0.023, stalled). Certify anything that converges.

### 3b-outcome (2026-07-11) — see `RESULTS_RUNG01_RUNG2.md`

Executed. Summary: (0) no clean easy gate exists — 1.12x (k=10) is the smallest
compatible system. (1) `ifs_solve2` (scaled truncated-SVD + adaptive truncation
continuation) **descends** the cold 1.12x seed 1.96 -> ~0.43 (vs frozen/1.8-stall
of earlier cuts) but does **not** converge — the cold 40-rev shooting basin wall,
as predicted. (2) Rung 2: the **min-time k=0 anchor converges** (IFS's first
end-to-end rendezvous solve, ||R||~2e-7, S<0 with margin), but naive t_f-stepping
**fails at the min-time fold** (0.1% step already breaks; vertical-tangent
degeneracy). (2b) **Pseudo-arclength** (`ifs_tf_arclength.m`) was built and is
mechanically correct (correctors converge, extended 8x9 Jacobian full rank), but
the min-time anchor carries a **near-null costate gauge** (scaled σ_min~8e-4) on
a vertical branch, so the march wanders the fixed-t_f min-time manifold (max S
non-monotone, factor frozen at 1.0000) — the reported "birth" is a gauge
artifact, not physical. **Next lever: gauge regularization** of the anchor
(phase/pinning on the near-null λ_0 direction), or a non-degenerate anchor above
the fold. Full record: `RESULTS_RUNG01_RUNG2.md`.

## 4. Odds and stopping rule

Rung 0+1 alone probably gets the 3-switch gate and maybe 1.12×; the
terminal-cluster / many-switch cases will need Rung 2. Stop short of Rung 3
unless t_f-continuation demonstrably stalls *between* converged points — it
is the biggest build for the least certain payoff, and PSR already covers
the "good enough switch times" use case.
