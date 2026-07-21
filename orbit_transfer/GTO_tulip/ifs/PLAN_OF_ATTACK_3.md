# IFS Plan of Attack 3 — the Zhang audit and the two-prong plan (2026-07-12)

> **STATUS UPDATE (2026-07-12, same day): user decision — go straight to
> Prong Z; Prong B is SHELVED (doc kept as fallback). Detailed Prong-Z
> execution plan: `PLAN_PRONG_Z.md`. Note its §1 correction to §4 below:
> Zhang ladders the SMOOTH families (min-time, min-energy) in thrust and
> runs the ε-march ONCE at target thrust — the increments in PLAN_PRONG_Z
> supersede the Z0–Z5 sketch here.**

Supersedes the rung sequence of `PLAN_OF_ATTACK_2.md` (whose Rung A was
executed and FALSIFIED; see `RESULTS_RUNG01_RUNG2.md` §"Rung A"). Rung B of
that plan survives here unchanged as Prong B (`PLAN_RUNG_B.md` is still the
execution doc). What is new is Prong Z: a faithful replication of the Zhang
2015 recipe, motivated by the audit in §2 — we have never actually run it,
only subsets of it.

## 1. Ground truth (what is now proven)

1. **PSR works and is the deliverable.** Certified bang-bang direct solutions
   across [1.12, 1.95]x with IPOPT native local-min certificates, exported
   data products (mesh/control/trajectory/dual-mapped costates), batch + movie
   infrastructure. Energy backbones exist exactly on [1.12, 1.95]x
   (`energy-seed-band` memory); the 1.01–1.11x band is hard even for the
   smooth energy problem.
2. **IFS machinery is validated; the solver is not the bottleneck.**
   `ifs_solve2` (equilibrated truncated-SVD GN) descends everywhere LM crawled.
   EOM/residual/Jacobian/seed unit tests green.
3. **The seed-quality axis is CLOSED (Rung A, falsified conclusively).**
   Mesh-accuracy (~1%) KKT duals are dynamically inconsistent with any exact
   adjoint trajectory at the 40-rev amplification level (1.2e12). No seed
   *derived from direct-solution data* can converge IFS: the 1.12x stall
   floor ~2–4e-2 is seed-independent.
4. **The min-time anchor is dead** (manifold tangency at t_f,min; gauge
   pinning fixes the wrong direction; pseudo-arclength crawls the degenerate
   costate family at factor=1.0000). Machinery (`ifs_tf_arclength`) is kept —
   it needs a non-degenerate anchor, which this plan supplies.
5. **Corollary (the load-bearing sentence):** the only source of dynamically
   consistent costates is **a converged indirect solve of an easier problem**.
   Exactly two "easier" axes exist for this problem:
   - **epsilon** (smoothing): the smoothed extremal has a wide basin → Prong B.
   - **thrust** (Tmax): at high thrust the transfer is a few revolutions, the
     shooting basin is large, and amplification is tame → Prong Z.

## 2. The Zhang audit — why Prong Z exists

Zhang, Topputo, Bernelli-Zazzera, Zhao, JGCD 38(8), 2015
(`../min_fuel_papers/2015-8-J-LowThrustMinimumFuelOptimizationRTBP.pdf`;
digest in `PLAN_OF_ATTACK.md` §2) is the **existence proof for exactly our
problem class**: Earth–Moon CR3BP, indirect min-fuel, GTO departure. Their
hardest converged case — 0.6 N, 140 days, **~150 bang-bang switches over ~150
revolutions, single shooting with 7 unknowns** — is ~4x beyond our target
(25 mN, ~40 revs, ~25 switches at 1.15x). Their robustness comes from four
ingredients used TOGETHER:

- **(a) Exact derivatives**: analytic variational-equation STMs (14+196 ODEs)
  with a saltation matrix Ψ = I + (ẏ⁺−ẏ⁻)(∂S/∂y)/Ṡ composed at every switch.
- **(b) Exact switch handling**: machine-precision hybrid event detection
  (Newton on S with analytic Ṡ, 1e-15; bisection fallback), integration as a
  3-regime automaton so no step smears a switching layer; RelTol=AbsTol=1e-14.
- **(c) Never solve cold**: double continuation — thrust ladder 10 N → 0.3 N
  (min-time and min-energy families), then Bertrand–Épénoy ε-march 1 → 0 at
  target thrust with the quadratically clustered schedule ε_j=(j²−1)/(N²−1).
- **(d) The escape knob**: "varying c_tf is a practical method to overcome
  convergence problems" — move t_f when a solve resists.

Known weakness (theirs and ours): **grazing switches** (Ψ's 1/Ṡ) — the exact
pathology of our 1.01–1.11x band edge.

**Audit of every failed attempt in this campaign against (a)–(d):**

| attempt | (a) exact STM+saltation | (b) exact switches | (c) continuation | (d) t_f knob | outcome |
|---|---|---|---|---|---|
| old thrust ladder (`../../lowThrust_GTO_tulip/thrust_continuation_minfuel_indirect.m`) | no (CS through ode113) | no (smeared; switch counts 24→674 erratic) | ladder yes, but 37% steps and the TOP ANCHOR itself never converged (res 6.4e83); fine 7.4% retry also 0-for-5 | fixed 1.15 | failed every rung |
| ms_band ε-march (Gate D) | no | smoothed (finite-ε layers) | ε only, marched into the 1/ε crawl; no ladder, no anchor | no | crawl, blocked |
| IFS cold/dual/adjoint seeds | CS Jacobian (accurate but 10–20x costlier) | yes (explicit nodes — the one ingredient we do have) | **none** | no | floors 0.02–0.4 |
| pseudo-arclength (Rung 2b) | — | — | yes, but off the degenerate min-time anchor | t_f is the parameter | crawls the gauge family |

**Conclusion: the recipe has never been run whole.** Each failure used a
strict subset of the ingredients. Zhang says our target sits comfortably
inside the solvable envelope; the audit says our failures are explained by
the missing ingredients, not by the problem being harder than theirs. That is
a testable claim, and Prong Z tests it.

## 3. Prong B — smoothed-indirect anchor + early lock (run FIRST, cheapest)

Execute `PLAN_RUNG_B.md` as written: B0 (mint energy direct solution at 1.12x
with lamAll) → **B1, the gating experiment** (converge sms at ε=1 from
energy duals — the structure-matched cell of the matrix that Gate D never
tried) → B2 (guarded ε-march to ε_lock, NOT to ε→0) → B3 (hard-throttle IFS
finish) → B4 (replicate at 1.14x).

Two Zhang-lens amendments to B2:
1. Use the **quadratically clustered ε schedule** ε_j=(j²−1)/(N²−1) between
   ε=1 and the lock check, rather than ad-hoc geometric steps.
2. Adopt knob (d): if the march stalls at some ε with the displacement guards
   healthy, try a small t_f sidestep (±0.005–0.01 in factor) at fixed ε before
   declaring the stall structural.

Cost: all code exists; B0+B1 is days. B1's stop rule stands: if a
structure-matched seed with the better solver can't converge even the SMOOTH
ε=1 indirect problem at 40 revs, record the trace and stop the prong.

## 4. Prong Z — ZTL, the Zhang-style thrust ladder done right

A single-shooting min-fuel solver with ALL FOUR ingredients, laddered in
thrust from an easy anchor down to 25 mN at fixed factor (1.15x to start —
PSR's certified solution there is the oracle). Suggested home: `../ztl/`.
Increments, each with a falsifiable gate:

- **Z0 — exact variational STM for hard-throttle arcs.** Time-domain
  Cartesian first (Zhang fidelity; Sundman variant only if high-thrust perigee
  sensitivity bites — at 200 mN the perigee problem that motivated Sundman is
  far milder). States+costates 14 (or 16 with (t, λ_t) if Sundman), plus the
  n² variational ODEs, hard u ∈ {0,1} per arc.
  *Gate Z0:* per-arc STM matches complex-step to rel ~1e-8, at 10–20x lower
  cost. **Side payoff regardless of Z's fate: drop-in cheaper Jacobians for
  `ifs_residual`** (the old plan's Rung-1 lever, never built).
- **Z1 — event-detected switching + saltation + 3-regime automaton.** Newton
  on S(t)=0 with analytic Ṡ to 1e-15 (bisection fallback); integrate
  On/Medium/Off so no step crosses a switch; compose Ψ at each event.
  *Gate Z1:* full-span STM across ≥1 switch matches CS to ~1e-7; switch times
  reproducible to 1e-13 under tolerance changes.
- **Z2 — single-shooting residual.** Unknowns λ0 ∈ R^7 (Cartesian: λ_rv 6 +
  λ_m; t_f fixed = factor × t_f,min(Tmax)). Residual: terminal rendezvous
  rv (6) + λ_m(t_f)=0 (1). Solver: `fsolve` with the exact Jacobian from
  Z0+Z1; fall back to the equilibrated truncated-SVD stepping of `ifs_solve2`
  if damping is needed.
  *Gate Z2:* residual is zero (to integration tol) on a manufactured
  ground-truth arc, exactly as `test_ifs_residual` does for IFS.
- **Z3 — the top anchor.** At high thrust (start ~200 mN ≈ 8x, where
  t_f,min ≈ 6.07 and the transfer is a few revs): min-time indirect (already
  converges — `run_gto_tulip_indirect` machinery) → Bertrand–Épénoy ε-march
  1→0 with the quadratic schedule → hard bang-bang min-fuel at factor 1.15.
  This is precisely where the old ladder died (6.4e83); with (a)+(b) it is
  Zhang's EASIEST case, so it is the decisive test of the audit's claim.
  *Gate Z3:* converged (‖R‖ ≤ 1e-10) + sign-law certified at top thrust.
  **If Z3 fails with the full machinery verified, the audit's explanation is
  wrong** — record and reassess (that would be genuine evidence for the
  regularized-coordinates thesis).
- **Z4 — ladder down 200 → 25 mN.** Warm-start λ0 from the previous rung;
  steps ≤ 10% in Tmax (Zhang: ~8%; the old attempt's 37% was reckless), halve
  on failure; knob (d) allowed. Switch births/deaths are AUTOMATIC — event
  detection re-discovers the structure each solve; no structural surgery
  (this is single shooting's decisive advantage over explicit-switch MS for
  continuation).
  *Gate Z4 (the campaign prize):* converged at 25 mN, factor 1.15x, and the
  solution matches PSR's certified 1.15x product (ΔV, ~25 switches, switch
  times) — the first fully converged indirect min-fuel at nominal thrust.
- **Z5 — the band.** From the converged 25 mN point, t_f-continuation down
  from 1.15x toward 1.01x: naive stepping first, `ifs_tf_arclength`'s
  pseudo-arclength (now with a NON-degenerate anchor) near folds. Grazing
  switches (1/Ṡ) are expected at structure changes — detour in t_f around
  graze points, and record any hard graze wall honestly: that wall would be
  the indirect twin of the direct side's 1.01–1.11x energy-seed wall and is a
  publishable characterization by itself.

**Where IFS fits (nothing is wasted):** a converged ZTL point hands IFS a
dynamically consistent (λ0, switch times, arc structure) — the seed class
Rung A proved cannot be built from direct data. IFS then serves as (i) the
independent multiple-shooting verifier/certifier (`ifs_certify`), and (ii)
the robustifier if single shooting turns sensitive deep in the band (short
arcs + explicit nodes, now warm). Z0's exact STMs also retrofit into
`ifs_residual` for cheap Jacobians.

## 5. Ordering, decision tree, stopping rule

```
now ──> B0+B1 (days; all code exists)          Z0+Z1 (unit-testable; useful to
        │                                       IFS regardless — can start in
        │                                       parallel or immediately after)
        ├─ B1 passes → B2 (locked march) → B3 (IFS finish) → B4
        │              │ (Z continues as the band-mapper either way)
        │              └─ crawl before lock → measure the gap, stop prong B
        └─ B1 fails at the diagnosed stop → Prong Z is the main line
Z3 anchor:
        ├─ converges → Z4 ladder → gate vs PSR 1.15x → Z5 band
        └─ fails with (a)+(b) verified → the Zhang-gap explanation is wrong →
           regularized coordinates (Leomanni 2021 digest, PLAN_OF_ATTACK.md §2)
           becomes the standing conclusion — stop and decide deliberately.
```

- **Run B first** because it is nearly free (two working codebases, detailed
  plan) and its B1 gate is decisive either way.
- **Z is the strategic build** (~3 sessions to Z3): it is the published
  recipe for exactly this problem, it is the only path that handles the
  band's structure changes automatically, and its increments pay IFS
  dividends even on failure.
- **Both prongs failing at their diagnosed stops** = the campaign thesis
  (sharp many-switch indirect at 40 revs needs regularized coordinates) is
  upheld from two independent directions; that rebuild is then a deliberate
  new project, not a patch.

## 6. Odds (honest)

- B1: a genuinely untried cell (structure-matched smooth seed + the better
  solver). Coin flip; cheap either way.
- Z3 (top anchor): should be Zhang's easiest case — if the audit is right,
  high confidence; its failure is the single most informative negative
  available.
- Z4 (full ladder to 25 mN): Zhang converged 4x beyond this with the same
  recipe; the honest unknown is our tulip rendezvous terminal (they targeted
  a halo insertion point) and the fixed-factor choice. Moderate-to-good, with
  knob (d) in reserve.
- Z5 (band): unknown — the graze wall is real physics both sides have now
  hit. Even partial penetration (e.g. 1.05–1.11x) plus a characterized wall
  is a strong result.

## 7. Pointers

- Prong B execution doc: `PLAN_RUNG_B.md` (read its §1–§3 traps first).
- Rung A falsification + solver history: `RESULTS_RUNG01_RUNG2.md`.
- Zhang digest: `PLAN_OF_ATTACK.md` §2; PDF in `../min_fuel_papers/`.
- Old ladder post-mortem data: `../../lowThrust_GTO_tulip/
  thrust_continuation_results.mat` (+ `_fine`): top anchor res 6.4e83,
  erratic switch counts — the (a)/(b)-less control experiment for Z3.
- PSR oracle products: `../PSR_data/psr_data_tf*_minEps*.mat` (certified
  bang-bang + IPOPT local-min certificates + dual-mapped costates).
- Consults: `CONSULT_GPT56_response.md` (note: its Q4 sweep prediction was
  refuted by Rung A — calibrate trust; it never evaluated the single-shooting
  route).
