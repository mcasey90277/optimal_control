# Rung B — smoothed-indirect anchor + early-lock handoff to IFS (2026-07-12)

Detailed, self-contained execution plan for Rung B of `PLAN_OF_ATTACK_2.md`.
Written for a future session (AI or human) to pick up cold. Read
`PLAN_OF_ATTACK_2.md` §1 (ground truth) first; this doc assumes it.

## 0. Why Rung B, in one paragraph

Rung A (adjoint sweep/smoother, `ifs_seed_adjoint.m`) proved that **no seed
built from the direct solution's mesh-accuracy data can be dynamically
consistent** — the ~1% dual costates are inconsistent with any exact adjoint
trajectory at the 40-rev amplification level (1.2e12), and the IFS solver
crawls from every such seed (best: 0.224 at 60 iters from the smoother-blended
seed, same cond~1e16 crawl). The only source of *dynamically consistent*
costates is **a converged indirect solve of an easier problem**. The smoothed
(finite-ε) min-fuel problem is that easier problem: its extremal is smooth, its
basin is wide, and as ε shrinks its solution approaches the bang-bang extremal.
Rung B: converge the smoothed indirect problem at a smooth anchor, march ε
down **only until the throttle saturates arc-wise** (`ε_lock`), then hand the
states+costates+switch structure to hard-throttle IFS to finish. The march
never approaches ε→0, so the 1/ε-layer crawl that killed ms_band is never
entered.

## 1. Critical history — ms_band Gate D is NOT this plan (but read it first)

`../ms_band/MS_BAND_CAMPAIGN.md` (2026-07-10 entries) records that an ε-march
at 1.12× was tried and BLOCKED. Rung B is a different experiment in three
load-bearing ways:

1. **Seed–problem structure match.** Gate D seeded the ε=1 (smooth) solve with
   the **min-FUEL dual seed** (near-bang costates). Its own diagnosis: "the
   eps=1 extremal is far from the near-bang seed... switches 0, bang 40.7%" —
   a structure mismatch, LM crawling across it at ~1e-4 nats/iter. Rung B
   seeds the ε=1 solve with **min-ENERGY duals** (smooth extremal ↔ smooth
   seed). This cell of the matrix has never been tried.
2. **Solver.** Gate D used `ms_solve` (lsqnonlin LM). Its trace showed a
   healthy system (GN-consistent to 6e-12, full rank, cond 6.3e8) with GN
   steps far outside the linear regime — the classic case where the
   equilibrated truncated-SVD stepping of `ifs_solve2` (built AFTER ms_band,
   Rung 1) does better than LM damping.
3. **Goal.** Gate D marched toward ε→1e-4 (the full crawl). Rung B stops at
   `ε_lock` (arc saturation) and lets hard-throttle IFS do the sharp part —
   IFS is exactly the machine that has no ε anywhere.

Also inherit from ms_band: the **switch-displacement guard discipline** of
`eps_march_adaptive.m` (discard steps that move S-crossings too far / change
their count; geometric bisection toward the last accepted ε). This is the
protection against the Rung-0 "mint cliff" (a free re-solve sliding into the
many-switch global basin, 6.35→5.4).

## 2. The smoothing family — precise definitions (do not confuse the two ε's)

**The indirect side (`sms_eom`) uses ENTROPY smoothing**, not the direct
solver's objective homotopy:

- Throttle law: `u = logistic(-S/ε) = (1 - tanh(S/(2ε)))/2`, S the min-fuel
  switching function `S = 1 - ||λv||c/m - λm`.
- Time-domain Hamiltonian includes the entropy running cost
  `Ht = λr'v + λv'(g+h) + (Tmax/c)(u·S + ε·Lear)`, `Lear = u·log u +
  (1-u)·log(1-u)` — required for H-conservation; implemented CS-safe
  (softplus identity) in `sms_eom.m`. Read its header before touching it.
- ε→∞: u→1/2 everywhere. ε→0: hard bang-bang. ε=1 is "smooth" at this
  problem's S-scale (Gate D's smooth end).

**The direct side (`casadi_minfuel_sundman`) uses the Bertrand–Épénoy
objective homotopy** `J(eps) = ∫s dt - eps·∫s(1-s) dt` (eps=1 energy → eps=0
fuel), whose interior throttle law is the LINEAR ramp `u* = clamp(1/2 -
S/(2·eps'), 0, 1)` — a different family with saturation kinks. Rung B uses
the direct eps only once (B0: produce the energy direct solution); all
indirect marching is in the entropy-ε family. Mapping between them is loose
(both interpolate energy-like ↔ bang-bang); do NOT assume the ε values
correspond numerically.

**Lock criterion (define before running, don't improvise):** at a converged
ε-step, sample u(σ) on the fixed 4000-point σ grid; identify S-crossings;
exclude a transition window of half-width `w = ε / min|dS/dτ|_crossings`
(the layer width) around each crossing; LOCKED when outside all windows
u ∈ [0, 0.05] ∪ [0.95, 1] uniformly, AND the crossing count and locations
(grid-node displacement < the eps_march_adaptive DISP_TOL) are stable across
two consecutive accepted ε-steps. Record `ε_lock` and the crossing set.

## 3. The β subtlety for smooth solutions (a trap for the implementer)

`sms_seed_duals` fits the dual→costate scale β via `beta_from_duals`, which is
a **switching-law fit assuming near-bang arcs** (classify burn/coast, fit S
signs). On a SMOOTH energy solution there are no bang arcs — that fit will be
garbage or error. For the B1 energy-dual seed, fit β instead from the smooth
throttle law: the direct energy solution's throttle s(σ) and the candidate
costates' S(σ; β) must satisfy the DIRECT family's ramp law
`s ≈ clamp(1/2 - S/(2·eps'), 0, 1)` with eps'=1 — β enters through S's scale,
so a 1-D fit of β (e.g. least squares of s_direct vs the ramp prediction on
interior nodes 0.05 < s < 0.95) is robust. ~20 lines; write it as
`beta_from_duals_smooth.m` next to `beta_from_duals.m`. (Alternatively verify
the mode-'d' theory value β ≈ 1 — the h-cancellation argument in
`sms_seed_duals`'s header applies to any integrand — but MEASURE it; at 1.85×
the fuel-side fit failed 4 orders, so never trust an unfitted β.)
Note: fuel-side β at 1.12× is 0.03102 (spread 0.45%) — if the smooth fit
lands orders away from that, something is wrong.

## 4. Assets on disk (verified 2026-07-12) and what must be built

| asset | status |
|---|---|
| Smooth ENERGY direct solution at any factor | **DOES NOT EXIST on disk.** The `results/minfuel/minfuel_f*_en.mat` files are NOT energy solutions — 'en' labels the energy-seeded *branch*; f1300_en has fracEdge=0.99 (bang-bang). B0 must create one. |
| Fuel direct solutions | `legacy_ms_f1120.mat` (k=10 certified; raw throttle crossings 12 = artifact — always use dual-S crossings), f1140, f1250, f1850; f13..f18 series. f1120/f1140 carry full `out.lamAll` (terminal-BC duals at tail: rv 6 + tf 1; layout `8N | nN | 24nN | 8 | 7` verified). |
| Energy re-solve machinery | `../sundman_minfuel/energy_step.m` + `casadi_minfuel_sundman(..., epsilon=1, warmTight=false)` — the loose-warm-start IPOPT path exists and is documented in the solver's comments ("the eps=1 energy re-solve"). |
| Smoothed indirect solver | `../ms_band/sms_{problem,eom,residual,pack,unpack,seed_duals}.m`, `ms_solve.m` (LM), `eps_march_adaptive.m` (guarded march + resume), `sms_jacobian_cs.m`. Campaign gates A–C green (EOM/Jacobian/min-time anchor). |
| Truncated-SVD solver | `ifs_solve2.m` — **hardcodes `ifs_residual`** (lines 69/105/139). B1 needs a small genericization: use `prob.resFun(Z, prob)` when the field exists (sms_problem already sets `prob.resFun = @sms_residual`), default to `ifs_residual`. Keep the API otherwise identical; rerun the ifs unit tests after. |
| Hard-throttle finisher | `ifs_residual/ifs_solve2/ifs_certify` + `ifs_seed`'s Z/prob assembly pattern (mirror it in B3). |

## 5. Increments and gates

### B0 — mint the energy anchor's direct solution (cheap, mechanical)
Produce a converged eps=1 (energy) direct solution at the target factor and
save it WITH `lamAll` (same save layout as `direct_build_minfuel`):
- Target factor: **1.12×** (the gate; certified structure to compare against).
  Warm start: `legacy_ms_f1120.mat`'s X/U, `casadi_minfuel_sundman(...,
  epsilon=1, warmTight=false)` (via `energy_step.m` if its interface fits).
- Gate B0: IPOPT optimal; throttle SMOOTH (fracEdge well below ~0.5;
  no bang structure); save `results/minfuel/energy_f1120.mat`.
- Sanity row: `lamMassEnd` ≈ 0 still (free final mass); `lamAll` tail present.

### B1 — THE gating experiment: converge sms at ε=1 from energy duals
- Seed: `sms_seed_duals('energy_f1120.mat', M, 1, 'd')` with β from
  `beta_from_duals_smooth` (§3). M: start 40 (Gate-D used 40/50).
- Solve at ε=1 with BOTH solvers (compare): (a) `ms_solve` (LM baseline),
  (b) genericized `ifs_solve2` driving `sms_residual`.
- **Gate B1: ||R|| ≤ 1e-9** (the campaign's certified-solve gate; odeOpts
  RelTol 1e-13/AbsTol 1e-15 — floor analysis says 1e-9 is honest there).
  Fallback gate: a documented floor < 1e-6 with a trace capture
  (diag_s1_gateD-style: is R in range(J)? cond? step size vs linear regime?).
- If B1 FAILS from a structure-matched seed with the better solver: capture
  the trace, compare to Gate D's, and STOP — that is strong evidence the
  smoothed-indirect basin at 40 revs is itself out of reach of mesh-accuracy
  seeds, and the campaign's "regularized coordinates" thesis takes over.
  Do not grind iterations past a diagnosed crawl.

### B2 — guarded ε-march down to ε_lock (not to ε→0)
- Reuse/extend `eps_march_adaptive` (it already has: geometric step targeting,
  switch-displacement guards, discard-and-bisect, resume-from-state). Change
  the TERMINATION: stop on the §2 lock criterion, not on reaching epsFloor.
- Track per accepted step: ε, ||R||, iterations, crossing count/locations,
  saturation fraction, min layer width. The interesting curve is
  "iterations-to-converge vs ε" — if it blows up before lock, record where:
  that measures the ε_lock-vs-ε_crawl gap, which is THE quantity this rung
  exists to discover.
- Gate B2: locked structure with crossing count == the direct solution's
  certified count (10 at 1.12× by dual-S; if it locks at a different count,
  adjudicate against `diag_verify_1120`-style recount before proceeding —
  do NOT assume the direct count is right).

### B3 — handoff to hard-throttle IFS (the finish)
- From the ε_lock converged solution: λ0 = its initial costates; switch times
  = its S-crossings (NOT the u=0.5 crossings — at finite ε they differ by
  O(ε)); node states+costates sampled AT those crossings; uArc from sign(S)
  per arc. Assemble Z/prob exactly as `ifs_seed.m` 'full' does (rendezvous
  terminal, tauParam sigmoid default). ~60 lines; suggest
  `ifs_seed_smoothed.m` mirroring `ifs_seed.m`'s interface.
- Run `ifs_solve2` (tolR 1e-8), then `ifs_certify`.
- **Gate B3: converged + certified = the first full IFS min-fuel convergence.**
  Watch iteration count: from a dynamically consistent seed the crawl should
  be GONE (5–20 iters per sol's estimate); if it still crawls, the wall was
  never seed consistency and the whole indirect program needs the
  regularized-coordinates rebuild (record and stop).

### B4 — replicate + band entry
- Repeat B0–B3 at 1.14× (second data point, k=24). Then Rung C
  (t_f-continuation between converged IFS points) attacks 1.01–1.11×.

## 6. Budgets, risks, stopping rule

- B0: one IPOPT solve (~minutes–1h). B1: hours (two solver runs + diagnosis).
  B2: the long pole — Gate-D steps ran ~10–70 min each; expect O(10) accepted
  steps; use `eps_march_adaptive`'s resume/state discipline (an external
  ~70-min watchdog killed three Gate-D runs; plan for background + resume).
  B3: ~1h. Total: a solid multi-session campaign; every increment has a
  falsifiable gate and independent diagnostic value.
- **Risk 1 — B1 fails**: pre-diagnosed stop (see B1). This is the cheapest
  possible test of the rung's core premise; run it FIRST, build nothing
  speculative before it passes.
- **Risk 2 — crawl arrives before lock (B2)**: the measured gap is itself the
  deliverable; if ε_lock unreachable, stop and pivot (regularized coords,
  `../LOW_THRUST_MINFUEL_CAMPAIGN.md` thesis; Leomanni 2021 digest in
  `PLAN_OF_ATTACK.md` §2).
- **Risk 3 — structure slide at lock** (the mint cliff): guarded by the
  displacement discipline; if the locked count ≠ certified count and the
  recount upholds the certified one, treat as Risk-2 evidence.
- **Risk 4 — β/scale bugs**: §3. Symptoms: seed ||R|| orders above ~5
  (Gate D's structure-mismatched seed was 4.285 — a structure-MATCHED seed
  should start lower), or S-scale inconsistent with u_direct.
- **Stopping rule**: if B1, or B2's lock, or B3's finish fails for the
  *diagnosed structural reason* (not a bug), Rung B is falsified; the honest
  conclusion is the standing campaign thesis (sharp many-switch indirect at
  40 revs needs regularized coordinates) — a new build, decided deliberately,
  not a patch.

## 7. Pointers

- Parent plan + Rung A falsification: `PLAN_OF_ATTACK_2.md`,
  `RESULTS_RUNG01_RUNG2.md` (§Rung A: why only a converged indirect solve can
  supply consistent costates).
- ms_band record (Gate D traces, guard discipline, dual-map adjudication):
  `../ms_band/MS_BAND_CAMPAIGN.md`; diagnostics `../ms_band/diag_s1_gateD*.m`.
- Consult that proposed the rung (and whose Q4 sweep prediction Rung A
  refuted — calibrate trust accordingly): `CONSULT_GPT56_response.md` Q2.
- Direct-side homotopy reference: Bertrand & Épénoy, OCAM 23(4), 2002;
  Zhang et al., JGCD 38(8), 2015 (`../min_fuel_papers/`).
- Campaign context / two-walls analysis: `../LOW_THRUST_MINFUEL_CAMPAIGN.md`.
