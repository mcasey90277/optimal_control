# IFS Rung 0+1 and Rung 2 — session results (2026-07-11)

Work done after `PLAN_OF_ATTACK.md`: built and ran the combined Rung 0 + Rung 1
increment, then Rung 2 (t_f-continuation). Honest record of what converged and
what didn't.

## New/changed code

| file | what |
|---|---|
| `ifs_solve2.m` | Rung-1 solver: two-sided Jacobian equilibration + **truncated / rank-revealing SVD Gauss-Newton step** with a Levenberg fallback and an alpha-floor line search; **adaptive truncation continuation** (start aggressive, relax `relTrunc` on plateau/stall). Same `(Z0,prob,opts)` API as `ifs_solve`. |
| `ifs_taus.m`,`ifs_gseed.m` | `mode='direct'` (switch times as unknowns) vs `'sigmoid'` (default). Retires the stick-breaking sigmoid whose `dtau/dg->0` compounds the crawl. |
| `ifs_residual.m` | thread `prob.tauParam`; **k=0 (all-burn) support** in `ifs_arcs` (single arc, no node/switch blocks). |
| `ifs_certify.m` | thread `prob.tauParam`. |
| `ifs_seed.m` | `opts.tauParam` -> builds the seed + prob in either parameterization. |
| `ifs_eom.m` | expose `Ht` (3rd output) for the `lamT0=-Ht(0)` anchor construction. |
| `ifs_seed_mintime.m` | **k=0 min-time anchor**: min-time indirect costates (`run_gto_tulip_indirect`) + `lamT0=-Ht(0)`, `tauf` by integrating hard-burn to `t=tf`. |
| `ifs_tf_continuation.m` | Rung-2 driver: march factor up from the anchor, k=0, track `max S` toward the first switch birth. |
| `mint_easy_gate.m` | attempt to mint a 3-switch easy gate (see finding below). |

Unit tests unchanged and still green (sigmoid default preserved); direct-mode
seed reproduces the sigmoid seed residual exactly (1.960465 both).

## Rung 0 — no clean easy gate exists (finding)

- The genuine 3-switch local optimum lives only in the old 7-state
  `minfuel_from_energy_seed.mat` (no costates/`factor`); minting it via an
  `eps=0 warmTight` Sundman re-solve **does not preserve 3 switches** — the
  objective slides 6.35 -> 5.4+ toward the many-switch global basin (the
  3-switch point is not a fixed point of the eps=0 solve, unlike the certified
  global solutions `prep_refine_seed` handles).
- Every compatible legacy solution is many-switch; **1.12x (k=10) is in fact the
  smallest** on disk (others 26/50/22). So there is no easier real gate — the
  1.12x case is both the target and the smallest system.

## Rung 1 — the solver descends the cold 1.12x seed but does not converge

Progression of `||R||` on the 1.12x cold direct seed (k=10, sigmoid), target 1e-8:

| solver variant | outcome |
|---|---|
| first cut (`relTrunc=1e-10`, weak Armijo) | **frozen** at 1.96 (alpha~1e-7 microsteps; two bugs: truncation too weak to bite; Armijo accepts microsteps) |
| `relTrunc=1e-6` (weak truncation) | crawl to **~1.8 stall** (kept subspace still cond~1e6, full step overshoots) |
| `relTrunc=1e-2` (aggressive, fixed) | descends 1.96 -> **0.52**, lam=0 (pure truncated GN), then floors (residual trapped in the ~28 truncated near-null directions) |
| **adaptive truncation continuation** (1e-2 -> 1e-11) | descends 1.96 -> **~0.43**, best of all; but crawls again once fully relaxed (the ill-conditioned lambda_r0 direction re-enters) |

**Read:** the scaled truncated-SVD step genuinely fixes the *step strategy* (the
solver now descends where lsqnonlin crawled and the first cut froze), but on the
**cold direct seed** it floors ~0.43 -- it does **not** converge 1.12x. The wall
is the tiny convergence basin of 40-rev shooting seeded cold, exactly what the
campaign and Zhang 2015 say needs *continuation*, not a better local step. This
matches the plan's own odds ("terminal-cluster/many-switch cases will need
Rung 2").

## Rung 1b — clean-band check: does IFS converge where PSR works? (No.)

The Rung-1 solver was first tested only on 1.12x, which sits at the band edge and
carries the terminal switch-cluster the original `RESULTS.md` blamed. To separate
"terminal cluster" from "cold-seed basin," `ifs_solve2` (adaptive truncation) was
run on the clean, well-separated cases where the direct/PSR solution works
(`test_clean_band.m`, `legacy_ms_f11{20,40}` / `f1250` / `f1850`; ode 1e-10):

| case | k | seed ||R|| | final ||R|| | outcome |
|---|---|---|---|---|
| 1.12x | 10 | 1.96 | ~0.43 | stall (flag -2) |
| **1.14x** (clean, well-sep) | 24 | **1.52** | **0.47** | **stall** |
| **1.25x** (clean, well-sep; GPT benchmark) | 47 | **1.31** | **~0.30** (crawl/stall) | **stall** |
| 1.85x | 22 | **1.18e5** | 1.2e4 | stall — *seed unusable* |

**Two findings:**
1. **The wall is the cold-seed basin, not the terminal cluster.** The clean,
   well-separated 1.14x (k=24) and 1.25x (k=47) cases -- good seeds (1.5, 1.3),
   no cluster -- floor at the SAME ~0.3-0.47 as 1.12x. `ifs_solve2` descends
   ~3-5x from the direct dual-map seed and stalls at *every* t_f tested. More
   switches (shorter arcs, less per-arc sensitivity) helps only marginally
   (1.25x floors a bit lower). So IFS does **not** converge where PSR works,
   for the same reason it fails at 1.12x: the 40-rev cold-seed shooting basin is
   small across the whole band. (Consistent with the campaign thesis: shooting
   over 40 revs needs a truly warm seed, not the direct solution.)
2. **Dual->costate SEED quality is erratic across t_f.** 1.12x/1.14x/1.25x seed
   at ~1.3-2.0 (usable), but 1.85x seeds at **1.18e5** -- a beta-scaling failure
   in the mode-'d' dual map (`sms_seed_duals`). So even where PSR gives a clean
   direct solution, the IFS seed derived from it can be unusable; at 1.85x the
   bottleneck is the seed map, not the solver.

**Consequence for the two levers:** both remaining paths are seed-quality paths.
(a) A **warm/continuation seed** from a converged neighbor -- blocked by the
min-time fold+gauge (Rung 2 below). (b) A **better costate seed** than the dual
map -- the erratic seed quality (1e5 at 1.85x) says the dual->costate recovery is
itself a weak link worth fixing independent of the solver.

## Rung 2 — min-time anchor converges (milestone); t_f-continuation hits the fold

**Anchor (factor=1.00, k=0 all-burn):** the FIRST end-to-end IFS rendezvous
convergence (previously only unit tests). Seed residual 2.4e-5 ->
`ifs_solve2` 2.3e-7; switching function **S<0 everywhere** (min -45.8, max
-1.48) -> all-burn confirmed, k=0 correct, margin to march up. `||R||` floors at
~2e-7 because the min-time point is genuinely near-singular (fsolve on the
min-time reference also reports "locally singular").

**Marching up in t_f FAILS by naive stepping — the min-time point is a fold.**
Fresh k=0 solves from the anchor at small t_f offsets:

| factor | seedRes | final `||R||` | converged | maxS |
|---|---|---|---|---|
| 1.001 | 0.303 | 0.021 | no | -1.53 |
| 1.002 | 0.551 | 0.031 | no | -1.58 |
| 1.005 | 0.950 | 0.086 | no | -1.77 |

Even a **0.1%** t_f step throws the seed residual to 0.30 (terminal is
hyper-sensitive to t_f near min-time), and the near-singular Jacobian floors the
solve at 0.02-0.09 without converging; `maxS` moves *away* from 0, so these are
not on the true min-fuel branch. This is the classic **vertical-tangent
degeneracy at min-time** (d(terminal)/d(lambda) singular there) — naive t_f
parameter-stepping cannot cross it.

## Where it stands / next lever

- **Solver infrastructure (Rung 1): built, validated, working** — descends,
  truncation-continuation active, k=0 supported, direct-tau option. It is the
  solver a working continuation will drive.
- **Anchor (Rung 2): working** — IFS converges the min-time all-burn rendezvous.
- **Blocker:** crossing the min-time fold needs **pseudo-arclength continuation**
  (parameterize by arclength along the branch; treat t_f as an unknown with an
  arclength constraint), the standard turning-point tool GPT-5.6-sol flagged.
  This is the recommended next build. Alternatives: free-`tau_f` reparameterization
  (sigma in [0,1] + unknown Sundman length), or anchoring the continuation above
  the fold via a different (non-min-time) converged k=0 point if one can be found.

IFS remains OPEN for the full 1.12x gate. But it now has (a) a solver that
descends, and (b) its first real rendezvous convergence (the anchor) — the two
pieces continuation needs.

## Rung 2b — pseudo-arclength continuation across the fold (`ifs_tf_arclength.m`)

Built the standard turning-point tool: state `x=[lambda0(8); factor]`, an
arclength constraint, predictor (tangent) + Newton corrector on the 9x9
`[R; arclength]` system, complex-step `dR/dx`, scaled metric, `tauf` fixed per
corrector and updated between steps.

Two bugs found and fixed while bringing it up:
1. **Corrector acceptance** — the min-time anchor R-floors at ~3e-7 (the
   min-time reference costates are themselves only ~1e-6 accurate; `svd`/`fsolve`
   both flag near-singularity), so demanding `||G||<1e-8` made every corrector
   fail. Fixed: accept when `||R||` reaches its achievable floor (`corrAcceptR`)
   with the arclength constraint satisfied.
2. **Tangent from `svd(...,'econ')`** — for the 8x9 extended Jacobian, `'econ'`
   returns V as 9x8 and DROPS the 9th (null) column, so `V(:,end)` picked the
   smallest non-null singular direction (nearly pure lambda0, factor-component
   4e-4) → factor frozen. Fixed to full `svd` and `V(:,9)` (the true tangent).

**Spectrum at the anchor** (scaled `[dR/dlambda0 | dR/dfactor]`, 8x9): singular
values 1.7e6 … 8.3e-4 — **full rank 8** (not singular), cond ~2e9. So the null
space is a clean 1-D tangent; the fold is well-posed for pseudo-arclength.

**Status (march running):** with both fixes the corrector converges (rn ~1e-7)
and the continuation moves ALONG the branch — **max S rises monotonically**
(-1.48 -> -1.12 over 8 steps), i.e. toward the first switch birth. BUT the
min-time fold is nearly **vertical**: `factor` stays 1.0000 to 4 decimals while
lambda0 changes a lot (a degenerate ~1-parameter family at exactly min-time).
`factor` should advance only quadratically in arclength, so it needs many steps
(ds auto-grows). A longer march (`ifs_tf_arclength(1.00,1.15,0.10,...)`,
`arclen_march.log`) is running to reach either the birth (max S -> 0) or the
point where `factor` starts advancing. Open question being resolved by that run:
whether max S reaches 0 essentially AT min-time (first switch born right at the
band's lower edge) or after `factor` climbs into the 1.01-1.11x band.

**March outcome (2026-07-11) — the min-time anchor is too degenerate; the
"birth" is a gauge artifact.** The 44-step march (`arclen_march.log`) never
advanced `factor` off 1.0000 (to 4+ decimals) and `max S` **wandered
non-monotonically**: -1.48 -> -0.88 -> -0.95 -> ... -> **-3.23** -> -2.80 ->
... -> -0.50 -> **+0.037**, at which point it declared a switch birth. A physical
t_f-branch would have `max S` rising *monotonically* as t_f grows. Instead the
continuation crawled around the **degenerate 1-parameter costate family at fixed
t_f = min-time** (the near-null gauge, smallest scaled singular value 8.3e-4),
and the "birth at factor=1.0000" is merely where `max S` happened to cross 0
while wandering that manifold — NOT a physical transition-band result.

**Diagnosis:** at min-time the all-burn branch is effectively *vertical* in
(t_f, lambda0) AND carries a near-null costate gauge, so the arclength tangent's
factor-component is ~0 and a loose corrector (accept `||R||`<1e-5, forced by the
~1e-6 min-time seed accuracy) lets lambda0 drift along the gauge at fixed t_f.
Pseudo-arclength is mechanically correct (correctors converge, extended Jacobian
full-rank) but cannot advance t_f from this pathological anchor.

**What it would take to cross for real (next levers, not built):**
1. **Regularize the gauge** — add a phase/pinning condition on the near-null
   lambda0 direction (or truncate it), so the corrector lands on unique branch
   points and the tangent's factor-component is recoverable.
2. **A non-degenerate anchor** — anchor the continuation slightly ABOVE min-time
   where the branch is not vertical; needs one converged non-min-time k=0 point
   (chicken-and-egg without a gauge-regularized solve).
3. **A tight corrector** — the ~1e-6 min-time seed limits corrector accuracy;
   the min-time reference costates would need re-converging tighter first.

**Net for Rung 2b:** the pseudo-arclength machinery is built and validated, and
it precisely *characterized* why the min-time anchor resists t_f-continuation
(vertical branch + costate gauge) — the same degeneracy class the direct
campaign hit in the 1.01-1.11x band, now seen from the indirect side. IFS
remains OPEN; the clean next step is gauge regularization of the anchor.

## External consult — GPT-5.6-sol (2026-07-12)

Sent a self-contained review of the above to GPT-5.6-sol (via crush/OpenRouter).
Prompt: `CONSULT_GPT56_prompt.md`. Full reply: `CONSULT_GPT56_response.md`.
Summary of its verdict + our cross-check against what we actually ran:

**It corrects our "gauge" diagnosis (and this matters).** The near-null
direction at the min-time anchor is *not* a scale gauge of the min-fuel
problem — it is the **separation between two distinct BVP manifolds**. At
t_f=t_f,min the min-time (H(τ_f)=0, free costates) and min-fuel (λ_m(τ_f)=0, H
not pinned) problems share the all-burn *trajectory* but have different
transversality, so their costate manifolds are near-tangent-but-distinct.
Consequence: our proposed **lever #1 (pin ‖λ0‖=1 / phase condition) would fix
the wrong null direction** — the min-fuel problem has no scale gauge (the cost
breaks it), so gauge pinning cannot recover the continuation. This retires
"gauge regularization of the anchor" as the next step.

**It says drop the min-time anchor entirely.** Confirms the inconsistency we
flagged: the *direct* side homotopes from min-energy, but we anchored the
*indirect* side at min-time. Min-energy anchor is strictly better (no gauge, S
bounded away from 0, dλ0/dt_f nondegenerate). If a t_f gap remains after
re-seeding, use a Bertrand-Épénoy cost homotopy `J_s=(1-s)½∫(T/Tmax)²+s∫(T/Tmax)`,
s:0→1, with **switch-structure locking** (promote each arc to burn/coast when
`T*` saturates >0.95 / <0.05, freeze structure at `s_lock`<1, extract costates
there). This is a *disciplined* version of the Rung-0 mint that failed — the
difference is it never runs a free ε=0 solve, so it should not slide to the
many-switch global basin.

**Its headline recommendation — the actual next build: a backward adjoint
sweep to replace the KKT-dual→costate seed.** Instead of *sampling* the direct
NLP's discrete duals (`sms_seed_duals` mode-'d'), *integrate the adjoint ODE
backward* along the frozen direct state trajectory (dense-output interp, tight
rtol), initialized at τ_f from terminal transversality (λ_m(τ_f)=0 exact; λ_r,
λ_v from the terminal rendezvous multipliers). Seeds then carry
integration-tolerance accuracy, not O(h) mesh accuracy. Its Lyapunov argument
(0.01% error in λ0 → O(1) endpoint residual over 40 revs) is exactly the
mechanism behind our observed seed pathology (‖R‖ floor 0.43 at 1.12×; blowup to
1.2e5 at 1.85×). **This directly attacks the wall Rung 1b identified (cold-seed
basin, not terminal cluster), and it is categorically different from anything we
tried — we only ever sampled duals.**

**Cross-check caveats we hold (things sol glosses):**
1. Backward-adjoint conditioning is *not* free — the backward adjoint of an
   unstable 40-rev symplectic flow can still amplify error. The real edge is the
   *initialization point* (τ_f: λ_m=0 exact, r/v are boundary duals), not a
   guarantee of accuracy through the sweep. Test empirically; don't assume.
2. sol's `λ_t(τ_f)=0` is wrong for our formulation — τ_f is fixed and t_f is a
   terminal constraint, so λ_t is a constant set by that constraint's multiplier
   (our anchor already used `λ_t0=-H_t(0)`). Init λ_t from the direct dual / H,
   not to zero.
3. The s-homotopy switch-lock is the same cliff as the Rung-0 mint (ε=0 slid
   6.35→5.4 into the global basin). The locking discipline is plausibly the fix,
   but eyes open.

**Revised next lever (supersedes "gauge regularization"):**
1. Build the **backward adjoint sweep** seed → re-seed IFS on the existing
   direct solutions (1.12×/1.14×/1.25×). Falsifiable: if the wall is seed
   quality, the ~0.43 floor should become real convergence; if not, the wall is
   more than seeding.
2. Only if a t_f gap remains: **min-energy anchor + s-homotopy** with switch
   locking. Drop min-time and pseudo-arclength.

## Rung A (plan 2) — adjoint sweep AND adjoint smoother: FALSIFIED (2026-07-12)

Built `ifs_seed_adjoint.m` (+ `test_ifs_seed_adjoint.m`) per `PLAN_OF_ATTACK_2.md`
Rung A and ran the gates on 1.12×. Three variants, all fail GA2. This closes
the "integration-accuracy seed" idea conclusively — and the *reason* is now
measured, not conjectured.

**A-1, pure backward sweep (`method='sweep'`).** Terminal init verified good
(`lamAll` terminal duals; q = 0.98 at the switch nearest τ_f, i.e. S≈0 to 2%)
— but the adjoint flow amplifies **1.2e12** over the 40-rev backward span, so
the 1–2% terminal-dual error rides to a seed residual of **2.98e10** (vs
dual-map 1.96). sol's Q4 accuracy prediction ("~1e-8–1e-10") assumed the sweep
error is set by integration tolerance; it is actually set by
(terminal-init error) × (backward amplification). Our held caveat, confirmed.
Corollary (struck from the plan): the A2 "linear collocation" insurance solves
the SAME IVP — same solution, same blowup; it guarded against marching error,
but the disease is genuine flow amplification.

**A-2, adjoint smoother, free data scale (`fitScale=true`).** Fit the terminal
costate to the whole dual history (truncated-SVD GN, complex-step sweep
sensitivities, 256 data nodes + k switch rows). Fit descends 7.0e11 → 9.2e3
and floors; **best-fit trajectory stays ~48 scaled units from the dual data**
(blend trusts 0% of nodes). Implementation artifact found: the free scale `s`
drifted 5.6× (10 switch rows too weak against 2048 data rows), corrupting the
dual fallback → seed 17.0. Fixed by pinning.

**A-3, adjoint smoother, scale pinned (`fitScale=false`, now default) — the
clean control.** Fit floors at 9.7e3; misfit median **54.6** scaled units;
only **2.2%** of nodes (the tail near τ_f) trust the sweep within 5%. The
blended seed has an essentially perfect *switching* structure — S-sign
agreement 100%, 10/10 crossings, **q at switches = 1.001 ± 0.008** — but the
MS residual is dominated by arc-continuity rows at data accuracy:
**seedRes 2.13 vs dual-map 1.96** (0.92×, no gain).

**The measured mechanism.** The dual-map costate history (O(h) mesh accuracy,
~1%) is **dynamically inconsistent at the 1e12-amplification level**: no exact
adjoint trajectory passes within ~50 scaled units of it over the full span.
Single-trajectory costate reconstruction from mesh-accuracy data is therefore
impossible at 40 revs — the information is not in the data. (This is also
*why multiple shooting exists*: the dual seed only ever propagates short arcs,
which is why its residual is 1.96 and not 1e10.)

**What survives of the seed-quality axis:** only reducing the node-wise dual
error itself (finer/PSR-refined direct mesh, tighter direct tolerance — plan
Rung D substrate), with bounded expected gains: per-arc residual ≈ (node
error) × (per-arc amplification ~e^2.7≈15), consistent with the observed 1.96.

**GA3 addendum (run, 2026-07-12):** `ifs_solve2` from the A-3 blended seed,
1.12×, 60 iters (ode 1e-10): **2.13 → 0.224** — roughly **half** the dual-seed
floor at a comparable budget (1.14× clean-band hit 0.47 at 60 iters; 1.12×
dual-seed runs floored ~0.43) — but with the SAME crawl signature (α=0.03,
~0.3%/iter steady, equilibrated cond ~1e16, truncation relaxing on plateau).
At that rate 1e-8 needs thousands of iterations. **Read:** the better switch
structure (q=1.001) buys a real but bounded improvement; the continuity rows
and the conditioning crawl still dominate. Falsification of "seed quality
alone closes the gate" stands; Rung B is next.

**GA3-long (300-iter due diligence, `ga3_long_results.mat`):** 2.128 →
**0.0398 at iter 272, then HARD STALL** (flag -2: no acceptable step at any
damping, λ→6.9e11, relTrunc fully relaxed to 1e-11). ~10× below the dual-seed
floor and comparable to the best number ever reached from any seed (the
original 400-iter lsqnonlin run's 0.0232) — but the same wall, now hit
decisively rather than asymptotically. The 1.12× stall level ~2–4e-2 is
seed-independent; conclusion unchanged.
