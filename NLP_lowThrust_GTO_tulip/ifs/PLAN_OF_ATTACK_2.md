# IFS Plan of Attack 2 (2026-07-12)

Supersedes the rung sequence of `PLAN_OF_ATTACK.md` (2026-07-11), whose Rungs
0/1/2/2b were executed and recorded in `RESULTS_RUNG01_RUNG2.md`. This plan
folds in what that campaign *proved* and what the GPT-5.6-sol consult
(`CONSULT_GPT56_prompt.md` / `CONSULT_GPT56_response.md`) corrected. The two
reference-paper digests (Zhang 2015, Leomanni 2021) in the old plan §2 remain
valid and are not repeated here.

## 1. What we now know (the ground truth this plan is built on)

1. **The solver is no longer the bottleneck.** `ifs_solve2` (equilibrated
   truncated-SVD GN + adaptive truncation) descends everywhere the old LM
   crawled or froze. Machinery (EOM/residual/Jacobian/seed) is unit-tested and
   validated.
2. **The wall is SEED QUALITY, uniformly across t_f.** The cold KKT-dual seed
   floors at the same ~0.3–0.47 at 1.12× (k=10), 1.14× (k=24), 1.25× (k=47) —
   with and without a terminal switch cluster. At 1.85× the dual map itself
   fails (seed 1.2e5). Mechanism (sol, and consistent with everything we
   measured): collocation KKT duals are O(h)-accurate costates, and 40-rev
   shooting amplifies a ~0.1–1% λ0 error to O(1) endpoint residuals.
3. **The min-time anchor is a dead end for continuation.** It converges
   (2.3e-7 — IFS's first rendezvous convergence) but t_f-continuation off it
   fails: the near-null direction there is NOT a scale gauge to be pinned; it
   is the **separation between two distinct BVP manifolds** (min-time
   H(τf)=0-free-costates vs min-fuel λm(τf)=0-H-unpinned, tangent at t_f,min).
   Gauge pinning (‖λ0‖=1) fixes the wrong direction. Drop this anchor.
4. **Hard-throttle (no ε) is the right formulation** given a warm seed; ε is a
   cold-start crutch whose 1/ε layers caused the ms_band crawl. Keep IFS as-is.
5. **New disk fact (verified 2026-07-12):** the direct `.mat`s carry what a
   terminal-transversality costate init needs:
   - `legacy_ms_f1120/f1140`: full `out.lamAll` present; layout confirmed
     *exactly* (`len = 8N + nN + 24·nN + 8 + 7 = 132040`), so the **terminal
     rendezvous multipliers** (last 7 entries: rv 6 + t_f 1) are directly
     readable.
   - `legacy_ms_f1250/f1850`: no `lamAll`, but the terminal duals agree with
     `−lamDef(:,end)` to ~1–3% (checked on f1120/f1140), so `−lamDef(:,end)`
     is a validated fallback.
   - λ_t convention pinned by data: t_f-dual = `−lamDef(8,end)` to 4 digits,
     and λ_t is a constant of the motion (τ-autonomous H) — init once, hold.

## 2. Central hypothesis (falsifiable)

> The IFS basin is reachable from a costate seed whose error is set by
> *integration tolerance*, not by *mesh spacing*. Build that seed by
> integrating the adjoint equations backward along the frozen direct state
> trajectory, initialized at τ_f from terminal transversality + terminal KKT
> multipliers. If this seed still floors at ~0.3, the wall is not seed quality
> and the plan pivots to Rung B.

## 3. The rungs

### Rung A — backward adjoint sweep seed (`ifs_seed_adjoint.m`) — BUILD FIRST

Replace the dual-*sampling* seed (`sms_seed_duals` mode-'d') with an adjoint
*integration* seed:

1. **Terminal init at τ_f** (this is the whole trick — the best-known point):
   - λ_m(τ_f) = 0 **exact** (min-fuel transversality),
   - λ_rv(τ_f) = β · (terminal rendezvous duals) — from `lamAll(end-6:end-1)`
     where present, else `−lamDef(1:6,end)`,
   - λ_t = β · (t_f dual) (constant; equals `−lamDef(8,end)`),
   - β: the min-fuel costates have a definite scale (the "1" in
     S = 1 − ‖λv‖c/m − λm breaks scale invariance). For this transcription the
     mode-'d' analysis says β≈1; `beta_from_duals` gives the empirical fit.
     Because the sweep is *affine* in (β·λf), a small 1-D β search (minimize
     seed residual / S-sign violations over the sweep) is cheap. This also
     retires the erratic-β failure (1.85×) — β becomes a fitted scalar against
     the whole sweep, not a switching-law fit on noisy interval duals.
2. **Sweep**: integrate the 8 costate ODEs **backward** τ_f→0, with the state
   frozen to the direct solution (dense interpolation of `out.X` on the σ
   grid, 4001 nodes) and the throttle from the direct arc structure
   (hard u∈{0,1} between `diag.tauCr` switch times). Implementation: call
   `ifs_eom` with Y = [x_dir(τ); λ] and keep only dY(9:16) — **exactly
   consistent** with the residual's EOM by construction, zero new derivation.
   Tight tolerances (RelTol 1e-11 / AbsTol 1e-13).
3. **Node extraction**: sample the sweep at the IFS shooting nodes → replaces
   only the costate half of the current seed; states/switch times unchanged.
4. **A2 variant (stability insurance, build only if A1 blows up):** the adjoint
   system is *linear* in λ, so instead of backward *marching* (which can
   amplify over 40 revs — sol glosses this; we don't), solve it as one
   trapezoidal **linear collocation system** on the direct σ-mesh with the same
   terminal values — a banded linear solve, unconditionally stable, no
   marching error accumulation.

**Gates (in order, each falsifiable):**
- **GA1** — sweep sanity: λ(τ) finite over the whole sweep; S(τ) sign pattern
  reproduces the direct burn/coast arcs (report % agreement; dual-map seeds
  should be beaten decisively).
- **GA2** — seed residual: `ifs_residual` norm at the sweep seed vs the dual
  map's 1.96/1.52/1.31. sol predicts 1e-8–1e-10; anything ≤1e-2 already
  changes the game. If it *degrades from τ_f backward* → backward instability
  → build A2.
- **GA3** — the real test: `ifs_solve2` from the sweep seed on
  **1.12×/1.14×/1.25×** to tolR 1e-8, then `ifs_certify` green. Success on any
  one of these = IFS's first full min-fuel convergence.

Effort: small (one new seed function + a test + reruns). Highest
information-per-hour in the whole plan.

> **RUNG A OUTCOME (2026-07-12): FALSIFIED — all variants.** Built
> (`ifs_seed_adjoint.m`, methods 'sweep' and 'smooth') and gated on 1.12×:
> (1) pure sweep: terminal init good (q=0.98 at the last switch) but the
> adjoint flow amplifies **1.2e12** backward over 40 revs → seed 2.98e10;
> (2) A2-as-specified is struck — linear collocation solves the same IVP,
> same blowup (it guarded against marching error; the disease is flow
> amplification); (3) adjoint SMOOTHER (terminal fit to the whole dual
> history, truncated-SVD GN — built as the rescue): best-fit trajectory stays
> ~50 scaled units from the dual data; blended seed 2.13 vs dual 1.96. The
> measured mechanism: mesh-accuracy (~1%) dual costates are **dynamically
> inconsistent at the 1e12 level** — no single 40-rev adjoint trajectory fits
> them; the information is not in the data. The central hypothesis (§2) is
> falsified in its strong form. Seed-quality gains are bounded to reducing
> node-wise dual error (Rung D substrate). → proceed to **Rung B**.
> Full numbers: `RESULTS_RUNG01_RUNG2.md` §"Rung A".

### Rung B — min-energy anchor + s-homotopy with early structure lock (only if A floors)

> **Detailed execution plan (2026-07-12): `PLAN_RUNG_B.md`** — grounded in the
> ms_band Gate-D history (why this is not a repeat), the entropy-vs-objective
> smoothing distinction, the smooth-β trap, verified disk assets (the `_en`
> files are NOT energy solutions), and increments B0–B4 with gates and a
> stopping rule. Start there.

sol's fallback, and a disciplined version of the ε-march that failed in
ms_band. Composition of two *existing* codebases:

1. Converge the **smoothed indirect MS** (`ms_band/sms_*`) at s=0 (pure energy,
   ε=1): smooth throttle, no layers, no gauge — should be the easy indirect
   solve, seeded from the direct energy solution's duals.
2. March s upward (Bertrand–Épénoy J_s = (1−s)·½∫u² + s·∫u dt) **only until
   arcs saturate** — monitor T*(τ) per arc; when every arc is uniformly
   ∈[0,0.05]∪[0.95,1] for two consecutive steps, **freeze the switch
   structure at s_lock < 1** and stop the march *before* the 1/(1−s) layers
   form (this early stop is exactly what ms_band's ε-march didn't do — it
   marched to the crawl).
3. Hand off states + costates + switch times at s_lock to hard-throttle IFS.

Known cliff (eyes open): the Rung-0 mint slid 6.35→5.4 into the many-switch
global basin during a free ε=0 re-solve. Rung B never runs a free ε→0 solve —
the lock is the fix — but arc-count changes during the s-march must be watched
(promote/demote arcs only by the saturation rule, never re-derive structure).

### Rung C — t_f-continuation *between converged IFS points* (the band prize)

Once ANY t_f converges (Rung A or B), the 1.01–1.11× band gets attacked by
marching t_f with **warm IFS restarts from the converged neighbor** — not from
min-time (manifold-tangency, retired) and not from cold seeds. Switch
birth/death across the band handled structurally: insert a switch where S(τ)
develops a tangency (S=0, dS/dτ=0), delete when an arc length → 0. This is
Zhang 2015's medicine, now with a legitimate anchor.

### Rung D — structural fallbacks (only if A and B both stall)

- Interior non-switch shooting nodes (shrink per-arc sensitivity; we know from
  1.25× k=47 that arc count alone is marginal *with bad seeds* — retest with
  good seeds before judging).
- PSR-refined substrate for the sweep (tighter direct defects → cleaner
  x_dir(τ) forcing, directly reduces the sweep's error floor).
- Tighten the direct collocation tolerance before sweeping (sol's Q5 note: if
  a sweep seed stalls at ~1e-4, blame the direct defect, not ε).

## 4. Decision tree

```
Rung A (adjoint sweep)
├── GA2 seed ≤ ~1e-2 and GA3 converges → DONE at that t_f → Rung C (band)
├── GA2 good, GA3 floors lower but > tol → Rung D substrate/tolerance, retry
├── GA2 bad, error grows backward       → A2 linear-collocation sweep, retry
└── A2 also floors ~0.3                 → hypothesis falsified: wall is not
                                          seed quality → Rung B
Rung B (s-homotopy + lock)
├── converges → hand off to IFS → Rung C
└── slides to many-switch basin at lock → structure-freeze discipline failed:
    revisit lock criterion; if robust failure → problem needs new coordinates
    (campaign doc's "regularized coords" thesis) — stop and reassess.
```

## 5. Odds and stopping rule (honest)

- Rung A is cheap and directly falsifiable; even its failure is decisive
  information (it isolates the wall to something other than seed accuracy).
  Odds of full convergence at ≥1 clean t_f: genuinely unknown — the Lyapunov
  argument cuts both ways (backward sweep can amplify too; that is what A2 is
  for). Call it even.
- Rung B reuses two working codebases; its risk is concentrated in one place
  (the lock). 
- If A, A2, AND B all fail, the honest conclusion is the campaign doc's
  standing thesis: sharp many-switch indirect solves at 40 revs need
  **regularized coordinates** (orbital elements / Leomanni-style), which is a
  new build, not a patch — stop and decide deliberately.

## 6. Pointers

- Consult that reshaped this plan: `CONSULT_GPT56_prompt.md`,
  `CONSULT_GPT56_response.md` (cross-checked summary in
  `RESULTS_RUNG01_RUNG2.md` §"External consult").
- Executed prior plan + results: `PLAN_OF_ATTACK.md`, `RESULTS.md`,
  `RESULTS_RUNG01_RUNG2.md`.
- Paper digests (Zhang 2015 = continuation existence proof; Leomanni 2021 =
  regularized-coordinate alternative): `PLAN_OF_ATTACK.md` §2, PDFs in
  `../min_fuel_papers/`.
- Dual map being replaced: `../ms_band/sms_seed_duals.m` (mode 'd'),
  β fit: `../ms_band/beta_from_duals.m`.
- Terminal-dual disk layout (verified): defects 8N | unit nN | bounds 24·nN |
  init BC 8 | **terminal BC 7 (rv 6, t_f 1)** — tail of `out.lamAll`.
