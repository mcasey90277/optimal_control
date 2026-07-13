# ZTL campaign results (Zhang-style thrust ladder)

Honest running record. Plan: `../ifs/PLAN_PRONG_Z.md` (execution),
`../ifs/PLAN_OF_ATTACK_3.md` (strategy + Zhang audit).

## Lessons learned (transferable) — the distilled findings

These transcend this problem and apply to any multi-revolution indirect
optimal-control / shooting work. Detailed evidence is in the dated sections
below.

1. **Differencing-based Jacobians are structurally inadequate for multi-rev
   CR3BP shooting.** BOTH complex-step (through an adaptive integrator) and
   central finite-difference fail: the derivative scales span ~11 orders, CS
   is corrupted at O(1) by adaptive-step coupling, FD takes secants across a
   curved valley. Only the VARIATIONAL STM (exact derivative of the continuous
   system, integrated alongside the state) is usable. This was measured, not
   assumed (P0h/P0i), and it confirms Zhang 2015's ingredient (a) as
   load-bearing. It likely explains a large part of the whole indirect
   campaign's historical crawl.

2. **The indirect shooting Jacobian is ARTIFICIALLY ill-conditioned by the
   disparate state-vs-costate scales — equilibrate it.** Raw cond ~3e11;
   two-sided (Ruiz) equilibration -> 5e6, column-only -> 8e8. This is a large,
   cheap, likely campaign-wide-missed lever. CAVEAT: row-scaling the RESIDUAL
   changes the objective and produces ascent steps far from the solution
   (instant stall); COLUMN scaling is a pure reparametrization (objective-safe)
   -- use it, or use two-sided only to compute an accurate Newton DIRECTION
   with a true-objective globalization.

3. **Node placement does NOT fix conditioning; equilibration does.** Perigee
   amplification is intrinsic and un-splittable -- cond(J_MS) plateaus ~1e11
   for uniform M=26..208, and amplification-equidistributing placement does
   not help (int|A|dt placement makes it WORSE). Placement is a LINEARITY lever
   (shorter stiff arcs), a separate axis from conditioning.

4. **Multiple shooting is NECESSARY, not merely nicer, for multi-rev
   indirect.** The decisive test: at a converged-ish MS iterate (||R||~1e-5,
   continuity gaps ~5e-6), the SINGLE-shooting residual at the SAME lam0 is
   2.56 -- enormous. The tiny continuity gaps are LOAD-BEARING: they absorb the
   13-rev amplification that single shooting cannot represent. This is why SS
   floors dead and MS descends 80x past it. Corollary: "re-chop to zero the
   gaps" is a dead end (throws away the amplification the gaps encode).

5. **Cold-seed landscape asymmetry: "wide basin at high thrust" really means
   ARRIVES-WARM easy.** High thrust -> a wrong seed crashes into Earth/escapes
   within a rev (explosive residual 1e10-1e21); low thrust -> tame landscape
   but lethal 40-rev amplification. Cold is hard everywhere; only the failure
   mode changes. Zhang's thrust ladder works because every solve after the
   first arrives warm.

6. **Min-time is a bad continuation substrate here** (near-singularity at
   factor 1.0; machinery-independent). It is only ever a one-shot seed source.

7. **A proper trust region converts STALL into CONVERGENCE, but nonlinearity
   caps the rate.** The SVD-LM trust region (More) descends monotonically where
   diagonal-LM and backtracking-Newton crawl/stall -- but the stiff
   perigee-continuity nonlinearity keeps the trust radius small, so it
   asymptotes slowly. The remaining lever is LINEARITY (perigee-concentrated
   nodes), not conditioning or step strategy.

## P0 — preflight (2026-07-12)

### P0a — graze margin at the PSR bang-bang solutions: PASS, favorable

From the dual-mapped switching function S(sigma) of every
`PSR_data/psr_data_tf*_minEps0.mat` (slope = LSQ line through the ±3 nodes
bracketing each crossing; `p0a_graze_margin.m`):

| factor | crossings | min \|dS/dtau\| | med \|dS/dtau\| | min/med | min \|dS/dt\| |
|---|---|---|---|---|---|
| 1.12x | 10 | 1.130e-1 | 4.359e-1 | 0.26 | 4.35e-1 |
| 1.15x | 26 | 5.348e-2 | 3.511e-1 | 0.15 | 7.84e-2 |

The shallowest crossing at the Z4 target (1.15x) carries 15% of the typical
slope — no graze; the saltation 1/Sdot is well-bounded. Zhang's known
weakness is NOT active at the target (expect it to activate somewhere in the
Z5 band march below 1.12x — guarded there).

### P0b — reproduce the 2025 min-time ladder: FAIL, and the old table is an artifact

`p0b_mintime_ladder.m` re-ran the old Phase-1 up-march (solve_tfmin_indirect
= complex-step + trust-region-dogleg), warm from the validated 25 mN
min-time point:

| f (thrust) | tf_min | ||R|| | flag |
|---|---|---|---|
| 1.0 (25 mN) | 6.290694 | 1.25e-6 | -3 (near-singular — known, benign) |
| 1.2 | 6.107257 | 3.16e0 | 0 (iter limit) |
| 1.5 | 6.084190 | 2.62e0 | -3 |
| 2.1 | 6.084190 | 3.96e3 | -2 |
| 3.2 | 6.079370 | 1.02e1 | -3 |
| 5.0 | 6.065106 | 1.32e3 | -2 |
| 8.0 | 6.065106 | 7.71e14 | -2 |

**Finding: the old 2025 tfMin table was never a converged backbone.** Its
duplicate entries (6.0842 at f=1.5 AND 2.1; 6.0651 at f=5 AND 8) are stalled
solves returning an unmoved guess — the same duplicates reproduce here with
residuals O(1)–O(1e14) and ode113 step-underflow crashes (trajectories
hitting the Earth singularity mid-iteration). PLAN_PRONG_Z §1's "Phase 1
validated asset" claim is hereby CORRECTED: only the 25 mN rung was ever
converged. This strengthens the audit — the old machinery (CS-through-ode113
+ dogleg, no events, no STM) fails even for MIN-TIME above ~1.2x thrust.

### P0c — energy probe at 200 mN: INVALIDATED (garbage seed)

`p0c_energy_top.m` was seeded from the P0b f=8 "min-time costates," which
P0b showed are garbage (res 7.7e14). Best of the beta grid reached
||R||=1.67 with a plausible interior throttle (min 0.074 / sat-hi 2.2%) —
suggestive of a nearby smooth solution, but not evidence. Superseded by P0d.

### The fixed-tf simplification (kills the Z2 ladder requirement)

tf_min is **provably nonincreasing in Tmax** (a larger control set cannot
lengthen a min-time), and tf_min(25 mN) = 6.2907 IS validated. Therefore
holding **tf = 1.15 x 6.2907 = 7.2344 ND fixed across the entire thrust
ladder** guarantees >= 15% margin over tf_min at every rung, and the
effective factor drifts only ~1.15 -> ~1.19 across 25–200 mN (tf_min is
phasing-limited, near-flat in thrust). Consequence: **Z2 (min-time ladder)
is NOT needed for Z3** — the energy ladder runs at fixed tf with thrust as
the only continuation parameter. Min-time at the top rung is still useful
once, as the Stage-1 seed for the top-rung energy solve.

### P0d — min-time march with pumpkyn's analytic-STM solver: FAIL (3 runs) — min-time seeding route retired

`pumpkyn.cr3bp.tfMin` (Koblick, Oct 2025) has genuinely Zhang-grade
machinery: full 14+196 variational STM, analytic Jacobian, ode45 1e-10/1e-12
WITH throttle-switch event-restart. Three runs of `p0d_top_anchor.m`:

1. **Direct 25→200 mN jump: hopeless.** lsqnonlin ground through an
   Earth-crash region, 2.65M integration warnings, no progress. Killed.
2. **1.25x thrust march via pumpkyn's wrapper: every rung stalls at its
   internal cap** — the wrapper hardcodes MaxFunctionEvaluations =
   MaxIterations = 100 (each rung = exactly ~100 evals then quit,
   ||R|| 0.9→18 growing up the ladder). Useful side-finding: tf_min
   genuinely DROPS with thrust (6.29 → 5.75 → 5.55 → 5.14, unconverged but
   directionally clear) — the old table's plateau was pure artifact. The
   fixed-tf argument survives (needs only monotonicity), but the top-rung
   effective factor is ~1.4, not ~1.19.
3. **Same march via `ztl_mintime_solve.m`** (same residual/J rebuilt from
   the public propagator, real budgets: 1500 evals, TRR, tf bounds):
   STILL fails — rung 1 quits at step-tolerance without polishing the
   legacy seed past 1.1e-3 (formulation offset + the known min-time
   near-singularity), later rungs burn the full budget crawling
   (||R|| ~0.5–30, 900–1800 s/rung, one tf excursion to 8.21).

**Verdict:** the min-time family is a bad continuation substrate here — the
near-singularity that killed IFS Rung 2/2b t_f-continuation ALSO makes the
thrust march ill-conditioned, machinery-independent. Min-time was only ever
the SEED source for the top-rung energy solve, and the fixed-tf argument
already removed every other need for it. **Route retired; the top anchor
gets attacked directly (P0e).** `ztl_mintime_solve.m` is kept (analytic-STM
residual pattern, reusable for Z0/Z1 reference).

### P0e — energy multistart at 200 mN: FAIL 0/20 — the landscape is explosive cold

`p0e_energy_multistart.m`: `solve_energy_indirect` (LM, CS Jacobian) at
200 mN, tf = 7.2344 fixed, 20 seeds (4 throttle-rescales of the validated
25 mN min-time costates, 8 perturbed, 8 cold). **All 20 fail**, almost all
instantly: ||R|| ~ 1e10–1e21 with flag 4, two seeds integrate to Inf/NaN
outright; best attempt 5.4e3 (a near-all-burn point, sat-hi 99.7%).

**Mechanism (the P0 discovery):** at 8x thrust a wrong costate seed does not
produce a mildly wrong trajectory — it crashes into the Earth singularity or
escapes within about a revolution, so the GLOBAL residual landscape is
explosive even where the local basin around the true solution may be wide.
Mirror-image of 25 mN, where garbage seeds still give ||R|| ~ 2 (tame
landscape) but 40-rev amplification kills shooting. **"Wide basin at high
thrust" must be read as: arrives-warm easy — cold anywhere is hard.** Zhang's
ladder works because every solve after the first arrives warm; the open
question is only where the FIRST cold solve is cheapest.

### P0f — cold-convergence thrust sweep (5x, 3x, 2x): no convergence, but the sweet-spot signal is real

Same multistart at 125/75/50 mN (high->low), tf = 7.2344 fixed, 8 seeds per
level, Inf-guard, deliberately tight LM caps (150 it / 900 evals). Best of 8
per level:

| thrust | best ||R|| | note |
|---|---|---|
| 200 mN (P0e) | 5.4e3 | explosive landscape, most 1e14+ |
| 125 mN | 9.36 | first sub-1e3 attempts appear |
| 75 mN | **1.42, flag 0** | rescaled-MT seed; BUDGET-CAPPED WHILE STILL DESCENDING |
| 50 mN | 4.04 (flag 4), 6.05 (flag 0) | mixed stall/capped |

No level converged cold, but best-residuals improve monotonically toward
intermediate thrust — consistent with the sweet-spot picture (~13 revs at
75 mN: landscape no longer explosive, amplification not yet lethal). The
two flag-0 (capped-while-descending) attempts are live leads.

### P0g — warm restarts: still descending, still budget-capped (LM crawl)

`p0g_warm_restart.m` (600 it / 5000 evals each): 75 mN reached **0.379**,
50 mN reached **0.99** — both flag 0 AGAIN (budget exhausted while
descending, ~0.2%/iter). Neither stalled; no explosions. The 75 mN iterate's
throttle is fully interior (0.087–0.980, zero saturation) — a genuinely
smooth energy extremal taking shape. Verdict: LM damping is the rate
limiter on this 7-unknown system, the campaign's familiar crawl — a step-
strategy problem, not a wall (yet).

### P0h — tsvd-GN broke the crawl, then floored at 1.6e-2 @ 75 mN

`p0h_gn_finish.m` (column-equilibrated truncated-SVD GN + Armijo, the
ifs_solve2 step recipe on the 7x7 J): **0.379 -> 0.0162 in 400 iterations**
(23x; LM had needed 1200 iterations for 1.42 -> 0.379), including a genuine
acceleration burst (0.21 -> 0.023 between iters 180-250) — then a hard floor
at ~1.6e-2 with cond(Js) climbing 1e7 -> 5e8 and steps rejected at every
truncation level. Final iterate: fully interior smooth throttle
(0.077-0.812, zero saturation), prop 2.72 kg, dV 4.12 km/s. The floor value
echoes IFS's seed-independent 2-4e-2 floor at 1.12x.

### P0h diagnostic — THE SMOKING GUN: the complex-step Jacobian is corrupted

`p0h_diag_floor.m` at the floor iterate:
- **CS vs central-difference Jacobian relative error = 7.77** — they
  disagree at O(1); cond(J_CS) = 9.5e10 vs cond(J_FD) = 1.4e8. Every legacy
  solve in this campaign steered with a CS-through-adaptive-ode113 Jacobian;
  at this iterate that J is unusable beyond its leading subspace.
- Residual is diffuse (all 7 components ~1e-3-1e-2; lam_m smallest) — NOT a
  single bad constraint, so not a formulation defect.
- Weakest singular direction lives in lam_r (rx 0.82 / ry 0.17) — the same
  weakly-determined initial position-costate the IFS scaled-SVD localization
  found at mid-crawl. Full Newton step on the corrupted J is huge (9.2) and
  fails (1.6e-2 -> 2.3).

**Read: derivative quality through the adaptive integrator is the
campaign-wide rate limiter, now measured — the experimental confirmation of
Zhang ingredient (a) (exact variational STM) as load-bearing, not doctrine.**

### P0i — GN on the central-difference Jacobian: the floor survives — NO differencing scheme works here

`p0i_fd_finish.m` (same GN loop, central-difference J, 14 integrations/it):
1.617e-2 -> 1.562e-2 then HARD STALL at iteration 19 — steps rejected at
every alpha and every truncation, despite the equilibrated J being benignly
conditioned (1.7e4). Iterate banked in `results/p0i_fd_finish.mat` (fully
interior throttle 0.077-0.812, prop 2.72 kg, dV 4.12 km/s).

**The sharpened mechanism.** A square nonsingular J with no descent
direction is impossible for a TRUE Jacobian (J'R = 0 with J nonsingular
forces R = 0) — so the FD J is also not the true derivative. The raw J's
scales explain it: singular values span 1.5e8 -> 1.6e-3 (11 orders). The
residual valley is so steep (|dR/dlam| ~ 1e8) that ||R|| = 1.6e-2 means the
steep components are satisfied to ~1e-10, and so curved that a central
difference at h ~ 1e-6 takes a secant ACROSS the valley rather than a
tangent along it. Complex-step fails differently (adaptive-step coupling +
real() guards, P0h diag: O(1) error). **No differencing scheme — CS or FD —
can produce a usable J for multi-rev CR3BP shooting. Only the variational
STM (the exact derivative of the continuous system, no differencing, no
adaptive-step coupling) can.** This is the campaign's crawl mechanism,
measured end-to-end.

## P0 CONCLUSION (2026-07-12)

Preflight complete. No anchor converged yet, but P0 delivered exactly what a
preflight should — it moved the build order and armed the build:

1. **Graze margin healthy at the Z4 target** (P0a) — saltation safe at 1.15x.
2. **Min-time retired** as seed/continuation substrate (P0b/P0d, three
   machineries) — and the fixed-tf argument removed every need for it.
3. **Anchor thrust moved 200 -> ~75 mN** (P0e/P0f): the cold landscape is
   explosive at high thrust ("arrives-warm easy"); cold multistart signal
   peaks at intermediate thrust.
4. **The campaign-wide crawl mechanism identified and measured** (P0h/P0i):
   differencing-based Jacobians (CS AND FD) are structurally inadequate for
   multi-rev CR3BP shooting (derivative scales span 11 orders; CS corrupted
   at O(1), FD secants across the curved valley). Zhang ingredient (a) is
   confirmed load-bearing by experiment, not doctrine.
5. **A warm 75 mN iterate at ||R|| = 1.56e-2 is banked**
   (`p0i_fd_finish.mat`) as the Z0 acceptance test: the falsifiable
   prediction is that a variational-STM Newton step punches through this
   exact floor in O(10) iterations. If it does, the ladder anchor falls out
   immediately and Z3 begins from it.

**Next: build Z0** (`ztl_eom` ramp family + `ztl_A` + variational STM +
`ztl_flow` event automaton, per PLAN_PRONG_Z §4), gated first on the banked
iterate above, THEN the standard Z0 unit gates. The direct-side dual seed
remains the backup anchor route if the variational J somehow also floors.

## Z0 — variational-STM machinery (built 2026-07-12)

Build spec: `Z0_BUILD.md` (self-contained, written for cold pickup). Files:
`ztl_eom.m` (BE ramp family, regime-explicit, CS-safe), `ztl_A.m` (exact
14x14 field Jacobian by complex step OF THE FIELD), `ztl_flow.m` (3-regime
event automaton + variational STM + saltation at eps=0 + graze guard),
`test_ztl_z0.m` (gates), `z0_accept_75mN.m` (pre-registered acceptance).

**Unit gates: 5/5 PASS.**
- G1 legacy equivalence at eps=1 (ode89 oracle): u0 to 1.1e-16, terminal
  state to 3.1e-8 over the full 13-rev arc (cross-integrator floor
  calibrated at 9.6e-7 by `diag_g1_integrator_floor.m`; like-for-like floor
  ~3e-8 from AbsTol on rescaled costates). Validates EOM + automaton + the
  costate map lam_BE = (2 Tmax/c) lam_legacy in one shot.
- G2 ztl_A vs FD-of-field: 1.4e-10 (worst regime).
- G3 variational STM vs FD-of-flow (0.05 arc): 3.7e-9.
- G4 events at eps=0.5: 55/55 reproducible to 3.7e-9 under tolerance change,
  S at boundaries to 4.5e-13, u continuous to 4.5e-13.
- G5 saltation STM across one eps=0 switch: 3.1e-10 (|Sdot| = 74 — healthy).

**Acceptance test (the campaign fork): running** — Newton with the exact
variational-STM J from the banked P0i floor iterate (75 mN, eps=1,
||R|| = 1.56e-2 where CS and FD Jacobians both failed). Gate: 1e-8 in <= 30
iterations. PASS = first converged indirect multi-rev solve of the campaign
= the Z3 ladder anchor.

### Z0 acceptance test result (2026-07-12): FAIL at the pre-registered gate — and the verdict is sharp

Part 1 (`z0_accept_75mN.m`, plain Newton + Armijo, exact STM J): sanity gate
PASSED perfectly (mapped seed reproduces the P0i floor: 1.5601e-2 vs
1.562e-2 — EOM + automaton + costate map validated end-to-end). But the
solve crawled: 1.56e-2 -> 1.526e-2 in 30 iters, alpha pinned at 0.002,
cond(Jeq) ~ 7e8 (the FD J's benign 1.7e4 was secant-smoothing artifact).

Part 2 (`z0_accept2_lm.m`, LM + exact STM J — the Zhang configuration;
amendment: step strategy only, same J source, same 1e-8 gate): 1500
iterations, 100 min: **1.526e-2 -> 4.95e-3, then capped, decelerating, no
quadratic burst.** The telling detail: LM's damping stayed ~1e-13 (near-full
GN steps, ACCEPTED — with the corrupted CS/FD J's they were rejected), yet
each step bought only ~0.1%. Exact linear model, trusted and followed,
still only inches along the valley.

**VERDICT (the campaign fork, honestly settled): exact derivatives are
NECESSARY (they turned rejected steps into accepted ones and 4x'd the
distance covered) but NOT SUFFICIENT. The remaining wall is the narrow
curved valley geometry of 13-rev SINGLE shooting itself.** Two pre-scoped
escapes:
1. **Multiple shooting with Z0's exact per-arc STMs** (the IFS architecture
   with the derivative quality it always lacked), seeded by CHOPPING the
   current best iterate (a genuinely integrated trajectory, so the MS seed
   is dynamically consistent by construction). Arcs of ~1 rev collapse the
   per-arc amplification and straighten the valley. This is the natural Z1.
2. Direct-side dual seed (start inside the quadratic basin).

Assets banked: `z0_accept_trace.mat` / `z0_accept2_trace.mat` (best 75 mN
iterate ||R|| = 4.95e-3 in BE convention — the Z1-MS seed).

## Z1 — multiple shooting (built 2026-07-13)

Files: `ztl_ms_residual.m` (MS residual + exact block-bidiagonal Jacobian
from Z0's per-arc STMs), `ztl_ms_seed.m` (dynamically-consistent seed by
chopping a single-shooting trajectory -> continuity residual EXACTLY 0 at the
seed), `ztl_ms_nodes.m` (uniform / amplification node placement),
`ztl_ms_solve.m` (column-scaled augmented-QR LM), `ztl_ms_solve_trr.m`
(lsqnonlin TRR), `test_ztl_z1.m` (gates), `z1_run_75mN.m` (driver).

**Unit gates 3/3 PASS:** H1 seed on the flow (maxCont = 0, termErr = SS floor
4.947e-3 exactly); H2 block Jacobian vs FD shrinks with arc length
(M=6 2.8e-4 -> M=26 2.3e-7, confirming the block placement is exact and the
M=6 gap was long-arc FD-reference degradation); H3 cond(J_MS) < cond(J_SS)/10.

**Three findings from bringing it up:**

1. **Node placement is NOT the conditioning lever.** cond(J_MS) plateaus at
   ~1e11 for uniform M=26..208 (does not fall with more nodes). Per-arc STM
   norms are wildly non-uniform (median 15, max 1.85e5 at an early-perigee
   arc -- 12000x). But amplification-EQUIDISTRIBUTING placement (equidistribute
   log||Phi(t,0)||) did NOT help (cond stayed ~1e12, ratio still ~5000x): the
   perigee amplification is INTRINSIC and un-splittable -- you cannot divide a
   perigee pass into equal-amplification sub-arcs. (int||A||dt placement was
   worse still, ~1e16: it packs nodes into the perigee instant.) Kept as an
   option; default uniform.

2. **THE conditioning lever is Jacobian EQUILIBRATION, and it is large.** The
   MS Jacobian is artificially ill-conditioned by the disparate scales of the
   state vs costate columns (and the residual-component rows). Two-sided
   (Ruiz) equilibration drops cond 3.4e11 -> 5.0e6 at M=104 (FOUR-plus orders);
   column-only drops it to 8.4e8. This likely explains a large part of the
   whole indirect campaign's crawl -- the shooting Jacobian was never
   equilibrated.
   - CAVEAT (a real bug, then fixed): ROW-scaling the residual changes the
     objective, so its LM step is an ascent direction for the true ||R|| far
     from the solution -- the two-sided solver stalled instantly (mu -> 1e14,
     every step rejected). The safe form is COLUMN-only reparametrization
     dz = Dc w (objective-preserving, textbook Marquardt), cond ~8e8,
     guaranteed descent.

3. **With column-scaled LM at M=104, MS BREAKS the single-shooting wall
   (the headline Z1 result).** From the banked Z0 seed (4.947e-3, where SS
   floored dead and TRR stalled at 4.78e-3), the objective-preserving
   column-scaled LM descends monotonically: 4.95e-3 -> 9.97e-4 (it 15) ->
   9.2e-5 (it 67) -> 6.2e-5 (it 79) -- an ~80x break past the ~1e-3 indirect
   floor that stopped EVERY prior method in the whole campaign
   (ms_band, IFS, Z0 single-shooting). The terminal BC converges well
   (termErr 3.3e-6 at it 79); the remaining residual is the CONTINUITY block,
   which plateaus at maxCont ~3e-5.

4. **The residual floor is the stiff costate-continuity directions, and the
   solver-side fix is still open.** maxCont ~3e-5 sits far above integration
   tol (1e-13), so it is a real un-reduced residual: the near-null (stiff)
   directions of the perigee-continuity rows. Four solver variants were
   characterized:
   - column-scaled LM (objective-preserving): BEST -- monotone to 6.2e-5,
     floors on continuity. Default.
   - two-sided-equilibrated LM (row+col): row-scaling the residual changes
     the objective -> ascent steps, instant stall (mu->1e14). Rejected.
   - equilibrated-Newton + Armijo: exact Newton direction (cond 3.5e7) but the
     full step overshoots the seed (13-rev sensitivity), backtracks to tiny
     alpha -> crawls (4.95e-3 -> 4.80e-3 / 13 it). Rejected.
   - weighted LM, fixed Dr = 1/rownorm(J_seed) (Deuflhard): the weighted norm
     converges (5.3e-3 -> 1.5e-4) but Dr DOWN-weights the stiff perigee rows
     (rownorm ~1e5), so true continuity BALLOONS to 2e-2. Wrong weighting.
   Candidate next levers (not yet built): (a) a genuine trust-region on the
   TRUE objective with the equilibrated Newton direction (dogleg/2-D subspace,
   controlling step length not just backtracking a fixed direction);
   (b) node refinement targeting the stiff directions specifically;
   (c) a physically-scaled residual DEFINITION (scale costate-continuity vs
   state-continuity by characteristic magnitudes) rather than rownorm.

**Z1 STATUS: multiple shooting is built, unit-gated (3/3), and demonstrably
BREAKS the single-shooting wall (80x past the campaign-wide ~1e-3 indirect
floor, monotone) -- the qualitative thesis is PROVEN. Full 1e-8 convergence
is not yet reached; it is now a bounded, well-characterized solver-tuning
problem on the stiff costate-continuity directions (candidate levers above),
not a wall. The machinery (exact block-STM residual, dynamically-consistent
seed, equilibration) is the reusable substrate for Z3/Z4.**

### Z1 trust-region globalization (2026-07-13)

Built `ztl_ms_solve_tr.m` -- SVD-based Levenberg-Marquardt trust region
(More-style) on the column-scaled Jacobian: one SVD per accepted iterate
makes the boundary mu-search closed-form/free; Delta adapts by the true-
objective reduction ratio, interpolating between short globalizing steps and
the full Newton step. `z1_resume_tr.m` warm-restarts it (resumable).

**Result (M=104, from the banked Z0 seed):** the TR is the best method yet --
descends monotonically 4.95e-3 -> 1.44e-5 (200 it) -> <1e-5 (resume),
**flag=0 (iteration cap, NOT stalled)**, with terminal BC essentially
converged (termErr ~2e-7) and the residual now dominated by continuity
(maxCont ~4.4e-6). It is the FIRST method in the whole campaign to converge
monotonically past 1e-5 with no plateau. BUT the convergence rate DEGRADES as
||R|| falls (1%/it -> 0.5%/it), and Delta stays pinned small (~0.01): the
stiff perigee-continuity nonlinearity caps the productive step, so it
asymptotes slowly and does NOT reach 1e-8 at M=104 in practical iteration
counts (extrapolates ~1e3+ iters, rate degrading).

**Re-chop is a DEAD END, and the reason is illuminating.** At the TR iterate
(MS ||R||~1e-5, gaps ~4.6e-6), the SINGLE-SHOOTING residual at the extracted
lam0 is **2.56** -- enormous. The tiny continuity gaps are LOAD-BEARING: they
absorb the 13-rev amplification that single shooting cannot represent.
Re-integrating lam0 gap-free lands back at ||R||~2.6 (worse than the seed).
This is the cleanest proof yet of WHY multiple shooting is necessary here (and
why single shooting floored dead) -- and why closing the gaps to zero is
intrinsically the slow part.

**Z1 status after the TR push:** multiple shooting + SVD trust region is the
first method to CONVERGE (monotone, no stall) on this 13-rev indirect problem,
reaching ~1e-5 with the terminal BC solved and only the stiff continuity
remaining. It asymptotes short of 1e-8 at M=104 (nonlinearity caps Delta).
The clear next lever (flagged from the start as insurance) is PERIGEE-
CONCENTRATED nodes: shorter arcs through the perigees reduce the per-arc
amplification/nonlinearity, which should let Delta grow and the gaps close --
targeting the floor exactly where the re-chop test localized it. (Uniform-node
scaling of the continuity floor: M=26 ~3.9e-4, M=104 ~1.2e-5 -> ~7e-6; uniform
alone won't reach 1e-8 at feasible M, but perigee-concentrated should be far
more efficient.)

### Z1 perigee-concentrated nodes (2026-07-13)

Added a `perigee` placement to `ztl_ms_nodes.m` (equidistribute w = 1/r1^p,
r1 = distance to Earth, on the integrator's native perigee-dense output;
`probe_perigee.m`). p=1 is the sweet spot: it flattens the per-arc
amplification profile from max/med **2918 (uniform) -> 42**, worst-arc ||Phi||
1.77e4 -> 4.44e3 (higher p over-concentrates into the perigee instant and is
worse). Combined with the trust region (M=104):
- **Helped the EARLY descent:** the trust radius grew to Delta=0.125 (vs 0.016
  pinned for uniform) and held through ~it 37, descending 4.9e-3 -> 1.2e-3
  much faster than uniform.
- **Did NOT crack the late-stage asymptote:** near ~2e-5 Delta shrinks again
  (0.008) and the rate degrades to ~1%/it, same as uniform. Reason: the WORST
  perigee arc's amplification (~4e3) is INTRINSIC and barely responds to node
  count -- `probe_perigee_M.m`: worst ||Phi|| = 4.44e3 (M=104) -> 3.96e3 (156)
  -> 3.71e3 (208) -> 3.18e3 (M=312). The single closest-approach passage
  cannot be split by time-based nodes (the same "perigee amplification is
  un-splittable" lesson as the conditioning, now for LINEARITY/curvature).

**Read:** perigee nodes are a real but partial win (faster early descent,
lower unknown-count needed than uniform for the same progress). The remaining
barrier is the intrinsic CURVATURE of the closest-perigee arc in PHYSICAL
time -- which is exactly why the DIRECT side (PSR) uses SUNDMAN regularization
(dt/dtau = r^1.5, stretching the perigee passage so per-tau dynamics are
gentle). The principled structural fix for Z1's last two orders is
**Sundman-regularized arcs** (PLAN_PRONG_Z's named fallback for perigee
stiffness), consistent with the whole campaign's experience. Meanwhile the
perigee+TR run continues toward its floor as the practical anchor.
