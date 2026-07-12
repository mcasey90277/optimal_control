# IFS Results — Rung 2 (full 1.12× gate)

## 1. Status

IFS machinery is **built and validated**. Full multi-switch convergence at the
Rung-2 gate (the 1.12×, k=10-switch, rendezvous-mode solve) is **BLOCKED** by
terminal-cluster shooting conditioning — a characterized negative result, not
a code bug. Two robustness guards were added post-mortem (see §3d and the
commit) so a future re-run cannot crash mid-solve or lose a completed solve to
a certify-time error; they do not change the conditioning verdict below.

## 2. Machinery validated

Unit tests (`test_ifs_residual`, `test_ifs_jacobian`, `test_ifs_seed`) — **ALL
PASS**:
- EOM matches `sms_eom` (gap 0.0, H_sigma 8.88e-15) — Task 1.
- Residual is zero on a continuous ground-truth trajectory: `||R_fixedState|| =
  2.35e-15` — Task 2.
- Complex-step vs finite-difference Jacobian, covering `fixedState` (k=1),
  `rendezvous` (k=1), and `rendezvous` k=2 (middle-arc dependency): rel errors
  7.71e-7 / 7.71e-7 / 7.64e-6 — Task 3.
- Seed builder (`ifs_seed`) constructs a valid Z from the direct/PSR `.mat` for
  both modes: window k=1 seed residual **3.07e-3**, full k=10 seed residual
  **1.96e+00** — Task 4.

The window seed residual (3.07e-3) being nearly self-consistent at the seed was
a good rung-1 omen for the mode-'d' dual→costate map; the full problem's larger
seed residual (1.96) is expected — driving it down is exactly IFS's job.

## 3. The four failure modes diagnosed & fixed

This is the core finding of the build. Four independent failure modes were hit,
diagnosed with `systematic-debugging`, and fixed (or in one case, correctly
ruled out of scope) in order:

**(a) λ_m gauge rank-deficiency in the `fixedState` WINDOW.**
The originally-planned Rung-1 ground-truth test (a one-switch interior window
with `fixedState` terminal BC) is structurally rank-deficient: mass dynamics
are costate-independent (`dm/dτ = -κ·u·Tmax/c`), λ_m decouples identically on
coast arcs, and `fixedState` fixes the end mass rather than imposing a λ_m(τ_f)
= 0 transversality — leaving λ_m an unconstrained gauge freedom. Diagnosed via
SVD: window smallest singular value **1.5e-16**, isolated ~9 orders of
magnitude below its neighbor (2.53e-7) — the textbook signature of an exact
rank-1 null space, not a smooth conditioning tail. The τ-column complex-step
hypothesis was tested and refuted (CS Jacobian accurate to 2.5e-9). The FULL
problem (`rendezvous` terminal, which does carry λ_m(τ_f)=0) was checked
separately and found full-rank: all 178 singular values healthy, smallest/
largest ratio 2.85e-12 but no isolated near-zero outlier, `rank(J)=178/178` at
MATLAB's default tolerance, GN-consistency `resInRange = 5.35e-11`. **User-
approved decision:** drop the window as the gate; repoint the gate directly to
the full 1.12× solve. `ifs_solve`/`ifs_certify` code was correct throughout —
the window was a bad test input, not a code defect. Full diagnostics:
`.superpowers/sdd/ifs-diag-floor-report.md`, `ifs-diag-fullrank-report.md`.

**(b) Last switch crossing the fixed τ_f → negative final arc.**
The first full-gate attempt (unbounded switch-time parameterization) stalled
at ||R|| = 1.960 → 0.0266 over 400 iterations with a decelerating LM crawl
(cond(J) ≈ 3.5e11). Root-cause diagnosis (Phase 4, `ifs-diag-stall-report.md`)
found switch 10 had migrated past the fixed terminal boundary τ_f
(151.791 > 151.684), inverting the sign law on the final arc (`cert.signViol =
1.12`, order-1) and producing a negative-length terminal arc
(`cert.minArcLen = -0.108`). This is anchored to a genuine structural
degeneracy (a terminal switch cluster: 5 switches packed into ~1.87 τ-units),
not pure unstructured ill-conditioning.

**(c) LM reordering switches → arc-reversal `ode113` failures.**
Fix attempt 1 (an independent sigmoid bound per switch time, commit
`87310c1`) removed the τ_f-crossing but let LM reorder switches during a step,
producing reversed-span `ode113` calls that fail outright. Confirmed that
*ordering*, not just bounding, is essential.

**(d) Collapsed switch gap → zero-length arc.**
Fix attempt 2 (monotone stick-breaking reparameterization of switch times —
`τ_0 < τ_1 < … < τ_k < τ_f` by construction, commit `1a999a9`) fixed (b) and
(c): the re-run descended 1.96 → 0.04 over 244 iterations with zero arc-
reversal failures, but crashed at iteration 244 when an LM step collapsed a
switch gap to exactly zero, and `ode113` errors on a zero-length tspan. Fixed
in `ifs_residual.m`'s per-arc integrator helper (`ifs_int_arc`, commit
`16a506e`): spans with `real(span(2)-span(1)) <= 1e-13` short-circuit to
endpoint = start instead of calling `ode113`. **This same class of bug was
still present in `ifs_certify.m`'s independent post-solve sign-law sampling
loop** (it re-integrates each arc separately to check S(τ) at the midpoint) —
confirmed by the actual crash log (`ifs_1p12.log`, tail):
```
Error using odearguments
The last entry in tspan must be different from the first entry.
Error in ode113 ... Error in ifs_certify (line 41) ... Error in run_ifs_1p12 (line 13)
```
Closed out in this pass: `ifs_certify.m` now skips (near-)zero-length spans
in that loop with the same 1e-13 guard, and `run_ifs_1p12.m` now saves the
solve (`out`) to `ifs_1p12_results.mat` immediately after `ifs_solve` returns,
before calling `ifs_certify`, and wraps the certify call in `try/catch` so a
certify-time error cannot lose a completed solve. The final summary `fprintf`
tolerates `cert` being empty.

## 4. Full 1.12× gate result

With all four fixes in place, the final gate run (`ifs_1p12.log`) ran the full
400-iteration LM budget cleanly — **zero integration failures during the
solve itself**:

- seed `||R||` = **1.960465** → final `||R||` = **0.0232** (`ifs_solve`
  report: `k=10 ||R0||=1.96e+00 ||R||=2.32e-02 iters=400 flag=0`)
- `lsqnonlin` stopped on the iteration-limit (`flag=0`), not convergence
  (`success=0`); target is `<1e-8` — **6 orders of magnitude short**.
- Descent had decelerated to ~0.2%/iter by the end of the budget; first-order
  optimality had frozen around 0.65 (not decreasing further).
- Root barrier, diagnosed independently (`ifs-diag-stall-report.md`,
  earlier checkpoint of the same crawl): a **terminal switch cluster** — 5 of
  the 10 switches packed into only ~1.87 τ-units near the end of the transfer,
  with tiny, highly sensitive inter-switch gaps (min gap ~8.6e-3 τ at one
  checkpoint). `cond(J) ≈ 5.9e9` there; numerical rank 170/178 (8 near-null
  directions living in that cluster).
- GN-consistency at the stall: `resInRange ≈ 8.3e-12` to `5.4e-11` (both
  checkpoints) — the residual is fully reachable from the Jacobian's range;
  the Gauss-Newton step is well-defined. This means the crawl is
  **conditioning-limited, not a singular/structural dead end**: LM's damped
  steps are simply too small, relative to the terminal cluster's sensitivity,
  to close 6 orders of magnitude in 400 iterations.
- This crawl signature — well-posed Jacobian, GN-consistent residual, slow
  ill-conditioned descent anchored to a tightly packed switch/perigee cluster
  — **echoes the `ms_band` conditioning wall** (`msband-indirect-campaign`
  memory) that motivated building IFS's switch-explicit, saltation-free
  formulation in the first place. IFS removed the *smoothing* wall (no ε, no
  1/ε layers) and the single-shooting perigee-sensitivity wall (short arcs via
  multiple shooting) — but a third, narrower conditioning wall specific to
  *this* switch pattern's terminal cluster remains.

## 5. Next lever (not yet done)

**Primary:** tighter multiple shooting — add interior shooting nodes *within*
the long/sensitive arcs (perigee crossings, and especially inside the terminal
switch cluster), not only at existing switches. This directly shortens the
sensitive arcs that drive `cond(J)` up, the same cure that made switches-as-
nodes work for the rest of the trajectory.

**Lighter alternatives, cheaper to try first:**
- Row/column scaling of the residual/Jacobian by state-vs-costate magnitude
  (LM already uses `ScaleProblem='jacobian'`, but a manual physical scaling of
  the state block vs. the costate block may separate the terminal cluster's
  sensitivity from the well-behaved perigee arcs).
- A `t_f` / smoothing homotopy seed: walk the terminal-approach cluster in
  gradually (e.g., from a slightly shorter transfer, or a lightly smoothed
  switch structure) before coupling it to the full k=10 rendezvous solve, so
  LM starts closer to the cluster's true structure rather than crawling into
  it cold from the direct-solution seed.

No further gate re-runs were performed in this pass (the 400-iteration solve
takes ~30 minutes and its numbers are already recorded above and in
`.superpowers/sdd/progress.md`).

---

## Post-merge diagnostic investigation (2026-07-11)

After banking IFS with the conditioning stall characterized, we ran a focused
investigation into the stall's true cause, including an external GPT-5.6-sol
methodology review (`reviews/gpt56sol_2026-07-11.md`). **Context: PSR (the
point-4 mesh-refinement method) already works and is merged; IFS is the open
one.** Summary of what we found — three candidate culprits were tested and the
first two were ruled out:

**GPT-5.6-sol review (validated the formulation, sharpened the diagnosis):**
- Explicit-switch-node hard-throttle multiple shooting is standard and sound;
  saltation matrices correctly unnecessary. Our vanishing-arc read is credible —
  it named it a "vanishing bang/coast arc near a switching-structure fold."
- Cautions that proved decisive: (a) `cond(J)` and `R in range(J)` are weak
  evidence — confirm from the **smallest singular VECTORS of a physically-scaled
  Jacobian**; (b) the stick-breaking sigmoid parameterization DEGENERATES as a
  gap -> 0 (dtau/dg -> 0), compounding the crawl; use center/width variables;
  (c) audit whether fixing both tau_f and t(tau_f)=t_f over-constrains.

**Scaled-SVD localization (the decisive test) — near-double-root REFUTED:**
Physically-equilibrated Jacobian, smallest singular vector, at two iterates:
- **Seed** (||R||=1.96): smallest SV 1.2e-11 (near-singular, survives scaling);
  null direction 76% on **switch #4's node** (the shallowest S-crossing,
  |dS/dtau|=0.11) — NOT the tight terminal pair (only 1.8% of the mass).
- **Mid-crawl** (||R||=0.15): smallest SV 1.9e-9; null direction has MOVED — now
  83% on the **initial costate lambda_0** (dominated by the position-costate
  components lambda_r,x / lambda_r,y).
- The near-double-root pair {switches 8,9 / 9,10}: 1.8% -> ~0%. **Not the
  culprit.** (The solve does compress that pair, tau-gap 9.95e-3 -> 3.17e-3 —
  the collapse symptom that the ifs_int_arc guard handles — but it is not the
  conditioning driver.)

**tau_f over-constraint — RULED OUT:** terminal-block residual split at both
iterates has |t - t_f| ~ 4e-3-5e-7 and |lambda_m| ~ 6e-5-2e-2, both negligible
vs the continuity/switch blocks. Fixing tau_f is not implicated.

**Conclusion:** the conditioning is **diffuse and shifting** (a shallow switch at
the seed, the initial position-costate lambda_r at mid-crawl), not a single
surgical target. This is the textbook weak spot of indirect shooting — lambda_r
is only indirectly coupled to the trajectory and is weakly determined by a short
first arc — and it confirms the "echoes the ms_band conditioning wall" reading.
There is no near-coincident pair to merge and no grazing switch to regularize.

**Remaining levers (all heavier; a strategy decision, not a surgical fix):**
1. Canonical physical scaling (reciprocal-costate) + a regularized / rank-
   revealing (truncated-SVD) Newton step to move through the near-null lambda_r
   direction instead of crawling against it. Cheapest to try.
2. A better-conditioned costate seed via t_f-homotopy (converge an easier t_f
   first, continue in t_f into 1.12x).
3. Tighter multiple shooting (interior non-switch nodes) so no single arc weakly
   determines its upstream costate. The principled structural cure, biggest build.

IFS remains OPEN. PSR is the working deliverable.
