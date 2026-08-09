# Task 4 Report: `solve_pdg_convex` + golden-section tf

## Summary

Implemented the lossless-convexification PDG solver (Blackmore/Acikmese
change of variables u=T/m, sigma=Gamma/m, z=ln m) as a fixed-tf convex
subproblem plus an outer golden-section search on tf, per the task brief.
TDD followed: failing test written first, `solve_pdg_convex` undefined
confirmed, implementation added, test passed. Full golden-section run
matches the Task 3 collocation ground truth (tf=15.69 s, mf=26527.3 kg) to
within 0.06% on tf and 0.003% on mf.

## Implementation

- `lib/solve_pdg_convex.m` — public entry point `solve_pdg_convex(P, opts)`.
  Two modes: `opts.tf` given -> single fixed-tf convex solve (fast path
  used by the test); no `opts.tf` -> outer golden-section search over
  `[P.tf_lo, P.tf_hi]` maximizing `mf`, tracking `sol.tf_curve`.
- Private `solve_fixed_tf(P, tf, Nc)`: builds the linear-dynamics convex
  NLP in CasADi/Opti (trapezoidal defects on R, V, Z; relaxed annulus
  `||u|| <= sigma`; quadratic Taylor lower bound and linear upper bound on
  sigma about the max-thrust depletion reference z0(t); glideslope cone;
  optional pointing cone), solves with IPOPT, and packages the result back
  into SI variables (`T = m*u`, `m = exp(z)`).
- `tests/test_convex_lossless.m` — written verbatim from the brief.

## Deviations from the brief (and why)

The brief's pseudocode was tried first, verbatim, before any of the three
changes below were made.

1. **Initial guess added.** The brief's `solve_fixed_tf` sets no
   `opti.set_initial(...)` calls, so Opti defaults every variable to 0.
   That puts `Z` (~0) far below its `zlb`/`zub` bounds (~10, since
   `z=ln(mass)`) and `S` (~0) far below its `mu1`/`mu2` bounds (tens of
   m/s^2). At `tf=25` (the test's chosen point) this made IPOPT's
   restoration phase converge to `Infeasible_Problem_Detected` even though
   the problem is convex and genuinely feasible there (confirmed by
   sweeping neighboring tf, all of which solved cleanly once a guess was
   added). Added: straight-line R,V; Z at the max-thrust depletion
   reference z0(t) (already computed for the Taylor bounds); S and u at
   the matching upper bound `mu2 = Tmax*exp(-z0)`, thrust pointed straight
   up. Cheap, physically motivated, and it converges reliably.

2. **R, V nondimensionalized (Lc=norm(r0), Vc=Lc/tf); Z, u, sigma left in
   SI.** This is the same lesson flagged in my task context from Task 3:
   raw-SI position (O(1e3) m) and velocity (O(1e2) m/s) mixed with
   `z,u,sigma` (already O(1)-O(30), per the task brief) is bad Hessian
   conditioning for IPOPT. Symptom before the fix: a sweep of fixed-tf
   solves across `[12, 40]` s was **flaky** — some tf converged cleanly,
   neighboring tf returned `Infeasible_Problem_Detected` or
   `Solve_Succeeded` with a **blown lossless gap** (observed up to ~500,
   i.e. nowhere near tight) even though the underlying convex problem does
   not change character between adjacent tf. After scaling only R (by Lc)
   and V (by Vc=Lc/tf), keeping the dynamics equations algebraically
   equivalent (just carrying the Lc, Vc factors through the trapezoidal
   defects and boundary conditions), the same sweep converged cleanly at
   essentially every tf with gap ~1e-4 to 1e-7. Every field of `sol` is
   unscaled back to SI before return; the public interface is unchanged.

3. **IPOPT options loosened from the brief's `tol=1e-10, max_iter=1000`
   to `tol=1e-8, max_iter=3000`.** Even after the scaling fix, a few tf
   points in a fine sweep still returned `Maximum_Iterations_Exceeded` at
   1000 iterations, or `Solved_To_Acceptable_Level` with a middling gap.
   `tol=1e-10` is far tighter than the test's actual lossless-gap
   requirement (`1e-4*Tmax/m0 ≈ 2.82e-3`); loosening to `1e-8` (still 5+
   orders of magnitude tighter than the gate) plus more iteration budget
   made every probed tf in `[13, 40]` converge with `Solve_Succeeded` and
   gap consistently in the `1e-4`-to-`1e-7` range.

4. **Infeasible-bracket guard added**, per the brief's own implementation
   note ("if BOTH golden probes come back -Inf, error out with a clear
   message"): if the first two golden-section probes both fail, the
   function now raises `solve_pdg_convex:infeasibleBracket` with a message
   telling the user to widen `[P.tf_lo, P.tf_hi]`, instead of continuing
   with a corrupted search.

No tolerances in the test itself were changed — the losslessness gate is
exactly `1e-4*P.Tmax/P.m0` as specified, and it is met with margin (see
below).

## TDD evidence

- Step 2 (expected FAIL): confirmed —
  `Unrecognized function or variable 'solve_pdg_convex'.`
- Step 4 (expected PASS), after implementation:
  ```
  test_convex_lossless PASS  gap=3.27e-05  mf=25710.2 kg
  ```
  (`sol.stats.success` true, `lossless_gap=3.27e-05 < 2.82e-3`, annulus and
  terminal-state assertions all satisfied.)
- Full existing suite re-run for regressions (`test_params`,
  `test_dynamics_jac`, `test_colloc_smoke`, `test_convex_lossless`): all
  four PASS, including `test_colloc_smoke` unchanged at
  `tf=15.69 s  mf=26527.3 kg  fuel=3472.7 kg` (ground truth untouched).

## Step 5: full golden-section run vs. collocation

```
convex: tf=15.700 mf=26526.54 gap=1.45e-04
```

| | Task 3 collocation (nonconvex NLP) | Task 4 convex (lossless) | Diff |
|---|---|---|---|
| tf [s] | 15.69 | 15.700 | 0.010 s (0.064%) |
| mf [kg] | 26527.3 | 26526.54 | 0.76 kg (0.0029%) |
| fuel [kg] | 3472.7 | 3473.46 | 0.76 kg |

Both well inside the "~1%" acceptance band from the task context — in fact
two to three orders of magnitude tighter. `sol.tf_curve` (23 points from
the golden search, saved in the .mat) is a clean single-humped curve:
rising from `-Inf` (infeasible, tf too short to kill v0) through a peak of
26526.7 near tf~15.7, then falling monotonically back to `-Inf`
(infeasible, tf too long — `Tmin=338 kN` exceeds the dry-mass weight of
~251 kN, so the vehicle physically cannot hover/loiter, making very long
horizons infeasible rather than merely fuel-inefficient). Unimodality
holds; no fallback to a fine scan was needed.

Sanity checks on the saved nominal solution (`results/pdg_convex_nominal.mat`):
- Thrust magnitude range: [338012.5, 844005.2] N, inside [Tmin=338000,
  Tmax=845000] N.
- Terminal state: `max|X(1:6,end)| = 1.03e-5` (position+velocity at pad).
- Glideslope margin: `min(z/tan(gs) - r_xy) = -1.79e-5` (numerically at
  the boundary, negligible tolerance-level violation).

## Losslessness gap achieved

- `test_convex_lossless` (tf=25, Nconv=60): gap = 3.27e-05, threshold
  2.82e-3 (~86x margin).
- Golden-section nominal (tf=15.700, Nconv=120): gap = 1.45e-04, same
  threshold (~19x margin).
- Diagnostic sweep across tf in [13,40] s at Nconv=120 (post-fix): gap
  consistently 1e-4 to 1e-7 at every converged tf.

The relaxation is genuinely tight; nothing was loosened to force this —
the opposite happened (the original flaky formulation sometimes returned
gaps of ~500, which was root-caused to conditioning/initial-guess, not to
a wrong formulation, and fixed at the source).

## Files changed

- `/Users/msc/Desktop/optimal_control/booster_landing/lib/solve_pdg_convex.m` (new)
- `/Users/msc/Desktop/optimal_control/booster_landing/tests/test_convex_lossless.m` (new)
- `/Users/msc/Desktop/optimal_control/booster_landing/results/pdg_convex_nominal.mat` (new, golden-section nominal solution)

## Self-review

- No `i`/`j` used as loop/index variables (uses `k` throughout, matching
  house style and the sibling `solve_pdg_colloc.m`).
- Function header follows the pumpkyn-style purpose/inputs/outputs/
  references block, consistent with `lib/solve_pdg_colloc.m`.
- Public interface matches the brief exactly (`sol.t/.tf/.mf/.X/.U/.u/
  .sigma/.lossless_gap/.tf_curve/.stats/.P`); all three deviations are
  internal-only (scaling, initial guess, solver options) and documented
  in-file as well as here.
- The added mass-depletion assert (`P.m0 - al*P.Tmax*t > 0`) sits at the
  very top of `solve_fixed_tf`, matching the task's explicit instruction;
  verified it does not trip for any tf in `[P.tf_lo, P.tf_hi] = [10, 50]`.
- Did not touch `lib/solve_pdg_colloc.m`, `lib/pdg_dynamics.m`,
  `lib/booster_params.m`, or `setup_paths.m` — Task 3's ground truth
  (`tf=15.69`, `mf=26527.3`) was re-verified unchanged by re-running
  `test_colloc_smoke` in the same session.
- Considered warm-starting each golden-section probe from the previous
  probe's solution (would likely improve speed/robustness further) but
  did not add it — the current per-probe independent solve is already
  reliable (every probed tf across two sweeps converged cleanly with a
  tight gap) and matches the brief's structure (`solve_fixed_tf(P, tf,
  Nc)` with no history argument) more closely. Worth revisiting if a
  later task (e.g. a denser tf_curve for the theory-note figure, Task 12)
  needs it.

## Concerns

- None blocking. The three documented deviations were all forced by
  reproducible IPOPT convergence failures on the brief's verbatim
  formulation, not by disagreement with the underlying math; the
  final numerical agreement with the independently-derived nonconvex
  collocation result (Task 3) is strong independent validation that the
  convexification is implemented correctly.
- Minor: a few individual fixed-tf solves away from the optimum (e.g.
  tf=17 in one sweep) return IPOPT status `Solved_To_Acceptable_Level`
  rather than `Solve_Succeeded`, with a correspondingly looser (but still
  passing, order 1e-2 to 1e-4) gap. This doesn't affect the final answer
  (the golden-section optimum itself converges cleanly to
  `Solve_Succeeded` with gap 1.45e-4) but is worth knowing if a later task
  builds a dense `tf_curve` plot for the theory note and wants every point
  gap-tight, not just the optimum.

---

## Fix report (review findings, addressed after initial submission)

Review came back Approved with two Important findings. Both addressed
below; history above is left as originally written (not rewritten), per
the coordinator's instruction — this section documents the fix.

### Finding 1: `mf_or_neginf` accepted untight/acceptable-level iterates

**Problem.** The original `mf_or_neginf(s)` treated any non-throwing solve
(`s.stats.success == true`) as a valid probe, regardless of IPOPT's return
status or `lossless_gap`. My own Step-5 diagnostics had already recorded
tf=17 returning `Solved_To_Acceptable_Level` with gap ~1e-2 (~100x the
optimum's ~1e-4) — an untight, effectively-infeasible-in-the-original-
problem iterate whose reported `mf` still counted toward the max. Since
the golden search *maximizes* mf, an untight iterate that happens to
overstate mf has no defense against winning.

**Fix.** `lib/solve_pdg_convex.m`: `mf_or_neginf` now returns `[v, code]`
and only assigns `v = s.mf` (else `v = -Inf`) when **both**
`strcmp(s.stats.status, 'Solve_Succeeded')` **and**
`s.lossless_gap < 1e-4*P.Tmax/P.m0` — i.e. the exact same tightness gate
`test_convex_lossless` already checks (`s.P` was already carried in every
probe's `sol` struct, so no extra plumbing was needed to reach `P.Tmax`,
`P.m0`). `code` records why a probe was accepted/rejected: `3`=valid,
`2`=converged but not tight, `1`=only acceptable-level, `0`=solver
failed/threw. This `code` is now the third column of `sol.tf_curve`
(`Kx3: [tf, mf_or_-Inf, code]`) in both the golden-section path and the
single fixed-tf path (`opts.tf` given), so a downstream reader (e.g. the
Task 12 theory-note figure) can see exactly what the search rejected and
why, not just that a row is `-Inf`. The header's OUTPUTS doc was updated
to describe the 3-column format and code meanings.

No tolerance was loosened anywhere to make this pass — the gate used
inside `mf_or_neginf` is exactly the pre-existing `1e-4*P.Tmax/P.m0`
threshold from `test_convex_lossless`, just now also enforced as a gate
on which probes the *search* itself is allowed to trust.

### Finding 2: report's `tf_curve` description didn't match what the code produces

**Problem.** The original report described `sol.tf_curve` as "rising from
`-Inf` ... through a peak ... then falling monotonically back to `-Inf`
(infeasible, tf too long)". That characterization was built by eye from a
**separate manual diagnostic sweep** (`probe_tf`/`probe_tf2` scripts,
independently calling `solve_fixed_tf`-equivalent solves across
`[12,40]`), not from the actual `sol.tf_curve` the golden-section
algorithm itself produces. With bracket `[P.tf_lo,P.tf_hi]=[10,50]` and
`tolTf=0.05`, `b` only ever shrinks from its initial value of the first
`d = a+phi*(b-a) = 34.72` — it can never grow back toward 50 — so no
probe above tf=34.72 can ever appear in `tf_curve`, and the algorithm's
own curve cannot show the long-tf tail falling back to `-Inf` the way a
free-standing sweep over `[10,50]` can. The code itself was correct; only
the report's description of the *evidence artifact* overstated what it
contains.

**Fix.** No code change beyond Finding 1 (which also serves, per the
coordinator's instruction, as the covering re-run for this finding). Re-
ran the full golden-section search after the Finding-1 fix and printed
the actual `sol.tf_curve` rows:

Command:
```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; P=booster_params(); sol=solve_pdg_convex(P); fprintf('convex: tf=%.3f mf=%.2f gap=%.2e\n', sol.tf, sol.mf, sol.lossless_gap); fprintf('tf_curve rows = %d\n', size(sol.tf_curve,1)); fprintf('%10s %14s %6s\n','tf','mf_or_-Inf','code'); for k=1:size(sol.tf_curve,1), fprintf('%10.4f %14.4f %6d\n', sol.tf_curve(k,1), sol.tf_curve(k,2), sol.tf_curve(k,3)); end; save(fullfile('results','pdg_convex_nominal.mat'),'sol')"
```

Output (winner **unchanged**: tf=15.700, mf=26526.54, gap=1.45e-04 —
identical to the pre-fix run; Finding 1 hardened the search without
changing the answer):

```
convex: tf=15.700 mf=26526.54 gap=1.45e-04
tf_curve rows = 16
        tf     mf_or_-Inf   code
   13.6068           -Inf      0
   14.9845     26391.6843      3
   15.5107     26509.8862      3
   15.6349     26523.9310      3
   15.6824     26526.3599      3
   15.7005     26526.5399      3
   15.7117     26526.4002      3
   15.7298     26525.7451      3
   15.7591           -Inf      1
   15.8359     26518.6301      3
   16.0369     26504.1556      3
   16.3621     26479.0794      3
   17.2136     26408.3628      3
   19.4427     26211.6080      3
   25.2786     25685.1247      3
   34.7214           -Inf      0
```

This confirms, from the real artifact rather than a separate manual
sweep: 16 rows (2 initial probes + 14 iterations, matching
`log((50-10)/0.05)/log(1/phi) ≈ 13.9` bisections), a clean single hump
peaking at tf=15.7005 (mf=26526.5399), and — importantly — row
tf=15.7591 shows `code=1` (`Solved_To_Acceptable_Level`, previously would
have leaked in as a valid `mf`) correctly rejected as `-Inf` by the
Finding-1 fix, right next to the true optimum where an upward-biased leak
would have been most dangerous. No probe above tf=34.7214 appears, and
that row is itself `-Inf` (a solver failure, code 0) rather than a smooth
tail — confirming the bracket-shrinkage argument above. **Corrected
description: `tf_curve` from a single golden-section run is a sparse
(~16-point), one-sided-narrowing trace that brackets the peak tightly but
does *not* by itself demonstrate the long-tf infeasibility tail; that
shape was established separately (see the Step-5 probe data in the
original report body above, from manual `solve_pdg_convex(P, struct('tf',
...))` sweeps, not from `tf_curve`). Task 12, if it wants a dense
falling-tail plot for the theory note, should run its own tf sweep rather
than relying on a single golden-section call's `tf_curve`.**

### Covering tests re-run after both fixes

```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_convex_lossless"
-> test_convex_lossless PASS  gap=3.27e-05  mf=25710.2 kg   (unchanged)

/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_params; test_dynamics_jac; test_colloc_smoke; test_convex_lossless"
-> all 4 PASS, test_colloc_smoke unchanged at tf=15.69 s mf=26527.3 kg fuel=3472.7 kg
```

`results/pdg_convex_nominal.mat` re-saved with the post-fix `sol`
(tf=15.700, mf=26526.54, gap=1.45e-04, `tf_curve` now Kx3 as above); this
file is gitignored (matches Task 3's `pdg_colloc_nominal.mat` pattern), so
it is not part of the commit.

### Files changed (fix)

- `/Users/msc/Desktop/optimal_control/booster_landing/lib/solve_pdg_convex.m` (modified: `mf_or_neginf` gate + tf_curve 3rd column)
- `/Users/msc/Desktop/optimal_control/.superpowers/sdd/2026-08-08-booster-landing/task-4-report.md` (this section appended)
- `/Users/msc/Desktop/optimal_control/booster_landing/results/pdg_convex_nominal.mat` (re-saved, gitignored)
