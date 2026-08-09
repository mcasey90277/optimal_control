# Task 3 Report: `solve_pdg_colloc`

## Status: DONE

## Implementation

Files created:
- `booster_landing/lib/solve_pdg_colloc.m`
- `booster_landing/tests/test_colloc_smoke.m`

Followed the brief's TDD sequence and formulation exactly (Hermite-Simpson
free-tf min-fuel NLP, thrust annulus handled directly/nonconvex, `opti.lam_g`
dual extraction per the house sign-bug lesson, `pdg_rhs_casadi` re-expressing
the RHS symbolically) **with one substantive adaptation**, documented below.

## Deviation from the brief: nondimensionalization

The brief's Step 3 code, run verbatim (raw SI units, cold start, N=15,
default `maxIter=3000`), does **not** converge. IPOPT hits
`Maximum_Iterations_Exceeded` with the iterate badly diverged (final mass
overshooting `m0`, `tf` well outside its `[10,50]` bracket, dual
infeasibility ~5.6e8). Per the brief's own contingency chain I first tried
`opts.maxIter=6000` (Step 3 implementation notes) — it got *worse*, not
better (mass moved further from physical range, more iterations spent in
what is evidently a restoration-phase struggle). The brief's next fallback
(warm-start from the Task 4 convex solver) is unavailable — Task 4 doesn't
exist yet — and the brief is explicit that "the smoke test must pass cold
(spec requirement)."

I diagnosed this empirically (not guessed): built a standalone probe with
the *identical* formulation but nondimensionalized decision variables, and
it converged to `Solve_Succeeded` in 32-51 iterations depending on the
scale choice, landing on the exact same physical optimum as the raw-SI
run's (still-improving) trend. Root cause: raw SI mixes O(1e3) m positions,
O(1e5) N thrust, and O(1e4) kg mass in one Newton step — a textbook-bad
Hessian conditioning case for a barrier-method NLP (this is exactly what
this repo's own top-level `optimal_control/CLAUDE.md` "Scaling
Considerations" section calls out, and standard doctrine in Kelly 2017 /
Betts 2010).

**Adaptation applied:** the NLP is built and solved in nondimensional
variables using characteristic scales derived from `P` itself (not
hardcoded magic numbers, so it generalizes if `P` changes):
- `Lc = norm(P.r0)` (position scale)
- `Tc = 0.5*(P.tf_lo + P.tf_hi)` (time scale, bracket midpoint)
- `Mc = P.m0` (mass scale)
- `Vc = Lc/Tc`, `Fc = Mc*Lc/Tc^2` (velocity/force scales, F=ma-consistent)

Dynamics, bounds, and the initial guess are all divided into this scaled
space before being handed to `opti`; `pdg_rhs_casadi` was extended to take
a scale struct `S` (`ghat`, `kMdot`, and the drag coefficients rescaled
consistently, though drag is off by default and untested here — that's
Task 11). Every field of the returned `sol` (`X`, `U`, `Um`, `tf`, `t`) is
unscaled back to SI before return, so **the public interface is byte-for-
byte what the brief specifies** — only the internal `Opti` construction
differs. `sol.lam_defect` is the dual of the *nondimensional* defect block
(documented in the header); no downstream consumer in this task's scope
reads it numerically.

The one other deviation is cosmetic and pre-authorized by the brief itself:
`aD = 0*v;` instead of `casadi.MX.zeros(3,1)*0` (brief's own note: "replace
with a plain `aD = 0*v;` if MX zeros is awkward").

`gDefStart`/`nDefRows` bookkeeping, constraint ordering, and the HS defect
math are unchanged from the brief.

## TDD evidence

**RED** (Step 2, before `solve_pdg_colloc.m` existed):
```
{Unrecognized function or variable 'solve_pdg_colloc'.
Error in test_colloc_smoke (line 10)
sol  = solve_pdg_colloc(P, struct('N', 15));
```

**GREEN** (Step 4, after implementation, cold start, default maxIter=3000,
51 IPOPT iterations, `Solve_Succeeded`):
```
test_colloc_smoke PASS  tf=15.69 s  mf=26527.3 kg  fuel=3472.7 kg
```
All four physics assertions passed: mass in `(mdry, m0)`, thrust annulus
`[338000.4, 844999.8]` N inside `[Tmin-1, Tmax+1]`, glideslope satisfied at
every node, terminal state `max(abs(X(1:6,end))) < 1e-3`.

## N=15 vs N=60 grid-convergence sanity (Step 5)

| grid | tf [s] | mf [kg] | fuel [kg] | IPOPT iters | status |
|------|--------|---------|-----------|-------------|--------|
| N=15 | 15.69  | 26527.3 | 3472.7    | 51 | Solve_Succeeded |
| N=60 | 15.691 | 26527.24| 3472.76   | 54 | Solve_Succeeded |

Δtf = 0.001 s (spec tolerance ~1 s), Δmf = 0.06 kg (spec tolerance ~5 kg) —
well inside spec. Saved to
`booster_landing/results/pdg_colloc_nominal.mat` (gitignored per
`results/*` in `.gitignore`, per the brief's own commit list which only
adds `lib/` and `tests/`).

Supplementary check on the saved N=60 solution (not required by the smoke
test, run for my own self-review): mass range `[26527.2, 30000.0]` kg
(bounds `mdry=25600`, `m0=30000`), thrust magnitude range
`[338000.4, 844999.8]` N (bounds `[338000, 845000]`), glideslope max
violation `-1.45e-42` (i.e. none), terminal `|X(1:6,end)| = 8.4e-43`.

## Numbers vs. the brief's sanity expectations

Brief expected tf in ~[15,40] s and fuel "a few hundred kg to ~2 t" with a
max-min-max throttle profile. tf=15.69 s lands right at the low edge of
that window (interior to the `[10,50]` bracket, not bound-active — a
genuine unconstrained optimum, not an artifact). Fuel = 3472.7-3472.76 kg
is somewhat above the "~2 t" upper end of the expectation note. I did not
force-fit this — it is the physical optimum for the given `P` (Tmin=338 kN
already exceeds the hover thrust needed at *any* mass in `[mdry,m0]`,
`P.m0*g0=294.2 kN` and `P.mdry*g0=251.1 kN`, so the vehicle cannot loiter
and the descent is necessarily brisk with real fuel cost). Reported
honestly rather than adjusted to match the expectation note; flagging for
awareness in case `booster_params` values are revisited later. I did not
inspect the throttle profile shape (bang-bang max-min-max) numerically —
brief defers that eyeball check to Task 9's viz.

## Files changed

- `/Users/msc/Desktop/optimal_control/booster_landing/lib/solve_pdg_colloc.m` (new)
- `/Users/msc/Desktop/optimal_control/booster_landing/tests/test_colloc_smoke.m` (new)
- `/Users/msc/Desktop/optimal_control/booster_landing/results/pdg_colloc_nominal.mat` (new, gitignored, not committed)

## Self-review

- Header block: full pumpkyn-style header (purpose, INPUTS with sizes,
  OUTPUTS with sizes, REFERENCES) present on `solve_pdg_colloc`; the
  private subfunction `pdg_rhs_casadi` also carries an explanatory header
  even though not strictly required for a local subfunction.
- No `i`/`j` used as loop or index variables anywhere (checked: all loops
  use `k`).
- Interface fields match the brief's spec exactly (`sol.t/.tf/.mf/.X/.U/.Um/
  .lam_defect/.stats/.P`); no fields renamed or dropped.
- `opti.lam_g` used for duals, never `opti.dual` (house lesson honored).
- Commit staged only the two brief-specified files; `results/*.mat` correctly
  excluded by `.gitignore` (`results/*`, `!results/.gitkeep`), consistent
  with the brief's own `git add` list.
- Commit message carries the actual measured numbers (tf, mf, fuel for both
  grids), not the `<fill from run>` placeholder.
- Verified the N=60 solution independently (not just N=15) against all four
  physics checks the smoke test exercises, even though the brief only
  requires the smoke-test assertions at N=15.

## Concerns

1. **Deviation from brief's literal formulation code.** The brief said "use
   verbatim except where genuine CasADi runtime issues force a documented
   adaptation." I judged the IPOPT divergence (not a CasADi API
   incompatibility, but a real numerical/conditioning failure at solve
   runtime, explicitly anticipated by the brief's own troubleshooting
   ladder) as within that allowance, especially since the brief's own two
   documented fallbacks (bump maxIter; warm-start from Task 4) were tried/
   unavailable and failed. Flagging this clearly for review since it's the
   largest deviation of the three tasks so far.
2. **Fuel number above the brief's "~2 t" sanity ceiling** (3.47 t). Traced
   to `Tmin` exceeding hover thrust at every feasible mass — a property of
   `booster_params`, not of the solver. Not adjusted; reported as-is.
3. **`sol.lam_defect` units.** These are duals of the *nondimensional*
   defect equations, not SI-unit adjoints. No consumer in Task 3's scope
   needs SI duals, but if a later task (certify, or costate-style physical
   interpretation) wants them in SI, a per-row-block rescale by
   `Mc/Tc, Lc/Tc^2*Mc`-style factors will be needed — flagging now so it
   isn't a surprise later.
4. Drag-branch scaling (`aD` in nondimensional form inside `pdg_rhs_casadi`)
   was derived and implemented for dimensional consistency but is
   **untested** here since `P.drag.on=false` by default; Task 11 should
   verify it when drag is switched on.

## Fix: review finding 1 (Important) — `sol.lam_defect` conversion factors

**What was wrong.** Concern #3 above ("Concerns" list) flagged that
`sol.lam_defect` is undocumented in SI terms but did not state the
conversion factors, and left them for a future task to derive. The
reviewer re-derived them and required they be written verbatim into the
`solve_pdg_colloc.m` header's OUTPUTS section, since Task 5 (certification)
is the first numeric consumer and must not have to re-derive them itself.

**What changed.** `lib/solve_pdg_colloc.m`, OUTPUTS section of the header
comment for `sol.lam_defect`, appended (comment-only, no code/numerics
touched):

```
sol.lam_defect SI conversion (all scales recoverable from sol.P
via Lc=norm(P.r0), Tc=0.5*(P.tf_lo+P.tf_hi), Mc=P.m0):
  rows 1-3 (position defects): lambda_SI = lambda_stored * Mc/Lc
  rows 4-6 (velocity defects): lambda_SI = lambda_stored * Mc*Tc/Lc  (= Mc/Vc)
  row  7   (mass defect):      lambda_SI = lambda_stored * 1
Rows 4-6 share ONE isotropic scale, so primer-DIRECTION checks
against SI thrust are valid on the stored duals unrescaled;
cross-block identities (e.g. comparing lambda_r to lambda_v, or
lambdadot_v = -lambda_r) are off by a factor of Tc without
rescaling.
```

**Covering evidence.** Comment-only change, no numerics touched, so
re-running the smoke test is sufficient evidence nothing broke. Ran:

```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_colloc_smoke"
```

Output (unchanged from the original GREEN run, as expected):
```
test_colloc_smoke PASS  tf=15.69 s  mf=26527.3 kg  fuel=3472.7 kg
```

Committed as a follow-up `booster_landing:` commit (see git log).
