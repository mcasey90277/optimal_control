# Task 6 Report: `tvlqr_design` (TVLQR gains via backward Riccati)

## Implementation

`lib/tvlqr_design.m` implements the brief's design verbatim:

- **Nominal interpolants**: `ctrl.xnom(t)` (pchip over `sol.t`/`sol.X`, 7x1) and
  `ctrl.Tnom(t)` (pchip over the interleaved node+midpoint grid `sol.t` ∪
  midpoints, values `sol.U`/`sol.Um`, 3x1). Both clamp `t` to `[0, sol.tf]`
  before interpolating.
- **Dense output grid**: `ctrl.tgrid`, `M = 4*(Nn+1)` points, `Nn = N` from
  `solve_pdg_colloc` (30 here → M=124).
- **Backward Riccati**: `ode45` integrated over `fliplr(ctrl.tgrid)` (tf→0)
  with terminal condition `vec(Qf)`, RHS `ricrhs` calling
  `pdg_dynamics(xnom(t), Tnom(t), P, nargout=3)` for analytic A,B at each
  query time, then `flipud(PV)` to reorder rows to `0→tf`.
- **Gains**: `K(t) = R^{-1} B(t)' P(t)`, stored `ctrl.K` (3x7xM) alongside
  `ctrl.Pt` (7x7xM), each `P` symmetrized (`(Pk+Pk.')/2`) before storage/use.
- Default weights exactly as specified: `Q = diag([1e-4 1e-4 1e-4 1e-2 1e-2
  1e-2 0])`, `R = 1e-10*eye(3)`, `Qf = diag([1e-2 1e-2 1e-2 1 1 1 0])` — mass
  row/column zero throughout (mass observable, not regulated).

No code deviations from the brief. Only addition beyond the brief: a
standalone check script (not committed — throwaway, in `/tmp`) to explicitly
confirm the ode45 row-order caution before trusting the implementation (see
below).

## Numerical health — no fallback needed

The brief flagged possible Riccati stiffness given `R=1e-10` against SI
thrust O(1e5-1e6 N). In practice `ode45` with `RelTol=AbsTol=1e-8` (as
specified) converged cleanly on the N=30 nominal trajectory — no stalling, no
blow-up in `P`. Fallbacks (tighter/looser tol, `ode15s`, nondimensionalize)
were **not required**; documenting per the task's instruction that they be
reported if unused.

## TDD evidence

1. Wrote `tests/test_tvlqr_riccati.m` verbatim from the brief.
2. Ran it before `tvlqr_design.m` existed:
   ```
   Unrecognized function or variable 'tvlqr_design'.
   Error in test_tvlqr_riccati (line 12)
   ```
   confirms the test fails for the right reason (undefined function), not a
   syntax/path error.
3. Wrote `lib/tvlqr_design.m` verbatim from the brief.
4. Reran: `test_tvlqr_riccati PASS` — first try, no code changes needed
   (P symmetric to machine precision every step, PSD with the stated
   tolerance including the unweighted mass direction, finite gains,
   `Pt(:,:,end)` matches `Qf` to `<1e-9`).

Command used both times:
```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_tvlqr_riccati"
```

## Explicit ode45 row-order verification (brief's own caution)

Reproduced the internal ODE call standalone and inspected `PV(1,:)` before
`flipud`:
```
max|PV(1,:) reshaped - Qf| = 0.000e+00
```
Confirms `ode45(ric, fliplr(tgrid), Qf(:), ...)` returns its first row at
`t=tf` (matching the terminal condition exactly, as it must — it's the
initial condition of the backward integration), so `flipud` correctly
re-orders rows to align with `ctrl.tgrid` (0→tf). This is exactly the
brief's requested sanity check, done outside the committed test (the
committed test checks the equivalent condition — `ctrl.Pt(:,:,end) ≈ Qf`
post-flip — directly).

## Gain magnitudes (N=30 nominal, tf=15.69 s, M=124)

| t (s)  | location | ‖P(t)‖_F     | ‖K(t)‖_F      |
|--------|----------|--------------|---------------|
| 0.000  | start    | 6.950e-02    | 2.310e+04     |
| 7.782  | mid      | 8.050e-02    | 2.762e+04     |
| 15.691 | end (tf) | 1.732e+00    | 6.529e+05     |

Gains are O(1e4)–O(1e5) N per unit state error, growing sharply near tf as
`P` grows toward the terminal weight (`Qf` diag entries O(1)–O(1e-2) vs the
mid-flight steady value ~0.07–0.08). Magnitudes are large because `R=1e-10`
barely penalizes control effort relative to SI thrust scale (O(1e5-1e6) N) —
this is a design choice in the brief's default weights, not a bug; a much
softer `R` (or nondimensionalized units) would be advisable before flying
this in `sim_closed_loop` (Task 7) if commanded corrections saturate
excessively. Flagged as a concern below, not fixed here (public interface
must stay exactly as specified).

## Files changed

- `booster_landing/lib/tvlqr_design.m` (new)
- `booster_landing/tests/test_tvlqr_riccati.m` (new)

## Self-review

- Header block present, pumpkyn-format (purpose/inputs/outputs/references).
- No `i`/`j` used as loop/index variables (`k` used throughout).
- `nargout=3` used exactly on `pdg_dynamics` calls, matching Task 2's
  certified analytic-Jacobian interface.
- `P` symmetrized before storage and before use in `K = Rinv*B'*P`, guarding
  against ODE integration asymmetry drift.
- Local functions (`clampt`, `ricrhs`) placed after the main function per
  MATLAB function-file convention; verified working (test passes).
- Public interface (`ctrl.tgrid/.K/.Pt/.xnom/.Tnom`, `opts.Q/.R/.Qf`) matches
  the brief exactly — Task 7 (`sim_closed_loop`) can consume this unchanged.
- No repo-wide regression sweep run (out of scope for this task; only the
  two new files were touched, no existing file edited).

## Concerns

1. **Gain scale**: `‖K‖` reaches ~6.5e5 near `tf` — a large control gain
   given the very small `R`. If `sim_closed_loop` (Task 7) feeds this raw
   into a tracking law without saturation, corrections near touchdown could
   dominate and possibly cause chattering against the brief's stated
   `T_cmd(t) = T*(t) - K(t)(x-x*)` saturation-by-the-sim plan. Not a defect
   in this task (interface followed exactly) — worth watching in Task 7.
2. **pchip smoothing of the thrust discontinuity**: as noted in the task
   context, `Tnom` is a pchip interpolant of a genuinely discontinuous
   min-max bang-bang profile (switch ~t=7.06s). This smooths the reference
   thrust used for both linearization (`B(t)` mildly softened near the
   switch) and tracking; accepted per the brief's explicit instruction that
   this is fine for gain design (Task 5 already established the honest
   reconstruction is a separate, HS-quadratic method used elsewhere).
3. Riccati stiffness fallbacks were not exercised (not needed on N=30); if a
   future run uses a much finer `sol.X`/`sol.U` grid, this should be
   rechecked.

## Commit

```
booster_landing: TVLQR gains via backward Riccati along guidance

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
```
