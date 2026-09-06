# MfMax code review: what to borrow for the Huber / continuation line

Read 2026-09-06, before starting "Huber on the cells eps cannot enter". Source:
`earth_elliptic_to_geo/indirect/mfmax/mfmax-v1/src/` (Gergaud group, Fortran,
HOMPACK90 path following; built and validated here 2026-08-15, see
`MFMAX_V1_RUNBOOK.md`). This is a read of the code, not of the papers; every
claim below has a file:line.

## 0. What MfMax is, in our vocabulary

Single shooting on the PMP BVP (8 unknowns: t_f + 7 costates at L0), the
shooting function `Sfun` = terminal residual after one full-arc `rkf45`
propagation at RelTol 1e-10 / AbsTol 1e-14 (`Defs.f90:33-34`), solved by
MINPACK `hybrd` (Powell dogleg, `ftir.f90:144-197`) -- i.e. our `ss_bvp_accept`
with a different root finder. Two homotopies, both on that shooting function:

- `lambda(1)`: initial-condition homotopy from the TARGET orbit (trivial
  problem) to the true departure orbit -- `homCI`, `path.f90:223-351`.
- `lambda(2)`: energy -> fuel regularization -- `Fullpath`, `path.f90:46-218`,
  followed by HOMPACK `fixpqf` as a zero CURVE in (lambda, z), not as a
  monotone parameter.

**The regularization is our eps family, not Huber** (`user.f90:236-244`):
`coeff = lambda - beta*Tmax*p_m`, `coeffc = 2(1-lambda)`; control is 0 below
`coeff`, unit above `coeff + coeffc`, a linear ramp between. That is the
Bertrand-Epenoy ramp with lambda_2 = 1 - p. So MfMax has nothing to say about
Huber as a family; everything it has to say is about HOW to walk a
continuation parameter. That is exactly the part we do crudely.

## 1. The five borrowable ideas, ranked for the Huber walls

### 1.1 Accept intermediate points LOOSELY; converge tightly only at the end

`homCI` accepts a continuation step when hybrd reports success AND
`|S| < 1e-3`, or when `|S| < 1e-5` even if hybrd reports failure
(`path.f90:293-310`). Along the HOMPACK curve the corrector tolerance is
`arcre = arcae = 1e-6` (`Defs.f90:50-51`); only the final answer at lambda = 1
is held to `ansre = ansae = 1e-10` (`Defs.f90:53-54`), and even then `Fullpath`
finishes with one more `ssolve` refinement (`path.f90:188-191`). The
philosophy: continuation points are PREDICTORS, not deliverables.

Our ladder demands `ms_bvp` `tolR = 1e-10` AND `Hdrift < 1e-6` at EVERY rung
(`run_minfuel_race.m:45`, `P.HdriftTol`). The Huber walls stalled at
`normR = 3e-5 .. 1e-3` with `Hdrift = 8e-6 .. 9e-4` (FINDINGS 23) -- MfMax
would have ACCEPTED every one of those points and kept walking, refining at
the floor. This is the cheapest idea to test and the most likely to move a
wall: a rung gate of `normR < 1e-6` (say) with the full gate only at `p_floor`.
Risk: a loose intermediate point that is not near the true branch seeds the
next rung badly; MfMax bounds that risk with the 1e-3 / 1e-5 pair and the
final refinement. We would bound it with `Hdrift` (a physics check MfMax does
not have) at a relaxed level plus the tight endpoint gate.

### 1.2 Arclength continuation: let p turn back

HOMPACK `fixpqf` follows the zero curve of rho(lambda, z) = S parameterized by
ARC LENGTH, with a Hermite-cubic predictor and a corrector in the hyperplane
normal to the tangent (`HOMPACK90.f:54ff`, header). It requires only
`rank [dS/dlambda, dS/dz] = N` along the curve -- NOT that `dS/dz` be
nonsingular. So a FOLD in lambda (the solution branch turns back: no solution
at the next smaller lambda on this branch, one further along the curve) is
passed routinely; a monotone ladder cannot pass it at all, and bisection
toward the fold produces exactly a stall with a slowly shrinking residual.
That is a candidate mechanism for the Huber walls that is DIFFERENT from the
grazing hypothesis in FINDINGS 23, and it is testable: the ms Jacobian's
condition number should blow up approaching a fold (`dS/dz` singular at the
turning point) and not at a grazing crossing.

We have every ingredient for a pseudo-arclength `ms_bvp`: the analytic
block Jacobian in `residual()` (`ms_bvp.m:221-268`), and `dF/dp` for free --
the smoothing parameter is already a CasADi Function INPUT in
`cr3bp_minfuel_pmp` (`ps`), so `jacobian(F14, ps)` is one line. Keller's
pseudo-arclength adds p to the unknowns and one equation
`(Y - Y0)'*Ydot0 + (p - p0)*pdot0 = ds`. A day's work on the engine, reusable
by every continuation in the repo (eps, Huber, gamma, thrust).

### 1.3 Dynamic rescaling of the unknowns at every refinement

`rescale` (`path.f90:356-392`) rescales each shooting unknown to a power of
ten so its magnitude sits in [10^-scalrange, 10^(1-scalrange)), skipping
anything below `scaltol = 1e-4`, with a special factor for the longitude
costate (`scalpL`, `Defs.f90:25-30`); called after every accepted step
(`path.f90:152,186,319,337`). `Sfun` evaluates the terminal conditions on
the DESCALED state (`ftir.f90:33`), so the 1e-3 acceptance above is in scaled
units. Our costates span ~1e-4 .. 30 within one lambda_0 (`lam_r` ~ 15-33,
`lam_m` ~ 1e-4); `ms_bvp` scales nothing and hands fsolve the raw vector.
hybrd's `MODE = 2` with unit `NLEW` (`ftir.f90:79-82`) means MfMax's root
finder sees the rescaled unknowns as O(1). Cheap to add to `ms_bvp` as an
option; likely helps Newton's residual floors regardless of cause.

### 1.4 Finite-difference Jacobian vs our variational STM -- and why they never needed saltation

`JACFUN` (`user.f90:342-419`) builds BOTH `dS/dlambda` and `dS/dz` by
CENTRAL FINITE DIFFERENCES with a relative step `jac_step*|x_i|`
(1e-4, floor 1e-4). A finite difference re-integrates the whole arc with the
perturbed initial condition, so the switching TIMES move with the
perturbation and the sensitivity across a bang-bang switch comes out right
automatically -- no saltation matrix, no event detection (`rkf45` handles
the jump by step rejection). That is why MfMax converges at lambda_2 = 1
exactly (bang-bang, `|u| in {0,1}` measured) with a crude Jacobian, and it
is the mirror image of our 09-05 bug: our variational STM is exact between
switches and WRONG across them unless saltation is applied. Two lessons:
(a) an FD Jacobian is the right CROSS-CHECK for any STM we propagate through
a discontinuity -- `test_huber_saltation` is exactly this, keep it as the
standing test for every new field; (b) `dS/dlambda = 0` is hard-coded for
lambda >= 1 (`user.f90:50-52`) and `Sfun` clamps `lambda_2 <= 1`
(`ftir.f90:29`) -- the bang-bang endpoint is handled as a boundary of the
homotopy, not as a point the smooth field must reach.

### 1.5 Initial-condition homotopy: a seed generator that needs no direct solve

`lambda(1)` blends the departure state from the target orbit (where the
transfer is trivial and the costates are ~0) to the true departure orbit
(`B2fun`/`Pfun`, `user.f90:1-31,254-299`; `homCI`). Step control: double the
step after two consecutive successes at the same size, halve on failure,
give up below `par(11)`, then fall back to HOMPACK on lambda_1
(`path.f90:283-348`). Our catalog line seeds every indirect solve from a
direct collocation solution -- fine where the direct solve exists, and the
exact gap on the HIGH-GAMMA band, where the direct min-energy solve is what
fails (FINDINGS 19). An IC homotopy on the min-energy field from a nearby
converged cell (or from the arrival orbit itself) is a seed route that does
not go through IPOPT at all. This is the one idea aimed at step 1's hardest
target rather than at the Huber walls.

## 2. What NOT to borrow

- Single shooting over 20-30 days of CR3BP. MfMax gets away with it on a
  6-day, 8-revolution two-body spiral; our catalog seeds miss by
  36,000-560,000 km when single-shot (STATUS_AND_ROADMAP). Multiple shooting
  stays.
- `Pfun` ignoring `rkf45 iflag = 6` (runbook gotcha 1): a truncated
  integration is used as if complete. Our propagators throw; keep that. (The
  review P2 item "check ode113 reached dt" is the same lesson from our side.)
- The silent clamps in `Phifun` (P >= 1e-3, dt = 0 on non-increasing L).
- The hard-coded target-orbit costate zeros as the lambda_1 = 0 start are
  specific to "start at the target"; for us the analogous trivial problem
  is "depart from the arrival orbit", which needs the arrival state to be a
  legal departure -- true for periodic orbits, so it transfers.

## 2b. Outcome of steps 1-2 (same day, FINDINGS 24)

Both walls diagnosed with `cond(J)` (idea 1.2's test) and re-walked with
MfMax's loose acceptance (idea 1.1): **neither fold nor gate**. cond(J) is
flat across every failed iterate near the branch; the loose gate bought one
rung on (6,8) that would not tighten (net shallower) and never fired on
(1,2). `huber_switch_diag` then found the actual mechanism: **grazing
bifurcations** of the Q = 1 switch structure -- a near-tangent crossing
(|dQ/dt| = 0.045 vs 0.10-0.80 on clean cells) on (6,8), a Q maximum of
0.9907 about to cross on (1,2). The solution curve has a corner there; no
parametrization passes it. Ideas 1.1 and 1.2 are closed for the Huber
walls. Idea 1.3 (rescaling) is moot for a corner. Idea 1.5 (IC homotopy)
stands, for the high-gamma seeds. Idea 1.4 (FD cross-check) stands as the
standing test. The cure is a FAMILY change at the corner: eps handoff, a
Huber-eps hybrid ramp, or stepping over the bifurcation.

## 3. Recommended order for step 1 ("Huber where eps cannot enter") -- as written before 2b

1. **Diagnose the two Huber walls first, cheaply (one afternoon):** re-run
   the (6,8)@1.2 and (1,2)@1.2 Huber walks with `ms_bvp` returning
   `cond(J)` at each accepted rung and at the stalled iterate. Blow-up ->
   fold (idea 1.2 is the cure); flat -> grazing/conditioning of the field
   (FINDINGS 23 hypothesis; idea 1.3 and a saltation-denominator log are
   the cure). This decides which engine change to make before making it.
2. **Loose-rung / tight-floor gate (idea 1.1)** as a `run_minfuel_race`
   option, tested on the same two walls. Cheapest possible test; if it
   passes the walls, the Huber story on the grid changes from 4/7 to 6/7 or
   7/7 before any engine work.
3. **Rescaling option in `ms_bvp` (idea 1.3)**, tested on the same walls
   and on the eps arms' bisection counts (eps fails 8-10 times per cell --
   rescaling may cut that too).
4. **IC homotopy for the high-gamma band (idea 1.5)** -- the seed route for
   the cells the direct solve cannot reach; run Huber AND eps from it.
5. **Pseudo-arclength `ms_bvp` (idea 1.2)** only if step 1 shows folds.
   Largest build, largest payoff if the walls are folds; wasted if they
   are not.

Each of 2-5 reuses `run_minfuel_race`/`ms_bvp` and the FINDINGS 23 records
as the baseline; none changes the shipped catalog.

## 4. Pointers

- `MFMAX_V1_RUNBOOK.md` -- build, `in.dat` decode, reference 10 N run.
- `mfmax_docs/mfmax.pdf`, `MfMaxmethod.pdf` -- the method papers.
- Watson, Sosonkina, Melville, Morgan, Walker, "Algorithm 777: HOMPACK90,"
  ACM TOMS 23(4), 1997 (`fixpqf` = quasi-Newton augmented-Jacobian curve
  tracker).
- Keller, "Numerical solution of bifurcation and nonlinear eigenvalue
  problems," 1977 (pseudo-arclength continuation).
- `DRO_tulip/FINDINGS.md` 22-23; `doc/extremal_and_local_min_survey.md`.
