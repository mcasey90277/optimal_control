# Where the interpolation work stands (2026-09-22)

The one-page state of the costate-interpolation line, for whoever picks it
up next. Everything here is measured; the numbers point to FINDINGS sections.

## The question

Can the costates stored at the library's grid points be blended into a
usable guess for a transfer BETWEEN grid points? (Mike, 2026-09-20: this is
the reason the library's resolution is being refined at all.)

## The answer so far

| library | grid | usable, 60 random off-grid queries | nearest entry alone | median miss before polish | FINDINGS |
|---|---|---|---|---|---|
| 96-phase spine (one departure phase) | s_D fixed, s_A 1/96 | **52 of 60 = 87%** | 42 of 60 | 629 km | 90 |
| 24 x 48 | 1/24, 1/48 | 1 of 60 | 0 of 60 | 39,800 km | 92 |
| 24 x 24 record | 1/24, 1/24 | 0 of 60 | 0 of 60 | 42,200 km | 91 |

- Along ARRIVAL phase, 1/96 is enough: the blend beats the nearest entry
  and every converged solve stays on its corners' branch. The 8 misses sit
  at folds, where a 6 km guess still stalls (a root-finder cannot cross a
  near-singular Jacobian; the arclength polish can).
- Halving the arrival step alone (24 x 48) bought nothing. So the DEPARTURE
  axis at 1/24 is the limit. Its costates change ~32% per step (FINDINGS 84);
  what step it needs has NOT been measured.

## What exists (all tested, all pushed)

| piece | where | tests |
|---|---|---|
| the blend + the rule for which corners may be blended | `costate_common/phase_catalog_interp.m` | 30 |
| fly + junctions + polish, one fenced call | `costate_common/polish_costate_guess.m` | 5 |
| setup shared by consumers (sensitivities, closures, caps) | `indirect/catalog_blend_setup.m`, `catalog_blend_guess.m` | in the scorer's test |
| an arrival sheet read as a one-row catalog | `indirect/arrival_sheet_as_catalog.m` | 7 |
| ONE query, scaffolding exposed (transfer_study style) | `indirect/interp_study.m` | mutation-checked (V1) |
| the hit rate over random queries | `indirect/score_interpolator.m`, job `batch/score_interpolator_job.m` | 12 |

The blend rule that works is TIME CONSISTENCY: neighbours on one branch have
flight times that agree with their own sensitivities dT/ds (from the costates)
to third order in the step; another branch misses by the jump in t_f. Family
labels and a costate-difference cap are both too blunt (FINDINGS 89).

## How to run things

```
% one query, in the interactive session (a fenced solve; needs a worker pool)
interp_study                      % edit `library` and (sD, sA) in sections 1-2

% score a library (detached; ~3 h for 60 queries; resumes)
SCORE_LIBRARY=<catalog .mat, or 'record' | 'spine96'> nohup \
  /Applications/MATLAB_R2026a.app/bin/matlab -batch \
  "run('.../DRO_tulip/indirect/batch/score_interpolator_job.m')" > ~/score.out 2>&1 &
% verdict: <catalog folder>/interp_score/SCORE_VERDICT.txt
```

## Data on disk (results/ is NOT git-tracked)

- `indirect/results/library_70mN_24x24_final/` -- the record (576/576), its
  audit, `phase_transversality.mat`, `interp_score/` (0 of 60).
- `indirect/results/sheet96_resolution_test/` -- the 96-phase spine sheet,
  `spine96_catalog.mat` (the one-row view), `sheet96_residuals.mat`,
  `interp_seed_experiment/`, `interp_score/` (52 of 60).
- `indirect/results/reproduce_70mN_24x48/` -- the 24 x 48 campaign; `final/`
  holds its catalog (1,095 of 1,152 cells, audit 1,095/0); `interim_round1/`
  is the copy that was scored (1 of 60).

## Open defects

1. **The hole filler's pool dies.** `fill_holes_direct` runs its direct solves
   in-process; the worker pool idles past its 30-minute timeout and shuts
   down, and every fenced call after that throws "The parallel pool has shut
   down". 47 of the 24 x 48 filler's 58 failures are this, not the physics.
   Fix: create the pool with `IdleTimeout = Inf`, revive a dead pool before
   each cell. Test-first; then rerun the filler on the 57 holes (it resumes).
2. `reproduce_library_70mN`'s verdict says REPRODUCED = FAIL for a build on a
   finer grid than the record, because its same-library test is not the
   right question for a different grid; the comparison on shared cells (551
   of 560 agree, 9 slower, 16 missing) is the useful part. Make the verdict
   say "finer grid: compared on shared cells" instead of FAIL.

## Next steps, in order

1. **Measure the departure axis** the way the arrival axis was measured: ONE
   rib at 96 departure phases on a single arrival column
   (`rib_from_crossing`, hours), wrap it as a one-column catalog
   (`arrival_sheet_as_catalog` has the pattern; the rib needs its own
   wrapper), score it. That says whether departure needs 1/48 or 1/96.
2. Fix defect 1, refill the 24 x 48 holes, decide whether to adopt the
   24 x 48 as the record (it is a superset of the record on 551 cells; the
   9 slower cells and 16 holes are the argument against).
3. Then a full build at the measured resolution (a 48 x 96 build is ~4 days
   at 8 workers; a 96 x 96 build about 8).
4. For the fold cells: an arclength polish from the blended guess instead of
   Newton (the 8 misses on the spine).
5. Junction-state interpolation (interpolate the neighbours' multiple-shooting
   junctions rather than re-flying a guessed z8; the August rule) -- the
   guesses that fly within 400 km still need 100-500 iterations because the
   re-flight amplifies the error.
