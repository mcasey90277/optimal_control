# Phase 3 results — continuous PMP residuals (2026-07-28)

> **RETRACTION 2026-07-28 (post-review).** The "0.233 kg production-mesh error"
> reported below is **WITHDRAWN as a measurement** and downgraded to a
> conditional model estimate. Both external reviewers returned NOT DEFENSIBLE
> on the same grounds, and they are right:
>
> - `E = k·Σh³` factors a single constant out of `Σ Cᵢhᵢ³`, i.e. assumes the
>   local coefficient is phase-independent — which **this study's own data
>   falsifies**, since redistribution made the worst cell 8x worse. The 0.233 kg
>   is then an artifact of the two particular grids compared.
> - `19 → 19` switches does **not** establish same-branch. Our own 2.5 N row held
>   75 switches ACROSS a real branch change.
> - The `R_x/h³` constancy is a **state-residual** scaling and does not
>   establish **objective-error** scaling; the only bridge between them is the
>   first-order estimator, which failed by 5.3x.
> - Gemini adds a stronger reading of that failure: with bang-bang control on an
>   unaligned grid the objective is **non-smooth** in node placement, so a smooth
>   adjoint sensitivity is not miscalibrated but *theoretically invalid* here.
>
> **What survives:** two same-node-count 10 N solutions differ in mass by
> **0.1456 kg**. Since mass rose toward the limit, that is a defensible **LOWER
> BOUND** on the production mesh's mass error at 10 N — not a point estimate.
> The "basin dominates by 4.4x" claim is withdrawn with the denominator.
> Full review: `review_2026-07-28_phase3_results.md`.
>
> **BOTTOM LINE (as originally written).** The per-interval residual of the earth MEE solutions is the
> **local truncation error of the trapezoidal rule** — `median R_x / h^3` is
> 0.171 / 0.171 / 0.163 across the 10 N, 2.5 N and 1 N rows, constant to 5%
> while revolutions, switches and node count all vary 9x. It is **not** driven
> by the bang-bang structure, and it peaks at **apogee** (16–19x the perigee
> decile) because the mesh is uniform in true longitude and apogee is where the
> time-steps are longest. Judged by error equidistribution rather than angular
> resolution, the mesh over-resolves perigee and under-resolves apogee.

Code: `verify_common/pmp/` (`pmp_residual`, `pmp_residual_mee`, `pmp_plot`).
Figures + data: `earth_elliptic_to_geo/direct/results/pmp/` (`.mat` git-ignored).

## Why this replaced the mesh ladder

Mesh convergence is a **comparative** measurement: two or more solutions, and
therefore a proof they are the same branch. That proof is what defeated Phase 2
— every ladder crossed a topology change, and both external reviewers returned
NO on the criteria for getting past it.

A PMP residual is computed from **one** solution. Nothing is compared, so the
branch problem does not arise. It also closes the gap the review called
blocking: `Solve_Succeeded` plus a nodal defect verifies only that IPOPT solved
the transcription, and nothing checked that the result approximates a
continuous extremal.

## P3.1 — toy validation (GATES EVERYTHING)

Double integrator, `min ∫u dt`, whose extremal is closed-form:
`ts = 2−√2`, `c1 = −1/√2`, `H = 1−√2`.

| check | result |
|---|---|
| self-gate, both arcs | 1.7e-17 / 8.3e-17 PASS |
| primal residual `R_x` | max 1.5e-16 |
| costate residual `R_λ` | max 1.1e-15 |
| Hamiltonian | −0.414213562373, matching analytic to 12 digits, CoV 1.5e-15 |
| unaligned grid: spike | **1.71e-03** at the straddling interval |
| unaligned grid: elsewhere | 3.3e-16 |

Both halves were required. (a) shows the checker is **correct**; (b) shows it
is **sensitive** — 13 orders of separation, peak at the right index. A checker
passing (a) but missing (b) would be blind to the only error that matters.

**Lesson found here, not assumed:** at a bang-bang switch the control is
genuinely double-valued, so a single shared node cannot carry both `u=1` and
`u=0`. A true retained breakpoint is a **duplicated** node — which is what a
multi-phase formulation provides. The first test version shared a node and the
self-gate refused the run rather than reporting a spurious 1.77e-02.

## P3.2 / P3.3 — the three earth rows

All self-gates pass against the solver's own reported defect:

| row | nodes | arcs | self-gate | solver reported | `λ_t` CoV |
|---|---|---|---|---|---|
| 10 N | 194 | 20 | 6.245e-15 | 6.273e-15 | 2.25e-08 |
| 2.5 N | 709 | 76 | — | — | — |
| 1 N | 1741 | 171 | 8.266e-14 | 8.266e-14 | 2.26e-08 |

`λ_t` is the conserved quantity for this transcription (`H_σ` is **not**
constant — the Gauss equations carry σ explicitly through `L = π + σ·ΔL`).
Constant to 8 digits in every row: the costate field is a genuine extremal's.

### The headline gap

| | median |
|---|---|
| discrete defect (what was solved) | 5.0e-16 |
| continuous residual (what we want) | 3.7e-04 |
| **ratio** | **~7e11** |

The solution satisfies the equations it was *solved against* about 7×10¹¹ times
better than the continuous equations we care about. That is discretization
error made directly visible — no ladder, no branch argument, no extrapolation.

### Finding 1 — the costate residual was a FALSE ALARM (a correction)

`R_λ = 13.9` was flagged as possibly alarming. It is an **absolute** norm, and
costates are sensitivities with |λ| of 24–264. Normalized:

| row | median relative `R_λ` | max |
|---|---|---|
| 10 N | 4.6e-06 | 2.2e-01 |
| 2.5 N | 7.7e-06 | 1.8e-02 |
| 1 N | 1.4e-06 | 7.4e-03 |

One interval anywhere reaches 22%; the rest are parts-per-million.

### Finding 2 — the peak is at APOGEE, not perigee

The stated guess was perigee ("where the dynamics are stiff"). Wrong.
Median `R_x`, apogee-most decile vs perigee-most decile:

| row | ratio |
|---|---|
| 10 N | **16.5x** |
| 2.5 N | **19.2x** |
| 1 N | **18.2x** |

**Mechanism.** The mesh is uniform in true LONGITUDE. Near apogee the angular
rate is slow, so a fixed ΔL spans a long Δt — the biggest time-steps sit at
apogee, and that is where a trapezoid hurts most. This is the same fact as the
positive slope in the residual-vs-step panel.

### Finding 3 — it tracks NEITHER switches NOR revolutions

| row | revs | switches | sw/rev | nodes/rev | median `R_x` | **median `R_x`/h³** |
|---|---|---|---|---|---|---|
| 10 N | 7.33 | 19 | 2.59 | 26.3 | 3.73e-04 | **0.171** |
| 2.5 N | 27.82 | 76 | 2.73 | 25.4 | 4.61e-04 | **0.171** |
| 1 N | 69.18 | 171 | 2.47 | 25.2 | 4.65e-04 | **0.163** |

Median `R_x` is flat while revolutions vary 9.4x, switches 9x and node count
9x. Normalized by `h^3` it is **constant to 5%**.

That is exactly the trapezoidal rule's local truncation error, `O(h^3)` per
step, with a coefficient set by the **orbital dynamics** — not by thrust level,
not by control structure. The campaign fixes nodes-per-revolution (26.3 / 25.4
/ 25.2), so the same steps recur at the same orbital phases in every row and
the same residual follows.

### Finding 3 re-tested with SUMS, after review (Gemini's objection)

Gemini objected that a MEDIAN dilutes sparse switch spikes into a large
population of smooth cells. Re-tested with sums, which is what an integrated
error actually cares about:

| row | switch cells | share of total ΣR_x | population share | mean sw/int |
|---|---|---|---|---|
| 10 N | 19 | **16.2%** | 9.8% | 1.77 |
| 2.5 N | 76 | **15.3%** | 10.7% | 1.50 |
| 1 N | 171 | **10.5%** | 9.8% | 1.08 |

The objection is fair but the conclusion **survives**: switch cells carry
10–16% of the total residual against a 10–11% population share. They are
mildly over-represented, not dominant, and the over-representation SHRINKS as
switches multiply. Finding 3 stands on sums as well as medians.

### Consequence: the earlier prediction is properly closed out

Phase 3 pre-registered: *"Expected: residuals concentrated at switch-crossing
intervals. If instead uniform, the bang-bang structure is not the dominant
error source and H1's mechanism is wrong."*

Measured at 10 N: switch cells median 1.25e-03 vs interior 3.48e-04 — ratio
**3.6** — and the single largest residual is an interior cell. The prediction
**failed**, and Findings 2 and 3 say what is happening instead.

This points **away** from switch-aligned refinement (the PSR claim, Phase 2's
Step 6) as the lever, and **toward** node redistribution by orbital phase.

## What this does NOT yet deliver

**A mass error bar.** The residual is a *local* per-interval quantity. It says
the mesh is mis-allocated; it does not say what that costs in kilograms.
Converting it needs a sensitivity argument — see next steps.

**Anything about basins.** PMP conditions are NECESSARY; every local extremal
satisfies them. The 10 N solution that is 1.015 kg worse than its neighbour
passes this check cleanly. The basin finding stands separate.

## Step 1 result — the residual converts to KILOGRAMS

`pmp_objective_error.m`. The discrete solution satisfies its own defect
constraints exactly; measured against the exact flow each interval carries a
signed residual `d_k`, so the discrete solution solves a problem whose defect
constraints are displaced by `d_k`. The multiplier of defect k IS
`dJ*/d(defect k)`, so to first order `ΔJ ≈ Σ_k lamDef(:,k)ᵀ d_k`. No extra
adjoint propagation — the NLP already computed the sensitivities, which also
avoids the backward instability that forced per-interval propagation.

| row | m_f (kg) | implied \|Δm\| (kg) | no-cancel bound | cancellation | relative |
|---|---|---|---|---|---|
| 10 N | 1377.101 | **1.237** | 3.517 | 2.8 | 9.0e-04 |
| 2.5 N | 1370.296 | **5.714** | 8.264 | 1.4 | 4.2e-03 |
| 1 N | 1370.798 | **4.046** | 6.634 | 1.6 | 3.0e-03 |

### Finding 1.1 — the magnitudes are comparable to the basin gains

Independently measured basin gains on the same rows were **+1.015**, **+5.971**
and **+1.454** kg. Against implied mesh errors of 1.237, 5.714 and 4.046 kg,
the two effects are the SAME ORDER — and at 2.5 N they agree to within 5%.

**This is an observation, not a claimed relationship.** The two quantities
measure different things (distance to a better local optimum on the same mesh,
versus distance to the continuous solution). Their similarity may reflect a
common underlying scale or may be coincidence. It is recorded because it is
striking and because it bears directly on the earlier, now-withdrawn claim that
basin error dominates discretization error: on these numbers they are
comparable, not separated by orders of magnitude.

### Finding 1.2 — cancellation is WEAK, which tests H2 directly

The cancellation ratio (no-cancel bound over actual) is **2.8 / 1.4 / 1.6**.
H2 argued that alternating-sign errors across many switches would partially
cancel in an integrated quantity, giving the objective a better effective rate
than the trajectory. Measured, the per-interval contributions largely
REINFORCE: cancellation removes only a factor of 1.4–2.8, and it is WEAKEST on
the rows with the most switches (2.5 N and 1 N), which is the opposite of what
the argument predicts. H2's cancellation mechanism does not operate here.

### THREE CAVEATS, none of them small

1. **The estimate overpredicts against the ladder.** Tier A observed the 1 N
   mass move 1.45 kg across the whole 1x→8x ladder with 0.10 kg of Richardson
   remainder — about 1.55 kg total. The first-order estimate says 4.05 kg,
   roughly 2.6x larger. Either higher-order terms reduce it, or the ladder
   never converged far enough to show the full error, or the sensitivity
   argument has a flaw. **Unresolved**, and the mesh-redistribution experiment
   is what discriminates.
2. **The SIGN is not yet established.** All three come out negative under the
   current multiplier convention, which would mean the discrete solution
   OVERSTATES achievable mass — the opposite direction from the basin effect,
   so the two would partly offset. That is a materially different story from
   the one in Section "basin", and it must not be asserted until the sign is
   fixed empirically.
3. **At 2.5 N and 1 N the residuals describe an IMPROVED point.** The dual
   refresh reported `m_f` rising (0.913194 → 0.913531 at 2.5 N; 0.913572 →
   0.913865 at 1 N), so the analysis characterizes that better extremal, not
   the published row. This is the seed-sensitivity effect appearing again.

## Step 2 result — mesh redistribution: the closed loop

Predictions were computed and printed BEFORE the re-solve. Earth 10 N,
**same node count** (193 intervals), nodes redistributed to be uniform in
physical TIME — which is what minimizes `Σh³` at fixed `N` and `T`, and is
therefore the optimal redistribution given the measured `E ∝ h³`.

| # | prediction | measured | verdict |
|---|---|---|---|
| 1 | `Σh³` drops **2.68x** | **2.66x** | **HIT** |
| 2 | residual falls by ~the same factor | median 1.47x, **max 0.13x (WORSE)** | **MISS** |
| 3 | implied \|Δm\| falls 1.237 → 0.461 kg | — | — |
| 4 | measured Δm ≈ **0.776 kg** | **+0.146 kg** | **MISS, 5.3x** |
| 5 | switch count stays 19 | **19 → 19** | **HIT** |

`h` range collapsed from 36.2x to 1.6x, so the redistribution did exactly what
it was designed to do. Switch count held, so the mass comparison is on ONE
BRANCH and is not confounded — the objection that defeated Phase 2 does not
apply here.

### The SIGN is now resolved, and it was inverted

Step 1's estimate came out negative, which would have meant the discrete
solution OVERSTATES achievable mass, so a better mesh should give a LOWER
`m_f`. **`m_f` rose**, 1377.101235 → 1377.246863. The multiplier-convention
sign in the first-order estimate was inverted: the discrete solution
UNDERSTATES achievable mass, and improving the mesh raises it. That direction
agrees with what Tier A saw, and it is now established by perturbation rather
than by convention-chasing.

### The first-order estimate OVERPREDICTS by ~5x

Predicted 0.776 kg, measured 0.146 kg. This is the same defect flagged in
Step 1's caveat (a), where the 1 N estimate ran 2.6x above the ladder's total
observed movement. It is now measured directly on a controlled, single-branch
perturbation: **the sensitivity estimate is good for order of magnitude and
not usable as a calibrated error bar.**

### The deliverable — a measured error bar for the production mesh

Two solutions on the same branch, same node count, differing only in node
placement, with `Σh³` of 2.668 and 1.003 and masses differing by 0.1456 kg.
Taking `E = k·Σh³`, which the measured constancy of `R_x/h³` supports:

| | value |
|---|---|
| `k` | 0.0874 kg per unit `Σh³` |
| **error of the production mesh** | **0.233 kg** (1.7e-04 relative) |
| error of the redistributed mesh | 0.088 kg |
| first-order estimate, for comparison | 1.237 kg (overpredicts 5.3x) |
| basin gain on this row | 1.015 kg — **4.4x the mesh error** |

**This is the number the study set out to obtain**, and unlike every earlier
attempt it rests on a same-branch comparison with the branch verified, not on
an extrapolation across a topology change. Caveats: two points and one
assumption (`E ∝ Σh³`), though that assumption is supported by the independent
`R_x/h³ = 0.171/0.171/0.163` measurement rather than fitted here.

### Prediction 2's miss is a real finding

The residual median improved only 1.47x, and the **maximum got 8x WORSE**
(1.05e-02 → 8.38e-02). Making `h` uniform in time removes the apogee penalty
but evidently creates a worse worst-cell somewhere else — most likely at
perigee, where a now-larger `h` meets rapidly varying dynamics, i.e. the
coefficient `C` in `E ≈ C·h³` is NOT phase-independent after all. Uniform-in-time
is therefore not the optimal mesh; equidistributing `C·h³` would need `C`
measured per phase. The mass still improved because the total fell even as the
worst cell rose.

## Experiments still to run

1. **Objective sensitivity: ΔJ ≈ Σ λᵀ·R.** The costate IS the sensitivity of
   the objective to a state perturbation, and we already have both the costates
   and the residual map. This converts the residual into kilograms and is the
   missing link to the original question.
2. **Mesh redistribution test.** Re-solve 10 N at the SAME node count with
   nodes redistributed to equidistribute `h^3`, and measure whether the
   residual drops and `m_f` moves by the amount step 1 predicts. Together, 1
   and 2 form a closed falsifiable loop.
3. **Tulip and ELFO adapters** — both need their own, since the tulip builds
   dynamics inline (no reusable RHS, no stored duals) and ELFO uses a 9-state
   two-primary solver.
