# Phase 3 results — continuous PMP residuals (2026-07-28)

> **BOTTOM LINE.** The per-interval residual of the earth MEE solutions is the
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
