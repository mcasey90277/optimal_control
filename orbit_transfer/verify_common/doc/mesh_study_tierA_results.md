# Tier A results — earth MEE order study (2026-07-27)

Run: `run_mesh_study_mee({'MEE_M2_10N','MEE_M2_2p5N','MEE_M2_1N'}, [1 2 4 8])`.
80.3 min wall. **All 12 levels converged** (`Solve_Succeeded`, defect 5.5e-15 to
9.7e-14). Artifacts under `earth_elliptic_to_geo/direct/results/mesh_study/`
(`.mat` files are git-ignored; this note is the tracked record).

## The measurement

| row | N (×1→×8) | m_f (kg) per level | deltas | p windows | verdict |
|---|---|---|---|---|---|
| 10 N | 193→1544 | 1377.101235 / 1377.109848 / 1377.117536 / **1377.874577** | +8.6e-3, +7.7e-3, **+7.6e-1** | 0.164, −6.62 | **NOT-CONVERGED** |
| 2.5 N | 708→5664 | 1370.156899 / 1371.612013 / 1372.610498 / **1375.793760** | +1.46, +1.00, **+3.18** | 0.543, −1.67 | **NOT-CONVERGED** |
| 1 N | 1740→13920 | 1370.603320 / 1371.587322 / 1371.912916 / 1372.051933 | +0.984, +0.326, +0.139 | **1.596, 1.228** | **CONVERGED** |

Switch counts: 10 N `19 20 20 18`; 2.5 N `76 75 75 75`; 1 N `173 171 171 171`.

## Finding 1 — the one clean order is p ≈ 1.2, and it supports H1

1 N is the only row that produced a usable order: deltas shrink monotonically,
the two sliding three-level windows are mutually consistent (1.596 → 1.228,
13% apart), and the switch count is stable at 171 from ×2 on. Observed
**p = 1.228**, trending downward toward 1.

This supports **H1** (p≈1: the switch-crossing interval dominates) over **H2**
(p≈2: sub-grid switch placement plus alternating-sign cancellation). It is one
row, so it is evidence, not a settled question — but it is the only row of the
three that was in a position to answer at all.

**Richardson error bar on the 1 N final mass: 1372.156 kg.** The finest level
(13,920 nodes) is still **0.10 kg** short of it, and the ×1 level — the
production resolution — is **1.55 kg** short.

## Finding 2 — 10 N and 2.5 N did not converge, for two different reasons

Both had their ×8 delta GROW rather than shrink, which is why no order exists
for them.

- **10 N changed topology**: switch count 20 → 18 at ×8, mass +0.757 kg, ΔV
  −0.011 km/s, and the rung took 179.5 s against 11.7 s at ×4 (15×, far above
  the N^1.88 scaling that predicts ~43 s). The solver worked much harder
  because it left the basin. This is a **branch finding**, not a solver
  failure, and the ×8 answer is the BETTER one on both metrics.
- **2.5 N did NOT change topology** — 75 switches at ×2, ×4 and ×8 — yet still
  jumped +3.18 kg and −0.045 km/s at ×8. A large improvement with the switch
  count unchanged.

## Finding 3 — refinement increases final mass, but the claim is WEAKER than it first looked

Every level of every row moved mass UP and ΔV DOWN, across 9 refinement steps
with no exceptions. The physical story is that a coarse mesh mislocates
switches and burns for slightly the wrong duration, wasting propellant.

**CAVEAT ADDED AFTER THE DOWN-PROJECTION TEST (see Finding 5): most of those
steps also crossed a topology change, so they are not clean measurements of
refinement at fixed structure.** 10 N went 19 -> 20 -> 20 -> 18, so EVERY one
of its steps spans a topology change. Within a FIXED 18-switch structure the
direction actually REVERSES: 1378.117 kg at 193 nodes against 1377.875 kg at
1544 nodes, i.e. refinement LOWERS mass by 0.24 kg.

The claim survives only where topology held: 1 N at x2/x4/x8 (171 switches
throughout, deltas +0.326, +0.139) and 2.5 N at x2/x4/x8 (75 switches). Those
are the rows the statement rests on; 10 N does not support it at all.

**Consequence for the campaign:** at 1 N the certified row understates
achievable mass by ~1.55 kg (Richardson). At 10 N the shortfall is real but
its cause is basin selection, not mesh resolution.

## Finding 5 — the 10 N x8 "branch change" is a BASIN effect, not a resolution effect

**The decisive test.** Take the x8 solution (18 switches, 1377.875 kg), project
it onto the x1 PRODUCTION grid (193 nodes), and re-solve there. If the
18-switch structure needs a fine mesh, it collapses back. It does not:

| on the SAME 193-node production mesh | m_f (kg) | sw | dV (km/s) |
|---|---|---|---|
| seeded from the row itself (certified) | 1377.101235 | 19 | 1.676631 |
| seeded from the x8 solution | **1378.116718** | **18** | **1.662173** |

Both `Solve_Succeeded`, both machine-tight (defect 4.8e-15 and 6.3e-15).
**The better structure exists at production resolution and is 1.015 kg better
there** (+7.4e-4 relative), with 0.014 km/s less dV.

So the x8 rung did not resolve something the coarse mesh could not represent.
It ESCAPED A BASIN, and once the better basin is known it is reachable at
coarse resolution. Framing this as "real or a solver artifact" is a false
choice: the x8 solution is a genuine converged optimum AND it is not a mesh
phenomenon.

**Two consequences, and the first is bigger than the mesh study.**

1. **The certified 10 N row is not the best solution of its own
   discretization** — it is 1.015 kg short at its own node count. This is a
   campaign finding independent of any mesh question, and it is 100x larger
   than 10 N's within-branch mesh error (0.009 kg).
2. **For 10 N, basin selection dominates mesh error by two orders of
   magnitude.** Refining the mesh is the wrong lever there; searching basins is
   the right one.

## Finding 4 — switch TIMES are not mesh-converged anywhere

The matched-switch counts FALL as the mesh refines:

| row | matched vs finest, ×1 / ×2 / ×4 |
|---|---|
| 10 N | 11 / 8 / 2 (of ~20) |
| 2.5 N | 51 / 32 / 21 (of 75) |
| 1 N | 125 / 160 / 162 (of 171) |

This looks backwards — a finer level should agree with the finest better — and
the explanation is the per-switch window. The window is 2× that switch's own
local step, so it TIGHTENS by ~2× per refinement. Matched counts falling means
switch positions are moving by more than their own local step between levels,
i.e. **switch times converge more slowly than the mesh refines**. The max
matched drift shrinks by roughly 2× per doubling (10 N: 0.486 → 0.270 → 0.109),
consistent with O(h) — H1 again, and far from H2's O(h²).

1 N is the exception that proves the mechanism: it is the row whose topology
stabilized, and it is the only one whose matched count rises with refinement.

## What this does NOT establish

- **A single p for the campaign.** One row converged. 10 N and 2.5 N need the
  plan's BRANCH CHECK (re-solve the finest level from an independently
  constructed seed) before their ×8 jumps can be attributed to discretization
  rather than to basin escape.
- **That the ×8 jumps are discretization error at all.** A better objective at
  a finer mesh is equally consistent with the solver finding a different, better
  optimum. The direction (always better) mildly favours "the coarse mesh was
  leaving mass on the table", but that is an inference, not a measurement.
- **Anything about the CR3BP / Sundman / PSR transcriptions.** That is Tier B.

## Open item found while reading the results

`mesh_order` reported p = −6.62 for 10 N. A negative observed order means the
differences GREW, which is not an order at all. `monotone = false` already
drives the verdict to NOT-CONVERGED so nothing is misreported, but the number
should be NaN rather than printed. Fix before Tier B.
