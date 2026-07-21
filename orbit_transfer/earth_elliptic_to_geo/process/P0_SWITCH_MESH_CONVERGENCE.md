# P0 — Deep-rung switch count: mesh-convergence certification (0.2 N)

**Date:** 2026-07-21  **Status:** DONE (0.2 N); 0.1 N open.
**Scope:** certifies whether the deep-rung bang-bang *switch count* is mesh-converged.

## The question

An external code review (GPT-5.6-terra + Gemini 3.1 Pro + host, 2026-07-20, on the
deep-rung recipe of `process/DEEP_THRUST_LESSONS.md`) converged, independently and
forcefully, on ONE high-value concern: the deep rungs run at **`nodesPerRev = 8`**
(`../direct/drivers/reproduce_deep_rung.m`), so 0.2 N packs ~823 switches into ~2773 nodes —
**~3 nodes per switch**. A machine-tight collocation defect (~2e-13) certifies the
*discrete transcription equations*, **not** that the continuous switch structure is
resolved. Reviewers' charge: the reported switch count could be a discretization
artifact, and defect-based certification of a bang-bang solution is not enough. The
campaign's own design notes had argued 15-40 nodes/rev with a convergence study; that
study had not been done for the deep rungs. This note does it.

## Method (two-pronged; the primal half is load-bearing)

**(A) PMP switching-function cross-check — INCONCLUSIVE here.** Both reviewers
suggested reconstructing the continuous PMP switching function from the node costates
(`out.lamDef`) to count switches independent of mesh density. We ran it
(`../direct/verify/verify_pmp_mee.m` on the banked 0.2 N solution): it is **blocked by the
campaign's open raw-dual/primer anomaly** (`process/DESIGN_dual_map.md`). At the 0.2 N
eccentricity the reconstructed costates are badly corrupted — primer **42.98 deg**
misaligned (gate <1 deg FAIL), S-sign agreement **71.9 %** (gate >=99 % FAIL) — so its
220 S=0 crossings cannot be trusted against the nodal count. The costate route is a
dead end until Campaign B closes.

**(B) Primal mesh-refinement — the real test.** Purely primal (throttle + node
longitudes only, `../direct/verify/switch_structure.m`), so it is **immune to the dual anomaly**.
Warm-refine the certified 8/rev solution through 16 -> 24 -> 40 nodes/rev (chained
`interp_warmstart` + a short warm eps-tail `[0.01 0.003 0.001 0]` via `homotopy_mee`,
deep-rung levers on), re-solving to eps=0 at each, and compare the switch structure.
Driver: `../direct/verify/meshstudy_switch.m`. All four densities reached eps=0 with defect
~2-3e-13; single run, zero crashes (~3.6 h wall).

## Result (0.2 N, c_tf=1.5, ~346.7 rev; all rows eps=0-certified)

| nodes/rev | nodes  | revs   | switches | sw/rev | m_f [kg]  | maxDefect |
|-----------|--------|--------|----------|--------|-----------|-----------|
| 8 (orig)  | 2 774  | 346.73 | **823**  | 2.374  | 1377.287  | 2.5e-13   |
| 16        | 5 549  | 346.70 | 865      | 2.495  | 1375.918  | 2.3e-13   |
| 24        | 8 322  | 346.67 | 871      | 2.512  | 1375.836  | 2.4e-13   |
| 40        | 13 868 | 346.68 | 863      | 2.489  | 1375.819  | 3.4e-13   |

Figure: `../direct/results/p0_switch_mesh_convergence.png` (`../direct/viz/fig_switch_convergence.m`).

## Findings

1. **Mass is mesh-converged: ~1375.8 kg.** Corrections -1.37, -0.08, -0.02 kg ->
   limit ~1375.81 kg. The original 8/rev value (1377.29) was **+1.5 kg / +0.11 %
   high** — a real under-resolution bias, not noise.
2. **Revolution count is mesh-invariant: 346.7 rev.** DeltaL is fully converged, so the
   deep-rung *geometry* claims (~345 rev at 0.2 N, ~690 at 0.1 N) are solid. Only the
   switch count was ever in question.
3. **Switch count converges to a BAND ~863-871 (~866 +/- 5, 2.49-2.51 sw/rev), not a
   single integer** — and it does NOT diverge (24->40 goes *down*, staying in-band). The
   original **"823" is a genuine ~5 % undercount** from coarse mesh. The residual
   +/-5-switch wobble across 16/24/40 is the expected non-invariance of exact switch
   *placement* on a uniform mesh, not a solver defect.

## Verdict

The reviewers were **right on the narrow point, wrong on the alarming one.** Right:
8/rev under-resolves the switch count, so quoting an exact integer is not defensible.
Wrong (refuted): this is NOT the direct method "fooling itself" — the physics (mass,
revs) is trustworthy and the switch structure converges to a bounded band rather than
climbing without limit. This **validates the campaign's own instinct** (README/CAMPAIGN:
"switch count should be read as a band, not a fixed integer") and now supplies the
number.

**Recommended reporting for the deep rungs.** State the converged quantities and flag
the coarse-mesh integer as a lower bound. For 0.2 N:

> m_f = 1375.8 kg (mesh-converged), 346.7 rev, **~866 +/- 5 switches (2.50 sw/rev)**.
> The 8-node/rev point estimate of 823 switches is a ~5 % undercount and should not be
> quoted as exact.

The mass ladder headline numbers are unaffected in substance (0.1 % shifts at most).

## Reproduce

```matlab
cd earth_elliptic_to_geo/direct; setup_paths
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'))
rows = meshstudy_switch('../direct/results/MEE_M2_0p2N.mat', [16 24 40]);  % ~3.6 h, resume-safe
fig_switch_convergence(rows)                                     % or fig_switch_convergence() for the baked result
```

## Open

- **0.1 N** not yet run. Same protocol confirms the pattern; cap at 24/rev (40/rev is
  ~28k nodes / ~330k vars at ~690 rev — very heavy). Expect the same conclusion: mass
  converged, switch count a band, 8/rev an undercount.
- **PMP switching-function verification** stays blocked by the raw-dual anomaly
  (`process/DESIGN_dual_map.md`, Campaign B). Until that closes, the primal
  mesh-refinement here is the only trustworthy switch-structure certification at deep
  rungs.
