# G1 continuous-residual sweep — both certified MEE ladders (2026-08-26)

First full rung-2 coverage of the two Earth-MEE campaigns: the shared gate
`verify_common/mee_residual` (on `oc.local_residual`) run over every
certified fuel row — 7 rungs x 2 campaigns, 14/14 rows, warm re-solves via
`refresh_duals_mee` / `cr3bp_residual_gate`. Raw table:
`g1_sweep_2026-08-26.mat` (variable `TT`); columns: campaign, row, ok, N,
defect, RPmax_km, RPmed_nd, intervals>1km, top10-switch-fraction, Rt_s, wall.

| rung | 2b RPmax km | 3b RPmax km | >1 km (2b/3b) | interior med (2b / 3b, ND) | Rt max s (2b/3b) |
|---|---|---|---|---|---|
| 10 N | 39.89 | 39.89 | 57 / 57 | 2.2e-13 / 6.2e-8 | 143 / 143 |
| 5 N | 78.35 | 30.81 | 102 / 92 | 3.4e-13 / 5.7e-8 | 179 / 197 |
| 2.5 N | 45.08 | 12.89 | 172 / 153 | 5.6e-14 / 5.8e-8 | 214 / 196 |
| 1 N | 20.03 | 5.03 | 259 / 253 | 2.7e-14 / 5.6e-8 | 221 / 196 |
| 0.5 N | 25.54 | 10.56 | 546 / 556 | 1.5e-13 / 5.7e-7 | 1554 / 1083 |
| 0.2 N | 9.20 | 4.90 | 974 / 1201 | 1.3e-13 / 1.1e-6 | 2069 / 1048 |
| 0.1 N | 4.60 | 4.61 | 1787 / 1974 | 6.4e-14 / 1.8e-6 | 2070 / 2032 |

## Findings

1. **Switch attribution is total, ladder-wide.** The top-10 residual
   intervals are switch-straddling on every 2-body rung (fraction 1.00; one
   0.80 at 5 N). The 10 N single-row finding (2026-08-25) is the general
   law: the transcription's real error lives at the throttle switches.
2. **Per-switch error SHRINKS with depth** (prediction refuted — the
   2026-08-25 expectation was growth): RPmax falls 39.9 → 4.6 km as the
   deep-rung meshes densify, while the >1 km interval COUNT grows
   57 → 1787 with switch count. Severity per event improves; burden count
   conserves.
3. **The 2-body interior is machine-exact at every rung** (median
   1e-14..3e-13 ND). Between switches the certified ladder has no
   discretization story at all.
4. **The CR3BP interior floor grows with depth**: 6e-8 → 1.8e-6 ND
   (2.6 m → 75 m), tracking the time-row error growth (143 → ~2000 s).
   HYPOTHESIS (not established): t-row error misphases the lunar term
   (phi = phi0 + nM·t), a feedback absent in 2-body. A targeted check
   would re-run one deep rung with the time row replaced by quadrature.
5. Time-row error is large in absolute terms at deep rungs (~2000 s) but
   the deep transfers are months long — relative ~1e-4, consistent with
   trapezoid O(h^2) on the 1/Ldot row. Do not read the seconds alone.

Consequence for certification language: "certified" MEE rows are
position-accurate to km-level AT SWITCHES and m-level elsewhere; anyone
consuming the trajectories at switch resolution needs PSR-style refinement,
which these numbers now justify quantitatively on BOTH campaigns.
