# earth_elliptic_to_geo/indirect — indirect (PMP shooting) side

## What is here

| where | what |
|---|---|
| `mfmax/` | **MfMax v0/v1** — the Gergaud group's own Fortran indirect solver for exactly this problem (single shooting + HOMPACK differential homotopy; ENSEEIHT-IRIT, 2004). Ported to gfortran and validated on this machine 2026-07-20: v0 at 10 N / c_tf=1.5 converges in ~1 s to m_f = 1378.37 kg vs our direct 1377.10 kg — an independent indirect cross-check of the direct campaign. v0 fixes both t_f and L_f; v1 (doc in `mfmax/mfmax_docs/`) frees t_f via the s=t/t_f + dtf/ds=0 device with L_f fixed. Build recipe: gfortran with `-std=legacy -fallow-argument-mismatch` (+ `SDKROOT` for the link); run `path` with an `in.dat`. **v1 built + validated 2026-08-15** (10 N, `par(1)=2.0`: m_f = 1381.41 kg, free t_f = 146.39 h, ‖F‖ = 2.2e-9, GEO hit to 8.9e-10 Mm) — full build/run/input/output reference in **`mfmax/MFMAX_V1_RUNBOOK.md`**, with a step-by-step how-to and slide deck at **`doc/mfmax_v1_howto.tex`** / **`doc/mfmax_v1_slides.tex`**. |

## Not yet built

Our own MATLAB indirect solver for this problem (PMP shooting on the MEE/
L-domain formulation, seeded from the certified direct solutions — the same
direct-seeded strategy as `../../GTO_tulip/indirect/ifs/`). Blocked-adjacent:
the raw-dual/primer anomaly (`../process/DESIGN_dual_map.md`) means direct-KKT
costate seeds are currently unreliable at high eccentricity; MfMax's converged
costates (written to its `out.dat`) are a candidate clean seed source.
