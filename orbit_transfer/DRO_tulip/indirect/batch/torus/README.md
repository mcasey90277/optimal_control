# The scratch jobs that built the 24 x 24 torus (rounds 2-8, 2026-09-13..15)

Copied verbatim from the session scratchpad so the campaign can be replayed;
`$SCRATCH` stands for the folder they were run from (any folder works --
they only write their own `.out` files there). Every MATLAB job is launched
as `nohup /Applications/MATLAB_R2026a.app/bin/matlab -batch "run('<job>.m')"
> <job>.out 2>&1 &` (R2026a: the Parallel Computing Toolbox is licensed
there; R2025b is not, and `certify_root` dies on `gcp`). Watchers are
`zsh watch_*.sh` launched in the background; each exits on an event or an
hourly tick and prints a summary -- re-arm it after reading.

| file | what it did | FINDINGS |
|---|---|---|
| `arc_fast2_{dn,up}.m`, `arc_direct18_*.m`, `arc_direct11_*.m` | pseudo-arclength arcs from a new anchor, both arrival directions, 24-level ladder, partial saves every 50 steps | 61-68 |
| `direct_gap_probe{,2,3,4}.m` | branch-blind direct solves at gap phases warm-started from the other family's certified root; `recert_one.m` re-certifies at the exact grid phase | 61-63, 68 |
| `v3_sheet_job.m`, `v4/v5/v6_campaign_job.m` | a ROUND: rebuild the sheet from every arc + seed, compare with the previous round, reuse ribs of unchanged columns, launch the supervised rib campaign with every earlier round's ribs offered to the packager | 62-69 |
| `v6_repackage_job.m` | re-package + audit + sweep without walking (rounds 7 and 8) | 70-72 |
| `fill_holes_job.m`, `improve_job.m`, `improve_job4.m` | `fill_holes_direct`: the holes, then the improve pass (`improveDays` 2), then one stuck cell with a 900 s cap | 70-71 |
| `watch_*.sh`, `chain_*.sh` | watchers (event or hourly tick) and the chains that launched the next stage when the previous one printed its DONE line | -- |

The runbook is `DRO_tulip/process/PHASE_TORUS_RUNBOOK.md`; the methods are
in `DRO_tulip/doc/phase_torus_methods.tex`.
