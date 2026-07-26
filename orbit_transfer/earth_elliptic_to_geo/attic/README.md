# earth_elliptic_to_geo/attic — superseded campaign scripts

**Off-path by design.** `direct/setup_paths.m` does not add this folder (it
once listed an `attic` entry, but that resolved to `direct/attic`, which has
never existed — the dead entry was removed 2026-07-26).

One-off drivers and watchdogs from the Task-9/Task-11 era of the MEE thrust
ladder, kept as a record of how the deep rungs were actually reached. The
production paths that replaced them are `drivers/`, `reproduce/`, and
`frontdoor/run_gergaud.m`.

`run_task9_rung.m` in particular is superseded by `reproduce/reproduce_row.m`
(see that file's own DEPRECATED header).
