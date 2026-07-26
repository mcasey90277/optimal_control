# elfo/attic — superseded ELFO experiments

**Off-path by design.** Nothing in `setup_paths` adds this folder; live code
must not depend on it. Files are kept for the record of what was tried and why
it was abandoned — several campaign notes cite them as provenance.

| file | superseded by |
|---|---|
| `gen_elfo_energy_backbone.m` | `../gen_elfo_energy_gravhom.m` (the fixed-t_f predecessor; its wall record is cited in that file's header) |
| `gen_elfo_energy_tangential.m` | same — an earlier seeding strategy |

Before reviving anything here, check it against the current
`casadi_energy_freetf` / `casadi_mintime_freetf` interfaces: these predate them.
