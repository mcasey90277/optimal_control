# run_phase_torus — design (approved 2026-09-15)

One entry script that builds a minimum-time costate library over a
user-specified set of departure and arrival phases, running the round
loop that was driven by hand for the 70 mN DRO -> tulip torus
(FINDINGS 59-72, `DRO_tulip/process/PHASE_TORUS_RUNBOOK.md`).

## Interface

`out = run_phase_torus(spec)` in `orbit_transfer/DRO_tulip/indirect/`.

| field | meaning | default |
|---|---|---|
| `.sD`, `.sA` | sorted phase lists in [0,1), any spacing; `.sD(1)` is the spine's departure phase | required |
| `.orbits`, `.engine` | as `run_costate_library` section 0 | 70 mN campaign |
| `.anchors` | `{name, file, sA, label; ...}`, at least one certified root at `sD(1)` | required |
| `.outDir` | campaign root; rounds in `<outDir>/round_NN`, final copy in `<outDir>/final` | required |
| `.arc` | `.nStep` 4000, `.deadlineSec` 6 h, `.span` 1.15 | |
| `.rib` | `.wallSec` 900 | |
| `.nWorkers` | rib workers | 4 |
| `.improveDays`, `.slowDays` | improve-pass and discovery thresholds (days) | 2, 2 |
| `.maxRounds` | | 6 |
| `.discover` | run the discovery step | true |
| `.plan` | print the next round's plan and stop; nothing launched | false |

## The loop (per round)

1. **Arcs**: for every anchor without both arcs on disk, spawn a batch
   MATLAB job per direction (partial saves on); wait on the files.
2. **Sheet**: `build_arrival_sheet` at the explicit `.sA` levels with every
   arc and every seed; columns whose spine changed vs the previous round
   are the ones to walk.
3. **Ribs**: `run_costate_library(launch = true)` in `round_NN` with the
   unchanged columns' ribs copied in and every earlier round's rib files
   offered to the packager; wait on `SUPERVISOR_VERDICT.txt`.
4. **Holes and improve**: `fill_holes_direct` (holes, then `.improveDays`);
   if it added points, re-package (package/audit/sweep without walking).
5. **Discover**: at every empty column and every column whose spine is
   more than `.slowDays` above a neighbour's, a direct solve seeded from the
   fastest certified root of a different family nearby (the family map
   names the families), certified at the exact phase; a root faster than
   the spine is appended to the certified-direct seed list and becomes a
   new anchor row.
6. Stop when a round adds no anchor and changes no spine (or at
   `.maxRounds`); copy the last round's catalog, receipt, sidecar, sheet
   and pictures to `<outDir>/final` with a README.

Every stage's product is on disk before the next starts; every wait has
a deadline and reads the stage's own verdict file, never stdout.

## Lattice removal

`sheet_from_arcs` (column = nearest list value mod 1), `build_arrival_sheet`
(levels = the list unwrapped by +-1), `sheet_to_catalog_file` (grid = the
lists), `run_costate_library` (non-uniform lists accepted; ribs get
explicit `targets` = the sD list as offsets from the spine),
`build_ribs`/`rib_validate` (targets through; points must land on the
list), `build_70mN_library` (arc levels = the list). A helper
`rib_targets(sD, sD0, direction)` owns the offset rule.

## Tests

- unit: `test_phase_lists` — sheet columns from a non-uniform list, an
  off-list crossing refused, catalog grid from lists, `rib_targets`,
  `rib_validate` on a list.
- plan mode on the finished 70 mN campaign: "nothing to do".
- acceptance: a 3 x 3 live torus with one anchor and small budgets.

## Out of scope

Physics, the gate stack, the catalog schema, running arcs in-process.
