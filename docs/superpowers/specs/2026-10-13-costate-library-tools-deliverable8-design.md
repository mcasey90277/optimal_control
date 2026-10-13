# Design: costate-library tools and deliverable 8 (70 mN DRO -> tulip, with the status layer)

**Date:** 2026-10-13
**Status:** approved in design review with Mike (sections 1-4); spec awaiting review
**Library:** `orbit_transfer/DRO_tulip/indirect/results/library_70mN_24x48_v2/costate_catalog_dro_tulip_70mN.mat`
(the library of record since 2026-10-05, FINDINGS 98): 1,152 primaries on a
24 x 48 phase grid, 656 alternatives, optimality status and junction states on
every entry, fail-closed audit CLEAN (1,808 ok / 0 bad).
**Predecessors:** deliverables 2-7 under `*/deliverables/` (the file set, README
shape and self-test rule they established; the library-deliverable skill).

## 1. Purpose

Put the tools for the current costate library in good shape and ship the
library with them as **deliverable 8**, in the style of deliverables 2-7.

**Recipient: Darin.** He picks a transfer by phases, gets a converged z8 that
`pumpkyn.cr3bp.tfMin` accepts unchanged, and flies it. Deliverable 8 adds
what the status layer makes possible: every transfer the builds found at a
phase pair, not only the fastest, each with what was established about its
optimality; the junction states that fly a 17-38 day arc accurately; and a
map of where the non-optimal (status-2) transfers sit.

**Approach A (decided):** the tools live in the repo (`DRO_tulip/indirect`
and `costate_common`), are built test-first and reviewed, and the deliverable
is assembled from them by the builder. One copy of each tool; our own studies
use the code Darin gets.

### Success criteria

1. The deliverable folder runs standalone: a fresh MATLAB with only pumpkyn +
   pumpkynPie and the folder on the path runs both example scripts and every
   helper's self-demo; a dependency scan finds no call outside the folder.
2. The exported catalog's two content keys equal the audited library's, so
   the README's verification numbers apply to the shipped file.
3. Darin's existing calls (`costate_catalog_pick(cat, tauDep, Np, depDays,
   arrDays, thrustN)`, `costate_lib_describe(cat)`,
   `costate_catalog_extremes(cat, filt, doPlot, metric)`) behave as before.
4. A primary flown from its junctions arrives within 1 km / 0.1 m/s of the
   catalog's own endpoint rule; `tfmin_accepts` reports |dz| < 1e-6 on it.
5. Every README number is computed from the exported catalog.

## 2. Decisions (from the design review)

| decision | choice |
|---|---|
| audience | Darin (pumpkyn tfMin); the mission-design use is served by the same tools |
| content | all of v2: primaries, alternatives, statuses, junctions |
| tool set | lookup (pick + all-at-phase), fly (endpoints + thrust program + tfMin check), survey (status-aware describe + extremes, status torus map), movies + two examples; **interpolation left out** until the junction-interpolation experiment is done |
| catalog file | a cleaned export: all data kept, machine paths and internal bookkeeping replaced by provenance text; same content keys |
| organisation | A: tools in the repo, packaged by the builder |

## 3. The deliverable folder and the catalog export

**Folder:** `DRO_tulip/deliverables/costate_catalog_70mN_status_deliverable8/`,
zipped from its parent so the archive unpacks into that name.

**`export_deliverable_catalog(recordMat, outMat)`** (new, `DRO_tulip/indirect`)
writes `costate_catalog_dro_tulip_70mN.mat` (variable
`costate_catalog_dro_tulip_70mN`):

- **Kept as is:** every sheet field (primaries, `status`, `status_reason`,
  `junctions`, conjugate verdicts, gate values, margins, `family_index`,
  `entry_notes`), the full `alternatives` table, `status_key`, `constants`,
  `thruster`, `rungs_N`, orbit recipes (`dep_family`/`dep_params`,
  `arr_family`/`arr_params`), `derive`, `conj_test`, `hyp_gates`,
  `second_order`, `notes_key`, `families`, `schema`, `usage`.
- **Replaced:** `merge` and `status_layer` collapse into one `provenance`
  struct (text + counts: how the library was built, refilled, merged,
  backfilled, re-certified, audited; no paths, no harvest key).
  Each alternative's `source` string loses file paths and keeps the kind
  (`sheet candidate, re-certified`, `rib refusal (stop point)`, `displaced
  primary (merge)`, ...), by a fixed mapping with a test.
- **Guaranteed, asserted at export:** `catalog_content_key` and
  `alternatives_content_key` equal the record's; `catalog_schema('validate')`
  passes; the file loads with `load` alone.

**Contents:**

| file | role |
|---|---|
| `costate_catalog_dro_tulip_70mN.mat` | the catalog (cleaned export) |
| `costate_catalog_pick.m` | five-coordinate lookup, honesty contract, unchanged |
| `transfers_at.m`, `status_key.m` | every transfer at a phase pair, ranked, with status |
| `catalog_endpoints.m`, `fly_catalog_entry.m`, `tfmin_accepts.m` | the fly layer (section 4) |
| `entry_thrust_program.m`, `fly_transfer.m` | the two flight engines the fly layer wraps |
| `costate_lib_describe.m`, `costate_catalog_extremes.m`, `plot_status_torus.m` | the survey layer (section 5) |
| `costate_catalog_extremes_movies.m`, sample `.mp4` + `.gif` | movies (section 6) |
| `print_transfer_summary.m`, `plot_transfer_3d.m` | the summary and the 3D figure the examples use |
| `costate_catalog_example.m`, `status_layer_example.m` | worked examples (section 6) |
| `get_family_orbit.m`, `phase_state.m`, `periodic_pp.m`, `validate_flight.m`, `scalar_verdict.m`, `catalog_schema.m` | dependencies |
| `README.md` | for someone with pumpkyn and nothing else |

The builder derives this list from one table (section 7); the self-test's
dependency scan is what proves it complete.

## 4. The fly layer

**`[B, orbits] = catalog_endpoints(cat, sheetIdx)`** -- the certified endpoint
rule from the catalog's own keys: `get_family_orbit(dep_family, dep_params)`
and `(arr_family, arr_params)` with the catalog's `muStar`, then
`phase_state` on each period. Returns closures `B.stateD(sD)`, `B.stateA(sA)`
(6 x 1, phases as FRACTIONS of a period, wrapping), `B.Tnd`, `B.cnd`, `B.mu`
(from `thruster`, `rungs_N`, `constants`), and `orbits.tD/rvD/tT/rvT` for
plotting. A test asserts bitwise-equal endpoints against the repo's
`catalog_setup_request` + `arclength_arrival('setup', physicsOnly)` on the
record, so the library's certification carries over to the shipped rule.

**`F = fly_catalog_entry(cat, entry, opts)`** -- `entry` is a `transfers_at`
row, or `struct('iD', iD, 'iA', iA)` for a primary by cell, or a z8 with
phases. `opts.method`: `'junctions'` (default when the entry has them;
`entry_thrust_program`, segment by segment) or `'z8'` (`fly_transfer`, one
flight, what tfMin does). Returns `t`, `Y` (n x 14), `u` (thrust direction),
`massKg`, `tfDays`, `dvKms`, `propellantKg`, `flyKm`, `flyVms`, `defectKm`/
`defectVms` (junction method), `status`, `statusName`, `reason`, `method`,
plus the fields `print_transfer_summary` and `plot_transfer_3d` read (a test
runs both on the output).

**`[ok, dz, z8out, info] = tfmin_accepts(cat, entry, opts)`** -- hands z8 and
the endpoints to `pumpkyn.cr3bp.tfMin`, reports |dz| against the catalog's
gate (`opts.tolDz` [1e-6]) and the witness's flown miss. Fenced by a time cap
when a pool is available. A status-1 entry (witness failed at certification)
is reported `ok = false` with `info.expected = true` and a sentence, never an
error.

**Conventions (README + every header):** phases are fractions of a period in
`transfers_at` and the fly layer, DAYS in `costate_catalog_pick` (its existing
contract); ND unless the field name says otherwise (`tfDays`, `dvKms`,
`massKg`, `flyKm`).

**Tests:** endpoints bitwise-equal to the record's rule; a primary flown both
ways arrives within 1 km / 0.1 m/s; thrust direction against the equations of
motion (the Task-4 oracle, several samples); junction defects below 0.05 km;
a status-1 alternative -> `ok = false`, `expected = true`, no error; summary
and 3D plot run on `F`; the z8-with-phases input form.

## 5. The survey layer

**`costate_lib_describe(cat)`** -- current printout kept; with the layer
present it adds the status table for primaries and alternatives (counts per
code, borderline rows named), the junction count and K, and the `status_key`
legend. Without the layer: today's output, byte for byte (test).

**`costate_catalog_extremes(cat, filt, doPlot, metric)`** -- new filter fields
`filt.includeAlternatives` [false] and `filt.status` [all]; defaults reproduce
today's answer exactly (test). With alternatives included, an extreme may be
an alternative and the report prints its status and reason beside the usual
coordinates; the returned extreme carries the row so the movie renderer can
film it.

**`plot_status_torus(cat, opts)`** -- the phase torus (departure x arrival), one
panel each: primaries' flight time; alternatives per cell; cells holding a
status-2 transfer (the first map of where the conjugate points sit).
Off-grid alternatives drawn as markers at their phases. `opts.out` saves a
PNG (used by the README). Test on a synthetic 2 x 3 catalog with known counts
+ a smoke run on v2.

## 6. Movies and examples

**`costate_catalog_extremes_movies`** -- as deliverables 3-7, rendered through
`fly_catalog_entry` so alternatives can be filmed. The sample pair shipped:
one cell's primary (sufficient) beside its fastest status-2 alternative, same
phases, same clock, titled with both statuses; house polish rules; 1280 x 720
frames (the H.264 streak rule).

**`costate_catalog_example.m`** -- the first run, deliverable-3 shape: load,
describe, pick by phases (days), derive units, rebuild orbits
(`catalog_endpoints`), fly (`fly_catalog_entry`), `tfmin_accepts`, plot.

**`status_layer_example.m`** -- study-script style: `transfers_at` at a cell
with several transfers, the ranked table with statuses, fly the primary and
one alternative from their junctions, both thrust programs plotted, the
status-2 row's reason shown, `plot_status_torus`.

Both load the catalog by bare filename ("sits beside this script"); no
absolute paths.

## 7. The builder and verification

**`build_dro_deliverable(cfg)`** extended: `cfg.catMat` defaults to v2;
`cfg.name` defaults to `costate_catalog_70mN_status_deliverable8`. Stages:
export (`export_deliverable_catalog`), copy helpers + dependencies from ONE
list, render the sample movie, write the README with every number computed
from the exported catalog (counts, ranges, status table, audit results read
from `audit_70mN.mat` of the merged library and the v2 `audit_alt_*.mat` /
verdict, each bound by content key), refuse to zip unless the verdict file
says CLEAN and both keys match the export.

**README sections** (deliverable-7 shape plus): what the entries are; contents
table; dependencies; quick start; what it covers (orbits with recipes, thrust,
grid, counts, t_f and Delta-V ranges); **the statuses** (the four codes, what
each establishes, the borderline rule and the five rows); **phase
conventions** (fractions vs days); entry-field reference (sheet, alternatives,
legends); how it was made + verification numbers (certification gates,
audits, worst cases); gotchas (interpolation between grid points NOT included
and why; status-1 rows and tfMin; off-grid alternatives; `source` wording).

**Verification order, before the zip:**
1. Unit tests for every new or changed function (test-first; deliberate-failure
   checks), run from the repo.
2. **Recipient self-test**, a fresh `matlab -batch` with only pumpkyn +
   pumpkynPie and the deliverable folder on the path: both examples end to
   end, every helper's bare self-demo, a dependency scan
   (`matlab.codetools.requiredFilesAndProducts` on every `.m`) proving no
   file outside the folder and pumpkyn is called. A failure fixes the
   deliverable, not the test.
3. The library-deliverable skill's checklist: paths localized, README
   sections present, `zip -qr` from the parent, `unzip -l` verified, no
   `.DS_Store`.

## 8. Testing summary

Per function, in `DRO_tulip/indirect/tests` (or `costate_common/tests` where
the function lives): `test_export_deliverable_catalog` (keys equal, schema
valid, paths gone, source mapping), `test_catalog_endpoints` (bitwise vs the
record's rule), `test_fly_catalog_entry` (both methods, oracle, status-1,
consumers), `test_tfmin_accepts` (a primary accepted; a status-1 expected
refusal), `test_describe_status` (with/without layer), `test_extremes_status`
(defaults unchanged; alternatives on), `test_plot_status_torus` (synthetic
counts), `test_build_dro_deliverable` (file list complete vs the dependency
scan; refuses without a CLEAN verdict), plus the recipient self-test script
`deliverable_selftest.m` kept in the repo (not shipped).

## 9. Out of scope

- Interpolation between grid points (`phase_catalog_interp` and the junction
  experiment): shipped later, once its hit rate is measured on this library.
- Status layers for the older catalogs (deliverables 3-7).
- A MATLAB package namespace (option C); the function names stay as Darin
  knows them.
