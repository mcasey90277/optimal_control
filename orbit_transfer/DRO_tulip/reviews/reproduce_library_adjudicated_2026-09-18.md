# reproduce_library_70mN and the P0 repairs: GPT-6 Astra (xhigh), adjudicated (2026-09-18)

**Bundle.** `reproduce_library_70mN.m`, `compare_phase_catalogs.m`, `adopt_walked_arcs.m`, the repaired `run_phase_torus.m` in full, the diff of the other files touched by the P0 and Track A commits, and the unchanged callees the script's assumptions rest on. Raw API, xhigh, 154 KB, 430 s, 48,716 prompt / 17,001 completion tokens, $1.46. Transcript: `reproduce_library_astra_2026-09-18.md`.

**Astra's verdict: NO-GO as supplied**, with a six-item GO list. The host agrees with the verdict as of the bundle, and splits the findings by what a FRESH, UNINTERRUPTED build at the shipped operating point needs.

## What Astra confirmed correct
Every spec field is one the driver accepts; ten adopted arcs make `need` empty and match the sheet's pattern; `discover = false` disables discovery but not filler registration or promotion; the campaign ends `done` at the first round that registers no spine root and adds no anchor (so not necessarily three rounds) and publishes `final` either way; the wrapper's quoting, `echo $!` / `wait $!` and the `.done` race are sound for a foreground `matlab -batch`; `roundExtras` + `unique(..., 'stable')` is valid; the closures-only changes and `commonOpts` repair the operating-point defect; the greedy family pairing reports zero differences if and only if the partitions are identical (greediness can only overstate the count for unequal partitions).

## Launch blockers: fixed today, test-first
| # | finding | verdict | repair |
|---|---|---|---|
| 1 | The H2/H3 bound fields were added to a certificate only on the success path, so a failed and a successful certificate had different field sets; the sheet's seed loop (`build_arrival_sheet.m:119`) and the rib walker (`rib_from_crossing.m:183`) append certificates raw, and MATLAB refuses that. With `librarySeeds = true` the full build would have thrown in its first sheet. | **CONFIRMED** (the host's own defect from Track A; RED shown live) | both fields born in the initializer; `test_certify_schema` appends a pass, an early refusal and an H2 refusal into one array |
| 2 | The audit's new validators were weaker than the certifier's: `Inf > floor` is true, `logical(2)` is true, and `logical(NaN)` THROWS instead of making a BAD row. | **CONFIRMED** (RED: "NaN values cannot be converted to logicals") | `realScalar` before every comparison; `h6Ok == 1`; five malformed-value cases in `test_audit_fail_closed` |
| 3 | The comparator could pass on bad data: `NaN > tol` is false, so a NaN flight time or costate compared EQUAL; a missing family map reported "0 differences"; two catalogs of different engines on one grid were compared; the "cells agreeing" arithmetic was wrong. | **CONFIRMED** (RED: a NaN flight time returned a match) | `.nNonfinite`, `.familiesCompared` (required for a match), problem identity asserted, `.nAgree` counted over cells, table columns reshaped |
| 4 | The z8 metric put t_f in the norm with the costates, diluting costate changes. | **AGREED** | costates compared on their own (1:7); t_f has its own line |
| 5 | `REPRODUCED: PASS` could pair the final catalog with an unrelated audit (the newest `audit_*.mat` found in any round; `nBad == 0` true of an empty audit). | **CONFIRMED** | `finalAudit`: the LAST round of the driver's state, and it must cover every entry of the final catalog; otherwise NOT ESTABLISHED |
| 6 | The adopted registry's fields must match the driver's append exactly. | **CHECKED: they do** (`sD sA tfDays z src`) | none needed |
| 7 | `compareOnly` demanded the anchors; `.adoptArcs = false` does not re-walk arcs already in the folder. | **CONFIRMED** | the anchor assert skips for a comparison; re-walk asserts an empty arc folder |
| 8 | "Run the new audit against the record before committing the rebuild budget": the stricter audit now ABORTS a round, and nothing showed the record passes it. | **AGREED** (also the host's own step 2) | running: `results/audit_failclosed_record_2026-09-18/`, chunked and resumable |

## Accepted as limitations of THIS run; recorded for the next rewrite
These matter when a run is interrupted, resumed or mixed with other inputs; a fresh folder and an uninterrupted run do not exercise them. None is dismissed.
- **Final publication is not transactional** and terminal status is committed before `final` exists (also finding P1-12 of the 09-17 adjudication). A crash in those few seconds strands the campaign.
- **The fixed-point decision uses this invocation's counters**; a resumed round can declare `done` on a sheet that predates a registration. **Registration precedes promotion**, so a crash between them suppresses the promotion; and in discovery a second seed returning the same root within the dedup tolerance replaces `best` with `added = false`. Discovery is off in this run; the filler path shares the first half of the defect.
- **PID existence is not job ownership**, legacy job records are dropped rather than migrated, and jobs are recorded after all are launched. **Supervisor adoption judges a verdict by mtime.**
- **`stageOutcomes` infers execution from switches and blocker prefixes**; NaN = not run is ignored by `assertPackaged` although these callers require all three stages.
- **Arc adoption checks bytes, not identity** (`A.anc.sD`, the anchor, `A.B.problem`), and any extra `arrival_arc_<tag>_*.mat` in the folder enters the sheet. The ten arcs here are the record's own, copied into a fresh folder.
- **Root dedup by t_f within 1e-3 d** (known, 09-17 P1-10).
- **The H2/H3 bounds are estimates, not enclosures.** Astra ranks this first. The code and the audit document already say so in those words; the certifier's refusal text now should too. Accepted: a validated enclosure is out of scope, and on this library the margin is five orders of magnitude.

## GO criteria, as the host holds them
1. blockers 1-7 fixed and green (done);
2. the 3 x 3 rehearsal reaches a fixed point under the new wiring;
3. the fail-closed audit of the record completes, and its BAD rows, if any, are understood before the build is launched.
