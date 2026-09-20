# Rebuilding the 70 mN DRO → tulip costate library: a guide to `reproduce_library_70mN`

Written 2026-09-19, the day after the first full rebuild, for whoever runs it
next (probably a future Claude session with none of this in context). It says
what the script does, how long it REALLY takes, how to run it unattended, how
to read its verdict, and every trap the first run fell into. The record of how
it came to be is `DRO_tulip/FINDINGS.md` sections 78–84; the hand-driven process
it replaces is `PHASE_TORUS_RUNBOOK.md`.

## 1. The thirty-second version

```matlab
cd orbit_transfer/DRO_tulip/indirect
reproduce_library_70mN                                        % PLAN: prints everything, creates nothing
out = reproduce_library_70mN(struct('go', true));             % BUILD + compare + verdict   (about ONE DAY)
out = reproduce_library_70mN(struct('compareOnly', true));    % re-run the comparison on a finished folder (seconds)
```

- It builds a complete, certified 24 × 24 library (576 minimum-time transfers)
  into `indirect/results/reproduce_70mN_24x24/`, then compares it cell by cell
  with the library of record, `indirect/results/library_70mN_24x24_final/`.
- **It takes about 24 hours on this machine with 4 rib workers.** Not 6–10.
  See section 4.
- It needs the shared MATLAB path (pumpkyn/pumpkynPie via `startup`), CasADi
  (`~/casadi-3.7.0`), and **R2026a** (the Parallel Computing Toolbox: the
  certifier refuses to run its external calls without a pool).
- Never run the `go` build in the shared interactive session: launch it
  detached (section 5).

## 2. What it is, and what it is not

It IS the build of THIS library as one reproducible chain: one engine
(70 mN, Isp 900 s, 150 kg), one orbit pair (DRO τ = 1.0 → 7-petal tulip, branch
−1), one grid (24 departure × 24 arrival phases).

It is NOT yet a generator for a new orbit pair or engine. It starts from five
KNOWN anchor roots, adopts ten arcs ALREADY WALKED and seven seed roots ALREADY
FOUND. For a new problem those must be made first (`anchor_by_direct_solve`,
then let `run_phase_torus` walk the arcs and discover families); the pipeline
is still wired to `tauDRO / NpTulip / pmTulip` (FINDINGS 78, "Track B").

`build_70mN_library.m` is NOT this. That is the older chain script; today it is
only the packager (package, audit, sweep, pictures) that `run_costate_library`
calls. Its own rib stage is wrong for 24 departure phases. Do not run it to
rebuild the library.

## 3. What happens, section by section

The script is a function with six numbered sections; read it top to bottom.

| § | what | detail |
|---|---|---|
| 0 | switches | `.go` (false = plan), `.outDir`, `.adoptArcs` (true), `.discover` (false), `.nWorkers` (4), `.maxRounds` (3), `.compareOnly` |
| 1 | the problem and grid | `sD = (0:23)/24`; `sA = sort(mod(0.0754 + (0:23)/24, 1))`; `librarySeeds = true` |
| 2 | the FIVE families | `anchor` 0.0754 (fast), `cell11` 0.9087 (A2, the 26 d family), `fast2` 0.8671, `direct18` 0.7837, `direct11` 0.4921; all certified roots at departure phase 0, files in `results/mintime_70mN_*.mat` |
| 3 | arcs and seed roots | `adopt_walked_arcs` copies `results/arrival_arc_<name>_{dn,up}_long.mat` to `<outDir>/arcs/arrival_arc_70mN_<name>_{dn,up}_long.mat` (never overwrites); `results/mintime_70mN_direct_certified.mat` (7 roots) becomes `<outDir>/direct_certified.mat` |
| 4 | the build | ONE call: `run_phase_torus(spec)` |
| 5 | the comparison | `compare_phase_catalogs(new, record)` |
| 6 | the verdict | campaign status, the FINAL catalog's fail-closed audit, same-library PASS/FAIL |

Why five anchors and no discovery: the hand build FOUND its families one round
at a time. They are known now, so discovery has nothing to find and the run is
the deterministic part. With all five given it reaches its fixed point in ONE
round.

What `run_phase_torus` does in that round (each stage leaves files in
`<outDir>/round_01/` before the next starts):

```
arcs (adopted: nothing to walk)
 → arrival sheet      build_arrival_sheet: every grid crossing of every arc, and every seed root,
                      through certify_root; per arrival phase the fastest CERTIFIED root wins   → arrival_sheet_70mN_nA24.mat
 → ribs               campaign_supervisor.sh + 4 workers, one arrival column per unit, each walks the
                      23 other departure phases off the column's spine (rib_from_crossing)      → fine_rib_colNN.mat
 → finalizer          run_costate_library again, packaging stages only → build_70mN_library stages 5–8:
                      package → audit (FAIL-CLOSED) → second-order sweep → pictures             → catalog, receipt, audit_70mN.mat, sidecar
 → holes + improve    fill_holes_direct: a direct collocation solve per empty cell, seeded from its
                      certified neighbours, harvested, certified                                 → fine_rib_direct_holes.mat
 → re-package         package + audit + sweep again, with the direct cells
 → fixed point?       no new spine root registered, no new anchor → status 'done'
 → final/             the last round's products, copied
```

## 4. How long it takes (measured 2026-09-18/19)

| stage | time |
|---|---|
| sheet (≈300 candidates, pool of 4) | 1 h 34 |
| ribs (24 columns × 23 points, 4 workers) | ≈ 9 h |
| package + audit 519 + sweep 519 | ≈ 4 h |
| filler (57 holes + 1 improve, ≈ 9 min each) | 7 h |
| re-package + audit 576 + sweep of the new cells | 2 h 31 |
| **total** | **24 h 12 min** |

The runbook's "ribs 3–5 h per round" is for rounds that walk only the CHANGED
columns; a rebuild walks all 24. With `.adoptArcs = false` add ≈ 20 h (ten
arcs, 2–6 h each). More rib workers is the first thing to raise.

## 5. Running it unattended (the pattern that worked)

Follow the `matlab-campaign` skill. Concretely:

1. **Rehearse first** if any code under the driver changed: the 3 × 3 acceptance
   torus (`sD = [0 1/3 2/3]`, `sA = 0.0754 + [0 8 16]/24`, one anchor, `arc.nStep
   400`, `nWorkers 2`, `maxRounds 3`) reaches a fixed point in 2 rounds, ≈ 90
   min, 7 cells, and its catalog should be IDENTICAL to `results/torus3d_p0/final`
   (compare with `compare_phase_catalogs`). Fast tests do not prove a live round.
2. **Launch detached**, from a job file, with a verdict FILE:
   ```bash
   nohup /Applications/MATLAB_R2026a.app/bin/matlab -batch \
       "run('.../DRO_tulip/indirect/batch/reproduce_library_job.m')" > ~/reproduce_job.out 2>&1 &
   ```
   `indirect/batch/reproduce_library_job.m` is the job that ran on 2026-09-18,
   kept as a template: it builds the path, calls the script inside a
   `try/catch`, saves `out`, and writes a one-line verdict to
   `~/REPRODUCE_VERDICT.txt` whether it succeeded or threw.
3. **Watch three things**: `<outDir>/torus.log` (one line per stage — it is
   SILENT for hours during the sheet, the ribs and the filler, which is normal);
   the job's `.out` (the sheet's per-candidate lines and the audit rows land
   there, unbuffered enough to read); the verdict file / process death.
   During ribs: `<outDir>/round_01/supervisor.log`, `hb/`, `worker_*.log`.
   During the filler: `<outDir>/round_01/fill_holes_direct.log`.
4. **It resumes.** Call it again with the same options: adoption skips files
   already copied, the driver reads `<outDir>/torus_state.mat`, the sheet and
   finished ribs are reused. But see section 8 for what resume does NOT handle.

## 6. Reading the verdict

`compare_phase_catalogs` reports, for the 576 cells matched BY PHASE:

- **coverage** — cells in one catalog and not the other;
- **t_f** — flight time within 1e-6 d; if not, which catalog is FASTER
  (`.nNewFaster`, `.nNewSlower`);
- **costates** — the seven initial costates within 1e-6 relative (t_f is left
  out of this norm on purpose);
- **finite** — a NaN is a difference, never a match;
- **families** — family labels compared as a PARTITION (index numbers do not
  matter); "not established" if a catalog has no family map;
- **VERDICT** `ok` = the same library; **`noWorse`** = covers the reference,
  nothing non-finite, slower nowhere.

`REPRODUCED: PASS` needs `ok` AND the final catalog's fail-closed audit with
zero BAD rows covering every entry.

**What the first rebuild gave, and how to read a similar result.** 565 cells
the same root; 11 cells a different root, the rebuild FASTER in all 11 (the
filler found better roots than eight hand rounds had: rows 14–22 of the column
at sA 0.6587, plus two); 82 cells the same root under another family LABEL.
Strict verdict FAIL, `noWorse` true. That was judged a BETTER library and
adopted as the record on 2026-09-19 (FINDINGS 82–83). Family-label differences
with identical roots are not differences in the library: `family_index` is
provenance, rib labels attach by flight time and are not reproducible.
A rebuild against the NEW record should come much closer to `ok`; expect label
differences to persist until rib attachment is fixed.

## 7. Traps this work fell into (each is fixed; know them anyway)

1. **Two files named `run_costate_library.m`**: `DRO_tulip/run_costate_library.m`
   (the August thrust-ladder script) and `DRO_tulip/indirect/run_costate_library.m`
   (the one the driver calls). A `find | head -1` picked the wrong one for a
   review bundle. List ALL copies of a name.
2. **The record lists its arrival phases from the anchor** (0.0754 … 0.9921,
   0.0337); **the driver sorts them** (0.0337 first). Same 24 phases, columns
   cyclically shifted by one. Address cells by phase. (The adopted record is
   sorted; the archived hand-built one is not.)
3. **The packaging chain names its audit `audit_70mN.mat` whatever tag the
   driver has**; the script finds the last round's audit by pattern.
4. **MATLAB struct arrays need one field set.** A certificate field added only on
   the success path made the sheet's seed loop throw on its first mixed
   pass/fail. Every `certify_root` field is born in the initializer
   (`test_certify_schema`).
5. **`NaN > tol` is false**: a comparator or audit written as `if x > tol, bad`
   passes NaN. Validate finite first (`realScalar`), and treat a missing result
   as BAD (the audit is fail-closed now; it was fail-open).
6. **A bulk substring rename** (`capped(pool,` → `cap(pool,`) also renamed a
   function definition and `run_capped(`. Use a word boundary, and diff against
   git afterwards.
7. **The shared session's pool times out**; `certify_root` then refuses to run.
   `capped_pool(2)` first.
8. **Estimates**: see section 4.

## 8. Known limitations (from the Astra review of 2026-09-18, not fixed)

They bite INTERRUPTED or RESUMED runs, not a fresh uninterrupted one
(`reviews/reproduce_library_adjudicated_2026-09-18.md`):
final publication is not transactional and the terminal status is saved before
`final/` exists; the fixed-point decision uses this invocation's counters, so a
resumed round can declare `done` on a stale sheet; registration precedes
promotion; a recorded PID is "alive", not "ours"; the supervisor's verdict is
judged by mtime; arc adoption checks bytes, not the arc's identity, and any
extra `arrival_arc_70mN_*.mat` in the arc folder enters the sheet. If a run is
interrupted, read those before trusting a resume — or start a fresh `outDir`.

## 9. The pieces, and their tests

| file | role | test |
|---|---|---|
| `indirect/reproduce_library_70mN.m` | this script | `tests/test_reproduce_library` (plan mode, `finalAudit`) |
| `indirect/compare_phase_catalogs.m` | the comparison | same test: self-match, planted differences, relabelled families, reordered grids, NaN, other engine, faster/slower |
| `indirect/adopt_walked_arcs.m` | arc adoption | same test |
| `indirect/run_phase_torus.m` | the round-loop driver | `tests/test_run_phase_torus_p0` (through `run_phase_torus('localfunctions')`) |
| `indirect/run_costate_library.m` | one round's front door | `tests/test_run_costate_library_seams` |
| `indirect/audit_phase_catalog.m` | the fail-closed audit | `tests/test_audit_fail_closed` |
| `indirect/certify_root.m` | the gate stack | `costate_common/tests/test_certify_enforcement`, `test_certify_schema` |
| `indirect/family_map.m` | families; root identity by costates | `tests/test_family_map`, `test_family_costate_identity` |
| `indirect/phase_transversality_check.m` | costates vs a finite-difference test; time-consistent edges and candidate components (a heuristic, not an interpolation licence) | `tests/test_phase_transversality_check` |

The `localfunctions` seam: `H = some_file('localfunctions')` returns handles to
that file's local functions, so a test calls the real helper.

## 10. Where things live

- Library of record: `indirect/results/library_70mN_24x24_final/` (README inside;
  adopted 2026-09-19). Its predecessor: `library_70mN_24x24_handbuilt_2026-09-15/`.
- The first rebuild's campaign folder: `indirect/results/reproduce_70mN_24x24/`.
- The 3 × 3 rehearsal: `indirect/results/torus3d_p0/`.
- **`results/` is not tracked by git.** The record exists on this disk only.
- Reviews: `DRO_tulip/reviews/` (the 2026-09-17 and 09-18 adjudications are the
  ones that matter for this script).
