# Phase-torus runbook: how the 70 mN DRO -> tulip 24 x 24 library was built

The reproducible procedure behind FINDINGS 59-72 (2026-09-13 .. 2026-09-15),
written so the next torus -- another engine, another orbit pair, another
grid -- follows the same steps without the trial and error. The methods
themselves are explained in `doc/phase_torus_methods.tex`; this file is the
operational recipe: what to run, in what order, what to watch, what each
failure looks like, and what it cost.

Library of record for this torus: `indirect/results/library_70mN_24x24_final/`
(README inside). The rounds' rib files: `indirect/results_fine_v2 .. v6`.
The exact jobs: `indirect/batch/torus/` (README inside).

## 0. Vocabulary

| term | meaning | file |
|---|---|---|
| **cell** | one (departure phase sD, arrival phase sA) pair on the grid, both in [0,1) as fractions of the orbits' periods | |
| **certified entry** | a root of the minimum-time PMP boundary-value problem that passed the whole gate stack (`certify_root`): multiple-shooting residual, flown arrival in position and velocity, the free-time conjugate test, the BCT hypothesis gates (min lambda_v, Q_mt, dim S = 1 with a 10x lift margin), transversality lambda_m(t_f) = 0 on a tight ode113 flight, lunar clearance 1900 km | `indirect/certify_root.m` |
| **anchor** | one certified root at sD = 0 from which arcs start (`.mat` with `best.z`, `best.it.Y`) | `indirect/results/mintime_70mN_anchor_*.mat` |
| **arc** | a pseudo-arclength continuation of the root in the ARRIVAL phase at sD = 0, one direction, recording every crossing of the grid's arrival levels and every fold | `costate_common/arclength_ms.m`, `indirect/arclength_arrival.m` |
| **sheet** | per arrival phase, every certified candidate from every arc's crossings and every seed, and the fastest one as the column's spine | `indirect/build_arrival_sheet.m`, `sheet_from_arcs.m` |
| **rib** | the departure-phase walk off one spine: sD in steps of 1/nD, each step re-solved and certified, bisection on failure | `indirect/rib_from_crossing.m`, `build_ribs.m` |
| **family** | one connected branch of extremals, named by its anchor; the arcs from that anchor trace it | `indirect/family_map.m` |
| **catalog** | the packaged library: fastest certified entry per cell, the verdicts, the family map | `indirect/package_phase_catalog.m` |
| **sidecar** | the second-order sweep's resumable record file, keyed by cell and z8 | `costate_common/second_order_pass.m` |
| **round** | one pass of: sheet -> ribs -> package -> audit -> sweep; a new round whenever a new family changes spines | `indirect/run_costate_library.m` |

## 1. Prerequisites

- MATLAB **R2026a** for every batch job (the Parallel Computing Toolbox is licensed there; under R2025b `certify_root` dies on `gcp`). `startup.m` builds the path; jobs `cd` to pumpkynPie and call `startup()` first.
- CasADi 3.7 on the path for the direct solves (`~/casadi-3.7.0`).
- 16 cores. A round with four rib workers, an arc pair and a finalizer runs at load 30-40 and everything slows; it still finishes.
- Never shell out to `git` from a job (`run_costate_library` reads `.git` itself); if `/usr/bin/git` says "You have not agreed to the Xcode license", run `sudo xcodebuild -license accept`, and until then commit with `/Library/Developer/CommandLineTools/usr/bin/git`.

## 1b. The one-call route: `run_phase_torus`

Since 2026-09-15 the whole loop below is one entry script:

```matlab
spec = struct('sD', [0 1/3 2/3], 'sA', sort(mod(0.0754 + (0:23)/24, 1)), ...   % any phases in [0,1)
    'orbits', struct('tauDRO', 1, 'NpTulip', 7, 'pmTulip', -1), ...
    'engine', struct('thrustN', 0.070, 'ispS', 900, 'm0kg', 150), ...
    'anchors', {{'anchor', 'results/mintime_70mN_anchor.mat', 0.0754, 'fast'}}, ...  % one certified root at sD(1)
    'outDir', 'results/my_torus', 'tag', 'mine', 'nWorkers', 4);
run_phase_torus(setfield(spec, 'plan', true));   % prints round 1's plan, launches nothing
out = run_phase_torus(spec);                     % rounds until nothing changes; out.final = the library
```

Each round is a `run_costate_library` campaign in `<outDir>/round_NN`
(arcs spawned as batch jobs, the sheet at the listed phases, ribs for
changed columns, finalizer, holes + improve, then discovery of new
families by direct solves at empty or slow columns); the driver resumes
from `<outDir>/torus_state.mat` and copies the last round to
`<outDir>/final`. Sections 2-9 describe what it does at each step and
how to watch it; they remain the manual route when a step needs a hand.

## 2. Declare the campaign (section 0 of `run_costate_library`)

`run_costate_library(struct(...))` is the entry point. Section 0 of the file
is the ONE place the orbits, engine and grid are chosen; the call can pass
`nD`, `nA`, `sD0`, `sA0`, `outDir`, `nWorkers`, `extraRibFiles`, `run`
(which stages), `launch`. The grid is a LATTICE: `sD = sD0 + (0:nD-1)/nD`,
`sA = sA0 + (0:nA-1)/nA`. The chain script `build_70mN_library.m` holds
the anchors table (name, seed file, sA0, family label) and the arc/rib
budgets; every family the campaign finds gets a row there.

A campaign directory (`outDir`) carries a manifest (written on first use,
checked on every call: same orbits, engine, grid or refusal), a
`campaign.lock`, the queue (`ribq/`), heartbeats (`hb/`), the generated
`rib_unit_job.m` and `finalize_job.m`, the sheet, the ribs, the catalog,
the receipt, the sidecar and the pictures.

## 3. The first anchor

You need one certified root at sD = 0. Routes used:
- `run_dro_tulip` (the front door) walks to one from the shipped catalog;
- a direct HS+Sundman solve (`casadi_mintime_dro`, N = 800, Sundman, lunar clearance enforced) warm-started from any nearby certified flight, then `harvest_ms_seed(o, 24)` -> `certify_root` (33 s when it converges).

Save it as `results/mintime_70mN_anchor_<name>.mat` with `best.z`,
`best.it.Y` (14 x K), `best.sA`, `best.sD`, `best.tfDays`, `best.origin`
and `Tnd`, `cnd`. The anchor's phase may be off the lattice (0.8671 vs
0.8670667); the arcs cross the lattice levels regardless.

## 4. Arcs

One batch job per (anchor, direction) -- templates `batch/torus/arc_*.m`:

```matlab
so = struct('thrustN', 0.070, 'ispS', 900, 'm0kg', 150, 'tauDRO', 1, 'NpTulip', 7, 'pmTulip', -1, 'sD', 0, ...
            'sA0', sA0, 'anchorMat', anchorFile);
[~, anc] = arclength_arrival('setup', so);
ao = so;  ao.direction = +1;  ao.sAStop = sA0 + 1.15;          % or -1 / -1.15
ao.levels = sA0grid + (-24:48)/24;  ao.nStep = 4000;  ao.deadlineSec = 6*3600;
ao.logFile = '...long.log';  ao.partialFile = '...long.partial.mat';  ao.saveEvery = 50;
A = arclength_arrival(anc, ao);  save('arrival_arc_<name>_<dn|up>_long.mat', 'A', '-v7.3');
```

- 4000 steps take 2-3 h at load 10, 5-6 h at load 40. `deadlineSec` ends a walk cleanly.
- **Partial saves are not optional**: the first direct18 dn walk died in a MATLAB segfault at step ~2160 with nothing on disk. `.partialFile` writes the arc so far every 50 steps, atomically; the partial is readable mid-walk (`A.crossings` already carry t_f) and the sheet builder ignores `*.partial.mat`.
- Read an arc before using it: span, t_f span, folds, min |rho| (normality), stop reason -- `family_map` prints all of it; a per-level table of `arc t_f vs sheet best` (see the session's `arrival_arc_*` reads in FINDINGS 67-68) says at once which columns will change.
- The arc file name is the family's name: `arrival_arc_<anchor>_<dn|up>_long.mat`; the sheet builder globs `results/arrival_arc_*.mat`, so a new pair is picked up by the next sheet rebuild automatically.

## 5. A round

Template `batch/torus/v5_campaign_job.m` (rounds 3-6 are the same script with the version bumped):

1. `run_costate_library(struct('outDir', Vn, 'nD', 24, 'nA', 24, 'sA0', 0.0754, 'launch', false, 'run', struct('sheet', true, ...)))` -- the sheet: every arc's crossings re-scanned at the grid levels, every library seed (`dro_tulip_library`, which includes `results/mintime_70mN_direct_certified.mat`), each candidate certified, duplicates merged only on equal t_f AND equal z8, the fastest certified per column. 1-3 h (about 300 candidates, a pool of 4).
2. Compare `S_new.TF` with the previous round's: columns whose spine changed are walked; the others reuse the previous round's rib file (copy it into `Vn`).
3. `run_costate_library(struct(..., 'nWorkers', 4, 'launch', true, 'extraRibFiles', {every earlier round's fine_rib_col*.mat}, 'run', struct('sheet', false, 'ribs', true, 'package', true, 'audit', true, 'sweep', true)))` -- launches `campaign_supervisor.sh` with four supervised workers on the queue; when every column is published the supervisor runs `finalize_job.m` once: package (every rib file of every round offered; the fastest certified point per cell wins), audit (every entry re-flown and re-gated, ~10 s each), sweep (second-order measurements written back), pictures, verdict.
4. Watch it (section 8). Read the verdict: `SUPERVISOR_VERDICT.txt`, `finalize_job.out` (CHAIN CLEAN / BLOCKERS), the receipt.

A rib is 23 points at 2-10 min each; the slowest column sets the round's length (3-5 h). A rib that stops early says why in `R.stop` ("normal-chart polish did not converge", "dense conjugate scan not clear"); the cells beyond it are holes for section 7.

## 6. The gap protocol: when a column is empty or slow

This is where the trial and error was. The rule that came out of it:
**continuation follows its own sheet; a faster family is found by a
branch-blind direct solve warm-started from the OTHER family at the phase
in question, never by walking further.**

1. At an empty or slow arrival phase, run a direct solve (`direct_gap_probe3.m` pattern) seeded from the nearest certified root of any family: fly its z with `tfMinProp`, interpolate states and directions onto N + 1 nodes, `casadi_mintime_dro` with `minAltKm` = 1900 - 1737.4. 1-5 min when it converges; 10 min of junk (100-400 d) when the seed's arrival geometry is wrong -- cap at 300-900 s.
2. Harvest and certify at the EXACT grid phase (`harvest_ms_seed` -> `certify_root`; a seed 3e-5 off the lattice fails the sheet's grid assert).
3. Append the certified root to `results/mintime_70mN_direct_certified.mat` (`direct(end+1) = struct('sD', 0, 'sA', ..., 'tfDays', ..., 'z', ..., 'src', ...)`); the next sheet rebuild seeds from it.
4. If it is faster than the sheet, make it an anchor (section 3) and walk both arcs (section 4); add the row to the chain's anchors table.
5. Rebuild the sheet (a new round, section 5): the arcs' crossings enter, the family map attaches every column, and the changed columns get ribs.
6. Repeat while a `family_index = -1` root (a certified root no arc passes through) is faster than its column's family.

What this found here: three families beyond the two the first campaign
knew (fast2 at 0.87, direct18 at 0.78, direct11 at 0.49), each 2-7 days
faster than the family the continuation was on, over 0.49-0.78 and
0.78-1.03. A refused direct result is also information: a solution whose
periselene sits ON the clearance floor is a constrained arc, not an
unconstrained extremal (0.5754 at 20.0 d; the 21.45 d root stands).

## 7. Holes and the improve pass

After the ribs, `fill_holes_direct(catMat, opts)` (template `batch/torus/fill_holes_job.m`):

- every empty cell is solved directly, warm-started from up to three certified neighbours -- **the same column's first** (they share the arrival geometry; a rib IS this continuation), then the same row's -- harvested, certified; the fastest certified root is kept and becomes a seed for the cells still to come, so a column chains like a rib; 300 s per direct solve; resumable, failed cells retried on resume. 38 of 38 holes certified here (4 min per cell).
- `opts.improveDays = 2` then re-solves every filled cell more than 2 d slower than a column neighbour, from faster neighbours only, keeping a root only if faster, and queues the next slower neighbour when a cell improves. Column 15 went from 24-25 d to 17.1-17.8 d over nine cells; nine cells above sD 0.54 resisted from both sides (a fold in sD, probably -- open).
- the output `fine_rib_direct_holes.mat` is a rib file; re-package (section 5 step 3 with `run.ribs` on an already-published queue = package/audit/sweep only, template `v6_repackage_job.m`).

## 8. Monitoring discipline (what "don't let it get stuck" meant in practice)

- Every long job logs to a file and writes its product after every item (arcs: partial saves; ribs: `walk_checkpoint` per point; filler: after every cell; sweep: after every chunk). `matlab -batch` buffers stdout: never judge a job by its `.out`.
- A watcher per job (`batch/torus/watch_*.sh`): exits on the job's DONE line, its process disappearing, an error line, a silent log (45 min for a filler cell, 50 min for a rib worker's heartbeat), or an hourly tick; re-arm after reading. Chains (`chain_*.sh`) launch the next stage on the previous DONE line -- put them in script files: an inline `pgrep -f` pattern matches the shell that carries it, and patterns are case-sensitive (`MATLAB_maca64`, not `matlab`).
- A MATLAB batch job at 0 % CPU after its banner: `pgrep -P` its `matlab_helper` child; a stuck `system()` command shows there (the git shim, 2026-09-15).
- Two writers on one file corrupt it: `pgrep -fl <job>` before every relaunch.
- Heartbeats are per solve stage; a 10-minute-old rib checkpoint with a fresh heartbeat is a long certification, not a stall.

## 9. Reading the family map

`cat_.families` (a `family_map` output) lists each family's anchor, arrival-phase span (unwrapped and mod 1), t_f span, folds, normality floor and the kind of each end (fold / walk budget / normality lost / unwalked anchor). `sheet.family_index` codes every entry: `1..n` a mapped family, `-1` a certified root at sD = 0 that no arc passes through (a branch nobody has walked -- a candidate anchor), `-2` a rib point whose spine root could not be identified (a rib from an earlier round whose spine is not among this sheet's roots, or a direct-solved cell with no spine), `0` no entry. Percentages quoted in FINDINGS are shares of the 576 cells.

## 10. The round ledger for this torus

| round | dir | what | result |
|---|---|---|---|
| 1 | `results_fine` | first 24 x 24 from the one known family, old launcher | 406 entries, 19 columns |
| 2 | `results_fine_v2` | corrected transversality gate (tight flight), supervised queue | 409, audit 409/0, sweep clean |
| 3 | `results_fine_v3` | fast2 (0.8671) and the direct roots at 0.7837/0.8254/0.9087/0.0337 seeded; ribs 18/19/24 | 432 |
| 4 | `results_fine_v4` | sheet from all six arcs + exact-phase seeds: 24/24 columns; ribs 19-23 | 515, audit 515/0 |
| 5 | `results_fine_v5` | direct18 arcs: cols 14-17, 20 take 16.8-19.0 d; ribs | 534, audit 534/0 |
| 6 | `results_fine_v6` | direct11 arcs: cols 11-13 take 18.0-18.8 d; ribs | 538, audit 538/0 |
| 7 | (v6 re-packaged) | + `fill_holes_direct` 38/38 | **576**, audit 576/0, sweep refused stale sidecar |
| 8 | (v6 re-packaged) | + improve pass (col 15), sidecar merge | **576/576, audit 576/0, sweep clean: 0 crossings, H6 4.44x, lift 10.6x** |

Wall time: 2026-09-13 09:40 to 2026-09-15 14:42, of which about 30 h were
compute (arcs 6 x 3-5 h, ribs 4 rounds x 3-5 h, finalizers 8 x 2-4 h,
filler 5 h) and the rest diagnosis and code.

## 11. Costs and rates (for planning the next one)

| step | rate |
|---|---|
| arc, 4000 steps | 2-6 h (load-dependent), ~10 s per step |
| sheet rebuild, ~300 candidates | 1-3 h with a pool of 4 |
| rib point | 2-10 min; a 23-point rib 1-4 h |
| audit | ~10 s per entry |
| sweep | ~15 s per entry (new ones only, with the sidecar) |
| direct cell solve + certify | 1-5 min when it converges |
| package | seconds |
