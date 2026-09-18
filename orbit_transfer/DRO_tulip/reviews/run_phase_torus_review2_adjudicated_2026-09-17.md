# run_phase_torus, second review pass: GPT-6 Astra (xhigh) + host, adjudicated (2026-09-17)

**Scope.** The driver `indirect/run_phase_torus.m` (743 lines, post-rewrite commits 1fb90dc3..26b08ea4) with every callee interface inlined: `indirect/run_costate_library.m`, `build_arrival_sheet`, `arclength_arrival`, `direct_cell_solve`, `fill_holes_direct`, `family_map`, `rib_targets`, `costate_common/capped_pool`, the header of `certify_root`; the first-pass review and its FINDINGS 73 adjudication; the runbook; both 3 x 3 acceptance logs (09-15 pre-rewrite, 09-16 current code).

**Astra run.** Raw OpenRouter API (`astra_raw.sh`), effort xhigh, 228 KB bundle, HTTP 200 in 423 s, 69,173 prompt / 19,951 completion tokens (14,420 reasoning), $1.86. Transcript: `run_phase_torus_astra_round2_2026-09-17.md`. A first launch was killed after 5 min: the bundle had inlined `DRO_tulip/run_costate_library.m` (the August thrust-ladder script) instead of `DRO_tulip/indirect/run_costate_library.m` (the struct-driven front door the driver calls). Two functions with that name sit one folder apart; the driver's `addpath(here, ...)` puts `indirect/` first, so MATLAB resolves it correctly, but a reader, a `which`, or a future path edit will not. **Rename or delete the August script.**

**Host pass.** The host read the driver in full and checked the disputed MATLAB semantics in the shared R2026a session: the sheet fields are rows (`S.TF` 1 x nA, `S.Z8` 8 x nA, `sheet_from_arcs.m:68`), so every `for j = find(...)` loop is safe; the MATLAB launcher script `exec`s its binary (`bin/matlab:1850`), so the recorded PID is the real process; the supervisor-adoption `pgrep -f` does not match its own shell (MATLAB's `system` spawns `zsh -c`, which execs a single command; empirically no match). Each Astra finding below carries the host's verdict and where it was verified.

## 1. The 17 "applied" claims from FINDINGS 73, re-verified

Astra: 5 CONFIRMED (arcs campaign-local; terminal campaigns remembered; physical defaults + `.rib.wallSec` forwarded; discovery try/catch; 1e-5 phase-list resolution) and 12 PARTIAL. The host agrees with every PARTIAL. The recurring shape is the one FINDINGS 36/41/47 named: the mechanism was built, the thing it was meant to guarantee is not what it checks. Three PARTIALs matter most:

| claim | what was built | what it does not do | where |
|---|---|---|---|
| manifest checked on resume | `sD, sA, orbits, engine, tag` compared exactly | `.arc.span`, `.arc.nStep`, the anchors table and anchor file contents are outside it; arcs are keyed by anchor NAME, so a resumed run reuses arcs walked under another span, and a changed `.anchors` input is silently ignored once `torus_state.mat` exists | `run_phase_torus.m:130-146` |
| PIDs recorded, jobs killed at deadline | `st.jobs` saved after launch | never read on resume: a killed driver restarted while its arc jobs still run respawns them for the same output files (two MATLABs, one `.part` -> one `movefile`) | `:134-141`, `:168-183`, `:421-444` |
| stale rib set aside when its spine changed | `ribMatchesSpine` | compares the rib's FIRST point to the spine within 0.5 d; the first point is the first departure TARGET (`rib_targets` removes the spine), so this is a heuristic, not identity, and only runs on `changed` columns | `:643-654`, `rib_targets.m:47` |

## 2. New findings, adjudicated

Verdicts: **CONFIRMED** = the host reproduced the mechanism at the cited lines (and in MATLAB where a semantic was in question); **PLAUSIBLE** = mechanism consistent with the code, consequence not demonstrated; **DISAGREE** stated with reason. Ranked by expected cost to a campaign.

### P0: would lose the campaign's products or hours of compute

1. **CONFIRMED. The sheet, the filler and discovery are set up at the SHIPPED operating point, not the campaign's.** `run_costate_library.m:239-243` passes `build_arrival_sheet` no `sD`, no `anchorMat`, no `sA0`; `arclength_arrival('setup')` then defaults to `sD = 0`, `results/mintime_70mN_anchor.mat`, `sA0 = 0.0754` (`arclength_arrival.m:113, 139`). A campaign with `sD(1) ~= 0` walks its arcs (hours) and then dies on the grid assert at `run_costate_library.m:270-272`; a campaign at another engine or orbit pair fails polishing the 70 mN anchor before its sheet exists. The filler (`run_phase_torus.m:222-223`) and discovery (`:560-561`) pass `st.anchors{1,2}` but not `st.anchors{1,3}`, so a first anchor not at 0.0754 is polished against the wrong arrival state. Latent today because every campaign has used `sD0 = 0` and the shipped anchor; fatal for `phase_torus_examples` example 3 (DRO tau 2 -> 8-petal) and for the any-orbit-pair goal. Both host forks flagged the same root in the pipeline (see the pipeline adjudication). **Fix:** thread `sD(1)`, `anchorMat`, `sA0` into `common` and from there into the sheet call; use `physicsOnly = true` for the filler and discovery (they need closures, not a polished anchor). Add an acceptance case with `sD(1) ~= 0` and a non-shipped anchor.

2. **CONFIRMED (host, independently). A failed arc can never be retried.** `runArcs` tests `.fail` before `.done` and `spawnArc` never removes a stale `.fail`, so a rerun errors on its first poll while the fresh job runs (`:430-435`). Astra generalises correctly: attempts share file names and verdicts, and a MATLAB that fails before the generated `try` (path, pool, license) writes no verdict at all, so the driver waits the full deadline + 30 min. **Fix:** delete both markers in `spawnArc`; write the verdict from a shell wrapper around `matlab -batch` so a launch failure is a `.fail` too.

3. **CONFIRMED (host, independently). A `budget` campaign cannot be resumed.** The stop reason says "raise .maxRounds and rerun" (`:253-255`); the status gate returns "nothing to do" for `budget` as well as `done` (`:148-152`). The only route is deleting `torus_state.mat`, which discards the anchors and rounds. **Fix:** treat `budget` as resumable when `maxRounds > roundsDone`.

4. **CONFIRMED. A resumed round drops the filler's roots.** `fill_holes_direct` skips cells already certified in its file and `fh.nCert` counts only this call's (`fill_holes_direct.m:148-156, 230`); the round's own `fine_rib_direct_holes.mat` reaches the packager only in the conditional re-package (`run_phase_torus.m:226-234`). A crash between the fill and the re-package therefore leaves a round whose catalog lacks cells that are on disk, and whose spine roots are never registered. **Fix:** offer the current round's holes file (if present) to the first `runRound`, register from `fh.file` unconditionally, and re-package when the file changed rather than when `nCert > 0`.

5. **CONFIRMED. Recorded PIDs give no resume ownership** (claim table above): duplicate writers on restart. **Fix:** on resume, for each missing arc check `st.jobs` liveness (`kill -0`) and adopt or retire before spawning.

### P1: wrong bookkeeping, wrong status, or a wasted round

6. **CONFIRMED. Discovery can anchor a root the filler just anchored.** The spine path promotes only when `added` (`:526`); the discovery path does not test `added` (`:594-604`), and it attaches against a family map built from `S.arcs` only, so a root promoted minutes earlier in the same round (arcs not yet walked) is "on no family" and gets a second anchor name and a second pair of arc walks. **Fix:** require `added` (or a registry/anchor lookup by root) before `promote` in discovery.

7. **CONFIRMED. Filler promotions are not counted.** `registerSpineRoots` promotes but returns only `nReg`; `newAnchors` comes from discovery alone (`:225-242`). The 09-16 log's "NEW ANCHOR d03_1" followed by "0 new anchor(s)" is this. The fixed-point rule still holds because `nRegistered >= 1`, so the consequence is a wrong README and round record, not a wrong stop.

8. **CONFIRMED. Last-round discovery is not in `final`.** Registry roots and anchors written by discovery in round `maxRounds` are never consumed by a sheet/package pass before `final` is copied (`:236-291`). The status is `budget`, correctly, but the library of record omits committed roots. **Fix:** one final sheet/package/audit/sweep generation from all committed roots before publication, no further arcs.

9. **CONFIRMED. `packaged` is not audit/sweep success.** `run_costate_library` sets `state = 'packaged'` when this invocation wrote a receipt, with audit/sweep blockers appended rather than fatal (`run_costate_library.m:518-531`); the driver asserts on `packaged` (`:217, :232`). **Fix:** return package/audit/sweep outcomes separately and require all three.

10. **CONFIRMED. Family attachment and registry dedup identify a root by t_f alone.** Registry: same `sA` within 1e-8 and t_f within 1e-3 d (`:635`); family map: t_f within 0.02 d of an interpolated arc (`family_map.m:240-267`). The sheet itself uses t_f AND z8 (`sheet_from_arcs.m:105-121`). Near a family crossing, exactly where discovery is interesting, a new root is "registered, not anchored" and never walked. Host and both forks flagged this. **Fix:** compare the normal-chart `lam0` direction too (the arcs store `p` at every point, so this is free).

11. **CONFIRMED. The sheet is rebuilt only when its file is absent** (`run_costate_library.m:238`), so a resume after a mid-round promotion walks the new arcs and then reuses the old sheet; the new family enters one round late (or never, in the last round).

12. **CONFIRMED. Terminal status is committed before `final` is built** (`:257` then `:264-291`), and `final` is removed before the required products are checked. A crash there leaves no library of record and a campaign that answers "nothing to do". **Fix:** build `final.part`, validate, rename, then save the terminal status.

13. **CONFIRMED. The fixed-point count ignores off-spine filler cells.** Only spine roots reach `nRegistered`; an off-spine cell certified this round can seed a hole that failed earlier, so the loop can stop one productive pass early. Astra's own caveat holds: `maxRounds` bounds every invocation and the 1e-3 d dedup means root jitter does not run the loop forever.

### P2: robustness, validation, hygiene

14. **CONFIRMED.** `withDefaults` keeps unknown NESTED fields, so `engine.thrustn` silently leaves 70 mN active; `maxRounds = 1.5` passes validation, runs round 1, never hits `k == maxRounds`, and publishes with status `running` (`:330-373, :704-706`). **Fix:** validate nested structs against their field tables; require integer counts and finite positive spans.
15. **CONFIRMED.** `run('%s')` at `:482` bypasses `mlq`; a relative `.startup` resolves relative to `jobs/` inside the child (`run` cds to the script's folder); `.startup`/`.matlab` are not forwarded to the rib and finalizer jobs, which hard-code pumpkynPie (`run_costate_library.m:665, 721`).
16. **PLAUSIBLE.** Supervisor adoption accepts a verdict by mtime, so an adopted supervisor that has already written its terminal verdict is waited on for 48 h. The mechanism reads correctly; not reproduced.
17. **PLAUSIBLE.** The anchor-name allocator counts prefixes rather than finding an unused suffix; an initial anchor named `d03_2` collides. Edge case.
18. **DISAGREE in part.** Astra lists the `pgrep -f` self-match as a live hazard. Verified in the shared session: no self-match under MATLAB's `zsh -c`. Keep the `[c]ampaign_supervisor` idiom anyway; it costs nothing.
19. **STYLE (host).** `j` is used as a loop index at `:195, :205, :438, :552, :565`; house rule is never `i`/`j`.

**Declined items (FINDINGS 73):** Astra reopened none. The host agrees.

## 3. What the acceptance runs did and did not exercise

Astra's reading of the logs is right and worth keeping: the 09-16 run's new anchor came from the FILLER, its discovery stage ran no probe (at three columns a third of a period apart, `seedRadius = 0.15` excludes every neighbour), and round 2 declared a fixed point with one freshly certified off-spine cell. Neither log exercises: a restart during arc publication, a late supervisor, a changed manifest, non-default physics or `sD(1) ~= 0`, or last-round discovery. Those are the acceptance cases to add before the driver is used for another orbit pair.

## 3b. Status, 2026-09-18

P0 1-5 and P1 6, 7 (filler anchors counted), 9 (`packaged` is not audit success) are FIXED, test-first; the pipeline adjudication's P0-3 (audit fails open) is fixed with them. Record and the list of tests: FINDINGS 79. Still open from this document: P1 8, 10-13 and all of P2 except the `j` loop variables (renamed) and the unquoted `run()` literal (now `mlq`-quoted inside the arc wrapper).

## 4. Fix order

P0 1-5 first (1 and 4 are small edits; 2, 3 are one-liners; 5 is a liveness check). Then 6, 8, 9, 10, 12. The rest with the next rewrite. The one-phase-path block (`run_costate_library.m:274-308`, blocked unconditionally when `nPts == 0`) belongs to the pipeline adjudication but blocks `phase_torus_examples` example 1 today.
