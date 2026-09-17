# campaign_common — job control for long MATLAB campaigns

How a campaign RUNS, not what it computes: a disk work queue that N worker
processes pull from, process-held locks, a supervisor that keeps the workers
alive and finalizes once, heartbeats, atomic publishing and resumable
checkpoints. Moved out of `../costate_common/` on 2026-09-16 (its cleanup
step 11), because none of it is optimal control.

The rules these files enforce, and the incidents that paid for each one, are
in `../doc/CAMPAIGN_DISCIPLINE.md`. The `matlab-campaign` skill carries the
same discipline for ad-hoc runs.

**Dependencies.** The nine MATLAB files call only each other and MATLAB /
Java; nothing here calls `costate_common` or a campaign. The two shell
scripts find each other through their own directory. The execution fences
(`run_capped`, `capped_pool`, `current_pool`) stayed in `costate_common`,
because the certifier needs them.

**Consumers (2026-09-16).** All in `DRO_tulip/indirect`:
`run_costate_library` (queue, campaign lock, generated worker and finalize
jobs, which receive this folder among their code roots; it launches
`campaign_supervisor.sh` from here), `build_ribs`, `build_arrival_sheet`,
`fill_holes_direct` (`publish_atomic`), `rib_from_crossing`
(`walk_checkpoint`, path supplied by `build_ribs`), `run_phase_torus` (drives
`run_costate_library`; finds a live supervisor by its script name), and four
archived `batch/torus/v*_job.m` scripts (`fmt_num`).

## Contents

| file | what | used by |
|---|---|---|
| `work_queue.m` | A disk work queue so N worker processes on one host pull the next unclaimed unit (a rib column, an entry, a rung) instead of static ranges. Ownership is a process-held kernel lock (`unit_lock`), never transferred on heartbeat age; attempts are recorded before work starts, and a unit that always fails cannot livelock the queue. `tests/test_work_queue`, `tests/test_campaign_processes`. | DRO |
| `campaign_worker.m` | One worker process: claim a unit, run it, VALIDATE its temporary output and publish it in one rename under the unit's lock, release, repeat. A unit function that returns without a valid output is a failed attempt, not a done unit. Launched N at a time by `run_campaign_workers.sh`. | DRO |
| `unit_lock.m` | A process-held file lock (java.nio `FileChannel.tryLock`, POSIX record locking) whose authority lives in a process-wide REGISTRY keyed by the file's device:inode, so a released or stale handle can never act on a newer holder's lock. Measured: refused while held, acquired on release and after `kill -9`. | DRO (+ internal) |
| `publish_atomic.m` | Move a finished file onto its final name in ONE rename(2) (java.nio `ATOMIC_MOVE`), so a reader sees the old file or the new one, never a partial or missing one. Not `movefile`: onto an existing directory it nests the source inside it. The only way a campaign artifact, queue record or heartbeat reaches its final name. | DRO (+ internal) |
| `walk_checkpoint.m` | A resumable checkpoint for a sequential walk (a rib column): the state after the last accepted point, saved atomically at `<output>.ckpt`, carrying the unit's identity and refused by name on a mismatch, so a column killed hours in resumes from its last point. | DRO |
| `campaign_heartbeat.m` | Heartbeats as the only liveness signal, and the five states a monitor must tell apart. No heartbeat is NOT evidence of health -- a missing log once read as "alive, age 0 s", and a `pgrep` on an environment-only tag reported live workers dead. | internal |
| `campaign_status.m` | What a campaign is actually doing, read from ARTIFACTS ONLY (queue outputs and heartbeats, no process matching): units not being worked on, and workers that are not working. Entry point. | entry point |
| `safe_report.m` | Run a REPORTING block so it can never fail the work it reports on: catch, say so loudly, return false. Twice a correctly saved stage was reported FAILED because its summary print threw. | DRO |
| `fmt_num.m` | Format a number at a fixed width, rendering NaN or empty as dashes of the same width instead of throwing -- the two shapes that broke a campaign print, and the reason a movie title or table column does not jump frame to frame. | DRO |
| `run_campaign_workers.sh` | Launch N workers against a queue, each under a supervising PARENT shell that waits for READY, kills on heartbeat silence, confirms the exit and records it. Usage: `run_campaign_workers.sh <jobScript> <nWorkers> <hangSec> <outDir>`. | `campaign_supervisor.sh`, tests |
| `campaign_supervisor.sh` | Keep N workers alive from the queue's own files, within a launch budget; when every expected output exists run the finalize job exactly once; BLOCKED / STALLED verdicts otherwise. Usage: `campaign_supervisor.sh <jobScript> <nWorkers> <hangSec> <outDir> <queueDir> <maxAtt> [finalizeJob]`. | `run_costate_library` |

## Tests

- `tests/test_work_queue` — the queue, locks, heartbeats, reports,
  checkpoints and atomic publish in one process (seconds). Its rib-validation
  checks use `DRO_tulip/indirect/rib_validate` as a real unit validator and
  put that folder on the path for that block only.
- `tests/test_campaign_processes` — the guarantees that only mean anything
  across PROCESSES, with real MATLAB workers launched by the scripts above
  (about six minutes).

After any change here, run both, then a small `run_phase_torus` acceptance
campaign (the 3 x 3 `torus3b` size) to its fixed point: the unit tests
cannot see a generated job that lost a code root.
