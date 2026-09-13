# Campaign discipline

Companion to `CERTIFICATION_DISCIPLINE.md`. That one is about what a result
must satisfy before it counts. This one is about running the machine that
produces results, and every rule below was paid for.

The reference incident throughout is the 24x24 DRO-to-tulip library
(2026-09-12/13), where the science was never wrong and roughly a day of
machine time was lost anyway.

## The primitives

| unit | kills |
|---|---|
| `work_queue` | static work assignment, lost units, livelock |
| `campaign_heartbeat` | monitors that report unverified health |
| `campaign_worker` | a worker holding a queue it cannot finish |
| `campaign_status` | liveness by process matching |
| `run_campaign_workers.sh` | launchers that silently do nothing |
| `safe_report`, `fmt_num` | reports that fail the work they report on |

## The rules

### 1. The unit of work is the unit of loss

A worker owns ONE unit at a time and claims the next from a shared queue.
Never hand a worker a list.

Units are wildly uneven: rib columns ran from 90 minutes to over five hours.
With static ranges, workers idled behind their own slow column while others
had nothing to do -- about five hours wasted in one afternoon -- and a
worker killed mid-list took its whole remainder with it.

### 2. Done is an artifact, never a flag

A unit is done when its output file exists. Nothing else records completion,
so nothing can disagree with what is on disk.

### 3. A dead worker must not take its unit with it

Claims go stale. A claim not touched for `staleSec` is reclaimable, so a
killed worker returns its unit to the queue automatically. Two workers were
killed by watchdogs holding columns 13 and 16 after nine hours each; nothing
could pick that work up, and both columns started again from zero.

### 4. Count attempts, and write the count BEFORE the work

A unit that always fails is released and immediately re-claimed, forever.
The first end-to-end test of the queue livelocked on exactly that, and took
the shared MATLAB session with it. After `maxAtt` attempts a unit is retired.

The count is written when the claim is taken, not when the work returns, so
a unit that hangs hard enough to kill its process still spends an attempt --
otherwise the next run retries it and hangs identically.

### 5. The watchdog must be longer than one unit

A watchdog exists to catch a hang, not to interrupt work. Sized shorter than
a unit it becomes the thing that destroys progress. `run_campaign_workers.sh`
derives it as 4x the measured unit time and REFUSES to launch if it is
shorter. The job's own budget stops cleanly BETWEEN units; the watchdog is
only the backstop.

### 6. Measure one unit before planning the campaign

The 24x24 plan was built on 25-40 seconds per certification. The real figure
was about 2.5 minutes. Every budget, watchdog and estimate downstream was
wrong by that factor. Run one unit, measure it, then size everything.

### 7. Liveness is an artifact's freshness, never a process match

Two monitor bugs in one campaign, both reporting health never verified:

* liveness was `pgrep` on the worker's tag, but the tag was in the
  ENVIRONMENT and never appeared in the command line, so four working
  processes were reported dead;
* a MISSING log read as "alive, age 0 s", because the age expression fell
  back to the current time. Three workers that never launched looked exactly
  like workers that had just checked in.

A monitor must distinguish five states and never collapse them: **never
started**, running, stalled, failed, done. `never` is the one that keeps
getting lost, and it is the one that means "go look".

### 8. A launcher must verify it launched

Three workers never started because a shell loop split a bracketed column
list and zsh globbed it. Nothing warned; the loss surfaced nine hours later.
The launcher now waits for each worker's heartbeat and reports loudly if it
does not appear. Related: pass arguments positionally and echo them back --
one mangled tag turned a column list into "all columns", and three workers
spent nine hours redoing each other's work.

### 9. Reporting must never be able to fail the work

Twice a completed, correctly saved stage was reported as a FAILED job because
its summary print threw: once on `char(string(NaN))` for a phase with no
certified transfer, once on `[c.ok]` where a column with no candidates holds
a plain double rather than an empty struct. Wrap reporting in `safe_report`
and format missing values with `fmt_num`. Both are ordinary states of a
sparse grid, not errors.

### 10. A gate blocks shipping, not measuring

One bad audit row used to abort the whole chain, so the second-order sweep --
the stage whose evidence a bad row needs interpreting WITH -- never ran.
Collect blockers, finish every stage that can run, and gate only the
deliverable.

### 11. Honour the contract of what you produce

A re-scanned crossing was a bracket, not a converged root, so 120 of 120 were
refused with "not converged" -- the producer's own flag read back, not a
solver failure. If a consumer expects a corrected root, correct it, using the
SAME corrector the original producer used.

### 12. Beat at every capped STAGE, and size the hang deadline from the largest stage

The first generated worker accepted a heartbeat callback and never called
it. The second called it once per certification -- which runs SEVEN
separately capped stages, 4500 s in all, so a healthy walk could outlast
the 2700 s deadline (pass 3). Now `certify_root` ticks the heartbeat after
every capped stage; the largest stage cap is 900 s and the deadline is
three of them. Derive the deadline from the longest possible silence
between beats, counted from the code, never from a column duration or a
guess. There is no calibration solve: it ran outside the queue's
ownership and could duplicate a live worker's column. A cancelled capped
call is VERIFIED to have stopped within 30 s, or the pool is deleted and
the attempt fails as infrastructure, not as a numerical wall.

### 13. Re-opening a campaign must not touch what workers own

"Re-run to package" is the normal path. If the entry point re-initialised
the queue, it freed columns being walked and reset the attempt counts that
stop the livelock. `open` is READ-ONLY and validates the unit set; `reset`
is a separate, explicit action that cannot touch a held unit.

### 14. Ownership is a process-held lock whose authority lives in a registry, not in the caller's struct

The first queue transferred a claim whose beat was 30 minutes old (an ABA
race). The second held a kernel lock but let the claim struct carry a
copied `held` flag: a released claim could still publish over the next
owner (pass 3, reproduced). Now `unit_lock` keeps a process-wide,
mlock'ed registry keyed by file IDENTITY (device:inode) that maps to the
live Java objects and the current holder's token; `holds`, `beat`,
`publish` and `release` ask the registry, so a released handle has no
authority whatever its copy says, and release is idempotent. POSIX frees
every lock a process holds on a file when ANY descriptor on it closes, so
the registry is consulted BEFORE any open, and an overlapping-lock report
parks the channel rather than closing it. Nobody can take a unit from a
live owner; a hung owner is killed by the supervisor, the kernel frees the
lock, and only then is the unit claimable. Measured: refused across
processes, released on exit and on SIGKILL, and (pass 3) dropped by a
second channel closed in the owner -- which is what the registry prevents.
`reset` holds the lock while it acts.

### 15. Publish under the lock, from an exclusive temporary, through one rename

The queue reads "artifact exists" as "unit done", so the artifact must
appear whole or not at all, and only the owner may put it there. The unit
writes to an attempt-specific temporary name the queue hands it; the WORKER
validates it and moves it onto the output while still holding the lock
(publish_atomic: java.nio ATOMIC_MOVE, which renames or throws -- movefile
onto an existing directory nests the source inside it). "Returned
normally" is not success; a validated publication is.

### 16. Say what you can prove: pending, launched, blocked, packaged, failed

The entry script returns a state. Packaging runs only when every column
validates as a rib and no lock is held, whether or not the rib stage ran in
this call. A short column (the walker stalled) is real data but not full
coverage, and is named as a blocker. Unreadable telemetry is UNKNOWN.
Unseeded columns of the requested grid are named, not dropped.

### 17. One controller at a time; one queue per campaign; drain foreign workers

The entry script holds a lock on the campaign directory for the whole
call, so two controllers cannot both create the manifest, both launch, or
both package (pass 3). Workers from an older launcher hold no unit lock,
so the queue cannot see them: the entry script refuses to launch OR
package while any such process is alive. A queue's unit set is fixed when
it is created; a different grid is a different directory. Existing
artifacts are validated against their unit (column, phase, lattice,
certification flags, problem identity) before the queue may call them
done; an invalid one is quarantined and the column walked again.

### 18. Concurrency is tested with processes, not with one session

Sequential tests prove ordinary-case behaviour; they cannot reach a race.
`test_campaign_processes` launches real MATLAB workers through the real
launcher and checks: a second process is still refused a lock after the
owner clears its functions and probes itself; units are computed once and
a unit whose owner is killed after its temporary write is published by
the REPLACEMENT; an always-failing unit retires after exactly three
attempts with its partial results kept as evidence; a kill on the LAST
attempt retires the unit; the SUPERVISOR kills a worker that goes silent,
records the reason and exit code, and the unit is abandoned and claimable.
The supervisor is the worker's PARENT: it reaps the exit and no recycled
pid can be mistaken for it.

## The shape of a campaign

    declare the problem ->  orbits, engine, phases, in ONE visible place
    drain foreign workers -> nothing outside the queue may be walking
    open the queue      ->  units and their artifacts, fixed at creation
    launch N workers    ->  each must report READY (queue opened); unique launch
    own by LOCK         ->  a live owner is never stolen from; the registry decides
    beat per STAGE      ->  silence beyond three stage caps = the parent supervisor kills
    publish under lock  ->  validated, exclusive temp, one rename
    blockers, not aborts ->  package only when complete and unheld
