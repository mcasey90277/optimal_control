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

### 12. Beat INSIDE the unit, and size the hang deadline from the solver, not the column

The first generated worker accepted a heartbeat callback and never called
it: claims were refreshed only between columns, columns take hours. The
unit function must call the beat at its natural inner cadence (here, after
every solve, whose wall cap is 900 s). The supervisor's inactivity deadline
is a multiple of THAT cap (hangSec = 2700 s), not of a measured column
duration -- a 9-hour column would otherwise have earned a 36-hour deadline.
There is no calibration solve: it ran outside the queue's ownership and
could duplicate a live worker's column.

### 13. Re-opening a campaign must not touch what workers own

"Re-run to package" is the normal path. If the entry point re-initialised
the queue, it freed columns being walked and reset the attempt counts that
stop the livelock. `open` is READ-ONLY and validates the unit set; `reset`
is a separate, explicit action that cannot touch a held unit.

### 14. Ownership is a process-held lock; heartbeat age is an alarm, never an authority

The first queue transferred a claim whose beat was 30 minutes old. Two
workers could both "reclaim" it (rename-takeover is an ABA race), and a
slow-but-alive owner would come back and publish over its replacement.
Tokens do not close that: a check followed by an action is the same stale
observation. Now a unit is owned by a kernel file lock (unit_lock) that the
process holds until it releases or DIES. Nobody can take it from a live
owner; a hung owner is killed by the supervisor, the kernel frees the lock,
and only then is the unit claimable. Measured on this host: the lock is
refused across processes, released on exit, released on SIGKILL.

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

### 17. One queue per campaign; drain foreign workers before it starts

Workers from an older launcher hold no lock, so the queue cannot see them
and would hand out a column one of them is walking. The entry script
refuses to launch while any such process is alive. A queue's unit set is
fixed when it is created; a different grid is a different directory.

### 18. Concurrency is tested with processes, not with one session

Sequential tests prove ordinary-case behaviour; they cannot reach a race.
`test_campaign_processes` launches three real MATLAB workers through the
real launcher, kills one mid-unit, and checks that every unit was computed
exactly once, the killed unit was taken over, the always-failing unit
retired after exactly three attempts, and the queue finished with no lock
held.

## The shape of a campaign

    declare the problem ->  orbits, engine, phases, in ONE visible place
    drain foreign workers -> nothing outside the queue may be walking
    open the queue      ->  units and their artifacts, fixed at creation
    launch N workers    ->  each must report READY (queue opened); unique launch
    own by LOCK         ->  a live owner is never stolen from
    beat per solve      ->  silence beyond the solver's cap = supervisor kills
    publish under lock  ->  validated, exclusive temp, one rename
    blockers, not aborts ->  package only when complete and unheld
