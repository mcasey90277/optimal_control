# Third-pass review

**Verdict: not yet fit for the advertised unattended, end-to-end multi-day campaign.** The replacement is a substantial improvement: the age-based ownership-transfer race is gone, ordinary queue opening is read-only, and the worker now controls publication. Those are real fixes.

However, four findings prevent approval:

1. **The claim’s authority is still just a copied Boolean.** A released claim can publish, beat, or delete a subsequent owner’s metadata.
2. **The lock registry can itself cause loss of a live lock.** It is neither protected against clearing nor keyed by file identity, and its exception path closes the very second descriptor that the design must avoid closing.
3. **The 2700-second healthy-silence premise is false.** One `certify_root` can traverse **4500 seconds of separately capped calls**, plus uncapped work, without calling `progress()`.
4. **Reset still has a check-then-act race.** It probes and releases the lock, then deletes attempt and owner records without exclusion.

There are also important controller, identity-validation, and finalization defects.

This is a static review of the supplied source. I have not executed MATLAB, traced this host’s JVM, or independently reproduced the reported measurements. I accept those measurements as evidence for the specific operations tested—not for stronger guarantees. Some significant dependencies remain absent, including `capped_pool`, `package_phase_catalog`, `guard_catalog_overwrite`, `audit_phase_catalog`, `second_order_pass`, and the full endpoint/sheet-building dependencies. Their behavior cannot be inferred from their names or comments.

All line references below refer to the **current** source unless explicitly stated otherwise.

---

# 1. Pass-2 findings, one by one

“FIXED” below means the particular earlier defect is fixed. It does not imply that every surrounding operation is now safe.

## P1-01 through P1-41

### Entry point, staging, and resumability

- **[CORRECTNESS]** `rib_validate.m:37–44; work_queue.m:134,209–210` — **P1-01: PARTIALLY FIXED.** Worker-side validation and short-column reporting exist, but existing files still mean queue completion, and the validator does not establish certification, unit identity, or requested coordinates. Validate against the expected unit contract before accepting either new or imported artifacts.

- **[CORRECTNESS]** `run_costate_library.m:179–199; rib_validate.m:43–44` — **P1-02: PARTIALLY FIXED.** A manifest exists, but code changes are merely noted, certification policy/input provenance are absent, and rib identities are not compared. Bind artifacts to a validated campaign identity; require explicit legacy import.

- **[CORRECTNESS]** `work_queue.m:98–108` — **P1-03: FIXED.** Ordinary `open` is read-only and checks units and output paths. Concurrent creation remains a separate defect.

- **[CORRECTNESS]** `work_queue.m:98–123` — **P1-04: FIXED.** Ordinary reopening no longer resets attempts. Reset’s exclusion defect remains.

- **[CORRECTNESS]** `run_costate_library.m:337–375` — **P1-05: PARTIALLY FIXED.** Packaging now has a barrier even with `.run.ribs=false`, and it checks `nHeld`. The invalid-artifact reporting branch can throw, foreign writers are not excluded here, and there is no finalization exclusion. Repair those paths; do not describe a probe sweep as a transactional quiescence barrier.

- **[GAP]** `run_costate_library.m:315–322,400–419; fine_library_autochain.sh:16–31` — **P1-06: PARTIALLY FIXED.** Chain exceptions become `failed`; general queue campaigns still require a later human call to finalize. The separate autochain is hardcoded to the old 24×24 run. Provide a campaign-specific finalizer/controller or narrow the “one call, end to end” contract.

- **[GAP]** `run_costate_library.m:384–390; build_70mN_library.m:95–126,201–208` — **P1-07: FIXED for the reported 24×24 handoff.** Exact grid, sheet, and rib paths now reach the chain. Non-default operating-point handoff still needs repair, discussed below.

- **[GAP]** `run_costate_library.m:391,411–419` — **P1-08: PARTIALLY FIXED.** An mtime threshold is better than unrestricted existence, but it is not a production receipt. Return a structured finalizer result tied to this invocation and these inputs.

- **[CORRECTNESS]** `run_costate_library.m:246–265` — **P1-09: FIXED for silent coverage omission.** `.onlyA` exclusions and missing seeds are now explicit blockers. A structured coverage result would be useful hardening.

- **[CORRECTNESS]** `run_costate_library.m:270–273` — **P1-10: FIXED.** Unowned calibration is removed.

### Timing and worker behavior

- **[CORRECTNESS]** `rib_from_crossing.m:88–90; certify_root.m:180–400` — **P1-11: PARTIALLY FIXED.** The callback reaches the walker, but its interval is not bounded by 900 seconds. Add progress at bounded stage boundaries and fix cancellation semantics.

- **[CORRECTNESS]** `run_campaign_workers.sh:71–74` — **P1-12: FIXED.** The watchdog measures heartbeat age, not worker lifetime.

- **[GAP]** `run_costate_library.m:270–273` — **P1-13: WITHDRAWN.** Calibration/fallback timing is no longer part of the safety mechanism. The replacement hang deadline has a different, demonstrable defect.

- **[GAP]** `run_costate_library.m:125,276,510–515; certify_root.m:146–149` — **P1-14: PARTIALLY FIXED.** Retry, stall, and hang settings are separated, but the supposed solve deadline is not an end-to-end deadline. Persist and validate the actual execution/certification policy, including stage caps.

- **[GAP]** `rib_from_crossing.m:70–111; build_ribs.m:92–94` — **P1-15: NOT FIXED.** No intra-column checkpoint. This remains an explicitly acceptable loss-budget choice once ownership and supervision are safe.

- **[GAP]** `campaign_worker.m:72,87–95; run_costate_library.m:511–515` — **P1-16: FIXED for premature idle departure.** Idle waiting now defaults to infinity, and the generated job does not override it. Replacement after losing all workers remains unsolved.

- **[CORRECTNESS]** `campaign_worker.m:104–109; work_queue.m:163–183` — **P1-17: FIXED.** Normal return without a publishable temporary artifact is no longer success.

- **[CORRECTNESS]** `campaign_worker.m:59–65,129–130` — **P1-18: FIXED for ordinary MATLAB lifecycle exceptions.** Failure telemetry is attempted and the error is rethrown. Resource cleanup on exceptional exit remains deficient.

- **[CORRECTNESS]** `campaign_worker.m:104–121; build_ribs.m:86–97; safe_report.m:38–41` — **P1-19: PARTIALLY FIXED.** Worker accounting follows publication, and the unsupported “output exists” reporting assertion is removed. Builder reporting can still throw after producing `tmpOut`, turning usable work into a failed attempt. Separate computation/commit from reporting throughout.

### Queue atomicity and ownership

- **[GAP]** `unit_lock.m:59–70` — **P1-20: WITHDRAWN.** Exclusive `mkdir` is no longer the ownership primitive. The replacement’s lifecycle must instead be tested.

- **[CORRECTNESS]** `work_queue.m:134–141` — **P1-21: FIXED.** There is no stale-path takeover or age-based transfer of ownership.

- **[CORRECTNESS]** `work_queue.m:158–191,252–256` — **P1-22: REGRESSED.** Sequential stale-claim rejection is now worse: `holds(c)` checks only the copied `.held` Boolean. After release, the caller’s old claim still passes. Use an active, identity-checked handle whose released state is shared by every reference.

- **[CORRECTNESS]** `work_queue.m:136–145` — **P1-23: FIXED.** Eligibility is rechecked after acquiring the unit lock. The lock-wrapper defects below must still be corrected.

- **[CORRECTNESS]** `work_queue.m:269–275,322–331,116–121` — **P1-24: PARTIALLY FIXED.** Text write/close results are checked and publication is atomic. Missing/non-file attempt paths still become zero, and reset can delete a newly written attempt under a live owner. Fix reset and distinguish absent records from malformed storage objects.

- **[CORRECTNESS]** `work_queue.m:205–230` — **P1-25: FIXED for exhausted-owner limbo and hidden ownership counts.** A free final-attempt unit is retired; `held` is separate from `done`. Status remains a sampled observation, not an atomic snapshot.

- **[GAP]** `run_costate_library.m:190–193; work_queue.m:128,196–197; campaign_status.m:33–36` — **P1-26: PARTIALLY FIXED.** The entry point adopts persisted policy, but the queue does not enforce it and monitors still default independently. Store policy in the queue definition and reject conflicting callers.

### Heartbeats and monitoring

- **[CORRECTNESS]** `campaign_heartbeat.m:97–114` — **P1-27: PARTIALLY FIXED.** Truncated records are rejected. The time is regex-shaped, not parsed; PID need only be non-NaN; campaign/instance identity is not validated. Validate a versioned structured record.

- **[GAP]** `campaign_heartbeat.m:82,113–114; campaign_status.m:44–46,76–77` — **P1-28: PARTIALLY FIXED.** “No check-in observed” improved, but status still prints “NEVER STARTED” and treats recent telemetry as live-process evidence. Separate telemetry freshness, observed process state, and progress.

- **[CORRECTNESS]** `run_campaign_workers.sh:42–54; campaign_heartbeat.m:66` — **P1-29: PARTIALLY FIXED.** Launch/tag collisions are addressed by `mktemp`; records still lack explicit queue/campaign identity and expected-instance validation. Include these fields.

- **[GAP]** `campaign_status.m:64–75` — **P1-30: PARTIALLY FIXED.** There is now an owner-tag join. It excludes done-but-held units and lacks authoritative expected-worker/PID identity. Reconcile all held units against a durable instance registry.

- **[GAP]** `run_campaign_workers.sh:59–78; campaign_worker.m:87–95; run_costate_library.m:315–316` — **P1-31: PARTIALLY FIXED.** Persistent idle recovery capacity and per-process watchdogs are improvements. There is still no controller that owns/reaps children, restores capacity, or finalizes general campaigns.

### Generated jobs and launcher

- **[CORRECTNESS]** `run_costate_library.m:504–506` — **P1-32: PARTIALLY FIXED.** Code roots are derived correctly; bootstrap remains hardcoded and unversioned. On this fixed host that is primarily deployment hardening: preflight and record the dependency explicitly.

- **[CORRECTNESS]** `run_costate_library.m:124,446–449` — **P1-33: FIXED.** Driver-derived output paths are canonical absolute paths.

- **[CORRECTNESS]** `run_campaign_workers.sh:51–52` — **P1-34: FIXED.** The job path crosses into MATLAB as environment data, not interpolated MATLAB syntax.

- **[CORRECTNESS]** `run_costate_library.m:289–290,516–520` — **P1-35: PARTIALLY FIXED.** Temporary job names are unique, but every invocation still overwrites one executable job pathname. Use immutable campaign/launch-specific jobs.

- **[GAP]** `campaign_worker.m:78–80; run_campaign_workers.sh:92–94` — **P1-36: PARTIALLY FIXED.** READY follows opening the queue and persists. The launcher accepts any `.ready` file before checking failure or exit, without validating its contents or queue identity. Parse READY and reconcile subsequent exit.

- **[GAP]** `run_campaign_workers.sh:35–42,49–54,81,87` — **P1-37: FIXED for zero-child loop success.** Bounded arithmetic loops, directory checks, and spawned-cardinality assertion replace `seq`.

- **[CORRECTNESS]** `run_campaign_workers.sh:42–54` — **P1-38: FIXED.** Concurrent launches receive separate directories, tags, and worker output files.

- **[CORRECTNESS]** `run_campaign_workers.sh:60,66,73` — **P1-39: PARTIALLY FIXED.** Escalation rechecks the command, but “contains `-batch`” is not instance identity. A recycled MATLAB PID can still be killed. Supervise owned children and reap them.

### Tests

- **[GAP]** `test_campaign_processes.m:68–129; test_work_queue.m:67–70,107–110` — **P1-40: PARTIALLY FIXED.** There is now a useful real-process test. It lacks controlled boundary races, and the stale-claim test manually falsifies the Boolean instead of testing an actually released claim.

- **[GAP]** `test_campaign_processes.m:49–65,84–129` — **P1-41: PARTIALLY FIXED.** Worker/launcher integration and one mid-unit kill are exercised. Real rib execution, supervisor killing, pool descendants, final-attempt death, simultaneous controllers, and crash-boundary publication remain untested.

## Grouped pass-1 robustness/style items

- **[CORRECTNESS]** `run_costate_library.m:399–435; build_70mN_library.m:304–310` — **Workspace isolation: PARTIALLY FIXED.** Base variables and cwd are protected, but root graphics defaults and the caller’s open figures are not. Remove `close all`; save/restore graphics defaults or use a separate process.

- **[ROBUSTNESS]** `campaign_worker.m:104–121; build_ribs.m:39–77` — **Shared-failure retry policy: NOT FIXED.** Deterministic setup/storage failures can consume every unit’s retry budget. Preflight once and distinguish infrastructure failure from unit failure.

- **[ROBUSTNESS]** `work_queue.m:136–150` — **Claim initialization recovery: FIXED for ownership/attempt reservation.** The lock spans initialization and attempts precede work. Missing owner telemetry still needs an explicit alarm state.

- **[ROBUSTNESS]** `work_queue.m:205–214; run_campaign_workers.sh:63–74` — **Clock/storage assumptions: PARTIALLY FIXED.** Age no longer authorizes takeover, but wall-clock jumps and host sleep can still trigger unnecessary kills. Disable sleep for production and use monotonic elapsed-time observation where feasible.

- **[ROBUSTNESS]** `run_campaign_workers.sh:41; run_costate_library.m:511–515` — **Heartbeat-directory agreement: NOT FIXED as an API guarantee.** The normal generated invocation agrees, but the launcher can still be supplied a job writing elsewhere. Bind READY to the expected configuration.

- **[ROBUSTNESS]** `run_campaign_workers.sh:64–75,86–103` — **Startup deadlines: PARTIALLY FIXED.** No-heartbeat startup has a deadline; heartbeat-without-READY does not. Enforce a separate READY deadline regardless of heartbeat activity.

- **[ROBUSTNESS]** `run_campaign_workers.sh:51–52,59–78` — **Detachment: PARTIALLY FIXED.** stdin/output redirection improved. The watchdog is still an unowned sibling process without an explicit survival/reaping contract. Put it under a durable parent/service.

- **[ROBUSTNESS]** `run_costate_library.m:117–156; run_campaign_workers.sh:36–37` — **Configuration/resources: PARTIALLY FIXED.** Launcher counts and supported lattice spacing are checked. Engine/policy domains, total campaign process count, pool children, threads, and memory are not. Validate and enforce aggregate limits.

- **[ROBUSTNESS]** `test_work_queue.m:153–160,185–186; test_campaign_processes.m:131–132` — **Tests failing CI: FIXED.** Tests throw; reporting exceptions now have a deliberate test.

- **[ROBUSTNESS]** `fmt_num.m:36–47` — **Fixed-width overflow: PARTIALLY FIXED.** Overflow rendering is implemented. Width/precision arguments remain unvalidated; validate them if this is intended as a general primitive.

- **[STYLE]** `run_costate_library.m:25–29,73–76; work_queue.m:24–29; run_capped.m:4–6` — **Guarantee language: NOT FIXED.** Calibration, hard-cap, and publication-authority claims exceed the implementation. Update the headers and reconcile the contradictory old/new rules in `CAMPAIGN_DISCIPLINE.md`.

## Pass-2 sections A–N

These status lines cover each separate finding in those sections, including findings that overlap the P1 list.

### A. Rename takeover

- **[CORRECTNESS]** `work_queue.m:134–145` — **A: FIXED.** The stale-directory ABA protocol is removed. Keep lock pathnames permanent and fix the new wrapper’s authority defects.

### B. Token checks and stale owners

- **[CORRECTNESS]** `work_queue.m:185–191,252–256` — **B, stale beat/release: REGRESSED.** Released claim copies still satisfy `holds`. Replace value-Boolean authority with an active claim handle.

- **[CORRECTNESS]** `campaign_worker.m:138–139; work_queue.m:167–179` — **B, loss of publication authority: NOT FIXED.** Lost authority is still not detected reliably; telemetry failure remains swallowed. Make invalid authority fatal to publication, not merely best-effort telemetry.

### C. Publication primitive

- **[GAP]** `publish_atomic.m:32–45` — **C: FIXED for same-directory atomic namespace publication on the stated host.** `ATOMIC_MOVE` fails rather than copying; existing-file replacement was measured. This is not a power-failure durability guarantee.

### D. Temporary names

- **[CORRECTNESS]** `work_queue.m:143–153; build_ribs.m:92–94; run_costate_library.m:516–520` — **D, shared `.part` paths: FIXED.** Attempts/builders/job writers no longer intentionally share a temporary pathname.

- **[ROBUSTNESS]** `work_queue.m:143,325,339` — **D, weak random suffixes: FIXED.** UUIDs replace the short custom PRNG suffixes. UUID naming is still not exclusive file creation.

### E. Final-attempt death

- **[CORRECTNESS]** `work_queue.m:211–230` — **E: FIXED in the state logic.** Final-attempt death yields retired once the lock is free. The claimed real-process test of this case does not exist.

### F. Completion and coverage

- **[CORRECTNESS]** `work_queue.m:205–230; run_costate_library.m:366–375` — **F, hidden ownership: PARTIALLY FIXED.** Ownership is independently exposed and checked. Probes are not a simultaneous quiescence snapshot; correctness must rest on immutable committed artifacts plus controller exclusion.

- **[CORRECTNESS]** `rib_validate.m:39–44; run_costate_library.m:346–364` — **F, short columns: PARTIALLY FIXED.** They are now named, but only counts and stop strings are inspected. Validate exact per-target coordinates and certification fields.

### G. Idle tail

- **[GAP]** `campaign_worker.m:72,87–95` — **G: FIXED for the finite idle-timeout defect.** Surviving idle workers remain available. Losing all workers still requires external recovery.

### H. Timing

- **[GAP]** `certify_root.m:180,210,233,258,276,304,348,357,400; rib_from_crossing.m:90` — **H, silence bound: NOT FIXED.** The supplied implementation positively disproves the 900-second premise. Instrument/cap stages and size the watchdog from the resulting inter-progress contract.

- **[CORRECTNESS]** `work_queue.m:136–141; run_campaign_workers.sh:71–74` — **H, lease expiry before termination: WITHDRAWN.** There is no lease expiry or live-owner takeover.

### I. Resume and reset

- **[CORRECTNESS]** `work_queue.m:89–108; run_costate_library.m:183–199` — **I, configuration mutation: PARTIALLY FIXED.** `open` is fixed; creation remains check-then-replace without controller exclusion. Serialize manifest/queue creation and validate after acquiring that exclusion.

- **[CORRECTNESS]** `work_queue.m:115–121` — **I, reset race: NOT FIXED.** Probe-then-delete is the same race in a different form. Hold the unit lock while resetting.

### J. Grid and identity

- **[CORRECTNESS]** `run_costate_library.m:153–156,205–207,234–236` — **J, discarded phase vectors: PARTIALLY FIXED.** Unsupported spacing is rejected, but requested departure origin is not passed when building a new sheet. Pass `sD(1)` or reject unsupported origins before expensive construction.

- **[CORRECTNESS]** `run_costate_library.m:214–236; build_ribs.m:61–68; rib_validate.m:43–44` — **J, incomplete identity: PARTIALLY FIXED.** Sheet identity is required, but only selected fields are checked; rib identity is not compared at all. Validate materialized identities and source provenance explicitly.

### K. Finalization

- **[CORRECTNESS]** `run_costate_library.m:384–390; build_70mN_library.m:95–126,201–208` — **K, wrong 12×12 input contract: FIXED for the current 24×24 campaign.**

- **[CORRECTNESS]** `run_costate_library.m:391,415–419` — **K, old catalog reported as newly packaged: PARTIALLY FIXED.** The mtime heuristic remains insufficient. Return an invocation-bound receipt.

- **[CORRECTNESS]** `run_costate_library.m:399–435; build_70mN_library.m:304–310` — **K, workspace destruction/escaping chain error: PARTIALLY FIXED.** Variables/cwd and chain exceptions are addressed; graphics destruction remains. Isolate process-global side effects too.

### L. Launcher

- **[CORRECTNESS]** `run_campaign_workers.sh:64–75` — **L, missing first heartbeat: PARTIALLY FIXED.** Truly absent heartbeat is timed out. A heartbeat without READY defeats the independent startup deadline. Track READY separately.

- **[GAP]** `run_campaign_workers.sh:92–94` — **L, failed worker accepted as launched: PARTIALLY FIXED.** Persistent READY proves a historical queue attachment, not current service or clean exit. Reconcile READY, failure, and exact child exit.

- **[CORRECTNESS]** `run_campaign_workers.sh:42–54,60–77; run_costate_library.m:294–308` — **L, identity/capacity: PARTIALLY FIXED.** Launch uniqueness is fixed; exact-child ownership and campaign-wide concurrency limits are not.

- **[GAP]** `run_campaign_workers.sh:35–42,49–54,81,87` — **L, positive count but no children: FIXED.**

### M. Interrupted non-rib saves

- **[ROBUSTNESS]** `build_arrival_sheet.m:104–108; run_costate_library.m:270–273; work_queue.m:328–331` — **M: FIXED for the identified sheet/calibration/attempt-write paths.** Sheet publication is atomic, calibration is removed, and attempt writes are checked. Downstream catalog/sidecar crash behavior is still not auditable from this bundle.

### N. Old launcher cutover

- **[CORRECTNESS]** `run_costate_library.m:298–305,343–375; fine_library_autochain.sh:17–28` — **N: PARTIALLY FIXED.** The normal launch path checks for foreign workers; the live autochain waits for them. Manual launcher invocation and direct packaging bypass that exclusion, and `pgrep` failure is not distinguished from no matches. Make cutover an explicit, enforced campaign operation.

## Pass-2 rulings on the pushbacks

- **[GAP]** `run_costate_library.m:179–199` — **Immutable identity: PARTIALLY FIXED.** A small manifest is appropriate, but this one neither freezes inputs/code nor validates imported ribs. Expand and enforce it; no distributed metadata service is needed.

- **[CORRECTNESS]** `unit_lock.m:53–86; work_queue.m:252–256` — **Declining generations: PARTIALLY FIXED.** Lifetime locks are a valid alternative, but this wrapper does not reliably preserve or verify that lifetime. Fix the handle/registry model rather than adding superficial token checks.

- **[GAP]** `work_queue.m:205–214` — **Declining a lease-clock coordinator: FIXED/ACCEPTED.** Age is now an alarm, not transfer authority.

- **[GAP]** `test_campaign_processes.m:68–129` — **Separate-process verification: PARTIALLY FIXED.** Useful evidence exists; the critical new boundaries are not tested.

- **[GAP]** `rib_from_crossing.m:70–111` — **Deferring checkpoints: NOT FIXED, CONDITIONALLY ACCEPTED.** State the whole-column loss budget explicitly.

- **[CORRECTNESS]** `run_costate_library.m:270–273` — **Calibration ownership: FIXED.** Removing calibration is the simplest correct resolution.

---

# 2. The lock model, adversarially

## 2.1 What primitive is this on macOS?

For the conventional OpenJDK macOS implementation, `FileChannel.tryLock()` reaches the native file dispatcher and uses **POSIX record locking through `fcntl(..., F_SETLK, ...)`**, with an exclusive write lock for this call. Blocking acquisition uses the corresponding blocking operation; unlocking uses `F_UNLCK`.

It is **not normally `flock(2)`**, and it is not Linux open-file-description locking.

However, **Java’s public API does not mandate that syscall**. The supplied three-process measurements distinguish useful behavior, but do not distinguish `fcntl` from every alternative.

- **[GAP]** `unit_lock.m:59–62` — **The installed JVM’s lock implementation has not been identified by the supplied evidence.** Record `version -java` and the JVM vendor/build, then inspect that implementation or trace one acquisition/release. The expected result is traditional POSIX `fcntl` semantics; do not claim syscall verification merely from exclusion tests.

For traditional POSIX record locks on this deployment:

| Event | Consequence |
|---|---|
| Another process opens/closes the same file | Does **not** release this process’s lock. |
| This process closes **any descriptor referring to the same file** | Can release **all this process’s record locks on that file**, even if the descriptor did not acquire them. |
| `FileLock.release()` | Releases the lock. |
| Closing the owning channel or `RandomAccessFile` | Releases the lock. |
| Process exit/SIGKILL | Releases its locks. |
| Sending SIGTERM | Does not itself prove exit or release. |
| Forking a child | Traditional process-associated record locks are not inherited as the child’s locks. |
| Unlinking/replacing the pathname | The lock remains associated with the old file object; a new pathname target can be locked independently. |
| Truncating the same file without closing a descriptor | Does not itself create a new lock object or revoke the record lock. |
| Changing permissions | Does not retrospectively revoke an already held lock. |

The “any descriptor” rule includes descriptors opened through MATLAB `fopen`, a MEX routine, Java, or another library **inside the same MATLAB process**. A Java-only registry cannot govern arbitrary native opens.

## 2.2 The immediate authority defect: a released claim still owns everything

- **[CORRECTNESS]** `work_queue.m:252–256,158–191; unit_lock.m:73–80` — **`holds(c)` is not a lock check. It is a stale value check.** `release` closes Java objects, but does not mutate the caller’s struct. Consequently:

  ```matlab
  c = work_queue('claim', q, 'A');
  work_queue('release', q, c);

  % c.lock.held is still true.
  % Write a valid artifact to c.tmpOut:
  P = work_queue('publish', q, c, c.tmpOut, validator);
  ```

  `P.ok` can be true, although A no longer holds the lock.

  If B has meanwhile claimed the unit:

  - A’s old `beat` overwrites B’s owner metadata.
  - A’s old `release` deletes B’s owner metadata.
  - A’s old `publish` overwrites the final result without ownership.

  **Concrete fix:** represent the claim as a mutable handle with an active/released state shared by all references. Validate the claim’s registered token, process identity, queue/unit/output identity, open channel, and `FileLock.isValid()`. Release must be idempotent and invalidate authority before exposing the unit for reuse.

`FileLock.isValid()` is **necessary but insufficient**: it cannot necessarily discover a kernel lock silently released by closing an unrelated descriptor. That failure must be prevented structurally.

The synchronous current rib walker does not intentionally call an old callback after release. That narrows the ordinary execution path; it does not make the advertised ownership API true. The regression is directly testable without racing anything.

## 2.3 Registry clearing and pathname aliases can drop the actual kernel lock

- **[CORRECTNESS]** `unit_lock.m:53–67` — **The registry is not durable for the lifetime of the locks it represents.** There is no `mlock`, and the map stores only Booleans.

  A deterministic failure sequence is:

  1. A claim holds Java channel X.
  2. `clear unit_lock`, an applicable `clear functions`, or `clear all` clears the persistent registry while another reference still retains X.
  3. A subsequent probe opens channel Y on the same file.
  4. Java rejects the overlapping lock in this JVM.
  5. The catch converts that to ordinary contention.
  6. `ch.close()` closes Y.
  7. Under the expected POSIX semantics, closing Y releases X’s kernel lock too.

  The original claim still says `.held=true`, and Java may still consider its `FileLock` valid.

  **Concrete fix:** use one process-wide, strongly retaining lock manager; protect the MATLAB entry point against ordinary clearing, and reject duplicate manager instances. An overlapping-lock exception must be treated as an internal invariant failure, not as routine contention followed by unsafe cleanup.

  Prevention is important here: after an accidental second descriptor has already been opened, simply “close it carefully” does not solve traditional `fcntl` semantics.

- **[CORRECTNESS]** `unit_lock.m:57–58,71,83; work_queue.m:259` — **Registry keys are raw path strings, not file identities.** Relative/absolute paths, `..`, symlinks, hard links, and case aliases on a case-insensitive APFS volume can refer to the same file with different keys. The same destructive second-channel sequence then occurs without clearing any function.

  The entry point’s canonical `outDir` helps the normal generated path. It does not make `unit_lock` or manual queue operations safe.

  **Concrete fix:** normalize paths at the locking API boundary, prohibit alternate aliases, and ideally key/check by the filesystem’s file identity. Canonical path normalization alone does not resolve hard-link aliases. Keep lock files in a dedicated namespace that no other code opens.

- **[CORRECTNESS]** `unit_lock.m:75–78` — **A repeated release can remove a newer registry entry.** An old `L` retains `.held=true`; after a new same-process acquisition, releasing old `L` again removes `heldHere(L.file)` without checking which lock generation it represents. The next probe can then open the destructive second channel.

  **Concrete fix:** store the actual active lock object/token in the registry and remove an entry only when it is that exact object. Make release idempotent through shared handle state.

A second function implementation, package-qualified copy, or separately loaded manager can have independent persistent state. Merely adding the same function directory twice does not automatically instantiate a second registry; **resolving operations through different implementations can**. A supported worker should use one immutable code path for its lifetime.

## 2.4 Object copying and garbage collection

A MATLAB struct copy containing Java objects copies references to those Java objects; it does not clone the underlying file descriptor or acquire another lock.

In the ordinary worker:

- `c` remains live through publication and release.
- `beat` also captures `c`.
- Therefore there is no basis for claiming that a reachable claim is spontaneously garbage-collected midway through `unitFcn`.

That part of the author’s design is reasonable.

- **[ROBUSTNESS]** `unit_lock.m:53–54,70–71; campaign_worker.m:86–121` — **The registry does not retain the Java objects or clean up abandoned handles.** If the last real reference to a claim disappears without release—for example during an exception unwind—resource cleanup/garbage collection can eventually close the descriptor while MATLAB remains alive. The Boolean map entry can remain forever, reporting a lock that no longer exists.

  Conversely, relying on GC for prompt release can leave a lock held unpredictably long.

  **Concrete fix:** strongly retain active lock handles in the manager and use an explicit per-claim cleanup guard. Remove registry entries and close resources together, exactly once.

Other JVM-related release paths include explicit channel closure, `RandomAccessFile` closure, and closure resulting from an interrupted operation on an interruptible channel. There is no evidence that MATLAB normally performs such an operation on this idle lock channel. The correct conclusion is **not** “Java randomly drops live locks”; it is that the wrapper must exclusively own and govern the channel lifecycle.

## 2.5 Lock-file deletion, rename, replacement

- **[CORRECTNESS]** `unit_lock.m:59; work_queue.m:259` — **The protocol requires permanent lock-file identity, but does not establish an enforced maintenance contract.** Unlinking `7.lock`, replacing it, restoring the queue directory from another copy, or moving/recreating the queue lets another process lock a different file while the old owner still locks the original inode.

  Renaming a locked file does not make the old lock vanish; it makes pathname-based coordination diverge.

  **Concrete fix:** never delete or replace lock files during a campaign, including during reset, cleanup, migration, or backup restoration. Create/protect the lock namespace under controller exclusion. Record/check file identity if the wrapper is expected to detect replacement.

The current `reset` fortunately does **not** delete `.lock` files. Keep that property.

## 2.6 Child pool process

- **[GAP]** `run_costate_library.m:506; build_ribs.m:42,81–83; run_capped.m:44–50` — **The client’s pool child is outside the supervisor’s lifecycle contract.**

  For the expected POSIX locking implementation:

  - The pool child does **not** inherit ownership of the client’s process-associated record lock.
  - Here the pool is also created before the unit claim.
  - The submitted functions receive numerical inputs, not `c`, `tmpOut`, or the publication callback.
  - Nothing in the supplied submission path authorizes the pool child to publish a rib.

  Therefore I do **not** identify a normal path by which that child retains the unit lock or publishes after client death.

  What is not established is whether it stops computing. A local pool worker or its job infrastructure may survive client SIGKILL, particularly while executing non-interruptible native code. It can continue consuming a core, memory, and licenses while a replacement client starts another pool.

  **Concrete fix:** record the actual pool/job/process relationship, test client SIGKILL during a running task, and make supervisor cleanup cover owned descendants/jobs—not a broad name-based `pkill`.

This is worse than a cosmetic orphan: repeated replacement can degrade the host until otherwise healthy workers miss their own deadlines.

## 2.7 Is publication actually protected?

**On the intended happy path, yes:** acquisition precedes work; the parent produces/validates the attempt output; publication occurs before release. There is no deliberate release-before-rename interval.

**As implemented, the advertised guarantee is nevertheless false**, because:

- `holds(c)` accepts released claims.
- Registry failure can release the kernel lock while `holds(c)` remains true.
- Lock pathname replacement creates two lock domains.
- Foreign/direct builders do not participate in the protocol.

Atomic rename solves publication visibility. It does not establish publication authority.

---

# 3. Defects in and exposed by the rewrite

## 3.1 Reset is still a race

- **[CORRECTNESS]** `work_queue.m:115–121; unit_lock.m:84–86` — **Reset acquires a lock only long enough to release it before doing the protected operation.**

  Interleaving:

  1. Reset probes unit 8: free.
  2. Probe releases its temporary lock.
  3. Worker acquires unit 8, increments `.att`, writes `.owner`.
  4. Reset deletes both files.
  5. Worker continues holding the unit, but its reserved attempt has disappeared.

  Ownership of the computation is not stolen, but the retry budget is reset under a live owner and monitoring loses its owner record.

  **Concrete fix:** use `unit_lock('try', ...)` in reset and retain the acquired lock through all mutations. Validate the requested unit IDs against the queue. Release with cleanup protection. A probe is not exclusion.

## 3.2 Manifest and queue creation are not exclusive

- **[CORRECTNESS]** `run_costate_library.m:183–199; work_queue.m:89–95` — **Two first callers can both create incompatible “immutable” campaign definitions.** Both observe absence, both write unique temporary files, and both atomically replace the final definition. Atomic last-writer-wins is not create-once.

  This can occur with different `.onlyA` selections even when orbit/grid requests agree: `.onlyA` is not part of the manifest, while it controls queue units.

  **Concrete fix:** acquire a permanent campaign/controller lock before checking or creating the manifest and queue. Recheck after acquisition. Publish once, then validate the existing definition. Serialize sheet creation and finalization through the same controller contract.

- **[CORRECTNESS]** `run_costate_library.m:289–290,294–308,380–415` — **Concurrent/repeated entry calls are not controlled at the campaign level.** They can overwrite the shared job, launch additional batches, or run two packagers/sweeps against the same catalog and sidecar.

  Per-unit locks prevent cooperating workers from computing the same unfinished unit. They do not limit total workers or serialize catalog writers.

  **Concrete fix:** maintain a campaign instance registry and reconcile to a requested total worker count. Use immutable jobs, and permit only one finalizer for a catalog generation.

For this deployment, a single local controller lock and a modest registry are sufficient. No distributed coordinator is required.

## 3.3 The manifest silently adopts too much

- **[CORRECTNESS]** `run_costate_library.m:179–199,211–236; rib_validate.m:43–44` — **Creating a manifest in a populated directory does not validate the artifacts it adopts.**

  The current manifest:

  - compares orbit/engine request structs and phase vectors;
  - adopts previous scheduling policy after printing a note;
  - merely prints a note about code revision changes;
  - records only a Git HEAD hash, not dirty working-tree content;
  - omits sheet/anchor/arc hashes and materialized endpoint identity;
  - omits certification tolerances, gate configuration, stage caps, artifact schema, and selected units;
  - does not bind workers to the manifest;
  - does not compare any rib’s `problem` to the campaign.

  A same-named legacy rib from another grid or certification policy can be treated as done, and the new manifest gives that adoption a misleading appearance of provenance.

  **Concrete fix:** require an explicit legacy-import operation that validates each artifact against a frozen unit contract. Freeze the relevant code/dependencies or reject unapproved changes. Record resolved policy—not just scheduling knobs—and identify the exact sheet used by workers.

The note about scheduling-policy adoption is not silent, and using persisted scheduling policy is reasonable. What is missing is enforcement outside this entry-point invocation.

## 3.4 `rib_validate` does not validate a certified rib

- **[CORRECTNESS]** `rib_validate.m:37–44` — **The validator accepts an identity variable’s name, not its value, and counts arbitrary `.pts` elements as certified points.** It does not establish:

  - scalar/expected rib structure;
  - `R.j` equals the claimed column;
  - `R.sA` equals that column’s phase;
  - exactly one intended rib is present;
  - each point has `ok=true` and the expected certification/schema markers;
  - points are unique and lie on the requested departure grid;
  - the `problem` variable is a valid identity, let alone the expected one;
  - per-rib coverage, rather than a total count across several ribs.

  For example, `.pts` can contain unrelated or uncertified data while a variable named `problem` contains an empty value.

  **Concrete fix:** change the validator contract to accept the expected campaign/unit specification. Validate structure, identity, coordinates, and certification flags. Return completed/unresolved target sets, not just a count.

- **[ROBUSTNESS]** `rib_validate.m:32–40` — **Only `whos` is caught; subsequent load/schema errors escape.** The worker catches validator exceptions, but the entry barrier calls the validator without such protection. A corrupt MAT body or missing `.stop` can abort the entry point rather than return a named blocker.

  **Concrete fix:** encompass load and structural validation in the validator’s error-to-result boundary.

## 3.5 The new barrier has an indexing bug

- **[CORRECTNESS]** `run_costate_library.m:344–356` — **`short` and `good` do not share an indexing domain.** `good` has one element per selected column; `short` contains only short or invalid columns. This expression is invalid:

  ```matlab
  short(~good(ismember(cols, cols)))
  ```

  `ismember(cols, cols)` is all true, so it reduces to indexing the shorter `short` cell array by `~good`.

  Example: three columns, first two complete, third invalid. `short` has one element; the logical index has its true element at position three. The intended blocker path throws.

  **Concrete fix:** maintain a per-column `reason` cell array parallel to `good`, plus a separate `isShort` mask. Use `reason(~good)` and `reason(isShort)`.

This is a small deterministic defect, not a concurrency edge case.

- **[GAP]** `work_queue.m:134,209–210; run_costate_library.m:353–358` — **An invalid existing output is permanently “done” to the queue but “pending” to the packager.** Workers skip it forever, while repeated entry calls do not repair it.

  **Concrete fix:** validate imported/existing outputs before queue activation. Under unit/controller exclusion, quarantine invalid outputs and define whether repair consumes a normal attempt. Return a blocked-invalid-artifact state rather than an apparently progressing pending state.

## 3.6 The worker commit boundary: what actually happens

The supplied code does **not** implement “publish, then release, then accounting.” It implements:

```text
unit function
publish
account + log
release
```

at `campaign_worker.m:105–121`.

That ordering is not intrinsically wrong. In particular:

| Crash point | Appropriate interpretation |
|---|---|
| Before final rename | Attempt spent; no committed result; eligible for retry/retirement. |
| After rename, before local accounting | Unit is committed; worker-local counters may be lost. |
| After accounting, before release | Unit is committed and temporarily still held. |
| During owner-record cleanup | Unit remains committed; cleanup/worker failure must not reclassify the artifact. |

The output-as-commit design handles the first two important cases well, subject to real ownership and validation.

- **[ROBUSTNESS]** `campaign_worker.m:99–121; work_queue.m:185–191` — **There is no exception-safe claim guard, and cleanup can throw before releasing the kernel lock.** Initial unit heartbeat, metadata deletion, or temporary deletion can throw outside the per-unit try. The outer lifecycle catch writes `fail` and rethrows, but does not explicitly release the acquired lock.

  A `-batch` process will ordinarily exit and free it eventually. An interactive caller catching the error can instead retain a lock or stale registry until GC/exit.

  Logging is best-effort against exceptions, but it can still block while the lock is held.

  **Concrete fix:** install per-claim cleanup immediately after acquisition. Make release close/invalidate the lock in a `finally`-equivalent path regardless of metadata cleanup errors. Report commit status separately from cleanup status.

- **[ROBUSTNESS]** `build_ribs.m:92–94; work_queue.m:189; campaign_worker.m:117–118` — **Killed writes can leave uncollected attempt debris.** The builder writes an additional nested temporary pathname before publishing to `c.tmpOut`; neither release nor recovery knows that nested pathname. SIGKILL during either save leaves debris. `.failed` movement also uses `movefile` without checking its returned status.

  **Concrete fix:** use an attempt directory or registered attempt-file set and garbage-collect only files belonging to confirmed inactive attempts. Check evidence-preservation operations. Do not claim that arbitrary kill points leave no `.part`.

## 3.7 Status probes: correct scope, actual costs, actual side effects

- **[GAP]** `work_queue.m:202–230; unit_lock.m:82–86` — **Status is not a read-only snapshot.** It opens/possibly creates every lock file, briefly acquires free locks, then releases them. Another observer can momentarily report the probing process as an owner; workers can skip a unit while a monitor briefly owns its lock. Attempts, lock state, owner metadata, and output existence are sampled at different times.

  Thus “owner died,” “claimable now,” and “nothing will change” are stronger than the observations establish.

  **Concrete fix:** label status as sampled/advisory, distinguish unknown owner metadata, and never use a released probe as authorization for mutation.

- **[ROBUSTNESS]** `work_queue.m:204–206; campaign_worker.m:88–94` — **Every idle worker performs a full lock-probe sweep every 15 seconds.** At 24 units and four idle workers, that is approximately **553,000 unit probes per day**, plus metadata reads and heartbeat writes.

  On local APFS with only 19–24 units, this is not a performance blocker. It is unnecessary metadata traffic and an amplifier for the same-process alias/registry bugs.

  **Concrete fix:** first make probing safe. Then centralize monitoring or use modest jitter/backoff and avoid repeated probes of immutable committed units where the application does not need an ownership census.

There is an important distinction here: once every output has been correctly committed and no normal queue worker can replace a committed output, packaging can safely consume those immutable outputs without a perfect simultaneous “zero locks” snapshot. Late claimers can briefly acquire a lock and fail their done recheck without modifying data.

The current code should adopt that precise contract—or acquire real finalization exclusion—not claim that sequential probes prove global quiescence.

## 3.8 Supervisor identity and lifecycle are still inadequate

- **[CORRECTNESS]** `run_campaign_workers.sh:60,66,73` — **`isours` does not mean ours, or even MATLAB.** It means that the current command string for a PID contains `-batch`. Another MATLAB batch process—or another command containing that string—passes after PID reuse.

  Rechecking before SIGKILL reduces exposure but does not close it. There is also a check-to-signal interval.

  **Concrete fix:** make the supervisor the worker’s actual parent, retain ownership until `wait`/reaping, and persist the exit result. Keep descendants in a deliberately owned process/job group. Command/start-time checks are useful diagnostics or supplementary guards, not a replacement for child ownership.

- **[GAP]** `run_campaign_workers.sh:59–78` — **These supervisors are worker siblings, not parents, and do not collect exit status.** When `isours` becomes false, they log “no longer running.” They cannot distinguish clean completion, MATLAB error, signal death, failed `ps`, or a process whose command representation changed.

  After requesting termination they break without confirming exit. If kill fails or the process remains alive, idle peers can wait forever on its lock.

  **Concrete fix:** use an owned parent supervisor that waits/reaps, records exit cause, confirms termination, and escalates unresolved termination as a campaign failure.

### What happens on normal worker exit?

If the worker really exits and its PID is not reused, the next `isours` check ends the supervisor, possibly after one final 15-second sleep. **The code does not inherently kill a normally exited worker.**

But it also does not obtain its exit code. A `done` heartbeat is written before MATLAB shutdown; shutdown can subsequently hang or fail. The supervisor does not parse that distinction and will eventually apply the inactivity timeout if the process remains.

That is a legitimate shutdown watchdog only if it is documented and reported as such.

- **[GAP]** `run_campaign_workers.sh:64–75,92–103` — **The startup deadline is independent of heartbeats only while no heartbeat exists.** A job that writes `.hb` but never `.ready` escapes `STARTUP_SEC`; an actively updated heartbeat can keep it alive indefinitely. The verification message promising that the supervisor will kill it at `STARTUP_SEC` is false.

  **Concrete fix:** track an explicit startup state. Until validated READY, apply the spawn-to-READY deadline regardless of heartbeat files. After READY, apply progress inactivity.

- **[GAP]** `run_campaign_workers.sh:92–94; campaign_heartbeat.m:58–68` — **READY is only existence-checked, and wins over an already written failure.** A worker that attaches and immediately fails can still produce `LAUNCH OK`.

  **Concrete fix:** validate READY’s schema, PID/instance, queue identity, and configuration digest. Report “attached then failed/exited” separately from “serving” and “cleanly completed.”

- **[ROBUSTNESS]** `run_campaign_workers.sh:60–78` — **A single process-inspection/storage failure can disable supervision or produce a false timeout.** `ps`, `stat`, arithmetic, and signal outcomes are not robustly handled; supervisor stderr is discarded.

  **Concrete fix:** make observation failure an explicit supervisor alarm, retry transient failures, and preserve supervisor diagnostics. Do not equate inability to inspect a process with confirmed exit.

## 3.9 Finalization still has side effects and incomplete contracts

- **[CORRECTNESS]** `build_70mN_library.m:304–310` — **The local-function wrapper does not isolate process-global graphics state.** Default `pictures=true` changes `DefaultFigureVisible` and calls `close all`, closing figures belonging to the caller. Exceptions can leave the default changed too.

  **Concrete fix:** save/restore defaults with a cleanup guard outside the script’s `clearvars`, and close only figure handles created by the chain. A separate finalizer MATLAB process is an even cleaner boundary.

- **[CORRECTNESS]** `run_costate_library.m:384–390; build_70mN_library.m:60–63,128–131,233–235,305–306` — **Non-default operating points still meet default finalizer configuration.** The entry accepts arbitrary declared engine/orbits, but passes neither to the chain. The chain prints default values, supplies default engine options to packaging, and constructs a default-orbit picture pathname.

  The missing `package_phase_catalog` implementation prevents determining which defaults it overrides from the sheet. The handoff itself is contradictory.

  **Concrete fix:** pass the resolved operating-point identity, or have the finalizer derive and validate it entirely from the accepted sheet/manifest. Do not retain two competing parameter sources.

- **[GAP]** `build_70mN_library.m:141–147,175` — **Package-only finalization still requires default anchors, an arc file, and a pool.** These are checked even though the supplied exact sheet and ribs are the material inputs and arc/sheet/rib execution is disabled.

  **Concrete fix:** gate prerequisites by the stages that actually need them. A validated catalog packaging operation should not fail because an unrelated default anchor is absent or pool startup fails.

- **[GAP]** `build_70mN_library.m:255–262,280–294; run_costate_library.m:400–406` — **“A stage that throws does not cost later measurements” remains false for audit.** An audit exception escapes the script and prevents sweep/pictures, although a sweep exception is contained.

  **Concrete fix:** give audit the same stage-outcome boundary, and return structured per-stage results. Distinguish a bad audit row from an audit infrastructure exception.

- **[CORRECTNESS]** `run_costate_library.m:391,415–419` — **Catalog mtime cannot prove this invocation produced this campaign’s catalog.** A recent pre-existing catalog, a concurrent writer, or later sweep writeback can satisfy the timestamp check. The two-second subtraction explicitly admits a pre-call interval.

  **Concrete fix:** have packaging return a receipt containing invocation ID, manifest/input digest, output identity, and produced/reused status. Publish the receipt with the corresponding catalog generation.

Audit-only/sweep-only successful calls also still return an unhelpful default `pending` state and blank catalog path. That is an outcome-model gap, not evidence that those stages did no work.

## 3.10 The live finish job/autochain

- **[CORRECTNESS]** `run_costate_library.m:298–305,343–375; fine_library_autochain.sh:18–28` — **Foreign-worker exclusion is an operational observation, not a coordinated cutover.** The entry checks it only before launching new queue workers, not before packaging. A direct launcher call bypasses it. The autochain’s “no old workers” check can race a later old-worker launch.

  Also, `pgrep` errors are collapsed into no matches.

  **Concrete fix:** administratively disable the old launcher, drain and verify its exact instances, record cutover, and enforce that state before queue activation/finalization. Treat process-query errors as unknown/blocking.

- **[ROBUSTNESS]** `fine_library_autochain.sh:13–22` — **The preliminary barrier counts glob matches, not the 19 expected validated units.** Unexpected matching files can satisfy the count; missing/corrupt expected files then fail later.

  The entry’s validator should prevent packaging those missing units, so this is not independently proof of wrong packaging. It is a poor launch trigger and combines badly with the barrier indexing bug.

  **Concrete fix:** derive expected unit paths from the sheet/manifest and check those exact paths. Let the finish job issue the authoritative validation outcome.

- **[GAP]** `fine_library_autochain.sh:14,28–31; fine_library_finish_job.m:20–41` — **Finalization has no single-instance ownership, robust supervision, or reliable shell failure result.** Two autochains can run concurrent finalizers and overwrite logs/verdicts. A hung finalizer blocks the shell indefinitely. MATLAB exceptions are caught into text, so the MATLAB process can exit zero after packaging or movie failure. The shell’s final command can also leave an unhelpful exit status.

  **Concrete fix:** give finalization its own controller lock, progress/deadline supervision, unique logs, atomic structured verdict, and explicit exit status. Preserve existing verdicts as history rather than deleting them at arming time.

Rendering a movie despite coverage/audit blockers is acceptable **as a diagnostic action**. It must not be mistaken for a shipping approval. `state='packaged'` is not a clean-library verdict.

---

# 4. The `hangSec` premise

## 4.1 The nominal capped work between beats is 4500 seconds, not 900

`progress()` is called only after `certify_root` returns:

```matlab
Ct = certify_root(...);
...
try progress(); catch, end
```

at `rib_from_crossing.m:88–90`.

On a successful full-stack certification, these separately capped calls can all execute:

| Stage | Call site | Default cap |
|---|---:|---:|
| Polish plus conjugate test | `certify_root.m:180` | 900 s |
| Accepted-solution flight | `certify_root.m:210` | 300 s |
| Foreign witness solve | `certify_root.m:258` | 300 s |
| Witness flight | `certify_root.m:276` | 300 s |
| First hypothesis-gates build | `certify_root.m:304` | 900 s |
| Second hypothesis-gates build | `certify_root.m:348` | 900 s |
| Dense conjugate scan | `certify_root.m:400` | 900 s |
| **Total** | | **4500 s = 75 min** |

A cap failure returns early, but a healthy call can finish each stage just below its cap and proceed to the next. There is no 900-second aggregate budget.

The `.wallSec=900` passed by the generated job reaches the in-process polish option at `certify_root.m:180–181`. It is not an end-to-end `certify_root` cap.

- **[CORRECTNESS]** `run_costate_library.m:73–76,510; certify_root.m:146–149,180–400` — **A healthy certification can outlast `hangSec=2700` without a heartbeat.** The supervisor can destroy a productive column after tens of minutes of legitimate multi-stage work.

  **Concrete fix:** either:
  
  1. beat at meaningful bounded stage boundaries within `certify_root`, while supervising each stage correctly; or
  2. enforce an actual aggregate certification deadline and choose a longer watchdog with measured overhead margin.

  Simply multiplying the polish cap by three is not a valid derivation.

## 4.2 Work outside the capped calls

Between progress callbacks, the following also occurs outside those per-call waits:

- `build_ribs` loading and pool acquisition: `build_ribs.m:39–42`;
- setup, endpoint generation, interpolant construction, anchor loading: `build_ribs.m:56`; `arclength_arrival.m:71–101`;
- endpoint evaluation and seed construction: `rib_from_crossing.m:71,83–86`;
- accepted and witness flight validation: `certify_root.m:218,280`;
- pointwise PMP checks: `certify_root.m:233`;
- `lift_margin`: `certify_root.m:357`;
- validation/bookkeeping and returned-data handling;
- logging;
- final save, validation, and publication: `build_ribs.m:86–97`; `campaign_worker.m:106`;
- `parfeval` submission before the timed wait and output retrieval afterward: `run_capped.m:44,50`.

Before the first worker heartbeat, bootstrap and pool startup run at `run_costate_library.m:504–506`. They are governed only by the launcher’s startup policy, not the unit hang policy.

**The strongest bound supportable from this bundle is therefore:**

> Nominally up to 4500 seconds of separately capped execution per full certification, plus uncapped stages and orchestration overhead; no proven finite worst-case inter-progress bound.

## 4.3 `cancel(fut)` is not evidence that a worker was killed/restarted

- **[GAP]** `run_capped.m:4–6,44–50` — **The implementation requests future cancellation; it does not implement or verify process termination/restart.** There is no worker PID handling, pool replacement, cancellation-completion check, or post-cancellation health check.

  A cancellation request is not a general guarantee that arbitrary native code has stopped. Even where MATLAB interrupts ordinary computations promptly, that does not justify the header’s unconditional “worker is killed and restarted.”

  **Concrete fix:** test cancellation of both interruptible MATLAB work and a representative non-interruptible/native workload on this installation. Define a bounded escalation to verified pool-worker/job replacement when cancellation does not stop execution.

This has a particularly bad interaction with a one-worker pool: a stuck task can monopolize the pool, subsequent tasks can time out waiting for that worker, and the client can keep returning named certification failures and beating between them. **The outer heartbeat watchdog then sees activity while useful computation has ceased.**

- **[GAP]** `certify_root.m:160–163,454–458; run_capped.m:44–48` — **Production can run unfenced or repeatedly fail against an unusable pool without a campaign-level infrastructure verdict.** An empty pool merely emits a warning and executes directly.

  **Concrete fix:** require a healthy nonempty process pool for this campaign. Distinguish timeout/pool failure from a terminal numerical inability to advance. Restart or block infrastructure failures rather than silently converting them into short-column coverage results.

### Is 2700 seconds safe?

**No, not under the current callback placement.**

It might become a reasonable **stage-inactivity** deadline after stage-level progress, uncapped-stage treatment, and real pool-cancellation supervision. It is not presently a justified **whole-certification** inactivity deadline.

---

# 5. Test adequacy

## 5.1 What `test_campaign_processes` actually proves

For the reported successful run, it establishes useful evidence that:

- the real launcher started three MATLAB processes that wrote persistent READY markers;
- ordinary shared-queue claiming worked for this schedule;
- four cheap successful units had one logged start each;
- unit 3 was restarted after its owner was explicitly SIGKILLed;
- ordinary exceptions exhausted unit 5’s three-attempt budget;
- the queue reached a sampled terminal state with no observed held lock;
- the killed worker did not write `done`;
- the two survivors wrote `done`.

That is materially better than the pass-2 test evidence.

It does **not** prove exactly-once execution in general. Recovery deliberately provides at-least-once attempts: unit 3 runs twice. The useful desired property is **exclusive active ownership with one accepted committed result**, not universal exactly-once computation.

## 5.2 Specific overclaims and test defects

- **[GAP]** `test_campaign_processes.m:52–55,112–114` — **The `.failed` evidence check is vacuous for the failing unit.** Unit 5 throws before writing `tmpOut`; the assertion explicitly allows no `.failed` files.

  **Concrete fix:** have one failing unit write a partial attempt output and then throw. Require its preserved failure artifact.

- **[GAP]** `test_campaign_processes.m:53–55,88–90,112–113` — **The killed unit is killed before it creates a temporary result.** Consequently, “no stray `.part`” does not test interrupted save cleanup, and it does not prove that publication required ownership.

  **Concrete fix:** introduce barriers after temporary-file creation, during a staged write, and immediately before/after publication.

- **[GAP]** `test_work_queue.m:67–70` — **The authority test bypasses the actual bug.** It explicitly writes `fake.lock.held=false`. A real released claim still has `true`.

  **Concrete fix:** release a claim normally, retain its unchanged copy, then attempt beat/release/publish after a replacement claim is acquired.

- **[GAP]** `test_campaign_processes.m:84–109; test_work_queue.m:128–134` — **No test kills a final attempt.** Unit 3 dies on its initial attempt; unit 5 fails normally; the sequential simulated death is not the final-attempt case claimed in the test header.

  **Concrete fix:** use `maxAtt=1`, claim one unit, SIGKILL its owner, and assert retired/finished without another unit start.

- **[GAP]** `test_campaign_processes.m:123–129` — **A `done` heartbeat is not verified process exit.** The test says survivors “exited cleanly” but checks only their heartbeat records.

  **Concrete fix:** collect exact child exit codes and confirm no client/pool/supervisor processes remain.

- **[CORRECTNESS]** `test_campaign_processes.m:138–140` — **Failure cleanup still matches the job path that lives in the environment, not argv.** The comment says to match a tag prefix, but the implementation does neither that nor exact PID cleanup. Both `pkill -f` expressions ordinarily miss the actual MATLAB command:

  ```text
  -batch run(getenv('CAMPAIGN_JOB'))
  ```

  **Concrete fix:** retain exact launched process identities and terminate/reap those instances and owned descendants. Do not use environment strings as ordinary `pgrep -f` patterns.

- **[ROBUSTNESS]** `test_campaign_processes.m:70,84–90` — **The kill schedule depends on launcher return timing.** By the time all workers are READY and `system` returns, unit 3 can already have advanced far into or completed its 40-second run. The later eight-second pause is not a controlled “mid-unit” barrier.

  **Concrete fix:** make the test unit wait on an explicit release file after announcing the targeted phase.

## 5.3 Smallest additional cases with the highest payoff

These are small, deterministic tests; they do not require another multi-day campaign.

### 1. Released-handle and registry lifetime test

**Two MATLAB processes, one unit, no solver.**

- A acquires then releases a claim, retaining the original struct.
- B acquires the unit and waits.
- A attempts old beat, old release, and old publish.
- B’s lock and metadata must remain intact; A must be refused.
- Separately, A holds a lock, clears the registry function or accesses a path alias, and a third observer verifies that exclusion remains intact or the API fails closed.

This directly catches the most serious current defect and the `fcntl` wrapper hazard.

### 2. Reset at the exact vulnerable boundary

**Two processes, one explicit barrier.**

- Reset pauses after finding/acquiring the candidate unit.
- A worker attempts to claim.
- Resume reset.
- Assert that either reset completed before the claim or reset refused/waited; it must never erase a live reserved attempt.

No sleeps pretending to be races.

### 3. Commit/death matrix

**One unit, cheap artifact, four kill points.**

- after attempt reservation;
- after temporary artifact creation;
- immediately before final rename;
- immediately after final rename, before release.

Include `maxAtt=1`.

Assert correct retirement versus committed completion, no duplicate commit, preserved attempt evidence, and no permanent held state.

### 4. Real supervisor startup/hang/shutdown cases

Make test deadlines configurable and short.

- no heartbeat and no READY;
- continuously updated heartbeat but no READY;
- valid READY followed by immediate failure;
- READY followed by silence;
- normal rapid completion;
- TERM-resistant worker requiring SIGKILL.

Require actual supervisor-issued termination, verified exit/reaping, and a durable exit reason. Testing a harness-issued SIGKILL is not testing the supervisor.

### 5. Real pool lifecycle test

**Use the actual one-worker pool.**

- cancel a running task;
- verify whether the task/process really stops;
- kill the client during a task;
- verify descendant/job cleanup;
- start a replacement and verify available capacity;
- run a healthy multi-stage synthetic certification whose total exceeds the old 2700-second relationship, scaled down for the test.

Then run a small real rib workload through the **generated job**, under the intended 3–5-client contention. That verifies bootstrap, pool topology, validation, save/publication, and actual progress gaps.

### 6. Two-controller and finalization contract test

**Two entry/controller processes with prebuilt tiny valid inputs.**

- simultaneous creation with equal and unequal manifests;
- simultaneous `.launch=true` against one queue;
- simultaneous finalize requests;
- wrong-column/wrong-identity rib;
- two complete ribs plus one invalid rib, exercising the barrier indexing branch;
- pre-existing recent catalog;
- audit exception;
- caller figure/default preservation.

Assert one accepted campaign definition, enforced total capacity, one finalizer, and an invocation-bound result.

If minimizing work further, tests **1–4** are the smallest mandatory fault suite. Tests **5–6** close risks that the generic toy queue cannot address.

---

# 6. Five changes before the next unattended run

Ranked by risk removed:

1. **Make lock ownership an actual lifecycle-controlled capability.**  
   Replace copied-Boolean authority; use idempotent active handles, a strongly retaining protected registry, normalized/permanent lock identity, and exception-safe cleanup. Hold the unit lock throughout reset. Test stale release/publish and registry reset/alias cases.

2. **Fix the progress/deadline and pool-cancellation contract.**  
   The present 2700-second premise is wrong. Add bounded stage-level progress or a real aggregate deadline, fail closed without a healthy pool, and verify cancellation/escalation and pool-child cleanup.

3. **Serialize campaign definition and finalization; validate adopted artifacts.**  
   Introduce a local controller lock, immutable jobs, frozen/resolved manifest inputs, strict unit-aware rib validation, invalid-artifact quarantine, and a real packaging receipt. Fix the barrier indexing bug immediately.

4. **Turn the launcher into an owned supervisor/controller.**  
   Parent and reap exact children, enforce READY independently of heartbeat, verify termination, enforce total capacity including pool workers, record exits, and retain responsibility through recovery and finalization. Give the live autochain the same single-finalizer and explicit-result discipline.

5. **Run the deterministic fault suite, then a real generated rib job under production contention.**  
   Especially: released-claim reuse, reset contention, final-attempt SIGKILL, post-rename death, supervisor hang-kill, client/pool death, and simultaneous controllers. Keep evidence on failure and verify that cleanup actually removes processes.

---

# VERDICT

**No—not yet for an unattended multi-day campaign with the stated guarantees.**

The rewrite **does** remove the pass-2 architecture’s central ABA takeover problem. I would not ask the author to revert to leases, add a distributed service, or introduce elaborate fencing generations merely because multiple MATLAB processes are involved. **A one-host lifetime-lock model is appropriate.**

What stands between this implementation and approval is concrete:

## Blockers

- Released/stale claim objects still retain publication and metadata authority.
- The persistent registry can lose synchronization with live Java/kernel locks, including through a destructive second-channel close.
- Reset mutates records without retaining exclusion.
- The watchdog can kill a healthy full-stack certification; cancellation/child cleanup is not established.
- Campaign creation, launch capacity, and finalization are not serialized.
- Existing/imported ribs are not validated against the campaign/unit identity.
- The new invalid-artifact barrier path can throw.
- The actual supervisor kill/recovery path and real pooled rib path have not been exercised.
- General campaigns still have no unattended completion/finalization owner.

## Hardening, not independent blockers once the above is resolved

- Intra-column checkpoints: highly valuable, but a declared whole-column loss budget is acceptable.
- Hardcoded bootstrap path on this one host: acceptable if preflighted and version-recorded.
- Converting the chain script to a function: desirable, but process/local-workspace isolation can suffice once global graphics effects are removed.
- Status-probe traffic, richer monitoring, formatting validation, and exclusive temporary creation rather than UUID-only naming.
- Power-failure durability beyond atomic rename, if host-crash recovery is part of the required contract.

**The remaining work is no longer “invent a sound ownership architecture.” It is “implement the chosen architecture faithfully, correct the false timing premise, and test its failure boundaries.”** That is a much smaller distance than after pass 2—but it is still a necessary distance before leaving this chain unattended for days.