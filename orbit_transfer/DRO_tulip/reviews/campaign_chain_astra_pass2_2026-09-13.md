# Second-pass review

**Verdict: still not fit for an unattended multi-day campaign.** There are meaningful fixes, but the ownership protocol still admits concurrent owners, publication is not fenced, a killed final attempt becomes permanently “dispatchable” without being claimable, and finalization still does not consume the library this driver built.

This is a static review of the supplied commit. I have **not** run MATLAB R2026a or inspected its filesystem syscalls. In particular:

- I cannot certify that MATLAB’s `movefile` implements atomic replacement on this installation merely because APFS supports atomic `rename(2)`.
- `certify_root` is not supplied, so its `wallSec` enforcement cannot be audited.
- The actual packaging section of `build_70mN_library`—between the supplied excerpts—is missing. Its visible input configuration already establishes a handoff mismatch; additional internal behavior cannot be assumed.

The sequential tests establish useful ordinary-case behavior. They do not establish the concurrency guarantees claimed in the adjudication.

## 1. Pass-1 findings, one by one

The identifiers below enumerate the original **CORRECTNESS/GAP** findings in their original order.

### Entry point, staging, and resumability

- **[CORRECTNESS]** `build_ribs.m:84–95; work_queue.m:142,200–220` — **P1-01: PARTIALLY FIXED — existence as completion.** Rib saves now use a temporary pathname, but a stalled, short rib is published as DONE; schema, coverage, identity, and exclusive publication remain unchecked. Validate a terminal result before committing it, and distinguish complete coverage from terminal unresolved work.

- **[CORRECTNESS]** `run_costate_library.m:145,170–181,210` — **P1-02: NOT FIXED — artifact identity.** The new guard checks six fields of an optional sheet identity. It does not identify the rib computation, either phase grid, certification policy, or source version. Use an immutable campaign identity and reject incompatible existing artifacts.

- **[CORRECTNESS]** `work_queue.m:94–115` — **P1-03: PARTIALLY FIXED — destructive re-run.** Ordinary `init` no longer removes claims. However, `open` still rewrites shared configuration, replaces `outFcn`, and merges units without serialization or compatibility checks. Normal resume should validate and read an immutable definition, not mutate it.

- **[CORRECTNESS]** `work_queue.m:94–115,118–130` — **P1-04: FIXED — ordinary resume resetting attempts.** Attempt deletion is now an explicit `reset` operation. Reset’s concurrency safety remains a separate problem below.

- **[CORRECTNESS]** `run_costate_library.m:244–290; work_queue.m:200–221` — **P1-05: PARTIALLY FIXED — launch-to-package barrier.** Launcher errors are checked, and successful launch returns before packaging. But `.run.ribs=false` bypasses the barrier, and status counts an artifact-bearing unit as done even while its claim exists. Establish the completion barrier independently of the execution switch and track outstanding ownership separately.

- **[GAP]** `run_costate_library.m:260–288,312–318` — **P1-06: PARTIALLY FIXED — lifecycle.** Explicit launch/pending states now exist. There is still no unattended finalization, and completed-but-unpackaged and packaging-failed paths can remain `pending`. Return stage outcomes with defined semantics; either provide automatic finalization or explicitly define the campaign as worker-only until a separate finalize command runs.

- **[GAP]** `run_costate_library.m:294–306; build_70mN_library.m:58–60,81–109` — **P1-07: NOT FIXED — downstream input handoff.** The packager still receives only output directory and switches. It selects `arrival_sheet_70mN.mat`, not this driver’s `arrival_sheet_70mN_nA24.mat`, and retains its own 12×12 configuration. `pictures` is now configurable, but the material handoff remains wrong. Pass exact inputs and configuration.

- **[GAP]** `run_costate_library.m:307–318` — **P1-08: PARTIALLY FIXED — returned outcome and blockers.** Some blockers are retrieved, and a nonexistent catalog is no longer returned as present. However, existence of an old catalog still authorizes `packaged`; no production receipt or input identity is checked. Return the finalizer’s actual result, including whether it produced, validated, or merely reused an artifact.

- **[CORRECTNESS]** `run_costate_library.m:182–206,247` — **P1-09: PARTIALLY FIXED — missing arrival seeds.** Empty eligibility is handled and missing seeds are reported. The queue denominator still contains only eligible columns, and `onlyA` can exclude additional requested columns without making coverage explicit. Preserve requested, selected, seeded, completed, and unresolved coverage separately.

- **[CORRECTNESS]** `run_costate_library.m:212–235` — **P1-10: NOT FIXED — unowned calibration.** Calibration still executes before queue ownership and attempt reservation. Persisting its timing does not prevent duplicate computation. Run calibration through the queue, or derive timing from the first ordinary claimed attempt.

### Timing and worker behavior

- **[CORRECTNESS]** `run_costate_library.m:396–401; build_ribs.m:81–83; rib_from_crossing.m:88–90` — **P1-11: PARTIALLY FIXED — dropped heartbeat.** The callback now reaches the walker and runs after each returned solve. A bounded inter-beat interval is not established, and ownership-loss/error information is discarded. Audit the complete solve deadline and make loss of publication authority actionable.

- **[CORRECTNESS]** `run_campaign_workers.sh:55–65` — **P1-12: FIXED — process-lifetime watchdog.** The timeout now measures heartbeat age rather than elapsed process lifetime. The new watchdog’s missing-heartbeat, identity, and lease-ordering defects are separate findings below.

- **[GAP]** `run_costate_library.m:213–239` — **P1-13: PARTIALLY FIXED — calibration fallback.** Timing is persisted, but its stored `nD` is never checked, its file is not committed atomically, and a short/stalled column can supply the measurement. The 3600-second fallback remains. Treat timing as compatible planning data, not a safety bound.

- **[GAP]** `run_costate_library.m:242,396–401; run_campaign_workers.sh:38–39` — **P1-14: PARTIALLY FIXED — budget policy.** Worker lease/retry settings are now forwarded. `wallSec` is documented as per point, but the “every budget is measured” claim remains false; the watchdog does not receive the configured lease threshold. Separate and explicitly configure solve deadline, telemetry interval, hang deadline, and drain policy.

- **[GAP]** `rib_from_crossing.m:70–111; build_ribs.m:81–95` — **P1-15: NOT FIXED — intra-column recovery.** No checkpoint survives until the entire walker returns. Reclamation restarts the column. Either checkpoint accepted grid points or explicitly accept a loss budget of up to a whole column per failed attempt.

- **[GAP]** `campaign_worker.m:88–98; run_costate_library.m:399–401` — **P1-16: PARTIALLY FIXED — idle workers leaving the tail.** Workers now wait, but only for one continuous hour by default. They can all leave while another worker continues a nine-hour column, then be absent when that worker dies. Keep a dispatcher alive until campaign termination, or supervise replacement.

- **[CORRECTNESS]** `campaign_worker.m:108–119; build_ribs.m:75–77` — **P1-17: NOT FIXED — normal return counted as successful work.** `okUnit=true` still means only that the function returned. A missing candidate can produce no artifact and still increment `nDone`. Validate the committed result at the work boundary.

- **[CORRECTNESS]** `campaign_worker.m:67–72,133–135` — **P1-18: FIXED — uncaught lifecycle errors leaving ordinary running telemetry.** MATLAB exceptions from the lifecycle now attempt a failure heartbeat and are rethrown. Clean worker exit is explicitly distinguished from campaign completion.

- **[CORRECTNESS]** `campaign_worker.m:108–127,156–168; safe_report.m:38–41` — **P1-19: PARTIALLY FIXED — reporting corrupting accounting.** Worker logging can no longer double-count an attempt or throw through release. But success still follows function return rather than commit, builder reporting is not uniformly isolated, and `safe_report` still asserts that an output exists without knowing that. Base accounting on a commit result and remove that unsupported assertion.

### Queue atomicity and ownership

- **[GAP]** `work_queue.m:149–152; tests/test_work_queue.m:51–55` — **P1-20: NOT FIXED — evidence for exclusive `mkdir`.** The losing-existing-directory case is tested sequentially; actual simultaneous creation is not. Use a documented exclusive operation and test it with separate processes on the target installation.

- **[CORRECTNESS]** `work_queue.m:145–147,229–236` — **P1-21: NOT FIXED — stale reclamation races.** Renaming a pathname does not condition takeover on the owner originally observed. A delayed reclaimer can rename a replacement owner’s fresh directory. Serialize ownership transitions or stop stealing from potentially live owners.

- **[CORRECTNESS]** `work_queue.m:171–183; build_ribs.m:92–95` — **P1-22: PARTIALLY FIXED — old owners acting on new ownership.** Tokens reject sequential late beat/release calls. They do not protect the check-to-action interval, and publication does not check ownership at all. Protect transitions and commit under the same ownership mechanism.

- **[CORRECTNESS]** `work_queue.m:153–167` — **P1-23: PARTIALLY FIXED — post-acquisition eligibility checks.** The checks were added. They are only useful if the acquisition remains exclusive, which the takeover protocol does not ensure. Repair ownership; retain these rechecks.

- **[CORRECTNESS]** `work_queue.m:245–269,278–283,329–334` — **P1-24: PARTIALLY FIXED — persistence failing open.** Malformed readable attempt files now block, and opening/publication errors receive more checking. Write/close success is not checked, atomic replacement is unverified, and a missing attempt file still means zero. Implement a checked persistence primitive and protect counter updates with actual exclusion.

- **[CORRECTNESS]** `work_queue.m:143,202–221` — **P1-25: REGRESSED — state accounting.** Counts now partition the units, but a stale claim on its exhausted final attempt is permanently classified as open while `claim` always skips it. Also, `done` hides outstanding claims. Correct terminal/exhausted-owner handling and expose ownership independently of completion.

- **[GAP]** `work_queue.m:104,136–137,187–188; campaign_status.m:33–37` — **P1-26: PARTIALLY FIXED — policy disagreement.** Options are forwarded consistently within an invocation. No authoritative policy is persisted: a re-run, old job, and default monitor invocation can still apply different policies to the same queue. Store and validate policy with the campaign.

### Heartbeats and monitoring

- **[CORRECTNESS]** `campaign_heartbeat.m:57–62,79–101` — **P1-27: PARTIALLY FIXED — malformed heartbeat accepted as healthy.** Empty/unknown first tokens become `unknown`. A truncated record containing only `beat` still becomes running; `done` alone becomes terminal. The stat/read race remains. Validate the entire record and its expected instance identity.

- **[GAP]** `campaign_heartbeat.m:74–76,93–98; campaign_worker.m:133–135` — **P1-28: PARTIALLY FIXED — labels exceeding evidence.** Worker `done` is now explained correctly. Missing telemetry is still asserted to mean “never started,” and a recent beat is still presented as running without process/progress qualification. Report “no check-in observed,” recent telemetry, and process state separately.

- **[CORRECTNESS]** `run_campaign_workers.sh:40,47; campaign_heartbeat.m:53,60` — **P1-29: PARTIALLY FIXED — terminal identity.** Tags usually differ across launches, but second-resolution timestamps can collide, and records contain no campaign/queue identity. Use an exclusively created launch ID and structured instance-bound records.

- **[GAP]** `campaign_status.m:36–60` — **P1-30: PARTIALLY FIXED — ownership/worker reconciliation.** Stale ownership is no longer called proven death. Claims and worker records are still counted independently rather than joined. Reconcile expected instances, claims, units, and tokens.

- **[GAP]** `run_campaign_workers.sh:55–92; campaign_worker.m:93–97` — **P1-31: NOT FIXED — unattended supervision.** There is no durable controller that collects exits, replaces lost capacity, resolves exhausted claims, detects lack of campaign progress, or finalizes. Idle polling and detached kill loops are not that controller.

### Generated jobs and launcher

- **[CORRECTNESS]** `run_costate_library.m:243,390–391` — **P1-32: PARTIALLY FIXED — source bootstrap.** Solver/common roots now come from the driver. The hardcoded `/Users/msc/Desktop/proj7/external/pumpkynPie` bootstrap remains, visibly in the generated job. Make that dependency explicit, validated, and version-recorded.

- **[CORRECTNESS]** `run_costate_library.m:112–116,333–336,400` — **P1-33: FIXED — generated relative output paths.** The entry point canonicalizes its output directory before deriving generated artifact paths.

- **[CORRECTNESS]** `run_costate_library.m:258–259,340–349; run_campaign_workers.sh:49` — **P1-34: PARTIALLY FIXED — quoting boundaries.** Driver shell arguments and generated MATLAB literals are quoted. The launcher still interpolates `JOB` directly into `run('$JOB')`; an apostrophe in the path breaks or changes MATLAB code. Pass the path as data through an environment variable or manifest.

- **[CORRECTNESS]** `run_costate_library.m:255–256,402–407` — **P1-35: PARTIALLY FIXED — mutable executable job.** Publication replaces a temporary file rather than truncating the live job directly. The job remains mutable shared configuration, and all controllers use the same `.part` name. Use an immutable job/configuration per campaign or launch.

- **[GAP]** `campaign_worker.m:81–87; run_campaign_workers.sh:76–84` — **P1-36: PARTIALLY FIXED — readiness verification.** Exit before the first heartbeat is detected. Any heartbeat file still passes, including `fail`, a malformed file, or an initial beat immediately followed by queue-load failure. Require an instance- and queue-bound READY acknowledgment after preflight.

- **[GAP]** `run_campaign_workers.sh:31–39,46,73,92` — **P1-37: PARTIALLY FIXED — validation and launching zero workers.** Positive counts, positive timing, readable job, writable output directory, and the MATLAB executable are checked. Failed `seq` commands can still leave both loops empty and print success; counts are unbounded and spawned cardinality is never checked. Use a bounded native shell loop and verify child count.

- **[CORRECTNESS]** `run_campaign_workers.sh:40,47–49; run_costate_library.m:393–394` — **P1-38: PARTIALLY FIXED — reused identities and logs.** The fallback worker identity is removed and ordinary relaunches get different names. Concurrent launches within one second still collide and truncate the same stdout logs. Exclusively allocate a launch directory/nonce.

- **[CORRECTNESS]** `run_campaign_workers.sh:55–66` — **P1-39: PARTIALLY FIXED — stale watchdog killing a reused PID.** Polling reduces the old long blind interval. `kill -0` still proves only PID existence, and the unconditional delayed `kill -9` is especially unsafe after the original process exits. Supervise exact children and reap them; do not control an instance solely by a recycled number.

### Tests

- **[GAP]** `tests/test_work_queue.m:51–79,142–148` — **P1-40: NOT FIXED — race tests that never race.** Every operation remains sequential. The purported post-mkdir recheck test creates the output before `claim`, so the *first* existence check passes the test without exercising the new recheck. Add process barriers at the actual vulnerable boundaries.

- **[GAP]** `tests/test_work_queue.m:81–165,180–181` — **P1-41: PARTIALLY FIXED — missing assembled-chain tests.** Persistence, retirement, malformed telemetry, and failure exit status gained useful tests. There are still no worker-process, launcher, generated-job, identity, publication-crash, or finalization tests. Generic one-byte artifacts are fine for a generic queue test; they do not substitute for rib validation tests.

### Grouped pass-1 robustness/style items

- **[CORRECTNESS]** `build_70mN_library.m:48–51,269–275; run_costate_library.m:302–306` — **Workspace isolation: NOT FIXED; severity upgraded.** Restoring `chainOverrides` does not restore cleared variables, `pwd`, figure defaults, or closed figures. This is actual data loss in the stated shared session. Use a function or separate MATLAB process.

- **[ROBUSTNESS]** `campaign_worker.m:108–125` — **Shared-failure retry policy: NOT FIXED.** Deterministic infrastructure errors can exhaust every unit immediately. Preflight dependencies and storage; distinguish campaign-wide failures from unit-local retryable failures.

- **[ROBUSTNESS]** `work_queue.m:159–165,322–325` — **Claim initialization recovery: PARTIALLY FIXED.** Attempts precede work, but missing-beat age still comes from a directory listing and the minimum of its entries, not a defined reservation timestamp. Use explicit reservation metadata and a defined initialization state.

- **[ROBUSTNESS]** `work_queue.m:325; run_campaign_workers.sh:59` — **Clock/storage assumptions: NOT FIXED.** NFS support is irrelevant to this deployment, but local wall-clock changes and sleep still affect age-based decisions. Specify local-filesystem support and avoid age alone authorizing ownership transfer.

- **[ROBUSTNESS]** `run_campaign_workers.sh:28,37; run_costate_library.m:400` — **Heartbeat-directory agreement: NOT FIXED.** Requiring argument four removes the conflicting default, but the launcher and embedded job can still name different directories. Both should consume one authoritative configuration.

- **[ROBUSTNESS]** `run_campaign_workers.sh:73–90` — **Startup deadlines: PARTIALLY FIXED.** Early exit is recognized; sequential per-worker waits and abandoned late starters remain. Supervise concurrent per-child deadlines and retain responsibility for unready children.

- **[ROBUSTNESS]** `run_campaign_workers.sh:49,55–66` — **Detachment: PARTIALLY FIXED.** Worker output and watchdog output no longer retain the launcher’s output pipe. There is still no explicit watchdog lifetime, signal, stdin, or child-reaping contract. Use one owned supervisor/service.

- **[ROBUSTNESS]** `run_costate_library.m:109–140; run_campaign_workers.sh:31–39` — **Configuration/resources: PARTIALLY FIXED.** Executable selection and basic launcher checks improved. Engine/grid/policy validation and aggregate worker/thread/memory limits remain absent. Validate before expensive stages and cap campaign-wide concurrency.

- **[ROBUSTNESS]** `tests/test_work_queue.m:168–181` — **Tests failing CI: FIXED.** Failure now throws. The reporting regression still needs a deliberate exception test in addition to the version-specific conversion case.

- **[ROBUSTNESS]** `fmt_num.m:34–39` — **Fixed-width overflow: NOT FIXED.** Define overflow behavior and validate width/precision.

- **[STYLE]** `work_queue.m:19–42; run_costate_library.m:25–27,74–85` — **Guarantee language is ahead of evidence.** “Exactly one,” “safe repeatedly,” “every budget,” and “never while claimed” are not current contracts. Correct the documentation alongside the implementation.

## 2. Defects in, or exposed by, the fixes

### A. Rename takeover is still an ABA race—even assuming perfect atomic rename

- **[CORRECTNESS]** `work_queue.m:145–167,229–236` — **Two reclaimers can both succeed because the source pathname is reused.** This execution requires neither a random-name collision nor a non-atomic filesystem:

  1. A and C both observe B’s claim as stale.
  2. A renames B’s claim to A’s unique tomb and removes it.
  3. A creates `j.claim`, writes its own token and attempt, and starts work.
  4. C resumes its already-authorized `takeover(j.claim)`.
  5. C renames **A’s fresh claim** to C’s different tomb.
  6. C creates `j.claim` and starts the same unit.

  The post-`mkdir` checks do not help: no result exists yet, and both acquisitions can be within the attempt budget. A single delayed reclaimer can likewise act after the original owner refreshes its beat.

  **Concrete fix:** for this deployment, hold a kernel-released per-unit lock throughout computation and publication. Never steal it based on heartbeat age. A supervisor may terminate a hung holder and wait for its confirmed exit; only then can another worker acquire the lock. Alternatively, serialize all lease transitions and compare the exact observed owner under that serialization, with publication fencing.

### B. Token checks do not make beat/release conditional operations

- **[CORRECTNESS]** `work_queue.m:174–182` — **`ownedBy()` followed by a pathname operation has the same stale-observation problem.**

  - Old owner A passes `ownedBy`.
  - A pauses.
  - C takes over; a replacement claim is created.
  - A resumes and overwrites the replacement’s `beat`, or recursively removes the replacement claim.

  One host does not eliminate this interleaving. Process scheduling, memory pressure, filesystem calls, and intentional recovery are sufficient. These windows are short in healthy operation, but the system’s recovery path deliberately brings delayed owners and reclaimers together.

  **Concrete fix:** ownership checking and mutation must share a lock/transaction. Do not delete or replace the pathname used for a kernel lock; doing so creates a new lock inode and defeats the lock.

- **[CORRECTNESS]** `campaign_worker.m:140–152; rib_from_crossing.m:90; build_ribs.m:92–95` — **Loss of ownership is explicitly ignored, and old owners retain publication authority.** `work_queue('beat')` can return false; `beatBoth` discards the result and still updates the worker heartbeat. Exceptions are also swallowed twice. The superseded worker can remain apparently healthy while continuing to overwrite the final rib.

  **Concrete fix:** distinguish telemetry failure from loss of ownership. On lost/uncertain publication authority, preserve work to an attempt-specific checkpoint, but do not publish the final artifact. Prefer a lifetime lock that remains the publication authority even if telemetry fails.

### C. What `movefile` does—and what has not been established

- **[GAP]** `work_queue.m:235,268,283,334; campaign_heartbeat.m:62; build_ribs.m:94; run_costate_library.m:406` — **APFS atomic rename is not proof that these MATLAB operations are atomic publication.**

  The distinctions matter:

  | Operation | What can be concluded |
  |---|---|
  | Native same-filesystem `rename(2)` to an absent pathname | Atomic namespace change. |
  | Native `rename(2)` replacing an existing regular file | Atomic replacement: a pathname lookup sees the old or new file, not an intentional absent interval. |
  | MATLAB `movefile` to an existing directory | It can move the source **inside** that directory; this is not rename-to-an-exact-nonexistent-name semantics. |
  | MATLAB `movefile` in this R2026a installation | The supplied code and sequential tests do not reveal whether each path uses rename, pre-removal, or copy/delete. |

  Same-directory temporary files make a same-filesystem rename possible. They do not force the wrapper to use it. Nor does successful ordinary overwriting prove crash atomicity. Adding `'f'` would not establish that property.

  **What does a losing takeover contender see?** If the source remains absent, failure; if a replacement claim has appeared, it can successfully move that replacement. If its tomb already exists as a directory, directory nesting semantics apply. The actual ABA defect above survives the best possible `movefile` implementation.

  **Can readers see no file where an old file existed?** Not from a correctly executed atomic regular-file replacement. They can if the wrapper removes the old destination before replacement or implements another multi-step overwrite. An interrupted such replacement is particularly serious for `.att`: missing becomes zero. Claim takeover itself also deliberately creates a missing-claim/beat interval.

  **Concrete fix:** provide one tested publication primitive with explicit same-filesystem atomic-replacement semantics and checked errors. A small native wrapper is reasonable; an atomic-move API must fail rather than silently downgrade. Inspect/trace the actual R2026a behavior on APFS, including replacement of an existing file and injected interruption. Do not base the ownership protocol on an undocumented wrapper detail.

### D. The temporary names are not exclusive

- **[CORRECTNESS]** `build_ribs.m:92–94; run_costate_library.m:402–406` — **Concurrent writers share the same `.part` pathname.** The queue currently permits duplicates, calibration bypasses it, old workers do not participate, and two driver calls can generate jobs concurrently.

  One process can move a temporary file while another is still writing that pathname. Depending on how the save opens the file, the first process can publish the other process’s incomplete bytes, or the second process can continue writing through a handle to the now-published file. At minimum, one process can publish the wrong producer’s temporary result.

  **Concrete fix:** create attempt-specific temporary files exclusively. Publish only while holding current ownership. Atomic final rename does not compensate for sharing the temporary file.

- **[ROBUSTNESS]** `work_queue.m:234,264,279,332,338–343` — **`randomHex` has substantially less uniqueness than its use suggests.**

  - Eight hex characters provide at most **32 bits** of suffix space, not 64.
  - Tomb suffixes have 24 bits; metadata temporary suffixes have only 16.
  - The seed space is only 31 bits.
  - The time component wraps modulo \(2^{31}\) microseconds—about **35.8 minutes**.
  - Adding a PID to time is not a uniqueness construction; processes can generate the same seed, and clock resolution can correlate calls.

  Collision relevance is scoped to the same base pathname and overlapping/orphaned names—not every random value generated in the campaign. Different worker tags also distinguish claim tokens despite equal random suffixes. Thus random collision is **not** the primary ownership failure here.

  It is nevertheless unnecessary risk, particularly because an existing tomb changes `movefile` semantics and an existing temporary file is opened with truncation.

  **Concrete fix:** use an OS-generated UUID and exclusive temporary creation. Treat existing names as collisions to retry, never as destinations to overwrite.

### E. Final-attempt crash is now a nonterminal limbo

- **[CORRECTNESS]** `work_queue.m:143,202–221; campaign_worker.m:92–97` — **A killed third attempt can never retire and can never be reclaimed.**

  Reproduction:

  1. Claim a unit on attempt three.
  2. Kill its worker without releasing.
  3. Let the claim become stale.
  4. `status`: `stale`, `nOpen=1`, `finished=false`.
  5. `claim`: attempt count is three, so it skips the unit before considering takeover.

  Every re-run can announce open work and launch workers that perform none. Idle heartbeats make those workers look active until their idle timeout.

  **Concrete fix:** represent exhausted-held ownership explicitly. Once the owner is confirmed dead and its lock is free, retire the unit without granting another attempt. `nOpen` must describe genuinely claimable work; terminal status must not depend on a release from a dead process.

### F. Completion hides ownership and numerical incompleteness

- **[CORRECTNESS]** `work_queue.m:200–221; run_costate_library.m:279–285` — **“No running units” does not mean “no claims.”** The `done` branch precedes claim inspection. Immediately after saving, but before release, a unit vanishes from running ownership. An old superseded worker can also remain active after another worker publishes.

  **Concrete fix:** maintain separate facts for committed result validity and outstanding ownership. Under immutable, fenced publication, finalization can safely consume committed artifacts even if cleanup metadata remains—but then document that contract rather than claiming no owner exists. Until publication is fenced, require a real quiescence barrier.

- **[CORRECTNESS]** `rib_from_crossing.m:98–101; build_ribs.m:84–95` — **An ordinarily returned stalled rib satisfies the queue’s complete-library predicate.** `R.stop` records failure to advance, but nothing consuming the artifact inspects it or verifies the requested `nD-1` targets.

  A sparse terminal numerical result can be a valid campaign artifact. It is not automatically complete grid coverage.

  **Concrete fix:** validate the exact requested coordinates and record per-target completed/unresolved status. A terminal short rib should produce explicit coverage blockers, not silently count as a fully completed column. For the generic multi-column `build_ribs` interface, saving after each rib also means file existence cannot imply the entire requested list finished.

### G. Idle polling does not supervise the tail

- **[GAP]** `campaign_worker.m:93–97; run_costate_library.m:399–401` — **The one-hour idle timeout starts too early for the claimed recovery guarantee.** It is not reset by another worker’s useful progress. With only one long column left, all peers leave after an hour even while that column is productive. A subsequent failure has no reclaimer.

  **Concrete fix:** keep at least one dispatcher alive until the campaign is complete or explicitly blocked, or have a supervisor replace dispatchers. An inactivity watchdog must distinguish useful campaign progress from a worker repeatedly saying “idle.”

### H. The new timing policy does not establish a safe silence bound

- **[GAP]** `rib_from_crossing.m:71–90; build_ribs.m:39–56; run_costate_library.m:397` — **“After every solve” is a cadence location, not a time bound.** If `certify_root` enforces a hard 900-second deadline around its entire call, then the nominal largest solve-related silence is approximately **900 seconds plus surrounding work**, not two or three minutes.

  The supplied code does not establish that hard bound. A cooperative deadline can be exceeded by a native call, pool wait, or an expensive gate between checks. Initial setup, endpoint construction, loading, logging, and final save are also outside the callback.

  Consequently, the provable worst silence from this bundle is **unbounded**. This does not prove ordinary solves exceed 1800 seconds; it means the required premise has not been verified.

  **Concrete fix:** inspect `certify_root` and all deadline-bearing stages, measure maximum observed callback gaps under production contention, and put progress/deadline checks around substantial non-solve stages. Define what happens when a native call cannot be interrupted safely.

- **[CORRECTNESS]** `run_campaign_workers.sh:38–39,59–61; work_queue.m:146–147` — **Lease expiry precedes watchdog termination by hours.** With fallback `unitSec=3600`, takeover begins after 30 minutes of silence while termination waits four hours. A measured nine-hour column produces a 36-hour watchdog. Thus a hung or merely delayed worker can remain alive long after replacements start.

  The comment “never shorter than the stale lease” is also only true for the default: the launcher never receives a customized `staleSec`.

  **Concrete fix:** do not transfer ownership to another process while the old process can still publish. For one host, use lock lifetime plus supervised kill-and-confirm-exit. Size hang detection from bounded progress expectations; use column duration only for planning.

### I. Resume is still configuration mutation

- **[CORRECTNESS]** `work_queue.m:98–115; run_costate_library.m:247,255–256` — **Concurrent `create/open` calls can lose changes or change live workers’ interpretation.** Two opens can load the same unit set, add different units, and last-writer-win away one addition. `outFcn` changes without checking compatibility. Already running generated closures may continue writing one output definition while future claims consult another.

  The union behavior also means `.onlyA` does not restrict an existing queue: old units remain. A later run at lower resolution can retain out-of-range old unit IDs. `reset` does not remove those units, despite the driver suggesting units can be “dropped.”

  **Concrete fix:** normal open must be read-only and compatibility-checking. Put intentional expansion/migration behind a controller lock and an explicit operation; reject contraction or identity change. Prefer explicit artifact-path data over a mutable saved closure.

- **[CORRECTNESS]** `work_queue.m:121–130` — **Reset’s freshness check is not exclusion.** An idle worker can claim after the check; reset then deletes its attempt record, and possibly its claim. A stale owner may also still be alive. `.force` explicitly permits destruction but does not prevent that process subsequently publishing.

  **Concrete fix:** require a quiesced campaign under controller exclusion, with all worker instances stopped or blocked from new claims. Reset individual retired units explicitly rather than globally clearing unrelated retry history.

### J. Section 0 does not actually control the declared phase grid

- **[CORRECTNESS]** `run_costate_library.m:135–141,161–163,229–230,353–361; rib_from_crossing.m:61,75` — **Accepted explicit phase vectors are discarded and replaced by different grids.**

  Examples:

  - `.sD = 0.1 + (0:23)/24` passes and is printed, but only `nD` reaches the builder. Departure still starts at the sheet’s anchor.
  - A four-element `.sA = [0.1 0.2 0.3 0.4]` is uniform and passes. The sheet builder instead uses `0.1 + (0:3)/4`.
  - Repeated phases pass because their differences are uniformly zero.
  - Descending uniform grids pass, but the builders use their own ordering.

  **Concrete fix:** either forward exact target coordinates through the builders, or reject unsupported vectors. The smallest safe restriction is to accept only the full supported \(1/n\) periodic lattice, check its origin against the actual sheet anchor, and compare the resulting sheet phases against the request.

- **[CORRECTNESS]** `run_costate_library.m:170–181; build_ribs.m:48–68,88` — **The identity guard is optional and incomplete.** A sheet with no `problem` bypasses it. A sheet with the right six scalars but the wrong arrival origin or departure anchor passes. Rib setup compares seven fields, not the entire stamped identity; it then writes the sheet’s identity into the result.

  `calibration.mat`, existing ribs, the queue definition, the generated job, and downstream artifacts are not checked against this guard.

  **Concrete fix:** require identity for accepted artifacts, with an explicit legacy-import validation path. Include exact grids, materialized operating-point identity, schema/certification policy, and relevant code/input provenance. Validate legacy outputs once before allowing their existence to imply completion.

  The source-arc path is also fixed to `here/results` in `build_arrival_sheet.m:49`. Reusing foreign-problem roots as **recertified guesses** can be legitimate, but the supplied `crossings_from_arc`/`sheet_from_arcs` implementations are absent, so this bundle does not establish that every such reuse is recertified. Do not confuse a newly stamped sheet identity with proof of source compatibility.

### K. Finalization still packages a different input contract

- **[CORRECTNESS]** `run_costate_library.m:145,294–306; build_70mN_library.m:52,58–60,100–109` — **The finalizer is not connected to the newly completed grid.** Its selected sheet basename differs from the driver’s; its configuration remains default 12×12; its documented input rib location is `results/`, not the campaign’s per-column files.

  If those old/default inputs exist, this is potentially a wrong-library operation rather than a clean missing-file failure.

  **Concrete fix:** make finalization accept `sheetMat`, exact validated rib paths, grids, identity, and output path. Do not repair this with copies or symlinks to default basenames: that preserves ambiguous provenance and shared mutable aliases.

- **[CORRECTNESS]** `run_costate_library.m:312–318` — **`packaged` means “a suitably named file exists after the chain returned,” not “this call produced this campaign’s catalog.”** A pre-existing catalog passes the predicate whenever the invoked chain returns normally, regardless of whether it replaced or merely retained that file. The supplied packaging excerpt is insufficient to enumerate its internal skip paths.

  Conversely, a normal return without the catalog adds a blocker but leaves state `pending`. Audit-only/sweep-only runs also end with no catalog path and no meaningful completed-stage state.

  **Concrete fix:** return a structured finalizer result containing input identity, output identity, produced/reused status, validation result, and blockers. Use that result—not a final `isfile`—to select state.

- **[CORRECTNESS]** `build_70mN_library.m:48–51,121–124,219–224; run_costate_library.m:306` — **Finalization can erase the user’s workspace and then throw without returning blockers.** Missing prerequisites or an audit exception escape before `chainBlockers` is assigned. This is worse than merely failing to restore a variable: the function can destroy interactive state and provide no `out`.

  **Concrete fix:** execute a function or isolated child process. Catch stage failures into structured outcomes while preserving a hard failure for invalid inputs; do not use the base workspace as an RPC mechanism.

### L. The launcher can still misverify, leak, or launch no useful work

- **[CORRECTNESS]** `run_campaign_workers.sh:58–65,76–90` — **A process that never writes its first heartbeat is never killed by the watchdog.** The watchdog only computes age when the file exists. Verification eventually exits failure, but a license wait, startup hang, or wrong-heartbeat-directory process can survive indefinitely.

  **Concrete fix:** track a startup deadline from spawn independently of heartbeat existence. On timeout, deliberately terminate/reap or keep supervising the unready instance.

- **[GAP]** `run_campaign_workers.sh:77–81; campaign_worker.m:81–87` — **A failed worker can still produce `LAUNCH OK`.** Initial beat happens before queue load; a subsequent fatal failure writes `fail`. The shell checks only file existence, before checking PID state.

  **Concrete fix:** require READY after preflight and inspect both structured readiness and child exit. A worker that already completed an empty queue is a different startup outcome, not proof that the intended campaign is being served.

- **[CORRECTNESS]** `run_campaign_workers.sh:40,47–49,61` — **Instance uniqueness and process identity remain insufficient.** Same-second launches share tags/logs; the five-second escalation can target a reused PID; the launcher does not collect final exit statuses. Repeated `.launch=true` calls can also add batches without any campaign-wide concurrency ceiling.

  **Concrete fix:** allocate instance directories exclusively, keep an expected-worker registry, supervise owned children/process groups, and enforce the desired total worker count rather than adding four on every request.

- **[GAP]** `run_campaign_workers.sh:46,73,88–92` — **Positive count validation does not prove children were created.** If `seq` is missing or fails, both loops can execute zero times and the empty failure array produces success.

  **Concrete fix:** use zsh arithmetic loops with validated bounded integers and assert the exact number of created children. Verify the job is a regular readable file and preflight the actual heartbeat/log directories, not just their parent.

### M. Interrupted saves still poison resume outside the rib publisher

- **[ROBUSTNESS]** `build_arrival_sheet.m:104; run_costate_library.m:215–218,235` — **Sheet and calibration are still saved directly to their final names.** An interrupted write leaves a file that suppresses reconstruction or repeatedly fails loading. A malformed calibration can prevent an otherwise completed campaign from reaching finalization.

  **Concrete fix:** use checked atomic publication for these artifacts too, validate on load, and quarantine/rebuild invalid artifacts under controller exclusion. Check write/close results for text metadata before moving it into place.

### N. The live old launcher is outside every new guarantee

- **[CORRECTNESS]** `work_queue.m:142–149; run_costate_library.m:247` — **Starting this queue before old explicit-range workers finish can reproduce today’s duplicate-column incident immediately.** Those workers hold no queue claims. An absent final file makes their in-progress column claimable to the new queue. New tokens and heartbeats cannot coordinate with a nonparticipant.

  **Concrete fix:** perform a one-time cutover: drain or stop all old workers, verify their exits, validate imported completed artifacts, then activate the queue. Do not infer old-worker quiescence from file existence.

## 3. Rulings on the pushbacks

- **[GAP]** `work_queue.m:104–115; run_costate_library.m:170–181,210` — **Declining a heavyweight manifest system: acceptable. Declining immutable campaign identity: not acceptable.** Five local workers do not require a distributed metadata service. They do require agreement about what unit 14 means and which outputs count.

  **Smallest closure:** one immutable configuration file in a campaign-specific directory, containing exact grids, operating-point identity, input hashes/provenance, artifact schema, policy, code revision, and explicit paths. Resume validates it. No mutable generational database is necessary for this deployment.

- **[CORRECTNESS]** `work_queue.m:145–182; build_ribs.m:92–95` — **Declining fencing generations: acceptable only if live-owner takeover is eliminated.** With expiring leases and concurrent old publishers, rejecting fencing is not justified by one-host scope.

  **Smallest closure:** a per-unit kernel-held lock retained through final publication, with an exact-child supervisor that kills and confirms termination before recovery. All writers—including calibration—must participate. This avoids needing fencing generations because an old live process cannot coexist with a new lock owner.

- **[GAP]** `work_queue.m:316–325; run_campaign_workers.sh:55–66` — **Declining a lease-clock coordinator: acceptable with a non-lease local ownership mechanism.** A coordinator-owned clock is not intrinsically required here. The current alternative—wall-clock age authorizing takeover of a potentially live publisher—is not adequate.

  **Smallest closure:** use heartbeat age for alarms/supervised cancellation, not ownership transfer. Let process-bound locks decide whether another worker may start.

- **[GAP]** `tests/test_work_queue.m:51–79,142–148` — **Deferring separate-process tests: not acceptable before unattended deployment.** Concurrency is the core feature being introduced, and the current implementation contains races the sequential tests cannot reach.

  **Smallest closure:** a miniature campaign using real MATLAB processes, controlled barriers, cheap unit functions, and short test deadlines. Required cases: simultaneous fresh claimers; simultaneous reclaimers; old owner delayed across beat/release/commit; final-attempt kill; live resume; interrupted replacement; malformed READY; missing first heartbeat; and end-to-end finalization of the exact intended inputs.

- **[GAP]** `rib_from_crossing.m:70–111` — **Deferring intra-column checkpoint/resume: conditionally acceptable.** I would not make exact bisection-state recovery a correctness prerequisite once ownership, publication, and supervision are safe. The operator must explicitly accept losing up to nine hours per failed attempt, with potentially several such losses before retirement.

  **Smallest improvement:** checkpoint after each completed grid target, including the accepted root and remaining target list. Resume from the last accepted target; preserving every intermediate bisection is not required initially. The documentation must stop implying reclamation preserves already-computed work.

- **[CORRECTNESS]** `run_costate_library.m:226–235` — **Deferring calibration ownership: not acceptable while the driver can run it on resume.** The easiest safe deployment option is to disable driver-side calibration and persist timings from ordinary queue attempts. There is no need to solve an unowned column to size an inactivity watchdog.

## 4. Five changes before the next unattended run

Ranked by the immediate deployment risk removed, with cheap deterministic protections ahead of larger changes of comparable risk:

1. **Remove the unsafe driver/finalizer shortcuts.** Disable unowned calibration; enforce completion independently of `.run.ribs`; isolate finalization from the base workspace; pass the exact sheet and rib inputs and return an actual production result.

2. **Freeze and validate campaign identity.** Use a new campaign directory or validated immutable configuration; reject unsupported explicit grids; validate legacy artifacts before adoption. Drain the old explicit-range workers before enabling queue claims.

3. **Replace stale-path ownership with a real local exclusion/publication protocol.** Prefer lifetime kernel locks through commit, checked atomic publication, exclusive temporary files, and correct exhausted-owner retirement. Do not try to repair this by adding another token recheck.

4. **Install one actual supervisor.** Bound startup, require structured READY, retain exact-child ownership, enforce total worker count, keep recovery capacity through the tail, collect exits, and schedule finalization. Base hang detection on progress bounds rather than measured column duration.

5. **Run the assembled miniature campaign and multiprocess fault tests on this exact MATLAB/APFS host.** Test overwrite/crash behavior of the chosen publication primitive explicitly. Then run a small real workload under the intended 4–5-process resource configuration before committing to the multi-day campaign.

## VERDICT

**No—not yet.**

The biggest pass-1 failures were addressed in recognizable ways: ordinary resume preserves claims and attempts, the callback reaches the solver loop, the watchdog no longer measures lifetime, and launch normally returns before packaging.

But three present defects alone are disqualifying:

- **A delayed reclaimer or old owner can still destroy new ownership, and superseded workers can still publish.**
- **A killed final attempt can remain permanently open but unclaimable.**
- **The finalizer does not receive the campaign’s actual input contract and can erase the interactive workspace.**

Between this code and “yes” stand: **exclusive ownership through commit, validated immutable computation identity and completion, correct terminal/recovery behavior, isolated exact-input finalization, and real-process verification of the assembled chain.**

A distributed coordinator and elaborate fencing manifest are not mandatory on one host. A concurrency protocol that is actually exclusive—and tested as such—is.