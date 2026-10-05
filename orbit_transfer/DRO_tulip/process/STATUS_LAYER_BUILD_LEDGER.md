# SDD ledger — plan: docs/superpowers/plans/2026-10-04-library-optimality-status.md

Spec: docs/superpowers/specs/2026-10-04-library-optimality-status-design.md (approved; deviations 1-6 of the plan accepted by Mike 2026-10-04)
Worktree: /Users/msc/Desktop/optimal_control-wt-status (branch status-layer, from main b58b6a4). Untracked result data symlinked in from the main checkout (202 links).
MERGE_BASE: b58b6a4

## Setup rulings
- Ruling: implementers do NOT commit; the controller commits each task after its review passes, staging only the task's files — Mike's standing rule "subagents never commit, one commit per reviewed unit" overrides the skill's implementer-commits step — cost if wrong: none (commits identical, just made by the host).
- Ruling (revised before Task 1 review): the controller commits each task's files on branch status-layer as a CANDIDATE commit before review (fix rounds add commits), so the standard review-package script works; still no subagent commits — cost if wrong: a few extra commits on the feature branch (squashable at merge).
- Ruling: work on branch status-layer in a worktree, merge to main at the finish (Mike's standing commit+push rule covers the push) — cost if wrong: a merge step Mike would not have needed.
- Ruling: tests that solve run with `matlab -batch` from the worktree's tests folder; pure tests may use the MCP session (tests addpath their own folder first, so worktree code shadows any main-checkout path) — cost if wrong: a test exercising main-checkout code (would show as unchanged behaviour; reviewers check the paths in reports).

## Pre-flight scan
| pair / task | produces -> consumes | found |
|---|---|---|
| T1 -> T2 | classifier reads .stage .conjFound .hypAfterConj .flyKm .ok; T2 creates them | consistent; T2's reason formats ('conjugate test verdict %g (%s)', dense 'dense conjugate scan not clear: ...') match T1's tests and regexes |
| T1 -> T3 | make_alternative uses C.status when .stage exists, else optimality_status | consistent |
| T2 -> T5 | rib reads Ct.status; copts.override.conj forces status 2 | consistent; rib's copts already carries pool/sA/sD/progress |
| T3 -> T6 | fill_holes keeps candidates; uses same_root (T3) | consistent |
| T3 -> T7, T8, T11 | make_alternative row fields, dedup_alternatives(A, sheet) needs .has_solution .entry_index .z8 | consistent; T7 builds a 'view' struct with those fields |
| T4 -> T9 | entry_thrust_program(z8, junctions, rv0, rvf, phys) .defectKm .flyKm .startErr | consistent |
| T5 -> T7 | R.refused with .walked, .sD, .sA | consistent |
| T7 -> T8 | sheet fields status/status_reason/junctions, top-level alternatives/status_key | consistent; merge treats status as a grid field and the two cells as per-entry fields (its generic detection) |
| T7 (internal) | package_phase_catalog asserts ribs share the sheet's problem (assertSameProblem) | test saves rib .mat files with R only -> Ruling below |
| T9 -> T11/T12 | audit_status_layer(catMat, opts) .nBad .altContentKey | consistent |
| T10 | standalone | consistent |
| T11 | consumes T1-T10 | consistent; harvest reads S.cand, R.pts/R.refused, rec(k).others, checkpoint pts/refused |
| each task self-consistency | tests vs code | T2 Step 3 wording "Rename line 1 ... to stay as is" is awkward: meaning = keep the public signature, add the wrapper; T9 test's fake certifier fixed in plan commit b58b6a4 |

- Ruling: Task 7's test saves `problem` (the sheet's S.problem) beside `R` in each synthetic rib file, as build_ribs does, so package_phase_catalog's identity check is satisfied — the plan's test omitted it — cost if wrong: none (test fixture only).

## Progress
Task 1: complete (commits b58b6a4..a210877, review clean; test re-run by controller 19/19)
Task 1: minor (deferred): equivalent mutant (pointwise/gates order) survives; add a legacy 'tfMin witness threw' case to cover the swap
Task 1: minor (deferred): optimality_status reads C.hypAfterConj unguarded when conjFound (use fieldd)
Task 1: minor (deferred): legacy 'dense scan: ... malformed' / 'conj_spectrum threw' fall to tier 1 though gates+H6 had passed (arguably 3)
Task 1: minor (deferred): legacy multiplicity-only / sign-change-only dense reasons -> 2 untested
Task 1: minor (deferred): legacy DIAGNOSTIC override reason -> 3 (unreachable in the library)
- Ruling (Task 2): test_certify_enforcement and test_certify_crossing need a pool opened first in a -batch session (their unchanged no-pool assert); every later dispatch runs such suites as `matlab -batch "addpath(<worktree>/orbit_transfer/costate_common); capped_pool(); cd(...); exit(~test())"` (capped_pool lives in costate_common, not on the startup path) — the plan's Step-8 commands omitted it — cost if wrong: none (harness only).
Task 2: complete (commits a210877..430db72, review clean; controller re-run of test_certify_status 8/8, MATLAB exit 0)
Task 2: minor (deferred): dead line `if C.conjFound, C.conjReason = ''; end` (certify_root.m ~552, plan-mandated)
Task 2: minor (deferred): header should say stage 5 is reached only when the coarse conjugate test passes (FAIL path goes 4->6->7)
Task 2: minor (deferred): certify_core repeats a dead nargin<5 line; short header
Task 2: carry -> Task 5 and Task 11: old certificate structs (pre-Task-2 field set) may be concatenated with new ones; use append_point-style field harmonisation, never plain [a, b]
- Ruling (Task 3): make_alternative applies the floor via optimality_status even to stamped results, then uses the stamped status/reason — the plan's code trusted a stamped status below the floor, contradicting its own test and the spec's floor — cost if wrong: none (the floor is the spec's rule).
Task 3: complete (commits 430db72..8dba468, review clean; controller re-run 9/9 on worktree code)
Task 3: minor (deferred): make_alternative reads C.status_reason unguarded for stamped results (use fieldd)
Task 3: minor (deferred): fill_holes_direct physicsFromCatalog keeps an unused `sh` argument (signature kept for callers/test)
Task 3: minor (deferred): dedup_alternatives trusts entry_index > 0 wherever has_solution (data invariant)
- Ruling (Task 4, pre-review): the test fixture is a PHASE-CATALOG entry (library_70mN_24x24_final, cell sD 0 / sA 0.0754, audited to arrive within 1 km of its keyed phase), not dro_tulip_library's anchor (its z8 belongs to sA 0.07537794 but is keyed 0.0754 — a ~50 km mislabel in that list, logged for later); restore the brief's 1 km / 0.1 m/s arrival tolerance — a 120 km gate would hide real errors — cost if wrong: one more test iteration.
Task 4: note (for later): dro_tulip_library's anchor entry is keyed sA 0.0754 but its costates arrive at 0.07537794 (49.5 km miss at the key) — check other entries of that list before trusting their phases
Task 4: complete (commits 8dba468..b2d7a35, review clean; controller re-run 6/6, arrival 0.000/0.001 km)
Task 4: minor (deferred): thrust-direction oracle sampled at one point of the K=12 flight only (check several, incl. a boundary)
Task 4: minor (deferred): entry_thrust_program has no size asserts on z8/rv0/rvf
Task 4: minor (deferred): no guard if tfMinProp returns non-monotonic/short time stamps (interp1)
Task 4: minor (deferred): test file uses the older % NAME header form
Task 5: review — Important: status-1/3 refusal retention untested; Minor: "status 0" wording (true value -1). Fix round 1 dispatched (resume implementer). FIX_BASE 8d39bbe.
Task 5: minor (deferred): R.pts has no entry for a target reached on a status-2 point (with walkPast) — consumers must not index R.pts by step
Task 5: fix round 1/5 (2 addressed, 0 open — status-1/3 retention tests, -1 wording; commits 8d39bbe..1aab427)
- Ruling (Task 5): test fixture nD = 48 (sD 47/48) instead of the plan's nD = 12 — the 1/12 and 1/24 steps out of the anchor fail the polish (the existing rib test needs 8 bisections for the same step) — cost if wrong: none (fixture only).
Task 5: complete (commits b2d7a35..1aab427, review clean after 1 fix round)
- Ruling (Task 6): keepCandidate skips a result whose z equals one already kept (the plan's push points could keep a winner twice: when solved and when displaced); dedup_alternatives repeats this at packaging — cost if wrong: none.
- Ruling (Task 6 review): clearance-floor and failed-direct-solve candidates are correctly not kept — a root under the Moon clearance fails validate_flight (below the spec floor), a failed direct solve has no root — cost if wrong: a few lost non-flyable candidates.
Task 6: complete (commits 1aab427..ec80bbb, review clean)
Task 6: minor (deferred): resume guard `~isempty(rec)` skips harmonising an empty old-field-set rec (then rec(end+1)=cellRec fails); drop the guard
Task 6: minor (deferred): no behavioural test of the winner exclusion or of resuming an old rec without `others`
- Ruling (Task 7): the packager tolerates a rib refusal without `.walked` ("walked: unrecorded") — the plan's own test refusals lack the field; real ribs always set it — cost if wrong: none.
- Ruling (Task 7): test_package_status strengthened with a primary twin and a repeated refusal so the no-dedup mutant is killed — the plan's data had no duplicates — cost if wrong: none.
- Ruling (Task 7 review): sheet candidates of columns with NO certified winner must be collected (the plan placed the loop after the no-winner skip; spec 5 says the packager collects sheet candidates) — fix round 1 dispatched — cost if wrong: none.
Task 7: minor (deferred): a displaced legacy primary becomes an alternative with inferred=false (no INFERRED grid)
Task 7: review — fix round 1 dispatched (FIX_BASE e1f4483): no-winner columns, usableEntry-failed winner, view->sheetView
Task 7: fix round 1/5 (3 addressed, 0 open — winnerless columns collected, unusable winner routed, sheetView; commits e1f4483..c576c6b)
Task 7: complete (commits ec80bbb..c576c6b, review clean after 1 fix round)
Task 7: minor (deferred): the unusable-winner branch of sheet_to_catalog_file has no test
- Ruling (Task 8 review): a merge where only one catalog carries the status layer is REFUSED by name in both directions (layered base + plain donor already refused by the field check; plain base + layered donor gets an explicit error) — silently dropping a layer is never acceptable — cost if wrong: a caller must backfill before merging.
Task 8: review — fix round 1 dispatched (FIX_BASE 5380622)
Task 8: minor (deferred): displaced primary row built with flyKm=flyVms=0 (base flight diagnostics not stored per entry)
Task 8: minor (deferred): a donor primary that loses (within tolDays, different root) is dropped, not kept as an alternative (spec only mandates displaced base primaries)
Task 8: fix round 1/5 (1 addressed, 0 open — one-sided status layer refused both ways; commits 5380622..06cd32e)
Task 8: complete (commits c576c6b..06cd32e, review clean after 1 fix round)
- Ruling (Task 9 review): alternatives are ALSO flown from their junctions and must clear the floor before re-certification; primaries gain velocity checks (junctionVms 1 m/s, flyVms 0.1 m/s) — spec 7 required both, the plan omitted them — cost if wrong: a few seconds of flight per entry.
- Ruling (Task 9 review): full certify_root per alternative is kept as the "targeted" check (it returns at the first failure), documented — spec 7's per-status recompute is satisfied in effect — cost if wrong: longer audits for status-2/4 alternatives.
- Ruling (Task 9 review): unfenced re-certification only by explicit opts.allowUnfenced (default false), never inferred from an empty pool — cost if wrong: none.
- Ruling (Task 9 review): .tolMove removed — same_root is the library-wide identity rule — cost if wrong: none.
Task 9: review — fix round 1 dispatched (FIX_BASE 1f54cb8): fence, alt junction flight + floor, primary velocity, NaN fail-open, setup try, crash-salvage keys, tolMove removed
Task 9: fix round 1/5 (8 addressed, 0 open; commits 1f54cb8..c0fbdad)
Task 9: complete (commits 06cd32e..c0fbdad, review clean after 1 fix round)
Task 9: minor (deferred): with skipPrimaries false, setup runs before the noLayer assert (a layer-less catalog may report a setup failure first)
Task 10: complete (commits c0fbdad..54d1c13, review clean; controller re-run 4/4)
Task 10: minor (deferred): no test for the empty result, the multi-sheet error, or propellantKg; an alternative identical to the primary would be listed twice (dedup prevents it upstream)
- Ruling (Task 11): the backfill job audits in its own chunked stage (run_recertify.sh N audit; primaries in chunk 1) and a separate verdict stage, not inside assemble — 650+ full re-certifications do not fit one process — cost if wrong: one more launch step.
Task 11: real harvest (no solve, 7 s): 2,944 candidates; needRecert 367 + 4 rib stop points; 0/1,152 primaries unmatched; trial assemble 653 alternatives, content key equal
- Ruling (Task 11 review): a 'cand' re-certification that MOVED enters as a new row (a genuine transfer at that phase); the legacy row keeps its conservative inferred status — cost if wrong: one extra alternatives row per moved candidate.
Task 11: review — fix round 1 dispatched (FIX_BASE 9725d7b): junction match by phase, harvest key binding, verdict coverage; minors nRecertErrors, recert dedup, NaN phases, dedup-order comment
Task 11: fix round 1/5 (7 addressed, 0 open; commits 9725d7b..6988427)
Task 11: complete (commits 54d1c13..6988427, review clean after 1 fix round)
Task 11: minor (deferred): duplicate-copy (item 7), NaN-phase skip and recert-overlap dedup are untested; nRecertOverlap not printed in the verdict; onCell 1e-8 assumes identical phase rounding (safe failure: unmatched -> re-polish)
- Ruling (Task 12): the controller runs the backfill stages itself (launch + monitor, no code changes); a subagent writes the docs afterwards — long unattended runs do not fit an implementer's turn — cost if wrong: none.
Final review (opus, b58b6a4..6988427): Ready with fixes — C1 conjugate label on floor-level zeros (CONFIRMED in conj_resolve: nZero = nZeroSign + nZeroFloor); C2 harmonised legacy records read as stamped, audit fails open on empty status; I1 only verdict-0 rows re-certified (other inferred rows cannot reproduce in the audit); I2 dedup keeps first copy regardless of status; I3 rib with no certified points loses refusals; I4 merge remaps donor alternatives by index not phase; minors M1-M7.
- Ruling (final review): the running re-certification was STOPPED at ~25 min (31 items) — C1 would have written wrong status-2 labels and conjDense lacks nZeroSign to correct them; partial recert files to be moved aside, harvest re-run after the fixes — cost if wrong: 25 min of compute.
- Ruling (final review): status 2 needs a DEFINITE refutation: coarse FAIL (trusted sign change) or dense nZeroSign + nInterior > 0; floor-level zeros and corank 'multiplicity' alone give 3; legacy dense reasons (no sign/floor split) give 3 unless coarse sign changes > 0 — cost if wrong: some true conjugate points labelled 3 (conservative, per the established-only principle).
- Ruling (final review): every INFERRED alternative candidate (not the cell primary's root) is re-certified, deduplicated by (phase, root) first — success criterion 3 requires it — cost if wrong: ~3-5 h more compute.
- Ruling (final review): M3 (status-layer key), M5 (report faster status-4 alternatives), M6 (fly convenience wrapper) DEFERRED; M7 already handled (v2 folder is a symlink into the main checkout) — cost if wrong: follow-up work.
Final fix wave: commit a1e01ac (6988427..a1e01ac); scoped re-review: all 9 addressed, no new Critical/Important.
Final: minor (deferred): a candidate whose re-certification moved/errored keeps its inferred row -> audit BAD -> verdict NOT CLEAN (fails closed)
Final: minor (deferred): coverage failure prints FAILED rather than NOT CLEAN; CLEAN rule does not read nOnlyNew (content key covers it)
Task 12: harvest (fixed code) 2,944 candidates, 657 items; recert 8 chunks 15:40-17:24: 0 threw, 0 moved; status 4: 268, 3: 110, 2: 269, 1: 9, below floor: 1. Assemble launched 17:25.
Task 12: assemble 17:25 — 1,152 entries, alternatives 656 (4: 268, 3: 110, 2: 269, 1: 9), re-certified 653, moved 0, content key equal, vs record 1,152/1,152 agree. Audit launched 8 chunks 17:25 (primaries in chunk 1).
Task 12: audit 17:25-19:10 (8 chunks): primaries 1,152 ok; alternatives 651 ok / 5 BAD. VERDICT: NOT CLEAN (audit nBad ~= 0); everything else clean (content key equal, 1,152/1,152 agree, coverage complete, moved 0, errors 0).
Task 12: the 5 BAD are run-to-run instability at numerical thresholds, not code defects: #419/#491/#544 lift margin 9.1/8.7/5.3 vs gate 10 (stored 3 -> audit 2); #21 stored 2 -> audit 3 (gate not recorded; likely cap under load); #219 dim S = 6 (stored 3) -> audit witness cap timeout (1). Decision handed to Mike: how to label threshold-sensitive entries. Docs + merge on hold until then.
- Mike's decision 2026-10-05: OPTION 1 for threshold-sensitive entries; then re-audit -> CLEAN -> adopt v2; merge; docs; backup.
- Ruling (option 1 design): BAD rows that did not move are relabelled to min(stored, re-audit) with `borderline = true` and a reason naming both runs and the check; the audit treats a borderline status as a LOWER BOUND (OK iff same root and re-certified status >= stored); moved rows are never relabelled. A relabel changes the alternatives key, so the old audit chunks are set aside and the full audit re-runs (~1.7 h) — cost if wrong: a borderline row could be under-labelled (conservative by construction).
Borderline change: commits a1e01ac..0d53a43 (1 fix round: relabel only re-established roots that disagree; statuses validated 1..4); re-review clean.
- CONTROLLER ERROR (my option-1 design): "lower of the two" by NUMERIC code is wrong — 2 (conjugate point found) is a stronger claim than 3 (it needs 3 + a reproducible refutation). Relabel attempt 1 produced 3->2 on three rows. UNDONE via the stage's recovery steps (catalog restored from the pre-relabel copy, audit files moved back; the pre-relabel copy kept in library_70mN_24x48_v2_relabel_attempt1/).
- Ruling (corrected): statuses form a lattice by what was ESTABLISHED: 1 < 3 < {2, 4} (2 and 4 both extend 3, incomparable). Conservative label = MEET of the two runs (meet(2,3)=3, meet(2,4)=3, meet(3,4)=3, meet(x,1)=1, meet(x,x)=x). Audit lower bound: a borderline stored status s is reached iff meet(s, now) == s. Expected relabel: k219 -> 1, k419/k491/k544 stay 3, k21 -> 3 — cost if wrong: none beyond rerun.
Borderline fix round 2 (lattice meet): commits 0d53a43..75d65ac, re-review clean (table verified entry by entry).
Task 12: relabel (meet) k219 3->1, k21 2->3, k419/491/544 stay 3 (borderline); round-2 audit 13:08-14:50: 1,808 ok / 0 bad. VERDICT: CLEAN. Adoption next (Mike's step 2).
