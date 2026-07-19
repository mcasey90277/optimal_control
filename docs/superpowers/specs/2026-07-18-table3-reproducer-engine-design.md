# DESIGN — Table-3 reproducer: an extensible thrust-ladder descent engine

**Goal.** A click-go, unattended reproducer that re-solves our version of Gergaud
Table 3 **from scratch** (not by loading the campaign caches), capturing every
per-rung "trick" as runnable code — built so the *same* engine extends to the
0.2/0.1 N push later. This build reproduces the **five attained rows**
(10/5/2.5/1/0.5 N); it does not attempt to crack the 0.5 N min-time wall or
0.2/0.1 N (explicitly next-session work).

Date: 2026-07-18. Status: DESIGN (design agreed conversationally; recipes
harvested from the actual campaign hand-scripts). Depends on the certified MEE
ladder (`earth_elliptic_to_geo/CAMPAIGN.md`) and its drivers.

---

## 1. The key decision: rungs are data, not scripts

The campaign's magic currently lives in scattered hand-scripts
(`results/task7c_*.m` for 1 N; `run_task9_rung.m` for the deep-rung body) and in
driver defaults. We consolidate it into **one parameterized descent engine +
a recipe registry**, so reproducing today's table and pushing to 0.2/0.1 N are
the same tool with different registry entries.

- **`reproduce_row(T, recipe)`** — solves one rung's exact recipe from scratch,
  into its own namespace, resumably; returns the certified row (or errors loudly).
- **recipe registry** (`table3_recipes.m`) — a struct-per-rung of the *proven*
  knobs for 10/5/2.5/1/0.5 N, plus *seeded* (tunable, not-yet-run) entries for
  0.2/0.1 N.
- **anchor strategy is pluggable** — the single biggest lever for the deep push:
  `'coldB' | 'chain' | 'smallN_first' | 'R0law'` (see §4).
- **watchdog orchestrator** (`reproduce_table3.sh`) — runs each rung in its own
  `matlab -batch` process, auto-relaunching on the uncatchable MEX/SIGBUS
  crashes, resuming from caches.
- **verify + probe harness** — each row asserts against the certified campaign
  numbers; deep-rung attempts log their config + outcome systematically.

## 2. Namespace and "from scratch"

Reproduction re-solves; it must not silently load `MEE_M2_*.mat`. All engine
output goes to a fresh namespace **`results/repro/`** (own tags), so a run
genuinely re-derives each row. The warm-chain is internal to the run:
10 N (cold) → feeds 5 N → 2.5 N → 1 N → 0.5 N, each rung consuming the
previous rung's freshly-converged solution from `results/repro/`. Resumability
is automatic (per-round/per-step/per-rung caches under `results/repro/`); a
killed run re-launches and skips completed work. A `reuseCampaignCache` flag
(default false) can instead point at the existing certified `.mat` for a fast
"present the numbers" mode.

## 3. `reproduce_row` contract

```
row = reproduce_row(T, recipe)   % recipe defaults to table3_recipes(T)
```
- Reads `prev` (the previous rung's converged anchor + fuel) from
  `results/repro/` when `recipe.anchor.strategy` needs a warm chain; errors with
  a clear "run rung <prevT> first" if absent.
- Runs the recipe's three stages (anchor → fuel → optional PSR), each resumable.
- Writes `results/repro/REPRO_row_T<10T>.mat` with the full solution + the row.
- Calls `gergaud_row`/`gergaud_row_str` for the printed row; runs the verify
  assert (§5).
- Returns the row struct; throws on any stage that does not certify (so the
  watchdog sees a nonzero exit and the row is never faked).

## 4. Anchor strategies (the pluggable part)

| strategy | used by | what it does |
|---|---|---|
| `coldB` | 10 N | `run_mintime_mee(T, npr)` — cold Stage-B 3-rev seed + keep-best (multi-basin) |
| `chain` | 5, 2.5 N | `run_mintime_mee(T, npr, warmStartAnchor=prev)` — C-law ΔL rescale from the rung above |
| `smallN_first` | 1 N | anchor at `nprLo` (15) warm from prev, **manual relaxed-stall continuation** (maxIter 75→150, no decadeMin floor, always-advance, per-round save — harvested from `task7c_step1_manual.m`), then `interp_warmstart` up to `nprHi` (25) + one `warmTight` solve (`task7c_step1b_refine.m`) |
| `R0law` | 0.5 N (+ future deep) | no min-time solve: `tfMinAnchor = R0const/T` with `R0const = 223.14` ND; records `anchorSource='R0law'` and the footnote |

## 5. Verification — ONE-SIDED (find-best, 2026-07-18 user decision)

**The reproducer is a keep-best-mass OPTIMIZER, not a bit-exact reproducer.**
Minimum-fuel = maximize final mass, so the goal is to find the highest-mass
solution — at least equalling, and where possible BEATING, the campaign. (This
is not hypothetical: at 10 N the keep-best-mass multi-start finds
18 sw / 7.56 rev / **1378.46 kg**, +1.36 kg over the campaign's 1377.10 and
structurally closer to the paper — the campaign under-optimized. See memory
`tenN-minfuel-razor-basin`.)

`table3_certified.m` holds the campaign numbers as the **FLOOR**. `verify_row(row,
cert, tol)` is **one-sided**: it throws only if `row.m_f_kg < cert.m_f_kg -
tol.m_f_kg` (a regression below the floor). A higher mass always passes and is
flagged `info.improved`. **Structure (switches/revs) is REPORTED, not gated** —
the best min-fuel optimum can have a different bang-bang structure than the
campaign row, so gating on the campaign's exact structure would wrongly reject a
better solution. As the full ladder is re-run, the best-found numbers become the
**updated Table 3** (collected in Task 4/6; the documented Table 3 in README/
CAMPAIGN is updated once the improved ladder is validated).

## 6. The recipe registry (harvested, proven)

`table3_recipes(T)` returns:

| T | anchor | fuel | PSR | verify (m_f/sw/revs) |
|---|---|---|---|---|
| 10  | `coldB`, npr 25 | seedThr 0.4, npr 25, tfMinAnchor from anchor | none | 1377.10 / 19 / 7.33 |
| 5   | `chain` from 10, npr 25, mtMaxIter 300 | warm from 10, npr 25 | none | 1364.54 / 32 / 14.16 |
| 2.5 | `chain` from 5, npr 25, mtMaxIter 300 | warm from 5, npr 25 | none | 1369.79 / 76 / 27.84 |
| 1   | `smallN_first` nprLo 15 → nprHi 25, warm from 2.5 | warm from 2.5 (C-law), npr 25 | rounds=2, nbr 2 | 1371.44 / 171 / 69.15 |
| 0.5 | `R0law` (446.28 ND) | coarse npr 12, warm from 1 (C-law) | maxRounds 5, globalEvery 3, globalFactor 1.3 | 1375.28 / 362 / 138.60 |
| 0.2 | `chain`/`R0law` (seed, TBD) | coarse npr 10–12 | maxRounds ≥5 | — (not run) |
| 0.1 | `chain`/`R0law` (seed, TBD) | coarse npr 8–12 | maxRounds ≥6 | — (not run) |

The 0.2/0.1 rows are present so the engine is complete, but are **seeded, not
executed** in this build (§9).

## 7. Watchdog orchestrator

`reproduce_table3.sh` (modeled on `PSR/psr_batch.sh` / `run_task9_watchdog.sh`):
```
for T in 10 5 2.5 1 0.5; do
  until matlab -batch "reproduce_row($T)"  # own process per rung
  do  classify crash via fresh ~/matlab_crash_dump.*; relaunch (caches resume)
      cap relaunches per rung; log each attempt
  done
done
reproduce_table3_collect   # assemble + print the full table + R0 spread
```
A thin `reproduce_table3.m` convenience wrapper runs the crash-free top rungs
in-process for interactive use; the shell path is the robust one.

## 8. Files

```
earth_elliptic_to_geo/
  reproduce_row.m            NEW  per-rung engine (anchor→fuel→PSR, resumable, verify)
  table3_recipes.m           NEW  recipe registry (proven 10..0.5, seeded 0.2/0.1)
  table3_certified.m         NEW  certified numbers for the verify assert
  anchor_smallN_first.m      NEW  the 1 N manual-continuation + refine strategy
                                  (harvested from results/task7c_step1_manual.m + _step1b_refine.m)
  reproduce_table3_collect.m NEW  assemble + print the reproduced table
  reproduce_table3.sh        NEW  watchdog orchestrator (per-process, relaunch, resume)
  reproduce_table3.m         NEW  thin in-process wrapper (top rungs)
  test_table3_recipes.m      NEW  registry + certified-table + verify_row unit tests (no solve)
  test_reproduce_row_smoke.m NEW  10 N live re-solve smoke (proves the engine end-to-end)
  results/repro/             NEW  fresh-namespace outputs
```
`run_task9_rung.m` is refactored into / subsumed by `reproduce_row` (the deep-rung
body it already implements becomes the `chain`+coarse-fuel+PSR path).

## 9. Validation plan (this build)

- **Unit (no solve):** registry returns the right knobs; `verify_row` passes on
  the certified numbers and fails on a perturbed row; certified table matches
  `CAMPAIGN.md`.
- **Live end-to-end:** re-solve **10 N from scratch** in `results/repro/`
  (~minutes) and assert it reproduces 1377.10/19/7.33. Re-solve **5 N and 2.5 N**
  via the `chain` strategy if the build window allows (~tens of minutes).
- **1 N and 0.5 N:** wired faithfully from the harvested hand-scripts and
  checkpoint-resumable; **not required to run to completion in the build** (hours,
  crash-prone). Their recipe code is unit-covered and the anchor `smallN_first`
  strategy is validated against the cached 15/rev anchor where possible. The full
  multi-hour reproduction is the user's click-go via the watchdog.

## 9b. Keep-best-mass multi-start fuel stage (razor-basin finding, 2026-07-18)

**Finding (memory `tenN-minfuel-razor-basin`):** the fuel basin is razor-sensitive
to `t_f` — at 10 N, the *exact* min-time anchor (22.220578) lands the fuel solve
in a slightly-worse local optimum (24 sw / 8.118 rev / 1377.086 kg) while a value
2.2e-5 larger, or a different *seed* at the exact tf, reaches the better basin
(19 sw / ~7.3 rev / up to 1377.19 kg). A single cold solve at the exact anchor is
not guaranteed to find the best optimum.

**Fix — the fuel stage does a keep-best-mass multi-start (A→B hybrid)**, since
min-fuel means *maximize final mass*:
- **A (seeds at exact tf):** solve the fuel problem at the exact `tf = c_tf·tfMinAnchor`
  over a small candidate set — for a cold-seed rung, vary `seedThr ∈ {0.40,0.45,0.50}`
  × `betaMode ∈ {tangential,transverse}` (the driver's own seed-revs window guard
  skips out-of-window seeds; a hung/failed seed is skipped, never fatal); for a
  warm-started rung, the warm trajectory is the seed (optionally with beta
  variants). Keep the highest-mass **certified** candidate.
- **B (tiny tf-bracket, fallback):** only if A's best certified mass does not reach
  the campaign reference (`>= table3_certified(T).m_f_kg - tol.m_f_kg`), also sweep
  `tfMinAnchor` over a tiny bracket (a few ×1e-5, i.e. within `c_tf=1.5`'s
  numerical precision) and keep the best over A∪B. For 10 N, A suffices; B is the
  safety net for any rung where seeds alone can't escape a worse basin.
- **Verify** against the campaign number with tolerance (structure exact/banded,
  mass within tol): the reproducer reproduces the STRUCTURE and returns mass
  `>=` the campaign's — sometimes marginally better (10 N: 1377.19 vs 1377.10).

Confirmed at 10 N (seedThr 0.45 transverse → 19 sw / 7.324 rev / 1377.19 kg at the
exact anchor). This makes the engine a **find-the-best-min-fuel-optimum** tool.
The recurring intermittent MUMPS-init hang is absorbed by skipping a hung seed and
by the watchdog's per-process relaunch.

## 10. Non-goals (this build)

- Cracking the 0.5 N min-time anchor wall, or attempting 0.2/0.1 N to success —
  next session. The registry seeds them; the engine supports them; we do not run
  them here.
- Changing any certified solver behavior. The engine composes the existing
  drivers (`run_mintime_mee`, `run_transfer_mee`, `psr_mee_refine`,
  `casadi_lt_mee`, `interp_warmstart`) — no solver-core edits.
- A GUI. "Click go" = one shell command (or one `.m`).

## 11. Open decisions for the plan

1. Refactor `run_task9_rung.m` in place vs wrap it — lean: subsume its body into
   `reproduce_row`'s `chain` path, keep the file as a thin deprecated shim.
2. Verify tolerances for the PSR-sensitive switch counts (1 N/0.5 N) — lean:
   exact for 10/5/2.5, ±few for 1/0.5 with a logged note.
3. Whether `reproduce_table3.m` (in-process wrapper) is worth building alongside
   the shell watchdog — lean: yes, thin, for the crash-free top rungs.
