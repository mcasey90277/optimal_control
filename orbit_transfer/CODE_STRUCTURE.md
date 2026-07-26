# Code structure — survey synthesis and extraction plan

Consolidates the four campaign surveys of 2026-07-26 into one ranked, gated
plan. Per-campaign detail lives in each campaign's `doc/EXECUTION_PATHS.md`;
this file owns the cross-campaign picture and the decisions.

Survey was **read-only** — no code changed. Tool: `docs/tools/callgraph.py`.

| campaign | `.m` | sharing strategy | measured outcome |
|---|---|---|---|
| `earth_elliptic_to_geo` | 136 | is the source | already well-layered; the reference pattern |
| `earth_elliptic_to_geo_CR3BP` | 18 | **extend** — path delegation + opt-in `par.pert` branch | **zero duplication**; the `opti.dual` fix arrived for free |
| `GTO_tulip` | 266 | **vendor** — `PSR/lib` frozen snapshot | 21 duplicate basenames, 14 identical; a live 4.7 kB drift that silently broke a feature |
| `GTO_ELFO` | 24 | **fork** — own `*_freetf` solvers | 76% identical to the parent; 34-file path dependency for one function |

---

## 1. The finding that should drive the library work

**The repo already ran the experiment.** Three campaigns each faced the same
problem — *same solver, different transfer* — and each answered it differently.
The outcomes are not close:

- **Extend won.** CR3BP added lunar gravity as an opt-in branch inside the
  *shared* `lt_mee_rhs`, and delegates its path to the 2-body campaign. It
  duplicates nothing, and when the `opti.dual` bug was fixed upstream the fix
  reached it automatically.
- **Vendoring silently failed.** `PSR/lib` froze a snapshot for
  self-containment. The isolation works exactly as designed — and that is the
  problem: upstream fixes never arrive. `casadi_minfuel_sundman` drifted 4.7 kB
  apart, and `run_psr` stage 5c called a function that could not resolve,
  degrading to a caught warning on every run, undetected for a day.
- **Forking left 76% duplication.** ELFO's `casadi_energy_freetf` shares 122
  code lines with tulip's solver, including a 41-line identical IPOPT options
  block — while the *substantive* 24% (the `cScale` state, free `t_f`, the
  two-primary clock) genuinely justifies a separate file.

**Standing rules derived from this evidence:**

1. **Prefer an opt-in branch in the shared file** over vendoring or forking.
   The winning move was not "extract a common core" — it was "add a switch to
   the existing one."
2. **Share problem *definitions* freely; share *algorithms* only behind a
   gate.** `cr3bp_common` (6 functions, up to 18 callers) is frictionless
   because it holds parameters and endpoints. Every painful case in this repo
   is a shared *algorithm*.
3. **If you must vendor, add a drift test.** A frozen copy with no
   byte-identity check is a silent-divergence generator.
4. **Every consolidation needs a gate** — byte-identity of results, or
   reproduce-the-certified-result. This generalises the rule already recorded
   for `PSR/lib`.

## 2. Extraction plan, ranked by risk

### Tier 0 — no behaviour can change — **DONE 2026-07-26**

| item | scale | gate | outcome |
|---|---|---|---|
| **`psr_vendor_check` drift test** | 15 identical + 4 known-divergent | the test *is* the gate | **DONE** — `PSR/tests/test_psr_vendor_drift.m` |
| **IPOPT options block** → `cr3bp_ipopt_opts` | **3 CR3BP solvers** (not 5 — see below) | byte-identity of the options struct | **DONE** — `cr3bp_common/cr3bp_ipopt_opts.m` |
| **CasADi path bootstrap** → one helper | 12 copies repo-wide | — | **REFUSED — see below** |

The drift test is deliberately listed *instead of* consolidating `PSR/lib`.
PSR's self-containment is a deliberate, verified design; the defect was never
the vendoring, it was that nothing noticed the drift. It found 15 byte-identical
files (the survey's count of 14 was one low) and was negative-tested by
injecting a line into `sms_pack.m`.

#### The IPOPT block was 3 solvers, not 5

The survey said "all 5 solvers; 41 identical lines." Measured before extracting:
the **three CR3BP-family Sundman solvers** (`casadi_minfuel_sundman`,
`casadi_energy_freetf`, `casadi_mintime_freetf`) share **20 identical option
assignments** — the only difference was one *redundant* `mu_strategy='monotone'`
in tulip, overridden by both arms of its own `warmTight` branch.

`earth/casadi_lt_mee` shares **zero**. It carries `mumps_pivot_order = 0`, the
AMD-ordering workaround for a hard METIS crash at its problem sizes. Forcing one
helper across both families would either impose that workaround where it is not
wanted or dilute it into a flag. The "41 lines" figure came from an
ELFO-vs-tulip comparison and was never true of the earth solver.

Gate: `cr3bp_common/tests/test_cr3bp_ipopt_opts.m` reconstructs the
pre-extraction inline struct verbatim and asserts `isequal` in both warm-start
regimes — which is why **no certified row had to be re-solved** to justify the
change. Also verified: the helper resolves under each campaign's own
`setup_paths` from `restoredefaultpath`, and IPOPT accepts every option name in
both regimes at a real `opti.solver` call.

#### The CasADi bootstrap consolidation is REFUSED

Listed as Tier 0 on the assumption it was free. It is not, and the survey did
not check the one thing that decides it: **no folder is on all four campaigns'
paths.** `cr3bp_common` serves tulip and ELFO only; earth reaches only its own
subfolders and cannot see it; PSR reaches nothing outside itself by design.

So a shared bootstrap helper would have to live somewhere new and be added to
every campaign's path — **creating exactly the kind of cross-cutting dependency
this survey concluded is the expensive mistake** — to remove three lines per
solver. Worse, the bootstrap is *inside each solver* rather than in
`setup_paths` for a reason: it makes a solver callable without campaign setup.
Hoisting it into a helper would move that dependency from "reads an env var" to
"needs another file on the path," which is strictly worse for the one property
the current form buys.

Left as is, deliberately. If it ever moves, it should move by the winning
pattern — into a file every campaign *already* loads — not into a new one.

### Tier 1 — small, one gate each

| item | why | gate |
|---|---|---|
| ~~**`insertion_states` → `cr3bp_common`**~~ **DONE 2026-07-26** | ELFO dropped the 34-file path dependency; **verified with `requiredFilesAndProducts`: all 7 ELFO entry points now reach 0 files in `GTO_tulip`** | `cr3bp_common/tests/test_insertion_states.m` — six endpoints pinned *and* checked against the certified seeds |

#### What the move turned up

The function was already a pure `cr3bp_common` client — all three of its callees
(`cr3bp_lt_params`, `gto_tulip_endpoints`, `gto_elfo_endpoints`) live there — so
both of its own `addpath` calls were wrong, and the 2026-07-21 reorg had made one
of them point at a directory that no longer exists. **Every ELFO endpoint call
had been emitting `Name is nonexistent or not a directory` since then.** Both
lines are gone; the gate asserts no path warning returns.

The old `sundman_minfuel/test_insertion_states.m` had been **broken since the
same reorg** — it loaded `'../elfo/results/energy_elfo_freetf.mat'`, a path that
stopped existing when ELFO moved — so it failed at load and had not run for five
days. Merged into the new gate with corrected paths; its seed comparison is the
half with real authority and is preserved.

#### OPEN DECISION — PSR runs the upstream solver, not its vendored copy

Found while doing the above, and **not** resolved here because it changes
behaviour on a certified pipeline.

`PSR/lib` was vendored 2026-07-12 and verified self-contained. On **2026-07-15**
the insertion-points feature added `addpath(../sundman_minfuel)` to all three PSR
entry points to reach `insertion_states`. `addpath` prepends, so that one line
put the upstream folder ahead of `PSR/lib`. Measured:

| path state | `casadi_minfuel_sundman` resolves to |
|---|---|
| `PSR/setup_paths` alone | `PSR/lib/` — as designed |
| after the entry-point `addpath` (what actually runs) | **upstream** |

So `PSR/lib/casadi_minfuel_sundman.m` and `PSR/lib/minfuel_at_tf.m` are **dead
code**, and PSR has been running the upstream solver for eleven days. The
vendoring did not fail loudly — an unrelated later feature silently overrode it.

`insertion_states` is now vendored into `PSR/lib` and the `addpath` is retained,
so **PSR's behaviour is unchanged by this work**. What remains is a genuine
choice:

- **(a) Restore isolation** — drop the `addpath`. The two dead files go live,
  i.e. PSR's solver silently reverts to the 2026-07-12 snapshot. Requires
  re-running the certified PSR result.
- **(b) Accept reality** — delete the two dead vendored files and document PSR
  as depending on upstream for the solver. No behaviour change, but PSR is no
  longer self-contained.

Do **not** take either as cleanup. `test_psr_vendor_drift` now pins the
resolution of all six load-bearing names, so this cannot flip again unnoticed
in either direction (negative-tested).
| ~~**certified-quantity guard** → one function~~ **DONE 2026-07-26** | **4 files / 5 instances**, not "8+" — `psr_mee_refine`'s hits were a refinement-loop stop-reason string | `verify_common/tests/test_certified_guard.m` (all 3 refusals, both directions) + 8/8 `verify_common` tests + 4 campaign smokes |
| **`optdef` vs per-campaign `getdef`/`gd`** | **REFUSED — see below** | — |
| ~~**Fix stale headers**~~ **DONE 2026-07-26** | ELFO `setup_paths` rewritten (the claimed `casadi_minfuel_sundman`/`minfuel_at_tf` reuse never happened); CR3BP README's "pipeline stages" corrected | none needed |

#### The certified-quantity guard: the divergence was not real

Two sites tested `out.success` alone, two also whitelisted the IPOPT status —
which looked like exactly the silent divergence this survey warns about. It was
not. `casadi_lt_mee` folds the whitelist *into* `out.success` itself, while the
Sundman solvers set `success = true` on any non-throwing solve, so their callers
had to apply it. All five sites were behaviourally equivalent; the shared guard
applies the whitelist unconditionally (no-op for the former, required for the
latter).

The one **real** difference it had to generalize over is which way is better:
min-fuel maximizes `m_f`, min-time minimizes `t_f`. `spec.better` is validated,
never defaulted — a mislabelled direction inverts the whole gate, refusing good
points and accepting degraded ones.

Fallout worth noting: this extraction is what surfaced that
`test_foc_check_10N` had been throwing `Unrecognized field name` since
`b3c6eaf` rather than gating anything. A consolidation that forces you to
actually run the affected tests finds things a read-through does not.

#### `optdef`/`getdef` consolidation is REFUSED

Measured before deciding: **27 definitions** across 4 names (`getdef` ×17,
`getfield_default` ×5, `local_default` ×2, `fcdef` ×1, plus the shared
`optdef`). Normalizing whitespace and comments, they collapse to **4 distinct
bodies that differ only in the parameter name** (`d` vs `dflt`):

```matlab
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
```

So there is **no divergence and no hazard** — nothing to fix. Consolidating
would delete ~54 trivial lines while adding a path dependency to 27 files,
including `casadi_energy_freetf` and `casadi_mintime_freetf`, which today have
none. That is the same trade refused for the CasADi bootstrap in Tier 0, and it
is refused here for the same reason: a local one-line subfunction buys
self-containment for free, and self-containment is the property these files
should keep.

The four *names* for one function are a readability wart, not a defect. If it
ever bothers anyone, rename in place — do not centralize.

### Tier 2 — real work, real gates

| item | why it is harder |
|---|---|
| **Two-pass seed protocol** — 3 copies, one carrying an explicit *"mirrored VERBATIM from run_transfer_mee.m lines 132-161"* comment | crosses the earth/CR3BP boundary; needs a reproduce gate on both 10 N rows |
| **CR3BP's double pipeline** — `run_cr3bp_geo` inlines stages A–D while `bridge_mu_continuation` + `solve_cr3bp_minfuel` implement the same sequence as callable stages | must decide which is authoritative before either can delegate |
| **`warmstart_on_mesh` (tulip) vs `interp_warmstart` (earth)** | two implementations of one idea in *different transcriptions* — needs genuine comparison, not a move |

### Never

- **Do not merge the solver cores.** `casadi_lt_mee` (MEE/L-domain) and
  `casadi_minfuel_sundman` / `casadi_*_freetf` (Cartesian/Sundman) differ in
  state dimension, independent variable, clock and horizon treatment. Merging
  buys tidiness and risks four campaigns' certified results.
- **Do not restructure `earth_elliptic_to_geo`.** It is the reference pattern.
- **Do not prune the 130 walled indirect files** (`ms_band` 57, `ztl` 47,
  `ifs` 26). Recorded status is "did not work," but three campaign records cite
  them as the evidence base. Understand the folder as an archive; leave it.

## 3. What this means for the other goals

- **Execution paths — DONE.** Four `doc/EXECUTION_PATHS.md` files.
- **Entry scripts + ladder scripts.** Do not design these: **three working
  precedents already exist.** `run_gergaud` (parameter block, run modes) and
  `run_cr3bp_geo` (thrust, lunar gain 0→1, ε, `t_f` in one editable section),
  with ladder shells in `run_cr3bp_ladder.sh`, `reproduce_table3.sh` and
  `elfo_{batch,energy_sweep,movies}.sh`. **Only `GTO_tulip` lacks a front
  door** — that is the whole of the remaining work for this goal.
- **LaTeX per campaign — last.** Template: `cr3bp_geo_phase1_note.tex`
  (reviewed, near-publishable). The `EXECUTION_PATHS.md` files are the raw
  material for the "how to run it" sections. Writing these before the tier-0/1
  cleanup means writing them twice.

## 4. Two patterns worth knowing before reading any survey number

- **75–87% of non-test files are unreached by any test** in every campaign.
  Mostly viz, one-off drivers and analysis scripts, which is defensible — but
  it means "it has a test" is not a safe assumption anywhere here.
- **Folder-local orphan detection over-reported in all four campaigns.** Every
  candidate turned out to be a doc-referenced figure script or a cross-campaign
  callee (`hamiltonian_const_check` is called by CR3BP; `cart_to_elements` by
  tests). Never delete on a single-folder graph; confirm repo-wide.
