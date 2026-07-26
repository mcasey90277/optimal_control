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

### Tier 0 — no behaviour can change (do these first)

| item | scale | gate |
|---|---|---|
| **`psr_vendor_check` drift test** — assert each of the 14 byte-identical `PSR/lib` files still matches its origin; fail loudly on drift | 14 files | the test *is* the gate |
| **CasADi path bootstrap** → one helper | **12 copies** repo-wide | byte-identity of results |
| **IPOPT options block** → `ipopt_opts(profile)` | **all 5 solvers**; 41 identical lines | byte-identity of results |

The drift test is deliberately listed *instead of* consolidating `PSR/lib`.
PSR's self-containment is a deliberate, verified design; the defect was never
the vendoring, it was that nothing noticed the drift.

### Tier 1 — small, one gate each

| item | why | gate |
|---|---|---|
| **`insertion_states` → `cr3bp_common`** | lets ELFO drop a 34-file path dependency it uses for one function, and removes a `casadi_minfuel_sundman` shadowing surface | ELFO + tulip smoke |
| **certified-quantity guard** → one function | **8+ copies** (`refresh_duals_mee`, `refresh_duals_cr3bp`, `run_foc_tulip`, `run_foc_elfo`, three earth drivers …) | byte-identity |
| **`optdef` vs per-campaign `getdef`/`gd`** | one utility, four local reimplementations | byte-identity |
| **Fix stale headers** — ELFO's `setup_paths` claims a reuse that does not happen; CR3BP's README calls functions "pipeline stages" the front door never calls | misleads the next reader | none needed |

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
