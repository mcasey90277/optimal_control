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

### Tier 2 — real work, real gates — **resolved 2026-07-26**

All three turned out to be **refusals with a smaller real fix inside them**. That
is the useful outcome of a Tier-2 item: the comparison is the deliverable, and
in two of three cases it found a defect the consolidation would have papered
over.

| item | outcome |
|---|---|
| Two-pass seed protocol | **not extracted**; shared the physics constant instead; false provenance claim corrected |
| CR3BP double pipeline | **front door declared authoritative**; not restructured |
| `warmstart_on_mesh` vs `interp_warmstart` | **not merged** — different requirements; **ported a robustness guard that fixes a NaN** |

#### Two-pass seed: share the constant, not the protocol

Three copies, but only two are close — `run_cr3bp_geo` deliberately uses a
different guard (relative error against the registry's certified rev count,
probe re-run every call). And the pair that claimed to be identical is not:
`bridge_mu_continuation`'s *"mirrored VERBATIM from run_transfer_mee.m lines
132-161"* is false. It has since gained a `table3_recipes` knob lookup and
`resume`-gated cache reads, while the earth driver takes knobs from `cfg` and
reads its cache unconditionally. (It also pointed at line numbers in another
file — a reference guaranteed to rot.) Comment corrected.

What they genuinely share is **four lines**: an `N = 50` probe, the rev-window
assert, and `N = round(nodesPerRev * nRev)`. Extracting those would need a
seven-argument helper — longer than the body, and it would drag each caller's
caching policy into a shared signature. **Refused.**

But the window and the N formula are *physics*, and silently diverging them
would give the two campaigns different admissible seeds. So the window alone
moved to `earth/direct/lib/mee_seed_rev_window.m`, which CR3BP already reaches.
A shared **definition**, not a shared algorithm — standing rule 2.

#### CR3BP double pipeline: the front door is authoritative

Both paths have certification history, which is why this was deferred. The
deciding evidence is provenance, not code:

- `run_cr3bp_ladder.sh` invokes **`run_cr3bp_geo`** — so the front door produced
  the certified 10 → 0.1 N ladder, the campaign's headline result.
- `bridge_mu_continuation` + `solve_cr3bp_minfuel` produced the Phase-1
  development artifacts (`gate1`, `gate2`, `gainwalk`) and are reached today
  only by `verify_cr3bp_pmp` and `compare_vs_2body`.

An equivalence check was attempted and is **not cheaply available**: the staged
path's banked artifacts stop at the energy/gain stages (`cr3bp_bridge_T10N_phi0_gate2`
carries `gain: 0` and m_f = 0.906449 as an *energy* solve) while the front door's
`cr3bp_T10N_phi0_fuel` is min-fuel (m_f = 0.918103). Different problems — the
difference is not a discrepancy, and comparing them would require a real
`solve_cr3bp_minfuel` run.

**Decision: `run_cr3bp_geo` is the authoritative pipeline**; the staged pair is
the Phase-1 development and analysis path. Making the front door *delegate* to
it is deliberately NOT done: the gate would be a full re-run of every certified
rung, and the benefit is removing a duplicated sequence in one campaign. The
option is preserved, and the cost is now written down.

#### The two warm-starts are not the same function — and comparing them found a bug

| | `interp_warmstart` (earth) | `warmstart_on_mesh` (tulip) |
|---|---|---|
| state | 7 (MEE) | 8 (Cartesian/Sundman) |
| original nodes | **resampled** | **copied verbatim** (no-resample discipline) |
| state / direction | linear | pchip |
| throttle | `nearest` | step-hold to the **left** |
| degenerate bracket | *(none — see below)* | norm floor + fallback |

Different requirements, not two solutions to one problem. **Not merged.**

The comparison earned its keep anyway. Earth's renormalization — itself added as
a "LATENT BUG FIX" — divides by a norm that collapses between nearly
**antipodal** directions, which is exactly what `beta` does across a throttle
switch on a coarse source mesh. Measured on the shipped function:

- exactly antipodal → `beta = [NaN NaN NaN]`, handed straight to IPOPT;
- near-antipodal → a confident-looking **unit vector whose direction is pure
  round-off**.

Both silent. Tulip guards both. The guard is now ported, with a regression test,
and is provably inert on well-conditioned input (`max(n,1e-12) == n` whenever
`n > 1e-12`).

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
  `elfo_{batch,energy_sweep,movies}.sh`. **DONE 2026-07-26** — `GTO_tulip` was
  the last campaign without one and now has `direct/sundman_minfuel/run_gto_tulip.m`
  plus `run_tulip_front.sh`.

  The tulip front door departs from the others in one respect, deliberately:
  its sweep script walks **t_f, not thrust**. Thrust is the dimension that works
  for the earth and CR3BP campaigns; for the tulip it is an open problem (the
  20 mN pilot failed against a fixed-τ_f topology wall), so `run_gto_tulip`
  *refuses* off-nominal thrust with that explanation rather than running a solve
  known to fail. It likewise refuses factors below ~1.12, where the ε=1 energy
  backbone itself will not converge. A front door whose knobs silently fail is
  worse than one that says which knobs are real.
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

## 5. Design checklist for a homotopy min-fuel OCP

Read this before building or seriously modifying a campaign. It is deliberately
*not* a code generator (see the note at the end): the transcription was never
the hard part, and every attempt to share the skeleton across campaigns
measured worse than sharing nothing.

### The four recurring decisions

**(a) Independent variable — put the mesh where the stiffness is.**

| campaign | choice | why |
|---|---|---|
| earth, CR3BP-GEO | true longitude `L` (MEE) | separates slow shape from the fast angle; survives hundreds of revolutions |
| tulip, ELFO | Sundman pseudo-time `τ`, `dt/dτ = r₁^1.5` | auto-concentrates nodes at perigee; tames the `1/r₁³` Hessian terms |

Uniform-in-time discretization fails both problems for the same reason: the
mesh you need at perigee is ruinous everywhere else.

**(b) Horizon — and the dense-column trap that decides how you get it.**

This is the decision that has cost the most, and there is one general rule
behind all three answers:

> **A free scalar coupled to every node produces a dense column in the KKT
> matrix, and the sparse factorization dies at the node counts these problems
> need.** Both working answers are the same trick: replicate the scalar per
> node and tie the copies with *local* continuity constraints. Same
> formulation, banded instead of arrowhead.

| campaign | answer | what it bought / cost |
|---|---|---|
| tulip | `τ_f` **fixed**; `t(τ_f) = t_f` as a terminal condition | banded KKT — but the winding number cannot grow, which is exactly why its thrust ladder hits a topology wall |
| ELFO | free time via slack **state** `cScale`, `dc/dτ = 0` | banded *and* revolutions can grow — its 20 mN rung certified where tulip's failed |
| earth, CR3BP-GEO | free `ΔL` scalar, then **`liftDL`** (per-node replication) | the naive scalar segfaulted MUMPS at ≥2.5 N; `liftDL` fixed it with no change of formulation |

**(c) Where the ε=1 seed comes from.** The homotopy needs a smooth root, and
getting one is problem-specific:

- **two-pass probe** (earth, CR3BP-GEO): cheap `N=50` solve reads the
  revolution count, then sample at full density.
- **banked energy backbone** (tulip): pre-solved min-energy solutions per `t_f`
  (available 1.12–1.95; below that the *backbone itself* will not converge).
- **gravity homotopy** (ELFO): four legs bringing lunar gravity up, because the
  target sits deep in the Moon's well.
- **μ-continuation** (CR3BP-GEO): solve two-body, then walk `gain: 0 → 1`.
- **neighbour seeding** (all): rescale a solved *bang-bang* solution to a nearby
  `t_f` and re-sharpen lightly. Far cheaper — this is how a front is walked.

**(d) Which dimension you can ladder in — decided by (b).**

You can walk a **thrust** ladder only if your formulation lets the revolution
count grow: a lower-thrust rung needs more revolutions.

- **earth** and **CR3BP-GEO** have walked full ladders, 10 N → 0.1 N (free `ΔL`).
- **ELFO** has one certified off-nominal rung (20 mN, `cScale`) — the
  formulation permits it, but a full ladder has not been demonstrated.
- **tulip** cannot: its 20 mN pilot was an honest failure against the fixed-`τ_f`
  topology wall, so it ladders in `t_f` instead. If a thrust ladder is a goal, decide (b) with that in mind — it
is not a knob you can add later without reformulating.

### Traps, each with the campaign that paid for it

1. **`opti.dual()` returns canonicalized-orientation duals.** Index
   `opti.lam_g` by the constraint's own row range instead. Cost: a 10°–60°
   primer misalignment that resisted explanation for nine days across two
   campaigns; the fix took it to 0.000° and the sign law to 100%.
2. **`maxIter` under-iteration silently collapses the deep-ε tail.** The run
   *looks* converged; the switching structure is simply unresolved. (earth deep
   rungs; CR3BP needs ≥4000.)
3. **`scaleNLP` is harmful** — it fights IPOPT's own automatic scaling. (earth)
4. **Renormalize interpolated thrust directions, but floor the norm.**
   Linear interpolation between near-antipodal unit vectors collapses toward
   zero: exactly antipodal gives `NaN` straight into IPOPT, near-antipodal gives
   a confident unit vector made of round-off. (earth `interp_warmstart`; guard
   ported from tulip's `warmstart_on_mesh`.)
5. **Certify only at the *requested* ε.** A clean intermediate step is a clean
   intermediate step, not a certified result. (tulip, triage C2)
6. **Never save an uncertified iterate** — it poisons neighbour-seed lookups
   later. (tulip `minfuel_at_tf`)
7. **Expect multiple basins.** Distinct seed routes at *identical* `t_f`
   converge to genuinely different extremals: measured ΔV spreads 0.02% to
   **9.84%**, in one case 17 switches apart, all ε=0 and machine-tight. Build
   the front as an envelope (`min` per factor), and treat a single-route sweep
   as an upper bound. (tulip)
8. **Fingerprint every cache boundary.** A nominal-thrust seed must not be
   served to an off-nominal request. (ELFO)
9. **Endpoints from exactly one function.** `insertion_states` exists so every
   pipeline starts and ends at literally the same numbers; a transfer to
   slightly-wrong endpoints still converges and still looks valid.
10. **The switch count is not a reproducible quantity; mass and ΔV are.** Three
    tulip re-solves of the certified recipe give 24 switches against the
    published 25, with mass agreeing to ~0.1%. Report it as a band.

### Why this is a checklist and not a skill

A code-generating skill would emit the easy 60% — the collocation loop and the
`s − ε·s(1−s)` objective — and leave every decision above to the user. It could
not be gated (what is the byte-identity check on a generated solver?), and it
would go stale silently: a generator written a week before this section would
today still be emitting the `opti.dual` bug from trap 1. Decisions and traps
transfer between campaigns; skeletons do not.
