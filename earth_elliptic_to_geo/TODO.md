# earth_elliptic_to_geo — Outstanding Work

Companion to `README.md`. Same structure as `proj7/pipelines/crlb/TODO.md`
(Done / Open ranked by priority / Not-a-goal). Source of record for the numbers
and provenance below: `process/CAMPAIGN.md` and `.superpowers/sdd/progress.md`.

---

## Reproducing Table 3

**"Reproduce" has two meanings — pick the right tool:**

**(A) Present the numbers we already certified** (the table/figure as reported).
Both the row cache and the campaign's own generators read the same certified
`.mat` files; use the generators:
```matlab
run_ladder([10 5 2.5 1])   % prints the ladder summary from the rung caches (no re-solve)
fig_table3                 % the Table-3 figure (switches/revs/R0 panels)
```
`run_gergaud(struct('thrustN',T,'runMode','auto'))` gives the same numbers one
row at a time (plus a plot/movie), and is the right tool for a single row or a
custom-endpoint row — not for the whole table.

**(B) Re-solve from scratch.** Do **what the campaign did** — a warm-chained
thrust continuation with per-rung recipes — **not** `run_gergaud solve`. The
per-rung recipe and the reasons `run_gergaud` will not reproduce the deep ladder
are documented in **`doc/campaign_reproduction_runbook.tex`**. Summary:

- **10 N** — `run_transfer_mee` (seed δ≈0.4, ε:1→0 homotopy) + `run_mintime_mee(10,25)`
  (multi-basin: cold 3-rev seed wins over the fuel-warm spurious basin; keep-best).
- **5 / 2.5 N** — `run_ladder` warm-chain (each rung C-law-rescaled from the one
  above); 2.5 N needs the raised per-round budget `cfg.mtMaxIter≥150–300` (the
  "wall" was a checkpoint-budget artifact, not conditioning).
- **1 N** — the hand recipe: min-time anchor **small-N-first** (15/rev, manual
  continuation) → **warm mesh-refine to 25/rev** (the size that crashed cold);
  fuel needs the tfTarget-relative time box; then PSR round 2 for the headline
  1371.44 kg.
- **0.5 N** — **anchor-free R0-law** target `t_f=c_tf·(223.14/T)` (its min-time
  anchor hit a conditioning wall) + coarse fuel solve + 5 PSR rounds.
- **0.2 / 0.1 N** — not attained; blocked on the 0.5 N min-time wall (item 1 below).

**Why `run_gergaud solve` will not reproduce (B):** it is single-rung (no
warm-chain across thrusts), it does not encode the 1 N small-N-first or 0.5 N
anchor-free recipes (it would call `run_mintime_mee(T,25)` and walk into the
crash wall), and its solve/probe path was wired but never run-to-convergence on
the deep rungs. See the runbook for the full argument.

---

## Done

### 2026-07-16..18 — MEE thrust-ladder campaign (commits `5f839dd..567801b`)
Rebuilt the solver in Modified Equinoctial Elements with `ΔL` a scalar decision
variable and `L` the independent variable; certified the full 10/5/2.5/1 N
min-fuel ladder (the Cartesian/Sundman stack died at 5 N) plus 0.5 N via an
anchor-free R0-law path. Cross-formulation gate passed (MEE 10 N m_f = 1377.10
vs Cartesian 1376.74 kg). PSR switch-refinement ported. PMP verifier delivered.
Full record + 6 binding footnotes in `process/CAMPAIGN.md`.

### 2026-07-18 — `run_gergaud` front door (commits `c04a057..a5b568b`)
Single PARAMETERS-block entry with user-definable initial AND final orbits,
three run modes (auto/solve/probe), Table-3 row printout, trajectory plot, and
movie. Endpoint parameterization (`opts.xf` terminal, `opts.initElems` initial)
threaded through the solver stack, default-preserving (existing certified caches
still load). Four rendered movies (`results/movie_MEE_{10N,5N,2p5N,1N}`). Feature
test suite 6/6. Built via subagent-driven development; final whole-branch review
Ready-to-merge after two custom-path fixes (I-1 PSR-skip guard, I-2 note
visibility).

### 2026-07-18 — Table-3 reproducer engine (best-found optimizer)
Built `reproduce_row(T)` — a Table-3 reproducer **engine**, not a bit-exact
replay: it re-solves one rung's min-fuel transfer entirely from scratch
(`REPRO_`-tag isolation so it can never load or clobber the campaign's own
production caches), composing the certified drivers (`run_mintime_mee`,
`run_transfer_mee`, `psr_mee_refine`, `anchor_smallN_first`) via a per-rung
recipe (`table3_recipes.m`: `coldB`/`chain`/`smallN_first`/`R0law` anchor
strategies). Because minimum-fuel means maximize final mass, the fuel stage
is a **keep-best-mass multi-start** (razor-sensitive fuel basin — a ~2e-5
change in `t_f` can flip the switch structure); `verify_row.m` is
**one-sided**, throwing only if the reproduced mass falls below the
campaign floor (`table3_certified.m`) — a higher mass always passes and is
flagged an improvement, with structure (switches/revs) reported but not
gated. Three entry points: `reproduce_table3.sh` (per-process watchdog,
crash/hang relaunch), `reproduce_table3.m` (thin in-process wrapper),
`reproduce_row(T)` (single rung); `reproduce_table3_collect.m` prints the
updated (best-found) Table 3 with a per-rung floor comparison. **Live proof
at 10 N: the engine found 18 switches / 7.56 rev / 1378.46 kg, beating the
campaign's certified 1377.10 kg by 1.36 kg** (and closer to the paper's
~18-switch structure) — the campaign had under-optimized that rung; Table 3
is updated with best-found numbers as the ladder is re-run. 1 N harvested a
`smallN_first` recipe from the hand campaign; 0.2/0.1 N recipes are
registered (`table3_recipes.m`, `.seeded=true`) but not yet run. Full detail
in README.md "Reproducing from scratch (best-found)". Supersedes
`run_task9_rung.m` (now a deprecated shim — see Open item 7 below and the
file's own DEPRECATED header).

---

## Open — ranked by priority

### 1. Certify a 0.5 N min-time anchor (conditioning wall)

**Files:** `run_mintime_mee.m`, `casadi_lt_mee.m`.

**What:** the 0.5 N free-longitude min-time solve never converged — 7 configs
tried, best defect 0.0545 (12/rev × 75 nodes; 15/rev was worse), bit-identical
stall on retry, 4+1 reproducible `libcoinmumps` MEX/SIGBUS crashes. So the
0.5 N fuel row (`MEE_M2_0p5N_PSR_psr_final.mat`, m_f = 1375.28) is currently
built against an **R0-law tfmin ESTIMATE** (`t_{f,min} ≈ 446.27 ND` from
`T·t_{f,min} ≈ 223.14 ND`), not a certified anchor.

**Fix path:** small-N-first then mesh-refine (the recipe that cracked the
crash-prone 1 N anchor: 15/rev ~660 nodes → one warm refine to 25/rev), and/or
warm-chain the 0.5 N anchor from the certified 1 N anchor via the C-law
`ΔL`-rescale. If a certified anchor differs from 446.27 ND by >~1%, re-solve the
0.5 N fuel row against the new target.

### 2. Descend to 0.2 N and 0.1 N (never attained)

**Files:** `run_gergaud.m` (probe mode), `run_ladder.m`.

**What:** these Table-3 rungs were honestly never solved — the deep-descent
effort stopped at the 0.5 N min-time wall (item 1). `run_gergaud probe` wires
the live attempt and reports `certified=false` rather than faking a row.

**Fix path:** blocked on item 1 (chain a certified 0.5 N anchor down), or extend
the anchor-free R0-law path to 0.2/0.1 N with PSR at very large N. Expect ~300
and ~600+ switches, N in the tens of thousands, and the crash class of item 1.
`table3_recipes.m` already registers `chain` recipes for both rungs
(`.seeded=true`) so `reproduce_row(0.2)`/`reproduce_row(0.1)` are wired and
ready to run once item 1 unblocks the chain — see item 7 below.

### 3. Thread a custom terminal target through PSR

**File:** `psr_mee_refine.m` (`solve_psr_round`), `run_gergaud.m`.

**What:** `psr_mee_refine` and its internal `solve_psr_round` call
`casadi_lt_mee` without an `xf` field, so every refinement round re-terminates
at the default GEO target. `run_gergaud` therefore **skips PSR for custom
endpoints** (I-1 guard) and reports the un-refined fuel solve. That is honest and
correctly-targeted, but a custom-endpoint run at T ≤ 1 N gets no switch-time
sharpening.

**Fix:** add an `xf` field to `psr_mee_refine`'s options and forward it into
both `casadi_lt_mee` calls in `solve_psr_round`; then drop the `~isDefaultEndpoints`
guard in `run_gergaud`. Re-validate the certified 1 N / 0.5 N PSR results are
byte-unchanged for the default (GEO) target before removing the guard.

### 4. PMP dual/primer anomaly — the escalate-branch probe

**Files:** `verify_pmp_mee.m`, `mee_dual_to_costate.m`.

**What:** the first-order PMP gates fail (primer misalignment 10–60°,
eccentricity-correlated) because the **raw IPOPT duals** fail cone-elided KKT
stationarity at high eccentricity — proven not a verifier bug by an independent
KKT re-derivation. Primal certifications are unaffected (they never use the duals).

**Fix path (Campaign B):** recover the raw `lam_g` via `nlpsol` bypassing
`opti.dual` (suspected incomplete `opti.dual` un-scaling). See `process/DESIGN_dual_map.md`.

### 5. Map the full Fig-23 front (multiple c_tf per thrust)

**Files:** `run_gergaud.m`, `run_transfer_mee.m`.

**What:** the campaign only ever solved **one** `c_tf = 1.5` per thrust level.
The paper's Fig 23 overlays several `c_tf` curves. `run_gergaud`/`run_transfer_mee`
already accept `ctf`, so this is a sweep, not new machinery.

**Fix:** loop `c_tf ∈ {1.2,1.5,2.0,2.5,3.0}` per thrust, collect `m_f(c_tf)`, and
plot the multi-curve front. Watch the basin scatter documented in `process/CAMPAIGN.md`
(take the best certified point per `(T,c_tf)`).

### 6. Housekeeping minors (deferred from the front-door review)

Low-risk items logged during the subagent-driven build, none affecting certified
numbers:

- `casadi_lt_mee.m`: `assert(numel(xf)==5)` checks count not shape; a permanent
  `selftest` early-return hook sits in the production NLP builder. Consider
  factoring `xf`-resolution into a small helper if more test hooks appear.
- `mee_seed.m`: the explicit `initElems=[]` (empty-but-present) case is correct
  by inspection but untested — add a one-line assertion.
- `run_transfer_mee.m` / `run_mintime_mee.m`: `initElems` is fingerprinted as an
  `initElems_isset` boolean only (two different custom `initElems` collide under
  one tag). Mitigated for the supported entry point by `run_gergaud`'s
  endpoint-hash tag suffix; only a risk for direct driver calls with a fixed tag.
- `gergaud_plot.m` duplicates `transfer_movie.m`'s ring/axis/burn styling (forced
  by the no-touch-renderer constraint). Consolidate into a shared helper if a
  third consumer appears.
- `run_gergaud.m`: a dead `if ~anchorOut.certified` branch (`run_mintime_mee`
  throws instead of returning uncertified); the probe-mode warning text is not
  rung-conditional (fires the 0.2/0.1 N wall message even for a 10 N probe).

### 7. Reproducer-engine follow-ups

**Files:** whole `earth_elliptic_to_geo/` tree (a); `reproduce_row.m` /
`table3_recipes.m` (b, c); `reproduce_table3.sh` (d).

- **(a) Coming: code-tidy + generic `lib/` refactor.** The directory has grown
  flat and dense (100+ `.m` files at the top level: solver core, per-rung
  drivers, PSR machinery, viz, tests, and now the reproducer engine all
  side-by-side). A refactor into a shared, generic `lib/` (mirroring the
  `mlabTools` layout used elsewhere in this repo) is planned so the pipeline
  stays navigable — flagged here so it is not lost as a priority. Scope: pure
  reorganization + path fixes, not a rewrite of certified numerics.
- **(b) Deeper warm-rung optimization.** The 10 N cold rung's multi-start
  spans a real seed-diversity sweep (multiple `seedThr`/`betaMode` pairs) that
  is how it found the better 18-switch basin. The warm-chained rungs
  (5/2.5/1/0.5 N) currently only explore the inherited warm-start plus
  `fuel_multistart`'s tiny `t_f`-bracket (`fuel_seed_set` returns a single
  inert candidate on the warm path) — they do not yet get the cold rung's
  seed-diversity treatment. Worth exploring once the full ladder validation
  (Task 6) shows whether the warm rungs already meet the floor or need a
  wider search themselves.
- **(c) 0.2/0.1 N recipes are seeded, not run.** `table3_recipes.m` carries
  `chain` recipes for both (`.seeded=true`, warm-started 0.5→0.2→0.1 N), but
  neither has been executed to a certified `REPRO_row_T*.mat` in this build —
  see item 2 above (blocked on the 0.5 N min-time anchor wall, item 1).
- **(d) Calibrate the watchdog wall cap from real timing.** `reproduce_table3.sh`
  uses a flat, unscaled 21600 s (6 h) per-attempt wall cap for every rung
  because no rung has yet been run start-to-finish under this orchestrator
  (see the script's own header). Once Task 6's controller-run produces real
  per-rung wall-clock numbers, replace the flat guess with a thrust-scaled (or
  at least rung-specific) cap via `WALLCAP_S`.

---

## Not a goal — intentional scope boundaries

- **CR3BP / third-body / lunar gravity.** This pipeline is strictly two-body
  Earth-centered — that is the paper's problem. The CR3BP + Moon low-thrust work
  lives in `../NLP_lowThrust_GTO_tulip/` and is a separate campaign; do not add
  third-body terms to `lt_mee_rhs.m`.
- **Full SPICE ephemeris / perturbations (J2, drag, SRP).** The paper is a clean
  two-body + thrust benchmark; keep it that way unless a study demands otherwise.
- **An indirect (PMP shooting) solver here.** The paper solves indirectly; this
  project's contribution is the *direct* reproduction. The PMP machinery present
  (`verify_pmp_mee`) is a *verifier* of the direct solution, not a second solver.
- **Reviving the Cartesian/Sundman stack for the deep ladder.** It is retained
  only for the cross-formulation gate; the MEE + ΔL formulation is the production
  path for all thrust levels.
- **Routine live certification of 0.2/0.1 N.** Until item 1 is solved these are
  research probes, not a supported `auto`-mode row.
