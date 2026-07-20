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

### 2026-07-20 — Hamiltonian optimality verification + render densification
Independent **first-order (PMP)** verification of the certified rows, complementary to
SOSC. `verify/hamiltonian_const_check.m`: the time-costate `λ_t` (row-7 defect dual) is
**constant** along the trajectory (⇔ the time-domain Hamiltonian `H_t` is conserved),
because the L-domain problem is autonomous in the time state — measured **CoV ~1e-9
(0.1 N), 2e-8 (10 N)**; a mass-costate control varies **56%**, so the check
discriminates (not trivially satisfied). `verify/hamiltonian_along_traj.m`: reconstructs
`H_L(t)` and `H_t(t) = −λ_t` node-by-node. **Tricks worth keeping:** the dual **sign** is
pinned by transversality `H_t = dJ*/dt_f < 0` (verified on the textbook `min∫u²`: `J*=1/T`,
`H=−1/T²`) — the extremal identity does *not* discriminate sign; and that extremal identity
`dH_L/dσ = ∂H_L/∂σ` holds only to **collocation order** (~1.6% @26 nodes/rev — a
reconstruction check, not a precision claim). Movie `viz/hamiltonian_movie.m` (`H_L`
breathes / `H_t` flat, one shared scale). **Polygonal-trajectory fix:**
`viz/mee_res_to_cart_res.m` gained `nDense` (default 1 = byte-identical; `>1` pchip-densifies
the render) — the polygon was **8 nodes/rev + linear segments**, not physics. Also
re-confirmed **verbatim** from the HMG-2004 preprint (p.6–7,
`min_fuel_papers/Gergaud-Haberkorn-Martinon-JournalGuidance2004-preprint.pdf`) that
`t_f = c_tf·t_{fMin}, c_tf>1` is a **fixed**-time transfer ("not obvious that the minimum
fuel problem with free final time has a solution"; min-time `TfMin` solved first for the
feasibility floor) — validates our `c_tf = 1.5` convention.

### 2026-07-19 — Thrust-ladder external review + deep-thrust fixes (branch `ladder-deep-thrust`)
Three-way deep review (GPT-5.6-terra + Gemini 3.1 Pro + host) of the ladder core
(`kepler_lt_params`, `lt_mee_rhs`, `mee_seed`, `casadi_lt_mee`, `homotopy_mee`,
`run_mintime_mee`, `run_transfer_mee`, `run_ladder`, `psr_mee_refine`). Both
reviewers independently pinned the **single scalar `ΔL` as the root cause of the
MUMPS/METIS conditioning wall** (a scalar in every defect = a dense Jacobian
column). Verified + FIXED (all commits on `ladder-deep-thrust`, each inert at
feasible points — 10 N reproduces `m_f=1377.1012`, |Δ|=0):
- **`dL` bound rung-adaptive** (`dbc17e3`): the fixed `dL≤2000` (~318 rev) made
  0.2 N (ΔL~2168) and 0.1 N (~4335) STRUCTURALLY INFEASIBLE; now
  `max(2000, 5·dL0)`.
- **Guarded `Ldot` division** (`dbc17e3`): `1/Ldot` (rhs + objective) was
  unguarded → NaN on an IPOPT trial step with `Ldot≤0`; now `fmax(Ldot,1e-6)`
  (far below `LdotMin`, inert feasibly).
- **`opts.liftDL`** (`8a3c78c`, plumbed `9933383`): lift scalar `ΔL` → tied `N+1`
  sequence → block-banded KKT. Verified equivalent (byte-identical 10 N); max
  Jacobian column nnz **1546 → 19** (81× at N=193; ~1700× projected at 0.2 N).
Full reviews: `scratchpad/{gpt,gemini}_review.md`. Remaining reviewer
recommendations still open — see the deep-thrust item below.

### 2026-07-19 — SOSC certificate (branch `sosc-certificate`, `verify/sosc/`)
NLP-level second-order local-minimum certificate: warm-re-solve KKT recovery
(bounds from `opti.lbg/ubg`, per-kind dual feasibility), active-set classification,
and a **direct reduced-Hessian `eig(Z'HZ)` + `zt`-sensitivity** verdict (PASS /
WEAK_MIN / FAIL / INCONCLUSIVE). Opt-in gate in driver+reproducer (`cfg.certifySosc`,
default off; only FAIL demotes). **Finding:** 10 N certifies as **WEAK_MIN (270 flat
directions)** — min-fuel bang-bang extremals are weak, non-strict minima; strict SOSC
is generically unreachable. Batch (`recertify_table3`): 10 N WEAK_MIN; 5 N/2.5 N
INCONCLUSIVE (near-flat directions of unresolvable sign — corrected from spurious
FAILs of the first-cut `ldl`/Gould method); 1 N/0.5 N ERROR/uncomputable (recovery
fails at scale). Built via subagent-driven development (10 tasks + amendments, all
reviewed; whole-branch review Ready-to-merge). Design + method evolution:
`process/DESIGN_sosc.md` (§11–12), `process/PLAN_sosc.md`.

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

### 0. SOSC certificate — deep-rung scalability + deferred minors

**Files:** `verify/sosc/sosc_inertia.m`, `sosc_recover_kkt.m`, `recertify_table3.m`.

**What:** the certificate certifies 10 N cleanly (WEAK_MIN) but stops at the deep
rungs. Two distinct scale walls: (a) **recovery** — the warm re-solve itself fails
at 1 N/0.5 N (`n≈16.5k`); (b) **inertia** — `Z=null(full(A))` is a dense
allocation, guarded at `n≤maxNullDim=10000` (so 5/2.5 N run but 1/0.5 N
scale-skip). Also: 5 N/2.5 N return INCONCLUSIVE (near-flat reduced-Hessian
directions of numerically-unresolvable sign) — not a bug, a genuine precision
limit.

**Fix paths:** a scalable warm-resolve recovery for the deep rungs; a sparse
null-space / `eigs`-on-the-near-zero-cluster reduced-eig so 1/0.5 N need not
scale-skip. **Minors (non-blocking, from the whole-branch review):**
`recertify_table3` should persist an ERROR-verdict sidecar in its `catch` path (so a
hard MEX crash still leaves a per-rung record); `sosc_kkt_residual` two-sided-range
row hardening (fail-safe, provably absent today); drop the unused `K` param from
`sosc_active_set`; `datestr(now)`→`datetime` in `verify_sosc_mee`; header-style
consistency across `verify/sosc/`.

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

**0.2 N CERTIFIED (2026-07-20)** — `results/MEE_M2_0p2N.mat` (gitignored cache):
`m_f=1377.29 kg`, 823 switches, 346.7 rev, defect **2.5e-13**, termErr 7.5e-36,
incДeg 0, IPOPT Solve_Succeeded, ε=0 bang-bang (edge 99.9%). Never attained
before. **Recipe:** warm-chain from certified 0.5 N → rung-adaptive `dL` bound
(feasibility) + `opts.liftDL` (block-banded KKT, no crash at n~30k) +
**phase-correct β warm start** (`warmstart_phase_beta`, kills the σ-interp
aliasing) + **ε-continuation 1→0 with `maxIter=3000`** (`homotopy_mee`
`adaptiveEps`) — every ε-step converged tight (ok=1, defect ~2e-13). The
critical fix vs the first attempt was `maxIter` 1500→3000 (under-iterated steps
collapsed the tail); `adaptiveEps` bisection was armed but not triggered.
`scaleNLP` was tried and DROPPED (it fought IPOPT's gradient-based auto-scaling
→ ε=1 restoration failure; a proper *complete* user-scaling is a separate item).

**0.1 N CERTIFIED (2026-07-20)** — `results/MEE_M2_0p1N.mat`: `m_f=1377.29 kg`,
1644 switches, 693.6 rev, defect **5.0e-13**, termErr 0.00, incl 0°,
Solve_Succeeded, ε=0 (edge 99.9%). Reproduced by `reproduce_deep_rung(0.1,
'results/MEE_M2_0p2N.mat')` (warm-chained from 0.2 N, `maxIter=5000`, all 17
ε-steps ok=1) — the driver + recipe validated at the last rung. **THE FULL
10 → 0.1 N THRUST LADDER IS NOW CERTIFIED** (10/5/2.5/1/0.5/0.2/0.1 N; the deep
two were never attained before the external review). Only remaining open item on
this front is the proper complete user-scaling (deferred) and, optionally, PSR
switch-time refinement of the deep rungs.

**Progress (2026-07-19, external review):** the `dL≤2000` infeasibility and the
conditioning root cause are FIXED (rung-adaptive `dL` bound + `opts.liftDL`
block-banded KKT + guarded `Ldot`; see the Done entry). **Remaining
reviewer-recommended levers, in priority order:**
1. **Fix the β warm-start aliasing** (`interp_warmstart.m`): linear interp of β in
   σ keeps the SOURCE oscillation frequency; onto a finer-rev rung it phase-aliases
   (0.5 N→0.2 N is 2.5× revs). Reconstruct β by L-phase (evaluate at the target
   node's true longitude mod 2π), and rescale the warm time-row consistently
   (currently unrescaled on the warm/mintime path). *Likely the next wall after
   liftDL.*
2. **Explicit NLP scaling** (`casadi_lt_mee`): nondimensionalize `t` by `tfTarget`,
   `ΔL` by the C-law prediction, defect rows by expected state increments — both
   reviewers call gradient-based auto-scaling inadequate at these scales.
3. **Bypass the fragile min-time anchor**: continue jointly in (thrust, ε, t_f)
   from certified 1 N with smaller thrust ratios (1→0.8→0.65…) + adaptive-bisection
   homotopy (`homotopy_mee` currently jumps to the next ε on a failed step instead
   of bisecting).
4. **Architectural, if factorization stays the limit**: Sundman/time-domain
   transcription (precedent in `cartesian_legacy/`); HSL MA57/MA97 or Pardiso vs
   the AMD-workaround MUMPS; Schur-condensing/multiple-shooting; PSR `maxAdd`
   scaled to problem size + a genuinely PMP-steered (switching-function) refinement.
A first 0.2 N solve with the fixes + liftDL is being trialed (warm from 0.5 N).

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

- **(a) DONE (2026-07-19): code-tidy + subfolder/`lib/` reorg.** The formerly
  flat directory (~87 `.m` at top level) was split into functional subfolders
  (`core/ drivers/ psr/ verify/ frontdoor/ reproduce/ viz/ coords/
  cartesian_legacy/ lib/ tests/ attic/`) with generic helpers extracted to
  `lib/` and `setup_paths.m`/`module_root.m` handling the path + `results/`
  resolution; the working process docs moved to `process/` (see README "Code
  layout"). Pure reorganization + path fixes, certified numerics unchanged;
  verified by the 15-test no-solve sweep, a live IPOPT smoke, and a full
  from-scratch 10 N reproduce through the new shell path.
- **(b) Deeper warm-rung optimization.** The 10 N cold rung's multi-start
  spans a real seed-diversity sweep (multiple `seedThr`/`betaMode` pairs) that
  is how it found the better 18-switch basin. The warm-chained rungs
  (5/2.5/1/0.5 N) currently only explore the inherited warm-start plus
  `fuel_multistart`'s tiny `t_f`-bracket (`fuel_seed_set` returns a single
  inert candidate on the warm path) — they do not yet get the cold rung's
  seed-diversity treatment.
  **Empirical finding (2026-07-19, deep-ladder run):** this limitation now has
  teeth. The 2.5 N reproduce (warm-chained from the *improved* 18-switch 10→5 N
  basin) solved cleanly but topped out at **1368.13 kg — below the 1369.79 kg
  campaign floor**, so one-sided `verify_row` correctly rejected it. A better
  upper rung warm-starts the lower rung into a *different, worse* basin, and the
  warm-rung multi-start can't climb out — improvement does **not** monotonically
  chain down the ladder. Fix options under discussion: give the warm rungs a
  cold-seed diversity sweep, warm-start 2.5 N from a different source (e.g. the
  campaign's 5 N rather than the improved one), or widen the `t_f`-bracket.
  Blocks the 2.5/1/0.5 N deep-ladder validation until a recipe decision is made.
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
