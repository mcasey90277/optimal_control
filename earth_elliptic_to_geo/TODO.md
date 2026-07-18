# earth_elliptic_to_geo — Outstanding Work

Companion to `README.md`. Same structure as `proj7/pipelines/crlb/TODO.md`
(Done / Open ranked by priority / Not-a-goal). Source of record for the numbers
and provenance below: `CAMPAIGN.md` and `.superpowers/sdd/progress.md`.

---

## Done

### 2026-07-16..18 — MEE thrust-ladder campaign (commits `5f839dd..567801b`)
Rebuilt the solver in Modified Equinoctial Elements with `ΔL` a scalar decision
variable and `L` the independent variable; certified the full 10/5/2.5/1 N
min-fuel ladder (the Cartesian/Sundman stack died at 5 N) plus 0.5 N via an
anchor-free R0-law path. Cross-formulation gate passed (MEE 10 N m_f = 1377.10
vs Cartesian 1376.74 kg). PSR switch-refinement ported. PMP verifier delivered.
Full record + 6 binding footnotes in `CAMPAIGN.md`.

### 2026-07-18 — `run_gergaud` front door (commits `c04a057..a5b568b`)
Single PARAMETERS-block entry with user-definable initial AND final orbits,
three run modes (auto/solve/probe), Table-3 row printout, trajectory plot, and
movie. Endpoint parameterization (`opts.xf` terminal, `opts.initElems` initial)
threaded through the solver stack, default-preserving (existing certified caches
still load). Four rendered movies (`results/movie_MEE_{10N,5N,2p5N,1N}`). Feature
test suite 6/6. Built via subagent-driven development; final whole-branch review
Ready-to-merge after two custom-path fixes (I-1 PSR-skip guard, I-2 note
visibility).

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
`opti.dual` (suspected incomplete `opti.dual` un-scaling). See `DESIGN_dual_map.md`.

### 5. Map the full Fig-23 front (multiple c_tf per thrust)

**Files:** `run_gergaud.m`, `run_transfer_mee.m`.

**What:** the campaign only ever solved **one** `c_tf = 1.5` per thrust level.
The paper's Fig 23 overlays several `c_tf` curves. `run_gergaud`/`run_transfer_mee`
already accept `ctf`, so this is a sweep, not new machinery.

**Fix:** loop `c_tf ∈ {1.2,1.5,2.0,2.5,3.0}` per thrust, collect `m_f(c_tf)`, and
plot the multi-curve front. Watch the basin scatter documented in `CAMPAIGN.md`
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
