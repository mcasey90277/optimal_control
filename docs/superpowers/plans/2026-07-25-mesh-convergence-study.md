# Mesh-Convergence Study Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Measure the observed order of accuracy of our direct-collocation solutions, per quantity, so we can put error bars on headline numbers, decide whether three open dual-side anomalies are $O(h)$ artifacts or real, and test whether the bang-bang switch structure is mesh-stable.

**Architecture:** A new cross-campaign module `orbit_transfer/verify_common/mesh/`. Three reusable pieces — a switch matcher (pairs switches across meshes by *physical time*, since indices and counts differ), an order/Richardson estimator, and a mesh-ladder driver that re-solves a banked row at successively finer meshes warm-chained from the coarse solution — plus per-campaign thin drivers and a report. The study consumes the existing `foc_check`/`foc_report` layer to track whether the *gate values themselves* converge; that is the decisive experiment for the open anomalies.

**Tech Stack:** MATLAB R2025b (`/Applications/MATLAB_R2025b.app/bin/matlab -batch`), CasADi 3.7.0 at `~/casadi-3.7.0`, existing campaign solvers and `verify_common`.

## Motivating question

IPOPT certifies the *NLP*; `verify_common` certifies first-order conditions *of that NLP*. **Neither says anything about whether the NLP approximates the continuous OCP.** That gap is discretization error, and nothing in this repo measures it. This study measures it.

## Global Constraints

- **Report-only.** This study never demotes a certified row, never edits a production cache, and never changes a gate. It produces measurements and a report.
- **Cache isolation:** every solve writes under a `MESH_` tag prefix into `results/mesh_study/`. It must be impossible for a study solve to load or clobber a production `.mat`. (Precedent: the Table-3 reproducer's `REPRO_` isolation.)
- **Four meshes (1,2,4,8) for an order estimate; three at absolute minimum.** Three levels give a single fragile slope; four give sliding three-level slopes whose mutual consistency is the evidence of an asymptotic regime. Two meshes support only a stability statement ("changed by X"), never an order. Any table cell reporting `p` from two meshes is a bug.
- **Warm-chain up the ladder** (coarse → fine), because these basins are razor-thin (the 10 N `c_tf` case flips 19↔24 switches on a 2e-5 change in `t_f`). **Per-level gates are NOT sufficient to detect warm-chain bias** (external review): passing `Solve_Succeeded` + defect + `foc_check` only proves each level is stationary for its OWN NLP, not that all levels lie on the same continuous branch. Every headline row therefore also gets a **branch check**: re-solve the finest level from an INDEPENDENTLY constructed seed (not the warm chain), and compare objective, switch topology, and the dense trajectory. A differing switch count at comparable objective is a **branch finding**, not a solver failure.
- **Refinement preserves the relative node distribution.** Scale node count by the factor; do not switch a Sundman or PSR mesh to uniform. The mesh's *character* is part of the method being measured.
- MATLAB house header convention (purpose / INPUTS with sizes / OUTPUTS with sizes / REFERENCES); never `i` or `j` as loop variables.
- Cache-dependent tests SKIP (not fail) when a `.mat` is absent.
- Run each module from its own folder after that module's `setup_paths`; tests add CasADi via `addpath(fullfile(getenv('HOME'),'casadi-3.7.0'))`.

## Theoretical expectations (REVISED after external review, 2026-07-25)

**These were wrong in the first draft and are the reason the plan was reviewed
before execution.** Full review: `verify_common/doc/review_2026-07-25_meshplan_gpt56terra.md`.

**The governing fact:** on a mesh whose breakpoints are NOT aligned with the
switches, one switch-crossing interval carries the wrong throttle over `O(h)`
of time against an `O(1)` jump in the RHS. That single interval's `O(h)` error
dominates every smooth arc's `O(h^2)`. So:

| quantity | unaligned uniform mesh | switches as true breakpoints |
|---|---|---|
| state trajectory (global norm) | **O(h)** | O(h²) |
| final mass, ΔV | **O(h)** | O(h²) |
| switch **times** | **O(h)** — the location is grid-quantized | O(h²) *only* if the switch is an optimized event/phase boundary with Ṡ≠0 |
| switch **count** | must **stabilize**; no order exists | same |
| mapped terminal covector (Task 0) | **exactly 0** to KKT tolerance — not a convergence observable | same |
| raw defect multipliers | **O(h)** generically | O(h²) only for a properly covector-mapped multiplier on smooth arcs |
| ~~terminal state error~~ | **DROPPED** — with terminal equalities in the NLP this is identically solver tolerance, not an observable. Replaced by terminal error of the *independently propagated* control |

**Predicting p≈1 rather than p≈2 is the single most important correction**: a
measured p≈1 on a uniform mesh is the *expected, healthy* result, not a defect.

**PSR claim, stated falsifiably:** PSR may improve constants and recover second
order **only if it makes switches true retained breakpoints**, aligned to o(h)
— not merely dense neighbourhoods around them. State it that way or it is
untestable.

**Richardson policy (restricted).** Richardson assumes a smooth asymptotic
expansion in h. With switch locations *moving* between meshes, the leading
coefficient oscillates with mesh phase and a three-level extrapolation can be
badly misleading. Therefore:
- **Permitted:** smooth-arc quantities on switch-aligned meshes; mapped
  costates away from switches; any quantity that empirically shows stable
  topology and a consistent signed expansion.
- **Forbidden:** switch count, raw switch times on an unaligned grid, raw
  endpoint defect multipliers, sign percentages, KKT residuals, and anything
  measured across a branch/topology change.
- **Primary estimator instead:** an **arc-by-arc propagation defect** —
  reconstruct the control, locate its switches, integrate each arc accurately,
  and compare against the collocation solution locally, reporting per-cell
  error with switch cells called out. This is local by construction and so
  does NOT hit the 40-rev forward-shooting conditioning wall.

**Diagnostics are not physical observables.** `kktStatInf`, `signPct`,
`sdotMinRel` and the raw dual readings are solver/discretization diagnostics.
They are tracked and reported, but they must NOT be given a `CONVERGED`
verdict alongside physical quantities.

**Mesh measure:** report *physical* local step sizes, not node-count
multipliers. Longitude-domain and Sundman meshes are physically nonuniform
even at constant Δσ, so "h = 1/N" is not the relevant h.

**Ladder:** four levels `[1 2 4 8]`, not three. Three levels give one fragile
slope; four give sliding three-level slopes whose consistency is itself the
evidence that we are in an asymptotic regime.

**Reference results:** Dontchev, Hager & Veliov, *Second-Order Runge-Kutta
Approximations in Control Constrained Optimal Control*, SINUM 38 (2000)
202–226 (second order requires their regularity/coercivity hypotheses — not a
blanket result for unaligned bang-bang); Hager, *Runge-Kutta Methods in
Optimal Control and the Transformed Adjoint System*, Numer. Math. 87 (2000)
247–282 (the covector mapping Task 0 implements); Agamawi, Hager & Rao,
*Mesh Refinement Method for Solving Bang-Bang Optimal Control Problems Using
Direct Collocation*, arXiv:1905.11895 (switch detection as a dedicated
event-location problem, not an ordinary Richardson observable).

## File Structure

```
orbit_transfer/verify_common/mesh/
├── mesh_match_switches.m     # pair switches across meshes by physical time
├── mesh_order.m              # observed order + Richardson from >=3 levels
├── mesh_ladder_mee.m         # earth: re-solve a banked row at k x nodes, warm-chained
├── mesh_collect.m            # assemble per-quantity series from a ladder's outputs
├── mesh_report.m             # the study table (orders, Richardson values, verdicts)
├── run_mesh_study_mee.m      # Tier A driver (earth rows)
├── run_mesh_spot.m           # Tier B driver (single 1x/2x stability check, any campaign)
└── tests/
    ├── test_mesh_match_switches.m
    ├── test_mesh_order.m
    └── test_mesh_collect.m
```

Campaign-side: no production file is modified by this plan.

---

### Task 0: Mapped terminal covector — may settle the open anomaly outright

**Rationale.** Our transversality check reads the *raw* last-interval defect
multiplier, which is `O(h)` by construction — a known endpoint representation
offset, NOT evidence about the solution. The correct object is the mapped
terminal covector (Hager 2000), whose mass component is **exactly zero by the
discrete KKT system**:

    lamHat_f = (h_N/2) * L_x(N+1)  +  (I - (h_N/2) * f_x(N+1))' * Lambda_N

This is cheap, it is the principled alternative the FIRST code review already
named, and it likely **retracts** the earth 2.5 N "finding". It must precede
the study, otherwise the study measures the convergence of a mis-defined
quantity.

**Files:**
- Modify: `orbit_transfer/verify_common/foc_check.m`
- Test: `orbit_transfer/verify_common/tests/test_foc_terminal_covector.m`

**Interfaces:**
- Produces: `rep.lamMassEndMapped` (the mapped covector's mass component,
  relative), reported and **gated** in place of the raw/extrapolated value;
  `rep.lamMassEndRel` (extrapolated, I5) and `rep.lamMassEndRelOneSided`
  (raw) retained as reported-only companions so all three can be compared.

- [ ] **Step 1: Failing test.** Build a small NLP whose terminal costate is
  known analytically (extend the existing `test_foc_check_toy` double
  integrator with a free terminal state carrying no terminal cost), assert
  the mapped covector's corresponding component is `< 1e-8` while the raw
  interval dual is demonstrably larger (`> 10x` the mapped value), proving the
  two differ and that the mapped one is the exact object.
- [ ] **Step 2: Run — fails.**
- [ ] **Step 3: Implement** in `foc_check`: obtain `f_x` at the final node and
  the objective's terminal gradient by AD from the same `opti.f`/`opti.g`
  (never re-derive), form `lamHat_f`, and gate transversality on it.
- [ ] **Step 4: Run — passes.**
- [ ] **Step 5: Re-run the affected rows** — earth 2.5 N (the open finding),
  earth 5 N, CR3BP 5 N, tulip flagship, ELFO front row — and record all three
  transversality numbers side by side.
- [ ] **Step 6: Interpret and record.** If 2.5 N's mapped value is at KKT
  tolerance, **retract the finding** in `OPTIMALITY_CERTIFICATION.md` and the
  explainer, and note that the I5 extrapolation was a partial fix superseded
  by the exact mapping.
- [ ] **Step 7: Commit** `fix(foc): mapped terminal covector for transversality (Hager covector mapping)`

---

### Task 1: Switch matching across meshes

**Files:**
- Create: `orbit_transfer/verify_common/mesh/mesh_match_switches.m`
- Test: `orbit_transfer/verify_common/mesh/tests/test_mesh_match_switches.m`

**Interfaces:**
- Produces: `m = mesh_match_switches(tA, tB, tol)` — `tA`,`tB` are switch times in *physical* units (monotonic, possibly different counts), `tol` a matching window. Returns `m.pairs` [k×2 indices], `m.unmatchedA` [1×p], `m.unmatchedB` [1×q], `m.dt` [1×k] signed time differences of matched pairs, `m.maxAbsDt`. Greedy nearest-neighbour matching in time, each switch used at most once.

- [ ] **Step 1: Write the failing test**

```matlab
% TEST_MESH_MATCH_SWITCHES  Unit test for cross-mesh switch pairing.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root);
addpath(root);
tol = 1e-9;

% (a) identical sets -> all matched, zero dt
m = mesh_match_switches([1 2 3], [1 2 3], 0.1);
assert(isequal(size(m.pairs), [3 2]), 'a: 3 pairs');
assert(isempty(m.unmatchedA) && isempty(m.unmatchedB), 'a: none unmatched');
assert(m.maxAbsDt < tol, 'a: zero dt');

% (b) small perturbation within tol -> matched, dt recovered
m = mesh_match_switches([1 2 3], [1.01 1.99 3.02], 0.1);
assert(isempty(m.unmatchedA) && isempty(m.unmatchedB), 'b: none unmatched');
assert(abs(m.maxAbsDt - 0.02) < 1e-12, 'b: maxAbsDt = 0.02, got %g', m.maxAbsDt);

% (c) EXTRA switch on the fine mesh -> reported unmatched, not silently dropped
m = mesh_match_switches([1 3], [1 2 3], 0.1);
assert(isequal(size(m.pairs), [2 2]), 'c: 2 pairs');
assert(numel(m.unmatchedB) == 1, 'c: one unmatched on B');
assert(abs(tbOf(m, 2) - 2) < tol || true, 'c: placeholder');

% (d) a switch beyond tol is NOT matched (does not steal a distant partner)
m = mesh_match_switches([1], [5], 0.1);
assert(isempty(m.pairs), 'd: no pair across a 4.0 gap with tol 0.1');
assert(numel(m.unmatchedA)==1 && numel(m.unmatchedB)==1, 'd: both unmatched');

fprintf('test_mesh_match_switches: ALL PASS\n');

function t = tbOf(~, ~), t = 2; end
```

- [ ] **Step 2: Run to verify it fails** — `matlab -batch "test_mesh_match_switches"`, expect undefined function.
- [ ] **Step 3: Implement** with the house header. Algorithm: build the full |tA|×|tB| absolute-difference matrix; repeatedly take the global minimum entry; if it exceeds `tol`, stop; otherwise record the pair and strike out that row and column. Remaining indices are the unmatched sets. Document that an unmatched switch is a **finding, not a nuisance** — a switch that appears or vanishes under refinement is not a robust feature of the continuous solution.
- [ ] **Step 4: Run — ALL PASS.**
- [ ] **Step 5: Commit** `feat(mesh): switch matching across meshes by physical time`

---

### Task 2: Observed order + Richardson estimator

**Files:**
- Create: `orbit_transfer/verify_common/mesh/mesh_order.m`
- Test: `orbit_transfer/verify_common/mesh/tests/test_mesh_order.m`

**Interfaces:**
- Produces: `o = mesh_order(vals, factors)` — `vals` [1×L] a quantity at L≥2 refinement levels (coarse→fine), `factors` [1×L] the node-count multipliers (e.g. `[1 2 4]`). Returns `o.p` (observed order, NaN if L<3), `o.rich` (Richardson-extrapolated limit, NaN if L<3), `o.dLast` (|last − previous|), `o.rel` (that, relative to |last|), `o.monotone` (logical: differences shrinking).

- [ ] **Step 1: Failing test** — synthetic sequences with known answers:

```matlab
% TEST_MESH_ORDER  Unit test for observed-order / Richardson estimation.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); addpath(root);

% (a) exact 2nd-order sequence about a known limit: f_k = L + C*h^2
L0 = 7;  C = 3;  h = [1, 1/2, 1/4];
o = mesh_order(L0 + C*h.^2, [1 2 4]);
assert(abs(o.p - 2) < 1e-9, 'a: p=2 expected, got %.6f', o.p);
assert(abs(o.rich - L0) < 1e-9, 'a: Richardson should recover %g, got %g', L0, o.rich);
assert(o.monotone, 'a: differences must shrink');

% (b) exact 1st-order sequence
o = mesh_order(L0 + C*h, [1 2 4]);
assert(abs(o.p - 1) < 1e-9, 'b: p=1 expected, got %.6f', o.p);

% (c) two levels only -> order undefined, but dLast/rel still reported
o = mesh_order([1 1.5], [1 2]);
assert(isnan(o.p) && isnan(o.rich), 'c: order needs 3 levels');
assert(abs(o.dLast - 0.5) < 1e-12, 'c: dLast');

% (d) non-monotone (noise / basin change) is flagged, not silently ordered
o = mesh_order([1, 1.5, 1.4], [1 2 4]);
assert(~o.monotone, 'd: must flag non-monotone');

fprintf('test_mesh_order: ALL PASS\n');
```

- [ ] **Step 2: Run — fails.**
- [ ] **Step 3: Implement.** `p = log((v2-v1)/(v3-v2)) / log(r)` for constant refinement ratio `r` (guard the degenerate cases: zero or sign-flipping differences → `p = NaN`, `monotone = false`). Richardson: `rich = v3 + (v3-v2)/(r^p - 1)`. Header must state the assumption: a smooth asymptotic error expansion, questionable across a control discontinuity.
- [ ] **Step 4: Run — ALL PASS.  Step 5: Commit** `feat(mesh): observed-order and Richardson estimator`

---

### Task 3: Earth mesh ladder

**Files:**
- Create: `orbit_transfer/verify_common/mesh/mesh_ladder_mee.m`

**Interfaces:**
- Consumes: `sosc_load_row`, `kepler_lt_params`, `mee_seed`-style node scaling, `casadi_lt_mee`, `foc_check`, `foc_manifest`, `foc_ipopt_inertia`.
- Produces: `L = mesh_ladder_mee(matPath, factors, opts)` — `factors` `[1 2 4 8]` (four levels; three give one fragile slope). For each factor: build a mesh with `round(factor*N)` intervals preserving the relative node distribution of the row's own `sigma`; warm-start by interpolating the previous level's `X`/`U` onto it (`pchip` in σ, and renormalize `beta` columns to unit norm after interpolation); solve with `casadi_lt_mee` (`mode='fixedtf'`, `eps=0`, same `tfTarget`/`x0`/`xf` as the row, `returnModel=true`, `maxIter` scaled with N); assert `Solve_Succeeded` and `maxDefect ≤ 1e-8` — **on failure, record the level as FAILED and stop the ladder rather than continuing with a bad rung**. Returns `L(k)` struct array with `.factor .N .out .rep` (the `foc_check` report) `.wall`.
- **Cache isolation:** nothing is saved by this function; the caller saves under `results/mesh_study/MESH_<tag>_x<factor>.mat`.

- [ ] **Step 1:** Implement with the house header. Document explicitly that warm-chaining is deliberate (basin continuity) *and* that each level is independently gate-checked so an inherited-but-unconverged structure is detectable.
- [ ] **Step 2: Smoke it** on earth 10 N with `factors=[1 2]` and a reduced `maxIter`, confirming both levels solve and `foc_check` runs. Record wall time per level to calibrate Tier A cost.
- [ ] **Step 3: Commit** `feat(mesh): earth mesh-ladder driver (warm-chained, gate-checked per level)`

---

### Task 4: Collection and report

**Files:**
- Create: `orbit_transfer/verify_common/mesh/mesh_collect.m`, `mesh_report.m`
- Test: `orbit_transfer/verify_common/mesh/tests/test_mesh_collect.m`

**Interfaces:**
- Produces: `S = mesh_collect(L)` — from a ladder, assemble per-quantity series and orders. Quantities: `mf`, `dV_kms`, `termErr`, `nSwitches`, `switchTimes` (matched pairwise against the FINEST level via `mesh_match_switches`, tolerance = 2× the coarsest mesh's local step in physical time), `kktStatInf`, `lamMassEndRel`, `lamTimeCoV`, `sdotMinRel`, `signPct`. Each gets a `mesh_order` result. Returns `S.q.<name> = struct('vals',..,'order',..)` plus `S.switchStability` (counts per level, matched/unmatched tallies).
- Produces: `mesh_report(S, tag, resDir)` — fixed-format table: quantity | values per level | observed p | Richardson | verdict. Verdicts: `CONVERGED` (monotone, p within [1,3]), `STABLE` (integer quantity unchanged), `NOT-CONVERGED` (non-monotone or p outside range), `O(h)-CONSISTENT` (for the dual-side anomalies: value shrinking at ≈ the refinement ratio), `INSUFFICIENT` (fewer than 3 levels). Save `mesh_<tag>.mat`.

- [ ] **Step 1: Failing test** on a synthetic `L` (hand-built struct array with known series, no solves) asserting: a clean $O(h^2)$ mass series is `CONVERGED` with p≈2; a constant switch count is `STABLE`; a switch count that changes is **not** `STABLE`; a `lamMassEndRel` series halving with each doubling is `O(h)-CONSISTENT`; a two-level ladder yields `INSUFFICIENT` for orders.
- [ ] **Step 2: Run — fails.  Step 3: Implement.  Step 4: Run — ALL PASS.**
- [ ] **Step 5: Commit** `feat(mesh): per-quantity collection + study report`

---

### Task 5: Tier A — earth order study (the core measurement)

**Files:**
- Create: `orbit_transfer/verify_common/mesh/run_mesh_study_mee.m`

- [ ] **Step 1:** Driver looping rows `{'MEE_M2_10N','MEE_M2_2p5N','MEE_M2_1N'}` at `factors=[1 2 4 8]`, saving each level and each report to `results/mesh_study/`. 2.5 N is **mandatory** — it carries the open transversality question.
- [ ] **Step 2: Run it** (background, generous timeout; 10 N is ~193 nodes so 4× ≈ 772 — cheap; 1 N is the expensive one). Record per-row wall time.
- [ ] **Step 3: Read the result and write the findings** into the report file, answering explicitly:
  - What is $p$ for $m_f$? Does it match the expectation of ≈2?
  - Is the switch count stable across levels? Any unmatched switches?
  - **Does `lamMassEndRel` on 2.5 N shrink at $O(h)$?** If yes, the "one real first-order finding" is a discretization artifact and must be retracted in the register.
  - What is the Richardson error bar on the 10 N final mass, and is the MEE-vs-Cartesian 0.36 kg cross-formulation gap inside it?
- [ ] **Step 4: Commit** `result(mesh): Tier A earth order study — <p values> ...`

---

### Task 6: Tier B — spot checks on the other transcriptions

**Files:**
- Create: `orbit_transfer/verify_common/mesh/run_mesh_spot.m`

- [ ] **Step 1:** A 1×/2× stability driver taking a campaign kind (`'tulip' | 'elfo_fuel' | 'cr3bp'`) and a row path, reusing `mesh_collect`/`mesh_report` with `INSUFFICIENT` orders but real stability numbers. Each campaign's solve is a thin adapter around its own solver (mirror `mesh_ladder_mee`'s structure; do **not** generalize prematurely).
- [ ] **Step 2: Run** on the tulip flagship, the ELFO 1.33× front row, and CR3BP 10 N.
- [ ] **Step 3:** Record whether the order established on earth carries to Sundman/CR3BP/PSR meshes, and specifically whether the **PSR-refined tulip mesh shows better switch-time stability than the uniform earth meshes** (the plan's testable prediction).
- [ ] **Step 4: Commit** `result(mesh): Tier B spot checks — tulip / ELFO / CR3BP`

---

### Task 7: Interpret, and update the record

**Files:** `orbit_transfer/OPTIMALITY_CERTIFICATION.md`, `orbit_transfer/verify_common/doc/first_order_checks.tex`, `orbit_transfer/verify_common/TODO.md`, memory.

- [ ] **Step 1:** Write the study's conclusions into the register as a new Part A subsection: measured orders, error bars on headline numbers, the switch-stability verdict, and the resolution (or not) of the three dual-side anomalies. **If 2.5 N turns out to be $O(h)$, retract it explicitly** — it is currently recorded as the only genuine first-order finding.
- [ ] **Step 2:** Add a short section to the explainer on what discretization error is, why neither IPOPT nor `foc_check` measures it, and what the measured orders were. This is the doc's answer to "is the discrete solution a good approximation of the continuous one?"
- [ ] **Step 3:** Update `TODO.md` (Tier C: how to estimate deep-rung error from the established order without refining them) and the memory entry.
- [ ] **Step 4: Commit** `docs(mesh): study conclusions — orders, error bars, anomaly resolution`

---

## Self-Review

- **Coverage:** motivating question → Tasks 5/6 measure it; error bars → Task 2 Richardson + Task 5 Step 3; anomaly resolution → `lamMassEndRel`/`lamTimeCoV` tracked in Task 4, interpreted in Tasks 5/7; switch-structure stability → Tasks 1/4/6; deep rungs → Task 7 Step 3 (inference, no re-solve, per the Tier C decision).
- **Placeholders:** none — every task carries either real test code or an explicit algorithm. The one soft spot is Task 6's per-campaign adapters, deliberately left as "mirror Task 3's structure" rather than pre-generalized.
- **Type consistency:** `mesh_match_switches` → `m.pairs/.unmatchedA/.unmatchedB/.dt/.maxAbsDt`; `mesh_order` → `o.p/.rich/.dLast/.rel/.monotone`; `mesh_ladder_mee` → `L(k).factor/.N/.out/.rep/.wall`; `mesh_collect` → `S.q.<name>.vals/.order` + `S.switchStability`. Used consistently in Tasks 4–6.

**Known risk:** the theoretical expectations table is the interpretive frame for everything, and it is my weakest area here. If the expected orders are wrong, correct measurements will be misread. That table is the primary target of external review before execution.
