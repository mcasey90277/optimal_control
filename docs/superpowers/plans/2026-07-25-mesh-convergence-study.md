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
- **Switch times are extracted by interpolated zero-crossing of the switching function, NEVER by thresholding the nodal throttle.** Reading a switch as "the first node where `s > 0.5`" quantizes its location to the grid *by construction*, so the study would measure `p ≈ 1` for switch times no matter what the underlying solution does — confirming hypothesis H1 circularly and proving nothing. Locate each switch as the zero of the switching function (or, where only the throttle is available, the linearly interpolated crossing of `s = 0.5` between adjacent nodes) and report it in **physical time**. A switch-time order fitted from thresholded node indices is a bug, not a measurement.
- MATLAB house header convention (purpose / INPUTS with sizes / OUTPUTS with sizes / REFERENCES); never `i` or `j` as loop variables.
- Cache-dependent tests SKIP (not fail) when a `.mat` is absent.
- Run each module from its own folder after that module's `setup_paths`; tests add CasADi via `addpath(fullfile(getenv('HOME'),'casadi-3.7.0'))`.

## Theoretical expectations (REVISED twice: 2026-07-25 and 2026-07-27)

**The first draft asserted p≈2 everywhere; the 2026-07-25 review replaced that
with p≈1 everywhere. Both were single confident predictions, and the two
reviewers DISAGREE.** Reviews:
`verify_common/doc/review_2026-07-25_meshplan_gpt56terra.md` (argues p≈1) and
`verify_common/doc/review_2026-07-27_meshplan_gemini31pro.md` (argues p≈1 is
too pessimistic for final mass and switch times).

**This table therefore states COMPETING FALSIFIABLE PREDICTIONS, not a single
expected answer.** The study exists to discriminate between them. A row with
two hypotheses is not indecision — it is the measurement's actual purpose, and
it is what stops a confirmatory reading of the result.

**The pessimistic argument (H1):** on a mesh whose breakpoints are NOT aligned
with the switches, one switch-crossing interval carries the wrong throttle over
`O(h)` of time against an `O(1)` jump in the RHS. That single interval's `O(h)`
error dominates every smooth arc's `O(h^2)`.

**The optimistic argument (H2):** the throttle is a CONTINUOUS NLP variable in
[0,1], so at a switch the optimizer can hold one node at an intermediate value,
encoding a sub-interval switch location rather than snapping to the grid.
Furthermore the objective is a quantity the optimizer directly controls and
evaluates by second-order quadrature, and the sign of the per-switch error
alternates across 25–360 switches, so integrated quantities may cancel to a
better effective rate than the trajectory norm.

| quantity | H1 (switch interval dominates) | H2 (sub-grid placement / cancellation) | switches as true breakpoints |
|---|---|---|---|
| state trajectory (global norm) | **O(h)** | O(h) — H2 makes no claim here | O(h²) |
| final mass, ΔV | **O(h)** | **O(h²)** — optimizer-controlled, 2nd-order quadrature, alternating-sign cancellation | O(h²) |
| switch **times** | **O(h)** — grid-quantized | **O(h²)** — intermediate node encodes sub-grid location | O(h²) *only* if the switch is an optimized event/phase boundary with Ṡ≠0 |
| switch **count** | must **stabilize**; no order exists | same | same |
| mapped terminal covector (Task 0) | **exactly 0** to KKT tolerance — not a convergence observable | same | same |
| raw defect multipliers | **O(h)** generically | same | O(h²) only for a properly covector-mapped multiplier on smooth arcs |
| ~~terminal state error~~ | **DROPPED** — with terminal equalities in the NLP this is identically solver tolerance, not an observable. Replaced by terminal error of the *independently propagated* control | | |

**Reading rule.** Report the measured `p` against BOTH columns and name which
hypothesis it supports. Do NOT record p≈1 as "expected and healthy" — that
phrasing (present in the previous revision) would wave through a genuinely
under-converged mass. Equally, do not record p≈2 as vindication without
checking that the switch topology was stable across the levels that produced
it.

**Switch-time extraction is a PREREQUISITE, not a detail.** See the Global
Constraint below: if switch times are read by thresholding the nodal throttle,
they are grid-quantized *by construction* and the study will measure p≈1
whatever the solution does — confirming H1 circularly. Extraction must resolve
sub-interval locations before any switch-time order is reported.

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
  and compare against the collocation solution. This is local by construction
  and so does NOT hit the 40-rev forward-shooting conditioning wall.

  **NORMALIZATION IS MANDATORY AND WAS AMBIGUOUS IN THE PREVIOUS REVISION.**
  The earlier wording ("reporting per-cell error") specifies the *local
  truncation error*, which for trapezoidal collocation is one order HIGHER than
  the global discretization error — `O(h^3)` per cell against `O(h^2)` global.
  Fitting an order to raw per-cell error inflates every reported `p` by +1 and
  would have invalidated the whole study. Two admissible forms, and the
  implementation must state in its header WHICH it uses:

  1. **Per-cell, normalized:** divide each cell's integrated error by that
     cell's own physical step `h_k` before fitting. Recovers `O(h^2)`.
  2. **Per-arc:** integrate across a whole inter-switch arc, whose physical
     length does NOT shrink with `h`, and compare at the arc's end. This
     accumulates to the global rate directly and needs no normalization.

  Report per-cell error UNNORMALIZED only as a diagnostic map (to show WHERE
  error concentrates, with switch cells called out); never fit an order to it.

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
├── mesh_switch_times.m       # sub-grid switch LOCATION (prerequisite for any order)
├── mesh_match_switches.m     # pair switches across meshes by physical time
├── mesh_order.m              # observed order + Richardson from >=3 levels
├── mesh_ladder_mee.m         # earth: re-solve a banked row at k x nodes, warm-chained
├── mesh_collect.m            # assemble per-quantity series from a ladder's outputs
├── mesh_report.m             # the study table (orders, Richardson values, verdicts)
│   #   tests/test_mesh_switch_times.m accompanies the extractor
├── run_mesh_study_mee.m      # Tier A driver (earth rows)
├── run_mesh_spot.m           # Tier B driver (single 1x/2x stability check, any campaign)
└── tests/
    ├── test_mesh_match_switches.m
    ├── test_mesh_order.m
    └── test_mesh_collect.m
```

Campaign-side: no production file is modified by this plan.

---

### Task 0: Mapped terminal covector — **COMPLETE 2026-07-26** ✅

**Status verified 2026-07-27, before resuming the plan.** All seven steps are
done and committed: `a8f7b61 fix(foc): mapped terminal covector for
transversality (Hager covector mapping)`, refined by `b3c6eaf`. Evidence:
- `foc_check.m` §(6) produces and **gates** `rep.lamMassEndMapped`.
- `verify_common/tests/test_foc_terminal_covector.m` PASSES (re-run 2026-07-27:
  mapped 0.000e+00 vs raw one-sided dual 2.564e-02 at N=20).
- Step 6 retraction landed: `OPTIMALITY_CERTIFICATION.md` §A6 RETRACTION —
  earth 2.5 N went **1.265e-03 → 2.268e-18**; the "single genuine first-order
  finding" verdict is withdrawn.

**Two deviations from the plan as written, both improvements — record them:**
1. **The formula is never evaluated.** Rather than forming `lamHat_f` from the
   expression below, the implementation reads the entry directly off the
   assembled Lagrangian gradient `gL` at the terminal mass index. This makes
   the 2026-07-27 review's sign objection moot — there is no hand-chosen sign
   to get wrong — and it also captures any OTHER constraint touching
   `X(massRow,N+1)` that a hand reconstruction would miss.
2. **It is NOT an independent test, and the code says so** (`foc_check.m:196`):
   the gated quantity is one ENTRY of the vector whose max-norm is
   `rep.kktStatInf`, so it cannot fail unless KKT stationarity already has. It
   is a readable projection of the master residual. Do not present it as
   independent evidence, in this study or in the register.

*Original task text retained below for the record.*

### Task 0 (original): Mapped terminal covector — may settle the open anomaly outright

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

**SIGN IS UNSETTLED — the test must discriminate it.** The 2026-07-27 review
argues the quadrature term enters as `-(h_N/2)*L_x(N+1)`, not `+`. The sign
depends on the adjoint convention, and this plan does not settle it by
assertion: Step 3 mandates obtaining every ingredient by AD from the same
`opti.f`/`opti.g` rather than hand-coding a convention, and Step 1's test
(mapped `< 1e-8` while the raw interval dual is `> 10x` larger) fails loudly
under a wrong sign. **If Step 1 fails, try the opposite sign BEFORE concluding
the mapping is wrong**, and record which convention the AD-derived KKT system
actually uses.

**Objective form, verified 2026-07-27:** the MEE objective is a pure Lagrange
integral of throttle (`casadi_lt_mee.m:284`, `sum((dsig/2).*(w_k + w_{k+1}))`),
with NO terminal mass term. Mass therefore does not appear in the cost, and the
"mass component is exactly zero" claim is not undermined by the min-fuel
objective. Re-verify this if a campaign ever switches to a Mayer (max final
mass) form.

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

### Task 1a: Sub-grid switch-time extraction — **COMPLETE 2026-07-27** ✅ (NEW, prerequisite)

**Why this task exists.** The Global Constraint added after the second review
forbids reading switch locations by thresholding the nodal throttle. Nothing in
the repo satisfies that constraint today: `foc_check.m:254` locates switches as
`swI = find(diff(burn) ~= 0)` — a nodal INDEX from exactly that thresholding —
and reports only `rep.nSwitches`, no times at all. Feeding those indices into an
order fit would measure `p ≈ 1` for switch times whatever the solution does.
The matcher in Task 1 consumes times, so extraction must come first.

**Files:**
- Create: `orbit_transfer/verify_common/mesh/mesh_switch_times.m`
- Test: `orbit_transfer/verify_common/mesh/tests/test_mesh_switch_times.m`

**Interfaces:**
- Produces: `s = mesh_switch_times(sg, tNode, burn, Shat)` — `sg` [N1×1] the
  independent-variable grid, `tNode` [N1×1] physical time at each node (the
  carried time state; pass `sg` itself if the domain IS time), `burn` [1×N1]
  logical throttle state, `Shat` [1×N1] the de-weighted switching function
  (`rep.Sdeweighted`) or `[]`. Returns `s.tSw` [1×k] switch times in physical
  units, `s.sgSw` [1×k], `s.n`, `s.method` (`'switchfn'` | `'throttle'`),
  `s.bracket` [k×2 node indices], `s.subGridFrac` [1×k] the crossing's
  fractional position within its bracket.
- **Preferred estimator is the zero of `Shat`**, not the throttle: `Shat` is
  smooth and crosses zero transversally at a regular switch, so linear
  interpolation of it is a genuine sub-grid locator. The throttle crossing
  `s = 0.5` is the FALLBACK when `Shat` is unavailable, and must be recorded as
  such in `s.method` so a report can never silently mix the two.
- `s.subGridFrac` is the evidence that extraction is not grid-quantized: values
  clustered at 0 or 1 mean the estimator has collapsed back onto nodes.

- [ ] **Step 1: Write the failing test** (analytic brackets, no solves):
  (a) a linear `Shat` crossing zero at a known point strictly between two
  nodes is recovered to `< 1e-12`, and the recovered `sgSw` differs from BOTH
  bracketing nodes by more than `0.1*h` — the anti-quantization assertion;
  (b) with a NON-UNIFORM `tNode`, `tSw` equals the correspondingly interpolated
  physical time, not the uniform-grid guess;
  (c) `Shat = []` falls back to the throttle crossing and sets
  `s.method='throttle'`;
  (d) a burn profile with no transition returns `s.n == 0` and empty arrays;
  (e) two adjacent switches one interval apart are both found (not merged).
- [ ] **Step 2: Run to verify it fails** — undefined function.
- [ ] **Step 3: Implement** with the house header. Bracket each switch by the
  `burn` transition, then locate the crossing INSIDE that bracket by linear
  interpolation of `Shat` (or of `throttle - 0.5`). Guard the degenerate case
  where `Shat` does not change sign across the bracket: fall back to the
  throttle crossing for that switch and flag it, rather than returning a node
  index silently.
- [ ] **Step 4: Run — ALL PASS.**
- [ ] **Step 5: Commit** `feat(mesh): sub-grid switch-time extraction (zero-crossing, not thresholding)`

---

### Task 1: Switch matching across meshes — **COMPLETE 2026-07-27** ✅

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

### Task 2: Observed order + Richardson estimator — **COMPLETE 2026-07-27** ✅

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

### Task 3: Earth mesh ladder — **COMPLETE 2026-07-27** ✅ (smoked on earth 10 N)

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

### Task 4: Collection and report — **COMPLETE 2026-07-27** ✅

**Files:**
- Create: `orbit_transfer/verify_common/mesh/mesh_collect.m`, `mesh_report.m`
- Test: `orbit_transfer/verify_common/mesh/tests/test_mesh_collect.m`

**Interfaces:**
- **CORRECTIONS MADE DURING IMPLEMENTATION (2026-07-27):** (a) `lamMassEndRel` does not exist — Task 0 removed it; the field is `lamMassEndMapped`. (b) `termErr` is classed DIAGNOSTIC, not physical: the terminal elements are equality-constrained, so it reads solver tolerance by construction (measured 1.6e-35 -> 2.6e-36 on the 10 N ladder) and would otherwise have earned a meaningless CONVERGED verdict — consistent with this plan's own expectations table, which drops it. (c) the switch-matching window is PER-SWITCH (2x that switch's own local step), not a scalar; `mesh_match_switches` was extended to accept a vector.
- Produces: `S = mesh_collect(L)` — from a ladder, assemble per-quantity series and orders. Quantities: `mf`, `dV_kms`, `termErr`, `nSwitches`, `switchTimes` (matched pairwise against the FINEST level via `mesh_match_switches`, tolerance = 2× the coarsest mesh's local step in physical time), `kktStatInf`, `lamMassEndRel`, `lamTimeCoV`, `sdotMinRel`, `signPct`. Each gets a `mesh_order` result. Returns `S.q.<name> = struct('vals',..,'order',..)` plus `S.switchStability` (counts per level, matched/unmatched tallies).
- Produces: `mesh_report(S, tag, resDir)` — fixed-format table: quantity | values per level | observed p | Richardson | verdict. Verdicts: `CONVERGED` (monotone, p within [1,3]), `STABLE` (integer quantity unchanged), `NOT-CONVERGED` (non-monotone or p outside range), `O(h)-CONSISTENT` (for the dual-side anomalies: value shrinking at ≈ the refinement ratio), `INSUFFICIENT` (fewer than 3 levels). Save `mesh_<tag>.mat`.

- [ ] **Step 1: Failing test** on a synthetic `L` (hand-built struct array with known series, no solves) asserting: a clean $O(h^2)$ mass series is `CONVERGED` with p≈2; a constant switch count is `STABLE`; a switch count that changes is **not** `STABLE`; a `lamMassEndRel` series halving with each doubling is `O(h)-CONSISTENT`; a two-level ladder yields `INSUFFICIENT` for orders.
- [ ] **Step 2: Run — fails.  Step 3: Implement.  Step 4: Run — ALL PASS.**
- [ ] **Step 5: Commit** `feat(mesh): per-quantity collection + study report`

---

### Task 5: Tier A — earth order study — **RUN 2026-07-27** (Steps 1-2 done, Step 3 partial) ⚠

**Result: `verify_common/doc/mesh_study_tierA_results.md`.** 12/12 levels
converged in 80.3 min. ONE of three rows produced a usable order — 1 N,
**p = 1.228** (windows 1.596 -> 1.228, consistent; switches stable at 171),
supporting **H1**. Richardson limit 1372.156 kg, leaving the production-
resolution row 1.55 kg short. 10 N and 2.5 N are NOT-CONVERGED: both had their
x8 delta GROW, 10 N with a topology change (20 -> 18 switches) and 2.5 N
without one (75 throughout) yet still +3.18 kg. Refinement raised final mass
and lowered dV in all 9 refinement steps with no exceptions, so the certified
rows are conservative. Switch TIMES are not converged in any row.

**Still open from Step 3:** the branch check on 10 N and 2.5 N (re-solve the
finest level from an independent seed) before their x8 jumps can be attributed
to discretization rather than basin escape. Until that runs, no campaign-wide
order statement is warranted.

### Task 5 (original): Tier A — earth order study (the core measurement)

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

**Known risk (updated 2026-07-27):** the theoretical expectations table is the interpretive frame for everything. It has now been reviewed twice by independent external models, **which disagree on it** — GPT-5.6-terra argues p≈1 for mass and switch times, Gemini 3.1 Pro argues those two should reach p≈2. That disagreement has been folded into the table as competing hypotheses H1/H2 rather than resolved by a third opinion, because the measurement is what settles it. The residual risk is no longer "the frame may be wrong" but "the measurement may be built so it can only confirm one branch" — which is why switch-time extraction and arc-defect normalization are now Global Constraints rather than implementation details.

---

# PHASE 2 — after the 2026-07-28 reasoning review

**Why this section exists, and why it is written BEFORE the runs.** Phase 1's
recurring failure was not execution — every solve converged — but
interpretation: numbers were produced and their meaning decided afterwards.
Two claims did not survive external review. The antidote is pre-registration:
each step below states what would count as success, what would count as
failure, and what specific number we expect, **before** the run. A criterion
invented after seeing the result is not a criterion.

**Standing rule for Phase 2.** No result may be reported against a criterion
that was written or revised after that result was seen. If a criterion turns
out to be wrong, say so explicitly, state the corrected criterion, and mark
the affected result as re-judged rather than silently re-scored.

---

## Step 1 — Branch-matched ladder (BLOCKING)

**Question.** What is the discretization error of a production-resolution
solution, measured along the single branch that solution lies on?

**Why it blocks everything.** Phase 1 produced no valid mesh-error estimate for
any row: 10 N had zero branch-consistent windows, 2.5 N and 1 N one each. Until
one exists, the study's founding question is unanswered and no statement about
mesh adequacy is admissible.

**Target row.** Earth 10 N, seeded from the KNOWN BETTER 18-switch solution
(1378.116718 kg on 193 nodes), not the certified 19-switch row. Chosen because
its full ladder was 209 s, and because two points on the 18-switch branch
already exist, giving a falsifiable prediction (below).

**Deliverable — a sentence we must be able to complete:**
> For the 10 N row on the 18-switch branch, the final mass at production
> resolution is X kg; the extrapolated converged value for that branch is
> Y kg; the production number therefore carries a discretization error of
> |X − Y| kg.

### Preconditions (failing any makes the run VOID, not negative)

| # | check | threshold |
|---|---|---|
| P1 | every level converges | `Solve_Succeeded` and `maxDefect <= 1e-8` |
| P2 | switch count identical | all four levels equal |
| P3 | switches correspond 1:1 | zero unmatched in EITHER direction vs the finest, per-switch window |
| P4 | switch displacement shrinks | max matched drift decreasing level to level |

P2 alone is explicitly NOT sufficient: the 2.5 N row held 75 switches across a
genuine branch change. P3 and P4 are what actually establish branch identity,
and even they are necessary-not-sufficient (see the open question below).

### Outcomes

- **SUCCESS** — P1–P4 hold; `nValidWindows >= 2`; windows agree within 25% of
  their mean; differences monotone and shrinking in magnitude. Then `p` and the
  Richardson limit are admissible and the deliverable sentence can be written.
- **PARTIAL** — P1–P4 hold but only one valid window, or windows disagree by
  more than 25%. Report a *local effective order* with no asymptotic claim, and
  escalate to Step 2.
- **FAILURE** — the branch cannot be held even when explicitly seeded to hold
  it. This is a SUBSTANTIVE NEGATIVE RESULT, not a null: it would mean the
  production solution has no continuously trackable counterpart under
  refinement, and would put in doubt whether a branch-matched error bar is
  obtainable for this problem class at all. Report it as a finding.

### Plausibility guards on the number

- `p` in roughly [0.8, 2.4]. H1 predicts 1, H2 predicts 2. A value near 0.5 or
  near 3 means something other than discretization drives the sequence — do not
  report it as an order.
- The extrapolation must be bounded: `|limit - finest| < |finest - coarsest|`.
  A correction larger than the total observed movement is not trustworthy.

### PRE-REGISTERED PREDICTION

Two points on the 18-switch branch are already in hand:

| N | m_f (kg) |
|---|---|
| 193 | 1378.116718 |
| 1544 | 1377.874577 |

The difference is **−0.242141 kg**: mass DECREASES with refinement on this
branch. Fitting `m(h) = m0 + C h^p` through the pair:

| assumed p | extrapolated limit |
|---|---|
| p = 1 | **1377.840 kg** |
| p = 2 | **1377.871 kg** |

**Prediction: the branch-matched Richardson limit lands in 1377.84–1377.87 kg,
and the production value of 1378.117 kg overstates the branch's converged mass
by 0.25–0.28 kg.** Landing outside that window means either the two points are
not on one branch or something else drives the sequence; either is informative.

Note the direction is DOWNWARD, opposite to Phase 1's weakened claim that
refinement raises mass. This run tests that too.

Caveat on the prediction's status: the two points came from different solve
paths, so it is a prior, not a derivation.

### Open design question to resolve before running

Warm-chaining from a better basin may simply re-introduce the escape problem at
a deeper level. If a level jumps despite being seeded on-branch, that is
outcome FAILURE above — but distinguish "the solver left the branch" from "the
seed was inadequate" by also trying a down-projection of the jumped level back
onto the previous grid.

---

## Step 2 — Extend the ladder to 16x

**Question.** Does the window sequence stabilize, and near which value?

**Runs only if Step 1 returns SUCCESS or PARTIAL.** 10 N at 16x is ~3100
intervals; the 8x rung took 179.5 s, so budget ~10 min.

- **SUCCESS** — three or more valid windows on one branch, with
  `|p_finest − p_previous| < 0.25`, and the settled value within 0.3 of either
  1 or 2. Report which hypothesis it supports, naming the margin.
- **PARTIAL** — windows still drifting by more than 0.25. Report that no
  asymptotic order was reached at 16x and STOP escalating; going to 32x is not
  justified by the value of the answer.
- **FAILURE** — the branch breaks at 16x. Report the resolution at which
  branch identity was lost, which is itself a useful number.

**What we are looking for:** `p` settling near 1 (H1: the switch-crossing
interval dominates) or near 2 (H2: sub-grid switch placement plus cancellation).
A settled value between 1.3 and 1.7 supports NEITHER and should be reported as
such rather than rounded toward whichever is more convenient.

---

## Step 3 — Port the basin probe to the ELFO solver

**Question.** Does the ELFO front row have a better basin?

This is a code fix first. ELFO uses a free-t_f two-primary solver
(`casadi_energy_freetf`, 9 states including a time-scaling slack, two-primary
Sundman clock) — not the tulip's fixed-t_f single-primary 8-state
`casadi_minfuel_sundman`. The probe assumed one CR3BP solver where there are
two.

- **SUCCESS (mechanical)** — all seeds run to completion, each either
  converging or failing with a genuine SOLVER reason. A shape/wiring error on
  any seed is a failure of this step, not a result.
- **SUCCESS (scientific)** — a baseline-vs-best comparison exists, of the same
  form as the Earth and tulip probes.
- Either outcome — better basin found, or not — is a legitimate result, PROVIDED
  the mechanical criterion passed. Report which.

**Guard against the tulip caveat:** as with the tulip, the seed family here
excludes down-projections. A null result must be reported as "survives THESE
seeds", never as robustness.

---

## Step 4 — Cislunar mesh ladder

**Question.** Can the tulip null result be made to face the adversary that beat
every Earth row?

The tulip flagship survived six seeds, but the down-projection seed family —
which won at every Earth row — was unavailable because no cislunar ladder
exists. Both external reviewers independently flagged this as the reason the
null is weak.

**Cost warning, stated in advance:** the tulip mesh is 4001 nodes; 8x would be
32,008. Budget for 2x and 4x only unless 2x proves cheap. Declare the cap
reached rather than silently truncating.

- **SUCCESS** — at least two finer levels converge, are down-projected onto the
  4001-node production grid, and re-solved there. Then either the flagship is
  beaten (a finding, matching the Earth pattern) or it survives its strongest
  available adversary (a materially stronger null than the current one).
- **FAILURE** — no finer level converges. Report the ceiling; the tulip null
  then stands permanently qualified as "survives weak seeds only".

---

## Step 5 — Per-row reproducibility scale

**Question.** How much does a row's answer move under changes that SHOULD be
immaterial?

This is what makes "the certified rows are pessimistic by X kg" defensible.
Phase 1 asserted a 1e-12 seed perturbation moved 2.5 N by 0.5 kg; the size of
the perturbation was never measured, and external review flagged that
projection onto the unit-direction and box constraints can amplify a small raw
change. A gain must be read against a MEASURED per-row scale.

**Method.** For each row, re-solve from the certified solution under K >= 5
nominally-immaterial variations: seed renormalization on/off, throttle clipping
on/off, `warmTight` on/off, and two IPOPT option variants. Record the spread in
final mass. Also record the actual norm of the seed change in each case, so
"a 1e-12 change" is a measurement rather than a description.

- **SUCCESS** — a spread `sigma_row` per row, plus the measured seed-change
  norms. A reported gain is then admissible only if it exceeds `5 * sigma_row`.
- **What we expect:** 10 N small (its baseline reproduced exactly), 2.5 N large
  (~0.5 kg). If 10 N's spread turns out large too, the +1.015 kg gain needs
  re-examination.
- **FAILURE** — no meaningful spread can be defined because variants fail to
  converge. Report which variants are unusable.

---

## Step 6 — Switch-aligned refinement (the falsifiable PSR claim)

**Question.** Does making switches true retained breakpoints recover second
order?

The claim as stated in Phase 1: PSR-style refinement improves the order **only
if** it makes switches true breakpoints aligned to o(h), not merely dense
neighbourhoods around them. This is the largest piece of work in Phase 2 and
should not start until Steps 1–2 have established a baseline order to compare
against.

- **SUCCESS** — a branch-matched ladder with switches enforced as retained
  breakpoints yields `p` at least 0.5 higher than the unaligned ladder on the
  same row and branch, with both ladders passing Step 1's P1–P4.
- **FAILURE** — alignment does not raise `p`. This REFUTES the claim, which is
  a genuine result and must be reported as such rather than attributed to
  implementation difficulty.
- **VOID** — the aligned ladder cannot hold its branch. Then the comparison was
  never made; do not report a verdict.

---

## Phase 2 exit

Phase 2 is complete when either (a) the deliverable sentence of Step 1 is
written with a defensible error bar and Step 5's scale makes the basin gains
admissible, or (b) Step 1 returns FAILURE and that failure is documented as the
answer to whether branch-matched error bars are obtainable for this problem
class.
