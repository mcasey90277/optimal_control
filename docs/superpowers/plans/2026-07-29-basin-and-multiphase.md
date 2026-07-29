# Basin Sweep and Multi-Phase Build — specification

**Written 2026-07-29, BEFORE any run.** This study's two retractions both came
from deciding what a number meant after seeing it. Success and failure criteria
are therefore fixed here, in advance, and the standing rule from the mesh study
carries over: *no result may be reported against a criterion written or revised
after that result was seen.*

---

## Why these two, and why in this order

Three papers are planned: GTO→tulip (direct method + PSR), GTO→ELFO, and
elliptic→GEO in the CR3BP. Assessing the first against what exists produced one
**blocking** gap and one important one:

- **BLOCKING.** The tulip campaign has *two certified optima at
  t_f = 1.150×*, both machine-tight, both 25 switches, 1.43% apart in ΔV
  (m_f 0.847086 vs 0.849066). A paper cannot report one without knowing why the
  other exists. If a referee reproduces and lands on the other, the result is
  worse than not publishing.
- **IMPORTANT.** There is no accuracy statement beyond the defect — and we now
  know a 1e-14 defect does not imply an accurate trajectory, because the
  continuous residual is ~7e11 larger.

Part A addresses the blocking gap and is cheap. Part B addresses the accuracy
gap, the branch-identity failure that invalidated two measurement attempts, and
the seeds the stalled indirect campaign lacked — all with one construction.

Part A does **not** depend on Part B. Run A first; it may change what B needs.

---

# PART A — Basin sweep

## What we are doing

Systematically searching for better local optima of already-certified rows, on
each row's **own unmodified mesh**, and recording what is found.

## Why

Every certified earth row examined is beaten on its own mesh — by 1.015, 5.971
and 1.454 kg — with the winner also carrying lower ΔV. At 10 N a *one-node
throttle shift* suffices. The tulip independently records the two-optima problem
above. "Certified" currently means converged, feasible and first-order optimal.
It does not mean best-found, and nothing in the certification path looks.

This is an operational gap, not a research question.

## Scope

| row | why |
|---|---|
| tulip flagship at t_f = 1.150× | **the blocking paper gap** — resolve which optimum is right |
| tulip, both routes at 1.150× | the two known optima, probed from each other's solution |
| earth 10 N / 2.5 N / 1 N | already probed; re-run under the final protocol for the record |
| ELFO 1.33× front row | needs the 9-state adapter first (see A4) |

## Success

- **A1 (the blocking gap).** For t_f = 1.150× we can state which of the two
  optima is better and whether either is beaten by a third. Concretely: a table
  of every seed tried, its converged m_f, ΔV and switch count, with a single
  best identified. **The paper then reports that one, with the sweep as
  evidence it was sought.**
- **A2.** A documented, repeatable probe protocol — the seed families, the
  acceptance test (converged + machine-tight defect), and the comparison rule
  (higher m_f wins, one-sided, per the campaign convention).
- **A3.** Every row's best-found result recorded with the seed that produced it,
  so any of it is reproducible.
- **A4.** The ELFO adapter exists and runs (currently it calls the wrong
  solver — the campaign uses a free-t_f two-primary nine-state formulation, not
  the tulip's fixed-t_f single-primary eight-state one).

## Failure, and what each failure means

- **The sweep finds nothing better anywhere.** That is a *good* outcome, not a
  null: it means the certified rows are already best-found under these seed
  families, and the papers can say so. It must be reported with the seed
  families listed, because a null under weak seeds is not robustness — the
  tulip flagship already survived six seeds that excluded the family which beat
  every earth row.
- **The sweep keeps finding better optima without converging on a best.** Then
  the basin structure is richer than a fixed seed set can survey, and the papers
  must report a *range* rather than a value. This is a real possible outcome and
  would be a finding in its own right.
- **A1 cannot be resolved** — the two 1.150× optima persist with neither
  dominating and no third found. Then the paper must present both and explain
  the bifurcation, which is honest but weaker.

## Explicitly NOT claimed by Part A

Global optimality. Every one of these is a local optimum; a sweep raises the
floor and never proves a ceiling. Language in the papers must say
**"best found under the seed families listed"**, never "the optimum".

## Pre-registered expectations

1. The 1.150× pair resolves in favour of **m_f = 0.849066** (the higher mass),
   since minimum fuel maximizes final mass and both are machine-tight.
2. A seed derived from one optimum will *not* reach the other, since they are
   distinct basins — if it does, they are the same basin and the "two optima"
   record is wrong.
3. At least one earth row will be beaten again by a seed not yet tried, since
   the current best came from a family (down-projection) discovered by accident.

---

# PART B — Fixed-structure multi-phase formulation

## What we are doing

Reformulating the min-fuel problem as a multi-phase NLP in which the bang-bang
structure is **fixed by construction** and the switch times are **explicit
decision variables**, rather than emerging from a throttle that the solver may
place anywhere.

## Why — four problems, one construction

1. **Branch identity becomes structural.** Both retractions in the mesh study
   trace to being unable to prove two solutions are the same solution. With the
   switch sequence fixed, they are the same by definition, and mesh convergence
   becomes measurable for the first time.
2. **It supplies a genuinely switch-aligned mesh.** The toy validation
   established the mechanism: a true retained breakpoint is a *duplicated* node,
   which is exactly what a phase boundary provides. This is what the PSR /
   second-order claim has been trying to test indirectly throughout.
3. **It is the bridge to indirect methods, not a detour.** Switch times as
   explicit variables *is* the indirect structure — the same switching
   conditions PMP imposes — but solved by a well-conditioned NLP instead of by
   shooting. It produces exactly the seeds the stalled indirect campaign
   lacked (it stalled at ‖R‖~0.023 on terminal-cluster conditioning).
4. **Every external reviewer pointed at it independently**, across two review
   rounds and two different questions.

## Build order

- **B1.** Toy first: the double integrator, whose extremal is closed-form
  (t_s = 2−√2). A two-phase formulation must recover t_s to integrator
  tolerance. This gates everything, exactly as P3.1 did for the residual
  checker.
- **B2.** Earth 10 N, 19 switches → 20 phases. Chosen because it is the
  cheapest real row (~3 s/solve), it has the most surrounding data, and its
  basin structure is already characterized.
- **B3.** Mesh convergence ON the multi-phase formulation — the measurement
  that has failed twice, now with branch identity structural.
- **B4.** Tulip flagship, 25 switches → 26 phases. This is the paper row.

## Success

- **B1.** Switch time recovered to <1e-10 against the analytic value; final
  mass matching the single-phase solution to solver tolerance.
- **B2.** A converged 20-phase solution whose switch times agree with the
  single-phase solution's sub-grid switch times to within the single-phase
  mesh's local step, and whose m_f is **not worse** (it should be equal or
  better — the phase formulation can place switches exactly).
- **B3.** Three or more mesh levels on a *verified identical* structure, giving
  an observed order with mutually consistent sliding windows. **This is the
  first defensible convergence result the study would have.**
- **B4.** The tulip flagship re-solved in phase form, giving the paper an exact
  switch-time table rather than grid-quantized estimates.

## Failure, and what each means

- **B1 fails** → the formulation or its implementation is wrong; nothing
  downstream is interpretable. Stop and fix.
- **B2 fails to converge** → the multi-phase NLP is harder to solve than the
  single-phase one at this scale. Plausible: 20 phases means 20 extra
  variables and 20 interior-point boundary conditions. Report the failure mode
  (which phase, what residual) rather than concluding the approach is unusable.
- **B2 converges but to a worse m_f** → the fixed structure is not the right
  structure, i.e. the single-phase solution's apparent switch count is not the
  true one. That would be a significant finding about the certified rows.
- **B3 still cannot hold a structure across meshes** → then branch identity was
  never the obstacle, and the difficulty lies elsewhere. This would refute the
  reasoning behind Part B, which is worth knowing.

## Explicitly NOT claimed by Part B

- That the fixed structure is the *optimal* structure. It is fixed to what the
  single-phase solution found; if that structure is wrong, the phase solution
  inherits the error. Part A is what searches structures; Part B refines one.
- Second-order sufficiency. Everything here remains first-order.

## Pre-registered expectations

1. B2's m_f will be **equal or slightly better** than the single-phase value,
   because exact switch placement removes the one-interval error the residual
   map localizes.
2. B3 will give an order **closer to 2 than to 1**, since switch-aligned
   meshes remove the O(h) switch-crossing term — this is the H1/H2 question,
   asked properly for the first time.
3. B2 will be **harder to converge** than the single-phase problem, and may
   need a continuation in the number of freed switch times.

---

## What connects this to the three papers

| paper | needs from here |
|---|---|
| GTO→tulip | **A1** (resolve the two optima — blocking), **B4** (exact switch times), and a tulip PMP-residual adapter for the accuracy section |
| GTO→ELFO | **A4** (the adapter), then a basin statement of the same form |
| elliptic→GEO CR3BP | already externally reviewed; a basin statement would strengthen it |

The GTO→tulip paper's remaining gaps beyond this plan: PSR's value is not yet
quantified (a before/after on a stated metric is needed if PSR is a claimed
contribution), and the 1.01–1.11× band is unsolved — publishable as an open
finding, since it fails for the smooth energy problem too and is therefore
conditioning rather than bang-bang structure.
