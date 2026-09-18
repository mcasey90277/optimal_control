# Cart-pole chains and the shared library — design (2026-09-17, revised)

This revises the 2026-09-17 three-objectives spec
(`2026-09-17-cartpole-three-objectives-design.md`) after sub-project A landed
and the goal widened. **The three problems, the bang-bang treatment (events plus
saltation), and the seeding rules in that spec stand unchanged.** This document
supersedes its *Layout*, *Study scripts*, *Sub-projects* and *Verification*
sections, and adds what the widened goal requires: chain scripts, the library
promotion program, a normality gate, and the orbit-transfer lessons as binding
constraints.

## The goal, as widened

Mike, 2026-09-17: *"One of the goals is to get good direct and indirect code for
this problem and good study scripts and to augment our oclib; I think that there
may be more generic library functions that we can pull out of the orbit transfer
code and put in the oclib folder. I want good cart-pole optimal control code and
good optimal control code for orbit transfers."*

So this is **one library serving two plants**, and both bodies of code are
deliverables. The cart-pole is not an example that borrows from the orbit work;
it is the **second plant that forces the library to be generic.** A function used
only by CR3BP campaigns can quietly assume seven states, a thrust/Isp/mass triple,
a `muStar` or a pumpkyn call, and no test will object. One that must also serve a
4-state plant with a scalar force cannot. The cart-pole work is therefore the
instrument that *finds* the couplings — which is also what makes `oclib`'s
admission rule (a second top-level consumer, an equivalence gate, a delegate left
behind) satisfiable at scale for the first time.

Two consequences run through everything below:

- **Both trees end green, always.** Every promotion touches code five orbit
  campaigns solve through. Its gates run *before* the move as well as after.
- **Promotions are pulled, never pushed.** A function moves when a cart-pole
  chain needs it, so it arrives with a live consumer and a live gate. Nothing is
  promoted speculatively. Some candidates will be refused after measurement, and
  a refusal is a result.

## Layout

```
optimal_control_examples/
  cartpole_common/            the plant (field, Jacobian + generator, constants),
                              the power-balance ORACLE, and -- new -- the
                              bang-bang propagator bangbang_prop (events +
                              saltation) once the min-time chain needs it
  ex3_cart_pole_minenergy/    was ex3_cart_pole_pmp; renamed in the docs
                              sub-project so the three siblings read as ONE
                              PLANT, THREE COSTS. Holds the min-energy chain,
                              the existing study script, its solver front doors
  ex4_cart_pole_mintime/      min-time chain + study + front doors
  ex5_cart_pole_minfuel/      min-fuel chain + study + front doors
  run_all_tests.m             one runner, one exit code (exists)
  teaching_docs/              revised once all three objectives have code
oclib/+oc/                    grows only by promotion; tests/ per function
orbit_transfer/               keeps a delegate at every old path; every campaign
                              keeps running bit-identically
```

The rename of ex3 is deferred to the last sub-project deliberately: it touches
READMEs, the teaching document and CLAUDE.md, and it is a `git mv` gated on the
suite reproducing its numbers — the same gate the plant extraction passed — so
it costs least when the tree is otherwise quiet.

## Chain scripts

One skeleton, three instances, in the discipline of
`orbit_transfer/DRO_tulip/indirect/build_70mN_library.m`: parameter block and
stage switches first; every stage's output is a file the next stage reads; a
stage that is off is skipped and its saved output reused, so a chain resumes
from any point; **no new machinery in the chain — each stage calls one front
door.** Each chain spans direct through indirect, because the direct solve's
multipliers *are* the costate seed and a chain cut at that seam would have to
define, version and test a file contract three times over.

```
0  parameters + stage switches   plant, objective, bound, horizon rule, budgets
1  prerequisites                 the plant oracle passes; metadata asserted
2  anchor                        a direct solve at an EASY point, warm-started
                                 from the min-energy answer (L5)
3  walk                          a direct homotopy toward the target problem
4  refine                        switch-aware mesh refinement (bang problems)
5  verify (direct)               first variation, true residual, mesh-density
                                 basin check, bound activity, normalised
                                 switch statistics
6  harvest                       duals -> costates; switch structure; then the
                                 SEED PROBE: fly the harvested costates and
                                 measure the miss before anyone shoots (L15)
7  shoot                         oc.ms_bvp, with bangbang_prop for bang problems
8  verify (indirect)             pointwise PMP checks, conjugate point,
                                 normality; the study script's gates, run
                                 headless
9  export + pictures             the .mat that downstream consumers read, with
                                 its metadata; the figures
```

Per objective:

| stage | min-energy | min-time | min-fuel |
|---|---|---|---|
| 2 anchor | the problem itself (it converges cold) | large `u_max`, seeded from the min-energy solution | `t_f = 1.5·t_min`, seeded from min-energy |
| 3 walk | dark | `u_max` ladder down to 40 N, each rung warm-seeding the next and banked; the `t_min(u_max)` front falls out | energy→fuel ladder: ε-relaxation AND a Huber-style law raced head to head (L6, L-race); the winner is a measurement |
| 4 refine | dark | detect sign changes in `u`, subdivide, re-solve until the switch times stop moving | same, on the three-branch control |
| 7 shoot | `oc.ms_bvp`, fixed `t_f` (exists) | 5 unknowns (`λ(0)`, `t_f`), closed by `H(t_f) = 0`; `bangbang_prop` | fixed `t_f`; `bangbang_prop` with the `|σ| = 1` event pair |

**Every homotopy step is band-gated** (L10): a rung whose cost or `t_f` lands
outside `[1/3, 3]×` what *that step's own predictor* expected is refused, however
clean its residual. **Every rung converges tight or is discarded** (L11); a
partially converged rung never seeds the next.

**Transcription:** trapezoid, plus switch-aware refinement (decided 2026-09-17).
The orbit work found the solver's Hessian approximation — not collocation
order — to be the wall for min-fuel, and that Hermite–Simpson made it worse
under the same solver (L6). So `fmincon` is the default for all three, and
**CasADi+IPOPT's exact Hessian is the named escalation for min-fuel** if it
smears or stalls. That is a decision point recorded in the min-fuel chain's
parameter block, not a surprise.

## Study scripts

`cartpole_minenergy_study.m` exists and is the template. `cartpole_mintime_study.m`
and `cartpole_minfuel_study.m` instantiate it with their objective's sections
6–7 from the earlier spec's table, with these revisions:

- **N7 (switch structure) REPORTS the switch count and times; it does not gate on
  the count.** The count is a mesh-sensitive, basin-sensitive integer (L3). N7
  gates on the switching function crossing its threshold *at* the realised event
  times, and on the sign law holding between them. A count that differs from the
  direct seed's is a finding to look at, printed beside a band from two or three
  mesh densities.
- **V3 (normality)** is added to every study: the extremal is normal — see the
  gate below. Without it the necessary conditions can hold vacuously.
- **X2 (perturbed seed) passes the duals** (L13). A re-solve that receives only
  the primal is a dual restart and proves nothing.
- Three gates in the template are annotated *carried for the bounded problems*
  (N6's finite-`d` probe, `tol.u`, V2b's time axis). Each study must re-test
  those claims once the control is bounded rather than inherit the comment
  (`mutation-sweep-finds-dead-checks`).

The study script keeps its explicit fail state: a failed N, X or V gate throws;
a failed or unresolved S gate is a finding, thrown only under `strict`.

## The library: what moves, and what pulls it

From the 2026-09-17 inventory of every optimality check in `orbit_transfer/`.
**A** = generic now; **B** = one named coupling to cut; **C** = campaign
physics, a cart-pole version is a rewrite not a promotion.

| candidate | list | coupling to cut | pulled by |
|---|---|---|---|
| `certified_guard`, `scalar_verdict` | A | — | min-time chain, stage 5 |
| `verify_common/pmp/pmp_residual` | A | — (already closure-based on `fFun`/`LFun`) | min-time chain, stage 8 |
| `ss_bvp_accept` | A | — (rides `ms_bvp`'s `prob` contract) | min-time chain, stage 7 |
| `conj_resolve` | A | — (tested on synthetic matrices) | min-time chain, stage 8 |
| `pmp_pointwise_checks` | **B** | hardcoded 14-state `[r;v;m;λ]` layout and the thrust/exhaust formula for recovering the applied control. Cut: externalise "recover the control from the field" as a caller closure, the move already made for `oc.fly_control` | min-time chain, stage 8 — **this is the instrument Mike named** (H constancy, transversality, adjoint residual, min-principle gap) |
| `conj_spectrum` | B | hardwired `mintime_prop_seg` and `muStar`/`Tmax`/`c`. Cut: a caller-supplied STM-providing `prop` closure, mirroring `ms_bvp` | min-fuel chain, stage 8 |
| `lift_margin`, `lift_space_dim` | A | — (pure linear algebra) | the normality gate |
| `foc_ipopt_inertia` | A | — (interprets a `regHistory` vector) | min-fuel chain, only if IPOPT is escalated to |
| `foc_check` | B | needs a CasADi `opti` + a labelled constraint-row registry and a `foc_manifest` entry | **refused for now**: our direct solves are `fmincon`, so the KKT/first-variation check is written against `fmincon`'s `lambda` in stage 5, with its sign convention verified by a minimal test first (L23) |
| `h6_margin`, `mintime_hypothesis_gates`, `validate_flight`, `certify_root`, `report_optimality` | C | the theory is min-time CR3BP (all-burn switching, reduced-Hamiltonian spurious zero, pumpkyn witness) | not promoted; `report_optimality`'s three-section scorecard *pattern* is adopted for the cart-pole certificate printer |

Rules for every promotion, all of which the `ms_conjugate_test` move followed:

1. The consumer exists and calls `oc.<fn>` before the README says it does.
2. Baseline gates run **before** the move: `golden_cells`, the four-campaign
   ladder, the function's own tests. Then the move. Then the same gates, compared
   line by line. A delegate that changes a number is not a delegate.
3. Executable lines of the promoted function are byte-identical unless the task
   is explicitly a *generalisation* (list B), in which case the orbit path
   through the generalised function must reproduce its numbers exactly and the
   cart-pole path is the new test.
4. Measure before extracting (L22). A candidate that turns out to carry a hidden
   coupling is refused and the refusal recorded with the coupling named.

## The normality gate

No check anywhere in the repo excludes an abnormal extremal — `(λ₀, λ) ≢ 0`
with `λ₀ = 0`, on which the Pontryagin conditions hold vacuously. The register
says so in italics (*"NOT CHECKED anywhere"*). The only instance is
`mintime_hypothesis_gates`' H1′, built for min-time CR3BP. For the cart-pole this
matters most in min-time, where an abnormal extremal is plausible rather than
academic, and it is exactly the class of gap the widened goal exists to close for
both plants.

**What it checks.** For a candidate trajectory, the multipliers `λ(·)` that
satisfy the adjoint equation and every condition an abnormal extremal would have
to meet form a linear space `S` (the conditions are linear in `λ` once `x(t)` and
`u(t)` are fixed). Normality holds when `S` contains no non-zero abnormal
multiplier — operationally, when the numeric rank of the constraint matrix built
from those conditions, evaluated through the state-transition matrices the chain
already has, leaves `dim S = 0` with a margin above the measured matrix error.
The per-objective conditions differ: for min-energy, `λ'g ≡ 0` along the arc;
for min-time, `σ(t_s) = 0` at each switch together with `H ≡ 0` on the free
horizon; for min-fuel, the corresponding conditions on burn and coast arcs.

**How it is built.** The rank-with-margin machinery is `lift_space_dim` and
`lift_margin`, promoted from `costate_common` as-is (list A). The **per-objective
constraint matrix is derived in the normality sub-project's plan, not here**, and
that derivation gets an adversarial review before it is coded — the orbit work
found three steps of a written second-order proof wrong under exactly such a
review while the numeric verdicts stayed correct (L21). The gate is **V3** in
every study script and runs in stage 8 of every chain.

## Lessons from the orbit-transfer work — binding constraints

Each of these was paid for. Each becomes a Global Constraint in every plan
written from this spec. The number in each is what it cost.

**Seeds and basins**
- **L1.** A machine-tight defect is a statement about the discretisation, not the
  trajectory: defect 1.4e-14 against a true re-integrated error of 1441 km, and
  refinement was *not* monotone (N 400→800 got worse). Every headline solution
  is re-integrated with an independent integrator, and refinement is checked
  for monotonicity.
- **L2.** Mesh density selects the basin: same solver, same guess, N = 400/800/1600
  gave 4.68/4.99/6.41 against a true 4.02 — 37% spread, none right. Every cold
  solve runs at two or three densities before it is called the extremal.
- **L5.** Min-energy warm-starts both bang problems — "by far the most effective
  single move in this campaign." Min-fuel shooting from a covector seed stalled at
  0.236 where min-energy reached 6e-12.
- **L15.** Direct duals are not shooting seeds until flown: PMP flight from catalog
  `λ₀` missed by 36,000–400,000 km while the recorded control landed within 10 km.
  The seed probe in stage 6 is not optional.
- **L-branch.** Continuation maps only the branch connected to its anchor; a
  branch-blind cold direct re-solve found a family 28% faster over a range the
  continuation structurally could not reach. Cold-solve periodically at the
  operating point and compare.

**Switch structure**
- **L3.** The switch count is a mesh-sensitive, basin-sensitive integer: 823 at
  8 nodes/rev against a converged band of 866 ± 5; three re-solves gave 24 vs the
  published 25; the 10 N basin flips 19 ↔ 24 switches on a `t_f` change of 2e-5.
  Report a band from two or three densities; judge on cost and `t_f`; never gate
  on the integer.
- **L4.** Every switch statistic is mesh-normalised before it is believed: a raw
  regular-switching value of 1.4e-7 read as "node-grazing" became 27 once divided
  by the local trapezoid weight, and the finding was retracted.
- **L7.** The saltation matrix is mandatory in any STM through a switch, and the
  STM is tested against finite differences of the state map: 1.35e-3 without it
  against 1.35e-7 for the continuous case. Omitting it produced a published
  "structurally unsuited" verdict that was retracted.
- **L8.** A grazing switch (`nᵀF⁻ → 0`) is a corner no continuation passes; the
  propagator refuses by name, and a distance-to-grazing statistic (min |dσ̇/dt|
  over active switches, or the largest sub-threshold extremum of σ) is printed
  because it predicts the wall one or two rungs ahead.
- **L6.** For min-fuel the solver's Hessian approximation, not the transcription,
  was the wall: L-BFGS interior point plateaued at defect ~1e-3, Hermite–Simpson
  made it worse, an exact Hessian fixed it. Suspect the solver × Hessian pair
  before blaming mesh order.
- **L-weak.** Bang-bang extremals are weak minima at the NLP level (270 flat
  directions in one reduced Hessian, curvature zero by construction). A native
  inertia check certifies a weak minimum; a strict claim needs the switching-time
  Hessian, which no campaign has fully built.

**Homotopy and continuation**
- **L9.** A failed rung proves nothing about feasibility: a "wall" at 75 mN was a
  stepping artefact pseudo-arclength walked straight through; a "fold" at 72 mN
  was a loss of normality (`|λ₀|` 46 → 3409) and the branch continued on a
  homogeneous chart. Distinguish fold from normality loss by costate boundedness;
  cold-solve at the stalled point before writing down a mechanism.
- **L10.** Band-gate every step: a rung returned defect 3.9e-14, feasible and
  clean, at `t_f` 16× its own predictor. Refuse anything outside `[1/3, 3]×` the
  step's prediction.
- **L11.** Under-iteration compounds silently: at `maxIter` 1500 a small residual
  carried forward until the tail collapsed into a restoration spiral; at 3000 the
  ladder walked clean. Every rung converges tight or is discarded.
- **L-sharpen.** Sharpening a smoothing parameter past the point the objective
  converged cost 8–12 failures per cell for a 1e-6 change: the endgame is a
  switching-time problem, not a homotopy. Stop sharpening at a stated rule.
- **L-race.** No smoothing family dominates: under one budget across 21 cases the
  ε-ladder failed 1, a Huber jump law failed 5 (but ran 5× faster where it
  worked), a continuous-ramp hybrid failed 3, and each had cases only it solved,
  unpredictable from any feature in advance. Race families; do not commit to one
  before the real problem has been tried.

**Solver and numerics**
- **L12.** Gates sit at the measured floor of the hardest case with margin, not at
  a round number: `tolR` 1e-10 retired an arm whose floor was 1.0e-10; 3e-10
  recovered it.
- **L13.** A warm-start verification that passes no multipliers is a dual restart:
  7 of 17 rows wandered up to 800 iterations before anyone noticed.
- **L14.** Say which evaluation produced a residual: 7.84e-12 in the solver's
  propagation, 1.41e-09 re-propagated plainly — 180×, both legitimate.
- **L23.** The multiplier sign convention is verified by a minimal test before
  anything differentiates against it: `opti.dual()` canonicalisation corrupted
  44–60% of costate signs behind a physically plausible eccentricity correlation
  and cost nine days and a published note. ex3's correlation vote with the
  amplitude ratio beside it is the defensive pattern; keep it.
- **L24.** Do not stack manual scaling on a solver's auto-scaling: a single hand-
  scaled row sent an easy convex solve into restoration failure.

**Verification discipline**
- **L17.** A gate that is computed is not a gate that is enforced — the most
  recurring failure in the project, named as "one habit, three instances" in a
  two-day span and recurring four more times after being written up. Every gate
  ships with a test that feeds it a known-bad input and asserts a named refusal.
- **L-dead.** Nine checks on this very demo could not fail; reading them found
  none, mutating them found most. Every new gate gets a deliberate-failure
  experiment before it counts (`mutation-sweep-finds-dead-checks`).
- **L19.** First-order agreement is not local optimality: 12 of 14 candidates that
  passed residual, flown-arrival and witness checks to 1e-10 were refuted by the
  conjugate test alone.
- **L16.** Metadata rides with the data and is asserted at every consumer: a
  catalog labelled Isp 1710 s beside physics at 900 s made every downstream ΔV
  1.9× too large. Every saved cart-pole artifact carries `p`, the bound, the
  horizon rule and the mesh, and every reader asserts them.
- **L18.** A chain's live path is executed under test before it is trusted: one
  chain had never run its non-cached path and carried three defects, including
  checkpoints matched positionally rather than by key. Checkpoints are keyed by
  content.
- **L20.** Language: "numerically certified under a stated policy", never
  "certified"; "fastest found", never "minimum". A sampled necessary condition is
  not a proof of anything between the samples.
- **L21.** Adversarial review before shipping — one pass cost $1.49 and found 41
  findings including three release blockers on code with nine green suites — and
  when a review names the decisive test, that test is built first: the one that
  was skipped stayed open nine days and took under an hour once built.
- **L22.** Measure before extracting; expect refusals: the 2026-07-26 pass refused
  5 of 8 candidates, and the comparison found a NaN a merge would have hidden.
- **L-mech.** "A satisfying mechanism that fits wrong measurements is more
  dangerous than an unexplained anomaly, because it stops the search." An anomaly
  that correlates with a physical quantity says where to look in the code first.

## Sub-projects, in order

Each gets its own plan. Order is forced by data flow and by where the skeleton
is cheapest to establish.

**1 — the min-energy chain.** `cartpole_minenergy_chain.m` over the pieces ex3
already has (`gen_direct_ref` becomes a front door with metadata; the harvest,
seed probe, shoot and verify stages are ex3's existing steps, staged and
banked). Stages 3–4 dark. This defines the whole 0–9 skeleton on the case that
already converges, and its stage-9 export is the warm-start artifact both bang
problems read (L5). Pulls nothing new. Smallest; a day.

**2 — the min-time chain, stages 0–6.** The direct half: `mintime_direct_solve`
(free `t_f` on a normalised mesh), the `u_max` ladder with band gates, switch-
aware refinement, the true-residual and mesh-density checks, harvest and seed
probe. Produces `t_min`, the switch band, the `t_min(u_max)` front, and the
banked seed. Pulls `certified_guard`, `scalar_verdict`. The `fmincon` first-
variation check is written here with its sign convention tested first (L23).

**3 — the min-time chain, stages 7–9, and `cartpole_mintime_study.m`.**
`bangbang_prop` into `cartpole_common` with its FD-of-the-state-map oracle (L7)
and grazing refusal (L8); the 5-unknown BVP; the free-time conjugate form (which
needs `info.Yend` and `spec.flow`, or the instrument silently returns
UNDETERMINED). Pulls `pmp_residual`, `ss_bvp_accept`, `conj_resolve`, and the
generalised `pmp_pointwise_checks` — the promotion with the most reach into both
trees, gated on the orbit path reproducing its numbers exactly.

**4 — the normality gate.** Derivation of the per-objective constraint matrices,
adversarially reviewed before coding; `lift_space_dim` and `lift_margin`
promoted; V3 added to all three studies and stage 8 of all three chains; a
deliberate abnormal fixture that the gate must refuse by name.

**5 — the min-fuel chain and `cartpole_minfuel_study.m`.** Needs `t_min`. Both
smoothing families raced (L-race); the sharpening stop rule (L-sharpen); the
three-branch control with `|σ| = 1` events; singular-arc detection reported, not
constructed; IPOPT escalation available. Pulls the generalised `conj_spectrum`
and, if IPOPT is used, `foc_ipopt_inertia`.

**6 — documentation and the rename.** ex3 → `ex3_cart_pole_minenergy` by `git
mv` gated on the suite; the teaching document revised to cover all three
objectives with the formula-vs-switch contrast now backed by code; `oclib`'s
README and CLAUDE.md brought to reality; `orbit_transfer` delegates and TODO
reconciled.

**S — stretch: the switching-time second-order test.** The strict-minimality
instrument for bang-bang extremals (Maurer–Osmolovskii) that no orbit campaign
has fully built (L-weak). The cart-pole, with a handful of switches, is the right
testbed: cheap, checkable by mutation, and a second consumer for the orbit work
the day it is generic. Queued, not promised.

## Verification — what "done" means

Every sub-project, before it is called done:
- `run_all_tests` green in the example tree **and** `golden_cells` 20/20 with
  unchanged numbers **and** the four-campaign ladder at its recorded figures, for
  any sub-project that touches `oclib` or `costate_common`.
- Every new gate has a refusal test (L17) and a recorded deliberate-failure
  experiment (L-dead).
- Every headline solution is re-integrated independently (L1), was solved at two
  or three mesh densities (L2), and quotes its switch count as a band (L3).
- The chain's live path was executed under test, not only its cached path (L18).
- Every saved artifact carries its metadata and every reader asserts it (L16).
- For a promotion: the consumer calls `oc.<fn>` before the README says so; gates
  ran before and after; executable lines are byte-identical or the orbit path
  reproduces its numbers exactly.
- An adversarial review pass on the verification code before the sub-project is
  merged (L21), with GPT-6 Astra reserved for the normality derivation and the
  min-fuel switch structure.

## Out of scope

- State constraints (a track limit on `q1`).
- Singular arcs as a *solved* case: detected and reported, never constructed.
- A min-fuel front over many `t_f`; one `t_f = 1.5·t_min` here.
- The switching-time Hessian (sub-project S) as a promised deliverable.
- Rewriting ex2; it remains the direct-collocation reference.
- Generalising list-C functions whose theory is CR3BP-specific.

## Risks

| risk | consequence | response |
|---|---|---|
| a promotion moves an orbit number | five campaigns silently changed | gates before and after, line by line; a delegate that changes a number is reverted, not tuned |
| a "generic" candidate carries a hidden coupling | the cart-pole path is wrong while the orbit path stays green | measure before extracting; refuse and record (L22) |
| the normality constraint matrix is derived wrong | a gate that passes everything | adversarial review of the derivation before code; a deliberate abnormal fixture it must refuse |
| grazing switch on the `u_max` ladder | saltation denominator → 0 | refuse by name; the distance-to-grazing statistic; cold-solve past the rung (L8, L9) |
| min-fuel smears under `fmincon` | no bang-off-bang structure resolved | IPOPT escalation recorded in the parameter block (L6) |
| chattering near 40 N | event loop stalls | cap the event count; report the band; the count is not a gate (L3) |
| mesh density selects a wrong basin | a confident wrong `t_min` | two or three densities per cold solve (L2) |
| a chain script that has only ever run cached | a defect ships on first live run | live-path test (L18) |
| the switch count from the seed does not match the shoot | mistaken for a failure | reported, not gated; the band and the cold re-solve decide (L3, L-branch) |
