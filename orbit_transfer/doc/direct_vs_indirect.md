# Direct and indirect methods across min-time, min-energy, min-fuel

**Goal:** understand what each method is good and bad at for each of the three
objectives, and how to use them together. Written 2026-07-31, grounded in what
this repo has actually measured rather than in general claims.

---

## 1. Why the three problems are different: it is all about the control structure

The objective determines the structure of the optimal control, and that
structure is what makes one method easy and the other hard.

| objective | optimal control | structure to discover |
|---|---|---|
| **min-time** | thrust at its **upper bound always**; direction = primer vector | **none** — there is nothing to switch |
| **min-energy** | throttle a **smooth** function of the switching function | none — control is continuous |
| **min-fuel** | **bang-bang**, `s ∈ {0,1}`, switched by the sign of `S` | **the number and location of every switch** |

That single column explains most of what follows. Min-time and min-energy have
no combinatorial content; min-fuel has a discrete structure that must be found
before anything can be refined.

---

## 2. What we have measured

### Min-time

- **Indirect is the natural fit and it works.** `pumpkyn.cr3bp.tfMin` solves
  DRO→tulip in **8 unknowns** (7 costates + `t_f`) and converges in ~3 s from a
  good seed. There is no switching structure to guess, and free final time
  supplies the clean transversality condition `H(t_f) = 0`.
- **But only at low revolution count.** DRO→tulip is ~2.7 ND, a few
  revolutions. The GTO→tulip spiral is ~**40 revolutions**, and there the
  homogeneous costate map amplifies by ~**1e11** even with Sundman scaling —
  single shooting is hopeless regardless of the guess.
- **Direct handles the many-revolution case.** Min-time anchors were computed by
  direct collocation for tulip, ELFO and earth. With thrust pinned at maximum
  there is no bang-bang structure, so the direct method's main weakness does not
  apply.

**Rule of thumb: revolution count, not the objective, decides.** Few revs →
indirect. Tens of revs → direct.

### Min-energy

- **Both work, and this is the easy problem.** The control is smooth, so
  neither method faces a discrete structure or a discontinuous right-hand side.
- **Its real role is as a bridge.** Both campaigns reach min-fuel by continuing
  from min-energy (Bertrand–Épénoy, `J(ε) = ∫[s − ε·s(1−s)]dt`, ε: 1 → 0). The
  energy solution is the seed that makes the fuel solution reachable at all.
- **The band is limited.** The energy backbone exists only for
  t_f ≈ 1.12×–1.95× of minimum time; below that, even the *smooth* problem will
  not converge. That is a conditioning wall, not a bang-bang artifact.

### Min-fuel

This is where the methods genuinely diverge.

**Direct: finds the structure, blurs the details.**
- Discovers switch count and ordering with no prior assumption — the certified
  tulip flagship has 25 switches, none of them specified in advance.
- Scales: 13,920 nodes solved on the 1 N earth row.
- **Weakness 1 — switch times are grid-quantized.** Reading a switch as "the
  first node where `s > 0.5`" ties its location to the mesh. Measured across
  three earth rows: switch times never converged, and matched-switch counts
  *fall* under refinement.
- **Weakness 2 — basin multiplicity, and it is severe.** At one t_f the tulip
  admits **at least six local optima spanning 24, 25 and 26 switches**. Every
  certified earth row is beaten on its own mesh by 1.0–6.0 kg. The certified
  tulip flagship was beaten by 15.3 g. A single converged direct solution is an
  **upper bound on propellant**, not the answer.
- **Weakness 3 — no error bar.** Two attempts to measure discretization error
  were retracted under review. What survives is a lower bound.

**Indirect: exact where it converges, and it often does not.**
- Switch times are **exact**, not grid-quantized: `tfMinProp` uses event
  detection to stop and restart integration at each switch.
- Satisfies the necessary conditions by construction rather than approximately.
- **Weakness — conditioning, quantified.** The IFS attempt at 1.12× stalls with
  the Levenberg–Marquardt residual crawling 1.96 → **0.023 over 400
  iterations**, `cond(J) ≈ 5.9e9`, six orders short of the 1e-8 target.

---

## 3. The coordination you propose has been tried here, and the diagnosis is not what one would guess

**The proposal — use the direct solve to generate a costate guess for the
indirect solve — is already built.** `foc_check` maps NLP duals to costates via
the Hager covector mapping; `verify_pmp_mee` reconstructs the primer vector and
switching function; PSR supplies KKT-dual costates. The IFS ("Indirect
Finishing Solve") campaign is exactly this handoff.

**It stalls, and the reason is not guess quality.** The scaled-SVD
investigation found the near-null direction of the Jacobian sits:

- **76% on switch #4's node** — the shallowest zero-crossing of the switching
  function, `|dS/dτ| = 0.11` — at the seed, then
- **83% on the initial position costate `λ_r0`** at mid-crawl.

Two candidate explanations were tested and *refuted*: the near-double-root
terminal pair, and `τ_f` over-constraint. What remains is the textbook weakness
of indirect shooting — `λ_r` couples to the trajectory only indirectly and is
weakly determined over a short arc.

**So the correction to the hypothesis is:** a good costate guess is
**necessary but not sufficient**. The obstruction is the conditioning of the
*shooting operator itself*, which a better starting point does not change. Two
things do: shortening the arcs (multiple shooting), and removing the shallow
switch from the residual.

---

## 4. The connection that only became visible today

The IFS null direction sits on the switch with `|dS/dτ| = 0.11`.

Independently, the front sweep found that the four t_f values where basin
multiplicity appears have `sdotMinRel` of **0.11 to 1.57**, against **27–28**
for the well-behaved 1.150× solutions — and the worst, **0.1147** at 1.120×, is
the row with the largest basin gain (+53 g).

**Near-tangential switch crossings are simultaneously where the indirect
Jacobian goes singular and where the direct method's basins proliferate.** One
degeneracy, two symptoms. When `Ṡ → 0` at a crossing the switch location stops
being well determined by the problem, so the indirect residual loses rank there
and the direct solver finds many nearly-equivalent structures.

This suggests `sdotMinRel` is a **predictor**: it should be computable from a
direct solution *before* attempting an indirect finish, and used to choose which
t_f to attempt.

---

## 5. A coordinated program

**Each method should do what it is good at, and the handoff should carry
structure, not just numbers.**

```
   min-energy solve  ──►  min-fuel DIRECT  ──►  BASIN SWEEP  ──►  INDIRECT finish
   (smooth, robust)       (finds structure)     (best of many)     (exact switches)
                                │                      │
                                └── costates ──────────┘
                                    switch count + ordering
                                    switch times (sub-grid)
                                    sdotMinRel  ──► go / no-go
```

1. **Min-energy seeds min-fuel.** Already how both campaigns work.
2. **Direct finds the structure.** Its combinatorial strength; nothing else
   discovers 25 switches unaided.
3. **Basin sweep before any refinement.** Refining the wrong local minimum to
   machine precision is wasted work, and we now know the first solution found is
   routinely not the best.
4. **Hand off four things, not one:** costates (Hager-mapped), the switch
   *count and ordering*, sub-grid switch *times*, and `sdotMinRel` as a
   go/no-go.
5. **Finish indirect — but with multiple shooting, not single**, and ideally
   with switch times as explicit variables. Single shooting over 40 revolutions
   is 1e11-amplified regardless of seed. The `ms_band` campaign reaching 1.12×
   with ten certified switches is evidence this is the right direction.

### Why this converges with the multi-phase build

A **fixed-structure multi-phase formulation with switch times as explicit
decision variables** is where the two methods meet:

- it is a *direct* NLP, so it inherits the direct method's conditioning;
- switch times are *variables*, so they are exact, as in the indirect method;
- the structure is *fixed*, which removes the shallow-switch rank deficiency
  from the residual — the thing that stalled IFS;
- it also fixes branch identity, which is what invalidated two mesh-study
  measurements.

That is the same construction three separate lines of work have now pointed at:
the mesh reviewers, the indirect stall, and this analysis.

---

## 6. What is not yet known

- **Whether `sdotMinRel` really predicts indirect convergence.** One
  correlation, two data points. Testable: compute it for the t_f values IFS
  attempted and see whether the failures line up.
- **Whether the direct-derived costates are accurate enough** once conditioning
  is handled. Never isolated, because conditioning always dominated first.
- **Where the revolution-count crossover sits.** Indirect works at ~3 revs and
  fails at ~40. The boundary is unmeasured, and it decides which method leads
  on a new problem.
- **Whether `minDeltaV`'s path-independent lower bound** can arbitrate between
  basins — untried, and the only method-independent check available.
