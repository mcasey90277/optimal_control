# DESIGN — Campaign A: Thrust ladder 10 N → 0.1 N

**Goal.** Complete the blocked leg of the HMG-2004 reproduction: the thrust
ladder (paper Table 3) from 10 N down to 0.1 N, delivering (i) the Fig-23
overlay — m_f(c_tf) curves at multiple thrusts and the near-independence test,
(ii) empirical law R0 (T_max·t_f,min ≈ C) across two decades of thrust, and
(iii) the low-thrust structure counts (revs, switches) vs Table 3.

Date: 2026-07-17. **Status (2026-07-18, Task 11 close-out): Phase 2 (MEE +
ΔL) IMPLEMENTED and the ladder run 10 N → 0.5 N.** Phase 0/1 (Cartesian
5 N/2.5 N patch attempts) superseded — Phase 2 subsumed the problem entirely,
per the Phase-2 gate note below. Dependency: inherits the verification
front-end from Campaign B (DESIGN_dual_map.md) — Campaign B delivered a
correct MEE PMP verifier (Task 10) whose primal gates are unaffected by the
open raw-dual finding (see Results summary).

**Results summary (all reviewer-verified, `.superpowers/sdd/progress.md`
Tasks 4-11):**
- **Gate P2a (10 N cross-formulation) PASSED**: MEE m_f=1377.10 kg vs
  Cartesian 1376.74 kg (diff 0.36 kg < 0.5 kg gate) — the linchpin (Task 4).
- **Ladder 10→0.5 N complete** (fuel, c_tf=1.5): m_f = 1377.10 / 1364.54 /
  1369.79 / 1371.44 / 1375.28 kg at 10/5/2.5/1/0.5 N; switches 19/32/76/171/362;
  revs 7.33/14.16/27.84/69.15/138.60 vs paper 7.5/15/30/74.5/149.
- **R0 law**: 4 independently certified min-time anchors (10/5/2.5/1 N) give
  T·t_f,min spread **0.72%** (mean 850.0 N·h vs paper ≈850 N·h). The 0.5 N
  row is an **R0-law estimate**, not an independent anchor (see footnote 1
  below) — deliberately excluded from the fit.
- **PSR (switch-aware refinement) ported and validated**: 1 N m_f went
  171-switch uniform 1370.36 kg → PSR-refined **1371.44 kg** (supersedes);
  0.5 N reached m_f=1375.28 kg/sw=362 after 4 PSR rounds (budget-limited,
  stopReason=`maxRounds`).
- **PMP verifier (Task 10) delivered and proven correct**; gates FAIL on raw
  IPOPT duals at high eccentricity (characterized, not a verifier bug — see
  Campaign B escalate branch). Primal certifications (m_f/switches/revs
  above) are defect/terminal-gated and unaffected.
- Deliverable figures: `fig_table3.m` → `results/fig_table3.png` (Table-3
  analog + R0-law panel), `fig_front_mee.m` → `results/fig_front_mee.png`
  (single-c_tf Fig-23-adjacent overlay, honestly not a multi-c_tf overlay —
  only one c_tf per thrust was ever solved).

**Open items (not closed by Task 11, carried forward as future work):**
1. **0.5 N min-time anchor wall** — 7 configurations attempted, best defect
   0.0545, reproducible MEX crashes (same `libcoinmumps` signature); the
   0.5 N row's anchor stays an R0-law estimate until this is resolved.
2. **0.2 N and 0.1 N** — honestly not attempted; the 0.5 N wall is where the
   deep-descent effort stopped (footnote 6, README.md).
3. **Dual/PMP escalate branch** (Campaign B, `DESIGN_dual_map.md`) — probe
   raw `lam_g` via `nlpsol` bypassing `opti.dual` to root-cause the raw-dual
   KKT-stationarity failure at high eccentricity; primal work does not block
   on this.
4. **Full Table 3 (down to 0.1 N)** — not reached; this campaign delivers
   10→0.5 N (two decades minus the last quarter-decade), not the paper's
   full 10→0.1 N span.
5. Minor/inherited: decadeMin stall-guard not yet rung-size-aware; the
   171-vs-179 switch-count question at 1 N remains open per Task 8's
   reviewer insight (windowed PSR cannot discover switch pairs outside its
   own neighbor windows — a hybrid periodic uniform-sweep round is the
   suggested fix, not yet built).
6. **T7c cross-process hang, uninvestigated** (Task 7c, `task-7-report.md`
   item 3): round 14 of the 1 N manual continuation hung for an unexplained
   42+ minutes on a nominally bit-identical computation. If CasADi/MUMPS's
   underlying BLAS is multi-threaded, the exact Newton path through an
   ill-conditioned region may not be bit-reproducible run-to-run across
   fresh process launches (unlike the within-session determinism confirmed
   at rounds 6-13) — an important caveat for anyone relying on
   retry-determinism reasoning for this solver stack. Not root-caused; not
   blocking primal results above, but worth investigating before the next
   deep-descent attempt (Task 9's 0.2/0.1 N rungs).

---

## 0. Front door + endpoint parameterization (post-Task-11 SDD, `.superpowers/sdd/task-1..8`)

On top of the ladder work above, a separate small SDD pass generalized the
fixed GEO-from-paper-GTO pipeline into a parameterized one and wrapped it in
a single front door, `run_gergaud.m` (README.md, "Front door: `run_gergaud`"
section, has the user-facing usage; this note is the design-side pointer).

**Core change.** `casadi_lt_mee.m` gained `opts.xf` — the 5-element MEE
terminal target `[P;ex;ey;hx;hy]`, defaulting to `[1;0;0;0;0]` (GEO in the
solver's own normalized units) — replacing five previously-hardcoded GEO
equality constraints with a loop over `xf`. `mee_seed.m` gained
`opts.initElems` — the 7-element seed initial state
`[P;ex;ey;hx;hy;m;t]`, which when absent or empty falls through to the
exact pre-existing literal (no `tan(3.5°)` recompute at seed time) rather
than a numerically-close-but-not-identical reconstruction. Both were then
threaded through the downstream drivers — `homotopy_mee.m`,
`run_transfer_mee.m`, `run_mintime_mee.m` — including their cache
fingerprints, so a custom-endpoint run cannot silently reuse or corrupt a
default-endpoint cache.

**Default-preservation guarantee.** At the paper's own endpoints
(P0=11625 km, e0=0.75, i0=7°, GTO-like start; Pf=42165 km, ef=0, if=0°,
GEO), `initElems` resolves to `[]` and `xf` resolves to `[1;0;0;0;0]` —
byte-identical to the pre-parameterization code path. Every certified
number in this document (the 10/5/2.5/1/0.5 N ladder, the R0 law, the PSR
results) was produced before this change and reproduced after it at the
same values (Task 1/2/4 regression checks in `.superpowers/sdd/progress.md`,
e.g. MEE 10 N m_f=1377.10 kg / tfmin=22.2206 ND reproduced with zero
fingerprint mismatch) — the parameterization is additive, not a rewrite of
the ladder's numerics.

**Front door.** `run_gergaud.m` is the user-facing consumer of both new
option fields: it resolves `(P0_km,e0,i0_deg)`/`(Pf_km,ef,if_deg)` into
`initElems`/`xf` (falling back to the defaults above whenever the inputs
match the paper/GEO endpoints exactly), selects cache-hit vs. live-solve
(`'auto'`/`'solve'`/`'probe'` modes), and assembles/prints one Table-3 row
(`gergaud_row.m`/`gergaud_row_str.m`) plus an optional plot/movie via the
`mee_res_to_cart_res.m` → `transfer_movie.m`/`gergaud_plot.m` chain. It adds
no new solver physics — it is a thin front end over the ladder machinery
documented in the rest of this file. A significantly non-GEO-like custom
final orbit is explicitly out of the validated scope (research-probe
territory), since the solver/seed were only ever validated against the
paper's GTO→GEO case and the certified ladder rungs.

---

## 1. Where we actually stand (corrected diagnosis)

The 5 N anchor blocked after six strategies (Task 14 report). Post-mortem
CORRECTION to the ledgered hypothesis: the failed anchor ran at N=1200 over ~9
min-time revs ≈ 133 nodes/rev — the SAME density as the successful 10 N anchors,
so nodes-per-rev parity does NOT explain the anchor failure (it does bind for
FUEL stages: 5 N fuel ≈ 15 revs wants N≈2000+, and Table 3's counts grow ∝ 1/T).
Better-supported suspects, both untested:

- **S1 — warmTight on loose iterates.** Continuation rounds 2+ ran warmTight
  (mu=1e-4, monotone) from defect ~5e-3 iterates. Tight warm starts are for
  re-solving AT a converged point; on an infeasible iterate they push IPOPT into
  restoration — and "restoration + false Infeasible_Problem_Detected at a
  defect floor" is exactly the observed signature (and the tulip campaign's
  documented pre-no-resample signature). The 10 N anchors survived only because
  round 1 landed close enough.
- **S2 — iteration starvation.** 3000 iters/round on a 2×-bigger NLP; the
  Route-B energy stage was given NO continuation rounds at all (single shot,
  errored at defect 3.0e-2).

Also standing (from Task 14): stretch-seeding across thrust is topology-flawed
(stretching adds no revolutions) — cold tangential seeds wind the correct revs
by construction; per-thrust pipelines stay independent.

**Honest-claim update (2026-07-17, per the three-way core review).** The
paper's own 2004 verdict that a direct/collocation attack on the low-thrust
ladder is "predictably unsuitable" below ~1 N assumed dense-Jacobian SQP with
finite-difference derivatives — a different machine than what we are running.
The two reviewers pushed on this from different, both-correct axes (host
adjudication, triage §Q4-disagreement-2): Gemini argues the **dimensional**
objection is defeated outright — a ~1e5-variable **sparse** NLP with exact-AD
Hessians into IPOPT/MUMPS (Phase 2's MEE+ΔL build) is categorically not the
regime the 2004 paper was warning about, and our own campaign routinely solves
problems of that size. Terra agrees on the dimensional point but holds that
the paper's REAL difficulty was never purely dimensional — it is
**combinatorial**: locating ~1500 bang-bang switches and holding the KKT
basin over 754 revolutions is untouched by a change of coordinates. Both are
right about their own axis. The honest claim for the paper writeup: *"MEE+L
defeats the dimensional objection to a direct attack; the residual difficulty
is combinatorial (switch/basin structure), not dimensional — and Phase 2's
thrust-continuation backbone (ΔL free) plus PSR-ported switch-aware
refinement is our answer to the combinatorial half, not to a dimensional one."*
Do not claim MEE+L alone "solves" the paper's objection; claim it correctly,
on the axis it actually addresses.

## 2. Phased plan with decision gates

### Phase 0 — cheap diagnosis fixes at 5 N (Cartesian stack, ~1 day)

**Amended 2026-07-17 per the three-way core review** (`doc/reviews/2026-07-17_triage.md`).
Both reviewers endorse testing the μ-staged fix first (terra: "cheapest and
most discriminating"; Gemini ranks it #1: "exactly right"), but terra reframes
it and warns it may not be sufficient alone.

1. **μ-staged continuation, selected by FEASIBILITY not call number (terra's
   framing, supersedes the original one-line "warmTight=false while defect
   ≥1e-6" rule).** In `run_mintime`'s continuation rounds: pick the
   barrier/warm-start policy from the current iterate's feasibility, not from
   which round number it is. Concretely — warmTight=false (adaptive μ) while
   maxDefect ≥ 1e-6, tight only below *that measured value*, regardless of
   round index; **never reuse multipliers from a restoration-phase or failed
   solve**; after each round, **retain the new primal iterate only if
   feasibility improved** over the incoming one (otherwise keep the prior
   iterate and retry with adaptive μ). Tests S1. Terra's warning: this cures
   the policy bug for certain but is **not guaranteed to cure a basin
   problem** — terra ranks a doubled-rev-count winding-basin issue #1 overall
   (Gemini ranks warmTight #1). Host adjudication: run warmTight first (both
   agree it's cheapest) but do not declare S1 closed on a pass alone; if the
   5 N anchor still stalls after (1)+(2), treat it as basin evidence, not a
   policy-bug residual. Note also: the cold tangential seed already winds the
   correct rev count by construction, which weakens the pure-basin story for
   the *cold-seed* leg specifically — but not for the continuation *rounds*,
   where each re-solve can still drift basins.
2. **maxIter 6000/round**; guard constants unchanged. (Tests S2.)
3. **Route-B energy stage gets the same patient rounds** (it is a fixedtf ε=1
   solve; rounds are legal there too).
4. **Prograde guard (NEW — required, not optional; see also Phase 1 below).**
   `casadi_lt_2body.m:103-109`'s free-longitude insertion manifold admits the
   retrograde GEO circle (h_z<0) as a legitimate solution branch — flagged
   independently by BOTH reviewers (terra #4, Gemini #2), which the host
   triage treats as high-confidence convergent evidence. `DESIGN.md` §2
   explicitly anticipated this ("the set also admits the retrograde orbit...
   add an h_z>0 guard only if the solver ever drifts") and deferred the fix;
   **that deferral is overturned** by the review. Add `opti.subject_to(h_z(end)
   >= h_min > 0)` (angular momentum z-component, h_min a small positive
   constant) to the manifold terminal, or a branch homotopy if the hard
   inequality fights convergence. Apply before or alongside (1)-(3) — a 5 N
   run that "succeeds" onto the wrong branch is not evidence for S1/S2 either
   way.
5. Order of attempts at 5 N: cold-seed N=1200 two-stage with (1)+(2)+(4); if
   stalled, Route-B with (3)+(4).

**Gate P0:** 5 N anchor converged (success, defect <1e-8, tfmin ≈ 44.4 ND
±30%, prograde — now *enforced*, not merely checked post hoc). PASS ⇒ Phase 1.
FAIL after both attempts ⇒ skip to Phase 2 (the MEE build subsumes the
problem; do not grind Cartesian further).

### Phase 1 — Cartesian ceiling: 5 N and 2.5 N fuel points (+1 N stretch)

The prograde guard (Phase 0 item 4) lives on the shared manifold terminal in
`casadi_lt_2body.m`, so it is inherited automatically by every Phase 1 run —
no separate action needed here, just confirm it stayed on for the 2.5 N and
1 N points too (they reuse the same terminal builder).

- 5 N fuel: c_tf=1.5, manifold, fresh seed, N=2400 (~160/rev at 15 revs).
  Expect Table 3: ~15 revs, ~36 switches.
- 2.5 N: anchor N≈2400 (18 min-time revs), fuel N≈4000 (30 revs; ~36k vars —
  hours/solve; budget accordingly). Expect ~30 revs, ~73 switches.
- 1 N (STRETCH, only if 2.5 N is smooth): fuel ~75 revs ⇒ N≈10k, ~90k vars —
  the realistic Cartesian limit; one attempt, honest report either way.
- Per-point resume + certified-only caching throughout (existing machinery);
  the MEX init-crash protocol stands (relaunch once; log the pattern).

**Deliverables at Phase 1 exit:** Fig-23 overlay with 2–3 curves; law R0 spread
across {10, 5, 2.5} N (gate: <10%, paper C ≈ 850 N·h); m_f near-independence
check (gate: within ~5 kg at c_tf=1.5). This alone closes the paper's own
Fig-23 scope (their overlay is 5/2.5/1 N).

**Gate P1:** at least the 5 N fuel point certified. The ladder BELOW 1 N is
declared out of Cartesian scope regardless of P1 (150–750 revs ⇒ 10⁵–10⁶
variables — a formulation problem, not a tuning problem).

### Phase 2 — the paper's own move: equinoctial (MEE) formulation

**Substantially rewritten 2026-07-17** to incorporate Gemini's structural
insight from the three-way core review (`doc/reviews/2026-07-17_triage.md`)
— a point neither the host nor terra had. The original framing below (fixed
L_f span, ~10–20 nodes/rev, node budget 10–20) is superseded by the three
structural revisions and the node-budget correction that follow.

Coordinates (P, eₓ, e_y, hₓ, h_y, m) slow + true longitude L fast, **L as the
independent variable** (the role Sundman τ plays now; the paper's conclusion
flags longitude-as-variable as their own next step). Both reviewers endorse
this: **L̇ > 0 strictly** for this regime (thrust ~1e-4 g cannot stall h/r²),
and **L subsumes the Sundman clock outright — no separate τ or cScale in the
MEE core**, which also sidesteps Campaign B's whole dual/primer problem class
by construction. Gemini adds: dividing by L̇ automatically concentrates mesh
density at apogee, exactly where the burns are.

**Three structural revisions (Gemini, all three recorded — this is the core
of the rewrite):**

1. **ΔL (total longitude span) becomes a scalar DECISION VARIABLE**, not a
   fixed target computed up front from the empirical laws. Mesh on a FIXED
   unit grid σ∈[0,1]; longitude is parametrized `L(σ) = L₀ + σ·ΔL`; physical
   time is carried as a state with clock relation `dt/dσ = ΔL/L̇` (ΔL enters
   the clock exactly where τ_f used to sit in the Sundman formulation, except
   now it is being solved for, not fixed).
2. **The homotopy objective is rescaled by the spatial measure**: since the
   integration variable is now σ (fixed grid) and not L or t, the ε-continuation
   objective must integrate against the σ-measure with the ΔL/L̇ Jacobian
   folded in: `J(ε) = ∫₀¹ Φ_ε(s) · (ΔL/L̇) dσ`. Getting this Jacobian right is
   load-bearing — an unscaled σ-quadrature would silently misweight the
   objective toward wherever L̇ happens to be large.
3. **Consequence — the revolution count becomes a decision variable.** This
   is precisely what Cartesian+Sundman could NOT do: there, the rev count is
   frozen into the seed's geometry (fixed τ_f, tangential-thrust propagation
   winds a specific number of revs by construction), which is exactly WHY the
   rejected stretch-seeding-across-thrust approach was topology-flawed —
   stretching τ_f adds duration, not revolutions, so a seed built for one
   thrust's rev count can't be warm-started into another's. In MEE with ΔL
   free, thrust continuation works because **the optimizer grows ΔL itself**
   as T_max is stepped down — the rev count is no longer a seed-time
   commitment, it is an outcome of the solve.

**Thrust continuation returns as the Phase-2 BACKBONE (reverses the
Cartesian-era decision to abandon it).** Geometric ladder
10 → 5 → 2.5 → 1 → 0.5 → 0.2 → 0.1 N, each point warm-starting the next.
Gemini's argument for why this is now mandatory rather than merely helpful: a
COLD 0.1 N solve will trap against neighbouring 699/701-rev local optima (the
754-rev basin is one of many near-degenerate rev-count basins at that thrust
level); only a geometric Tmax-continuation path, walking ΔL up gradually as
thrust steps down, reaches the correct 754-rev basin reliably. This is also
exactly the paper's own "discrete continuation" device (their thrust ladder,
Table 3) — we are not inventing a new idea, we are recovering the paper's
device in a coordinate system where it actually works. **Explicitly note the
reversal:** Cartesian-era thrust continuation via time-stretching was
abandoned for a real reason (topology-flawed, per Task 14 and the Phase-1
Cartesian ceiling doc above) — that reason was the *time-stretch seeding
mechanism*, not the *continuation idea* itself. MEE+ΔL-free removes the flaw
by making the rev count optimizable, so the continuation idea is rehabilitated
without reviving the flawed mechanism.

**Node budget corrected.** The reviewers disagreed and the host split the
difference: terra argues **25–40 nodes/rev** ("15 is aggressive, not a safe
production number"), citing ~22.6k intervals at 754 revs as still tractable;
Gemini argues a **12–16/rev absolute floor** (~9–12k points at 754 revs). Host
adjudication (triage): plan for **25/rev nominal with a 15/rev floor probe**,
settled empirically by terra's prescribed convergence study at **20/30/40
nodes/rev** before committing the ladder's bottom end. At 754 revs and 25/rev
that is ~19k nodes (~O(1.5–2)×10⁵ variables with the 8-state-plus-control MEE
vector) — both reviewers agree this order of magnitude is tractable; this
supersedes the original doc's "~10–20 nodes/rev" placeholder.

**Switch smearing — terra's prescription adopted, Gemini's rejected (with
rationale).** Both reviewers agree Hermite-Simpson/pseudospectral methods
don't remove the Gibbs-ringing/smearing at bang-bang switches on a fixed mesh.
They prescribe different fixes:
- **terra: hp/adaptive mesh refinement steered by the switching-function
  zeros**, or switch times as explicit decision variables in a multi-phase
  formulation. **Adopted.** We already own this exact mechanism —
  **`NLP_lowThrust_GTO_tulip/PSR/`** (PMP-Steered Refinement) is switching-
  function-steered mesh refinement, already built and validated on the tulip
  problem. Campaign A Phase 2 should **port PSR to the MEE solver** rather
  than build a new refinement scheme from scratch.
- **Gemini: stop the ε-continuation at ε≈1e-4** instead of driving to 0,
  absorbing the smearing into physical fuel bias rather than KKT violation.
  **Rejected as the default.** Stopping short of ε=0 biases m_f, and our own
  tulip campaign explicitly *fixed* a legacy schedule that stopped at ε≈1e-3
  for exactly this reason (see `NLP_lowThrust_GTO_tulip/LOW_THRUST_MINFUEL_CAMPAIGN.md`).
  Recording the disagreement for the record: Gemini's option remains a
  fallback if PSR-on-MEE proves harder to port than expected, but it is not
  the plan of record.

**[ROBUSTNESS]** (Gemini) — at ΔL ≈ 4700 rad (the 754-rev end of the ladder),
raw `cos(L)`/`sin(L)` evaluated on large operands risks precision loss from
argument-reduction cancellation, and can damage AD derivative accuracy through
the trig chain. Wrap trig evaluations to `mod(L, 2π)` inside the dynamics (or
otherwise confirm CasADi's AD path produces analytically exact derivatives
through the wrap) before running the bottom-of-ladder points.

Build (reusing the entire campaign architecture — only the dynamics/clock
change):

| module | role |
|---|---|
| `casadi_lt_mee.m` | NEW solver core: Gauss variational equations (paper p.6, already transcribed in DESIGN.md §2's source), independent variable σ∈[0,1] with L(σ)=L₀+σ·ΔL and ΔL a decision variable, dt/dσ = ΔL/L̇ carried like the Sundman clock (t as a state), objective rescaled `∫Φ_ε·(ΔL/L̇)dσ`, `mod(L,2π)` trig guard, L subsumes cScale (no separate slack state), modes fixedtf/mintime, ε-homotopy objective, cone-eliminated [α;s] in the RTN frame (q,s,w — the paper's control frame), `L̇ ≥ L̇_min > 0` enforced/logged |
| `mee_seed.m` | seed by element-space propagation of the tangential-thrust law (elements drift smoothly — seeds are nearly defect-free by construction); feeds the geometric Tmax-continuation ladder, not a per-thrust cold start |
| `elements_to_cart.m` / `cart_to_elements.m` | already built + roundtrip-tested to 1e-10 (Task 3) — used for BCs and cross-checks |
| `run_mintime` / `run_transfer` / `run_ctf_sweep` | gain a `formulation` switch ('cartesian'|'mee') or thin MEE siblings — decide at plan time; prefer siblings to avoid destabilizing certified drivers |
| `verify_pmp_2body` | MEE variant of the switching-function check (Campaign B's front-end pattern); per triage action 9, do NOT gate Phase-2 validation on the raw-dual primer metric until Campaign B (`DESIGN_dual_map.md`) closes — m_f/structure gates stand on their own here |
| `NLP_lowThrust_GTO_tulip/PSR/` (ported, not new) | switching-function-steered mesh refinement, adapted from the tulip solver to the MEE solver — the adopted switch-smearing fix |

**Validation before use (non-negotiable).** Solve the 10 N / c_tf=1.5 case in
MEE and require m_f within ~0.5 kg of the Cartesian 1376.74/1377.05 kg and the
same burn structure — a cross-FORMULATION check stronger than any mesh check.
Then (if Phase 1 delivered them) cross-check 5 N and 2.5 N too. This gate is
independent of Campaign B and does not wait on it (see verifier row above).

**The ladder:** {1, 0.5, 0.2, 0.1} N via geometric Tmax-continuation
(10 → 5 → 2.5 → 1 → 0.5 → 0.2 → 0.1 N, each warm-starting the next, per the
backbone above — NOT independent cold solves per point). At each rung: hold
c_tf=1.5 fuel solve → structure counts, with ΔL emerging from the solve rather
than being prescribed. N ∝ revs at 25/rev nominal; iteration budgets ∝ N;
overnight jobs at the bottom end. Per-point resume, certified-only caching,
refinement spot-check (N→1.5N, plus the 20/30/40-nodes/rev convergence study
above) at 1 N and 0.1 N.

**Gate P2a:** MEE-vs-Cartesian agreement at 10 N. **Gate P2b:** each ladder
point certified + refinement-stable. **Exit:** full Table-3 analog (revs,
switches, t_fmin per thrust), complete Fig-23 overlay, law R0 across 10→0.1 N.

## 3. Milestones

- **L0** = Gate P0 (5 N anchor, Phase-0 fixes) — also feeds back into
  `run_mintime` for everyone.
- **L1** = Phase 1 exit (5/2.5 N fuel points, first overlay + law spread).
- **L2** = Gate P2a (MEE core validated at 10 N).
- **L3** = ladder complete (1/0.5/0.2/0.1 N) + final figures/tables.

## 4. Honesty & stop rules

- Table 3's switch counts are reported as BANDS with refinement checks (house
  rule); expect count growth ∝ 1/T (paper: 18→1786).
- Every phase has an explicit gate; a failed gate produces a documented block
  (attempt trajectories, not vibes) and a controller decision — no unbounded
  grinding. Phase 0 is time-boxed to ~1 day, Phase 1 to ~3 days of compute,
  Phase 2 build to ~2 days + ladder compute as scheduled background work.
- Isp remains 2000 s (family cross-check caveat carries into all new numbers).
- Per the 2026-07-17 triage (action 9): do NOT gate Phase-2 MEE validation on
  Campaign B's raw-dual primer metric until Campaign B closes — the m_f-match
  and burn-structure gates (Gate P2a/P2b) stand on their own.

## 5. Estimate

Phase 0: ~1 day. Phase 1: ~2–3 days (mostly background solves). Phase 2:
~1–2 days build/validation + ladder solves over several nights. Total to L3:
~1.5–2 weeks calendar with the machine working nights.

## 6. References

- Task 14 report (six-strategy 5 N record; topology finding; guard history).
- Paper: Table 3, Figs 19/20/23, conclusion (longitude-as-variable remark);
  p.6 Gauss equations (MEE dynamics source).
- DESIGN.md §5 (original out-of-scope declaration this campaign lifts).
- Campaign B: DESIGN_dual_map.md (verification front-end dependency).
- **`doc/reviews/2026-07-17_triage.md`** (source of this file's Phase 0
  prograde-guard/feasibility-policy amendments and the Phase 2 rewrite) + the
  two full reviews, `2026-07-17_core_review_gpt56terra.md` and
  `2026-07-17_core_review_gemini.md`.
- `NLP_lowThrust_GTO_tulip/PSR/` (PMP-Steered Refinement — ported into Phase 2
  as the switch-smearing fix) and `NLP_lowThrust_GTO_tulip/LOW_THRUST_MINFUEL_CAMPAIGN.md`
  (the legacy ε-schedule fuel-bias precedent cited against Gemini's ε≈1e-4
  stopping-point alternative).
