# DESIGN — Campaign A: Thrust ladder 10 N → 0.1 N

**Goal.** Complete the blocked leg of the HMG-2004 reproduction: the thrust
ladder (paper Table 3) from 10 N down to 0.1 N, delivering (i) the Fig-23
overlay — m_f(c_tf) curves at multiple thrusts and the near-independence test,
(ii) empirical law R0 (T_max·t_f,min ≈ C) across two decades of thrust, and
(iii) the low-thrust structure counts (revs, switches) vs Table 3.

Date: 2026-07-17. Status: approved design, pre-implementation.
Dependency: inherits the verification front-end from Campaign B (DESIGN_dual_map.md)
— run B first or in parallel; A's milestone gates are primal and do not block on B.

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

## 2. Phased plan with decision gates

### Phase 0 — cheap diagnosis fixes at 5 N (Cartesian stack, ~1 day)

1. **μ-staged continuation** in `run_mintime`'s rounds: warmTight=false
   (adaptive μ) while defect ≥ 1e-6; warmTight=true only below. (Tests S1;
   one-line change.)
2. **maxIter 6000/round**; guard constants unchanged. (Tests S2.)
3. **Route-B energy stage gets the same patient rounds** (it is a fixedtf ε=1
   solve; rounds are legal there too).
4. Order of attempts at 5 N: cold-seed N=1200 two-stage with (1)+(2); if
   stalled, Route-B with (3).

**Gate P0:** 5 N anchor converged (success, defect <1e-8, tfmin ≈ 44.4 ND
±30%, prograde). PASS ⇒ Phase 1. FAIL after both attempts ⇒ skip to Phase 2
(the MEE build subsumes the problem; do not grind Cartesian further).

### Phase 1 — Cartesian ceiling: 5 N and 2.5 N fuel points (+1 N stretch)

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

Coordinates (P, eₓ, e_y, hₓ, h_y, m) slow + true longitude L fast, **L as the
independent variable** (the role Sundman τ plays now; the paper's conclusion
flags longitude-as-variable as their own next step). Slowly-varying states need
only ~10–20 nodes/rev vs ~130 Cartesian:
- 1 N: 74.5 revs × 15 ≈ 1.1k nodes. 0.1 N: 754 revs × 15 ≈ **11k nodes —
  tractable** (~90k vars, same order as the Cartesian 1 N case).

Build (reusing the entire campaign architecture — only the dynamics/clock
change):

| module | role |
|---|---|
| `casadi_lt_mee.m` | NEW solver core: Gauss variational equations (paper p.6, already transcribed in DESIGN.md §2's source), independent variable L, dt/dL = 1/L̇ carried like the Sundman clock, t as a state, fixed L_f span with a cScale-style slack ONLY if needed (prefer none — Campaign B's lesson), modes fixedtf/mintime, ε-homotopy objective, cone-eliminated [α;s] in the RTN frame (q,s,w — the paper's control frame) |
| `mee_seed.m` | seed by element-space propagation of the tangential-thrust law (elements drift smoothly — seeds are nearly defect-free by construction) |
| `elements_to_cart.m` / `cart_to_elements.m` | already built + roundtrip-tested to 1e-10 (Task 3) — used for BCs and cross-checks |
| `run_mintime` / `run_transfer` / `run_ctf_sweep` | gain a `formulation` switch ('cartesian'|'mee') or thin MEE siblings — decide at plan time; prefer siblings to avoid destabilizing certified drivers |
| `verify_pmp_2body` | MEE variant of the switching-function check (Campaign B's front-end pattern) |

**Validation before use (non-negotiable):** solve the 10 N / c_tf=1.5 case in
MEE and require m_f within ~0.5 kg of the Cartesian 1376.74/1377.05 kg and the
same burn structure — a cross-FORMULATION check stronger than any mesh check.
Then (if Phase 1 delivered them) cross-check 5 N and 2.5 N too.

**The ladder:** {1, 0.5, 0.2, 0.1} N, each: MEE min-time anchor → c_tf=1.5 fuel
solve → structure counts. N ∝ revs; iteration budgets ∝ N; overnight jobs at
the bottom end. Per-point resume, certified-only caching, refinement spot-check
(N→1.5N) at 1 N and 0.1 N.

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
