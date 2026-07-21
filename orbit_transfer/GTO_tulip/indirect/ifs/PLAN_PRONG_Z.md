# Prong Z execution plan — ZTL: the Zhang recipe, run whole (2026-07-12)

Detailed, self-contained execution plan for Prong Z of `PLAN_OF_ATTACK_3.md`.
**User decision 2026-07-12: go straight to Prong Z; Prong B is SHELVED (not
deleted — its doc `PLAN_RUNG_B.md` stands if Z fails).** Written for a future
session to pick up cold. Read `PLAN_OF_ATTACK_3.md` §1–2 first.

Suggested home for the build: `../ztl/` (sibling of `ifs/`, `PSR/`).

## 0. The recipe, precisely (what Zhang 2015 actually does)

Zhang converges 0.6 N / ~150 revs / ~150 switches with single shooting
(7 unknowns) by composing THREE phases — note carefully which family is
laddered and where the ε-march happens:

1. **Min-TIME thrust ladder** (down from high thrust): robust family; gives
   t_f,min(T) and a costate backbone at every rung.
2. **Min-ENERGY thrust ladder** (down from high thrust, fixed
   t_f = c_tf × t_f,min(T) per rung): SMOOTH strictly-convex family, wide
   basin, no switches to track — each rung warm-starts λ0 from the previous.
3. **Bertrand–Épénoy ε-march 1 → 0, ONCE, at the TARGET thrust**, seeded by
   the min-energy solution at that same thrust; quadratically clustered
   schedule ε_j = (j²−1)/(N²−1), N=10 (dense near ε=0). The hard bang-bang
   solve happens exactly once, at the bottom, warm.

Robustness substrate under all three: exact variational-equation STMs (with
saltation only at ε=0), machine-precision event-detected regime boundaries
(no step ever smears a switching layer), RelTol=AbsTol≈1e-14, and t_f as the
escape knob ("varying c_tf is a practical method to overcome convergence
problems").

## 1. Post-mortem: the 2025 ladder as the control experiment

`../../lowThrust_GTO_tulip/thrust_continuation_minfuel_indirect.m` failed
0-for-7 (top anchor res 6.4e83; fine 7.4%-step retry 0-for-5). Against the
recipe it made FOUR deviations, each now correctable:

| deviation | old ladder | Zhang / this plan |
|---|---|---|
| family laddered | **min-FUEL costates laddered down** (Phase 3 warm-started the near-bang problem across thrust) | ladder the SMOOTH families only; ε-march once at 25 mN |
| smoothing family | logistic/tanh `u=(1−tanh(S/2ε))/2` — never saturates exactly, no regime boundaries, stiff layers integrated blind | BE ramp `u*=clamp(1/2−S/(2ε),0,1)` — EXACT saturation events at S=±ε, integrable as a 3-regime automaton |
| derivatives | complex-step **through ode113** (differentiates the adaptive stepper) | variational ODE Φ̇=A(y)Φ integrated alongside y; saltation Ψ at ε=0 switches |
| march discipline | thrust steps 37%, ad-hoc ε schedule | ≤10% thrust steps; quadratic ε clustering; halve-on-fail; t_f knob |

The old Phase 1 (min-time ladder) DID converge at every rung — the tfMin
table [6.0651, 6.0651, 6.0794, 6.0842, 6.0842, 6.1073, 6.2907] for thrust
multipliers [8, 5, 3.2, 2.1, 1.5, 1.2, 1]× lives in
`thrust_continuation_results.mat`. That is a validated asset (Z2 is mostly a
port), and it localizes the old failure to Phases 2–3, exactly where the
deviations sit.

## 2. Formulation decisions (locked before coding)

1. **Time-domain Cartesian** (Zhang fidelity), NOT Sundman. Augmented state
   y = [r; v; m; λr; λv; λm] ∈ R¹⁴, autonomous CR3BP. The perigee
   sensitivity that motivated Sundman is mild at high thrust and is handled
   at low thrust by exact STMs + tight tolerances (Zhang did 150 revs this
   way). Sundman-ZTL is the named fallback if Z3's bottom rungs choke on
   perigee stiffness — do not build it preemptively.
   Reuse: `lt_pmp_eom_minfuel.m` already has the analytic G (gravity
   gradient) and Hc (Coriolis) blocks and the min-fuel costate equations —
   the ZTL EOM is that file with the throttle law swapped (see 2) and an
   A(y) output added (see 3).
2. **The BE ramp family for everything smooth.** J_ε = ∫[u − ε u(1−u)]
   (Tmax/c) dt; interior minimizer u* = 1/2 − S/(2ε), clamped to [0,1];
   S = 1 − ‖λv‖c/m − λm unchanged. Three regimes with EXACT boundaries:
   - On: S ≤ −ε → u=1 (hard-burn EOM)
   - Medium: |S| < ε → u = 1/2 − S/(2ε) (u depends on y through S)
   - Off: S ≥ +ε → u=0 (coast EOM)
   ε=1 IS min-energy (J = ∫u²·(Tmax/c)dt up to an affine shift). ε=0 has no
   Medium regime — pure bang-bang. u is continuous for ε>0 ⇒ NO saltation
   during the march (Φ is continuous; only A jumps at regime events);
   saltation appears ONLY in the final ε=0 solve.
   **Do not reuse the tanh law anywhere in ZTL.** (Two-ε-family warning of
   `PLAN_RUNG_B.md` §2 applies; the direct side's eps and this eps are the
   same BE family, so PSR's epsMin>0 products are direct cross-checks.)
3. **A(y) = ∂f/∂y via complex step OF THE EOM FIELD, per evaluation** —
   machine-precision, ~15 lines, no tensor algebra. The variational system
   Φ̇ = A(y(t))Φ (14+196 = 210 ODEs) is integrated alongside y; this is
   "exact derivatives" in the sense that matters (the STM solves the true
   variational ODE to integrator tolerance; nothing differentiates through
   an adaptive stepper). Analytic A (CR3BP third-derivative tensor
   contracted with λv) is a later optimization ONLY if profiling demands.
   Within Medium, ∂u/∂y = −(1/2ε)∂S/∂y enters A automatically via the CS.
4. **Events and saltation.**
   - Event functions: S(y)∓ε (marches), S(y) (ε=0). Newton on the event
     time with ANALYTIC Ṡ (chain rule off the EOM:
     Ṡ = −(c/m)(λvᵀλ̇v)/‖λv‖ + (c‖λv‖/m²)ṁ − λ̇m), bisection fallback;
     locate to |S∓ε| ≤ 1e-13. Integrate regime-by-regime — never step
     across an event (the 3-regime automaton).
   - ∂S/∂y analytic: ∂S/∂m = c‖λv‖/m², ∂S/∂λv = −c·λv/(m‖λv‖),
     ∂S/∂λm = −1, zeros elsewhere.
   - Saltation at an ε=0 switch (f jumps): Ψ = I + (f⁺−f⁻)(∂S/∂y)/Ṡ⁻,
     composed into Φ across the event. GRAZING GUARD: if |Ṡ⁻| at any event
     falls below a floor (tune ~1e-3 of typical |Ṡ|), flag the solve, do
     not trust Ψ, and use knob (d) (t_f sidestep) — this is Zhang's known
     weakness and ours.
5. **Unknowns and residual (single shooting).** z = λ0 ∈ R⁷.
   R(z) = [rv(t_f) − rv_tgt (6); λm(t_f) (1)] ∈ R⁷; t_f fixed
   = factor × t_f,min(Tmax); m(t_f) free (consistent with λm(t_f)=0).
   Exact Jacobian: rows of the composed Φ(t_f, 0) = Φ_n Ψ_n … Ψ_1 Φ_1,
   columns 8–14 (∂/∂λ0).
6. **Solver.** `lsqnonlin` Levenberg–Marquardt with
   `SpecifyObjectiveGradient=true` (the old campaign measured LM 4× better
   than dogleg here); port the equilibrated truncated-SVD stepper from
   `ifs_solve2` as the fallback when LM damping crawls. Convergence gate
   ‖R‖ ≤ 1e-10 (honest against 1e-13/1e-14 integration).
7. **Integrator.** `ode89` (or ode113), RelTol 1e-13, AbsTol 1e-14, per
   regime segment with terminal events. Cost estimate at 25 mN: 210 ODEs
   over t_f ≈ 7.2 ND with ~40 perigee passes — expect minutes per
   flow+STM; LM warm from a rung neighbor should need O(10) iterations.
   High rungs are seconds. Run long solves via nohup/background per
   `matlab-headless` discipline.
8. **Endpoints/scaling**: lift verbatim from
   `thrust_continuation_minfuel_indirect.m` (muStar/lStar/tStar, GTO rv0,
   tulip max-ydot rvf, m0=15 kg, c; Tmax(mN) → ND accel conversion).
   factor = 1.15 throughout the ladder (PSR's certified 1.15× is the
   bottom oracle).

## 3. Preflight P0 — one session, before any building

Three cheap measurements that de-risk the build order:

- **P0a — graze margin at the target.** From `../PSR_data/
  psr_data_tf1p150_sw*_minEps0*.mat` (costate + pmp products): compute
  min over the ~25 crossings of |dS/dσ|. This pre-registers the Z4/Z5
  saltation risk (R3 below) with a number, before we build Ψ.
- **P0b — confirm the min-time backbone reproduces.** Load
  `thrust_continuation_results.mat`, re-run `solve_tfmin_indirect` at the
  top rung (200 mN) from the stored costates; confirm ‖R‖ and tfMin match.
  Validates the Z2 port target.
- **P0c — free look at the top anchor.** Run the EXISTING
  `solve_energy_indirect` (CS-Jacobian, no STM) at 200 mN,
  tf = 1.15 × 6.0651, seeded from the rescaled min-time costates. At ~2
  revs the old machinery may simply converge — if it does, Z3's top rung
  is done before ZTL exists, and it hands Z0–Z1 a converged ground-truth
  arc for unit tests (the same trick `test_ifs_residual` used).

## 4. Increments and gates

### Z0 — regime-explicit EOM + variational STM (`ztl_eom.m`, `ztl_flow.m`)
- `ztl_eom(t, y, P, regime)`: hard-coded u per regime (1 / ramp / 0);
  outputs yDot, and (nargout>1) S, Ṡ, u, Ht. P carries Tmax, c, muStar, ε.
- `ztl_A(y, P, regime)`: 14×14 ∂f/∂y via complex step of `ztl_eom`.
- `ztl_flow(y0, tspan, P)`: the 3-regime automaton — integrate y (+ Φ,
  optional) segment-by-segment with S∓ε (or S) terminal events, Newton
  event polish, Ψ composition at ε=0 switches; returns yf, Φ, event log
  (times, regimes, Ṡ at each event), dense output handle.
- **Gate Z0:** (i) with ε=1 and no events, Φ from the variational system
  matches complex-step-through-integrator to rel ≤1e-6 on a 1-rev arc
  (tolerance-limited); (ii) hard-throttle arc EOM matches
  `lt_pmp_eom_minfuel` at ε→0 limits (u pinned) to machine precision;
  (iii) event times reproducible to 1e-12 under RelTol change 1e-12→1e-13.
- **Side deliverable:** `ztl_A` + per-arc variational flow retrofit into
  `ifs_residual` later (the never-built cheap-Jacobian lever) — do NOT do
  it now, just keep the interface compatible (14-state layout matches).

### Z1 — shooting residual + exact Jacobian (`ztl_residual.m`, `ztl_solve.m`)
- `ztl_residual(lam0, P)` → [R (7×1), J (7×7), info(events, m_f, dV)];
  J from composed Φ/Ψ.
- `ztl_solve(lam0, P)`: LM w/ exact J; tsvd fallback port.
- **Gate Z1:** (i) on the P0c converged arc (or any converged smooth
  solve): R at the solution ≤ integration tol; (ii) J vs complex-step of
  the full residual: rel ≤1e-6 WITH at least one saturation event inside
  the span (this is the test that catches a wrong Ψ/event derivative);
  (iii) at ε=0 on a manufactured 1-switch arc: J vs CS across the switch.

### Z2 — min-time ladder port (`ztl_mintime_ladder.m`)
- Port old Phase 1 (validated): up/down continuation of
  `solve_tfmin_indirect` across the rung set; store tfMin(T) + costates.
- Rung set: 200 → 25 mN in ≤10% steps: [200, 180, 162, 146, 131, 118,
  106, 96, 86, 78, 70, 63, 57, 51, 46, 41, 37, 33, 30, 27.5, 25] mN
  (geometric ~0.9; ~21 rungs), halve a step on failure.
- **Gate Z2:** converged min-time at every rung; tfMin curve smooth and
  matching the old 7-point table where they coincide.

### Z3 — min-energy ladder (`ztl_energy_ladder.m`) — the anchor chain
- At each rung k: t_f = 1.15 × tfMin(T_k); solve ε=1 with `ztl_solve`;
  seed = previous rung's λ0 (top rung: rescaled min-time costates, or the
  P0c solution). Save per-rung .mat (λ0, R, events, dV, throttle profile).
- Saturated arcs will appear as thrust drops (u hitting 0/1 at ε=1) —
  handled automatically by the automaton; no structure bookkeeping.
- **Gate Z3 (per rung):** ‖R‖ ≤ 1e-10. **Bottom gate (25 mN): compare
  against the DIRECT energy backbone `energy_f1150.mat`** — dV, m_f, and
  throttle profile u(t) vs the direct s(σ) (map σ→t via the direct κ);
  agreement to direct-mesh accuracy (~1e-3 rel) certifies the whole chain.
  **This is the first falsifiable "Zhang works here" checkpoint: a
  converged 40-rev indirect solve, smooth family, at nominal thrust.**
- If a rung fails: halve the thrust step; if still stuck, knob (d) — nudge
  factor ±0.01 at fixed thrust, re-ladder factor back at the lower thrust
  (2-D continuation, Zhang's exact escape).

### Z4 — ε-march at 25 mN (`ztl_eps_march.m`) — the campaign prize
- From the Z3 bottom anchor: ε_j = (j²−1)/(99), j = 10…1
  (1, 0.798, 0.616, 0.455, 0.313, 0.192, 0.091, 0.030, 0) — insert
  midpoints on failure (never more than halving in ε-space); each step
  seeded from the last converged λ0.
- The ε=0 endpoint uses saltation-composed Jacobians and the S-event
  automaton; watch the graze guard (P0a told us the expected margin).
- **Gate Z4:** ε=0 converged (‖R‖ ≤ 1e-10) + certified (sign law S<0 ⇔
  u=1 along the whole trajectory — reuse `ifs_certify` logic) + **oracle
  match vs PSR certified 1.15×**: switch count (25), ΔV and prop mass to
  ~1e-3 rel, switch times to direct-mesh accuracy. Success = **the first
  fully converged indirect min-fuel at nominal thrust — the goal IFS was
  built for.** Update `LOW_THRUST_MINFUEL_CAMPAIGN.md` + memory
  immediately.
- Optional (cheap, high value): also stop the march at ε = 0.5 and 0.1
  and compare against PSR's epsMin=0.5/0.1 products — a family-level
  cross-validation of both pipelines at three points, not one.

### Z5 — the band (`ztl_tf_march.m`)
- From the converged 25 mN / 1.15× bang-bang point: march factor DOWN
  1.15 → 1.01 in Δ ≤ 0.01 steps at ε=0, warm λ0, structure re-discovered
  by event detection each solve (no surgery). Near a fold / structure
  change: expect a graze (Ψ blowup) — the guard flags it; response is
  (i) smaller Δfactor, (ii) a temporary ε lift (solve at ε=0.03, step
  factor, re-drop — a local homotopy detour), (iii) pseudo-arclength in
  (λ0, factor) with `ifs_tf_arclength`'s machinery, now off a
  NON-degenerate anchor.
- **Gate Z5:** each converged factor extends the certified set below
  1.12×; the deliverable is the reachable set + a characterized wall if
  one exists (its location and mechanism — graze density vs conditioning —
  is publishable either way, and directly comparable to the direct side's
  1.01–1.11× energy-seed wall).
- Bonus symmetric run: march factor UP toward 1.95× and beyond 25-switch
  structures to cross-check the whole PSR atlas indirectly.

## 5. Where IFS fits (unchanged from PLAN_OF_ATTACK_3 §4)

Any converged ZTL point hands IFS a dynamically consistent
(λ0, switch times, arc structure) — the seed class Rung A proved cannot come
from direct data. Use `ifs_seed`-layout packing + `ifs_solve2`/`ifs_certify`
as the independent multiple-shooting verifier at Z4/Z5 milestones, and as
the robustifier if single shooting gets sensitive deep in the band. Z0's
variational flow is the future cheap-Jacobian retrofit for `ifs_residual`.

## 6. Budget, risks, stopping rule

Budget: P0 one session. Z0+Z1 one to two sessions (the automaton + Ψ tests
are the care points). Z2 half (a port). Z3 one session + ladder compute
(background). Z4 one session + march compute. Z5 open-ended research with
per-step payoff. Total to the Z4 prize: ~4–6 focused sessions.

- **R1 — Z3 bottom rungs choke on 40-rev perigee stiffness** (time-domain
  integration cost or event misses). Mitigation: tighter tolerances first;
  then the named fallback: Sundman-ZTL (swap independent variable; the
  automaton and STM machinery carry over; `sms_eom` has the Sundman
  Hamiltonian pattern).
- **R2 — ε-march crawls near ε=0 despite events.** Distinguish from
  ms_band's crawl: there the Jacobian was mesh-noisy CS-through-integrator
  and the seed structure-mismatched; here both are clean, so a persistent
  crawl is REAL information. Capture the trace (cond(J), step-vs-linear
  regime, which rows dominate) before reacting; knob (d) first.
- **R3 — grazing switches at Z4/Z5** (the 1/Ṡ wall). P0a measures the
  1.15× margin up front. Response ladder: Δ-step reduction → ε-lift detour
  → pseudo-arclength → record the wall honestly.
- **R4 — single-shooting basin too small even warm at 40 revs.** The
  Zhang-faithful counter is that rung spacing controls seed distance;
  if a rung fails at ≤2% thrust steps AND knob (d) fails, hand the rung to
  IFS (multiple shooting) seeded from the last converged neighbor — the
  hybrid that neither method could attempt alone.
- **STOPPING RULE:** if Z3's bottom gate (indirect energy at 25 mN) fails
  with Z0/Z1 unit tests green and ≤2% rungs — i.e., the SMOOTH problem
  resists a warm, exact-derivative, event-exact single shooter at 40 revs —
  that falsifies the audit's explanation and is the strongest evidence yet
  for the regularized-coordinates thesis (Leomanni 2021 digest,
  `PLAN_OF_ATTACK.md` §2). Stop, record, and decide the rebuild
  deliberately. Do not grind past a diagnosed structural failure.

## 7. Pointers

- Parent plan + Zhang audit: `PLAN_OF_ATTACK_3.md` §2.
- Zhang digest: `PLAN_OF_ATTACK.md` §2; PDF
  `../min_fuel_papers/2015-8-J-LowThrustMinimumFuelOptimizationRTBP.pdf`.
- Reusable old-campaign code (`../../lowThrust_GTO_tulip/`):
  `lt_pmp_eom_minfuel.m` (analytic G/Hc + costate eqs; throttle law to be
  REPLACED), `solve_tfmin_indirect.m`, `solve_energy_indirect.m`,
  `shoot_residual_*.m` (CS pattern), `thrust_continuation_minfuel_indirect.m`
  (endpoints/scaling block + Phase-1 min-time loop; its Phases 2–3 are the
  control experiment, not a template).
- Solver fallback to port: `../ifs/ifs_solve2.m` (equilibrated tsvd step).
- Certification pattern: `../ifs/ifs_certify.m` (sign law), PSR's
  `psr_ipopt_certify.m` (direct-side local-min certificates for the oracle
  files).
- Oracles on disk: `../PSR_data/psr_data_tf1p150_*` (bang-bang + epsMin
  variants), `../sundman_minfuel/results/energy/energy_f1150.mat` (direct
  energy at 25 mN), `thrust_continuation_results.mat` (min-time backbone).
