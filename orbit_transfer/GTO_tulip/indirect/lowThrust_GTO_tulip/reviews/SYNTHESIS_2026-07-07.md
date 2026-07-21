# Three-way review synthesis — 2026-07-07

Reviewers: host Claude (Fable 5) + GPT-5.5 + Gemini 3.1 Pro (docs and NLP code).

## Documents (theory note + tutorial)

Convergent (all three reviewers) — all FIXED:
- **[critical] Theory §4 described the 4-control throttle transcription while the
  implemented/taught solver uses the 3-control fixed-throttle reduction** (nZ,
  cone, Jacobian blocks all inconsistent). Fix: §4 now presents the general
  (w, s) formulation with the ballast pitfall, then a "min-time specialization
  actually implemented" subsection (s ≡ 1 justified by S < 0, sphere constraint,
  corrected counts, warm-start bound-activity rationale). The general form is
  retained deliberately — min-fuel needs it.
- **[major] "five significant figures" agreement overclaim** → "mesh accuracy,
  ~3e-4 relative at N = 12000, improving under refinement"; "30,000 variables"
  → "up to 1e5".

GPT-5.5 unique — FIXED:
- ND exhaust velocity written dimensionally (c̃ = Isp·g0) → full formula
  Isp·g0·t*/(1000·l*) ≈ 20.24, with the 20.59 km/s dimensional value shown.
- Theory defect/tf-column formulas were uniform-mesh (h = tf/N, 1/(2N)) →
  generalized to h_k = tf·Δσ_k with the uniform special case noted; mesh
  paragraph now names density-matching.
- Terminal-manifold remark: "pin the normal components" → costate LIES IN the
  normal space; tangent (mass) component vanishes.
- "rendezvous with the tulip" → "fixed sampled state"; free-phase transversality
  noted. Normal-multiplier (lambda_0 = 1) convention stated. Covector-mapping
  caveat (quadrature weights, Δσ scaling) added. O(h³) local vs O(h²) global
  defect-rate wording. `if real(S) > 0` made explicit (hint + reference EoM).
  Reproducible rng(7) gradcheck recipe added to the Phase E hint.

Gemini unique — FIXED (amended):
- Tutorial Phase A signature had undefined `aThrust` → defined in the exercise.
- Gemini's proposed ballast-pitfall rewrite was REJECTED as wrong physics (with
  s fixed, mass flow is constant — no exploit exists); the pitfall stays with
  the general formulation where it is real.

Host unique — FIXED: GTO eccentricity 0.7245 → 0.7248.

## NLP code

No correctness findings (GPT-5.5 ran independent FD checks on A/B and the
gradient layout). Applied:
- Inlined triplet assembly in nlp_constraints (was per-block meshgrid helper);
  nargout guard so value-only calls skip Jacobian building (GPT+Gemini).
- Vectorized lt_dynamics Jacobian assembly (3D implicit expansion) (Gemini).
- build_guess: unique() guard on tangential abscissae; density_matched_mesh
  hardened (dedupe, nonzero-start normalization, monotonicity check); zero-norm
  guards on unit-vector normalizations via unit_columns() (GPT+Gemini).
- solve_tfmin_nlp: Z0/sigma assertions (GPT).
- Drivers: removed dead one-period GTO propagation (rv0 = rv0ND), fixed header
  field names (.devPos_km/.W), removed unused Z (GPT+Gemini).
- Declined: making the pumpkyn path configurable (deliberate, documented).

Re-verification: all five harvest checkpoints identical after the rewrite
(gradcheck 1.58e-9, CP-D exact 0). Note: N=3000 fmincon stall point moved
6.264259 → 6.264901 (~6e-4 scatter under numerically-equivalent
reformulations — flat valley at coarse-mesh accuracy; refinement re-run
numbers quoted in the docs with this caveat).

## Addendum: min-fuel variant (added post-review, same day)

Built after the review round. Direct side CONVERGES (arrival-leg replan,
burn+coast warm start, single-switch bang-bang, defects 2e-15, propellant
1.0622 vs 1.0650 kg reference); indirect polish OPEN (four-seed progression
||R|| 1.55/0.83/0.33/0.14 -- raw, rescaled, switch-anchored, LSQ costate
reconstruction). Failure modes measured and documented in theory S6 +
tutorial Phases G-H: costate scale anchor ("+1" in S), structural seed
mismatch, time-stretched warm starts, defect cliffs at pinned nodes,
lbfgs-stall multiplier quality. The min-fuel sections have NOT been through
an external review round yet. [superseded -- see round 2 below]

## Round 2: min-fuel code + doc review (same day, later)

Reviewers: host Claude + GPT-5.5 + Gemini 3.1 Pro (docs three-way; code
host + both externals). Raw outputs: minfuel_{code,doc}_{gemini,gpt55}_*.md.

### Code (9 min-fuel files) -- all FIXED

Convergent (host + GPT + Gemini, independently):
- **[critical] run_gto_tulip_minfuel.m: `solve_minfuel_indirect` called
  without the `m0` argument** -- every argument shifted one slot, runtime
  failure. Never surfaced because the pipeline driver
  (NLP_lowThrust_GTO_Tulip_minfuel.m) makes the call correctly and the
  standalone driver was never executed end-to-end. Fixed: pass `1`.
- costate_seed_from_nlp.m empty-index crashes (no node with s > 0.9;
  no throttle switch) -- guards added: descriptive error for no-burn,
  last-burn-node gauge fallback for no-switch (GPT + Gemini).

GPT-5.5 unique -- FIXED:
- solve_minfuel_indirect.m `nargin < 8` guard off by one (epsSchedule is
  the 9th arg; documented 8-arg call would error) -> `nargin < 9`.
- Fractional `nCoasts` (half-integer when the arc ends mid-coast) ->
  count maximal S > 0 runs, both drivers.
- "EXACTLY feasible" warm-start overclaim in the NLP driver header (the
  throttle clip makes it ~1e-4) -> honest wording.
- m0 (kg) vs mass-fraction naming in run_gto_tulip_minfuel -> m0kg.
- Declined: zero-norm guard on the primer direction in
  lt_pmp_eom_minfuel.m (same accepted pattern as min-time lt_pmp_eom.m);
  replacing the clipped warm start with exact s = 1/0 controls
  (deliberate -- interior-point needs strictly interior guesses).

Host unique -- FIXED: dead original-target block (getTulip + prop) in the
NLP driver (rvf overwritten by the phase-shifted target; lint-flagged).

Both externals explicitly verified correct: min-fuel PMP equations +
switching function + costate ODEs, complex-step safety, throttle-dynamics
Jacobians A/B, 154/segment + 4/node triplet layout, cone gradient,
costate-reconstruction algorithm. Post-fix: lint clean, nCoasts unit
tests pass, shoot_residual_minfuel smoke test (corrected arg order)
finite 7x1/7x7, both seed-guard paths exercised in MATLAB.

### Docs (theory S6 + tutorial Phases G-H) -- no CRITICAL findings

All quoted numbers verified independently by both externals (Gemini
re-derived the N=5 gradcheck dimensions 41/66/304 and the mass/time/dV
chain; GPT confirmed the results table against ground truth).

GPT-5.5 (5 MAJOR, precision/pedagogy) -- FIXED:
- "+1 anchors the scale" keyfact conflicted with S3's normal convention
  (H = 1 + lambda'f, lambda_0 = 1) -> scoped: min-time CONTROL LAW is
  scale-invariant, H(tf) = 0 pins the scale without feeding back into
  controls; min-fuel makes the throttle law itself scale-sensitive.
- "we measured both methods failing" on the full spiral lacked the
  direct-side evidence -> added (NLP grinds in feasibility mode;
  covector-seeded shooting breaks the integrator at a perigee dive).
- Phase G checkpoint used full-transfer shooting with no setup ->
  concrete inputs added (Phase C endpoints, m0 = 1, tf = 1.2 x 6.2907,
  z(1:7) seed); solve_minfuel_indirect.m added as exercise item 3.
- costate_seed_from_nlp.m used in Checkpoint H2 but never assigned ->
  added as Phase H exercise item 4 (signature + the rows to stack);
  H2 now references it instead of re-describing the algorithm.
- What-comes-next item 1 overpromised convergence -> "standard
  escalations, neither verified here; closing ||R|| = 0.14 is the open
  exercise".
- MINOR: intro-table driver row split and named; control 4-vector no
  longer called u (collides with PMP throttle); "lighter is better" ->
  dumping mass buys acceleration though the objective charges for it.

Gemini (1 MINOR, good catch) -- FIXED: Phase G said "ADD the running
cost to Ht" -- a literal learner copying lt_pmp_eom.m would keep the
min-time bare 1 -> "REPLACE the 1 with (Tmax/c)u".

Both PDFs recompiled clean (triple-pass, no undefined refs, no >40pt
overfulls), aux cleaned. Checkpoint numbers unchanged (no fix touched a
verified quantity).
