# Elliptic→GEO with Lunar Gravity — Phase-0 Design Spec

**Date:** 2026-07-22  **Status:** Phase 1 implemented through the 10 N rung
(plan 2026-07-22-elliptic-geo-cr3bp-phase1.md); deep rungs + phi0 sweep
pending.
**Campaign:** `orbit_transfer/earth_elliptic_to_geo_CR3BP/` (Goal 2 of the
2026-07-21 orbit-transfer goals discussion).
**Companion theory:** Bonnard–Caillau–Picot 2010,
`earth_elliptic_to_geo_CR3BP/papers/Geometric_And_Numerical_Techniques_In_3_Body_Low_Thrust_Transfers.pdf`
(μ-continuation provenance + convergence theory; conjugate-point
certification; σ-bound rev-count law).

## 1. Goal

Re-solve the certified HMG-2004 elliptic→GEO min-fuel campaign (1500 kg,
P⁰ = 11625 km, e⁰ = 0.75, i⁰ = 7° → equatorial GEO; thrust ladder
10 → 0.1 N) **with the Moon's gravity included**, and answer quantitatively:
**how much does lunar gravity move the certified 2-body answers** — final
mass, switch structure (as mesh-bands), and the ladder laws (R0,
T·t_f,min)? Production driver: a reusable solver + documented drivers that
reproduce across the thrust range, not a one-off solve.

## 2. Decisions (locked, with rationale)

- **D1 — Formulation: Earth-centered MEE/L-domain + lunar third-body
  perturbation** (option (a), user-approved 2026-07-22). Keeps the certified
  2-body machinery (free ΔL, L-domain mesh, trivial GEO terminal) and adds
  the Moon as a perturbing acceleration. Rotating-frame Cartesian (option b)
  rejected: loses free ΔL (the tulip topology wall), makes GEO a moving
  target, and suits Moon-bound targets — ours is Earth-bound throughout.
- **D2 — Moon plane = reference plane** (user-approved 2026-07-22). Faithful
  CR3BP idealization: the Moon's circular orbit defines the equatorial
  plane; the GTO is inclined i⁰ = 7° to it; GEO lies in it. The real
  Moon-equator obliquity (~18–29°) is explicitly NOT modeled — recorded as
  the idealization; an inclined-Moon fidelity step is out of scope (§10).
- **D3 — Solver reuse: extend the shared 2-body core, opt-in**
  (user-approved 2026-07-22). `casadi_lt_mee`/`lt_mee_rhs` in
  `earth_elliptic_to_geo/direct/core/` gain an OPTIONAL third-body
  perturbation input; absent ⇒ byte-identical nominal (the ladder-prep
  back-compat invariant, with the same nominal-regression gate). One solver
  serves both campaigns.
- **D4 — Comparison convention: same physical t_f as the certified 2-body
  rungs** (t_f = c_tf·t_fMin with the 2-BODY anchors, c_tf = 1.5
  convention). This makes Δm_f a controlled experiment (identical boundary
  conditions and time budget; only the dynamics differ). CR3BP-native
  min-time anchors are a later, separate result — not Phase 1.
- **D5 — Bridge strategy: μ_M-continuation on the ENERGY solution, then
  ε-sharpen to fuel at full Moon.** Dial the lunar mass 0 → μ_M in warm-
  started steps from the certified 2-body energy solution (the
  `gen_elfo_energy_gravhom` pattern = Bonnard–Caillau–Picot μ-continuation,
  §4.4), THEN run the standard energy→fuel ε-homotopy with the Moon fully
  on. Rationale: continuation lives on the provably benign energy landscape
  (their Thm 3.5); the bang-bang structure is resolved only once, at the
  final physics.
- **D6 — Lunar phase φ₀ is a first-class parameter, fixed for the
  baseline.** Baseline: φ₀ = 0 (Moon on the +x axis of the reference frame
  at t = 0). It enters every artifact fingerprint. A φ₀-sweep ("does lunar
  phasing measurably help?") is a recorded later experiment, not Phase 1.
- **D7 — Params home: campaign-local `lunar_params.m`** in
  `earth_elliptic_to_geo_CR3BP/direct/`, expressing the Moon constants in
  the 2-body campaign's canonical units (read from `kepler_lt_params` —
  conversion inside the function, physical values in the header).
  `cr3bp_common/` is NOT used: its normalization is the Earth–Moon
  rotating-frame scale set for the 25 mN tulip craft; this campaign is
  "2-body campaign + perturbation" and inherits kepler units.

## 3. Dynamics specification

**Perturbing acceleration (Earth-centered frame, direct + indirect terms):**

    a_M(r, t) = μ_M · [ (r_M(t) − r)/|r_M(t) − r|³  −  r_M(t)/|r_M(t)|³ ]

The **indirect term is mandatory** (the Earth-centered frame is
non-inertial because the Moon accelerates the Earth); dropping it is the
classic third-body bug. This is the exact Earth-centered restatement of
CR3BP gravity for a spacecraft of negligible mass.

**Moon ephemeris (circular, per D2):**

    r_M(t) = D_EM · [cos(n_M t + φ₀); sin(n_M t + φ₀); 0]
    n_M    = sqrt( (μ_E + μ_M) / D_EM³ )      (sidereal rate, ~27.32 d period)

**Physical constants** (converted to campaign canonical units inside
`lunar_params`): μ_E = 398600.4418 km³/s², μ_M = 4902.800 km³/s²
(ratio μ_M/μ_E = 0.012300, consistent with the CR3BP μ* = 0.0121506 via
μ*/(1−μ*)), D_EM = 384400 km.

**Coupling rules:**
- a_M is a **pure acceleration**: it enters the Gauss VOP equations exactly
  as the thrust acceleration does (same RTN resolution path in
  `lt_mee_rhs`), but contributes **nothing to ṁ** — no mass flow from
  gravity. This asymmetry must be explicit in the RHS.
- Time is already a state in the L-domain formulation (dt/dσ = ΔL/L̇,
  t(1) = t_f pinned), so the time-dependent r_M(t) needs **no structural
  change** — the RHS reads the t-state.
- **Hook contract (D3):** the core accepts an optional perturbation spec
  `pert = struct('muM', μ_M_canonical, 'DM', D_canonical, 'nM',
  n_canonical, 'phi0', φ₀)`; empty/absent ⇒ term compiled out entirely
  (byte-identical nominal). A continuation scale `pert.gain ∈ [0,1]`
  multiplies μ_M — this single knob IS the D5 bridge parameter.

## 4. Architecture / file plan

- `earth_elliptic_to_geo/direct/core/lt_mee_rhs.m`, `casadi_lt_mee.m` —
  opt-in `pert` input threaded through (D3); nominal regression gate.
- `earth_elliptic_to_geo_CR3BP/direct/`:
  - `lunar_params.m` — constants per §3 (+ header with physical values).
  - `bridge_mu_continuation.m` — D5 stage 1: load certified 2-body ENERGY
    solution at a rung, walk `pert.gain` 0 → 1 (adaptive step, halve on
    fail, the gravhom/tfsweep step pattern incl. `Solve_Succeeded` gates,
    rF fallback, checkpoint/resume).
  - `solve_cr3bp_minfuel.m` — D5 stage 2: ε-homotopy energy→fuel at
    `gain = 1` (reuses `homotopy_mee` with the pert threaded).
  - `compare_vs_2body.m` — the §1 deliverable: Δm_f, switch mesh-bands,
    per-rung table vs the certified 2-body values.
  - `sanity_bound.m` — §7 numbers from code, not prose.
- **Fingerprints:** reuse the earth campaign's `check_cache_fp` pattern;
  fp gains `muM, phi0, gain` fields so 2-body and CR3BP caches can never
  cross-seed silently.
- Results under `earth_elliptic_to_geo_CR3BP/direct/results/` (gitignored
  .mat, committed figures), README/TODO kept current.

## 5. Terminal set and conventions

Unchanged from the 2-body campaign (that is the point of D1/D4): GEO
terminal via the existing `geo_terminal` MEE conditions; same c_tf = 1.5;
same certified t_f values per rung; mass/Isp per `kepler_lt_params`
(1500 kg, Isp ≈ 2000 s class — read from the params file, do not restate).

## 6. Phase-1 plan preview (for the implementation plan)

1. Sanity bound (§7) — tabulated from certified data, BEFORE any solve.
2. Core opt-in extension + nominal byte-path regression.
3. Bridge at **10 N** (fast, cheap, Moon effect predicted ~0.1% of control
   authority): μ-continuation → ε-sharpen → certify (defect + PMP checks as
   in the 2-body campaign) → first Δm_f data point. Expect Δ ≈ tiny; that
   is the null-hypothesis validation of the machinery, not a disappointment.
4. Walk the ladder down (1 N, 0.2 N, …) with per-rung fingerprints; the
   Moon effect grows as t_f stretches toward and past the lunar month.
5. Comparison table + figures; README/TODO close-out.

## 7. Sanity bound (null hypothesis the solves must beat)

Lunar **tidal** acceleration at radius r (the indirect term cancels the
uniform part): a_tide ≈ 2 μ_M r / D_EM³. At GEO radius (42164 km):
**≈ 7.3×10⁻⁶ m/s²**. Against thrust authority T/m₀:

| rung | T/m₀ [m/s²] | lunar tide / authority |
|---|---|---|
| 10 N | 6.7×10⁻³ | ~0.11% |
| 1 N | 6.7×10⁻⁴ | ~1.1% |
| 0.2 N | 1.3×10⁻⁴ | ~5.5% |
| 0.1 N | 6.7×10⁻⁵ | ~11% |

Transfer durations follow t_f = 1.5·t_fMin(T) with T·t_fMin ≈ R0 (campaign
README: ≈ 850 N·h — `sanity_bound.m` must recompute from the certified
data, not trust this prose): 10 N ⇒ t_f ≪ lunar month (Moon nearly static,
phase matters); 0.1 N ⇒ t_f ≫ lunar month (Moon phase-averages; oscillatory
part largely cancels, secular part survives). **Prediction:** Δm_f
negligible at 10 N, sub-percent but resolvable at deep rungs; switch-count
bands may shift near thresholds. The campaign's result is the measured
deviation from this null model.

## 8. Validation gates

1. `pert` absent ⇒ 2-body solver byte-identical (regression on a certified
   rung reproduction, ladder-prep style).
2. `gain = 0` with pert PRESENT ⇒ matches the 2-body solution to solver
   tolerance (the hook itself introduces no drift).
3. Bridge at 10 N reaches `gain = 1` with certified (Solve_Succeeded +
   defect + PMP-consistency) steps throughout; ε-sharpen certifies at ε = 0.
4. Comparison table entries carry mesh-band switch counts (P0 protocol),
   never bare integers.

## 9. Risks

- **Certified-core edits (D3)** — mitigated by the opt-in invariant +
  gates 1–2; the ladder-prep package proved this pattern twice.
- **Deep-rung bridge cost** — 0.1 N solves are hours-long; μ-continuation
  multiplies solves. Mitigation: bridge coarse (few gain steps — the
  perturbation is ≤11% of authority), sharpen once; per-rung checkpoints.
- **Phase-dependence at fast rungs** — at 10 N the answer depends on φ₀
  (Moon effectively static). D6 pins φ₀ = 0; the fingerprint prevents
  silent mixing; the φ₀-sweep is the recorded follow-up.
- **Averaging optimism** — the §7 cancellation argument is heuristic; if
  deep-rung Δm_f comes out larger than predicted, that is a finding, not an
  error (house honesty rule).

## 10. Out of scope (recorded, not promised)

Inclined/ephemeris Moon, solar third body, SRP; CR3BP-native min-time
anchors; the φ₀ optimization sweep; **Phase 2 (indirect counterpart +
conjugate-point certification)** — already scoped in the campaign TODO,
gets its own spec after Phase 1 lands.
