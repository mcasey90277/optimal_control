# earth_elliptic_to_geo — Design Spec

**Goal.** Reproduce the low-thrust **minimum-fuel** Earth orbit transfer of Haberkorn,
Martinon & Gergaud, *"Low thrust minimum-fuel orbital transfer: a homotopic
approach,"* JGCD 27(6), 2004 (`orbit_transfer/min_fuel_papers/Gergaud-Haberkorn-Martinon-JournalGuidance2004-preprint.pdf`)
— transfer a 1500 kg satellite from a low, elliptic, inclined orbit to equatorial
GEO — using **our direct (collocation) machinery**, whereas the paper solves it
**indirectly** (single shooting + homotopy).

Date: 2026-07-16. Status: **implemented (2026-07-17), T1–T16 complete.** Isp
pinned at 2000 s (default, §7 item 1; validated by the M2 mass match, not an
independent citation — the intended ref [6] PDF was unreadable). Headline
result (M2, free-longitude insertion manifold, 10 N, c_tf=1.5): m_f = 1376.74 kg
(N=600) / 1377.05 kg (N=1200), within the plan's [1355,1385] kg gate band, ~2 kg above the paper's Fig-23 reading (1370–1375 kg); Δmesh
0.31 kg. M3 front (c_tf 1.2→3.0, 10 N) monotone 1360.37→1386.81 kg matching
Fig 23's shape; thrust law C≈846.6 N·h at 10 N vs paper's ≈850 N·h, but the
5 N/2.5 N law leg is BLOCKED (defect floor ~5e-3, false-infeasibility
signature, nodes-per-rev hypothesis — future work). Open, non-gating finding:
a reproducible ~20° PMP primer/dual misalignment on M0/M1 (primal
certification unaffected). Full record: `README.md`, `PLAN.md`,
`.superpowers/sdd/progress.md`.

---

## 1. Why this is a natural fit

The paper's method **is our method, done indirectly**:

| paper (indirect) | ours (direct) |
|---|---|
| min ∫‖u‖dt, fixed t_f = c_tf·t_{f,min} | same |
| solve min-time first for t_{f,min} | same (our min-time anchor) |
| energy→mass homotopy ∫λ‖u‖+(1−λ)‖u‖² (λ:0→1) | energy→fuel homotopy ∫s−ε∫s(1−s) (ε:1→0) — **identical**, opposite parametrization |
| bang-bang via switching function ψ, no a-priori switch count | bang-bang via ε→0, no a-priori switch count |
| PMP costates from shooting | PMP costates from KKT duals (verification) |
| empirical law T_max·t_{f,min} ≈ C | we check the same law |

So the work is: **keep our recipe, swap the dynamics (CR3BP → 2-body Kepler) and the
terminal condition (fixed point → orbit-insertion manifold).**

**Decisions locked (brainstorm):** (1) *2-body physics via our CR3BP machinery* —
reuse the Sundman/cone/ε/fixed-t_f architecture, Earth-only 1/r² gravity, no Moon,
no rotating frame; (2) *scope Path 1* — validate + modest front, 10 N → ~1 N; (3)
*fixed terminal first, then free-L_f*.

---

## 2. Problem statement (as we will pose it)

**Inertial, Earth-centered, 2-body + low thrust**, nondimensionalized so μ=1 and
states are O(1). State `x = [r(3); v(3); m; t]` (8), control `u = [α(3); s]` (4),
cone-eliminated: thrust = s·T_max·α/m, ‖α‖=1, s∈[0,1].

- ṙ = v
- v̇ = −μ r/‖r‖³ + (T_max/m) s α
- ṁ = −(T_max/c) s,  c = Isp·g₀   (exhaust velocity)

**Sundman regularization** (the classic Kepler transform, ideal here): dt/dτ = κ = ‖r‖^p,
p = 1.5; t carried as the 8th state; τ_f **fixed** (= seed value); transfer time enforced
by the terminal time-state condition t(τ_f) = t_f. (Same trick as `casadi_minfuel_sundman`;
avoids the dense-KKT-column blow-up of a free τ_f.)

**Objective (energy→fuel homotopy):** J(ε) = ∫ s dt − ε ∫ s(1−s) dt, swept ε:1→0.
ε=1 is min-energy (∫s² dt, strictly convex, smooth ramp, big basin); ε=0 is min-fuel
(∫s dt, bang-bang). Δm and thus m_f come from the flown mass state.

### Boundary conditions

Initial elements (paper): P=11625 km, (eₓ,e_y)=(0.75,0), (hₓ,h_y)=(0.0612,0), L=π,
m=1500 kg. With eₓ=0.75,e_y=0 ⇒ Ω+ω=0; hₓ=0.0612 ⇒ i≈7°, Ω=0 ⇒ ω=0; L=π ⇒ θ=π ⇒
**start at apogee** (r_a≈46 500 km, r_p≈6 640 km). Converted to inertial (r,v) by
`elements_to_cart`.

Terminal = GEO: **P=42165 km**, (eₓ,e_y,hₓ,h_y)=0, L free, m free.
*NB the paper is internally inconsistent: the p.5 BC table says 42125 km but the
p.6 formulation (P_mf) says 42165 km; GEO radius is 42164.17 km, so we adopt
42165 (the p.5 value is a typo; ref [6]'s benchmark also uses 42165).*
- **M1 (fixed):** full r,v pinned at a **prescribed** L_f (reuses our fixed-rvf
  terminal). Choose L_f from the paper's own empirical laws so M1 is directly
  comparable: R1 gives (L_{f,min}−L₀)·T_max ≈ 264 rad·N ⇒ at 10 N,
  L_{f,min}−L₀ ≈ 26.4 rad; R2 gives c_{Lf,opt} ≈ 1.12·c_tf + 0.09 ≈ 1.78 at
  c_tf=1.5 ⇒ **L_f − L₀ ≈ 47 rad ≈ 7.5 revolutions** — reproducing Table 3's
  rev count by construction, so the switch count (~18) and apogee-burn
  structure are the genuine tests.
- **M2 (free-L_f insertion manifold):** 5 scalar constraints, phase free —
  equatorial (r_z=0, v_z=0), circular radius ‖r‖=a, circular speed ‖v‖=√(μ/a),
  zero radial rate r·v=0. (Equivalent to P=42165, eₓ=e_y=hₓ=h_y=0.) The set
  also admits the retrograde orbit (h_z<0); the prograde seed selects the right
  branch — add an h_z>0 guard only if the solver ever drifts.

### Constants

μ_⊕ = 398600.47 km³/s² (standard); m₀ = 1500 kg; T_max ∈ {10, 5, 2.5, (1)} N;
c_tf = 1.5 (primary). Control frame: **inertial α** (physically identical to the
paper's radial/transverse/normal q,s,w — just a representation choice; inertial is
simpler for a Cartesian solve). **Isp: OPEN** — not stated numerically in the paper;
pin from ref [6] (Caillau & Noailles 2001, same benchmark) before coding the M2
mass-match; default 2000 s. Trajectory *structure* (revs, switches, apogee-burns) is
Isp-independent, so M0/M1 don't need it; only M2's m_f number does. Sanity check on
the default: m_f≈1370 kg at 10 N/1.5 ⇒ c = ΔV/ln(m₀/m_f) with ΔV ≈ 1.8–2 km/s
(perigee raise + 7° plane change) ⇒ Isp ≈ 1900–2250 s, so 2000 s is consistent
with the paper's own numbers.

---

## 3. Module architecture

New top-level study folder `optimal_control/earth_elliptic_to_geo/` (sibling of
`lambert/`, `NLP_lowThrust_GTO_tulip/`). Each module is a small, single-purpose unit.

| module | role | depends on |
|---|---|---|
| `kepler_lt_params.m` | constants (μ, m₀, Isp, g₀, T_max), canonical nondim units (LU, TU so μ=1), Sundman p | — |
| `elements_to_cart.m` | (P,eₓ,e_y,hₓ,h_y,L) → inertial (r,v); + inverse for reporting osculating elements | params |
| `geo_terminal.m` | terminal builders: fixed rvf at L_f (M1) and the 5-constraint insertion manifold (M2) | params |
| `casadi_lt_2body.m` | **the solver core** (single mode-switched file — planning-time DRY consolidation of the two spec'd siblings): Sundman clock via cScale slack state (dt/dτ = cScale·κ, κ=‖r‖^p, Betts sparse free-time trick, τ_f fixed), cone-eliminated [α;s], inertial 2-body dynamics, exact-AD Hessian + IPOPT. Modes: `'fixedtf'` (pins t(τ_f)=t_f, ε energy→fuel objective) and `'mintime'` (s≡1, min t(τ_f)). Terminal pluggable (fixed rvf / free-L manifold). **Min-time runs on the free-L manifold from M0 on** — the paper's TfMin is free-longitude and every c_tf scale hangs off it; the fixed-terminal-first de-risking applies to the energy→fuel stage | params, terminal, dynamics |
| `seed_2body.m` | dynamically-exact tangential-thrust warm start: propagate α=v̂ at constant throttle s̄ (s̄=1 → min-time cold seed; s̄≈1/c_tf **bisected on arrival longitude** → fuel-stage seed in the right rev-topology — a stretched min-time seed lands ~3 revs short of the law-prescribed L_f), dense-output sampled at uniform Sundman-τ nodes (defect-free by construction) | dynamics, params |
| `homotopy_2body.m` | ε:1→0 sweep with tight re-clean + guard (the `sundman_homotopy` pattern) | solver core |
| `run_transfer.m` | single-case driver: mintime → seed → homotopy → verify → export | all above |
| `run_ctf_sweep.m` | loop c_tf (and T_max) → m_f-vs-c_tf front; check law C | run_transfer |
| `verify_pmp_2body.m` | PMP checks (primer/switching/transversality) on the 2-body Hamiltonian | solver output |
| `transfer_movie.m` | trajectory + throttle animation (reuse `psr_movie` pattern) | run_transfer |

We do **not** edit `casadi_minfuel_sundman` — `casadi_lt_2body` is a clean sibling
that shares the *architecture*, not the file. (A small `lt2b_rhs_time.m` dynamics
helper is shared by the solver, the seed propagation, and the tests.)

---

## 4. Solve strategy (mirrors our campaign)

1. **min-time anchor (free-L manifold):** solve s≡1, free t_f, terminal = insertion
   manifold → t_{f,min}(T_max) and the min-time longitude span ΔL_mt. Cross-check
   empirical law **T_max·t_{f,min} ≈ C** (≈ 850 N·h). L_f for the fixed-terminal fuel
   stage then comes from the paper's law: L_f = L₀ + c_Lf·ΔL_mt, c_Lf ≈ 1.12·c_tf+0.09.
2. **energy seed:** at t_f = c_tf·t_{f,min}, build a tangential constant-throttle
   (s̄≈1/c_tf) propagation seed, **bisected on s̄** so its arrival longitude lands within
   ~π of the prescribed L_f (right rev-topology), sampled at uniform Sundman-τ nodes
   (`seed_2body`), then solve ε=1 (energy) → machine-tight seed.
3. **sharpen:** ε:1→0 with tight re-clean per step → bang-bang min-fuel; read m_f, revs,
   switches from the flown solution.
4. **verify:** PMP first-order (primer/switching/transversality) + structure vs the paper.

---

## 5. Milestones & success criteria

Numbers below are the paper's (Table 3, Figs 21–23), at c_tf=1.5 unless noted.

- **M0 — coplanar, fixed terminal, 10 N.** Coplanar data (i=0) run through the 3D
  machinery — z stays quiescent, no separate 2D code — fuel-stage terminal = fixed rvf.
  *Success:* homotopy ε→0 converges machine-tight (defect ≤1e-8), throttle goes
  bang-bang, sane structure. (Smallest solve; validates dynamics+Sundman+homotopy.
  No paper number — internal.)
- **M1 — full 3D (7° incl.), fixed terminal at the law-prescribed L_f (≈ L₀+47 rad,
  see §2), 10 N.** *Success:* 7.5 revolutions hold by construction; the genuine
  tests are ~**18 switches** and **burns concentrated near apogee / coasts near
  perigee** (Fig 22; last-perigee corrective burns allowed).
- **M2 — free-L_f insertion manifold, 10 N, c_tf=1.5.** *Success (headline match):*
  **m_f ≈ 1370–1375 kg** (Fig 23), 7.5 revs, ~18 switches, apogee-burn structure
  (Figs 21/22). Requires the pinned Isp.
- **M3 — modest front.** m_f-vs-c_tf curve at 10 N (Fig 23 shape: m_f ≈ 1350→1388 kg
  over c_tf 1.05→3), and the empirical law **T_max·t_{f,min}≈C** across {10, 5, 2.5} N.
  *Comfortable:* 10 N (7.5 rev) → 2.5 N (30 rev). *Stretch:* 1 N (74.5 rev, ~2× our
  demonstrated Cartesian capacity — attempt, report honestly if the mesh wall bites).

**Out of scope (Path 1):** T_max ≤ 0.5 N (149–754 revs). These need an equinoctial
formulation — noted as the future escalation, not attempted here.

---

## 6. Verification & testing

- **Paper-match:** m_f, rev count, switch count vs Table 3; apogee-burn geometry vs
  Fig 22; law C.
- **PMP first-order** (`verify_pmp_2body`): primer alignment (thrust ∥ −p_v), switching-
  function sign law (burn where ψ<0, coast where ψ>0), transversality (p_m(t_f)=0 free
  mass; in M2 also λ·(∂x/∂L)=0 along the free phase direction of the insertion
  manifold — the Cartesian equivalent of the paper's p_L(t_f)=0). Reuses our
  KKT-dual approach.
- **Unit tests:** the paper's toy example **P₂** (double-integrator min-fuel, min∫|u|dt,
  x(0)=0→x(2)=(0.5,0)) as a cheap known-answer check of the energy→fuel homotopy;
  plus standard defect/edge/`‖α‖=1` guards. A no-solve smoke test for path/params.

---

## 7. Open items to resolve during implementation

1. **Isp** — ~~pin from ref [6] Caillau & Noailles 2001 (same benchmark); default 2000 s;
   validate by the M2 m_f match. (Only M2 depends on it.)~~ **CLOSED 2026-07-19.**
   Ref [6] obtained and read (`orbit_transfer/min_fuel_papers/COCV_2001__6__239_0.pdf`, p.255): the
   benchmark's mass-flow coefficient is δ = 0.05112 km⁻¹·s in ṁ = −δ·|thrust|, and
   δ = 1/(Isp·g₀) ⇒ c = 1/δ = 19.562 km/s ⇒ **Isp = c/g₀ = 1994.8 s** (Caillau &
   Noailles' exact value). Our default 2000 s was 0.27% high (masses ~0.3 kg high;
   validated by the M2 match, as planned). All other benchmark constants confirmed
   identical (P⁰=11625 km, Pᶠ=42165 km, e⁰=0.75, L⁰=π, m⁰=1500 kg, μ⁰=398600.47).
   Default kept at 2000 s with the exact value documented in `kepler_lt_params.m`
   (pass ispS=1994.8 for the source value).
2. **Canonical units** — choose LU (candidate: initial P or GEO radius), TU=√(LU³/μ),
   T_max → nondim acceleration. Fix in `kepler_lt_params`.
3. **Free-L_f local minima** (Fig 18 non-monotone at 10 N) — the paper resolves L_f by a
   min-longitude sub-anchor + c_Lf sweep (law R2: c_Lf,opt≈1.12·c_tf+0.09). M2 may need
   the same; M3's clean front uses the optimal-L_f (rendezvous) branch (Fig 23, smooth).
4. **Three-way core review (2026-07-17)** validated the Cartesian baseline as clean
   (1376.74 kg stands as the MEE validation gate) and produced two code findings now
   tracked in the follow-on specs: the primer/dual-map correction (`DESIGN_dual_map.md`)
   and the free-longitude manifold's missing prograde guard (`DESIGN_thrust_ladder.md`
   Phase 0). **[ROBUSTNESS]** nit also raised: `casadi_lt_2body.m`'s post-solve numeric
   defect re-check uses `norm(r)^pSund` while the NLP itself uses the softened
   `(r'r+1e-12)^(pSund/2)` — negligible at our radii, but a verifier should certify the
   actual NLP expression, not a slightly different one. See `doc/reviews/2026-07-17_triage.md`.

---

## 8. References

- Haberkorn, Martinon, Gergaud, JGCD 27(6), 2004 (the target paper).
- ref [6] Caillau & Noailles, *Coplanar control of a satellite around the Earth*, ESAIM
  COCV 6, 2001 (benchmark constants incl. Isp).
- Our machinery: `NLP_lowThrust_GTO_tulip/sundman_minfuel/` (Sundman min-fuel engine),
  `PSR/` (energy→fuel homotopy + PMP verification), `LOW_THRUST_MINFUEL_CAMPAIGN.md`.
