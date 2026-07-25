# earth_elliptic_to_geo_CR3BP — TODO

**Phase 1 COMPLETE (2026-07-24):** full certified ladder 10→0.1 N + 10 N φ₀
sweep + note + movies. What remains: deep-rung φ₀ sweep, stronger optimality
certificates, and the indirect Phase 2.

## Phase 0 — formulation decisions (design before code)

- [x] **Dynamics representation.** Decided: (a) MEE/L-domain solver from
  `../earth_elliptic_to_geo/` + lunar perturbation acceleration in the Gauss
  equations (D1, keeps the winning ΔL-free formulation and mesh behavior).
  Option (b) rotating-frame Cartesian rejected — see spec sec 2 D1.
- [x] **Terminal set in the chosen frame.** Decided: equatorial GEO stays the
  trivial MEE target ([1;0;0;0;0]); t_f = c_tf·t_f,min convention carries
  over unchanged, anchored to the 2-body t_f,min (spec sec 2 D4).
- [x] **Params home.** Decided: own `lunar_params.m` file in `direct/`
  (spec sec 2 D2/D3), not a generalized `cr3bp_common/cr3bp_lt_params`.
- [x] **Moon model sanity bound.** Built: `direct/sanity_bound.m` +
  `direct/results/sanity_bound.md` — tide/authority ratio 0.11% at 10 N
  growing to 10.9% at 0.1 N (spec sec 7 null model).

## Phase 1 — direct

- [x] Gravity-homotopy bridge at 10 N: warm-started from the certified
  2-body 10 N solution, Moon mass dialed 0 → μ* (`bridge_mu_continuation.m`,
  T4), then energy → fuel eps-sharpened to a CERTIFIED CR3BP min-fuel
  solution (`solve_cr3bp_minfuel.m`, T5;
  `direct/results/minfuel_cr3bp_T10N_phi0.mat`).
- [x] Compare vs 2-body: `direct/compare_vs_2body.m` (T6) —
  Δm_f = +0.0545 kg (+0.00396%), 19/19 switches (nodal, mesh-band caveat),
  maxDefect = 4.19e-15. Full table:
  `direct/results/compare_vs_2body.md`. (The analytic
  Bonnard–Caillau–Picot rev-count/thrust bound cross-check is still open —
  folded into the deep-rung walk below.)
- [x] **Complete thrust ladder 10→0.1 N** (2026-07-24, commit `f3692b0`):
  all 7 rungs certified via the front door + same-chain gain=0 controls;
  Δm_f(T) = +52.0…+31.3 g. The predicted deep-rung ε-wall never appeared —
  `liftDL=true` + `maxIter≥4000` carried 0.2 N (890 sw, 265 d) and 0.1 N
  (1724 sw, 531 d) straight through (see the MUMPS/liftDL lesson, note §
  thrust ladder). Interestingly the Moon effect did NOT grow toward the ~11%
  sanity-bound ratio — it stays ~50 g then gently *declines*, the
  phase-averaging picture.
- [x] **φ₀ sweep at 10 N** (2026-07-24, commit `f3692b0` + note §phase,
  fig `fig_phi_sweep.m`): +52.0/−56.3/+68.9/−52.8 g at 0/π/2/π/3π/2 —
  π-periodic tidal quadrupole, the sign FLIPS with phase, so "the Moon
  helps" is phase-specific not generic. Control-law sensitivity quantified
  (`phase_control_sensitivity.m`): ≤1.3 min switch shifts, envelope
  argument. Movies: `phase_quad_movie.m` (4-panel synced).
- [ ] **Deep-rung φ₀ sweep** — the direct test of the phase-averaging
  interpretation: at 0.2/0.1 N (t_f ≫ lunar month) the quadrupole amplitude
  should collapse, leaving only the ~30 g secular floor. Recorded in the
  note as the testable corollary.
- [ ] **Switch-count band marker in `compare_vs_2body.m`** — now that
  0.2/0.1 N are in the table (890/1724 nodal switches), the swStr column
  should carry an explicit band marker, not bare "N/N" (spec sec 8 gate 4;
  final-review-report.md gate-4 finding).
- [x] **CR3BP-aware PMP/primer verification** (task B, 2026-07-23,
  `feat(verify): lunar-aware PMP verification + CR3BP campaign driver`):
  `mee_primer_switch.m` now subtracts the zero-throttle ballistic/lunar
  bracket out of its B(X)/pel extraction and its S-formula G0 term before
  forming the primer vector and switching function (opt-in via `par.pert`,
  byte-identical when absent); `verify_cr3bp_pmp.m` drives it on the
  front-door artifacts. Ran it on the certified 10 N and 5 N CR3BP fuel
  solutions: primer/sign GATE STILL FAILS (10 N: 32.4 deg / 78.4%; 5 N:
  26.7 deg / 76.6%) — but this is the SAME pre-existing, already-tracked
  eccentricity-correlated raw-`lam_g`/KKT-dual anomaly that already fails
  the PURE 2-body certified solutions identically (`../earth_elliptic_to_geo/
  TODO.md` sec "eccentricity-correlated" item, `.superpowers/sdd/
  task-10-report.md`) — NOT a new lunar-specific defect. The lunar
  coupling itself (`A0/Ldot0` diagnostic) is confirmed tiny (~1e-5) and its
  removal changes the gate numbers by <0.1 deg. Certification still rests
  on the four NLP metrics + bound-saturation check, not on primer
  agreement, pending the raw-`lam_g` fix landing campaign-wide.

## Phase 2 — indirect

- [ ] PMP shooting counterpart (costate dynamics gain the lunar-gradient
  terms), seeded from the direct solutions — same direct-seeded strategy as
  `../GTO_tulip/indirect/ifs/`.
- [ ] **Second-order certification: conjugate-point (Jacobi-field) check**
  along each converged extremal, verifying t_conj > t_f — upgrades "shooting
  converged" to "locally optimal" (Bonnard–Caillau–Picot 2010 §2.3–2.4 +
  their cotcot code ref [12]; the variational integration reuses the shooting
  Jacobian machinery). Our current indirect checks are all first-order; this
  closes that gap and should also back-port to `../GTO_tulip/indirect/ifs/`.

## Stronger optimality evidence (note sec 9-10)

- [x] **CR3BP-aware primer** half DONE (task B, 2026-07-23, commit
  `feat(verify): lunar-aware PMP verification + CR3BP campaign driver`):
  ballistic rate subtracted from the primer reconstruction, opt-in via
  `par.pert`, regression-clean on the 2-body suite; see the Phase-1 item
  above for the honest gate numbers (still FAIL, pre-existing dual anomaly,
  not lunar-specific).
- [ ] **PSR half remains open.** PMP-steered switch-time refinement (PSR)
  per the tulip pattern — needed before any published switch-time claim,
  decisive at deep rungs. (Blocked in practice on the same raw-`lam_g`
  dual anomaly above, since PSR consumes the same primer/switching
  reconstruction this task made lunar-aware but did not make agree with
  the solver's own switches.)
- [ ] **SOSC tier 1 (cheapest, machinery half-built):** finish the
  2-body campaign's PLAN_sosc reduced-Hessian check on casadi_lt_mee's
  returnModel registries; run it on the certified 10 N and 5 N solutions
  (2-body and CR3BP).
- [ ] **SOSC tier 2:** induced switching-time-problem Hessian
  (Maurer-style bang-bang sufficient conditions) on the arc-parameterized
  representation shared with PSR.
- [ ] **SOSC tier 3:** Jacobi-field conjugate-point test with the Phase-2
  indirect solver (BCP 2010 pattern).
