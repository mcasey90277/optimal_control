# earth_elliptic_to_geo_CR3BP — TODO

Nothing built yet. Phases in intended order:

## Phase 0 — formulation decisions (design before code)

- [ ] **Dynamics representation.** Choose: (a) MEE/L-domain solver from
  `../earth_elliptic_to_geo/` + lunar perturbation acceleration in the Gauss
  equations (keeps the winning ΔL-free formulation and mesh behavior), vs
  (b) Earth–Moon rotating-frame Cartesian (the `../GTO_tulip/` engine's
  native form, Sundman-regularized). Decide with a short spec, not in code.
- [ ] **Terminal set in the chosen frame.** Equatorial GEO: trivial in MEE
  ([1;0;0;0;0]); a time-dependent circle in the rotating frame. Also decide
  whether t_f conventions (t_f = c_tf·t_f,min) carry over unchanged.
- [ ] **Params home.** New craft-specific CR3BP params (1500 kg, 10..0.1 N,
  Isp ~2000 s, Earth–Moon μ*) — own file here, or generalize
  `../cr3bp_common/cr3bp_lt_params` to take the craft as input.
- [ ] **Moon model sanity bound.** Back-of-envelope first: lunar perturbing
  acceleration along the spiral vs thrust accel (10 N/1500 kg ≈ 6.7e-3 m/s²
  down to 0.1 N ≈ 6.7e-5 m/s²) — predicts where in the ladder the Moon
  starts to matter, and gives the null hypothesis the solves must beat.

## Phase 1 — direct

- [ ] Gravity-homotopy bridge at 10 N: warm-start from the certified 2-body
  10 N solution, dial Moon mass 0 → μ* (reuse the
  `gen_elfo_energy_gravhom` two-primary ladder pattern), then energy → fuel
  to a certified CR3BP min-fuel solution.
- [ ] Compare vs 2-body: Δm_f, switch structure (as a mesh-band), R0-law drift
  — and against the ANALYTIC rev-count/thrust bound of Bonnard–Caillau–Picot
  2010 Prop 3.4/3.7 (σ²(n,e) dominating the control quadratic form: an
  l_f-vs-thrust-bound estimate, the R0 law's rigorous cousin).
- [ ] Walk the thrust ladder down while it stays interesting (the Moon effect
  should grow as thrust drops and transfer time stretches).

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

## Housekeeping

- [ ] Create `direct/`/`indirect/` when the first code lands; keep this
  README/TODO current.
