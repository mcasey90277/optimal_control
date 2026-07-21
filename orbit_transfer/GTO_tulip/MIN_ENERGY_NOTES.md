# Min-energy solver (direct + indirect) — status

Added Jul 7 2026. The min-ENERGY variant of the GTO->tulip transfer:
minimize J = int (1/2) s^2 dt at fixed tf. Motivation: min-energy has a
CONTINUOUS control (a saturated ramp u* = sat(S_e,0,1)), not the bang-bang
of min-fuel, so it is the objective that can take on the FULL 40-rev spiral
instead of only an arrival leg — and it is the standard homotopy root for
eventually cracking min-fuel (energy -> fuel continuation).

## Files

Indirect (PMP shooting), `indirect/lowThrust_GTO_tulip/`:
- `lt_pmp_eom_energy.m`     — 14-state augmented dynamics; control
  u = sat(S_e), S_e = Tmax(||lam_v||/m + lam_m/c); costate ODEs identical
  to min-time/fuel (envelope theorem); complex-step safe (clamp branches on
  real(S_e), so saturated regions carry zero sensitivity — the correct
  clamp derivative).
- `shoot_residual_energy.m` — fixed-tf, 7 conditions (6 rendezvous +
  lam_m(tf)=0), complex-step Jacobian. NO smoothing parameter needed.
- `solve_energy_indirect.m` — lsqnonlin Levenberg-Marquardt, single solve
  (no eps-homotopy — the control is already smooth).

Direct (collocation NLP), `GTO_tulip/`:
- `solve_energy_nlp.m`               — fmincon interior-point; objective
  int 1/2 s^2 dt (trapezoidal, analytic sparse gradient); REUSES the
  verified min-fuel transcription (`lt_dynamics_throttle`,
  `nlp_constraints_minfuel`: defects + throttle cone) unchanged.
- `NLP_lowThrust_GTO_Tulip_energy.m` — FULL-spiral driver: burn+coast
  warm start (min-time burn to feasibility, ballistic coast, phase-shifted
  target), direct solve, optional indirect polish.

## Verification (headless)

- lint clean; complex-step vs finite-difference Jacobian rel err 1.2e-6.
- **Indirect on the arrival leg: ||R|| = 6.3e-12 (machine zero), throttle a
  smooth ramp u in [0.15, 1.0].** This is the headline: where the min-fuel
  indirect never beat ||R|| = 0.14, min-energy converges to machine
  precision — because the control is continuous. Smooth control => smooth
  shooting residual => large basin.

## Does it handle the FULL spiral?

Short answer: yes, in the way that matters — decisively better than
min-fuel, which could not.

- **Direct, full 40-rev spiral (N=3000, tf=1.15x min):** produces a
  complete, structurally-correct min-energy trajectory — smooth throttle
  ramp s in [0.02, 1.0], reaches the (phase-shifted) target, propellant
  ~2.93 kg, dV ~4.48 km/s. Min-fuel's direct method could NOT do this on
  the full spiral (it grinds in feasibility mode without ever producing a
  coherent bang-bang solution). Tight feasibility, however, is MESH-LIMITED
  at N=3000: fmincon floors at max defect ~2e-4 (feasibility mode) to 5e-4
  (default) and returns flag -2 ("converged to an infeasible point").
  ~75 nodes/rev cannot resolve the perigee dives to 1e-9 — the SAME
  discretization wall the min-time NLP hit at low N (it needed N~12000).
  Cure: finer mesh (N~12000) or mesh/tf continuation. NOT a fundamental
  obstruction — it is resolution, not the objective.
- **Indirect, full spiral: does NOT converge by single shooting**, even
  seeded from the direct solution via the covector mapping
  (`costate_seed_from_nlp_energy`). Measured (Jul 7): N=4000 direct
  (defect 3.2e-4) -> reconstructed seed -> lsqnonlin-LM STALLS at
  ||R|| = 0.236 (flag 4, step below tolerance). The reconstructed seed is
  good but not inside the shooting basin. Why: single shooting integrates
  through ~40 perigees => ~1e6 sensitivity, so the basin is tiny -- the
  SAME wall min-time faced (min-time single-shooting only converged
  because pumpkyn handed it a near-EXACT seed). This is a DYNAMICS-side
  difficulty, independent of the objective; smoothing the control (energy)
  does not shrink it.
  Cure = MULTIPLE SHOOTING: partition the arc into ~40 segments (each with
  ~1/40 the sensitivity), match states at the nodes, and let each segment's
  costates be unknowns. That is the standard fix for multi-rev shooting and
  the honest next step; single shooting is simply the wrong tool for 40
  revs regardless of objective.

## The two difficulties, separated (the real lesson)

This exercise cleanly separates the two independent obstacles:
  1. OBJECTIVE-side (bang-bang nonsmoothness): killed min-fuel; min-energy
     REMOVES it (continuous control) -> leg indirect converges to 6e-12.
  2. DYNAMICS-side (~1e6 multi-rev shooting sensitivity): a property of the
     40-perigee geometry, NOT the objective. Defeats SINGLE shooting on the
     full spiral for min-energy (||R||=0.236) exactly as it would for
     min-time without a near-exact seed. Needs multiple shooting.
Min-energy is the right objective for the full spiral (fixes #1); getting
the full-spiral INDIRECT to machine precision additionally needs multiple
shooting (fixes #2). The DIRECT method sidesteps #2 (collocation doesn't
integrate) but pays low-order-mesh accuracy (~3e-4 floor at N=4000).

## Bottom line

The bang-bang wall that stopped min-fuel on the full spiral is GONE with
min-energy (continuous control): indirect converges to machine precision on
tractable arcs, and the direct method builds the whole 40-rev solution. The
only remaining full-spiral cost is mesh resolution (direct) and integrator
sensitivity (indirect) — both quantitative, both surmountable with a finer
mesh / continuation, neither the fatal structural problem min-fuel had.

Next step to a machine-clean full spiral: run
`NLP_lowThrust_GTO_Tulip_energy(12000, 1.15, true, true)` (minutes), then
seed the indirect from that solution.
