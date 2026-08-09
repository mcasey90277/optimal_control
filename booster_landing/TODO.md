# booster_landing — TODO

Curated from the theory note's future-work section, the external code-review
parked items (rulings in `process/CAMPAIGN_LEDGER.md`), and post-campaign
discussion (2026-08-09). The campaign itself is COMPLETE and certified —
everything here is extension or hardening, nothing blocks current use.

## P1 — natural next builds (each is a real campaign)

- [ ] **State estimation in the loop.** The tracker currently receives the
      true state (perfect navigation). Feed it an estimated state from an
      INS/GNSS filter instead — the `~/Desktop/navigation` repo's GPS-aided
      error-state EKF is the natural block — and add sensor-grade
      dispersions to the Monte Carlo. Closes the loop in the full GNC sense
      and marries two existing projects.
- [ ] **Engine-off coast + ignition-time optimization ("Phase 3").** The
      landing burn matches how a real hoverslam flies; what's missing is
      the free-fall before ignition. No hybrid OCP needed: wrap the
      existing (unchanged, still-convex) solver in a 1-D outer search that
      slides the ignition point along an incoming ballistic arc and picks
      the cheapest. Recovers most of the Tmin-arc propellant (~750 kg
      class) at near-zero formulation cost.
- [ ] **6-DOF rigid body.** Engine gimbal + grid fins, attitude dynamics,
      aerodynamic moments — and a terminal-verticality constraint, which is
      what makes a legged rocket land straight (the 3-DOF optimum genuinely
      lands ~20 deg tilted; the cinematic movie eases it cosmetically).
      Pointing-cone machinery (`P.theta_max_deg`) is already in place and
      G5's cone-aware primer semantics are ready, but never exercised
      end-to-end — add a cone-active test when this starts.
- [ ] **Onboard replanning (MPC / re-solved guidance).** The measured limit
      of the fixed-reference tracker (a -5% engine needs the braking point
      re-timed ~90 m higher; feedback around a stored plan cannot do that)
      is the motivation. This is what the flight-proven G-FOLD lineage does
      onboard.
- [ ] **SCvx for the drag phase.** Restores the two-solver cross-validation
      (G3/G4, currently 'skipped' under drag because the convex route is
      vacuum-only). Successive convexification re-linearizing about the
      converged trajectory.
- [ ] **Full return profile.** Boostback burn, entry burn, aero descent —
      each its own phase; the landing burn slots in as the final leg.

## P2 — hardening (parked external-review items, ruled 2026-08-09)

- [ ] G1 per-row scaled tolerance (retire the drag-path `tol=1e-11`
      workaround; the mass row is effectively gated ~30000x tighter than
      position in the current absolute-SI max).
- [ ] `tvlqr_design`: validate `tBlend > 0` and R symmetric positive
      definite; replace `inv(R)` with a factorization (GPT#12).
- [ ] `annulus_switch`: reject or properly schedule unsupported throttle
      structures (first-upward-crossing picks the wrong switch on a
      max-min-max profile — live the moment drag/cone re-solves produce
      one; TODO comment in-file) (GPT#15).
- [ ] `run_monte_carlo`: validate `Nrun`/`sig`/`P.seed` finiteness; use an
      isolated RandStream instead of mutating the global one (GPT#13).
- [ ] `hs_quad_ctrl`: the bound-touch step condition can still misfire on
      an arc that tangentially touches a bound (documented limitation);
      tie switch detection to a certified discrete switching function.
- [ ] Truncate/bound the Gaussian engine dispersions (negative draws are
      ~66-sigma non-events at current sigmas, but the sampler allows them).
- [ ] Golden-section validity-gap mapping (a rejected-probe region inside
      the bracket could steer it; map feasible intervals first).
- [ ] `doc/figs/` drift check: the note's committed PNGs are hand-copied
      from `results/` — add a checksum/regeneration script so a re-run with
      different numbers can't silently diverge from the figures.

## P3 — polish / nice-to-have

- [ ] Title card: label the Monte Carlo line "vacuum" (Phase 2 wind MC is
      99.0% and lives in the note/README).
- [ ] Ghost-fleet Monte Carlo scene for the cinematic movie (20 translucent
      dispersed trajectories converging on the pad) — the "99.5% as one
      image" shot.
- [ ] Side-by-side vacuum-vs-drag race cut (Phase 2 story: drag lengthens
      the flight but saves 434.7 kg).
- [ ] `solve_pdg_colloc`: delete the dead `S.ve` field; harden the failure
      path (`opti.debug` dual/stats fetch can throw before packaging).
- [ ] Conic-solver drop-in (ECOS/CVX) for the convex route — replaces the
      convex-NLP-via-IPOPT argument with a certificate-grade solve.

## Done markers

Campaign record: `process/CAMPAIGN_LEDGER.md`. Certified baseline:
`README.md` (expected flagship results). Adjudication history: the design
spec (`docs/superpowers/specs/2026-08-08-booster-landing-design.md`) is the
tiebreaker document.
