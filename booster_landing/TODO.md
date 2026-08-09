# booster_landing — TODO

Curated from the theory note's future-work section, the external code-review
parked items (rulings in `process/CAMPAIGN_LEDGER.md`), and post-campaign
discussion (2026-08-09). The campaign itself is COMPLETE and certified —
everything here is extension or hardening, nothing blocks current use.

## P0 — COMMITTED ROADMAP (user-adjudicated 2026-08-09, in this order)

Three campaigns Mike has committed to, each with the full spec → plan →
certify discipline of the base campaign. Order chosen deliberately:
controller upgrades first on the certified 3-DOF plant (algorithm bugs
isolated from plant bugs, near-total machinery reuse), then the heavy
plant campaign ports two proven controllers.

- [ ] **R1 — On-the-fly relinearization for K.** Evaluate the (existing,
      complex-step-verified) Jacobians about the CURRENT state each cycle
      instead of the stored nominal, extending the tracker's basin of
      validity beyond the small-deviation regime. Smallest campaign;
      changes where A(t),B(t) are evaluated + a gain-update strategy
      (recompute cadence, Riccati re-sweep vs SDRE-style algebraic solve).
- [ ] **R2 — Full MPC (onboard replanning).** Re-solve the guidance from
      the current state receding-horizon — `solve_pdg_convex` warm-started
      at each replan is the natural engine (sub-second convex solves; the
      G-FOLD onboard story). Directly addresses the measured
      fixed-reference-tracker limit (a -5% engine needs the braking point
      re-timed ~90 m higher). Decisions at spec time: replan cadence,
      solve-latency modeling, fallback to TVLQR between/on failed solves.
- [ ] **R3 — 6-DOF rigid body.** Attitude on SO(3)/quaternions (lieFiltering
      home turf), inertia, engine gimbal + grid fins, aero moments,
      terminal-verticality constraint (fixes the tilted landing
      physically). Both controllers port here. Pointing-cone machinery
      (`P.theta_max_deg`) and G5's cone-aware primer semantics are ready
      but never exercised end-to-end — add a cone-active test when this
      starts.

## P1 — natural next builds (unordered, after or alongside the roadmap)

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
