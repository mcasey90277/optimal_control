# Astra review #2 (2026-09-06): the Huber line, code + theory -- adjudicated

Two raw-API GPT-6 Astra passes (`huber_code_astra_2026-09-06.md`, 85 KB
bundle, 179 s, $0.73; `huber_theory_astra_2026-09-06.md`, 31 KB, 155 s,
$0.45) on everything built since the 09-05 review: the corrected
conjugate instrument, the event-split Huber propagator, the huberc
family, the (p, delta) race, the switch diagnostic, the catalog binding,
and the five theory claims of FINDINGS 24-26 and the certification
survey. Host verified the checkable claims in MATLAB before accepting.

## Verified by measurement (accept as fact)

- **Q-dot is continuous across the PMP switch.** Astra: for this
  Hamiltonian the jump in F is along the symplectic gradient of Q, so
  `n'(F+ - F-) = 0` and `Qdot = -T lam_v'lam_r / (m rho)` exactly.
  Measured at a Huber event state: `n'F- = n'F+ = 6.3986117192`,
  difference 8.9e-16; closed form matches to 1e-9. Consequences: the
  grazing assert in `propHuber` compares a number to itself (its scale is
  vacuous); `huber_switch_diag` should use the closed form, not
  `gradient(Q,T)`; the saltation numerator structure is confirmed.
- **Never-full-rank input returns PASS.** `kFull = nS+1` leaves `live`
  empty, `nIn = 0`, `pass = true`. Confirmed on a fixture where the state
  row never sees the costate. Must be UNDETERMINED.
- **Root exactly at t_f is roundoff-signed.** LQ with T = pi: det at t_f
  is -1.5e-16; whether the last bracket "crosses" is decided by noise. A
  genuine endpoint root (weak, non-strict minimum) must not become a
  refutation by luck.

## Accepted on argument (no measurement needed)

Code:
- `dets(k) = sign(det(Me))*abs(dets(k))` feeds a possibly-underflowed raw
  magnitude back into the zero test; keep sign and log-magnitude apart.
- `[+,0,-]` counts two crossings; runs of zeros count once per sample;
  merge zeros with adjacent brackets.
- `propHuber` at Q(y0) == 1 always picks 'lo'; choose by sign of the
  closed-form Qdot. `dt < 0` is broken by `te(end) < dt - 1e-14`: reject
  it explicitly. A switch exactly at a shooting boundary gets the
  incoming STM silently: flag it.
- huberc: validate `delta > 0` finite; explicit `p = 1` dispatch instead
  of the `1e-300` trick; `aux.sRaw` reports the core law only (rename or
  extend to the piecewise inverse law).
- `run_minfuel_race`: stage 2 starts with `xGood = []` (first stage-2
  failure has no bisection anchor); `retired` reset only because
  `deltaSched` is nonempty; loose-floor fallback trims every per-rung
  array EXCEPT `.delta` and leaves `.z`/`seed` at the discarded rung, so
  the acceptance gate re-solves from a stale seed and binds its verdict
  to a different solution; after a successful refinement `coastFrac/
  iters/wall` describe the loose solve. None of the reported tables came
  through the loose/fallback path except the 09-06 wall-diag loose runs,
  whose conclusion (loose gate net worse) does not change.
- `run_conj_fixedtf_sweep`: `.mf = it.Y(7,end)` is the last junction
  START mass, not m_f (use `Yend`); verdicts recorded regardless of
  `converged` (all 15 were converged, so no live effect); binding should
  also check `p` and the spec.
- `huber_switch_diag`: `[-,0,-]` becomes two fictitious crossings after
  `s(s==0)=1`; `tCross` is the left sample, not the root; extremum
  "within 10% of 1" is a screening window, not a bifurcation verdict.
- Test fixtures: the coast fixture `Z = [0 0; 0 1]` is not a fundamental
  matrix (use identity: frozen subsystem); add underflow, even-
  multiplicity, endpoint-root and no-testable-sample cases.

Theory:
- **Claim 2 was undersold.** For a CONTINUOUS clipped law (eps, huberc)
  `F+ = F-` at a transversal branch boundary, so the saltation matrix is
  the identity and the branchwise-AD STM IS the derivative of the full
  flow -- not of an active-set-frozen surrogate. Our case-B Jacobi test
  is therefore the ordinary Jacobi necessary condition under regular
  crossings, stronger than the survey said. A separate switching-time
  Hessian is not automatically required for continuous saturation.
- **Claim 3 overstated in two places.** (i) `|Qdot| = 0.045` and
  `Q_max = 0.9907` are grazing RISK INDICATORS; a grazing bifurcation
  needs `h = 0, hdot = 0, hddot != 0` on the branch plus a nonzero
  unfolding in p. (ii) "No parametrization passes a corner" is false:
  the local behaviour is typically square-root (`t_switch ~ mu^{-1/2}`
  sensitivity, endpoint correction ~ sqrt(mu)), which a reparametrization
  `w = sqrt(mu)`, explicit switching times, or a structure-change step
  can traverse. The correct statement: the current smooth Newton on p is
  unsuitable near a structure change. "Flat cond(J) excludes a fold" is
  weaker evidence than stated (finite-distance samples at nonzero
  residuals). The family change (huberc) is a valid cure, not the only
  one -- and it worked, which is the empirical fact that stands.
- **Claim 4 over-generalized.** Three cells on one solver configuration
  give a configuration-dependent empirical floor, not a family-
  independent one; ramp width in Q is not a temporal resolution
  (`dt_ramp ~ delta/|Qdot|`); segment count alone does not set switch
  resolution; and `delta -> 0` at fixed `p = 0.001` recovers the Huber
  JUMP, not exact bang-bang fuel. m_f convergence to 1e-6 is one scalar;
  switch times may still differ. The switching-time endgame is
  sensible; its convergence certifies stationarity, not minimality, and
  the direction control alpha(t) on burns must be in the accessory
  problem.
- **Claim 1 needs its hypotheses audited, not gates added.** The BCT
  sufficiency theorem must be applied with our actual endpoint manifold
  (mass free; for all-burn min-time mass is eliminable, `m = m0 - Tt/c`,
  which is why rows 1:6 / cols 8:13 is the right block); "no conjugate
  time at integrator resolution" is a numerical assessment, not a proof
  (needs between-sample bounds, even-multiplicity handling, endpoint
  focal degeneracy); normality: tfMin acceptance gives a normal LIFT, it
  does not exclude an abnormal lift of the same trajectory; the costate-
  scaling quotient is an invariance of the CONTROL (so `Phi_xl*lam0 = 0`
  exactly), but with lambda_0 = 1 and H = 0 it is not a symmetry of the
  BVP -- the equivalence to the theorem's construction must be shown.
  The catalog's "min-time" is the ALL-BURN problem; optimality against
  throttle variations is assumed, not tested.
- **Claim 5 correct** for 0 < p < 1, delta > 0 on s in [0,1]; treat the
  saturated branch as `L + indicator_[0,1]` (the constant extension is
  not convex at 1); explicit p = 1 and p -> 0 limits; at delta = 0 the
  argmin at Q = 1 is the whole interval [p, 1] and the jump law is a
  selection.
- **Register M1 disputed.** "Fuel linearity makes every NLP Hessian
  incapable of detecting strictness" is not a valid impossibility
  argument: the Lagrangian's reduced Hessian on the kernel of the
  endpoint-constraint Jacobian need not vanish because `L(s) = s`. The
  measured 270 flat directions at earth 10 N stand as a measurement; the
  theoretical ceiling claimed from it does not. Recorded for the
  register; not this line's problem.

## Rejected / out of scope

- "`maxBisect` is per-gap" -- known P2 from the first review, unchanged.
- "Coast fixture is not a fundamental matrix" -- true and harmless for a
  synthetic indexing test; will switch to identity anyway.
- `seedScale` on all junction columns is CORRECT for a multiple-shooting
  seed (Astra concurs); making it family-aware (eps needs 1) is a
  convenience.
- Higher m_f as the rewalk tie-break is a legitimate selection
  preference, not a minimality claim (Astra concurs); the catalog note
  already says "no detected conjugate point", not "certified".

## What this changes in the record

FINDINGS 24: "grazing bifurcation" -> "grazing-risk signature; a
structure change in the switch set"; delete "no parametrization of p
passes it"; keep "the smooth p-Newton stalls there and the continuous
family passes it" (measured). FINDINGS 26: "family-independent floor" ->
"configuration-dependent empirical floor on three cells; factors to vary
listed"; add "delta -> 0 at fixed p is the Huber jump, not bang-bang".
Survey 2.1 case B: branchwise-AD STM is the true flow derivative under
regular crossings; case A: the hypotheses-audit list above replaces
"two cheap additions"; note M1 disputed.

## Fix plan (proposed, not yet executed)

P0 (instrument correctness, re-sweep after):
  1. `ms_conjugate_test`: UNDETERMINED when no sample is testable
     (`pass = false`, `.tested = false`); sign from the equilibrated LU,
     magnitude via log-det, never fed back; merge exact zeros with
     adjacent brackets; endpoint-only crossing -> `.endpointOnly`, not
     a refutation; `.sampledThrough` when `Yend` is absent.
  2. `propHuber`: branch at Q == 1 from closed-form Qdot; grazing guard
     on `|Qdot|` against an error scale, state-only path included;
     assert `dt > 0`; flag switch-at-boundary.
  3. `huberc`: validate delta; explicit p = 1 dispatch; fix `sRaw`.
  4. `run_conj_fixedtf_sweep`: `mf` from `Yend`; verdict only if
     converged; bind p + spec in the catalog builder.
P1 (race bookkeeping; no reported number changes): stage-2 `xGood`,
  `retired`, full-array trimming incl. `.delta` and `.z`, acceptance
  seed rebuilt from the retained rung, metrics refreshed after
  refinement, diagnostic propagations inside try/catch.
P1 (diagnostic): `huber_switch_diag` on closed-form Qdot and the
  propagator's event records; touch vs crossing classification.
Docs: the three record edits above.
Tests: identity coast fixture; underflow / even-multiplicity /
  endpoint-root / no-testable cases; huberc knees and limits;
  ascending + descending + boundary switches for the propagator.
Then: re-run the 15-verdict sweep and golden cells on the fixed
instrument (expected unchanged, must be shown).
