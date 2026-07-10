# ms_band — indirect multiple shooting for the 1.01–1.11x transition band

Campaign record. Spec: `../../docs/superpowers/specs/2026-07-09-ms-band-design.md`.
Plan: `../../docs/superpowers/plans/2026-07-09-ms-band.md`.

## Status
- [ ] M1 machine validated (Jacobian + min-time + 1.12x reproduction)
- [ ] M2 anchors (1.01x up, 1.11x down)
- [ ] M3 band march
- [ ] M4 assembly (front figure + docs)

## Log
(append dated entries here as work proceeds)

**2026-07-09 — Task 5: LM wrapper + M1(a) min-time anchor; residual-floor fix.**
`ms_solve.m` (guarded LM via lsqnonlin, analytic sparse Jacobian) + `test_ms_reproduce_mintime.m`.
First run at RelTol 1e-12: physics exact (dV 4.4665 km/s, prop 2.9247 kg, bang 100.0%) but
||R|| stalled at 1.206e-9 vs the 1e-9 gate. Row-class diagnostic at the converged Z: costate
continuity rows dominate (block norm 9.5e-10; worst row = joint-1 lambda_ry defect 7.9e-10 on a
~2e2-magnitude costate = ~4e-12 relative) — i.e. the stall is the integrator's own evaluation-noise
floor (~costate scale x RelTol x sqrt(Nrows)), not a solver deficiency. Review-verified fix applied:
`ms_problem` odeOpts tightened to RelTol 1e-13 / AbsTol 1e-15 (gates unchanged, tolR still 1e-9).
Re-run: ||R|| = 2.253e-10, dV 4.4665, prop 2.9247, bang 100.0% — PASS test_ms_reproduce_mintime
(10 LM iters, ~1 min; seed residual is ~2e-5 at 1e-13 because the tighter integrator exposes the
min-time reference's own ~1e-6 gate; LM absorbs it in one step).
**Known regression (needs gate-owner decision):** `test_ms_jacobian` now FAILS Gate 1 on col 49
(arc-4 lambda_m): |J-FD| 1.086e-5 > tol 4.241e-6. The CS backstop (Gate 2) is green (h-independence
1.9e-13), and an h-sweep shows the FD error is noise-dominated ~1/h at the test's hFD=1e-7
(falls to 2.8e-9 at hFD=1e-4, confirming the CS column) — the h-vs-h/2 self-consistency tolerance
underestimates correlated integration noise at RelTol 1e-13. The Jacobian is correct; the FD gate
heuristic is miscalibrated. Candidate fixes: larger hFD (1e-5..1e-4) for Gate 1, or defer to the
green CS backstop when Gate 1 is noise-limited. Not changed here (gate discipline).
Logs: `test_mintime_reltol13.log`, `test_jacobian_reltol13.log`, diagnostics `task5_diag*.log`.

### 2026-07-10 — M1(b) BLOCKED, adjudicated; pivot to up-march-first
The 1.12x dual-seed reproduction failed to converge under every recipe
(time-domain MS): best case (M=48, eps=1, 600 LM iters) decelerates
asymptotically at ||R||~1e-2. Seed machinery itself VALIDATED (beta=0.03102,
spread 0.5%, burnAgree 100%, arcCheckErr 2.6e-3 at M=48) — committed
(beta_from_duals.m, seed_from_duals.m). Root-cause candidate: KKT duals are
SUNDMAN-domain costates; time-domain conversion carries a discrete
(dkappa/dr)L offset on burn arcs (largest at perigee burns), and the
unscaled mixed-magnitude system cripples LM. Full traces + analysis:
.superpowers/sdd/task-6-report.md. test_ms_reproduce_112.m left uncommitted
(fails; revisit at the endpoint cross-check).
Rulings (user): (1) proceed up-march-first (Task 7; clean integration seeds);
1.12x reproduction deferred to an endpoint cross-check where the up-march
meets the band edge; Sundman-domain MS rebuild is the fallback if the
up-march also crawls. (2) Guard amended: within-step LM relay allowed
(iteration-capped solve with monotonically decreasing resnorm may continue
from its own iterate at the SAME eps; never across steps; worse iterates
still discarded). Also adopting ScaleProblem='jacobian' in ms_solve
(probe: ~10x early-phase speedup).

### 2026-07-10 — Task 7 step 1: ms_solve ScaleProblem='jacobian' + regression
Added `'ScaleProblem', 'jacobian'` to the lsqnonlin options in `ms_solve.m`
(Marquardt scaling; probe-verified ~10x early-phase speedup on the crawling
1.12x dual-seed system — `probe_scaled_lm.m`/`.log`). Re-ran
`test_ms_reproduce_mintime` as a regression gate: PASS, ||R||=2.384e-10
(15 iters; comparable to the pre-change 2.253e-10 at 10 iters — this
easy near-converged problem doesn't exercise the scaling benefit the probe
measured on the hard dual-seed case, but the gate stays green). physics
exact (dV 4.4665 km/s, prop 2.9247 kg, bang 100.0%). Log:
`regress_mintime_scaled.log`.

### 2026-07-10 — Task 7: M2(a) 1.01x up-anchor BLOCKED (adjudicated); Sundman pivot triggered
Machinery landed and validated: `eps_march.m` (guarded schedule march +
within-step relay per the 2026-07-10 guard amendment — relay accounting
worked exactly as designed), `seed_from_mintime.m` (clean integration seed:
||R(seed)|| = 0.72 at eps=1, entirely the 7-row terminal miss; joint
defects ~0 by construction), `run_anchor_up.m` (M=24 with one-shot M=48
escalation).
**Result: the eps=1 solve at 1.01x plateaus and never converges.** M=24
relay sequence (200 LM iters each, ScaleProblem=jacobian):
0.72 -> 4.887e-3 -> 2.566e-3 -> 1.975e-3 -> 1.826e-3 (per-relay improvement
99.3% -> 47% -> 23% -> 7.5%; <10% => step failed per guard). M=48
escalation: 0.72 -> 3.308e-3 -> 3.264e-3 (1.3% => abandoned).
`FAIL run_anchor_up (eps floor Inf)`. Total wall ~95 min. Log: `anchor_up.log`.
**Plateau mechanism (from the LM traces):** at the M=48 plateau the
first-order optimality fell to 1.2e-2 (from 9.3e+2 at the seed, a ~8e4
drop) while ||R|| froze at 3.26e-3 — J'R -> 0 with R != 0. For a SQUARE
MS system that means the Jacobian is (numerically) singular at the plateau:
a genuine nonzero local minimum of ||R||^2, not a conditioning crawl that
more iterations would fix. M=24 shows the softer version (optimality
0.5–17, asymptotic deceleration to a ~1.8e-3 floor).
**Conclusion:** time-domain MS is ill-conditioned at ANY in-band tf — the
clean integration seed did NOT dodge the Task-6 pathology; at 1.01x even
eps=1 (maximal smoothing) stalls. Together with the 1.12x dual-seed block
this is decisive: the pre-authorized fallback triggers — Sundman-domain
MS rebuild (16-dim augmented state, tau as independent variable).
ScaleProblem=jacobian stays in (mintime regression green).

### 2026-07-10 — Task S1 Gates A+B: 16-dim Sundman EOM + Jacobian GREEN
Sundman-domain MS machinery built (sms_ file set: eom/problem/pack/unpack/
residual/jacobian_cs/traj/seed_mintime + tests). System: sigma with
dt/dsigma = kappa = r1^1.5, 16-dim Y = [r;v;m;t;lamR;lamV;lamM;lamT],
lamRdot = kap*(-G'lamV) - dkapdr*(Ht+lamT), Ht ENTROPY-smoothed
(Lear via CS-safe softplus identity; Lear(S=0) = -log(2) verified to
1e-15). **Gate A PASS** (test_sms_eom): H_sigma conservation 1.9e-12;
cross-domain state match — short-span 6.9e-10/2.0e-9, full-span 8.5e-8/
2.5e-7 vs measured time-domain self-noise 3.7e-8/1.9e-7 (flat 1e-9
full-span unattainable: STM through ~30 perigee passes amplifies
1e-13-tolerance noise; diag_s1_gateA.m attribution — matched-span
agreement 1.6e-11 proves the EOM). Brief's conservation check is BLIND to
the dkapdr sign on Hval=0 trajectories (drift rate prop. to Hval); added
offset-lamT short-span check with in-test flipped-sign assertion:
detection power 1.1e-14 vs 1.8e-1. **Gate B PASS** (test_sms_jacobian,
v3 design at 16 dims): all 12 columns pass FD Gate 1 at the FIRST h
(1e-3) — no col-49-style noise limitation; CS h-independence <=6e-13, B2
<=5e-14, structure exact. The regularized domain also fixes the FD
conditioning. sigf(1.00x) = 149.75.
