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

### 2026-07-10 — Task S1 Gate C: sms min-time anchor GREEN
prob.resFun wired into ms_solve + eps_march (default @ms_residual;
sms_problem sets @sms_residual); time-domain regression
test_ms_reproduce_mintime PASS after the edit (||R|| = 2.384e-10, 15
iters — identical record; regress_mintime_resfun.log). **Gate C PASS**
(test_sms_reproduce gateSel='C', M=24, eps=1e-3, seed sms_seed_mintime
1.00x): ||R|| = 9.711e-11 (19 iters, wall 1.6 min), dV 4.4665 km/s, prop
2.9247 kg, bang 100.0%, t(sigf)-tf = -7.8e-13. The sigma-domain system
even lowers the residual floor below the time-domain's 2.4e-10.

### 2026-07-10 — Task S1 Gate D: 1.12x dual-seed eps-march BLOCKED (healthy system, strategy mismatch)
Native sigma-domain dual seed validated (sms_seed_duals, M=50: beta =
0.03102, spread 0.45%, burn/coastAgree 100%, node1Err 0, arcCheckErr
2.25e-3, lamT dual relStd 8.7e-3; ||R(eps=1)|| = 4.285, costate joints
dominant). eps-march per the S1 brief ([1 0.3 ... 1e-4], 200 iters/step,
relays): the eps=1 step capped 3x — 4.285 -> 9.082e-3 -> 7.307e-3
(19.5%) -> 6.607e-3 (9.6% < 10% guard) -> abandoned. FAIL
test_sms_reproduce(D). **Trace capture (diag_s1_gateD.m/.mat): this is
NOT the time-domain pathology.** At the plateau iterate the GN system is
CONSISTENT to 6.4e-12 (R entirely in range(J); J full numerical rank) —
no singular residual minimum. Instead: (1) the eps=1 extremal is far
from the near-bang seed (terminal rv rows 4.4e-3 dominate; dV 4.22 vs
3.83; switches 0, bang 40.7% — the smooth extremal has no bang
structure); (2) cond(J) = 6.3e8 with GN step norm 129 (solution-scale,
outside the linear regime) — LM crawls at ~1e-4 nats/iter (~1.6e5 iters
to 1e-9). The Sundman machinery itself is certified (Gates A-C green,
lower residual floor, clean FD conditioning); the block is the
smooth-first continuation schedule applied to a SHARP seed.
Recommendation for the gate owner: march start near the source's
effective sharpness (e.g. eps0 ~ 1e-2..1e-3, sharp-end warm start)
instead of eps=1; alternatives: M=40, larger per-step budget (likely
insufficient alone). Logs: test_sms_gateD.log, diag_s1_gateD.log.

### 2026-07-10 — Task S1 Gate D retry (sharp-start) also plateaus; M=40 runs killed by external watchdog
Adjudicated retry (sharp schedule [1e-2 3e-3 1e-3 3e-4 1e-4], native dual
seed): eps=1e-2 step capped 4.184 -> 1.435e-2 -> relay 1.301e-2 (9.3% <
10% guard) -> abandoned. Trace capture (diag_s1_gateD2.m/.mat): again NO
singular minimum (GN consistent 5.7e-11) but cond(J) 2.6e9, ||dGN|| =
1131 (13x seed norm), LM moved only 9.6 in 400 iters; at the iterate the
trajectory is a ZERO-SWITCH 99.4%-saturated full burn (dV 4.54) — the
discrete-dual costate noise collapses the 12-switch structure (see
2026-07-10 switch-count adjudication: 10 certified + 1 near-graze dip) on
propagation. M=40 escalation (authorized): relay cuts 15.6%/14.6%,
relay 3 trending ~9.1e-3 — three successive runs killed by an external
~70-min watchdog at the identical wall point (run_s1_gateD_m40*.log);
kill-robust resume driver written (run_s1_gateD_m40b.m, state saved per
solve) but superseded by the external-review directive before relay 3
completed.

### 2026-07-10 — Task S1: external review + weighted dual-map adjudication (mode d wins)
GPT-5.6 external review (.superpowers/sdd/gpt56_review_S1.md): EOM/
Jacobian/residual confirmed correct; prime suspect = the interval-dual ->
node-costate map in sms_seed_duals. Four candidate maps implemented
(sms_seed_duals mode 'a' baseline left-interval, 'b' h-weighted, 'c'
adjacent-h averaged, 'd' midpoint-principled) + full-trajectory
validation harness test_sms_dualmap.m (one-arc propagation at 5 arcs
incl. perigee + switch arcs; |Ht+lamT| along the trajectory; FD-vs-RHS
adjoint defects). KKT derivation of the casadi trapezoid transcription:
interval multipliers are MIDPOINT costates with NO h-weights (objective
and defect quadrature weights cancel) — mode 'd'. Table
(legacy_ms_f1120, M=40): 'd' wins everything mode-sensitive — one-arc
costate err 0.497 -> 0.023 (21x vs baseline), |Ht+lamT| rms 2.40 ->
0.40, max 13.9 -> 1.41, adjV rms 2.4x; 'b' catastrophically falsified
(adj defects ~1e12): duals carry NO interval weights. Caveat: the mesh
has near-duplicate nodes (min h = 4.2e-12, 356 intervals < mean/10)
whose duals are noise (adjacent lamT-dual jumps to 0.89) — inflates the
FD adjoint metric for ALL modes (adjM/adjT floors ~1e2) and is a
residual seed-noise source no node map can remove. Robustness fixes
folded in: beta_from_duals error()s on no-switch data; sms_eom rejects
||lamV|| < 1e-8 via identified error (sms_residual catches it, returns a
large finite miss — documented threshold rejection, no regularization).
Logs: test_sms_dualmap.log, dualmap_table.mat, probe_s1_meshmin.m.

### 2026-07-10 — Task S1 close-out: Gate-D reproduction CLOSED (BLOCKED); verifier reframe delivers the Tier-1 adjoint-ODE certificate
Adjudication (user-approved): the 1.12x LM reproduction crawl is intrinsic
(full-rank GN-consistent J at every plateau; trust-region collapse on
1/eps switch layers per the GPT-5.6 methods verdict) — reproduction line
closed. Reframe: the Sundman 16-dim machinery + mode-'d' dual map now
serve as a VERIFIER. `verify_direct_pmp.m` (no LM anywhere): per-arc
propagation of the direct solution's own (state; costate) through
sms_eom, per-block defects vs the direct joints, along-arc |Ht+lamT|,
dual-S vs direct-throttle switch structure, primer alignment, terminal
transversality; table + verify_pmp_<name>.mat + two-panel png.
**1.12x verdict (legacy_ms_f1120, M=40, mode d, eps 1e-4):** arcs 1-39
at/below the dual-map floor (state blocks <= 1.03e-2, 37/39 under 1e-2);
primer 0.0971 deg mean (PASS); |lamM(sigf)| = 1.7e-6 (PASS); switches
matched 10/12 (primary stat; the 2 unmatched adjudicated below).

> **Switch structure: 10 PMP-certified switches (10/12 of the historical count
> matched to dual-S crossings at 0-1 nodes).** The remaining "pair" (#3/#4,
> tau 19.46-19.51, t = 0.277) is not a switch pair: the direct throttle only
> dips to u = 0.43 (never reaches the coast bound; 2 intermediate-throttle
> nodes of 4001), and the un-interpolated interval-dual switching function
> stays burn-side (S = -7.3e-4 at its local max toward zero) while separating
> all 542 true coast intervals (S >= +2.3e-3) from all 3446 true burn
> intervals (S <= -1.6e-3) at 100%. The historical "12 switches" was a
> counting artifact (s>0.5 threshold crossing a shallow dip over a
> near-grazing S). Whether the underlying extremal strictly burns
> (S_max ~ -8e-4) or marginally grazes there is unresolvable at this mesh;
> either way the certified switch count is **10 (+1 near-graze throttle dip,
> not a certified switch)**.

The arc-40 defect (0.68) is switch-crossing AMPLIFICATION, not costate
error: the last arc contains the terminal switch complex (5 threshold
crossings in the last 0.5 tau); defect <= 4.2e-3 up to the first in-arc
switch (tau 149.763, crossing matched at 0 nodes), then the u-vs-s
disagreement (25.9% of samples) integrates to O(1). The earlier "6
crossings" observation was seed-propagation degradation at eps=1e-2, not
grid resolution. Verdict framing: legacy_ms_f1120 is CONSISTENT WITH A
CONTINUOUS PMP EXTREMAL at its transcription's O(h^2) resolution. Review
verdict on the tool + adjudication: APPROVED/CONFIRMED (2026-07-10);
v2 gates the certificate text on all checks passing-or-adjudicated,
reports coast-bound dwell vs threshold crossings, and adds a p95 primer
row. Logs: run_verify_1120.log, diag_verify_1120.log; figure
verify_pmp_legacy_ms_f1120.png.
