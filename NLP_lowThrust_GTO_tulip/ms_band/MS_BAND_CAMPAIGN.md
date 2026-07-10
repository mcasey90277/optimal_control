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
