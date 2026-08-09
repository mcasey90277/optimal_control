# External code-review fix wave — report

Campaign: `booster_landing`. Date: 2026-08-09.
Inputs: `process/external-review-gpt56terra.md` (GPT-5.6-terra),
`process/external-review-gemini31pro.md` (Gemini 3.1 Pro), and the SDD's own
code-versus-documentation list (`doc/booster_landing_sdd.tex` §9).

Everything below was applied, re-certified, and documented. MATLAB R2025b,
synchronous `-batch` throughout.

---

## Headline

The single biggest finding was real, and larger than it looked. **G5's
primer-vector check was comparing a segment's Hermite–Simpson defect dual
against the wrong control.** Fixing the comparison collapsed the campaign's
two longest-standing "unexplained artifacts" — a vacuum discretization
effect (0.574°) and a drag-model degradation (2.605°) — to the
double-precision floor of `acos` (1.207e-6°) at *every* grid, in vacuum and
under drag alike. Both had been documented at length as open questions.
Neither was about physics; they were one bug.

The consequence is that a gate this campaign had *loosened* (1° → 10° under
drag) is now **1000× tighter than the vacuum threshold it replaced**:
0.01°, one number for both configurations.

---

## F1 — G5 primer time-shift (both reviewers; Gemini gave the fix)

**What was wrong.** `certify_pdg.m` computed
`Tdir = solC.U(:,1:end-1)./|...|` — the *node* controls — and compared them
against `solC.lam_defect(4:6,:)`, the duals of the *segment* defect
constraints.

**Why it is wrong.** `lam_defect(:,k)` is the multiplier of segment `k`'s
defect. The only control whose discrete stationarity condition involves that
multiplier *alone* is the segment's midpoint control `Um(:,k)`: it enters
the NLP through exactly one constraint block (segment `k`'s defect, via the
Simpson weight `4h/6`) plus its own annulus row, whose gradient is parallel
to `Um` itself. So

```
0 = dL/dUm_k = -(4h/6) * B(x_m)' * lam_k + mu_k * Um_k
```

and since `B`'s thrust-touching rows are `I/m` (velocity) and
`-T'/(|T| Isp g0)` (mass, itself radial in `T`), this forces `Um_k` exactly
parallel to `lam_v,k` — with no discretization error at all. A node control
appears in **two** adjacent defect blocks *and* inside both neighbours'
Hermite interpolants `x_m = ½(x_k+x_{k±1}) ± (h/8)(f_k - f_{k±1})`, so its
stationarity condition mixes `lam_{k-1}` and `lam_k`. That is GPT's finding
#3 stated precisely — and it turns out the transcription *does* have a place
where the naive comparison is exact; it is just not the node.

**Fix.** `Tdir = solC.Um ./ sqrt(sum(solC.Um.^2,1))`. `Um` is 3×N, matching
`lam_defect`'s N columns exactly, so the old truncation of the last node
disappears too.

**Measurement — grid sweep** (same duals, node vs midpoint comparison):

| N | case | node_deg | mid_deg |
|---|---|---|---|
| 20 | vac | 1.7238 | 8.54e-07 |
| 30 | vac | 1.1424 | 8.54e-07 |
| 40 | vac | 0.8534 | 8.54e-07 |
| 60 | vac | **0.5741** | 8.54e-07 |
| 20 | drag | 8.6275 | 1.21e-06 |
| 30 | drag | 5.1122 | 1.21e-06 |
| 40 | drag | 3.6264 | 1.21e-06 |
| 60 | drag | **2.6051** | 1.21e-06 |

The node column halves as N doubles — the signature of an O(h) time-base
error. The midpoint column is *constant*, and it is not a physical number:
`acosd(1-eps) = 1.20742e-06`. Alignment is exact to floating point.

**Gate decision — RE-TIGHTENED to 0.01°, one threshold for vacuum and drag;
the drag branch is deleted.** The brief asked for "~3× margin over the
measurement". That is not applicable here: the measurement *is* machine
epsilon, so 3× it would be 3× numerical noise. What the angle now measures
is the solve's KKT stationarity residual, so the honest calibration is
against **solver convergence depth**, measured at N=60 by varying `opts.tol`
alone:

| tol | vac_deg | drag_deg |
|---|---|---|
| 1e-04 | 0.03704 | 92.113 |
| 1e-06 | 5.09e-05 | 0.00758 |
| 1e-08 | 1.21e-06 | 4.62e-04 |
| 1e-09 (vac default) | 8.54e-07 | 2.40e-05 |
| 1e-11 (drag default) | 1.21e-06 | 1.21e-06 |

0.01° carries ~4 orders of margin over both shipped defaults, tolerates 3
orders of convergence degradation, and still **fails** a genuinely
unconverged solve (the `tol=1e-4` row — which the old 10° drag gate would
have passed in the vacuum column at 0.037°). A sign flip or gross
misalignment lands at 90–180°, nowhere near the boundary.

**Cone guard re-measured.** With `P.theta_max_deg = 15`, the corrected
midpoint primer reads **7.322°** — essentially identical to the node value
(7.322°). Unlike the vacuum/drag artifact, the cone gap is a real
constrained-vs-unconstrained direction offset that survives the time-base
fix, so the `'skipped-cone'` guard is still required — and matters *more*
now that the bound is 0.01° rather than 1°. (Midpoint tilt on that solution:
max 15.000°, min 2.262° — it really does ride the cone boundary.)

**Docs.** The whole "G5 primer LOOSENING under drag" header note in
`certify_pdg.m` is replaced by a "G5 primer TIME BASE" note carrying both
tables and the superseded explanation. The theory note's §Certification now
derives the midpoint stationarity condition (new eq. `midstat`) and its
Phase-2 section records the withdrawal explicitly, including the
methodological lesson: *an analysis that establishes there is no mechanism
for an observed effect is evidence that the measurement is wrong.* The old
note had correctly proved drag cannot change the direction condition (drag
enters `A`, never `B`) and then attributed the residual to the duals anyway.

---

## F2 — `annulus_switch` no-switch fallback (Gemini)

**What was wrong.** The fallback returned `ts = 0`, intending "phase-B
weights over the whole flight". The blend is
`Q(t) = (1-w)QA + w*QB`, `w = ½(1+tanh((t-ts)/tb))`, so `ts=0` gives
`w(0)=½` — a 50/50 QA/QB blend at t=0, injecting half the aggressive
phase-A lateral gains into exactly the opening transient the fallback exists
to keep conservative.

**Fix.** `ts = -10*tBlend` (saturates the tanh: `w(0) = 1 - 2.1e-9`).
Expressed in units of `tBlend` so it stays saturated for any caller-supplied
blend width; `ctrl.tBlend` is therefore now resolved *before*
`annulus_switch` is called, and the local function takes it as a third
argument.

**Unit check** (forged no-switch profile: real solution with all `U`/`Um`
rescaled to `etaT*Tmax`, so no mid-annulus crossing exists):

```
tSwitch (new fallback)  = -4.0000 s   (expect -10*tBlend = -4.0000)
||Qfun(0)-QB||/||QB||  NEW (ts=-10*tb) = 2.242e-11   <- pure QB
||Qfun(0)-QB||/||QB||  OLD (ts=0)      = 5.439e-03   <- 50/50 blend
Qfun(0)(1,1): NEW 0.0001   OLD 0.108095   QB(1,1) 0.0001
F2 unit check PASS
```

`QB` is taken as `Qfun(1e4)` (blend fully saturated). The (1,1) entry makes
the failure concrete: the old fallback put a lateral position weight of
0.108 where 1e-4 was intended — three orders too aggressive.

---

## F3 — mass depletion misclassified as `arrest` (GPT)

**What was wrong.** `sim_closed_loop` mapped terminal events {2: `vz=0`
arrest, 3: `m=mdry` depletion} both to `stop='arrest'`, so
`run_monte_carlo` reported propellant exhaustion as a tracking failure.

**Fix.** A distinct `'depleted'` class; `landed` stays false.
Classification order is touchdown → depleted → arrest → horizon: if both
[2] and [3] fire, the arrest is a *consequence* of the dry tank, so the
cause is the more useful label. Added `mc.n_depleted`, the print paths in
`run_booster_landing` (both phases), `plot_footprint`'s annotation box, and
a new exhaustiveness assertion
(`run_monte_carlo:stopClassesIncomplete`) that the four counts partition the
run set — without it, a fifth class added later would vanish from the
breakdown, which is exactly how `depleted` hid inside `arrest`.
`test_monte_carlo_small` also now asserts every `mc.stop` string is one of
the four documented classes (a sum check alone cannot catch a miscount that
cancels).

**200-run vacuum MC re-run** (`rng(P.seed)=42`, N=60, deterministic):

```
success_rate = 0.9950  (199/200)
landed 199 / arrest 0 / depleted 1 / horizon 0
  run  84  stop=depleted   mprop=   0.000 kg  miss=   4.031 m  vtd= 8.290 m/s  ok=0
not-ok runs: 84
```

**Run 84's true class is `depleted`** — `mprop = 0.000 kg`, i.e. it burned
to `mdry`. Success rate unchanged at 99.5%, as expected: this is a relabel,
not a re-scoring. The theory note's Table `mc` row, its run-84 prose ("and
arrested rather than touching down" → "the simulation terminated on the
propellant-depletion event"), the README flagship row, and the SDD were all
updated.

The Phase-2 wind MC was re-run too (to regenerate its footprint annotation)
and reproduced exactly: 99.0% (198/200), landed 200 / arrest 0 / depleted 0
/ horizon 0, miss mean 6.02 / p95 12.38 / max 15.52 m, vtd max 1.125 m/s,
mprop min 477.3 kg — every published Phase-2 number unchanged. Both
`footprint.png` and `phase2_footprint.png` regenerated (in `results/` and
`doc/figs/`) so their breakdown boxes carry the new class.

---

## F4 — cold-start `tf` guess outside its own bracket (GPT)

**What was wrong.** `opti.set_initial(tf, 30/Tc)` — a hardcoded 30 s, stale
from before the bracket narrowed to `[P.tf_lo, P.tf_hi] = [10, 22]` s. Every
cold solve therefore began outside a simple scalar bound and IPOPT spent
restoration iterations pushing it back in.

**Fix.** `opti.set_initial(tf, 0.5*(P.tf_lo + P.tf_hi)/Tc)`. Since
`Tc = 0.5*(tf_lo+tf_hi)` this is exactly 1 in scaled units; written as the
ratio so a future bracket change carries the guess with it rather than going
stale a second time.

**Verification — same optimum.** `test_colloc_smoke` (cold, N=60):
`tf = 16.60 s, mf = 26464.3 kg`. Nominal certification block:
`tf = 16.595, mf = 26464.6 kg`, matching the published flagship to the
digits the README quotes.

---

## F5 — `certify_pdg`'s "never throws" contract (GPT; also SDD-2)

**What was wrong.** The header promised "report-only and never throws"; the
file carried two `assert(...,'certify_pdg:timebase',...)` calls. The failure
direction is the bad one: a malformed or failed-solver `solC` — precisely
the input a report is most wanted for — threw instead of producing one,
taking every other gate's diagnostics with it.

**Fix.** The two checks are now **gate G0**, scored like any other:
`rep.G0_tf_err`, `rep.G0_dt_err`, `rep.G0_tf_tol`, `rep.G0_dt_tol`,
`rep.G0_pass`, `rep.G0_msg`. The assertion messages are preserved *verbatim*
in `G0_msg` and printed under the G0 rows on failure, so nothing readable
from the old exception text is lost — it is delivered as data instead.
`rep.all_pass` includes G0, so no caller can certify on an inconsistent time
base. Thresholds are relative (`1e-9*max(1,|·|)`) and are stored in the
report rather than left for the printer to guess — the same discipline the
`tolScale`/effective-threshold machinery exists for.

---

## F6 — division guards in `pdg_dynamics` analytic Jacobians (both)

**What was wrong.** `A(4:6,4:6)` divides by `vmag` and `B(7,:)` divides by
`Tmag`; both are `0/0 → NaN` at zero. Neither is hypothetical:
`sim_closed_loop`'s arrest mode passes through `vz=0` by construction, and
`tvlqr_design` requests `A` and `B` along the whole trajectory.

**Fix.** Both magnitudes are now `sqrt(sum(.^2) + 1e-12)` — the
`sqrt(sum+eps)` form specified, matching the CasADi twin
`pdg_rhs_casadi`, which has always carried the same `1e-12`. No `abs`,
`norm` or `max` on any complex-step path. Two properties make this the right
guard rather than a patch:

* the smoothed forms have the correct limits (drag Jacobian → 0, `B(7,:)`
  → 0);
* the analytic Jacobians remain **exact** derivatives of the smoothed
  `xdot` — `d/dv [-kD/m · vmag · v]` with `vmag = sqrt(v'v+e)` is exactly
  `-(kD/m)(vmag·I + vv'/vmag)` — so complex-step still validates them to
  machine precision.

Magnitudes here are SI (`|T| ~ 1e5` N, `|v| ~ 1e2` m/s), so the guard is
~1e-22 relative: invisible away from the singular point. A side benefit is
that the MATLAB and CasADi RHS now smooth *identically*, removing a small
inconsistency G1 was silently absorbing.

**Verification.** `test_dynamics_jac PASS (vacuum + drag)` — complex-step
`h = 1e-30i` against `A` and `B` at the file's own `1e-12` tolerance.

---

## F7 — documentation corrections in code

| Item | File | What changed |
|---|---|---|
| G3 trajectory row scored | `certify/certify_pdg.m`, `print_certify_report.m` | `G3_traj_Linf` now enters `G3_pass` at a 5.0 m gate and prints as a scored row (was `(info)`). Threshold dimensioned from the mission — one third of `P.pad_radius = 15` m — not from the measurement. Measured 0.373 m, a 13× margin. |
| SDD-3 | `viz/plot_pdg_solution.m` | "max-min-max money plot" → **min-max**, with a note on why (Tmin exceeds weight at every mass, so there is no hover-capable opening arc). |
| SDD-4 | `viz/movie_landing.m` | `sol` documented as **unused**; the throttle axis comes from `out.t(end)` (the closed-loop duration, which is the right choice). Argument kept — it is third of five and positional, so dropping it would silently re-bind `outfile`/`opts` in every caller. |
| SDD-5 | `lib/solve_pdg_colloc.m` | `.tol` INPUTS line now states both defaults: 1e-9 vacuum, 1e-11 when `P.drag.on`. |
| SDD-1 | `run_booster_landing.m` | `OUTPUTS` header now gates `.mc` behind `cfg.doMC`, matching how it already documented `.mcD`. |
| SDD-7 | `README.md` | Test bootstrap: the real cwd-independent form is given in the folder map, with an explicit note that the shortened `addpath('..')` only works from `tests/` and must not be copied. |
| SDD-9 | `README.md` | `results/` row now states the `cfg.doMC` condition on `phase2_footprint.png` (and on `footprint.png`, under the same gate). |
| SDD-8 | `tests/test_monte_carlo_small.m` | Stale pointer to "`run_monte_carlo.m`'s own probe at N=P.N" corrected — the 50-run probe is a one-off recorded in the task-8 report; that function contains no probe. |
| Gemini style | `lib/tvlqr_design.m` | Header ADAPTATION 5(b) corrected: lateral **velocity** weights (`qV`) are phase-scheduled too, not only position. The code was right; the comment was wrong. `qP` alone sets the undamped frequency — `qV` is what buys the specified damping ratio `zetaA`. |

---

## Parked items (ruled by the controller, not implemented)

Recorded here so the next reader knows they were considered:

* **Unbounded Gaussian draws** (GPT #5) — a sign reversal needs 66σ
  (thrust) or 100σ (Isp); unreachable at any feasible campaign size.
  Truncating would change every reported success rate to defend against an
  event that cannot occur. Noted in `run_monte_carlo`'s header and the SDD.
* **IPOPT globality wording** (GPT #11) — softened in one comment in
  `solve_pdg_convex.m`: globality follows from convexity of the objective
  and feasible set, *verified after the fact via G4 and G3*, not from any
  property of IPOPT (a local method) and not from handing an SOC to a cone
  solver (the annulus is in squared form here).
* **Golden-section gap mapping** (GPT #10), **rng stream isolation** (GPT
  #13), **G1 per-row scaling** (GPT #14, already SDD future work),
  **`hs_quad_ctrl` tangential-touch condition** (GPT #7, a documented
  limitation), **retry jitter** (Gemini) — not implemented.
* **`hs_quad_ctrl` zero-direction edge** (GPT #6) — defensive comment only.
  `Traw = 0` requires the quadratic through three vectors of magnitude ≥
  `Tmin = 338` kN to pass through the origin, i.e. a ~180° direction
  reversal inside one segment; the glideslope cone and min-fuel structure
  forbid it, and G2ff (20 samples/segment) would catch it. Handling it would
  add an untested branch to the function every component's control
  reconstruction routes through.

---

## Verification

### Full suite — 11/11 PASS

Each test in its own `-batch` call (per the README's own warning about
`run()` clobbering script workspaces):

```
test_params                PASS (min-throttle dry T/W = 1.346)
test_dynamics_jac          PASS (vacuum + drag)
test_colloc_smoke          PASS  tf=16.60 s  mf=26464.3 kg  fuel=3535.7 kg
test_convex_lossless       PASS  gap=4.81e-05  mf=26435.3 kg
test_tvlqr_riccati         PASS
test_closed_loop_nominal   PASS  (nom miss 0.007 m vtd 0.978 m/s, disp miss 1.58 m vtd 0.968 m/s)
                                 + 5/5 dispersion battery cases PASS
test_monte_carlo_small     PASS  success=65%  landed=20 arrest=0 depleted=0 horizon=0
test_viz_smoke             PASS
test_phase2_drag           PASS  fuel vac=3535.4 drag=3100.9 kg
test_run_front_door        PASS  (ALL GATES PASS)
test_certify_nominal       PASS  (coarse) + PASS (nominal, all_pass)
```

### Nominal-grid gate table (Phase 1, N=60 / Nconv=120, tolScale=1)

```
Gate                                    value              threshold   verdict
------------------------------------------------------------------------------------
G0 t(end) vs tf err [s]                     0          < 1.65952e-08   PASS
G0 grid spacing err [s]           2.77556e-15                < 1e-09   PASS
  -> G0 gate                                                           PASS
G1 max HS defect                  1.39238e-07                < 1e-06   PASS
G2 pos residual [m]                 0.0088486                    < 1   PASS
G2 vel residual [m/s]             0.000299075                  < 0.1   PASS
G2 mass residual [kg]             1.02479e-05                  < 0.5   PASS
  -> G2 gate                                                           PASS
G2ff below-Tmin frac              1.72212e-16                < 1e-06   PASS
G2ff above-Tmax frac              1.58356e-16                < 1e-06   PASS
  -> G2ff gate                                                         PASS
G3 |dmf| [kg]                         0.42792                    < 1   PASS
G3 |dtf| [s]                       0.00340864                  < 0.2   PASS
G3 traj Linf [m]                     0.372616                    < 5   PASS
  -> G3 gate                                                           PASS
G4 lossless gap [m/s^2]           0.000135597        < 1e-4*TmaxG/m0   PASS
G5 bound fraction                    0.991736                >= 0.95   PASS
G5 interior switches                        1                   <= 2   PASS
G5 structure (bound+switch+max-last)                                         PASS
G5 primer angle [deg]             1.20742e-06                 < 0.01   PASS
  -> G5 gate                                                           PASS
------------------------------------------------------------------------------------
ALL GATES                                                              PASS
```

Coarse block (N=40 / Nconv=90) also ALL PASS, with the same primer value
(1.20742e-06) — the grid dependence is gone.

### Nominal-grid gate table (Phase 2, drag, N=60)

`tf = 17.5308 s, mf = 26899.28 kg`:

```
G0 t(end) vs tf err [s]                     0          < 1.75308e-08   PASS
G0 grid spacing err [s]           2.94209e-15                < 1e-09   PASS
  -> G0 gate                                                           PASS
G1 max HS defect                  2.28368e-08                < 1e-06   PASS
G2 pos residual [m]                 0.0122795                    < 1   PASS
G2 vel residual [m/s]              0.00066087                  < 0.1   PASS
G2 mass residual [kg]             9.04675e-05                  < 0.5   PASS
  -> G2ff gate                                                         PASS
G3 cross-method (dmf,dtf,traj)             --                     --   skipped
G4 lossless gap                            --                     --   skipped
G5 bound fraction                           1                >= 0.95   PASS
G5 interior switches                        1                   <= 2   PASS
G5 primer angle [deg]             1.20742e-06                 < 0.01   PASS
  -> G5 gate                                                           PASS
ALL GATES                                                              PASS
```

Drag now passes the *same* 0.01° primer gate as vacuum, with four orders of
margin.

### Documents

Both recompiled twice, aux cleaned, no undefined references:

* `doc/booster_landing_note.tex` — §Certification gains the midpoint
  stationarity derivation and a G0 paragraph; flagship gate table replaced;
  §Phase-2 records the withdrawn loosening; MC table + run-84 prose
  relabelled; future-work item "the drag primer question" replaced by
  "higher-order costate recovery" (the old one is closed).
  `verify-paper`: **10 verified, 1 to review, 0 failed, 1 skipped**; 4/4
  figures OK. The single WARN (`betts2010`, a book with no DOI) is
  pre-existing.
* `doc/booster_landing_sdd.tex` — report struct, gate table, `tolScale`
  semantics, G5 primer semantics, G0 (replacing "Internal assertions"),
  stop classes, MC counts + invariants, `pdg_dynamics` guards, cold-start
  guess, `annulus_switch` fallback, test-grid rationale, extension points;
  and §9 rewritten with eight of nine items marked **FIXED** and their
  entries updated to describe what they now say. Item 6 (the shipped
  `booster_run.mat` predating `.ctrlD`/`.Pd`) stands — no full flagship was
  re-run.

### Commands

```
# probes (scratchpad)
matlab -batch "run('.../probe_primer.m')"    # node vs midpoint on saved solutions
matlab -batch "run('.../probe_grid.m')"      # N = 20/30/40/60, vacuum + drag
matlab -batch "run('.../probe_tol.m')"       # tol = 1e-4 .. 1e-11
matlab -batch "run('.../probe_cone.m')"      # theta_max = 15 deg
matlab -batch "run('.../verify_a.m')"        # F2 unit check + drag gate table
matlab -batch "run('.../verify_mc.m')"       # 200-run vacuum MC
matlab -batch "run('.../regen_figs.m')"      # both footprints + Phase-2 MC

# suite
for t in test_params test_dynamics_jac test_colloc_smoke \
         test_convex_lossless test_certify_nominal test_tvlqr_riccati \
         test_closed_loop_nominal test_monte_carlo_small test_viz_smoke \
         test_phase2_drag test_run_front_door; do
  /Applications/MATLAB_R2025b.app/bin/matlab -batch \
    "cd('<campaign>'); setup_paths; $t" || echo "FAILED: $t"
done
```

---

## Residual concerns

1. **No full flagship re-run.** Per the brief. `results/booster_run.mat`
   therefore still predates this fix wave (SDD §9 item 6 stands, and its
   stored `rep` has the old primer semantics). Everything the README and
   note quote was re-measured directly; the next no-args
   `run_booster_landing` will refresh the `.mat`.
2. **`hs_quad_ctrl`'s tangential-touch condition** (GPT #7) remains a real,
   documented limitation: a smooth interior arc that merely grazes a bound
   is rewritten as a full Tmin→Tmax step. Vacuous on the certified bang-bang
   solution, live the moment a genuine interior arc appears — which a
   pointing cone or a more aggressive drag case could produce. Parked, not
   solved.
3. **The `N=40` coarse test grid is now justified by history, not by a
   measurement.** Its original reason (the primer needed N≥40 to clear 1°)
   is gone. It is kept because it is the pair whose margins on the *other*
   gates are characterized and asserted; re-picking it would mean
   re-characterizing all of them. Documented in both tests and the SDD.
4. **G0 reports but does not halt.** Downstream gates still run when G0
   fails, so their numbers are computed on a time base known to be
   inconsistent. That is the intended behaviour (report everything, gate on
   `all_pass`), but a reader skimming a failed table should start at G0.
