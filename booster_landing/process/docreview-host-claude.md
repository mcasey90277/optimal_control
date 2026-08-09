# Three-way document review — Anthropic slot (host Claude, fresh eyes)

**Reviewer:** host Claude (Opus 5), no involvement in authoring either document
**Date:** 2026-08-09
**Documents reviewed:**
- `doc/booster_landing_note.tex` (theory note, 1644 lines)
- `doc/booster_landing_sdd.tex` (software description document, 2258 lines)

**Method.** Both documents read in full. Because this slot has repository access
(the two external slots do not), load-bearing and checkable claims were verified
against the MATLAB source (`lib/`, `certify/`) and against
`results/booster_run.mat` via one read-only `matlab -batch` probe. Nothing in the
repository was modified. Where a document and the code disagree, the code is
reported as measured, not assumed correct — but in every case below the code was
also cross-read for context.

**Adjudicated, not re-litigated:** `etaT = 0.87`; `vf = [0;0;-1.5]`; the
G3 1.0 kg gate as a Taylor model-error measurement; the min–max structure; the
recently re-derived G5 midpoint primer comparison.

**Verified correct (worth stating, because these are the load-bearing pieces).**
- The G5 midpoint-stationarity derivation (note Eq. 20) is **right**, and I
  checked the algebra end to end. `B^T λ = λ_v/m − (λ_m/(I_sp g_0)) T̂`; both
  terms on the constraint side of Eq. (20) are radial in **T**, so
  `(4h/6) λ_v/m = ((4h/6)c + μ) T` and `λ_v ∥ T` exactly, independent of `μ_k`
  and of the nondimensional scaling. I confirmed in `solve_pdg_colloc.m` that
  midpoint states carry no glideslope / `z≥0` / `m≥mdry` rows (nodes only) and
  that the annulus is written as `Tmin² ≤ ||T||² ≤ (ηT·Tmax)²`, whose gradient is
  `2T` — so the stationarity condition really is as clean as the note claims.
  This is the strongest single argument in the note and it survives scrutiny.
- The Taylor-bound conservatism argument is correct. On `δz ≥ 0`,
  `1 − δz + δz²/2 ≥ e^{−δz}` (f(0)=f'(0)=0, f''≥0) and `1 − δz ≤ e^{−δz}`, so
  the imposed lower bound is above and the imposed upper bound is below the true
  band — a genuine inner restriction on both sides. The observation that the
  *lower* bound was already convex and that the quadratic buys representability
  rather than convexity is a real contribution and is stated precisely.
- The ARE gain formulas (Eq. 15), the pole ceiling `ω_r = √(q_r/(2m√(q_r r)+q_v))`
  (Eq. 16), and the inversion `q_r = m²rω⁴`, `q_v = m²rω²(4ζ²−2)` all reproduce
  by hand from the 2×2 Riccati equations. The Butterworth floor `ζ ≥ 1/√2`
  falling out of `q_v ≥ 0` is correct and is a nice observation.
- The lossless-convexification proof sketch is honest about the exact place it is
  loose (the degenerate `λ_v ≡ 0` branch is genuinely *not* closed off by
  transversality here, and the normality caveat is placed correctly).
- Every entry of the flagship gate table in the note reproduces bit-for-bit from
  `rep` in `booster_run.mat`. Monte-Carlo statistics reproduce exactly
  (99.5%, 199/0/1/0, miss 3.720/3.665/6.723/9.238, vtd 1.0535/1.0001/1.2055/8.2898,
  landed-only span 0.879–1.476, mprop 0.00/537.26/882.66). `out0` reproduces
  (miss 0.0072, vtd 0.9777, 854.98 kg, sat_frac 0.3506). `tf_curve` has exactly
  14 rows, matching "the flagship search evaluated 14 fixed-`t_f` convex
  subproblems". Bound fraction is 120/121. `Qf(6,6) = 9.2906`. Hoverslam ratios,
  `|r0|`, `|v0|`, `Isp·g0`, `ηT·Tmax`, the de-rate percentages, and the Phase-2
  arc arithmetic all check.
- The Phase-2 result is physically coherent. Backing out ΔV: vacuum
  `2765.5·ln(30000/26464.6) = 346.9` m/s against a `Δ|v| + g₀t_f` budget of
  343.7 m/s (≈3 m/s steering); drag `301.7` m/s against a budget of 352.9 m/s,
  implying ≈51 m/s of free drag impulse. At `z=2 km, v=180 m/s` the model gives
  `a_D = 5.6 m/s²`, so ≈51 m/s over an 11 s margin arc is the right order. The
  counterintuitive headline survives an independent budget check.

---

## Findings

### CORRECTNESS

- **[CORRECTNESS]** `note:§5.2 (The lateral pole ceiling)` — **The quoted
  old-weight closed-loop eigenvalues do not reproduce, and are internally
  inconsistent with the note's own mass.** The note states
  `−0.0929 ± 0.0497i` with a "10.8 s time constant". Computed from
  `booster_run.mat` with `q_r=1e-4, q_v=1e-2, r=1e-9` and `m_A` = mean mass over
  the margin arc (29 628 kg by the node-grid switch, 29 567 kg by the shipped
  `ctrl.tSwitch = 6.1552` — the latter is what `tvlqr_design` uses, confirmed by
  its `q_r = 0.2099`), the roots of `s² + (K_v/m)s + K_r/m` are
  **`−0.09047 ± 0.04989i`, τ = 11.05 s**. The two statements in the same
  paragraph pin different masses: `2m√(q_r r) = 0.0187` requires `m ≈ 29.6 t`
  (correct, measured 0.018738), while a real part of `0.0929` requires
  `m ≈ 28.7 t`. Fix: recompute the pair, quote it as `−0.0905 ± 0.0499i`
  (τ ≈ 11.0 s), and state `m_A` explicitly in the sentence so the number is
  checkable. The argument is unaffected — the envelope time constant is still
  two-thirds of the 16.6 s flight — but a wrong number in the one paragraph that
  overturns a prior campaign verdict is the worst place to have one.

- **[CORRECTNESS]** `note:§5.4 (The guidance thrust de-rate)` — **The
  reachability table is quoted without its configuration, and the "107 m" figure
  is not a comparison against the campaign's nominal.** The table comes from
  `process/task-7-report.md §5`, which was measured with the *then-shipped*
  `ηT = 0.93`, and whose own narrative reads "806.6 m instead of **the nominal
  718 m**" — i.e. 88.6 m, not 107 m. The note's 107 m is the difference between
  two rows of the model's own bisection column (806.6 − 699.6). Meanwhile the
  *shipped* `ηT = 0.87` feed-forward brakes at **914.7 m** (the note's own
  §4 table). Three distinct "nominal brake altitudes" (699.6 / 718 / 914.7) are
  therefore live in a passage whose whole force is "a −5% vehicle must brake
  *this much* higher than nominal, and a fixed-reference tracker cannot".
  Fix, in two parts: (a) state the table's configuration and model assumptions
  (vertical-only, all thrust vertical, zero tracking error, `ηT = 0.93`) so a
  reader stops trying to reconcile 699.6 with 914.7; (b) **close the loop, which
  the note currently leaves implicit and which is the strongest version of its
  own argument**: at the shipped de-rate the guidance brakes at 914.7 m, above
  the 806.6 m the −5% vehicle needs, so the de-rate's sufficiency can be stated
  in the reachability table's own currency rather than only as "reserved thrust
  headroom". As written, a sharp reader can construct that inequality and wonder
  why the note didn't.

- **[CORRECTNESS]** `note:§5.4` — **Two margin metrics are compared as if
  commensurate.** "A −5% thrust scale costs 7.2–7.7% of the net deceleration
  `Tmax/m − g0` … so any de-rate that reserves less than that is buying nothing"
  puts a *fraction of net deceleration* on one side and a *fraction of thrust*
  (`1 − ηT = 13%`) on the other. At these masses 13% of thrust is ≈19% of net
  deceleration; and the `ηT = 0.93` rebuttal two paragraphs later
  ("0.95/0.93 = 2.1% of margin") is back in thrust units again. The conclusion is
  right in every currency, but the inference as literally written does not follow
  from the numbers as literally given. Fix: convert once and stay there — e.g.
  "reserving `1 − ηT` of thrust returns `(1−ηT)(Tmax/m)/(Tmax/m − g0)` of net
  deceleration, so `ηT = 0.93` returns 10% against a 7.2–7.7% loss (marginal) and
  `ηT = 0.87` returns 19% (2.5× cover)".

- **[CORRECTNESS]** `note:§Reproduction` — **"Every number quoted here was read
  from that file" is not true of the shipped artifact.** The top-level variables
  of `results/booster_run.mat` are exactly `{solC, solV, rep, ctrl, out0, mc, P,
  when}` (measured). There is no `solD`, `repD`, `ctrlD`, `Pd` or `mcD`, so
  *none* of Table 8 (vacuum vs drag), the §7.2 drag gate table, or Table 10
  (wind Monte Carlo) can have come from it; nor can the `ηT = 1.00 / 0.93` rows
  of Table 6, nor the §5 before/after tracking tables, nor the reachability
  table. Fix: scope the sentence — "every Phase-1 number quoted here was read
  from that file; Phase-2 numbers come from a `cfg.phase2 = true` run, and the
  sweep tables in §5 from the experiments recorded in `process/task-{5,7}-report.md`".
  This matters more than it looks: the note's credibility rests on
  claims-vs-evidence discipline, and this is the one sentence that over-claims.

### CONSISTENCY

- **[CONSISTENCY]** `note:§2.3 vs sdd:§2.3 and certify_pdg.m` — **The Taylor
  model-error band is quoted two different ways.** The note says the refinement
  sweep "held the gap in a **0.70–0.74** kg band"; the SDD says
  "**0.70–0.94** kg, non-shrinking under refinement"; `certify/certify_pdg.m`'s
  header (the source of both) says "held the gap in a **0.70–0.94** kg band".
  The note is the outlier. Fix: 0.70–0.94 in the note. (The change also slightly
  strengthens the note's point — a wider non-shrinking band is even less like a
  discretization tail.)

- **[CONSISTENCY]** `sdd:§6.2 (Why the grids differ) vs sdd:§5.3 and certify_pdg.m`
  — **The SDD quotes two different primer-angle sweeps.** §5.3 gives
  `1.72°, 1.14°, 0.85°, 0.57°` at `N = 20/30/40/60` (matching `certify_pdg.m`'s
  measurement table: 1.7238, 1.1424, 0.8534, 0.5741, and matching the note).
  §6.2 gives `≈1.7°, 1.19°, 0.92°, 0.57°` — stale numbers from an earlier
  measurement. Fix: make §6.2 quote the same four values as §5.3, or drop the
  numbers there and cross-reference.

- **[CONSISTENCY]** `note:§3 (G5) and sdd:§5.3` — **The arccos floor is asserted
  as one number where the code measures two.** Both documents say the corrected
  midpoint angle is `1.2×10⁻⁶` deg "at every grid, in vacuum and under drag
  alike". `certify_pdg.m`'s own table records **8.5e-07 deg for all four vacuum
  grids** and 1.2e-06 for all four drag grids; the flagship `rep.G5_primer_deg`
  is 1.20742e-06. Both are on the floor, but the floor is two ulps wide:
  `primer_deg` is a `max` of `acosd`, and near unity `acosd` is quantized —
  `acosd(1−eps/2) = 8.54e-07` deg (the smallest attainable nonzero value) and
  `acosd(1−eps) = 1.207e-06` deg. So the note's parenthetical
  "`arccos(1−ε) = 1.207×10⁻⁶` degrees" identifies the wrong ulp as the floor.
  Fix: "one to two ulps of `arccos` near unity (8.5×10⁻⁷ to 1.2×10⁻⁶ deg,
  depending on which representable double the dot product lands on)". This is a
  small correction that makes the claim *more* convincing, because it explains
  why the vacuum and drag columns differ by a factor of 1.4 and neither is
  physical.

- **[CONSISTENCY]** `sdd:§9 item 6` — **The stated variable list of the shipped
  `booster_run.mat` is wrong.** The SDD says its top-level variables are
  "`P, ctrl, mc, mcD, out0, rep, repD, solC, solD, solV, when` — no `ctrlD`, no
  `Pd`". Measured: `{solC, solV, rep, ctrl, out0, mc, P, when}` — no `mcD`, no
  `repD`, no `solD` either. The shipped artifact is Phase-1 only. Fix: correct
  the list, and say the file is from a default (`phase2 = false`) run, which is
  also what §4.1's "may see a Phase-1-only schema" anticipates.

- **[CONSISTENCY]** `sdd:§1.4 (Repository layout)` — **"Seventeen `.m` files
  ship" does not match its own decomposition.** `1 + 1 + 8 + 2 + 4 = 16`
  (measured: 16 non-test `.m` files, 11 tests, 27 total). Fix: "Sixteen `.m`
  files ship … plus the 11 test scripts (27 in all)".

- **[CONSISTENCY]** `note:§2.1 (Why nondimensionalize)` — "the production solve
  at `N = 60` converges in **30** IPOPT iterations"; `solC.stats.iter = 32`
  (measured). Fix: 32, or "about 30".

- **[CONSISTENCY]** `note:§6 (Monte Carlo)` — two small mismatches against
  `mc`: run 84's drawn `Δr0` is `[−34.2; 425.3; 121.6]` m (note: `−34.3`), and
  the drawn thrust range is `[0.9561, 1.0447]`, i.e. `−4.4%/+4.5%` rather than
  "±4.5%". Trivial individually; worth fixing because the surrounding paragraph
  is explicitly an audit trail.

- **[CONSISTENCY]** `note:§5.2 vs §5.5` — **The same dispersion has two different
  outcomes with no time-stamp.** §5.2's before/after table reports
  `Δr0 = [50;−30;0]` landing at miss 0.944 m / `v_td` 1.363 m/s "after"; Table 7
  (the acceptance battery, shipped defaults) reports 1.582 m / 0.968 m/s for the
  same case. Both are presumably right for their configuration (the §5.2 table
  predates altitude indexing and the 0.87 de-rate), but nothing says so, and a
  reader checking the note against itself will find a contradiction. Fix: label
  the §5.2 (and §5.3) tables as mid-campaign snapshots at the then-current
  configuration, superseded by Table 7. Same treatment for §5.3's
  `0.977 m/s / 0.395 m`, which *does* match Table 7 and so currently makes the
  §5.2 mismatch look like an error rather than a chronology.

### CLARITY

- **[CLARITY]** `note:§2.2 (Route B)` — **The note never states Route B's
  discretization.** §2.2 develops the continuous convexified problem and the
  Taylor bounds, then jumps to "the fixed-`t_f` subproblem is a genuine convex
  program"; the only mention of *how* it is discretized is `N_conv` = "trapezoid
  nodes, convex solver" in Table 1. From the source it is trapezoidal on the
  linear dynamics with `u` and `σ` piecewise-linear, second order, against Route
  A's Hermite–Simpson. That order gap is load-bearing for §2.3 and deserves two
  sentences in §2.2. Fix: state the transcription (trapezoid, `N_conv = 120`
  nodes, `u`/`σ` as node decision variables) where the convex program is
  introduced.

- **[CLARITY]** `note:§2.3` — **The restriction argument is applied
  asymmetrically.** "The convex problem is a restriction of the true
  continuous-time problem, so its optimum is bounded above by the true optimum"
  is a statement about the *continuous* convex program; what is solved is that
  program *plus* a second-order trapezoid transcription, whose error has unknown
  sign — exactly the caveat the very next sentence applies to Route A. The
  evidence that closes this is already in the paragraph (the `N_conv ∈
  {120,180,240,300}` sweep at fixed `t_f`), but it is presented as excluding
  "discretization" generically. Fix: say explicitly that the `N_conv` sweep is
  what rules out Route B's own trapezoid error, and that the residual is then
  the only remaining candidate.

- **[CLARITY]** `note:§2.2 (Tightness)` — Two half-sentences would close the
  gap between the theorem stated and the problem solved. (i) The theorem is
  stated for the exact relaxation (Eq. 8–10); what is solved is the Taylor
  restriction (Eq. 12). The argument transfers unchanged — the Taylor bounds
  alter only the `σ`-admissible set, never the pointwise minimization of
  `λ_v·u` over `‖u‖ ≤ σ` — and saying so in one clause pre-empts the obvious
  objection. (ii) Eq. (12)'s upper bound `μ₂(1 − δz)` requires `δz < 1` to be
  non-vacuous; that holds here with enormous room (measured max `δz = 0.034` on
  the flagship), and one parenthetical makes the restriction non-vacuous by
  measurement rather than by assumption.

- **[CLARITY]** `note:§4.2` — "Pontryagin's principle **mandates** bang–bang with
  a **bounded number of switches** when there is no singular arc" is stronger
  than PMP gives. PMP gives bang–bang off the singular set; a *bound* on the
  switch count is an additional result requiring extra structure. And the absence
  of a singular arc here is measured (bound fraction 0.9917, one switch, at
  every grid) rather than proved. Fix: "PMP forces bang–bang wherever the
  switching function is nonzero; it does not mandate which bound comes first,
  and it does not by itself bound the number of switches — that the count is one
  here is measured, not derived."

- **[CLARITY]** `note:§5.4` — **The reachability table's first column is
  non-monotone and unexplained**: minimum achievable `|v_td|` reads
  1.150 / 0.104 / 0.635 / 0.755 m/s at thrust scales 1.00 / 0.98 / 0.95 / 0.90,
  with the *nominal* vehicle the worst of the four. That is almost certainly a
  bisection-residual artifact (the bisection stops when it can hit the target,
  so the column is a convergence residual, not a physical minimum), but as
  printed it invites the reader to conclude something false. Fix: one sentence
  saying what the column is, or drop it and keep only the brake-altitude column,
  which is the column the argument actually uses.

- **[CLARITY]** `note:§2.3` — "0.43 kg … or `1.6×10⁻⁵` of the vehicle mass".
  `0.428/30000 = 1.4×10⁻⁵`; `1.6×10⁻⁵` is against `m_f`. Fix: say which mass, or
  use 1.4×10⁻⁵ against `m₀`.

- **[CLARITY]** `note:§5.2 (pole ceiling paragraph)` — "the initial lateral
  command is `ω²|Δr| = 0.49 × 58 = 28 m/s² ≈ 840 kN`, which together with the
  ~338 kN vertical feed-forward **sits just outside the 845 kN rim**". The
  commanded magnitude is `√(840² + 338²) = 905 kN`, 7% outside, not "just".
  As written the 840-vs-845 juxtaposition reads as if the lateral command alone
  is the binding quantity. Fix: quote 905 kN and keep the conclusion.

- **[CLARITY]** `note:§1.1 / Table 1 / §7` — **The Phase-2 headline result has no
  aerodynamic-model caveat.** `C_d = 1.0` with `A = 10.75 m²` is the body frontal
  area (`π·1.85²`) with no base-drag, plume-interaction, or grid-fin
  contribution, and no Mach or attitude dependence. The 434.7 kg saving is very
  nearly linear in `C_d·A` over this range, so the headline number is a direct
  function of a single unvalidated coefficient. Fix: one sentence in §7.1 stating
  that the saving scales with `C_d·A` and that a retro-propulsive base flow makes
  the effective coefficient uncertain — this *protects* the result, since the
  mechanism (coast extension, switch moves later) is robust to the coefficient
  even if the 434.7 kg is not.

- **[CLARITY]** `note:§5.5` — "spending 35.1% of the flight on a thrust bound"
  is `out0.sat_frac` (measured 0.3506), a time-weighted duty cycle on *either*
  bound of the closed-loop command. Placed one page after a guidance solution
  that is 100% bang–bang by construction, it reads as a contradiction. Fix:
  define it inline ("the closed-loop command sits on an annulus bound 35.1% of
  the time, versus 100% for the open-loop reference — the feedback correction
  spends most of the flight strictly inside the annulus").

- **[CLARITY]** `note:§2.2 (Free final time)` — the shipped `tf_curve` shows the
  first golden probe (`t_f ≈ 14.58 s`) returning validity code 0 and entering the
  search as `−∞`. So the function the golden search actually optimizes is not
  `m(t_f)` but a `−∞`-substituted surrogate, and unimodality of the former does
  not by itself license golden section on the latter. It happens to be fine here
  because the infeasible set is a left tail. Fix: one clause — "infeasible probes
  enter the search as `−∞`; since the infeasible region is a left tail (too
  little time to arrest), the surrogate remains unimodal".

- **[CLARITY]** `sdd:` (document level) — **The SDD has no as-built stamp.** It
  is a description of software at a moment, it quotes measured numbers, and it
  documents a `.mat` whose schema has already drifted from what §9 item 6 claims.
  Fix: a one-line "as-built against commit `<sha>`, `booster_run.mat` of
  `<R.when>`" under the title. This is cheap and it is exactly what keeps §9 from
  going stale silently again.

### STYLE

- **[STYLE]** `note:abstract` — "a **five-gate** certification layer" while the
  layer reports six blocks (G0, G1, G2, G2ff, G3, G4, G5). §3 explains this
  carefully; the abstract does not. Fix: "a five-gate certification layer plus a
  time-base gate on the basis the other five are computed on" (the SDD abstract
  and §1.1 already phrase it this way — the note is the one out of step).

- **[STYLE]** `note:§5` — the section is the most valuable in the document and
  the least navigable, because five of its tables are snapshots at different
  configurations (pre-/post-altitude-indexing, `ηT` 0.93 vs 0.87). A single
  two-column "configuration at which each table was measured" line under each,
  or one summary sentence at the head of §5, would remove most of the
  cross-checking friction a reader currently hits (and would subsume three of
  the CONSISTENCY findings above).

- **[STYLE]** `sdd:§3` — the module catalog is excellent and would be more useful
  still with a one-line "stability" marker per interface (frozen / expected to
  change), since the document's stated audience is someone extending the code.

---

## Per-document assessment

### Theory note — `booster_landing_note.tex`

This is a strong document and, unusually, an honest one. Its distinguishing
feature is that it reports the *history* of each result — what was measured,
what was concluded, what was wrong, and what retracted it — and it does so
without softening. §3's account of the G2 reconstruction (a gate that was
failing a physically accurate solution by a factor of 60 700, i.e. a false
negative), §5.2's overturning of the campaign's own "actuator-authority limit"
verdict, and §7.2's withdrawal of the 10° drag loosening with the accompanying
lesson ("an analysis that establishes there is *no mechanism* for an observed
effect is evidence that the measurement is wrong") are the parts of the note
that would survive being read by a skeptical reviewer, and they are the parts
that make the certified numbers believable.

The mathematics is sound where I checked it, which is everywhere it is
load-bearing. The G5 midpoint-stationarity argument — the newest content — is
correct, and correct for the right reason; I verified the transcription supports
it (midpoint states carry no path constraints, the annulus gradient is radial).
The lossless-convexification sketch is honest about its one loose branch and
places the normality caveat correctly. The Taylor-bound conservatism argument is
right on both bounds and the observation that the lower bound was already convex
is a genuine improvement on the usual presentation. The ARE derivation and the
pole-ceiling result reproduce by hand.

The weaknesses are concentrated in one place: **§5.4, the de-rate section, is
the least defensible part of the note**, not because its conclusion is wrong but
because its supporting numbers arrive without provenance, mix two configurations
(`ηT = 0.93` for the reachability table, 0.87 everywhere else), compare
incommensurate margin metrics, and stop one inequality short of the strongest
form of their own argument. §5.2's eigenvalue pair is simply wrong and
self-inconsistent. And the "every number was read from that file" claim is not
supportable against the shipped artifact. Everything else is small.

**Assessment: strong. The mathematics holds; the fixes needed are numerical and
presentational, not structural.**

### Software description document — `booster_landing_sdd.tex`

As a standalone extension reference this is close to complete and is better than
most SDDs I have read for a research campaign. The parts that do the real work
are the ones an extender would otherwise have to rediscover: the `cfg.P`
derived-field re-derivation contract with its explicit "invariant to preserve";
the `lam_defect` nondimensional-scaling table with the note about which
comparisons are valid unrescaled; the `anyOn` conjunct in `hs_quad_ctrl`'s
transition test (with the correction of the file's own former comment); the
`landed ⇔ stop` invariant and why every consumer must gate on it; the
negative-`T_z` guard in `allocate_thrust`; the `Qf(6,6)` units invariant; the
figure-ownership contract; and §6.3's "no suite runner, on purpose" with the
MATLAB-script-scope reason. Each of those is a defect waiting to be reintroduced
by someone who doesn't know, and each is documented at the point where it would
be reintroduced. §5.5 ("what the gates can and cannot catch", with the *blind to*
column) and §9 (code-vs-documentation discrepancies, kept rather than deleted)
are both unusual and both correct instincts.

Cross-consistency with the note is good: interface names, gate semantics,
thresholds, adjudicated values and the skip patterns all match, and the SDD
consistently defers derivations rather than re-deriving them. The failures are
three stale numbers (the Taylor band, the §6.2 primer sweep, the `.mat` variable
list) and one arithmetic slip ("Seventeen"), all of which are the same failure
mode: a number transcribed once and not re-checked when its source moved. The
as-built stamp recommended above is the systemic fix.

Gaps that remain, none large: no per-interface stability marker; no worked
"how do I add a gate / a dispersion / a viz product" walkthrough (the invariants
are all present, but assembled by the reader); and §9's own drift shows the
document needs a refresh trigger tied to the flagship run.

**Assessment: strong, and genuinely usable standalone. Fix the four transcription
errors and add the as-built stamp.**

---

## Submit-readiness verdict (theory note)

**Not yet — but close: one focused revision pass, not a rewrite.**

The note is technically sound and the mathematics survives independent checking.
Four things should be fixed before it goes out, and all four are bounded:

1. **§5.2 eigenvalue pair** — recompute (`−0.0905 ± 0.0499i`, τ ≈ 11.0 s) and
   state `m_A`.
2. **§5.4 reachability table** — give it its configuration (`ηT = 0.93`,
   vertical-only model), correct or drop "107 m higher than nominal", and add
   the 914.7 m > 806.6 m inequality that completes the argument.
3. **§5.4 margin metrics** — put the de-rate reserve and the dispersion loss in
   the same units.
4. **§Reproduction** — scope "every number was read from that file" to Phase 1.

Plus the low-cost consistency fixes (0.70–0.94 kg; 32 iterations; the arccos
floor as a two-ulp range; the mid-campaign table labels) and, if there is
appetite for one substantive addition, the two sentences on Route B's
discretization in §2.2 — currently the only place where a reader cannot tell
what was actually solved from the note alone.

With those, this is a document I would be comfortable defending. The min–max
result, the G5 midpoint stationarity argument, the pole-ceiling inversion, and
the drag coast-extension mechanism are each worth publishing, and the campaign's
willingness to record its own retractions is what makes the rest credible.
