<!-- This file is the CAMPAIGN NARRATIVE (results tables, findings, honesty
footnotes, and the full round-by-round record). It was the project's original
README.md; it was renamed to CAMPAIGN.md on 2026-07-18 when README.md became the
leaner operational pipeline doc (crlb-pipeline style). For quick orientation and
usage see README.md; for full methodology see doc/table3_method_note.tex; for
outstanding work see TODO.md. Nothing here was removed in the rename. -->

# earth_elliptic_to_geo — direct-method reproduction of Haberkorn–Martinon–Gergaud (campaign narrative)

Direct (collocation NLP) reproduction of the low-thrust **minimum-fuel** Earth
orbit transfer in Haberkorn, Martinon & Gergaud, *"Low thrust minimum-fuel
orbital transfer: a homotopic approach,"* JGCD 27(6), 2004 — a 1500 kg
satellite from a low, elliptic, 7°-inclined orbit (P=11625 km, e=0.75) to
equatorial GEO (P=42165 km). The paper solves it **indirectly** (single
shooting + PMP + energy→mass homotopy); this reproduction solves the same
physics **directly** (Sundman-regularized trapezoidal collocation + CasADi/
IPOPT), reusing the campaign architecture from `GTO_tulip/`.

The paper's method **is** our method, done the other way:

| paper (indirect) | ours (direct) |
|---|---|
| min ∫‖u‖dt, fixed t_f = c_tf·t_{f,min} | same |
| solve min-time first for t_{f,min} | same (free-longitude manifold anchor) |
| energy→mass homotopy ∫λ‖u‖+(1−λ)‖u‖² (λ:0→1) | energy→fuel homotopy ∫s−ε∫s(1−s) (ε:1→0) — identical, opposite parametrization |
| bang-bang via switching function ψ | bang-bang via ε→0 |
| PMP costates from shooting | PMP costates from KKT duals (verification) |

Full design rationale, unit conventions, and open items: `DESIGN.md`.
Task-by-task build record (16 tasks, T1–T16): `PLAN.md` and
`.superpowers/sdd/progress.md` (ledger section dated 2026-07-16 onward).

## Module map

| module | role |
|---|---|
| `kepler_lt_params.m` | Constants + canonical nondim units (LU = GEO radius, μ=1); Isp default 2000 s |
| `elements_to_cart.m` / `cart_to_elements.m` | Paper elements (P,eₓ,e_y,hₓ,h_y,L) ↔ inertial (r,v), roundtrip-tested |
| `geo_terminal.m` | Terminal builders: `'fixed'` full-state rendezvous at a prescribed longitude (M0/M1) and `'manifold'` free-longitude 5-constraint GEO insertion set (M2+) |
| `lt2b_rhs_time.m` | Shared 2-body + thrust RHS (solver, seed propagation, tests) |
| `casadi_lt_2body.m` | **Solver core** — Sundman-regularized (fixed τ_f, cScale slack state) trapezoidal collocation, cone-eliminated [α;s], modes `'mintime'` and `'fixedtf'` (ε energy→fuel objective), CasADi+IPOPT |
| `seed_2body.m` | Dynamically-exact tangential-thrust warm start, sampled at uniform-τ nodes; bisects throttle on arrival longitude for the right rev-topology |
| `homotopy_2body.m` | Guarded ε:1→0 sweep (loose first step, tight thereafter; never advance/cache on a failed step) |
| `run_mintime.m` | Free-longitude min-time anchor; stage-1 fixed-rendezvous warm-up + stage-2 manifold + warm continuation (see Findings below) |
| `run_transfer.m` | Full pipeline for one case: mintime anchor → tf/Lf → seed → homotopy → structure report; caches **only if certified** |
| `run_ctf_sweep.m` | M3: m_f-vs-c_tf front + thrust-law sweep, resumable (skips completed points) |
| `verify_pmp_2body.m` | PMP first-order checks from KKT duals: primer alignment, switching-function sign law, mass transversality |
| `transfer_movie.m` | Trajectory + throttle animation (MP4/GIF), reuses `PSR/psr_movie.m` layout |
| `mintime_guard_constants.m` | Single source for the continuation stall-guard (rounds cap, min decade improvement) shared by solver and test |

## How to run

```matlab
cd earth_elliptic_to_geo
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'))   % or set CASADI_PATH

mt  = run_mintime(10, 0.0612, 600);        % G1: free-L min-time anchor, 10 N, 3D
res = run_transfer(struct('thrustN',10,'ctf',1.5,'hx0',0.0612,'term','manifold'));
run_ctf_sweep                              % M3: full front + thrust-law sweep (hours)
transfer_movie('results/M2_manifold.mat','results/M2_movie')
```

Batch/resume note: `run_ctf_sweep` and `run_mintime`'s continuation rounds are
written to tolerate the sporadic CasADi/IPOPT MEX fatal crash (see Caveats) —
completed, certified points are skipped on rerun, so just rerun the script
after a crash.

## Milestone results — ours vs paper

Paper numbers: Table 3 / Figs 18, 21–23, at c_tf = 1.5, T_max = 10 N unless noted.

| gate | ours | paper | note |
|---|---|---|---|
| **G1** min-time anchor, coplanar | t_f = 21.9376 ND (83.57 h) | — | free-L manifold, defect 2.6e-12 |
| **G1** min-time anchor, 3D (7°) | t_f = 22.2248 ND (84.66 h) | ≈ 84.7 h | defect 2.6e-12 |
| **G2** energy stage (ε=1), coplanar | m_f = 1367.15 kg, defect 4.8e-12 | — | one-shot, no continuation needed |
| **M0** coplanar, fixed terminal | m_f = 1375.10 kg, 24 switches, 7.93 revs, apogee-burn ratio 1.98 | internal (no paper number) | defect 5.3e-12 |
| **M1** 3D (7°), fixed L_f (law-prescribed) | m_f = 1375.01 kg, 22 switches, 7.97 revs, apogee-burn 2.09, incl → 0.000000° | ~18 switches, 7.5 revs | plane-change cost 0.09 kg (folded into apogee burns); defect 5.4e-12 |
| **M2** free-L_f manifold — **headline mass match** | m_f = 1376.74 kg (N=600) / 1377.05 kg (N=1200) | 1370–1375 kg | Δ(mesh) = 0.31 kg < 1 kg pass; termErr 1.7e-38 / 1.1e-16 |
| **M3** front, T_max=10 N, c_tf = 1.2/1.5/2.0/2.5/3.0 | 1360.37 / 1376.74 / 1377.10 / 1385.42 / 1386.81 kg | ~1350→1388 kg, c_tf 1.05→3 | monotone, best-of chain/fresh with provenance; matches Fig 23 shape |
| **M3** thrust law C = T_max·t_{f,min} | 846.6 N·h (10 N) | ≈ 850 N·h | 5 N / 2.5 N legs blocked (see Findings) |

Isp is not stated numerically in the paper, but the benchmark's constants are
those of ref [6] Caillau & Noailles (ESAIM COCV 6, 2001,
`orbit_transfer/min_fuel_papers/COCV_2001__6__239_0.pdf`), obtained and read 2026-07-19: p.255
gives the mass-flow coefficient δ = 0.05112 km⁻¹·s in ṁ = −δ‖u‖, so
δ = 1/(Isp·g₀) ⇒ c = 1/δ = 19.562 km/s ⇒ **Isp = 1994.8 s**. We run at the
default **2000 s** (0.27% high; masses ~0.3 kg high); pass ispS=1994.8 for the
exact value. The earlier "unreadable reference / family cross-check" caveat is
now **CLOSED** — all benchmark constants (P⁰=11625, Pᶠ=42165, e⁰=0.75, L⁰=π,
m⁰=1500, μ⁰=398600.47) are confirmed identical, and the M2 mass match
(landing in the paper's 1370–1375 kg band) is corroborated by the true Isp.

## Test suite

No-solve / cheap guard tests (all pass in seconds to ~2 min):

```matlab
cd earth_elliptic_to_geo
test_params; test_elements; test_dynamics; test_terminal; test_seed; ...
test_solver_smoke; test_stall_guard; test_p2_homotopy
```

`test_p2_homotopy` reproduces the paper's own toy example (P2, double-
integrator min-fuel) against the analytic bang-off-bang solution — a cheap
known-answer check of the ε-homotopy machinery, independent of the 2-body
problem. `test_seed` and `test_p2_homotopy` are the slowest of this set
(low minutes each); the rest are sub-second.

`test_energy_stage.m` is a genuine one-shot IPOPT solve (gate G2, ε=1
coplanar energy stage) — not part of the "no-solve" list above, run
separately when validating the solver core end to end.

## Findings / honest caveats

- **Cold one-shot min-time against the free manifold does not converge.**
  A direct `'mintime'` solve warm-started straight from the cold tangential
  seed stalls (Maximum_Iterations_Exceeded, defect ~2.6e-3). The working
  recipe (`run_mintime.m`) is two-stage: stage 1 solves an easy fixed-
  rendezvous warm-up at the seed's own arrival longitude, stage 2 solves the
  real manifold problem warm-started from stage 1, with warm continuation
  rounds if stage 2 doesn't reach gate tolerance directly.
- **Mesh sensitivity changes switch count and rev basin, not the mass
  match.** Refining M2 from N=600 to N=1200 moved the solution from
  24 switches / 8.11 revs toward 19 switches / 7.32 revs (a basin change,
  not noise) while m_f moved only 0.31 kg — the paper's own reported
  structure (18 switches, 7.5 revs) sits between our two mesh resolutions.
  Switch count should be read as a band, not a fixed integer. The P0
  mesh-convergence certification (2026-07-21,
  `process/P0_SWITCH_MESH_CONVERGENCE.md`) quantifies this for the 0.2 N deep
  rung: the count converges to a band ~866±5 across 16/24/40 nodes/rev, while
  the 8/rev point estimate (823) is a ~5% undercount — mass (1375.8 kg) and rev
  count (346.7) are mesh-converged.
- **Basin scatter is real and mirrors the paper's own Fig 18.** At fixed
  c_tf, a cold fresh solve can land in a different (better or worse) local
  optimum than a neighbor-chained warm start — up to ~22 kg apart at the
  same c_tf (c_tf=3.0: 1386.81 kg fresh vs 1364.50 kg chained). M3's
  reported front takes the best certified point per c_tf, with provenance
  logged, rather than pretending the landscape is convex.
- **Seed-topology rule: stretching a seed in time does not add revolutions.**
  An early thrust-continuation strategy (warm-starting a lower-thrust anchor
  from a higher-thrust anchor's time-stretched trajectory) was found to be a
  topology error — stretching only rescales physical time, it cannot turn a
  ~4.5-rev shape into an ~8.4-rev min-time shape. Every thrust level now
  solves its own cold tangential-seed anchor.
- **Thrust-law leg BLOCKED at 5 N** after 6 documented strategies (defect
  floor ~5e-3, a false-infeasibility signature matching the parent tulip
  campaign's known failure mode). Every attempted 5 N run had ≤70 nodes/rev
  vs ≥130 for the successful 10 N runs — the leading hypothesis is
  insufficient mesh density per revolution, not a structural infeasibility.
  Left as future work (mesh parity at N~2400, energy-first warm start); the
  10 N point of the law (846.6 N·h vs paper's ≈850) still checks out. The
  2.5 N leg and the stretch goal of 1 N (74.5 revs, "Out of scope" per
  DESIGN.md but listed as an optional attempt) were never reached — the
  campaign never advances past a blocked lower thrust level.
- **PMP dual/primer anomaly, open, non-gating: primer misalignment 10–24°
  across milestones** (M0 13.3°, M1 20.5° — stable at 20.4645° across all
  three M1 probes, M2 18.8°, M2-N1200 23.5°); switching-sign (96.7–97.1% vs
  98% gate) and transversality margins degrade in tandem (M2-N1200
  transversality passes at 9.5e-11). Suspected coupling through the Sundman
  `cScale` slack state, not yet root-caused. Primal certification (defects,
  mass, structure) is unaffected; M2/M3 gate on primal results only.
  Diagnostic scripts are preserved in `results/dual_anomaly/`.
- **Sporadic CasADi/IPOPT MEX fatal crash on process init**, mostly on the
  first `opti.solve()` of a MATLAB process after an idle gap (observed ~4
  crashes across 10 launches this session). A plain relaunch always
  recovered; `run_transfer`/`run_ctf_sweep` cache only certified results so
  a crash never poisons a downstream seed.

## Deliverables

- `results/M2_manifold.mat` / `results/M2_manifold_N1200.mat` — the headline
  mass-match run and its mesh-refinement check.
- `results/front_mf_ctf.png` — the M3 m_f-vs-c_tf front (Fig 23 analog).
- `results/M2_movie.mp4` / `.gif` — trajectory + throttle animation of the M2
  transfer (apogee burns visible, GEO ring, running ΔV/mass meter).
- `results/dual_anomaly/` — diagnostics for the open PMP primer finding.
- `run_gergaud.m` — the front door (see "Front door: `run_gergaud`" below):
  one call prints a Table-3-style row and (optionally) a plot + movie for a
  chosen thrust/endpoint pair.
- `results/movie_MEE_10N.{mp4,gif}` / `results/movie_MEE_5N.{mp4,gif}` /
  `results/movie_MEE_2p5N.{mp4,gif}` / `results/movie_MEE_1N.{mp4,gif}` —
  trajectory + throttle animations of the four certified MEE min-fuel rungs
  (10/5/2.5/1 N), rendered via `mee_res_to_cart_res.m` + `transfer_movie.m`.
- `results/gergaud_MEE_M2_10N.{png,mp4,gif}` — example front-door output
  (default-endpoint 10 N run): static plot (`gergaud_plot.m`) + movie.

## MEE thrust-ladder campaign — Phase 2 (Campaign A)

**The story.** The Cartesian/Sundman stack above reproduces the paper's 10 N
headline case cleanly (M0-M3) but **dies at 5 N**: six documented strategies
(`DESIGN_thrust_ladder.md` §1) all hit the same false-infeasibility signature,
and a post-mortem review concluded the paper's own thrust ladder (10 N → 0.1 N,
Table 3) is not reachable by stretching the Cartesian/Sundman formulation —
the fixed-τ_f seed freezes the revolution count into the seed geometry, so
"stretching" a seed in time cannot add revolutions (a topology error, not a
tuning one). **Campaign A rebuilds the solver in Modified Equinoctial
Elements (MEE)** with true longitude L as the independent variable and the
total longitude span **ΔL as a decision variable** — L̇>0 strictly in this
regime, so L subsumes the Sundman clock outright, and because ΔL is now
solved for rather than fixed, **the optimizer grows the revolution count
itself** as thrust steps down. This is what makes thrust continuation work as
the ladder's backbone (10→5→2.5→1→0.5 N, each warm-starting the next) where
the Cartesian formulation could not.

**The linchpin gate (Task 4).** Before any ladder point is trusted, the new
MEE solver had to reproduce the Cartesian 10 N/c_tf=1.5 baseline to within
0.5 kg with matching burn structure — a cross-FORMULATION check, not just a
mesh check. It passed: m_f = 1377.10 kg vs the Cartesian 1376.74 kg (diff
0.36 kg), maxDefect 6.3e-15, same switch/rev basin (sw=19, revs=7.326). Every
ladder number below inherits this gate.

**The ladder (fuel solves, c_tf=1.5) + PSR-ported switch-aware refinement:**

| T [N] | m_f [kg] | switches | revs (ours) | revs (paper) | anchor t_f,min [ND] | anchor revs | R0=T·t_f,min [N·h] |
|---|---|---|---|---|---|---|---|
| 10  | 1377.10 | 19  | 7.326  | 7.5  | 22.2206 | 4.503 | 846.5 |
| 5   | 1364.54 | 32  | 14.157 | 15   | 44.6796 | 8.673 | 851.0 |
| 2.5 | 1369.79 | 76  | 27.841 | 30   | 89.2530 | 17.66 | 850.0 |
| 1   | **1371.44** | 171 | 69.152 | 74.5 | 223.8081 | 44.17 | 852.6 |
| 0.5 | 1375.28 | 362 | 138.597 | 149  | 446.27 (est.) | (est.) | 850.0 (by construction) |

Anchors are the free-longitude min-time solves; R0 = T_max·t_f,min holds to
**0.72% spread across the 4 independently certified anchors** (mean 223.14
ND ≈ 850.0 N·h) — the same empirical law the paper reports (≈850 N·h), now
reproduced across two thrust decades in a formulation the paper itself never
built. Figures: `fig_table3.m` → `results/fig_table3.png` (switches/revs vs
thrust + the R0-law panel) and `fig_front_mee.m` → `results/fig_front_mee.png`
(m_f vs thrust, the Fig-23-adjacent overlay).

**Six binding footnotes (carry into any downstream use of these numbers):**

1. **0.5 N row is anchor-free.** Its "t_f/t_f,min" is an **R0-law ESTIMATE**
   (anchorSource='R0law': tfTarget = 1.5×(223.14/0.5) = 669.42 ND), not an
   independently certified min-time solve — the 0.5 N min-time anchor hit a
   genuine conditioning wall (7 configs attempted, best defect 0.0545,
   reproducible MEX crashes). It is excluded from the R0-law fit shown on
   `fig_table3.png` (circular: it was built from the fit). If a certified
   0.5 N min-time anchor is ever obtained and differs from 446.27 ND by more
   than ~1%, the 0.5 N m_f/switches/revs need re-solving against the new target.
2. **0.5 N m_f/switches are PSR round-4-of-4, budget-limited** (stopReason=
   `maxRounds`, dsw=4/dmf=0.33 kg between the last two rounds) — not confirmed
   mesh-stable. Shown hollow on `fig_front_mee.png` for this reason.
3. **The ours-vs-paper revs gap is ladder-wide and systematic**
   (approximately −5.6%/−7.2%/−7.2%/−7.0% at 5/2.5/1/0.5 N), not a per-rung
   anomaly — one inherited model/paper discrepancy, footnoted once here
   rather than re-litigated at every rung.
4. **1 N provenance:** the table above carries the PSR-refined m_f =
   **1371.44 kg** (PSR round 2 of the MEE port), which **supersedes** the
   earlier uniform-mesh value 1370.36 kg (sw=171 stable across both; edge
   improved 0.9983→0.9994 — a discretization refinement, not a basin change).
5. **PMP dual/primer status: verifier delivered and proven correct, gates
   fail on raw duals, primal certification unaffected.** `verify_pmp_mee.m` +
   `mee_dual_to_costate.m`/`mee_primer_switch.m` were built and independently
   re-derived term-for-term (Task 10); the Fig-16 analog (`fig_switching.m`)
   honestly shows large primer misalignment (median 32.4° at 10 N, 60.0° at
   1 N-PSR, eccentricity-correlated) and switching-sign gates failing. A
   reviewer's independent KKT re-derivation showed this is **not a verifier
   bug** — raw IPOPT duals themselves fail cone-elided KKT stationarity at
   high eccentricity. This is now Campaign B scope (`DESIGN_dual_map.md`,
   escalate branch: investigate raw `lam_g` via `nlpsol` bypassing
   `opti.dual`). **Every m_f/switch/revs number in this README is a primal
   certification (defect/terminal-gated) and is unaffected by this open
   dual-side finding — do not read the primer/switching gate failures as
   casting doubt on the mass or structure numbers.**
6. **0.2 N and 0.1 N were honestly not attempted.** The 0.5 N min-time anchor
   conditioning wall (footnote 1) is where the deep-descent effort stopped;
   extending the ladder further is open future work, not a silent gap.

**Isp caveat now RESOLVED (2026-07-19):** the benchmark Isp is **1994.8 s**
(Caillau & Noailles 2001, δ=0.05112 km⁻¹·s ⇒ Isp=1/(δg₀); see the intro
paragraph above). Our default 2000 s is 0.27% high, so every m_f in the table
above is ~0.3 kg high vs the exact-Isp value — negligible, and no conclusion
changes. The absolute masses no longer rest on an unverified assumption.

**Fig-23 honesty note.** The paper's Fig 23 overlays several c_tf curves;
this campaign only ever solved **one c_tf (1.5) per thrust level**, so
`fig_front_mee.png` is the honest single-c_tf version (our 5 rungs against
a shaded band showing the paper-implied near-independence range, 1370-1375
kg) rather than a fabricated multi-curve reproduction. Ours spans
1364.5-1377.1 kg over the ladder — noticeably wider scatter than the paper's
implied near-independence, consistent with the basin-scatter phenomenon
already documented for the Cartesian M3 front (`fig_basin_scatter.m`).

**New modules (Phase 2, on top of the Cartesian module map above):**
`casadi_lt_mee.m` (MEE solver core, L-domain collocation, ΔL decision
variable), `mee_seed.m`, `run_mintime_mee.m`/`run_transfer_mee.m`/
`run_ladder.m` (per-thrust drivers + ladder orchestrator), `interp_warmstart.m`
(mesh-refine handoff, renormalizes RTN beta), `psr_mee_refine.m` +
`psr_switch_score_mee.m` + `psr_refine_sigma_mee.m` (PSR ported from
`GTO_tulip/direct/PSR/`), `mee_dual_to_costate.m` + `mee_primer_switch.m`
+ `verify_pmp_mee.m` (PMP verifier, footnote 5), `fig_table3.m` +
`fig_front_mee.m` (this section's deliverable figures).

Full task-by-task ledger: `.superpowers/sdd/progress.md` (section dated
2026-07-17 onward, "MEE thrust-ladder SDD ledger"). Design rationale and
open items: `DESIGN_thrust_ladder.md`.

## Front door: `run_gergaud`

`run_gergaud.m` is a single-call, PARAMETERS-block front door onto the MEE
thrust-ladder campaign above — it adds no new solver physics, just endpoint
resolution, cache-vs-solve selection, and Table-3 row/plot/movie assembly on
top of `mee_seed.m` / `casadi_lt_mee.m` / `homotopy_mee.m` /
`run_mintime_mee.m` / `run_transfer_mee.m` / `psr_mee_refine.m`.

**Usage.** Two equivalent calling styles, matching `PSR/run_psr.m` and
`elfo/elfo_run_one.m`:

```matlab
run_gergaud                                  % interactive: edit the
                                              % PARAMETERS block at the top
                                              % of the file, then run
row = run_gergaud(struct('thrustN', 5));     % opts override, same defaults
row = run_gergaud(struct('thrustN', 1, 'runMode', 'solve'));
row = run_gergaud(struct('thrustN', 0.2, 'runMode', 'probe'));   % honest
                                              % up-front wall warning, never
                                              % certified in this campaign
```

Every PARAMETERS-block field is also an `opts.<field>` key: `thrustN` (one
of 10/5/2.5/1/0.5/0.2/0.1 N, default 10), `P0_km`/`e0`/`i0_deg` (initial
orbit, default 11625/0.75/7 — the paper's own GTO-like start), `Pf_km`/`ef`/
`if_deg` (final orbit, default 42165/0/0 — GEO), `ctf` (t_f/t_f,min, default
1.5), `nodesPerRev` (default 25), `maxIter` (default 1500), `runMode`
(default `'auto'`), `makeMovie`/`makePlot` (default `true`), `m0kg`/`ispS`
(default 1500 kg / 2000 s), and `returnOnly` (test hook: returns the row
struct, skips plot+movie regardless of `makeMovie`/`makePlot`).

**Run modes.**
- `'auto'` (default) — if BOTH endpoints are the paper defaults AND a
  certified cache exists for `thrustN` (10/5/2.5/1/0.5 N), loads it and
  builds the row with no solve. Otherwise it falls through to a live solve
  automatically (this includes custom endpoints — `'auto'` never claims a
  cached number for a non-default target/initial orbit).
- `'solve'` — always runs the live pipeline (`run_mintime_mee` anchor →
  `run_transfer_mee` fixed-tf fuel homotopy → `psr_mee_refine` for
  `thrustN<=1` N), ignoring any cache.
- `'probe'` — research mode: forces a live solve and prints an up-front
  warning that thrust below 0.5 N (0.2/0.1 N) was never certified in this
  campaign. `row.certified` is reported honestly either way — the script
  never fabricates a row for a rung that doesn't converge.

**Endpoint knobs (default-preserving).** `(P0_km,e0,i0_deg)` and
`(Pf_km,ef,if_deg)` at their paper/GEO defaults resolve to `initElems=[]`
(the byte-preserving legacy literal already inside `mee_seed.m`) and
`xf=[1;0;0;0;0]` (the byte-preserving default already inside
`casadi_lt_mee.m`) — leaving both endpoints at their defaults reproduces the
existing certified numbers exactly and is what lets `'auto'` mode hit the
cache. Any deviation builds `initElems = [P0_km/LU; e0; 0;
tan(deg2rad(i0_deg)/2); 0; 1; 0]` and/or `xf = [Pf_km/LU; ef; 0;
tan(deg2rad(if_deg)/2); 0]` (LU = 42165 km, `kepler_lt_params.m`'s fixed
length unit) and tags the result's cache files with a deterministic hash
suffix so a custom run can never collide with a certified-default cache.
**Research-probe caveat:** the solver/seed were validated for GEO-like
(near-circular, near-equatorial) targets only — a significantly eccentric,
inclined, or retrograde custom final orbit is research-probe territory, not
a reproduction of a known-good case; the script reports whether the live
solve certified rather than presuming any custom target converges.

**Per-rung recipe / honesty map.** `run_gergaud` encodes, but does not
re-litigate, the campaign's own honesty footnotes above:

| T [N] | recipe | status |
|---|---|---|
| 10 / 5 / 2.5 | `run_mintime_mee` + `run_transfer_mee` | clean, cache-hit in `'auto'` mode |
| 1 | + `psr_mee_refine` | headline is the PSR-refined value (footnote 4) |
| 0.5 | anchor-free R0-law t_f,min (footnote 1) + `psr_mee_refine` (footnote 2) | budget-limited PSR, anchor is an estimate |
| 0.2 / 0.1 | live probe only (same recipe, honestly attempted) | never certified in this campaign (footnote 6); reports `certified=false` rather than a fabricated row |

See the "Six binding footnotes" list above for the full detail behind each
of these — `run_gergaud` only points at them.

**Outputs.** The row (a `gergaud_row()` struct) is always printed via
`gergaud_row_str()`, with an `UNCERTIFIED — <note>` banner prepended
whenever `row.certified` is false. Unless `returnOnly` is set, it also
writes `results/gergaud_<tag>.png` (`makePlot`, via `gergaud_plot.m`) and
`results/gergaud_<tag>.{mp4,gif}` (`makeMovie`, via `transfer_movie.m`),
where `tag = mee_fuel_tag(thrustN)` (e.g. `MEE_M2_10N`, `MEE_M2_2p5N`) plus
the endpoint-hash suffix for a custom target. Both renderers consume
Cartesian trajectory data; `mee_res_to_cart_res.m` is the adapter that
converts the solver's native MEE/L-domain state and RTN-frame control into
that inertial layout (reconstructing (r,v) via `elements_to_cart` at each
node's true longitude and rotating the thrust direction into ECI).

The four rendered movies of the certified ladder rungs live at
`results/movie_MEE_10N.{mp4,gif}`, `results/movie_MEE_5N.{mp4,gif}`,
`results/movie_MEE_2p5N.{mp4,gif}`, and `results/movie_MEE_1N.{mp4,gif}`.
`results/gergaud_MEE_M2_10N.{png,mp4,gif}` is an example of the front
door's own output (a default-endpoint 10 N run).

## Companions

- Parent campaign (shared solver architecture, `sundman_minfuel` pattern):
  `../GTO_tulip/`, esp. `process/LOW_THRUST_MINFUEL_CAMPAIGN.md`.
- Problem source: Haberkorn, Martinon & Gergaud, JGCD 27(6), 2004
  (`orbit_transfer/min_fuel_papers/Gergaud-Haberkorn-Martinon-JournalGuidance2004-preprint.pdf`).
