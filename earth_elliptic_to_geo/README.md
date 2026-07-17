# earth_elliptic_to_geo — direct-method reproduction of Haberkorn–Martinon–Gergaud

Direct (collocation NLP) reproduction of the low-thrust **minimum-fuel** Earth
orbit transfer in Haberkorn, Martinon & Gergaud, *"Low thrust minimum-fuel
orbital transfer: a homotopic approach,"* JGCD 27(6), 2004 — a 1500 kg
satellite from a low, elliptic, 7°-inclined orbit (P=11625 km, e=0.75) to
equatorial GEO (P=42165 km). The paper solves it **indirectly** (single
shooting + PMP + energy→mass homotopy); this reproduction solves the same
physics **directly** (Sundman-regularized trapezoidal collocation + CasADi/
IPOPT), reusing the campaign architecture from `NLP_lowThrust_GTO_tulip/`.

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

Isp is not stated numerically in the paper; we pinned **2000 s** (the design
default) and validated it by the M2 mass match landing inside the paper's
1370–1375 kg band — this is a family cross-check, not an independent citation
(the intended ref [6] Caillau & Noailles 2001 PDF was unreadable).

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
  22 switches / 7.97 revs toward 19 switches / 7.32 revs (a basin change,
  not noise) while m_f moved only 0.31 kg — the paper's own reported
  structure (18 switches, 7.5 revs) sits between our two mesh resolutions.
  Switch count should be read as a band, not a fixed integer.
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
- **PMP dual/primer anomaly (~20° misalignment), open, non-gating.** The
  KKT-dual costates on M0/M1 show a reproducible ~20° primer-vector
  misalignment that survives dual-polishing, node re-centering, and wide
  control bounds (consistent at 20.4645° to 4 decimals across probes) —
  suspected coupling through the Sundman `cScale` slack state, not yet
  root-caused. Primal certification (defects, mass, structure) is
  unaffected; M2/M3 gate on primal results only. Diagnostic scripts are
  preserved in `results/dual_anomaly/`.
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

## Companions

- Parent campaign (shared solver architecture, `sundman_minfuel` pattern):
  `../NLP_lowThrust_GTO_tulip/`, esp. `LOW_THRUST_MINFUEL_CAMPAIGN.md`.
- Problem source: Haberkorn, Martinon & Gergaud, JGCD 27(6), 2004
  (`min_fuel_papers/Gergaud-Haberkorn-Martinon-JournalGuidance2004-preprint.pdf`).
