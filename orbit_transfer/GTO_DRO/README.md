# GTO_DRO -- GTO -> DRO minimum-time low-thrust transfers

The Earth -> DRO leg the orbital network lacks (see `orbit_network/`),
solved in the CATALOG REGIME (150 kg, Isp 1710 s, thrust ladder 15 N down)
through the shared ladder engine in `costate_common/` -- the same
pipeline every costate catalog used: direct solve -> fly the control ->
`ms_tfmin` multiple shooting -> pumpkyn `tfMin` acceptance unchanged. No
solver code lives here. Design: `docs/superpowers/specs/2026-11-09-gto-dro-mintime-design.md`.

Stage 1 (this folder today): ONE transfer with its parameters exposed.
Stage 2 (planned): the catalog chain script over 16 sheets. Stage 3
(planned): the GTO catalogs enter the orbital network.

## Run

```matlab
cd orbit_transfer/GTO_DRO
gto_dro_study              % one transfer: edit the knobs block, press Run (~minutes)
cd tests; run_tests        % clearance instrument, driver planning, the study headless
```

Paths come from `setup_gto_dro_paths` (this dir, lib/, costate_common,
oclib, casadi, pumpkyn + pumpkynPie). Outputs go to `results/<run>/`
(git-ignored): three PNGs and `gto_dro_study.mat` with `run_meta`.

## Knobs (section 0 of the study script)

| knob | default | meaning |
|---|---|---|
| `thrustN` | 15 | returned rung [N]; the ladder walks 15 N down to it; [1, 15] only |
| `m0kg`, `ispS` | 150, 1710 | spacecraft |
| `gto.orientDeg` | 0 | Earth->Moon line to perigee [deg], rotating sense |
| `gto.sma_km`, `gto.ecc` | 24446, 0.7248 | the 350 x 35786 km GTO |
| `gto.depFrac` | 0.10 | departure phase, TIME fraction since perigee (0.5 = apogee); 0.10 by measurement, see Findings |
| `dro.tau` | 1 | DRO period [ND]; 0.5, 1, 2, 3 are the catalog periods (1 = 4.43 d) |
| `dro.arrFrac` | 0 | arrival phase along the DRO from its +x crossing |
| `recipe.tf0` | 0.30 | top-rung t_f guess [ND], inherited from the tulip campaign; the script prints t_f / tf0 |
| `recipe.floorEarthKm` | 200 | Earth altitude floor imposed in the direct solve [km] (engine option added 2026-11-09) |
| `out.reuse` | '' | gate a saved `gto_dro_study.mat` instead of solving (knobs must match); the reuse run saves to `gto_dro_study_reuse.mat`, never over the original |
| `out.conjTest` | true | run the conjugate re-solve (S4); false = NOT APPLICABLE |

`study_opts` in the workspace overrides any knob (the tests use it).

## What the study script checks

| id | condition |
|---|---|
| V1 | GTO perigee and apogee radii match a(1 -+ e) to 1 km |
| V2 | DRO closes to 1e-7 ND after cont_np |
| V3 | t_f ratio between consecutive rungs in [1/3, 3] (steps only) |
| N0 | every walked rung accepted |
| N1 | multiple-shooting residual < 1e-8 |
| N2 | max |H| along the flown arc < 1e-6 (autonomous, free time) |
| N3 | flown arrival within 100 km and 10 m/s |
| N4 | |lambda_m(t_f)| < 1e-6 |
| N5 | adjoint equations hold (finite differences of H) |
| N6 | minimum principle for the applied control, throttle 1, Q_mt >= 0 |
| N7a | Earth altitude >= 200 km on the flown arc (post-solve check) |
| N7b | lunar altitude >= 500 km on the flown arc |
| S1-S3 | strong Legendre, strict bang, dim S = 1 (hypothesis gates + lift margin) |
| S4 | no conjugate point: junction test and dense spectrum scan; UNRESOLVED when either cannot decide |
| X1 | pumpkyn tfMin at the accepted z8 returns within 1e-6 (NOT CHECKED if it throws) |
| X2 | pumpkyn's field agrees with an independent field row by row |

V, N and X are code gates: the script throws when one fails. S1-S4 are
results and are reported either way (strict mode throws on them too).

## Findings (stage 1, 2026-11-09)

- The perigee departure (phase 0) at 15 N does not converge. Measured 2026-11-09 before the Earth floor: IPOPT hits the 6000-iteration cap with the defect at 1e-15 while t_f creeps near 0.41 ND and the all-burn arc ends with about 5 percent of the mass (15 N on 150 kg burns 77 kg per day). Measured 2026-11-09 WITH the 200 km floor: still Maximum_Iterations_Exceeded at 6000 iterations, defect 1.1e-08, t_f 0.41701 ND (1.848 d), the Earth altitude riding the 200 km floor, m_f/m0 0.0476, 516 s. Phase 0.10 converges in about 300 iterations (about 100 s with the floor). The study's default departure phase is therefore 0.10; the perigee cell belongs to stage 2 (phase continuation).
- Without an Earth floor the direct solver returned a defect-clean arc passing 1722 km from Earth's CENTRE (inside the planet) whose flown control missed the DRO by 880,000 km. The shared engine now has an Earth-altitude floor option (`casadi_mintime_dro` `minAltEarthKm`, `thrust_ladder_library` `floorEarthKm`, default off, byte-identical for every other campaign); this driver sets it to 200 km.
- With the floor, the cold start from tf0 = 0.30 lands on t_f 0.3040 ND (1.348 d), m_f 45.84 kg, dV 19.88 km/s, flown miss 1e-4 km, pumpkyn tfMin acceptance |dz| 7e-14, lowest Earth altitude 13,888 km at departure. The tf0 probes 0.24, 0.26 and 0.27 ended at the same t_f with status Solved_To_Acceptable_Level (they were not flown or accepted). A second floor-satisfying solution exists at t_f 0.2668 ND (reached once by an unfloored cold start, warm-startable, never reached by a floored cold start).
- The 0.3040 ND solution passes every necessary condition and S1 to S3, but the conjugate test REFUTES it: a conjugate point near t/t_f 0.71 where |lambda_v| falls to 6e-4 (junction test and dense scan agree). It is a certified extremal, not a certified local minimum. The verdict line says so.
- Stage 2's pilot therefore owns basin selection at the top rung (a direction multistart or a warm start from the 0.2668 solution) and the conjugate test on whichever basin it picks; the 0.2668 ND solution is a KKT point of the direct NLP recorded as a number here; it was never flown, refined or accepted.
- The 12 N rung warm-started from the 15 N solution (tests/test_gto_dro_study) lands at t_f 0.3070 ND (1.361 d), m_f 65.86 kg, dV 13.80 km/s, flown miss 1e-4 km, and PASSES the conjugate test (junction test and dense scan, 0 zeros). Scaling the two 15 N basins by the ladder's own rule t_f ~ T^-0.6 gives 0.305 ND for the 0.2668 branch and 0.347 ND for the 0.3040 branch, so the 12 N solution is consistent with the continuation of the faster branch. Measured 2026-11-09 (reverse warm start): 15 N cold 0.30404 ND, 12 N warm 0.30700 ND, then 15 N warm-started FROM the 12 N solution 0.26682 ND (Solve_Succeeded) -- it does not return to 0.3040, so the 12 N solution continues the faster 15 N branch. The V3 ratio 1.010 (printed 1.010 against 3) is far below the 1.14 the scaling predicts, which is the same statement. So the 12 N transfer is a certified local minimum-time extremal with every gate PASS in the recorded test run (the test itself requires only the V, N and X gates to pass), while the 15 N cold start is not.

## What it does not claim

- No new solver and no change to the physics of `costate_common`; a violated Earth floor is now a CONSTRAINT in the direct solve (default 200 km) and still a post-solve CHECK (N7a) on the flown arc, since the constraint binds at the nodes only.
- No free arrival or departure phase: rendezvous at the chosen states.
- No thrust below 1 N (many-rev regime: junction states, Sundman mesh).
- No sheet writing; the catalog is stage 2.
- A full PASS table says "a certified local minimum-time extremal that
  pumpkyn's tfMin accepts unchanged", never global optimality.
- The GTO is an algebraic Kepler locus frozen at epoch zero in the
  rotating frame: a set of departure states keyed by orientation and
  time-since-perigee, not a trajectory.

## Stage 2: the catalog (built 2026-11-10; pilot and fleet status below)

| file | role |
|---|---|
| `catalog/run_gto_dro_catalog.m` | fleet front door: 16 sheets (orientDeg {0,90,180,270} x DRO tau {0.5,1,2,3}), 12 x 6 phasing, rungs [15 12 10 7 5], Earth floor 200 km, two-point top-rung multistart; `planOnly` returns the engine options without solving |
| `catalog/run_gto_dro_batched.sh` | nohup'd OS-timeout driver (matlab-campaign discipline) |
| `lib/earth_clearance_sweep.m` | flown Earth altitude of every accepted entry, sidecar + writeback (demotes below the floor) |
| `lib/seed_top_from_entry.m` | a flown entry as the engine's top-rung seed (reverse warm start, phase continuation) |
| `lib/repair_gto_dro_cells.m` | 'unsolved' (phase continuation from the nearest solved departure phase) and 'basin' (S4 junction test on every top rung; refuted ones re-solved from their 12 N entry; backup + restore) |
| `catalog/build_gto_dro_catalog.m` | the post-fleet chain: census, earth, repairs, earth re-sweep of repaired sheets, package (clean staging copy of the official sheets), conj, gates, audit, pictures, deliverable (refused on blockers; `allowConjBlockers` waives only refuted-entry blockers) |
| `catalog/audit_gto_dro_entries.m` | up to 3-entry replay: flown miss, tfMin acceptance, Earth floor (gated); the tulip pick rule when the sheets allow it, a reduced rule otherwise |
| `catalog/costate_catalog_pick.m`, `costate_lib_describe.m` | vendored pickers (deliverable convention); the describe copy prints a DRO arrival key instead of Np |

The fleet front door and the repair tool run with `directGateKm = inf` (no direct gate; the refined gate stays 100 km, with tfMin acceptance and the ms hard timeout) because the open-loop flight of the N = 400 direct control is not a filter on these transfers. Measured 2026-11-10: with `directGateKm = 1000` the pilot sheet solved, per rung [15 12 10 7 5] N, only 57 / 50 / 33 / 26 / 10 of 72 cells, with 54 ladders stopped at the direct gate; a probe re-solving six gate-stopped cells with no direct gate went from 14 to 26 of 30 OK rungs, refinement closing direct misses of 1,778 km and 42,964 km to under 0.1 km, and a failed rung no longer stops the ladder. Four rungs stayed open (126, 488, 1,082, 30,285 km): the size of the direct miss does not predict whether refinement succeeds.

Run the pilot (sheet 2), then the chain on it, then the fleet:

    cd orbit_transfer/GTO_DRO/catalog
    nohup ./run_gto_dro_batched.sh 2 6 1800 > /dev/null 2>&1 &
    # MATLAB, after the pilot's completion marker:
    build_gto_dro_catalog(struct('expectSheets', 1, 'deliverable', false))
    nohup ./run_gto_dro_batched.sh "[9:12 1 3:8 13:16]" 6 1800 > /dev/null 2>&1 &

Chain safeguards (fix round 1): a call with a sheet-writing stage on takes the
fleet driver lock (`results/.gtodrocat_driver.lock`) and refuses with
`build_gto_dro_catalog:fleetRunning` if it exists (`ignoreLock` skips it). The
deliverable stage runs a preflight on what is on disk whatever stages ran:
catalog present, packaged entries equal to the sheets' (`stale:`), conj and
hypothesis-gate records complete and a current Earth sweep (`unverified:`),
and an audit file whose fingerprint still matches the sheets; none of these
can be waived. `allowConjBlockers` and `allowEarthDemotions` waive refuted
entries (`conj:`) and Earth demotions (`earth-demoted:`) from the ship
decision only; both stay in the returned blocker list.

Chain notes: the packager is fed `results/package_stage/` (the official
sheets only), because the catalog folder also holds `_earth`, `_repairs`,
`_repair_backup_*`, `_deep` and `_direct` files that match the sheet glob. Re-packaging
deletes the conj and gates sidecars (their entry indices go stale).
Engine options added for this stage (costate_common, default off): `seedTop`,
`redo`, `directGateKm`, `msHardCapSec`. The front door passes `msHardCapSec =
300` (override `optsIn.msHardCapSec`, 0 = off), opens the fence pool once per
call, writes `ms fence ACTIVE (pool of N)` (or `UNAVAILABLE`) to
`combined_progress.txt`, refuses to run unfenced unless `optsIn.requireFence
= false` (`run_gto_dro_catalog:noFence`), and silences the ode45 tolerance
warning while the cap is on, because the pilot's cell (8,4) of
`gto_dro_o000_tau1` hung 55 min inside `ms_tfmin` at 15 N (an iterate's
junction mass runs to zero mid-segment and pumpkyn's `tfMinProp` restarts the
failed `ode45` forever). `mintime_prop_seg` now refuses such a segment at the
source (`massDepleted`): reproduced with the guard and a 60 s cap as safety net, the
15 N rung of that cell now fails refinement in-process (37 s for the whole rung,
direct solve included; the cap never fired), and the 12 N rung is accepted
(t_f 0.36518 ND, flown 4e-4 km). The batch driver trims `gtodrocat_matlab.log`
to its first and last 2 MB whenever a pass leaves it above 200 MB.
The fleet also keeps each rung's direct solution in `<tag>_direct.mat` next
to the sheet (engine `keepDirect`, front-door knob default true, about 12 MB
per sheet), because 46 of the pilot's 280 refined entries fly below the Earth
floor while the direct solution measured on cell (11,5) holds it, and the
sheets alone cannot answer that question without re-solving; a repair that
restores a cell marks its sidecar slices not recorded.

## Pilot findings (stage 2, 2026-11-11; sheet 0 deg / tau 1, `gto_dro_o000_tau1`)

Recipe: no direct gate, refined flown gate 100 km, tfMin acceptance 1e-6, ms
hard cap 300 s. 13 batch passes, 4 h 54 min, no hang, the cap never fired.
Each line names the stage that measured it.

- **Fleet engine (raw sheet, kept in `results/pilot_attempt3_raw/`).** Cells
  solved per rung [15 12 10 7 5] N: 56 / 58 / 55 / 57 / 54 of 72 (280
  entries; the 1000 km direct gate gave 57 / 50 / 33 / 26 / 10). 34 cells hold
  the full ladder; 8 cells hold nothing (both attempts spent), four of them
  at arrival phase 1. A failed rung no longer ends a ladder.
- **tf0 calibration (raw sheet).** Accepted 15 N t_f: median 0.3098 ND, range
  0.2589 to 0.3238 (n = 56). The median is within 0.02 of 0.30 (the plan's
  rule), so the fleet keeps `tf0 = 0.30`.
- **Earth sweep (chain stage 2).** 46 of the 280 entries fly below the 200 km
  Earth floor and are demoted; 41 of them pass below the surface (deepest
  -3,855 km). Per rung: 5 / 5 / 12 / 17 / 7. The plan expected none. 43 of
  the 46 lie in a band of departure phases that moves with thrust (phases
  11-12 of 12 at 15 and 12 N, 10-11 at 10 N, 8-10 at 7 N, 3-5 at 5 N), and the
  dip is at the first perigee pass (t / t_f 0.03 to 0.12). Re-flown and
  sampled at the 401 uniform node times, 45 of the 46 are below the floor at
  a node, and 43 stay below it longer than one node spacing (up to 3.8): the
  refined trajectory is not the direct solution, which holds the floor at
  its nodes. The refinement (`ms_tfmin`) carries no Earth constraint.
  Measured on one demoted cell, (11,5) at 15 N, re-solved cold: the direct
  solution touches the floor (243.9 km at the nodes, 199.9 km at a
  Hermite-Simpson midpoint), its control flown open loop stays within 7 km
  of the node states through perigee and bottoms at 195.5 km, and the
  refined entry of the same cell bottoms at -1,541 km for a t_f shorter by
  0.3 percent. So there the floor constrains the optimum and the
  unconstrained refinement leaves it. A demotion condemns that trajectory,
  not the cell: the fleet keeps every direct solution (`<tag>_direct.mat`)
  so floor-active cells can be re-seeded or shipped as constrained direct
  solutions afterwards. The same sweep on the 1000 km-gate sheet finds 6 of 176 below
  the floor, including the same two 15 N entries at -1,540 and -1,564 km, so
  the 15 N violations do not come from dropping the direct gate.
- **Unsolved cells (chain stage 3, phase continuation).** 8 of 8 repaired
  (15 N t_f 0.2616 to 0.3274). The perigee row (departure phase 0), where the
  engine had retired 2 of 6 cells, is solved at every rung afterwards (30 of
  30 entries).
- **Top-rung basin census (chain stage 4).** 18 refuted 15 N entries
  re-solved from their 12 N entry: 10 landed on a faster branch (t_f shorter
  by 0.034 to 0.040 ND, about 12 percent) and all 10 pass the conjugate test
  afterwards; 5 returned the same solution and stay refuted; 3 were restored
  (2 not faster, 1 faster but lost lower rungs). 15 N t_f after repairs:
  median 0.3024, range 0.2589 to 0.3274 (n = 59).
- **After the chain.** 59 / 60 / 51 / 47 / 54 of 72 per rung, 271 entries,
  minimum flown Earth altitude 240 km.
- **Conjugate sweep (stage 6).** 211 pass, 60 refuted; per rung 10 of 59, 7
  of 60, 11 of 51, 12 of 47, 20 of 54. Only top rungs are repaired, so the
  lower-rung refutations are a census, not a repair queue.
- **Hypothesis gates (stage 7).** 271 of 271 (H2, H3, no abnormal lift).
  **Audit (stage 8, reduced rule).** fly 3 / 3, accept 3 / 3, Earth 3 / 3.
- **Blockers returned.** `earth-demoted` (46) and `conj` (60), both of the
  waivable kind; nothing stale or unverified.
- **Deep probe (six full-ladder cells, 9 rungs, separate sheet in
  `results/deep_probe/`, cold start, not Earth-swept).** Solved per rung
  [15 12 10 7 5 3 2 1.5 1] N: 5 5 5 5 5 5 1 0 0 of 6. One cell failed its
  15 N cold start although the pilot had solved it. All five others reach
  3 N; one reaches 2 N (t_f 1.027 ND); below that the refined flight misses
  by 300,000 to 425,000 km or the direct solve fails. The wall on this sheet
  is between 3 N and 2 N.
- **Margins on the accepted entries (sheet after the chain, 271).** Flown
  miss at most 0.23 km (gate 100 km), tfMin |dz| at most 7.7e-9 (gate 1e-6),
  refinement residual at most 8.0e-11. Earth and Moon clearance re-measured
  on a dense spline of every flight: Earth 240.4 km, within 0.1 km of the
  sample minimum on every entry; Moon 533.5 km, none below the 500 km floor
  (the chain does not sweep the Moon yet).
- **Review.** `reviews/astra_stage2_2026-11-11.md` and its dispositions.
- **Not answered here.** The 180-degree orientation; the fleet runs those
  sheets first (`"[9:12 1 3:8 13:16]"`).

## Fleet findings (stage 2, 2026-11-15; all 16 sheets)

Fleet: 141 batch passes, 78.5 h (2026-11-11 04:16 to 2026-11-14 10:49), the
solver fence fired 5 times, nothing hung. Chain on the 16 sheets: 13 h 47 min.
The raw sheets are kept in `results/fleet_raw/`. The four counts Astra asked
to be kept apart are reported apart: direct solutions on record, accepted by
the fleet, accepted after the chain (floor-clear), and not refuted. Each cell
of the table is out of 360 (72 cells x 5 rungs).

| GTO orientation | DRO size | direct on record | accepted (fleet) | after the chain | not refuted | refuted | primer flag |
|---|---|---|---|---|---|---|---|
| 0 | 0.5 | 323 | 240 | 262 | 189 | 73 | 3 |
| 0 | 1 (pilot) | not kept | 280 | 271 | 211 | 60 | 0 |
| 0 | 2 | 321 | 262 | 275 | 209 | 66 | 1 |
| 0 | 3 | 343 | 280 | 289 | 260 | 29 | 0 |
| 90 | 0.5 | 343 | 221 | 187 | 185 | 2 | 1 |
| 90 | 1 | 337 | 229 | 199 | 193 | 6 | 4 |
| 90 | 2 | 313 | 222 | 180 | 174 | 6 | 2 |
| 90 | 3 | 344 | 249 | 214 | 209 | 5 | 2 |
| 180 | 0.5 | 252 | 116 | 130 | 130 | 0 | 4 |
| 180 | 1 | 273 | 114 | 126 | 126 | 0 | 1 |
| 180 | 2 | 271 | 131 | 128 | 127 | 1 | 1 |
| 180 | 3 | 268 | 137 | 131 | 131 | 0 | 1 |
| 270 | 0.5 | 347 | 251 | 215 | 199 | 16 | 0 |
| 270 | 1 | 344 | 255 | 210 | 194 | 16 | 1 |
| 270 | 2 | 353 | 268 | 225 | 207 | 18 | 1 |
| 270 | 3 | 336 | 243 | 223 | 210 | 13 | 1 |
| all | | 4,768 of 5,400 | 3,498 | 3,265 | 2,954 | 311 | 23 |

"Direct on record" is the sidecar after the chain (a restored repair marks
its slices not recorded), so it slightly undercounts what the fleet solved.

- **The catalog.** 3,265 entries of 5,760 possible (57 percent); per rung
  [15 12 10 7 5] N: 706 / 681 / 625 / 612 / 641. 2,954 are not refuted by the
  conjugate test. Coverage is by orientation, not by DRO size: the four DRO
  sizes of one orientation differ by at most 34 entries.
- **Coverage by orientation.** Cells holding at least one entry, of 72 per
  sheet: 0 deg 69 to 72, 90 deg 62 to 68, 270 deg 68 to 70, 180 deg 48 to 50.
  Cells holding the full ladder: 0 deg 18 to 27, 270 deg 8 to 14, 90 deg 4
  to 11, 180 deg 1 to 5.
- **The 180-degree corner.** The fleet accepted about a third of the
  entries there (114 to 137 of 360) and the repairs did not change that:
  phase continuation recovered 25 of 98 unsolved cells on these four sheets
  (73 of 99 on the other eleven fleet sheets). The entries that exist are
  the cleanest in the catalog: 1 refuted of 515. So the 180-degree
  orientation is hard to solve from this recipe, not poor in good
  extremals. The fastest transfers of the whole catalog are there (15 N
  t_f 0.2210 ND, 0.96 d, on the size-3 sheet).
- **Earth floor (chain stages 2 and 2b).** 707 entries fly below 200 km
  after refinement and are demoted, on every sheet (16 to 64 per sheet); the
  deepest passes 5,760 km below the surface. By orientation and rung:

  | orientation | 15 N | 12 N | 10 N | 7 N | 5 N |
  |---|---|---|---|---|---|
  | 0 | 21 | 20 | 34 | 63 | 26 |
  | 90 | 33 | 46 | 56 | 49 | 29 |
  | 180 | 3 | 16 | 23 | 23 | 20 |
  | 270 | 68 | 66 | 60 | 45 | 6 |

  The pilot's one-cell measurement holds at fleet scale. Of the 661 demoted
  entries on the 15 sheets that kept their direct solutions, 640 (97
  percent) have a direct solution whose minimum altitude is within 1 km of
  the 200 km floor; 12 have no direct solution on record and 9 have one
  clear of the floor. So nearly every demotion is a cell where the floor
  constrains the optimum and the unconstrained refinement leaves it. Those
  640 direct solutions are on disk (`<tag>_direct.mat`). Lowest accepted
  entry: 201.4 km.
- **Unsolved cells (stage 3, phase continuation).** 98 of 197 repaired.
- **Top-rung basin census (stage 4).** 54 refuted 15 N entries had a 12 N
  entry to re-solve from: 41 landed on a faster branch, 13 were restored (11
  not faster, 2 faster but lost lower rungs, the case of disposition C5).
  All 54 are at orientation 0 or 90. 60 top-rung entries remain refuted (35
  at 0 deg, 2 at 90, 23 at 270).
- **Conjugate sweep (stage 6).** 2,954 pass, 311 refuted, none undecided.
  Refuted share per rung: 8 / 10 / 8 / 10 / 12 percent. 228 of the 311 are
  at orientation 0 and 63 at 270; orientations 90 and 180 hold 20 between
  them. These are different populations per rung, so no trend down the
  ladder is claimed.
- **Hypothesis gates (stage 7).** 3,265 gated; no abnormal lift, no H3
  failure; 23 entries fail H2 (minimum |lambda_v| on the arc below 1e-6;
  values 8e-8 to 9.8e-7). The minimum is a continuous quantity across the
  catalog (median 2.2e-4, 98 more entries between 1e-6 and 1e-5), so the 23
  are the tail under a fixed threshold: arcs where the thrust direction is
  close to a reversal. None of the 23 is refuted by the conjugate test. The
  chain treats this blocker as never waivable.
- **Margins on the 3,265.** Refinement residual at most 1.0e-10 (median
  6.3e-14); tfMin |dz| at most 1.5e-7 (gate 1e-6); flown miss median 0.0014
  km, 28 entries above 1 km, largest 33.1 km (90 deg / size 3, cell (1,4),
  5 N; gate 100 km).
- **Perigee row (departure phase 0), entries of 24 per orientation and
  rung.** 0 deg 20 to 21; 90 deg 23 to 24 down to 7 N and 20 at 5 N; 180
  deg 18 to 19 down to 7 N and 8 at 5 N; 270 deg 23 down to 7 N and none at
  5 N.
- **15 N t_f per sheet.** Medians 0.273 to 0.303 ND (1.19 to 1.32 d);
  catalog range 0.2210 to 0.3548 ND (0.96 to 1.54 d).
- **Audit (stage 8, reduced rule).** fly 3 / 3, accept 3 / 3, Earth 3 / 3.
- **Blockers returned.** `earth-demoted` (16, one per sheet, waivable),
  `conj` (311 refuted, waivable), `gates` (23, not waivable). The
  deliverable is not built.
- **Not done yet.** The Moon floor is not swept at catalog scale (the pilot
  measurement found none below 500 km); lower-rung refutations are a census,
  not a repair queue; the floor-active cells have no admissible entry; the
  package field `tauDRO` holds the GTO orientation (the DRO size is
  `arr_params.tau`).
