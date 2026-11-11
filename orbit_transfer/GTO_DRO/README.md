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
    nohup ./run_gto_dro_batched.sh "[1 3:16]" 6 1800 "" <tf0 from the pilot> > /dev/null 2>&1 &

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
`_repair_backup_*` and `_deep` files that match the sheet glob. Re-packaging
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
