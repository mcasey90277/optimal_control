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
| `out.reuse` | '' | gate a saved `gto_dro_study.mat` instead of solving (knobs must match) |
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

- The perigee departure (phase 0) at 15 N does not converge: IPOPT hits the 6000-iteration cap with the defect at 1e-15 while t_f creeps near 0.41 ND and the all-burn arc ends with about 5 percent of the mass (15 N on 150 kg burns 77 kg per day). Phase 0.10 converges in about 300 iterations. The study's default departure phase is therefore 0.10; the perigee cell belongs to stage 2 (phase continuation).
- Without an Earth floor the direct solver returned a defect-clean arc passing 1722 km from Earth's CENTRE (inside the planet) whose flown control missed the DRO by 880,000 km. The shared engine now has an Earth-altitude floor option (`casadi_mintime_dro` `minAltEarthKm`, `thrust_ladder_library` `floorEarthKm`, default off, byte-identical for every other campaign); this driver sets it to 200 km.
- With the floor, every cold start from tf0 in {0.24, 0.26, 0.27, 0.30} lands on the same solution: t_f 0.3040 ND (1.348 d), m_f 45.84 kg, dV 19.88 km/s, flown miss 1e-4 km, pumpkyn tfMin acceptance |dz| 7e-14, lowest Earth altitude 13,888 km at departure. A second floor-satisfying solution exists at t_f 0.2668 ND (reached once by an unfloored cold start, warm-startable, never reached by a floored cold start).
- The 0.3040 ND solution passes every necessary condition and S1 to S3, but the conjugate test REFUTES it: a conjugate point near t/t_f 0.71 where |lambda_v| falls to 6e-4 (junction test and dense scan agree). It is a certified extremal, not a certified local minimum. The verdict line says so.
- Stage 2's pilot therefore owns basin selection at the top rung (a direction multistart or a warm start from the 0.2668 solution) and the conjugate test on whichever basin it picks; both basins are recorded here.
- The 12 N rung warm-started from the 15 N solution (tests/test_gto_dro_study) lands at t_f 0.3070 ND (1.361 d), m_f 65.86 kg, dV 13.80 km/s, flown miss 1e-4 km, and PASSES the conjugate test (junction test and dense scan, 0 zeros). Scaling the two 15 N basins by the ladder's own rule t_f ~ T^-0.6 gives 0.305 ND for the 0.2668 branch and 0.347 ND for the 0.3040 branch, so the 12 N solution is consistent with the continuation of the faster branch: the first ladder step moved off the refuted branch. The V3 ratio 1.010 (printed 1.010 against 3) is far below the 1.14 the scaling predicts, which is the same statement. So the 12 N transfer is a certified local minimum-time extremal with every gate PASS in the recorded test run (the test itself requires only the V, N and X gates to pass), while the 15 N cold start is not.

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
