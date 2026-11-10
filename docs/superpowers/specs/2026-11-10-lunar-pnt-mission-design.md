# Lunar PNT mission simulation (DRO -> six-satellite tulip constellation -> a month of south-pole GDOP): design

Date: 2026-11-10. Design approved by Mike 2026-11-10 (five sections presented in the brainstorm session; decisions 1 to 4 in section 2 are his, made by question and answer).

## 1. Purpose

Show, end to end, how a lunar PNT constellation gets put in place: six spacecraft parked in one distant retrograde orbit (DRO) each fly a certified minimum-time low-thrust transfer into the tulip orbit, timed so that they land on the six optimally phased constellation slots, and the assembled constellation is then scored for a month of lunar south-pole navigation geometry (GDOP) against the same LunaNet requirements the tulip paper used.

Nothing in the repository composes transfers with constellation geometry today. The transfer side (the 70 mN DRO -> tulip costate library, the tfMin acceptance instruments) lives in `optimal_control/orbit_transfer`; the geometry side (phasing, GDOP grid, R1 to R3 scorers) lives in `proj7`. This study is the seam between them, and its product is a mission-design demonstration: a schedule, six certified transfers, the position histories, the GDOP series with the arrivals marked, and the requirement verdicts.

Success is: one study script that runs from a fresh MATLAB session, prints PASS/FAIL for every condition it relies on, and produces the figures; a validation gate that proves the GDOP half reproduces the paper's numbers on the paper's phasing; every transfer accepted by the same instruments that certify the library.

### Goals

- A schedule that assembles the constellation as early as possible (the last arrival is minimised), with propellant reported per spacecraft.
- Each spacecraft arriving on its exact slot state, not merely near it, by an accepted minimum-time extremal.
- Thirty days of south-pole GDOP for the assembled constellation, plus the ramp-up window while it assembles, scored with proj7's own engine and requirement scorers.
- A result a reader can check: every gate computed in the script, every number traceable to a named function.

### Non-goals (recorded refusals, section 9)

- No ephemeris model: CR3BP throughout, as in the library and the paper.
- No station-keeping, no navigation filter, no signal model: GDOP is geometry only.
- No propellant-optimal schedule: propellant is reported, not optimised.
- No GTO leg: the GTO -> DRO catalog is a different regime (1710 s, 15 N to 5 N) and is not composed here.
- Nothing in `proj7`, `pumpkyn`, `pumpkynPie` or the deliverable-9 kit is modified.

### Anti-goals

- A schedule whose transfers are interpolated or unpolished catalog entries presented as rendezvous. A cell's arrival phase is within a tolerance of a slot, not on it; only the polished transfer counts.
- A GDOP number that has not been checked against the paper's engine on the paper's phasing.
- A stale phasing vector. `run_gdop_analysis_6sat.m` and `opt_GDOP_Tulip_Demo.m` carry old `dtau` values; the only source is `tulip_phasing()`.

## 2. Decisions (Mike, 2026-11-10)

1. **Regime.** The 70 mN / Isp 900 s / 150 kg DRO -> tulip phase library of record: `orbit_transfer/DRO_tulip/deliverables/costate_library_70mN_user_deliverable9/costate_catalog_dro_tulip_70mN.mat`, 24 departure phases x 48 arrival phases, 1,152 of 1,152 cells certified (status 4), t_f 16 to 23 d, dV 0.68 to 0.99 km/s. Kit functions: `transfers_at(c, sD, sA)` selects, `catalog_transfer(cat_, entry, opts)` flies an entry.
2. **DRO start.** Six slots evenly spread in phase: spacecraft k starts at DRO phase (k-1)/6 of the DRO period at the clock epoch (pumpkynPie `phaseConstellation` convention). DRO: tau 1 ND = 4.433 d.
3. **Schedule objective.** Earliest complete constellation: minimise the LAST arrival over the 720 assignments of spacecraft to tulip slots, the catalog cell each spacecraft uses, and a whole number of DRO revolutions of waiting before each departure. Propellant is reported.
4. **Rendezvous, approach B.** Enumerate on catalog cells with a residual tulip-phase tolerance of 0.005 of a tulip period (about 2.8 h), then polish each spacecraft to the EXACT slot state with pumpkyn `tfMin` seeded by the cell's z8, with the arrival epoch pinned and the departure state updated until t_f converges; then the stage 1 instruments on each polished transfer (flown miss < 100 km, pointwise PMP checks, min-time hypothesis gates, conjugate test).
5. **Scheduler candidates (2026-11-10, spec review).** Primaries only: the one certified status-4 transfer per (departure, arrival) cell, 1,152 candidates. The 656 alternatives are not scheduled; the study reports as a sensitivity whether any status-4 alternative would have given an earlier last arrival.
6. **GDOP grid (2026-11-10, spec review).** The paper's Table 1 grid, i.e. `run_requirements_table`'s: the five latitude chunks `run_requirements_table` passes (`lat_range` pairs from `latv = [-89.9, -89:0.5:-75]` in groups of six with `lat_res` 0.5, so the first chunk evaluates -89.9:0.5:-87.4; 30 latitudes), longitude 0 to 358 every 1 degree (10,770 points), elevation mask 0 degrees, Nt = 5760 per tulip period (about 348 s). The validation gate and the 30-day R1 to R3 table both use it. A coarse grid (latitude -89.5 to -75 in 2 degree steps plus the pole, longitude 0 to 358) is used only for time-series figures and the assembly window.

Constellation facts the design rests on (gathered 2026-11-10, to be re-verified by the plan's first task):

- Phasing: `[dtau6, dtau5, Np, tau0, pm] = tulip_phasing()` in `proj7/pipelines/gdop/`; six-satellite `dtau6 = [0.7783 1.1456 2.5174 2.8930 4.2668 4.6316]` (dimensionless time offsets in [0, tau0), tau0 = 5.236; Koblick's GA optimum, adopted 2026-07-06).
- Tulip: Np 7, pm -1, period 23.209 d. Tulip phase fraction x tau0 equals dtau.
- pumpkyn `tulipConstellation` convention: satellite j sits at orbit time `mod(tau + dtau_j, tau0)`.
- One orbit definition for everything this study flies or targets: `get_family_orbit` (continued with `cont_np`) evaluated through the periodic `phase_state`, which is what the library's endpoints use. `tulipConstellation` builds its orbit without the continuation step and differs by up to 6.7e-6 ND (about 2.6 km); that is the paper's reference geometry and the validation gate's tolerance absorbs it.
- Catalog grids: departure phases (k-1)/24; arrival phases 0.0129 + (k-1)/48 (not zero-based).
- The catalog keys only phases; it carries no epoch. The mission clock is this study's own.

## 3. Folder

```
orbit_transfer/lunar_pnt_mission/
  lunar_pnt_mission_study.m     the study script (numbered sections, PASS/FAIL gates)
  setup_lunar_pnt_paths.m       path setup: proj7/pipelines/gdop (addpath, never genpath),
                                the deliverable-9 kit folder, costate_common, pumpkynPie startup
  lib/
    schedule_rendezvous.m       the enumeration (section 4)
    polish_rendezvous.m         catalog cell -> exact-slot accepted transfer (section 5)
    assemble_histories.m        DRO coast + transfer + tulip slot on one time grid (section 6)
    mission_gdop.m              histories -> proj7 GDOP grid + R1 to R3 over a window (section 6)
  tests/
    run_tests.m
    test_schedule_rendezvous.m  test_polish_rendezvous.m
    test_mission_gdop.m         test_lunar_pnt_mission_study.m
  results/                      git-ignored; study .mat with run_meta, PNGs, optional movie
  README.md
```

House rules as everywhere in `orbit_transfer`: pumpkyn-style headers with sizes, `if nargin == 0` self-demo on library functions, never `i`/`j` as indices, never `norm`, no `%#ok`, ASCII README without a signature line, results never committed. Shadowing (checked 2026-11-10): the kit's `get_family_orbit` is byte-identical to `costate_common`'s; `transfers_at` and `catalog_transfer` exist identically in the kit and in `DRO_tulip/indirect`; `costate_catalog_pick` has a DIFFERENT copy in `GTO_DRO/catalog` but the study does not use it. The path setup adds the kit folder and the study asserts with `which` that the catalog functions it calls resolve to the kit and the instruments to `costate_common`.

## 4. `lib/schedule_rendezvous.m`

Clock: epoch 0 is when spacecraft k sits at DRO phase (k-1)/6. All times ND unless a name says days.

A candidate for spacecraft k and tulip slot j is a primary catalog cell (iD, iA) (decision 5) and a wait of n whole DRO revolutions:

- Departure time: `tauD*mod(sD(iD) - (k-1)/6, 1) + n*tauD`, DRO phase 0 being `get_family_orbit`'s first state (the six slots fall on the 24-point departure grid, since 1/6 = 4/24).
- Arrival time: departure time + the cell's t_f.
- Slot j's tulip phase at the arrival time, from the `tulipConstellation` convention.
- Residual: the wrapped gap between the cell's arrival phase sA(iA) and slot j's phase at arrival, as a fraction of the tulip period. The candidate is admissible when |residual| <= tol (default 0.005).

For each (k, j) keep the admissible candidate with the earliest arrival (and, as tie-break, the smaller |residual|). The objective is then a bottleneck assignment: over the 720 permutations (or an equivalent bottleneck-assignment solve), minimise the maximum arrival; ties broken by total propellant.

Inputs: the catalog struct, the DRO and tulip periods, the six DRO phases, the six slot offsets, `tol`, `nMaxRevs` (default 10, about 44 d of waiting; a knob, not a derived bound: the period ratio 5.236 is irrational, so there is no relative-phase cycle to cover, and each extra revolution shifts the slot's phase by 0.191 of a tulip period. My call; the sensitivity row reports 5 and 20 beside it). Outputs: the assignment, per-spacecraft cell, wait, departure and arrival times, residual, t_f, propellant (kg, from the regime constants), the last-arrival time, and the full 6 x 6 table of best arrivals so the choice can be audited. If any (k, j) has no admissible candidate within `nMaxRevs`, the function says so in the table and the assignment is over what exists; if no complete assignment exists it refuses with a named error.

Refused here: optimising over non-integer waits or over cells outside the tolerance. The tolerance exists because the polish (section 5) closes a small phase gap; the study reports how the last-arrival time moves with the tolerance (0.0025, 0.005, 0.01) as a sensitivity row, not as a tuning knob.

## 5. `lib/polish_rendezvous.m`

Turns one scheduled candidate into a transfer that arrives ON slot j's state at a pinned arrival epoch.

- The arrival epoch t_arr is pinned to the candidate's arrival time; the target state is slot j's exact state at t_arr.
- Unknowns: z8 = [lambda_r; lambda_v; lambda_m; t_f]. The departure time is t_arr - t_f, so the departure state is spacecraft k's DRO state at that time, which moves as t_f changes.
- Iteration: this is a scalar fixed point, t_f = T(x_DRO(t_arr - t_f), slot_j(t_arr)), where T is the min-time solve. Plain iteration contracts only where |dT/dt_dep| < 1; on the catalog's t_f grid the median slope along departure phase is 0.40 but 6 percent of edges exceed 1 (family hand-overs), so the polish is a bracketed scalar root-find (secant with a bracketing fallback) on g(t_f) = t_f - T(...), each evaluation one `tfMin` solve warm-started from the previous z8. Stop at |g| below a stated tolerance (default 1e-9 ND) or an evaluation cap; the count is reported. (My call, replacing the plain iteration of the presented design.)
- `tfMin` returns z with no exit flag, so convergence of each inner solve is judged by the instruments below (flown miss), not by the solver.
- The polished departure time t_arr - t_f must be >= 0 (the spacecraft cannot leave before the clock starts); a candidate whose polish needs an earlier departure is refused and the scheduler's next-best candidate for that pair is polished instead. A polished departure that differs from the scheduled one by the phase-gap amount is expected and reported.
- If the root lies on a different transfer family than the cell's (t_f jumps by more than a stated fraction, default 2 percent, from the cell's), the transfer still counts if it passes every acceptance gate, and is flagged `family changed` in the table.
- Acceptance on the polished z8, each a PASS/FAIL the study prints: `fly_transfer` flown arrival miss < 100 km (and the velocity miss it reports); `pmp_pointwise_checks`; `mintime_hypothesis_gates`; the conjugate test (`ms_tfmin` with `conjTest`, reproduced-solution rule `converged && max|dz| < 1e-6` as in `gto_dro_entry`). The model for how these four are called and thresholded on a 70 mN, 16 to 23 d transfer is `DRO_tulip/indirect/transfer_study.m` (the standard study script for this library); `gto_dro_entry` is the model only for the conjugate re-solve.
- A spacecraft whose polish does not converge or fails a gate is reported as such with its numbers; the study does not substitute the unpolished cell. The schedule's verdict is then FAIL for that spacecraft, and the README records it.

Builder's calls: the inner solver wiring (direct `pumpkyn.cr3bp.tfMin` call vs the kit's `catalog_transfer` options) and whether `ms_tfmin` is needed as a fallback when single shooting loses the basin. The pinned epoch and the departure state moving with t_f are fixed by this design.

## 6. `lib/assemble_histories.m` and `lib/mission_gdop.m`

**Histories.** One time grid, 300 s step, from epoch 0 to the last arrival + 30 d. For each spacecraft: its DRO coast until departure, the flown polished transfer, then its tulip slot from arrival onward. Output `rTgt [N x 3 x 6]` in the dimensionless barycentric rotating CR3BP frame (the same frame and units the transfers use; `gdop_grid_from_states` takes exactly that, checked 2026-11-10), a per-spacecraft phase label per sample (DRO / transfer / on-station), and the event times. Continuity is a gate: the position jump at each hand-over (DRO -> transfer, transfer -> slot) is below a stated tolerance.

**Geometry.** `gdop_grid_from_states(tau, rTgt, mu, lStar, tStar, params)` from `proj7/pipelines/gdop` on the two grids of decision 6: the paper's Table 1 grid for every requirement number, the coarse grid for the time-series figures and the assembly window. Elevation mask 0 degrees on both. The requirement scoring samples at the paper's cadence (about 348 s), not the 300 s history grid, so block lengths near the 4 h and 9 h thresholds quantise as the paper's do; the paper grid over 30 d is about 7,450 epochs x 10,770 points and needs the latitude chunking `run_requirements_table` already uses. Only on-station satellites count as navigation sources (a spacecraft in transfer is not broadcasting). `gdop_grid_from_states` has no satellite mask, so `mission_gdop` calls it once per interval over which the on-station set is constant, passing only that subset (`rTgt [N x 3 x M]`, M = satellites on station), and stitches the intervals; fewer than four on station is reported as no fix. NaN positions are forbidden (the function's spline resampling would poison the column). The function resamples on a half-open grid, so the stitch is tested for a duplicated or dropped epoch at each seam. (My call.)

**Two windows.**

- Assembly: from the first arrival to the last. Reported as the pole GDOP series and the number of satellites on station versus time; no requirement verdict is claimed while fewer than six are up, but the first time GDOP at the pole drops below 6 is reported.
- Service: the 30 days after the last arrival, scored for R1, R2 and R3 over a finite window with no periodic wrap. The proj7 requirement files (`verify_req1_link`, `verify_req2_gdop`, `verify_req3_eva`, `run_requirements_table`) hardcode their constellations and take no state histories, and proj7 is not modified; so `mission_gdop` applies `run_requirements_table`'s recipe on the proj7 utilities underneath it (`max_min_access_bounds`, `max_min_access_block_bounds`): R1 = visibility count >= 2 at 90 percent; R2 = GDOP < 6, best day-start, worst-day fraction; R3 = blocks of 4 h or more summing to 9 h, scored on GDOP < 6 exactly as `run_requirements_table` does for the paper's Table 1 (it uses the same GDOP < 6 coverage as R2; only `verify_req3_eva` uses GDOP < 7), with the GDOP < 7 version reported as an extra row. A ported recipe needs its own equivalence test: see the validation gate. (My call on how; that the mission is scored against the paper's three requirements is the design.)

**Validation gate.** With the six satellites placed directly on their slots (no transfers), `mission_gdop`'s ported recipe, on the paper grid over one tulip period, reproduces the six-satellite tulip column of `run_requirements_table` (the paper's Table 1: R1, R2, R3 numbers) to a stated tolerance; the reference is the saved result `proj7/figures/requirements_table_results.mat` (`S2(1).r`, current phasing, proj7 HEAD 3a7f898): R1 floor 1.000, R2 floor 0.754032, R3 floor 17.987 h. The script itself is not re-run (it computes all five columns and overwrites that file). `analyze_tulip_gdop_complete` is NOT the reference (it returns no R1 to R3 numbers, uses a closed grid and a 10 degree default mask); it may be used only as a second look at the pole GDOP series with its mask set to 0. The 30-day service numbers are reported beside the one-period ones (a longer window changes worst-day statistics; that difference is a result, not an error). If this gate fails, the study stops before claiming any mission result.

## 7. `lunar_pnt_mission_study.m`

Numbered sections in the study-script style (each condition computed in the script, PASS/FAIL printed, results saved with `run_meta`: date, git HEAD, MATLAB version, knobs):

0. Knobs: regime constants read from the catalog, tolerance, nMaxRevs, grid step, window lengths, figure/movie switches.
1. Inputs verified: catalog loads and its regime matches decision 1; `tulip_phasing()` returns the six offsets; DRO and tulip periods from `get_family_orbit` match the catalog's; the tulip-phase-times-tau0 identity holds.
2. Schedule: the 6 x 6 arrival table, the chosen assignment, the tolerance sensitivity row.
3. Polish: six transfers, iteration counts, the acceptance gates per spacecraft.
4. Histories: continuity gates.
5. GDOP validation gate against the paper's engine.
6. Assembly window and service window results; R1 to R3 table.
7. Figures: schedule chart (per-spacecraft bars: wait, transfer, on station), 3-D assembly view in the rotating frame, pole GDOP series with arrivals marked, a service-volume GDOP snapshot, propellant table. Slide-quality fonts; legends on every plot.
8. Movie: off by default; when on, in the style of `pumpkynPie/Demos/lowThrust_constellation_animation.m` (which has the launch-delay logic for one hardcoded transfer and no GDOP), frames sized to a multiple of 16.
9. VERDICT line and the saved `.mat`.

## 8. Tests

- `test_schedule_rendezvous`: a hand-built t_f table and periods small enough to enumerate by hand; the function returns the known bottleneck assignment, respects the tolerance, and refuses when no complete assignment exists.
- `test_polish_rendezvous`: one spacecraft from the real catalog; the polished transfer arrives on the slot state (flown miss gate), t_f converged, gates present.
- `test_mission_gdop`: the validation seam (on-slot satellites versus the recorded `run_requirements_table` six-satellite column, on a reduced grid if the full one is too slow for a test, with the full-grid comparison in the study), and the seam stitch (no duplicated or dropped epoch when the on-station set changes).
- `test_lunar_pnt_mission_study`: headless full run through `study_opts`; every gate PASS or reported N.A. with the reason; the VERDICT line, the PNGs and the `.mat` exist.
- `run_tests`: one runner, one logical verdict, throws on failure.

Every gate gets one deliberate-failure experiment recorded in the test or the README (a gate that cannot fail is not a gate).

## 9. Refusals, recorded in the README

CR3BP only (no ephemeris); no station-keeping after arrival (the tulip slot is the ideal periodic orbit); no navigation filter or signal model (geometry-only GDOP); the transfers use the 70 mN / 900 s / 150 kg regime, which is NOT the GTO -> DRO catalog's regime, so this is not a continuous GTO -> DRO -> tulip mission; all six spacecraft share one DRO and one clock; collision and proximity between spacecraft are not checked.

## 10. Risks carried into the plan

- **Shadowed functions** between the kit folder and `costate_common` (and the Aerospace Toolbox shadowing mlabTools coordinate functions on this machine): resolved and asserted in the path setup, never by bare-name luck.
- **Frame and unit seam:** checked at spec review (both sides are dimensionless, barycentric, rotating); the validation gate remains the proof.
- **Compute for the paper grid:** about 7,450 epochs x 10,770 points x 6 satellites for the service window; measured on a one-day slice before the full run, chunked by latitude.
- **The polish may leave the certified basin** for a spacecraft whose cell sits near a family hand-over (17 of 96 arrival edges were hand-overs in the sheet96 measurement): reported per spacecraft; the fallback order is the builder's call within section 5.
- **No admissible candidate** for some (spacecraft, slot) pair at the default tolerance within the revolution bound: the table shows it; the sensitivity row shows whether a looser tolerance changes the assignment.
- **Machine load:** the GTO -> DRO fleet runs in a background `matlab -batch` driver at the same time; this study runs in the shared session. The study's heavy steps (GDOP grid over ~50 d at 300 s) must be measured once and sized before a full run.

## 11. Acceptance

The study script runs from a fresh session; the validation gate passes; six polished transfers each pass the four acceptance instruments (or the failures are reported as the result); the R1 to R3 table for the 30-day service window is produced by proj7's scorers; the figures exist; `run_tests` passes; the README states the measured schedule (assignment, last-arrival day, propellant), the requirement verdicts and the refusals.
