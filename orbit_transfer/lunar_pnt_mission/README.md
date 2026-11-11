# lunar_pnt_mission -- six spacecraft DRO -> tulip slots, then the south-pole GDOP time series

Six spacecraft parked in one 4.43 d distant retrograde orbit each fly a
certified 70 mN minimum-time low-thrust transfer onto one of the six slots of
the 7-petal tulip constellation (the paper's phasing), timed so that each
lands on its exact slot state at a pinned arrival epoch. The product is ONE
time series: GDOP at the lunar south pole (latitude -90, 0 deg mask) from
clock 0 -- nothing on station yet -- through the assembly and for
`serviceDays` after the last arrival, computed with proj7's own GDOP engine
after a validation gate checks the pole series against the paper's own
constellation construction. The LunaNet R1-R3 requirement scoring
(validated against the paper's Table 1 floors) is built and tested but sits
behind `scoreRequirements = true` and is not run by default. Design:
`docs/superpowers/specs/2026-11-10-lunar-pnt-mission-design.md`.

## Run

```matlab
cd orbit_transfer/lunar_pnt_mission
lunar_pnt_mission_study    % edit the knobs block, press Run (CLEARS the workspace it runs in)
cd tests; run_tests        % every test, one verdict, throws on failure (~10 min)
```

The script clears its workspace (`clearvars -except study_opts`). Run it
from a function to keep the base workspace (the test does this).
Paths come from `setup_lunar_pnt_paths` (this dir, lib/, costate_common, the
deliverable-9 kit, proj7/pipelines/gdop, oclib, casadi, pumpkyn); it asserts
with `which` that every function resolves in the folder it is meant to.
Outputs go to `results/` (git-ignored): four PNGs (five with
`scoreRequirements`) and `lunar_pnt_mission_study.mat`:
`pole` (tDays, gdop, nOn, nVis, arrDays, stats), `S S0 P Pattempted
polishCache H Vpole V Gservice gates rounds I2 run_meta`; positions are
stored as single.

## Knobs (section 0 of the study script)

| knob | default | meaning |
|---|---|---|
| `tol` | 0.005 | arrival phase tolerance [fraction of a tulip period] |
| `nMaxRevs` | 10 | longest DRO wait [whole revolutions] |
| `dtSec` | 300 | history grid step [s] |
| `serviceDays` | 30 | service window after the last arrival [d] |
| `validate` | true | run the validation gate(s); false gives STOPPED AT VALIDATION (G0 NOT RUN) |
| `maxRounds` | 12 | schedule -> polish -> exclude rounds before giving up |
| `polishGates` | true | full certification on every polish (false: fly gate only) |
| `i2Km` | 1 | I2 limit: orbit handle vs direct propagation [km] |
| `i2Phases` | 2400 | I2 uniform phases (plus 200 around the tulip perilune) |
| `scoreRequirements` | false | also score R1-R3 on the paper grid (G0R validation vs Table 1, R gates, service_map.png) |
| `scoreGrid` | 'paper' | grid R1-R3 are scored on (scoreRequirements only) |
| `chunkSel` | [] | score-grid chunks ([] = all five; scoreRequirements only) |
| `reuse` | '' | an earlier `lunar_pnt_mission_study.mat` whose polishes seed the cache, guarded (see "Reusing polishes") |
| `g0SlotOverride` | [] | TEST SEAM: slot offsets for G0's own series (the deliberate-failure experiment); [] = none |
| `makeFigures` | true | section 7 |
| `makeMovie` | false | section 8: after the save, `mission_movie` renders `results/lunar_pnt_mission.mp4` from the saved results |
| `headless` | false | figures invisible and closed after export |
| `outDir` | '' | results folder ('' = results/) |

`study_opts` in the workspace overrides any knob; unknown names are refused.
A run with any result-changing knob off its default says so in its VERDICT
line, e.g. `VERDICT: PASS (reduced: serviceDays 3)`.

## How the schedule is found

`schedule_rendezvous` assigns spacecraft to slots so the last arrival is as
early as possible, each (spacecraft, slot) pair represented by its earliest
catalog cell whose arrival phase lies within `tol` of the slot at the arrival
epoch (whole-revolution DRO waits allowed). `polish_rendezvous` then moves
each transfer onto the slot's exact state at that pinned epoch, sliding the
departure along the DRO; it can refuse (a root that needs a departure before
clock 0, a root not found, a failed or erroring certification gate). A refusal
changes the bottleneck, so the study re-optimises the whole assignment rather
than walking one spacecraft's list: every refused candidate `[k j iD iA nRev]`
joins the scheduler's `opts.exclude` and stays excluded in every later round
(the list is cumulative), the schedule is recomputed, and only spacecraft
whose chosen candidate changed are polished again (polishes are cached by
candidate key). The loop ends when all six chosen candidates are accepted,
after `maxRounds`, or when no complete assignment remains (C1 FAIL). Every
round is printed and saved in `rounds`; on a stop inside the loop the last
attempted round's polishes are saved in `Pattempted` and its table is
labelled as such.

## Gates

| id | condition |
|---|---|
| I1 | every `mission_inputs` check (catalog, periods, phasing, regime, two flown catalog cells land on the tulip) |
| I2 | the periodic phase_state handles (tulip and DRO) agree with a direct ode113 propagation (tol 1e-13) to < 1 km at >= 2,400 phases incl. the perilune neighbourhood |
| C1 | a complete assignment exists (in every round) |
| P1-P6 | spacecraft k's final candidate accepted by the polish |
| PF | flown miss < 100 km and < 10 m/s on all six |
| PR | ms_tfmin reproduces the polished root, residual < 1e-8 |
| PP | pmp_pointwise_checks (H, lambda_m, adjoint, minimum principle, throttle, Q_mt) on all six |
| PH | hypothesis gates: min abs(lambda_v) > 0, min Q_mt > 0, dim S = 1 |
| PC | conjugate test passed and covered the arc (PASS / FAIL / UNRESOLVED) |
| PA | pumpkyn tfMin reproduces each polished z8 to 1e-6 |
| H1 | DRO -> transfer position jump < 1 km on all six |
| H2 | transfer -> slot position jump < 100 km on all six |
| G0 | `validate_pole_series`: on-slots pole series over one tulip period vs a direct engine call on `tulipConstellation` states; relative GDOP difference < 1e-3 (GDOP < 20, same visibility count), visibility and finite/NaN mismatches only next to a visibility edge |
| G0R | scoreRequirements only: R1-R3 floors within 0.005 / 0.005 / 0.10 h of the paper's Table 1 column on the same points |
| R1-R3 | scoreRequirements only: the LunaNet floors over the service window on the paper grid |

Anything but PASS from G0 (or G0R), including NOT RUN, stops the study
before the GDOP section. `VERDICT: PASS` only when every gate is PASS;
otherwise the line names the failing gates, and a stop (I2, C1, polish
loop, validation) is named as such. The pole series itself carries
statistics, not verdicts: first finite and first GDOP < 6 day, and over the
service window the fraction of time GDOP < 6, the fraction finite, median
and 95th percentile GDOP, the longest outage (GDOP >= 6 or NaN) and the
count of samples with 6 <= GDOP < 7.

That last count is zero in every run so far, and the R3 rows on GDOP < 6
and GDOP < 7 come out identical. This is MEASURED for this constellation,
not a general property: on the pole series and on one paper-grid latitude
chunk (chunk 5, -77.5 to -75 deg, five days on slots), GDOP was at most
5.81 whenever 5 or 6 satellites were in view and at least 16.7 whenever
exactly 4 were, so no sample fell in [6, 7). It has not been checked on the
full grid.

## Reusing polishes

`reuse` points at an earlier results file; its polishes seed the cache so
the GDOP part can be re-run without re-solving. It is guarded:

- The whole reuse is refused (`lunar_pnt_mission_study:reuse`) unless the
  file's `fingerprint` equals this run's: catalog name and created date,
  MD5 of `sh.z8` and `sh.tf_nd`, thrust, Isp, m0, muStar, lStar, tStar,
  `dtau6`, the polish thresholds and `polishGates`. A file without a
  fingerprint (made before this guard) is refused.
- Each cached entry, when taken for the current candidate, must have the
  candidate's arrival epoch (1e-12 ND), `tDep + tf = tArr` (1e-12), the
  current slot state as its target (1e-12) and the current thresholds; an
  entry that fails is ignored, printed as IGNORED, and polished afresh.
- Reused polishes are NOT re-certified in the run that reuses them:
  `run_meta.polishSource` records fresh / reused per spacecraft with the
  certifying file, date and git HEAD, the P1..P6 / PF / PR / PP / PH / PC /
  PA rows are tagged REUSED, and the VERDICT reads `PASS (REUSED polishes
  ... certified <date> at <HEAD>; not re-certified in this run)`.

## Movie

```matlab
info = mission_movie('');                              % results/lunar_pnt_mission.mp4 (~2 min)
info = mission_movie('', struct('preview', true));     % 14 story frames -> results/movie_preview/
info = mission_movie('', struct('dayList', [0 10 30])); % PNGs at chosen mission days
```

`lib/mission_movie.m` reads `results/lunar_pnt_mission_study.mat` (H, P, pole)
and re-solves nothing. 1280 x 720, 30 fps, 70 s, dark background. Top: the
Earth-Moon rotating frame (orthographic camera, no axes; a km scale bar and an
Earth-direction pointer, since Earth is 389,703 km away), the textured Moon,
the DRO and tulip outlines with labels (one CR3BP period each, closure
asserted; drawn in every frame at constant brightness), six spacecraft with
fixed colours and 2-day fading trails, marked parked (hollow), thrusting
(filled, with plume) or on station (filled with a ring). Bottom: the
south-pole GDOP series drawn up to the current day on a fixed 2-8 axis
(values above 8, and infinite GDOP, are clipped and marked; undefined GDOP is
a gap), the GDOP = 6 line, coloured arrival ticks and the number of
satellites on station on a right axis. Counters (day, on station, current
GDOP) and the story caption (hard cuts; the "parked" caption ends at the first
departure) sit at fixed positions.

Camera: the view direction (azimuth 12 deg, elevation 18 deg) and the target
never change; only the zoom does. The first and last frames are the end
framing, where the DRO and the tulip fill the scene. The zoom comes from the
data: for every frame, the view size that keeps all six spacecraft (12 %
margin while transferring) and their trails inside the scene area, never
smaller than the end framing, made smooth and single-peaked (running max,
unimodal hull, raised-cosine smoothing over several seconds, re-inflated
where smoothing cut below the requirement). Measured on the saved run: one
zoom-out and one zoom-in, widest 1.75 x the end framing at day 12.5, fastest
change 6.6 % of the view size per second of screen time, within 2 % of the
end framing from day 19.7 and exactly on it from day 25.5. The renderer
asserts that every spacecraft is at least 4 px inside the scene area in every
frame (smallest margin 16.5 px). The scale bar shows the largest round length
(1,000 to 100,000 km) whose bar fits in 160 px and is redrawn every frame.
Markers, trails and plumes are still clipped to the scene area as a safety
net; no trail reaches its edge in the rendered movie.

The Moon follows `pumpkyn.util.moon3D` (used by `showMoon` in the other orbit
movies): an ellipsoid mesh with x, y, z negated and pumpkyn's `Moon_2k.jpg`
as a texture map, loaded by path from `~/Desktop/proj7/external/pumpkyn/src/
+pumpkyn/+util/` (`moonTexture` option); north is +z and the prime meridian
faces the Earth (-x), with no rotation in this frame (tidal lock). If the
file is missing it falls back to a shaded grey sphere and prints one line.
Its lighting is a fixed light near the camera, chosen so the texture reads;
it is not the physical Sun direction. The Moon is drawn to scale: about 33 px
across at the end framing and 19 px at the widest zoom, so the maria are
visible but fine detail is not.

The time map is non-uniform and smooth: a 2 s hold at day 0, about 5.5 s for
the first DRO revolution, 29 s for the transfers to the last arrival, 30.6 s
for the 30-day service window, a 3 s hold; the day counter and a days-per-
second readout show the changing speed.

The thrust plume points along the exhaust, `+lam_v/|lam_v|` (the flown thrust
is `u = -lam_v/|lam_v|`; checked against the velocity change not explained by
CR3BP gravity), and switches off at arrival. Its length is a fixed number of
pixels (`plumePx`) at every zoom: it is illustrative, not to scale.

The MP4 is written to `<name>.partial.mp4` and moved onto `outFile` only after
VideoReader confirms its frame count and size, so a failed render never
truncates the movie. A `maxFrames` render needs an explicit `outFile` that is
not the default movie (`mission_movie:partialOutFile`).

Options (`opts`): `outFile`, `fps` [30], `durationSec` [70], `frameSize`
[1280 720], `visible` [false], `maxFrames` [inf], `frameList` / `dayList`
(PNG previews; days must lie in the mission), `preview` [false],
`previewDir`, `tailDays` [2], `plumePx` [62], `moonTexture`, `probeKm` (test
markers), `quality` [95]. `info` returns the frame days and drawn states,
plume counts and directions (3-D and on screen), the camera and zoom per
frame with the containment margins, outline visibility, the Moon texture
state and the time map.

## Tests

`tests/run_tests.m` runs all eight tests. `test_lunar_pnt_mission_study`
runs the script in throwaway function workspaces: the default pole-series
path (3 service days, the polish loop for real), `validate = false`,
`maxRounds = 1`, `tol = 1e-7`, a reduced `scoreRequirements = true` run
(chunk 5, 3 days, checking the REUSED disclosure), reuse refused for a
different fingerprint and for a file without one, a tampered cached entry
ignored and re-polished, and the G0-FAIL stop through `g0SlotOverride`;
the later runs reuse the first run's polishes.
`test_mission_movie` renders 12 frames to a temporary MP4 and three preview
PNGs from the saved results (it skips, passing, when no results file exists)
and checks the size, the frame count, plume on mid-transfer and off after the
last arrival, and the drawn plume direction against `+lam_v/|lam_v|`.

## Refusals

- CR3BP only, no ephemeris model.
- No station-keeping after arrival: the tulip slot is the ideal periodic orbit.
- No navigation filter or signal model: GDOP is geometry only.
- By default no requirement verdicts (R1-R3): the product is the south-pole
  GDOP series. The scoring exists behind `scoreRequirements`.
- The transfers use the 70 mN / 900 s / 150 kg regime, which is NOT the
  GTO -> DRO catalog's regime, so this is not a continuous GTO -> DRO -> tulip
  mission.
- All six spacecraft share one DRO and one clock.
- Collision and proximity between spacecraft are not checked.
- No propellant-optimal schedule: propellant is reported, not optimised.

## Measured results
