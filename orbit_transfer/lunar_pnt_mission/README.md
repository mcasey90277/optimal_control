# lunar_pnt_mission -- six spacecraft DRO -> tulip slots, then the south-pole GDOP time series

Six spacecraft parked in one 4.43 d distant retrograde orbit each fly a
certified 70 mN minimum-time low-thrust transfer (Measured results says what
"certified" covers) onto one of the six slots of
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

The study has two legs. Leg 1 (section "Leg 1: GTO to DRO", `lunar_pnt_leg1_study`)
flies six 5 N transfers from one inertial GTO onto the six DRO slots and hands
over at epoch T0; leg 2 (everything else below, `lunar_pnt_mission_study`)
starts there. Leg 1's regime is a catalog regime, not a vehicle design.

## Run

```matlab
cd orbit_transfer/lunar_pnt_mission
lunar_pnt_mission_study    % edit the knobs block, press Run (CLEARS the workspace it runs in)
cd tests; run_tests        % every test, one verdict, throws on failure (~15 min)
```

The script clears its workspace (`clearvars -except study_opts`). Run it
from a function to keep the base workspace (the test does this).
Paths come from `setup_lunar_pnt_paths` (this dir, lib/, costate_common, the
deliverable-9 kit, proj7/pipelines/gdop, oclib, casadi, pumpkyn); it asserts
with `which` that every function resolves in the folder it is meant to.
Outputs go to `results/` (git-ignored): four PNGs (five with
`scoreRequirements`) and `lunar_pnt_mission_study.mat`:
`pole` (tDays, gdop, nOn, nVis, arrDays, stats), `S S0 P Pattempted
polishCache H Vpole V Gservice gates rounds I2 fingerprint run_meta`,
and (with `exportTransfers`) the sidecar `mission_transfers_70mN.mat` (see
"Captured transfers (sidecar)").
`H.r` is saved in double precision (18.0 MB file, about 1 MB more than
single), so the saved history passes `mission_gdop`'s 1e-9 ND slot check
and can be fed back into it. `run_meta` records the date, MATLAB version,
knobs, the git HEAD of this repository (found from the script's folder) and
of proj7, each with a dirty flag (`optimal_controlDirty`, `proj7Dirty`:
tracked files modified), and `verdictSuffix`, the VERDICT's disclosure
text. The run header prints the same git state, e.g. `optimal_control
1d43572 (DIRTY: 105 tracked path(s) modified)`.

## Knobs (section 0 of the study script)

| knob | default | meaning |
|---|---|---|
| `tol` | 0.005 | arrival phase tolerance [fraction of a tulip period] |
| `nMaxRevs` | 10 | longest DRO wait [whole revolutions] |
| `dtSec` | 300 | history grid step and south-pole series cadence [s] (G0's own one-period comparison keeps its 300 s grid) |
| `serviceDays` | 30 | service window after the last arrival [d] |
| `validate` | true | run the validation gate(s); false gives STOPPED AT VALIDATION (G0 NOT RUN) |
| `maxRounds` | 12 | schedule -> polish -> exclude rounds before giving up |
| `polishGates` | true | full certification on every polish (false: fly gate only) |
| `i2Km` | 1 | I2 limit: orbit handle vs direct propagation, position [km] |
| `i2Mps` | 0.1 | I2 limit: the same, velocity [m/s] (measured max 0.0016 m/s) |
| `i2Phases` | 2400 | I2 uniform phases (plus 200 around the tulip perilune) |
| `scoreRequirements` | false | also score R1-R3 on the paper grid (G0R validation vs Table 1, R gates, service_map.png) |
| `scoreGrid` | 'paper' | grid R1-R3 are scored on (scoreRequirements only) |
| `chunkSel` | [] | score-grid chunks ([] = all five; scoreRequirements only) |
| `reuse` | '' | an earlier `lunar_pnt_mission_study.mat` whose polishes seed the cache, guarded (see "Reusing polishes") |
| `g0SlotOverride` | [] | TEST SEAM: slot offsets for G0's own series (the deliberate-failure experiment); [] = none |
| `makeFigures` | true | section 7 |
| `makeMovie` | false | section 8: after the save, `mission_movie` renders `results/lunar_pnt_mission.mp4` from the saved results |
| `exportTransfers` | true | section 9: after the save, `export_mission_transfers` writes the six transfers as costate-library entries to `mission_transfers_70mN.mat` beside the results; NOT written (printed line) when the run stopped, i.e. before six accepted transfers or without a PASSED validation; an export error is printed as `EXPORT FAILED: <id>: <message>`, recorded in `run_meta.export` (appended to the results file) and the run goes on (movie included) |
| `exportFullStatus` | true | the library's certifier on each exported transfer (about 4 min for six); false: the fallback status (status 1, seconds) |
| `headless` | false | figures invisible and closed after export |
| `outDir` | '' | results folder ('' = results/) |

`study_opts` in the workspace overrides any knob; unknown names are refused.
A run with any result-changing knob off its default says so in its VERDICT
line, e.g. `VERDICT: PASS (reduced: serviceDays 3)`.

## How the schedule is found

`schedule_rendezvous` assigns spacecraft to slots so as to minimise the last
arrival OVER ITS CANDIDATE SET: each (spacecraft, slot) pair is represented by
its earliest primary catalog cell whose arrival phase lies within `tol` of the
slot at the arrival epoch (whole-revolution DRO waits allowed), and that
cell's arrival epoch is then PINNED. `polish_rendezvous` moves
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

Every gate has a deliberate-failure experiment (spec section 8): a test that
makes it fail on purpose, so the quantity behind each gate is known to move. For I1, H1 and H2 the
experiment shows the measured quantity crossing the threshold (the oracle or
the jump is sensitive); it does not run the study's own gate row to FAIL, so
those three tests would not notice the gate row being deleted. For PC the
experiment exercises the reproduction prerequisite, not a conjugate FAIL.

| id | condition | shown to fail by |
|---|---|---|
| I1 | every `mission_inputs` check (catalog, periods, phasing, regime, two flown catalog cells land on the tulip) | `test_mission_inputs` (quantity only): the same cell flown with 1.01 x Tnd misses by 1.37e4 km |
| I2 | the periodic phase_state handles (tulip and DRO) agree with a direct ode113 propagation (tol 1e-13) to < 1 km and < 0.1 m/s at >= 2,400 phases incl. the perilune neighbourhood | `test_lunar_pnt_mission_study` (J): `i2Km = 1e-6`, then `i2Mps = 1e-6` -> I2 FAIL, STOPPED AT INPUTS |
| C1 | a complete assignment exists (in every round) | `test_lunar_pnt_mission_study` (D): `tol = 1e-7`; `test_schedule_rendezvous` (4, 6): `noAssignment` |
| P1-P6 | spacecraft k's final candidate accepted by the polish | `test_polish_rendezvous` (3a, 3b, 4-8): refusals, failed and erroring gates give ok = false |
| PF | flown miss < 100 km and < 10 m/s on all six | `test_polish_rendezvous` (4): `flyKmMax = 1e-9` |
| PR | ms_tfmin reproduces the polished root, residual < 1e-8 | `test_polish_rendezvous` (7): `residMax = 0` |
| PP | pmp_pointwise_checks (H, lambda_m, adjoint, minimum principle, throttle, Q_mt) on all six | `test_polish_rendezvous` (7): `HmaxMax = 0` |
| PH | hypothesis gates: min abs(lambda_v) > 0, min Q_mt > 0, dim S = 1 | `test_polish_rendezvous` (6): `dimS = 2` |
| PC | conjugate test passed and covered the arc on a reproduced root (PASS / FAIL / UNRESOLVED) | `test_polish_rendezvous` (8): `reproTol = 0` (the verdict only counts for a reproduced root) |
| PA | pumpkyn tfMin reproduces each polished z8 to 1e-6 | `test_polish_rendezvous` (7): `accDzMax = 0` |
| H1 | DRO -> transfer position jump < 1 km on all six | `test_assemble_histories` (3) (quantity only): tDep shifted 0.05 ND -> 9,090 km |
| H2 | transfer -> slot position jump < 100 km on all six | `test_assemble_histories` (3) (quantity only): tArr shifted 0.05 ND -> 16,354 km |
| G0 | `validate_pole_series`: on-slots pole series over one tulip period vs a direct engine call on `tulipConstellation` states; relative GDOP difference < 1e-3 (GDOP < 20, same visibility count), visibility and finite/NaN mismatches only next to a visibility edge | `test_validate_pole_series` (3) and `test_lunar_pnt_mission_study` (I): slot 6 moved 0.3 ND -> FAIL, STOPPED AT VALIDATION |
| G0R | scoreRequirements only: R1-R3 floors within 0.005 / 0.005 / 0.10 h of the paper's Table 1 column on the same points | `test_mission_gdop` (7): satellite 6 moved 0.3 ND -> FAIL |
| R1-R3 | scoreRequirements only: the LunaNet floors over the service window on the paper grid | covered by the G0R experiment (the same scorers) |

Checks inside the library have their own experiments: `mission_gdop:slot`
(an on-station position moved 1e-6 ND, `test_mission_gdop`),
`mission_gdop:phase`, `assemble_histories:flightGrid` (NaN / Inf / reversed
flight times, `test_assemble_histories`), `mission_movie:containment` (zoom
envelope forced to half through the `zoomScale` test seam,
`test_mission_movie`) and the reuse guards (`test_lunar_pnt_mission_study`
F, G).

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
  `dtau6`, the spacecraft DRO phases `droPhase0`, MD5 of the DRO and tulip
  orbit handles sampled at 64 phases (`droMd5`, `tulipMd5`), the polish
  thresholds and `polishGates`. A file without a fingerprint, or whose
  fingerprint lacks one of these fields, is refused: every results file
  written before 2026-11-10 21:00 (no `droPhase0` / `droMd5` / `tulipMd5`)
  is therefore no longer reusable.
- Each cached entry, when taken for the current candidate, must have the
  candidate's arrival epoch (1e-12 ND), `tDep + tf = tArr` (1e-12), the
  current slot state as its target (1e-12), the spacecraft's current DRO
  state at its cached `tDep` as its departure state `rv0` (1e-12) and the
  current thresholds; an entry that fails is ignored, printed as IGNORED,
  and polished afresh.
- Reused polishes are NOT re-certified in the run that reuses them:
  `run_meta.polishSource` records fresh / reused per spacecraft with the
  certifying file, date and git HEAD (with `-dirty` appended when that run's
  working tree had modified tracked files), the P1..P6 / PF / PR / PP / PH /
  PC / PA rows are tagged REUSED, and the VERDICT reads `PASS (REUSED
  polishes ... certified <date> at <HEAD>; not re-certified in this run)`.
- Reused REFUSALS decide the exclusions, so they are disclosed too, even
  when all six final polishes are fresh: `(exclusions reused: 3 refused
  polish(es) from <file>, certified <date> at <HEAD>; not re-polished in
  this run)`.

## Movie

```matlab
info = mission_movie('');                              % results/lunar_pnt_mission.mp4 (~2 min)
info = mission_movie('', struct('preview', true));     % 14 story frames -> results/movie_preview/
info = mission_movie('', struct('dayList', [0 10 30])); % PNGs at chosen mission days
```

`lib/mission_movie.m` reads `results/lunar_pnt_mission_study.mat` (H, P, pole)
and re-solves nothing. It takes muStar, lStar and tStar from the file's
`fingerprint` (the catalog's values; `mission_inputs()` when a file has none)
and reports in `info.provenance` (never on screen) whether the run was a fresh
full-knob run or carried a REUSED / reduced note (`run_meta.verdictSuffix`).
1280 x 720, 30 fps, 70 s, dark background. Top: the
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
markers), `quality` [95], `zoomScale` [1] (test seam: scales the zoom
envelope before the containment check). `info` returns the frame days and
drawn states, plume counts and directions (3-D and on screen), the camera and
zoom per frame with the containment margins, outline visibility, the Moon
texture state, the time map and the provenance.

## Captured transfers (sidecar)

Every run that accepts all six transfers AND passes validation (a REUSED
run too) writes `mission_transfers_70mN.mat` beside its results: ONE
variable, `mission_transfers`, holding the six new transfers as
costate-library entries. A stopped or unvalidated run writes none. New
transfers are an opportunity to add to the library for free, but the shipped
library of record
(`DRO_tulip/deliverables/costate_library_70mN_user_deliverable9/costate_catalog_dro_tulip_70mN.mat`)
is NOT modified: `export_mission_transfers` refuses an output inside the kit
folder (the target is resolved to a canonical absolute path first, so a
relative path, a `..` or a not-yet-existing subfolder cannot slip through),
and its test checks the library file's MD5 before and after.

- `.entries` [1 x 6]: EXACTLY the fields of the library's
  `cat_.alternatives` (same names, order, class and shape), off grid
  (`iD = iA = 0`): `sD` = the DRO phase at departure, `sA` = the tulip phase
  at arrival, `source` = e.g. `lunar_pnt_mission 2026-11-10: sc 3 -> slot 6,
  polished from primary (16,31)`. A CERTIFIED entry stores what the library
  itself would store (`make_alternative`: `z8 = C.z`, `junctions = C.Y`,
  `flyKm = C.flyKm`): the certifier's root, its [14 x 24] multiple-shooting
  junctions (column 1 = `[rv0; 1; z8(1:7)]`) and its flown miss. z8 and
  junctions come from the same solve because the kit flies the junctions
  over `h = z8(8)/24` (`catalog_transfer.m` 171-184). The root the mission
  flew (the polish's) is `.mission.polishZ8`; the two differ by
  `.mission.certifyDz` (at most 4.8e-11 here; more than 1e-6 means the
  status is not attributed and the entry falls back). A FALLBACK entry
  stores the polish's z8, the K = 24 cut of its single flight as junctions
  (matching library primary (16,31)'s stored junctions in column 1 exactly
  and elsewhere to 3.3e-7 relative) and the polish's flown miss.
- `.mission`: spacecraft, slot, seed cell, wait, departure and arrival
  epochs (ND and days), propellant, dV, the polish's gates and thresholds,
  nEval, `reused` and the polish's origin, `statusRun` ('certify_root
  (full)' or 'fallback: <why>'), `polishZ8`, `certifyDz`, and the
  certifier's numbers.
- `.catalog`: the results file's catalog fingerprint (checked against the
  catalog on the path with the same recipe, `lib/catalog_identity.m`; a
  mismatch is refused), the regime, the orbit families. `.provenance`: the
  results' `run_meta`, the export date, `fenced`, and a disclosure note.
  `.usage`: how to fly one.

Fly an entry with the kit:

```matlab
L = load('costate_catalog_dro_tulip_70mN.mat');  cat_ = L.costate_catalog_dro_tulip_70mN;
load('results/mission_transfers_70mN.mat');      % mission_transfers
F = catalog_transfer(cat_, mission_transfers.entries(1));                         % junctions
F = catalog_transfer(cat_, mission_transfers.entries(1), struct('method', 'z8')); % one flight
```

**Status.** Each entry is run through the library's OWN certifier,
`certify_root` (`DRO_tulip/indirect`; appended to the path with `-end` only
when it does not already resolve, so nothing the campaign asserts can be
shadowed; the session path is restored exactly afterwards), and the
library's own packager `make_alternative` stamps the status with the
library's rule, `optimality_status`: 4 only when `C.ok` is true, i.e. every
certify_root stage passed (polish, flight, tight-flight pointwise checks,
tfMin witness and its flight, conjugate test, hypothesis gates with X2 and
lift_margin, H6, dense conjugate scan); 3 or 2 for a stack stopped at or
after the witness stage (2 only for a definite conjugate refutation with
gates and H6 held); 1 before it. The certifier re-polishes; its root must
reproduce the stored z8 to 1e-6 or the status is not attributed (tested: a
z8 costate moved by 1e-5 comes back as a fallback entry, status 1, "certify_root
converged to a root 1.0e-05 from the stored z8"). Measured cost 19-63 s per
entry (mean 40 s, about 4 min for six), so `fullStatus` defaults to true.
With `fullStatus = false` the polish's own evidence reaches certify_root's
stage 2 only (its pointwise checks run on the loose flight and it does not
fly the witness's solution), so by the same rule every entry is status 1 and
its `status_reason` lists the checks not run. All six 2026-11-10 transfers
are status 4, "full stack passed".

**Unfenced.** The certifier runs without a worker pool (`opts.pool = []`,
recorded as `provenance.fenced = false`), as the polish does: an unattended
run can hang inside it, the same exposure as the library's own unfenced
certification. Pass a connected `parallel.Pool` as `opts.pool` to fence it
(`fenced = true`; the workers need the campaign's paths).

## Leg 1: GTO to DRO

Leg 1 puts the six spacecraft into the DRO that leg 2 starts from. One launch
leaves six stacks on ONE geostationary transfer orbit (350 x 35,786 km) that is
fixed in inertial space (in the Earth-Moon rotating frame its orientation turns
at the frame rate, sgn = -1, measured). Each stack leaves the GTO at a catalog
departure phase and flies a 5 N per 150 kg, Isp 1710 s minimum-time transfer
with the engine on throughout, landing exactly on its DRO slot. Every departure
is polished: the GTO -> DRO catalog cell's costates are walked to the GTO's
true orientation at that departure, then the arrival is root-found onto the
exact slot. The six slots are 1/6 of a DRO period apart, and the hand-off epoch
T0 is the first instant after the last arrival at which the satellites sit at
the DRO phases (k-1)/6 leg 2 starts from; leg 2 then runs unchanged from T0.
Design: `docs/superpowers/specs/2026-11-11-lunar-pnt-leg1-gto-dro-design.md`.

**THE CAVEAT. These are minimum-time, engine-always-on transfers in the
catalog's numerical regime (5 N per 150 kg, Isp 1710 s), not a vehicle
design.** The numbers that make this plain: each transfer burns 41.8 to 45.0 %
of the stack's launch mass (mf/m0 = 0.550 to 0.582), the six stacks total
1,774 kg wet and 760 kg of stage propellant to deliver 900 kg, the delta-V is
8.64 to 10.02 km/s, and each stack needs 9.3 to 10.4 N and 78 to 87 kW of jet
power. A propellant-realistic leg 1 would be a minimum-fuel or much lower-thrust
transfer, and this catalog does not contain one. The stage dry fraction 0.15
(dry mass per kg of stage propellant) is an assumption, a knob (`dryFrac`), not
a design value. Stack mass, thrust and power are scaled from the catalog
trajectory by the closed form below and say only what the catalog regime would
cost.

### Run leg 1

```matlab
cd orbit_transfer/lunar_pnt_mission
lunar_pnt_leg1_study                        % CLEARS its workspace (clearvars -except leg1_opts); ~100 s
% optional, BEFORE the run: leg1_opts = struct('doSecond', false);  overrides any knob
leg1Mat = fullfile(pwd, 'results', 'lunar_pnt_leg1_study.mat');
leg2Mat = fullfile(pwd, 'results', 'lunar_pnt_mission_study.mat');
[T, info] = combined_timeline(leg1Mat, leg2Mat, struct('outPng', fullfile(pwd, 'results', 'mission_timeline.png')));
leg1_movie(leg1Mat, struct());              % results/lunar_pnt_leg1.mp4
full_mission_movie(leg1Mat, leg2Mat, struct());   % results/lunar_pnt_full_mission.mp4
```

| call | needs | writes (all under `results/`, git-ignored) |
|---|---|---|
| `lunar_pnt_leg1_study` (edit the knobs block, or set `leg1_opts`) | the GTO -> DRO package, the leg-2 results for the timeline | `lunar_pnt_leg1_study.mat` (`S S0 P1 T0 T0sched knobs run_meta failed gates exclude massTot slotOutcome`), `leg1_transfers_gto_dro_5N.mat` (sidecar), `leg1_schedule.png`, `leg1_3d.png`, `leg1_mass.png`, `mission_timeline.png` |
| `[T, info] = combined_timeline(leg1Mat, leg2Mat, opts)` | the two results files, a package file with the leg-1 run's SHA-256 (`opts.catalogMat`, else a snapshot; `leg1_package_for`) | the figure when `opts.outPng` is set; `T` is a six-row table of [start end] pairs in days on the leg-1 clock (`gtoWait leg1 droCoast droWait leg2 onStation`); `info.T0Days`, `info.endDays`, `info.packageMatch`, `info.packageFile`, `info.packageSource`, `info.checks` |
| `info = leg1_movie(leg1Mat, opts)` (`''` = default results file) | `lunar_pnt_leg1_study.mat` with verdict PASS, six accepted with `G.pass`, no test seam (else `leg1_movie:verdict`; `opts.allowUnvetted` stamps UNVETTED RESULTS) | `lunar_pnt_leg1.mp4` (1280 x 720, 30 fps, 30 s); `opts.preview = true` writes `lunar_pnt_leg1_preview.mp4`, `opts.dayList` / `opts.frameList` write PNG stills to `movie_preview/`; never `lunar_pnt_mission.mp4` (`leg1_movie:protected`) |
| `info = full_mission_movie(leg1Mat, leg2Mat, opts)` (`''` = default files) | both results files, compatible by `combined_timeline`'s checks (run first); with `render = false`, both MP4s newer than their results (else `full_mission_movie:stale`) | `lunar_pnt_mission_dayoffset.mp4` (leg 2 re-rendered with its day counter continuing from T0 in leg 1's two-decimal format, so it never steps back at the join; leg 1's end framing from the same `leg2Mat`) and `lunar_pnt_full_mission.mp4` (leg 1 + leg 2, 15-frame cross-fade, 2,985 frames, 99.5 s); no output may resolve to `lunar_pnt_mission.mp4` (`full_mission_movie:protected`) |

`leg1Mat` and `leg2Mat` are the paths of `results/lunar_pnt_leg1_study.mat` and
`results/lunar_pnt_mission_study.mat`. The leg-1 study is a numbered script like
the leg-2 one: 0 knobs, 1 inputs and identity, 2 schedule, 3 polish loop, 4
gates table and VERDICT, 5 mass and power, 6 save, 7 capture, 8 figures. Run it
from a function to keep the base workspace.

### Measured results

Run on 2026-11-11 (97 s, of which 91.8 s is the six polishes; every blocking
gate PASS, VERDICT PASS). The same inputs rerun bit-identically (three runs
compared field by field: dT0 = 0, max |dtf| = 0, max |dz8| = 0, identical
propellant, clearances and misses). Provenance: run from a working tree with
uncommitted changes on top of 1f02bf6 (`run_meta` records `1f02bf6` with the
dirty flag set).

Schedule: orient0 0 deg, GTO phase at t = 0 of 5/12, DRO slot offset
delta = 5/72, chosen from 2,974 feasible triples; maximum orientation offset
from the catalog sheet 3.319 deg. All six closed on their FIRST candidate (0
exclusions, 6 polishes). Satellites are in leg-2 numbering k; slot j is the DRO
slot the stage lands on.

| k | slot | cell (iD,iA,iS) n | departs (d) | flight (d) | arrives (d) | GTO orientation at departure (deg) | arrival miss (km, m/s) |
|---|---|---|---|---|---|---|---|
| 1 | 3 | (11,1,1) 0 | 0.1834 | 2.4342 | 2.6176 | -2.371 | 0.000413, 6.2e-06 |
| 2 | 4 | (1,2,1) 0  | 0.2568 | 2.3434 | 2.6002 | -3.319 | 0.0473, 0.000663 |
| 3 | 5 | (12,3,1) 0 | 0.2201 | 2.3869 | 2.6070 | -2.845 | 0.00196, 3.6e-05 |
| 4 | 6 | (6,4,1) 0  | 0.0000 | 2.5717 | 2.5717 | 0.000 | 0.000616, 1.2e-05 |
| 5 | 1 | (6,5,1) 0  | 0.0000 | 2.5898 | 2.5898 | 0.000 | 0.00398, 5.9e-05 |
| 6 | 2 | (6,6,1) 0  | 0.0000 | 2.6188 | 2.6188 | 0.000 | 0.000884, 1.4e-05 |

| k | min Earth (km) | min Moon (km) | second order | stage prop (kg) | stack wet (kg) | thrust (N) | jet power (kW) | dV (km/s) |
|---|---|---|---|---|---|---|---|---|
| 1 | 1068.0 | 26839.7 | pass | 120.765 | 288.879 | 9.629 | 80.74 | 9.078 |
| 2 | 349.9  | 27803.2 | pass | 112.381 | 279.238 | 9.308 | 78.04 | 8.635 |
| 3 | 520.2  | 27931.0 | pass | 116.329 | 283.779 | 9.459 | 79.31 | 8.846 |
| 4 | 6306.7 | 15837.0 | pass | 134.628 | 304.822 | 10.161 | 85.19 | 9.773 |
| 5 | 6427.6 | 841.3   | pass | 136.577 | 307.063 | 10.235 | 85.82 | 9.867 |
| 6 | 6630.5 | 27936.6 | pass | 139.739 | 310.700 | 10.357 | 86.84 | 10.018 |

Totals: stage propellant 760.418 kg, stage dry 114.063 kg, stack wet (launch
mass of the six stacks) 1,774.481 kg, largest thrust 10.357 N, largest jet power
86.84 kW. Minimum Earth altitude 349.9 km is satellite 2 at its own perigee
departure (the GTO's perigee, floor 200 km); minimum Moon altitude 841.3 km is
satellite 5 (floor 500 km). The three satellites with departure day 0.0000 leave
at the GTO phase the schedule starts on; the others wait 0.18 to 0.26 d on the
GTO. Min Earth and Moon are LOCATED minima (see Gates).

Hand-off: T0 = 2.6473 d (last polished arrival 2.6188 d, 0.0285 d earlier; the
polished T0 equals the table schedule's T0 to the printed digits, gate T1).
The flown arrival states sit on the DRO phases (k-1)/6 at T0 to within
2.04e-07 of a period (gate GH, tolerance 1e-5). Combined mission, on the
leg-1 clock: leg-2 departures 2.682 to 3.941 d, last tulip arrival day 21.446,
end of the 30-day service window day 51.446 (leg 2's own clock is shifted by
T0).

Mass and power bookkeeping (closed form, `leg1_gates` mode `'mass'`), with
f = 1 - mf/m0 of the transfer:

    mp = 150 f / (1 - f - dryFrac f),   dry = dryFrac mp,   wet = 150 + dry + mp
    thrust = 5 N wet / 150,             jet power = thrust Isp g0 / 2   (g0 = 9.80665 m/s^2)

The catalog trajectory holds for any vehicle with the same thrust-to-mass ratio
and Isp, so each stack is sized to arrive at the 150 kg leg 2 assumes. The jet
power is the ideal thrust power; no efficiency is applied.

Jacobi constants (CR3BP, same frame as the mission): DRO (size 1, period 1 ND =
4.4327 d) 3.0576; the tulip 3.1866; for reference L1 3.1883 and L2 3.1722 (the
tulip's constant lies between those two).

### Gates

The gate ids printed by the study are I1 C1 T1 P6 GA GF GO GH GE GM GW G2 GT
GD. Every gate is computed on an INDEPENDENT re-flight of the polished costates
(`fly_transfer`) and on the scheduler's context; none reads the polish's own
bookkeeping of the quantity it checks. A satellite is accepted only if the
final gates pass at the final T0 and its label k.

BLOCKING (a failure excludes the candidate, re-schedules that satellite alone
with the triple locked, and tries the next cell):

- arrival: flown miss < 100 km and < 10 m/s (GA; measured max 0.0473 km, 0.000663 m/s);
- flight admissibility (`validate_flight`) (GF);
- one GTO: the departure state recomputed from orient0, the GTO phase and tDep
  alone equals the stage's departure state to < 1e-9 ND (GO; measured 0);
- hand-off phase: the DRO phase of the flown arrival state, advanced to T0,
  equals (k-1)/6 within 1e-5 of a period, about 1.9 km along track (GH; the gate
  exists to catch a clock or a 1/6 relabelling error, and the arrival gate
  bounds the miss itself), and the flown arrival epoch matches the schedule's to 1e-9 ND;
- clearance on the LOCATED minimum: Earth altitude >= 200 km (GE), Moon altitude
  >= 500 km (GM); the range-rate roots are refined on the true dynamics, so a
  perigee that falls between samples is caught;
- stage feasibility: the closed-form mass relation has a solution (GW);
- the stored thrust direction alpha = -lambda_v/|lambda_v| re-flies the flight
  below 1 km, and the opposite sign misses (GD).

REPORTED, not blocking: second-order status (G2; pass, refuted, notReproduced or
error: `ms_tfmin` reproduces the root to 1e-6, conjugate-point test passed and
covered, min |lambda_v| > 0, min Q_mt > 0, dim S = 1; all six pass); T1 (the
hand-off epoch from the polished arrivals equals the table's); GT (T0 is >= every
arrival and a slot instant). A refuted transfer is not blocking but the next
candidate is tried first.

NOT covered:

- The second-order verdict and "minimum-time" refer to the fixed-endpoint
  problem to the slot's arrival point (a local test at the root); the
  moving-target rendezvous condition is not checked.
- The refined solutions carry no path constraint in the optimisation: Earth and
  Moon clearance is checked after the fact, on the located minimum.
- The GTO wait is Keplerian, not propagated in the CR3BP.
- No launch windows, eclipses, navigation or dispersions.
- The off-grid fallback solver (`lib/solve_gto_dro_offgrid.m`) is not built; the
  study has an empty hook (`offgridSolver`) because no satellite needed it. The
  exclusion and failed-satellite branches ran only through the test seams.
- The second-order check omits the pumpkyn `tfMin` acceptance, the lift margin,
  the H6 margin and the dense conjugate scan; the leg-2 polish also runs
  `pmp_pointwise_checks`, which leg 1 does not.
- No spacecraft-to-spacecraft proximity check.

Deliberate-failure experiment per gate (each prints and asserts; measured
values):

| gate | experiment | what happens |
|---|---|---|
| GA arrival | `z8(8)` (t_f) + 0.01 ND | miss 2017.2 km, 189.4 m/s: FAIL |
| GF flight | flight truncated at half its time, re-validated by `validate_flight` | not admissible ("did not reach t_f"), mass block infeasible, G.pass 0 |
| GO one GTO | `orient0` + 1 deg | departure state differs by 0.17 ND: FAIL (arrival and hand-off still pass) |
| GH hand-off | T0 + 0.01 ND | phase error +0.01 of a period: FAIL (oneGto still passes) |
| GH relabelling | spacecraft number 2 -> 3 | phase error -1/6 of a period: FAIL |
| GH epoch | `P.tArr` + 1e-8 ND | "flown arrival epoch differs from P.tArr by -1e-08 ND": FAIL |
| GE Earth, located | synthetic CR3BP arc through a 150 km perigee between samples 226.20 km apart | sample minimum 226.199 km PASSES the 200 km floor; located minimum 150.0000 km FAILS |
| GE / GM, re-integration stops early | propagator stopped at half | "minimum ... not located": both FAIL |
| GM Moon | floor set to measured + 1 km / - 1 km (26,886.7 / 26,884.7) | FAIL / pass |
| GW mass | `dryFrac` 1.5 (> 1.41 allowed) | denominator -0.03687, stage cannot exist: FAIL; arrival and hand-off still pass |
| G2 second order | injected error | status error, secondUnknown 1, G.pass 1 (reported, not blocking) |
| GD thrust direction | -alpha used for the re-flight (every run) | misses by 391,955 km |
| C1 schedule | a table with a slot no cell can fill | `schedule_gto_dro:noAssignment`, per-slot candidate counts, no partial assignment |
| P6 closure | seam: every candidate of slot 4 refused | "SLOT 4 FAILED: candidates exhausted (3 tried of 3; no off-grid solver configured)", P6 FAIL 5 of 6, no export, no timeline |
| final T0 | seam: T0 + 0.01 ND at the final gates | 0 of 6 accepted ("final gates failed: handoff at final T0"), no sidecar |
| exclusion | seam: slot 1 refused | locked re-schedule to another cell; T0 jumped 2.647 -> 3.386 d and every label changed, the other five transfers unchanged |
| timeline refusals | package replaced; leg-2 DRO differs; leg-2 phases differ; DRO period 0.5; time units; satellite not accepted; T0 before the last arrival; leg-1 verdict not PASS; G.pass false; leg-2 gate table failed | each throws its own `combined_timeline:<check>` identifier |
| export refusals | outFile outside `results/`; accepted entry with G.pass false; none accepted | `export_leg1_transfers:outFile` / `:notPassed` / `:noTransfers`, nothing written |

The test seams (`seamRefuse`, `seamRefute`, `seamFinalT0`) are refused when
`outDir` is the default `results/`; with a seam set the export is skipped and the
VERDICT names the seam.

### Method notes

- The GTO orientation drift sign is MEASURED, not assumed: the inertial GTO's
  orientation at leg-1 time t is `orient0 + sgn (180/pi) t` deg with sgn = -1
  (rotating-frame residual 8.9e-16 ND with sgn, 0.699 with -sgn).
  GTO period 0.0993132465 ND (10.5653 h).
- Departure pinned, arrival sliding. For a trial duration t_f the target is the
  slot's state at tDep + t_f, which slides slowly along the DRO, so the root of
  t_f - T(t_f) is well conditioned. Pinning the arrival instead would slide the
  departure around the GTO, a tenth of a revolution per 0.01 time units. This is
  the reverse of leg 2.
- Orientation walk (probe on a perigee cell and an apogee cell): 2 deg steps with
  step halving reach +10 and -25 deg at both ends of the GTO with no halving ever
  needed, about 0.3 s per degree; a perigee cell walks to +60 deg in 30 accepted
  steps (t_f +10.7 %, miss <= 0.515 km). A single 25 deg leap also converges (36
  iterations); a single +60 deg leap on the perigee cell does not. A failed or
  branch-jumping step returns `ok = false` with the orientation reached, never a
  transfer from a different departure. The real schedule needs at most 3.3 deg.
- The scheduler's tie-break: after earliest T0 and least propellant, the smallest
  maximum orientation offset from the catalog sheet. Without it the lexicographic
  rule chose orient0 = -20 deg and pushed departures to -23.3 deg against the
  25 deg bound for no benefit; with it the maximum offset is 3.3 deg and the T0
  and propellant (385.27 kg on the unpolished table) are unchanged.
- Results are deterministic: the study reproduced to the last printed digit
  across runs, so the test pins allow only cross-session solver noise (T0 1e-9 d,
  t_f 1e-6 d, propellant and wet 1e-3 kg, jet power 1e-4 kW, clearances 0.01 km,
  maximum miss 1e-4 km and 1e-6 m/s).
- Closest margins: the satellite-2 perigee departure at 349.9 km (floor 200 km)
  and the satellite-5 Moon pass at 841.3 km (floor 500 km).
- In the leg-1 movie the plume points opposite the thrust, which from about day
  1.4 (the stages brake) is ahead of the motion in the rotating frame; that is the
  physics, not a sign error. The GTO is drawn rotating about the Earth, not the
  barycentre (a barycentre rotation misses `gtoState` by 2,370 km).

### Which catalog the numbers come from

The PILOT package `orbit_transfer/GTO_DRO/results/costate_catalog_gto_dro.mat`
(git-ignored): one sheet, GTO orientation 0 deg, DRO size 1, 271 entries, 34
clean 5 N cells (`has_solution` and `conj_pass == 1` at rung 5), SHA-256 beginning
`e57d9236` (full `e57d9236f8264d026f3f5a7bcd1fb6ce9e3c698958d04ed0b78cbbb0669574d3`),
created 2026-11-11. A snapshot is kept at
`results/inputs/costate_catalog_gto_dro_pilot_e57d9236.mat`. Every leg-1 run now
keeps its package automatically (`leg1_snapshot_package`): any file in
`<outDir>/inputs` whose own SHA-256 equals the package's is the snapshot (the name
is not used), else the package is copied there as
`costate_catalog_gto_dro_<first 8 hex>.mat`, never overwriting; the path is
recorded in `run_meta.packageSnapshot`. Everything that reads saved leg-1 results
(`leg1_movie`, `combined_timeline`, the tests) finds the package with
`leg1_package_for`: the default path if its SHA-256 matches, else the recorded
snapshot, else a matching file in `inputs/`, HASHED BEFORE anything loads it.
When the fleet's package replaces the pilot, the saved results keep working
from the snapshot; `combined_timeline` refuses with `combined_timeline:identity`
only when no file anywhere matches (unless `allowStalePackage`), and the study
test applies its headline pins only to results from the pilot package (SHA-256
`e57d9236`), otherwise it reports PASS* with a SKIP note and still runs every
structural check. In this package `sheets(k).tauDRO` holds the
GTO orientation in degrees, not a DRO period; sheets are selected by
`arr_params.tau == 1`.

### Captured transfers (leg-1 sidecar)

`results/leg1_transfers_gto_dro_5N.mat` holds ONE variable,
`leg1_transfers` (six entries, one per ACCEPTED transfer), written by
`export_leg1_transfers` (atomic, refused outside this folder's `results/`).
Each entry stores: the satellite and slot, the costates `z8` = [lambda0(7); t_f],
`tf_nd`, the endpoint states `rv0`, `rvf` (the orientation and arrival phase are
off the package grid, so the phases alone do not rebuild them), the regime
(`thrustN ispS m0kg`), `orientDeg`, `sD`, `arrPhase`, the epochs (`tDep tArr`, ND
and days), `seededFrom` [iD iA iS], `nWait`, `route`, `flyKm flyVms`, and
`verdicts` and `notChecked` below. Top level: `regime package provenance note
usage`; the note is built from the stored verdicts.

CHECKED (stored in `verdicts`): the flown arrival miss, flight admissibility, one
GTO, hand-off error, LOCATED Earth and Moon minima with epochs and floors, the
second-order status with conjugate verdict, minLamV, minQmt, dimS and the
reproduction dz. The test re-flies every entry and recovers the stored
misses (0.0004133 / 0.0473 / 0.001959 / 0.0006155 / 0.003984 / 0.0008843 km).
NOT CHECKED (listed in `notChecked`): the pumpkyn `tfMin` catalog acceptance
beyond the polish's own `ms_tfmin` solve; an Earth floor inside the
optimisation; tight-flight pointwise checks, lift margin, H6 and the dense
conjugate scan; admission to any library or package. Nothing is written into
the GTO -> DRO package or the shipped library.

### Leg-1 tests

Nine tests belong to leg 1 (`run_tests` now runs 19 in all, about 17 min):
`test_leg1_inputs`, `test_leg1_package_for`, `test_polish_gto_dro`,
`test_schedule_gto_dro`, `test_leg1_gates`, `test_lunar_pnt_leg1_study`,
`test_export_leg1_transfers`, `test_combined_timeline`, `test_leg1_movie`. The
study, export and timeline tests read the saved results (and skip, passing, when
no results file exists) and pin the headline numbers above when the results come
from the pilot package. A test may return a third output, a note: `run_tests`
then reports PASS* (passed, a section degraded, the note printed), never a
failure; a false result or an error is always FAIL.

## Tests

`tests/run_tests.m` runs all 19 tests, ten of leg 2 listed here and nine of
leg 1 (see Leg-1 tests; about 17 min in all). The ten leg-2 tests (about 15 min; measured 923 s
on 2026-11-10, the six transfers certified ONCE, in the export test; the
study test's export runs with `exportFullStatus = false`): `test_mission_inputs`,
`test_schedule_rendezvous`, `test_polish_rendezvous`,
`test_assemble_histories`, `test_mission_gdop`, `test_validate_pole_series`,
`test_pole_series_stats`, `test_lunar_pnt_mission_study`,
`test_mission_movie`, `test_export_mission_transfers`. Each gate's deliberate-failure test is in the Gates
table.

The headline numbers are pinned: `test_schedule_rendezvous` (9) asserts the
unpolished optimum 18.537 d with assignment `[2 4 6 1 5 3]` and, with the
three refused candidates excluded (`[1 2 1 1 0; 2 4 5 15 0; 6 4 21 16 0]`,
rows `[k j iD iA nRev]` from the saved run's `rounds`), 18.799 d with
`[2 5 6 1 3 4]`; `test_lunar_pnt_mission_study` (A) asserts on the real
polish loop the last arrival day 18.799, the assignment `[2 5 6 1 3 4]`,
71.261 kg, exactly three rounds and no family change.
`test_pole_series_stats` checks the service statistics on hand-built series
(all NaN, an outage of known length at the start / middle / end, values
exactly 6, an Inf).

`test_lunar_pnt_mission_study` runs the script in throwaway function
workspaces: (A) the default pole-series path (3 service days, the polish
loop for real, the pins, the saved double-precision H fed back into
`mission_gdop`, the new fingerprint fields and git dirty flags), (B)
`validate = false`, (C) `maxRounds = 1`, (D) `tol = 1e-7`, (E) a reduced
`scoreRequirements = true` run (chunk 5, 3 days, checking the REUSED and
reused-exclusions disclosure), (F) reuse refused for a different
fingerprint and for a file without one, (G) two tampered cached entries
(an arrival epoch, a departure state) ignored and re-polished, (I) the
G0-FAIL stop through `g0SlotOverride`, (J) the I2 failures; the later runs
reuse (A)'s polishes. A run that reuses ONLY refusals (all six final
polishes fresh) goes through the same disclosure code as (E) but is not run
on its own (six fresh polishes, about 5 min).

`test_mission_movie` reads the saved results (it skips, passing, when no
results file exists) and writes only to a temporary folder: a 12-frame MP4
(size, frame count, no partial file left, partial renders aimed at the
default movie refused, the default movie untouched), six preview PNGs (camera
fixed in direction and target, end framing first and last, containment
margins, plume on mid-transfer and off after the last arrival, plume
direction against `+lam_v/|lam_v|` at a node and between nodes, the drawn
plume on screen), two probe PNGs at two zoom levels (the km -> pixel map
against where MATLAB draws), two PNGs of a file with a spacecraft on station
from t = 0, the provenance, and the containment failure through `zoomScale`.
It uses plain MATLAB only (the probe blobs are found with `find` and `mean`;
no Image Processing Toolbox).

`test_export_mission_transfers` reads the saved results (it skips, passing,
when no results file exists) and writes only to a temporary folder:
`catalog_identity` reproduces the saved fingerprint; the refusals of a
non-accepted transfer, a catalog-fingerprint mismatch (tampered copies) and
four output paths inside the kit (absolute new subfolder, relative, `..`);
one variable with the library's alternatives fields matched by name, class
and shape; the HEADLINE PIN: all six certified by the library's certifier
(`certifyDz` < 1e-6) and status 4, with the certifier's z8, junctions and
miss stored; `sD` / `sA` give back the polish's `rv0` / `rvf` (1e-9); the
fallback cut against a library primary; the `fullStatus = false` statuses
with the session path restored exactly; the deliberate failure of the
root-reproduction gate (one certifier run, `only = 3`); the kit round trip
by both methods below 0.05 km (the library's gate is 100 km; measured max
2.2e-3 km); the library file's MD5 unchanged.
`test_lunar_pnt_mission_study` (E) checks the fallback sidecar of a REUSED
run (entries marked reused, six printed rows, the status-run line), (B),
(C) and (I) the skip line, and (G) an export error that is printed after
the VERDICT, recorded in `run_meta.export`, and does not stop the run.

## Refusals

- CR3BP only, no ephemeris model.
- No station-keeping after arrival: the tulip slot is the ideal periodic orbit.
- No navigation filter or signal model: GDOP is geometry only.
- By default no requirement verdicts (R1-R3): the product is the south-pole
  GDOP series. The scoring exists behind `scoreRequirements`.
- The transfers use the 70 mN / 900 s / 150 kg regime, which is NOT the
  GTO -> DRO catalog's regime, so leg 2 alone is not a continuous GTO -> DRO
  -> tulip mission; leg 1 (below) flies the GTO -> DRO part in its own 5 N /
  1710 s regime and hands over at T0.
- All six spacecraft share one DRO and one clock.
- Collision and proximity between spacecraft are not checked.
- No propellant-optimal schedule: propellant is reported, not optimised.

## Measured results

Full default run, 2026-11-10 21:23 (342 s, of which 337 s is the certified
polishes; results file fingerprint z8 MD5 1cbd1390...). Every gate PASS.
Provenance: the run was made from a working tree with UNCOMMITTED changes
(this fix wave) on top of 1d43572; `run_meta` records `optimal_control
1d43572` with `optimal_controlDirty = true` (105 modified tracked paths in
the repository, most of them unrelated) and `proj7 3a7f898` with
`proj7Dirty = true`. The first run that produced these numbers (20:06) was
likewise made from uncommitted changes on top of a20d563, not from commit
a20d563 itself. Both runs give identical numbers.

**What "certified" covers for a polished transfer.** Each accepted transfer
passed, on its final z8 and departure state: the flown miss (< 100 km,
< 10 m/s, admissible flight); a separate multiple-shooting `ms_tfmin` solve (24 segments,
tolR 3e-11) that REPRODUCES the root (max |dz| < 1e-6) with residual
< 1e-8; the pointwise PMP checks (`pmp_pointwise_checks`: H, lambda_m,
adjoint, minimum principle, throttle, Q_mt >= 0); the three hypothesis gates
(min |lambda_v| > 0, min Q_mt > 0, dim S = 1); the conjugate-point junction
test of `ms_tfmin` (pass and covered, on the reproduced root); and the
pumpkyn `tfMin` acceptance (tfMin from the polished z8 returns it to 1e-6).
It does NOT include, relative to the library's own study script
(`DRO_tulip/indirect/transfer_study.m`): the lift margin, the null-residual
and H-residual thresholds, the H6 margin, or the dense conjugate-spectrum
scan. Every gate integrates through pumpkyn's `tfMinProp` (the same
integrator as the solve): the checks are integrator-consistent, not an
independent re-flight in a different integrator.

**Schedule.** Three rounds. Round 1 (unpolished optimum, last arrival day
18.537): spacecraft 1 and 2 refused, their exact-slot solution needs a
departure before clock 0. Round 2: spacecraft 6 refused for the same reason.
Round 3: all six accepted. Constellation complete on day 18.799.

Pinned-epoch caveat. Each candidate's arrival epoch is pinned to its catalog
cell's value and is never searched. All three refusals were "departure
before clock 0", which a slightly later arrival epoch could cure. So day
18.799 is the optimum of THIS formulation (primary catalog cells, whole-
revolution waits, pinned arrival epochs, refused candidates excluded), and
an UPPER bound on what a free-epoch polish could reach; the unpolished
18.537 is not a proven lower bound (a polish moves t_f, and cells off the
primary set were not searched).

| sc | slot | cell (iD,iA) | departs (d) | arrives (d) | t_f cell (d) | t_f polished (d) | flown miss (km) | propellant (kg) | dV (km/s) |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 2 | (6,2)   | 0.931 | 18.799 | 17.87558 | 17.86760 | 8.4e-05 | 12.244 | 0.7515 |
| 2 | 5 | (6,29)  | 0.112 | 18.036 | 17.85087 | 17.92312 | 9.0e-05 | 12.282 | 0.7540 |
| 3 | 6 | (16,31) | 1.293 | 17.496 | 16.20336 | 16.20276 | 6.6e-05 | 11.103 | 0.6787 |
| 4 | 1 | (15,43) | 0.366 | 17.161 | 16.79148 | 16.79444 | 2.2e-03 | 11.508 | 0.7045 |
| 5 | 3 | (17,11) | 0.035 | 17.140 | 17.13972 | 17.10510 | 2.5e-04 | 11.721 | 0.7181 |
| 6 | 4 | (22,17) | 0.227 | 18.327 | 18.14271 | 18.09995 | 4.0e-04 | 12.403 | 0.7617 |

No spacecraft waits a DRO revolution; total propellant 71.261 kg. Each
polished transfer passes every gate listed above under What "certified" covers;
none changed family.

**Captured transfers (sidecar), export of 2026-11-10.** The library's own
certifier, run on each polished transfer by `export_mission_transfers`,
gives all six status 4 (so the checks listed above as NOT included in the
polish -- lift margin, null and H residual thresholds, H6, dense conjugate
scan -- were run on these six by the export, and passed); its roots match
the polished ones to at most 4.8e-11. The kit's own flight
(`catalog_transfer`, endpoints from the kit's spline at `sD`, `sA`) lands on
the tulip by both methods; the junction flights fly the certifier's
multiple-shooting junctions and land to below 0.1 mm.

| entry (sc) | slot | seed cell | sD | sA | t_f (d) | status | certifyDz | H6 margin | lift margin | kit miss z8 (km) | kit miss junctions (km) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2 | (6,2)   | 0.210133 | 0.028763 | 17.8676 | 4 | 4.8e-11 | 9.05  | 1412  | 8.3e-05 | 2.1e-08 |
| 2 | 5 | (6,29)  | 0.192034 | 0.591976 | 17.9231 | 4 | 3.8e-11 | 9.74  | 4602  | 9.0e-05 | 4.2e-08 |
| 3 | 6 | (16,31) | 0.625135 | 0.638424 | 16.2028 | 4 | 3.4e-12 | 9.69  | 7041  | 6.6e-05 | 4.3e-08 |
| 4 | 1 | (15,43) | 0.582667 | 0.888043 | 16.7944 | 4 | 2.3e-13 | 25.42 | 40491 | 2.2e-03 | 1.9e-08 |
| 5 | 3 | (17,11) | 0.674476 | 0.219276 | 17.1051 | 4 | 2.1e-11 | 13.20 | 829   | 2.4e-04 | 5.4e-08 |
| 6 | 4 | (22,17) | 0.884647 | 0.342178 | 18.0999 | 4 | 6.2e-12 | 18.39 | 646   | 3.2e-04 | 3.4e-08 |

**Orbit handles (I2).** Against a direct propagation at 2,582 phases: tulip
handle max 0.015 km and 0.0016 m/s, DRO handle max 0.006 km and 0.00015 m/s.
Against pumpkyn's `tulipConstellation` at its own samples the slot states
differ by a median of 0.48 km and at most 2.96 km (it applies no
continuation). A 33.8 km figure seen during development is not a slot-state
error: the slot states agree with direct propagation to 0.015 km and with
`tulipConstellation` at its own samples to 2.96 km; the larger difference
appeared only when `tulipConstellation`'s non-uniform samples were
spline-resampled near perilune (the likely explanation, not proven), and
gate G0 bounds its effect on the pole GDOP (max relative difference
1.5e-4).

**Validation (G0).** On-slots south-pole series against a direct engine call
on `tulipConstellation` states over one tulip period: max relative
difference 1.5e-04 (tolerance 1e-03), no visibility mismatch away from an
edge. Moving slot 6 by 0.3 ND fails it.

**South-pole GDOP series** (300 s samples from day 0).

| measure | value |
|---|---|
| first finite GDOP (four on station) | day 18.247 |
| first GDOP < 6 | day 18.247 |
| service window | day 18.799 to 48.799 |
| time with GDOP < 6 | 99.41 percent |
| time with a fix | 100 percent |
| median / 95th percentile GDOP | 4.514 / 5.593 |
| longest stretch at GDOP >= 6 | 0.67 h |
| samples with 6 <= GDOP < 7 | 0 |

With five or six satellites above the horizon the pole GDOP is at most 5.81;
with exactly four it is at least 16.7 (measured, this constellation).

**Requirement scoring** (`scoreRequirements = true`, not the default; the run
of 2026-11-10 before the default changed): validation against the saved
Table 1 column, full paper grid, on slots, one period: the three FLOORS are
identical (1.000000 / 0.754032 / 17.987224 h). Per point the results are not
identical: 74 R2 points and 75 R3 points differ (R3 by at most 1.35 h, far
from the floor), and those per-point differences vanish when the
reference's own orbit construction (`tulipConstellation`) supplies the
states, i.e. they come from the two orbit constructions, not from the
scoring. 30-day service floors R1 1.000, R2 0.754, R3 17.695 h.

**Movie.** `results/lunar_pnt_mission.mp4`, 70 s, 1280 x 720 (revision 3:
zoom widens to 1.75 times the end framing at day 12.5 so all six spacecraft
stay in view; smallest margin to the scene edge 16.5 px).

**What these numbers do not show.** The validation shares proj7's GDOP
engine with what it checks (it catches wiring on this side, not an engine
error). The schedule optimum is over primary catalog cells with whole-
revolution waits, pinned arrival epochs and the exclusions made by refused
polishes; alternative cells and free arrival epochs were not searched (see
the pinned-epoch caveat). The certification is integrator-consistent
(pumpkyn `tfMinProp` throughout) and omits the study script's lift, null-
and H-residual, H6 and dense conjugate-spectrum checks.
