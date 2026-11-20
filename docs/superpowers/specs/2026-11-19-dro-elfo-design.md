# DRO to ELFO minimum-time transfers: design

Date: 2026-11-19. Owner: Mike. Status: for approval.

## Why

We have a mature DRO to tulip line (catalog down to 0.5 N, a 70 mN phase
library, study and chain scripts). The question is how much harder or easier
the same transfer is when the target is an elliptical lunar frozen orbit
(ELFO), the orbit the IM lunar data network uses. Answering it needs the same
machinery pointed at a new target, a measured comparison, and then a costate
library for the new pair. Along the way, anything written for this pair that
is not specific to it goes into the shared library so the next pair is cheaper.

## Decisions Mike made (2026-11-19)

1. The default ELFO is the IM shape from `proj7/pipelines/gdop/im_elfo_optimum`
   (semi-major axis 12,000 km, eccentricity 0.69, inclination 56.5 deg,
   argument of perilune 90 deg). The script also accepts raw elements.
2. The catalog thrust ladder is run first (150 kg, Isp 1710 s, 15 N down
   toward 0.5 N), because a matched DRO to tulip data set exists there. 70 mN
   is tried afterwards on one anchor cell.
3. The ELFO arrival is a fixed curve of states at a given plane orientation,
   the way the GTO departure is treated in the GTO catalogs.
4. Minimum time first.
5. pumpkyn and pumpkynPie supply the astrodynamics; the ELFO definition comes
   from proj7; the user script is written in the style of
   `SDA/run/run_ekf_tulipVelfo_cr3bp_resacc.m`.

## What the ELFO is here

A Keplerian orbit about the Moon, converted to the rotating frame with
`pumpkyn.cr3bp.fromOrb` exactly as `im_elfo_states` does. It is not a periodic
orbit of the CR3BP: propagated one Kepler period from perilune it misses its
start by 1,158 km in the rotating frame (measured), and its plane turns once a
month relative to the Earth-Moon line.

So the arrival locus is algebraic (no propagation): the osculating states over
one revolution at a fixed plane orientation `raanDeg` (node angle measured in
the rotating frame at the arrival instant), parametrized by time since perilune
over the Kepler period. "Arrive on the ELFO" means "match that osculating
state"; what the satellite does afterwards under Earth's perturbation is not
part of this problem. A library for the pair is a stack of tori: (DRO phase x
ELFO phase), one sheet per (DRO size, `raanDeg`).

Measured: period 0.3080 ND (1.365 d), radius from the Moon's centre 3,720 to
20,280 km.

## Stage 0 result (done, scratch only)

A scratch copy of `get_family_orbit` with the ELFO added, run through the
unchanged ladder engine with the GTO to DRO recipe (`thrLock`, `tf0 = 0.3`,
fence 300 s). DRO size 1, `raanDeg` 0, DRO phase 0, ELFO phase 0 (perilune)
and 0.5 (apolune), rungs 15 to 1 N: 18 of 18 accepted in 4.4 minutes, flown
miss below 1e-5 km. Transfer times 0.090 ND at 15 N to 0.52 ND at 1 N (0.40 to
2.3 d), all well under one revolution. The shipped DRO to tulip catalog at DRO
size 1 has medians 0.093 to 0.106 ND at 15 N and 0.43 to 0.49 ND at 1 N. So in
the catalog regime the two targets look alike on two cells. Two cells are not
a comparison; stage 2 is.

## Goals

1. A user script: set which DRO, which ELFO, the two phases and the thrust;
   it solves the minimum-time transfer, prints every verdict, saves and plots.
2. A measured difficulty comparison, DRO to ELFO against DRO to tulip, on a
   matched experiment.
3. Fleet and chain scripts that build the phase tori and package a DRO to ELFO
   costate catalog in the existing schema.
4. Generic pieces promoted to the shared library, each with a test and an
   equivalence check against the campaign code it generalizes.

## Non-goals

- Minimum fuel or minimum energy.
- A 70 mN phase library for this pair (one anchor attempt only, as a measurement).
- An ELFO that moves during the transfer (rendezvous with a drifting plane).
- A multi-satellite mission or constellation deployment.
- Shipping a deliverable to Darin (separate, gated step).
- Restyling or restructuring DRO_tulip.

## Anti-goals

- A second ELFO definition. One builder, checked against proj7.
- A number reported without its verdict.
- A shared-library function that only this campaign can call.
- Changing the default behaviour of the ladder engine or of any shipped
  catalog, deliverable kit, pumpkyn, pumpkynPie or proj7 file.

## Parts

### A. Shared library (`orbit_transfer/costate_common/`)

- **`'elfo'` family in `get_family_orbit`.** Parameters `.oe` ([sma_km ecc
  inc_deg argp_deg], default the IM shape), `.raanDeg` (default 0). Algebraic
  locus like `'gto'`: sampled uniform in eccentric anomaly, `tau` is time
  since perilune in ND over one Kepler period, a locus parameter and never an
  epoch offset. Oracle test: for each of the five IM satellites the locus
  state at that satellite's mean anomaly and RAAN equals the `t = 0` state
  `im_elfo_states` builds, and the radius limits are 3,720 and 20,280 km.
- **`pair_transfer`.** One transfer for any orbit pair the family provider
  knows: `T = pair_transfer(dep, arr, sD, sA, thrustN, opts)` with `dep` and
  `arr` structs of `.family` and `.params`. It runs the ladder engine on one
  cell from the top rung down to `thrustN`, re-flies the accepted entry, and
  attaches the verdicts: flown miss, tfMin acceptance, flight admissibility,
  located minimum Moon and Earth altitude, conjugate test, hypothesis gates.
  It returns the flight (states, costates, thrust direction, mass) for
  plotting and the entry in library format, always with its junction
  states (below 1 N a bare `z8` is not a usable entry). The rung list is an
  option (default: the standard rungs above `thrustN`, then `thrustN`), with
  no hard thrust limits; the engine's recipe options (`tf0`, `thrLock`,
  `directGateKm`, floors, fence) pass through. It lives in the shared
  library and puts no campaign folder on the path. It generalizes
  `GTO_DRO/lib/gto_dro_entry`; equivalence check: with that wrapper's recipe
  it reproduces the wrapper's demo cell (same `z8` to the acceptance
  tolerance). On one DRO to tulip cell its result must pass every gate, and
  its transfer time is reported beside the shipped catalog's value, not
  asserted equal (the shipped catalog came from an older recipe and basin
  selection is not reproducible). `gto_dro_entry` itself is left alone.
- **`plot_pair_transfer`.** One figure from a `pair_transfer` result: the
  transfer and both orbits in the rotating frame, a Moon-centred view, thrust
  direction, distances to Moon and Earth, mass. Any pair.
- **`transfer_difficulty`.** A scorecard computed from any sheet the ladder
  engine wrote (and its direct sidecar when present): per rung, share of cells
  accepted, attempts, wall time, transfer time statistics, revolutions about
  the target's primary, delta-V, phase sensitivity of transfer time, share
  refuted when verdicts are present. One table per sheet; same code for every
  pair.

`gto_elfo_endpoints` (older GTO to ELFO campaign) builds the same ELFO its own
way. Pointing it at the new builder is done only if its outputs are bit-equal
afterwards; otherwise it is left alone and the difference is recorded.

### B. Campaign folder (`orbit_transfer/DRO_ELFO/`)

- **`dro_elfo_transfer.m`**, the user script. Style of the SDA example: a
  header paragraph; a `%% knobs` block with units and commented-out
  alternates (DRO size; ELFO from `im_elfo_optimum` + satellite number, or raw
  elements; `raanDeg`; DRO phase; ELFO phase; thrust; mass and Isp; plot and
  movie switches); then one library call per `%%` section with the output
  fields listed in comments; one-line `fprintf` summaries; a verdict table; a
  `.mat` in `results/`; the plot. A commented-out knob line switches the
  target to the tulip, since `pair_transfer` does not care.
- **`compare_dro_targets.m`**, the stage 2 script: the matched experiment and
  the two scorecards side by side.
- **`run_dro_elfo_torus.m`** and **`run_dro_elfo_batched.sh`**, the fleet
  front door and its resumable batched driver, copied from the GTO to DRO
  pair with the parameter block changed.
- **`build_dro_elfo_catalog.m`**, the chain: census, Moon and Earth clearance
  on the located minimum, repairs, package, conjugate test, hypothesis gates,
  audit, pictures. Copied from `build_gto_dro_catalog` with the parameter
  block changed. The deliverable step is not built.
- `setup_dro_elfo_paths.m`, `README.md`, `tests/run_tests.m`.

### C. The matched experiment (stage 2)

Same DRO (size 1), same 6 x 6 phase grid, same cold recipe, same rungs
[15 12 10 7 5 3 2 1.5 1 0.75 0.5] N, fresh solves for both targets (the
shipped tulip catalog was built with an older recipe, so it is a reference,
not the control). ELFO at `raanDeg` 0 and 90. Then a deep probe on a few
full-ladder cells of each target below 0.5 N, down to 70 mN, to locate where
each ladder stops. Reported per target: the scorecard, where the ladder
stops, and the share refuted by the conjugate test. The finding is stated as
measured, including "no difference" if that is what comes out.

### D. The library (stage 3)

Pilot: one sheet (DRO size 1, `raanDeg` 0), 12 x 12 phases, rungs 15 to 0.5 N,
then the chain on it, then an Astra review of the recipe before the fleet.
Fleet (my call, to be revisited after the pilot): DRO size 1 only, `raanDeg`
every 30 deg, 12 sheets. Other DRO sizes are a later extension.

### E. Sheet identity in the catalog

The packager keys a sheet by departure period, arrival period and `Np`. Every
ELFO sheet has the same Kepler period, so those keys cannot tell plane
orientations apart. The GTO catalogs solved the same problem by storing the
orientation in the field named `tauDRO`, which later caused a mislabelling
bug. Here nothing is overloaded: `tau_dep` is the DRO period, `tau_arr` the
Kepler period, the orientation is `arr_params.raanDeg`, and a new additive
sheet field `sheet_key` (struct: `depSize`, `raanDeg`) is written by the
campaign's packaging step. A campaign picker `dro_elfo_pick` selects by
(DRO size, `raanDeg`, phases, thrust). The shared schema validator and the
shipped pickers are not changed; the tulip-named schema fields (`Np`, `pm`,
`period_tulip_nd`) are filled as the packager already does for non-tulip
arrivals (NaN, NaN, arrival period).

## Constraints

- House MATLAB style for functions (header with sizes, no `i`/`j` indices, no
  `norm`, no `%#ok`, ASCII); the user script follows the SDA example's style.
- No dependence on the SDA folder. proj7 is read (ELFO shape), never written.
- Every solved transfer is captured in library-entry format (standing rule).
- Long runs follow the campaign discipline (file logs, nohup, resumable,
  attempt counters on disk, verdict files).
- Every gate gets a deliberate-failure experiment.

## Left to the builder

- Internal structure of `pair_transfer` (how much of `gto_dro_entry` is
  lifted versus called).
- How revolutions are counted in the scorecard.
- Figure layout.
- Whether the chain is a copy or a parametrized shared function; a copy is
  the standing convention.

## My calls (not Mike's)

- `raanDeg` is measured in the rotating frame at the arrival instant.
- The fleet is DRO size 1 x 12 orientations.
- Stage 2 uses fresh tulip solves as the control.
- `gto_elfo_endpoints` is touched only if bit-equal.
- Sheet identity by an additive `sheet_key` field and a campaign picker, no
  overloaded fields (section E).
- The tulip equivalence check reports the shipped value instead of asserting it.

## Builder check (2026-11-19)

A reviewer read this against the code. Settled here: sheet identity (E),
`pair_transfer` thrust range and junction states, the tulip check, no
campaign dependency in the shared function. Confirmed by the stage 0 run: the
unchanged engine accepts an algebraic arrival locus (`phase_state` and the
period fields work with the Kepler period). To check while building: the
units of the anomaly passed to `fromOrb` (the oracle test against
`im_elfo_states` settles it) and whether `phase_sensitivity` assumes a
uniformly sampled periodic orbit.
