# Lunar PNT mission, leg 1: six transfer stages from one GTO to the six DRO slots

Date: 2026-11-11. Extends `2026-11-10-lunar-pnt-mission-design.md` (leg 2:
DRO to tulip, built, in `orbit_transfer/lunar_pnt_mission/`).

## Intent

The mission sim starts with six satellites already spread at sixths of a
DRO. Mike wants the mission analysis to start one step earlier: the six
satellites leave a GTO, reach the DRO, and then fly the existing 70 mN leg
to the tulip constellation that serves the lunar south pole. The result is
one timeline from GTO to a working PNT constellation, and a movie of the
six GTO to DRO transfers.

Decided by Mike, 2026-11-11:

1. **Two-stage story.** Leg 1 is flown by a transfer stage in the GTO to
   DRO catalog's regime; leg 2 is the satellite's own 70 mN thruster, as
   built.
2. **Six stages, one per satellite, from one inertially fixed GTO, with
   each departure polished** (seeded from the catalog, re-solved at its
   true departure state). Not "catalog entries as they are".
3. **The 5 N rung** (the cheapest in the catalog), with the propellant and
   power caveat stated plainly wherever numbers are shown.

## What is already true (measured 2026-11-11)

- The mission DRO (`get_family_orbit('dro', tau = 1)`, period 4.433 d,
  Jacobi constant 3.0576) is the catalog's DRO size 1. The tulip's Jacobi
  constant is 3.1866.
- The catalog's six arrival phases (0, 1/6, ..., 5/6) equal
  `mission_inputs().droPhase0`.
- Pilot package `orbit_transfer/GTO_DRO/results/costate_catalog_gto_dro.mat`
  (git-ignored): one sheet, orientation 0 deg, DRO size 1, 12 departure
  phases by 6 arrival phases by rungs [15 12 10 7 5] N, 150 kg, Isp 1710 s.
  At 5 N, 34 cells are clean (`has_solution` and `conj_pass == 1`; packaged
  entries are Earth-clear). Clean 5 N flight times 2.35 to 3.07 d,
  propellant 61 to 79 kg of 150 kg (delta-V 8.7 to 12.6 km/s). The
  departure rows at phases 0, 5/12, 1/2, 5/6 and 11/12 are well covered;
  phases 1/6 to 1/3 and 7/12 to 3/4 are sparse.
- A fleet run is filling the other orientations (90, 180, 270 deg) for
  every DRO size; it takes about 40 h from 2026-11-11 04:16 PST.

## Formulation

**Clock.** Leg-1 clock time t = 0 is the earliest instant any stage may
depart. T0 is the leg-1 time at which the DRO constellation is declared
complete; leg 2's own clock (`t = 0` in the existing study) is leg-1 time
T0. Mission time for the combined products is leg-1 time.

**The GTO.** One Kepler orbit (350 x 35,786 km, the catalog's default
elements), fixed in inertial space, treated as the catalog treats it: an
algebraic locus, not a CR3BP propagation. In the rotating frame its
perigee direction turns at the frame rate (1 rad per ND time unit, 12.9 deg per day, about 5.7 deg per GTO revolution), so

    orient(t)   = orient0 + sgn * (180/pi) * t        [deg, t in ND]
    gtoPhase(t) = mod(gtoPhase0 + t/tauG, 1)

where tauG is the Kepler period in ND (about 10.6 h) and all six stages
ride the same orbit at the same phase until each departs (one launch). The
sign `sgn` is fixed by a test against `pumpkyn.cr3bp.fromPCI` at two
epochs, not assumed. `orient0` and `gtoPhase0` are mission knobs (launch
date and launch time).

**The DRO slots.** Slot j at leg-1 time t is at DRO phase
`mod(droPhase0(j) + (t - T0)/tauD, 1)`. T0 is free (a satellite may coast
in the DRO), so the six slots are six points one sixth of a period apart
with a free common offset.

**One stage's transfer.** Stage k leaves the GTO at tDep_k from the state
the GTO has there (phase and orientation both follow from tDep_k), flies a
minimum-time, thrust-always-on, 5 N per 150 kg transfer, and arrives at
tArr_k = tDep_k + tf_k exactly on its slot's state. After arrival the
satellite coasts on the DRO until T0.

**Schedule, then polish** (the pattern of leg 2):

- *Schedule* on the catalog table: choose for each satellite a clean 5 N
  cell (iD, iA), a number of whole GTO revolutions to wait, and a slot, so
  that every cell's arrival phase is within a tolerance of its slot and the
  constellation is complete as early as possible. The objective is T0,
  taken equal to the last arrival (ties to less propellant); the slot
  offset follows from T0 and is not a separate choice. Departures are
  taken at catalog departure phases (tDep is chosen so the GTO phase is
  exactly a grid phase), so only the orientation is off-grid. `orient0`, `gtoPhase0` and the slot
  offset are searched on a grid as part of the schedule, with every
  departure's orientation kept within a bound of an available catalog
  orientation. `schedule_rendezvous` is the model (candidate rule,
  assignment, tie-breaks) but does not fit unchanged (its waits are DRO
  revolutions and its table is two-dimensional): leg 1 gets a sibling.
- *Polish* each scheduled transfer with the DEPARTURE epoch pinned (the
  reverse of leg 2, my call): the departure state is fixed on the GTO, and
  for a trial duration tf the target is the slot's state at tDep + tf, so
  the target slides slowly along the DRO and the root of tf - T(tf) is
  well conditioned. Pinning the arrival instead would slide the departure
  around the GTO, a tenth of a revolution per 0.01 ND, and the root-find
  would be steep and multi-rooted. Two steps: first walk the catalog
  cell's solution from the sheet orientation to the true orientation at
  the cell's own arrival phase by continuation, then root-find the
  arrival. Seeds are single flights of the stored costates (the package
  holds no junction states; its entries re-fly to 0.23 km). The plan's
  first task is a probe of this polish on one cell before anything else
  is built on it.
- *Fallback* when continuation from the catalog fails for a satellite:
  solve that satellite's cell at its true orientation through the shared
  ladder engine with the FLEET's recipe (no direct gate, the 300 s fence
  on the refinement, 15 N down to 5 N), from a new wrapper in this folder.
  The existing `gto_dro_entry` is not that recipe (it keeps the direct
  gate and has no fence) and is not used. Then polish from that entry. A satellite that neither route closes is reported as failed with
  its reason; the study does not substitute an unpolished catalog entry
  silently.

## Gates (each with a deliberate-failure experiment, as in leg 2)

1. **Arrival.** Each polished transfer, re-flown from its stored costates,
   meets its slot state at tArr within a tolerance the builder sets from
   what is measured (leg 2's polish accepts at 100 km and 10 m/s and
   achieves far less); the achieved maxima are reported, not the
   threshold.
2. **Slot phases at T0.** The six satellites, coasted on the DRO to T0,
   sit at `droPhase0` within the same tolerance: this is the hand-off to
   leg 2 and the reason leg 2's saved results can be reused unchanged.
3. **One GTO.** Every departure state equals the state of the single
   inertial GTO at that stage's departure epoch (orientation and phase
   recomputed independently of the polish).
4. **Clearance.** Minimum Earth altitude at least 200 km and minimum Moon
   altitude at least 500 km over each flown transfer, using the located
   minimum, not node samples. The refined solutions carry no path
   constraint and the pilot showed they can dive through the Earth at the
   first perigee pass, so a failure here is a real outcome: the satellite
   is re-scheduled onto another cell or reported failed.
5. **Second-order status.** Conjugate test and hypothesis gates on each
   polished transfer with the catalog's own instruments; a refuted
   transfer is reported as refuted (it still flies) and the scheduler
   prefers another cell if one exists.
6. **Leg-2 compatibility.** The saved leg-2 results' DRO, `droPhase0` and
   catalog fingerprint match what leg 1 delivers; otherwise the combined
   products refuse to build.

## Mass and power bookkeeping (the caveat, made quantitative)

The catalog trajectory holds for any vehicle with the same thrust-to-mass
ratio and Isp. Each stack is sized so the satellite arrives at 150 kg, the
mass leg 2 assumes:

    payload        = 150 kg
    stage dry mass = dryFrac * stage propellant        (knob, default 0.15)
    stack wet mass = (payload + stage dry) / (mf/m0 of the transfer)
                     (the two relations solved together, closed form)
    stack thrust   = 5 N * stack wet / 150 kg
    jet power      = thrust * Isp * g0 / 2

The study prints, per satellite and in total: flight time, propellant,
stack wet mass, thrust, jet power and delta-V, followed by one sentence
that these are minimum-time, engine-always-on transfers in the catalog's
numerical regime and not a vehicle design. The README says the same.
`dryFrac = 0.15` is my assumption, not Mike's.

## Products

All new code under `orbit_transfer/lunar_pnt_mission/`; results under its
git-ignored `results/`.

- **Study script** for leg 1 in the house study-script style (numbered
  sections, each gate computed in the script with PASS or FAIL), saving a
  results file with, per satellite: cell and slot, tDep, tArr, tf,
  orientation and GTO phase at departure, costates, flown history,
  clearances, verdicts, mass and power; plus T0, the knobs, and run
  metadata with the catalog fingerprint.
- **Combined timeline**: one figure (per satellite: GTO wait, leg 1, DRO
  coast, DRO wait, leg 2, on station) and a table, built from the leg-1
  results and the saved leg-2 results with leg-2 times shifted by T0.
- **Captured transfers**: the six polished transfers written in the
  fields a GTO to DRO sheet entry has (costates, t_f, endpoint states,
  thrust regime, verdicts) plus the off-grid orientation, phases and
  epochs, to a sidecar in `results/`, with the status each one actually
  earned and the identity of the package it was seeded from (file hash,
  creation stamp, entry count). No new certifier or schema is built
  (standing rule, 2026-11-10). Nothing is written into the GTO
  to DRO catalog.
- **Movie of leg 1**, rendered from the saved results, in the style of the
  existing mission movie: dark frame, textured Moon, Earth drawn, GTO and
  DRO outlines on in every frame (the GTO outline turns with the frame as
  time advances), yellow-orange thrust plume behind each thrusting stage,
  day counter, the camera zooming out and in so every satellite stays in
  view (wide for the departures near Earth, closing on the Moon as the
  slots fill), ending on the framing where the existing leg-2 movie
  begins. A bottom strip shows satellites in the DRO versus mission day.
- **Full-mission movie**: the leg-1 movie followed by the existing leg-2
  movie in one file, with the day counter continuous across the join
  (`mission_movie` gains a day-offset option; its default output is
  unchanged).
- README section for leg 1 with the measured numbers and the caveat.

## Non-goals

- A propellant-realistic leg 1 (minimum-fuel or 70 mN GTO to DRO).
- Propagating the GTO wait in the CR3BP, launch windows, eclipses,
  navigation or dispersions.
- A different DRO size (the near-equal-Jacobi size 0.5 would need a new
  DRO to tulip library).
- Re-solving or re-rendering leg 2.

## Anti-goals

- Presenting leg 1's mass or power as a vehicle design.
- Six satellites leaving six different GTOs, or a catalog entry shown as
  flown from a state it was not solved for.
- A transfer in the movie that passes below the Earth floor.

## Constraints

- While the fleet runs: no edits under `orbit_transfer/costate_common`,
  `orbit_transfer/GTO_DRO/lib` or `orbit_transfer/GTO_DRO/catalog`, and
  nothing written under `orbit_transfer/GTO_DRO/results`. Calling those
  functions and reading those results is fine.
- Done means done on the pilot package: the study runs today on the one
  0 deg sheet. If fewer than six satellites close on it, the study
  reports which and why, and is re-run when the fleet's package lands; it
  does not wait for the fleet.
- The catalog path is a knob. The study uses every DRO size 1 sheet the
  package holds, so it runs on the pilot package today and gains the other
  orientations when the fleet's package replaces it; the run metadata
  records which.
- Existing leg-2 tests stay green; leg-2 saved results are read, never
  rewritten.
- House rules as in the leg-2 plan's rules file (MATLAB style, tests
  first, a deliberate-failure experiment per gate, no commits by
  subagents, MATLAB through the shared session).

## Left to the builder

- Whether `polish_rendezvous` is generalised (endpoint handles) or leg 1
  gets a sibling; either way the leg-2 tests stay green.
- The continuation path in (orientation, phase) and its step rule.
- The schedule's grid over `orient0`, `gtoPhase0` and the slot offset, the
  phase tolerance and the orientation bound (start from 25 deg).
- The arrival tolerance in gate 1, set from what is measured.
- The movie's time map and zoom schedule, following `mission_movie`.
- How the two movies are joined (one renderer or a concatenation).

## Known risks

- Continuation over tens of degrees of orientation moves the departure
  state by thousands of km near Earth; the fallback exists for that.
- With one orientation sheet the schedule's flight-time table is exact
  only at 0 deg; the polish corrects it, and a polished flight time far
  from the table can break a slot assignment (then re-schedule).
- Earth-floor failures after the polish (gate 4) may remove candidate
  cells; 34 clean cells and the free knobs give room, but a satellite can
  end as failed and the study must say so.
