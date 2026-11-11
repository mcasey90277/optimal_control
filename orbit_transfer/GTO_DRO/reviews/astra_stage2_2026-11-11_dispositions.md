# Dispositions: Astra review of the stage 2 recipe, 2026-11-11

Review: `astra_stage2_2026-11-11.md` (prompt beside it). Astra saw a 10 KB
summary, not the code, so several findings ask for things that exist. Each
finding is marked ACCEPTED (and when), CHANGED (accepted in another form),
ALREADY COVERED, or REFUSED, with the reason. Numbers are measured on the
pilot sheet `gto_dro_o000_tau1` unless stated.

## Measurements made in response (controller, shared session, 2026-11-11)

- **Is the Earth floor active in the direct solution of a demoted cell?**
  Cell (11,5), 15 N, re-solved cold exactly as the engine does: direct t_f
  0.32254 ND; minimum Earth altitude 243.9 km at the nodes, 199.9 km at the
  Hermite-Simpson midpoints (the floor is 200 km, so it IS active, at a
  midpoint: an isolated contact, not an arc). The direct control flown open
  loop with ode113 (1e-11) through the first perigee stays within 7 km of the
  node states and bottoms at 195.5 km. The refined entry of the same cell has
  t_f 0.32144 ND and bottoms at -1,541 km. So for this cell the direct
  solution is accurate and admissible, the floor constrains it, and the
  unconstrained refinement leaves it by 1,740 km at perigee for a t_f gain of
  0.3 percent. One cell; the fleet will keep every direct solution so the
  rest can be measured (see C1).
- **Side finding.** The direct solver's `altEarthActive` flag looks at the
  nodes only (it read false above, with the midpoint on the floor). Deferred
  minor in `costate_common`; the new sidecar records the node-and-midpoint
  minimum.
- **Is the clearance instrument (minimum over the flight's output samples)
  good enough, and does the Moon floor hold on the refined entries?** All
  271 accepted entries re-flown and their positions spline-interpolated on
  400,001 points (clamped with the flown end velocities): Earth minimum
  240.4 km by samples and by dense interpolation, largest difference over
  all entries 0.1 km, none below 200 km. Moon minimum 533.5 km (next 1,780
  km), none below the 500 km floor. Flights carry a median of 2,337 samples.

## Critical

**C1. The direct solve and the refinement solve different problems.**
ACCEPTED, in two parts.
- Before the fleet: the engine gains `keepDirect` (default off; the fleet
  front door turns it on) and writes every successful direct solution to
  `<tag>_direct.mat` with its node-and-midpoint Earth and Moon minima. The
  40 h of compute then supports either architecture Astra names without
  re-solving: (A) the constrained direct solution as the product for
  floor-active cells, or re-seeding another admissible branch.
- Wording: a demotion condemns that trajectory, not the cell. The README
  bullet that said "demotion is the final verdict" is corrected.
- After the fleet (not before): the decision between A and B for the
  floor-active cells, mesh refinement at perigee and floor-perturbation
  sensitivity. A constrained indirect formulation (boundary conditions of a
  second-order state constraint) is out of scope for this catalog.

**C2. Continuous path feasibility is not established; the Moon floor is not
swept.** CHANGED. Measured above: the sample minimum agrees with a dense
interpolation to 0.1 km on every entry and the Moon floor holds on all 271,
so no pilot verdict changes. Accepted as a chain addition AFTER the fleet
(post-processing, no recompute): a Moon sweep beside the Earth sweep, and
the located minimum (root of (r - r_b) . v) in place of the sample minimum.
The chain already re-sweeps every entry a repair changed (stage 2b, keyed
to each entry's z8), which Astra could not see.

**C3. The acceptance gates do not establish a rendezvous PMP solution.**
ALREADY COVERED in substance; Astra is right about what the acceptance
number is. Acceptance is pumpkyn `tfMin` started at z8, whose residual is
the full terminal state (position and velocity) plus the free-time
condition. The pass criterion is |dz| <= 1e-6 (how far tfMin moves z8), a
reproducibility measure as Astra says, but the engine also records the
multiple-shooting residual norm of every entry (`RES`). On the 271 pilot
entries: RES at most 8.0e-11 (median 3.6e-14), |dz| at most 7.7e-9, flown
miss at most 0.23 km (median 0.001 km). So the 100 km and 1e-6 thresholds
are screens that no accepted entry comes near; the fleet findings will
report the three maxima rather than the thresholds. The mass-costate and
Hamiltonian conditions are checked per entry by the study script's N-gates
on demand, not at catalog scale. The arrival state is a fixed point in the
rotating frame (stage 1 README). No code change.

**C4. Conjugate-refuted entries cannot ship as minimum-time solutions.**
REFUSED as a removal, ACCEPTED as a labelling rule. The catalog's standing
convention since the 0.5 N extension is that refuted entries ship with
`conj_pass = 0` and the README tells the consumer to consult it: they are
certified extremals and useful seeds, and removing them hides where the
weak cells are. The count that may be called "not refuted" is 211, and the
pilot findings say "pass / refuted", never "optimal". Making every rung a
repair queue is accepted for after the fleet (I2). The hypothesis gates
(what `dim S = 1` means, the |lambda_v| scale) are documented in
`DRO_tulip/doc/mintime_second_order_audit.tex`.

**C5. "Faster and loses no lower rung" can keep an inferior top rung.**
ACCEPTED, after the fleet. One pilot case: (3,6), 0.31064 -> 0.27232 ND
rejected because the re-solved ladder lost rungs. The fix is in the repair
tool only (accept the better top rung, keep the previous lower rungs as
entries of their own branch), so it does not touch fleet compute. Recorded
as a repair-tool task.

## Important

**I1. A large direct miss should be a flag, not ignored.** ACCEPTED via
`keepDirect`: with the direct t_f and trajectory kept, the direct-to-refined
change in t_f and in perigee altitude is measurable per entry after the
fleet. No pre-gate is reinstated.

**I2. Basin search too narrow; lower rungs not repaired.** ACCEPTED, after
the fleet: lower-rung refutations become a repair queue (re-solve from the
rung below, then from neighbouring phases), and the five same-solution
returns get neighbouring-phase seeds. Astra's caution that the growth of
the refuted share down the ladder is on unmatched populations is correct
and the findings do not claim a trend beyond the counts.

**I3. Easy phase-continuation repairs.** ALREADY COVERED. The repair tool
changes the endpoint (the cell's own departure state from the sheet's meta)
and the engine re-solves, refines and accepts against that state; the
repaired perigee entries passed the Earth re-sweep and the audit replays
entries from the sheet's own keys. Approaching a phase from both sides is a
good branch check and is noted for the fleet findings.

**I4. The deep-rung wall is a solver observation.** ACCEPTED as wording.
The pilot text says where the ladder stops on six cells, not that 2 N is
unreachable. Junction-state storage for deep rungs is the known route (the
70 mN library ships junctions for this reason) and is out of scope here.

**I5. Mass model and always-on premise.** ALREADY COVERED. Thrust magnitude
is fixed by construction (`thrLock`), stated in the stage 1 README as the
min-time model. Mass depletion is guarded at the source
(`mintime_prop_seg:massDepleted`), and final mass is derivable per entry
from t_f (constant mass flow).

**I6. Exhaustive invalidation, not a three-entry audit.** ALREADY COVERED.
Flown miss and acceptance are computed by the engine on every entry; Earth
clearance, conjugate test and hypothesis gates run on every packaged entry;
the deliverable preflight refuses a catalog whose packaged z8 differs from
the sheets' or whose Earth sidecar does not cover the current z8. The
three-entry audit is an end-to-end replay of the shipped file through the
user's pickers, on top of those.

## Minor

**M1. tf0 statement.** CHANGED: the "within 0.02" is about the median
(0.3098 vs 0.30), which is the plan's rule; the README sentence is reworded
so it cannot be read as the range. A multi-tf0 comparison on matched cells
is not run before the fleet: the front door already uses a two-point top-rung
multistart, and basin selection is the repair stage's job.

**M2. "Solved" conflates categories.** ACCEPTED as wording in the fleet
findings (direct converged / accepted / floor-clear / not refuted are
reported separately; the sidecar adds the first).

**M3. One sheet does not calibrate the fleet.** ACCEPTED in the form the
plan already has: the 180-degree sheets run first, and the fleet is
resumable, so a recipe failure there costs the first sheets, not 40 h. A
separate stratified preflight is refused: it would spend the same hours on
cells the fleet solves anyway.

## What changes before the fleet

1. `keepDirect` in the engine, on in the fleet front door (reviewed task).
2. README wording (C1, M1).

Everything else is post-processing on the sheets and waits for them.
