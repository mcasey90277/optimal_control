You are reviewing the recipe for a catalog of minimum-time low-thrust transfers in the Earth-Moon CR3BP, GTO -> distant retrograde orbit (DRO), before a ~40-hour compute fleet is launched. Be a skeptical optimal-control referee: say what is wrong, what is unproven, and what you would change before spending the compute. Rank findings (Critical / Important / Minor), each with the concrete failure it causes and a concrete fix or measurement. Do not praise.

## Problem and pipeline

- Dynamics: CR3BP rotating frame, state (r, v, m), thrust magnitude T fixed (min-time, always-on), direction free. 150 kg, Isp 1710 s, thrust rungs [15 12 10 7 5] N. Rendezvous: fixed departure state on a GTO locus (keyed by orientation and time since perigee; 12 departure phases) to a fixed state on a planar DRO (6 arrival phases). 16 sheets (4 GTO orientations x 4 DRO sizes), 72 cells x 5 rungs each.
- Per cell, a thrust ladder. Top rung (15 N): cold-start direct transcription (CasADi/IPOPT, N = 400 uniform-in-time nodes, free t_f, path constraints AT THE NODES ONLY: Moon altitude >= 500 km, Earth altitude >= 200 km). Lower rungs warm-start from the rung above.
- Each direct solution is turned into costates (covector harvest from the NLP duals) and REFINED by indirect multiple shooting on the min-time PMP two-point boundary value problem (`ms_tfmin`: unknowns are junction states/costates and t_f; NO path constraint of any kind). The result is an 8-vector z8 = [lambda_r0; lambda_v0; lambda_m0; t_f].
- Gates on the refined z8: single-shooting flight from z8 must hit the target within 100 km; a third-party single-shooting min-time solver (pumpkyn tfMin) started at z8 must converge with |dz| <= 1e-6 (acceptance). There is NO gate between the direct solution and the refinement (the open-loop flight of the N = 400 direct control misses by 1e3 to 4e4 km on these multi-rev arcs, yet the refinement converges from it; with a 1000 km pre-gate the pilot solved 57/50/33/26/10 of 72 cells per rung, without it 56/58/55/57/54).
- Post-fleet chain: (a) Earth sweep: fly every accepted z8, take the minimum Earth altitude over the integrator's output samples, demote entries below 200 km; (b) repair 'unsolved' cells by phase continuation: seed the top-rung direct solve with the flown solution of the nearest solved departure phase (resampled to the nodes); (c) repair 'basin': run the conjugate-point (Jacobi) test on every top-rung entry; each refuted one is re-solved at 15 N seeded from its own 12 N entry (reverse warm start); keep the new solution only if it is accepted, strictly faster and loses no lower rung, else restore the backup; (d) conjugate test and hypothesis gates (min |lambda_v| > 1e-6 i.e. no singular arc, abnormal-lift dimension 1) on every entry; refuted entries ship flagged conj_pass = 0, they are not removed.

## Stage 1 finding that motivated the basin repair

One cell, 15 N: the floored cold start lands at t_f 0.3040 ND and is REFUTED by the conjugate test (conjugate point near t/t_f 0.71 where |lambda_v| dips to 6e-4). The 12 N warm start from it gives 0.3070 ND and PASSES. 15 N warm-started from that 12 N solution gives 0.2668 ND (faster, passes). So the first ladder step leaves the refuted branch.

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
  0.2589 to 0.3238 (n = 56). Within 0.02 of 0.30, so the fleet keeps
  `tf0 = 0.30`.
- **Earth sweep (chain stage 2).** 46 of the 280 entries fly below the 200 km
  Earth floor and are demoted; 41 of them pass below the surface (deepest
  -3,855 km). Per rung: 5 / 5 / 12 / 17 / 7. The plan expected none. 43 of
  the 46 lie in a band of departure phases that moves with thrust (phases
  11-12 of 12 at 15 and 12 N, 10-11 at 10 N, 8-10 at 7 N, 3-5 at 5 N), and the
  dip is at the first perigee pass (t / t_f 0.03 to 0.12). Re-flown and
  sampled at the 401 uniform node times, 45 of the 46 are below the floor at
  a node, and 43 stay below it longer than one node spacing (up to 3.8): the
  refined trajectory is not the direct solution, which holds the floor at
  its nodes. The refinement (`ms_tfmin`) carries no Earth constraint. NOT
  measured: whether the floor is active in the direct solution of those
  cells (the sheet does not keep it); if it is, an unconstrained costate
  entry cannot represent the optimum there and demotion is the final
  verdict. The same sweep on the 1000 km-gate sheet finds 6 of 176 below
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
- **Not answered here.** The 180-degree orientation; the fleet runs those
  sheets first (`"[9:12 1 3:8 13:16]"`).

## Questions

1. THE EARTH-FLOOR DEMOTIONS (lead question). 46 of 280 refined entries fly through or below the Earth floor at the first perigee pass, although the direct solutions they were refined from satisfy the floor at the nodes. (a) Is my reading right: where the state constraint is active in the constrained optimum, the unconstrained PMP BVP has no admissible solution near it, the refinement converges to the unconstrained extremal of the same homotopy class, and the entry must be demoted? What measurement would confirm it cheaply (e.g. the floor's multiplier in the direct solve, which I can read from IPOPT)? (b) The band moves with thrust (departure phases 11-12 at 15 N down to 3-5 at 5 N). Does that pattern fit the explanation, or does it suggest something else (e.g. the node spacing of 292 s at 15 N versus perigee speed ~10 km/s, i.e. a between-node violation in the DIRECT solution that the refinement merely inherits)? (c) What should the fleet do about these cells: accept the loss (about 16 percent of a sheet), look for another admissible basin (extra revolution, longer t_f) by seeding, densify/replace the constraint in the direct solve (perigee-radius constraint, mesh refinement near perigee), or add the state constraint to the indirect refinement (boundary arc / penalty)? Which is worth doing BEFORE the fleet and which after?
2. No gate between direct and refined. Refinement closes direct open-loop misses of up to 43,000 km. What can go wrong when the refined solution is that far from the direct one, beyond the Earth-floor case above (different homotopy class, different local minimum than the NLP's, entries in neighbouring cells on different branches)? Is the set of post-gates (flown miss, tfMin acceptance, Earth sweep, conjugate test, hypothesis gates) sufficient, and what is missing?
3. Basin repair logic. First-success-wins plus reverse warm start: can a repaired 15 N entry be a worse choice that happens to be faster? What exactly does "passes the conjugate test after repair" prove and not prove? 5 of 18 re-solves returned the identical refuted solution: what would you try next for those? The share of refuted entries GROWS down the ladder (10/59 at 15 N to 20/54 at 5 N) and lower rungs are not repaired: is shipping them flagged acceptable, and is there a cheap repair (re-solve from the rung below, as was done at the top)?
4. Phase continuation for the perigee row: a neighbour at 1/12 of the GTO period seeds phase 0, where the all-burn arc departs from perigee. All 8 unsolved cells repaired. Anything suspicious in that being so easy when the cold start fails outright?
5. Earth sweep as an instrument: minimum over the ODE solver's output samples of a single-shooting flight. Is that enough near perigee (speed ~10 km/s, the dip lasts minutes)? What would you require (event function, dense output, refine factor)?
6. Deep rungs: the wall between 3 N and 2 N shows as flown misses of 300,000+ km after refinement or as direct failures. Is that consistent with a sensitivity wall of single-shooting flight (the gate, not the solution, failing) and how would you tell the two apart?
7. Anything in the fleet recipe that the pilot (one orientation, one DRO size) did not measure and that could invalidate 40 hours of compute.
