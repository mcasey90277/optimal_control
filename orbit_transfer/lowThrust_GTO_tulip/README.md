# lowThrust_GTO_tulip

Guided-tutorial home for the minimum-time low-thrust GTO -> tulip transfer
in the Earth-Moon CR3BP (the problem behind pumpkynPie's
`Demos/LunaNet Analysis/lowThrust_GTO_Tulip.m`). Contains the theory note,
the build-it-yourself coding tutorial, and the reference INDIRECT (PMP
shooting) solver. The reference DIRECT (collocation NLP) solver lives in
`../GTO_tulip/`.

## Documents

- `gto_tulip_mintime_theory.pdf` — problem + vocabulary (GTO, CR3BP, tulip,
  primer, switching function), OCP formulation, indirect method (all 14
  ODEs, 8 boundary conditions, target manifold + transversality), direct
  NLP transcription, and the results table (both methods vs pumpkyn).
- `building_the_gto_tulip_solvers.pdf` — guided exercises, Phases A-F:
  augmented PMP dynamics -> complex-step shooting residual -> indirect
  solve -> control-explicit dynamics -> transcription + density-matched
  mesh -> fmincon solve. All checkpoint numbers verified by running the
  reference implementations headless.

## Reference solution (indirect)

| file | purpose |
|------|---------|
| `lt_pmp_eom.m` | 14-state augmented PMP dynamics (complex-safe) |
| `shoot_residual_tf.m` | terminal residual + Jacobian (complex step, analytic tf column) |
| `solve_tfmin_indirect.m` | fsolve trust-region-dogleg driver |
| `run_gto_tulip_indirect.m` | end-to-end driver: endpoints via pumpkyn, solve, dV, figure |
| `setup_paths.m` | adds pumpkyn (proj7/external/pumpkyn) to the path |

Verified result: tf = 6.290694 ND = 27.8845 days (agrees with
`pumpkyn.cr3bp.tfMin` to 8 significant figures), propellant 2.9247 kg of
15 kg, dV = 4.4665 km/s, ~25 s wall.

## Min-fuel extension (indirect side)

`lt_pmp_eom_minfuel.m` (Bertrand-Epenoy smoothed throttle),
`shoot_residual_minfuel.m` (fixed-tf, 7 conditions, complex-step Jacobian),
`solve_minfuel_indirect.m` (Levenberg-Marquardt + eps-continuation),
`run_gto_tulip_minfuel.m` (indirect driver; requires a supplied costate
seed). The min-fuel TPBVP itself remains unconverged from the best seed
available (see `../GTO_tulip/README.md` and the theory note
S6 for the honest accounting); the direct solver in the sibling folder is
the converged min-fuel reference.

Learner workspace: `mytry/` (build your own versions there; consult the
reference only after attempting each exercise).
