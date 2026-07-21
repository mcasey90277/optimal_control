# NLP_lowThrust_GTO_tulip

Direct-NLP (collocation) solver for the minimum-time low-thrust GTO -> tulip
transfer in the Earth-Moon CR3BP — the direct-method twin of pumpkynPie's
`Demos/LunaNet Analysis/lowThrust_GTO_Tulip.m` (which solves the same problem
indirectly via PMP costate shooting). Produces the same outputs: transfer
time, propellant / delta-V accounting, and the two-panel trajectory figure.

Problem instance: 15 kg spacecraft, 25 mN, Isp 2100 s, GTO (350 x 35,786 km,
argp -25 deg) to the 7-petal southern tulip (ND period 5/6 * 2*pi), rendezvous
at the max-ydot tulip point. Reference answer (indirect):
tf = 6.290694 ND = 27.8845 days, propellant 2.9247 kg, dV = 4.4665 km/s.

## Files

| file | purpose |
|------|---------|
| `NLP_lowThrust_GTO_Tulip.m` | driver: endpoints via pumpkyn, guess, solve, report, figure |
| `solve_tfmin_nlp.m` | fmincon interior-point wrapper (sparse analytic gradients, lbfgs, warm-start barrier) |
| `nlp_constraints.m` | trapezoidal defects + unit-sphere control equalities, sparse Jacobian |
| `lt_dynamics.m` | 7-state max-thrust dynamics f(x,w) + analytic A, B |
| `build_guess.m` | density-matched nonuniform mesh + warm start ('indirect' or 'tangential') |
| `unpack_z.m` | decision-vector bookkeeping Z = [X(:); W(:); tf] |
| `setup_paths.m` | adds pumpkyn (proj7/external/pumpkyn) to the path |

## Usage

```matlab
cd NLP_lowThrust_GTO_tulip
out = NLP_lowThrust_GTO_Tulip();            % N = 3000, indirect warm start, plot
out = NLP_lowThrust_GTO_Tulip(6000);        % finer mesh
```

## Formulation notes (hard-won)

- **Throttle fixed at 1** (control = unit direction vector, `w'*w = 1`
  equalities). Justified post-hoc: the indirect switching function stays in
  [-46.8, -2.5] (never coasts). Carrying a throttle variable puts the warm
  start exactly on its bound and interior-point fmincon shoves it off the
  wall, coasts, and stalls ~2% suboptimal.
- **Density-matched nonuniform mesh** (node fractions = quantiles of the
  reference integrator's adaptive grid). A uniform N=3000 mesh leaves
  warm-start defects ~0.9 (perigee passes unresolved); density-matched is
  6.4e-4 at the same N.
- **`InitBarrierParam` 1e-6**: the default barrier (0.1) drags a warm-started
  iterate far off the solution (observed tf -> 4.0 at feasibility 1.1).
- **No open-loop replay validation**: ~40 perigee passes amplify control
  interpolation error by ~1e6; replay diverges for any control, including
  the true optimum. Validate by defect norm + node-wise deviation from the
  indirect arc + mesh-refinement convergence of tf.
- **Mesh-refinement convergence** (indirect warm start, fmincon flag 2,
  ~0.5 min each):

  | N | tf (ND) | error vs indirect | max node dev | dV (km/s) |
  |---|---------|-------------------|--------------|-----------|
  | 3000  | ~6.2658  | -2.5e-2 | 11,000 km | 4.4468 |
  | 6000  | 6.281802 | -8.9e-3 | 3,970 km  | 4.4595 |
  | 12000 | 6.288574 | -2.1e-3 | 1,010 km  | 4.4648 |

  (N=3000 stalls in a flat valley: its tf scatters ~1e-3 across
  numerically equivalent code variants; N >= 6000 reproduces to 7 digits.)

  Driver default is N = 12000 (~1 min): tf = 27.8751 days and
  dV = 4.4648 km/s vs the indirect 27.8845 days / 4.4665 km/s.

  The discrete optimum exploits trapezoidal discretization error near
  perigee (defects at machine zero, tf slightly below truth) and walks back
  to the continuous answer as N grows.

## Min-fuel variant (arrival-leg replan)

`NLP_lowThrust_GTO_Tulip_minfuel.m` + `lt_dynamics_throttle.m` +
`nlp_constraints_minfuel.m` + `solve_minfuel_nlp.m` + `costate_seed_from_nlp.m`.
Fixed leg time (1.3x the leg minimum from tau = 4.0 on the min-time arc,
arrival at the tulip point one coast downstream), objective = final mass,
throttle in [0,1] via the general 4-control cone transcription.

- Direct solve CONVERGES (defects ~2e-15): propellant 1.0622 kg vs the
  burn+coast reference 1.0650 kg, clean single-switch bang-bang at
  t = tfMinLeg. Warm start must be the burn+coast construction (feasible
  to 1e-4); time-stretched or original-rvf-pinned guesses fail (measured).
- Indirect polish is OPEN: seed progression (raw / rescaled / anchored /
  LSQ-reconstructed costates) stalls at ||R|| = 1.55 / 0.83 / 0.33 / 0.14.
  The full 25-mN min-fuel spiral (~80 switches) is research-grade.
  Details: theory note S6, tutorial Phases G-H.

## Companion documents

- Theory: `../lowThrust_GTO_tulip/gto_tulip_mintime_theory.pdf`
- Guided build-it-yourself tutorial (this solver is its Phase D-F answer
  key): `../lowThrust_GTO_tulip/building_the_gto_tulip_solvers.pdf`
- Indirect counterpart: `../lowThrust_GTO_tulip/` (complex-step shooting)
