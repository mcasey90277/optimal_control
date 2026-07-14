# PSR → ELFO retarget — build record

Adapting the working direct (PSR) pipeline to fly GTO → **elliptical lunar
frozen orbit (ELFO)** instead of GTO → tulip. ELFO shape from proj7
(`im_elfo_optimum.m`): sma 12000 km, ecc 0.69, inc 56.5°, argp 90° — the
south-pole nav orbit from the lunar-nav paper. Started 2026-07-13.

## Key architectural finding — the PSR pipeline is target-agnostic

`minfuel_at_tf` threads `rv0`/`rvf` **from the seed file** (`E.rvf`) straight
into `casadi_minfuel_sundman`; it never calls `gto_tulip_endpoints`. So the
ONLY tulip-specific thing in the entire direct chain (solve → refine → export →
verify → movie) is the `rvf` baked into the energy-backbone `.mat` files.
**Retargeting = produce an energy backbone whose `rvf` is the ELFO. Nothing
else changes.**

## ELFO rendezvous geometry (probe_elfo_target.m)

Rendezvous at **apolune** (highest, slowest, over the south pole — the natural,
gentle low-thrust capture point):

| state | dist Earth | dist Moon | speed (ND) |
|---|---|---|---|
| GTO rv0            |   6,728 km |   383,616 km | 9.92 |
| ELFO **apolune**   | 384,376 km |    20,388 km | 0.23 |
| ELFO perilune      | (1,921 km alt) | 3,659 km | (fast) |

The tulip endgame was **already near the Moon**, so ELFO is not a new region;
apolune is farther out and 5× slower than the tulip target → a **gentler**
terminal, not a stiffer one. The two-primary Sundman `kappa` (r1·r2) I expected
to need for lunar capture is likely unnecessary for an apolune rendezvous.

## Endpoint-helper discrepancy (flagged, not a blocker)

The energy backbones target a DIFFERENT tulip point than `gto_tulip_endpoints`
returns:

| tulip target | dist Earth | dist Moon | speed (ND) |
|---|---|---|---|
| backbone `E.rvf` (what the solves actually hit) | 397,944 km | 28,303 km | 0.31 |
| fresh `gto_tulip_endpoints` (max-ẏ) | 393,023 km | 6,040 km | 1.15 |

`E.rv0` matches the fresh GTO exactly (Δ=0) — only the *target* differs, by 1.13
ND. Because the pipeline threads `rvf` from the seed, **every certified PSR
tulip solution rendezvouses at `E.rvf`, not the max-ẏ point.**
`gto_tulip_endpoints` is stale relative to the backbones (different `getTulip`
args when they were built). Consequence for us, both favorable: the backbone's
real target (slow, dMoon 28k km) is a **near-neighbor of the ELFO apolune**
(slow, dMoon 20k km), so the target homotopy is only **0.48 ND** — the ELFO
retarget is a small perturbation of what the campaign already solved.

## Two seed routes, fired in parallel

### (A) Homotopy from tulip — PRIMARY  [`gen_elfo_energy_backbone.m`]
Reuse a converged tulip energy backbone (a real ~40-rev GTO→lunar-vicinity
spiral) and slide `rvf` from the tulip target to the ELFO apolune in small
steps: loose continuation + tight re-clean per step (the `gen_energy_seed`
recipe, but stepping the TARGET instead of t_f). Adaptive step, **per-step
checkpoint/resume** (`energy_elfo_ckpt.mat`) added after seeing the campaign's
~1-in-10 CasADi MEX crash risk over ~18 solves.
- Mechanism validated; but at **f=1.50 the walk stalled at s≈0.2**: the loose
  continuation could not satisfy the moved terminal constraint (returned the
  warm start; smaller steps then diverged, inf_du→1e12).

**Root cause = throttle saturation (edge), and the fix is COUNTERINTUITIVE.**
The min-energy solution's saturated-node fraction (edge) is NOT monotone in t_f:

| factor | edge | tf |
|---|---|---|
| **1.20** | **7.8%** | 33 d |
| 1.50 | 54.7% | 42 d |
| 1.70 | 70.0% | 47 d |
| 1.85 | 51.9% | 51 d |

More time does NOT mean smoother control — higher t_f buys long coast + full-burn
arcs (both count as edge), so the SMOOTHEST energy solution is at the LOW end.
At f=1.50's 55% edge, a small terminal move must restructure the bang arcs and
the continuation fails; at **f=1.20's 8% edge the intermediate throttle just
deforms**, and the walk sails: step 1 defect **3.8e-14** (edge 12%), step 2
(the f=1.50 killer) defect **3.5e-14** (edge 12%). Walk in progress at f=1.20.

**Lesson:** run a target/parameter homotopy on the LEAST-saturated member of the
family (lowest edge), not the one with the most time. Check edge before picking
the continuation root.

### (B) Tangential from scratch — INDEPENDENT CHECK  [`gen_elfo_energy_tangential.m`]
Target-agnostic bootstrap: propagate max-thrust velocity-aligned steering from
the GTO, map to the Sundman mesh targeting the ELFO apolune, solve energy with
NO tulip seed. **VERDICT: does not converge.** The pure tangential spiral over
41 days ends **3.8 M km from the Moon** (mass frac 0.71) — an enormous cold
rendezvous gap. The energy NLP stalled: inf_pr frozen ~0.57, inf_du diverged to
1.2e9 over 440+ iterations. Confirms `build_guess.m`'s own warning and the
campaign's "cold multi-rev seed is explosive" lesson. → **The homotopy route is
the right one; a from-scratch spiral cannot find the ELFO.**

## Lesson (transferable)

Retargeting a converged multi-rev low-thrust transfer is **cheap via target
homotopy, hopeless via cold re-bootstrap.** The expensive, fragile part of
these problems is obtaining ANY converged N-rev spiral into the lunar vicinity;
once you have one, sliding the terminal state to a nearby target is a short,
robust continuation. Do not rebuild the spiral from a heuristic guess.

---

# Where the energy-target homotopy actually stalled (2026-07-13)

The direct energy-target homotopy is the FURTHEST-progress route but does not
finish: it is clean to **s=0.45** (machine precision, edge<1%) then the back
half fails. Two causes, both diagnosed:
1. **Moon-ward terminal stiffness** — the nearest-insertion ELFO target sits
   deeper in the Moon's gravity well; at fixed t_f the collocation "can't reach
   terminal" (both loose and tight return the warm start, defect ~3e-10, ok=0).
2. **A CasADi MEX bus-error crash** mid-walk (the ~1-in-10 fatal). Checkpointing
   (`energy_elfo_ckpt.mat`, s=0.45) made this survivable.

Also characterized: the **apolune** target is un-homotopable at fixed t_f — the
tulip-terminal-to-apolune velocity angle is **121 deg**, so a linear velocity
interpolation COLLAPSES the terminal speed (0.31 -> 0.13 mid-path), a degenerate
rendezvous; the bifurcation sits at s~0.3 regardless of factor. Nearest-
insertion (velAngle 63 deg) avoids that but hits cause (1) above.

# Min-time as a root (min_time/ module, 2026-07-13)

Since min-time was the tulip's own root (min-time -> energy -> fuel) and is
always-burn (tf floats, no edge/saturation), we tried it as the ELFO root:

| route | result |
|---|---|
| indirect min-time **single shooting** (`mintime_solve`) | floors ~1e-3 (13-rev STM-product sensitivity); ELFO homotopy stalls at s=0.05 |
| indirect min-time **multiple shooting** (`mintime_ms_*`) | **tulip: 4e-9 (validated, beats the wall).** ELFO target homotopy fights min-time's shooting sensitivity even with predictor-corrector -- impractically slow |
| direct min-time (attic `solve_tfmin_nlp`, fmincon) | does not scale: t_f plunges / infeasible at usable N |

Net: the indirect min-time MS is a real, validated result for the tulip, but
**retargeting to the ELFO resists both fixed-t_f direct continuation (Moon-ward
stiffness) and shooting continuation (sensitivity).** See `min_time/README.md`.

# External design review (GPT-5.6-terra + Gemini 3.1 Pro, 2026-07-13)

Fired a two-model design review (not a code audit) on how to manufacture the
GTO->ELFO energy seed. Strong, independent convergence:

- **KILL the direct-min-time-collocation plan.** Both models reject it as a
  detour: min-time is s==1 bang-everywhere, structurally far from a smooth energy
  ramp, so min-time->energy re-introduces a hard restructuring. Both also flagged
  the SAME trap in a naive free-t_f min-time: minimizing t(tau_f)=Int[kappa]dtau
  rewards the optimizer for shrinking r (diving at a primary to slow the clock)
  -- which retroactively explains route (5)'s "t_f plunges."
- **Fix Route (1) instead**, via three changes (all now built into
  `sundman_minfuel/casadi_energy_freetf.m`):
  1. **Two-primary Sundman clock** kappa=(r1^-q + (r2/D)^-q)^(-p/q), p=1.5, q~4,
     D=moonZone~0.15 (~lunar SOI). Recovers r1^p near Earth, (r2/D)^p near Moon;
     redistributes mesh into the lunar-capture arc where the single-primary clock
     starved. (terra's soft-min form + gemini's D-scaling; D<=0 -> original.)
  2. **Free physical t_f with a BANDED KKT** via a constant slack STATE cScale
     (Betts): dt/dtau=cScale*kappa, dcScale/dtau=0. One number tied by LOCAL
     continuity constraints -> Jacobian stays banded (a free SCALAR tau_f would
     make one dense column -> OOM). Every intermediate target is now reachable at
     some t_f, killing the fixed-t_f "can't-reach-terminal" wall.
  3. **Moon-gravity homotopy**: hold rvf FIXED, continue muGain: 0->1 scaling
     ONLY the Moon-gravity term -muGain*muStar*rr/r3 (NOT muStar in the frame/
     Coriolis -- that would move the barycenter under the BCs). muGain=0 is a
     well-less near-2-body transfer; the linear Cartesian retarget that was toxic
     at muGain=1 is benign with the well off.
- Where they diverged: terra preferred a free ELFO-insertion-phase terminal
  manifold; gemini preferred the gravity homotopy. Taken as complementary
  (gravity homotopy first -- cheaper, fixed true target; phase-freedom is the
  fallback structural upgrade). Both ranked meet-in-the-middle two-phase
  collocation as the #2 fallback (as DIRECT matching, NOT backward shooting).

Solver `casadi_energy_freetf.m` smoke-validated 2026-07-13 (N=4000 f1.20
backbone): with the clock matching the backbone (moonZone=0) the free-t_f slack
formulation reproduces it to **9.5e-8** in 40 iters (t_f floats 7.55->8.30,
cScale~1.0, no runaway) -- proving the KKT stays solvable and t_f is tame. The
two-primary + gravity-off cases construct and step (defect high pre-continuation,
as expected; that is the driver's job).

# GTO->ELFO min-ENERGY seed: DONE (2026-07-13)

**SOLVED.** `sundman_minfuel/gen_elfo_energy_gravhom.m` produced a
machine-precision GTO->ELFO min-energy solution: `results/energy_elfo_freetf.mat`
(9-row free-t_f, tf=7.5488 ND=33.46 d, mf=0.8430 [15.7% prop], terminal dMoon
16799 km at full CR3BP gravity). Independently verified (solver-free MATLAB defect
recompute): defect **1.77e-15**, unit-norm 1.2e-12, endpoints exact (GTO 8e-33,
ELFO 0.00). `verify_elfo_seed.m`.

**Two fixes were needed beyond the three review changes:**
1. **Pin t_f** (opts.tfTarget) -- free-t_f min-ENERGY is ill-posed (energy optimum
   drifts t_f -> longer time, thinner thrust, lower Int[s^2]dt; IPOPT wandered off
   the warm start, loose->3.7e-3, tight->0.21). Constraining t(tau_f)=tfTarget and
   letting cScale float to satisfy it makes it well-posed (machine precision).
   [The slack state still earns its keep: a clean single-DOF way to hold a fixed
   t_f under the *changing* two-primary clock.]
2. **Leg ORDER** (found empirically, each step observation-driven):
   LEG 0  free-t_f convert (mu=1, single-primary, tulip)   -- machine precision
   LEG A  gravity OFF muGain 1->0 (single-primary, tulip)  -- cleanest leg
   LEG B  clock ON moonZone 0->0.15 (mu=0, tulip)          -- BENIGN with well off
   LEG C  retarget tulip->ELFO (mu=0, clock on)            -- THE crux, walked clean
   LEG D  gravity ON muGain 0->1 (ELFO, two-primary)       -- dissolves s=0.45 wall
   Two dead orders were tried and rejected: (i) retarget BEFORE clock-on ->
   mesh-starved near-Moon terminal -> dual stalls (iter 973, tiny steps); (ii)
   clock-on at mu=1 -> concentrating nodes in the full-gravity well stiffens at
   moonZone~0.09. Lesson: do the clock-on with gravity OFF (benign re-mesh), and
   the retarget with the clock already on (terminal resolved). Every leg then
   converges to ~1e-14; edge stayed ~30-38% throughout.

**Min-fuel GTO->ELFO at tf=1.20x: DONE (2026-07-13).** `gen_elfo_minfuel.m`
ramped epsilon 1->0 from the energy seed (fixed tf=7.5488 ND, two-primary clock,
full gravity): reached eps=0 (pure FUEL) at machine precision. Result
`minfuel_elfo.mat`: **34-switch bang-bang, edge 99.6%, mf=0.8545 (14.5% prop)**,
defect 5.7e-14, independently verified. The epsilon ladder: eps=0.8 -> 8 sw,
0.6 -> 8, 0.4 -> 18, 0.2 -> 40, 0.0 -> 34 sw (count settles as near-switches
resolve). The sharpening wall that left the tulip many-switch case open did NOT
stop this -- the two-primary clock keeps the lunar-arc mesh resolved through
bang-bang.

**Next: the min-fuel tf-GRID (per [[minfuel-tf-grid-strategy]]).** Min-fuel is a
tf-grid convergence map, not a single-tf solve (energy band is wider than the
eps=0-convergent band; some tf's stall before fuel). tf=1.20x is grid-point #1
(reaches eps=0). Plan: tf-continuation on `casadi_energy_freetf` (step opts.tfTarget)
to spread ELFO energy seeds across a tf band -- cheaper than re-running the full
gravity ladder per tf -- then `gen_elfo_minfuel` at each grid tf, recording
eps-reached / switches / mf. Output = the ELFO min-fuel convergence map.
