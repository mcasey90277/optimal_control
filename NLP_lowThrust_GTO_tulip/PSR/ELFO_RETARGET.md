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

# Where we are going (open)

**Goal: min-fuel GTO->ELFO** (as we have for GTO->tulip). That needs a
**GTO->ELFO min-ENERGY seed** (the homotopy root the PSR fuel pipeline consumes),
which is exactly the missing piece. The one path that would have BOTH
direct-collocation robustness AND floating t_f -- and moot the edge issue via
always-burn -- is a **direct min-time collocation** (modify `casadi_minfuel_
sundman` to minimize t_f, free final time). That is the leading candidate for
manufacturing the GTO->ELFO energy seed; not yet built.
