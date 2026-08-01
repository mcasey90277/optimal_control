# pumpkyn / pumpkynPie — what we can leverage

Hand-written companion to **`pumpkyn_catalog.md`**, which is generated and holds
every signature. This file holds the judgement: what is relevant to our four
campaigns, what we have been duplicating, what conventions bite, and what is
worth trying next.

Verified against pumpkyn `5f5ca31` / pumpkynPie `47be599` (2026-07-30).
Signatures move — regenerate the catalog and re-check before relying on one.

---

## The two things this survey found that change what we should do

### 1. `minDeltaV` is an independent lower bound we are not using

```matlab
dV = pumpkyn.cr3bp.minDeltaV(x0, x1, muStar)   % x0, x1: [Nx6] ND states
```

> "the theoretical minimum cost associated with changing the energy of the
> spacecraft from an initial orbit to a final orbit **without considering its
> transfer path**"

This is exactly the kind of check the basin work showed we lack. Every ΔV we
report — 3.35 km/s for GTO→tulip, 0.72 km/s for DRO→tulip — is the cost of *a
particular local minimum*, and we have no way to say how far off the floor it
sits. `minDeltaV` is path-independent, so it does not care which basin we
landed in. **Any transfer ΔV should be quotable against it.**

*Not yet run.* If a reported ΔV came out below this bound, something is wrong;
if it sits far above, that is a quantified argument that better solutions exist
— which the basin sweeps have been finding by brute force.

### 2. `toModEq` exists, and the earth campaigns hand-rolled MEE

```matlab
[meOEV, theta, psi] = pumpkyn.cr3bp.toModEq(tau, rvND, aT, M, muStar, tStar, lStar, P)
```

CR3BP rotating-barycentric → modified equinoctial elements about either
primary, citing the JPL modified-equinoctial reference. **It also converts the
thrust acceleration vector `aT` into the RTN frame and returns pitch and yaw
angles** — i.e. it handles the control, not just the state.

`earth_elliptic_to_geo` built its own MEE formulation (`lt_mee_rhs`, the
L-domain Gauss equations). That work is done and certified, so this is not a
call to rewrite it, but any *new* MEE work should start here, and it is worth a
cross-check of the two conventions before the campaigns are compared.

---

## By campaign

| campaign | what pumpkyn offers |
|---|---|
| GTO→tulip | `getTulip` (already used via `gto_tulip_endpoints`), `minDeltaV` bound, `manifolds` for the unsolved near-min-time band, `tulipConstellation` |
| GTO→ELFO | `minDeltaV`, `orbitProperties` for the target |
| elliptic→GEO (2-body) | `toModEq` cross-check, `toOrb`/`fromOrb` |
| elliptic→GEO (CR3BP) | `threeBody` (Cartesian 3-body with perturbation outputs `aP1`, `aP2`) as an independent check on the lunar-perturbation branch |
| all | `tfMin`/`tfMinEoM`/`tfMinProp` — a working indirect min-time solver |

---

## The indirect min-time solver, which we now use

```matlab
sol = pumpkyn.cr3bp.tfMin(rv0, rvf, lambda0_tf, Tmax, c, muStar)   % 8 unknowns
[tau, y, Hf, dHdy] = pumpkyn.cr3bp.tfMinProp(tf, y0, Tmax, c, muStar)
[yDot, Ht, dHdy, aThrust] = pumpkyn.cr3bp.tfMinEoM(tau, y, Tmax, c, muStar)
```

Decision vector `[lambda_r(3); lambda_v(3); lambda_m; tf]`. Residuals:
`r(tf)=r_f`, `v(tf)=v_f`, `lambda_m(tf)=0`, `H(tf)=0`. Analytic Jacobian from
the propagated STM.

**Two things worth knowing:**

- **`tfMin` returns no success flag.** Verify independently by re-propagating
  and checking those four conditions — see
  `abstracts/data/bht1500_continuation.m`. Note that an `ode45` re-propagation
  at default tolerances is *stricter* than the solver's internal residual: the
  published 70 mN case re-propagates to 5.0e-06 where `fsolve` reports ~1e-11.
- **`tfMinProp` uses event detection to stop and restart integration at
  switches.** That is switch-aligned integration done properly, and it is the
  thing our direct campaigns do *not* do — the mesh study spent two phases
  discovering that an unaligned mesh is the problem. Worth studying before
  building the multi-phase formulation.

**Its own header says thrust continuation is the way to get convergence**, which
matches what we measured: see `bht1500_continuation.m` for a working example
(101 mN / 1710 s reached in four steps, zero failures).

---

## Frames and conventions

Most routines are **nondimensional CR3BP rotating-barycentric**, scaled by
`muStar`, `tStar`, `lStar` from `pumpkynPie.cr3bp.getCR3BParams()`.

- `toOrb` / `fromOrb` — classical elements about a chosen primary
- `toModEq` — modified equinoctial (+ control angles)
- `toPCI` / `fromPCI` — planet-centred inertial
- `toJ2K` / `fromJ2K` — J2000 inertial
- `toLLA` / `fromLLA` — lat/lon/alt
- `eci2orb` / `orb2eci` — two-body element conversions

**Traps:**
- A trailing **`P`** argument selects the primary body. Get it wrong and the
  elements are about the other body, silently.
- A trailing **`dim3`** flag appears on several routines (`eom`, `threeBody`,
  `elAng`, `bsxAng`, …) and controls which array dimension is the vector axis.
- **Dimensional vs nondimensional is not uniform.** `minDeltaV` takes ND
  states; `toOrb`/`toModEq` take ND states *and* the scale factors and return
  dimensional elements. Read the header.
- `getTulip` returns **catalogued** initial states, then `cont_np` refines to a
  closed periodic orbit. Both steps are needed — this is how
  `cr3bp_common/gto_tulip_endpoints` does it.

---

## Orbit characterization we are not using

`jacobi` (Jacobi constant), `monodromy`, `stabilityIndex`, `orbitProperties`
(aggregated diagnostics), `lagrangePts`, `manifolds` (stable/unstable manifolds
of an unstable periodic orbit).

**`manifolds` is the interesting one.** The 1.01–1.11× near-minimum-time band
is unsolved and fails for the *smooth* energy problem too, so it is
conditioning rather than bang-bang structure. Manifold-guided transfers are the
standard tool there, and the Zhang et al. Note we cite describes prior work
patching a tangential-thrust spiral onto a stable-manifold leg. Untried.

## Constellation and coverage

`tulipConstellation` (equally spaced satellites on a tulip),
`occultationCalc` (occultation by either primary — the arXiv paper's R5
Earth-occultation metric), `dop`, `maxGaps`, `elAng`, `pointSphere`. These are
the proj7 coverage tools; they live in pumpkyn, not proj7.

## Other transfer machinery

`directLambert`, `minLambert` (single-revolution Lambert in CR3BP — compare
against our own `orbit_transfer/lambert/`), `stationKeeping_deltaV` (bounds,
relevant to the paper's ΔV budget), `cont`/`cont_np` (pseudo-arclength and
natural-parameter continuation — for *periodic orbit families*, not transfers,
but the predictor-corrector pattern is what the mesh-study reviewers said our
branch tracking needs).

---

## Demos as worked examples

`demos/` in pumpkynPie:

- **`lowThrustDRO2Tulip.m`** — indirect min-time DRO→tulip. **Its three
  commented-out costate blocks are converged solutions at its three commented
  thrust tiers** (final times 2.709 / 4.015 / 6.189 against 0.1 / 0.07 /
  0.055 N). Not stated anywhere; found by matching. Seeding from the right tier
  is what made the BHT-1500 continuation work.
- **`lowThrust_constellation_animation.m`** — staggered deployment into a phased
  tulip constellation.
- **`constellation_Animation.m`** — the scene our movie recorders mirror.

## Operational rules

1. **Our code lives in `optimal_control/`.** Call into pumpkyn; never write into
   it. Pattern: `cd(pieRoot); startup();` as in `bht1500_continuation.m`.
2. **Pull pumpkyn and pumpkynPie together.** The newer `SatelliteAnimator`
   passes `'AddLighting'` to pumpkyn's `showMoon`; updating one alone breaks the
   other — this cost a debugging cycle on 2026-07-30.
3. **When a campaign solver builds dynamics inline** (the tulip's
   `casadi_minfuel_sundman` does), any rebuild must be checked against the
   original by a self-validating gate. A missing `tauf0` factor was caught
   exactly this way and would otherwise have read as a physics finding.
4. **Re-diff the movie recorders after a pumpkynPie pull.** They mirror the
   demo rather than calling it; `SatelliteAnimator` changed by 1145 lines on
   2026-07-30 and the recorders have not been re-checked since.
