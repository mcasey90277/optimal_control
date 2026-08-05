# The DRO → Tulip Costate Libraries

*Reference for `costate_lib_dro_tulip.mat` (v1) and
`costate_lib_dro_tulip_v2.mat` (v2). Written 2026-08-05.*

Both libraries answer the same question — *what costates fly a minimum-time
low-thrust transfer from a DRO to a tulip orbit?* — for a grid of cases. Each
entry is a **converged solution of the indirect (PMP) problem**, not a hint:
handed to `pumpkyn.cr3bp.tfMin`, an entry comes back unchanged.

| | **v1** | **v2** |
|---|---|---|
| File | `costate_lib_dro_tulip.mat` | `costate_lib_dro_tulip_v2.mat` |
| Structure name | `costate_lib_dro_tulip` | `costate_lib_dro_tulip_v2` |
| Axes | departure phase × arrival phase | departure phase × arrival phase × **thrust** |
| Entries | 126 | 990 |
| Coverage | 126 of 144 phase pairs | 130 of 144 phase pairs; 49 with all 9 rungs |
| Thrust | 0.07 N (single) | 1, 1.5, 2, 3, 5, 7, 10, 12, 15 N |
| Isp / initial mass | 900 s / 150 kg | 1710 s / 150 kg |
| Flight times | 16.21 – 33.31 days | 0.269 – 2.560 days |
| Size on disk | 22 KB | 145 KB |
| Helper | `costate_lib_nearest.m` | `costate_lib_pick.m` |
| Example | `costate_lib_example.m` | `costate_lib_v2_example.m` |

Files live in `orbit_transfer/DRO_tulip/`: libraries under `direct/results/`,
code under `indirect/`.

---

## 1. The problem each entry solves

Minimum-time low-thrust transfer in the Earth–Moon circular restricted
three-body problem (CR3BP), rotating frame, nondimensional units. Pontryagin's
Minimum Principle reduces it to finding eight unknowns

```
z8 = [ λ_r(3) ; λ_v(3) ; λ_m ; t_f ]
```

such that the trajectory propagated from the departure state under the PMP
control law satisfies

```
r(t_f) = r_f ,   v(t_f) = v_f ,   λ_m(t_f) = 0 ,   H(t_f) = 0
```

the last being the free-final-time transversality condition. This is exactly
the system `pumpkyn.cr3bp.tfMin` solves; the libraries store its roots.

**Why a library exists:** single shooting on this system converges in seconds
*from a good guess* and not at all otherwise. The guess is the entire
difficulty. A library of converged costates turns every future solve into a
lookup plus a seconds-long polish, and provides continuation anchors for
neighboring cases.

## 2. The two orbits

Identical in both libraries, and recorded inside each file so it is
self-describing.

**Departure — DRO**, selected *by period*: `pumpkynPie.cr3bp.getDRO(tau)`
returns the distant retrograde orbit whose period equals `tau` in
nondimensional time. Here **τ = 1.0**, i.e. a period of **4.4327 days**.

**Arrival — Tulip**, `pumpkyn.cr3bp.getTulip(tau, Np, pm)` with
**τ = 5·2π/6, Np = 7, pm = −1**, giving a period of **23.2093 days**. The
tulip's period is exactly 5/6 of the Moon's revolution, so its seven-lobed
pattern closes after 6 tulip periods = 5 lunar revolutions.

Both `departure_params` and `arrival_params` carry a `reconstruction` string
with the literal pumpkyn calls needed to rebuild the orbit and evaluate its
state at any phase.

**Phase convention.** A phase is the fraction of an orbit period elapsed past
the reference point returned by `getDRO`/`getTulip`. Every entry stores it
three ways: `_frac` (0–1), `_nd` (nondimensional time), `_days`. Both axes
wrap, so the domain of phase pairs is a **torus**. The arrival axis is offset
by the anchor phase 0.075378 (the tulip's maximum-velocity-angle point); both
libraries use the same offset grid, so their cells correspond one-to-one.

**Grid.** 12 departure phases × 12 arrival phases = 144 pairs. Grid spacing is
8.9 hours along the DRO and 1.93 days along the tulip.

## 3. Physical constants and thruster

```
muStar = 0.012150585609624      lStar = 389703.264829278 km
tStar  = 382981.289129055 s     1 ND time = 4.432654 days
```

`tfMin` takes nondimensional thrust and exhaust velocity:

```
Tmax_nd = (thrust_N / m0_kg) * tStar^2 / (lStar_km * 1000)
c_nd    = (isp_s / tStar) * g0 ,   g0 = 9.80665 * tStar^2 / (1000 * lStar)
```

Mass is normalized so m(0) = 1; multiply by `m0_kg` for kilograms. `c_nd` is
stored directly in `lib.thruster`; `Tmax_nd` follows from each entry's own
`thrust_N` (v2 stores the formula in `thruster.thrust_nd_formula`).

**v1** is a single operating point: 0.07 N, Isp 900 s, 150 kg — an initial
acceleration of 0.467 mm/s², giving transfers of 16–33 days spanning several
DRO revolutions. **v2** covers 1–15 N at Isp 1710 s and 150 kg — 6.7 to
100 mm/s², giving sub-revolution transfers of 0.27–2.56 days.

The two libraries are **neighboring sheets, not points on one continuum**. The
gap between 0.07 N and 1 N contains several revolution-count transitions,
which is why v2 was built by anchoring at high thrust rather than by
continuing upward from v1.

## 4. Structure layout

Both files hold one struct, identical in shape apart from v2's thrust axis.

**Top level:** `name`, `description`, `created`, `provenance`, `constants`,
`thruster`, `departure_orbit`, `departure_params`, `arrival_orbit`,
`arrival_params`, `entries`, `n_entries`, `grid`, `usage`.

**`entries(k)`** — one converged solution each:

| Field | Meaning |
|---|---|
| `departure_phase_frac / _nd / _days` | where on the DRO the transfer starts |
| `arrival_phase_frac / _nd / _days` | where on the tulip it ends |
| `thrust_N`, `isp_s`, `m0_kg` | the propulsion case for *this* entry |
| `lambda0` [7×1] | initial costates [λ_r(3); λ_v(3); λ_m] |
| `tf_nd`, `tf_days` | minimum flight time |
| `z8` [8×1] | `[lambda0; tf_nd]` — the tfMin-ready vector |
| `ms_residual` | multiple-shooting residual at convergence |
| `flight_miss_km` | arrival error when the PMP law is flown end to end |

**`grid`** — the same data as lookup arrays: phase axes (fractions and days),
`thrust_N` (v2), `has_solution` logical mask, `tf_days`, and `entry_index`
pointing into `entries` (0 = no solution). v1's arrays are 12×12; v2's are
12×12×9.

## 5. How each was produced

Both follow the three-step pipeline in
`process/COSTATE_LIBRARY_PIPELINE.md`:

1. **Direct solve** — Hermite–Simpson collocation produces a trajectory *and*
   a costate trajectory (Hager covector mapping of the NLP defect
   multipliers). Verified by flying the solved control end to end.
2. **Refinement** — `ms_tfmin` (multiple shooting, K segments, analytic block
   Jacobian built from pumpkyn's own segment STMs) drives the PMP residual to
   ~1e-12. Collocation duals are only ~1e-3 accurate and single shooting
   amplifies that ~1000× over a multi-revolution arc, so this step is
   mandatory, not cosmetic.
3. **Acceptance** — `pumpkyn.cr3bp.tfMin` must return the entry unchanged.

They differ in how the grid was covered:

- **v1** swept the torus with a two-wave strategy — cold straight-line seeds
  opened isolated cells, then solutions flooded outward to neighbors — because
  at 0.07 N a cold solve is expensive and often fails.
- **v2** anchored each phase pair at **15 N**, where the transfer is nearly
  impulsive and converges cold in seconds, then walked thrust **down** the
  rung set with each rung warm-started from the one above. Continuation keeps
  a cell on one solution family; independent cold solves do not (two meshes at
  15 N were measured landing 50% apart in t_f, both feasible).

## 6. Verification, in numbers

| Check | v1 (126 entries) | v2 (990 entries) |
|---|---|---|
| Multiple-shooting residual | max 9.5e-11 | max 9.3e-11, median 8.7e-13 |
| Flown arrival miss | max 0.064 km | max 1.5e-5 km |
| `tfMin` acceptance ‖Δz‖ | max 2.8e-8, **126/126 pass** | max 1.2e-10, **990/990 pass** |

For v1 there is also a cross-method check: direct and indirect flight times
agree to a median of 2e-8 ND across all 126 phase pairs.

## 7. Coverage and known gaps

**v1** solves 126 of 144 phase pairs. The 18 missing are dominated by the
arrival row at phase 0.0754 — the tulip's fastest-moving point. Reaching it
after 4–7 revolutions is genuinely hard, and the sweep's seeds never got
there.

**v2** covers 130 of 144 pairs, with 49 carrying the complete nine-rung
ladder. Coverage by rung:

| Thrust | Pairs | t_f range (days) |
|---|---|---|
| 15 N | 109 | 0.269 – 0.661 |
| 12 N | 128 | 0.301 – 0.583 |
| 10 N | 128 | 0.330 – 0.642 |
| 7 N | 121 | 0.398 – 0.773 |
| 5 N | 125 | 0.483 – 0.919 |
| 3 N | 119 | 0.704 – 1.375 |
| 2 N | 106 | 0.923 – 1.553 |
| 1.5 N | 91 | 1.122 – 2.054 |
| 1 N | 63 | 1.546 – 2.560 |

Two patterns, both expected:

- **The 15 N anchor is cold-started**, so it fails more often (109) than the
  rungs below it (128). A failed anchor does not stop the ladder — the next
  rung down recovers.
- **The bottom rung is thin.** At 1 N the transfer is longest and hardest, and
  only 63 pairs make it. Filling this is the main outstanding item; the plan
  is an intermediate 1.25 N rung to soften the continuation step.

Notably, v2 *does* cover phase pairs v1 could not: the hard arrival row is
reachable at sub-revolution thrust levels, which confirms that its difficulty
was many-revolution geometry rather than anything intrinsic to the endpoints.

## 8. Using a library

```matlab
% v2: choose phases (days) and any thrust in range
L = load('costate_lib_dro_tulip_v2.mat');  lib = L.costate_lib_dro_tulip_v2;
[tf_days, e, bracket] = costate_lib_pick(lib, depDays, arrDays, thrustN);

% e.z8 is a converged seed at the nearest rung; tf_days is interpolated
% linearly between the bracketing rungs. Rebuild the orbits from
% lib.departure_params / lib.arrival_params, then either fly it:
[t, y] = pumpkyn.cr3bp.tfMinProp(e.tf_nd, [rv0(1:6)'; 1; e.lambda0], ...
                                 Tnd, lib.thruster.c_nd, mu);
% ...or hand the seed to the single-shooting solver (returns it unchanged):
z = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), e.z8, Tnd, lib.thruster.c_nd, mu);
```

v1 is identical with `costate_lib_nearest(lib, depDays, arrDays)` — no thrust
argument, since v1 is a single operating point.

For a thrust between rungs, `e.z8` is a warm start a fraction of a rung away;
`tfMin` closes the difference in about a second. For a phase between grid
points, likewise — grid spacing is 8.9 hours (DRO) and 1.93 days (tulip),
comfortably inside the convergence basin.

Browse rather than query with `lib.grid.has_solution` (what solves) and
`lib.grid.tf_days` (how long it takes); `lib.grid.entry_index` maps a grid
cell straight to its row in `lib.entries`.

## 9. Growth axes

The pipeline is parameterized, so further libraries are reruns rather than
rebuilds:

- **DRO period** — τ *is* the period, so each τ produces one library sheet;
  together they form the ΔV–period–phasing map.
- **Thrust** — additional rungs or ranges, as v2 demonstrates.
- **Orbit pair** — endpoints come from any of pumpkyn's catalogued families.
