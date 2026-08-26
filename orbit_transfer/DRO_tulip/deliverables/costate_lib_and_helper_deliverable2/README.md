# DRO → Tulip Costate Library with a Thrust Axis (deliverable 2)

A catalog of **converged** minimum-time low-thrust transfer solutions between a
DRO and a tulip orbit in the Earth–Moon CR3BP, indexed by departure phase,
arrival phase, and thrust level. Each entry is a root of the indirect (PMP)
boundary-value problem, not a hint: hand one to `pumpkyn.cr3bp.tfMin` and it
comes back unchanged in about a second.

## Contents

| File | What it is |
|---|---|
| `costate_lib_dro_tulip_v2.mat` | the library — 1,105 entries, 130 phase pairs × up to 14 thrust rungs |
| `costate_lib_pick.m` | look up by phases (days) + thrust (N); interpolates flight time between rungs |
| `costate_lib_v2_example.m` | worked example: pick → fly → plot → phasing torus |
| `plot_costate_torus.m` | flat and rotatable 3-D phasing torus, colored by coverage |
| `costate_lib_extremes.m` | cheapest / most expensive (or fastest / slowest) transfer, side-by-side plot |
| `costate_lib_extremes_example.m` | worked example for the above, plus a per-rung spread sweep |

Requires **pumpkyn / pumpkynPie** on the MATLAB path (for orbit families,
propagation, and `tfMin`). The library itself is plain data — no dependencies.

## Quick start

```matlab
% from this folder, with pumpkyn on the path:
costate_lib_v2_example        % pick a phasing + thrust, fly it, plot it
costate_lib_extremes_example  % cheapest vs dearest transfers
```

Or directly:

```matlab
L   = load('costate_lib_dro_tulip_v2.mat');
lib = L.costate_lib_dro_tulip_v2;
[tf_days, e] = costate_lib_pick(lib, 1.0, 10.0, 4.0);   % dep d, arr d, thrust N
e.z8          % [lambda_r(3); lambda_v(3); lambda_m; tf] -- tfMin-ready
e.deltaV_kms  % exact Delta-V for this transfer
```

## What the library covers

- **Departure:** DRO with τ = 1.0 — τ *is* the orbit period in nondimensional
  time, so this is the 4.4327-day DRO.
- **Arrival:** tulip with τ = 5·2π/6, Np = 7, pm = −1; period 23.2093 days.
- **Propulsion:** Isp 1710 s, initial mass 150 kg.
- **Thrust rungs (N):** 15, 12, 10, 7, 5, 3, 2, 1.5, 1, 0.9, 0.8, 0.7, 0.6, 0.5
- **Phasing grid:** 12 × 12 = 144 phase pairs; 130 have at least one entry.
- **Flight times:** 0.269 – 5.088 days.  **ΔV:** 0.915 – 6.979 km/s.

Coverage thins at the low-thrust end (128 phase pairs at 10 N, 63 at 1 N, 20 at
0.5 N) because the transfer approaches one full DRO revolution there and the
solution family reorganizes. `lib.grid.has_solution(:,:,k)` is the
authoritative answer for what exists at rung *k*; the torus plot shows it at a
glance.

## Structure of an entry

```
departure_phase_frac / _nd / _days    where on the DRO the transfer starts
arrival_phase_frac   / _nd / _days    where on the tulip it ends
thrust_N, isp_s, m0_kg                the propulsion case for THIS entry
lambda0 [7x1]                         initial costates
tf_nd, tf_days                        minimum flight time
z8 [8x1]                              [lambda0; tf] -- pass straight to tfMin
deltaV_kms, m_final_kg, propellant_kg exact (all-burn) mass budget
ms_residual, flight_miss_km           how it was verified
```

`lib.grid` repeats the same information as lookup arrays over
(departure × arrival × thrust): `has_solution`, `tf_days`, `deltaV_kms`, and
`entry_index` into `lib.entries`.

Both orbits are reconstructible from the file alone — `lib.departure_params`
and `lib.arrival_params` carry the parameters *and* the literal pumpkyn calls
needed to rebuild them.

## How the entries were made, and how they were checked

Three steps per case: a direct collocation solve supplies a trajectory and a
costate trajectory; multiple shooting (`ms_tfmin`, analytic Jacobian from
pumpkyn's own segment STMs) refines the costates to ~1e-12; and
`pumpkyn.cr3bp.tfMin` must accept the result unchanged. Thrust levels are
covered by a ladder — anchored at high thrust, where the transfer is nearly
impulsive, then walked downward with each rung warm-starting the next — so a
phase pair stays on one solution family.

Across all 1,105 entries: multiple-shooting residual ≤ 9.3e-11, flown arrival
miss ≤ 7.3e-5 km, and `tfMin` acceptance ‖Δz‖ ≤ 6.8e-10 (1,105 of 1,105).

## Two things to know before relying on it

1. **At a fixed thrust, minimum flight time and minimum ΔV are the same
   transfer.** These solutions burn continuously, so
   ΔV = c·ln(1/(1 − T·t_f/c)) increases monotonically with t_f. The metrics
   differ only when comparing across thrust levels.
2. **Phasing matters far more at high thrust than at low.** The spread between
   best and worst phasing is 4.47 km/s at 15 N but under 0.4 km/s near 0.9 N —
   a slow transfer has time for the geometry to come to it.

M. Casey / D. Koblick, Coorbital Inc. — August 2026
