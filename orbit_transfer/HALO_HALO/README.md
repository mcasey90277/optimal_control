# HALO_HALO — L1 ↔ L2 halo-to-halo transfers (probe stage)

Darin's ask (2026-08-07): transfers between L1 and L2 halos — "visually it
will look really nice." Feasibility is PROVEN: the shared pipeline needed
**zero changes** (the ladder engine takes departure and arrival families
independently), and both probe pairs solved through all three gates.

## Files

| file | what |
|---|---|
| `probe_l1_l2_halo.m` | One-cell probe through the full pipeline; parameterized by (τ_L1, τ_L2); outputs keyed by the pair. |
| `render_l1_l2_movie.m` | Pumpkyn-style movie (SatelliteAnimator dark scene, exact 1280×720) of a probe rung; parameterized by rung + probe file. |
| `direct/results/l1l2/` | Probe sheets + movies. |

## Probes solved

| pair | character | rungs passed | ladder |
|---|---|---|---|
| L1 τ=1.8037 → L2 τ=2.2 | near/NRHO-like, dives past the Moon | **7/7** (15→1 N) | 0.52 d / 5.18 km/s @15 N → 2.05 d / 1.23 km/s @1 N |
| L1 τ=2.7433 → L2 τ=3.40 | far/circular, hugs the libration points (Darin's preferred look) | 6/7 | 0.76 d @15 N → 3.11 d @1 N; lone 10 N tfMin-acceptance NaN (flown 0.08 km — watch item) |

Movies: `l1_l2_halo_1N.mp4/.gif` (near pair), `l1_l2_halo_2p743_3p400_1N.mp4/.gif`
(far pair).

## Constraints for the full catalog

Only **two** admissible L1 members exist under the standing criteria
(τ = 1.8037 at 2,555 km periselene; τ = 2.7433 at 12.16 d) — the L1 family
mostly dips below the 500 km floor. The box is therefore 2 L1 × 11 L2
members. Blocking schema item: the catalog format keys arrivals by petal
count `Np`; a halo arrival needs an **arrival-period axis** (the SDD's
"new parameter axis" extension step).
