# HALO_tulip — L2-southern halo → tulip costate catalog

Second catalog campaign (2026-08-07/08), and the proof that the pipeline is
family-agnostic: the front door differs from the DRO catalog's ONLY in the
departure-family block — the engines (`DRO_tulip/indirect/
thrust_ladder_library`, `ms_tfmin`) and the shared library
(`costate_common`) run unmodified.

## Files

| file | what |
|---|---|
| `run_halo_catalog.m` | Front door: pilot / all / report stages; τ ∈ {1.75, 2.2, 2.8, 3.4} (L2 southern, pm = −1) × Np ∈ {5, 7, 9, 12}, 6×6 phasing, rungs 15→1 N. |
| `run_halo_batched.sh` | Hang-proof unattended driver (OS-timeout batches, lock, resume, abort-on-error, 600-pass limit). |
| `build_halo_catalog.m` | Packages the sheets via `costate_common/build_costate_catalog_family` → `costate_catalog_halo_tulip.mat`. |
| `direct/results/halo_bounds.mat` | Measured admissible box: L2 southern continuous 7.1–15.1 d (τ 1.61–3.42); L1 only two members; L3 outside. Northern branch = exact z-mirror. |
| `direct/results/catalog/` | 16 sheet .mats + progress logs. |
| `deliverables/costate_catalog_halo_deliverable4/` (+ .zip) | Shipped: catalog + catalog-agnostic pickers + worked example + extremes movie + README. Recipient self-tested. |

## Results

**3,980 entries, 530/576 pairs (92%), 242 full 9-rung ladders.**
ΔV 0.65–7.28 km/s; flight 0.68 d – ~3 d. Cheapest entry in ANY catalog so
far: τ=1.75 halo → Np=12 tulip at 1.5 N = **0.6546 km/s** (5.74 kg) —
halos sit energetically closer to the tulips than DROs (best DRO 0.984).
Solvability **improves** with halo period: τ=3.4 sheets run 35–36/36 pairs,
τ=1.75 runs 25–34/36. All entries three-gate verified; conjugate-point
(Jacobi necessary) check passes on the golden cell.

Ops note: the original run stopped at 84% by exhausting the driver's
200-pass loop limit (not an error); the limit is now 600 and resume is free.
