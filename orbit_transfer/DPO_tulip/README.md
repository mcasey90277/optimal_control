# DPO_tulip — distant-prograde-orbit → tulip costate catalog

Third catalog campaign (2026-08-08), completing Darin's roadmap (richer DRO
→ Halo → DPO). Same shared engines and pipeline as the halo campaign; the
front door differs only in the departure-family block.

## Files

| file | what |
|---|---|
| `run_dpo_bounds.m` | Admissibility survey over pumpkyn's NATIVE getDPO members (1,035 seeds): **583 admissible, τ 0.048–4.695 ND (0.21–20.8 d)**; the upper edge is the 100 Mm criterion (DPOs grow with period); the τ ≳ 9 tail fails the periodicity guard. → `direct/results/dpo_bounds.mat`. |
| `run_dpo_catalog.m` | Front door: τ ∈ {1, 2, 3, 4} × Np ∈ {5, 7, 9, 12}, 6×6 phasing, rungs 15→1 N; pilot / all / report. |
| `run_dpo_batched.sh` | Hang-proof unattended driver (this campaign's runs used it exclusively after the one-shot pilot process died silently — buffered-stdout crash mode). |
| `build_dpo_catalog.m` | Packages sheets → `costate_catalog_dpo_tulip.mat`. |
| `deliverables/costate_catalog_dpo_deliverable5/` (+ .zip) | Shipped: catalog + catalog-agnostic pickers + worked example + extremes movie + README. Recipient self-tested (flown miss 0.000 km, tfMin |Δz| = 0). |

## Results

**3,932 entries, 511/576 pairs (89%), 238 full ladders.**
ΔV 0.76–7.10 km/s; flight 0.24–3.58 d. Cheapest: τ=1 DPO → Np=12 tulip at
1 N = 0.7600 km/s (6.65 kg). The hard corner is the smallest DPO against
the longest-period tulip (τ=1 × Np=12: 27/36), mirroring the halo pattern
(short-period departures are hard). Some low rungs are rejected by the
500 km altitude floor as spirals dip toward the Moon — a designed stop and
one source of partial ladders. Integer-ND periods are incommensurate with
the lunar and tulip periods, so no exact-resonance degeneracies (review
observation, Gemini).

## Conjugate-point sweep (2026-08-23)

All 3,932 entries carry a Jacobi (conjugate-point) verdict in the catalog's
`conj_pass` grids: **3,931 pass, 1 FAIL** — the τ=2 → Np=7 sheet at phases
(2/3, 2/3), 15 N (tf = 0.115 ND, ΔV 5.12 km/s) has one interior sign change
(re-solve fidelity 1.6e-12) and is **refuted as a local minimum**: a cheaper
neighboring transfer exists for that cell. Its lower-thrust rungs pass. It
still flies and tfMin still accepts it, but continuation seeded from it walks
a non-minimizing branch. Sweep engine: `costate_common/conj_catalog_pass`.
