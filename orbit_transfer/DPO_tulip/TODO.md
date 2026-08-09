# DPO_tulip — TODO

- [ ] **Densify the 65 unsolved pairs** (worst sheet τ=1 × Np=12 at
  27/36); neighbor-seeded continuation first, as for halo.
- [ ] **Period-axis growth:** the admissible box is τ 0.048–4.695 with 583
  native members; the catalog samples 4. The small-τ end (fast, tiny DPOs
  hugging the Moon) is unexplored below τ=1 and likely hard (winding);
  the τ≈4.5 end is nearly free.
- [ ] **Floor-rejected low rungs:** several cells stop at 1.5–2 N because
  the 1 N spiral dips under the 500 km floor. If Darin wants those rungs,
  the direct solve needs the altitude floor as an active path constraint
  (the machinery exists — `floorKm` — but the ladder currently stops
  rather than re-solving constrained).
- [ ] Diagnose (only if it recurs) the isolated tfMin-acceptance NaN
  pattern first seen in the HALO_HALO far-pair probe at 10 N.
