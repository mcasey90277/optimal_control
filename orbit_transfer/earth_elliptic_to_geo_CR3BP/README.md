# earth_elliptic_to_geo_CR3BP — elliptic → GEO with lunar gravity

**Goal.** Re-solve the low-thrust minimum-fuel elliptic-orbit → GEO transfer
(the Haberkorn–Martinon–Gergaud benchmark: 1500 kg, P⁰=11625 km, e⁰=0.75,
i⁰=7° → equatorial GEO, thrust ladder 10 → 0.1 N) **with the Moon's gravity
incorporated** — Earth–Moon CR3BP dynamics instead of the 2-body `1/r²` model
used in `../earth_elliptic_to_geo/`.

**The question:** how much does lunar gravity move the certified 2-body
answers — final mass, switch structure, and the thrust-ladder laws?
Comparison baselines (from `../earth_elliptic_to_geo/`, all certified):

| baseline (2-body) | value |
|---|---|
| 10 N, c_tf=1.5 | m_f = 1377.10 kg, 19 switches, 7.33 rev |
| 0.2 N (mesh-converged, P0 study) | m_f = 1375.8 kg, ~866±5 switches, 346.7 rev |
| R0 law | T·t_f,min ≈ 850 N·h across two thrust decades |

## Status

**Not started** — this folder holds only README + TODO. Structure below is the
plan, mirroring the sibling campaigns.

## Planned structure

```
earth_elliptic_to_geo_CR3BP/
├── README.md, TODO.md
├── direct/          (collocation NLP with CR3BP dynamics)
└── indirect/        (PMP shooting counterpart)
```

## Design starting points (see TODO for the open decisions)

- **Warm-start bridge:** the certified 2-body solutions in
  `../earth_elliptic_to_geo/` are the natural seeds. The proven pattern is the
  **two-primary gravity homotopy** in
  `../GTO_ELFO/direct/elfo/gen_elfo_energy_gravhom.m` — dial the Moon's mass
  from 0 to μ*, warm-starting each step — which bridges exactly this gap
  (Earth-only solution → CR3BP solution). This is precisely the
  **μ-continuation** of Bonnard–Caillau–Picot 2010 (`papers/Geometric_And_
  Numerical_Techniques_In_3_Body_Low_Thrust_Transfers.pdf`, §4.4; traced to
  Poincaré), which also supplies its convergence theory: the continuation
  branch stays smooth exactly while the extremal stays clear of conjugate
  points (their Prop 2.3).
- **Dynamics/formulation:** the 2-body campaign's winning MEE/L-domain
  formulation assumes Earth-centered Keplerian structure; with the Moon it
  either gains a perturbation term (MEE Gauss + lunar acceleration) or is
  replaced by rotating-frame Cartesian like `../GTO_tulip/`'s engine. Open
  trade — see TODO.
- **Terminal set:** equatorial GEO is a fixed 5-element target in MEE, but in
  the Earth–Moon rotating frame it is a moving circle; the terminal manifold
  must be re-posed for whichever frame wins.
- **Shared code:** `../cr3bp_common/` holds the CR3BP constants/pattern
  (note: its `cr3bp_lt_params` is the tulip craft — 15 kg / 25 mN; this
  problem's craft is 1500 kg / 10..0.1 N, so it needs its own params file,
  possibly promoted into `cr3bp_common` with craft as an argument).
