# orbit_transfer/doc — cross-campaign reference documents

Documents that span campaigns (campaign-specific notes live in each
campaign's own `doc/`). The live cross-campaign status + roadmap is
**`../STATUS_AND_ROADMAP.md`** (one level up, beside the program TODO).

| file | what |
|---|---|
| `transfer_problem_space.md` | The official goal map (2026-07-31): all pumpkyn-reachable orbit pairs × {direct, indirect} × {min-time, min-energy, min-fuel}; the nine catalogued families. |
| `library_catalog.md` (+ `gen_library_catalog.py`) | **Generated function reference for OUR library** (oclib/+oc, costate_common, verify_common, cr3bp_common): signature + purpose + I/O per function, 44 functions. Regenerate after any library change — never edit by hand. |
| `pumpkyn_catalog.md` (+ `gen_pumpkyn_catalog.py`) | Generated catalog of pumpkyn/pumpkynPie orbit families and getters. |
| `pumpkyn_reference.md` | Working reference for the pumpkyn API surface this repo uses. |
| `direct_vs_indirect.md` | The sensitivity-vs-combinatorics axis: when each method wins (min-time low-thrust favors direct; indirect for precision once seeded). |
| `campaign_status.tex/.pdf` | Point-in-time rollup, **SUPERSEDED 2026-08-23 by `../STATUS_AND_ROADMAP.md`** (it predates the entire costate-catalog program). Kept as a historical snapshot only. |
