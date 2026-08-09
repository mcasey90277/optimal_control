# orbit_transfer/doc — cross-campaign reference documents

Documents that span campaigns (campaign-specific notes live in each
campaign's own `doc/`).

| file | what |
|---|---|
| `transfer_problem_space.md` | The official goal map (2026-07-31): all pumpkyn-reachable orbit pairs × {direct, indirect} × {min-time, min-energy, min-fuel}; the nine catalogued families. |
| `pumpkyn_catalog.md` (+ `gen_pumpkyn_catalog.py`) | Generated catalog of pumpkyn/pumpkynPie orbit families and getters. |
| `pumpkyn_reference.md` | Working reference for the pumpkyn API surface this repo uses. |
| `direct_vs_indirect.md` | The sensitivity-vs-combinatorics axis: when each method wins (min-time low-thrust favors direct; indirect for precision once seeded). |
| `campaign_status.tex/.pdf` | Point-in-time campaign status rollup (check the date inside before trusting it). |
