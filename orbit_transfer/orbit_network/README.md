# orbit_network -- stability versus cost, and the low-thrust orbital network

Two study lines lifted from Braik and Ross, "Orbital Networks in the
Three-Body Problem" (arXiv 2605.31543, 2026) onto the certified min-time
costate catalogs this repo already holds. No proxies, no new solves: every
number comes from `DRO_tulip`, `HALO_tulip`, `DPO_tulip` and `HALO_HALO`
(56 sheets, 6 x 6 phasing grids, thrust ladder 15 N to 0.5 N). Design:
`docs/superpowers/specs/2026-11-05-orbit-network-design.md`.

## Run

```matlab
cd orbit_transfer/orbit_network
stability_cost_study        % topic 1: does orbit stability predict transfer cost?
cd tests; run_tests         % loader, node table, lib helpers, scripts end to end
```

The shared MATLAB session already has pumpkyn and pumpkynPie on the path
(the startup adds them); a clean `-batch` run calls pumpkynPie's startup
first. Scripts add `lib/` and `../costate_common` themselves. Outputs go to
`results/` (git-ignored): figures and the node-table cache
`orbit_nodes.mat`.

## What is here

| file | role |
|---|---|
| `load_catalog_sheets.m` | the four catalogs -> one 56-row sheet table (v1 DRO catalog mapped onto the v2 names; refuses a missing file, mismatched rungs, a wrong grid shape) |
| `orbit_node_table.m` | the 18 distinct orbits: period, Jacobi constant, monodromy spectrum, sigma = ln|lambda|max / T (0 when stable), degrees; cached |
| `stability_cost_study.m` | topic 1: per sheet, dCJ, sigma_dep, sigma_arr against min t_f at 1 N and survival of the phasing grid to 0.5 N; nested regression with a partial F test and per-coefficient t tests; P1-P7 code gates, H1-H2 hypothesis verdicts; one figure |
| `lib/partial_f_test.m` | nested OLS, partial F, p-values through betainc (no toolbox), per-coefficient t tests |
| `tests/` | `test_load_catalog_sheets`, `test_orbit_node_table`, `test_lib`, `test_pipeline` (scripts headless, parse check), `run_tests` |

Topic 2, the network itself (`orbit_network_study.m`, the accessibility
plane, relay closure, the note in `doc/`), follows the spec's sections 7-10
and is not yet built.

## What it does not claim

- Relay-versus-direct cost is not testable on this data: no pair has both a
  direct edge and a two-leg route. In a genuine optimal control problem the
  direct feasible set contains every relay route, so direct never loses.
- Closure is monotone between 1 N and 0.5 N (the two rungs the study uses)
  but not along the whole ladder: the 15 N starting rung closes fewer cells
  than 12 N. That is solver coverage, printed as a fact, never gated.
- sigma_arr varies only on the 8 halo-to-halo sheets, so it is confounded
  with catalog; the verdicts rest on the departure term.
- pumpkynPie's `orbitProperties` errors on a stable orbit, so sigma is taken
  from `pumpkyn.cr3bp.monodromy` directly.
