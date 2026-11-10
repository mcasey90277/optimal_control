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
orbit_network_study         % topic 2: the network, the coverage and envelope planes, relay closure
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
| `load_catalog_sheets.m` | the four catalogs -> one 56-row sheet table (v1 DRO catalog mapped onto the v2 names; the `accepted` mask = solved, conjugate-pass (== 1; the stored -1 "not run" is not a pass) and the three hypothesis gates; refuses a missing file, mismatched rungs, a wrong grid shape) |
| `orbit_node_table.m` | the 18 distinct orbits: period, Jacobi constant, monodromy spectrum, sigma = ln\|lambda\|max / T (exactly 0 when \|lambda\|max < 1 + 1e-3), degrees; cached with provenance (keys, muStar, parameters) |
| `stability_cost_study.m` | topic 1: per sheet, dCJ, sigma_dep, sigma_arr against min t_f at 1 N and retention of the accepted phasing grid to 0.5 N; nested regression with a partial F test, per-coefficient t tests and 95% CIs; P1-P7 code gates, H1-H2 verdicts (DETECTED / NOT DETECTED, nominal OLS); one figure |
| `lib/partial_f_test.m` | nested OLS, partial F, p-values through betainc (no toolbox), per-coefficient t tests and CIs; refuses rank-deficient or non-nested designs; both-exact gives p = 1 |
| `orbit_network_study.m` | topic 2: digraph of 18 orbits and 56 edges over ACCEPTED cells with per-rung accepted fraction, min/median t_f and the lowest accepted rung; raw-coverage and known-feasible-envelope planes over (thrust bound, time cap) with the directed-reachability overlay; in/out strengths at 1 N with the topology caveat; L1 halo -> L2 halo -> tulip relay closure with PHASE-COMPATIBLE assembled times; the hypothetical reverse graph (off by default, geometry only); G1-G5 gates; four figures and `results/network_summary.txt` |
| `lib/edge_tables.m`, `lib/network_graph.m`, `lib/accessibility_plane.m`, `lib/relay_closure.m`, `lib/mirror_state.m`, `lib/symmetry_residual.m` | the network's computations, each with a demo and a hand-checked test in `tests/test_lib.m` |
| `doc/orbit_network_note.tex` | the topic-2 note (pdflatex; the PDF is git-ignored) |
| `slides/orbit_network_slides.tex` | 18-slide beamer deck over both topics (pdflatex; the PDF is git-ignored): six plain-language slides first (one analogy throughout, a train timetable: stations = orbits, trips = transfers, changing trains = relay, platform wait = phasing; TikZ diagrams and a plain-words-to-jargon table), then the 12 technical slides. `slides/make_slide_figures.m` runs both study scripts headless (outputs to a scratch folder, never `results/`) and redraws the figures at slide size from the same workspaces, plus a 3-D plot of the 18 node orbits; PNGs in `slides/figures/` are committed |
| `reviews/` | the GPT-6 Astra review of the spec and plan (2026-11-05) and its dispositions: what was accepted and what it changed |
| `tests/` | `test_load_catalog_sheets`, `test_orbit_node_table`, `test_lib`, `test_pipeline` (scripts headless, parse check), `run_tests` |

Both scripts print every gate PASS/FAIL and throw on a failed code gate;
hypothesis verdicts (H1-H2) are results either way. `results/network_summary.txt`
holds every network number the note cites; the stability numbers come from stability_cost_study's output.

## What it does not claim

- Relay-versus-direct cost is not testable on this data: no pair has both a
  direct edge and a two-leg route. The lemma in the note is about global
  infima: the feasible set of the unrestricted problem contains every relay
  route with compatible full states (mass included), so its infimum is never
  above the relay infimum. That says nothing about what a particular solver
  returned, and it does not dismiss Braik and Ross's result, which is
  relative to their restricted direct-maneuver class.
- Closure is monotone between 1 N and 0.5 N (the two rungs the study uses)
  but not along the whole ladder: the 15 N starting rung closes fewer cells
  than 12 N. That is solver coverage, printed as a fact, never gated. An
  unaccepted cell or rung is unresolved, never a demonstrated infeasibility.
- sigma_arr varies only on the 8 halo-to-halo sheets, so it is confounded
  with catalog; the verdicts rest on the departure term. There are 14
  distinct departure orbits and 8 arrivals, so sheets share nodes and the
  OLS p-values are nominal; NOT DETECTED is not "no effect", and 1 N and
  0.5 N sit far above the 70 mN regime of the DRO wall.
- pumpkynPie's `orbitProperties` errors on a stable orbit, so sigma is taken
  from `pumpkyn.cr3bp.monodromy` directly.
- An unaccepted cell or rung is unresolved, never a demonstrated infeasibility;
  the "lowest accepted rung" is not a thrust floor.
- Weak connectivity and betweenness are not reported: only L2 halos have both
  in- and out-edges, so those statistics would restate the catalog's topology.
- The hypothetical reverse graph is geometry: a mirrored thrusting arc belongs
  to a mass-GAIN model and carries no certificate.
