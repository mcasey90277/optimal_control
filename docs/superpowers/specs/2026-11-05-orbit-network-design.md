# Orbit network over the costate catalogs: design

Date: 2026-11-05. Status: approved in conversation (Mike, 2026-11-05);
revised the same day after the GPT-6 Astra review
(`orbit_transfer/orbit_network/reviews/astra_review_2026-11-05.md`, with
dispositions beside it). The plan argues from this document.

## 1. Purpose

Braik and Ross (arXiv 2605.31543, "Orbital Networks in the Three-Body
Problem") organize periodic orbits of the planar Earth-Moon CR3BP into a
weighted graph whose edges come from reachable-set overlap under
heading-change impulses on one Jacobi level, then read hubs, relays and
budget regimes off the graph. Their maneuver model is not ours. Their
multileg-cheaper headline is relative to their restricted direct-maneuver
class: in an optimal control problem whose feasible set already contains
every relay route, the global direct infimum is never above the relay
infimum (the lemma in section 10). What survives for us is the framing.
This folder lifts two of their ideas onto the certified min-time catalogs
we already hold, using NO proxies and NO new solves:

1. **Stability versus cost** (their section 6 reading: stable orbits are
   the hard nodes, unstable ones share a chaotic sea). Tested at SHEET
   level across the four catalogs.
2. **The low-thrust orbital network** (their sections 5-6). Nodes are
   pumpkyn orbits, edges are accepted catalog sheets, and the edge
   attributes are catalog coverage across thrust: the lowest accepted rung
   and a phase-resolved time of flight. Their (dV cap, time cap) budget
   plane becomes our (thrust capability, time cap) plane.

Decisions already made (do not reopen):

- Sheet level, not family level: every (departure orbit, arrival orbit)
  sheet is one point, 56 in all.
- Data-only: the 56 existing edges, thrust 15 N down to 0.5 N. The legs
  that would give DRO-to-tulip both a direct edge and a relay route
  (DRO-to-DPO, DRO-to-L2 halo, and their deep rungs) are a LATER campaign;
  this folder tells us which relay orbit is worth that ladder time.
- Deliverable form: a study script per topic plus one LaTeX note. No query
  tool/API until the missing legs exist.
- Relay-versus-direct cost is NOT testable on this data (no pair has both a
  direct edge and a two-leg route). The note carries the lemma; the script
  computes relay CLOSURE, and the phase-compatible assembled route time.
- **Accepted cells.** Every count, minimum, plane and route uses
  `accepted = has_solution & conj_pass == 1 & gate_min_lamv > 0 &
  gate_min_qmt > 0 & gate_dimS == 1` (the certified witnesses). The solved
  count is reported beside it. `conj_pass` stores -1 for "not run", so a
  `logical()` cast is wrong and never used.
- **The hypothetical reverse graph.** The x-axis symmetry
  (x, y, z, vx, vy, vz, t) -> (x, -y, z, -vx, vy, -vz, -t) maps a
  thrusting arc to a thrusting arc of a mass-GAIN model (thrust direction
  reflected by diag(1, -1, 1); mass runs m_f -> m_0). Reversed edges are
  therefore not transfers of the physical reverse problem at all, not even
  approximately. The option is off by default, named HYPOTHETICAL, carries
  no certificate, and its numbers never enter the headline sections.
- sigma = ln|lambda|_max / T from `pumpkyn.cr3bp.monodromy`, set to
  exactly 0 when |lambda|_max < 1 + 1e-3 ("no resolved exponential
  instability in the computed multipliers"; the tolerance is stated
  wherever sigma is reported). `pumpkynPie.cr3bp.orbitProperties` errors
  on stable orbits and is not used.
- What a closed or unclosed cell means: an accepted cell is a feasible
  witness; an unclosed cell or rung is unresolved, never a demonstrated
  infeasibility. The lowest accepted rung is reported as such, not as a
  floor or a bracket.

## 2. What the data is

Four min-time catalogs (schema v1/v2), all with `rungs_N = [15 12 10 7 5 3
2 1.5 1 0.75 0.5]` and 6 x 6 phasing grids per sheet, phases
`sD_frac = sA_frac = (0:5)/6` of the orbit period from the orbit's
`get_family_orbit` initial state (the same origin in every catalog):

| catalog | file | sheets | departure orbits | arrivals |
|---|---|---|---|---|
| DRO_tulip | `DRO_tulip/direct/results/costate_catalog_dro_tulip.mat` | 16 | DRO tau 0.5, 1, 2, 3 | tulip Np 5, 7, 9, 12 (pm -1) |
| HALO_tulip | `HALO_tulip/direct/results/costate_catalog_halo_tulip.mat` | 16 | L2 south halo tau 1.75, 2.2, 2.8, 3.4 | same tulips |
| DPO_tulip | `DPO_tulip/direct/results/costate_catalog_dpo_tulip.mat` | 16 | DPO tau 1, 2, 3, 4 | same tulips |
| HALO_HALO | `HALO_HALO/direct/results/costate_catalog_halo_halo.mat` | 8 | L1 south halo tau 1.8037, 2.7433 | the four L2 halos |

Sheet fields used: `has_solution`, `tf_nd`, `conj_pass` (1 / 0 / -1),
`gate_min_lamv`, `gate_min_qmt` (margins, NaN where unsolved),
`gate_dimS` (int8, 1 = pass, -1 where unsolved), each 6 x 6 x 11;
`sD_frac`, `sA_frac`, `period_tulip_nd` (= arrival period), and the keys
`tauDRO`/`Np`/`pm` (v1, DRO catalog) or `dep_family`/`dep_params`/
`arr_family`/`arr_params`/`tau_dep`/`tau_arr` (v2). Catalog-level
`constants` carries muStar, lStar_km, tStar_s.

Facts established 2026-11-05 that the gates pin:

- Closure is monotone between the two rungs the study uses: zero cells
  closed at 0.5 N but open at 1 N, in every catalog. Across the WHOLE
  ladder it is not: the 15 N starting rung closes fewer cells than 12 N
  (462 vs 494 in DRO_tulip; 593 non-monotone adjacent-rung cells over all
  four catalogs). That is solver coverage at the ladder's top, not a data
  error; scripts print it as a fact and never gate on whole-ladder
  monotonicity.
- Accepted counts per rung over all 56 sheets are [1578 1688 1674 1640
  1622 1570 1512 1381 1118 1069 898] against solved [1579 1688 1674 1640
  1623 1572 1516 1384 1124 1090 940]: a few cells are refuted or gated at
  every rung, most at 0.75 and 0.5 N. Four cells are accepted at 0.5 N
  but not at 1 N (refuted at the shallower rung), so the ACCEPTED sets are
  not nested between rungs even where the SOLVED sets are; retention is
  therefore defined as accepted-at-both over accepted-at-1 N, and the
  exceptions are printed.
- Every sheet has at least 9 cells closed at 1 N (range 9-34 of 36).
- Node stability (|lambda|_max of the monodromy matrix, sigma in 1/ND):
  all four DROs, all four tulips, DPO tau 1 and L2 halo tau 2.2 are
  linearly stable (|lambda|_max within 3e-4 of 1); L2 halos tau 1.75 / 2.8
  / 3.4 have sigma 0.63 / 1.18 / 2.05; L1 halos tau 1.8037 / 2.7433 have
  sigma 0.88 / 2.83; DPOs tau 2 / 3 / 4 have sigma 1.22 / 1.71 / 1.67.
  So sigma_arr is zero on 48 of 56 sheets; the live contrast is departure
  stability: 24 sheets leave a stable orbit, 32 an unstable one. There are
  only 14 distinct departure orbits and 8 distinct arrivals; sheets share
  nodes and are not independent observations.
- Every stored initial state lies on the symmetry plane (y = vx = vz = 0),
  so a mirror check at the initial state is vacuous; symmetry is checked
  at an off-plane point of the orbit.
- MATLAB `digraph` sorts edges by (source, target); attributes handed in
  beside the edge table do not follow. Every attribute is read from
  `G.Edges`.

## 3. Layout

```
orbit_transfer/orbit_network/
  README.md                     what the folder is, how to run, what it does not claim
  load_catalog_sheets.m         the four catalogs -> one 56-row sheet table (with `accepted`)
  orbit_node_table.m            the 18 orbits with period, C_J, |lambda|max, sigma, x0; cached with provenance
  stability_cost_study.m        topic 1 study script (sections 0-5, gates, figure)
  orbit_network_study.m         topic 2 study script (sections 0-7, gates, figures, summary file)
  lib/
    partial_f_test.m            nested linear models -> F, p, per-coefficient t, CI (betainc; no toolbox)
    mirror_state.m              x-axis symmetry R(x) = (x, -y, z, -vx, vy, -vz)
    symmetry_residual.m         distance of the mirror image of a state from the orbit through it
    edge_tables.m               per-edge, per-rung attributes over accepted cells
    network_graph.m             the digraph with a sheet_index edge variable (edge order safe)
    accessibility_plane.m       (thrust rung, time cap) -> pair and cell counts, raw and envelope
    relay_closure.m             two-leg routes with the phase-compatible assembled time
  tests/
    test_load_catalog_sheets.m
    test_orbit_node_table.m
    test_lib.m                  every lib helper against hand-computed answers
    test_pipeline.m             both scripts headless (plus the hypothetical-graph run); parse check
    run_tests.m                 the folder's runner, one logical verdict
  reviews/                      the Astra review and its dispositions
  results/                      git-ignored: figures, network_summary.txt, the node-table cache
  doc/
    orbit_network_note.tex      the note (+ .pdf, figures/)
```

Paths: scripts add `here`, `here/lib`, and `../costate_common` (for
`get_family_orbit`); pumpkyn/pumpkynPie come from the shared session path
(the startup adds them; a clean `-batch` run calls pumpkynPie's startup
first, as `transfer_study.m` does). No campaign folder goes on the path.

House style throughout: header quartet, `if nargin == 0` demo on every
function, no `i`/`j`, no `norm`, no `%#ok`, no emojis, ASCII README; study
scripts in the numbered-section style with every condition computed in the
script and reported PASS/FAIL; `\author{Michael Casey}` only in the note.

## 4. `load_catalog_sheets(root)` -> `T`

`root` is the `orbit_transfer` folder (`[]` for the parent of this file;
a bare call runs the demo). Returns a 56 x 1 struct array, one row per
sheet, catalog order DRO_tulip, HALO_tulip, DPO_tulip, HALO_HALO, sheet
order as stored. Fields:

| field | type | meaning |
|---|---|---|
| catalog | char | 'dro_tulip' \| 'halo_tulip' \| 'dpo_tulip' \| 'halo_halo' |
| sheet_index | double | position in that catalog's `.sheets` |
| dep_family, arr_family | char | 'dro' \| 'dpo' \| 'halo' \| 'tulip' |
| dep_params, arr_params | struct | exactly the `get_family_orbit` parameter struct |
| dep_key, arr_key | char | canonical node key, e.g. `dro_tau1`, `halo_L2_tau2.2`, `tulip_Np7`, `halo_L1_tau1.8037` |
| T_dep, T_arr | double | periods (ND) as the catalog stores them |
| rungs_N | 1 x 11 | copied from the catalog |
| has_solution | 6 x 6 x 11 logical | solved |
| conj_pass | 6 x 6 x 11 logical | `conj_pass == 1` (the -1 "not run" code is NOT a pass) |
| accepted | 6 x 6 x 11 logical | `has_solution & conj_pass & gate_min_lamv > 0 & gate_min_qmt > 0 & gate_dimS == 1` |
| tf_nd | 6 x 6 x 11 | copied |
| sD_frac, sA_frac | 1 x 6 | copied |
| constants | struct | muStar, lStar_km, tStar_s from the catalog |

The v1 DRO catalog is mapped onto the v2 field names here so no consumer
ever branches on schema version. Refusals: a missing catalog file
(`:missingCatalog`, names the path), rung vectors that differ between
catalogs (`:rungs`), and a sheet whose arrays are not 6 x 6 x numel(rungs)
(`:shape`).

## 5. `orbit_node_table(T, opts)` -> `N`

The distinct `dep_key`/`arr_key` values of `T`, 18 of them, each built once
with `get_family_orbit(family, params)` (first output the time grid,
second the 325 x 6 state history, third `info` with `periodND`; the initial
state is row 1). Fields per node: `key`, `family`, `params`, `period_nd`,
`x0` (6 x 1), `CJ` (`pumpkyn.cr3bp.jacobi`), `lam_max`, `lam_min`, `detM`
(from `pumpkyn.cr3bp.monodromy(x0.', period_nd, muStar)`), `stable` =
lam_max < 1 + 1e-3, `sigma` = 0 when stable else log(lam_max)/period_nd,
`in_degree`, `out_degree` (counts over `T`). The loader's continuation
prints a log; wrap in `evalc`.

`opts.cache` (default `results/orbit_nodes.mat` beside this file): the
file stores `N` and `meta` (muStar, the keys, the parameter structs, the
date). It is used only when the keys, muStar and every parameter struct
equal the request; `opts.rebuild = true` forces a rebuild, and the note's
numbers are produced from a rebuilt table. Gate G1 compares `period_nd`
with the catalog's `T_dep`/`T_arr` for every row of `T` to 1e-9.

## 6. `stability_cost_study.m` (topic 1)

Sections:

0. **settings**: paths, root, the two rungs (1 N as the shared cost rung,
   0.5 N as the deepest rung), the nominal significance level 0.05, figure
   folder.
1. **load and census**: `T`, `N`; print the 56-row table (catalog, dep,
   arr, accepted cells at 1 N and 0.5 N, solved cells beside them, min t_f
   at 1 N in days). Gates: P1 numel(T) == 56; P2 identical rungs across
   catalogs; P3 every cell SOLVED at 0.5 N is solved at 1 N (the ladder
   descends from 1 N), with the whole-ladder non-monotone count and the
   accepted-at-0.5-N-but-not-at-1-N count printed as facts; P4 every
   sheet has at least one accepted cell at 1 N.
2. **orbit properties**: for each node print period, C_J, |lambda|max,
   sigma, stable (tolerance 1e-3 stated). Gates: P5 symplectic pairing
   lam_max*lam_min = 1 to 1e-6 and |det M - 1| < 1e-6; P6 node periods
   match the catalogs to 1e-9 (G1); P7 the arrival nodes of every tulip
   sheet are stable (a fact check).
3. **covariates and responses** per sheet: dCJ = CJ_arr - CJ_dep;
   sigma_dep, sigma_arr; tf1 = min over cells ACCEPTED at 1 N of tf_nd, in
   days; surv = (accepted at both 1 N and 0.5 N) / (accepted at 1 N), a
   nested ratio by construction. Print the number of distinct departure
   (14) and arrival (8) orbits.
4. **regression**: for each response y in {tf1, surv}: model M0 y ~ 1 +
   dCJ; model M1 y ~ 1 + dCJ + sigma_dep + sigma_arr; least squares; the
   partial F for M1 over M0, and the per-coefficient t statistics, two-
   sided p-values and 95% confidence intervals from `partial_f_test`;
   group means of y for stable-departure versus unstable-departure sheets,
   labelled descriptive. Verdicts, printed as DETECTED / NOT DETECTED
   (nominal OLS), rest on the sigma_dep coefficient's own t-test and sign:
   H1 for tf1 (paper's sign negative, unstable departures faster); H2 for
   surv (sign positive). The joint F and the sigma_arr test print beside
   each verdict with the caveat that sigma_arr varies only on the 8
   halo-to-halo sheets. Both outcomes are results; NOT DETECTED is not a
   failure and is never written as "no effect".
5. **figure**: one PNG, two panels, each response against dCJ, markers
   colored by sigma_dep (stable group in one color), with the M0 fit line;
   written to `results/stability_cost.png`.

The summary block lists P1-P7 PASS/FAIL and H1-H2 verdicts and throws if
any P gate fails.

The note's wording for the result (section 10) says what was and was not
detected, that sheets share nodes and the p-values are nominal, that the
1 N and 0.5 N rungs sit far above the 70 mN regime of the DRO wall, and
that nonsignificance is not equivalence.

## 7. `orbit_network_study.m` (topic 2)

Sections:

0. **settings**: paths, the reference rung (1 N, every sheet has an
   accepted cell there), the deepest rung (0.5 N), the number of time
   caps (25), `hypothetical = false` (the reverse graph), the relay key
   prefixes (`halo_L1` -> `halo_L2` -> `tulip`), figure folder.
1. **load**: `T`, `N`; the per-rung solved and accepted counts, the
   whole-ladder non-monotone count and the accepted-at-0.5-N-but-not-at-1-N
   count printed as facts. Gates: G1 node
   periods match the catalogs to 1e-9; G2 56 edges, every node key found
   in the graph by `findnode`, zero cells SOLVED at 0.5 N but not at 1 N
   (the accepted-but-not-at-1-N count is printed, not gated).
2. **graph and edge tables**: `edge_tables(T)` gives, per edge and per
   rung, `n_solved`, `n_accepted`, `closed_frac` (= n_accepted / 36),
   `tf_min_days` and `tf_med_days` over accepted cells (NaN when none),
   and the lowest accepted rung `F_low_N` with `F_low_str` ("0.5 (lowest
   tested)" when it is the last rung). `network_graph(E, N)` builds the
   `digraph` with node variables (family, period, C_J, sigma, stable) and
   edge variables (sheet_index and the per-rung tables); every attribute
   used later is read from `G.Edges`. Print the edge table at 1 N and
   0.5 N.
3. **accessibility plane**: caps span the smallest to the largest
   accepted cell time over all rungs, extended to cover the assembled
   relay times from section 5. `accessibility_plane(T, caps)` returns, at
   each (rung, cap), `pairs` and `cells` (accepted cells at that rung with
   t_f under the cap: the raw CATALOG COVERAGE at that thrust bound) and
   `pairs_env`, `cells_env` (the known-feasible ENVELOPE: a trajectory
   accepted at a lower thrust bound is feasible under a higher one, so
   the envelope at rung r unions every rung at or below r's thrust). Gate
   G3: all four arrays non-decreasing along the cap axis, and the envelope
   non-decreasing from the deepest rung upward, with zero violations. The
   rung axis of the raw arrays is not gated (section 2).
4. **strengths**: in-strength and out-strength with weights 1 / tf_min
   at the reference rung, read from `G.Edges`; one table, with the stated
   caveat that between-family rankings reflect the catalog's topology
   (each tulip has 12 in-edges, each L2 halo 2 in and 4 out). No
   betweenness: with phase-compatible route times a per-edge additive cost
   does not exist, and the only candidates would be the L2 halos anyway.
5. **relay closure**: `relay_closure(T, N, src, via, dst)` enumerates
   the 32 routes L1 halo -> L2 halo -> tulip and, per rung, computes the
   phase-compatible assembled time
   t_route = min over accepted cells (a, b) of leg 1 and (c, d) of leg 2
   of t1(a, b) + T_via * mod(sD2(c) - sA1(b), 1) + t2(c, d), the wait of
   that minimizer, and the optimistic lower bound t_sum = min t1 + min t2.
   Per (L1 halo, tulip) pair: the best relay by t_route at each rung, the
   deepest rung with any closed route. The overlay for the plane is
   `reach` = pairs + the number of (src, dst) pairs with t_route under the
   cap (directed reachable ordered pairs beyond the direct edges), and
   the relay-usage count per L2 halo (how many pairs it serves as the best
   relay at the reference rung). Gate G4: exactly 32 routes over 8 pairs,
   each route's two sheets share the relay key, each leg's phase grid on
   the relay is the same set (sA_frac of leg 1 equals sD_frac of leg 2),
   `closed` equals the AND of the legs' accepted masks, t_route and t_sum
   NaN exactly when open, and t_route >= t_sum wherever closed. The note
   states the mass composition (a leg solved at 150 kg flown at the
   arrival mass of the previous leg by scaling thrust and mass by q =
   m1/m0 <= 1 keeps the same path and satisfies the thrust bound).
6. **hypothetical reverse graph**: G5 = for every node, periodic closure
   |phi(x0, T) - x0| < 1e-8 AND the mirror image of the off-plane state
   phi(x0, 0.37 T) lies on the orbit (`symmetry_residual` < 1e-8). G5 is
   printed always and throws only when `hypothetical` is on. When on, the
   56 reversed edges are added with `hypothetical = true`, `n_accepted =
   0`, the same `tf_*` numbers carried as geometry only; the script
   reports the number of strongly connected components and the number of
   ordered node pairs reachable by a directed path, every line prefixed
   HYPOTHETICAL. Those numbers never enter sections 3-5 or the figures.
7. **figures**: F1 the graph at 1 N and at 0.5 N, built from the
   edge-table subset of accepted edges (width = closed_frac, color =
   tf_min_days, labels with Interpreter none); F2 the raw coverage plane
   (pairs) with two boolean-mask contours, all 56 direct pairs under the
   cap and at least one relay route under the cap, and the envelope plane
   beside it; F3 catalog closure down the ladder (mean closed_frac per
   rung by departure family, the 15 N dip labelled as ladder coverage);
   F4 relay closure, per (L1 halo, tulip) pair the assembled t_route (solid)
   and the optimistic t_sum (dotted) against thrust. All PNG into
   `results/`.

Summary block: G1-G5 PASS/FAIL, throws on a FAIL of G1-G4 (G5 only when
`hypothetical`); `results/network_summary.txt` carries every printed
fact, the strength table, the relay table and the usage counts, so the
note cites one source.

## 8. `lib/`

- `mirror_state(x)`: 6 x 1 or 6 x n; returns (x, -y, z, -vx, vy, -vz).
  Refusal `:shape`.
- `symmetry_residual(x1, T, mu, opts)`: min over t in [0, T] of
  max|phi(x1, t) - mirror_state(x1)| by a coarse grid and `fminbnd`;
  returns the residual and the minimizing time. Meaningful only for an
  off-plane x1 (the caller passes phi(x0, 0.37 T)).
- `edge_tables(T)`: as in section 7.2.
- `network_graph(E, N)`: the digraph; node names are the keys; edge
  variable `sheet_index` maps each edge back to `T`/`E`.
- `accessibility_plane(T, caps_days)`: as in section 7.3; also returns the
  per-(rung, cap) list of sheet indices counted in `pairs`.
- `relay_closure(T, N, src_prefix, relay_prefix, dst_prefix)`: as in
  section 7.5.
- `partial_f_test(y, X0, X1)`: nested OLS; refuses a non-nested or
  rank-deficient design (`:nested`, `:rank`); F and upper-tail p through
  betainc; per-coefficient se, t, two-sided p and 95% CI through the
  inverse incomplete beta; RSS0 = RSS1 = 0 gives p = 1; a negative
  RSS0 - RSS1 within 1e-12 of RSS0 is clamped, a larger one refused.

## 9. Tests (`tests/`)

- `test_load_catalog_sheets`: 56 rows; catalog order; v1 rows carry
  dep_family 'dro' and a tau in dep_params; every array 6 x 6 x 11; keys
  unique per (dep, arr); `accepted` implies `has_solution`; accepted
  differs from solved only at 0.75 and 0.5 N; a -1 in the raw conj_pass is
  not a pass; refusal on a bad root.
- `test_orbit_node_table`: on a 2-row synthetic T: 3 nodes, periods equal
  T_dep/T_arr, symplectic pairing to 1e-6, DRO stable with sigma exactly
  0, cache round-trip, a changed muStar misses the cache, `rebuild`
  forces.
- `test_lib`: `partial_f_test` (hand case, exact fit, both-exact p = 1,
  CI contains the coefficient, `:nested`, `:rank`, `:shape`);
  `mirror_state` (hand vector, involution, `:shape`);
  `symmetry_residual` (DRO off-plane point < 1e-8; a non-periodic state
  > 1e-6 so the gate can fail); `edge_tables` (synthetic sheets: counts
  over accepted, NaN on an empty rung, the fastest solved-but-refuted
  cell is NOT selected, lowest-rung strings); `network_graph` (scrambled
  input order: every attribute found by (dep, arr) is the right one);
  `accessibility_plane` (counts, envelope includes a deeper-rung witness
  under a shallower cap, caps below the data give zeros, monotone);
  `relay_closure` (toy legs: closed iff both, a case where the smallest
  sum is not the smallest assembled time and the assembled one is chosen,
  mismatched phase grids refused).
- `test_pipeline`: runs both scripts headless (figures invisible, results
  into a temp folder via the base-workspace `study_opts`), then the
  network script again with `study_opts.hypothetical = true`; asserts the
  summary blocks, every gate PASS, the H verdict lines present, the
  HYPOTHETICAL prefix and the reversed-edge count in the second run, every
  expected PNG written, and a parse check on both scripts. Under 3 minutes
  with a warm node cache.
- `run_tests`: the four, in that order, one logical verdict.

## 10. The note (`doc/orbit_network_note.tex`, about 6 pages)

Sections: (1) motivation, Braik and Ross in one paragraph; (2) their
problem against ours: single Jacobi level, energy-preserving impulses and
a symmetric half-budget split, versus a mass-varying min-time OCP with
certified extremals; the lemma, stated as set inclusion: if F_AB is the
feasible set of the A -> B problem with all its state, mass, boundary and
path constraints and R_ACB is the subset of trajectories that pass
through C with compatible full states at the joins and no relay-only
operation, then inf over F_AB of J <= inf over R_ACB of J for the same
cost functional, a statement about global infima, not about what a
particular solver returned; a relay witness below the termination thrust
of the computed direct branch would be another feasible unrestricted
transfer and diagnoses the branch or the search, not a global barrier;
(3) data: the four catalogs, the sheet table, the accepted mask, the
coverage facts; (4) the graph, the lowest accepted rung, the strength
table with its topology caveat; (5) the plane: raw coverage against the
known-feasible envelope, the diagnostics (direct availability, directed
reachability with relays) named as diagnostics, not regimes; (6) relay
closure with the phase-compatible assembled times, the mass composition
remark, and the statement that relay-versus-direct cost is untestable
here; (7) stability versus cost in the review's wording: not detected at
the nominal threshold in this exploratory sheet-level model, the joint
time term carried by arrival instability on the halo-to-halo sheets,
shared nodes, nominal p-values, no claim about 70 mN, nonsignificance is
not equivalence; (8) the hypothetical reverse graph in one paragraph
(what it is and is not); (9) what the legs campaign would add, with the
present data supporting candidate prioritization among the L2 halos, not
identification of the best DRO relay. Figures from `results/`. Compile
with pdflatex, clean aux files, `\author{Michael Casey}`.

## 11. Out of scope

New solves of any kind; a query API; min-fuel edges; the 70 mN DRO phase
library (a single cell, not a sheet); physical reverse transfers; any
claim about relay-versus-direct cost, global thrust thresholds, or the
absence of a stability-transport relationship.
