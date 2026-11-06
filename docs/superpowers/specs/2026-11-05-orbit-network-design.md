# Orbit network over the costate catalogs: design

Date: 2026-11-05. Status: approved in conversation (Mike, 2026-11-05); the
plan argues from this document.

## 1. Purpose

Braik and Ross (arXiv 2605.31543, "Orbital Networks in the Three-Body
Problem") organize periodic orbits of the planar Earth-Moon CR3BP into a
weighted graph whose edges come from reachable-set overlap under
heading-change impulses on one Jacobi level, then read hubs, relays and
budget regimes off the graph. Their maneuver model is not ours and their
headline (multileg routes cheaper than direct) is an artifact of it: in a
genuine optimal-control formulation the direct feasible set contains every
relay route, so direct never loses on cost. What survives for us is the
framing. This folder lifts two of their ideas onto the certified min-time
catalogs we already hold, using NO proxies and NO new solves:

1. **Stability versus cost** (their section 6 claim: stable orbits are the
   hard nodes, unstable ones share a chaotic sea). Tested at SHEET level
   across the four catalogs.
2. **The low-thrust orbital network** (their sections 5-6). Nodes are
   pumpkyn orbits, edges are certified catalog sheets, and the edge
   attribute that matters for an OCP is feasibility across thrust: a thrust
   floor and a phase-resolved t_f distribution. Their (dV cap, time cap)
   budget plane becomes our (thrust, time cap) accessibility plane.

Decisions already made (do not reopen):

- Sheet level, not family level: every (departure orbit, arrival orbit)
  sheet is one point, 56 in all.
- Data-only: the 56 existing edges, thrust 15 N down to 0.5 N. The legs that
  would let a relay compete with the direct DRO-to-tulip wall (DRO-to-DPO,
  DRO-to-L2 halo, and their deep rungs) are a LATER campaign; this folder
  tells us which relay orbit is worth that ladder time.
- Deliverable form: a study script per topic plus one LaTeX note. No query
  tool/API until the missing legs exist.
- Relay-versus-direct cost is NOT testable on this data (no pair has both a
  direct edge and a two-leg route). The note carries the optimal-control
  argument; the script tests relay CLOSURE only.
- Mirrored (reverse) edges by the x-axis symmetry are an option, off by
  default, labelled approximate (mass flow breaks the exactness: the
  reversed trajectory runs m_f -> m_0). Headline numbers never use them.
- sigma is computed from `pumpkyn.cr3bp.monodromy` directly as
  max(0, ln|lambda|_max / T). `pumpkynPie.cr3bp.orbitProperties` errors on
  stable orbits (empty unstable-eigenvalue list) and is not used.

## 2. What the data is

Four min-time catalogs (schema v1/v2), all with `rungs_N = [15 12 10 7 5 3
2 1.5 1 0.75 0.5]` and 6 x 6 phasing grids per sheet:

| catalog | file | sheets | departure orbits | arrivals |
|---|---|---|---|---|
| DRO_tulip | `DRO_tulip/direct/results/costate_catalog_dro_tulip.mat` | 16 | DRO tau 0.5, 1, 2, 3 | tulip Np 5, 7, 9, 12 (pm -1) |
| HALO_tulip | `HALO_tulip/direct/results/costate_catalog_halo_tulip.mat` | 16 | L2 south halo tau 1.75, 2.2, 2.8, 3.4 | same tulips |
| DPO_tulip | `DPO_tulip/direct/results/costate_catalog_dpo_tulip.mat` | 16 | DPO tau 1, 2, 3, 4 | same tulips |
| HALO_HALO | `HALO_HALO/direct/results/costate_catalog_halo_halo.mat` | 8 | L1 south halo tau 1.8037, 2.7433 | the four L2 halos |

Sheet fields used: `has_solution`, `tf_nd`, `conj_pass` (each 6 x 6 x 11),
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
- Every sheet has at least 9 cells closed at 1 N (range 9-34 of 36).
- Node stability (|lambda|_max of the monodromy matrix, sigma in 1/ND):
  all four DROs, all four tulips, DPO tau 1 and L2 halo tau 2.2 are
  linearly stable (|lambda|_max within 3e-4 of 1); L2 halos tau 1.75 / 2.8
  / 3.4 have sigma 0.63 / 1.18 / 2.05; L1 halos tau 1.8037 / 2.7433 have
  sigma 0.88 / 2.83; DPOs tau 2 / 3 / 4 have sigma 1.22 / 1.71 / 1.67.
  So sigma_arr is zero on 48 of 56 sheets; the live contrast is departure
  stability: 24 sheets leave a stable orbit, 32 an unstable one.

## 3. Layout

```
orbit_transfer/orbit_network/
  README.md                     what the folder is, how to run, what it does not claim
  load_catalog_sheets.m         the four catalogs -> one 56-row sheet table
  orbit_node_table.m            the 18 orbits with period, C_J, |lambda|max, sigma, x0
  stability_cost_study.m        topic 1 study script (sections 0-5, gates, figure)
  orbit_network_study.m         topic 2 study script (sections 0-7, gates, figures)
  lib/
    mirror_state.m              x-axis symmetry R(x) = (x, -y, z, -vx, vy, -vz)
    accessibility_plane.m       (thrust rung, time cap) -> pair and cell counts
    partial_f_test.m            nested linear models -> F, p (betainc; no toolbox)
  tests/
    test_load_catalog_sheets.m
    test_orbit_node_table.m
    test_lib.m                  mirror_state, accessibility_plane, partial_f_test
    test_pipeline.m             both scripts headless into a temp folder; parse check
    run_tests.m                 the folder's runner, one logical verdict
  results/                      git-ignored: figures, the node-table cache
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

`root` is the `orbit_transfer` folder (default: resolved from this file).
Returns a 56 x 1 struct array, one row per sheet, catalog order DRO_tulip,
HALO_tulip, DPO_tulip, HALO_HALO, sheet order as stored. Fields:

| field | type | meaning |
|---|---|---|
| catalog | char | 'dro_tulip' \| 'halo_tulip' \| 'dpo_tulip' \| 'halo_halo' |
| sheet_index | double | position in that catalog's `.sheets` |
| dep_family, arr_family | char | 'dro' \| 'dpo' \| 'halo' \| 'tulip' |
| dep_params, arr_params | struct | exactly the `get_family_orbit` parameter struct |
| dep_key, arr_key | char | canonical node key, e.g. `dro_tau1`, `halo_L2_tau2.2`, `tulip_Np7`, `halo_L1_tau1.8037` |
| T_dep, T_arr | double | periods (ND) as the catalog stores them |
| rungs_N | 1 x 11 | copied from the catalog |
| has_solution, tf_nd, conj_pass | 6 x 6 x 11 | copied from the sheet |
| sD_frac, sA_frac | 1 x 6 | copied |
| constants | struct | muStar, lStar_km, tStar_s from the catalog |

The v1 DRO catalog is mapped onto the v2 field names here (dep_family
'dro', dep_params struct('tau', tauDRO), arr_params struct('Np', Np, 'pm',
pm), T_dep = tauDRO, T_arr = period_tulip_nd), so no consumer ever branches
on schema version. Refusals: a missing catalog file (`:missingCatalog`,
names the path), rung vectors that differ between catalogs (`:rungs`), and
a sheet whose arrays are not 6 x 6 x numel(rungs) (`:shape`).

## 5. `orbit_node_table(T, opts)` -> `N`

The distinct `dep_key`/`arr_key` values of `T`, 18 of them, each built once
with `get_family_orbit(family, params)` (its first output is the time grid,
second the 325 x 6 state history, third `info` with `periodND`; the initial
state is row 1). Fields per node: `key`, `family`, `params`, `period_nd`
(= `info.periodND`), `x0` (6 x 1), `CJ` (`pumpkyn.cr3bp.jacobi`),
`lam_max`, `lam_min` (largest and smallest |eig| of
`pumpkyn.cr3bp.monodromy(x0.', period_nd, muStar)`), `sigma` = max(0,
log(lam_max)/period_nd), `stable` = lam_max < 1 + 1e-3, `in_degree`,
`out_degree` (counts over `T`). The loader's continuation prints a log;
wrap in `evalc`.

`opts.cache` (default `results/orbit_nodes.mat` beside this file): when the
file exists and its keys equal the requested keys it is loaded; otherwise
the table is rebuilt and saved. `opts.rebuild = true` forces a rebuild.
Gate G1 (below) compares `period_nd` with the catalog's `T_dep`/`T_arr`
for every row of `T`, to 1e-9, so a stale cache cannot pass silently.

## 6. `stability_cost_study.m` (topic 1)

Sections:

0. **settings**: paths, root, the two rungs (1 N as the shared cost rung,
   0.5 N as the wall rung), the significance level 0.05, figure folder.
1. **load and census**: `T`, `N`; print the 56-row table (catalog, dep,
   arr, cells closed at 1 N and 0.5 N, min t_f at 1 N in days).
   Gates: P1 numel(T) == 56; P2 identical rungs across catalogs (the loader
   already refuses, the gate re-states it); P3 every cell closed at 0.5 N
   is closed at 1 N (so surv is a ratio of nested sets), with the
   whole-ladder non-monotone count printed as a fact; P4 every sheet has
   at least one cell closed at 1 N.
2. **orbit properties**: for each node print period, C_J, |lambda|max,
   sigma, stable. Gates: P5 symplectic pairing lam_max*lam_min = 1 to 1e-6
   and |det M - 1| < 1e-6 for every node; P6 node periods match the
   catalogs to 1e-9 (G1); P7 the arrival nodes of every tulip sheet are
   stable (the fact section 2 states), reported as a fact check.
3. **covariates and responses** per sheet: dCJ = CJ_arr - CJ_dep;
   sigma_dep, sigma_arr; tf1 = min over cells closed at 1 N of tf_nd,
   converted to days with tStar_s; surv = (cells closed at 0.5 N) / (cells
   closed at 1 N).
4. **regression**: for each response y in {tf1, surv}: model M0 y ~ 1 +
   dCJ; model M1 y ~ 1 + dCJ + sigma_dep + sigma_arr; least squares by
   backslash; partial F for M1 over M0 and its p-value from
   `partial_f_test`, which also returns per-coefficient standard errors,
   t statistics and two-sided p-values for the larger model; group means
   of y for stable-departure versus unstable-departure sheets. Verdicts,
   printed as SUPPORTED / NOT SUPPORTED (these are hypotheses, not code
   gates), rest on the sigma_dep coefficient's OWN t-test, because the
   paper's claim is about departure stability beyond energy and sigma_arr
   varies only on the 8 halo-to-halo sheets (confounded with catalog):
   H1 the sigma_dep coefficient for tf1 has two-sided p < 0.05 and the
   paper's sign (negative, unstable departures are faster); H2 the same
   for surv (sign positive, unstable departures survive the thrust wall
   better). The joint F p-value and the sigma_arr test are printed beside
   each verdict. Both outcomes are results; the script must not treat NOT
   SUPPORTED as a failure.
5. **figure**: one PNG, two panels, each response against dCJ, markers
   colored by sigma_dep (stable group in one color), with the M0 fit line;
   written to `results/stability_cost.png`.

The summary block at the end lists P1-P7 PASS/FAIL and H1-H2 verdicts and
throws if any P gate fails (so the pipeline test can assert on it).

## 7. `orbit_network_study.m` (topic 2)

Sections:

0. **settings**: paths, rungs, the time-cap grid (25 caps from the
   smallest closed t_f to the largest, in days), `useMirror = false`,
   figure folder.
1. **load**: `T`, `N`; gates G1 (periods, 1e-9) and G2 (56 edges, zero
   cells closed at 0.5 N but open at 1 N; the whole-ladder non-monotone
   count printed as a fact, see section 2).
2. **graph**: a `digraph` with 18 nodes (keys as names, node table columns
   as node variables) and 56 edges; per edge and per rung the edge tables
   `closed_frac` (of 36), `tf_min_days`, `tf_med_days` (over closed cells),
   `n_certified` (has_solution & conj_pass), and the thrust floor
   `F_floor_N` = the smallest rung with closed_frac > 0, reported as a
   bracket (`<= 0.5` when the 0.5 N rung closes). Print the edge table at
   1 N and 0.5 N.
3. **accessibility plane**: `accessibility_plane(T, rungs, caps)` returns
   two 11 x 25 integer arrays: pairs whose tf_min <= cap at that rung, and
   cells whose tf <= cap. Regimes overlaid: direct access (number of pairs
   closed), weak connectivity of the 18-node graph at that (rung, cap), and
   multileg closure (number of two-leg L1-halo-to-tulip routes closed, from
   section 5). Gate G3: both arrays non-decreasing along the cap axis with
   ZERO violations (a violation here is a coding error, not a finding).
   Along the rung axis the counts are NOT gated: the 15 N starting rung
   closes fewer cells than 12 N (section 2), and the plane shows that as a
   dip at the top, which the note explains as ladder coverage.
4. **centrality**: in-strength and out-strength with weights 1/tf_min_days
   at 1 N, and betweenness on the directed graph; one table. The text
   states that only the L2 halos can have nonzero betweenness on this
   graph.
5. **relay closure**: the 32 two-leg routes L1 halo -> L2 halo -> tulip.
   Per route and rung: closed iff both legs have a closed cell; t_total =
   tf_min(leg 1) + tf_min(leg 2); t_wait = t_total + period of the L2 halo
   (the phase hand-over bound: leg 2's min is over departure phases, so the
   wait at the relay is at most one relay period). Report the deepest rung
   at which each L1 halo reaches each tulip, and the best relay for each.
   Gate G4: a route is closed at a rung iff both of its edges are (checked
   against the edge table, not recomputed).
6. **mirror option**: for every node, `mirror_state(x0)` must lie on the
   same orbit: minimize |phi(x0, t) - R x0| over t in [0, T] (coarse grid
   then `fminbnd`); gate G5: residual < 1e-8 for all 18 nodes. If
   `useMirror`, add the 56 reversed edges with the same per-rung tables and
   an `approximate = true` flag, and report how many node pairs become
   mutually reachable and how many two-leg routes appear; this section's
   numbers are labelled APPROXIMATE in the printout and never enter
   sections 3-5.
7. **figures**: F1 the graph at 1 N and at 0.5 N (layered layout, edge
   width = closed_frac, edge color = tf_min_days, mirrored edges dashed if
   on); F2 the accessibility plane (pairs) with the connectivity and
   multileg-closure contours; F3 per-rung edge survival (closed_frac
   averaged) by departure family; F4 relay t_total and t_wait against
   thrust for the best route per (L1 halo, tulip). All PNG into `results/`.

Summary block: G1-G5 PASS/FAIL, throws on any FAIL.

## 8. `lib/`

- `mirror_state(x)`: 6 x 1 or 6 x n; returns (x, -y, z, -vx, vy, -vz).
  Demo: a DRO initial state. Refusal `:shape`.
- `accessibility_plane(T, rungs, caps)`: as in section 7.3; also returns
  the per-(rung, cap) list of closed pair indices so the connectivity
  overlay reuses it. Demo: a synthetic 2-sheet table.
- `partial_f_test(y, X0, X1)`: nested OLS; returns F, p, b0, b1, RSS0,
  RSS1, df. p = betainc(df2/(df2 + df1 F), df2/2, df1/2) (upper tail).
  Refusal `:nested` if X0's columns are not a subset of X1's. Demo: a
  known case with p printed against the textbook value.

## 9. Tests (`tests/`)

- `test_load_catalog_sheets`: 56 rows; catalog order; v1 rows carry
  dep_family 'dro' and a tau in dep_params; every array 6 x 6 x 11; keys
  unique per (dep, arr); refusal on a bad root.
- `test_orbit_node_table`: on a 2-row synthetic T (one DRO sheet, one halo
  sheet) with a temp cache: 3 nodes, periods equal T_dep/T_arr, symplectic
  pairing to 1e-6, DRO stable, cache round-trip equal, `rebuild` forces.
- `test_lib`: mirror_state on a hand vector and the refusal;
  accessibility_plane on a 2-sheet synthetic table with known counts and
  the monotone property; partial_f_test against a hand-computed 6-point
  case and the `:nested` refusal.
- `test_pipeline`: runs both scripts headless (figures invisible, results
  into a temp folder via a variable the scripts honour when it exists in
  the base workspace, as `test_sleigh_pipeline` does), asserts the
  summary blocks exist in the captured output with no FAIL, and that
  every expected PNG was written; plus a parse check (`checkcode`) on both
  scripts. Runtime budget: under 3 minutes with a warm node cache.
- `run_tests`: the four, in that order, one logical verdict, registered in
  the folder README (the orbit_transfer tree has no single runner; this
  folder's runner is its own).

## 10. The note (`doc/orbit_network_note.tex`, about 6 pages)

Sections: (1) motivation, Braik and Ross in one paragraph; (2) their
problem against ours: single Jacobi level, energy-preserving impulses and
a symmetric half-budget split, versus a mass-varying min-time OCP with
certified extremals; the direct-never-loses argument stated as a lemma
(feasible-set inclusion); (3) data: the four catalogs, the sheet table;
(4) the graph and the thrust-floor edge attribute; (5) the accessibility
plane and the three regimes; (6) relay closure, L1 halo to tulip through
L2 halos, and what cannot be tested yet; (7) stability versus cost, the
sheet-level regression and its verdict; (8) what the legs campaign
(DRO-to-DPO, DRO-to-L2 halo, deep rungs toward 70 mN) would add. Figures
from `results/`. Compile with pdflatex, clean aux files, `\author{Michael
Casey}`.

## 11. Out of scope

New solves of any kind; a query API; min-fuel edges; the 70 mN DRO phase
library (a single cell, not a sheet); reverse-transfer t_f beyond the
mirrored approximation; any claim about relay-versus-direct cost from
data.
