# root_origins_study.m -- operational notes

Pedagogy lives in `doc/root_origins_study_guide.tex`. This file is the
operator's view: what the script does, what it writes, how to rerun it, and
what has bitten before. Script: `DRO_tulip/indirect/root_origins_study.m`.
Record of runs: FINDINGS 75-77.

## What it answers

Where does the FIRST root of a min-time DRO->tulip transfer come from when no
neighbouring solution exists to seed from? Three mechanisms, each measured:
(A) a cold direct solve at the operating point across mesh densities (the
"lottery"); (B) a ladder down the engine axis -- direct collocation 15 N ->
0.5 N, a covector harvest, then multiple shooting 0.5 N -> the operating
point; (C) a branch-blind basin hunt at other arrival phases from the root
reached.

## Sections and what each produces

| # | section | writes to `S.` | gate / outcome |
|---|---|---|---|
| 0 | user inputs (problem, `op`, `ref` WITH its problem, `cold`, `dladder`, `iladder`, `hunt`, `num`, `tol`) | `config` | -- |
| 1 | orbits from parameters, endpoints, `refApplies` | -- | assert closure/seam |
| 2 | cold lottery: one cold solve per mesh, `directOK` each | `lottery` (per mesh) | outcome A |
| 3 | direct ladder, rung by rung, band screen, throttle check | `direct` (every rung, refused too) | B1 (fatal), B1b |
| 4 | Sundman re-solve of the bottom rung, harvest, shoot | `handoff` | H1 (fatal), H2, H3 (3-state) |
| 5 | indirect ladder from banked starts, exponent sweep, per-attempt records | `indirect` (every rung, `attempts` per exponent), `reached` | B2 |
| 6 | certify_root with forwarded options; adopt polished root WITH rescaled grid | `certificate`, `certOpts` (no pool), `reached` (polished) | C1; outcome B decided here |
| 7 | basin hunt from `flyFromJunctions` warm start; candidates banked as seeds | `hunt` (per probe: solverOK / accepted / harvested / seed / direct) | outcome C |
| 8 | self-check | `gateStatus`, `outcomes`, `h3applicable` | asserts B1 && H1 only |

Outcomes are SUPPORTED / PARTIAL / REACHED NUMERICALLY / CANDIDATES /
NOT OBSERVED / INCONCLUSIVE. A stopped ladder, a stalled walk and an
uncertified root are FINDINGS, not failures; only B1 and H1 throw.

## Running it

- **MATLAB R2026a** (`/Applications/MATLAB_R2026a.app/bin/matlab`): it has the
  Parallel Computing Toolbox, so the fences are real. R2025b has no PCT
  licence here; the script then runs UNFENCED and, unless
  `num.allowUnfenced = true`, records "not certified: no pool" at section 6.
- Launch as a batch job under an OS watchdog, log to a file (matlab -batch
  buffers stdout). Pattern used: `scratchpad/run_origins.sh` (4 h kill).
- Needs pumpkynPie on the path (`startup()` in
  `~/Desktop/proj7/external/pumpkynPie`) and `~/casadi-3.7.0`.
- Default run: ~55 min. Lottery N=1600 alone is ~3 min; the 0.09 N indirect
  rung can burn its full 600 s budget. Reduced settings for a smoke:
  `meshes [200 400]`, three direct rungs, two indirect rungs, one hunt
  phase, `N 400` -> ~6 min.
- Record: `indirect/results/root_origins_study.mat` (results/ is gitignored),
  saved after every rung/probe (`saveq` warns on failure).

## What the default run shows (2026-09-16)

Lottery 4.6809 / 4.9909 (iterate) / 6.4126 ND vs 4.0152 reference. Direct
11/11, worst throttle slack 1.3e-6. Handoff lambda_t = 1.000000, |R| 4.6e-12
in one iteration. Indirect 10/12: stalls at 0.11 N with every exponent
`unconverged` (a search-policy stall on THIS cell; the campaign reached
0.09 N from the sheet's fastest 0.5 N cell). Certified at 0.12 N / 900 s,
10.6060 d. Hunt: 2 of 3 accepted, one 14.4% faster at sA 0.5754, both banked.

## Knobs that matter

- `dladder.rungs`: ratios ~0.75. A 1 -> 0.5 N step jumped to a 54.7 d / 46.8
  km/s solution that passed every feasibility check; the band screen
  (`tfBand`, steps only) is what catches that.
- `iladder.rungs`: put the Isp stage where the walk still converges; finer
  ratios bought 0.16 -> 0.12 N; the Isp stage did not move the 0.11 N stall.
- `iladder.tfExps`: the fastest accepted candidate across exponents is kept
  -- a search, not branch-preserving continuation.
- `iladder.maxPropFrac` (0.6): a seed-policy cap, not exhaustion (which is
  at c/T).
- `tol.thrMin` (1e-3): the throttle is FREE by design; barrier slack is
  1e-7..1e-6; a real switch goes to 0. Checked on both sides.
- `hunt.cpuSec` (300): a capped probe is a capped probe; run to run it lands
  on either side (sA 0.7837 converged in 61 s once, capped at 308 s once).

## Contracts it depends on (verified 2026-09-16)

- `ms_tfmin`/`ms_bvp` return `info.Y` as 14 x K junction STARTS with a 1 x
  (K+1) grid; consumers read columns 1..K only. Feeding `info.Y` back is
  lossless. The library HEADERS still say K+1 -- open item.
- `certify_root` reads `m0kg`, `gateKm`, `gateVms`, `moonKmMin`,
  `allowUnfenced` from ITS options; returns `C.z`, `C.Y` but NO grid: after
  a polish, `tGrid = sig*C.z(8)`.
- `run_capped` `ok = false` means timeout OR worker error.
- `casadi_mintime_dro` needs `returnModel = true` for multipliers; `thrMin`,
  `maxInterp`, `tfSpread` are separate from `maxDefect`; `altMinKm` is an
  ALTITUDE over NODES.
- `flyFromJunctions` (local): banked start wins each seam; returns the
  worst seam mismatch (2.7e-8 km on the default run).

## Things that bit

1. Comment-only lines inside a `...`-continued struct literal are a parse
   error (twice).
2. The direct ladder's first rung has a cold guess; the band must not apply
   to it (41x off at 15 N).
3. A struct array declared without a field later assigned -> "dissimilar
   structures" at the first save.
4. `gcp('nocreate')` in a default expression is evaluated eagerly and throws
   without PCT (fixed library-wide via `current_pool`).
5. "dV rises down the ladder" -- it falls; T*t_f ~ T^0.4. Read the table
   before writing the sentence under it.
