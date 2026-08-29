# costate_catalog_gto_tulip — GTO → Tulip Minimum-Time Costate Catalog

Every entry in this catalog is a **converged minimum-time PMP solution**, not a
hint: `pumpkyn.cr3bp.tfMin` accepts each one **unchanged** (|Δz| < 1e-6). Give
the picker five coordinates — GTO departure orientation, tulip petal count,
departure phase, arrival phase, thrust — and you get the flight time and the
full costate/final-time vector `z8 = [λr(3); λv(3); λm; tf]` ready to fly or
to seed `tfMin` at your own endpoints.

This is the fifth costate-catalog case (after DRO/HALO/DPO/L1↔L2-halo → tulip)
and the first with a **non-periodic departure orbit** — a GTO. There is no
orbital period to key sheets on, so the sheet axis is a departure
**orientation angle** instead (see Schema notes below). It ships at the same
propulsion standard as its siblings (150 kg / 1710 s) walked down to 5 N, not
1 N — see Scope below for why.

## Contents

| File | Role |
|---|---|
| `results/costate_catalog_gto_tulip.mat` | The catalog: 16 sheets, 2,625 entries (schema v2) |
| `costate_catalog_pick.m` | **This catalog's own** five-coordinate lookup (NOT the shared/generic picker — see Schema notes) |
| `build_gto_catalog.m` | Packager: sheet .mats → the single schema-v2 catalog variable |
| `run_gto_catalog.m` | THE SWATH DRIVER — solves/resumes/regenerates any subset of the 16 sheets |
| `gto_entry.m` | THE SINGLE-ENTRY DRIVER — one (orientation, Np, phase, phase, thrust) cell through the same engine, with `nargin==0` self-demo |
| `audit_gto_entries.m` | 3-entry spot-audit replay (fly z8 → tfMin round-trip, sheet-local endpoint reconstruction) |
| `smoke_gto_cell.m` | 1×1 cell smoke test through the ladder engine (development artifact, kept for reference) |
| `setup_paths.m` | This module's path setup (costate_common, DRO_tulip/indirect engine, pumpkynPie) |
| `results/catalog/gto_o*_Np*.mat` | The 16 engine-native sheet files (`OK/TF/Z8/ATT/rungs/sD/sA/meta`) the packager reads |

**The 16 per-sheet `.mat` files are gitignored (only their `_progress.txt` companions are tracked); the packaged catalog `.mat` IS committed (force-added), so a fresh clone carries the catalog itself but not the raw sheets**
(`*.mat` in the repo `.gitignore`) — a fresh clone of this repo does not
carry them. To reproduce them, re-run the fleet (`run_gto_catalog()`, ~40 h
wall per the Provenance notes below) or obtain the data files directly
from whoever ran the campaign.

## Dependencies

pumpkyn + pumpkynPie on the MATLAB path (`setup_paths` handles this). The
packaged `.mat` itself is dependency-free and self-describing; only
`costate_catalog_pick.m` needs to be on the path to look it up (it has no
external dependency either — the true-Kepler-period formula it needs is
inlined, per deliverable-picker convention).

## Quick start

**Look something up in the built catalog:**

```matlab
L = load('results/costate_catalog_gto_tulip.mat');
cat_ = L.costate_catalog_gto_tulip;
[tf_nd, z8, info] = costate_catalog_pick(cat_, 90, 7, 0.0, 0.0, 10.0);
% tauDep = 90  -> orientDeg = 90 (NOT a period -- see Schema notes)
% arrKey = 7   -> Np = 7 tulip petals
% depPhaseDays/arrPhaseDays -> snapped to the nearest grid phase
% thrustN = 10 -> interpolated tf between the bracketing rungs, z8 from
%                 the nearest stored rung
```

**Solve or refresh one entry directly** (the single-entry driver — same
engine, same gates, same warm recipe as the fleet, never a bespoke solver):

```matlab
setup_paths
[entry, gates] = gto_entry(0, 7, 0.0, 0.0, 15);   % orient 0 deg, Np 7, perigee, 15 N
% entry.z8, entry.tf_nd, entry.dV_kms, entry.mf_kg
% gates.msNormR, gates.flownKm, gates.acceptDz

gto_entry(270, 5, 0.10, 0.30, 5, struct('writeSheet', true, 'conjTest', true));
% banks this entry (and every stepping-stone rung solved on the way down
% to 5 N) into results/catalog/gto_o270_Np5.mat, and runs the conjugate
% test on the accepted z8
```

**Regenerate or extend the fleet** (resumable — revisits only no-OK cells):

```matlab
setup_paths
run_gto_catalog();              % all 16 sheets, default rungs [15 12 10 7 5]
```

**Extending a shipped sheet's rung ladder is a deliberate act, not a
drop-in override.** `run_gto_catalog`'s rung-mismatch guard refuses any
call whose `rungsIn` does not match a sheet's on-disk `Q.rungs` exactly —
this is what catches an accidental silent re-solve-from-scratch that would
wipe the sheet's accumulated entries (task-4 fix, 2026-08-27). So passing
a longer ladder (e.g. `[15 12 10 7 5 3]`) against a shipped 5-rung sheet
**errors**, by design — it does not append the new rung. To extend a rung
deliberately, move the sheet aside first so the call starts a fresh sheet
at the new rung count:

```matlab
setup_paths
movefile('results/catalog/gto_o000_Np5.mat', ...
         'results/catalog/gto_o000_Np5_5rung.bak.mat');   % move aside
run_gto_catalog(1, [], [], [15 12 10 7 5 3]);              % sheet 1 only,
                                                            % fresh 6-rung ladder
```
(`run_gto_batched.sh` wraps this under `nohup` with per-pass OS-kill and
resume for unattended multi-hour/day runs — see the matlab-campaign skill.)

## Orientation-axis convention (for Darin)

`orientDeg` is **the angle from the Earth→Moon line to the ellipse's perigee
direction, measured in the rotating frame's rotation sense.** It is the
*only* planar orientation angle a GTO needs in this frame — there is no
separate argument of perigee, because the rotating frame already fixes the
Earth–Moon line as the reference. The physical GTO itself never changes
(350 × 35,786 km, `sma_km = 24446`, `ecc = 0.724781150290436`, fixed for
every sheet); only its orientation relative to the Moon varies, sampled at
four values: `{0°, 90°, 180°, 270°}`.

Departure **phase** (`sD_frac` / `depFrac`) is a separate axis: the **time
fraction since perigee** over one Kepler period (mean-anomaly fraction, not
true-anomaly fraction) — 0 = perigee, 0.5 = apogee. Every entry's departure
state is the algebraic ellipse state at that phase, at that orientation; no
propagation is needed to reconstruct it (`get_family_orbit('gto', p)`, the
`'gto'` pseudo-family, Task 1 of this campaign).

## What it covers

- **Departure:** fixed 350 × 35,786 km GTO at 4 orientations `{0°, 90°, 180°,
  270°}` (spec's flagship geometry, orientDeg = −25°, is reachable via the
  picker's nearest-sheet warning, not a stored sheet).
- **Arrival:** tulip orbits Np ∈ {5, 7, 9, 12}, pm = −1 branch (pm = +1 is
  the exact z-mirror, per catalog convention).
- **Phasing:** 12 departure (anomaly) × 6 arrival phase grid per sheet.
- **Thrust:** rungs **[15 12 10 7 5] N** at Isp 1710 s, m0 150 kg — a
  5-rung ladder, truncated above the siblings' 9-rung [...1] N ladder. See
  Scope below.
- **Totals:** 16 sheets (4 orientations × 4 Np), **2,625 entries**, **840 of
  1,152 phase pairs solved (73%)**, 136 full 5-rung ladders.
- **Sample range** (orientDeg=270, Np=5, 15 N cell): t_f = 0.2951 ND =
  1.308 days, m_f = 0.32612, ΔV = 18.79 km/s — all three derive-registry
  formulas (`days`, `m_final`, `deltav_kms`) evaluate cleanly against this
  catalog's `.constants`/`.thruster`.
- **Catalog-wide ranges** (all 2,625 entries, every sheet/rung, via the
  same three derive-registry formulas): t_f = 0.2343 – 0.7155 ND =
  1.04 – 3.17 days; ΔV = 8.11 – 26.48 km/s; m_f = 30.9 – 92.5 kg (mass
  fraction 0.2061 – 0.6164, m0 = 150 kg). Fastest/lightest-ΔV entries sit
  at the 15 N rung, slowest/heaviest-ΔV at the 5 N rung, as expected for a
  thrust-limited all-burn min-time family.

### Coverage by orientation (the headline physics finding)

| orientDeg | pairs solved (of 288) | coverage | entries | full ladders |
|---|---|---|---|---|
| 0°   | 240 | 83% | 726 | 27 |
| 90°  | 207 | 72% | 659 | 40 |
| 180° | 156 | 54% | 400 | 6  |
| 270° | 237 | 82% | 840 | 63 |

**The π-dip:** orientation 180° — apogee pointed toward the Moon — is the
hard corner across *every* petal count (53–56% per-sheet coverage, only 1–2
full ladders per sheet), while 0°/270° each clear 82–88%. This is not a
mesh artifact: Diagnostic A (see Provenance) traced the *initial* failure
mode to a cold-start defect that hit apogee departures hardest (a spurious
slow-coast manifold at `tf0=4.0`), and the warm recipe (`tf0=0.30`,
thrust-locked) fixed the *reachable* cells everywhere; the 180° gap that
remains after the warm recipe is the genuine cold-basin difficulty of that
geometry, not a leftover solver defect. Best/worst sheets: `orientDeg=0,
Np=7` (89%, 173 entries) vs `orientDeg=180, Np=12` (53%, 95 entries).

### Per-sheet coverage (from `costate_lib_describe`)

```
sheet                pairs solved   entries   full ladders
  orient=  0 Np= 5     62/72 ( 86%)     197      8
  orient=  0 Np= 7     64/72 ( 89%)     173      5
  orient=  0 Np= 9     52/72 ( 72%)     164      7
  orient=  0 Np=12     62/72 ( 86%)     192      7
  orient= 90 Np= 5     48/72 ( 67%)     152      9
  orient= 90 Np= 7     56/72 ( 78%)     180     13
  orient= 90 Np= 9     47/72 ( 65%)     150      9
  orient= 90 Np=12     56/72 ( 78%)     177      9
  orient=180 Np= 5     38/72 ( 53%)      97      1
  orient=180 Np= 7     40/72 ( 56%)     100      2
  orient=180 Np= 9     40/72 ( 56%)     108      2
  orient=180 Np=12     38/72 ( 53%)      95      1
  orient=270 Np= 5     61/72 ( 85%)     220     19
  orient=270 Np= 7     60/72 ( 83%)     216     19
  orient=270 Np= 9     63/72 ( 88%)     218     15
  orient=270 Np=12     53/72 ( 74%)     186     10
```
(Reproduced from `costate_lib_describe`'s coverage table, relabeled —
see the describe-tool caveat below for why the raw tool output should not
be quoted as-is for this catalog.)

## Conjugate census

**2,625 / 0 / 0** (pass / fail / notrun) — every accepted entry passed
`ms_conjugate_test` (K=24, free-time quotiented Jacobi), zero interior
sign changes measured on any of them (`conj_ncross = 0` everywhere).
Re-solve fidelity (the honesty gate, |z − z8| < 1e-6) ranged
9.6e-16 .. 1.1e-8 across the full sweep — comfortably inside gate. This is
a stronger result than the DPO→tulip catalog's 3,931/1/0 (one refuted
entry); two structural differences are recorded as *plausible, unconfirmed*
explanations (not proven mechanisms): the warm-recipe two-point multistart
here vs a cold-start ladder there, and this catalog's 5 N thrust floor vs
the siblings' 1 N floor (deep, many-revolution rungs are where extra
junction structure is more likely to appear). Verdicts are stored in the
catalog itself (`conj_pass`/`conj_ncross`/`conj_atfinal` per sheet,
`conj_test` provenance at the top level) — see
`../../.superpowers/sdd/2026-08-26-gto-tulip-catalog/task-7-report.md` for
the full sweep record.

## Schema notes (read this before writing a generic consumer)

- **`sheets(k).tauDRO` / `.tau_dep` holds `orientDeg` (degrees), NOT a
  period.** This is a deliberate review-mandated override (all three
  reviewers, Task 3): the picker selects sheets by nearest `(tau_dep, Np)`,
  and the Kepler period of this fixed GTO is *identical across all four
  orientations* — keying on it would make 3 of 4 sheets unreachable. The
  true Kepler period lives in `dep_params.sma_km`/`.ecc` (reconstruct via
  `get_family_orbit('gto', dep_params)`), never in `tauDRO`.
- **The shared `costate_lib_describe` tool's labels are WRONG for this
  catalog.** It prints `"DRO periods (tau = period): ND [0.00 90.00 180.00
  270.00]"` and a bogus `"days [0.00 398.94 797.88 1196.82]"` row (it
  multiplies the degree-valued key by `tStar_s/86400` as if it were an ND
  time). **Ignore both rows entirely** — the sheet labels in its coverage
  table (`"tau=90.00 Np=7"`, etc.) are `orientDeg` values mislabeled as
  `tau`, but the coverage/entries/full-ladder *numbers* themselves are
  correct and meaningful. This is a known, disclosed limitation of the
  shared tool (out of scope to fix generically here), not a defect in the
  catalog.
- **Use THIS folder's `costate_catalog_pick.m`, not the shared/generic
  copy** vendored into the other deliverables. The generic picker's
  departure-axis distance metric is `(log(tauDRO) - log(tauDep))^2` —
  safe for every period-keyed family (DRO/halo/DPO, whose periods are
  never zero) but **broken for `orientDeg = 0`**: `log(0) = -Inf` makes
  the *correct* sheet's distance `NaN` while every wrong sheet scores
  `+Inf`, and MATLAB's `min` ignores `NaN`, so it silently returns the
  wrong sheet for a full quarter of the sheet space. This catalog's own
  `costate_catalog_pick.m` uses a plain linear (zero-safe) distance on the
  departure axis instead, and reconstructs the true GTO Kepler period from
  `dep_params.sma_km` for phase-day conversion rather than dividing by
  `tauDRO`. Verified: all four orientations pick exactly (max|Δz8| = 0)
  through this fix, including two independent non-trivial phase requests.
- `dep_family = 'gto'`, `dep_params = {sma_km, ecc, orientDeg}` per sheet
  (sheet-local, not homogenized — the audit rebuilds each entry from its
  own sheet's `dep_params`).
- `catalog_schema('validate', cat_)` returns `{}` (0 problems).

## Scope (v1) and the identifiability rule for the future mN extension

**This catalog stops at 5 N.** The 3–1 N legs of the standard 9-rung ladder
were attempted and produced **zero entries** by both the warm
(thrLock+tf0=0.30) and cold-mop-up recipes across the full campaign — a
*measured* closure wall, not an attempt-budget artifact (every cell's `ATT`
counter shows it was tried; the densify-specific `ATTD` counter shows
activity only at the 15 N rung). Mike adjudicated (2026-08-27) to ship the
5-rung `[15 12 10 7 5]` fleet as deliverable-7 v1 and split the 3–1 N
closure off as its own deep-rung investigation — not attempted here.

**Binding rule for that future extension (spec, verbatim):** compact `z8`
is legitimate at *few-rev, standard-thrust* entries like this catalog's —
but **many-rev (mN-regime) entries must ship full multiple-shooting
junction states, not bare z8.** This was adjudicated the hard way on this
same campaign's flagship min-time anchor: two independently-converged
`z8` seeds at the 25 mN / ~40-rev regime agreed to only ~1e-4 ND in `t_f`
and ~100 km in arrival position from the *same* bare-z8 flight — not two
extremals (a head-to-head flight proved it is one, `OPTIMALITY_CERTIFICATION.md`
§6, 2026-08-26 two-root adjudication), but bare `z8` alone could not tell
the difference. Any extension of this catalog toward deep rungs or the
25 mN flagship regime must carry the full junction-state seed, not just
the terminal `z8`, or it will inherit the same non-identifiability.

## Promotion path

`gto_entry.m` is this campaign's single-entry driver (Mike's product
request, 2026-08-26): one phasing cell, the same engine + gates as the
fleet, never a bespoke solver. Per this repo's migration rule — code
promotes to a shared library once a *second* campaign wants it — a
generic `costate_common/catalog_entry(family_dep, family_arr, ...)` is the
natural promotion once another campaign asks for the same single-entry
capability `gto_entry` provides here. Not built yet; this file is the
reference implementation to generalize from.

## How it was made

Direct collocation solve (cold, at 15 N) → warm-recipe walk down the
5-rung ladder (two-point top-rung multistart: `tf0=0.30` thrust-locked as
the primary recipe, cold `tf0=4.0` as scoped mop-up on cells the warm
recipe itself attempted and missed) → multiple-shooting refinement
(`ms_tfmin`) → `pumpkyn.cr3bp.tfMin` acceptance (three gates: ms residual,
flown arrival <100 km, |Δz| < 1e-6) → conjugate test. Root cause of the
original cold-start failure (Diagnostic A, 2026-08-26): the DRO-family
default seed `tf0 = 4.0` ND (vs a true short-transfer time of ~0.29 ND)
combined with free throttle opened a spurious slow-coast manifold — **not**
a mesh-resolution problem (escalating mesh resolution made it *worse*, the
opposite of the initial hypothesis). Fleet run 2026-08-27 → 2026-08-29
(~40 h wall, two nohup driver restarts survived cleanly via resume);
independent 3-entry spot-audit (`audit_gto_entries.m`) replayed end-to-end
across 3 orientations and 2 rung classes, all pass both hard gates.

## Gotchas, honestly

- **27% of phase pairs have no solution at any rung** (840/1,152) —
  consult `sheets(k).has_solution` (or the picker's warnings) before
  trusting a request; orientDeg=180 is the weak quadrant.
- The generic/shared `costate_catalog_pick.m` and `costate_lib_describe.m`
  vendored elsewhere in this repo will silently mis-serve this catalog at
  `orientDeg = 0` and mislabel every coverage row — always use this
  folder's own picker, and read the describe-tool caveat above before
  quoting its printed labels.
- `z8` is never interpolated between thrust rungs — only `tf` is (linear);
  you always get a converged vector at the nearest stored rung.
- This is one physical GTO sampled at 4 orientations, not 4 different
  GTOs — `sma_km`/`ecc` are identical on every sheet; only `orientDeg`
  varies.

Casey / Koblick, 2026-08-29.
