# orbit_transfer — status and roadmap

**What this document is.** The single case-by-case answer to three questions:
*what have we solved, by which method, and how well do we know it is right* —
then *what remains*. Written 2026-08-23 by reading the campaigns, the
deliverable READMEs and the certification register, not from memory.

**What it is not.** It does not restate theory (`../OCP_UNIFYING_MATH.md`),
code layout (`README.md`, `CODE_STRUCTURE.md`), the certification detail
(`OPTIMALITY_CERTIFICATION.md`), or per-campaign narrative (`*/process/`). It
is the map above those.

| document | owns |
|---|---|
| **this file** | case-by-case status, coverage gaps, the roadmap |
| `README.md` | folder layout, what each campaign is |
| `TODO.md` | the actionable open items (short form of Part 4 here) |
| `OPTIMALITY_CERTIFICATION.md` | the verification register — instruments, per-row verdicts, blocking mechanisms |
| `doc/transfer_problem_space.md` | the *goal map*: which pairs are reachable at all |
| `../OCP_UNIFYING_MATH.md` | the one Bolza problem and the four solution routes |

---

## 1. The frame

Every campaign here instantiates the same problem — minimum-X low-thrust
transfer between two orbits — and differs only along four axes:

| axis | values |
|---|---|
| **orbit pair** | DRO/halo/DPO/L1-halo → tulip/halo; GTO → tulip/ELFO; elliptic → GEO |
| **objective** | min-time · min-energy (∫‖u‖²) · min-fuel (∫‖u‖) |
| **method** | direct (collocation NLP → IPOPT) · indirect (PMP shooting) |
| **thrust** | 25 mN … 15 N (two very different regimes — see §3.1) |

Two distinct product lines have grown out of this:

- **The costate catalogs** — libraries of converged PMP costates keyed by
  (orbit pair, phasing, thrust), consumed by Darin's `pumpkyn.cr3bp.tfMin`.
  Min-time, 1–15 N, high-revolution-count-free (0.2–20 day transfers).
- **The transfer campaigns** — deep single-problem studies reproducing
  published benchmarks, mostly min-fuel, mostly many-revolution spirals.

They share `costate_common/`, `cr3bp_common/`, `verify_common/` and the
`../oclib/+oc` package.

---

## 2. Part 1 — the case table

### 2.1 Costate-catalog cases (min-time, shipped to Darin)

All four use the **same pipeline**: direct collocation solve → covector harvest
→ `ms_bvp` multiple shooting → `tfMin` acceptance. So every one is
**direct AND indirect** — that is the point of the pipeline, not a coincidence.

| # | case | entries | pairs solved | thrust rungs [N] | grid | method | status |
|---|---|---|---|---|---|---|---|
| 1 | **DRO → tulip** (coarse) | 3,936 | 32–35 of 36 per sheet (89–97%) | 15 12 10 7 5 3 2 1.5 1 | 4 τ × 4 Np × 6×6 phasing | both | **shipped** (deliverable 3) |
| 2 | **DRO → tulip** (fine sheet) | 1,105 | 130 of 144 | 15 12 10 7 5 3 2 1.5 1 0.9 0.8 0.7 0.6 **0.5** | τ=1, Np=7, 12×12 phasing | both | **shipped** (deliverable 2) |
| 3 | **HALO(L2 S) → tulip** | 3,980 | 530 of 576 (92%) | 15 12 10 7 5 3 2 1.5 1 | 4 τ × 4 Np × 6×6 | both | **shipped** (deliverable 4) |
| 4 | **DPO → tulip** | 3,932 | 511 of 576 (89%) | 15 12 10 7 5 3 2 1.5 1 | 4 τ × 4 Np × 6×6 | both | **shipped** (deliverable 5) |
| 5 | **L1 halo → L2 halo** | 1,952 | 258 of 288 (90%) | 15 12 10 7 5 3 2 1.5 1 | 2 τ_dep × 4 τ_arr × 6×6 | both | **shipped** (deliverable 6, schema v2) |
| 6 | **L2 halo → L1 halo** | 2,096 | 277 of 288 (96%) | same | same | both | **shipped** (deliverable 6) |

**Totals:** 15,896 catalog entries + 1,105 fine-sheet entries = **17,001
accepted min-time PMP solutions.** Common propulsion for all of them:
**Isp 1710 s, m₀ 150 kg.** Every entry passes three gates — multiple-shooting
residual, flown arrival (<100 km), and **acceptance unchanged by `tfMin`**
(|Δz| < 1e-6, observed ~1e-9).

Physics headlines worth carrying: halo departures are cheapest (0.65 km/s best,
vs 0.76 DPO / 0.92–0.98 DRO); solvability improves with departure period; the
hard corner everywhere is shortest-departure × longest-arrival; L1→L2 and
L2→L1 are **measurably asymmetric** (90% vs 96% solved), so direction is a
real axis, not a symmetry.

### 2.2 Beyond min-time (the catalog line)

| case | objective | status |
|---|---|---|
| DRO → tulip, 3 flagship cells | **min-energy**, fixed t_f = γ·t_f^min, γ ∈ {1.1, 1.2, 1.4} | **pilot complete 2026-08-14, 5/5 cells through all seven gates.** Not a catalog — a proof the pipeline generalizes off min-time |
| any | min-fuel | **not started.** Explicitly deferred by Darin (Aug 2026) |

### 2.3 Transfer-campaign cases (deep single-problem studies)

| case | thrust | direct | indirect | objectives |
|---|---|---|---|---|
| **elliptic → GEO, 2-body** (HMG-2004 benchmark) | 10 → 0.1 N ladder | **certified** — MEE/L-domain, full Table-3 ladder, `run_gergaud` front door | **external only**: MfMax v0/v1 built + validated as cross-check. *Our own MATLAB indirect: not built* | min-fuel (+ min-time anchors) |
| **elliptic → GEO, CR3BP** (with lunar gravity) | 10 → 0.1 N ladder | **certified** — `run_cr3bp_geo`, reviewed technical note | **not started** (Phase 2) | min-fuel |
| **GTO → tulip** | 25 mN, ~40-rev spiral | **certified** — Sundman min-fuel engine, 25-switch flagship, ΔV–t_f front | **built but NOT certified** — `ms_band`, `ifs`, `ztl`, `min_time`; stalls short | all three (one homotopy chain) |
| **GTO → ELFO** | 25 mN | **working** — min-fuel front mapped, min-time anchor certified | **not started** (Route C) | min-energy, min-fuel, min-time |
| **DRO → tulip** (the bridge case) | 1–15 N | **certified** min-time | **certified** min-time | min-time, min-energy pilot |

**DRO → tulip is the only case where direct and indirect have been compared
head-to-head, and they agree**: t_f 4.0152501 vs 4.0152425 (5 significant
figures), primer directions to 1.2e-6 deg, costate scale factor 0.99999, state
histories to 0.396 km peak. That is the strongest single verification result in
the repository, and it is the template for everything else.

---

## 3. Part 2 — what the table shows

### 3.1 The thrust gap is the biggest structural hole

Two disjoint regimes, with nothing in between:

```
 25 mN ......................... 0.5 N ... 1 N ————————————— 15 N
   ^                               ^        ^                  ^
   GTO campaigns                   one      the four catalogs
   (15 kg, ~40 rev)                sheet    (150 kg, 0.2-20 d)
   NO catalog                      only
```

- The catalogs bottom out at **1 N** (one fine sheet reaches 0.5 N).
- The GTO campaigns live at **25 mN** — a factor of 40 below the catalog floor
  — and have **no catalog at all**.
- Sub-0.5 N rungs were deferred by Darin "until customers need it," but that
  deferral is about *catalog product*, not about our own need.

Why it matters beyond coverage: low thrust means **many revolutions**, and
revolution count — not the objective — is what decides whether indirect
shooting closes. The catalogs are easy precisely because they are 0.2–20-day,
few-revolution transfers. Extending them downward is not "more of the same";
it walks into the regime where the GTO campaigns already stall.

### 3.2 Objective coverage is thin off min-time

| objective | catalog cases | campaign cases |
|---|---|---|
| min-time | 6 of 6 shipped | anchors everywhere |
| min-energy | 1 pilot, 3 cells, 1 pair | GTO_tulip, GTO_ELFO (as homotopy roots) |
| min-fuel | **none** | 4 campaigns certified |

The two product lines are almost complementary — the catalogs are min-time and
the campaigns are min-fuel — and the bridge between them is the energy→fuel
homotopy, which exists in the campaigns and has never been run on a catalog
seed.

### 3.3 Method coverage is asymmetric

Direct leads almost everywhere. Indirect exists in force only in the catalog
pipeline (where it is *seeded by* the direct solve) and as external MfMax.

| | direct | indirect |
|---|---|---|
| certified | 5 campaigns + 6 catalogs | 6 catalogs (min-time), 1 campaign pair |
| built, not certified | — | GTO→tulip (ms_band/ifs/ztl) |
| not started | — | earth CR3BP, GTO→ELFO, our own earth 2-body |

### 3.4 Family-pair coverage

Nine pumpkyn families are catalogued; we have shipped **4 distinct pair
types** out of a much larger reachable set. Untouched and named as interesting:
GEO→DRO (Darin's own example), LPO→tulip, DRO→halo, Lyapunov↔halo (a published
benchmark, so externally checkable), cycler→anything.

---

## 4. Part 3 — verification: what we can actually claim

### 4.1 The ladder

Four rungs, in increasing strength. Knowing which rung a result sits on is the
whole game.

| rung | claim | instrument | where we are |
|---|---|---|---|
| 1 | **feasible on the transcription** | NLP defect ~1e-14 | everywhere — *and nearly worthless alone* (see 4.2) |
| 2 | **solves the continuous ODE** | G1/G1b re-integration, flown control | DRO→tulip only; **generalizing it to the other campaigns is the obvious next move** |
| 3 | **PMP-extremal** (first order) | `foc_check`/`foc_report`, primer + sign law + transversality + KKT + Ṡ≠0 | **all four campaigns**, coverage uneven |
| 4a | **weak local min** | IPOPT native inertia (δ_w = 0) | all four campaigns, mixed verdicts per row |
| 4b | **strict local min** | Maurer–Osmolovskii switching-time Hessian | **nowhere** |

### 4.2 The lesson that reframes rung 1

Measured on DRO→tulip: NLP defect **1.4e-14** while the true continuous local
error was **1.5** — a ratio of **1e7**. And refinement was **not monotone**:
N=400→800 made it worse and put a node 719.6 km inside the Moon; only N=1600
landed on the reference. *A machine-tight defect is a statement about the
discretization, not the trajectory.* Every campaign quotes 1e-14 defects; that
evidence is necessary and nowhere near sufficient.

### 4.3 Independence — where we have it and where we do not

The register's gap **G1** says the Earth and CR3BP campaigns verify PMP from
the NLP's *own* duals: self-consistency, not an independent witness. It catches
convergence to a non-extremal; it cannot catch an error shared with the
transcription.

**The costate pipeline has already solved this for its own entries, and that
deserves to be stated plainly.** Its chain is:

```
direct collocation solve  →  harvest duals as a SEED  →  ms_bvp multiple
shooting (converges on its own PMP residual, 1e-11)  →  tfMin acceptance
(third-party code, |Δz| < 1e-6)
```

The `ms_bvp` root satisfies the boundary-value problem regardless of where the
seed came from, and `tfMin` is Darin's independent implementation. So catalog
entries carry a genuinely independent witness that the campaign results do not.
**That chain is the template for closing G1 elsewhere** — it is the same thing
"seed our own MATLAB indirect solver from the certified direct solution" would
buy the Earth campaigns.

### 4.4 Second order — the honest position

- **Weak local minimality** is certified per-row on all four campaigns by
  IPOPT's inertia (δ_w = 0). Mixed: e.g. earth 2-body is LOCAL MIN at
  10 / 0.2 / 0.1 N but NOT CERTIFIED at 5 / 2.5 / 1 / 0.5 N; the tulip
  25-switch flagship is borderline NOT CERTIFIED (δ_w 1.8e-8 vs 1e-8 tol).
- **Strict local minimality is certified nowhere**, and this is structural,
  not lazy: the fuel objective is linear in throttle, so the reduced Hessian is
  genuinely flat (270 flat directions at earth 10 N). **Strictness for
  bang-bang lives in the switching times**, which no NLP-level test can reach.
- **A conjugate-point test exists** — `costate_common/ms_conjugate_test`,
  free-final-time Jacobi form, handling the two exact degeneracies (costate
  scaling, time reparametrization). It is the pipeline's first genuine
  second-order instrument.

> **Two staleness findings from writing this document.**
> 1. ~~`OPTIMALITY_CERTIFICATION.md` Part B §1 still lists the conjugate-point
>    test as "not built; gated on an indirect solver."~~ **FIXED 2026-08-23**:
>    the row now records it as built 2026-08-08 (migration #3) with its three
>    run sites, and an experiment-log entry documents the correction.
> 2. ~~None of the 17,001 shipped entries carries a conjugate-point verdict.~~
>    **CLOSED 2026-08-23**: the sweep ran on all 15,896 catalog entries
>    (`conj_catalog_pass`, ~33 min) — **15,895 pass / 1 fail / 0 unverified**,
>    verdicts stored in the catalog .mats under the schema's new optional
>    `conj_pass` fields. The one refuted entry (DPO τ=2 → Np=7, (2/3,2/3),
>    15 N) has a genuine interior conjugate point. Remaining: the deliverable
>    zips predate the verdicts, and the 1,105-entry deliverable-2 fine-sheet
>    library (older format) is not yet swept.

---

## 5. Part 4 — what remains

Grouped by the goals as stated, each with the concrete next action.

### A. More thrust levels

| | |
|---|---|
| **now** | catalogs 1–15 N; one fine sheet to 0.5 N; campaigns at 25 mN |
| **want** | continuous coverage down to the mN regime |
| **do** | (1) extend the four catalogs to 0.5 N using the existing `extend_thrust_ladder`/`densify_ladder` machinery; (2) pilot a single cell at 0.1 N and at 25 mN to find where `ms_bvp` stops closing — that boundary is the real deliverable, and it is currently unknown |
| **risk** | low thrust ⇒ many revolutions ⇒ the regime where GTO→tulip indirect already stalls. Expect a wall; locating it is the point |

### B. A GTO → tulip costate library — BUILT 2026-08-29

| | |
|---|---|
| **now** | **catalog SHIPPED at the Darin-standard propulsion regime** (150 kg / 1710 s, few-rev): `GTO_tulip/catalog/results/costate_catalog_gto_tulip.mat`, 16 sheets (4 orientations × 4 Np), **2,625 entries, 840/1,152 phase pairs (73%)**, rungs `[15 12 10 7 5]` N, schema v2 valid, conjugate census **2,625/0/0** (stronger than DPO's 3,931/1/0). Coverage by orientation: 240/207/156/237 of 288 (83/72/54/82%) — a real π-dip at orientDeg=180° (apogee toward the Moon): Diagnostic A traced the *initial* failure mode there to a cold-start defect, not a mesh artifact; the residual gap after the warm-recipe fix is *inferred* (not traced) to be genuine cold-basin difficulty of that geometry. |
| **how it got unblocked** | the sequence this row used to prescribe (§6 step 6) ran as planned: `ms_bvp` closed on the GTO→tulip min-time problem at K=60 (2026-08-25, Stage A), agreeing exactly with the certified direct anchor (t_f = 6.290694 ND); the catalog fleet (Stage B, this row) then wrapped the **standard 1–15 N regime first** (adjudicated 2026-08-25), not the 25 mN flagship — the flagship-regime catalog is still open, see below. |
| **the closure wall** | the ladder was attempted down to 1 N; **3–1 N produced zero entries** by both the warm (thrLock, tf0=0.30) and cold-mop-up recipes across the full fleet — a measured wall (every cell's attempt counter shows it was tried), not a budget artifact. Deliverable-7 v1 ships the 5-rung `[15 12 10 7 5]` fleet; the 3–1 N deep-rung leg is split off as **the open item** below. |
| **the 25 mN flagship regime** | still **no catalog** — this campaign deliberately stayed at the few-rev standard-thrust regime the pipeline is proven at (spec adjudication 2026-08-26); the flagship's ~40-rev/25 mN regime carries its own binding rule: **many-rev entries must ship full ms junction states, not bare z8** (measured on this same campaign's min-time anchor — two independently-converged bare-z8 seeds agreed to only ~1e-4 ND in t_f / ~100 km in arrival before a head-to-head flight proved they were one extremal, not two; `OPTIMALITY_CERTIFICATION.md` §6, 2026-08-26 two-root adjudication). A future mN extension of this catalog must apply that rule from the start. |
| **do next (open item)** | the **deep-rung investigation** (3–1 N closure wall) — root-cause it before attempting any deep-rung extension of this or any catalog; densify-style grinding at deep rungs already showed no per-attempt wall-time cap in this campaign (13 h churn, 8.3 MB log — needs a `wallSec`/attempt cap first) |
| **product + record** | catalog README `GTO_tulip/catalog/README.md`; full build/pilot/diagnostic/fleet/audit/conjugate-sweep record in `.superpowers/sdd/2026-08-26-gto-tulip-catalog/` (progress.md + task 1–8 reports) |

### C. Both direct and indirect everywhere

| case | missing leg | note |
|---|---|---|
| earth elliptic → GEO (2-body) | our own MATLAB indirect | MfMax's converged costates are the clean seed source; direct-KKT duals are unreliable at high eccentricity |
| earth elliptic → GEO (CR3BP) | indirect | Phase 2, not started |
| GTO → tulip | indirect **certification** | built, stalls; the 1.01–1.11× near-min-time band is the open wall |
| GTO → ELFO | indirect | Route C, not started |

Borrow from MfMax when this restarts: arclength continuation for energy→fuel
instead of our hand-scheduled ε ladder; the start-at-the-target
initial-condition homotopy; PMP-converged costates as seeds. See
`earth_elliptic_to_geo/indirect/mfmax/MFMAX_V1_RUNBOOK.md`.

### D. Clear, independent verification that we have an extremal

| | |
|---|---|
| **now** | first-order gate on all four campaigns, but **self-referential** (G1) except in the catalog pipeline |
| **do** | (1) **generalize the G1/G1b continuous-residual gate** from DRO→tulip to the other three campaigns — rung 2 of the ladder is currently one-case-only, and rung 1 is proven misleading; (2) replicate the catalog pipeline's independent chain on a campaign result: take a certified direct solution, harvest, multiple-shoot, and require an independent PMP root; (3) close the uneven `foc_check` coverage (earth CR3BP 4 of 7 rungs; ELFO 1 rung, not a ladder; tulip PSR front never exercised end-to-end) |
| **also** | adjudicate the recorded tulip disagreement: raw-dual PASS (0.058°) vs LS-reconstruction FAIL (2.00°) on the same flagship — recorded, never resolved |

### E. Methods for showing a solution is a local minimum

| | |
|---|---|
| **now** | weak local min on four campaigns; strict nowhere; a conjugate-point test built but unapplied at scale |
| **do (cheap, high value)** | **run `ms_conjugate_test` across the shipped catalogs and store the verdict in the schema.** The instrument exists, the STMs come free from the `ms_bvp` Jacobian machinery, and it would upgrade 17,001 entries from "extremal" to "extremal with no conjugate point in (0, t_f)" — the necessary second-order condition |
| **do (the real build)** | the **STM / multiple-shooting switching-time Hessian**, already specified in `switch_hessian.m`'s header: keep the collocation trajectory as the feasible base, get switch-time sensitivities from the STM, form the second variation with segment matching so nothing propagates across all revolutions. It is the only live path to strictness, it is campaign-agnostic, and it shares its variational core with the conjugate test |
| **target statement** | positive-definite projected switching-time Hessian + Ṡ≠0 at every switch + primer/Legendre direction optimality = the Maurer–Osmolovskii sufficient condition for a bang-bang local minimum. Ingredients (ii) and (iii) already exist and are trustworthy; only (i) is missing |
| **also open** | a **fixed-t_f** conjugate test — `ms_conjugate_test` is free-time-specific, so the min-energy entries have no second-order verdict at all |

### F. Objectives beyond min-time in the catalogs

1. γ grid per cell — J and m_f are **non-monotone in γ** (measured: m_f 0.9379
   → 0.9421 → 0.9374 for γ = 1.1 → 1.2 → 1.4), so the fixed-time problems at
   different γ are not nested and basin discipline applies.
2. A γ / t_f axis in `catalog_schema`.
3. Then the **energy→fuel ε-homotopy on the same seeds** — the route to
   min-fuel catalog entries, reusing the campaigns' machinery.

### G. More orbit pairs

Densify first (halo 46 unsolved pairs, DPO 65, DRO s_A=0.075 row + 6
stragglers), grow the period/petal axes toward the measured admissible box
(DRO τ 0.05–3.25, tulips Np 3–14 — both wider than what is shipped), then new
pair types: GEO→DRO, LPO→tulip, DRO→halo, Lyapunov↔halo.

---

## 6. Recommended sequence

Ordered by (value delivered) ÷ (work required), not by ambition:

1. ~~**Run the conjugate test across the shipped catalogs, store it in the
   schema.**~~ **DONE 2026-08-23**: 15,896 entries, 15,895 pass / 1 fail,
   ~33 min. Residue: sweep the deliverable-2 fine-sheet library (1,105
   entries, older format); re-ship deliverables when Darin wants the verdicts.
2. ~~Generalize the G1/G1b continuous-residual gate~~ — **PARTIAL
   2026-08-25**: engine promoted to `oc.local_residual` (TDD, injected-error
   locality), `dro_residual` rerouted (equivalence 1.7e-18), and the earth
   2-body MEE campaign got its first-ever gate (`mee_residual`, 10 N row:
   P median 2.2e-13 ND / max 39.9 km, **all large residuals on
   throttle-switching intervals** — confirms the PSR rationale). CR3BP-GEO
   wrapper DONE 2026-08-26 (`cr3bp_residual_gate`; same switch-interval
   concentration, small 2.6 m interior floor from the lunar term;
   `mee_residual` promoted to `verify_common` on its second consumer). Ladder sweep DONE
   2026-08-26 (14/14 rows: switch attribution total, per-switch error
   shrinks with depth — growth prediction refuted; CR3BP interior floor
   grows 2.6→75 m, t-row/lunar-phase hypothesis;
   `verify_common/doc/g1_sweep_results.md`). Open: the ELFO wrapper only.
3. ~~**Update the two stale records** found here (§4.4)~~ — **DONE
   2026-08-23** (register conjugate row + log entry; `transfer_problem_space`
   tables refreshed; `campaign_status.tex` marked superseded).
4. **Extend the catalogs to 0.5 N, then probe 0.1 N and 25 mN on one cell.**
   Locate the closure wall.
5. **γ grid → γ axis → energy→fuel homotopy** on catalog seeds.
6. ~~`ms_bvp` on GTO→tulip min-time~~ — **DONE 2026-08-25 (Stage A)**:
   shared engine closes at K=60 in 10.3 s, ‖R‖=2e-11, tfMin accepts, conjugate
   PASS, and **t_f = 6.290694 ND lands exactly on the certified direct
   reference** (the 2026-07-13 module root was 47 s off — shallower root).
   Direct↔indirect agreement is now closed on the flagship min-time anchor.
   Stage B (the GTO→tulip catalog, standard 150 kg / 1–15 N first per
   adjudication 2026-08-25) — **DONE 2026-08-29**: 2,625 entries / 16 sheets /
   73% pair coverage / conjugate 2,625/0/0, shipped at rungs `[15 12 10 7 5]`
   N (see §5.B). Dual capture in `gen_tulip_mintime` for harvest seeds is
   still open only for a future flagship-regime (25 mN) extension of this
   catalog, not for the shipped v1.
7. **The GTO deep-rung investigation** (3–1 N closure wall, split off from
   Stage B — see §5.B) — root-cause before any deep-rung catalog extension.
8. **The switching-time Hessian.** The largest piece, and the only path to a
   strictness claim worth putting in a paper.

---

## 7. Provenance

Compiled 2026-08-23 from: the four deliverable READMEs (entry counts, thrust
rungs, coverage percentages, ΔV/flight-time ranges); `OPTIMALITY_CERTIFICATION.md`
Parts A and B (instruments, coverage matrix, verdicts, blocking mechanisms);
`doc/transfer_problem_space.md` (goal map, direct↔indirect agreement numbers);
each campaign's `README.md` and `TODO.md`; `DRO_tulip/FINDINGS.md` (min-energy
pilot); and direct inspection of `costate_common/` (conjugate-test call sites,
schema fields). Numbers are as recorded in those sources; where this document
and an older one disagree, the disagreement is flagged in §4.4 rather than
silently resolved.
