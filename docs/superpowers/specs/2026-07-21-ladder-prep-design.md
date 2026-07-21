# Thrust-Ladder Prep (P2) — Design Spec

**Date:** 2026-07-21  **Status:** design approved (brainstorm); plan to follow.
**Scope:** both CR3BP campaigns — `GTO_tulip` (sundman engine) and `GTO_ELFO`
(freetf engine) — plus the shared `cr3bp_common/`. Source of requirements: the
2026-07-21 review triages (`GTO_tulip/doc/reviews/2026-07-21_triage.md` C4–C6,
`GTO_ELFO/doc/reviews/2026-07-21_triage.md` C4–C6): the "ladder-prep trio"
that gates the thrust-ladder goal in both TODOs.

## 1. Goal

Make both min-fuel pipelines safe to run a T_max thrust ladder (warm-chained
rungs around the nominal 25 mN): (1) per-rung thrust parameterization with
config fingerprints on every cached artifact, (2) adaptive bounds + a
bound-saturation diagnostic, (3) alias-free cross-rung warm-start chaining.
Validated by two live 20 mN pilot rungs.

**Honest scope limit (stated up front):** the tulip fixed-τf formulation
freezes the revolution count into the seed geometry (the earth campaign's
proven topology lesson), and ELFO's `cScale` gives only time-span slack, not
rev-count freedom. This prep therefore targets the **near-nominal band
(~15–40 mN)**. Deep rungs (rev-count changes) need a free-span reformulation
(the earth MEE+ΔL analog) — explicitly out of scope; recorded as the
escalation path in both TODOs.

**Back-compat invariant:** at nominal thrust with no new options, every
driver/solver path is byte-identical — same cache filenames, same NLP, same
results. All new behavior is parameter-gated or purely additive output.

## 2. Component 1 — fingerprints + per-rung provenance (`cr3bp_common`)

- **`cr3bp_fingerprint(p, extra)`** → struct with the config that determines a
  solution: `thrustN, m0kg, ispS, Tmax, cEx(=c), muStar, pSund` (+`qSund` when
  the caller has one), `rv0, rvf`, plus any `extra` fields (e.g. `tf`,
  `insertion`, `epsMin`). Built from the `cr3bp_lt_params` struct so it cannot
  drift from the physics.
- **`check_cr3bp_fp(Scached, fpNow, file, tag)`** — earth `check_cache_fp`
  semantics, generalized: cached struct lacking `.fp` → WARN and trust
  (legacy); field present in `fpNow` but absent in cached fp → WARN
  (schema-older); field present on both sides with different values →
  **ERROR** naming the field and file (fail-loud; stale/foreign cache under
  the same tag).
- **Applied at every cache boundary:**
  - tulip: `minfuel_at_tf` seed load + result save; `sundman_homotopy`
    saveFile; `gen_tulip_energy_2p` checkpoint + outFile; `gen_tulip_mintime`
    anchor save.
  - ELFO: energy-seed saves (`gravhom`, `tfsweep`), `gen_elfo_minfuel`
    checkpoint (today it checks only `tf0`) + result save, `elfo_run_one`
    rows, `gen_elfo_mintime` anchor.
  - **`elfo_find_energy_seed` gains an fp filter**: candidate seeds whose
    stored fp mismatches the current one (esp. `Tmax`) are skipped, with a
    legacy pass-through (no fp ⇒ eligible + warned) so the existing 25 mN
    seed bank keeps working at nominal.
- **Thrust as a parameter:** `minfuel_config(over)` accepts an optional
  override struct (e.g. `struct('thrustN',0.020)`) merged over the defaults;
  no args ⇒ byte-identical today's config. Drivers thread an optional
  `thrustN`/`cfgOver` option down to it.
- **Rung namespacing:** a shared tag helper `thrust_tag(thrustN)` returns ''
  at nominal (25 mN) and `'_T<mN>mN'` otherwise (e.g. `_T20mN`); appended to
  every artifact filename the drivers write. Nominal names unchanged.

## 3. Component 2 — adaptive bounds + saturation diagnostics

- **`casadi_minfuel_sundman` gains a trailing optional `opts` struct**
  (15th arg; absent ⇒ byte-identical): `opts.vBox` overrides the ±12 velocity
  box halves of `lbX/ubX` (line 114); `opts.rBox` likewise ±3 position.
  Defaults = current values.
- **`boundSat` output on BOTH solvers** (always on, output-only): per-row
  minimum slack to each nonphysical box (r, v, m-upper, cScale for freetf) at
  the returned solution; `out.boundSat = struct('minSlack',…,'worst',label)`
  plus a `warning(...)` when any nonphysical bound is within 1e-4 of binding.
  BCs/pinned rows excluded (they sit at bounds by construction).
- **ELFO `cBox` scaling in the drivers** (not the solver): the step-solve
  `base.cBox` becomes `cBoxNom .* [min(1,Tfac) max(1,Tfac)]` with
  `Tfac = T_nominal/T_max` — wider ceiling at lower thrust, wider floor at
  higher; at nominal, `Tfac=1` ⇒ unchanged `[0.15 6]`.
- **tulip `tauf0`:** never reused across rungs — recomputed by the chain
  helper (Component 3). No solver change.

## 4. Component 3 — cross-rung chaining helpers

- **`chain_rung_seed_tulip(solMat_or_struct, tfNew, pNew)`** → `[sigma, X0,
  U0, tauf0, fp]`: takes the previous rung's converged Sundman solution
  (8-state X with t = row 8), rescales ONLY the time row by `tfNew/t_end`
  (same-mesh scalar rescale — established, non-aliasing), converts to the
  time-mesh form and re-maps through **`sundman_seed_map`** (no-resample:
  controls stay attached to their own spatial nodes; fresh `tauf0`; endpoints
  pinned). No σ-interp, no t-interp of controls anywhere.
- **ELFO chaining is deliberately trivial:** `casadi_energy_freetf`'s
  `cScale` slack already decouples the clock, so a thrust rung warm-starts
  from the source `X,U,tauf0` unchanged (same nodes) with the new `Tmax`,
  scaled `cBox`, and the new fingerprint. A thin
  `chain_rung_seed_elfo(seedS, pNew)` packages that + fp for uniformity.
- Both helpers **refuse** (error) a source whose fp thrust equals the target
  (chaining to the same rung is a caller bug) and record
  `fp.chainedFrom = <source tag/thrust>`.

## 5. Component 4 — folded small fixes (review C4/C6, both campaigns)

1. **`rF` fallback cascades:** in `gen_tulip_energy_2p`>step_solve,
   `gen_elfo_minfuel`>step_solve, `gen_elfo_energy_gravhom`>leg stepper, and
   `gen_elfo_energy_tfsweep`: when the tight re-clean `rT` (after a
   successful loose probe) fails, fall back to tight-from-`Xk` (`rF`) before
   declaring the step failed.
2. **`gen_elfo_mintime`:** certification gate (`Solve_Succeeded` + defect)
   before writing the anchor file; failed solves are reported, not saved.
3. **`gen_elfo_energy_tfsweep`:** tight re-clean before `save_point` (today a
   500-iter loose solve can be banked directly).

## 6. Component 5 — validation gates

1. **No-solve unit tests** (fast, in each campaign's tests or alongside):
   fingerprint build/match/mismatch/schema-older/legacy paths;
   `thrust_tag` nominal-empty + off-nominal; `chain_rung_seed_tulip` output
   shapes, endpoint pinning, and small stencil defect against the saved
   certified solution; `boundSat` fields on a synthetic solution;
   `elfo_find_energy_seed` fp filter (legacy eligible, mismatched skipped).
2. **Nominal byte-path regression:** defaults untouched — assert
   `minfuel_config()` equals today's values, solver calls without `opts`
   reproduce the existing smoke results.
3. **Two live pilot rungs (background, resume-safe, the package's exit
   gate):**
   - **Tulip 20 mN:** `chain_rung_seed_tulip` from the certified 1.15×
     nominal solution (same t_f), energy re-clean at 20 mN, ε:1→0 sharpen via
     the hardened gates. PASS = `certified=1`, clean `boundSat`, fp with
     `thrustN=0.020` + `chainedFrom` recorded, artifacts under `_T20mN` tags.
   - **ELFO 20 mN:** `chain_rung_seed_elfo` from the f1.238 base seed, same
     recipe on the freetf engine, scaled `cBox`. Same PASS bar.
   - Pilots prove machinery, not optimality: per-rung min-time anchors (and
     factor labels vs a 20 mN anchor) belong to the ladder campaign itself.

## 7. Files (create/modify inventory)

| file | action |
|---|---|
| `cr3bp_common/cr3bp_fingerprint.m`, `check_cr3bp_fp.m`, `thrust_tag.m` | create |
| `cr3bp_common/minfuel_config.m` | modify: optional override arg |
| `GTO_tulip/direct/sundman_minfuel/casadi_minfuel_sundman.m` | modify: trailing opts (vBox/rBox), boundSat output |
| `GTO_tulip/direct/sundman_minfuel/{minfuel_at_tf,sundman_homotopy,gen_tulip_energy_2p,gen_tulip_mintime}.m` | modify: fp wiring (+ rF cascade in energy_2p) |
| `GTO_tulip/direct/sundman_minfuel/chain_rung_seed_tulip.m` | create |
| `GTO_ELFO/direct/elfo/casadi_energy_freetf.m` (+`casadi_mintime_freetf.m`) | modify: boundSat output |
| `GTO_ELFO/direct/elfo/{gen_elfo_minfuel,gen_elfo_energy_gravhom,gen_elfo_energy_tfsweep,elfo_run_one,gen_elfo_mintime,elfo_find_energy_seed}.m` | modify: fp wiring, cBox scaling, rF cascades, anchor gate, re-clean, fp filter |
| `GTO_ELFO/direct/elfo/chain_rung_seed_elfo.m` | create |
| pilot drivers `pilot_rung_20mN.m` (one per campaign) | create |
| tests (per campaign + cr3bp_common) | create |
| TODOs + triage docs | modify: mark C4–C6 items done/partial, record pilot results |

## 8. Risks

- **Solver-file edits on certified paths** — mitigated by the back-compat
  invariant (defaults byte-identical; boundSat output-only) + the nominal
  regression gate before pilots.
- **Legacy-cache friction** — the WARN-and-trust legacy path keeps the
  25 mN banks usable; only genuinely mismatched fp'd caches hard-error.
- **Pilot rungs may honestly fail** (20 mN may hit a wall the prep can't fix,
  e.g. sharpening stalls) — that is a *finding*, not a plan failure; the
  package still lands, and the wall gets recorded in the TODO with the
  attempt trajectory (house honesty rule).
