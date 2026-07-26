# earth_elliptic_to_geo_CR3BP — elliptic → GEO with lunar gravity

**Goal.** Re-solve the low-thrust minimum-fuel elliptic-orbit → GEO transfer
(the Haberkorn–Martinon–Gergaud benchmark: 1500 kg, P⁰=11625 km, e⁰=0.75,
i⁰=7° → equatorial GEO, thrust ladder 10 → 0.1 N) **with the Moon's gravity
incorporated** — Earth–Moon CR3BP dynamics instead of the 2-body `1/r²` model
used in `../earth_elliptic_to_geo/` — and answer quantitatively: *how much
does lunar gravity move the certified 2-body answers?*

## Status: Phase 1 COMPLETE (2026-07-24)

The full thrust ladder is solved, certified, compared, documented, and
visualized. The answer to the campaign question is a measured curve, not a
single number. What remains open (indirect solver, deeper optimality
certificates, deep-rung φ₀ sweep) is in `TODO.md`.

## Entry point

`direct/run_cr3bp_geo.m` — the front door. Set parameters (thrust, lunar
phase φ₀, gain, epsMin (1=energy/0=fuel), endpoints, run name), run, and get
the solved transfer + saved data products (states/controls/costates/mesh) +
plots + optional movie. Validated end-to-end: reproduces the certified 10 N
result (m_f=1377.1545 kg, 19 sw) exactly. The whole ladder is walked by
`run_cr3bp_ladder.sh` (per-rung process isolation).

## Headline results (same-chain baselines throughout)

**Moon effect vs thrust** (φ₀=0, complete ladder 10→0.1 N):
+52.0 / +48.8 / +49.7 / +49.3 / +45.9 / +36.8 / +31.3 g at
10 / 5 / 2.5 / 1 / 0.5 / 0.2 / 0.1 N (t_f = 0.19 → 19.4 lunar months). Flat
~50 g through ~2 lunar months, then a gentle decline to +31.3 g at 19 months
— consistent with a phase-averaging picture (untested until a deep-rung φ₀
sweep runs).

**Moon effect vs lunar phase** (10 N): +52.0 / −56.3 / +68.9 / −52.8 g at
φ₀ = 0, π/2, π, 3π/2 — **π-periodic tidal-quadrupole signature: the lunar
assist flips sign with phase.** Harmonic fit
Δm_f ≈ +3.0 − 8.5cosφ₀ − 1.75sinφ₀ + 57.5cos2φ₀ g; ~125 g peak-to-peak from
departure-epoch choice. The four control laws differ by ≤1.3 min of switch
timing (same 19-switch structure) — the mass effect is first-order along the
unperturbed extremal, the control response second-order (envelope argument).

Full table: `direct/results/compare_vs_2body.md`.

*Certification caveat:* "certified" = the four NLP metrics (defect /
unit-norm / terminal-error / IPOPT `Solve_Succeeded`) + bound-saturation
check — a first-order, solver-level certificate. The 2-body campaign's
primer/PMP check is now lunar-aware (`direct/verify_cr3bp_pmp.m`, task B) but
still fails on the pre-existing campaign-wide eccentricity-correlated
raw-dual anomaly (fails the pure 2-body solutions identically — not
lunar-specific). Genuine second-order evidence (SOSC / conjugate-point) is
future work — see `TODO.md` and the note §10.

## The note

`doc/cr3bp_geo_phase1_note.tex` — the campaign's technical note (16 pp, two
rounds of GPT-5.6 + Gemini review, "publishable-as-note"): the OCP, the
lunar incorporation (direct + indirect third-body terms), the fixed-t_f /
free-ΔL discussion, the pipeline (with a flow diagram), the μ-continuation
lineage (Bonnard–Caillau–Picot 2010), both result curves with figures, the
control-sensitivity/envelope discussion, the direct-vs-indirect working
position, and the three-tier minimality-certificate plan. Figures in
`doc/figs/`.

## Ladder solver recipe (the load-bearing lesson)

Rungs ≥2.5 N (N≳700 nodes) segfault in the bundled MUMPS factorization until
`liftDL=true` (the scalar ΔL's dense arrowhead column → block-banded KKT;
numerically identical) — then add `maxIter≥4000` (10000 at the deepest
rungs). With that recipe the whole ladder runs; the historically feared
deep-rung ε-wall never materialized. Details: note §"thrust ladder", lessons.

## Structure

```
earth_elliptic_to_geo_CR3BP/
├── README.md, TODO.md
├── doc/            cr3bp_geo_phase1_note.tex (+ pdf), figs/, reviews/
├── papers/         Bonnard-Caillau-Picot 2010 (μ-continuation theory), COCV 2001
└── direct/         MEE+lunar-pert campaign:
    ├── run_cr3bp_geo.m          front door; run_cr3bp_ladder.sh (ladder)
    ├── lunar_params.m           Moon constants in canonical units
    ├── bridge_mu_continuation.m, solve_cr3bp_minfuel.m   callable stages (see below)
    ├── compare_vs_2body.m, sanity_bound.m                analysis
    ├── verify_cr3bp_pmp.m, test_cr3bp_consistency.m      verification
    ├── viz/         fig_phi_sweep, fig_ladder_dmf, phase_quad_movie, render_ladder_outputs
    └── results/     .mat products (gitignored), figures, movies
    (indirect/ — Phase 2, not started)
```

**Two implementations of the same sequence.** `bridge_mu_continuation.m` and
`solve_cr3bp_minfuel.m` were labelled "pipeline stages" here, which reads as if
the front door calls them. It does not — measured 2026-07-26, `run_cr3bp_geo`
inlines stages A–D itself and references neither. They are a *callable-stage*
implementation of the same sequence, reached only by `verify_cr3bp_pmp.m` and
`compare_vs_2body.m`. Deciding which is authoritative (so the other can
delegate) is a Tier-2 item in `orbit_transfer/CODE_STRUCTURE.md`; until then,
treat `run_cr3bp_geo` as the campaign's front door and these two as the
analysis path.

The one shared-core edit lives upstream: `../earth_elliptic_to_geo/direct/
core/lt_mee_rhs.m` gained the opt-in `par.pert` lunar third-body branch
(absent ⇒ 2-body path character-identical). Provenance: spec
`docs/superpowers/specs/2026-07-22-elliptic-geo-cr3bp-phase0-design.md`, plan
`.../plans/2026-07-22-elliptic-geo-cr3bp-phase1.md`.
