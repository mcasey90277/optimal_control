# Three-way document review — fix wave, dispositions and verification

**Date:** 2026-08-09
**Implementer:** single fix-wave agent (Claude Fable 5)
**Targets:** `doc/booster_landing_note.tex`, `doc/booster_landing_sdd.tex`
**Reviews synthesized:**
- `process/docreview-gpt56terra.md` (GPT-5.6-terra, no repo access)
- `process/docreview-gemini31pro.md` (Gemini 3.1 Pro, no repo access)
- `process/docreview-host-claude.md` (host Claude, code-verified — authoritative on measured values)

---

## 0. The fact that resolved several findings

The full flagship (`cfg.phase2 = true`) was re-run **after** the reviews took
their snapshots. `results/booster_run.mat` is now current and complete.
Verified by direct probe (`fieldnames` on the loaded struct):

```
VARS: P, Pd, ctrl, ctrlD, mc, mcD, out0, rep, repD, solC, solD, solV, when
when: 09-Aug-2026 12:56:25
```

Every number the reviews questioned reproduces from it:

| quantity | measured |
|---|---|
| `solC.tf` / `mf` / fuel | 16.595192 s / 26464.5890 kg / 3535.4110 kg |
| `solC.stats.iter` | **32** (note said 30 — fixed) |
| `solD.tf` / `mf` / fuel | 17.530799 s / 26899.2760 kg / 3100.7240 kg |
| MC phase 1 | 0.9950; landed 199 / arrest 0 / depleted 1 / horizon 0 |
| MC phase 2 (wind) | 0.9900; landed 200 / arrest 0 / depleted 0 / horizon 0 |
| `rep` G3/G4/G5 | dmf 0.427920, dtf 0.00340864, traj 0.372616, G4 1.356e-04, primer 1.207418e-06, bound frac 0.991736 |
| `repD` | primer 1.207418e-06, bound frac 1.000000, G1 2.2837e-08 |
| `out0` | miss 0.0072 m, vtd 0.9777 m/s, sat_frac 0.3506, mprop 854.98 kg |
| `tf_curve` rows | 14 |

**Consequences.** Host finding "the `.mat` is Phase-1 only" is **RESOLVED**
(it no longer is). The SDD's stale-`.mat` item (`§9` item 6) is **CLOSED**.
The note's "every number read from that file" claim is true again and was
kept — but scoped, because a handful of *mid-campaign sweep tables* still
come from `process/task-{5,7}-report.md` (see N-Repro below).

---

## 1. Note — CORRECTNESS tier

### N1 (GPT #1 + #2; the submission blocker) — APPLIED
Separate the losslessness **theorem** from the **implemented** program.

- Theorem 1 restated for the *exact exponential relaxation* (Eqs.
  `eq:linear`–`eq:relax`), with its assumptions now enumerated **in the
  theorem statement**: (A1) normality/controllability, which is what
  excludes the degenerate `λ_v ≡ 0` branch; (A2) a constraint qualification
  giving multipliers of bounded variation for the state path constraints;
  (A3) existence.
- A lead-in sentence says explicitly that the theorem is *not* a statement
  about the Taylor-restricted program that is solved.
- The proof is now framed as **intuition**, not a proof: "offered as
  intuition for why the relaxation is tight … it is not a self-contained
  proof, and the place it is loose is identified explicitly."
- New paragraph **"What the theorem covers, and what it does not"** draws
  the boundary in two bullets: (i) losslessness *transfers* to the Taylor
  program because the Taylor bounds alter only the σ-admissible set and
  never the pointwise minimization of `λ_v·u` over `‖u‖ ≤ σ`, and **G4
  certifies `‖u‖ = σ` for the program that was actually solved**; (ii) the
  Taylor restriction is a **separately measured model error**, reported as
  G3.
- Non-vacuity of Eq. (11)'s (the Taylor band's) upper bound added by measurement: `max δz =
  0.034` on the flagship Route-B solution (measured this session).

**Glideslope activity — checked in the `.mat` as instructed.** The margin
`z − ‖(x,y)‖ tan γ_gs` over the 61 nodes is 1705.6 m at t=0, falls
monotonically, reads **0.9027 m at the last interior node** (t = 16.319 s),
and is exactly 0 only at the pad, where *both sides vanish* because
`r(t_f) = 0`. So the cone is **inactive on the whole interior**; what binds
at the endpoint is the terminal boundary condition, not the cone. This is
now stated in the proof sketch where the "glideslope inactive" assumption is
taken (the reviews expected possible activity at t=0 — the measurement says
the opposite end, and degenerately).

### N2 (host; most important) — APPLIED
§5.4 reachability passage rewritten. The table is now `Table 3` with a
caption carrying **its configuration**: measured at the *then-shipped*
`ηT = 0.93`, vertical-only model, all thrust vertical, zero tracking error,
brake altitude by bisection; labelled a mid-campaign snapshot. Added:

- a caution that column 1 is a **bisection residual**, not a physical
  minimum (which is why it is non-monotone and why nothing should be
  inferred from its ordering);
- explicit disambiguation of the **three live altitudes** — 699.6 m (the
  model's own unit-thrust-scale row), 718 m (where the `ηT = 0.93`
  feed-forward actually braked), 914.7 m (where the shipped `ηT = 0.87`
  feed-forward brakes) — each tagged to its configuration;
- **the strongest-form inequality**, which the note previously stopped one
  step short of: the `−5%` vehicle needs 806.6 m; the shipped 0.87
  feed-forward brakes at **914.7 m > 806.6 m** (and within 11 m of even the
  `−10%` vehicle's 925.2 m). "The de-rate works because it moves the nominal
  brake altitude up past the dispersed vehicle's requirement."
- "107 m higher than nominal" **deleted** (it was a difference of two rows of
  the model's own column). Replaced by the correct 88.6 m = 806.6 − 718,
  both here and in the Future-work MPC item that quoted it.

### N3 (host; measured) — APPLIED
§5.2 eigenvalues corrected from `−0.0929 ± 0.0497i` / τ = 10.8 s to
**`−0.0905 ± 0.0499i` / τ = 11.0 s**, and `m_A` is now stated in the
sentence so the number is checkable. Measured this session:

```
m_A (mean mass over [0, tSwitch=6.1552 s]) = 29623.7 kg
2*m_A*sqrt(q_r*r)                          = 0.018736
roots of s^2 + (K_v/m)s + K_r/m            = -0.090478 +/- 0.049886i,  tau = 11.05 s
```

The paragraph is now internally consistent: the same `m_A` produces both the
0.0187 figure and the eigenvalue pair.

### N4 (host) — APPLIED
The de-rate margin comparison is restated in **one** currency (net
decelerating acceleration `Tmax/m − g0`), with the conversion shown:

```
m = 30000.0 kg:  5% loss = 7.67% of net | eta=0.93 reserve = 10.74% | eta=0.87 reserve = 19.94%
m = 26464.6 kg:  5% loss = 7.22% of net | eta=0.93 reserve = 10.10% | eta=0.87 reserve = 18.76%
```

So `ηT = 0.93` covers a −5% dispersion by ≈1.4× (thin) and `ηT = 0.87` by
≈2.5×. The `0.95/0.93 = 2.1%` figure two paragraphs later is now explicitly
labelled as the same comparison **in thrust units**, kept only because it is
how the 0.93 decision was recorded, with a cross-reference to the
net-deceleration form and to the 88.6 m form.

### N5 (GPT #4) — APPLIED
"Coast-extension" and "the free part of the deceleration" renamed and
corrected. The paragraph is now **"The minimum-throttle-arc extension
mechanism"** and opens by refusing the old framing: by Eq. (5) the engine is
never off, so the vehicle is on a continuous `Tmin` burn and is spending
propellant the whole time. What is free is **drag's incremental
deceleration**. Wording corrected throughout the paragraph and in the
Conclusion. (The two surviving occurrences of "free coast" in the note are
the *assertions that there is none*, and are correct.)

### N6 (GPT #5, both externals) — APPLIED, both docs
Note §Future-work and SDD §5.3 both made conditional: a finite cone can be
inactive and the primer law still holds while it is; **when the cone
constraint is active** the constrained optimum lies on its boundary and
thrust need not be parallel to the costate. Both docs now also state that
the *software* policy (never score the primer under a finite cone) is
deliberately more conservative than the mathematics, because nothing tells
the gate in advance which instants bind — and the note cites the measured
7.32° at θ_max = 15° as evidence the false-fail is real.

### N7 (Gemini + host) — APPLIED
The false "all velocity weights are untouched" claim is gone. The note now
matches `tvlqr_design`'s `QA`/`QB`: the **lateral** channel is rescheduled in
*both* position and velocity — the derived pair `(q_pos, q_vel)` occupies
entries (1,1),(2,2) and (4,4),(5,5) of `QA`, and `QB` replaces both with
1e-4 and 1e-2 — while the **vertical** channel `Q(3,3)`, `Q(6,6)` is
identical in `QA` and `QB`. Consistent with SDD §3.6 and the source.

### N8 (Gemini) — APPLIED
New `\section{Conclusion}` inserted before Future work (~1 page): the
certified min–max profile and why it follows from `Tmin > mg0`; the
cross-route agreement and gate margins; the tracking discoveries with their
common pattern (every fix followed from a derivation or a reachability
bound, and in every case a sweep had searched the wrong axis); the
counter-intuitive drag finding and its mechanism; and the methodological
residue (a gate can produce a false negative that looks like physics).

### N9 (GPT #6) — APPLIED
"Table 1 is the complete parameter set" scoped to the **physical and PDG**
parameter set, with explicit cross-references to the three other
configuration tables (dispersions/seed, tracker weights and schedule,
certification thresholds). A `\label{sec:battery}` was added to the
acceptance-battery subsection to support the cross-reference.

### N10 (GPT #3) — APPLIED
The G5 / Eq. (25) prose no longer says "forces exactly parallel" without
qualification. It now shows the collection step — every term other than
`λ_v` is radial in `T_{k+1/2}`, so the condition reads
`(4h/6) λ_v/m = c·T_{k+1/2}` — and names the two conditions that make the
collection possible and that both hold here: **midpoint states carry no path
constraint row** (glideslope, `z ≥ 0`, `m ≥ m_dry` are node-only) and **no
active pointing cone**. A finite active cone would add a non-radial term,
which is exactly why the gate stops scoring under a cone.

### N-Repro (host, resolved-then-scoped) — APPLIED
The Reproduction section now says the note is written against a
`phase2 = true` run (the shipped file, timestamped 2026-08-09), that every
result number — Phase 1 and Phase 2 alike — was read from it, and lists the
**exceptions**: the `ηT` cost table, the reachability table, the before/after
tracking pairs in §5.2–§5.3, and the Taylor-gap refinement sweeps, all
mid-campaign sweeps at superseded configurations recorded in
`process/task-{5,7}-report.md`.

---

## 2. Both docs — CONSISTENCY tier

### B1 (GPT #9 + host) — APPLIED, identically in both docs
The Taylor-gap band is now stated **with its conditions** in both documents.
Traced to `process/task-5-report.md`, which holds both halves:

- fixed-`t_f` sweep, `N_conv = 120/180/240`: 0.7369 → 0.7337 → 0.7310 kg
  (the note's old "0.70–0.74", correctly a **subset**);
- with the golden search allowed to reselect `t_f` at `tolTf = 0.01`:
  0.699 kg at `N_conv=120` and 0.943 kg at `N_conv=240`.

So **0.70–0.94 kg** is the envelope of every refinement tried (which is what
`certify_pdg.m`'s header documents and what the SDD said), and 0.73–0.74 kg
is the fixed-`t_f` subset. Both are now stated, with conditions, in both
docs, and both are labelled as mid-campaign measurements at the
pre-de-rate configuration. The flagship **0.428 kg** is tagged in both docs
as the same quantity measured at the shipped `ηT = 0.87`.

### B2 (GPT #7 + host-verified) — APPLIED
File counts measured: root 2, `lib/` 8, `certify/` 2, `viz/` 4 → **16
production modules**; `tests/` 11 → **27 `.m` files in all**. The SDD's
"Seventeen" sentence is corrected and now agrees with its own inventory
table (which was already right).

### B3 (GPT #8 + host-verified) — APPLIED, three places
Verified in source:

```
lib/pdg_dynamics.m:46   magEps = 1e-12;
lib/pdg_dynamics.m:48   Tmag = sqrt(sum(T.^2) + magEps);
lib/pdg_dynamics.m:53   vmag = sqrt(sum(v.^2) + magEps);
lib/solve_pdg_colloc.m:234  Tmag = sqrt(sum(T.^2) + 1e-12);
lib/solve_pdg_colloc.m:239  vmag = sqrt(sum(v.^2) + 1e-12);
```

Both paths use `sqrt(sum(.^2) + 1e-12)`. Corrected in all three places that
described it: the **note** §1.1 (was an unqualified `sqrt(Σ(·)²)`), the SDD
**§3.3** (now names `magEps`), and the SDD **Conventions** entry (which had
said MATLAB was unregularized and only CasADi regularized — the direct
contradiction GPT flagged). All three now state the shipped truth and say
the two implementations smooth identically.

### B4 (Gemini, SDD §3.3) — APPLIED
The false claim that `‖T‖ = 0` is reachable is removed. The two guards are
now given **different** justifications: for `v` it is **operational**
(`sim_closed_loop`'s arrest event is `v_z = 0` by construction), and for `T`
it is **formal smoothness only** — the control constraints hold
`‖T‖ ≥ Tmin = 338.0 kN` at all times, so a zero-thrust state is not
reachable. An explicit "*Do not* justify it by claiming a zero-thrust state
occurs" is left in place so the error is not reintroduced.

### B5 (host) — APPLIED
- SDD §6.2's stale primer sweep (`≈1.7°, 1.19°, 0.92°, 0.57°`) replaced by
  the same four values §5.3 and `certify_pdg.m` carry
  (`1.72°, 1.14°, 0.85°, 0.57°`), with a cross-reference so they cannot
  drift apart again.
- SDD §4.1 now records the shipped artifact explicitly: a completed
  `phase2 = true`, `doMC = true` run of 09-Aug-2026 12:56:25 carrying all
  thirteen variables, verified by `fieldnames`. The old (and wrong) variable
  list in §9 item 6 is gone with the ledger.

---

## 3. SDD structure (controller ruling: compress, don't delete)

**Ledger compression.** §9 is retained as a *section*, not as a per-item
changelog: one paragraph records that nine doc-vs-code discrepancies were
found while transcribing interfaces, that the code is what §3 documents,
that none changed a certified number, and that eight were fixed in code on
2026-08-09. The two that were substantive enough to change the software are
documented **where they now matter** (`certify_pdg`'s "never throws" →
gate G0 §5.6; `plot_pdg_solution`'s "max–min–max" label → the certified
min–max structure §2.3). The ninth (stale `.mat`) is recorded as closed. A
second paragraph points at the archive files, including this fix wave's own
review set. Section retitled "Code-versus-documentation discrepancies: where
the ledger lives".

This is the middle path between Gemini's "delete §8 entirely" and GPT's
"dedupe": the SDD stays as-built, the history stays reachable.

**Dedupe (GPT #12).** One authoritative telling per lesson, cross-references
elsewhere:

| lesson | authoritative telling | shortened to a cross-ref |
|---|---|---|
| pre-fix reconstruction inadmissible between samples (18% below `Tmin` over 11.7%) | §3.5 `hs_quad_ctrl` | §6.2 `N=60` paragraph |
| `annulus_switch` first-crossing limitation | §3.6 `tvlqr_design` | §8.2 |
| IPOPT on an SOCP / conic drop-in | §3.4 `solve_pdg_convex` | §8.5 bullet |

**As-built stamp.** Added as a boxed block under the title: commit
`77bf8e9` (2026-08-09), `booster_run.mat` of 09-Aug-2026 12:56:25, plus an
explicit **refresh trigger** naming the two things to re-check after any run
that replaces the `.mat`. This is the systemic fix the host review asked for
(§9 went stale silently because no trigger existed).

---

## 4. Additional host-measured corrections applied

These are one-line numeric fixes the host review measured; all re-verified
this session before applying.

| item | was | now |
|---|---|---|
| §2.1 IPOPT iterations | 30 | **32** (`solC.stats.iter`) |
| §2.3 gap as a mass fraction | `1.6e-5` of vehicle mass | **`1.4e-5` of `m0`** (0.428/30000) |
| §3 `arccos` floor | single value `1.207e-6` deg | **one-to-two-ulp range 8.5e-7 – 1.2e-6 deg**, with why the vacuum and drag columns differ by 1.4× |
| §5.2 lateral command | "840 kN … sits just outside the 845 kN rim" | **842 kN lateral, 907 kN commanded magnitude, 7% outside** |
| §6 run 84 draw | `Δr0 = [−34.3; 425.3; 121.6]` | **`[−34.2; 425.3; 121.6]`** |
| §6 thrust draw range | `±4.5%` | **`−4.4% / +4.5%`** (`[0.9561, 1.0447]`) |
| §5.2/§5.3 tables | unlabelled | labelled **mid-campaign snapshots**, superseded by Table 5 (the acceptance battery), with the one apparent self-contradiction (`Δr0 = [50;−30;0]`) explained as a configuration difference |

---

## 5. Deliberately NOT applied

| finding | why |
|---|---|
| GPT #10 — fully specify G3's `L∞` metric (components, interpolation, clock, sampling) | Outside the controller's apply list. It is a real gap in the SDD's interface definition and is worth a follow-up, but it needs a source read of `certify_pdg.m`'s G3 block and a decision about what to *promise*, not just a transcription. |
| Host — state Route B's discretization in §2.2 (trapezoid, `N_conv = 120`, `u`/σ as node variables) | Outside the apply list. It is the one place a reader cannot tell what was solved from the note alone; recommended for a follow-up. |
| Host — §4.2 "PMP mandates bang–bang with a bounded number of switches" is stronger than PMP gives | Outside the apply list. Still true as a criticism. |
| Host — `C_d·A` caveat on the 434.7 kg Phase-2 headline | Outside the apply list. |
| Host/GPT — SDD per-interface stability markers; a "how do I add a gate" walkthrough | Additive scope, not a correction. |

---

## 6. Verification

**Compilation.** Both documents compiled **twice** with
`/Library/TeX/texbin/pdflatex -interaction=nonstopmode`, exit 0, **zero `!`
errors**, **zero undefined references or citations** on the second pass
(the only residual warning is the harmless `OMS/cmtt/m/n` font-shape
substitution, pre-existing). Aux files cleaned
(`*.aux *.log *.out *.toc *.fls *.fdb_latexmk *.synctex.gz`).

| document | pages after fixes |
|---|---|
| `booster_landing_note.pdf` | **29** |
| `booster_landing_sdd.pdf` | **37** |

**verify-paper — note** (`--deep-figs`): **GREEN.**

```
CITATIONS  (12 cited, 12 in inline thebibliography)
  Summary: 10 verified, 1 to review, 0 failed, 1 skipped.
FIGURES  (4 \includegraphics)
  OK  pdg_solution.png / footprint.png / phase2_vac_vs_drag.png / phase2_footprint.png
DEEP FIGURE CHECK
  OK  no two figures share an identical plot interior
```

The single WARN is `betts2010` (Betts, *Practical Methods for Optimal
Control and Estimation Using Nonlinear Programming*, 2nd ed., SIAM, 2010) —
a monograph with no DOI, whose best Crossref hit is one of its own chapter
titles. Reviewed and correct as cited; pre-existing, untouched by this wave.
The SKIP is `casadi2019` (software paper). **0 FAILED on citations and 0 on
figures.**

**verify-paper — SDD**: **GREEN.** 0 citations, 0 figures (it is a TikZ-only
document), 0 failures.

**MATLAB probes.** Three read-only `matlab -batch` probes, run synchronously
against `results/booster_run.mat` and the source. Nothing in the repository
was modified by a probe.

---

## 7. Micro-pass: controller rulings on the four logged residuals

The four items logged in §5 as "deliberately NOT applied" were subsequently
authorized by the controller as a final micro-pass. All four are now applied.
Everything below was transcribed from source, not invented.

### R1 — GPT #10: G3's `L∞` metric. Ruling: **document as-implemented, no interface change.** APPLIED

`certify_pdg.m`'s G3 block reads, verbatim:

```matlab
tq  = linspace(0, min(solC.tf, solV.tf), 200);
rC  = interp1(solC.t.', solC.X(1:3,:).', tq.', 'pchip');
rV  = interp1(solV.t.', solV.X(1:3,:).', tq.', 'pchip');
rep.G3_traj_Linf = max(sqrt(sum((rC - rV).^2, 2)));
```

That is now transcribed into the SDD's certification section as a new
paragraph, **"G3's trajectory metric, exactly as implemented"**, carrying the
code block plus a seven-row table answering each of GPT's questions:

| aspect | as implemented |
|---|---|
| state components | **position only** (rows 1–3); velocity and mass do not enter |
| pointwise norm | Euclidean 3-vector norm of the position difference |
| reduction | `max` over samples — `L∞` **in time** of an `ℓ2`-in-space difference |
| comparison clock | **physical time in seconds** (not normalized time, not a shared node index) |
| sampling grid | **200** uniform samples on `[0, min(tf_C, tf_V)]` — a common window, so the longer solution's tail is never sampled |
| interpolation | **`pchip`** for both routes, each onto its own **node** grid (`solC.t` = N+1 nodes, `solV.t` = Nconv nodes); Route A's `Um` / `hs_quad_ctrl` are **not** used |
| terminal states | included only for whichever route owns `min(tf)`; the two terminal states are **not** compared to each other — `|Δtf|` is the row that scores endpoint disagreement |

Two consequences a maintainer needs are stated with it: the window
truncation means the row cannot see a disagreement confined to the last
`|Δtf|` of flight (3.4 ms on the flagship), and the 200-sample grid is fixed
and independent of `N`/`Nconv`, so the row deliberately does not sharpen
under refinement — it is dimensioned from the mission (5.0 m = one third of
`P.pad_radius`), not from the measurement. Measured value at the production
grid: 0.373 m. The gates longtable row now says "position-trajectory
difference … (defined precisely below)". **No code was changed.**

### R2 — Route B's discretization in note §2.2. APPLIED

New **Transcription** paragraph, transcribed from `solve_pdg_convex.m`
(`solve_fixed_tf`, lines ~200–255):

- `t = linspace(0, tf, Nc)` with `Nc = Nconv = 120` **uniform nodes**
  (hence 119 intervals), spacing `h`;
- **trapezoidal** on the linear dynamics — new displayed Eq. (12) giving the
  three defect rows for `r`, `v`, `z` in the paper's own notation, matching
  the source's `Rh/Vh/Z` updates;
- decision variables `(r, v, z, u, σ)` **at nodes only** — no midpoint
  variables, so `u` and `σ` are node values and the trapezoid rule treats
  them as piecewise linear;
- all convex constraints imposed **per node**: SOC `‖u_k‖ ≤ σ_k`, the Taylor
  band, `σ_k ≥ 0`, the `z` bracket, the glideslope in homogeneous quadratic
  form, `z_k ≥ 0`, and the pointing cone `u_z,k ≥ cos(θ_max)·σ_k` only when
  `θ_max` is finite;
- boundary conditions at first and last node; objective `min −z(tf)`;
  CasADi + **IPOPT**.

Two order facts are drawn out because §2.3 leans on them: Route B is
**second** order against Route A's third-order Hermite–Simpson, and its
control representation is **piecewise linear** against Route A's per-segment
quadratic — so the two routes share no discretization order and no control
basis, and the `Nconv` sweep is what rules out Route B's own trapezoid error.

The golden-section outer loop's **validity-code gating** is also now
described where the transcription is: codes 3/2 contribute `mf`, codes 1/0
contribute `−∞`, and the note records that the shipped `tf_curve` contains
exactly one such probe. Verified against the artifact:

```
tf_curve (14 rows): tf | mf | code
   14.5836          -Inf  0     <-- the only invalid probe
   15.6656    26037.9840  3
   ... (11 more, codes 3 and 2) ...
   19.1672    26249.2084  3
Nconv = 120 ; numel(solV.t) = 120
```

so the search optimizes a `−∞`-substituted surrogate, not `m(tf)` — which is
legitimate here only because the infeasible region is a **left** tail (too
little time to arrest), and the note now says exactly that.

### R3 — §4.2 "PMP mandates bang–bang with a bounded number of switches". APPLIED

Softened to what is actually true, in three separated claims:

- PMP gives the bang–bang **form** — because the Hamiltonian is linear in
  thrust magnitude, the optimal magnitude sits on a bound wherever the
  switching function is nonzero. **That is all it gives here.**
- It does **not** say which bound comes first (scenario-dependent — the
  content of §4.2).
- It does **not** bound the switch **count**: a count bound is an additional
  result requiring extra structure and is **not a theorem being invoked**.
  That the count is one on this problem family is an **empirical and
  structural observation** (Hamiltonian linear in the control, no singular
  arc detected at any grid) and is **measured, not derived** — G5 reports
  bound fraction 0.9917 and a single interior switch at every grid from
  `N = 15` to `N = 240`.

### R4 — `Cd·A` caveat on the 434.7 kg headline. APPLIED

One paragraph inserted in §7.1 immediately after the headline sentence,
where the number is first developed:

> the 434.7 kg is conditional on the drag model of Eq. (2) — `Cd = 1.0`,
> `A = 10.75 m²` (the body frontal area `π·1.85²`, with no base-drag,
> plume-interaction or grid-fin contribution and no Mach or attitude
> dependence), and an exponential atmosphere — and over this range the saving
> scales roughly with the product `Cd·A`. The **mechanism**
> (minimum-throttle-arc extension, switch moves later) is robust to that
> coefficient; the **number** is model-dependent, and a retro-propulsive base
> flow is exactly the regime in which an effective `Cd·A` is hardest to pin
> down.

This protects the result rather than weakening it: the finding that survives
is the mechanism, which is what §4 predicts independently.

### Micro-pass verification

Both documents recompiled **twice**, exit 0, **0 `!` errors**, **0 undefined
references or citations**. Aux cleaned.

| document | pages before micro-pass | pages after |
|---|---|---|
| `booster_landing_note.pdf` | 29 | **30** |
| `booster_landing_sdd.pdf` | 37 | **38** |

`verify-paper` re-run on both and **still green**: note (`--deep-figs`) 10
verified / 1 to review (the pre-existing `betts2010` no-DOI monograph) / **0
failed** / 1 skipped, 4 figures OK, deep-figs OK; SDD 0 citations, 0 figures,
**0 failed**.

One additional read-only `matlab -batch` probe was run to verify `tf_curve`'s
validity codes and `Nconv` before the Route-B paragraph asserted them.
Nothing in the repository was modified by a probe, and **no source file was
changed by this micro-pass** — all four items are documentation-only.
