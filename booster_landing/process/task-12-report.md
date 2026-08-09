# Task 12 Report: Theory note `doc/booster_landing_note.tex`

**Status: DONE.** 22-page LaTeX note, compiles clean, `verify-paper
--deep-figs` green (0 FAIL), committed as `fce8617`.

## Files changed

- **Created:** `booster_landing/doc/booster_landing_note.tex` (1426 lines,
  66 KB). The `doc/` directory itself did not exist and was created.
- **Generated, not committed:** `booster_landing/doc/booster_landing_note.pdf`
  (22 pp, 1.16 MB). `*.pdf` is git-ignored at the `optimal_control` repo
  root (checked: `optimal_control/.gitignore` line "# PDF course materials
  (reference only, not source)" / `*.pdf`), so the repo convention is
  source-only. Only the `.tex` was staged, per the brief's fallback rule.

No campaign code was touched. No `.md` report files other than this one.

## Section map

| § | Title | Content |
|---|---|---|
| — | Abstract | The whole campaign in one paragraph, with the six headline numbers |
| 1 | Problem | 3-DOF dynamics + drag term, thrust annulus, glideslope, BCs; full F9-class parameter table; the hoverslam T/W argument with the three downstream consequences named |
| 2 | Two solution routes | 2.1 free-$t_f$ Hermite--Simpson NLP + why nondimensionalize + the ND-dual caveat; 2.2 lossless convexification: change of variables, the relaxation, the Taylor mass bound (incl. *why* the max-thrust depletion reference), Theorem + PMP proof sketch, free-$t_f$ golden section; 2.3 cross-validation and what the residual **is** |
| 3 | Certification | G1/G2/G2ff/G3/G4/G5 architecture with the "defect is not accuracy" motivation and both G2 reconstruction stories; flagship gate table **verbatim**; Fig. `pdg_solution.png` |
| 4 | What the optimum actually is | min--max arc table, the physical argument, when max--min--max *would* be right, what the gate change did and did not weaken |
| 5 | Closing the loop: five things the tracking problem forced | 5.1 singular $v(t_f)=0$ BC; 5.2 the ARE pole ceiling + $Q(t)$ by ARE inversion; 5.3 altitude-indexed tracking; 5.4 the $\eta_T=0.87$ de-rate + reachability table; 5.5 acceptance battery 7/7 |
| 6 | Monte Carlo | dispersion table (with the thrust-scale vs Isp-scale distinction), results table, single-failure forensics, Fig. `footprint.png` |
| 7 | Phase 2 | drag model, vac-vs-drag table, the coast-extension mechanism, what certification can/cannot say under drag (G5 loosening + G1 tol), wind MC, Figs. `phase2_vac_vs_drag.png` + `phase2_footprint.png` |
| 8 | Future work | 6-DOF, full return profile, SCvx, conic solver, MPC/replanning (motivated by the de-rate's reachability argument), per-row-scaled G1 tolerance, the open drag-primer question |
| — | Reproduction | the one-line front-door command |
| — | Bibliography | 12 entries, all real |

The G2ff gate was added to the brief's five-gate outline because it exists
in the shipped code and is one of the campaign's genuine findings; the note
presents it as G2ff (a sub-gate of the continuous-accuracy family), matching
`certify_pdg.m` and `print_certify_report.m`.

## Provenance: where every number came from

All `.mat` numbers were extracted with a synchronous
`/Applications/MATLAB_R2025b.app/bin/matlab -batch` run against
`results/booster_run.mat` (saved 09-Aug-2026 09:50:28), via a scratchpad
script; nothing was typed from memory or from a report where the `.mat`
carries it.

### From `results/booster_run.mat` directly

| Quantity | Field | Value used |
|---|---|---|
| Parameter table | `P.*` | all rows of Table 1 |
| Colloc $t_f$, $m_f$ | `solC.tf`, `solC.mf` | 16.595192 s, 26464.589 kg |
| Convex $t_f$, $m_f$ | `solV.tf`, `solV.mf` | 16.598601 s, 26464.161 kg |
| Golden-section probe count | `size(solV.tf_curve,1)` | 14 |
| Flagship gate table | `rep` via `print_certify_report(rep)` | verbatim block in §3 |
| Phase-2 gate table | `repD` via `print_certify_report(repD)` | §7 table |
| Drag solve | `solD.tf`, `solD.mf` | 17.530799 s, 26899.276 kg |
| Nominal closed loop | `out0.td.*`, `out0.sat_frac` | miss 0.0072 m, vtd 0.9777, alt $-1.1\times10^{-15}$, landed, sat_frac 0.3506, 854.98 kg remaining |
| Vacuum MC | `mc.*` | 99.5%, 199/1/0, miss mean 3.721 / med 3.665 / p95 6.723 / max 9.238, vtd mean 1.054 / med 1.000 / p95 1.206 / max 8.290, mprop 0.0/537.26/882.66 |
| Vacuum MC failure forensics | `mc.dr0(84,:)`, `mc.dv0(84,:)`, `mc.thrust_scale(84)`, `mc.stop{84}` | $[-34.25;425.32;121.57]$ m, $[2.09;19.55;-7.37]$ m/s, 0.9805, `arrest` |
| Wind MC | `mcD.*` | 99.0%, 200/0/0, miss mean 6.021 / med 5.623 / p95 12.381 / max 15.517, vtd mean 0.919 / max 1.125, mprop min 477.26 / mean 994.12, max wind 34.72 m/s; the two failures are misses at 15.517 m and 15.145 m, both `touchdown` |
| TVLQR | `ctrl.tSwitch`, `ctrl.tBlend`, `diag(ctrl.Qf)` | 6.1552 s, 0.4 s, $Q_f(6,6)=9.2906$ |

### Computed in MATLAB from `.mat` fields (script in scratchpad)

- $m_0g_0 = 294.20$ kN, $m_{dry}g_0 = 251.05$ kN, $\Tmin/(m_0g_0)=1.1489$,
  $\Tmin/(m_{dry}g_0)=1.3463$, $\Tmax/(m_0g_0)=2.8722$, $I_{sp}g_0 = 2765.5$ m/s,
  $|\rvec_0| = 2063.98$ m, $|\vvec_0| = 182.48$ m/s.
- Fuel: vac 3535.411 kg, drag 3100.724 kg, $\Delta = 434.687$ kg (12.30% of
  vacuum fuel); $\Delta t_f = +0.93561$ s.
- Switch times by mid-annulus crossing of $\|T\|$ on the node grid with
  linear interpolation (the same detector `tvlqr_design`'s `annulus_switch`
  uses; cross-checked — it reproduces `ctrl.tSwitch = 6.1552` exactly for
  the vacuum case): vac 6.1552 s / 914.65 m / 175.44 m/s; drag 10.9567 s /
  378.85 m / 119.11 m/s. Brake arcs 10.4399 s (62.9%) and 6.5741 s (37.5%).
  **The method is footnoted in the note**, because the README quotes the
  drag switch as "10.8 s" from a different detector; 10.957 is what this
  measurement gives and it is the number in the note.
- Net accel at $\Tmin$: 1.460 m/s² at $m_0$, 3.397 m/s² at $m_{dry}$.
- G4 threshold $10^{-4}\eta_T\Tmax/m_0 = 2.4505\times10^{-3}$ m/s² (margin 18x);
  G1 margin 7.2x; G3 uses 57% of ceiling; bound fraction $120/121 = 0.99174$.
- 4.25$\sigma$ attribution for MC run 84 (425.3 m against a 100 m 1$\sigma$).

### From task reports (measurements not stored in the `.mat`)

| Claim | Source |
|---|---|
| pchip vs per-segment-quadratic G2 collapse ($0.7897 \to 1.30\times10^{-5}$ kg, 60,700x; holds at $N=240$); the $N=15$--$240$ plateau | task-5 report, Fix Report "Important 4+8" |
| $N_{conv}$ 120/180/240/300 refinement sweep holding 0.70--0.74 kg; tolTf 0.05→0.01 no help | task-5 report, "Important 7" |
| G1 defect / annulus / `onHi(end)` behaviour across $N=15$--240 | task-5 report N-sweep table |
| The $(\eta_T,N)$ G2-position grid lottery ($7.33\times10^{-5}$ to 1.287 m → 0.002--0.051 m after the switch-step fix); Eq. (switch step) | task-7 report §7b-1 §2 |
| G2ff root cause: 18% below $\Tmin$ over 11.7% of flight, +0.6--1 m open-loop bias | `certify_pdg.m` G2b comment block + task-7 round 4 |
| Singular BC: 0.53 m arrest, 30-config sweep, retraction of the "removes the arrest entirely" claim | task-7 rounds 3--4 + `booster_params.m` `P.vf` comment |
| ARE pole ceiling: Eqs. for $K_r/m$, $K_v/m$, pole formula, $\sqrt{q_r/q_v}=0.1$ rad/s at old weights, 0.059 rad/s at shipped $r$; 18 kN asked vs 507 kN spare; 32 m at switch, 10.9 m at touchdown, 16° tilt; before/after 10.851→0.944 m and 6.819→1.363 m/s, sat_frac 0.653→0.13 | `tvlqr_design.m` ADAPTATION 5 (verbatim measurements) + task-7 round 5 |
| ARE inversion $q_r = m^2r\omega^4$, $q_v = m^2r\omega^2(4\zeta^2-2)$, $\zeta\ge1/\sqrt2$; $\omega=0.70$, $\zeta=1.00$ → $q_r=0.2099$, $q_v=0.8569$; validation $K_{11}/m=0.4808$ vs 0.49, $K_{14}/m=1.377$ vs 1.40; $\omega$ feasibility ceiling 28 m/s² ≈ 840 kN | `tvlqr_design.m` ADAPTATION 5(a) |
| Altitude indexing: $-169.38$ vs $-169.84$ m/s at $t=6.76$; 50 m low → $-31.8$ m/s; 9.0→43.4 m/s under tightening; $v^*(1950)=-179.83$ vs $-180.00$; the three structural deletions | `sim_closed_loop.m` ALTITUDE-INDEXED note |
| De-rate: 67.4 m/s under two controllers; reachability table (1.00/0.98/0.95/0.90 → 699.6/741.1/806.6/925.2 m); bandwidth sweep 5→20 rad/s gives 3.91→2.66 m/s; 0.93 gives 2.1% margin and 5/7 | task-7 §5 + §7b + `booster_params.m` `P.etaT` comment |
| De-rate cost table ($\eta_T$ 1.00/0.93/0.87 → 3453.108/3493.075/3535.411 kg) | task-7 §7b FINAL §4 |
| Acceptance battery 7/7 final (post-polish) numbers | task-7 §7b polish "Final battery (post-I1)" |
| MC dispersion 1σ magnitudes and the thrust-scale vs Isp-scale distinction | `run_monte_carlo.m` + `sim_closed_loop.m` headers |
| Drag G5 primer: 1.14° vacuum vs 5.11° drag at coarse grid, 2.605° at production; the $B$-independence argument | `certify_pdg.m` G5-drag note |
| Drag G1 tol: $3.21\times10^{-6}$ FAIL → $8.87\times10^{-7}$ → $2.28\times10^{-8}$; $1.07\times10^{-10}\times M_c$ five-figure match; $m_f$ unchanged to 0.001 kg | `solve_pdg_colloc.m` "G1 tol under drag" note |
| Independent two-arc reproduction of the drag result to <1% | SDD ledger, Task 11 review line |
| Nondimensionalization necessity (divergence at maxIter 6000); convex-side blown gap up to ~500; 52 kg cost of status-string-first classification; $t_f\approx27$ s infeasibility wall and the 6.5e-5/0.448 nondeterminism | `solve_pdg_colloc.m` / `solve_pdg_convex.m` / `booster_params.m` `P.tf_hi` headers |

### Derived analytically in the note (my own work, checkable)

- The Theorem/proof sketch (costate $\lamr$ constant, $\lamv$ affine, ball
  minimizer, exclusion of $\lamv\equiv0$ by transversality + nontriviality).
- **The sign prediction for G3.** The note argues that because
  $z\ge z_0 \Rightarrow \delta z\ge0$, the Taylor pair
  $1-\delta z+\tfrac12\delta z^2 \ge e^{-\delta z} \ge 1-\delta z$ makes
  \eqref{eq:taylor} a *conservative inner approximation on both sides*, so
  the convex route must report $m_f$ **less than or equal to** the
  collocation route's. Measured: $26464.161 < 26464.589$. This is new
  framing (no report states it) but it is a two-line calculus check and it
  strengthens the "the residual IS the model error" claim from a magnitude
  argument to a magnitude-**and-sign** argument.
- Bound fraction $0.9917 = 120/121$ ($2N+1$ samples, one midpoint inside the
  annulus at the switch).
- The gravity-loss framing of why min--max beats max--min--max, and why drag
  moving the switch later corroborates it.

## Compile output

```
$ /Library/TeX/texbin/pdflatex -interaction=nonstopmode booster_landing_note.tex   (x2)
exit 0
Output written on booster_landing_note.pdf (22 pages, 1164238 bytes).
Overfull \hbox (2.27057pt too wide) in paragraph at lines 834--837
<no LaTeX Warnings, no undefined references, no undefined citations>
Figures resolved: ../results/pdg_solution.png, ../results/footprint.png,
                  ../results/phase2_vac_vs_drag.png, ../results/phase2_footprint.png
```

One residual 2.3 pt overfull box (a math-heavy line in §5.2) — cosmetically
invisible; the three larger ones found on the first compile (18.7, 26.9,
38.5, 72.8 pt) were all fixed by rewording/`sloppypar`, not suppressed.

Aux files cleaned after the final compile (house rule):
`rm -f *.aux *.log *.out *.toc *.lof *.lot *.fls *.fdb_latexmk *.synctex.gz`.
`doc/` now contains only the `.tex` and the `.pdf`.

## verify-paper output (mandatory house gate)

```
$ ~/ai_council/venv/bin/python ~/Documents/myLatex/tools/paper_verify.py \
    booster_landing/doc/booster_landing_note.tex --deep-figs

CITATIONS  (12 cited, 12 in inline thebibliography)
   OK   acikmese2007 / acikmese2011 / anderson1990 / blackmore2010 / boyd2004 /
        ecos2013 / ipopt2006 / kelly2017 / lawden1963 / mao2016   (Crossref title match)
  WARN  betts2010: no DOI; best Crossref match (score 54) — review
  SKIP  casadi2019: non-DOI resource — verify by hand
  Summary: 10 verified, 1 to review, 0 failed, 1 skipped.

FIGURES  (4 \includegraphics)
   OK   pdg_solution.png / footprint.png / phase2_vac_vs_drag.png / phase2_footprint.png

DEEP FIGURE CHECK
   OK   no two figures share an identical plot interior
```

**Verdict: green — 0 FAIL on citations, 0 FAIL on figures.** Both WARNs
reviewed:

- `betts2010` — Betts, J.T., *Practical Methods for Optimal Control and
  Estimation Using Nonlinear Programming*, 2nd ed., SIAM, 2010
  (ISBN 978-0-898716-88-7). Real; Crossref indexes SIAM monographs by
  chapter, which is why the best match is "1. Introduction to Nonlinear
  Programming" from that same book. Kept.
- `casadi2019` — Andersson et al., *Math. Prog. Computation* 11(1), 2019,
  pp. 1--36. Real, the canonical CasADi paper; skipped by the tool as a
  software resource. Kept.

The first verify run had two additional WARNs (`blackmore2010` and
`boyd2004` present but never cited). Rather than delete two real and
relevant references, both were given genuine citation sites in §2.2
(Blackmore et al. as the minimum-landing-error extension of the
Acikmese--Ploen construction; Boyd & Vandenberghe for local-equals-global
under convexity). Re-verified clean.

**No invented citations.** Every entry is a paper or book I can identify
by author, venue, year and page range; the brief's four required references
are all present and cited.

## Commit

```
fce8617  booster_landing: theory note (PDG, convexification, gates, TVLQR, MC, drag)
```
Staged: `booster_landing/doc/booster_landing_note.tex` only (PDF is
git-ignored repo-wide).

## Self-review

**What I checked before writing a number down.** I loaded the `.mat` and
printed the gate reports through the campaign's own
`print_certify_report`, rather than transcribing tables out of the reports
— this caught that the flagship `rep` in the shipped `.mat` differs from
several tables in the task reports (e.g. G2 pos is 0.00885 m, not the
0.000156 m of the pre-de-rate task-5 tables, and G3 `|dmf|` is 0.4279 kg,
not 0.7039). The note quotes the shipped `.mat` throughout and quotes
report numbers only for measurements the `.mat` does not carry (sweeps,
batteries, before/after comparisons), each flagged as such in the
provenance table above.

**Deviations from the brief's outline, all deliberate:**
1. The brief lists five gates; the shipped code has six (G2ff). Documented
   as G2ff, with the bug it was built to catch.
2. §5 is longer than the brief's "Tracking" section, and is written as
   engineering narrative with the measurements inline — per the dispatch's
   instruction that this is the note's most original content.
3. A short "Reproduction" section was added at the end; a note whose whole
   claim is "every number traces to one file" should say how to regenerate
   the file.
4. The de-rate's cost is quoted as 1.87% of *usable propellant* (82.30 /
   4400) and 2.38% of *fuel burned* (82.30 / 3453.108), both stated, because
   the reports use the two denominators in different places and the
   ambiguity is worth closing rather than inheriting.

**Things I deliberately did not claim:**
- I did not say drag "makes the problem easier" or that Phase 2 is
  certified equivalently — §7.2 states plainly that G3/G4 are skipped
  because there is no convex twin, and that reporting a G3 number under
  drag would be reporting a fiction.
- I did not present the G5 max-first relaxation as costless; the note says
  which half was retained and gives the evidence that the retained half is
  non-trivial (`onHi(end)` true at every $N$ from 15 to 240).
- I did not claim the MC vtd max of 8.29 m/s is a touchdown speed; the note
  explicitly attributes it to the single arrest run.
- The drag primer degradation is left as an open question with the analytic
  argument for why it *should not* happen, rather than force-closed.

**MATLAB/house conventions:** no emoji; LaTeX for the technical writing;
`pdflatex` twice then aux cleaned; `verify-paper --deep-figs` run and green
before declaring done; commit prefixed `booster_landing:` with the
Co-Authored-By footer; `.pdf` not staged (repo-wide ignore).

## Concerns for the controller

1. **The figures are git-ignored.** `booster_landing/.gitignore` excludes
   `results/*`, so the four PNGs the note `\includegraphics` are *not* in
   the repository. The `.tex` compiles for anyone who has run
   `run_booster_landing(struct('phase2', true))` and fails for anyone who
   has not. That matches the brief's instruction and the campaign's
   "products are regenerable" convention, but if the note is ever meant to
   travel (a collaborator, a paper submission, an arXiv drop), the figures
   will need to be either committed under `doc/figs/` or regenerated at the
   destination. Worth a decision, not a fix I made unilaterally.
2. **One switch-time number disagrees with the README.** The README says the
   drag coast extends "6.1 -> 10.8 s"; my measurement (mid-annulus crossing
   on the node grid, the same detector the tracker uses, cross-validated
   against `ctrl.tSwitch` for the vacuum case) gives 10.957 s. The note uses
   10.957 and footnotes the detector. If the README's 10.8 came from a
   different, preferred detector, the two should be reconciled — but the
   mechanism and its magnitude are unaffected either way.
3. **`P.theta_max_deg = Inf`.** The thrust-pointing cone is implemented but
   disabled for the whole campaign. The note says so in §1 and lists
   enabling it as the first step of the 6-DOF work item; nobody should read
   the certified result as covering a pointing-constrained vehicle.
4. **Offer standing:** the spec-time decision was to run a `doc-review` pass
   on the note rather than on the spec. I have not run it — that is Mike's
   call. The note is at a state where a three-way review (host Claude +
   GPT-5.6 + Gemini 3.1 Pro) would be productive: the proof sketch in §2.2,
   the sign argument in §2.3, and the ARE derivation in §5.2 are the three
   places where an independent math check would be worth the most.

---

# Fix Report (Review Round 1)

Review verdict: **Approved with fixes** — 1 Critical, 5 Important, 12 minors.
All addressed. Recompiled twice, aux cleaned, `verify-paper --deep-figs`
re-run and still green. **The note is now fully clean: 23 pages, zero
overfull boxes, zero LaTeX warnings** (the previous round shipped with one
2.3 pt box; that is gone too).

## Critical 1 — false claim against the .mat (§6). FIXED.

The reviewer is right and this was my error. I wrote "among the 199 genuine
touchdowns the worst lies inside the p95 figure." Re-measured:

```
p95 over ALL 200 runs       = 1.2055 m/s
worst vtd among the 199 landed = 1.4760 m/s
n landed exceeding p95(all)  = 9 of 199
landed vtd span              = 0.879 .. 1.476 m/s
```

So the claim was false twice over: the worst genuine touchdown is *outside*
p95, and nine runs are. Root cause of my error: p95 is computed over all 200
runs including the arrest, and I reasoned about the population without
re-slicing it by `mc.landed` — exactly the arrest-contaminates-the-statistic
trap this campaign built the `landed` classifier to avoid, which makes it a
worse slip than a typo. Replaced with the measured statement (worst genuine
touchdown $1.476$ m/s, population span $0.879$--$1.476$ m/s, comfortably
inside the $2.0$ m/s gate), keeping the correct point that the $8.29$ m/s
maximum belongs to the arrest run.

## Important 2 — the sign claim was an over-deduction (§2.2, §2.3, abstract). FIXED.

Accepted in full. $m_V \le m_{true}$ (restriction) and $m_C \approx
m_{true}$ (discretized *local* optimum) are two statements about the true
optimum; they do not order each other, and §4 of the note explicitly
declines to claim Route A finds the global optimum. Rewritten in three
places:
- §2.2 now says the restriction "does not order the two *computed* answers"
  but "gives a direction the measured residual can either corroborate or
  contradict."
- §2.3 says the measurement "is consistent with that and corroborates it,"
  then states the two extra assumptions a rigorous ordering would need
  (Route A at the continuous-time global optimum, plus a bound on its own
  discretization error) and cross-refs §4 for why the first is not claimed.
- Abstract: "shown to be, in both magnitude and *sign*, the model error" →
  "identified as the model error ... it does not shrink under refinement,
  and its sign is the one that bound's conservatism points to."

## Important 3 — the $\lamv \equiv 0$ exclusion (§2.2). FIXED.

Also accepted. With $\rvec(t_f)$ and $\vvec(t_f)$ fully constrained,
transversality leaves $\lamr(t_f), \lamv(t_f)$ free and nontriviality is
satisfied by $(\lambda_z, \lambda_0)$ alone — so my sentence proved nothing.
The proof sketch now says so explicitly and attributes the closure to the
controllability/normality assumption in Acikmese & Ploen and Acikmese &
Blackmore, noting that their theorem "is what makes the statement above a
theorem rather than a heuristic." Also flagged $\dot\lambda_z = 0$ as a
sketch simplification (unconstrained-path case; an active glideslope or
throttle-bound multiplier adds a term), with a note that it does not touch
the argument, which turns only on $\lamv$.

## Important 4 — convexity mis-attributed (§2.2). FIXED.

Correct and it does strengthen the section. Added a paragraph before the
Taylor construction stating that $\Tmin e^{-z} \le \sigma$ is *already*
convex (sublevel set of a convex function) and that only the upper bound is
nonconvex; the construction therefore does two different jobs — convexifies
the upper bound, and on the lower one trades an exactly-convex exponential
for a quadratic, buying SOC/quadratic representability and a solver-friendly
Hessian, not convexity. The follow-on sentence was corrected to match
("upper bound is now affine ... lower bound is convex quadratic (it was
already convex; it is now representable)").

## Important 5 — primer sign convention inconsistent. FIXED.

Added a dedicated "A sign convention, stated once" paragraph after
Eq. (primer): Eq. (primer) carries the minus sign of the PMP *minimum*
principle; G5 measures against $+\lamv$ because it reads the NLP *defect
duals*, whose sign is set by how the constraint was written. Records that
the literal form measured $179.7$--$179.9^\circ$ at every grid (a convention
flip, not a misaligned solution), that the sign was flipped once at
extraction, and that no `abs()` was introduced, so a genuine $180^\circ$
regression is still caught. Explicitly maps the $+\lamv$ statements in §3
and §7.2 onto it.

## Important 6 — "Five Certification Gates" vs six. FIXED (recommended option).

Title kept. §3's opening now reads "five gates, of which the second
(continuous accuracy) has two halves: G2 asks whether the reconstructed
control *flies* to the right place, and G2ff asks whether that same
reconstruction is *admissible* while doing it. They share a reconstruction
and a failure history, which is why they are one gate and not two." The
G2ff paragraph heading became "G2ff — feed-forward annulus feasibility
(second half of G2)."

## Minors — all folded in

| # | Item | Resolution |
|---|---|---|
| 1 | §5.2 "the actual pole" | Eq. relabelled $\omega_r \equiv K_r/K_v$ and described as the overdamped-limit corner frequency / dominant-pole surrogate. Added the measured exact eigenvalues: **$-0.0929 \pm 0.0497\mathrm{i}$** at the mean mass over the margin arc (I computed these rather than take the review's $-0.0906\pm0.0498$i; the small difference is the mass point — at $m_0$ they are $-0.0897\pm0.0499$i, at $m_{dry}$ $-0.1000\pm0.0486$i). Envelope time constant $10.8$ s on a $16.6$ s flight. Bound argument unaffected and now explicitly carries the conclusion. |
| 2 | §7.3 crosswind unattributed | Now "the worst miss (run 43, $15.52$ m) drew $\mathbf{w} = [25.74; -9.41; 0]$ m/s, a $27.4$ m/s steady wind." Measured from `mcD.wind(43,:)`. |
| 3 | §5 "each fix was derived rather than tuned" | Softened to the *need* for each fix and its form being derived, with an explicit concession that some constants were still swept for, "$\eta_T = 0.87$ most obviously." |
| 4 | §2.1 "separated (compressed)" | Now "the *compressed* form ... (compressed and separated are the two alternative formulations there; the midpoint state is eliminated here rather than carried as a decision variable)." |
| 5 | §5.4 "~2% net-decel to spare" | **Could not be traced to the .mat** — my own computation gives different framings (a $-5\%$ thrust costs $7.2$–$7.7\%$ of $\Tmax/m - g_0$ on the braking arc; $\eta_T=0.87$ reserves $23$–$25\%$). The figure came from a `booster_params.m` comment I could not reproduce. **Deleted** and replaced with two statements I can source: at $\eta_T=1$ the guidance ceiling *is* the engine ceiling so the braking arc has no reserve by construction, and the measured $7.2$–$7.7\%$ deceleration cost of a $-5\%$ engine. Also removed the resulting duplicate of the $0.95/0.93 = 2.1\%$ figure (it already appears in its proper place). |
| 6 | 4 ASCII quote pairs | Fixed to ``…'' at all four prose sites (lines ~417/745/880/908) plus two more introduced by the Important-5 paragraph. The three remaining `"` are inside `verbatim` (shell command) and in `\"a`/`\"u` umlaut escapes — correctly left alone. |
| 7 | In-text `\ref`s | Was 3 of 12 referenced; now **12 of 12**. Added natural-position references to `tab:params`, `fig:solution`, `tab:derate`, `tab:battery`, `tab:mc`, `tab:phase2`, `fig:phase2`, `tab:mcd`, `fig:phase2footprint`. |
| 8 | Eq. (2) missing $\vvec_{rel}$ | Now $\mathbf{a}_D = -\tfrac12 \rho_0 e^{-z/H} \frac{C_dA}{m}\lVert\vvec_{rel}\rVert\vvec_{rel}$ with $\vvec_{rel} = \vvec - \mathbf{w}$, plus a sentence that the guidance solves with $\mathbf{w}=\mathbf{0}$ and wind enters only in the truth model — "the guidance does not know the wind." |
| 9 | §5 forward-pointer to §3 | Added to the §5 preamble. |
| 10 | Note portability | Four PNGs copied to `doc/figs/` and **committed**; `\graphicspath{{figs/}{../results/}}` so a clean clone compiles from the committed copies while anyone with a fresh campaign run still picks up `../results/` second. Verified from the log: figures now resolve as `./figs/*.png`. No `.gitignore` change needed — `results/*` is ignored, `doc/figs/` is not (`git check-ignore` confirmed). 744 KB total. |

## Compile + verify (post-fix)

```
$ pdflatex -interaction=nonstopmode booster_landing_note.tex   (x2)
exit 0
Output written on booster_landing_note.pdf (23 pages, 1170780 bytes).
Overfull \hbox : (none)
LaTeX Warning : (none)
Figures: ./figs/pdg_solution.png  ./figs/footprint.png
         ./figs/phase2_vac_vs_drag.png  ./figs/phase2_footprint.png
```
Aux cleaned. `doc/` holds the `.tex`, the `.pdf` (git-ignored) and `figs/`.

```
$ paper_verify.py booster_landing/doc/booster_landing_note.tex --deep-figs
CITATIONS (12 cited, 12 in inline thebibliography)
  Summary: 10 verified, 1 to review, 0 failed, 1 skipped.
FIGURES (4 \includegraphics)  — all 4 OK
DEEP FIGURE CHECK — OK, no two figures share an identical plot interior
```
**Still green: 0 FAIL on cites, 0 FAIL on figures.** The one WARN
(`betts2010`, a SIAM monograph Crossref indexes by chapter) and one SKIP
(`casadi2019`, software resource) are the same two reviewed and kept last
round.

## Self-review of this round

- Critical 1 was a real false claim about data I had in hand, not a
  reporting slip. I re-derived the corrected numbers from `mc.landed` rather
  than trusting the review's figures, and they matched.
- Important 2 and 3 were both cases of me writing a *stronger* argument than
  the evidence supports — a restriction argument dressed as an ordering, and
  a hand-wave dressed as a transversality proof. Both are now weakened to
  what is actually true and, in the case of Important 3, attributed to the
  literature that does the real work. I consider these the most valuable
  findings in the review.
- Minor 5 is the one item where I deviated from "fix it" toward "delete it":
  the $\sim2\%$ figure could not be traced, and the note's whole contract is
  that every number traces. I replaced it with two sourced statements rather
  than keeping an untraceable one with a hedge.
- No numbers were changed anywhere else; the flagship results, gate table
  and battery are untouched.

## Commit (this round)

```
booster_landing: theory-note review fixes -- touchdown-population claim, Taylor sign framing, PMP sketch, G2ff as G2's second half
```
Files: `booster_landing/doc/booster_landing_note.tex`,
`booster_landing/doc/figs/{pdg_solution,footprint,phase2_vac_vs_drag,phase2_footprint}.png` (new).
