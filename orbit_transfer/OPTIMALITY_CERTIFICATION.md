# Optimality certification — cross-campaign register

**Purpose.** One place that answers "what do we actually certify, on which
transfer, and what is still missing" — for **both** first- and second-order
optimality (goal 1 of the standing program).

**Why one file and not two.** The two orders are not independent deliverables.
The Maurer–Osmolovskii sufficient condition for a bang-bang local minimum is
assembled from three ingredients that straddle the boundary:

    (i)   projected switching-time Hessian positive definite      [2nd order]
    (ii)  regular switching, Sdot != 0 at every switch            [1st order]
    (iii) primer / strengthened Legendre-Clebsch direction optimality  [1st order]

Splitting the record would scatter one certificate's ingredients across two
documents — the exact fragmentation that let one campaign re-plan work another
had already built and killed. Part A is a **coverage matrix** (first order is
largely solved; what matters is *where it runs*). Part B is a **research
register** (second order is open; what matters is *what broke and why*).

**Scope discipline — read this before adding to it.**
- Function headers stay authoritative for *how an instrument works* and for
  its own derivation. Do not restate derivations here; link them.
- This file owns what no header can: the **status matrix**, the **blocking
  mechanisms** (which are general, not per-file), what is **ruled out**, open
  **leads**, and the **decision record**.
- It is a register of experiments. It is only worth keeping if it is updated
  when one runs. Append to §6 every time.

Why it exists at all: in July 2026 the CR3BP campaign's TODO independently
planned a "tier 1 / tier 2 / tier 3" second-order program that the tulip
campaign had already built and found blocked. That is the same failure mode
diagnosed in `earth_elliptic_to_geo/process/LESSONS_DUAL_EXTRACTION.md` (L1:
the sibling implementation was invisible).

---

# PART A — First order (PMP extremality)

**State: largely solved.** The instruments exist and pass; the gaps are
*coverage* and one *structural* caveat, not correctness. (They were not
trustworthy before 2026-07-25 — see
`earth_elliptic_to_geo/process/LESSONS_DUAL_EXTRACTION.md` for why, which is
not restated here.)

## A1. Instruments

| # | condition tested | implementation |
|---|---|---|
| 1 | Hamiltonian minimization over thrust direction (β ∥ −primer) | `mee_primer_switch` + `verify_pmp_mee`; tulip `certify_minfuel_pmp` |
| 2 | switching-function sign law (S<0 ⇔ thr=1) | same |
| 3 | switch-time alignment (S zero-crossings vs actual throttle switches) | same |
| 4 | transversality λ_m(t_f)=0 (free final mass) | same, relative gate |
| 5 | Hamiltonian constancy (λ_t const ⇔ H_t conserved; autonomous only) | `hamiltonian_const_check`, `hamiltonian_along_traj` |
| 6 | full KKT certificate — stationarity, dual feasibility, complementarity | `sosc_recover_kkt` + `sosc_kkt_residual`; `results/dual_anomaly/diag_t1_beta.m` |
| 7 | no singular arc (no run of \|S\|≈0) | `verify_pmp_mee` H2 check |
| 8 | **regular switching, Sdot != 0** — *Maurer ingredient (ii)* | `second_order_nlp.m:278` (see gap G3) |
| 9 | **generic AD-based full first-order gate — one instrument, all four campaigns**: KKT stationarity, min-condition tangential residual, sign law, transversality, singular-arc count, regular-switching Sdot != 0, all from the NLP's own `opti.lam_g` via manifest-driven AD Lagrangian, printed+saved by a fixed-format standard report | `verify_common/foc_check.m` + `foc_report.m` + `foc_manifest.m` (built 2026-07-25, Tasks 1-5); wired into all four campaigns (Tasks 6-9): `run_foc_mee.m` (earth 2-body), `verify_cr3bp_pmp.m` (CR3BP), `run_foc_tulip.m` (tulip), `run_foc_elfo.m` (ELFO) |
| 10 | **continuous-time accuracy gate (G1) and method-agreement gate (G2)** — re-integrate each interval at 1e-11/1e-13 from its own left node and compare at the right node; POSITION and VELOCITY errors resolved separately and each given its own dimensional tolerance (km, m/s) — do NOT norm the full state and call it km, which is the bug external review caught here. G1b additionally flies the control end-to-end for a global figure. G2 compares t_f against an independently converged indirect solution where one exists | `DRO_tulip/direct/certify/certify_dro_mintime.m` + `dro_residual.m` + `tests/test_dro_residual.m` (built 2026-08-02, DRO→tulip only; **2026-08-25: engine generalized to `oc.local_residual`** — DRO rerouted with 1.7e-18 equivalence, earth-MEE gate `mee_residual` added; CR3BP-GEO + ELFO wrappers still open) |
| — | primal, dual-free cross-checks (switch structure, edge fraction, bound saturation) | `switch_structure`, solver `boundSat` |

**Note — two different costate sources.** The earth/CR3BP path reads the NLP's
own multipliers (`opti.lam_g`). The tulip path (`certify_minfuel_pmp`) instead
*reconstructs* costates by sparse least-squares on the discrete adjoint
recursion plus primer directions, with the scale pinned by S=0 at switches.
The second is methodologically stronger (it does not consume the solver's
duals) and is the partial answer to gap G1.

**2026-07-25 update.** `foc_check` generalizes the first (raw-dual) approach
into one manifest-driven instrument that now runs on all four campaigns,
including tulip and ELFO. On the tulip flagship the two sources were run
side by side on the *same* re-solved row (Task 8): raw-dual `foc_check` gives
primerAlignDeg = 0.058°, signPct = 100.00% (PASS); the independent
LS-reconstruction (`certify_minfuel_pmp`) gives primerDirErr = 2.000 (~180°),
signMatchFrac = 40.44% (FAIL/REVIEW). The two **disagree** — recorded honestly
as `rep.crossCheck.agree = false`, not adjudicated by suppressing either
number. `certify_minfuel_pmp`'s own header already documents itself as
`STATUS: PARTIAL` over the ~40-rev spiral (continuous-adjoint discretization
error amplified ~1e11), so the disagreement is read as the LS side being the
weaker of the two witnesses here, not as doubt on the raw-dual PASS — but this
does not close G1: it is still one family of duals (the NLP's own) checked
against a second family with a known precision limit, not an indirect solve.

## A2. The PMP necessary conditions — which do we actually check?

The textbook list, with our status. Conditions marked ★ are ones easy to omit
from an informal list but which matter here.

| condition | status |
|---|---|
| **(a) state dynamics** satisfied along the trajectory | YES — machine-tight, but **at the collocation nodes, to collocation order**. Not the same as satisfying the continuous ODE; see the caveat below. **MEASURED 2026-08-02 on DRO→tulip and the gap is not academic: NLP defect 1.4e-14 against a TRUE continuous local POSITION error of 1441 km.** The residual scales as altitude^-4.24, matching trapezoidal LTE. No second-order solve on that transfer passed the accuracy gate at any density tried. Fourth-order Hermite-Simpson at N=1600 reaches 3.3 m worst-interval position error (42 m end-to-end) and CERTIFIES. Instrument: `DRO_tulip/direct/certify/dro_residual.m` (self-tested 3/3). **Every campaign in this register quotes 1e-14 defects as evidence of a good solve; none has measured this.** |
| **(b) costate dynamics** (adjoint equation) | **Implicitly, and DISCRETELY only** — see the note below the table on what this does and does not establish. **Implicitly** — KKT stationarity w.r.t. the state variables *is* the discrete adjoint recursion, now verified generically via `foc_check` on all four campaigns (2026-07-25): ‖∇ₓL‖∞ ranges 1.5e-14 .. 6.9e-10 across earth's 9 rows, 1.2e-14 .. 2.5e-13 across CR3BP's 4 rows, 5.7e-13 on the tulip flagship, 5.9e-13 .. 1.1e-11 across ELFO's 3 rows — all PASS. Tulip additionally enforces it explicitly by least-squares (independent of the solver's own duals, but disagrees with the raw-dual read on this same flagship — see A1's Task 8 note). **Nobody checks a continuous-time costate ODE residual.** **EXTERNALLY VALIDATED 2026-08-02 on DRO→tulip, the only problem with an independent indirect answer: the defect multipliers reproduce tfMin's converged costates to 6e-04 deg median direction (λ_v AND λ_r), scale factor 0.99999 (same normalization, CoV 4e-05), Hager terminal covector to 5e-04 relative. First time the dual→costate mapping has been checked against a second opinion rather than against itself. `DRO_tulip/direct/certify/costate_compare.m`** |
| **(c) transversality** | λ_m(t_f)=0 checked via `foc_check` on all four campaigns, now on the **mapped terminal covector** (Hager 2000; Task 0, 2026-07-25 mesh-convergence-study, see A6/RETRACTION below) rather than the raw/extrapolated interval dual. Terminal elements are fully pinned, so no condition there. In the L-domain, free ΔL supplies a transversality-type stationarity row — enforced by KKT, never *reported* as a PMP quantity. **Finding I5 SUPERSEDED, formerly-real earth-2.5N miss RETRACTED** — see A6 |
| **(d) Hamiltonian constancy / H≡0 for free t_f** | see the correction below — partially; the free-t_f case is now partially checked (G4): `foc_check`'s dual-form on the ELFO min-time anchor confirms λ_t(t_f)=−1.000, exactly the theoretical ±1 value, but leaves a genuine open anomaly (`lamTimeCoV`=5.7%, expected ~0) unresolved — see A4/LEAD-4. Still not the literal H(t_f)=0 value test |
| ★ **(e) the minimum condition itself**: u\* = argmin over the admissible set, not merely ∂H/∂u = 0 | **YES, and it is our strongest check.** For control-affine dynamics with ‖β‖=1 and thr∈[0,1] the pointwise minimizer is exactly β = −p/‖p‖ with thr from sign(S) — so the primer-alignment + sign-law gates *are* the minimum condition. Note ∂H/∂u = 0 would be the wrong test: the throttle optimum is at a bound, not a stationary point |
| ★ **(f) non-triviality / normality**: (λ₀,λ) ≢ 0, λ₀ ≥ 0; no abnormal extremal | **NOT CHECKED anywhere.** In the direct/NLP form λ₀ ≡ 1 by construction, so abnormality is not even representable — it becomes a real question only for the indirect solver (goal 6) |
| ★ **(g) no singular arc** (else generalized Legendre–Clebsch / Kelley) | YES — H2 run-of-\|S\|≈0 check (earth) |
| ★ **(h) regular switching, Ṡ ≠ 0** | computed and, as of 2026-07-25, a **standing per-row line in every `foc_report`** (G3 closed as a standing line — advisory, not yet a hard gate). Live findings from that line: tulip flagship 1.4e-7 vs tol 1e-3 (25 switches), ELFO 1.33x front row 4.5e-5 vs tol 1e-3 (50 switches) — attributed to node-grazing, not costate defects — **but see finding I1 (final review, 2026-07-25): `Sd` carries the local trapezoid mesh weight and `sdotMinRel` is NOT mesh-normalized, so a small value on a refined/PSR-split mesh may be a DISCRETIZATION CONFOUND rather than physical grazing; not yet distinguished on either row** |
| ★ **(i) Weierstrass–Erdmann corner conditions** (λ, H continuous across a switch) | satisfied *by construction* — control-affine, no active state constraints; worth stating, not worth testing |
| ★ **(j) weak Legendre–Clebsch** ∂²H/∂u² ≥ 0 | vacuous at ε=0 (linear in throttle); meaningful only on the ε>0 homotopy legs |
| ★ **(k) state-constraint jump conditions** | not applicable — the state boxes are inactive; the solver's `boundSat` diagnostic is what confirms that, and it is checked |

**On (b): what "verifying the costate dynamics" actually means here.** KKT
stationarity in the state block *is* the discrete adjoint recursion — not an
approximation of it, the same equations. But that makes the check
near-tautological rather than powerful: we DEFINE the costate to be the
multiplier, so the multipliers satisfy the recursion by construction. What it
genuinely buys is (i) the solver reached a KKT point and (ii) the multipliers
we hold belong to the model we think we solved — which is precisely what
caught the `opti.dual` corruption. What it does NOT establish: the
*continuous* adjoint ODE (that is a consistency/order question governed by the
mesh, invisible to the residual); a clean recursion at nodes where other
constraints act; or any independence from the transcription (G1).

**The primal/dual asymmetry.** On the primal side the discrete-vs-continuous
gap has been *measured*: forward-integrating the exact direct control from x₀
over ~40 revs diverges by ‖r‖~3, ‖v‖~5 while the defects read 1e-14
(`GTO_tulip/direct/certify/switch_hessian.m`). **Read that correctly** — it is
a statement about IVP conditioning over 40 revs (the same wall that defeats
indirect single shooting), NOT evidence that the collocation solution poorly
approximates the continuous BVP solution. The affirmative evidence runs the
other way: two entirely independent transcriptions agree at 10 N to ~0.03%
(MEE 1377.10 vs Cartesian 1376.74 kg), both inside HMG-2004's published band.
On the dual side, the analogous gap has **never been measured**. Recorded as a
lead in `verify_common/TODO.md` §4.

**Correction worth being precise about, on (d).** Three distinct cases, and we
conflate them at our peril:
- *autonomous + **fixed** t_f* → H = **const, generally ≠ 0**. This is our
  min-fuel setup (t_f = c_tf·t_f,min, fixed, per HMG-2004 verbatim). Measured:
  H_t = −λ_t = −0.0277, constant to CoV ~1e-9. Correct and checked.
- *free t_f* → H(t_f) = 0, and if autonomous then **H ≡ 0**. This is our
  min-time anchors. **Partially checked as of 2026-07-25 (G4)** — `foc_check`
  reports the *dual form* of this condition rather than the literal H(t_f)=0
  value test; on the ELFO min-time anchor it confirms λ_t(t_f)=−1.000
  exactly, but leaves the trajectory-wide constancy of λ_t open
  (`lamTimeCoV`=5.7%, LEAD-4). The tulip min-time anchor has not been run
  through this gate at all. The literal H(t_f)=0 *value* test — the sharpest
  available, because it is a value test, not a constancy test — is still not
  built.
- *non-autonomous* (earth CR3BP: the Moon's position is an explicit function
  of t) → H is **not** constant, and we currently just label the check
  "informational" and stop. That discards a real test: the correct
  generalization is **dH/dt = ∂H/∂t**, which is computable here. Upgrading
  `hamiltonian_const_check` to test that residual would restore a genuine
  first-order check to the whole CR3BP campaign. (Recorded as G6.)

**Caveat on (a) that we should not soften.** We satisfy the *discrete* defects
to ~1e-14, not the continuous ODE. `switch_hessian`'s blocked finding is
the hard evidence: re-integrating the exact direct control from x₀ over ~40
revs diverges by ‖r‖~3, ‖v‖~5. So "dynamics satisfied at every time point" is
true of the transcription, and only true of the continuous problem up to
collocation order on the given mesh. Any claim aimed at a paper should say so.

## A3. Coverage matrix

**Updated 2026-07-25 (Tasks 6-9): `foc_check`/`foc_report` now run on all four
campaigns.** Columns "full KKT" and "Sdot ≠ 0" are the new green cells; both
were previously the `verify_pmp_mee`-only physical layer, now backed by the
generic AD instrument (A1 row 9) everywhere. "focPass" below is `foc_check`'s
own advisory `rep.pass` — a strictly stricter bar than the campaigns' existing
physical `pass` column (see A5); a FAIL here does not demote a certified row.

| campaign / target | primer + sign law | switch alignment | transversality | H constancy | full KKT | Sdot ≠ 0 |
|---|---|---|---|---|---|---|
| earth 2-body (10 → 0.1 N, + 2 PSR) | **9/9 @ 0.000° / 100%** | yes | **9/9 PASS on the mapped terminal covector** (Task 0, 2026-07-25 mesh-convergence-study). 10N/5N/2.5N re-run: `lamMassEndMapped` = 1.242e-25 / 1.405e-17 / **2.268e-18** — all at KKT-stationarity floor, `tolTrans`=1e-3 cleared by 14-18 orders of magnitude. **The former "2.5N is the single genuine first-order finding" verdict is RETRACTED** — its raw/extrapolated readings (1.265e-3, unchanged by I5) were themselves the representation artifact, not the solution; see A6/RETRACTION | yes (CoV ~1e-9) | **yes, 9/9 kktStat PASS** (1.5e-14 .. 6.9e-10 across the 9 rows — max is the 1N-PSR row) | **yes, 9/9 PASS** — mesh-normalized (I1); values 2.1 .. 1.0e2 across the ladder |
| earth CR3BP (lunar) | **4/4 @ 100% sign** (10/5/1/0.5 N; 2.5/0.2/0.1 N **NOT RUN** — deferred, large-N warm re-solves too slow for this pipeline) | yes (4/4) | **4/4 PASS** (re-run 2026-07-25): T5N moved 3.108e-3 → **6.437e-4** under the I5 endpoint correction — same artifact as earth 5N, now resolved | reported, informational (non-autonomous — correct; G6 still open, no new instrument) | **yes, 4/4 kktStat PASS** (1.2e-14 to 2.5e-13) | **yes, 4/4 PASS** |
| GTO→tulip | raw-dual **PASS** (0.058°, 100%) on flagship; LS-reconstruction **FAIL** (2.00, 40.44%) — **DISAGREE**, recorded not adjudicated (see A1 note) | yes (raw-dual) | **PASS** (5.2e-9) | no (autonomous H_t checked separately elsewhere, not via foc_check) | **yes, PASS** (5.746e-13) on the 25-switch flagship | **PASS** (re-run 2026-07-25): 1.376e-7 → **2.701e+01** under the I1 mesh normalization. The earlier reading was **entirely a mesh artifact** — the 25 switches are transversal, and the "node-grazing" attribution recorded here is **RETRACTED**. Flagship `rep.pass` is now PASS |
| GTO→ELFO | **energy seed PASS** (100%, eps=1 report-only caveat); **front row (1.33x, 50sw) PASS** (100%, exact at eps=0); min-time n/a (all-burn) | yes (fuel rows) | **3/3 PASS** — energy seed 1.00e-05, min-time 9.62e-06, front row 6.00e-09 | no (autonomous CR3BP-rotating; G4/lamTimeCoV below is the closer analog for the free-tf leg) | **yes, 3/3 kktStat PASS** (5.89e-13, 7.25e-14, 1.075e-11) — **first-ever FOC gate on this campaign** | energy seed **PASS** (8 sw); front row **FAILS** (sdotMinRel 4.485e-05, 50 sw, near-graze — same pattern as tulip — **see finding I1**); min-time n/a (all-burn) |
| min-time anchors (free t_f) | direction only, PASS on ELFO (3.30e-18 tan-max) | n/a (all-burn) | n/a (no free-mass transversality; see G4 instead) | n/a | **ELFO: yes, PASS** (7.25e-14); tulip min-time anchor not run this round | n/a |

**G4 dual-form note (ELFO min-time anchor only):** `lamTimeEnd = -1.000` —
lands exactly on the theoretical free-t_f transversality value
λ_t(t_f)=±1, a genuine confirmation, not just a threshold pass. But
`lamTimeCoV = 5.69e-2` (constancy of λ_t along an *autonomous* trajectory,
where it should be ~0) is a real, non-trivial departure, twelve orders above
the KKT-stationarity floor — **open observation**, not resolved (see Part B
§4 LEAD-4).

## A4. Gaps — what we lack

- **G1 — independence.** Earth and CR3BP verify PMP from the NLP's *own*
  duals: self-consistency, not an independent witness. It catches convergence
  to a non-extremal; it cannot catch an error shared with the transcription.
  Real fix: costates from an indirect solve (goal 6). Cheap partial fix: port
  tulip's least-squares reconstruction as a second opinion.
- **G2 — coverage. PARTIAL (2026-07-25).** ELFO now has a first-order gate
  (3 artifacts on the nominal rung: energy seed, min-time anchor, 1.33x/50sw
  front row) — the "none" cell is gone, but only one thrust rung is covered,
  not a ladder. CR3BP is verified at 10/5/1/0.5 N (4 of 7 rungs); 2.5 N was
  skipped (optional) and 0.2/0.1 N were explicitly **DEFERRED, NOT RUN**
  (large-N warm re-solves too slow for this pipeline at the time). Tulip's
  25-switch flagship is covered; the PSR front (stage 5c wiring exists but was
  never exercised end-to-end — see Task 8 concern 1) and the min-time anchor
  are not. Tulip's certifier still carries the known false negative on
  node-grazing switches (see the tulip TODO C3) — now independently confirmed
  by `foc_check`'s own `sdotMinRel` FAIL on the same flagship row.
- **G3 — Sdot != 0. CLOSED as a standing report line (2026-07-25).**
  `foc_check`/`foc_report` now compute and print `sdotMinRel` as a per-row
  line on every campaign, not buried inside `second_order_nlp.m`'s
  NOT-APPLICABLE verdict — the *reporting* gap is closed. It is advisory
  (folded into `rep.pass`, not a hard gate), and it already caught two real,
  honestly-reported findings: the tulip flagship (1.4e-7 vs tol 1e-3, 25 sw)
  and the ELFO 1.33x front row (4.485e-5, 50 sw) both attributed to
  node-grazing switches — consistent with M2/near-degenerate bang-bang
  minima, not a new defect. **Finding I1 (final-review, 2026-07-25):** that
  attribution is unverified. `Sd` (the discrete switching-function residual
  the gate differences) carries the local trapezoid mesh weight, and
  `sdotMinRel` is **not mesh-normalized** — it is not divided by the nodal
  weight or differenced against `dsigma`. On a refined or PSR-split mesh a
  small `sdotMinRel` value can therefore be a **DISCRETIZATION CONFOUND**
  (finer mesh → smaller raw residual by construction) rather than physical
  node-grazing. Neither the tulip flagship nor the ELFO front row has been
  re-checked with a mesh-normalized `Sd`. Mesh-normalizing `Sd` (divide by
  the nodal weight, difference by `dsigma`) is **REQUIRED work before G3 is
  promoted from an advisory standing line to a hard gate.**
- **G4 — free-final-time Hamiltonian condition. PARTIAL, dual-form confirmed
  on one anchor (2026-07-25).** `foc_check` reports the dual-form of this
  condition (λ_t(t_f) vs its constancy) rather than the classical value test
  H(t_f)=0; on the ELFO min-time anchor `lamTimeEnd = -1.000` lands exactly on
  the theoretical λ_t(t_f)=±1 value — a genuine confirmation. But
  `lamTimeCoV = 5.69e-2` (expected ~0 under this autonomous CR3BP model) is
  open and unresolved (Part B §4 LEAD-4); the tulip min-time anchor
  (`tfMin_tulip`) was not run this round. Still not the literal H(t_f)=0 value
  test A2 calls for — G4 stays open as a derivation question even where the
  dual-form check passes.
- **G5 — no mesh bands on the gates.** Switch counts are known mesh-sensitive;
  the first-order gate values are reported as single numbers, not bands.
- **G6 — the non-autonomous Hamiltonian test is abandoned rather than
  generalized. Unchanged by Tasks 6-9.** Under lunar gravity H is genuinely
  not constant, and the CR3BP driver correctly refuses to gate on
  constancy — but it then checks nothing. `foc_check` does not build the
  `dH/dt = ∂H/∂t` residual either (out of scope for Tasks 1-9); test that
  instead (see A2's correction on (d)). Still open, no new instrument.
- **G7 — normality / abnormal extremals never examined.** Invisible in the
  direct formulation (λ₀ ≡ 1 by construction); becomes a live question the
  moment an indirect solver exists.

## A5. What can and cannot be claimed today

**Can:** machine-tight *on the transcription*, PMP-extremal, with the pointwise
minimum condition verified (primer + sign law), bang-bang structure confirmed,
no singular arcs, and switch times consistent with the reconstructed switching
function.

**Cannot:** (i) that the continuous ODE is satisfied to the same accuracy — see
A2's caveat on (a); (ii) that any individual switch is *proven* transversal —
the mesh-normalized Ṡ statistic shows no evidence of grazing, but a
single-mesh threshold cannot establish regularity (two-mesh confirmation is
pre-promotion work, §A6); (iii) *strict* local minimality — Part B gives weak local
minimality on part of the tulip set and nothing elsewhere; (iii) anything
global.

## A6. Before hard-gate promotion

Findings from the final whole-branch review of the foc-gate-layer branch
(2026-07-25) that are deliberately **not** coded yet — checker-semantics
issues in `foc_check.m` itself, not campaign findings. G3 (Sdot != 0) and the
direction-check row currently run as *advisory* lines feeding `rep.pass`;
promoting either to a hard gate should wait on these:

- **I1 — `sdotMinRel` is not mesh-normalized.** `Sd` (the discrete
  switching-function residual) carries the local trapezoid mesh weight
  uncorrected. On a refined or PSR-split mesh, a small `sdotMinRel` can be a
  **discretization artifact** (finer mesh → smaller raw residual) rather than
  genuine regular switching. Required fix: mesh-normalize `Sd` — divide by
  the nodal weight, difference by `dsigma` — before this can be trusted as a
  hard gate. Currently contaminates the tulip-flagship and ELFO-front-row
  "node-grazing" reads (A3, A4/G3).
- **I2 — the direction-check (`dirTanMax`/`dirTanMed`) is not an independent
  check, and is sign-blind.** The tangential-residual computation is bounded
  by `kktStat` *by construction* (it is the same stationarity condition
  projected onto the tangent space of `‖β‖=1`), so a PASS here is not
  independent evidence beyond the KKT-stationarity row — it is currently
  reported as if it were a separate instrument. It is also **sign-blind to
  the ±p/‖p‖ branch**: tangential-residual-only cannot distinguish the primer
  direction from its antipode, so it can pass on a solution pointing the
  thrust the wrong way. Needs (a) a `betaNorm`-dual exclusion so the
  computation is not trivially guaranteed by the constraint it is meant to
  check, and (b) a multiplier-sign test to close the ±p/‖p‖ ambiguity.

Both are recorded here rather than fixed now — see the branch's final-review
report for why they were scoped out of this fix wave.


### FIX PASS + ATTRIBUTION (2026-07-25, commit `6b9c4a7`)

Three of the review's items are now **shipped**, and — because each corrected
check reports its legacy value beside it — the outstanding advisory FAILs are
now **attributed** rather than merely moved.

| row | check | legacy | corrected | attribution |
|---|---|---|---|---|
| tulip flagship (25 sw) | Ṡ regularity | 1.376e-07 | **2.701e+01** | **mesh artifact** → PASS |
| ELFO front 1.33× (50 sw) | Ṡ regularity | 4.485e-05 | **1.193e+02** | **mesh artifact** → PASS |
| earth 5 N | transversality | 3.189e-03 | **5.215e-04** | **endpoint bias** → PASS |
| CR3BP 5 N | transversality | 3.108e-03 | **6.437e-04** | **endpoint bias** → PASS |
| earth 2.5 N | transversality | 1.265e-03 | 1.265e-03 | **REAL** — unchanged (**RETRACTED 2026-07-25, Task 0 — see the RETRACTION block below: this was itself a representation artifact, not a solution defect**) |
| ELFO energy seed (ε=1) | sign law etc. | PASS incl. "100%" | PASS, 3 checks `--` | **false PASS removed** |

**Of five outstanding advisory FAILs, four were discretization artifacts and
one is real.** *(Superseded 2026-07-25 by Task 0 of the mesh-convergence plan:
**all five were artifacts.** See the correction under item 2 and the Task-0
block below.)* Two corrections to things previously recorded here as findings:

1. **The tulip "node-grazing switches" reading was an artifact.** It had been
   recorded as the FOC layer quantifying a behaviour the tulip campaign already
   suspected. Properly mesh-normalized the flagship reads 27 — comfortably
   transversal. The suspicion may hold on other grounds; this layer was never
   evidence for it.
2. ~~**The transversality misses split.**~~ **THIS ITEM IS RETRACTED IN FULL
   (Task 0, 2026-07-25).** It read: *the endpoint bias only appears where λ_m
   is still moving at t_f — i.e. where the final arc burns; both 5 N rows moved
   and pass, 2.5 N (final arc coasting) did not move at all, so its miss is a
   property of the solution … the only genuine first-order finding in the whole
   set.*
   **Every part of that is wrong**, and instructively so. Under the correct
   statistic — the Hager-mapped terminal covector — **all** rows land at the
   KKT floor regardless of final-arc type: 10 N 1.24e-25, 5 N 1.41e-17,
   **2.5 N 2.27e-18**, tulip flagship 1.19e-24. There is no burn/coast split
   because there was never a transversality violation to split.
   **The methodological lesson is the point:** the burn-vs-coast story was a
   physically plausible mechanism constructed on top of two *mis-defined*
   statistics (the raw interval dual and its extrapolation). It explained the
   numbers, so it was recorded as confirmation. A satisfying mechanism that
   fits wrong measurements is more dangerous than an unexplained anomaly,
   because it stops the search. The extrapolation "fixing" 5 N was not
   evidence for the mechanism — it was one wrong number moving nearer zero. **This conclusion does not survive Task 0 (2026-07-25) — the
   I5 extrapolation, though itself an improvement, was still an
   *approximation* to λ_m(t_f), and approximating the wrong thing more
   carefully does not make it the right thing: the exact discrete object was
   never Lambda_N or any function of it alone. See the RETRACTION block below.**

**Caveat that survives (do not over-read the new PASSes).** GPT's point stands:
no threshold on a *single* mesh establishes transversality of a switch.
Regularity should be asserted only if the normalized statistic stays above the
floor on **two** meshes, and the 1e-3 gate is provisional until recalibrated on
deweighted data. The Ṡ PASSes therefore mean "no evidence of grazing," not
"proven regular" — which matters, since this is Maurer ingredient (ii) and any
second-order certificate will lean on it.

**Regression lock:** `verify_common/tests/test_foc_mesh_invariance.m` asserts
the property in both directions on the toy problem — across a 4× refinement the
normalized statistic holds (ratio 1.02) while the legacy one falls by 4.1
(exactly the refinement factor, i.e. it was reporting `h`).

**Still open from the review:** I2 (direction-check independence + sign
blindness — note the *previously recorded fix was wrong*, see above), the
sign-resolution ambiguity guard, the δ_w relabel, and three robustness asserts.
Tracked in `verify_common/TODO.md` §1.2, §1.5, §3.

### RETRACTION — mapped terminal covector supersedes I5, earth 2.5 N finding withdrawn (Task 0, 2026-07-25 mesh-convergence-study)

The I5 fix above (endpoint extrapolation) was recorded, correctly at the time,
as a *partial* fix — both external reviewers also named a more principled
alternative: the exact discrete **mapped/transformed terminal covector**
(Hager, "Runge-Kutta Methods in Optimal Control and the Transformed Adjoint
System," Numer. Math. 87, 247-282, 2000),

    lamHat_f = (h_N/2)*L_x(N+1) + (I - (h_N/2)*f_x(N+1))' * Lambda_N

whose mass component is **exactly zero by the discrete KKT system** whenever
the final mass is genuinely free with no terminal cost/constraint on it — not
an estimate of λ_m(t_f), but the same full-Lagrangian-gradient identity
`kktStatInf` already certifies elsewhere, read at one specific row/node. It is
now implemented in `foc_check.m` §(6) as `rep.lamMassEndMapped` (reusing the
already-assembled `gL`, no new differentiation — see the file for the
derivation) and **gates** transversality; `rep.lamMassEndRel` (I5
extrapolation) and `rep.lamMassEndRelOneSided` (raw one-sided dual) are
retained as reported-only attribution companions. Regression test:
`verify_common/tests/test_foc_terminal_covector.m`.

**Re-run of the affected rows** (warm re-solve, live model, same procedure as
every other `foc_check` row):

| row | mapped (gated) | extrapolated (I5) | raw one-sided (legacy) |
|---|---|---|---|
| earth 10 N | 1.242e-25 | 1.635e-10 | 1.939e-10 |
| earth 5 N | 1.405e-17 | 5.215e-04 | 3.189e-03 |
| **earth 2.5 N (the open finding)** | **2.268e-18** | 1.265e-03 | 1.265e-03 |
| tulip flagship (25 sw) | 1.192e-24 | 3.229e-09 | 5.239e-09 |

**The decisive number: earth 2.5 N's mapped transversality is 2.268e-18** —
14 orders of magnitude inside `tolTrans`=1e-3, and at the same machine-KKT
floor as every other row in the table (`kktStatInf` for this row is 5.0e-10).
**The "single genuine first-order finding" recorded above (I5 fix-pass block,
A3 earth-2-body row) is RETRACTED.** It was never a property of the 2.5 N
solution — it was a property of the *statistic*: Lambda_N (and its I5 linear
extrapolant, exact only for a linear costate locally) is not λ_m(t_f), it is
a sample near sigma_end − h_N/2, and at 2.5 N the final arc's λ_m is moving
fast enough there that neither the raw dual nor its linear extrapolation gets
within 3 orders of the true (zero) endpoint value. The I5 extrapolation was a
real improvement over the raw dual in general (it resolved the earth-5N and
CR3BP-5N misses correctly, for the correct reason — see the FIX PASS table
above, unchanged) — it was simply the wrong tool for a row where the local
extrapolation itself carries residual curvature error. The mapped covector
has no such failure mode: it is exact by construction regardless of how fast
λ_m is moving at t_f, which is why all four rows above — burning or coasting,
tight or marginal — land at floor together.

**Consequence for A3/A6 above:** every "PASS (I5 endpoint bias)" and the one
"REAL (unchanged)" verdict recorded in the FIX PASS + ATTRIBUTION table and
the A3 earth-2-body coverage-matrix cell predate this fix and are superseded
by the mapped-covector re-run in this block, not by further edits to the
historical text (kept as a record of what was believed and why, per this
file's own append-don't-rewrite convention, §"Scope discipline").

**What this changes and does not change.** It closes the one open
first-order finding in Part A with no primal/solution implications — the 2.5
N certified row itself is untouched, only the *instrument that reads its
transversality* changed. It does NOT retroactively validate M1/M2/M3 (Part
B, second order) or any claim beyond transversality; those remain as
recorded. It also sharpens the general lesson already half-stated by I5: a
one-sided or linearly-extrapolated dual is never legitimately called
λ(t_f) — only the mapped covector is.

### External review of the checker itself (2026-07-25)

The generic layer and its explainer were sent for three-way review
(GPT-5.6-terra + Gemini 3.1 Pro + host Claude). Verbatim reviews:
`verify_common/doc/review_2026-07-25_{gpt56terra,gemini31pro}.md`.
Concrete fixes now specified in `verify_common/TODO.md` §1.

**Validated** (do not re-litigate): the switching-function recovery is exact
(the active bound multiplier cancels `Sd` identically, so zeroing its row
recovers magnitude and sign); the dual→costate map is the correct
*stationarity* combination, re-derived independently by both reviewers.

**I1 sharpened.** Both reviewers judged the current readings (tulip 1.4e-7,
ELFO front row 4.5e-5) **unable to diagnose grazing**, and converged on the
same fix (deweight by nodal trapezoid weight → divided difference in σ →
non-dimensionalize). GPT adds: no threshold on a *single* mesh is defensible —
require stability across two meshes. The same deweighting is needed for the
singular-arc test, which otherwise false-positives in densified regions.

**I2 CORRECTED.** The previously recorded fix — "exclude the betaNorm dual
before projecting" — does **not** work: the projector annihilates the radial
term anyway, `P_b(q+2μb) = P_b q`. Independence requires building the gradient
from a *different* source (defect duals only, or a separate Hamiltonian
reconstruction). The branch test is the cone multiplier's sign (`μ>0` ⟺
minimizer, under `L = L0 + μ(b'b−1)`).

**I5 confirmed by both**, with an agreed extrapolation fix and a more
principled discrete-endpoint-covector alternative whose mass component is zero
by terminal-node stationarity. **DONE (Task 0, 2026-07-25 mesh-convergence-
study): the principled alternative is now built and gated — see the
RETRACTION block below, which also withdraws the earth-2.5N finding this
extrapolation fix could not resolve.**

**New (not previously recorded):** the advisory verdict is **not ε-aware** —
`signPct` / singular-arc / `sdotMinRel` are folded into `rep.pass` even on ε>0
homotopy legs where the bang-bang law does not apply; and the `certLocalMin` /
"LOCAL MIN" naming of the δ_w line overstates what zero inertia correction
proves (relabel campaign-wide, including `ipopt_certify.m`).

---

# PART B — Second order (local minimality)

**State: partially delivered.** Weak local minimality is certified on 12 of 17
ε=0 tulip rows by IPOPT's native inertia; **strict** local minimality is
certified nowhere. **Update 2026-07-25 (LEAD-0 delivered):** the "three other
campaigns have no second-order verdict at all" clause above is now stale —
the native-inertia port (`foc_ipopt_inertia`) landed on earth 2-body, earth
CR3BP, and ELFO in Tasks 6, 7 and 9. All three now carry a weak-local-min
verdict per row wherever the first-order gate ran; none of them reach strict.

## 1. Status matrix

| instrument | campaign | verdict | where |
|---|---|---|---|
| **IPOPT native inertia (δ_w)** | GTO_tulip | **DELIVERS** — wired into production PSR; 137/137 rows carry a verdict. **ε=0 (bang-bang): 12 certified / 5 not.** ε>0: 100 / 20. Flagship 25-switch row re-checked via the new generic port (`foc_ipopt_inertia`, Task 8): **NOT CERTIFIED** (δ_w tail max 1.8e-8, just above the 1e-8 tol) — borderline, consistent with the campaign's known weak/near-degenerate Hessian at fuel-optimal rungs | `GTO_tulip/direct/certify/ipopt_certify.m`, called from `run_gto_tulip.m stage 5c` + `run_one.m`; verdicts stored as `ipoptCert` in `direct/data/psr_data_*.mat`. Generic port cross-check: `run_foc_tulip.m` |
| **IPOPT native inertia (δ_w) — ported (LEAD-0, 2026-07-25)** | earth 2-body | **9-row ladder, mixed: LOCAL MIN at 10N / 0.2N / 0.1N (δ_w=0); NOT CERTIFIED at 5N (4.51e-05) / 2.5N (7.19e-05) / 1N (2.93e-04) / 0.5N (4.01e-05) / 1N-PSR (5.31e-03) / 0.5N-PSR (1.34e-05→1.77e-04)** | `verify_common/foc_ipopt_inertia.m`, wired via `run_foc_mee.m` + `run_verify_pmp_all.m` (Task 6) |
| **IPOPT native inertia (δ_w) — ported (LEAD-0, 2026-07-25)** | earth CR3BP | **4/4 LOCAL MIN, δ_w=0** at all four covered rungs (10N/5N/1N/0.5N) — a clean weak-local-min certificate everywhere the first-order gate currently runs; 2.5/0.2/0.1 N not run (coverage gap, G2) | `verify_cr3bp_pmp.m` (Task 7) |
| **IPOPT native inertia (δ_w) — ported (LEAD-0, 2026-07-25)** | GTO_ELFO | **energy seed: LOCAL MIN** (δ_w=0 tail); **min-time anchor: LOCAL MIN** (δ_w=0 tail); **front row (1.33x, 50sw): NOT CERTIFIED** (δ_w max 1.98e-06) — same near-degenerate pattern as the tulip flagship, first-ever second-order verdict on this campaign | `run_foc_elfo.m` (Task 9) |
| NLP reduced-Hessian SOSC | earth 2-body | **WEAK_MIN** @10 N (270 flat directions); INCONCLUSIVE @5 / 2.5 N; ERROR @1 / 0.5 N; scale-skip above `maxNullDim=10000` | `earth_elliptic_to_geo/direct/verify/sosc/`, `process/DESIGN_sosc.md` §11–12 |
| NLP SSOSC via KKT inertia | GTO_tulip | **NOT APPLICABLE (structural)** | `GTO_tulip/direct/certify/second_order_nlp.m` (FINDING, 2026-07-12) |
| Maurer–Osmolovskii switching-time Hessian | GTO_tulip | **BLOCKED** (forward-flow conditioning) | `GTO_tulip/direct/certify/switch_hessian.m` (FINDING, 2026-07-12) |
| conjugate point (Jacobi field) | catalog pipeline (any campaign with an `ms_bvp` solution) | **BUILT 2026-08-08** (migration #3) — free-final-time quotiented form, the two exact degeneracies (costate scaling, time reparametrization) handled, λ_m excluded; **runs**: golden-cells regression (3 engine cells, verdicts pinned), opt-in in `ms_tfmin` (`.conjTest`), advisory hook in `foc_check` (`opts.msInfo`, migration #5). **RUN AT CATALOG SCALE 2026-08-23** (`conj_catalog_pass`, K=24): all 15,896 entries of the five catalogs decided — **15,895 PASS / 1 FAIL / 0 unverified**; verdicts stored per sheet (`conj_pass`/`conj_ncross`/`conj_atfinal` + `.conj_test` provenance, schema-validated). The one refuted entry: DPO τ=2 → tulip Np=7, phases (2/3, 2/3), 15 N (nCross=1 interior, re-solve dz=1.6e-12). Still open: no Earth-campaign ms solution has exercised the `foc_check` hook; **free-time only** — the min-energy (fixed-t_f) entries have no applicable form (open item, `costate_common/TODO.md`). Necessary condition only — a pass does not certify, a fail refutes | `costate_common/ms_conjugate_test.m` (header = derivation); `golden_cells.m`; `DRO_tulip/indirect/ms_tfmin.m`; `verify_common/foc_check.m:318` |

**Net (corrected 2026-07-25, twice now):** second order is **not** at zero,
and it is no longer tulip-only. The IPOPT native-inertia certificate delivers
on all four campaigns as of Tasks 6-9 (LEAD-0 closed) — weak local
minimality per row wherever the first-order gate runs, mixed verdicts
everywhere (some rows LOCAL MIN, some borderline NOT CERTIFIED). An earlier
version of this register said "zero delivered certificates" on the tulip
alone — that was wrong, and the miss is instructive: the instrument had **no
FINDING block** (unlike the two blocked ones), so it was invisible to a
survey that read headers. *An instrument that quietly works is easier to
lose than one that loudly fails.* None of this reaches **strict** local
minimality anywhere — that is still exclusively §5's open item.

**What δ_w = 0 does and does not prove.** IPOPT's inertia-controlled linear
solver adds Hessian regularization only when the reduced Hessian is indefinite;
converging with δ_w = 0 over the final iterations therefore means no negative
curvature *on the well-scaled barrier system at the final μ*. Three caveats,
all in `ipopt_certify.m`'s own header and none of them fatal:
- it is a statement about the **barrier subproblem** at small μ, not exactly
  the original NLP's critical cone;
- at **ε > 0** the problem is strictly convex in throttle, so this is a strict
  local minimum of the NLP and is the easy case;
- at **ε = 0** the barrier supplies the throttle-direction curvature, so it
  certifies a **weak (non-strict)** local minimum — fully consistent with M1,
  not a contradiction of it. **Strictness for bang-bang lives in the switching
  times**, which is precisely what §5 is about.

So the honest ladder is: *extremal* (Part A) → *no negative curvature / weak
local min* (this instrument) → *strict local min* (§5, still missing).

First order, by contrast, is in good shape — see Part A, and the 9/9 sweep in
`earth_elliptic_to_geo/direct/results/verify_pmp_all.mat`.

Note that Maurer ingredients (ii) and (iii) — regular switching and primer /
Legendre direction optimality — are **first-order** checks and already exist
(A1 rows 1 and 8). What is missing is only ingredient (i), the strict
curvature test.

---

## 2. The three blocking mechanisms

These are the transferable content. Each is a property of the *problem class*,
not of an implementation.

### M1 — bang-bang extremals are weak, non-strict minima
The fuel objective is linear in throttle, so naive control-space curvature is
zero and the reduced Hessian is genuinely flat in many directions (270 of them
at earth 10 N, with a clean spectral gap: ~270 eigenvalues at 1e-10, nothing
until ~1e-4). **Strict SOSC is generically unreachable for min-fuel.** A
WEAK_MIN verdict is the honest ceiling of the NLP-level test, not a
tolerance-tuning failure.

### M2 — strict complementarity fails at the throttle bounds
The tulip's KKT-inertia test was supposed to escape M1: at a bang-bang
solution the throttle sits at a bound on every arc, so with *strict*
complementarity (switching function `S != 0` ⇒ nonzero bound multiplier) those
directions leave the critical cone and the remaining curvature is real. In
practice the solver parks the throttle **near but not at** its bounds
(`s ~ 0.02 / 0.98`, edge ~99%), the bound constraints are not active, the
degeneracy is never removed, and the test is inapplicable. This is structural
for an ε→0 homotopy that approaches bang-bang rather than imposing it.

### M3 — forward flow across many revs destroys the base point
The correct instrument for a bang-bang extremal is the Maurer–Osmolovskii
reduced problem over the *k switch times*. The tulip built it as a
**forward-flow** formulation and it is blocked: re-integrating the exact
direct control from `x0` across ~40 revs diverges by `||r||~3`, `||v||~5`,
`|t|~10`. The collocation solution satisfies the *discrete* defects but does
not correspond to a nearby *continuous* trajectory under forward integration,
so the reduced problem's base point is infeasible (`||c|| ~ 20`), the reduced
gradient is not ~0, and any Hessian there is meaningless.

To its credit `switch_hessian.m` **detects** this via a `baseFeas`
self-check and returns BLOCKED rather than a bogus certificate. This is the
same conditioning wall that defeats indirect shooting (`ifs`, `ms_band`) — it
is a property of long multi-rev flows, not of the second-order test.

---

## 3. Ruled out — do not redo

- **Tuning `inertiaZero` to rescue the NLP test.** `ldl`/Gould pivot-sign
  inertia is unreliable on this near-singular KKT (mis-signs 56 pivots at
  zt=1e-9); its correct window is *disjoint* from `eig`'s, so no single
  threshold serves both. `eig` is the gold standard and `inertiaZero=1e-9` is
  settled. Reverting this has already been tried and undone once.
- **Expecting strict SOSC from the NLP for a min-fuel row.** See M1.
- **Any forward-shooting formulation of the switching-time problem at
  multi-rev scale.** See M3.
- **Planning a fresh "tier 1/2/3" program per campaign.** Tiers 1 and 2 are
  already built (earth, tulip) and their verdicts are in §1.

---

## 4. Open leads

**LEAD-1 (high value, cheap) — the SOSC "recovery wall" at 1 / 0.5 N may not
exist.** `earth_elliptic_to_geo/TODO.md` item 0 records two distinct scale
walls, the first being "*recovery* — the warm re-solve itself fails at 1 N /
0.5 N (n≈16.5k)", and proposes building a scalable warm-resolve recovery.
**Contradicting evidence (2026-07-25):** `refresh_duals_mee.m` performs
essentially that same warm re-solve (same `warmTight`, same `maxIter`, same
`returnModel` rebuild) and it *succeeded* at both rungs — 1 N `Solve_Succeeded`
defect **8.27e-14**, 0.5 N defect **1.83e-13** — comfortably inside
`sosc_recover_kkt`'s `recoverOK` bar of `maxDefect < tol.feas = 1e-8`.
If that reproduces under `sosc_recover_kkt` itself, the ERROR verdicts come
from *downstream* KKT assembly / dense null-space at scale (the second,
acknowledged wall), and the proposed recovery fix is aimed at a non-problem.
**Test first, build second.**

**LEAD-2 — certificates must be evaluated at tight optima.** Six of nine earth
rows improved on warm re-solve (5 N +0.235 kg, 2.5 N +0.51, 1 N +0.44, 1 N PSR
+0.47, 0.5 N +0.02, 0.5 N PSR +0.02); only 10 N / 0.2 N / 0.1 N were already
tight. A second-order certificate computed at an under-converged point
certifies nothing. Re-solve the six before certifying them.

**LEAD-0 — DONE (2026-07-25).** Ported the IPOPT native-inertia certificate
(`foc_ipopt_inertia`) to the other three campaigns via Tasks 6, 7, 9. It
delivered exactly as predicted: a weak-local-min verdict per row for a few
hours of plumbing, no new instrument needed — see Part B §1 for the per-row
verdicts (earth 2-body mixed 3 LOCAL MIN / 6 NOT CERTIFIED across 9 rows;
CR3BP clean 4/4 LOCAL MIN; ELFO 2 LOCAL MIN + 1 NOT CERTIFIED). It sidesteps
the conditioning wall that defeats our own factorization (M3 /
`second_order_nlp`'s tolEig noise floor), because MUMPS tests inertia on the
well-scaled system at every iteration — confirmed working identically on all
three new campaigns, not just tulip.

**LEAD-3 — M2 might be removable by construction.** If the throttle were
*pinned* to its bounds on a frozen arc structure (rather than approached via
ε→0), strict complementarity would hold and the KKT-inertia test might apply.
This is close to what the switching-time reduction does anyway, so it may not
be worth a separate build — but it has never been tried and it would reuse
existing machinery.

**LEAD-4 (new, 2026-07-25) — G4 dual-form partial: the ELFO min-time anchor's
free-t_f condition is half-confirmed, half-open.** `foc_check` on
`mintime_elfo.mat` (Task 9) gives `lamTimeEnd = -1.000e+00`, landing exactly
on the theoretical free-t_f transversality value λ_t(t_f)=±1 — a genuine,
notable confirmation that the raw-dual costate extraction is correct on this
anchor. But `lamTimeCoV = 5.69e-02` — the coefficient of variation of λ_t
along the trajectory, which should be ~0 under this autonomous CR3BP
rotating-frame model — is a real, non-trivial departure, twelve orders above
the kktStat floor (7.25e-14). Two candidate explanations, neither adjudicated:
(1) the two-primary Sundman-clock + `cScale` free-t_f slack-state structure
may distribute the horizon condition across the `cScale` adjoint rows
differently than `foc_check`'s plain `lamTimeCoV` computation assumes —
`foc_check`'s own `horizonNote` for `'freetf-cscale'` already flags this as
"value-form H(tf) check ... informational, derivation pending"; (2) the
duals may simply be less converged than the primal at a free-t_f slack
state. Open — feeds directly into G4 and the eventual freetf-cscale
horizon-condition derivation; should travel with the `-1.0` confirmation
above, not be flattened to a bare PASS.

---

## 5. Decision — what to build next, and why

**Build the STM / multiple-shooting switching-time Hessian.**

`switch_hessian.m`'s own header already specifies it: keep the collocation
trajectory as the (feasible) base, obtain switch-time sensitivities from the
state-transition matrix,

    dPsi/dsigma_i = Phi(sigma_f, sigma_i) * [ f_burn - f_coast ]

integrated along that base, and form the second variation with **segment
matching**, so no quantity is forward-propagated across all revs. Multiple
shooting defeating a forward-shooting conditioning wall, exactly as it does in
the indirect campaigns.

Rationale:
1. It is the **only live path to STRICTNESS** for a bang-bang solution. The
   IPOPT native-inertia certificate (§1) already gives *weak* local minimality
   on 12/17 ε=0 tulip rows, and LEAD-0 extends that cheaply to the other
   campaigns — but a weak minimum is exactly what M1 predicts and it cannot be
   sharpened by any NLP-level test (M1 kills the reduced-Hessian route, M2
   kills its repair). Strictness lives in the switching times.
2. It is a **build, not a research problem**: the formulation is already
   written down, and the blocked forward-flow version supplies the arc
   bookkeeping, the `baseFeas` self-check, and the output contract.
3. It is **campaign-agnostic**: the same second variation serves earth 2-body,
   CR3BP, tulip and ELFO — one instrument closes the Tier-2 row for all four.
4. Its STM/variational core is **shared with the conjugate-point test**, so it
   partly pre-builds Tier 3 and dovetails with the indirect goal (6).

Sequencing note: do LEAD-2 (re-solve the loose rows) before certifying
anything, and LEAD-1 (one cheap re-test) before investing in SOSC recovery.

**Sufficiency claim to aim for.** A positive-definite projected switching-time
Hessian, *together with* (i) strict bang-bang `Sdot != 0` at every switch and
(ii) primer / strengthened Legendre–Clebsch optimality of the thrust direction
— both already reported by the first-order verifier, and now trustworthy since
the dual fix — is the Maurer–Osmolovskii sufficient condition for a bang-bang
local minimum. That is the target statement for the methods paper.

---

## 6. Experiment log (both orders)

Append one entry per experiment — first- or second-order: date, instrument,
campaign/row, verdict, and what it changed in Part A or §1–§5.

| date | instrument | target | verdict | consequence |
|---|---|---|---|---|
| 2026-07-12 | NLP SSOSC (KKT inertia) | tulip 1.15× refined | NOT APPLICABLE | established M2 |
| 2026-07-12 | switching-time Hessian (forward flow) | tulip 1.15× refined | BLOCKED (`baseFeas` ~20) | established M3; specified the STM fix |
| 2026-07-19 | reduced-Hessian SOSC | earth 10 N | WEAK_MIN (270 flat) | established M1 |
| 2026-07-19 | reduced-Hessian SOSC | earth 5 / 2.5 N | INCONCLUSIVE | genuine precision limit, not a bug |
| 2026-07-19 | reduced-Hessian SOSC | earth 1 / 0.5 N | ERROR / scale-skip | recorded as a recovery wall — **disputed, see LEAD-1** |
| (ongoing) | **IPOPT native inertia (δ_w)** | tulip PSR, 137 rows | **12/17 certified at ε=0**; 100/120 at ε>0 | second order is not at zero; weak local min for bang-bang |
| 2026-07-25 | first-order PMP (primer/sign/alignment) | earth 9 rows + CR3BP 10/5 N | **PASS 0.000° / 100%** | Part A row 1–3 green for earth; makes Maurer ingredients (ii)/(iii) trustworthy |
| 2026-07-25 | full KKT certificate from raw `lam_g` | earth 10 N | ‖∇ₓL‖∞ 1.5e-14, tangential ∂L/∂β 8e-17 | Part A row 6; also root-caused the dual bug |
| 2026-07-25 | **external three-way review of the FOC layer itself** | `verify_common` core + explainer | recovery + dual map VALIDATED; I1/I2/I5 sharpened; 2 new findings (ε-awareness, δ_w naming) | reviews archived in `verify_common/doc/`; drove the fix pass below |
| 2026-07-25 | **FOC fix pass** (ε-aware verdict, mesh-normalized Ṡ, endpoint-corrected transversality) | tulip flagship, ELFO front + energy seed, earth 10/5/2.5 N, CR3BP 5 N | **4 of 5 advisory FAILs were artifacts; 1 real (earth 2.5 N transversality)** | tulip grazing claim RETRACTED; §A6 fix-pass block |
| 2026-07-25 (Task 6) | `foc_check`/`foc_report` full ladder | earth 2-body, 9 rows (10→0.1 N + 2 PSR) | **7/9 focPass** — 5N (3.189e-3) and 2.5N (1.265e-3) FAIL, transversality only, marginal (~3x tol); kktStat 1.5e-14 .. 6.9e-10 across all 9 rows, all PASS; `foc_ipopt_inertia` mixed (3 LOCAL MIN / 6 NOT CERTIFIED) | A3 earth 2-body row; A1 row 9 first real catch of the generic layer; Part B §1 earth row; LEAD-0 (earth) done; **final review: 5N/2.5N misses now also flagged under finding I5 (one-sided endpoint-dual bias, unadjudicated)** |
| 2026-07-25 (Task 7) | `foc_check`/`foc_report` subset | earth CR3BP, 4 rows (10/5/1/0.5 N, φ0=0; 2.5N skipped, 0.2/0.1N DEFERRED not run) | **3/4 focPass** — T5N FAILs transversality only (3.108e-3), same marginal pattern as earth; kktStat 1.2e-14 to 2.5e-13 all PASS; `foc_ipopt_inertia` **4/4 LOCAL MIN** (δ_w=0) | A3 CR3BP row updated; G2 partially closed (4/7 rungs); Part B §1 CR3BP row; LEAD-0 (CR3BP) done; **final review: T5N miss now also flagged under finding I5** |
| 2026-07-25 (Task 8) | `foc_check`/`foc_report` + LS cross-check | GTO→tulip, 25-switch flagship (`sundman_minfuel_certified.mat`) | KKT stationarity **PASS** (5.746e-13); primer/sign **PASS** (raw dual, 0.058°/100%) vs LS-reconstruction **FAIL** (2.00/40.44%) — **DISAGREE**, recorded not adjudicated; overall `rep.pass` **FAIL** (sdotMinRel 1.4e-7 node-graze + δ_w tail 1.8e-8 borderline) | A1 note on the two-source disagreement; A3 tulip row; G1 gets new data, not closed; G3 standing-line closed; Part B §1 tulip row cross-checked with the generic port; **final review: the "node-graze" read on sdotMinRel is now finding I1 — unverified without mesh-normalizing `Sd`** |
| 2026-07-25 (Task 9) | `foc_check`/`foc_report`, first-ever ELFO gate | GTO→ELFO, 3 artifacts (energy seed, min-time anchor, 1.33x/50sw front row) | energy seed **ADVISORY PASS** (KKT 5.89e-13, eps=1 sign-law caveat); min-time anchor **ADVISORY PASS** (KKT 7.25e-14, `lamTimeEnd=-1.000` CONFIRMS G4 dual-form, `lamTimeCoV=5.69e-2` OPEN); front row **ADVISORY FAIL** (sdotMinRel 4.485e-5 near-graze + δ_w 1.98e-6 NOT CERTIFIED) | A3 ELFO row goes from "none" to 3/3 covered on the nominal rung; G2 partial; G4 partial-dual-form-confirmed; LEAD-4 opened; Part B §1 ELFO row; LEAD-0 (ELFO) done; **final review: the front-row "near-graze" read is now finding I1 — unverified without mesh-normalizing `Sd`** |
| 2026-07-25 (final review) | register amendments (this fix wave) | `OPTIMALITY_CERTIFICATION.md` itself | recorded findings **I1** (`sdotMinRel` not mesh-normalized — discretization confound possible on refined/PSR-split meshes), **I2** (direction-check tangential residual not independent of `kktStat` by construction, and sign-blind to the ±p/‖p‖ branch), **I5** (`lamMassEndRel` one-sided endpoint-dual O(h) bias plausibly explains all three transversality misses); fixed M7 wording | none of I1/I2/I5 changes any PASS/FAIL verdict above; they are pre-promotion leads, see new §A6 |
| 2026-07-25 (Task 0, mesh-convergence-study) | **mapped terminal covector** (Hager 2000) — `foc_check.m` §(6), `rep.lamMassEndMapped`, now gated | earth 10/5/2.5 N, GTO→tulip flagship (warm re-solve, live model) | mapped transversality: 10N 1.242e-25, 5N 1.405e-17, **2.5N 2.268e-18**, tulip 1.192e-24 — all at KKT-stationarity floor, 14+ orders inside `tolTrans`=1e-3 | **RETRACTS** the "earth 2.5N is the single genuine first-order finding" verdict (I5 fix-pass block, A3 earth-2-body row) — see new "RETRACTION" block in §A6; I5 (`lamMassEndRel`) and the raw one-sided dual retained as reported-only companions; regression test `verify_common/tests/test_foc_terminal_covector.m` |
| 2026-08-26 (two-root ADJUDICATION, same day) | head-to-head flight of both roots on a common time base + symmetric cross-fly across environments | flagship min-time, both stored z8 artifacts | **ONE extremal, not two**: trajectory separation max 23.5 km / 28 d (median 1.3 km), periselene identical, λ_v0 angle 0.000°; cross-fly symmetric — each z8 lands 4.3 km only in its home arithmetic, 123–132 km in the other | the razor-basin reading is RETRACTED for this case: the 94 s spread is **z8 identifiability** (~1e-4 ND in t_f, ~100 km in arrival at ~40 revs), not a second basin; **certified direct t_f = 6.290694 STANDS**; Stage B rule: GTO-regime catalog entries must ship full ms junction states, not bare z8 (the refined-library lesson at mN scale); conjugate test on the alt is moot (same extremal) |
| 2026-08-26 (two-root discovery, flagship min-time) | library-move verification (`seed_from_z8`/`ms_tfmin` promotion; seed bit-identical, file identical to HEAD, in-session determinism 0.0) + `-batch` cross-run | GTO→tulip min-time K=60 probe | the SAME seed converges to two distinct extremals depending on MATLAB runtime: 6.290693961 (=certified direct) vs **6.290449990 — 94 s faster, periselene 6,039 km, ΔV 4.4663 vs 4.4665 km/s, tfMin-accepted** | razor-basin family grows (tenN-minfuel, DRO faster-basin, now flagship min-time); artifacts `mintime_ms_bvp_probe{,_altroot}.mat`; OPEN: adjudicate whether t_f,min moves — conjugate test + trajectory diff on the alt root; standing rule implied: near-min-time solves need basin disambiguation, one environment's root is not "the" root |
| 2026-08-26 (G1 ladder sweep, both MEE campaigns) | `verify_common/mee_residual` over every certified fuel row (7 rungs × 2 campaigns; warm re-solves; per-row try/catch batch) | earth 2-body MEE_M2_{10..0.1N} + CR3BP cr3bp_T{10..0.1}N_phi0_fuel — **14/14 rows ran**, incl. the previously-deferred 0.2/0.1 N (43–300 s each) | (1) switch attribution TOTAL ladder-wide (top-10 residual intervals are switch-straddling, fraction 1.00, one 0.80); (2) per-switch error SHRINKS with depth — RPmax 39.9→4.6 km — while >1 km interval count grows 57→1787 with switch count (**the 2026-08-25 growth prediction is REFUTED**); (3) 2-body interior machine-exact at every rung (1e-14..3e-13 ND); (4) CR3BP interior floor grows 6e-8→1.8e-6 ND (2.6→75 m) tracking t-row error — hypothesis: t-drift misphases the lunar term via phi0+nM·t (unestablished) | full table + raw .mat: `verify_common/doc/g1_sweep_results.md` / `g1_sweep_2026-08-26.mat`; rung-2 coverage now ladder-complete on both MEE campaigns (G2-gap: only DRO→tulip has an indirect cross-check); certification language updated: km-accurate AT switches, m-accurate elsewhere — quantitative justification for PSR-style refinement on both campaigns |
| 2026-08-26 (first CR3BP-GEO continuous-residual gate) | shared gate `verify_common/mee_residual` (MOVED from the 2-body campaign on its second consumer, migration rule; equivalence: 10 N 2-body numbers reproduce from the new home) + thin driver `cr3bp_residual_gate` (par+pert rebuilt from the artifact fingerprint, refresh_duals_cr3bp warm re-solve) | CR3BP-GEO certified lunar-aware 10 N row (`cr3bp_T10N_phi0_fuel`), drift 1.5e-11, defect 4.1e-15 | P max 9.46e-4 ND = **39.89 km at the SAME switch-straddling intervals as the 2-body row** (56/109/83/129/103) — the switch-interval concentration finding generalizes across the lunar perturbation; interior median 6.2e-8 ND (≈2.6 m) vs 2-body 2.2e-13 — a small nonzero floor from the richer dynamics; time-row 143 s unchanged | rung-2 coverage now: DRO→tulip, earth 2-body (10 N), earth CR3BP (10 N). Open: ELFO wrapper (Sundman/CasADi format), remaining rungs both MEE ladders |
| 2026-08-25 (flagship min-time via shared MS engine) | `costate_common/ms_bvp` (via `ms_tfmin`, K=60) seeded by flying the module's stored z8; gates: residual, flown arrival, tfMin acceptance, conjugate test | GTO→tulip min-time, 25 mN / 15 kg / ~40 rev (`GTO_tulip/indirect/min_time/mintime_ms_bvp_probe`) | **PASS**: ‖R‖=1.97e-11 (16 iters, 10.3 s), flown 4.3 km, accept |Δz|=5.1e-7 (ss floor at 6.3 ND), **conj PASS 0 crossings**; t_f = 6.290693961 ND — exactly the certified direct reference 6.290694 (the 2026-07-13 module root ‖R‖=4e-9 sat 47 s off at 6.290815: shallower root, same basin) | the shared engine closes in the flagship regime → GTO catalog (Stage B) unblocked; direct↔indirect agreement now closed on the flagship min-time anchor; finding: `gen_tulip_mintime` saves no duals — harvest leg needs dual capture (Stage B item) |
| 2026-08-25 (first MEE continuous-residual gate) | **`oc.local_residual`** (new oclib engine, TDD 6/6; `dro_residual` rerouted through it, equivalence 1.7e-18 on both stored refs, instrument harness 3/3) + new campaign wrapper `earth_elliptic_to_geo/direct/verify/mee_residual` | earth 2-body MEE 10 N certified row (warm re-solve, drift 3.6e-15, defect 6.3e-15) | P error: **median 2.2e-13 ND, max 9.5e-4 ND = 39.9 km**; 57/193 intervals >1 km — **every large residual sits on a throttle-switching interval** (linear trapezoid control vs interior switch); time-row max 1.0e-2 ND = 143 s on early full-burn intervals (O(h²) on 1/L̇, mesh widest in time); defect-to-truth gap ~1e11 | A1 row 10's "generalize to the other campaigns" now PARTIAL: engine in the library, DRO rerouted, MEE gated; CR3BP-GEO + ELFO wrappers open; the switch-interval concentration independently confirms the PSR switch-aware-refinement rationale |
| 2026-08-23 (catalog conjugate sweep) | `ms_conjugate_test` at catalog scale via `costate_common/conj_catalog_pass` (fly stored z8 → K=24 junction seed → 1-iteration ms re-solve with STMs → quotiented Jacobi dets; verdict only when re-solve reproduces the entry, measured dz 1e-12..1e-10) | ALL five shipped catalogs — 15,896 min-time entries (DRO 3,936 / HALO→tulip 3,980 / DPO 3,932 / L1→L2 1,952 / L2→L1 2,096) | **15,895 PASS / 1 FAIL / 0 unverified**, ~33 min wall (3 batched passes, ~0.14 s/entry). The FAIL: DPO τ=2 → Np=7, (sD,sA)=(2/3,2/3), 15 N, tf=0.115 ND, ΔV 5.12 km/s — one interior sign change, atFinal=0: **refuted as a local minimum** (an expensive branch; its own lower-thrust rungs pass) | every shipped entry now carries a second-order-necessary verdict (rung 3½ of the ladder); Part B §1 conjugate row updated; verdicts written into the canonical catalog .mats (backups `.bak_conj`); deliverable copies NOT yet re-shipped — Darin's zips predate the verdicts |
| 2026-08-23 (register correction, no new experiment) | Part B §1 conjugate-point row | `OPTIMALITY_CERTIFICATION.md` itself | the row claimed "not built; gated on an indirect solver" — **stale since 2026-08-08**: `costate_common/ms_conjugate_test` was built in migration #3 and runs in the golden-cells regression, `ms_tfmin` opt-in, and the `foc_check` advisory hook. The miss repeats the documented Part B failure mode ("an instrument that quietly works is easier to lose than one that loudly fails") — the working instrument had no FINDING block and the register row was never updated | Part B §1 row rewritten; the actionable gap is now stated correctly: catalog-scale application + schema storage, an Earth-campaign ms solution through the hook, and a fixed-t_f form (see `STATUS_AND_ROADMAP.md` §5E/§6) |
| 2026-08-14 (min-energy pilot) | first-order, indirect route: fixed-t_f `ms_minenergy` residual, generic single-shooting acceptance `ss_bvp_accept`, Hamiltonian first-integral check, direct-vs-indirect cost agreement, flown arrival | DRO→tulip 12×12 cells (2,5)/(6,8)/(1,2), fixed t_f = γ·t_f^min, γ ∈ {1.1,1.2,1.4}, MIN-ENERGY J=∫s² | **5/5 PASS all seven gates**: ms ‖R‖ ≤ 7e-11 (2–3 iters, K=12), accept |Δz| ≤ 3e-10, |ΔH| 3e-9..4e-8, indirect flight lands 0.000–0.053 km, J agree ~1e-5, m_f agree to 6 digits | first non-min-time coverage; **no second-order verdict for min-energy entries** (`ms_conjugate_test` is free-time-specific — a fixed-tf Jacobi test is an open item); single-shooting residual floor measured at ~6e-7 (with/without-STM flights) — acceptance must be judged by |Δz| at that floor, never at the ms 1e-10 (record in `DRO_tulip/FINDINGS.md`) |
| 2026-08-27..29 (GTO→tulip catalog fleet, Stage B) | the standing three-gate pipeline at catalog scale (multiple-shooting `ms_tfmin` residual, flown-arrival <100 km, `pumpkyn.cr3bp.tfMin` acceptance |Δz|<1e-6), run over 16 sheets via the warm-recipe two-point top-rung multistart (thrLock, tf0=0.30 primary + scoped cold tf0=4.0 mop-up), following the 2026-08-25 flagship-anchor row above onto the Darin-standard 1–15 N (few-rev) regime instead of the 25 mN flagship | GTO (4 orientations {0,90,180,270}° × 4 Np {5,7,9,12}) → tulip, 12×6 phasing grid, rungs `[15 12 10 7 5]` N | **2,625 entries accepted, 840/1,152 phase pairs (73%)**, 136 full 5-rung ladders, all three-gate certified (independently re-verified by a 3-entry sheet-local replay audit, `audit_gto_entries.m`, 2/2 hard gates on all 3 rows); coverage by orientation 240/207/156/237 of 288 (83/72/54/82%) — orientDeg=180° (apogee toward the Moon) is a genuine cold-basin hard corner (53–56% per sheet, 1–2 full ladders), not a mesh artifact (Diagnostic A: escalating mesh resolution made the original failure *worse*) | Stage B (§6 step 6 of `STATUS_AND_ROADMAP.md`) **DONE**; first costate catalog built on a non-periodic departure family (`'gto'` pseudo-family, Task 1); the 3–1 N legs of the standard ladder produced **zero** entries by both recipes — a measured closure wall, split off as the deep-rung investigation (open item, `STATUS_AND_ROADMAP.md` §5.B); catalog README `GTO_tulip/catalog/README.md` |
| 2026-08-29 (GTO→tulip catalog conjugate sweep) | `ms_conjugate_test` at catalog scale via `costate_common/conj_catalog_pass` (same instrument as the 2026-08-23 catalog sweep row: fly stored z8 → K=24 junction seed → 1-iteration ms re-solve with STMs → quotiented Jacobi dets; verdict only when re-solve reproduces the entry) | the GTO→tulip catalog above — 2,625 min-time entries, 16 sheets | **2,625 PASS / 0 FAIL / 0 unverified** (100%) — every accepted entry has zero interior conjugate-point crossings (`conj_ncross = 0` on all 2,625); re-solve fidelity 9.6e-16..1.1e-8, comfortably inside the 1e-6 honesty gate; 610 s wall (0.232 s/entry). **Stronger than the DPO precedent** (3,931 PASS / 1 FAIL, 99.97%) — two plausible, unconfirmed structural differences recorded: this fleet's warm-recipe multistart (vs DPO's cold-start ladder) and its 5 N thrust floor (vs DPO's 1 N floor, where deep many-rev rungs carry more junction structure) | every GTO→tulip catalog entry now carries a second-order-necessary verdict (rung 3½ of the ladder), matching the other five shipped catalogs; verdicts stored in the catalog itself (`conj_pass`/`conj_ncross`/`conj_atfinal` per sheet + `conj_test` provenance, schema-validated `{}`); the zero-refutation contrast with DPO is recorded as an open hypothesis, not proof — a genuine test needs a comparable deep-rung GTO sweep, which does not exist yet |

---

## 7. Pointers

| file | role |
|---|---|
| `verify_common/doc/first_order_checks.pdf` (`.tex`) | **the explainer for Part A's generic layer**: every first-order check stated formally + intuitively, with per-campaign run commands, the annotated report block, thresholds table, and the I1/I2/I5 caveats carried verbatim |
| `earth_elliptic_to_geo/process/DESIGN_sosc.md` §11–12 | NLP SOSC method evolution, `eig` vs `ldl`, threshold rationale |
| `earth_elliptic_to_geo/process/PLAN_sosc.md` | original build plan |
| `earth_elliptic_to_geo/direct/verify/sosc/` | the NLP-level implementation |
| `GTO_tulip/direct/certify/second_order_nlp.m` | KKT-inertia test + M2 finding |
| `GTO_tulip/direct/certify/switch_hessian.m` | Maurer test + M3 finding + the STM fix spec |
| `earth_elliptic_to_geo/process/LESSONS_DUAL_EXTRACTION.md` | why the first-order side is now trustworthy |
| `earth_elliptic_to_geo_CR3BP/TODO.md` | CR3BP tier plan (supersede with §1 before acting on it) |
