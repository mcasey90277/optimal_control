# Survey: what it takes to PROVE an extremal is a local minimum -- and what we actually test

Written 2026-09-05 after the step-5 review. Companion to
`../OPTIMALITY_CERTIFICATION.md` (the running register of instruments and
verdicts). That document records *what ran*; this one asks the question
the other way round: **for each problem we solve, what is the complete
chain of conditions a rigorous certificate needs, which link does each
existing test supply, and where is the chain broken?**

Scope: the CR3BP low-thrust transfers (DRO/halo/DPO -> tulip, halo <-> halo,
GTO -> tulip/ELFO) and the Earth elliptic -> GEO campaigns, for the three
objectives min-time, min-energy (fixed t_f) and min-fuel.

---

## 0. The problem class, stated once

Autonomous, control-affine Bolza problem on the 7-state x = (r, v, m):

    xdot = f0(x) + (s T/m) alpha,   mdot = -s T/c,
    u = (alpha, s),  alpha in S^2 (unit sphere),  s in [0, 1]
    J = int_0^tf L(s) dt,   L = 1 (min-time, tf free) | s (min-fuel) | s^2 (min-energy)
                            | (1-p)s + p s^2 (Bertrand-Epenoy, p in (0,1]) | huber_p(s)
    r(0), v(0), m(0) fixed;  r(tf), v(tf) fixed;  m(tf) FREE.

Normalized Hamiltonian (lambda_0 = 1): H = L(s) + lambda_r'v + lambda_v'(g + h + s T alpha/m) - lambda_m s T/c.
Switch quantity Q = T(|lambda_v|/m + lambda_m/c); the fuel law is s = 1 if Q > 1, 0 if Q < 1.

Three structurally different cases fall out, and the certificate is
different for each:

| case | control structure | what "regular" means | second-order theory that applies |
|---|---|---|---|
| **A. min-time** | s = 1 always; alpha = -lambda_v/\|lambda_v\| smooth on S^2 | \|lambda_v\| > 0 on [0,tf] (no direction singularity) | smooth case: Legendre + Jacobi (conjugate/focal points), free-time quotient -- Bonnard-Caillau-Trelat 2007, Agrachev-Sachkov ch. 20-21 |
| **B. smoothed fuel, p > 0** | interior arcs (Q in (1-p, 1+p)) where s is C^0 in Q, plus **bound-active arcs** s = 0 or 1 | strict complementarity: Q stays off the band edges except at isolated transversal crossings | smooth theory on interior arcs; bound-active arcs behave as bang arcs -> **mixed** bang / interior second-order conditions (Maurer-Osmolovskii-Buskens) |
| **C. bang-bang fuel, p = 0** | s in {0,1}, finitely many switches at Sdot != 0 | regular switching: Sdot(t_i) != 0, no singular arc | Maurer-Osmolovskii-Osmolovskii: positive-definite switching-time Hessian + regular switching + primer Legendre |

Every claim below is indexed to one of A/B/C.

---

## 1. Level 1 -- the extremal (Pontryagin necessary conditions)

The full list, not the usual short one. "IND" = the indirect pipeline
(`ms_bvp`/`ms_tfmin`/`ms_minenergy`/`ms_minfuel` + `tfMin` acceptance);
"DIR" = direct collocation + `foc_check`.

| # | condition | what a proof needs | IND: what we have | DIR: what we have | gap |
|---|---|---|---|---|---|
| 1.1 | state ODE holds on [0,tf] | continuous residual, not discrete defects | by construction: ode113 RelTol 1e-10; ms junction residual < 1e-10 | discrete defects ~1e-14; continuous LTE measured ONCE (`dro_residual`: 1441 km at N=800 trapezoid, 3 m at HS N=1600) | DIR: G1 accuracy gate exists but is not routine; **defect != accuracy** |
| 1.2 | costate ODE holds | same | by construction (the costates are integrated) | KKT stationarity = discrete adjoint (tautological); validated against IND once (`costate_compare`: 6e-4 deg) | DIR-only campaigns (earth, ELFO) have **no independent costate witness** (G1) |
| 1.3 | pointwise minimum condition u* = argmin_U H (not dH/du = 0) | alpha = -lambda_v/\|lambda_v\|; s from Q by the exact argmin of L(s) - sQ | by construction in the field (`cr3bp_*_pmp`; `test_minfuel_pmp` checks the argmin, both families) | primer alignment + sign law (`foc_check`), 0.000 deg / 100% on earth, PASS on tulip/ELFO fuel rows | DIR: the two tulip costate reads DISAGREE on the flagship (raw dual vs LS), never adjudicated |
| 1.4 | transversality, free mass: lambda_m(tf) = 0 | exact | terminal condition of the BVP (residual < 1e-10) | Hager-mapped terminal covector, PASS all campaigns | none |
| 1.5 | free-tf condition: H(tf) = 0 (autonomous => H == 0) | value test | min-time: terminal condition H(tf) = 0 in `tfMin`/`ms_tfmin` (exact) | dual form only (lambda_t(tf) = -1 confirmed on ELFO; lamTimeCoV 5.7% OPEN) | DIR: literal H(tf)=0 value test not built (G4) |
| 1.6 | fixed-tf: H == const (autonomous) | constancy along the arc | `Hdrift` at junctions (1e-11..1e-12) | earth: CoV 1e-9; CR3BP-earth non-autonomous: nothing (should test dH/dt = dH/dt_explicit, G6) | IND: junction-only sampling (interior not checked) |
| 1.7 | normality: lambda_0 != 0 | exclude abnormal extremals | normalized lambda_0 = 1 in every field; the min-time catalog's scale invariance quotient assumes normality | lambda_0 == 1 by construction | **never tested anywhere** (G7). For min-time an abnormal extremal is the tfMin BVP with the "1 +" dropped; a one-line check is whether the same lambda(0) direction solves it -- it will not if \|lambda\| is finite and H(tf) = 0 is met, so acceptance IMPLIES normality; state that, do not assume it |
| 1.8 | no singular arc | Q != 1 (or S != 0) except at isolated times | not computed on IND min-fuel walks at finite p (the band replaces the switch) | H2 run-of-\|S\|~0 check (earth) | IND: add a Q-band residence statistic per entry |
| 1.9 | regular switching Sdot != 0 (case C) / transversal band crossings (case B) | at each switch | not computed (finite p: the crossing of Q = 1-p, 1+p should be transversal; `cr3bp_minfuel_prop`'s grazing assert covers only huber) | `sdotMinRel` advisory line, mesh-normalized since 07-25 | IND: compute Qdot at band crossings, gate it |
| 1.10 | Weierstrass-Erdmann (lambda, H continuous at switches) | -- | by construction (control-affine, no state constraints) | same | none -- state it |
| 1.11 | **independence of the witness** | two methods agree | direct -> covector harvest -> ms -> **pumpkyn tfMin acceptance unchanged (\|dz\| < 1e-6)** -- a genuinely independent second solver | DIR-only: self-consistency | this is the strongest extremal evidence in the program and exists ONLY for min-time; min-energy/min-fuel use `ss_bvp_accept` (same field, K = 1) -- same code, not a second method |

**Verdict at level 1.** The catalog pipeline (case A, 18,360 entries)
establishes extremality about as well as it can be established
numerically: two methods, one of them foreign (pumpkyn), agree on the
costates to |dz| ~ 1e-9, with every terminal condition exact. The
fixed-tf line (cases B) is extremal by construction of the indirect solve
but has no foreign witness, no singular-arc statistic and no
crossing-transversality gate. The direct-only campaigns are extremal *on
their transcription* with one measured, large, continuous-accuracy gap.

---

## 2. Level 2 -- local minimality

### 2.1 What the theory actually gives

**Case A (smooth, min-time, free tf).** Bonnard-Caillau-Trelat (2007) /
Agrachev-Sachkov: for a *normal* extremal with the *strong Legendre
condition* (H_uu > 0 on the tangent space of the control set -- here
H_alpha alpha restricted to T S^2 = (sT/m)|lambda_v| I, positive iff
|lambda_v| > 0) and **no conjugate time on (0, tf]** (for free tf: the
quotient by the flow direction, exactly `ms_conjugate_test`'s free-time
form), the extremal is a **strict local minimum in the C^0 topology** on
trajectories with the same endpoints. Conversely a conjugate time in
(0, tf) refutes minimality (Jacobi necessary condition). So for min-time
the conjugate test is not just necessary -- with two cheap additions it is
the *sufficient* certificate:

  (i) |lambda_v(t)| bounded away from zero on [0, tf]  -- NOT GATED today;
  (ii) no conjugate time on the whole of (0, tf] at CONTINUOUS resolution,
       not at K junction samples  -- junction resolution today
       (12-48 samples; a conjugate PAIR inside one segment hides);
  (iii) normality  -- implied by tfMin acceptance (1.7), never stated.

  **Audit required before any of this is called a certificate (Astra
  review #2, 2026-09-06):** (1) the theorem must be applied with our
  ACTUAL endpoint manifold -- r, v fixed, mass FREE with lambda_m(tf) = 0
  and H(tf) = 0; for all-burn min-time mass is eliminable
  (m = m0 - Tt/c, delta m_f = -(T/c) delta t_f), which is why the
  6-state block rows 1:6 / cols 8:13 is the right object, but the
  reduced exponential map / accessory BVP whose kernel is the admissible
  null space of the second variation must be written down, not assumed
  from "the free-time quotient as in HamPath"; (2) tfMin acceptance
  gives a normal LIFT -- it does not exclude an abnormal lift of the same
  trajectory, and if the theorem needs strict/strong normality or a
  corank-one endpoint map, that is a separate hypothesis; (3) the
  costate-scaling quotient is an invariance of the CONTROL
  (Phi_xl * lam0 = 0 exactly), but with lambda_0 = 1 and H = 0 it is not
  a symmetry of the BVP -- its equivalence to the theorem's construction
  must be shown; (4) "no conjugate time at integrator resolution" is a
  numerical ASSESSMENT (needs between-sample bounds, even-multiplicity
  handling, endpoint focal degeneracy), not a proof; a rigorous
  certificate needs validated enclosures; (5) the catalog's "min-time" is
  the ALL-BURN problem -- optimality against throttle variations is
  assumed, not tested. Items (i)-(ii) remain the right cheap additions;
  they upgrade the assessment, they do not close the proof.

**Case B (smoothed fuel, p > 0).** On interior arcs H_ss = 2p > 0 (eps) or
1/p (huber) -- strong Legendre holds. **Corrected 2026-09-06 (Astra review
#2):** for a CONTINUOUS clipped law (eps, huberc) the minimized field is
continuous and piecewise smooth, so at a transversal branch boundary
F+ = F- and the saltation matrix is the identity: the branchwise-AD STM
`ms_conjugate_test` consumes IS the derivative of the full flow, not of an
active-set-frozen surrogate. Its verdict is the ordinary Jacobi necessary
condition under the regularity hypotheses (strict sign of Q - 1 inside
active arcs, isolated transversal crossings, positive interior throttle
curvature). "Frozen active set" is legitimate only as a statement about
the critical cone (delta s = 0 a.e. on bound-active arcs), not as a
restriction on perturbed trajectories. Sufficiency for continuous
saturation of a strictly convex control cost can be analysed with
critical-cone / Riccati methods; a separate bang-bang switching-time
Hessian is NOT automatically required. The original Huber family (jump at
Q = 1, nonunique minimizer there) is a genuinely different, discontinuous
case and does need the bang-bang machinery. The transversality gate (1.9)
still does not exist. What we can claim for the 7 shipped min-fuel entries after
the 09-05 fix: *no interior conjugate point at junction resolution for the
frozen-active-set accessory problem* -- a necessary condition, honestly
labelled in the catalog as such.

**Case C (bang-bang, p = 0).** The register's M1 argument -- "the fuel cost
is linear in s, so no NLP-level Hessian can reach strictness" -- is
**disputed (Astra review #2)**: the Lagrangian's reduced Hessian on the
kernel of the endpoint-constraint Jacobian includes the nonlinear dynamics
and need not vanish because L(s) = s. The 270 flat directions measured at
earth 10 N stand as a measurement of THAT transcription; the impossibility
theorem inferred from it does not. (Register M2 -- strict complementarity
fails at the throttle bounds in an eps -> 0 approach -- is a separate,
valid point.) Strictness lives in the switch times:
Maurer-Osmolovskii's sufficient condition is (i) positive-definite Hessian
of the reduced problem over the switch times, (ii) Sdot != 0 at every
switch, (iii) primer/Legendre optimality of alpha. (ii) and (iii) are
level-1 lines we already report; (i) is BLOCKED in its forward-flow form
(M3) and specified but unbuilt in its multiple-shooting/STM form
(register Part B section 5). IPOPT native inertia certifies a *weak* local
minimum on some rows -- the honest ceiling of that route.

### 2.2 What we have, mapped

| instrument | case | level | what it proves | coverage | status after 09-05 |
|---|---|---|---|---|---|
| `ms_conjugate_test` free-time form | A | 2-necessary (and, with (i)-(iii) above, sufficient) | no conjugate time at junction resolution; quotient handles costate scaling + time reparametrization | 18,249 catalog entries (K=24), 61 refuted | corrected (t_f sampled, equilibrated, coast skip); catalogs NOT yet re-swept with the corrected instrument; golden 20/20 identical |
| `ms_conjugate_test` fixed-time form, rows 1:7 | B | 2-necessary (frozen active set) | no interior conjugate point | 15 fixed-tf records | corrected 09-05 (was the terminal block; was counting coast zeros); 15/15 PASS |
| IPOPT native inertia delta_w | B, C (direct) | weak local min of the barrier NLP | no negative curvature on the barrier system | tulip 12/17 eps=0 rows; earth 3/9; CR3BP-earth 4/4; ELFO 2/3 | unchanged |
| NLP reduced-Hessian SOSC | C (direct) | weak min (270 flat directions) | -- | earth 10 N only | ceiling reached (M1); do not redo |
| KKT-inertia SSOSC | C (direct) | -- | inapplicable (M2) | tulip | ruled out |
| switching-time Hessian (forward flow) | C | 2-sufficient (Maurer) | -- | tulip | BLOCKED (M3); STM/multiple-shooting form is the register's stated next build |
| `ss_bvp_accept` / `tfMin` acceptance | A, B | level 1 (extremal), not level 2 | root of the single-shooting equations | all | tfMin: foreign witness (A); ss_bvp_accept: same field (B) |
| golden cells | -- | regression of the instruments | verdicts pinned on 4 cells | -- | 20/20 through every 09-05 change |

### 2.3 The honest one-line status per objective

- **Min-time (A):** extremal to the strongest standard we have (foreign
  witness); local minimality at the level of "no conjugate time at
  junction resolution" on 18,249 entries -- *two small additions away from
  the Bonnard-Caillau-Trelat sufficient certificate* (section 3, item 1).
- **Min-energy / smoothed min-fuel (B):** extremal by construction of the
  indirect solve, no foreign witness; frozen-active-set Jacobi necessary
  condition passes on 15/15; no sufficiency instrument.
- **Bang-bang min-fuel (C):** extremal on the transcription with a
  measured continuous-accuracy gap; weak local minimum on a subset of rows;
  strictness unreachable without the switching-time Hessian.

---

## 3. The gaps, ordered by (certificate gained) / (cost)

1. **Finish the min-time sufficiency certificate (A).** Three additions to
   the catalog pipeline, all cheap, all reusing `ms_bvp`'s STMs:
   (a) gate min |lambda_v| over the arc (strong Legendre) -- one number per
   entry, already computable from the junction states + flown arc;
   (b) sample the conjugate determinant at the integrator's own steps, not
   only at junctions -- propagate the variational equations per segment
   (already done for the STM) and evaluate det(Phi_xl(t,0) P, f) along the
   dense output; a conjugate pair inside a segment can then no longer hide;
   (c) state normality as implied by tfMin acceptance (1.7).
   With (a)-(c) and the corrected instrument, every PASS entry carries the
   BCT sufficient conditions for a strict C^0-local minimum. **This is the
   highest-value item in the program**: it turns 18k "necessary-condition
   pass" entries into certified local minima with no new theory.
   Prerequisite: re-sweep the catalogs on the corrected instrument (the
   t_f sample is new; it can only add refutations).

2. **Backward Riccati / focal-point test for the free-mass terminal
   manifold (B, interior arcs).** Integrate the Riccati equation of the
   accessory problem backward from tf with the terminal condition that
   encodes delta(r,v)(tf) = 0, delta m(tf) free; an escape time in [0, tf)
   is a focal point. Same STMs, one more 7x7 sweep. Together with 1(a)-(b)
   this gives the smooth-case sufficient condition on any p > 0 extremal
   whose active set is empty -- and, for the shipped entries, the interior-
   arc half of the mixed condition.

3. **Crossing-transversality and band-residence gates for B (level 1.8,
   1.9).** Qdot at every entry/exit of the band (1-p, 1+p), and the fraction
   of the arc spent inside the band. Without these, the frozen-active-set
   Jacobi verdict rests on an unstated assumption. Cheap: the flown arc is
   already computed for `mf`.

4. **The switching-time Hessian in multiple-shooting/STM form (C).**
   Register Part B section 5 already specifies it; the ms STMs and the new
   saltation machinery (`cr3bp_minfuel_prop`, huber branch) are exactly the
   ingredients: d Psi / d sigma_i = Phi(sigma_f, sigma_i) [f_burn - f_coast]
   is a saltation vector. This is the only route to a strict certificate
   for the GTO -> tulip / ELFO / earth bang-bang fronts.

5. **A foreign witness for B.** `ss_bvp_accept` re-solves with the same
   field. Options: (i) pumpkyn has no min-fuel twin; (ii) MfMax (built,
   `earth_elliptic_to_geo/indirect/mfmax`) is a genuinely foreign min-fuel
   indirect solver for the two-body case -- not CR3BP; (iii) the direct
   collocation solution IS the foreign witness for the fixed-tf line
   (`minenergy_pilot` compares direct vs ms m_f) -- promote that comparison
   to a stored, gated per-entry number.

6. **Direct-campaign hygiene (level 1):** routine continuous-accuracy gate
   (1.1), continuous adjoint residual (1.2), H(tf) = 0 value test (1.5),
   dH/dt = dH/dt_explicit for the non-autonomous Earth-CR3BP (1.6),
   adjudicate the tulip raw-dual vs LS disagreement (1.3).

---

## 4. What the 09-05 review taught about the instruments themselves

- A second-order instrument validated on a *scalar* analytic case
  (`test_conj_fixedtf`'s LQ pi-point) can be wrong on the block it monitors
  and still pass -- the fixture cannot distinguish blocks. Fixtures must be
  at least 2-state and must include the degenerate structures the real
  problem has (initial coast, saturated arc, final segment).
- "Structurally zero" and "focal point" both read as det = 0. Any
  determinant-based Jacobi test on a problem with bound-active arcs needs
  a rank diagnostic and a rule for where testing may begin.
- A discontinuous control law needs a saltation-correct STM before any
  claim -- first- or second-order -- is made from a Newton residual. The
  Huber "structural knockout" was a wrong Jacobian; a "conjugate point"
  in a Huber arc would have been too.
- A sufficiency claim is a claim about the *whole* interval, so the
  instrument's resolution (junction samples) is part of the claim.
  Item 3.1(b) exists because of this.

## 5. Pointers

- `../OPTIMALITY_CERTIFICATION.md` -- the register (instruments, verdicts,
  experiment log); Part B section 5 is the switching-time Hessian decision.
- `../costate_common/ms_conjugate_test.m` -- the instrument, header = math.
- `../DRO_tulip/reviews/minfuel_code_review_2026-09-05.md`, FINDINGS 22.
- Bonnard, Caillau, Trelat, "Second order optimality conditions in the
  smooth case and applications in optimal control," ESAIM: COCV 13 (2007).
- Maurer & Osmolovskii, "Second order sufficient conditions for time-optimal
  bang-bang control," SIAM J. Control Optim. 42 (2004); Osmolovskii &
  Maurer, *Applications to Regular and Bang-Bang Control*, SIAM 2012.
- Agrachev & Sachkov, *Control Theory from the Geometric Viewpoint*,
  ch. 20-21 (Jacobi/conjugate points, Legendre conditions).
- Bryson & Ho, *Applied Optimal Control*, sec. 6.3 (focal points, backward
  Riccati sweep with a terminal manifold).
