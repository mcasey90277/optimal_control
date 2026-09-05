# Step-5 min-fuel line -- external code review, adjudicated (2026-09-05)

**Scope.** The fixed-t_f objective line closed 2026-09-02 (FINDINGS 18-21):
`costate_common/{cr3bp_minfuel_pmp,cr3bp_minfuel_prop,cr3bp_minenergy_pmp,
cr3bp_minenergy_prop,ms_bvp,ms_conjugate_test,catalog_schema}.m`,
`DRO_tulip/{run_minfuel_race,run_minfuel_grid,build_minfuel_catalog}.m`,
`DRO_tulip/indirect/{ms_minfuel,ms_minenergy}.m`, and the three tests
(`test_minfuel_pmp`, `test_conj_fixedtf`, `test_catalog_schema_v3`).
15 files, 2,343 lines, inlined (no `-c`) per the GPT-slot rule.

**Reviewers.** GPT-5.6-sol (full 121 KB bundle via crush, 6 min,
`minfuel_code_sol.md`); GPT-6 Astra (the requested reviewer: 53 KB
math-core bundle via the raw OpenRouter API after crush wedged twice,
146 s, `minfuel_code_astra_core_rawapi.md`; plus a two-file crush probe,
`minfuel_code_astra_min.md`); host Claude (independent read + numerical
verification against the shipped `.mat` artifacts in MATLAB). The Astra
invocation problem is diagnosed in the last section.

**Prompt.** `minfuel_code_prompt.txt`, with FINDINGS 18-21 pasted in as
"hypotheses to attack". The prompt deliberately planted a wrong prior
(tanh/logistic smoothing) with an instruction to VERIFY; sol verified and
reported the actual polynomial form -- a good sign -- but then recommended
"implementing the advertised tanh law", which is not a defect. Discount that
item.

---

## Verdict in one paragraph

The first-order mathematics of the `eps` (Bertrand-Epenoy) line that
shipped is sound: costate equations, switch quantity, throttle law and
envelope-theorem Jacobian are all correct, and the catalog's own
reconstruction recipe reflies every entry to <= 0.15 km. Three defects
matter. (1) **The shipped catalog mislabels its propulsion**:
`thruster.isp_s = 1710` but `c_nd = 8.673746`, which is Isp = 900 s
exactly; any consumer recomputing c from `isp_s` (as every other catalog
builder does) gets delta-V wrong by 1.90x. (2) **The Huber arm's STM omits
the saltation matrix** across its Q = 1 jump, so FINDINGS 18's *mechanism*
("switch-jump residual floor is structural") is not established; the
*conclusion* (walk eps) stands since the eps law is continuous. (3) **The
fixed-t_f conjugate test monitors the terminal shooting block, not the
interior Jacobi block** (`[1:6 14]` instead of `1:7`), so the 14 verdicts
of FINDINGS 20 and the 7 `conj_pass` flags in the catalog are unmeasured
and must be re-run -- and when they are, the single refutation will very
probably turn out to be an **initial-coast artifact** (the refuted record
is the only one still coasting at the first det sample), which rewrites
the mechanism story of FINDINGS 20-21. A further defect, the `atFinal`
rule demoting a strictly interior crossing, is real but measured latent
(0 firings in 14 fixed-t_f verdicts and 4,405 min-time entries).

---

## P0 -- fix before the catalog is shipped to anyone

### P0.1 Propulsion label contradicts the physics (host finding; sol missed it)

`DRO_tulip/build_minfuel_catalog.m:112`
```matlab
'thruster', struct('isp_s', 1710, 'm0_kg', 150, 'c_nd', rec1.c, ...
```
`rec1.c = 8.67374595` comes from `minenergy_pilot.mat`, whose substrate is the
12x12 phasing torus built at **0.07 N / Isp 900 s / 150 kg**
(`DRO_tulip/direct/sweep_phasing_direct.m:55`,
`DRO_tulip/indirect/sweep_phasing.m:90`). The `1710` is the min-time
catalogs' Isp, copied by hand.

Measured on the shipped `.mat`:
```
isp_s=1710  c_nd=8.673746
c_nd EXPECTED from Isp=1710 : 16.480117   (stored 8.673746, ratio 0.5263)
=> stored c_nd corresponds to Isp = 900.0 s
```
Consequences: `catalog_schema('derive','deltav_from_mf')` is internally
right (it uses `c_nd`), but a consumer who trusts `isp_s` -- which is what
`build_costate_catalog_family.m:77` does -- gets delta-V 1.90x too large
(e.g. entry (1,2)@1.10: 0.919 km/s true vs 1.746 km/s). `m0_kg = 150` is
correct (Tmax_nd 0.175642 inverts to 0.070 N at 150 kg). `thrustN` is not
stored at all.

**Fix.** Set `isp_s` from the source (900), add `thrustN` (0.07), and port
the consistency assert from `build_costate_catalog_family.m:82` into the v3
builder AND into `catalog_schema('validate')`:
`abs(c_nd - (isp_s/tStar)*g0nd) < 1e-9*c_nd`. Then rebuild the catalog and
update FINDINGS 19/21 and the deliverable README wherever "Isp 1710" is
implied for the min-fuel set.

### P0.2 Huber STM omits the saltation matrix (sol; corroborated by host)

`cr3bp_minfuel_pmp.m:126-128` makes `s` jump from `p*Q` to `1` at Q = 1;
`cr3bp_minfuel_prop.m:84-86` propagates `Phi_dot = A*Phi` straight through
the jump. The correct STM across a state-dependent discontinuity surface
g(y) = Q(y) - 1 = 0 is
```
Phi+ = [ I + (F+ - F-) (grad g)^T / ((grad g)^T F-) ] Phi-
```
Without it the shooting Jacobian is wrong on every Huber segment that
crosses a switch, which is exactly where FINDINGS 18 measured the
"switch-jump residual floor". So the floor may be Newton fighting a wrong
Jacobian, not a structural property of the penalty. Host had reached the
same worry from a different direction: ode113 is a variable-order multistep
method and carries history across the jump, so its step rejection is the
weakest possible handling of a discontinuity.

**Scope.** Huber only. The `eps` throttle is continuous (clipped affine),
so Phi is continuous and no saltation is needed; the shipped eps records
are unaffected.

**Fix.** Either (a) rerun the race arm with event-located switches, split
propagation, and the saltation update -- and only then keep or retract the
"structurally unsuited" language in FINDINGS 18 -- or (b) soften FINDINGS 18
to "eps wins in continuation practice; the Huber residual floor is
unexplained pending a saltation-correct Jacobian". (b) is cheap and honest;
(a) is the experiment.

### P0.3 The fixed-t_f conjugate test monitors the wrong interior block (sol; settled by Astra; host reversed)

`ms_minfuel.m:70-72`, `ms_minenergy.m:93-95`, `build_minfuel_catalog.m:123`
use `stateRows = [1:6 14]`, i.e. the 7x7 block
`Phi([r v lam_m], [lam_r lam_v lam_m])`. That block is the **terminal
shooting Jacobian** and its nonsingularity at t_f is exactly what ms_bvp
convergence already certifies. It is NOT the interior Jacobi test.

Why (Astra's argument, host re-derived): a classical conjugate point tau
in (0, t_f) is a nontrivial Jacobi field with delta x(0) = 0 and
**delta x(tau) = 0 for the full state, mass included**. Such a field,
extended by zero on (tau, t_f], is an admissible variation for the
original problem (continuous state, satisfies delta(r,v)(t_f) = 0, and
delta m(t_f) = 0 is allowed since mass is free) with zero second
variation -- that is the refutation. The condition delta x(tau) = 0 is a
rank drop of **`Phi(1:7, 8:14)`**. The code's field instead has
delta(r,v)(tau) = 0 and delta lam_m(tau) = 0 with delta m(tau) free;
extending it by zero makes delta m discontinuous, so it is not admissible
and its singularity refutes nothing. Conversely the code's block is
generically nonsingular at a true conjugate point (delta lam_m(tau) != 0
there), so the current test can MISS real ones. The mixed block is
justified only AT t_f, where delta lam_m(t_f) = 0 is the natural condition
of the free-mass terminal manifold.

`test_conj_fixedtf` did not catch this because its LQ fixture is scalar
and cannot distinguish the two blocks (sol's point).

**Impact.** All 14 fixed-t_f verdicts in FINDINGS 20, and the 7
`conj_pass = 1` values in the shipped min-fuel catalog, were produced by
the wrong interior matrix. They are not known to be wrong; they are
**unmeasured**. The (2,5)@1.4 refutation and its gamma-continuation rescue
(FINDINGS 21) may or may not survive.

**Fix.** (1) `stateRows = 1:7, costateCols = 8:14` in both drivers.
(2) Expect small dets near t = 0: the state block is rank <= 3 at leading
order in t (only acc and mdot see lam), so scale/equilibrate before the
sign test and start sampling past the first junction. (3) Re-run the
14-verdict sweep and rebuild the catalog. (4) Add a 2-state fixture to
`test_conj_fixedtf` that distinguishes `Phi(1:7,8:14)` from the mixed
block. (5) Optionally keep the mixed block as a separate *terminal*
second-variation check, which is what it actually is. The min-time
catalogs are unaffected: there the control is lam_m-independent and mass
is a known function of time, so rows 1:6 / cols 8:13 is right.

### P0.4 The one conjugate refutation is probably an initial-coast artifact (Astra; measured by host)

Astra: on a strict coast s = 0 in a neighbourhood of the extremal, so the
state dynamics do not depend on the costates and `Phi(1:6, 8:14) == 0`
identically through the initial coast (and `Phi(7, 8:14) == 0` too, since
mdot = 0). Any det sampled inside that interval is structurally zero --
and `ms_conjugate_test.m:133` counts every exact-zero sample as a focal
point, while a roundoff-signed near-zero can flip against the first
burning sample. Neither is a conjugate point.

Measured on the shipped records (throttle at the first three junctions;
the first det sample is at t_2):
```
cell   gamma  p_floor    s(t_1)  s(t_2)  s(t_3)   fuel verdict
(2,5)  1.20   0.00168     0       1       0       pass=1 nCross=0
(6,8)  1.20   0.00128     1       1       0       pass=1 nCross=0
(1,2)  1.20   0.00100     1       1       1       pass=1 nCross=0
(2,5)  1.10   0.00179     1       1       0       pass=1 nCross=0
(2,5)  1.40   0.00500     0       0       1       pass=0 nCross=1   <--
(6,8)  1.10   0.00124     0       1       1       pass=1 nCross=0
(1,2)  1.10   0.00341     1       0       0       pass=1 nCross=0
```
Three records start on a coast; the ONLY refuted one is the only one
still coasting at t_2. The two other coast-starters had regained rank by
t_2 and passed. This is a strong hypothesis, not a proof (the det samples
were not stored); the P0.3 re-run settles it for free by logging `dets`.

**Implications.** FINDINGS 20's "the one refutation explains the
gamma-anomaly" and FINDINGS 21's rescue narrative are probably wrong in
mechanism: the non-monotone m_f in gamma was a basin effect (the rescue
did find a better optimum, 0.949005 vs 0.944025, so it was not wasted),
but the conjugate refutation that pointed at that record was most likely
a coast artifact. **Fix.** In the re-run: (a) skip samples where the
state block has not attained full rank (or start sampling after the first
burn junction); (b) return `indeterminate`, not a crossing, for an
exact-zero sample; (c) store `dets` with every verdict. Note this bites
HARDER after the P0.3 block fix (rows 1:7 add the mass row, which is
also zero on a coast).

---

## P1 -- real defects, measured latent on the shipped verdicts

### P1.1 `atFinal` demotes a strictly interior crossing (sol; measured by host)

`ms_conjugate_test.m:116-133`. `dets` has K-1 samples at
`tGrid(2:K)`; the last is at t_K, not t_f. A sign change on the bracket
(t_{K-1}, t_K) is strictly interior, yet line 132 classifies it `atFinal`
and line 133 subtracts it from `nCrossings`, so `pass` can be true on a
refuted extremal. Measured for the fixed-t_f runs (K = 12):
```
last sampled time = 4.022258   t_f = 4.387918
unmonitored tail  = 0.365660 ND = 8.33% of the transfer
```
Also: the final segment's STM is never chained (loop ends at K-1), so a
conjugate point anywhere in the last 1/K of the transfer is invisible. The
header (lines 33-40) documents "endpoint not sampled", which understates
both effects.

**Measured impact: none so far.**
```
fixed-tf verdict sweep : atFinal fired on 0 of 14
DRO->tulip conj sweep  : atFinal fired on 0 of 4,405 tested (8 refuted)
```
**Fix.** Chain all K STMs and sample `tGrid(2:K+1)`; count every bracket
inside (0, t_f); reserve `atFinal` for a sign change whose bracket ends at
t_f itself (a boundary case, not an interior refutation). Re-sweep once and
confirm the 0 count survives. The min-time catalogs share this instrument.

### P1.2 Conjugate test runs on unconverged solves (sol)

`ms_minfuel.m:67-73`, `ms_minenergy.m:87-96`: `info.conj` is computed
regardless of `info.converged`. Gate it. (The verdict sweep checks
`converged` downstream, so no shipped verdict is affected.)

### P1.3 Junction time grid is never persisted (sol + host)

`run_minfuel_race.m` arm output has no `tGrid`; `run_minfuel_grid.m:56`
stores `'tGrid', []` while the header advertises `.tGrid`;
`build_minfuel_catalog.m:76` reconstructs the last step as `rec.tf/K`. Host
verified the grid IS uniform (`linspace(0, tf, K+1)` at
`run_minfuel_race.m:76`, carried unchanged), so `rec.tf/K` is correct in
practice -- but the catalog's second reconstruction route ("reconverge
ms_minfuel from the Yj junctions") needs the grid, and nothing in the
deliverable says it is uniform. **Fix.** Store `tGrid` per entry (or state
`tGrid = linspace(0, tf_nd, K+1)` in `derive.reconstruction`), and use
`tGrid(end)-tGrid(end-1)` as the race already does at line 114.

### P1.4 Catalog builder trusts the rescue artifact and tolerates missing verdicts (sol)

`build_minfuel_catalog.m:29-43`: the (2,5)@1.4 rewalk is swapped in on
`Lr.R.conjPass` alone; `kBad` is not checked to exist; no residual /
endpoint / eps / fingerprint check binds the artifact to the entry.
`:59,80-89`: an entry with no verdict ships as `conj_pass = -1` and the v3
validator accepts it. Currently all 7 shipped verdicts are `1`, so latent.
**Fix.** Assert `~isempty(kBad)`, check `Lr.R.pDeepest`, `normR`, endpoint
miss; make `catalog_schema('validate')` reject `conj_pass(has_solution) == -1`
and enforce `int8` in {0,1,-1}.

### P1.5 Nothing enforces a minimum eps depth (sol)

`run_minfuel_grid.m:47-58` packages any arm with >= 1 converged rung. A walk
that retired at eps = 0.5 would ship labelled by its `p_floor`. All 7 shipped
records reached p_floor <= 0.005, so latent. **Fix.** A `pTarget` gate
(e.g. `A.p(end) <= cfg.pFloorMax`) before appending to `G`.

### P1.6 `coastFrac` is a junction count, not a time fraction (Astra; host had noted it)

`ms_minfuel.m:82`: `coastFrac = nnz(info.s < 1e-3) / K` counts junction
STARTS below threshold. With K = 12 its resolution is 8.3%, it aliases
short arcs, and it changes under mesh redistribution alone. FINDINGS 19's
"coast up to 58%" is 7/12 junctions. **Fix.** Integrate the indicator
over the flown trajectory (threshold-crossing times), and rename the
present statistic `coastJunctionFrac`. Same for `Hdrift`: it compares
junction Hamiltonians only; sample the propagated interior too.

---

## P2 -- robustness / hygiene

- `catalog_schema.m:186,192` -- `log(1/mf)` is matrix division; fine for
  the scalar calls made today, wrong for a grid. Use `1./mf`, and check
  `0 < mf <= 1`. (sol)
- `catalog_schema.m:73-80` -- unknown schema versions validate as v3;
  require version in {1,2,3}. (sol)
- `run_minfuel_grid.m:34-35` -- a `string` verdict falls into the non-char
  branch and can be accepted on stale gates; use
  `strcmpi(string(r.verdict),"PASS")`. (sol)
- `run_minfuel_race.m:129-139` -- `maxBisect` is documented per-gap but
  counted per-arm (up to 51 insertions under defaults). (sol)
- `run_minfuel_race.m` -- "save after every step" is not resumability; the
  queue/seed are never reloaded, so an interrupted arm restarts from the
  energy seed and may land in another basin. (sol; matches the
  matlab-campaign skill's resume rule)
- `cr3bp_minfuel_prop.m:65-75` -- no check that ode113 reached `dt`, no
  finite check, no m <= m_min / collision event. (sol)
- `cr3bp_minfuel_pmp.m:95-96,126-128` -- Huber accepts p > 1, where
  `s = min(max(pQ,0),p)` can exceed 1. Enforce `0 < p <= 1` for Huber. (sol)
- `build_minfuel_catalog.m:48,60,78` -- K taken from `G(1)` and imposed on
  every entry; mixed K would error. (sol)
- `cr3bp_minfuel_pmp.m:114-116` -- the 1e-300 primer regularisation makes
  `lam_v'*alpha = -|lam_v| + 1e-300/|lam_v|`; at lam_v == 0 the field can
  burn mass with zero thrust. Unreachable in practice; a comment or a
  `|lam_v| < tol` reject would close it. (sol; host rates this a nit)
- `run_minfuel_race.m:110-113,132`, `build_minfuel_catalog.m:73` -- growing
  arrays; preallocate. (sol)
- `run_minfuel_race.m:34-36,91` -- the race is not symmetric at its first
  rung. Both arms start at p = 1 from the energy seed; for `eps` that IS
  the energy field (s* = Q/2), but Huber kappa = 1 minimises s^2/2 - sQ at
  s* = Q, so its first step is a cold solve from an off-solution seed (the
  equivalent seed is lambda/2). Not a bug, but a confound to note beside
  the saltation issue when FINDINGS 18 is reworded. (Astra)
- Finite-eps entries are approximations of the bang-bang problem; the
  catalog says so but stores no error estimate. Store the neighbouring-eps
  m_f delta (the "~1e-5 by eps 0.005" figure) per entry. (Astra)
- REFERENCES sections missing on every `costate_common` math function in
  the bundle (house rule). Bertrand & Epenoy 2002 for the eps family;
  Caillau/Bonnard (hampath) for the conjugate test; the rocket equation for
  `derive`. (sol; host concurs)

---

## Refuted -- sol findings the code does NOT need to act on

**R1 -- WITHDRAWN, promoted to P0.3.** Host's first adjudication defended
the code's `Phi([1:6 14], 8:14)` block as "the focal-point condition for
this terminal manifold". Astra's Q2 answer (below) shows why that is wrong
for an INTERIOR test, and host concurs on re-derivation. Sol was right.

**R2. "The eps family is not the advertised tanh regularisation."** The
prompt planted tanh as a prior to be verified. The code implements the
polynomial Bertrand-Epenoy cost L = (1-p)s + p s^2 with
s* = clip((Q-(1-p))/(2p), 0, 1); host re-derived dH/ds = 0 and it matches
(lines 119-121). p = 1 reproduces the energy field exactly. Correct as is.

**R3. "Claim 19 not established -- acceptance-failed walks are published."**
Fact confirmed: `acceptOk` is true on 1 of 7 shipped records. But FINDINGS
19 already states that ss acceptance is not a valid gate at deep eps, and
host tested the deliverable's own recipe (single shot from `Yj(:,1,n)` over
`tf_nd` at `p_floor`):
```
cell   gamma  p_floor    miss_r[km]  miss_v[m/s]  acceptOk
(1,2)  1.10   0.00341        0.00        0.00       0
(2,5)  1.10   0.00179        0.09        0.00       0
(6,8)  1.10   0.00124        0.00        0.00       0
(1,2)  1.20   0.00100        0.01        0.00       0
(2,5)  1.20   0.00168        0.00        0.00       1
(6,8)  1.20   0.00128        0.02        0.00       0
(2,5)  1.40   0.00168        0.15        0.00       0
```
All 7 reconstruct to <= 0.15 km. The `acceptOk` flag measures the ss
gate's own `tolDz` on the costate vector, not endpoint accuracy. **Action:**
rename/document the field (sol's fair half) -- it invites exactly this
misreading -- but the record set is good.

**R4. "Claim 20 not established" -- UPHELD, not refuted.** With R1
withdrawn (P0.3), sol's blanket statement is correct: the 13/14 verdicts
were produced by the wrong interior block. They must be re-run.

---

## What each reviewer contributed

| | sol | Astra (raw API) | host |
|---|---|---|---|
| P0.1 Isp mislabel | missed | flagged the unchecked Isp/c requirement, no `.mat` | found + measured |
| P0.2 Huber saltation | found, sharp | found; explicit grad Q; eps exempt | corroborated via ode113 angle |
| P0.3 conjugate block | found | **settled it** (admissible-variation argument) | wrongly defended, reversed |
| P0.4 initial-coast zero | missed | **found** | measured: refuted record coasts at t_2 |
| P1.1 atFinal | found | found | measured 8.33% gap, 0 firings |
| P1.6 coastFrac | missed | found | had noted, dropped |
| Huber seed asymmetry | missed | found (lambda/2) | confirmed in race sched |
| R3 acceptOk | asserted unsafe | (out of scope) | refuted by reflying all 7 |
| hygiene | 10 items | ~8, overlapping | concur |

Sol: strong bundle reviewer, over-claims severity on what it cannot
measure; every claim against *code paths* held, its one claim against
shipped *results* dissolved. Astra: the best reviewer of the three on
theory -- it broke the block dispute, found the coast artifact nobody
else saw, and every one of its CORRECTNESS items survived measurement.
Host: the artifacts (P0.1, the reflight, the coast table, the 0-firing
counts) and one wrong theory call. Division of labour that worked:
external models on the source, host on the `.mat` files, Astra on any
theory dispute.

---

## GPT-6 Astra -- the invocation problem, diagnosed

Requested reviewer. Via **crush** it wedged on both real bundles -- 121 KB
(22 min, 0 bytes, 0% CPU) and 53 KB (15 min, 0 bytes) -- with no `-c`,
while answering a one-line smoke prompt instantly and a 15 KB two-file
probe in ~1 min. Sol digested the identical 121 KB bundle through crush in
6 min, so this is Astra-specific and NOT the documented
long-prompt-plus-`-c` failure.

**Raw OpenRouter API, same 53 KB bundle that wedged crush:**
```
HTTP 200   total 146 s   prompt 17,432 tok   completion 7,454 tok (3,624 reasoning)
finish_reason: stop   cost $0.59
```
16.7 KB of review, `minfuel_code_astra_core_rawapi.md`. **The wall is
crush, not Astra.** Crush's OpenRouter client evidently drops or never
surfaces Astra's response above some prompt size in the 15-53 KB range
(the request takes 2+ min of server-side reasoning; a client timeout that
is swallowed silently fits the symptoms). Rule recorded in the
opencode-headless skill: **for Astra, bypass crush and call the API
directly** (`references/astra_raw.sh`); crush remains fine for the
GPT-5.6 tiers and Gemini.

Quality: best of the three reviewers on theory (see the table above).
At $0.59 for a 53 KB bundle it is not expensive in absolute terms; the
"final gates only" positioning is about not spending it on iteration, not
about cost per call.

**Two-file crush probe (Q1 saltation, Q2 block): SUCCESS, ~1 min, 2.7 KB,
three findings, all correct.**

- Q1 (saltation): "Yes, generally" -- gives the saltation matrix with the
  explicit switch-surface gradient
  `n = grad Q = Tmax * [0; 0; -rho/m^2; 0; lam_v/(m rho); 1/c]`,
  one-sided fields on either branch, the grazing caveat, and confirms the
  eps family is exempt (F+ = F-, S = I). Matches sol; adds the concrete n.
- Q2 (block): `Phi(1:7, 8:14)` for the interior test; the mixed block is
  the terminal derivative, justified only at t_f; gives the
  admissible-variation argument that settled the dispute (P0.3).
- Q2 (endpoint): `atFinal` exclusion "No" -- every sample is at
  `tGrid(2:K)` < t_f; fix `atFinal = false; nIn = numel(flips) +
  nnz(dets == 0)`. Matches sol and host.

**Superseded lesson.** The earlier reading -- "Astra is bounded by prompt
size, use it only as a two-file tie-breaker" -- was wrong; the bound was
crush's. Through the raw API it is a full bundle reviewer. The tie-breaker
pattern (two files, pointed questions) remains a good, cheap way to use it
when a dispute is already localised.
