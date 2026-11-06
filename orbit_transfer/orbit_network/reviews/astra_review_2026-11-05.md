# GPT-6 Astra review of the orbit_network spec + plan, 2026-11-05

Raw OpenRouter API, reasoning effort high, 100.8 KB prompt (spec + plan + built partial_f_test + study sections 3-4), 550 s, 30,750 prompt / 17,994 completion tokens, cost $1.28. Host verification and dispositions: see reviews/astra_review_2026-11-05_dispositions.md.

## Overall verdict

**Do not ship the note as currently prescribed.** The feasible-set inequality and the incomplete-beta formulas are sound, but several proposed interpretations are not:

- A catalog edge is not automatically a globally optimal transfer or a physical feasibility boundary.
- Independently minimized leg times are **not** a feasible relay travel time under a time cap.
- Time-reversed, mass-increasing arcs are **not approximate certified reverse transfers** for the original spacecraft problem.
- The proposed “wall” narrative contradicts the feasible-set lemma.
- Several graph calculations can silently produce incorrectly attributed weights or incorrect summary statements.

Below, findings are ranked by the cost of shipping them: scientific claims first, numerical/reporting defects next, figures last.

---

## 1. The lemma needs full-state compatibility—and it contradicts the proposed “topology wall” explanation

### Quoted text

Task 9:

> “every concatenation of an admissible \(A\to C\) transfer, a coast on \(C\), and an admissible \(C\to B\) transfer is itself admissible for \(A\to B\).”

And:

> “the question is relay CLOSURE below the 67 mN direct wall, which the lemma does not settle because the wall is a basin/topology wall, not a cost.”

### What is wrong

The inequality is correct **when the relay trajectories are a subset of the same unrestricted \(A\to B\) feasible set**. But orbit labels do not establish that inclusion. Compatibility must cover:

- position **and velocity** at the joins;
- spacecraft mass;
- initial and terminal conditions;
- thrust, mass, and other path constraints;
- any absolute-time or phase restrictions;
- absence of a relay-only operation such as refueling.

Additivity is **not necessary** for the inequality. Set inclusion gives the inequality for any common cost functional. Additivity only explains how to evaluate the cost of concatenated pieces.

The “wall” sentence is a more serious contradiction. If an admissible relay trajectory exists below a claimed physical \(A\to B\) feasibility wall, it is also an admissible unrestricted \(A\to B\) trajectory. Thus the wall cannot be a global feasibility barrier under the same assumptions. It could instead be:

- termination of a particular continued extremal branch;
- a solver/basin limitation;
- a restriction on phases, topology class, or allowed trajectory structure.

Calling it a “basin/topology wall” does not evade feasible-set inclusion.

### Write instead

> **Lemma—unrestricted versus relay-constrained optimum.** Let \(\mathcal F_{AB}\) be the feasible trajectories of a specified \(A\to B\) optimal-control problem, including spacecraft mass and all boundary and path constraints. Let \(\mathcal R_{ACB}\subseteq\mathcal F_{AB}\) consist of trajectories passing through relay orbit \(C\), with compatible full states at every join and no additional relay-only operations. For the same cost functional \(J\),
> \[
> \inf_{x\in\mathcal F_{AB}}J[x]\leq
> \inf_{x\in\mathcal R_{ACB}}J[x].
> \]
> This statement concerns global infima, not the costs returned by particular direct or relay solution procedures.

Replace the wall sentence with:

> “A relay witness below the termination thrust of the computed direct branch would establish another feasible unrestricted transfer. It would diagnose a limitation of that branch or search, not passage through a proven global feasibility barrier.”

Also replace:

> “their headline … is an artifact of [their maneuver model]”

with:

> “A multileg advantage is relative to their restricted direct-maneuver class. It is not a comparison against the global optimum over a feasible set already containing those multileg trajectories.”

The lemma does not justify dismissing their result, nor the sentence that an OCP network is “never about relay savings”: decomposition can improve a **computed candidate**, even though it cannot beat the corresponding global unrestricted optimum.

### Important mass-model qualification: the catalog legs can often be composed exactly

Do **not** replace the current claim with the equally incorrect blanket claim that independently initialized 150 kg legs cannot be used.

For
\[
\dot m=-F/v_e,\qquad 0\leq F\leq F_{\max},
\]
suppose leg 1 arrives with mass \(m_1\), while the stored leg 2 starts at \(m_0=150\) kg. Set
\[
q=m_1/m_0,\qquad
\widehat m_2(t)=q\,m_2(t),\qquad
\widehat F_2(t)=q\,F_2(t).
\]
Then
\[
\widehat F_2/\widehat m_2=F_2/m_2,
\]
the mass equation still holds, and the thrust bound holds because \(q\leq1\). The stored leg-2 position–velocity trajectory and duration are preserved exactly.

Thus **with throttleability, free coast, and no violated dry-mass or other constraint**, this supplies a no-new-solves feasibility construction. It does not preserve leg-2 optimality or its certificate for the changed initial-mass problem.

If the stored terminal masses are \(m_{f1},m_{f2}\), the assembled terminal mass is
\[
m_{f,\mathrm{route}}=\frac{m_{f1}m_{f2}}{m_0}.
\]
Check it against any dry-mass requirement. If mass histories are unavailable, a conservative bound can be obtained from durations and \(F_{\max}/v_e\).

For minimum propellant, consumption is additive over the **actual assembled pieces**, not the two nominal 150 kg catalog consumptions:
\[
\Delta m_{\mathrm{route}}
=(m_0-m_{f1})+\frac{m_{f1}}{m_0}(m_0-m_{f2}).
\]

### Is the coast admissible?

Yes, in the stated ideal autonomous CR3BP: at an exact periodic-orbit position–velocity state, zero thrust follows that orbit and leaves mass constant. Linear instability does not destroy nominal admissibility; it affects sensitivity and operational robustness.

---

## 2. The relay time-cap contour currently declares access using a non-achievable lower bound

### Quoted code

Task 4:

```matlab
tt = E(e1).tf_min_days + E(e2).tf_min_days;
```

Task 6:

```matlab
P.multileg(rr, cc) = nnz(arrayfun(@(b) b.t_total_days(rr) <= caps(cc), best));
```

### What is wrong

The two minima can use different relay phases. Their sum omits the required hand-over coast. Consequently, the green contour can claim a route meets a cap when no assembled catalog route does.

Let the actual relay arrival phase of the chosen first leg be \(\alpha\), and the departure phase of the chosen second leg be \(\beta\), expressed as time fractions with a common origin and direction. Then
\[
w=T_C\,\operatorname{mod}(\beta-\alpha,1).
\]
The assembled duration is
\[
t_1+w+t_2.
\]

The sum \(L=\min t_1+\min t_2\) is a lower bound for the collection of these fixed-duration catalog-template assemblies. It is **not** a lower bound on the globally optimized relay OCP: unsearched phases and a reoptimized second leg at lower initial mass may be faster.

### Is the one-period bound correct and tight?

**Yes, conditionally.**

For arbitrary phases on the same periodic orbit,
\[
0\leq w<T_C.
\]
Thus \(L+T_C\) is a valid conservative upper bound on the duration of the assembled minimizing legs, subject to the mass and admissibility qualifications above. For arbitrary continuous phases, \(T_C\) is the tight supremum.

For **six distinct, equally spaced time phases on the same grid and with the same phase origin**:
\[
w\in\{0,T_C/6,\ldots,5T_C/6\}.
\]
Then \(5T_C/6\) is the tight uniform bound.

But “6-point grids” alone does not establish this. Inspect `sA_frac`, `sD_frac`, and their conventions:

- Are both grids equally spaced?
- Do they have the same phase origin?
- Are 0 and 1 both included, duplicating an endpoint?
- Is arrival phase the actual terminal phase, or the target’s phase at the transfer’s initial epoch?

In the last case, convert to the actual terminal phase before taking the difference. A \(T_C/6\) bound is **not** generally available: the next grid departure need not be closed or attain the leg minimum.

### Write instead

> “For each pair of catalog legs, compute the forward coast from the actual relay arrival state to the selected relay departure state. The independent-minima sum is an optimistic template bound; the phase-compatible assembled time is the quantity compared with a mission time cap.”

Change `relay_closure` to receive the phase-resolved data—or retain that metadata in `E`. With no new solves, compute
\[
\min_{\substack{\text{closed cell }(a,b)\text{ on leg 1}\\
                 \text{closed cell }(c,d)\text{ on leg 2}}}
\left[
t_1(a,b)+T_C\operatorname{mod}(\phi_D(c)-\phi_A(b),1)+t_2(c,d)
\right].
\]

If that computation is deferred, plot two explicitly labeled quantities:

- **optimistic template bound**: sum of independent minima;
- **conservative assembled-time bound**: sum plus the justified wait bound.

Only the second can establish access under a cap.

Also, Task 6 selects the relay minimizing `tt` and then copies its `tw`. That does **not** minimize the wait-inclusive time or upper bound. Optimize the criterion actually advertised.

---

## 3. The mirror option is not an approximate solution of the original reverse OCP

### Quoted text

Spec §1:

> “Mirrored (reverse) edges … labelled approximate (mass flow breaks the exactness: the reversed trajectory runs \(m_f\to m_0\)).”

Task 1:

> “With time reversed, the image of a trajectory is a trajectory…”

### What is correct

Let
\[
S=\operatorname{diag}(1,-1,1).
\]
For a forward arc, define
\[
\widetilde r(s)=S r(t_f-s),\qquad
\widetilde v(s)=-S v(t_f-s).
\]
This is the stated position–velocity reverser.

The required thrust acceleration transforms as
\[
\widetilde a_T(s)=S a_T(t_f-s).
\]
Therefore the thrust **direction** transforms as
\[
(u_x,u_y,u_z)\longmapsto(u_x,-u_y,u_z),
\]
not with the velocity signs \((-u_x,u_y,-u_z)\).

### What is wrong

With
\[
\widetilde m(s)=m(t_f-s),
\]
the transformed mass equation is
\[
\dot{\widetilde m}=+F(t_f-s)/v_e.
\]
That is a mass-gaining spacecraft. It violates the original physical dynamics, and usually the reverse problem’s initial condition of 150 kg.

This is not merely loss of an optimality certificate. **Feasibility in the stated reverse problem is lost.** No approximation error or near-optimality bound is provided.

For example, the reversed kinematics of an endpoint reached under full thrust at mass \(m_f<m_0\) demand acceleration \(F_{\max}/m_f\). Starting the reverse spacecraft at \(m_0\) would require thrust greater than \(F_{\max}\).

### Is it an extremal of a well-posed problem?

Yes: time reversal gives an exact correspondence to an artificial **mass-gain problem**, with appropriately swapped boundary conditions and transformed costates/transversality. Minimum-time local/global optimality transfers under that exact correspondence.

But the endpoint conditions matter: an original free terminal mass becomes a free initial mass in the transformed problem. This is not the usual reverse transfer initialized at 150 kg.

### Write instead

> “The CR3BP reverser maps the spatial trajectory exactly when thrust acceleration is reflected by \(\operatorname{diag}(1,-1,1)\). Reversing the mass history produces an exact trajectory of a mass-gain model, not of the physical mass-depleting reverse problem. Reverse edges shown here are hypothetical geometry-based edges; their physical feasibility, optimality, and certification are unestablished.”

Do not copy `n_certified` onto these edges as though certificates transfer. If retained, label the option **HYPOTHETICAL REVERSE GRAPH**, not “approximate certified reverse transfers.”

The node-symmetry gate only checks whether endpoint orbits map to themselves. It cannot repair the powered-trajectory or mass issue.

---

## 4. “Certified edges” are not what the primary calculations actually select

### Quoted code

Task 2:

```matlab
h = s.has_solution(:, :, rr);
e.closed_frac(rr) = nnz(h)/nCell;
e.n_certified(rr) = nnz(h & s.conj_pass(:, :, rr));
```

The minima, plane, graph, and relay closure all use `has_solution`, not the certificate mask.

### What is wrong

The synthetic test deliberately admits three closed cells but only two conjugate-pass cells. Thus the planned implementation explicitly permits an uncertified cell to determine an edge’s advertised minimum time.

Moreover, `conj_pass` alone may not encode all the hypothesis gates required by the certification procedure.

There is also a separate scope issue: a conjugate-point test with appropriate hypotheses generally certifies a local optimality statement. It does not, by itself, prove the global infimum over all controls, all branches, or continuous phases.

### Write instead

Choose and enforce one of these definitions:

1. **Certified-witness graph:** all counts and minima use
   ```matlab
   accepted = has_solution & conj_pass & all_required_hypothesis_gates;
   ```
   with finite positive times and the relevant feasibility checks.

2. **Feasible-solution graph:** counts use verified feasible solutions, and certified counts are separately reported. The graph and its minima must not be called certified minima unless the selecting witness is certified.

Use:

> “minimum duration among accepted catalog cells”

rather than an unqualified:

> “minimum time between the orbits.”

Add a test where the fastest `has_solution` cell fails certification. The certified graph must not select it.

---

## 5. The statistical formulas are correct; the scientific verdict needs narrower wording

### Quoted text

Spec §6:

> “Verdicts, printed as SUPPORTED / NOT SUPPORTED…”

Built code:

```matlab
sup = o.pt1(3) < alpha && signOK;
```

### What is correct

Both formulas are correct under their usual distributional assumptions:

\[
P(F_{d_1,d_2}\geq F)
=I_{d_2/(d_2+d_1F)}(d_2/2,d_1/2),
\]

\[
P(|T_\nu|\geq|t|)
=I_{\nu/(\nu+t^2)}(\nu/2,1/2).
\]

Using the departure coefficient’s own test, rather than attributing the joint F result to departure stability, is the right correction.

The supplied results do **not** detect the proposed departure effect:

- `tf1`: departure \(t=0.58\), with the wrong sign for the stated prediction;
- `surv`: departure \(t=1.20\), \(p=0.24\);
- the strong joint time result is associated principally with `sigma_arr`.

### What is wrong or unsupported

“NOT SUPPORTED” is defensible only as shorthand for **not detected by this particular analysis**, not as a rejection of the Braik–Ross mechanism.

There are several limitations:

- There are 56 sheets, but only **14 distinct departure orbits**, six stable and eight unstable. Repeated sheets sharing departure and arrival nodes are not automatically independent observations.
- Stability and family are partly confounded. The model does not include family effects.
- `sigma_arr` has support only in the halo-to-halo catalog. It varies within that catalog, so it is not literally identical to a catalog indicator, but its interpretation is strongly catalog-restricted.
- A linear coefficient on `sigma_dep` tests a particular monotonic association, not every stable-versus-unstable effect.
- A signed linear `dCJ` term does not control all energetic or geometric difficulty.
- `surv` is bounded, has varying denominators, and measures solver/catalog retention—not a known physical survival probability.
- `tf1` is an order statistic over different numbers of successful cells.
- These are selected deterministic catalogs, not an identified random sample from an orbit population.

The nominal OLS t/F reference distributions therefore need qualification.

### Does this design have power?

It is not correct to say it has **no** power. It can detect a sufficiently large association. But 56 alone does not establish useful power for the physical effect of interest.

Even an idealized independent two-group comparison with 24 versus 32 observations needs roughly a **0.76-standard-deviation** difference for 80% power at a two-sided 5% level. Confounding and shared-node dependence change that calculation, potentially substantially. Power for the conditional `sigma_dep` slope also depends on its residual variation after the other predictors.

Report coefficient confidence intervals and avoid interpreting nonsignificance as equivalence.

### Write instead

> “Departure instability was not associated with shorter sampled minimum durations or greater catalog retention at the nominal significance threshold in this exploratory sheet-level model. The time model’s significant joint stability term is associated primarily with arrival instability on the halo-to-halo sheets. These results neither establish absence of a departure-stability effect nor test the mechanism in the substantially lower-thrust regime. Family overlap is limited, sheets share orbit nodes, and the reported OLS p-values are nominal model-based summaries.”

Keep the raw means, but label them descriptive. Do not present the 1 N/0.5 N results as resolving behavior near 70 mN.

---

## 6. The proposed “thrust floor” is a lowest successful rung, not a measured feasibility threshold

### Quoted text/code

Spec §7:

> “the thrust floor … reported as a bracket”

Task 2:

```matlab
e.F_floor_str = sprintf('%g', s.rungs_N(deep));
```

### What is wrong

Failure to close a cell or sheet at a lower rung is not a proof of infeasibility there. This is already demonstrated by the nonmonotone coverage near 15 N.

If a feasible witness exists at \(F_{\mathrm{found}}\), the physical threshold—if defined—is bounded **above**:
\[
F_{\mathrm{physical}}^\star\leq F_{\mathrm{found}}.
\]
That reasoning applies at **every** successful rung, not only at 0.5 N.

An unsuccessful next rung does not provide the lower endpoint of a physical bracket.

### Write instead

> “Lowest successful tested thrust: 1 N; lower tested rungs unresolved.”

At the bottom:

> “Successful at the lowest tested thrust, 0.5 N; no lower-thrust conclusion.”

If displaying bounds on a physical threshold, use `<= F_found` for all witnesses, subject to their feasibility qualifications.

Likewise, “edge survival” should be described as **catalog closure/retention down the continuation ladder**, not demonstrated loss of physical accessibility.

---

## 7. Thrust is a legitimate axis here—but the current plane mixes capability with solver coverage

### Quoted text

Spec §1:

> “Their \((dV\text{ cap},time\text{ cap})\) budget plane becomes our \((thrust,time\text{ cap})\) accessibility plane.”

Task 3:

> “the gap between the two counts is phase robustness.”

### What is wrong

For the same 150 kg initial spacecraft and Isp, \(F_{\max}\) is a sensible common hardware-capability axis. Initial acceleration \(F_{\max}/m_0\) is just a constant rescaling. Changing the label to acceleration would not fix the substantive issues.

But thrust is not a consumable budget like total \(\Delta v\). Acceleration also varies with mass. For the stated propulsion model,
\[
\Delta v_{\mathrm{propulsive}}=\int F/m\,dt
=v_e\log(m_0/m_f),
\]
which generally cannot be recovered from duration alone without additional thrust/mass information.

Further, with a **bounded, throttleable** thrust model, a trajectory found at 12 N is feasible under a 15 N cap. Therefore true accessibility cannot decrease as the cap increases.

The raw rung-specific plane is a **catalog coverage plane**, not the full set of known feasible witnesses under a thrust cap.

### Write instead

> “We use a thrust-capability/time-cap plane, not a propellant-budget plane. Raw rung-specific counts describe catalog coverage. A separate known-feasible envelope includes every accepted trajectory computed at any thrust bound no greater than the displayed cap.”

That envelope can be constructed without new solves. It must be labeled a feasibility-witness envelope: a trajectory certified as locally minimum-time at 12 N is not thereby certified optimal at 15 N.

Also replace “the gap … is phase robustness” with:

> “Cell counts describe sampled phase-grid coverage. Report the covered fraction of the 36 cells, optionally conditional on the pair having at least one witness.”

`cells - pairs` is not a meaningful robustness measure, and six sampled phases do not establish robustness over continuous phase.

---

## 8. The three graph diagnostics are not three physical accessibility regimes

### Quoted text

Spec §7:

> “weak connectivity of the 18-node graph…”

Task 9:

> “The three regimes in our axes…”

### What is wrong

Weak connectivity ignores direction. In this graph:

- DROs, DPOs, and L1 halos have no incoming edges;
- tulips have no outgoing edges;
- only L2 halos can be internal vertices of directed paths;
- there are no directed cycles.

At 1 N, with all 56 edges present, the graph can be weakly connected while all strongly connected components remain singletons.

There are only eight additional directed endpoint pairs supplied by two-leg paths: L1 halo to tulip. Thus the full original topology contains at most **64 reachable ordered distinct-node pairs**, not all \(18\cdot17=306\).

Also, retaining each edge whose time is below a cap does not establish that a multiedge path meets that **total** cap. Weak connectivity establishes even less.

Finally, these diagnostics need not occur in the proposed order. A first relay route can appear before the graph spans all 18 nodes.

### Write instead

> “We display three distinct catalog diagnostics: direct-edge availability, weak structural connectivity after ignoring direction, and phase-compatible two-leg witness availability under a total time cap. They are not mutually exclusive or necessarily sequential regimes, and weak connectivity does not imply directed mission accessibility.”

Explicitly distinguish:

- **unmeasured edge**;
- **measured but unresolved edge**;
- **available witness**.

The missing 250 directed pair catalogs cannot be treated as physical prohibitions.

---

## 9. The centrality definitions are valid graph statistics, but their dynamical interpretation is severely constrained

### Quoted text

Spec §7:

> “in-strength and out-strength with weights \(1/t_f\)… and betweenness…”

Task 9:

> “Which relay orbit the present data points to…”

### What is misleading

Weighted degree with \(1/t_f\) is a legitimate statistic, but it is not throughput, transport capacity, or an invariant “hub” property.

The graph was designed with very different opportunities:

- each tulip has 12 incoming catalog edges;
- each L2 halo has two incoming and four outgoing edges;
- only L2 halos are eligible to have nonzero betweenness.

Thus much of the apparent hub/relay hierarchy is imposed by catalog selection.

Moreover, betweenness weighted by sums of independently minimized edge times uses phase-relaxed costs. It is not betweenness of physically timed assembled missions.

With the usual unnormalized directed definition, all betweenness is distributed among the L2 halos over the eight L1-to-tulip endpoint pairs. That is a useful structural fact, not evidence that L2 halos are globally privileged relays in CR3BP transport.

### Write instead

> “Centralities describe this selected, phase-relaxed catalog graph. Between-family rankings largely reflect the available catalog topology. L2-only betweenness is structural. Within-L2 rankings may suggest candidates for further computation, but do not predict the best DRO relay or performance near 70 mN.”

The existing L1-to-L2 legs do not tell you the missing DRO-to-L2 leg costs. Current data can support candidate prioritization, not identification of the best future DRO relay.

---

## 10. `digraph` edge reordering can silently assign the wrong costs and figure attributes

### Quoted code

Task 5:

```matlab
G = digraph(edgeTab, nodeTab);
cost = tfAll(:, rC);
bc = centrality(G, 'betweenness', 'Cost', cost);
```

Task 7:

```matlab
Gr = digraph(src(keep), dst(keep), ones(numel(keep), 1), nN);
h.LineWidth = 0.5 + 4*cfAll(keep, rr);
h.EdgeCData = tfMin(keep, rr);
```

### What is wrong

MATLAB graph constructors sort edges. Input `E` order is not guaranteed to equal `G.Edges` order.

The constructor carries table attributes with their edges, but the separately supplied `cost`, widths, and colors remain in catalog order. This can change centrality rankings and attach another edge’s data to a plotted edge.

### Write instead

For centrality:

```matlab
cost = G.Edges.tf_min_days(:, rC);
```

For filtered graphs, retain a `sheet_index` edge attribute and obtain widths/colors from the filtered graph’s own edge rows, or explicitly map its edges back to `E`.

Add a synthetic test with deliberately scrambled input edges and a known shortest route. The current tests will not catch this.

---

## 11. Several printed network results are wrong or misleading in boundary cases

### A. “Mutually reachable pairs” is hard-coded to the number of edges

Quoted code:

```matlab
fprintf(..., numedges(Gm), nS, nS, nS, nTwoLeg);
```

The fourth value, `nS`, is printed as mutually reachable pairs.

For the fully connected bidirected hypothetical graph, every two nodes are mutually reachable: **153 unordered pairs**, or **306 ordered pairs**, not 56.

**Replace with:** an SCC-based count, and state whether pairs are ordered. If the intended quantity is direct reciprocal links, print:

> “56 unordered directly reciprocal links.”

Likewise, `nTwoLeg` counts endpoint pairs, not distinct relay routes.

### B. The connectivity summary falsely reports success if no cap connects the graph

Quoted code:

```matlab
[~, cFull] = max(P.largest(rC, :) == nN);
```

If all entries are false, this returns index 1.

**Replace with:**

```matlab
cFull = find(P.largest(rC, :) == nN, 1, 'first');
```

and explicitly handle empty:

> “No sampled cap yields a spanning weak component.”

When nonempty, say **first sampled cap**, not the exact threshold.

### C. The cap range contradicts the spec and truncates phase coverage

Quoted code:

```matlab
tfAll = vertcat(E.tf_min_days);
caps = linspace(min(tfAll(:), [], 'omitnan'), ...
                max(tfAll(:), [], 'omitnan'), nCaps);
```

The spec requests the smallest to largest **closed-cell time**. The code uses the smallest to largest **sheet minimum**. Slow closed cells are consequently omitted from the plane’s range.

The range may also exclude relevant relay times entirely.

**Replace with:** a range derived from all accepted cell times, extended as needed for the advertised relay contour. Compute exact event thresholds separately from the 25-point plotting grid if quoting threshold values.

---

## 12. G4 and G5 do not validate what the note will imply they validate

### A. G4 repeats the construction and passes vacuously for no routes

Quoted code:

```matlab
g4 = true;
for kk = 1:numel(R)
    legs = E(R(kk).e1).closed_frac > 0 & E(R(kk).e2).closed_frac > 0;
```

This checks the same Boolean expression used to create `R.closed`. It is a useful indexing consistency check, but not a test of physical concatenability, phase hand-over, or mass feasibility. It also passes if route enumeration returns nothing.

**Write instead:**

> “G4 checks route enumeration and consistency with edge availability; separate checks establish phase-compatible assembly and applicable mass constraints.”

Require:

- 32 distinct route triples;
- eight endpoint pairs;
- correct edge endpoints;
- expected open/closed masks;
- correct phase waits and assembled times.

### B. G5 is effectively automatic at the stored symmetry-plane initial states

Quoted code:

```matlab
Rx = mirror_state(x0);
...
[~, k] = min(d);
```

For \(y=v_x=v_z=0\), \(Rx_0=x_0\), so the residual at \(t=0\) is zero regardless of the supplied period. It can pass for a nonperiodic trajectory initialized on the fixed set.

This identity is mathematically useful **if periodicity is independently established**, but it does not validate periodic closure, the claimed search near \(T\), or reverse powered trajectories.

The perturbed-state test does not fix this. A generic perturbation is not guaranteed to have a particular positive minimum residual.

**Write instead:**

> “G5 checks endpoint-orbit invariance under the reverser, conditional on independently verified periodic closure.”

Add independent period-closure checks and, if numerical reversibility is being tested, use off-symmetry-plane states. The current test does not implement the review-focus promise to verify behavior near \(T\).

There is also no reason for failure of an optional mirror-only check to invalidate the non-mirrored headline analysis.

G3’s monotonicity check, by contrast, is a reasonable code invariant—just do not present it as scientific validation.

---

## 13. `partial_f_test` has correctness gaps, despite having the right formulas

### Quoted code

```matlab
df1 = p1 - p0;
df2 = n - p1;
...
if RSS1 <= 0
    p = 0;
```

### What is wrong

- Column counts are valid degrees of freedom only when the designs have the required ranks.
- A rank-deficient `X1` can pass the column-space nesting test and still produce invalid F and coefficient tests.
- If both models fit exactly, `RSS0 = RSS1 = 0`, assigning `p = 0` falsely declares significant improvement.
- Tiny roundoff-induced negative `RSS0-RSS1` can put the incomplete-beta argument outside its valid domain.
- `inv(X1.'*X1)` amplifies numerical problems for poorly conditioned designs.

These do not establish that the supplied current results are wrong; they are defects in the general helper.

### Write instead

Require finite inputs, full-rank designs, positive residual degrees of freedom, and an explicit policy for zero residual variance. Clamp only tiny negative RSS differences attributable to roundoff; reject materially negative differences. Obtain coefficient covariance through a rank-aware factorization.

The built column-**space** nesting test is correct and more general than the spec’s literal “subset of columns” wording.

---

## 14. Two property-table claims are stronger than their checks

### A. “Sigma is zero when stable” does not match the stated implementation

Quoted definition:

```matlab
sigma = max(0, log(lam_max)/period_nd)
stable = lam_max < 1 + 1e-3
```

A numerically stable node with `lam_max = 1 + 2e-4` gets positive sigma.

**Write instead:** either explicitly set sigma to zero under the adopted numerical stability classification, or retain the raw estimate and stop claiming exact zeros on stable nodes. Report the tolerance and any sensitivity to it.

Also, “no resolved exponential instability in the computed multipliers” is more precise than implying nonlinear stability or invariant transport barriers. The planar stable-torus interpretation does not transfer automatically to spatial CR3BP.

### B. Period agreement does not make a cache trustworthy

Quoted text:

> “so a stale cache cannot pass silently.”

A cache with the same keys and periods but stale states, Jacobi constants, multipliers, constants, or algorithm versions can pass.

**Write instead:**

> “Period agreement detects period mismatches; cache validity additionally requires matching orbit parameters, physical constants, and numerical/provenance metadata.”

Include that provenance in the cache key or require rebuilding for publication runs.

---

## 15. Lower-cost plan/spec contradictions and figure defects

### A. The mirror figure option is not implemented

Spec:

> “mirrored edges dashed if on”

Task 7 always constructs `Gr` from original `src` and `dst`; `Gm` is unused.

**Replace with either:** a separate, clearly hypothetical reverse-graph figure, or remove the promised visual option. Keeping headline figures non-mirrored is preferable.

### B. The connectivity contour is drawn at the maximum integer value

```matlab
contour(..., P.largest, [nN nN], ...)
```

A contour at the maximum plateau can disappear or fail to mark the intended transition.

**Replace with:** a contour at 0.5 of the Boolean mask `P.largest == nN`, handling all-true/all-false masks explicitly.

### C. The summary file does not contain all numbers the note is instructed to cite

The centrality table and per-pair best-relay rows are printed, but not appended to `lines`. Task 9 nevertheless makes `network_summary.txt` the network-number source.

**Replace with:** persist those tables and the phase/mass assumptions alongside the summary, or explicitly capture and archive the complete network-script output.

### D. Figure-existence tests are insufficient for the highest-risk errors

Add hand-checked tests for:

- reordered graph edges;
- mismatched relay phases;
- a relay whose smallest independent sum is not its smallest wait-inclusive time;
- a fastest cell failing certification;
- no spanning component at any cap;
- mutually reachable versus directly reciprocal pairs.

Those tests protect published claims. Checking for `PASS`, `112 edges`, and PNG existence does not.

---

## Recommended note framing

A defensible central claim is:

> “We construct a directed graph of available catalog witnesses, with sampled phase coverage and duration attributes across thrust limits. The graph exposes the selected catalogs’ coverage and their sole measured relay structure, L1 halo–L2 halo–tulip. Phase-compatible, mass-consistent catalog assemblies can provide feasible multileg witnesses without new optimal-control solves. These data do not establish global thrust thresholds, reverse-transfer solutions, relay savings over unrestricted global optima, or absence of a general stability–transport relationship.”

That preserves the useful network study while removing the claims most likely to require a published correction.