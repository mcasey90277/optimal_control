Here is the expert review of the **GTO->ELFO retarget layer**, focused strictly on the post-triage lessons (L1–L6) and readiness for the upcoming $T_{max}$ thrust ladder. 

### 1. Findings

**[METHODOLOGY]** `gen_elfo_minfuel.m:101` (and `elfo_export_data.m:82`) -- **L1: Switch-Count Mesh Under-resolution**. The mesh size $N=4001$ is statically inherited from the backbone seed, but sweeping $t_f$ or $T_{max}$ alters the transfer duration. A static boolean edge-crossing count `sum(abs(diff(ss>0.5)))` is pseudo-exact and becomes corrupted by severe aliasing as nodes-per-rev dilutes. 
*Fix*: Parameterize $N$ dynamically based on $t_f$ (e.g., $N \approx t_f \times \text{nodes\_per\_ND}$) and do not hardcode integer switch counts (`sw%d`) directly into canonical data product filenames, which breaks parsers when resolution noise adds ghost switches.

**[ROBUSTNESS]** `gen_elfo_minfuel.m:123` (and `gen_elfo_energy_gravhom.m:202`) -- **L2: Under-Iteration Collapse**. Continuation acceptance demands `rT.success && rT.maxDefect < 1e-6`. CasADi registers IPOPT’s `Solved_To_Acceptable_Level` as `.success=true`, accepting states with sloppy 1e-5 dual inaccuracy. Cascading "acceptable" KKT costates into the warm-start of the next $\epsilon$-homotopy step accumulates dual error until the continuation tail inexplicably collapses. 
*Fix*: Explicitly demand full dual convergence for continuation acceptance: `ok = (strcmp(rT.ipoptStatus, 'Solve_Succeeded') && rT.maxDefect < 1e-6)`.

**[CORRECTNESS]** `casadi_energy_freetf.m:158` -- **L3: Partial Manual Scaling Trap**. The trap-rule collocation defects `D` are $O(\Delta\sigma) \approx O(10^{-4})$. IPOPT's `gradient-based` scaling establishes the Jacobian scale of these rows as 1.0 (from $\Delta X$), leaving the macroscopic $O(10^{-4})$ offset exposed directly to the `constr_viol_tol = 1e-7`. This permits a massive $0.1\%$ relative integration drift *per step* because of unmanaged $O(\Delta \sigma)$ vs $O(1)$ scaling disparities.
*Fix*: Manually multiply the defect equations by $N$ to balance their scale against boundary and unit constraints: `opti.subject_to( N * D(:) == 0 )`.

**[EFFICIENCY]** `gen_elfo_energy_tfsweep.m:138` -- **L4: Warm-Start Phase Aliasing**. The $t_f$ sweep propagates oscillatory commands `Uk` (which contain rev-locked apogee/perigee burns) directly by node index into a transfer $0.5$ ND longer ($\approx 1.5$ revolutions). Pure $\sigma$-indexing phase-aliases the controls wildly against the true anomaly of the prolonged orbit, stalling the $t_f$ sweep prematurely.
*Fix*: 1D-interpolate the warm-starts for $X$ and $U$ mapped to physical time (`t_f * sigma`) rather than fractional $\sigma$. 

**[ROBUSTNESS]** `casadi_energy_freetf.m:256` -- **L5: Bound-Saturation Blind Spot**. Explicit bounds enforce $v \in [-12, 12]$ (barely above Earth perigee speed of 9.9 ND) and `cScale` into a narrow span. IPOPT will comfortably declare `"Optimal_Solution_Found"` while state components are silently parked on these arbitrary ceilings, yielding unphysical answers with machine-zero defects. 
*Fix*: Add saturation diagnostics `max(lbX - X)` and `max(X - ubX)` into the output struct `out`. Throw an error or flag if any KKT solver limits $<1e-5$ dynamically binding constraints.

**[ROBUSTNESS]** `casadi_energy_freetf.m:101` & `gen_elfo_minfuel.m:120` -- **L6: Continuation-Adaptive Bounds**. `cBox` is rigidly hardcoded to `[0.10 8]` or `[0.15 6]`. During a $T_{max}$ thrust ladder, lowering thrust physically forces larger $t_f$. Because $\tau_f$ is fixed, the slack scalar `cScale` *must* stretch proportionately. The strict `cBox` will render deeper rungs of the $T_{max}$ ladder structurally infeasible. 
*Fix*: Widen `cBox` to a logarithmic span (`[0.01 100]`) or scale it to $\max(T_{base}/T_{max}, 1)$.

**[ROBUSTNESS]** `gen_elfo_minfuel.m:127` (Triage G2/G4 Deferred item) -- **Short-circuited Fallback Cascades**. The continuation `step_solve` tests the loose probe (`rL`) and skips to the tight target (`rT`). If `rT` fails, it halts `ok=false` instead of attempting the `rF` (tight-from-$X_k$) fallback strategy! 
*Fix*: Route `rT` failure directly into the `rF` routine so failure to clean a loose probe doesn't artificially terminate the ladder. Further, widen triage G3 velocity bounds (`[-25, 25]`).

---

### 2. LADDER-READINESS
If you run the $T_{max}$ thrust ladder on this specific pipeline configuration next week, **it will violently stall almost immediately.** As $T_{max}$ shrinks and $t_f$ balloons, three limiters trigger in cascade: (1) `cScale` attempts to stretch the time domain to accommodate more revolutions but artificially saturates against `cBox=6`, returning numerically false "optimal" traps; (2) the fixed $N=4001$ mesh collapses underneath the vastly longer physical duration, leaving bang-bang structures so drastically under-resolved that the KKT solver begins chattering the trajectory; and (3) under-iterated "acceptable" duals from the previous structural strain will prevent IPOPT from regaining warm-start traction. 
**Top 3 Changes required BEFORE starting the ladder:**
1. Dynamically decouple the grid size $N$ derived from $t_f$ continuously alongside physical 1D-time interpolations for all continuations.
2. Widen structural constraints deeply: $cScale \in [0.01, 100]$ and velocity spans to `[-25, 25]`.
3. Reject "acceptable" duals within `step_solve` (strictly check `Solve_Succeeded`) to protect the trajectory spine constraints.

---

### 3. OVERALL VERDICT
I fundamentally trust the published $T_{max}=25$ mN numbers in terms of *physics* because the geometry, two-primary collocation, and homotopy mechanics demonstrably solve the correct target. However, I question the *numerical precision* claims. Unscaled trap-rule defects ($O(\Delta\sigma)$ offset mapped to global $O(1)$) masked by a raw tolerance of `1e-7` mean the trajectory includes continuous fractional integration creep, resolving the dynamics well below the claimed capability limit of typical pseudospectral norms. 

**Single highest-value improvement:**
Explicitly force the collocation trap-rule alignment (`D(:)`) row-scaling with $N$ ($O(\Delta \sigma^{-1})$). This one line of code forces IPOPT's gradient routines to actually respect the raw $10^{-7}$ violation tolerance as a physical dynamic adherence floor, ensuring downstream primer alignment and PMP costates remain mechanically sound.
