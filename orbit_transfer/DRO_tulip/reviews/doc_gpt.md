## R1 — State-constrained Pontryagin theory

**Verdict:** `q=2` is correct for \(C(x)=\rho-\rho_{\min}\) (or its squared equivalent): \(C^{(1)}\) depends on position and velocity, while \(C^{(2)}\) first contains the thrust acceleration. The surrounding indirect-method description is materially oversimplified and the stated jump formula is not convention-complete.

- The multiplier treatment for an order-two state constraint is not generally “append \(\eta(t)C\) to the Hamiltonian.” A standard treatment differentiates to the mixed constraint \(C^{(2)}=0\) on the boundary arc (with tangency conditions \(C=C^{(1)}=0\)), or uses a measure-valued multiplier/direct adjoining formulation. State the convention and formulate it consistently.
- The displayed jump relation can be valid under one sign convention, but not as written without defining: the minimization/maximization Hamiltonian convention; whether feasibility is \(C\ge0\) or \(-C\le0\); whether the multiplier is attached as \(+\nu C\) or \(-\nu C\); and that \(\partial_x C^{(j)}\) is a gradient with respect to the **full state** \((r,v,m)\), not position alone.
- One cannot state that every \(\nu_j\) is nonnegative. The sign restriction belongs to the nonnegative state-constraint measure / appropriate boundary atom under a declared inequality convention; the derived higher-order junction multipliers do not in general each have independent nonnegativity restrictions.
- Costates may have jumps at junctions when the state-constraint measure has atoms. They need not jump at both entry and exit, and an equivalent direct-adjoining formulation can use continuous transformed costates. Say “the conventional costate may jump at junctions, depending on the adjoining convention and junction multipliers,” rather than asserting discontinuity categorically.
- The qualitative comparison is fair if narrowed: direct transcription makes it substantially easier to **pose and numerically explore** a state-constrained problem. It is not “nearly free” in a rigorous sense, because direct methods still require event/mesh refinement and continuous-time constraint verification; indirect methods need not always rely on manually guessed arcs if complementarity, continuation, multiple shooting, or boundary-arc detection is used.

## R2 — Altitude exponent

**Verdict:** The \(-3.5\) radial scaling is correct as a two-body, locally Keplerian heuristic:
\[
|\dddot r|\sim \mu |v|/r^3,\qquad |v|\sim(\mu/r)^{1/2},
\]
hence \(|\dddot r|\sim\mu^{3/2}r^{-7/2}\). It does **not** establish that a fitted residual should scale as altitude\(^{-3.5}\).

A measured exponent of \(-4.24\) differs by about 21% in magnitude from \(-3.5\); call it “consistent with a steep inverse-distance dependence expected from rapidly varying near-lunar dynamics,” not a match. The discrepancy can arise from fitting altitude \(h=r-R_M\), rather than lunar-center radius \(r\); varying step size/final time; velocity variation; CR3BP, Coriolis, and thrust terms; control reconstruction error; and fitting a non-asymptotic regime with one power law.

## R3 — Double reversal

**Verdict:** The final conclusions are visible in the abstract, and the document is admirably explicit about withdrawn evidence. Nevertheless, the body currently reads as a sequence of self-contradictions because the initial strong conclusions are stated before their retraction. More importantly, the purported “stronger” re-establishment of ill-posedness is not sufficient evidence: its \(N=800\) Hermite–Simpson trajectory has a reported worst continuous residual of \(1{,}017{,}917\) km.

Keep the history, but put the final status first in each affected section, then move the failed inference and its correction into a clearly titled “Investigation history” subsection or appendix. Add forward pointers at the first historical claim, not only at its later retraction.

## Findings

- **[CORRECTNESS]** `Section 6, lines 311–326` — “For minimum time one expects \(u^\star\equiv1\): there is no reason to coast” is not a general Pontryagin result. Endpoint geometry and dynamics can make coasts optimal even in a minimum-time problem. Concrete fix: say that full thrust is an empirical property of these computed extremals, supported by the reported \(\min u=1\), not a consequence of the objective alone.

- **[COMPLETENESS]** `Section 6, lines 311–326` — The indirect PMP summary omits the free-final-mass transversality condition \(\lambda_m(t_f)=0\) and the free-final-time condition \(H(t_f)=0\) for the stated autonomous convention. Concrete fix: include both, with the Hamiltonian convention used to obtain the displayed switching function.

- **[CORRECTNESS]** `Section 8.2, lines 518–540` — The state-constraint multiplier and junction discussion is not technically adequate for an order-two pure state constraint. Concrete fix: define \(C\), derive \(C^{(1)}\) and \(C^{(2)}\), state boundary-arc tangency conditions, declare the Hamiltonian/inequality convention, and replace the unconditional jump/discontinuity claim with convention-qualified measure-multiplier junction conditions.

- **[CLARITY]** `Section 8.2, lines 510–516` — “Millions of rows already” is incompatible with the documented NLP sizes, which are in the thousands to tens of thousands of constraints. Concrete fix: replace with “an additional sparse block of inequality rows,” avoiding a numerical claim that is false for this implementation.

- **[CORRECTNESS]** `Section 8.1, lines 502–504; Section 7.1, lines 431–434; Section 8.2, lines 545–556` — The document contradicts itself on enforcement points. It first says Hermite–Simpson enforces the floor at nodes **and midpoints**, then says the transcription adds only \(N+1\) rows and repeatedly says it is imposed “at the nodes only.” Concrete fix: distinguish schemes explicitly: trapezoid uses nodes; Hermite–Simpson uses nodes plus collocation midpoints, if that is implemented; neither establishes continuous-time feasibility.

- **[CORRECTNESS]** `Section 10, lines 712–777` — A machine-feasible Hermite–Simpson node inside the physical Moon does not prove that the continuous unconstrained OCP has no minimum. The same \(N=800\) solution is later reported to have a worst propagated interval discrepancy of \(1{,}017{,}917\) km. Concrete fix: withdraw “therefore has no minimum” unless supported by a converged continuous sequence of feasible trajectories with decreasing allowed periselene and decreasing \(t_f\), or by an analytical existence/nonexistence argument.

- **[CORRECTNESS]** `Section 10, lines 765–772` — “No quadrature argument produces a node inside the primary” is false. A collocation NLP can place nodes anywhere permitted by its constraints, including inside a body, while still having large inter-node dynamic error. Concrete fix: describe the result as evidence that the **discrete transcription** exploits the missing physical exclusion constraint; do not treat it as proof about the continuous trajectory without propagation.

- **[CLARITY]** `Section 10, lines 715–754` — The section heading and the strong assertions at lines 732–754 state the ill-posedness conclusion before warning that the evidence is immediately retracted at line 756. Concrete fix: open with a boxed final-status statement and label the table “historical trapezoidal evidence, later invalidated by continuous residual checks.”

- **[CORRECTNESS]** `Section 10, lines 785–794` — Adaptive integration does not make an indirect method structurally immune to lunar-collision or near-singularity failures, nor does it prove that a near-surface path would be “visibly” resolved. Shooting can still converge to incorrect/local extremals or fail near singular dynamics. Concrete fix: limit the claim to adaptive integration avoiding a fixed collocation mesh’s particular inter-node defect mechanism.

- **[CORRECTNESS]** `Section 11, lines 805–850` — \(R_k\) is the norm of a seven-component nondimensional state discrepancy, mixing position, velocity, and mass. It cannot be converted to kilometres by multiplying the total norm by \(L\), and “hundreds of thousands of kilometres” is therefore not justified by the displayed definition. Concrete fix: report separate \(R_r\) in km, \(R_v\) in km/s, and \(R_m\), or explicitly define a physically weighted norm and retain it as dimensionless.

- **[CORRECTNESS]** `Section 11, lines 841–850` — The text says the measured altitude exponent “matches” the \(-3.5\) theoretical prediction, but \(-4.2\) differs materially from \(-3.5\), and the derivation is for radius \(r\), not altitude \(r-R_M\). Concrete fix: call the result qualitatively consistent with steep near-lunar scaling, fit against lunar-center radius as a comparison, and discuss the observed discrepancy.

- **[CLARITY]** `Section 11, lines 805–816` — The diagnostic integrates linearly reconstructed controls, but the transcription itself does not uniquely prescribe continuous controls between nodes. Concrete fix: call \(R_k\) the residual for the specified **linear-control reconstruction**, not unqualified “the true local error” of the discrete solution.

- **[CORRECTNESS]** `Section 11, lines 861–865` — “Three routes are available and none has been tried here” contradicts lines 384–457 and lines 702–710, which document Hermite–Simpson as the tried, successful higher-order route. Concrete fix: change to “Before the Hermite–Simpson result, none had been tried,” or remove the sentence and point forward to Section 7.

- **[CORRECTNESS]** `Section 7, lines 466–470; Section 7, lines 438–450; Section 12.3, lines 943–945` — The \(N=400\) Hermite–Simpson final times conflict: \(4.0173381\) in the scheme comparison but \(4.6808938\) in both the warm-start discussion and the final ladder. Concrete fix: identify the differing initialization/problem variant, label both runs, or correct the erroneous table entry.

- **[CORRECTNESS]** `Section 12.2, lines 928–933` — \(4.0152501\) and \(4.0152425\) agree to five significant figures, not six. Also, agreement in final time and endpoints does not show agreement “on a trajectory.” Concrete fix: say “five significant figures” and report a trajectory/control comparison metric before claiming trajectory agreement.

- **[CLARITY]** `Section 12, lines 875–957` — The document calls an unconstrained result “certified” after concluding that the unconstrained problem has no minimum. The gate certifies local transcription accuracy and agreement with one indirect extremal; it does not certify global optimality or repair the OCP. Concrete fix: rename it “certified reproduction of the reference local extremal” until a constrained problem and appropriate reference/gate are established.

- **[COMPLETENESS]** `Section 8.3, lines 565–602; Section 14, lines 1028–1034` — No constrained solution is certified, and the only reported 100 km-floor runs are inactive and mesh-dependent. Concrete fix: state prominently that the constrained OCP is formulated but not yet solved/certified; do not imply that the 1600-node unconstrained result validates the mission-constrained problem.

- **[CORRECTNESS]** `Section 6.2, line 354` — \(m_k>0.05\) is a strict inequality, which a standard NLP cannot impose. Concrete fix: write the implemented non-strict bound, presumably \(m_k\ge0.05\).

- **[CLARITY]** `Section 9, lines 611–614` — The stated variable and defect counts describe trapezoidal transcription, not the preceding Hermite–Simpson transcription with midpoint states, controls, interpolation equations, and Simpson defects. Concrete fix: label the count as trapezoidal, or provide separate accurate counts for each scheme.

- **[CLARITY]** `Section 7, lines 387–413` — “Fourth-order scheme” needs a qualification. Hermite–Simpson has fourth-order state accuracy under suitable smoothness and control representation assumptions; independently optimized controls at nodes and midpoints do not automatically make every OCP quantity fourth-order accurate. Concrete fix: say “a fourth-order Hermite–Simpson state transcription for smooth arcs,” then report observed convergence separately.

- **[CLARITY]** `Section 9, lines 616–621; Section 12, lines 882–900` — The reported NLP defect, terminal error, and norm error are necessary solver checks but are not independent validation when recomputed from the same discrete variables. Concrete fix: explicitly call them discrete feasibility checks and reserve “continuous validation” for propagated residuals, constraint checks, and mesh-refinement evidence.

- **[STYLE]** `Section 7 and Section 10, lines 625–710 and 712–777` — The reversal narrative is honest but too long in the main technical flow. Concrete fix: retain a short “final finding” first, move the chronological failed hypotheses into an appendix or shaded “investigation record,” and leave one forward reference from each superseded claim.

- **[COMPLETENESS]** `Reproducing this, lines 1003–1049` — Reproduction instructions lack MATLAB release, IPOPT/MUMPS versions/build, exact IPOPT options, tolerances for all integrations, reference-solution provenance/version, and repository commit identifiers. Concrete fix: add a reproducibility table with versions, commit hashes, solver options, warm-start source, and expected outputs/checksums.


