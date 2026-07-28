# Mesh-convergence plan — Gemini 3.1 Pro review (2026-07-27)

Model: `openrouter/google/gemini-3.1-pro-preview` via crush, prompt inlined (no `-c`).
Reviewed: `docs/superpowers/plans/2026-07-25-mesh-convergence-study.md` (post-GPT-5.6-terra revision) + `foc_check.m`.
Second independent review; the first is `review_2026-07-25_meshplan_gpt56terra.md`.

Scope requested: (1) the theoretical expectations table, (2) Task 0's covector formula,
(3) the Richardson permitted/forbidden split, (4) the branch check. The review engaged
those four and nothing else, so it is silent on the rest of the plan.

- **[CORRECTNESS]** `Theoretical expectations` -- The prediction that switch times on unaligned meshes are "grid-quantized" to $O(h)$ is mathematically incorrect for numerical optimal control formulations using continuous throttles ($s \in [0,1]$). Direct collocation solvers natively formulate precise sub-grid switch placement by leveraging a single partially-throttled intermediate node, seamlessly interpolating the crossing and routinely recovering $O(h^2)$ convergence. Additionally, final mass benefits from integrated error cancellation across many internal arcs (unlike a hard 1-norm trajectory bound), meaning integrated quantities can perform asymptotically better than deterministic $O(h)$ worst-case bounds suggest.

- **[CORRECTNESS]** `Task 0` -- The prompt's concern about the min-fuel terminal cost is misplaced: the problem defines min-fuel purely as a Lagrange integral of throttle (`integral[ s... ]`), meaning mass does NOT appear in the objective and the plan's claim that its exactly mapped transversality is 0 correctly holds. However, the plan's textual definition for `lamHat_f` contains a relative sign error on the objective (quadrature) term: discrete KKT/Hamiltonian equivalence requires it to be $-(h_N/2)L_x(N+1)$ to properly isolate $\nabla \Phi$, rather than added (though this term safely zeros out for the mass variable anyway).

- **[CORRECTNESS]** `Richardson policy` -- Using an arc-by-arc propagation defect directly as an order estimator strictly measures *local truncation error*, which is one polynomial order higher (e.g., $O(h^3)$ for trapezoidal collocation) than the fundamental global discretization error ($O(h^2)$). Concrete correction: The local integrated cell errors must be explicitly divided by the local step size $h_k$ (or appropriately weighted and summed globally) before extracting $p$; otherwise, the reported measurement will falsely inflate the true global convergence rate by $+1$.

- **[ROBUSTNESS]** `Global Constraints` (Branch check) -- Re-solving only the finest level from an independently constructed seed successfully verifies that the chain's terminus has not settled into an artificial warm-chain trap, but it fails to prove that the sequence evaluated a continuous branch across the intermediate meshes. Concrete correction: Explicitly assert that the switch topology (specifically the count) of the independently seated fine mesh strictly matches the topology of the *coarsest* (production) level, failing the check if a masked bifurcation occurred mid-ladder.

VERDICT
PROCEED WITH THE CORRECTIONS ABOVE

The single most important thing to fix first is the `Richardson policy` arc-by-arc defect estimator. You must mathematically correct it by dividing local cell defects by the local step size; otherwise, the plan will inadvertently measure local truncation error instead of global discretization error, systematically inflating and invalidating every order reported by the study.
