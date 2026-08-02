## Q1 — Hermite–Simpson

**(a)** Correct: `casadi_mintime_dro.m:161,164` implements the separated Hermite–Simpson interpolation and defect with the correct signs and \(h/8\), \(h/6\) factors.

**(b)** Fourth-order *global state accuracy* holds for a fixed, sufficiently smooth control reconstruction and smooth solution. The per-interval endpoint local truncation error is \(O(h^5)\). Bang/coast switches, mesh-scale control variation, or nonsmooth reconstructed controls can reduce observed order.

**(c)** (iii) is the primary explanation: your reset-and-propagate interval endpoint mismatch is a local error, so an approximately \(2^5=32\) reduction is expected. The observed 37.6 is compatible with asymptotic-plus-pre-asymptotic effects, but only if both meshes represent the same limiting trajectory/control. Document it as “fourth-order globally; fifth-order local defect/propagation error.”

**(d)** A free `Um` is the standard and appropriate HS transcription. Replacing it with an arithmetic average generally supplies only an \(O(h^2)\) midpoint-control approximation and can reduce the state scheme’s order. It is not mathematically the *only* valid choice—e.g. a higher-order feasible spherical-direction reconstruction could work—but it is necessary for this formulation and avoids invalid averaged directions.

## Q2 — Control reconstruction

**(a)** No. HS constrains dynamics only at the three interval samples; it does not uniquely imply an intersample control. The quadratic is a reasonable reconstruction convention, not a consequence of the NLP.

**(b)** Correct: at \(w=0\), the basis is \((1,0,0)\); at \(w=1/2\), \((0,1,0)\); at \(w=1\), \((0,0,1)\).

**(c)** Linear control is likewise a reasonable trapezoidal reconstruction convention, but is not uniquely implied by endpoint-control trapezoidal collocation.

**(d)** It changes the off-grid control and therefore the measured propagation error. It agrees at collocation samples, where `alpha` is constrained unit norm, but the NLP did not solve the normalized quadratic control curve.

**(e)** Interval repropagation is a useful mesh-error indicator, but not a global trajectory-error certificate. A Betts-style polynomial residual is complementary: it diagnoses the reconstructed interpolant directly but still depends on a chosen control reconstruction. Use both, plus a continuously propagated full-trajectory endpoint/path-feasibility check.

## Q3 — Gates

**(a)** G1 is useful for local mesh adequacy, but max local error is not a global error bound: resetting every interval suppresses accumulated instability and error. Do not replace it with a sum; add a sequential propagation from the initial state under the reconstruction and compare propagated states with all nodes, especially terminal state and path margins.

**(b)** Make G2 a validation/convergence gate rather than an indirect-reference requirement: require mesh-refinement convergence of objective, trajectory, controls, propagated terminal error, and continuous path feasibility. Retain an independent reference only as optional benchmark evidence where available; it cannot certify global optimality or cover constrained cases.

**(c)** Missing: continuously propagated terminal feasibility, full propagated path feasibility, mesh-refinement convergence, solver-status/lifted-time-spread checks, and control-reconstruction validity. G6 is not universally valid: minimum-time problems can have legitimate coast arcs, particularly with path constraints.

- **[CORRECTNESS]** `dro_residual.m:71`; `certify_dro_mintime.m:69` — G1 claims a position error in km but computes a norm over position, velocity, and mass, then multiplies the combined nondimensional quantity by `lStar`. Velocity and mass-fraction errors cannot be converted to km this way. Compute `norm(Z(end,1:3).'-o.X(1:3,k+1))` for G1; report velocity and mass residuals under separate, dimensionally valid tolerances.
- **[CORRECTNESS]** `dro_residual.m:63-70` — The quadratic control is presented as “the Hermite–Simpson control,” but the transcription defines no intersample control. Normalizing `alpha` and clipping throttle produce a different nonlinear, nonsmooth off-grid control than the raw quadratic. Label this explicitly as a chosen reconstruction; ideally make the reconstruction configurable and assess sensitivity to it.
- **[CORRECTNESS]** `certify_dro_mintime.m:141-152` — G7 is not a continuous altitude guarantee. It resets propagation at every NLP node and only evaluates altitude at nine requested output times per interval; a periselene minimum between samples can be missed. Propagate sequentially from the initial state and locate the minimum with dense-output refinement or an event/root-search on lunar-radius margin.
- **[CORRECTNESS]** `certify_dro_mintime.m:106-108` — `C.pass` ignores G3–G7 by design, despite the report calling a passing solution “CERTIFIED.” A solution may pass G1/G2 while failing terminal feasibility, unit norm, NLP feasibility, or continuous altitude. Either require all mandatory feasibility gates in `pass`, or rename the present verdict to “accuracy/reference agreement passed.”
- **[ROBUSTNESS]** `dro_residual.m:95`; `certify_dro_mintime.m:168` — A quadratic interpolation of three unit directions can pass through or arbitrarily near zero. `al/max(norm(al),eps)` then yields a non-unit vector at zero, silently reducing thrust instead of representing a valid direction. Use a direction interpolation guaranteed nonzero/unit (e.g. a spherical/geodesic reconstruction), or detect and fail such intervals.
- **[ROBUSTNESS]** `dro_residual.m:95`; `certify_dro_mintime.m:168` — Quadratic throttle interpolation can overshoot `[0,1]`; clipping changes the reconstructed control and introduces nonsmoothness. Record reconstruction overshoot and treat it as a failed control-reconstruction validity gate rather than silently clipping it.
- **[ROBUSTNESS]** `casadi_mintime_dro.m:207-208` — Strict symbolic inequalities do not provide a strict feasible margin in an NLP; practical solvers enforce non-strict bounds. Use explicit lower bounds such as `TF >= tfMin` and `X(7,:) >= mMin` chosen for physical/numerical safety.
- **[CORRECTNESS]** `casadi_mintime_dro.m:245-254` — Supplying `X0=[]` with a valid user `U0` silently discards `U0`; supplying `X0` with `U0=[]` reaches `opti.set_initial(U, [])` and can fail. Initialize missing state and control guesses independently.
- **[ROBUSTNESS]** `casadi_mintime_dro.m:274-278` — Opposite adjacent seed directions yield a zero midpoint average; dividing by `eps` leaves a zero direction, violating the midpoint unit-norm initial constraint. Provide a deterministic unit fallback direction.
- **[CORRECTNESS]** `certify_dro_mintime.m:96-97` — Full throttle is not a universal minimum-time feasibility or optimality condition. It is a useful diagnostic for this unconstrained reference solution, but should not be a mandatory gate for constrained or alternative local optima.
- **[ROBUSTNESS]** `certify_dro_mintime.m:26-29`; `certify_dro_mintime.m:90-104` — The documentation says lifted-time spread is a supporting check, but no gate tests `o.tfSpread`. Add a tolerance-scaled lifted-time continuity gate.
- **[STYLE]** `casadi_mintime_dro.m:57-60` — “Hermite-Simpson is fourth order at the same … constraint count” is inaccurate: this separated formulation adds `7N` interpolation equalities, `4N` midpoint controls, and midpoint control constraints. Say it preserves the mesh interval count, not constraint count.
