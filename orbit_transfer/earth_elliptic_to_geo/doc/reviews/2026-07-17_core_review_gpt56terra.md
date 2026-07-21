## Q1 — dual → costate map

**Yes: the diagnostic’s naive nodal map is wrong.** `lamDef(:,k)` is an **interval-defect multiplier**, while `U(:,k)` is a **node control** shared by defects \(k-1\) and \(k\). `casadi_lt_2body.m:162-166` compares \(\alpha_k\) only with interval \(k\)'s \(\lambda_{v,k}\); `verify_pmp_2body.m:34-43` does the same.

Let \(h_k=\Delta\tau_k\), \(q_j=c_j\kappa_j\), \(F_j=[q_j f_j;0]\), and define defects
\[
D_k=x_{k+1}-x_k-\frac{h_k}{2}(F_k+F_{k+1})=0.
\]
Take the NLP Lagrangian as
\[
\mathcal L=J+\sum_{k=1}^N \Lambda_k^\top D_k+\cdots,
\qquad
\Lambda_k=[\Lambda_{r,k};\Lambda_{v,k};\Lambda_{m,k};\Lambda_{t,k};\Lambda_{c,k}].
\]
CasADi may return the global negative of this convention; that does not affect a direction check once one global sign is fitted.

For an interior node \(j=2,\ldots,N\), only defects \(D_{j-1}\) and \(D_j\) contain \(\alpha_j\). Since
\[
\frac{\partial F_{v,j}}{\partial\alpha_j}
=q_j\,\frac{T_{\max}}{m_j}s_jI,
\]
and the fixed-\(t_f\) objective has no direct \(\alpha\) dependence,
\[
0=
-\frac{T_{\max}s_j}{2m_j}
\left[
h_{j-1}q_j\Lambda_{v,j-1}
+h_jq_j\Lambda_{v,j}
\right]
+2\nu_j\alpha_j
+\text{active alpha-bound terms}.
\]
On a burn arc, away from artificial alpha box bounds,
\[
\boxed{
\alpha_j\ \parallel\
\pm\left(h_{j-1}\Lambda_{v,j-1}+h_j\Lambda_{v,j}\right).
}
\]
The **common node-local factor** \(q_j=c_j\kappa_j\) cancels from the direction. Thus the correct discrete primer is the **time-step-weighted adjacent interval average**, not `lamDef(:,j)`.

For a uniform \(\tau\)-mesh, \(h_{j-1}=h_j\), this reduces exactly to
\[
\boxed{\lambda_{v,j}^{\rm disc}=\tfrac12(\Lambda_{v,j-1}+\Lambda_{v,j})}
\]
for interior control nodes. Therefore, experiment (b) was the correct first test, and its tiny improvement means that **this indexing error alone does not explain 17–24°**. For a nonuniform mesh, it must instead be
\[
\lambda_{v,j}^{\rm disc}=
\frac{h_{j-1}\Lambda_{v,j-1}+h_j\Lambda_{v,j}}
{h_{j-1}+h_j}.
\]
There is no additional \(c\)-, \(\kappa\)-, or \(q\)-weight in that *nodal direction* because both adjacent contributions use the same \(q_j\). A \(q\)-weighted average would be wrong here: \(F_j\), not \(F_{j\pm1}\), is differentiated with respect to \(\alpha_j\).

For velocity, the internal-node stationarity is
\[
\begin{aligned}
0={}&\Lambda_{v,j-1}-\Lambda_{v,j}\\
&-\frac{h_{j-1}}2
\left(\frac{\partial F_j}{\partial v_j}\right)^\top\Lambda_{j-1}
-\frac{h_j}2
\left(\frac{\partial F_j}{\partial v_j}\right)^\top\Lambda_j
+\text{box/boundary terms}.
\end{aligned}
\]
Here
\[
\frac{\partial F_j}{\partial v_j}
=
\begin{bmatrix}
q_j I_3\\0_{6\times3}
\end{bmatrix},
\]
because \(\dot r=v\), while two-body acceleration has no velocity dependence. Hence
\[
\boxed{
0=\Lambda_{v,j-1}-\Lambda_{v,j}
-\frac{q_j}{2}\left(h_{j-1}\Lambda_{r,j-1}+h_j\Lambda_{r,j}\right)
+\cdots .
}
\]
This is the discrete adjoint coupling between \(\lambda_v\) and \(\lambda_r\).

The clock dependence matters in **\(r\)-stationarity**, not directly in alpha stationarity:
\[
\frac{\partial F_j}{\partial r_j}
=
c_j\left[
\kappa_j\frac{\partial f_j}{\partial r_j}
+f_j(\nabla\kappa_j)^\top
\right],
\qquad
\nabla\kappa_j=1.5\,\|r_j\|^{-1/2}r_j
\]
(up to the stated radius softening). The second term contracts with **all first eight defect multipliers**:
\[
\Lambda_j^\top c_j f_j(\nabla\kappa_j)^\top
=
c_j(\Lambda_{1:8,j}^\top f_j)\nabla\kappa_j^\top.
\]
So the Sundman clock does couple the discrete \(r\)- and \(v\)-adjoints. But the CR3BP sibling has this same \(\kappa(r)\) mechanism, so it cannot by itself explain why only this solver has bulk rotation.

The cScale defect row is
\[
D_{c,k}=c_{k+1}-c_k.
\]
It has **zero derivative with respect to \(\alpha_j\) and \(v_j\)**. Therefore \(\Lambda_{c,k}\) contributes neither to the alpha stationarity nor to the velocity stationarity above. It enters only cScale stationarity:
\[
0=
\Lambda_{c,j-1}-\Lambda_{c,j}
-\frac12\!\left[
h_{j-1}\kappa_j f_j^\top\Lambda_{1:8,j-1}
+h_j\kappa_j f_j^\top\Lambda_{1:8,j}
\right]
+\frac{\partial J}{\partial c_j}
+\cdots .
\]

**[CORRECTNESS]** `earth_elliptic_to_geo/casadi_lt_2body.m:162-166` and `earth_elliptic_to_geo/verify_pmp_2body.m:34-43` — the primer check assigns interval \(k\)'s defect dual directly to node \(k\)'s control. Use the adjacent-step weighted formula above, with one-sided values only at endpoints.

Because your simple adjacent average remains badly misaligned, the next decisive test is not another mapping variant:

1. Build the exact NLP Lagrangian gradient using the returned duals, including cone, terminal, equality, and bound multipliers.
2. At each interior burn node evaluate the **tangential** alpha residual
   \[
   (I-\alpha_j\alpha_j^\top)
   \frac{\partial\mathcal L}{\partial\alpha_j}.
   \]
3. If it is machine-small while the stated primer is rotated, the extracted `opti.dual(conDef{k})` is not the multiplier convention/order assumed by the diagnostic, or another multiplier/sign/scaling issue is present. If it is not small, IPOPT/CasADi’s exposed duals are not a usable KKT certificate under this Opti formulation.

Your cScale-elimination experiment is sound as a **formulation-equivalence test**. With converged \(c^\*\), set
\[
\tau'=c^\*\tau,\qquad \tau_f'=c^\*\tau_{f0},\qquad \frac{dt}{d\tau'}=\kappa(r),
\]
and remove the ninth state and its zero defect row. The primal defects for the eight retained states become numerically identical:
\[
D_k=x_{k+1}-x_k-\frac{\Delta\tau'_k}{2}
[\kappa_k f_k+\kappa_{k+1}f_{k+1}].
\]
At the same primal solution, the non-cScale defect Jacobian columns are also identical. Therefore its retained defect duals should agree, up to global sign/convention, with the 9-state formulation if both are reporting a well-conditioned KKT multiplier. If the reduced formulation produces an aligned primer, it strongly implicates numerical multiplier recovery or formulation conditioning caused by the cScale column; it does **not** make a single-interval map correct.

## Q2 — 5 N min-time stall

Your diagnosis is credible, but it is not the top physical/numerical cause.

| Rank | Candidate cause | Why 5 N is different | Cheapest decisive test |
|---|---|---|---|
| 1 | **Wrong/bad basin at doubled revolution count** | 5 N needs about twice the duration and revolutions. The cold tangential spiral reaches the right energy but is not close to the correct multi-revolution rendezvous/control topology. Cartesian defects accumulate phase error across about nine revolutions. | Freeze the final longitude to the cold seed’s longitude and inspect whether stage 1 reaches \(<10^{-8}\) defects with `warmTight=false`, 6000 iterations. If it cannot, the seed/basin is the primary failure before the manifold is introduced. |
| 2 | **Tight IPOPT restart from a \(5\times10^{-3}\)-infeasible point** | `warmTight=true` activates monotone barrier mode and near-zero bound pushes at `casadi_lt_2body.m:127-131`. Those settings are appropriate for a nearby feasible KKT point, not an infeasible restoration iterate. The observed restoration followed by `Infeasible_Problem_Detected` is consistent with this. | Repeat exactly round 2 with `warmTight=false`, adaptive barrier, `maxIter=6000`; preserve the same primal seed but do **not** warm-start multipliers. This is the cheapest and most discriminating test. |
| 3 | **Long-horizon Cartesian phase conditioning** | The same nodes/rev resolves the local orbit, but does not make the global KKT system equally conditioned: twice as many defect blocks, phase-sensitive terminal constraints, and more opportunities to converge to the wrong winding basin. | Solve the 5 N problem on the same physical seed with a fixed final rendezvous, then homotope terminal longitude/manifold release gradually. Compare reduced-Hessian/KKT conditioning or simply defect progress against the direct manifold solve. |
| 4 | **Free-longitude manifold includes unwanted branches** | The five conditions define both prograde and retrograde GEO circles; they do not impose \(h_z>0\). It also admits every longitude. A nearby bad orbital-phase branch is plausible. | Record \(h_z=(r\times v)_z\) at every iterate/final point. Add \(h_z(end)\ge h_{\min}>0\), or temporarily impose a fixed prograde GEO endpoint, and compare. |
| 5 | cScale range | A larger \(t_f\) can shift cScale, but `[0.05,20]` is generous unless `tauf0` is poorly normalized. | Print `min/max(out.X(9,:))` and IPOPT bound multipliers. If cScale is interior and constant, eliminate it as a cause. |
| 6 | `t<=300` or `|v_i|<=8` boxes | The reported 5 N \(t_f\sim44\) ND is far below 300. Initial-perigee speed is also well below 8 ND; reasonable transfer speeds should remain interior. | Print max absolute `t` and `v` and their bound multipliers. This should be nearly free. |

So: **apply your adaptive-barrier / 6000-iteration change first**, but do not expect it alone to cure a topology/phase basin issue. It is a proper continuation policy:

- `warmTight=false` while `maxDefect >= 1e-6`;
- use tight warm starts only after a feasible or nearly feasible iterate;
- do not reuse dual multipliers from failed/restoration solves;
- retain the new primal iterate only if it improves feasibility.

**[ROBUSTNESS]** `earth_elliptic_to_geo/run_mintime.m:196-197, 314-345` — every continuation call uses `warmTight=true`, including restarts from visibly infeasible iterates. Select barrier/warm-start policy by feasibility rather than call number.

**[ROBUSTNESS]** `earth_elliptic_to_geo/casadi_lt_2body.m:103-109` — the free-GEO manifold does not encode prograde angular momentum. Add a prograde terminal inequality or use a longitude/branch homotopy during anchor construction.

The rejected thrust-stretch seed is indeed topology-flawed: rescaling only `X(8,:)` at `run_mintime.m:229` changes neither orbital phase nor revolution count. It cannot turn a 4.5-revolution path into a nine-revolution one.

## Q3 — baseline integrity

Nothing shown invalidates the reported 10 N Cartesian baseline. The reported agreement with the independent min-time anchor and Haberkorn-scale fuel/revolution/switch structure is strong evidence that the **primal result** is valid; the dual anomaly is a PMP-certificate problem until proven otherwise.

- Dynamics are correct and MX-safe: explicit squared radius and power at `lt2b_rhs_time.m:20-23`; no symbolic `norm`, `abs`, or `max`.
- The Sundman equations are correct: `f=[cScale*kappa*f_time;0]` at `casadi_lt_2body.m:54-60`, with the cScale zero row included in every defect at `:64-69`.
- The physical-time objective measure is correct: \(dt=c\kappa\,d\tau\), trapezoid-integrated at `:117-118`.
- Initial conditions and both terminal modes are correctly imposed at `:95-109`.
- Nondimensional thrust, exhaust velocity, and postprocessed mass/\(\Delta v\) units are correct at `kepler_lt_params.m:24-31` and `casadi_lt_2body.m:175`.
- Free final mass is not terminal-pinned; its near-zero transversality is consistent with that design.

**[ROBUSTNESS]** `earth_elliptic_to_geo/casadi_lt_2body.m:147-152` — the post-solve defect check uses `norm(r)^pSund`, while the NLP uses \((r^\top r+10^{-12})^{pSund/2}` at `:56-59`. Recheck with the identical softened expression. This is negligible at your radii, but a verifier should certify the actual NLP.

The 1376.74 kg / 0.5 kg successor-validation gate is reasonable. Also require: terminal residual, independent Cartesian reconstruction from MEE, mass-flow quadrature agreement, and convergence under \(N=600\rightarrow1200\). Do not gate the successor on the current raw-dual primer metric until Q1 is resolved.

## Q4 — MEE with \(L\) as independent variable

**(a) Sound, with a monitored positivity condition.** In MEE,
\[
\dot L=
\sqrt{\frac{\mu}{p^3}}\,w^2+
\sqrt{\frac{p}{\mu}}\,
\frac{h_x\sin L-h_y\cos L}{w}\,a_n,
\qquad
w=1+e_x\cos L+e_y\sin L.
\]
Use \(d(\cdot)/dL=\dot{(\cdot)}/\dot L\), carrying \(t\) with \(dt/dL=1/\dot L\). For this prograde, low-thrust transfer, the Kepler term is strongly positive and the normal-thrust correction is small; \(\dot L\to0\) or reversal is not physically expected. Still enforce/log a margin such as \(\dot L\ge \dot L_{\min}>0\). If it approaches zero, the transcription becomes singular and must switch independent variable or regularization.

**(b) \(L\) is the clock.** A separate radial Sundman clock is normally unnecessary. The denominator \(1/\dot L\) is the state/control-dependent time dilation. You may still use a fixed normalized mesh \(\sigma=L/L_f\), with \(L_f\) free only through a sparse scaling-state device if necessary.

**(c) Main traps.**
- \(w\to0\): near-parabolic/very eccentric geometry; not expected for this initial ellipse but must be guarded.
- \(p\to0\), \(m\to0\), or \(\dot L\to0\): true transcription singularities.
- Equinoctial \(h_x,h_y\) remain nonsingular at \(i=0\), which is exactly why they are appropriate here; they become poorly conditioned only near \(i=\pi\).
- Longitude wrapping is not a state discontinuity if \(L\) is unwrapped and monotone; use many-revolution \(L_f\), not a modulo-\(2\pi\) terminal variable.
- The terminal circular-equatorial target leaves longitude free, so retain an explicit prograde branch condition.

**(d) Node floor.** Fifteen nodes/rev is an aggressive lower bound, not a safe production number, for trapezoid with two discontinuous switches/rev. A practical initial range is **25–40 nodes/rev**, then demonstrate mass/switch-time convergence at 20, 30, and 40 nodes/rev. At 754 revolutions, 30 nodes/rev is about 22.6k intervals: still far below Cartesian’s orbital-resolution burden, but not trivial.

Trapezoid smears a switch because a jump located between nodes must be represented by a long fractional-throttle interval. Mesh refinement or switch-aware mesh adaptation is more valuable than merely changing the state coordinates. Hermite–Simpson improves smooth-arc accuracy but does not remove a discontinuity-resolution floor. Global pseudospectral methods can be worse for discontinuities unless segmented at switches. Best practical upgrade: hp/adaptive mesh refinement with interval subdivision around switching-function zeroes, or a multi-phase formulation with switch locations as variables after an initial regularized solve.

**(e) The plan substantially weakens, but does not completely defeat, the paper’s direct-method objection.** MEE+\(L\) removes the need to resolve unforced Kepler oscillation in five states, so the problem becomes plausibly tractable at 0.1 N with sparse NLP and continuation. The remaining challenge is not raw state dimension; it is reliably locating roughly 1500 switching events and maintaining a favorable KKT basin across 754 revolutions. A direct Cartesian transcription remains predictably unsuitable; an MEE, longitude-domain, adaptive/switch-aware direct method is a materially different proposition.
