Reviewing code paths now.## Findings

- **[CORRECTNESS]** `run_sundman_from_seed.m:53` -- This file still resamples onto a uniform `sigma` mesh, contradicting the required NO-RESAMPLE strategy. Replace lines 53-58 with the homotopy pattern: use `sigma=tau/tauf0` after `unique`, `X0=[Xs(:,ku); tSeed(ku).']`, `U0=[alpha(:,ku); s_seed(ku)]`, and ignore `N`.
- **[ROBUSTNESS]** `casadi_minfuel_sundman.m:48` -- Earth gravity denominator is unguarded: `d3=(dd.'*dd)^1.5`. If the intent is bounded regularized perigee terms, use the guarded distance consistently, e.g. `d3 = r1^3`, or form the product as `kappa*dd/d3` using guarded denominators.
- **[ROBUSTNESS]** `run_sundman_from_seed.m:56` -- `aG = aG./sqrt(sum(aG.^2,1))` can divide by zero after PCHIP overshoot/cancellation. Guard with `max(norm,1e-9)` and fill tiny-norm columns from neighboring valid directions.
- **[ROBUSTNESS]** `run_sundman_homotopy.m:65` -- A failed/loose homotopy step always becomes the next warm start. Add the same `out.success && out.maxDefect<1e-6` guard or abort instead of poisoning later eps steps.
- **[STYLE]** `casadi_minfuel_sundman.m:13` -- Header says `tau_f` is free, but implementation fixes it at `tauf=tauf0` in line 70. Update the header to match the dense-KKT avoidance design.
- **[STYLE]** `casadi_minfuel_sundman.m:18` -- Header omits the `epsilon` input and all four reviewed functions lack a `REFERENCES` block. Add CasADi/IPOPT, Sundman, Bertrand-Epenoy, and CR3BP references.

## Confirmed correct

- **[CORRECTNESS]** `casadi_minfuel_sundman.m:54` -- Multiplying the whole RHS by `kappa`, including `mdot` and `dt/dtau=1`, correctly implements Sundman dynamics.
- **[CORRECTNESS]** `casadi_minfuel_sundman.m:49-50` -- CR3BP gravity and Coriolis signs match the existing project dynamics.
- **[CORRECTNESS]** `casadi_minfuel_sundman.m:74` -- Non-uniform `dsig` trapezoidal defect weighting is correct; with fixed `tauf=tauf0`, `tauf*dsig_k` is the interval `tau` step.
- **[CORRECTNESS]** `casadi_minfuel_sundman.m:95-98` -- Homotopy objective is correct: `eps=1` gives `∫s^2 dt`, `eps=0` gives `∫s dt`, equivalent to minimizing propellant up to a positive constant.
- **[CORRECTNESS]** `casadi_minfuel_sundman.m:88-90` -- Fixed `tauf` plus `t(0)=0`, `t(end)=tf` is well-posed for a compatible Sundman mesh; it enforces `∫kappa dτ=tf` without making `tau_f` dense.
- **[CORRECTNESS]** `run_sundman_homotopy.m:49-51` -- `dtau = dt.*0.5.*(1./kap_k+1./kap_{k+1})` is the correct trapezoidal quadrature of `dt/kappa`.
- **[CORRECTNESS]** `casadi_minfuel_sundman.m:78` -- Unit thrust-direction equality is correct for cone-eliminated `alpha`.
- **[ROBUSTNESS]** `casadi_minfuel_sundman.m:83-86` -- Bounds are explicit two-sided constraints; no malformed chained MATLAB bounds remain.
- **[ROBUSTNESS]** `run_sundman_tail.m:39-45` -- Tail guard correctly avoids advancing the warm start or overwriting the certified file on failed/loose steps.
