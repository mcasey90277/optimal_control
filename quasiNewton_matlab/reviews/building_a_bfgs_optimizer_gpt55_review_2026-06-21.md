Reviewing checkpoint numerics.Checking inexact-search claims.## Verdict

Strong tutorial. Core math, code alignment, and headline checkpoint are correct. I verified the reference demo in MATLAB R2025b: BFGS gives 20 iterations, `3.275e-14` final gradient, `6.081e-14` error, 0 skips.

## Issues

| Severity | Location | Problem | Fix |
|---|---|---|---|
| Major | `building_a_bfgs_optimizer.tex:412-418` | “This is exactly the Kanamori & Ohara result” overstates the Armijo experiment. The paper’s Table 3 is about line-search noise/accuracy sweeps, not necessarily this default Armijo implementation. | Say “qualitatively mirrors/foreshadows Table 3” unless reproducing their exact noise model. |
| Minor | `building_a_bfgs_optimizer.tex:221-230` | Random BFGS checkpoint can rarely choose very small positive `s'*y`, making `rho` huge and the `~1e-15` secant residual / PD check fragile. | Set `rng(0)` and construct `y = G*s` with SPD `G`, or enforce a margin on `s'*y`. |
| Minor | `building_a_bfgs_optimizer.tex:106-109` | Says references live in parent `quasiNewton_matlab/`, but this review package has them in the same directory. Learners may look in the wrong place. | Match the actual distribution path or say “in the supplied reference folder.” |
| Minor | `building_a_bfgs_optimizer.tex:92-103`, `398-403` | “What You Will Build” omits `dfp_update_H.m`, later introduced as stretch. | Add an optional/stretch row for `dfp_update_H.m`. |
| Minor | `building_a_bfgs_optimizer.tex:321-338` vs `qn_minimize_H.m:34-40,73` | Tutorial names helper `getdef`; reference uses `getfield_default` and also supports `c1`, `bt`. Not wrong, but drift. | Align helper name and mention optional `c1`, `bt` defaults if matching reference is desired. |

## Correct / verified

- BFGS inverse update `VHV' + rho ss'`, secant equation, and PD proof are correct.
- Exact quadratic line search `-(g'p)/(p'Ap)` is correct.
- `grad f = Ax-e`, Hessian `A`, `x* = A\e`, and SPD claim are correct.
- DFP inverse update is correct.
- Checkpoint magnitudes are plausible; headline BFGS numbers match the reference exactly under MATLAB.
