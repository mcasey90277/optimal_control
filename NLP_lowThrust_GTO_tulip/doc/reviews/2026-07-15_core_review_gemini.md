The mathematical and architectural redesign in this collocation pipeline is outstanding. The integration of Betts’ sparse `cScale` clock mapping, the continuous two-primary Sundman metric, and the `muGain=0` bypass for the Cartesian retargeting are exceptionally clever and address the physical/numerical walls you encountered perfectly. 

The implementation of `cScale`, the homotopy measure integrals, and the KKT costate (`lamAll(1:9*N)`) pointer logic correctly match the physics and the CasADi stack.

However, there is a lingering foundational bug in the old `minfuel` code that must be resolved, alongside a flaw in the continuation fallback logic that can artificially strand the solver.

Here are the concrete findings:

- **[CORRECTNESS]** `sundman_minfuel/casadi_minfuel_sundman.m:103` — Issue. The problem is physically overconstrained. `tauf` is fixed to `tauf0` (line 103) and terminal time is pinned via `X(8,nN) == tf` (line 123), but there is NO `cScale` slack state. This formulation removes the necessary timing degree-of-freedom, rigidly forcing the spatial trajectory to contort its average geometry simply to ensure the line integral $\int r^p d\tau$ hits the scalar limits precisely. This adds massive artificial stiffness and yields physically sub-optimal paths. Concrete fix: Implement the `cScale` state identical to `casadi_energy_freetf.m` (`dt/dtau = cScale*r1^pSund` with `X(9)` as a floating state) so physical time and $\tau_f$ can untangle appropriately.

- **[ROBUSTNESS]** `elfo/gen_elfo_energy_gravhom.m:187` — Issue. The `step_solve` continuation logic is flawed and discards its own fallback. If the loose probe `rL` succeeds, the tight re-clean `rT` runs. But if `rT` subsequently fails (which happens when bounds shift against the tight `1e-9` barrier), the step returns `FAIL` immediately (line 205). The solver never attempts the robust, fully-tight continuation from `Xk` (`rF`), unnecessarily halting the continuation. Concrete fix: Cascade the fallback safely:
```matlab
% (a) Loose probe
rL = casadi_energy_freetf(..., oL);
if rL.success && rL.maxDefect < 1e-6
    % Tight re-clean
    rT = casadi_energy_freetf(..., rL.X, rL.U, ..., oT);
    if rT.success && rT.maxDefect < 1e-6
        ok = true; Xn = rT.X; Un = rT.U; info = rT; return;
    end
end
% (b) Tight fallback ONLY if rL or rT failed
rF = casadi_energy_freetf(..., Xk, Uk, ..., oF);
if rF.success && rF.maxDefect < 1e-6
    ok = true; Xn = rF.X; Un = rF.U; info = rF; return;
else
    ok = false; Xn = Xk; Un = Uk; info = rF; return;
end
```

- **[ROBUSTNESS]** `elfo/casadi_energy_freetf.m:164` — Issue. State velocity bounds are statically clamped at `[-12, 12]`. GTO departure perigee relative to the Earth can exceed 10.2 ND velocity ($\sim10.2$ km/s). During loose IPOPT exploration (especially under `muGain=0`), the solver regularly inspects off-nominal dynamic bands that momentarily step outside these thresholds. Clamping velocity too closely to the perigee wall invites early search truncation or restoration-phase failure. Concrete fix: Widen velocity bounds to `[-25, 25]` (and conservatively extend position bounds to `[-5, 5]`) in all three CasADi scripts to guarantee bounds only block genuine diverges. 

- **[EFFICIENCY]** `elfo/gen_elfo_energy_gravhom.m:200` — Issue. In the current script, if the tight fallback `rF` is triggered and succeeds, its tight solution (`Xs = rF.X`) is then pointlessly shoved immediately into the `rT` tight re-clean block. Because `rF` and `rT` utilize identical `warmTight = true` specs evaluating at the identical `rvf_s` terminal, this re-computes an already-converged KKT station. Concrete fix: Applying the revised logical `if/else` cascade provided in the Robustness finding above avoids executing `rF` and `rT` back-to-back, saving an entire IPOPT sequence.
