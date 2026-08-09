# Task 2 Report: pdg_dynamics with Analytic Jacobians + Complex-Step Test

## Summary

Successfully implemented 3-DOF powered-descent dynamics (`pdg_dynamics.m`) with analytic Jacobians and complex-step verification, following strict TDD methodology and pumpkyn coding discipline.

## Files Created

1. **`lib/pdg_dynamics.m`** — Core dynamics function
   - Signature: `[xdot, A, B] = pdg_dynamics(x, T, P)`
   - State vector: `x = [r(3); v(3); m]` [7x1, SI units]
   - Outputs: `xdot` (7x1), `A = ∂f/∂x` (7x7), `B = ∂f/∂T` (7x3)
   - Physics: Vacuum + optional drag (`P.drag.on` flag)
   - Complex-step safe: all magnitudes computed as `sqrt(sum(.^2))`, no `norm/abs/max`
   - Full pumpkyn-style header with Purpose, INPUTS/OUTPUTS (with sizes), REFERENCES

2. **`tests/test_dynamics_jac.m`** — Complex-step verification harness
   - Perturbs each of 10 inputs (7 state + 3 thrust) with `h=1e-30i`
   - Compares `imag(f)/h` against analytic Jacobians A and B
   - Tests both vacuum and drag branches (vacuum default, drag opt-in)
   - Uses nonzero velocity test point to avoid singularity at `v=0`
   - Tolerance: 1e-12 (relative error)

## TDD Phases

### RED Phase
```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_dynamics_jac"
```
**Result:** ✓ Failed as expected — `Unrecognized function or variable 'pdg_dynamics'`

### GREEN Phase
**Same command after implementation:**
```
test_dynamics_jac PASS (vacuum + drag)
```
✓ All assertions passed (1e-12 tolerances on both A and B, both drag modes)

## Test Coverage

- **Vacuum (Phase 1 default):** Jacobians verified without drag effects
- **Drag (Phase 2, opt-in):** Full nonlinear drag model tested
  - Density model: `rho(z) = rho0 * exp(-z / H)`
  - Acceleration: `aD = -(kD / m) |v| v` where `kD = 0.5 * rho * Cd * A`
  - Jacobians: `∂aD/∂v`, `∂aD/∂z`, `∂aD/∂m` all verified

## Physics Implemented

**State derivative:**
- `r_dot = v`
- `v_dot = g + T/m + aD` (aD = drag acceleration, zero if `P.drag.on = false`)
- `m_dot = -|T| / (Isp * g0)` (mass consumption from thrust)

**Analytic Jacobians:**
- `A(1:3, 4:6) = I₃` (position-velocity coupling)
- `A(4:6, 7) = -T/m²` (thrust direction scaled by mass inverse-square)
- `A(4:6, 4:6)` = drag derivative (if enabled): `-(kD/m)(|v|I + vv'/|v|)`
- `A(4:6, 3)` = altitude dependence of drag density: `(kD|v|/H/m) v`
- `A(4:6, 7)` = mass dependence of drag: `+kD|v|v/m²`
- `B(4:6, :) = I₃/m` (thrust-to-acceleration)
- `B(7, :) = -T' / (|T| Isp g0)` (thrust direction scaled by thrust magnitude)

## Discipline Checklist

✓ Complex-step safe (no `norm`, `abs`, `max`)
✓ No `i`/`j` loop variables
✓ Full pumpkyn-style header
✓ MATLAB R2025b only
✓ TDD: RED → GREEN confirmed
✓ Staged only the two files created (verified via `git diff --cached`)
✓ Commit on main with proper prefix and co-author footer
✓ Test point has nonzero velocity (avoids `v=0` singularity in drag Jacobian)

## Commit Record

```
3d6fc23 booster_landing: 3-DOF dynamics + analytic Jacobians, complex-step verified
```

Files:
- `booster_landing/lib/pdg_dynamics.m` — 49 lines
- `booster_landing/tests/test_dynamics_jac.m` — 37 lines

## Quality Assurance

- **Numerical correctness:** Complex-step differences all <1e-12 (machine-precision level) relative error
- **Both physics branches:** Vacuum and drag both verified in single test run
- **Edge-case handling:** Test point at x(3)=1500 m altitude, v=[-25; 5; -140] m/s (nonzero, realistic descent velocity)
- **Code clarity:** Implementation exactly mirrors physics equations; no obfuscation

## Concerns

None. Function is complete, tested at tolerance 1e-12 across both physics branches, and ready for downstream tasks (ode45 integration and LQR linearization).

## Next Task

Task 3: `solve_pdg_colloc` (Hermite-Simpson NLP) can now consume `pdg_dynamics` as the forward model with verified Jacobians.
