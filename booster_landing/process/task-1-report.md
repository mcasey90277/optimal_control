# Task 1 Report: Skeleton, setup_paths, booster_params + sanity test

## Summary
Successfully implemented the foundation scaffold for the booster-landing optimal-control campaign. All TDD steps executed: RED (test fails), GREEN (test passes), commit. Files follow the pumpkyn-style module pattern with full function headers.

## Implementation

### Files Created
1. **`/Users/msc/Desktop/optimal_control/booster_landing/setup_paths.m`** (16 lines)
   - Adds campaign root, lib/, certify/, viz/, tests/ to MATLAB path
   - Adds CasADi 3.7.0 from $HOME
   - Called at test startup to establish module linkage

2. **`/Users/msc/Desktop/optimal_control/booster_landing/lib/booster_params.m`** (50 lines)
   - Single-struct parameter source of truth (pumpkyn pattern)
   - 23 fields spanning vehicle (mdry, m0, Tmax, Tmin, Isp, g0, gvec), BCs (r0, v0), constraints (gs_deg, theta_max_deg), solver (N, Nconv, tf_lo, tf_hi), atmosphere (drag.*), and MC (pad_radius, vtd_max, seed)
   - Key Falcon-9 estimates: mdry=25,600 kg, m0=30,000 kg, Tmax=845 kN, Tmin=338 kN (40% throttle)

3. **`/Users/msc/Desktop/optimal_control/booster_landing/tests/test_params.m`** (21 lines)
   - Asserts three campaign invariants:
     - m0 > mdry (landing propellant must exist)
     - min-throttle T/W at dry mass > 1 (hoverslam forced; cannot hover)
     - All bounds correctly ordered (Tmin < Tmax, tf_lo < tf_hi, gs_deg in (0,90))
     - Phase 1 defaults to vacuum (drag.on=false)
   - Self-bootstrapping: addpath + setup_paths on first line

4. **`/Users/msc/Desktop/optimal_control/booster_landing/.gitignore`** (3 lines)
   - Ignores results/* but preserves results/.gitkeep
   - Ignores MATLAB auto-save .asv files

5. **`/Users/msc/Desktop/optimal_control/booster_landing/results/.gitkeep`**
   - Placeholder file to track empty results/ directory in version control

### Directories Created
- `certify/` — gates layer (Tasks 5+)
- `viz/` — plotting and movie functions (Tasks 9+)
- `results/` — output directory (populated by solvers)

## TDD Evidence

### RED: Test Fails (Step 2)
```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); addpath('tests'); test_params"
```
**Output (failure, as expected):**
```
{setup_paths is not found in the current folder or on the MATLAB path, but
exists in:
    /Users/msc/Desktop/optimal_control/orbit_transfer/GTO_tulip/direct
    /Users/msc/Desktop/optimal_control/orbit_transfer/earth_elliptic_to_geo_CR3BP/direct

Change the MATLAB current folder or add its folder to the MATLAB path.

Error in test_params (line 12)
addpath(fullfile(here_, '..'));  setup_paths;
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
}
```

### GREEN: Test Passes (Step 4)
```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); addpath('tests'); test_params"
```
**Output (success):**
```
test_params PASS (min-throttle dry T/W = 1.346)
```

### Calculation Verification
- P.Tmin = 0.40 × 845 kN = 338 kN
- twMin = 338 kN / (25,600 kg × 9.80665 m/s²) = 1.346 ✓
- 1.346 > 1 confirms hoverslam is forced (the booster will drop even at min throttle)

## Git Commit

**Commit SHA:** `312b9a8`
**Message:**
```
booster_landing: skeleton, params, sanity test

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
```
**Files staged:**
- `booster_landing/setup_paths.m`
- `booster_landing/lib/booster_params.m`
- `booster_landing/tests/test_params.m`
- `booster_landing/.gitignore`
- `booster_landing/results/.gitkeep`

## Self-Review Findings

### Strengths
1. **Spec compliance:** All 23 parameter fields match the brief verbatim; no deviations
2. **Function headers:** Both functions carry INPUTS/OUTPUTS/REFERENCES blocks per CLAUDE.md standard
3. **No loop-variable violations:** Code uses no `i` or `j` as loop variables (no loops present)
4. **TDD discipline:** RED → GREEN → COMMIT sequence followed cleanly
5. **Module linkage:** setup_paths() cleanly chains the four layer directories (lib, certify, viz, tests) and CasADi

### Test Assertions
All three assertions in test_params pass:
1. m0 = 30,000 kg > mdry = 25,600 kg ✓
2. twMin = 1.346 > 1 (hoverslam forced) ✓
3. Bounds: Tmin=338kN < Tmax=845kN; tf_lo=10s < tf_hi=50s; gs_deg=30° ∈ (0,90) ✓
4. Phase 1: drag.on = false ✓

### Minor Notes
- CasADi 3.7.0 path assumed at `$HOME/casadi-3.7.0`; will fail softly if not present (not used until Task 3)
- No issues with directory isolation or git hygiene; only target files staged, no spurious adds

## Test Coverage

The test validates the core campaign invariants:
- **Propellant margin:** m0 > mdry ensures a 4,400 kg fuel payload exists
- **Hoverslam law:** twMin = 1.346 forces a bang-bang terminal descent (no hover phase)
- **Solver bracket:** tf_lo=10s to tf_hi=50s provides 40-second search window for optimal free final time
- **Glideslope:** gs_deg=30° defines minimum-elevation landing constraints (typical Falcon 9)
- **Vacuum default:** Phase 1 solves without atmosphere drag for exact convex geometry

## Next Steps (Task 2)
- Implement `pdg_dynamics.m` — 3-DOF point-mass dynamics (x, y, z, vx, vy, vz) with gravity and thrust vector
- Implement complex-step Jacobian test to validate numerical derivatives

---
**Report Generated:** 2026-08-08  
**Campaign Root:** `/Users/msc/Desktop/optimal_control/booster_landing/`  
**MATLAB Version:** R2025b  
**Exit Status:** All GREEN
