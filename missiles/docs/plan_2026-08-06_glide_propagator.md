# HGV Glide Propagator Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A validated 3-DOF hypersonic glide propagator that runs a prescribed-control trajectory from entry interface to impact and agrees with closed-form solutions.

**Architecture:** One shared `+coorbital` package under `missiles/`. Physics models (atmosphere, gravity, aerodynamics) are injected into the equations of motion as function handles carried in an `env` struct, so raising fidelity later is a one-line swap in an entry script rather than an edit to the EOM. A phase driver chains integrations with ODE events at the handoffs. Vehicle folders hold only parameter files and entry scripts.

**Tech Stack:** MATLAB R2025b (`/Applications/MATLAB_R2025b.app/bin/matlab`), `ode45`, pumpkyn (`~/Desktop/proj7/external/pumpkyn`) for frames, time, vector math, constants, and plotting.

## Global Constraints

- Style: pumpkyn house style per `~/Desktop/proj7/doc/pumpkyn_style_guide.md` — `%%` header quartet (`Purpose`/`Inputs`/`Outputs`/`Revision History`) with aligned three-column input and output blocks stating units, closed by the `%% ------------------------ Begin Code Sequence ---------------------------` divider; `=` signs column-aligned with right-justified names; colon-terminated `%%` step comments; `if nargin == 0` self-demo in every library function calling itself by full namespace.
- Add a `%% References:` block whenever a function's math has a source.
- Never use `i` or `j` as loop or index variables — they are the imaginary unit. Use `k`, `ii`, `idx`, or a meaningful name.
- Norms via `pumpkyn.util.vmag` or `sqrt(sum(x.^2,dim))`. **Never `norm`** — it does not vectorize and is not complex-step safe.
- Never modify `~/Desktop/proj7/external/pumpkyn`. Never consume `getConst`'s `deg2ArcSec` or `c` fields (both known bad).
- SI units throughout the library: metres, seconds, kilograms, radians, kelvin. Entry scripts may accept degrees and kilometres in the user block and convert immediately.
- All angles in the state vector are radians.
- Every library function lives in `missiles/+coorbital/+<subpkg>/`, one function per file, file name equal to function name, called by full namespace.
- Author line in every header: `Michael Casey` with the current date, then `Copyright 2026 Coorbital, Inc.`
- **On pumpkyn reuse:** the spec's §3 lists frame conversions (`LLA2ECEF`, `ECItoECEF`, …) as reusable. Phase 1 does not call them: the glide state carries geocentric latitude/longitude directly and never leaves that frame, so there is nothing to convert. They become load-bearing at geodetic altitude and ECI trajectory export, both out of scope here. `pumpkyn.util.vmag` is likewise unused because no 3-vector norms appear in the 6-state EOM. Do not add pumpkyn calls that serve no purpose just to satisfy §3.
- **No lint-suppression pragmas.** Never write `%#ok<NASGU>`, `%#ok<AGROW>`, `%#ok<*ANY>`, or any other `%#ok` directive — humans do not write them and they hide real problems. If a variable is unused, delete it or capture it with `~`. If a loop grows an array, preallocate it.
- Run MATLAB headlessly as: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/missiles'); <command>" 2>&1 | grep -vE "Home License|personal use|academic, research|organizational use"`

---

## File Structure

| File | Responsibility |
|---|---|
| `+coorbital/+util/missileConst.m` | Earth and air constants in one struct; the single source of truth |
| `+coorbital/+atmos/expAtmos.m` | Isothermal exponential atmosphere: density, pressure, temperature, sound speed |
| `+coorbital/+grav/sphereGrav.m` | Point-mass gravity, returns radial and latitudinal components |
| `+coorbital/+aero/constLD.m` | Constant-CL, constant-L/D aerodynamics |
| `+coorbital/+util/vehicleDefaults.m` | Vehicle parameter struct with validation |
| `+coorbital/+eom/glide3DOF.m` | The 6-state glide equations of motion |
| `+coorbital/+prop/eventAltitude.m` | ODE event: altitude crossing |
| `+coorbital/+prop/phaseRun.m` | Multi-phase integration driver |
| `+coorbital/+guide/prescribed.m` | Evaluate a control schedule at time `t` |
| `HGV/vehicle_hgv.m` | HGV parameter values |
| `HGV/run_glide.m` | Entry script with the user-parameter block |
| `tests/run_tests.m` | Dependency-free test runner |
| `tests/test_*.m` | One file per unit under test |
| `docs/LESSONS_LEARNED.md` | Running log, started in Task 1 |
| `docs/README.md` | Code organization and how to run |

---

### Task 1: Test harness and constants

**Files:**
- Create: `missiles/tests/run_tests.m`
- Create: `missiles/tests/test_missileConst.m`
- Create: `missiles/+coorbital/+util/missileConst.m`
- Create: `missiles/docs/LESSONS_LEARNED.md`

**Interfaces:**
- Consumes: nothing
- Produces: `c = coorbital.util.missileConst()` returning a struct with fields `muE` (m³/s²), `rE` (m), `omegaE` (rad/s), `rho0` (kg/m³), `Hscale` (m), `T0` (K), `Rair` (J/kg/K), `gamAir` (–), `g0` (m/s²)

- [ ] **Step 1: Write the test runner**

```matlab
function run_tests()
%% Purpose:
%
%  Run every test_*.m in this directory and report pass/fail. Each test file
%  is a script or function that throws on failure and returns silently on
%  success.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none                                         Prints a summary; sets a
%                                               non-zero exit code on failure
%                                               when run under -batch
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Locate this directory and the library root:
              here = fileparts(mfilename('fullpath'));
              root = fullfile(here,'..');
    addpath(root);

%% Vehicle folders are added only once they exist, so the suite stays quiet
%% while the library is still being built out:
           vehDirs = {'HGV','BM','ICBM'};
    for k = 1:numel(vehDirs)
              vDir = fullfile(root,vehDirs{k});
        if isfolder(vDir)
            addpath(vDir);
        end
    end

%% Discover the test files:
             files = dir(fullfile(here,'test_*.m'));
                nT = numel(files);
            passed = 0;
            failed = 0;
          failures = cell(nT,1);

%% Run each test in turn, recording the message of any that throws:
    for k = 1:nT
             name = erase(files(k).name,'.m');
        try
            feval(name);
            passed = passed + 1;
            fprintf('  PASS  %s\n',name);
        catch err
            failed = failed + 1;
            failures{k} = sprintf('%s: %s',name,err.message);
            fprintf('  FAIL  %s\n',name);
        end
    end
          failures = failures(~cellfun(@isempty,failures));

%% Report:
    fprintf('\n%d passed, %d failed\n',passed,failed);
    for k = 1:numel(failures)
        fprintf('\n----- %s\n',failures{k});
    end
    if failed > 0
        error('missiles:testsFailed','%d test(s) failed',failed);
    end
end
```

- [ ] **Step 2: Write the failing test**

```matlab
function test_missileConst()
%% Purpose:
%
%  Verify the constants struct exposes every field the library depends on,
%  with values in SI and of the right order of magnitude.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none                                         Throws on any failed assertion
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();

%% Every required field is present:
          required = {'muE','rE','omegaE','rho0','Hscale','T0','Rair','gamAir','g0'};
    for k = 1:numel(required)
        assert(isfield(c,required{k}), ...
            'missileConst is missing field %s',required{k});
    end

%% Values are physically sensible:
    assert(abs(c.muE  - 3.986004418e14) < 1e6,   'muE out of range');
    assert(abs(c.rE   - 6378137)        < 1,     'rE out of range');
    assert(abs(c.g0   - 9.80665)        < 1e-3,  'g0 out of range');
    assert(abs(c.rho0 - 1.225)          < 1e-3,  'rho0 out of range');
    assert(c.Hscale > 6000 && c.Hscale < 9000,   'Hscale out of range');
    assert(abs(c.gamAir - 1.4)          < 1e-6,  'gamAir out of range');
    assert(abs(c.T0   - 250)            < 1e-9,  'T0 out of range');
    assert(abs(c.Rair - 287.053)        < 1e-3,  'Rair out of range');

%% Surface gravity from muE and rE agrees with g0 to better than 0.5 percent:
              gSurf = c.muE/c.rE^2;
    assert(abs(gSurf - c.g0)/c.g0 < 5e-3, ...
        'mu/r^2 = %.4f disagrees with g0 = %.4f',gSurf,c.g0);
end
```

- [ ] **Step 3: Run it and confirm it fails**

Run:
```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/missiles'); run('tests/run_tests')" 2>&1 | grep -vE "Home License|personal use|academic, research|organizational use"
```
Expected: `FAIL test_missileConst` with an "Unrecognized function" message for `coorbital.util.missileConst`.

- [ ] **Step 4: Implement the constants**

```matlab
function c = missileConst()
%% Purpose:
%
%  Return the Earth and standard-air constants used across the missile
%  trajectory library. Single source of truth -- no routine may hard-code
%  these values.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  c                Struct                      Constants, SI units:
%                                               muE    (m^3/s^2)
%                                               rE     (m) equatorial radius
%                                               omegaE (rad/s) Earth rotation
%                                               rho0   (kg/m^3) sea-level density
%                                               Hscale (m) density scale height
%                                               T0     (K) isothermal temperature
%                                               Rair   (J/kg/K) specific gas constant
%                                               gamAir (-) ratio of specific heats
%                                               g0     (m/s^2) standard gravity
%
%% References:
%   [1] WGS-84, NIMA TR8350.2, 3rd ed., 1997. (muE, rE, omegaE)
%   [2] U.S. Standard Atmosphere, 1976. (rho0, Rair, gamAir, g0)
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Earth:
             c.muE = 3.986004418e14;    %m^3/s^2, WGS-84
              c.rE = 6378137;           %m, WGS-84 equatorial radius
          c.omegaE = 7.292115e-5;       %rad/s, WGS-84

%% Standard air:
            c.rho0 = 1.225;             %kg/m^3, sea level, US76
          c.Hscale = 7200;              %m, exponential density scale height
              c.T0 = 250;               %K, isothermal reference (see expAtmos)
            c.Rair = 287.0528;          %J/kg/K, specific gas constant for air
          c.gamAir = 1.4;               %ratio of specific heats
              c.g0 = 9.80665;           %m/s^2, standard gravity
end
```

- [ ] **Step 5: Run the tests and confirm they pass**

Run the same command as Step 3. Expected: `PASS test_missileConst`, `1 passed, 0 failed`.

- [ ] **Step 6: Start the lessons log**

Create `docs/LESSONS_LEARNED.md`:

```markdown
# Lessons Learned — Missile Trajectory Library

Running log. Newest entries at the top. Record what broke, what fixed it, and
what a future reader would otherwise rediscover the hard way.

## 2026-08-06 — Constants live in one place

`missileConst` is the single source of truth for Earth and air constants.
pumpkyn's `getConst` was not used for these: it carries `g` and `R` but no
`muE`, no Earth radius, and two fields known to be wrong (`deg2ArcSec`, `c`).
Mixing the two would make it ambiguous which constant a routine actually used.
```

- [ ] **Step 7: Commit**

```bash
cd ~/Desktop/optimal_control
git add missiles/tests missiles/+coorbital/+util/missileConst.m missiles/docs/LESSONS_LEARNED.md
git commit -m "missiles: test harness and constants"
```

---

### Task 2: Exponential atmosphere

**Files:**
- Create: `missiles/+coorbital/+atmos/expAtmos.m`
- Create: `missiles/tests/test_expAtmos.m`

**Interfaces:**
- Consumes: `coorbital.util.missileConst()`
- Produces: `[rho,P,T,a] = coorbital.atmos.expAtmos(h)` where `h` is geometric altitude in metres, accepts `[N x 1]`, returns `[N x 1]` each in kg/m³, Pa, K, m/s

- [ ] **Step 1: Write the failing test**

```matlab
function test_expAtmos()
%% Purpose:
%
%  Verify the exponential atmosphere reproduces sea-level density, decays by
%  exactly one e-fold per scale height, vectorizes, and returns pressure and
%  sound speed matching values computed independently of the implementation.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none                                         Throws on any failed assertion
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();

%% Sea level returns rho0 exactly:
    [rho,P,T,a] = coorbital.atmos.expAtmos(0);
    assert(abs(rho - c.rho0) < 1e-12,'sea-level density wrong');

%% One scale height is one e-fold:
       [rhoH,~,~,~] = coorbital.atmos.expAtmos(c.Hscale);
    assert(abs(rhoH - c.rho0/exp(1)) < 1e-12,'scale-height decay wrong');

%% Sea-level pressure and sound speed against values worked out by hand, so a
%% formula error cannot hide behind the implementation's own arithmetic:
    assert(abs(T - 250)      < 1e-9,  'isothermal temperature wrong');
    assert(abs(P - 87909.92) < 1e-2,  'sea-level pressure wrong: got %.2f Pa',P);
    assert(abs(a - 316.9677) < 1e-3,  'sea-level sound speed wrong: got %.4f m/s',a);

%% The ideal gas law then holds at an altitude the hand values do not cover:
    [rho2,P2,T2,~] = coorbital.atmos.expAtmos(40000);
    assert(abs(P2 - rho2*c.Rair*T2) < 1e-6,'P does not satisfy the ideal gas law');

%% Vectorizes over a column of altitudes:
              hVec = (0:10000:100000)';
    [rhoV,PV,TV,aV] = coorbital.atmos.expAtmos(hVec);
    assert(isequal(size(rhoV),size(hVec)),'rho not [N x 1]');
    assert(isequal(size(PV),  size(hVec)),'P not [N x 1]');
    assert(isequal(size(TV),  size(hVec)),'T not [N x 1]');
    assert(isequal(size(aV),  size(hVec)),'a not [N x 1]');
    assert(all(diff(rhoV) < 0),'density must decrease monotonically');

%% Negative altitude does not produce a non-finite value:
       [rhoNeg,~,~,~] = coorbital.atmos.expAtmos(-100);
    assert(isfinite(rhoNeg),'negative altitude gave a non-finite density');
end
```

- [ ] **Step 2: Run it and confirm it fails**

Run the Task 1 Step 3 command. Expected: `FAIL test_expAtmos`, unrecognized function.

- [ ] **Step 3: Implement**

```matlab
function [rho,P,T,a] = expAtmos(h)
%% Purpose:
%
%  Isothermal exponential atmosphere. Density decays with a single scale
%  height and temperature is held constant, which makes the model crude above
%  roughly 30 km but gives closed-form entry solutions to validate against.
%  Interface-compatible with richer models (us76Atmos, nrlmsise) so the
%  environment struct can swap them without touching the equations of motion.
%
%% Inputs:
%
%  h                [N x 1]                     Geometric altitude above the
%                                               reference sphere (m)
%
%% Outputs:
%
%  rho              [N x 1]                     Density (kg/m^3)
%
%  P                [N x 1]                     Pressure (Pa)
%
%  T                [N x 1]                     Temperature (K)
%
%  a                [N x 1]                     Speed of sound (m/s)
%
%% References:
%   [1] Vinh, N.X., Busemann, A., Culp, R.D., "Hypersonic and Planetary Entry
%       Flight Mechanics," Univ. Michigan Press, 1980, Ch. 2.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
              hVec = linspace(0,120e3,400)';
    [rhoV,~,~,~] = coorbital.atmos.expAtmos(hVec);
    figure('color',[1 1 1]);
    semilogx(rhoV,hVec./1000,'linewidth',1.5); grid on;
    xlabel('density (kg/m^3)'); ylabel('altitude (km)');
    title('Isothermal exponential atmosphere');
    [rho,P,T,a] = deal([]);
    return;
end

                 c = coorbital.util.missileConst();

%% Density decays one e-fold per scale height:
               rho = c.rho0.*exp(-h./c.Hscale);

%% Isothermal by construction:
                 T = c.T0.*ones(size(h));

%% Pressure from the ideal gas law, sound speed from the same temperature:
                 P = rho.*c.Rair.*T;
                 a = sqrt(c.gamAir.*c.Rair.*T);
end
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run the Task 1 Step 3 command. Expected: `2 passed, 0 failed`.

- [ ] **Step 5: Verify the self-demo runs**

Run:
```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/missiles'); coorbital.atmos.expAtmos()" 2>&1 | grep -vE "Home License|personal use|academic, research|organizational use"
```
Expected: no output and no error. The figure is created off-screen under `-batch`.

- [ ] **Step 6: Commit**

```bash
cd ~/Desktop/optimal_control
git add missiles/+coorbital/+atmos missiles/tests/test_expAtmos.m
git commit -m "missiles: isothermal exponential atmosphere"
```

---

### Task 3: Spherical gravity

**Files:**
- Create: `missiles/+coorbital/+grav/sphereGrav.m`
- Create: `missiles/tests/test_sphereGrav.m`

**Interfaces:**
- Consumes: `coorbital.util.missileConst()`
- Produces: `[gr,gLat] = coorbital.grav.sphereGrav(r,lat)` where `r` is geocentric radius in metres `[N x 1]` and `lat` is geocentric latitude in radians `[N x 1]`; `gr` is the downward (toward the centre) magnitude in m/s² and `gLat` is the northward component, zero for this model

- [ ] **Step 1: Write the failing test**

```matlab
function test_sphereGrav()
%% Purpose:
%
%  Verify point-mass gravity: correct surface value, inverse-square falloff,
%  no latitude dependence, and a zero latitudinal component.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();

%% Surface gravity is mu/r^2 and points down (positive by convention):
       [gr,gLat] = coorbital.grav.sphereGrav(c.rE,0);
    assert(abs(gr - 9.7983) < 1e-4,'surface gravity wrong: got %.6f m/s^2',gr);
    assert(gr > 0,'gr must be positive downward');
    assert(gLat == 0,'spherical gravity has no latitudinal component');

%% Inverse square: doubling the radius quarters the acceleration:
      [gr2,~] = coorbital.grav.sphereGrav(2*c.rE,0);
    assert(abs(gr2 - gr/4) < 1e-12,'not inverse-square');

%% No latitude dependence for a point mass:
      [grEq,~] = coorbital.grav.sphereGrav(c.rE,0);
      [grPo,~] = coorbital.grav.sphereGrav(c.rE,pi/2);
    assert(abs(grEq - grPo) < 1e-12,'point-mass gravity must not vary with latitude');

%% Vectorizes:
              rVec = c.rE + (0:10000:100000)';
      [grV,gLatV] = coorbital.grav.sphereGrav(rVec,zeros(size(rVec)));
    assert(isequal(size(grV),size(rVec)),'gr not [N x 1]');
    assert(isequal(size(gLatV),size(rVec)),'gLat not [N x 1]');
    assert(all(diff(grV) < 0),'gravity must decrease with radius');
end
```

- [ ] **Step 2: Run it and confirm it fails**

Run the Task 1 Step 3 command. Expected: `FAIL test_sphereGrav`.

- [ ] **Step 3: Implement**

```matlab
function [gr,gLat] = sphereGrav(r,lat)
%% Purpose:
%
%  Point-mass gravitational acceleration. Returns two components so that
%  oblate models (j2Grav, j24Grav) satisfy the same interface -- J2 produces a
%  latitudinal acceleration that a single "local up" scalar cannot represent.
%
%% Inputs:
%
%  r                [N x 1]                     Geocentric radius (m)
%
%  lat              [N x 1]                     Geocentric latitude (rad).
%                                               Unused here; present so the
%                                               signature matches j2Grav.
%
%% Outputs:
%
%  gr               [N x 1]                     Acceleration toward the centre,
%                                               positive downward (m/s^2)
%
%  gLat             [N x 1]                     Northward acceleration (m/s^2);
%                                               identically zero for a point mass
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
                 c = coorbital.util.missileConst();
              rVec = c.rE + linspace(0,200e3,200)';
       [grV,~] = coorbital.grav.sphereGrav(rVec,zeros(size(rVec)));
    figure('color',[1 1 1]);
    plot((rVec-c.rE)./1000,grV,'linewidth',1.5); grid on;
    xlabel('altitude (km)'); ylabel('g (m/s^2)');
    title('Point-mass gravity');
    [gr,gLat] = deal([]);
    return;
end

                 c = coorbital.util.missileConst();

%% Inverse-square magnitude, positive toward the centre:
                gr = c.muE./r.^2;

%% A point mass has no latitudinal component:
              gLat = zeros(size(r));
end
```

- [ ] **Step 4: Run the tests and confirm they pass**

Expected: `3 passed, 0 failed`.

- [ ] **Step 5: Commit**

```bash
cd ~/Desktop/optimal_control
git add missiles/+coorbital/+grav missiles/tests/test_sphereGrav.m
git commit -m "missiles: point-mass gravity with J2-compatible interface"
```

---

### Task 4: Aerodynamics and vehicle parameters

**Files:**
- Create: `missiles/+coorbital/+aero/constLD.m`
- Create: `missiles/+coorbital/+util/vehicleDefaults.m`
- Create: `missiles/tests/test_constLD.m`

**Interfaces:**
- Consumes: nothing
- Produces:
  - `veh = coorbital.util.vehicleDefaults()` returning a struct with `mass` (kg), `Sref` (m²), `CL` (–), `LD` (–), `noseRadius` (m)
  - `[CL,CD] = coorbital.aero.constLD(alpha,mach,veh)` — `alpha` and `mach` are accepted and ignored by this model, present so table-driven models satisfy the same signature

- [ ] **Step 1: Write the failing test**

```matlab
function test_constLD()
%% Purpose:
%
%  Verify the constant-L/D aerodynamic model returns the vehicle's CL, a CD
%  consistent with the requested lift-to-drag ratio, and ignores alpha and
%  Mach as documented.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               veh = coorbital.util.vehicleDefaults();

%% Required vehicle fields exist and are positive:
          required = {'mass','Sref','CL','LD','noseRadius'};
    for k = 1:numel(required)
        assert(isfield(veh,required{k}),'vehicleDefaults missing %s',required{k});
        assert(veh.(required{k}) > 0,'%s must be positive',required{k});
    end

%% CL is the vehicle CL and CD gives exactly the requested L/D:
          [CL,CD] = coorbital.aero.constLD(0.1,20,veh);
    assert(abs(CL - veh.CL) < 1e-12,'CL should equal veh.CL');
    assert(abs(CL/CD - veh.LD) < 1e-12,'CL/CD should equal veh.LD');

%% alpha and Mach are ignored by this model:
        [CL2,CD2] = coorbital.aero.constLD(0.4,6,veh);
    assert(CL2 == CL && CD2 == CD,'constLD must not depend on alpha or Mach');
end
```

- [ ] **Step 2: Run it and confirm it fails**

Expected: `FAIL test_constLD`.

- [ ] **Step 3: Implement the vehicle defaults**

```matlab
function veh = vehicleDefaults()
%% Purpose:
%
%  Default hypersonic glide vehicle parameters. PLACEHOLDER VALUES drawn from
%  the open literature for a generic lifting entry body -- they are of the
%  right order but are not any specific vehicle. Replace before any result
%  leaves this machine.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  veh              Struct                      Vehicle parameters:
%                                               mass       (kg)
%                                               Sref       (m^2) reference area
%                                               CL         (-) lift coefficient
%                                               LD         (-) lift-to-drag ratio
%                                               noseRadius (m) for heating
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

          veh.mass = 900;              %kg, PLACEHOLDER
          veh.Sref = 0.75;             %m^2, PLACEHOLDER
            veh.CL = 0.35;             %PLACEHOLDER
            veh.LD = 2.5;              %PLACEHOLDER, typical waverider range 2-4
    veh.noseRadius = 0.05;             %m, PLACEHOLDER
end
```

- [ ] **Step 4: Implement the aerodynamics**

```matlab
function [CL,CD] = constLD(alpha,mach,veh)
%% Purpose:
%
%  Constant lift coefficient and constant lift-to-drag ratio. The crudest
%  usable aerodynamic model, and the one for which equilibrium glide has a
%  closed-form solution. Signature matches table-driven models so the
%  environment struct can swap them.
%
%% Inputs:
%
%  alpha            [N x 1]                     Angle of attack (rad). Accepted
%                                               and ignored by this model.
%
%  mach             [N x 1]                     Mach number. Accepted and
%                                               ignored by this model.
%
%  veh              Struct                      Vehicle parameters; uses the
%                                               CL and LD fields
%
%% Outputs:
%
%  CL               scalar                      Lift coefficient (-)
%
%  CD               scalar                      Drag coefficient (-)
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
               veh = coorbital.util.vehicleDefaults();
       [CLd,CDd] = coorbital.aero.constLD(0,20,veh);
    fprintf('CL = %.4f, CD = %.4f, L/D = %.2f\n',CLd,CDd,CLd/CDd);
    [CL,CD] = deal([]);
    return;
end

%% Constant by construction; alpha and mach are deliberately unused:
                CL = veh.CL;
                CD = veh.CL/veh.LD;
end
```

- [ ] **Step 5: Run the tests and confirm they pass**

Expected: `4 passed, 0 failed`.

- [ ] **Step 6: Commit**

```bash
cd ~/Desktop/optimal_control
git add missiles/+coorbital/+aero missiles/+coorbital/+util/vehicleDefaults.m missiles/tests/test_constLD.m
git commit -m "missiles: constant-L/D aero and placeholder vehicle parameters"
```

---

### Task 5: Glide equations of motion

**Files:**
- Create: `missiles/+coorbital/+eom/glide3DOF.m`
- Create: `missiles/tests/test_glide3DOF.m`

**Interfaces:**
- Consumes: `coorbital.util.missileConst`, and through `env` the handles `env.atmos`, `env.grav`, `env.aero`, plus `env.omegaE`
- Produces: `xdot = coorbital.eom.glide3DOF(t,x,u,veh,env)` where `x = [r; lon; lat; V; gamma; psi]` is `[6 x 1]` in SI and radians, `u = [alpha; sigma]` is `[2 x 1]` in radians, and `xdot` is `[6 x 1]`

**Note on the two exact tests below:** they are the highest-value checks in the plan. The spec's §5 lists an "energy bookkeeping" check comparing drag work against mechanical energy lost; this plan uses the stronger instantaneous form instead — asserting `dV/dt` equals `-D/m - g sin(gamma)` computed independently — which catches the same class of error at a single point without integration tolerance in the way. Vacuum energy conservation catches sign errors and frame mistakes that the approximate analytic comparisons in Task 8 would hide.

- [ ] **Step 1: Write the failing test**

```matlab
function test_glide3DOF()
%% Purpose:
%
%  Two exact checks on the glide equations of motion. In vacuum with no lift
%  the specific orbital energy must be conserved to solver tolerance, and with
%  drag present the mechanical energy must decrease at exactly the rate the
%  drag force does work.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();
               veh = coorbital.util.vehicleDefaults();

%% Vacuum environment: zero density, point-mass gravity, no lift:
        envVac.atmos = @(h) deal(zeros(size(h)),zeros(size(h)), ...
                                 c.T0.*ones(size(h)),ones(size(h)));
         envVac.grav = @coorbital.grav.sphereGrav;
         envVac.aero = @(al,ma,vv) deal(0,0);
       envVac.omegaE = 0;

%% A lofted ballistic arc:
                x0 = [c.rE + 100e3; 0; 0; 5000; deg2rad(20); deg2rad(90)];
              odeF = @(t,x) coorbital.eom.glide3DOF(t,x,[0;0],veh,envVac);
              opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
           [~,X] = ode45(odeF,[0 300],x0,opts);

%% Specific energy V^2/2 - mu/r is conserved in vacuum:
              rArc = X(:,1);
              vArc = X(:,4);
              eArc = 0.5.*vArc.^2 - c.muE./rArc;
            eDrift = abs(eArc - eArc(1))./abs(eArc(1));
    assert(max(eDrift) < 1e-8, ...
        'vacuum energy drifted by %.3e relative; expected < 1e-8',max(eDrift));

%% With drag, energy must fall at exactly the rate drag does work:
          env.atmos = @coorbital.atmos.expAtmos;
           env.grav = @coorbital.grav.sphereGrav;
           env.aero = @coorbital.aero.constLD;
         env.omegaE = 0;

                xD = [c.rE + 40e3; 0; 0; 3000; deg2rad(-2); deg2rad(90)];
             xdotD = coorbital.eom.glide3DOF(0,xD,[0;0],veh,env);

%% Drag power per unit mass, computed independently of the EOM:
        [rhoD,~,~,~] = coorbital.atmos.expAtmos(xD(1) - c.rE);
          [~,CDd] = coorbital.aero.constLD(0,0,veh);
             dragA = 0.5*rhoD*xD(4)^2*veh.Sref*CDd/veh.mass;
          [grD,~] = coorbital.grav.sphereGrav(xD(1),xD(3));

%% dV/dt must equal -dragAccel - g*sin(gamma) exactly:
            vdotEx = -dragA - grD*sin(xD(5));
    assert(abs(xdotD(4) - vdotEx) < 1e-12, ...
        'dV/dt = %.9f does not match -D/m - g sin(gamma) = %.9f',xdotD(4),vdotEx);

%% Kinematics: with gamma = 0 the radius must be stationary:
                xF = [c.rE + 40e3; 0; 0; 3000; 0; deg2rad(90)];
             xdotF = coorbital.eom.glide3DOF(0,xF,[0;0],veh,env);
    assert(abs(xdotF(1)) < 1e-12,'dr/dt must vanish at gamma = 0');

%% Heading due east at the equator moves longitude, not latitude:
    assert(xdotF(2) > 0,'eastward flight must increase longitude');
    assert(abs(xdotF(3)) < 1e-12,'eastward flight at the equator must not change latitude');
end
```

- [ ] **Step 2: Run it and confirm it fails**

Expected: `FAIL test_glide3DOF`.

- [ ] **Step 3: Implement**

```matlab
function xdot = glide3DOF(t,x,u,veh,env)
%% Purpose:
%
%  Three-degree-of-freedom point-mass equations of motion for atmospheric
%  glide over a rotating spherical Earth. The rotation terms are gated on
%  env.omegaE, so setting it to zero recovers the non-rotating case exactly;
%  the latitudinal gravity component gLat is carried through so that an
%  oblate gravity model drops in without editing this file.
%
%% Inputs:
%
%  t                scalar                      Time since phase start (s).
%                                               Unused; present for ode45.
%
%  x                [6 x 1]                     State:
%                                               r     (m) geocentric radius
%                                               lon   (rad) longitude
%                                               lat   (rad) geocentric latitude
%                                               V     (m/s) planet-relative speed
%                                               gamma (rad) flight path angle
%                                               psi   (rad) heading, clockwise
%                                                     from north
%
%  u                [2 x 1]                     Control:
%                                               alpha (rad) angle of attack
%                                               sigma (rad) bank angle
%
%  veh              Struct                      Vehicle parameters from
%                                               coorbital.util.vehicleDefaults
%
%  env              Struct                      Environment model handles:
%                                               atmos, grav, aero, omegaE
%
%% Outputs:
%
%  xdot             [6 x 1]                     State derivative
%
%% References:
%   [1] Vinh, N.X., Busemann, A., Culp, R.D., "Hypersonic and Planetary Entry
%       Flight Mechanics," Univ. Michigan Press, 1980, Eqs. (2.28)-(2.33).
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
                 c = coorbital.util.missileConst();
               veh = coorbital.util.vehicleDefaults();
         env.atmos = @coorbital.atmos.expAtmos;
          env.grav = @coorbital.grav.sphereGrav;
          env.aero = @coorbital.aero.constLD;
        env.omegaE = 0;
                x0 = [c.rE + 60e3; 0; 0; 6000; deg2rad(-1); deg2rad(90)];
             xdot = coorbital.eom.glide3DOF(0,x0,[0;0],veh,env);
    fprintf('xdot = [%.4g %.4g %.4g %.4g %.4g %.4g]''\n',xdot);
    return;
end

%% Unpack the state and control:
                 r = x(1);
               lat = x(3);
                 V = x(4);
             gamma = x(5);
               psi = x(6);
             alpha = u(1);
             sigma = u(2);

                 c = coorbital.util.missileConst();
                 h = r - c.rE;

%% Guard the coordinate singularities so they fail loudly, not as silent NaN:
if abs(cos(lat)) < 1e-8
    error('coorbital:glide3DOF:polarSingularity', ...
        'Longitude rate is singular at the pole (lat = %.6f rad).',lat);
end
if abs(cos(gamma)) < 1e-8
    error('coorbital:glide3DOF:verticalFlight', ...
        'Heading rate is singular in vertical flight (gamma = %.6f rad).',gamma);
end
if V < 1
    error('coorbital:glide3DOF:zeroSpeed', ...
        'Equations are singular as V approaches zero (V = %.3f m/s).',V);
end

%% Environment: density, gravity, aerodynamic coefficients:
    [rho,~,~,aSnd] = env.atmos(h);
        [gr,gLat] = env.grav(r,lat);
          [CL,CD] = env.aero(alpha,V./aSnd,veh);

%% Aerodynamic accelerations:
              qbar = 0.5.*rho.*V.^2;
              aLift = qbar.*veh.Sref.*CL./veh.mass;
              aDrag = qbar.*veh.Sref.*CD./veh.mass;

%% Kinematics:
              rdot = V.*sin(gamma);
            londot = V.*cos(gamma).*sin(psi)./(r.*cos(lat));
            latdot = V.*cos(gamma).*cos(psi)./r;

%% Dynamics over a non-rotating sphere:
              Vdot = -aDrag - gr.*sin(gamma) ...
                     + gLat.*cos(gamma).*cos(psi);
          gammadot = aLift.*cos(sigma)./V - (gr./V - V./r).*cos(gamma) ...
                     - gLat.*sin(gamma).*cos(psi)./V;
            psidot = aLift.*sin(sigma)./(V.*cos(gamma)) ...
                     + V.*cos(gamma).*sin(psi).*tan(lat)./r ...
                     - gLat.*sin(psi)./(V.*cos(gamma));

%% Rotating-Earth Coriolis and centrifugal terms, gated on omegaE:
                om = env.omegaE;
if om ~= 0
              Vdot = Vdot + om.^2.*r.*cos(lat).* ...
                     (sin(gamma).*cos(lat) - cos(gamma).*sin(lat).*cos(psi));
          gammadot = gammadot + 2.*om.*cos(lat).*sin(psi) ...
                     + om.^2.*r.*cos(lat).* ...
                       (cos(gamma).*cos(lat) + sin(gamma).*cos(psi).*sin(lat))./V;
            psidot = psidot + 2.*om.*(sin(lat) - cos(lat).*cos(psi).*tan(gamma)) ...
                     + om.^2.*r.*sin(lat).*cos(lat).*sin(psi)./(V.*cos(gamma));
end

%% Assemble:
              xdot = [rdot; londot; latdot; Vdot; gammadot; psidot];
end
```

- [ ] **Step 4: Run the tests and confirm they pass**

Expected: `5 passed, 0 failed`. If the vacuum energy check fails, the fault is almost always a sign on `gr` in `Vdot` or a missing `V/r` in `gammadot` — check those two before anything else.

- [ ] **Step 5: Record the lesson**

Prepend to `docs/LESSONS_LEARNED.md`:

```markdown
## 2026-08-06 — Vacuum energy conservation is the test that matters

The glide EOM has six equations and many chances for a sign error. Propagating
a vacuum ballistic arc and checking that V^2/2 - mu/r holds constant to 1e-8
catches nearly all of them in one assertion, because any sign flip in gravity
or the centrifugal term breaks conservation immediately. Approximate checks
against analytic glide solutions do not — they tolerate a percent of error by
construction, which is enough to hide a wrong term.
```

- [ ] **Step 6: Commit**

```bash
cd ~/Desktop/optimal_control
git add missiles/+coorbital/+eom missiles/tests/test_glide3DOF.m missiles/docs/LESSONS_LEARNED.md
git commit -m "missiles: 3-DOF glide equations of motion"
```

---

### Task 6: Phase driver and events

**Files:**
- Create: `missiles/+coorbital/+prop/eventAltitude.m`
- Create: `missiles/+coorbital/+guide/prescribed.m`
- Create: `missiles/+coorbital/+prop/phaseRun.m`
- Create: `missiles/tests/test_phaseRun.m`

**Interfaces:**
- Consumes: `coorbital.eom.glide3DOF`, `coorbital.util.missileConst`
- Produces:
  - `[value,isterminal,direction] = coorbital.prop.eventAltitude(t,x,hStop)` — standard ODE event signature, terminal, decreasing
  - `u = coorbital.guide.prescribed(t,x,sched)` returning `[2 x 1]`, where `sched` has fields `tGrid` `[1 x K]`, `alpha` `[1 x K]`, `sigma` `[1 x K]`; values are linearly interpolated and held at the endpoints
  - `traj = coorbital.prop.phaseRun(phases,x0,veh,env)` returning a struct with fields `t` `[N x 1]`, `x` `[N x 6]`, `u` `[N x 2]`, `phaseIdx` `[N x 1]`, and `junction` `[(P-1) x 1]` struct array with fields `t` and `x`

- [ ] **Step 1: Write the failing test**

```matlab
function test_phaseRun()
%% Purpose:
%
%  Verify the phase driver terminates on an altitude event, returns a
%  monotonic time vector, records the junction state between phases, and
%  hands the terminal state of one phase to the next without a discontinuity.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();
               veh = coorbital.util.vehicleDefaults();
         env.atmos = @coorbital.atmos.expAtmos;
          env.grav = @coorbital.grav.sphereGrav;
          env.aero = @coorbital.aero.constLD;
        env.omegaE = 0;

%% A level-bank glide that ends when it reaches 20 km:
             sched = struct('tGrid',[0 5000],'alpha',[0 0],'sigma',[0 0]);
       phase1.eom = @coorbital.eom.glide3DOF;
     phase1.guide = @(t,x) coorbital.guide.prescribed(t,x,sched);
 phase1.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,20e3);
     phase1.tspan = [0 4000];

                x0 = [c.rE + 60e3; 0; 0; 6000; deg2rad(-1); deg2rad(90)];
              traj = coorbital.prop.phaseRun(phase1,x0,veh,env);

%% Time is monotonic and the run stopped on the event, not the horizon:
    assert(all(diff(traj.t) > 0),'time must increase monotonically');
    assert(traj.t(end) < 4000,'integration hit the horizon instead of the event');

%% It stopped at 20 km:
              hEnd = traj.x(end,1) - c.rE;
    assert(abs(hEnd - 20e3) < 1,'terminated at %.1f m, expected 20000 m',hEnd);

%% Shapes are as documented:
    assert(size(traj.x,2) == 6,'state must be [N x 6]');
    assert(size(traj.u,2) == 2,'control must be [N x 2]');
    assert(numel(traj.t) == size(traj.x,1),'t and x must have equal length');

%% Two phases: the second continues from the first without a jump:
       phase2 = phase1;
 phase2.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,5e3);
     phase2.tspan = [0 4000];
             traj2 = coorbital.prop.phaseRun([phase1 phase2],x0,veh,env);

    assert(numel(traj2.junction) == 1,'expected one junction between two phases');
              iJun = find(traj2.phaseIdx == 2,1,'first');
             jump = abs(traj2.x(iJun,:) - traj2.junction(1).x(:)');
    assert(max(jump) < 1e-6,'state jumped by %.3e across the junction',max(jump));
    assert(all(diff(traj2.t) >= 0),'time must not run backwards across a junction');
end
```

- [ ] **Step 2: Run it and confirm it fails**

Expected: `FAIL test_phaseRun`.

- [ ] **Step 3: Implement the event**

```matlab
function [value,isterminal,direction] = eventAltitude(t,x,hStop)
%% Purpose:
%
%  ODE event that fires when the vehicle descends through a given altitude.
%  Terminal and one-sided, so a climbing trajectory passing the same altitude
%  does not stop the integration.
%
%% Inputs:
%
%  t                scalar                      Time (s). Unused.
%
%  x                [6 x 1]                     Glide state; uses x(1) = r
%
%  hStop            scalar                      Altitude to stop at (m)
%
%% Outputs:
%
%  value            scalar                      Altitude above hStop (m)
%
%  isterminal       scalar                      1, always terminal
%
%  direction        scalar                      -1, descending crossings only
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();
             value = (x(1) - c.rE) - hStop;
        isterminal = 1;
         direction = -1;
end
```

- [ ] **Step 4: Implement the prescribed guidance**

```matlab
function u = prescribed(t,x,sched)
%% Purpose:
%
%  Evaluate a prescribed control schedule at the current time. Linear
%  interpolation between grid points, held constant outside the grid.
%
%% Inputs:
%
%  t                scalar                      Time since phase start (s)
%
%  x                [6 x 1]                     State. Unused; present so the
%                                               signature matches closed-loop
%                                               guidance laws.
%
%  sched            Struct                      Schedule:
%                                               tGrid [1 x K] (s)
%                                               alpha [1 x K] (rad)
%                                               sigma [1 x K] (rad)
%
%% Outputs:
%
%  u                [2 x 1]                     [alpha; sigma] (rad)
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Single-point schedules are constant; interp1 needs two points:
if numel(sched.tGrid) == 1
                 u = [sched.alpha(1); sched.sigma(1)];
    return;
end

%% Interpolate, clamping outside the grid:
             alpha = interp1(sched.tGrid,sched.alpha,t,'linear','extrap');
             sigma = interp1(sched.tGrid,sched.sigma,t,'linear','extrap');
    if t <= sched.tGrid(1)
             alpha = sched.alpha(1);
             sigma = sched.sigma(1);
    elseif t >= sched.tGrid(end)
             alpha = sched.alpha(end);
             sigma = sched.sigma(end);
    end
                 u = [alpha; sigma];
end
```

- [ ] **Step 5: Implement the phase driver**

```matlab
function traj = phaseRun(phases,x0,veh,env)
%% Purpose:
%
%  Integrate a sequence of trajectory phases, handing the terminal state of
%  each phase to the next and recording the junction states so a downstream
%  optimizer can enforce them as linkage constraints.
%
%% Inputs:
%
%  phases           [1 x P] struct              Each with fields:
%                                               eom       function handle
%                                               guide     function handle
%                                               terminate function handle
%                                               tspan     [1 x 2] (s)
%
%  x0               [6 x 1]                     Initial state
%
%  veh              Struct                      Vehicle parameters
%
%  env              Struct                      Environment model handles
%
%% Outputs:
%
%  traj             Struct                      Trajectory:
%                                               t        [N x 1] (s), cumulative
%                                               x        [N x 6] state
%                                               u        [N x 2] control (rad)
%                                               phaseIdx [N x 1] 1-based phase
%                                               junction [(P-1) x 1] struct
%                                                        with fields t and x
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
                 c = coorbital.util.missileConst();
               veh = coorbital.util.vehicleDefaults();
         env.atmos = @coorbital.atmos.expAtmos;
          env.grav = @coorbital.grav.sphereGrav;
          env.aero = @coorbital.aero.constLD;
        env.omegaE = 0;
             sched = struct('tGrid',[0 5000],'alpha',[0 0],'sigma',[0 0]);
            ph.eom = @coorbital.eom.glide3DOF;
          ph.guide = @(t,x) coorbital.guide.prescribed(t,x,sched);
      ph.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,0);
          ph.tspan = [0 4000];
                x0 = [c.rE + 60e3; 0; 0; 6000; deg2rad(-1); deg2rad(90)];
              traj = coorbital.prop.phaseRun(ph,x0,veh,env);
    figure('color',[1 1 1]);
    plot(traj.t,(traj.x(:,1)-c.rE)./1000,'linewidth',1.5); grid on;
    xlabel('time (s)'); ylabel('altitude (km)');
    return;
end

%% One cell per phase, concatenated once at the end. Segment lengths are not
%% known until ode45 returns, so the cells are the preallocation:
               nPh = numel(phases);
              tSeg = cell(nPh,1);
              xSeg = cell(nPh,1);
              uSeg = cell(nPh,1);
            idxSeg = cell(nPh,1);
          junction = repmat(struct('t',[],'x',[]),max(nPh-1,0),1);
              tOff = 0;
             xCurr = x0(:);

%% Integrate each phase in turn:
for kp = 1:nPh
                ph = phases(kp);
              odeF = @(t,x) ph.eom(t,x,ph.guide(t,x),veh,env);
              opts = odeset('RelTol',1e-10,'AbsTol',1e-10, ...
                            'Events',ph.terminate);
          [tk,xk] = ode45(odeF,ph.tspan,xCurr,opts);

%% Rebuild the control history on the returned grid:
                uk = zeros(numel(tk),2);
    for kt = 1:numel(tk)
         uk(kt,:) = ph.guide(tk(kt),xk(kt,:)').';
    end

%% Drop the duplicated first sample when continuing from a previous phase:
                k0 = 1;
    if kp > 1
                k0 = 2;
    end

         tSeg{kp} = tOff + tk(k0:end);
         xSeg{kp} = xk(k0:end,:);
         uSeg{kp} = uk(k0:end,:);
       idxSeg{kp} = kp*ones(numel(tk)-k0+1,1);

%% Record the junction and carry the state forward:
             xCurr = xk(end,:)';
              tOff = tOff + tk(end);
    if kp < nPh
        junction(kp).t = tOff;
        junction(kp).x = xCurr;
    end
end

%% Assemble:
            traj.t = vertcat(tSeg{:});
            traj.x = vertcat(xSeg{:});
            traj.u = vertcat(uSeg{:});
     traj.phaseIdx = vertcat(idxSeg{:});
     traj.junction = junction;
end
```

- [ ] **Step 6: Run the tests and confirm they pass**

Expected: `6 passed, 0 failed`.

- [ ] **Step 7: Commit**

```bash
cd ~/Desktop/optimal_control
git add missiles/+coorbital/+prop missiles/+coorbital/+guide missiles/tests/test_phaseRun.m
git commit -m "missiles: phase driver, altitude event, prescribed guidance"
```

---

### Task 7: Entry script

**Files:**
- Create: `missiles/HGV/vehicle_hgv.m`
- Create: `missiles/HGV/run_glide.m`

**Interfaces:**
- Consumes: everything from Tasks 1–6
- Produces: `traj = run_glide()` — also prints a summary and draws three figures

- [ ] **Step 1: Write the vehicle file**

```matlab
function veh = vehicle_hgv()
%% Purpose:
%
%  Hypersonic glide vehicle parameters for the run_glide entry script.
%  Starts from the library defaults so there is one place to change a value
%  that should apply everywhere, and overrides only what is vehicle-specific.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  veh              Struct                      Vehicle parameters; see
%                                               coorbital.util.vehicleDefaults
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               veh = coorbital.util.vehicleDefaults();

%% Vehicle-specific overrides (PLACEHOLDER values):
          veh.mass = 900;              %kg
          veh.Sref = 0.75;             %m^2
            veh.CL = 0.35;
            veh.LD = 2.5;
end
```

- [ ] **Step 2: Write the entry script**

```matlab
function traj = run_glide()
%% Purpose:
%
%  Propagate a prescribed-control hypersonic glide from entry interface to a
%  terminal altitude and report range, flight time, and peak conditions.
%  Everything a routine run needs to change lives in the USER PARAMETERS
%  block; nothing below it should require editing.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  traj             Struct                      Trajectory from
%                                               coorbital.prop.phaseRun
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% USER PARAMETERS:
           hEntry = 60;               %km, entry interface altitude
           vEntry = 6000;             %m/s, planet-relative speed at entry
       gammaEntry = -1.0;             %deg, flight path angle, negative is down
         psiEntry = 90;               %deg, heading clockwise from north
          latEntry = 0;               %deg, entry latitude
          lonEntry = 0;               %deg, entry longitude
            hStop = 5;                %km, terminal altitude
        bankAngle = 0;                %deg, constant bank
       alphaAngle = 0;                %deg, constant angle of attack
           tMaxSim = 4000;            %s, integration horizon
         earthSpin = false;           %true to enable Coriolis and centrifugal
%% END USER PARAMETERS

%% Resolve paths so the script runs from anywhere:
              here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here,'..'));

                 c = coorbital.util.missileConst();
               veh = vehicle_hgv();

%% Assemble the environment; swap any handle here to raise fidelity:
         env.atmos = @coorbital.atmos.expAtmos;
          env.grav = @coorbital.grav.sphereGrav;
          env.aero = @coorbital.aero.constLD;
        env.omegaE = 0;
    if earthSpin
        env.omegaE = c.omegaE;
    end

%% Build the control schedule and the single glide phase:
             sched = struct('tGrid',[0 tMaxSim], ...
                            'alpha',deg2rad([alphaAngle alphaAngle]), ...
                            'sigma',deg2rad([bankAngle  bankAngle]));
            ph.eom = @coorbital.eom.glide3DOF;
          ph.guide = @(t,x) coorbital.guide.prescribed(t,x,sched);
      ph.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,hStop*1000);
          ph.tspan = [0 tMaxSim];

                x0 = [c.rE + hEntry*1000; ...
                      deg2rad(lonEntry); ...
                      deg2rad(latEntry); ...
                      vEntry; ...
                      deg2rad(gammaEntry); ...
                      deg2rad(psiEntry)];

%% Propagate:
              traj = coorbital.prop.phaseRun(ph,x0,veh,env);

%% Derived quantities for the summary:
              hKm = (traj.x(:,1) - c.rE)./1000;
                V = traj.x(:,4);
    [rho,~,~,aSnd] = coorbital.atmos.expAtmos(traj.x(:,1) - c.rE);
              qbar = 0.5.*rho.*V.^2;
              mach = V./aSnd;
          rangeKm = c.rE.*coorbital.util.greatCircle(x0(3),x0(2), ...
                        traj.x(end,3),traj.x(end,2))./1000;

%% Report:
    fprintf('\nGlide summary\n');
    fprintf('  flight time      %8.1f s\n',traj.t(end));
    fprintf('  ground range     %8.1f km\n',rangeKm);
    fprintf('  terminal speed   %8.1f m/s\n',V(end));
    fprintf('  peak dyn press   %8.1f kPa\n',max(qbar)/1000);
    fprintf('  peak Mach        %8.1f\n',max(mach));

%% Plots:
    figure('color',[1 1 1]);
    subplot(3,1,1); plot(traj.t,hKm,'linewidth',1.5); grid on;
    ylabel('altitude (km)'); title('Prescribed-control glide');
    subplot(3,1,2); plot(traj.t,V,'linewidth',1.5); grid on;
    ylabel('speed (m/s)');
    subplot(3,1,3); plot(traj.t,qbar./1000,'linewidth',1.5); grid on;
    ylabel('q (kPa)'); xlabel('time (s)');
end
```

- [ ] **Step 3: Add the great-circle helper the entry script calls**

Create `missiles/+coorbital/+util/greatCircle.m`:

```matlab
function d = greatCircle(lat1,lon1,lat2,lon2)
%% Purpose:
%
%  Central angle between two points on a sphere, via the haversine form,
%  which stays accurate for small separations where the cosine form loses
%  precision.
%
%% Inputs:
%
%  lat1             [N x 1]                     Latitude of point 1 (rad)
%
%  lon1             [N x 1]                     Longitude of point 1 (rad)
%
%  lat2             [N x 1]                     Latitude of point 2 (rad)
%
%  lon2             [N x 1]                     Longitude of point 2 (rad)
%
%% Outputs:
%
%  d                [N x 1]                     Central angle (rad); multiply
%                                               by a radius for arc length
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

              dLat = lat2 - lat1;
              dLon = lon2 - lon1;
                 a = sin(dLat./2).^2 + cos(lat1).*cos(lat2).*sin(dLon./2).^2;
                 d = 2.*asin(min(1,sqrt(a)));
end
```

- [ ] **Step 4: Test the great-circle helper**

Create `missiles/tests/test_greatCircle.m`:

```matlab
function test_greatCircle()
%% Purpose:
%
%  Verify the haversine central angle against cases with exact answers, and
%  confirm it stays accurate at small separations where the cosine form
%  degrades. Covers the spec's great-circle range check.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% A quarter turn along the equator is pi/2:
                 d = coorbital.util.greatCircle(0,0,0,pi/2);
    assert(abs(d - pi/2) < 1e-12,'equatorial quarter turn wrong');

%% Equator to pole is pi/2:
                 d = coorbital.util.greatCircle(0,0,pi/2,0);
    assert(abs(d - pi/2) < 1e-12,'equator-to-pole wrong');

%% Coincident points give exactly zero:
                 d = coorbital.util.greatCircle(0.3,1.1,0.3,1.1);
    assert(d == 0,'coincident points must give zero');

%% Antipodal points give pi:
                 d = coorbital.util.greatCircle(0,0,0,pi);
    assert(abs(d - pi) < 1e-9,'antipodal separation wrong');

%% Small separations stay accurate: one arcsecond of latitude:
             dSmall = coorbital.util.greatCircle(0,0,arcsecToRad(1),0);
    assert(abs(dSmall - arcsecToRad(1)) < 1e-15,'small-angle accuracy lost');
end

function r = arcsecToRad(a)
%% Purpose:
%
%  Convert arcseconds to radians.
%
%% ------------------------ Begin Code Sequence ---------------------------
                 r = a.*pi./(180*3600);
end
```

Run the Task 1 Step 3 command. Expected: `PASS test_greatCircle`.

- [ ] **Step 5: Run the entry script**

Run:
```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/missiles/HGV'); run_glide" 2>&1 | grep -vE "Home License|personal use|academic, research|organizational use"
```
Expected: a summary block printing a flight time of order hundreds of seconds, a ground range of order 1000 km, and a peak Mach above 15. If the vehicle skips back out of the atmosphere and hits the 4000 s horizon instead of the 5 km floor, that is physical for `L/D = 2.5` at a shallow entry angle — steepen `gammaEntry` to -3 degrees and re-run.

- [ ] **Step 6: Run the full test suite to confirm nothing regressed**

Expected: `7 passed, 0 failed`.

- [ ] **Step 7: Commit**

```bash
cd ~/Desktop/optimal_control
git add missiles/HGV missiles/+coorbital/+util/greatCircle.m missiles/tests/test_greatCircle.m
git commit -m "missiles: HGV glide entry script with user parameter block"
```

---

### Task 8: Analytic validation

**Files:**
- Create: `missiles/tests/test_equilibriumGlide.m`
- Create: `missiles/tests/test_allenEggers.m`

**Interfaces:**
- Consumes: everything from Tasks 1–6
- Produces: no new interfaces; these are validation tests only

- [ ] **Step 1: Write the equilibrium glide test**

```matlab
function test_equilibriumGlide()
%% Purpose:
%
%  Compare the propagated glide against the equilibrium glide condition. With
%  the bank angle at zero and the flight path angle shallow, lift very nearly
%  balances weight less centrifugal relief:
%
%      0.5 rho V^2 S CL = (m g - m V^2 / r) cos(gamma)
%
%  which rearranges to a closed-form V at each altitude. The propagated speed
%  must track it once the trajectory has settled out of its entry transient.
%  The cos(gamma) factor falls out of setting dgamma/dt = 0; it is within 1e-4
%  of unity at the shallow angles used here, but is carried so the reference is
%  the equilibrium condition exactly rather than an approximation to it.
%
%% References:
%   [1] Vinh, N.X., et al., "Hypersonic and Planetary Entry Flight Mechanics,"
%       Univ. Michigan Press, 1980, Sec. 5.3.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();
               veh = coorbital.util.vehicleDefaults();
         env.atmos = @coorbital.atmos.expAtmos;
          env.grav = @coorbital.grav.sphereGrav;
          env.aero = @coorbital.aero.constLD;
        env.omegaE = 0;

%% Enter shallow and unbanked so the equilibrium assumption holds:
             sched = struct('tGrid',[0 4000],'alpha',[0 0],'sigma',[0 0]);
            ph.eom = @coorbital.eom.glide3DOF;
          ph.guide = @(t,x) coorbital.guide.prescribed(t,x,sched);
      ph.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,25e3);
          ph.tspan = [0 4000];
                x0 = [c.rE + 55e3; 0; 0; 5500; deg2rad(-0.5); deg2rad(90)];
              traj = coorbital.prop.phaseRun(ph,x0,veh,env);

%% Closed-form equilibrium speed at each propagated altitude:
              rArc = traj.x(:,1);
              vArc = traj.x(:,4);
      [rhoArc,~,~,~] = coorbital.atmos.expAtmos(rArc - c.rE);
        [grArc,~] = coorbital.grav.sphereGrav(rArc,traj.x(:,3));
          [CL,~] = coorbital.aero.constLD(0,0,veh);
             cgam = cos(traj.x(:,5));
             vEq2 = (veh.mass.*grArc.*cgam)./ ...
                    (0.5.*rhoArc.*veh.Sref.*CL + veh.mass.*cgam./rArc);
              vEq = sqrt(vEq2);

%% Skip the first 20 percent, which is the entry transient:
                k0 = max(2,round(0.2*numel(vArc)));
             relEr = abs(vArc(k0:end) - vEq(k0:end))./vEq(k0:end);

    assert(max(relEr) < 0.05, ...
        'equilibrium glide mismatch: max relative error %.3f, expected < 0.05', ...
        max(relEr));
end
```

- [ ] **Step 2: Run it and confirm it passes or reveals a real error**

Run the Task 1 Step 3 command. Expected: `PASS test_equilibriumGlide`.

If it fails with an error in the 5–15% band, the likely cause is that the entry transient has not damped by the 20% mark — print `relEr` against altitude and confirm the error decreases monotonically before loosening the window. **Do not loosen the tolerance to make a failing physics check pass.** If the error grows with time, the EOM is wrong and Task 5's tests were not sufficient.

- [ ] **Step 3: Write the Allen–Eggers test**

```matlab
function test_allenEggers()
%% Purpose:
%
%  Compare a zero-lift ballistic entry against the Allen-Eggers closed-form
%  solution. For a steep entry into an exponential atmosphere over a flat
%  gravity field, peak deceleration is
%
%      aMax = V_e^2 sin(gammaE) / (2 e H)
%
%  which is famously independent of the ballistic coefficient -- so this
%  checks the atmosphere and the drag term together, without the vehicle
%  parameters being able to mask an error.
%
%% References:
%   [1] Allen, H.J., Eggers, A.J., "A Study of the Motion and Aerodynamic
%       Heating of Ballistic Missiles Entering the Earth's Atmosphere at High
%       Supersonic Speeds," NACA TR-1381, 1958.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();
               veh = coorbital.util.vehicleDefaults();

%% Zero lift so the vehicle falls ballistically:
         env.atmos = @coorbital.atmos.expAtmos;
          env.grav = @coorbital.grav.sphereGrav;
          env.aero = @(al,ma,vv) deal(0,vv.CL/vv.LD);
        env.omegaE = 0;

%% Steep entry, where the Allen-Eggers assumptions are good:
            gammaE = deg2rad(-30);
                vE = 6000;
             sched = struct('tGrid',[0 600],'alpha',[0 0],'sigma',[0 0]);
            ph.eom = @coorbital.eom.glide3DOF;
          ph.guide = @(t,x) coorbital.guide.prescribed(t,x,sched);
      ph.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,0);
          ph.tspan = [0 600];
                x0 = [c.rE + 120e3; 0; 0; vE; gammaE; deg2rad(90)];
              traj = coorbital.prop.phaseRun(ph,x0,veh,env);

%% Peak AERODYNAMIC deceleration. Differencing the speed history gives the NET
%% rate, which carries gravity's along-track contribution; Allen-Eggers
%% predicts the drag term alone.
              vArc = traj.x(:,4);
              rArc = traj.x(:,1);
    [rhoArc,~,~,~] = coorbital.atmos.expAtmos(rArc - c.rE);
          [~,CDb] = env.aero(0,0,veh);
             accel = 0.5.*rhoArc.*vArc.^2.*veh.Sref.*CDb./veh.mass;
           aMaxNum = max(accel);

%% Allen-Eggers closed form:
           aMaxRef = vE.^2.*abs(sin(gammaE))./(2.*exp(1).*c.Hscale);

               rel = abs(aMaxNum - aMaxRef)./aMaxRef;
    assert(rel < 0.10, ...
        'peak deceleration %.1f m/s^2 vs Allen-Eggers %.1f m/s^2 (%.1f%% off)', ...
        aMaxNum,aMaxRef,100*rel);
end
```

- [ ] **Step 4: Run the full suite**

Expected: `9 passed, 0 failed`.

The Allen–Eggers tolerance is 10% rather than the 2% quoted in the design spec: the closed form neglects gravity's contribution to the speed history and assumes a flat Earth, both of which matter at a 30-degree entry. Ten percent confirms the drag model and atmosphere are right without demanding agreement the approximation cannot deliver.

- [ ] **Step 5: Record the lesson**

Prepend to `docs/LESSONS_LEARNED.md`:

```markdown
## 2026-08-06 — Allen-Eggers is a 10 percent check, not a 2 percent one

The design spec called for 2 percent agreement with the Allen-Eggers peak
deceleration. That is tighter than the approximation supports: the closed form
drops gravity from the speed equation and assumes a flat Earth, which at a
30-degree entry costs several percent on its own. The test is set at 10
percent. It still does its job, because the reference is independent of the
ballistic coefficient -- a wrong drag area or vehicle mass cannot move the
prediction to match a wrong simulation.
```

- [ ] **Step 6: Commit**

```bash
cd ~/Desktop/optimal_control
git add missiles/tests/test_equilibriumGlide.m missiles/tests/test_allenEggers.m missiles/docs/LESSONS_LEARNED.md
git commit -m "missiles: analytic validation against equilibrium glide and Allen-Eggers"
```

---

### Task 9: Organization documentation

**Files:**
- Create: `missiles/docs/README.md`

**Interfaces:**
- Consumes: nothing
- Produces: nothing

- [ ] **Step 1: Write the README**

```markdown
# Missile Trajectory Library

MATLAB 3-DOF trajectory generation for hypersonic glide vehicles, with
ballistic and ICBM variants sharing the same core. Design rationale is in
`DESIGN.md`; this file covers layout and use.

## Running

    cd ~/Desktop/optimal_control/missiles/HGV
    run_glide

Everything a routine run needs is in the `%% USER PARAMETERS:` block at the top
of `run_glide.m`. Nothing below that block should need editing.

Tests:

    cd ~/Desktop/optimal_control/missiles
    run('tests/run_tests')

Headless:

    /Applications/MATLAB_R2025b.app/bin/matlab -batch \
      "cd('/Users/msc/Desktop/optimal_control/missiles'); run('tests/run_tests')"

## Layout

    +coorbital/+atmos    atmosphere models
    +coorbital/+grav     gravity models
    +coorbital/+aero     aerodynamic models
    +coorbital/+eom      equations of motion
    +coorbital/+guide    guidance laws
    +coorbital/+prop     integration driver and ODE events
    +coorbital/+util     constants, vehicle parameters, shared helpers
    HGV/ BM/ ICBM/       vehicle parameters and entry scripts only
    tests/               one test file per unit
    docs/                design, plan, lessons, LaTeX notes

## Adding a fidelity level

Models are injected as function handles, so the equations of motion never name
one. To add a higher-fidelity atmosphere:

1. Write `+coorbital/+atmos/us76Atmos.m` with the signature
   `[rho,P,T,a] = us76Atmos(h)`, `h` in metres, SI out.
2. Give it an `if nargin == 0` self-demo and a test in `tests/`.
3. In the entry script, change one line:

       env.atmos = @coorbital.atmos.us76Atmos;

Gravity models take `[gr,gLat] = f(r,lat)` and aerodynamic models take
`[CL,CD] = f(alpha,mach,veh)`. Rotating-Earth terms are already present in
`glide3DOF`, gated on `env.omegaE` — set it to `c.omegaE` to enable them.

## Conventions

pumpkyn house style throughout: see `~/Desktop/proj7/doc/pumpkyn_style_guide.md`.
Header quartet with aligned input and output columns and units, column-aligned
`=` signs, colon-terminated `%%` step comments, `if nargin == 0` self-demo in
every library function. SI units inside the library; entry scripts may take
degrees and kilometres in the user block and convert immediately.
```

- [ ] **Step 2: Verify every documented command runs**

Run each of the three commands in the Running section. Expected: `run_glide` prints its summary, `run_tests` reports `9 passed, 0 failed`.

- [ ] **Step 3: Commit**

```bash
cd ~/Desktop/optimal_control
git add missiles/docs/README.md
git commit -m "missiles: code organization README"
```

---

## Out of scope for this plan

Deliberately excluded; each becomes its own plan once the glide propagator is validated:

- `+coorbital/+util/stateConvert.m` — the spec provides for phases with differing state definitions, but every Phase 1 phase uses the same 6-state glide vector, so `phaseRun` passes the state straight through. Needed when boost (Cartesian, with mass) hands off to glide.
- `boost3DOF` and a prescribed pitch program; `BM/run_ballistic.m`
- Terminal and descent phase; the full boost → glide → descent chain
- `+viz` package (ground tracks, 3-D trajectory over `pumpkyn.util.earth3D`)
- Fidelity increments: rotating Earth validation, `j2Grav`, geodetic altitude, tabulated aero, US76 atmosphere
- `hgv_dynamics_note.tex` and `software_design.tex` — written after the interfaces survive contact with working code
- Phase 2 optimization against `orbit_transfer/verify_common`
