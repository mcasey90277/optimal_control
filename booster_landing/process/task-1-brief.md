### Task 1: Skeleton, `setup_paths`, `booster_params` + sanity test

**Files:**
- Create: `booster_landing/setup_paths.m`
- Create: `booster_landing/lib/booster_params.m`
- Create: `booster_landing/tests/test_params.m`
- Create: `booster_landing/results/.gitkeep`, `booster_landing/.gitignore`

**Interfaces:**
- Consumes: nothing.
- Produces: `P = booster_params()` — single struct every later task reads:
  `P.mdry P.m0 P.Tmax P.Tmin P.Isp P.g0 P.gvec(3x1) P.r0(3x1) P.v0(3x1) P.gs_deg P.theta_max_deg P.N P.tf_lo P.tf_hi P.pad_radius P.vtd_max P.seed P.drag.{on,rho0,H,Cd,A}` — and `setup_paths()` (adds campaign dirs + CasADi).

- [ ] **Step 1: Write the failing test**

`tests/test_params.m` (house pattern: self-bootstrapping, throws on failure):

```matlab
% TEST_PARAMS  Physical sanity of the booster parameter set.
%
% Checks the three facts the whole campaign leans on:
%   1. m0 > mdry (there is landing propellant);
%   2. min-throttle thrust-to-weight at DRY mass > 1 (hoverslam is forced:
%      the booster cannot hover, so the optimizer must produce a terminal
%      max-thrust arc arriving at v=0 exactly at touchdown);
%   3. bounds ordered: 0 < Tmin < Tmax, tf_lo < tf_hi, gs in (0,90) deg.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P = booster_params();
assert(P.m0 > P.mdry, 'no landing propellant');
twMin = P.Tmin / (P.mdry * P.g0);
assert(twMin > 1, 'min-throttle T/W at dry mass = %.3f, expected > 1', twMin);
assert(P.Tmin > 0 && P.Tmin < P.Tmax, 'thrust bounds disordered');
assert(P.tf_lo < P.tf_hi, 'tf bracket disordered');
assert(P.gs_deg > 0 && P.gs_deg < 90, 'glideslope out of range');
assert(P.drag.on == false, 'Phase 1 default must be vacuum');
fprintf('test_params PASS (min-throttle dry T/W = %.3f)\n', twMin);
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); addpath('tests'); test_params"`
Expected: FAIL — `setup_paths` / `booster_params` undefined.

- [ ] **Step 3: Implement `setup_paths.m` and `booster_params.m`**

`setup_paths.m`:

```matlab
function setup_paths()
% SETUP_PATHS  Paths for the booster_landing campaign (one campaign, one setup).
%
% Adds: this folder (front door), lib (dynamics + solvers + tracking),
% certify (gates), viz (plots + movie), tests, and CasADi 3.7.
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, 'lib'));
addpath(fullfile(here, 'certify'));
addpath(fullfile(here, 'viz'));
addpath(fullfile(here, 'tests'));
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));
end
```

`lib/booster_params.m` — values verbatim from the spec (pumpkyn-style aligned block). `if nargin==0` self-demo prints the struct and the T/W numbers:

```matlab
function P = booster_params()
% BOOSTER_PARAMS  Falcon-9-class 3-DOF landing-burn parameter set.
%
% Single source of truth for physical constants, boundary conditions,
% solver grid and Monte-Carlo settings. Public F9 estimates; every number
% adjustable here and only here.
%
% INPUTS:  none
% OUTPUTS: P - parameter struct (fields documented inline below)
%
% REFERENCES:
%   [1] Blackmore, Acikmese, Scharf, "Minimum-Landing-Error Powered-Descent
%       Guidance for Mars Landing Using Convex Optimization," JGCD 2010.
%   [2] docs/superpowers/specs/2026-08-08-booster-landing-design.md

%% Vehicle (public Falcon 9 block-5 estimates):
P.mdry   = 25600;              % dry mass [kg]
P.m0     = 30000;              % mass at landing-burn start [kg]
P.Tmax   = 845e3;              % one Merlin 1D, sea level [N]
P.Tmin   = 0.40 * P.Tmax;      % ~40 percent min throttle [N]
P.Isp    = 282;                % sea-level Isp [s]
P.g0     = 9.80665;            % standard gravity [m/s^2]
P.gvec   = [0; 0; -P.g0];      % flat-Earth gravity, z up [m/s^2]

%% Boundary conditions (pad at origin, z up):
P.r0     = [500; 100; 2000];   % post-entry-burn position [m]
P.v0     = [-30; 0; -180];     % descending ~180 m/s [m/s]

%% Path constraints:
P.gs_deg        = 30;          % glideslope min elevation angle [deg]
P.theta_max_deg = Inf;         % thrust-pointing cone half-angle, Inf = off

%% Discretization / solver:
P.N      = 60;                 % Hermite-Simpson segments (collocation)
P.Nconv  = 120;                % trapezoid nodes for the convex solver
P.tf_lo  = 10;                 % free-final-time bracket [s]
P.tf_hi  = 50;

%% Atmosphere (Phase 2, OFF by default -- vacuum keeps convexification exact):
P.drag.on   = false;
P.drag.rho0 = 1.225;           % sea-level density [kg/m^3]
P.drag.H    = 8500;            % scale height [m]
P.drag.Cd   = 1.0;             % landing-leg config drag coefficient [-]
P.drag.A    = 10.75;           % reference area, 3.7 m diameter [m^2]

%% Success criteria / Monte Carlo:
P.pad_radius = 15;             % landing accuracy requirement [m]
P.vtd_max    = 2.0;            % max touchdown speed [m/s]
P.seed       = 42;             % rng seed, all random draws
end
```

- [ ] **Step 4: Run test to verify it passes**

Same command as Step 2. Expected: `test_params PASS (min-throttle dry T/W = 1.346)`.

- [ ] **Step 5: `.gitignore` + commit**

`booster_landing/.gitignore`:

```
results/*
!results/.gitkeep
*.asv
```

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/setup_paths.m booster_landing/lib/booster_params.m \
        booster_landing/tests/test_params.m booster_landing/.gitignore booster_landing/results/.gitkeep
git commit -m "booster_landing: skeleton, params, sanity test

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

