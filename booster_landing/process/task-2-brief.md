### Task 2: `pdg_dynamics` with analytic Jacobians + complex-step test

**Files:**
- Create: `lib/pdg_dynamics.m`
- Create: `tests/test_dynamics_jac.m`

**Interfaces:**
- Consumes: `P` from `booster_params`.
- Produces: `[xdot, A, B] = pdg_dynamics(x, T, P)` — `x = [r; v; m]` (7x1), `T` thrust vector [N] (3x1), `xdot` 7x1, `A = ∂f/∂x` 7x7, `B = ∂f/∂T` 7x3. Honors `P.drag.on` (Phase 2 branch written NOW, tested now, default off — the opt-in-flag house pattern). Must be complex-step safe: no `norm`, no `abs`, no `max`.

- [ ] **Step 1: Write the failing test**

`tests/test_dynamics_jac.m`:

```matlab
% TEST_DYNAMICS_JAC  Complex-step vs analytic Jacobians of pdg_dynamics.
%
% Perturbs each of the 10 inputs (7 state + 3 thrust) with h=1e-30i and
% compares imag(f)/h against the analytic A and B columns. Run twice:
% vacuum (Phase 1 default) and with the drag branch forced on, at a state
% with nonzero velocity (the ||v||v drag term is not differentiable at
% v=0; we never linearize there).
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P  = booster_params();
x0 = [400; 80; 1500; -25; 5; -140; 28000];
T0 = [2e4; -1e4; 5e5];
h  = 1e-30;

for dragOn = [false true]
    P.drag.on = dragOn;
    [~, A, B] = pdg_dynamics(x0, T0, P);
    Acs = zeros(7,7);  Bcs = zeros(7,3);
    for k = 1:7
        xp = complex(x0);  xp(k) = xp(k) + 1i*h;
        Acs(:,k) = imag(pdg_dynamics(xp, complex(T0), P)) / h;
    end
    for k = 1:3
        Tp = complex(T0);  Tp(k) = Tp(k) + 1i*h;
        Bcs(:,k) = imag(pdg_dynamics(complex(x0), Tp, P)) / h;
    end
    errA = max(abs(A(:) - Acs(:))) / max(1, max(abs(Acs(:))));
    errB = max(abs(B(:) - Bcs(:))) / max(1, max(abs(Bcs(:))));
    assert(errA < 1e-12, 'A mismatch (drag=%d): %.2e', dragOn, errA);
    assert(errB < 1e-12, 'B mismatch (drag=%d): %.2e', dragOn, errB);
end
fprintf('test_dynamics_jac PASS (vacuum + drag)\n');
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_dynamics_jac"`
Expected: FAIL — `pdg_dynamics` undefined.

- [ ] **Step 3: Implement `lib/pdg_dynamics.m`**

```matlab
function [xdot, A, B] = pdg_dynamics(x, T, P)
% PDG_DYNAMICS  3-DOF powered-descent dynamics with analytic Jacobians.
%
%   rdot = v
%   vdot = g + T/m + aD          aD = -(1/2) rho(z) Cd A |v| v / m  (opt-in)
%   mdot = -|T| / (Isp g0)
%
% Complex-step safe: magnitudes via sqrt(sum(.^2)), no norm/abs/max.
%
% INPUTS:
%   x - state [r(3); v(3); m] [7x1]  (SI: m, m/s, kg)
%   T - thrust vector [3x1] [N]
%   P - booster_params struct (uses gvec, Isp, g0, drag.*)
%
% OUTPUTS:
%   xdot - state derivative [7x1]
%   A    - d(xdot)/dx [7x7]   (only when requested)
%   B    - d(xdot)/dT [7x3]   (only when requested)
%
% REFERENCES:
%   [1] Acikmese & Ploen, "Convex Programming Approach to Powered Descent
%       Guidance for Mars Landing," JGCD 2007. (vacuum model)
v    = x(4:6);   m = x(7);
Tmag = sqrt(sum(T.^2));
aT   = T / m;
aD   = zeros(3,1);
if P.drag.on
    rho  = P.drag.rho0 * exp(-x(3) / P.drag.H);
    vmag = sqrt(sum(v.^2));
    kD   = 0.5 * rho * P.drag.Cd * P.drag.A;      % so aD = -kD |v| v / m
    aD   = -kD * vmag * v / m;
end
xdot = [v; P.gvec + aT + aD; -Tmag / (P.Isp * P.g0)];

if nargout > 1
    A          = zeros(7,7);
    A(1:3,4:6) = eye(3);
    A(4:6,7)   = -T / m^2;
    if P.drag.on
        % d(aD)/dv = -(kD/m)(|v| I + v v'/|v|); d/dz via drho/dz = -rho/H;
        % d/dm = +kD |v| v / m^2
        A(4:6,4:6) = -(kD/m) * (vmag*eye(3) + (v*v.')/vmag);
        A(4:6,3)   = A(4:6,3) + (kD*vmag/(P.drag.H*m)) * v;
        A(4:6,7)   = A(4:6,7) + kD * vmag * v / m^2;
    end
    B        = zeros(7,3);
    B(4:6,:) = eye(3) / m;
    B(7,:)   = -T.' / (Tmag * P.Isp * P.g0);
end
end
```

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: `test_dynamics_jac PASS (vacuum + drag)`.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/lib/pdg_dynamics.m booster_landing/tests/test_dynamics_jac.m
git commit -m "booster_landing: 3-DOF dynamics + analytic Jacobians, complex-step verified

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

