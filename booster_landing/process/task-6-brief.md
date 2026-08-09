### Task 6: TVLQR design (`tvlqr_design`)

**Files:**
- Create: `lib/tvlqr_design.m`
- Create: `tests/test_tvlqr_riccati.m`

**Interfaces:**
- Consumes: `sol` (colloc), `pdg_dynamics` analytic A,B.
- Produces: `ctrl = tvlqr_design(sol, P, opts)`:
  - `ctrl.tgrid` (1×M dense, M = 4·(N+1)), `ctrl.K` (3×7×M) gains, `ctrl.Pt` (7×7×M) Riccati solution
  - `ctrl.xnom` / `ctrl.Tnom` — pchip-interpolant handles `@(t)->7x1 / 3x1` of the nominal trajectory/control (what `sim_closed_loop` tracks)
  - `opts.Q` (default `diag([1e-4 1e-4 1e-4 1e-2 1e-2 1e-2 0])`), `opts.R` (default `1e-10*eye(3)`), `opts.Qf` (default `diag([1e-2 1e-2 1e-2 1 1 1 0])`) — SI units; mass column zero (mass is observable, not directly regulated)

- [ ] **Step 1: Write the failing test**

`tests/test_tvlqr_riccati.m`:

```matlab
% TEST_TVLQR_RICCATI  Riccati solution health along the nominal trajectory:
% P(t) symmetric positive semidefinite everywhere (PSD, not PD: the mass
% row/column is unweighted), gains finite, terminal condition P(tf)=Qf.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
sol  = solve_pdg_colloc(P, struct('N', 30));
assert(sol.stats.success);
ctrl = tvlqr_design(sol, P);
M    = numel(ctrl.tgrid);
for k = 1:M
    Pk = ctrl.Pt(:,:,k);
    assert(max(max(abs(Pk - Pk.'))) < 1e-6 * max(1,max(abs(Pk(:)))), ...
           'P not symmetric at k=%d', k);
    assert(min(eig((Pk+Pk.')/2)) > -1e-8 * max(1,max(abs(Pk(:)))), ...
           'P not PSD at k=%d', k);
end
assert(all(isfinite(ctrl.K(:))), 'non-finite gains');
Qf = diag([1e-2 1e-2 1e-2 1 1 1 0]);
assert(max(max(abs(ctrl.Pt(:,:,end) - Qf))) < 1e-9, 'terminal P ~= Qf');
fprintf('test_tvlqr_riccati PASS\n');
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_tvlqr_riccati"`
Expected: FAIL — `tvlqr_design` undefined.

- [ ] **Step 3: Implement `lib/tvlqr_design.m`**

```matlab
function ctrl = tvlqr_design(sol, P, opts)
% TVLQR_DESIGN  Time-varying LQR about a PDG guidance trajectory.
%
% Linearizes pdg_dynamics along (x*(t), T*(t)) and integrates the
% differential Riccati equation backward:
%   -Pdot = A'P + PA - PBR^{-1}B'P + Q,   P(tf) = Qf
% Gains K(t) = R^{-1} B(t)' P(t), stored on a dense grid; controller is
%   T_cmd(t) = T*(t) - K(t) (x - x*(t))  (saturation applied by the sim).
%
% INPUTS:
%   sol  - collocation solution (Task 3 interface)
%   P    - booster_params
%   opts - (optional) .Q .R .Qf weight overrides (7x7, 3x3, 7x7)
% OUTPUTS:
%   ctrl - .tgrid(1xM) .K(3x7xM) .Pt(7x7xM) .xnom(@t->7x1) .Tnom(@t->3x1)
%
% REFERENCES:
%   [1] Anderson & Moore, "Optimal Control: Linear Quadratic Methods."
if nargin < 3, opts = struct(); end
if ~isfield(opts,'Q'),  opts.Q  = diag([1e-4 1e-4 1e-4 1e-2 1e-2 1e-2 0]); end
if ~isfield(opts,'R'),  opts.R  = 1e-10*eye(3);                            end
if ~isfield(opts,'Qf'), opts.Qf = diag([1e-2 1e-2 1e-2 1 1 1 0]);          end

%% Nominal interpolants (nodes+midpoints for thrust, nodes for state):
Nn = size(sol.X,2) - 1;  h = sol.tf/Nn;
tU = sort([sol.t, sol.t(1:end-1) + h/2]);
TU = zeros(3, numel(tU));  TU(:,1:2:end) = sol.U;  TU(:,2:2:end) = sol.Um;
ctrl.xnom = @(t) interp1(sol.t.', sol.X.', clampt(t, sol.tf), 'pchip').';
ctrl.Tnom = @(t) interp1(tU.',   TU.',    clampt(t, sol.tf), 'pchip').';

%% Backward Riccati on vec(P), dense output grid:
M     = 4*(Nn+1);
ctrl.tgrid = linspace(0, sol.tf, M);
Rinv  = inv(opts.R);
ric   = @(t, pv) ricrhs(t, pv, ctrl, P, opts.Q, Rinv);
oo    = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[~, PV] = ode45(ric, fliplr(ctrl.tgrid), opts.Qf(:), oo);   % integrates tf->0
PV    = flipud(PV);                                          % re-order 0->tf
ctrl.Pt = zeros(7,7,M);  ctrl.K = zeros(3,7,M);
for k = 1:M
    Pk = reshape(PV(k,:), 7, 7);  Pk = (Pk + Pk.')/2;
    [~, ~, Bk] = pdg_dynamics(ctrl.xnom(ctrl.tgrid(k)), ...
                              ctrl.Tnom(ctrl.tgrid(k)), P);
    ctrl.Pt(:,:,k) = Pk;
    ctrl.K(:,:,k)  = Rinv * (Bk.' * Pk);
end
end

function t = clampt(t, tf), t = min(max(t, 0), tf); end

function pdot = ricrhs(t, pv, ctrl, P, Q, Rinv)
Pk = reshape(pv, 7, 7);
[~, A, B] = pdg_dynamics(ctrl.xnom(t), ctrl.Tnom(t), P);
Pd   = -(A.'*Pk + Pk*A - Pk*B*Rinv*B.'*Pk + Q);
pdot = Pd(:);
end
```

Note: `ode45(ric, fliplr(tgrid), ...)` with a decreasing tspan integrates backward directly — verify PV row order corresponds to the flipped grid before `flipud` (check `PV(1,:)` reshapes to ≈ Qf).

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: `test_tvlqr_riccati PASS`.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/lib/tvlqr_design.m booster_landing/tests/test_tvlqr_riccati.m
git commit -m "booster_landing: TVLQR gains via backward Riccati along guidance

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

