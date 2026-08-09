### Task 8: Monte Carlo (`run_monte_carlo`)

**Files:**
- Create: `lib/run_monte_carlo.m`
- Create: `tests/test_monte_carlo_small.m`

**Interfaces:**
- Consumes: `sol`, `ctrl`, `sim_closed_loop`.
- Produces: `mc = run_monte_carlo(sol, ctrl, P, opts)`:
  - `opts.Nrun` (default 200), `opts.sig` struct of 1σ values (defaults: `.r0 = [100;100;50]`, `.v0 = [10;10;10]`, `.thrust = 0.015`, `.isp = 0.01`, `.wind = [10;10;0]` — wind only drawn when `P.drag.on`)
  - `mc.land` (Nrun×2 touchdown xy), `mc.vtd` (Nrun×1), `mc.mprop` (Nrun×1 propellant remaining above dry), `mc.ok` (Nrun×1 logical: miss < pad_radius AND vtd < vtd_max AND m ≥ mdry), `mc.success_rate`
  - Deterministic: `rng(P.seed)` once at top; draws via `randn`.

- [ ] **Step 1: Write the failing test**

`tests/test_monte_carlo_small.m`:

```matlab
% TEST_MONTE_CARLO_SMALL  20-run MC sanity: executes, deterministic under
% the seed (two calls identical), success rate positive, outputs shaped.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
sol  = solve_pdg_colloc(P, struct('N', 30));
ctrl = tvlqr_design(sol, P);
mc1  = run_monte_carlo(sol, ctrl, P, struct('Nrun', 20));
mc2  = run_monte_carlo(sol, ctrl, P, struct('Nrun', 20));
assert(isequal(size(mc1.land), [20 2]) && numel(mc1.vtd) == 20);
assert(isequal(mc1.land, mc2.land), 'MC not deterministic under seed');
assert(mc1.success_rate > 0, 'zero successes at default dispersions');
fprintf('test_monte_carlo_small PASS  success=%.0f%%\n', 100*mc1.success_rate);
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_monte_carlo_small"`
Expected: FAIL — `run_monte_carlo` undefined.

- [ ] **Step 3: Implement `lib/run_monte_carlo.m`**

```matlab
function mc = run_monte_carlo(sol, ctrl, P, opts)
% RUN_MONTE_CARLO  Dispersed closed-loop landing campaign.
%
% Draws initial-state, thrust-bias, Isp-bias (and, with drag on, wind)
% dispersions and runs sim_closed_loop per sample. Deterministic under
% P.seed. Success = miss < pad_radius AND vtd < vtd_max AND m >= mdry.
%
% INPUTS:  sol, ctrl, P as usual; opts.Nrun [def 200], opts.sig 1-sigma
%          struct (see plan Task 8)
% OUTPUTS: mc - .land(Nx2) .vtd(Nx1) .mprop(Nx1) .ok(Nx1) .success_rate
if nargin < 4, opts = struct(); end
if ~isfield(opts,'Nrun'), opts.Nrun = 200; end
sig = struct('r0',[100;100;50], 'v0',[10;10;10], 'thrust',0.015, ...
             'isp',0.01, 'wind',[10;10;0]);
if isfield(opts,'sig')
    fn = fieldnames(opts.sig);
    for k = 1:numel(fn), sig.(fn{k}) = opts.sig.(fn{k}); end
end

rng(P.seed);
Nr = opts.Nrun;
mc.land = zeros(Nr,2);  mc.vtd = zeros(Nr,1);
mc.mprop = zeros(Nr,1); mc.ok  = false(Nr,1);
for krun = 1:Nr
    d = struct('dr0', sig.r0 .* randn(3,1), ...
               'dv0', sig.v0 .* randn(3,1), ...
               'thrust_scale', 1 + sig.thrust*randn(), ...
               'isp_scale',    1 + sig.isp*randn());
    if P.drag.on, d.wind = sig.wind .* randn(3,1); end
    out = sim_closed_loop(sol, ctrl, P, d);
    mc.land(krun,:) = out.td.r(1:2).';
    mc.vtd(krun)    = out.td.vtd;
    mc.mprop(krun)  = out.td.m - P.mdry;
    mc.ok(krun)     = out.td.miss < P.pad_radius && ...
                      out.td.vtd < P.vtd_max && out.td.m >= P.mdry;
    if mod(krun, 25) == 0, fprintf('  MC %d/%d\n', krun, Nr); end
end
mc.success_rate = mean(mc.ok);
end
```

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: PASS. Note the small-N success rate; the spec's ≥95% criterion is judged at Nrun=200 in the front door, not here.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/lib/run_monte_carlo.m booster_landing/tests/test_monte_carlo_small.m
git commit -m "booster_landing: dispersed Monte-Carlo landing campaign

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

