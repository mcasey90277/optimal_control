# PMP-Residual-Driven Adaptive Mesh Refinement (prototype) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a prototype loop that uses the direct min-fuel solution's own KKT-dual PMP switching function to localize mesh under-resolution at bang-bang switches, refine the σ collocation mesh there, and re-solve the direct Sundman solver — demonstrating that switch times stabilize (lose their node-boundary pinning) at essentially fixed propellant.

**Architecture:** Four focused MATLAB functions in a new `sundman_minfuel/refine/` folder — an indicator (reuses the validated `ms_band` mode-'d' dual→costate map + `sms_eom`), a σ-mesh refiner, a no-resample warm-start builder, and a driver loop — plus a seed-prep helper and a results note. The direct solver `casadi_minfuel_sundman.m` is reused unchanged (it already accepts a per-interval σ mesh). No LM/shooting anywhere; the indirect machinery is a *measurement* tool only.

**Tech Stack:** MATLAB R2025b, CasADi 3.7.0 + IPOPT (via `casadi_minfuel_sundman`), the `ms_band` Sundman-domain PMP machinery (`sms_seed_duals`, `sms_eom`, `sms_problem`, `beta_from_duals`).

## Global Constraints

- **MATLAB R2025b only:** `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('<dir>'); script"` (R2025a is license-broken). Multi-line `-batch` strings fail through the shell — always write a `.m` file and run it.
- **Filter the license banner** on every run: pipe through `grep -v -i "home license\|personal use\|academic, research\|organizational use"`.
- **CasADi** at `~/casadi-3.7.0` (`CASADI_PATH`); prebuilt arm64 loads, MEX compiler is broken — never build. `casadi_minfuel_sundman` adds it to the path itself.
- **No-resample discipline (hard requirement):** refining the mesh means INSERTING σ nodes and interpolating ONLY the inserts; every original node's state/control value is copied verbatim. Resampling a 40-rev trajectory onto a fresh mesh pins IPOPT in restoration (`OPTIMALITY_VERIFICATION_PLAN.md` §F.3).
- **Re-solve mode:** `casadi_minfuel_sundman(..., epsilon=0, warmTight=true)` — tight warm start for re-solving AT a bang-bang point. The wrong warm mode wedges IPOPT.
- **Sundman exponent** `pSund = 1.5` everywhere (campaign constant).
- **MATLAB function-header standard** (purpose / inputs w/ sizes / outputs w/ sizes / references) on every function. Never use `i`/`j` as loop indices (imaginary unit) — use `k`, `kk`, `q`, `ii`.
- **IPOPT writes a stray `fort.6`** in cwd — ignore/delete; never commit it.
- **Indicator source (design decision):** switching-function localization (option 1) DRIVES refinement; the Hamiltonian residual |H_σ| is computed as a PASSIVE diagnostic only (option 2 is a future escalation, not built here).

Spec: `docs/superpowers/specs/2026-07-11-pmp-mesh-refine-design.md`.

---

### Task 1: PMP switch-localization indicator

**Files:**
- Create: `NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/pmp_refine_indicator.m`
- Test: `NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/test_pmp_refine_indicator.m`

**Interfaces:**
- Consumes: `ms_band/sms_seed_duals(matFile, M, epsEval, mode) → [Zseed, prob, info]` where `info.tauN [1×nN]`, `info.Y16 [16×nN]` (rows 1:8 states `[r;v;m;t]`, 9:11 λ_r, 12:14 λ_v, 15 λ_m, 16 λ_t), `info.X [8×nN]`, `info.U [4×nN]`, `info.beta`, `info.spreadPct`; `prob.c`, `prob.Tmax`, `prob.muStar`, `prob.pSund`. Consumes `ms_band/sms_eom(σ, Y, Tmax, c, muStar, epsSmooth, pSund) → [dY, Ht, S, u]`.
- Produces: `[score, tauSwitch, diag] = pmp_refine_indicator(seedFile, opts)` — `score [1×N]` per-interval refinement score (≥0), `tauSwitch [1×nsw]` direct-throttle switch times (τ, sorted), `diag` struct with fields `Snode [1×nN]`, `tauN [1×nN]`, `tauCr [1×ncr]` (S=0 crossings), `swI [1×nsw]` (direct switch interval indices), `Hres [1×nN]`, `HresMax`, `nViol`, `betaSpread`, `beta`. Later tasks rely on these exact names.

- [ ] **Step 1: Write the failing test**

Create `test_pmp_refine_indicator.m`:

```matlab
function test_pmp_refine_indicator()
% TEST_PMP_REFINE_INDICATOR  Indicator sanity vs the known-good 1.12x file.
%
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', '..', 'ms_band'));   run_setup_paths(here);
seed = fullfile(here, '..', 'results', 'minfuel', 'legacy_ms_f1120.mat');
assert(isfile(seed), 'missing seed file %s', seed);

opts = struct('M', 40, 'epsEval', 1e-4, 'mode', 'd', 'nbr', 3);
[score, tauSwitch, diag] = pmp_refine_indicator(seed, opts);

nN = numel(diag.Snode);  N = nN - 1;
assert(isequal(size(score), [1 N]), 'score must be 1xN');
assert(all(score >= 0), 'score must be nonnegative');
assert(any(score > 0), 'some intervals must score > 0 (switches present)');
assert(~isempty(tauSwitch) && issorted(tauSwitch), 'tauSwitch nonempty & sorted');
assert(numel(tauSwitch) >= 8 && numel(tauSwitch) <= 14, ...
       'expect ~10-12 direct switches at 1.12x, got %d', numel(tauSwitch));
assert(diag.betaSpread < 5, 'beta spread should be small (mode-d), got %.2f%%', diag.betaSpread);
assert(isequal(size(diag.Hres), [1 nN]), 'Hres must be 1xnN');
assert(diag.nViol >= 0, 'nViol counts');
fprintf('ALL PASS (nSwitch=%d, betaSpread=%.2f%%, HresMax=%.2e, nViol=%d)\n', ...
        numel(tauSwitch), diag.betaSpread, diag.HresMax, diag.nViol);
end

function run_setup_paths(here)
% add ms_band + its dependency chain (sundman_minfuel, lowThrust, pumpkyn)
old = cd(fullfile(here, '..', '..', 'ms_band'));  c = onCleanup(@() cd(old));
setup_paths();
end
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/refine
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('$(pwd)'); test_pmp_refine_indicator" 2>&1 | grep -v -i "home license\|personal use\|academic, research\|organizational use"
```
Expected: FAIL with `Unrecognized function or variable 'pmp_refine_indicator'`.

- [ ] **Step 3: Write minimal implementation**

Create `pmp_refine_indicator.m`:

```matlab
function [score, tauSwitch, diag] = pmp_refine_indicator(seedFile, opts)
% PMP_REFINE_INDICATOR  Per-interval mesh-refinement score from the PMP
%   switching function, plus passive Hamiltonian-residual diagnostics.
%
% Recovers the discrete costates from a direct solution's KKT defect duals
% (mode-'d' midpoint map, ms_band/sms_seed_duals), forms the min-fuel
% switching function S(tau) = 1 - ||lamV||*c/m - lamM on the node grid,
% and scores each collocation interval by how poorly the mesh localizes an
% S=0 crossing (a switch): wide intervals with a central crossing score
% highest. Also returns, as PASSIVE diagnostics only, the Sundman-domain
% Hamiltonian residual |kappa*(Ht+lamT)| per node and the count of
% switching-law sign violations outside a switch deadband.
%
% INPUTS:
%   seedFile - .mat with out.X [8x(N+1)], out.U [4x(N+1)], out.lamDef [8xN],
%              factor, tauf0, sigma [(N+1)x1] (sms_seed_duals input layout)
%   opts     - struct: M arcs for dual map [default 40], epsEval smoothing
%              [default 1e-4], mode dual->costate map [default 'd'], nbr
%              neighbor half-window for scoring [default 3]
%
% OUTPUTS:
%   score     - per-interval refinement score [1xN], >= 0
%   tauSwitch - direct-throttle switch times in tau, sorted [1xnsw]
%   diag      - struct: Snode [1xnN], tauN [1xnN], tauCr [1xncr] (S=0
%               crossings), swI [1xnsw], Hres [1xnN], HresMax, nViol,
%               betaSpread, beta
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-11-pmp-mesh-refine-design.md
%   [2] ms_band/verify_direct_pmp.m (Snode + dual-map reuse)

if nargin < 2, opts = struct(); end
if ~isfield(opts, 'M'),       opts.M = 40;        end
if ~isfield(opts, 'epsEval'), opts.epsEval = 1e-4; end
if ~isfield(opts, 'mode'),    opts.mode = 'd';    end
if ~isfield(opts, 'nbr'),     opts.nbr = 3;       end

[~, prob, info] = sms_seed_duals(seedFile, opts.M, opts.epsEval, opts.mode);
tauN = info.tauN;   Y16 = info.Y16;   X = info.X;   U = info.U;
nN = size(X, 2);    N = nN - 1;

% switching function on the node grid (same construction as verify_direct_pmp)
lamV  = Y16(12:14, :);   lamM = Y16(15, :);   lamT = Y16(16, :);
Snode = 1 - sqrt(sum(lamV.^2, 1)).*prob.c./X(7, :) - lamM;

% S=0 crossings, localized by linear interpolation within each bracket
crossI = find(diff(sign(Snode)) ~= 0);        % node index just before crossing
tauCr  = zeros(1, numel(crossI));
for q = 1:numel(crossI)
    k = crossI(q);
    tauCr(q) = tauN(k) + (0 - Snode(k))*(tauN(k+1) - tauN(k))/(Snode(k+1) - Snode(k));
end

% per-interval score: for each crossing, spread weight over +-nbr intervals,
% weighted by (normalized local width) * (centrality of the crossing in [0,0.5])
score = zeros(1, N);
hInt  = diff(tauN);                            % [1xN] interval widths in tau
sigf  = tauN(end);
for q = 1:numel(crossI)
    kc  = crossI(q);
    off = min(tauCr(q) - tauN(kc), tauN(kc+1) - tauCr(q)) / (tauN(kc+1) - tauN(kc));
    for kk = max(1, kc-opts.nbr):min(N, kc+opts.nbr)
        score(kk) = score(kk) + (hInt(kk)/sigf)*off;
    end
end

% direct-throttle switch times
s     = U(4, :);
swI   = find(diff(double(s > 0.5)) ~= 0);      % switch interval indices [1xnsw]
tauSwitch = sort((tauN(swI) + tauN(swI+1))/2);

% ---- PASSIVE diagnostics ---------------------------------------------------
rE   = [-prob.muStar; 0; 0];
Hres = zeros(1, nN);
for k = 1:nN
    [~, Htk] = sms_eom(0, Y16(:, k), prob.Tmax, prob.c, prob.muStar, ...
                       opts.epsEval, prob.pSund);
    r1 = sqrt(sum((X(1:3, k) - rE).^2));
    Hres(k) = abs(r1^prob.pSund * (Htk + lamT(k)));
end
% switching-law sign violations outside a +-3-node deadband of a direct switch
viol = sign(Snode) ~= sign(0.5 - s);
dead = false(1, nN);
for w = -3:3
    idx = swI + w;  idx = idx(idx >= 1 & idx <= nN);  dead(idx) = true;
end
nViol = nnz(viol & ~dead);

diag = struct('Snode', Snode, 'tauN', tauN, 'tauCr', tauCr, 'swI', swI, ...
              'Hres', Hres, 'HresMax', max(Hres), 'nViol', nViol, ...
              'betaSpread', info.spreadPct, 'beta', info.beta);
end
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/refine
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('$(pwd)'); test_pmp_refine_indicator" 2>&1 | grep -v -i "home license\|personal use\|academic, research\|organizational use"
```
Expected: `ALL PASS (nSwitch=..., betaSpread=0.4...%, HresMax=..., nViol=...)`. (betaSpread ≈ 0.45% for `legacy_ms_f1120` mode-'d' per the campaign record.)

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/pmp_refine_indicator.m NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/test_pmp_refine_indicator.m
git commit -m "refine: PMP switch-localization indicator + passive Hamiltonian diagnostic"
```

---

### Task 2: σ-mesh refiner

**Files:**
- Create: `NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/refine_sigma.m`
- Test: `NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/test_refine_sigma.m`

**Interfaces:**
- Consumes: `score [1×N]` from Task 1.
- Produces: `[sigmaNew, isNew, nDropped] = refine_sigma(sigma, score, opts)` — `sigmaNew [(N'+1)×1]` sorted, contains every original node exactly; `isNew [(N'+1)×1] logical` marks inserted nodes; `nDropped` count of top-scored intervals skipped by guards. `opts` fields: `K` (max intervals to bisect, default `min(8, nnz(score>0))`), `hFloor` (min σ-interval width to bisect, default 1e-9), `maxAdd` (cap on inserted nodes per call, default 40).

- [ ] **Step 1: Write the failing test**

Create `test_refine_sigma.m`:

```matlab
function test_refine_sigma()
% TEST_REFINE_SIGMA  Bisection preserves originals, marks inserts, honors guards.
%
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));  addpath(here);
sigma = linspace(0, 1, 11).';        % 10 intervals
score = zeros(1, 10);  score(3) = 5;  score(7) = 9;   % two hot intervals
opts  = struct('K', 2, 'hFloor', 1e-9, 'maxAdd', 40);
[sigmaNew, isNew, nDropped] = refine_sigma(sigma, score, opts);

assert(numel(sigmaNew) == 13, 'expected 11+2 nodes, got %d', numel(sigmaNew));
assert(nnz(isNew) == 2, 'exactly 2 inserted');
assert(issorted(sigmaNew), 'sigmaNew sorted');
assert(nDropped == 0, 'no drops expected');
% every original node still present (exact)
for k = 1:numel(sigma)
    assert(any(abs(sigmaNew - sigma(k)) < 1e-15), 'original node %d lost', k);
end
% inserts are the midpoints of intervals 3 and 7
mids = sort([(sigma(3)+sigma(4))/2; (sigma(7)+sigma(8))/2]);
got  = sort(sigmaNew(isNew));
assert(max(abs(got - mids)) < 1e-15, 'inserts must be interval midpoints');

% guard: hFloor skips a hot but too-thin interval
sigma2 = [0; 1e-10; 0.5; 1];  score2 = [9 0 0];   % interval 1 is 1e-10 wide
opts2  = struct('K', 1, 'hFloor', 1e-9, 'maxAdd', 40);
[sn2, in2, nd2] = refine_sigma(sigma2, score2, opts2);
assert(nnz(in2) == 0 && nd2 == 1, 'sub-hFloor interval must be dropped');

fprintf('ALL PASS\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/refine
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('$(pwd)'); test_refine_sigma" 2>&1 | grep -v -i "home license\|personal use\|academic, research\|organizational use"
```
Expected: FAIL with `Unrecognized function or variable 'refine_sigma'`.

- [ ] **Step 3: Write minimal implementation**

Create `refine_sigma.m`:

```matlab
function [sigmaNew, isNew, nDropped] = refine_sigma(sigma, score, opts)
% REFINE_SIGMA  Bisect the top-scoring collocation intervals, preserving nodes.
%
% Inserts the midpoint of each selected sigma-interval, keeping every
% original node exactly (the no-resample discipline). Guards: never bisect
% an interval narrower than hFloor; cap inserted nodes per call at maxAdd.
% Dropped selections are counted (never silently truncated).
%
% INPUTS:
%   sigma - current normalized nodes, 0->1 [(N+1)x1]
%   score - per-interval refinement score [1xN], >= 0
%   opts  - struct: K max intervals to bisect [default min(8,nnz(score>0))],
%           hFloor min sigma-interval width to bisect [default 1e-9],
%           maxAdd cap on inserted nodes [default 40]
%
% OUTPUTS:
%   sigmaNew - refined nodes, sorted, originals preserved [(N'+1)x1]
%   isNew    - logical mask of inserted nodes [(N'+1)x1]
%   nDropped - number of top-scored intervals skipped by a guard [scalar]
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-11-pmp-mesh-refine-design.md

sigma = sigma(:);   N = numel(sigma) - 1;   score = score(:).';
if nargin < 3, opts = struct(); end
if ~isfield(opts, 'K'),      opts.K = min(8, nnz(score > 0)); end
if ~isfield(opts, 'hFloor'), opts.hFloor = 1e-9;             end
if ~isfield(opts, 'maxAdd'), opts.maxAdd = 40;               end

h = diff(sigma).';                             % [1xN] interval widths
[~, ord] = sort(score, 'descend');
sel = [];  nDropped = 0;
for q = 1:numel(ord)
    if numel(sel) >= min(opts.K, opts.maxAdd), break; end
    k = ord(q);
    if score(k) <= 0, break; end               % no more scored intervals
    if h(k) < opts.hFloor, nDropped = nDropped + 1; continue; end
    sel(end+1) = k; %#ok<AGROW>
end
add = false(1, N);  add(sel) = true;

sigmaNew = zeros(N + 1 + numel(sel), 1);  isNew = false(size(sigmaNew));
w = 1;
for k = 1:N
    sigmaNew(w) = sigma(k);  w = w + 1;
    if add(k)
        sigmaNew(w) = 0.5*(sigma(k) + sigma(k+1));  isNew(w) = true;  w = w + 1;
    end
end
sigmaNew(w) = sigma(N + 1);
end
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/refine
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('$(pwd)'); test_refine_sigma" 2>&1 | grep -v -i "home license\|personal use\|academic, research\|organizational use"
```
Expected: `ALL PASS`.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/refine_sigma.m NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/test_refine_sigma.m
git commit -m "refine: top-K interval bisection with hFloor/maxAdd guards"
```

---

### Task 3: No-resample warm-start builder

**Files:**
- Create: `NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/warmstart_on_mesh.m`
- Test: `NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/test_warmstart_on_mesh.m`

**Interfaces:**
- Consumes: `out` (with `out.X [8×nN]`, `out.U [4×nN]`), `sigma [nN×1]`, `sigmaNew [nN'×1]`, `isNew [nN'×1] logical` from Tasks 1–2.
- Produces: `[X0, U0] = warmstart_on_mesh(out, sigma, sigmaNew, isNew)` — `X0 [8×nN']`, `U0 [4×nN']`; original nodes copied verbatim, inserts pchip-interpolated for state and thrust direction (renormalized), throttle step-held from the pre-switch (left) original node.

- [ ] **Step 1: Write the failing test**

Create `test_warmstart_on_mesh.m`:

```matlab
function test_warmstart_on_mesh()
% TEST_WARMSTART_ON_MESH  Originals exact; insert straddling a switch holds left.
%
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));  addpath(here);
sigma = linspace(0, 1, 6).';                 % 5 intervals, 6 nodes
X = [ (1:6);  (11:16);  (21:26);  (31:36); ...
      (41:46); (51:56);  linspace(1,0.9,6);  linspace(0,1,6) ];   % [8x6]
% throttle switches between node 3 and node 4 (burn -> coast)
s   = [1 1 1 0 0 0];
al  = repmat([1;0;0], 1, 6);                 % unit +x direction
U   = [al; s];                                % [4x6]
out = struct('X', X, 'U', U);

% refine the interval that straddles the switch (interval 3: nodes 3-4)
sigmaNew = sort([sigma; 0.5*(sigma(3)+sigma(4))]);
isNew    = false(size(sigmaNew));  isNew(abs(sigmaNew - 0.5*(sigma(3)+sigma(4))) < 1e-15) = true;

[X0, U0] = warmstart_on_mesh(out, sigma, sigmaNew, isNew);

assert(isequal(size(X0), [8 7]) && isequal(size(U0), [4 7]), 'sizes');
% originals preserved exactly
origCols = find(~isNew);
assert(max(abs(X0(:, origCols) - X), [], 'all') < 1e-15, 'original X preserved');
assert(max(abs(U0(:, origCols) - U), [], 'all') < 1e-15, 'original U preserved');
% inserted throttle holds the LEFT (pre-switch, burn) value = 1, never averaged to 0.5
ins = find(isNew);
assert(abs(U0(4, ins) - 1) < 1e-15, 'insert throttle must step-hold left value 1, got %g', U0(4, ins));
% inserted direction is unit-norm
assert(abs(norm(U0(1:3, ins)) - 1) < 1e-12, 'insert direction must be unit norm');
% inserted state is between its neighbors (pchip monotone here)
assert(X0(1, ins) > X(1,3) && X0(1, ins) < X(1,4), 'insert state interpolated');
fprintf('ALL PASS\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/refine
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('$(pwd)'); test_warmstart_on_mesh" 2>&1 | grep -v -i "home license\|personal use\|academic, research\|organizational use"
```
Expected: FAIL with `Unrecognized function or variable 'warmstart_on_mesh'`.

- [ ] **Step 3: Write minimal implementation**

Create `warmstart_on_mesh.m`:

```matlab
function [X0, U0] = warmstart_on_mesh(out, sigma, sigmaNew, isNew)
% WARMSTART_ON_MESH  Build a no-resample warm start on a refined sigma mesh.
%
% Every original node's state/control is copied VERBATIM (the no-resample
% discipline); only inserted nodes are filled. States and thrust direction
% are pchip-interpolated (direction renormalized to unit norm); the throttle
% is STEP-held from the nearest original node to the LEFT (the pre-switch
% side), so an insert straddling a switch is never seeded with a smeared
% intermediate throttle -- the re-solve relocates the switch.
%
% INPUTS:
%   out      - struct with X [8xnN], U [4xnN] on the OLD mesh
%   sigma    - old normalized nodes [nNx1]
%   sigmaNew - refined normalized nodes [nN'x1], contains all of sigma
%   isNew    - logical mask of inserted nodes [nN'x1]
%
% OUTPUTS:
%   X0 - warm-start states on the refined mesh [8xnN']
%   U0 - warm-start controls on the refined mesh [4xnN']
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-11-pmp-mesh-refine-design.md
%   [2] OPTIMALITY_VERIFICATION_PLAN.md sec F.3 (no-resample requirement)

sigma = sigma(:);  sigmaNew = sigmaNew(:);  isNew = logical(isNew(:));
X = out.X;  U = out.U;
nNn = numel(sigmaNew);
X0 = zeros(8, nNn);  U0 = zeros(4, nNn);

% originals verbatim
X0(:, ~isNew) = X;
U0(:, ~isNew) = U;

ins = find(isNew).';
if ~isempty(ins)
    sv = sigmaNew(ins);
    % states + direction by pchip on the OLD (sigma, .) grid
    X0(:, ins)     = interp1(sigma, X.',        sv, 'pchip').';
    al             = interp1(sigma, U(1:3, :).', sv, 'pchip').';
    al             = al ./ sqrt(sum(al.^2, 1));
    U0(1:3, ins)   = al;
    % throttle: step-hold from the nearest original node to the left
    for q = ins
        kk = find(sigma <= sigmaNew(q), 1, 'last');
        kk = min(max(kk, 1), numel(sigma));
        U0(4, q) = U(4, kk);
    end
end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/refine
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('$(pwd)'); test_warmstart_on_mesh" 2>&1 | grep -v -i "home license\|personal use\|academic, research\|organizational use"
```
Expected: `ALL PASS`.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/warmstart_on_mesh.m NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/test_warmstart_on_mesh.m
git commit -m "refine: no-resample warm-start builder (step-throttle across switches)"
```

---

### Task 4: Seed preparation (regenerate duals)

**Files:**
- Create: `NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/prep_refine_seed.m`
- Test: `NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/test_prep_refine_seed.m`

**Interfaces:**
- Consumes: `sundman_minfuel/casadi_minfuel_sundman(sigma, tf, rv0, rvf, Tmax, c, muStar, X0, U0, tauf0, pSund, maxIter, epsilon, warmTight) → out` (with `out.lamDef [8×N]`); `sundman_minfuel/cr3bp_lt_params`.
- Produces: `outFile = prep_refine_seed(seedFile, outFile)` — reads a solution `.mat`, ensures it carries `out.lamDef` (regenerating via one ε=0 warmTight re-solve if absent), stamps `factor`, and writes `{out, factor, tauf0, sigma, rv0, rvf}` in the layout `sms_seed_duals`/`pmp_refine_indicator` expect. Returns the written path.

- [ ] **Step 1: Write the failing test**

Create `test_prep_refine_seed.m`:

```matlab
function test_prep_refine_seed()
% TEST_PREP_REFINE_SEED  Prepared 1.15x seed carries duals + required fields.
%
% NOTE: runs one eps=0 re-solve (~1-3 min).
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));  addpath(here);
addpath(fullfile(here, '..', '..', 'ms_band'));   % ms_band/setup_paths does NOT add ms_band itself
old = cd(fullfile(here, '..', '..', 'ms_band'));  c = onCleanup(@() cd(old));
setup_paths();  cd(old);

src = fullfile(here, '..', 'sundman_minfuel_certified.mat');
out = fullfile(tempdir, 'refine_seed_test.mat');
outFile = prep_refine_seed(src, out);

S = load(outFile);
assert(isfield(S, 'out') && isfield(S.out, 'lamDef') && ~isempty(S.out.lamDef), 'lamDef present');
assert(isequal(size(S.out.lamDef, 1), 8), 'lamDef is 8xN');
assert(isfield(S, 'factor') && abs(S.factor - 1.15) < 1e-9, 'factor = 1.15');
assert(all(isfield(S, {'tauf0', 'sigma', 'rv0', 'rvf'})), 'required fields present');
assert(S.out.maxDefect < 1e-6, 'seed re-solve converged tight, got %.2e', S.out.maxDefect);
fprintf('ALL PASS (switches=%d, defect=%.2e)\n', S.out.switches, S.out.maxDefect);
end
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/refine
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('$(pwd)'); test_prep_refine_seed" 2>&1 | grep -v -i "home license\|personal use\|academic, research\|organizational use"
```
Expected: FAIL with `Unrecognized function or variable 'prep_refine_seed'`.

- [ ] **Step 3: Write minimal implementation**

Create `prep_refine_seed.m`:

```matlab
function outFile = prep_refine_seed(seedFile, outFile)
% PREP_REFINE_SEED  Ensure a direct solution carries duals for refinement.
%
% The certified 1.15x .mat was saved before dual extraction and lacks
% out.lamDef and a factor field. This loads such a file, re-solves eps=0
% warmTight from its own (X,U) to regenerate the KKT-dual costates
% (out.lamDef), stamps factor = round(tf/tfMin, 2), and writes the layout
% pmp_refine_indicator / sms_seed_duals require. A file that already carries
% out.lamDef is passed through with fields normalized (no re-solve).
%
% INPUTS:
%   seedFile - source solution .mat (certified layout: out, sigma, tauf0,
%              rv0, rvf, tf; out may lack lamDef/factor)
%   outFile  - destination .mat path
%
% OUTPUTS:
%   outFile - the written path [char] (echoes the input)
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-11-pmp-mesh-refine-design.md
%   [2] sundman_minfuel/casadi_minfuel_sundman.m (solver + dual extraction)

S   = load(seedFile);
p   = cr3bp_lt_params(0.025, 15, 2100);
out = S.out;
tf  = S.out.X(8, end);                       % carried terminal time = tf
tfMin = 6.290694;                            % campaign constant (ND)
factor = round(tf/tfMin, 2);

if ~isfield(out, 'lamDef') || isempty(out.lamDef)
    fprintf('prep_refine_seed: regenerating duals (eps=0 warmTight re-solve)...\n');
    out = casadi_minfuel_sundman(S.sigma, tf, S.rv0, S.rvf, p.Tmax, p.c, ...
              p.muStar, S.out.X, S.out.U, S.tauf0, 1.5, 3000, 0, true);
    assert(out.success && out.maxDefect < 1e-6 && ~isempty(out.lamDef), ...
           'seed re-solve failed: success=%d defect=%.2e', out.success, out.maxDefect);
end

sigma = S.sigma;  tauf0 = S.tauf0;  rv0 = S.rv0;  rvf = S.rvf; %#ok<NASGU>
save(outFile, 'out', 'factor', 'tauf0', 'sigma', 'rv0', 'rvf');
fprintf('prep_refine_seed: wrote %s (factor=%.2f, switches=%d)\n', ...
        outFile, factor, out.switches);
end
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/refine
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('$(pwd)'); test_prep_refine_seed" 2>&1 | grep -v -i "home license\|personal use\|academic, research\|organizational use"
```
Expected: `ALL PASS (switches=25, defect=...e-14)` (or the certified count; a re-solve may report the near-graze-adjudicated 24–25). Delete any stray `fort.6`.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
rm -f NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/fort.6
git add NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/prep_refine_seed.m NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/test_prep_refine_seed.m
git commit -m "refine: seed-prep helper (regenerate KKT duals via eps=0 re-solve)"
```

---

### Task 5: Refinement-loop driver + smoke run

**Files:**
- Create: `NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/refine_loop.m`
- Test: `NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/test_refine_loop_smoke.m`

**Interfaces:**
- Consumes: `pmp_refine_indicator`, `refine_sigma`, `warmstart_on_mesh` (Tasks 1–3), `casadi_minfuel_sundman`, `cr3bp_lt_params`.
- Produces: `history = refine_loop(seedFile, opts)` — `history` is a `[1×R]` struct array (row per measured round, round 1 = the seed) with fields `nNodes`, `switches`, `tauSwitch`, `maxSwitchMove` (τ move vs previous round; `NaN` at round 1), `prop_kg`, `dProp` (vs previous; `NaN` at round 1), `nViol`, `HresMax`, `maxDefect`, `betaSpread`, `converged` (logical), `ipoptStatus`. Also saves `refine_history_<tag>.mat` and `refine_<tag>.png`. `opts` fields: `maxRounds` [default 4], `tag` [default from seed basename], plus indicator/refiner opts (`M`, `epsEval`, `mode`, `nbr`, `K`, `hFloor`, `maxAdd`), `propTol` [default 1e-4 kg].

- [ ] **Step 1: Write the failing test**

Create `test_refine_loop_smoke.m`:

```matlab
function test_refine_loop_smoke()
% TEST_REFINE_LOOP_SMOKE  2-round loop on the fast 1.12x file runs & records.
%
% NOTE: runs up to 2 eps=0 re-solves of the 10-switch 1.12x problem (~min each).
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));  addpath(here);
addpath(fullfile(here, '..', '..', 'ms_band'));   % ms_band/setup_paths does NOT add ms_band itself
old = cd(fullfile(here, '..', '..', 'ms_band'));  c = onCleanup(@() cd(old));
setup_paths();  cd(old);

seed = fullfile(here, '..', 'results', 'minfuel', 'legacy_ms_f1120.mat');
opts = struct('maxRounds', 2, 'tag', 'smoke_1p12', 'K', 6, 'maxAdd', 30);
history = refine_loop(seed, opts);

assert(numel(history) >= 2, 'expect >= 2 measured rounds, got %d', numel(history));
assert(isnan(history(1).maxSwitchMove), 'round 1 has no previous move');
assert(all([history.nNodes] == sort([history.nNodes])), 'nodes non-decreasing');
assert(history(end).maxDefect < 1e-6 || ~history(end).converged, 'tight or flagged');
assert(isfile(fullfile(here, 'refine_smoke_1p12.png')), 'figure written');
fprintf('ALL PASS (rounds=%d, final switches=%d, maxMove=%.2e)\n', ...
        numel(history), history(end).switches, history(end).maxSwitchMove);
end
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/refine
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('$(pwd)'); test_refine_loop_smoke" 2>&1 | grep -v -i "home license\|personal use\|academic, research\|organizational use"
```
Expected: FAIL with `Unrecognized function or variable 'refine_loop'`.

- [ ] **Step 3: Write minimal implementation**

Create `refine_loop.m`:

```matlab
function history = refine_loop(seedFile, opts)
% REFINE_LOOP  PMP-residual-driven adaptive mesh refinement (prototype driver).
%
% Each round: measure the current solution's PMP switch-localization score
% (pmp_refine_indicator), refine the sigma mesh where switches are worst
% localized (refine_sigma), build a no-resample warm start (warmstart_on_mesh),
% and re-solve the direct Sundman solver at eps=0 warmTight. Stops when the
% switch times stabilize (max move < local interval width AND |dProp| < propTol
% AND switch count unchanged) or at maxRounds, or if a re-solve fails to
% converge tight. History is persisted every round (crash-recoverable) and a
% summary figure is written. No LM/shooting anywhere -- the indirect machinery
% (via pmp_refine_indicator) is a measurement tool only.
%
% INPUTS:
%   seedFile - prepared seed .mat carrying out.lamDef, factor, tauf0, sigma,
%              rv0, rvf (see prep_refine_seed)
%   opts     - struct: maxRounds [default 4], tag [default seed basename],
%              propTol [default 1e-4 kg], and pass-through indicator/refiner
%              opts M/epsEval/mode/nbr/K/hFloor/maxAdd
%
% OUTPUTS:
%   history - [1xR] struct array, row per measured round (row 1 = seed):
%             nNodes, switches, tauSwitch, maxSwitchMove, prop_kg, dProp,
%             nViol, HresMax, maxDefect, betaSpread, converged, ipoptStatus.
%             Saved to refine_history_<tag>.mat; figure refine_<tag>.png.
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-11-pmp-mesh-refine-design.md

here = fileparts(mfilename('fullpath'));
if nargin < 2, opts = struct(); end
[~, base] = fileparts(seedFile);
if ~isfield(opts, 'maxRounds'), opts.maxRounds = 4;      end
if ~isfield(opts, 'tag'),       opts.tag = base;         end
if ~isfield(opts, 'propTol'),   opts.propTol = 1e-4;     end

p = cr3bp_lt_params(0.025, 15, 2100);
S = load(seedFile);
out = S.out;  sigma = S.sigma;  tauf0 = S.tauf0;
rv0 = S.rv0;  rvf = S.rvf;  factor = S.factor;
tf  = out.X(8, end);
tmpFile = fullfile(here, sprintf('.refine_tmp_%s.mat', opts.tag));

history = struct([]);  prevSwitch = [];  prevProp = NaN;
for r = 1:(opts.maxRounds + 1)
    % --- measure current solution ---
    write_seed(tmpFile, out, factor, tauf0, sigma, rv0, rvf);
    [score, tauSwitch, diag] = pmp_refine_indicator(tmpFile, opts);
    prop = p.m0kg*(1 - out.mf);

    [maxMove, converged, dProp] = deal(NaN, false, NaN);
    if ~isempty(prevSwitch)
        [maxMove, localH] = switch_move(prevSwitch, tauSwitch, diag.tauN);
        dProp = prop - prevProp;
        converged = numel(tauSwitch) == numel(prevSwitch) ...
                    && maxMove < localH && abs(dProp) < opts.propTol;
    end
    history(r).nNodes = numel(sigma); %#ok<*AGROW>
    history(r).switches = numel(tauSwitch);
    history(r).tauSwitch = tauSwitch;
    history(r).maxSwitchMove = maxMove;
    history(r).prop_kg = prop;
    history(r).dProp = dProp;
    history(r).nViol = diag.nViol;
    history(r).HresMax = diag.HresMax;
    history(r).maxDefect = out.maxDefect;
    history(r).betaSpread = diag.betaSpread;
    history(r).converged = converged;
    history(r).ipoptStatus = out.ipoptStatus;
    save(fullfile(here, sprintf('refine_history_%s.mat', opts.tag)), 'history');
    fprintf(['[round %d] nodes=%d sw=%d maxMove=%.2e dProp=%.2e nViol=%d ' ...
             'HresMax=%.2e defect=%.2e conv=%d\n'], r-1, numel(sigma), ...
            numel(tauSwitch), maxMove, dProp, diag.nViol, diag.HresMax, ...
            out.maxDefect, converged);
    if converged || r > opts.maxRounds, break; end

    % --- refine + re-solve ---
    [sigmaNew, isNew, nDropped] = refine_sigma(sigma, score, opts);
    if nDropped > 0
        fprintf('  refine_sigma dropped %d sub-hFloor interval(s)\n', nDropped);
    end
    if nnz(isNew) == 0
        fprintf('  no intervals refinable (all sub-hFloor); stopping.\n');  break;
    end
    [X0, U0] = warmstart_on_mesh(out, sigma, sigmaNew, isNew);
    o = casadi_minfuel_sundman(sigmaNew, tf, rv0, rvf, p.Tmax, p.c, p.muStar, ...
                               X0, U0, tauf0, 1.5, 3000, 0, true);
    if ~(o.success && o.maxDefect < 1e-6)
        fprintf('  re-solve did NOT converge tight (defect=%.2e, %s); stopping.\n', ...
                o.maxDefect, o.ipoptStatus);
        break;
    end
    prevSwitch = tauSwitch;  prevProp = prop;
    sigma = sigmaNew;  out = o;
end
if isfile(tmpFile), delete(tmpFile); end

make_figure(history, diag, here, opts.tag);
end

% -------------------------------------------------------------------------
function write_seed(f, out, factor, tauf0, sigma, rv0, rvf)
% Persist the current solution in sms_seed_duals input layout.
save(f, 'out', 'factor', 'tauf0', 'sigma', 'rv0', 'rvf');
end

function [maxMove, localH] = switch_move(prev, curr, tauN)
% Nearest-neighbor match of switch times between rounds; max move + the
% local interval width at the worst-moving switch (acceptance scale).
maxMove = 0;  localH = Inf;
n = min(numel(prev), numel(curr));
if numel(prev) ~= numel(curr)
    maxMove = Inf;  localH = 1;  return;   % count changed -> not converged
end
sp = sort(prev);  sc = sort(curr);
mv = abs(sp - sc);
[maxMove, w] = max(mv);
% local width near this switch time on the CURRENT node grid
[~, kk] = min(abs(tauN - sc(w)));
kk = min(max(kk, 1), numel(tauN) - 1);
localH = tauN(kk+1) - tauN(kk);
end

function make_figure(history, diag, here, tag)
% Three-panel summary: switch times vs round, S(tau)+throttle, escalation row.
fig = figure('Visible', 'off', 'Position', [50 50 1000 780]);
subplot(3, 1, 1);
hold on;
for r = 1:numel(history)
    ts = history(r).tauSwitch;
    plot(r*ones(size(ts)), ts, 'o');
end
grid on; xlabel('round'); ylabel('switch time \tau');
title(sprintf('%s: switch times vs refinement round', tag), 'Interpreter', 'none');
subplot(3, 1, 2);
yyaxis left;  plot(diag.tauN, diag.Snode, 'LineWidth', 1); ylabel('S (duals)');
hold on; yline(0, 'k-');
plot(diag.tauCr, zeros(size(diag.tauCr)), 'go');
yyaxis right; ylabel('|H_\sigma|'); plot(diag.tauN, diag.Hres, '-');
xlabel('\tau'); grid on; title('final: switching function, S=0 crossings, |H_\sigma|');
subplot(3, 1, 3);
rr = 0:(numel(history)-1);
yyaxis left;  plot(rr, [history.nViol], 's-'); ylabel('nViol');
yyaxis right; plot(rr, [history.HresMax], 'o-'); ylabel('HresMax');
xlabel('round'); grid on; title('escalation dashboard (drive option-2 decision)');
exportgraphics(fig, fullfile(here, sprintf('refine_%s.png', tag)), 'Resolution', 140);
close(fig);
end
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/refine
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('$(pwd)'); test_refine_loop_smoke" 2>&1 | grep -v -i "home license\|personal use\|academic, research\|organizational use"
```
Expected: `ALL PASS (rounds=..., final switches=..., maxMove=...)`. Delete any stray `fort.6`.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
rm -f NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/fort.6
git add NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/refine_loop.m NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/test_refine_loop_smoke.m
git commit -m "refine: adaptive mesh-refinement loop driver + 1.12x smoke test"
```

---

### Task 6: Headline 1.15× run + results note

**Files:**
- Create: `NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/run_headline_1p15.m`
- Create: `NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/RESULTS.md`
- Create: `NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/README.md`

**Interfaces:**
- Consumes: `prep_refine_seed`, `refine_loop` (Tasks 4–5).
- Produces: a saved `refine_history_headline_1p15.mat`, `refine_headline_1p15.png`, and a written `RESULTS.md` reporting whether switch times stabilized, the propellant drift, and the passive Hamiltonian-residual verdict on whether option 2 is warranted.

- [ ] **Step 1: Write the headline runner**

Create `run_headline_1p15.m`:

```matlab
function history = run_headline_1p15()
% RUN_HEADLINE_1P15  Prototype demonstration on the certified 1.15x solution.
%
% Prepares a duals-carrying seed from sundman_minfuel_certified.mat, then runs
% the refinement loop and prints the summary table for RESULTS.md.
%
% INPUTS:  none
% OUTPUTS: history - the refine_loop history struct array
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-11-pmp-mesh-refine-design.md

here = fileparts(mfilename('fullpath'));  addpath(here);
addpath(fullfile(here, '..', '..', 'ms_band'));   % ms_band/setup_paths does NOT add ms_band itself
old = cd(fullfile(here, '..', '..', 'ms_band'));  c = onCleanup(@() cd(old));
setup_paths();  cd(old);

src  = fullfile(here, '..', 'sundman_minfuel_certified.mat');
seed = fullfile(here, 'seed_1p15.mat');
if ~isfile(seed), prep_refine_seed(src, seed); end

opts = struct('maxRounds', 4, 'tag', 'headline_1p15', 'K', 8, 'maxAdd', 40);
history = refine_loop(seed, opts);

fprintf('\n=== HEADLINE 1.15x SUMMARY ===\n');
fprintf('%-6s %-7s %-4s %-11s %-11s %-7s %-11s\n', ...
        'round', 'nodes', 'sw', 'maxMove', 'dProp(kg)', 'nViol', 'HresMax');
for r = 1:numel(history)
    h = history(r);
    fprintf('%-6d %-7d %-4d %-11.2e %-11.2e %-7d %-11.2e\n', ...
            r-1, h.nNodes, h.switches, h.maxSwitchMove, h.dProp, h.nViol, h.HresMax);
end
end
```

- [ ] **Step 2: Run the headline demonstration**

Run:
```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/refine
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('$(pwd)'); run_headline_1p15" 2>&1 | grep -v -i "home license\|personal use\|academic, research\|organizational use" | tee headline_run.log
```
Expected: a per-round table; switch `maxMove` decreasing toward/below the round's local h, `dProp` within ~1e-4 kg, `nViol` non-increasing. This is the empirical result — record the actual numbers (do not fabricate).

- [ ] **Step 3: Write RESULTS.md from the actual run**

Create `RESULTS.md` capturing the observed per-round table (copied from `headline_run.log`), a one-paragraph verdict on whether switch times stabilized at fixed propellant, and the **option-2 decision**: did `HresMax` stay high in un-refined arcs (→ escalate to Hamiltonian-driven refinement) or drop with the switch fixes (→ option 1 sufficed)? Reference the figure `refine_headline_1p15.png`. Report honestly, including a negative result if switches moved without propellant/violation improvement or if a re-solve failed.

- [ ] **Step 4: Write README.md for the folder**

Create `README.md`: one-paragraph purpose (point-4 prototype), the file map (`pmp_refine_indicator`, `refine_sigma`, `warmstart_on_mesh`, `prep_refine_seed`, `refine_loop`, `run_headline_1p15`), how to run (the two commands above), the seed options (default 1.15× via `prep_refine_seed`; fast alternative `results/minfuel/legacy_ms_f1120.mat` at 1.12×), and a pointer to the spec and to `LOW_THRUST_MINFUEL_CAMPAIGN.md`.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
rm -f NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/fort.6
git add NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/run_headline_1p15.m NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/RESULTS.md NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/README.md NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/refine_headline_1p15.png NLP_lowThrust_GTO_tulip/sundman_minfuel/refine/refine_history_headline_1p15.mat
git commit -m "refine: headline 1.15x demonstration + RESULTS/README"
```

---

## Notes for the implementer

- **Run tasks in order.** Tasks 1–3 are pure/fast (seconds). Task 4 runs one ~1–3 min re-solve. Tasks 5–6 run 1–4 re-solves each (~min apiece for 1.12×, longer for the 4001-node 1.15×). Budget Bash timeouts to 600000 ms and prefer `run_in_background` for the headline run.
- **MEX-crash recoverability:** `refine_loop` saves `refine_history_<tag>.mat` every round before re-solving. If IPOPT MEX-crashes the MATLAB process (documented, uncatchable), the completed rounds are on disk — rerun and note it in RESULTS.md; do not treat a crash as success.
- **Do not resample.** If tempted to interpolate the whole trajectory onto a fresh uniform mesh, stop — that is the documented restoration-pinning failure. Only inserted nodes get interpolated.
- **`fort.6`** appears in the cwd after any IPOPT run — the commit steps delete it; never stage it.
- **Negative results are results.** If refinement does not stabilize switches, or propellant drifts past `propTol`, or `nViol` rises, RESULTS.md records that plainly.
