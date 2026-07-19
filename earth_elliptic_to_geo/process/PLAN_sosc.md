# SOSC Certificate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an NLP-level second-order-sufficient-conditions (SOSC) local-minimum certificate — independent KKT re-check + reduced-Hessian critical-cone inertia test — as the new bar for "certified", applied to reproducer and existing campaign rows.

**Architecture:** New `verify/sosc/` subsystem. For a saved certified fuel row, rebuild the exact CasADi Opti NLP from the saved primal + config, warm re-solve to recover the full multiplier set, assemble the Lagrangian Hessian `H` and constraint Jacobian `A` from Opti's native (unscaled) symbols, run a KKT residual re-check, classify the active set (strong/weak), and test the inertia of `[[H,Aᵀ],[A,0]]`. Verdict ∈ {PASS, FAIL, INCONCLUSIVE, ERROR} drives a tiered gate. One numerics-preserving `returnModel` hook is added to `core/casadi_lt_mee.m`; nothing else on the certified path changes.

**Tech Stack:** MATLAB R2025b, CasADi 3.7.0 (`~/casadi-3.7.0`) + IPOPT/MUMPS, sparse `ldl` for inertia.

**Spec:** `process/DESIGN_sosc.md` (authoritative — read it; this plan implements it).

## Global Constraints

- **MATLAB R2025b only** — invoke `/Applications/MATLAB_R2025b.app/bin/matlab`; R2025a license is broken.
- **CasADi on path:** `addpath(fullfile(getenv('HOME'),'casadi-3.7.0'))`; every test cd's to the module root and calls `setup_paths` first.
- **Work exclusively in Opti's native (unscaled) symbols** — `opti.x`, `opti.g`, `opti.f`, `opti.lam_g`. Never mix in IPOPT-internal scaled quantities. This is what sidesteps the Cartesian cScale dual-anomaly class.
- **`returnModel` is numerics-invariant** — with the flag off (all existing callers), `casadi_lt_mee` output is byte-for-byte unchanged.
- **Never clobber campaign caches** — existing-row verdicts go to sidecar `results/sosc/sosc_<tag>.mat`; `results/*.mat` campaign files are read-only here.
- **Verdict logic lives in the orchestrator** (`verify_sosc_mee`), not in `sosc_inertia` — `sosc_inertia` returns inertia counts + a subspace-OK bool only.
- **Canonical units** (`kepler_lt_params`): μ=1, LU=42165 km, mass unit=m0; residual magnitudes are O(1), so absolute tolerances in `sosc_defaults.m` are meaningful.
- **MATLAB function-header block required** on every new function (purpose / INPUTS with sizes / OUTPUTS with sizes / REFERENCES), per repo CLAUDE.md.
- **Reconstruction is the gate** — no certification is trusted until Task 1 proves the rebuilt NLP reproduces the saved primal to `tol.recon`.

---

## File Structure

| File | Responsibility |
|---|---|
| `verify/sosc/sosc_load_row.m` | Normalize a saved MEE_M2 `res` or PSR `out` into a common `saved` struct. |
| `verify/sosc/sosc_defaults.m` | Single source of the tolerance struct. |
| `verify/sosc/sosc_recover_kkt.m` | Warm re-solve at saved primal; return `x*`, `lam_g`, sparse `H`, sparse `A_all`, `grad_f`, `gval`, registries, drift. |
| `verify/sosc/sosc_kkt_residual.m` | Global-sign resolution + KKT residuals (stationarity/feas/dual/comp). |
| `verify/sosc/sosc_active_set.m` | Active/strong/weak classification; active Jacobian `A`; LICQ flag; weak-node labels. |
| `verify/sosc/sosc_inertia.m` | Sparse LDLᵀ inertia of the KKT matrix + subspace-SOSC bool + optional curvature margin. |
| `verify/sosc/verify_sosc_mee.m` | Orchestrator: verdict + tiered-gate status struct. |
| `verify/sosc/recertify_table3.m` | Batch re-cert of existing rows → sidecar + printed report. |
| `core/casadi_lt_mee.m` | **Modify:** add `opts.returnModel` (registry + model handles). Constraint block at lines 102–172. |
| `drivers/run_transfer_mee.m` | **Modify:** attach `res.sosc` + apply tiered gate (near line 250–258). |
| `reproduce/reproduce_row.m` | **Modify:** adopt keep-best candidate only if verdict ∈ {PASS, INCONCLUSIVE}. |
| `tests/test_sosc_*.m` | Test scripts (repo convention: cd root, `setup_paths`, assert, print PASSED). |

---

## Task 1: Row normalizer + reconstruction checkpoint (the gate)

**Files:**
- Create: `verify/sosc/sosc_load_row.m`
- Test: `tests/test_sosc_recon_10N.m`

**Interfaces:**
- Produces: `saved = sosc_load_row(matPath)` → struct with fields
  `sigma [(N+1)x1]`, `X [7x(N+1)]`, `U [4x(N+1)]`, `dL [scalar]`,
  `tfTarget [scalar]`, `xf [5x1]`, `thrustN`, `m0kg`, `ispS`, `maxIter`,
  `tag [char]`, `kind ['MEE_M2'|'PSR']`.
- Produces (for Task 4): the exact reconstruction recipe — call
  `casadi_lt_mee(saved.sigma, saved.X, saved.U, saved.dL, opts)` with
  `opts = struct('par', kepler_lt_params(saved.thrustN,saved.m0kg,saved.ispS), 'mode','fixedtf','eps',0, 'tfTarget',saved.tfTarget, 'x0',saved.X(:,1), 'xf',saved.xf, 'maxIter',saved.maxIter, 'warmTight',true, 'printLevel',0)`.

**Context:** MEE_M2 rows save a `res` struct (`run_transfer_mee.m:255-258`) with `res.sigma`, `res.fuel` (a `casadi_lt_mee` out struct: `.X/.U/.dL`), `res.tf` (the resolved fixed-tf target — use this directly as `tfTarget`), `res.fp` (`.thrustN/.m0kg/.ispS/.maxIter`; **`.xf` is present only on rows regenerated after the xf field was added — the older 5/2.5/1/0.5 N caches lack it**). PSR-final rows save `out` (`psr_mee_refine.m:294-305`) with `out.finalSigma`, `out.finalOut` (out struct), and `fpFinal` (`.baseTag/.thrustN/.m0kg/.ispS/.tf`; there `tfTarget = fpFinal.tf` — already resolved, no `ctf*tfMinAnchor` multiply; **`.xf` is absent from `fpFinal` and from every field of the PSR file**). **`xf` must therefore default to the GEO target `[1;0;0;0;0]` when absent** — correct for every campaign row (all target GEO; custom endpoints are a never-used research feature), via `optdef(fp,'xf',[1;0;0;0;0])`.

- [ ] **Step 1: Write the failing test** `tests/test_sosc_recon_10N.m`

```matlab
% TEST_SOSC_RECON_10N  Task-1 gate: sosc_load_row normalizes the 10 N MEE_M2
% row, and rebuilding+re-solving the NLP from it reproduces the saved primal
% to tol.recon. If this fails, res.fp is insufficient and must be fixed
% before any SOSC work is trusted.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));

matPath = fullfile(module_root(),'results','MEE_M2_10N.mat');
assert(isfile(matPath), 'need the certified 10 N cache MEE_M2_10N.mat');
saved = sosc_load_row(matPath);

assert(isequal(size(saved.X),[7, numel(saved.sigma)]), 'X shape');
assert(saved.thrustN==10 && saved.xf(1)==1, 'thrustN/xf');
assert(saved.tfTarget>0, 'tfTarget resolved');

par  = kepler_lt_params(saved.thrustN, saved.m0kg, saved.ispS);
opts = struct('par',par,'mode','fixedtf','eps',0,'tfTarget',saved.tfTarget, ...
    'x0',saved.X(:,1),'xf',saved.xf,'maxIter',saved.maxIter, ...
    'warmTight',true,'printLevel',0);
o = casadi_lt_mee(saved.sigma, saved.X, saved.U, saved.dL, opts);
drift = max(abs(o.X(:) - saved.X(:)));
fprintf('recon drift ||x_rebuilt - x_saved||_inf = %.3e\n', drift);
assert(o.success, 'rebuild re-solve did not converge');
assert(drift < 1e-6, sprintf('recon drift %.3e >= 1e-6 tol.recon', drift));

% Normalizer correctness on the real xf-less MEE_M2 row (5 N) and a PSR-final
% row -- both lack fp.xf; expect the GEO default. No re-solve (cheap checks).
s5 = sosc_load_row(fullfile(module_root(),'results','MEE_M2_5N.mat'));
assert(isequal(s5.xf,[1;0;0;0;0]), '5 N: xf must default to GEO when fp.xf absent');
assert(abs(s5.tfTarget - 67.0194) < 1e-3 && s5.thrustN==5, '5 N: tfTarget/thrustN');
sP = sosc_load_row(fullfile(module_root(),'results','MEE_M2_1N_PSR_psr_final.mat'));
assert(strcmp(sP.kind,'PSR') && sP.thrustN==1, 'PSR: kind/thrustN');
assert(isequal(sP.xf,[1;0;0;0;0]), 'PSR: xf defaults to GEO');
assert(abs(sP.tfTarget - 335.7122) < 1e-3, 'PSR: tfTarget = fpFinal.tf (already resolved)');
assert(isequal(size(sP.X),[7,numel(sP.sigma)]), 'PSR: X/sigma shapes consistent');
fprintf('test_sosc_recon_10N PASSED\n');
```

- [ ] **Step 2: Run it to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "run('$(pwd)/tests/test_sosc_recon_10N.m')"`
Expected: FAIL — `Unrecognized function 'sosc_load_row'`.

- [ ] **Step 3: Implement `verify/sosc/sosc_load_row.m`**

```matlab
function saved = sosc_load_row(matPath)
% SOSC_LOAD_ROW  Normalize a saved certified fuel row (MEE_M2 res-struct or PSR
% out-struct) into the common `saved` struct the SOSC recovery consumes.
%
% INPUTS:
%   matPath - path to a certified .mat (results/MEE_M2_*.mat or
%             results/*_PSR_psr_final.mat) [char]
% OUTPUTS:
%   saved - struct: .sigma[(N+1)x1] .X[7x(N+1)] .U[4x(N+1)] .dL[1] .tfTarget[1]
%           .xf[5x1] .thrustN .m0kg .ispS .maxIter .tag[char] .kind[char]
% REFERENCES:
%   [1] run_transfer_mee.m:255-258 (res-struct layout); psr_mee_refine.m:294-305
%       (PSR out layout); process/DESIGN_sosc.md sec 4.2.
S = load(matPath);
[~, base] = fileparts(matPath);
geoXf = [1;0;0;0;0];                             % GEO default when fp.xf absent
if isfield(S,'res')                              % MEE_M2 row
    r  = S.res;  fu = r.fuel;  fp = r.fp;
    saved = struct('sigma', r.sigma(:), 'X', fu.X, 'U', fu.U, 'dL', fu.dL, ...
        'tfTarget', r.tf, 'xf', optdef(fp,'xf',geoXf), 'thrustN', fp.thrustN, ...
        'm0kg', fp.m0kg, 'ispS', fp.ispS, 'maxIter', optdef(fp,'maxIter',1500), ...
        'tag', base, 'kind', 'MEE_M2');
elseif isfield(S,'out')                          % PSR-refined row
    o  = S.out;  fu = o.finalOut;  fp = S.fpFinal;
    saved = struct('sigma', o.finalSigma(:), 'X', fu.X, 'U', fu.U, 'dL', fu.dL, ...
        'tfTarget', fp.tf, 'xf', optdef(fp,'xf',geoXf), 'thrustN', fp.thrustN, ...
        'm0kg', fp.m0kg, 'ispS', fp.ispS, ...
        'maxIter', optdef(fp,'maxIter',1500), 'tag', base, 'kind', 'PSR');
else
    error('sosc_load_row:unknownShape', ...
        '%s has neither a res nor an out variable', matPath);
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "run('$(pwd)/tests/test_sosc_recon_10N.m')"`
Expected: PASS, `recon drift ... < 1e-6`.
**If drift ≥ tol:** `res.fp` is insufficient — STOP and escalate; do not proceed. This is the gate the spec (§10.1) flags.

- [ ] **Step 5: Commit**

```bash
git add verify/sosc/sosc_load_row.m tests/test_sosc_recon_10N.m
git commit -m "feat(sosc): row normalizer + reconstruction checkpoint (Task 1 gate)"
```

---

## Task 2: `returnModel` hook + constraint/variable registry (numerics-invariant)

**Files:**
- Modify: `core/casadi_lt_mee.m` (opts parse near line 74; constraint block 102–172; out struct 282–289)
- Test: `tests/test_sosc_returnmodel.m`

**Interfaces:**
- Produces: with `opts.returnModel=true`, `out.model = struct('opti',opti,'X',X,'U',U,'dL',dL,'creg',creg,'vreg',vreg)`.
  `creg`: struct array, one per constraint group, fields `label[char]`, `kind['eq'|'ineqLo'|'ineqHi']`, `rows[1xk]` (row range in `opti.g`), `bound[scalar|[]]`, `node[1xk|[]]`.
  `vreg`: `struct('Xrows',1:7,'Urows',1:4,'nNode',N+1)` — decision-block layout metadata. (The full per-index maps into `opti.x` are intentionally NOT built — YAGNI: no downstream task consumes them; `sosc_active_set`'s node labels come from `creg.node`, and `sosc_recover_kkt` assembles `H`/`A` from `opti.x`/`opti.g` wholesale.)

**Context:** Every constraint is added via `opti.subject_to` (lines 105,110,115,118,120,125-130,151-153,156,161,169). None are native `lbx/ubx`. Registry entries are recorded by bracketing each group: `r0 = size(opti.g,1)+1; <add>; r1 = size(opti.g,1); creg(end+1)=...`.

- [ ] **Step 1: Write the failing test** `tests/test_sosc_returnmodel.m`

```matlab
% TEST_SOSC_RETURNMODEL  (a) returnModel=false leaves the solve output
% unchanged (numerics invariant); (b) returnModel=true exposes opti + a
% registry whose row ranges partition opti.g exactly once.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));

saved = sosc_load_row(fullfile(module_root(),'results','MEE_M2_10N.mat'));
par = kepler_lt_params(saved.thrustN,saved.m0kg,saved.ispS);
base = struct('par',par,'mode','fixedtf','eps',0,'tfTarget',saved.tfTarget, ...
    'x0',saved.X(:,1),'xf',saved.xf,'maxIter',saved.maxIter,'warmTight',true,'printLevel',0);

oOff = casadi_lt_mee(saved.sigma,saved.X,saved.U,saved.dL, base);
optsOn = base; optsOn.returnModel = true;
oOn  = casadi_lt_mee(saved.sigma,saved.X,saved.U,saved.dL, optsOn);

% (a) numerics invariant
assert(isequal(oOff.X,oOn.X) && isequal(oOff.U,oOn.U) && oOff.dL==oOn.dL, ...
    'returnModel changed the numeric solution');
assert(~isfield(oOff,'model'), 'model leaked with flag off');

% (b) registry partitions opti.g
m = size(oOn.model.opti.g,1);
covered = [];
for i = 1:numel(oOn.model.creg), covered = [covered, oOn.model.creg(i).rows]; end %#ok<AGROW>
assert(isequal(sort(covered(:)'),1:m), 'creg rows must partition 1..m exactly once');
labels = {oOn.model.creg.label};
assert(any(strcmp(labels,'defect')) && any(strcmp(labels,'betaNorm')) && ...
       any(strcmp(labels,'termBC')), 'expected core labels present');
fprintf('test_sosc_returnmodel PASSED (m=%d constraint rows, %d groups)\n', m, numel(oOn.model.creg));
```

- [ ] **Step 2: Run it to verify it fails**

Expected: FAIL — `oOn.model` does not exist.

- [ ] **Step 3: Implement the hook.** In `casadi_lt_mee.m`:

(3a) After line 74, parse the flag:
```matlab
returnModel = d('returnModel', false);
creg = struct('label',{},'kind',{},'rows',{},'bound',{},'node',{});
addc = @(lab,kind,r0,bnd,nd) struct('label',lab,'kind',kind,'rows',r0:size(opti.g,1),'bound',bnd,'node',nd);
```
(Define `addc` AFTER `opti` is created at line 82; move the `addc`/`creg` init to just after line 82.)

(3b) Bracket each constraint group. Example for defects (replace loop at 102–106):
```matlab
r0 = size(opti.g,1)+1;
conDef = cell(1,N);
for k = 1:N
    conDef{k} = X(:,k+1) - X(:,k) - (dsig(k)/2)*dL*(dXdL(:,k)+dXdL(:,k+1)) == 0;
    opti.subject_to(conDef{k});
end
if returnModel, creg(end+1) = addc('defect','eq',r0,0,1:N); end
```
Apply the identical `r0 = size(opti.g,1)+1; ...; if returnModel, creg(end+1)=addc(<label>,<kind>,r0,<bound>,<node>); end` pattern to each remaining group with these labels/kinds/bounds:
- Ldot guard (109–111): `('ldotGuard','ineqLo', r0, par.LdotMin, 1:N+1)`
- beta norm (114–116): `('betaNorm','eq', r0, 1, 1:N+1)`
- throttle: mintime `('thrEq','eq',r0,1,1:N+1)`; fixedtf two groups `('thrLo','ineqLo',r0,0,1:N+1)` then `('thrHi','ineqHi',r0,1,1:N+1)`
- state boxes (125–130): one `ineqLo`+`ineqHi` pair per element, labels `boxP_lo/boxP_hi`, `boxEx_lo/…` etc., bounds as coded
- t box (151): `('tBox_lo','ineqLo',r0,0,1:N+1)`, `('tBox_hi','ineqHi',r0,tUB,1:N+1)`
- beta box (152): `('betaBox_lo','ineqLo',r0,-1.01,[])`, `('betaBox_hi','ineqHi',r0,1.01,[])`
- dL box (153): `('dLbox_lo','ineqLo',r0,0.1,[])`, `('dLbox_hi','ineqHi',r0,2000,[])`
- initial BC (156): `('initBC','eq',r0,0,1)`
- terminal BC (161): `('termBC','eq',r0,0,N+1)`
- fixedtf tf pin (169): `('tfPin','eq',r0,tfTarget,N+1)`

(3c) Build `vreg` and attach the model to `out` (after line 289, before `end`):
```matlab
if returnModel
    vreg = struct('Xrows',1:7,'Urows',1:4,'nNode',N+1);
    out.model = struct('opti',opti,'X',X,'U',U,'dL',dL,'creg',creg,'vreg',vreg);
end
```

- [ ] **Step 4: Run the test to verify it passes**

Expected: PASS, prints `m=... constraint rows, ... groups`.

- [ ] **Step 5: Regression guard.** Run the existing solver smoke to confirm no numeric drift:

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "run('$(pwd)/tests/test_mee_solver_smoke.m')"` (or the nearest existing solver smoke)
Expected: PASS unchanged.

- [ ] **Step 6: Commit**

```bash
git add core/casadi_lt_mee.m tests/test_sosc_returnmodel.m
git commit -m "feat(sosc): numerics-invariant returnModel hook + constraint registry (Task 2)"
```

---

## Task 3: Tolerance defaults

**Files:**
- Create: `verify/sosc/sosc_defaults.m`
- Test: `tests/test_sosc_defaults.m`

- [ ] **Step 1: Write the failing test** `tests/test_sosc_defaults.m`

```matlab
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
t = sosc_defaults();
for f = {'recon','drift','stat','feas','dual','comp','active','mu','inertiaZero'}
    assert(isfield(t,f{1}) && t.(f{1})>0, sprintf('missing/nonpos tol.%s',f{1}));
end
assert(t.feas==1e-8, 'feas must match the existing maxDefect<1e-8 gate');
fprintf('test_sosc_defaults PASSED\n');
```

- [ ] **Step 2: Run it — FAIL** (`sosc_defaults` undefined).

- [ ] **Step 3: Implement `verify/sosc/sosc_defaults.m`**

```matlab
function tol = sosc_defaults()
% SOSC_DEFAULTS  Single source of SOSC certificate tolerances (canonical units,
% magnitudes O(1)). See process/DESIGN_sosc.md sec 6.
% OUTPUTS: tol - struct of scalar thresholds.
tol = struct( ...
    'recon',       1e-6, ...  % rebuild reproduces saved primal
    'drift',       1e-6, ...  % warm-resolve drift (report/warn, not fail)
    'stat',        1e-6, ...  % stationarity ||grad L||_inf
    'feas',        1e-8, ...  % equality residual / inequality violation
    'dual',        1e-8, ...  % inequality dual-sign violation
    'comp',        1e-6, ...  % complementarity max|lam*slack|
    'active',      1e-7, ...  % inequality slack -> active
    'mu',          1e-6, ...  % relative multiplier -> strongly-active
    'inertiaZero', 1e-9);     % relative pivot magnitude -> zero eigenvalue
end
```

- [ ] **Step 4: Run — PASS.**
- [ ] **Step 5: Commit** `git add verify/sosc/sosc_defaults.m tests/test_sosc_defaults.m && git commit -m "feat(sosc): tolerance defaults (Task 3)"`

---

## Task 4: `sosc_recover_kkt` — warm re-solve + assemble H, A, duals

**Files:**
- Create: `verify/sosc/sosc_recover_kkt.m`
- Test: `tests/test_sosc_recover_10N.m`

**Interfaces:**
- Consumes: `sosc_load_row` output `saved`; `sosc_defaults`; the `returnModel` model.
- Produces: `R = sosc_recover_kkt(saved, tol)` → struct
  `.recoverOK[bool] .x[nx1] .lam_g[mx1] .gval[mx1] .grad_f[nx1] .H[nxn sparse] .A_all[mxn sparse] .creg .vreg .drift[scalar] .n .m .ipoptStatus`.

- [ ] **Step 1: Write the failing test** `tests/test_sosc_recover_10N.m`

```matlab
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));
saved = sosc_load_row(fullfile(module_root(),'results','MEE_M2_10N.mat'));
tol = sosc_defaults();
R = sosc_recover_kkt(saved, tol);
assert(R.recoverOK, 'recovery failed: %s', R.ipoptStatus);
assert(numel(R.x)==R.n && numel(R.lam_g)==R.m, 'x/lam_g dims');
assert(isequal(size(R.H),[R.n R.n]) && issparse(R.H), 'H shape/sparse');
assert(size(R.A_all,1)==R.m && size(R.A_all,2)==R.n && issparse(R.A_all), 'A_all shape/sparse');
assert(R.drift < tol.drift, sprintf('drift %.3e >= %.1e', R.drift, tol.drift));
fprintf('test_sosc_recover_10N PASSED (n=%d, m=%d, drift=%.2e)\n', R.n, R.m, R.drift);
```

- [ ] **Step 2: Run — FAIL** (`sosc_recover_kkt` undefined).

- [ ] **Step 3: Implement `verify/sosc/sosc_recover_kkt.m`**

```matlab
function R = sosc_recover_kkt(saved, tol)
% SOSC_RECOVER_KKT  Rebuild the NLP at a saved primal, warm re-solve to recover
% the full multiplier set, and assemble the KKT objects in Opti's native
% (unscaled) symbols.
%
% INPUTS:
%   saved - struct from sosc_load_row (primal + config)
%   tol   - struct from sosc_defaults
% OUTPUTS:
%   R - struct: .recoverOK .x[nx1] .lam_g[mx1] .gval[mx1] .grad_f[nx1]
%       .H[nxn sparse] .A_all[mxn sparse] .creg .vreg .drift .n .m .ipoptStatus
% REFERENCES:
%   [1] process/DESIGN_sosc.md sec 4.2. [2] casadi_lt_mee.m (returnModel hook).
import casadi.*
par  = kepler_lt_params(saved.thrustN, saved.m0kg, saved.ispS);
opts = struct('par',par,'mode','fixedtf','eps',0,'tfTarget',saved.tfTarget, ...
    'x0',saved.X(:,1),'xf',saved.xf,'maxIter',saved.maxIter, ...
    'warmTight',true,'printLevel',0,'returnModel',true);
o = casadi_lt_mee(saved.sigma, saved.X, saved.U, saved.dL, opts);
R.ipoptStatus = o.ipoptStatus;
R.recoverOK   = o.success && o.maxDefect < tol.feas;
R.drift = max(abs(o.X(:) - saved.X(:)));
if ~R.recoverOK, R.x=[]; R.lam_g=[]; R.gval=[]; R.grad_f=[]; R.H=[]; R.A_all=[];
    R.creg=[]; R.vreg=[]; R.n=0; R.m=0; return; end

opti = o.model.opti;  sol = opti.debug;   % the solved Opti (sol from last solve)
% Native symbols:
x   = opti.x;   g = opti.g;   f = opti.f;   lam = opti.lam_g;
R.x     = full(sol.value(x));
R.lam_g = full(sol.value(lam));
R.gval  = full(sol.value(g));
% Gradient, Jacobian, Lagrangian Hessian as CasADi Functions, evaluated at soln:
gradF = gradient(f, x);
Jg    = jacobian(g, x);
Hlag  = hessian(f + lam.'*g, x);           % returns [H, grad] -> take H
Fkkt  = Function('Fkkt', {x, lam}, {gradF, Jg, Hlag});
[gf, A, H] = Fkkt(R.x, R.lam_g);
R.grad_f = full(gf);
R.A_all  = sparse(A);
R.H      = sparse(H);
R.creg = o.model.creg;  R.vreg = o.model.vreg;
R.n = numel(R.x);  R.m = numel(R.gval);
end
```
Note: if `hessian(f+lam.'*g, x)` returns a 2-output `[H,g]` form in this CasADi version, capture only the first; the `Function` wrapper above already isolates `Hlag` as the Hessian expression. Confirm on first run and adjust the `hessian` call if the API differs.

- [ ] **Step 4: Run — PASS**, prints `n=..., m=..., drift=...`. Confirm the CasADi `hessian` API shape here; fix if the two-output form bites.

- [ ] **Step 5: Commit** `git add verify/sosc/sosc_recover_kkt.m tests/test_sosc_recover_10N.m && git commit -m "feat(sosc): warm-resolve KKT recovery (Task 4)"`

---

## Task 5: `sosc_kkt_residual` — sign resolution + residuals

**Files:**
- Create: `verify/sosc/sosc_kkt_residual.m`
- Test: `tests/test_sosc_kkt_residual.m`

**Interfaces:**
- Consumes: `R` from `sosc_recover_kkt`; `tol`.
- Produces: `K = sosc_kkt_residual(R, tol)` → struct `.sign[+1|-1] .signOK[bool] .stat .primalEq .primalIneq .dualFeas .comp .pass[bool]`.

- [ ] **Step 1: Write the failing test** `tests/test_sosc_kkt_residual.m`

```matlab
% Synthetic: n=2, one equality (row1, kind eq), one inequality g<=0 (row2,
% ineqHi bound 0). Choose grad_f, A_all, lam_g so stationarity is exactly 0
% under sign s=+1, and slack/comp are clean.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
R.n=2; R.m=2;
R.grad_f = [ -1; -1 ];
R.A_all  = sparse([1 0; 0 1]);        % d g1/dx, d g2/dx
R.lam_g  = [1; 1];                     % grad_f + A' * lam = 0  (s=+1)  -> lam forced nonzero
R.gval   = [0; 0];                     % both active: eq satisfied; ineq AT its bound
                                       % (a nonzero lam on row 2 requires slack 0 for
                                       %  complementarity -- an inactive row would violate it)
R.creg = struct('label',{'eqA','ineqB'},'kind',{'eq','ineqHi'}, ...
                'rows',{1,2},'bound',{0,0},'node',{[],[]});
K = sosc_kkt_residual(R, sosc_defaults());
assert(K.signOK && K.sign==1, 'sign should resolve to +1');
assert(K.stat < 1e-12, 'stationarity ~0'); assert(K.comp < 1e-9, 'comp ~0');
assert(K.dualFeas <= 0 || K.dualFeas < 1e-12, 'dual feasible');
assert(K.pass, 'overall KKT pass');
fprintf('test_sosc_kkt_residual PASSED\n');
```

- [ ] **Step 2: Run — FAIL.**

- [ ] **Step 3: Implement `verify/sosc/sosc_kkt_residual.m`**

```matlab
function K = sosc_kkt_residual(R, tol)
% SOSC_KKT_RESIDUAL  Resolve the one global Lagrangian-sign ambiguity, then
% re-check first-order KKT residuals at the recovered point.
%
% INPUTS: R - sosc_recover_kkt struct; tol - sosc_defaults struct
% OUTPUTS: K - struct .sign .signOK .stat .primalEq .primalIneq .dualFeas
%              .comp .pass
% REFERENCES: process/DESIGN_sosc.md sec 4.3; verify_pmp_mee.m:112-121 (sign trick).
isEq   = strcmp({R.creg.kind},'eq');
% expand per-group kind to per-row masks
kindRow = strings(R.m,1);
for i=1:numel(R.creg), kindRow(R.creg(i).rows) = R.creg(i).kind; end
eqRow   = kindRow=="eq";
ineqRow = ~eqRow;

% (1) global sign: choose s minimizing ||grad_f + s*A'*lam||_inf
rP = R.grad_f + (R.A_all.' * R.lam_g);
rM = R.grad_f - (R.A_all.' * R.lam_g);
if norm(rP,inf) <= norm(rM,inf), K.sign=+1; stat=norm(rP,inf);
else,                            K.sign=-1; stat=norm(rM,inf); end
K.stat   = stat;
K.signOK = stat < tol.stat;

% (2) primal feasibility (eq residual; ineq violation for g<=0 canonical form)
K.primalEq   = max(abs(R.gval(eqRow)), [], 'omitnan');
if isempty(K.primalEq), K.primalEq = 0; end
% ineqHi: g<=bound -> viol = max(0, g-bound); ineqLo: g>=bound -> viol=max(0,bound-g)
viol = zeros(R.m,1);
for i=1:numel(R.creg)
    c=R.creg(i); if strcmp(c.kind,'eq'), continue; end
    gv=R.gval(c.rows);
    if strcmp(c.kind,'ineqHi'), viol(c.rows)=max(0, gv - c.bound);
    else,                        viol(c.rows)=max(0, c.bound - gv); end
end
K.primalIneq = max(viol);

% (3) dual feasibility: s*lam_g >= 0 for all inequality rows (g<=0 convention)
lamSigned = K.sign * R.lam_g;
K.dualFeas = max(-lamSigned(ineqRow), [], 'omitnan');   % worst negative
if isempty(K.dualFeas), K.dualFeas = 0; end

% (4) complementarity: |lam * slack| over inequalities
slack = zeros(R.m,1);
for i=1:numel(R.creg)
    c=R.creg(i); if strcmp(c.kind,'eq'), continue; end
    gv=R.gval(c.rows);
    if strcmp(c.kind,'ineqHi'), slack(c.rows)=c.bound-gv; else, slack(c.rows)=gv-c.bound; end
end
K.comp = max(abs(R.lam_g(ineqRow).*slack(ineqRow)), [], 'omitnan');
if isempty(K.comp), K.comp = 0; end

K.pass = K.signOK && K.primalEq<tol.feas && K.primalIneq<tol.feas && ...
         K.dualFeas<tol.dual && K.comp<tol.comp;
end
```

- [ ] **Step 4: Run — PASS.**
- [ ] **Step 5: Commit** `git add verify/sosc/sosc_kkt_residual.m tests/test_sosc_kkt_residual.m && git commit -m "feat(sosc): KKT residual re-check + sign resolution (Task 5)"`

---

## Task 6: `sosc_active_set` — classification + active Jacobian + LICQ

**Files:**
- Create: `verify/sosc/sosc_active_set.m`
- Test: `tests/test_sosc_active_set.m`

**Interfaces:**
- Consumes: `R`, `K`, `tol`.
- Produces: `AS = sosc_active_set(R, K, tol)` → struct `.A[sparse] .m_active .nEq .nStrong .nWeak .weakLabels{cell} .licq[bool]`.

- [ ] **Step 1: Write the failing test** `tests/test_sosc_active_set.m`

```matlab
% n=3 vars; rows: 1 eq (active), 2 ineq. Ineq A active+strong (slack 0, |lam|
% big); ineq B active+WEAK (slack 0, |lam|~0). Expect nWeak==1, A has eq+strong.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
R.n=3; R.m=3;
R.A_all = sparse([1 0 0; 0 1 0; 0 0 1]);
R.gval  = [0; 0; 0];                  % all active (slack 0)
R.lam_g = [5; 3; 1e-12];              % eq; strong ineq; weak ineq
R.creg = struct('label',{'eqA','strongB','weakC'},'kind',{'eq','ineqHi','ineqHi'}, ...
    'rows',{1,2,3},'bound',{0,0,0},'node',{[],137,204});
K.sign = 1;
AS = sosc_active_set(R, K, sosc_defaults());
assert(AS.nEq==1 && AS.nStrong==1 && AS.nWeak==1, 'counts');
assert(AS.m_active==2, 'A = eq + strong only');
assert(any(contains(AS.weakLabels,'204')), 'weak label names the node');
assert(AS.licq, 'independent rows -> LICQ ok');
fprintf('test_sosc_active_set PASSED\n');
```

- [ ] **Step 2: Run — FAIL.**

- [ ] **Step 3: Implement `verify/sosc/sosc_active_set.m`**

```matlab
function AS = sosc_active_set(R, K, tol)
% SOSC_ACTIVE_SET  Classify inequality rows active/strong/weak, assemble the
% active Jacobian (equalities + strongly-active inequalities), flag LICQ.
%
% INPUTS: R - recover struct; K - kkt_residual struct (for sign); tol - defaults
% OUTPUTS: AS - struct .A[sparse m_active x n] .m_active .nEq .nStrong .nWeak
%               .weakLabels{cell} .licq
% REFERENCES: process/DESIGN_sosc.md sec 4.4.
muThresh = tol.mu * max(1, max(abs(R.lam_g)));
eqRows = []; strongRows = []; nWeak = 0; weakLabels = {};
for i = 1:numel(R.creg)
    c = R.creg(i);
    if strcmp(c.kind,'eq'), eqRows = [eqRows, c.rows]; continue; end %#ok<AGROW>
    for j = 1:numel(c.rows)
        r = c.rows(j);
        gv = R.gval(r);
        if strcmp(c.kind,'ineqHi'), slack = c.bound - gv; else, slack = gv - c.bound; end
        if slack < tol.active                       % active
            if abs(R.lam_g(r)) > muThresh
                strongRows = [strongRows, r]; %#ok<AGROW>
            else
                nWeak = nWeak + 1;
                nd = ''; if ~isempty(c.node), nd = sprintf(', node %d', c.node(min(j,numel(c.node)))); end
                weakLabels{end+1} = sprintf('%s%s (slack %.1e, lam %.1e)', c.label, nd, slack, R.lam_g(r)); %#ok<AGROW>
            end
        end
    end
end
actRows = sort([eqRows, strongRows]);
AS.A = R.A_all(actRows, :);
AS.m_active = numel(actRows);
AS.nEq = numel(eqRows); AS.nStrong = numel(strongRows); AS.nWeak = nWeak;
AS.weakLabels = weakLabels;
AS.licq = (sprank(AS.A) == AS.m_active);
end
```

- [ ] **Step 4: Run — PASS.**
- [ ] **Step 5: Commit** `git add verify/sosc/sosc_active_set.m tests/test_sosc_active_set.m && git commit -m "feat(sosc): active-set classification + active Jacobian (Task 6)"`

---

## Task 7: `sosc_inertia` — KKT inertia + subspace-SOSC bool

**Files:**
- Create: `verify/sosc/sosc_inertia.m`
- Test: `tests/test_sosc_inertia_qp.m`

**Interfaces:**
- Consumes: `H[nxn sparse]`, `A[m_a x n sparse]`, `tol`.
- Produces: `IN = sosc_inertia(H, A, tol)` → struct `.npos .nneg .nzero .expected[1x3] .subspaceOK[bool] .redMinEig`.

- [ ] **Step 1: Write the failing test** `tests/test_sosc_inertia_qp.m`

```matlab
% Hand-built KKT inertia cases (no NLP):
%   PD reduced Hessian -> inertia (n, m_a, 0), subspaceOK=true
%   Indefinite reduced Hessian -> wrong inertia, subspaceOK=false
%   Rank-deficient A -> nzero>0
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
tol = sosc_defaults();
% Case PD: H=I2, A=[1 1] -> reduced Hessian 1x1 = 1 > 0
IN = sosc_inertia(sparse(eye(2)), sparse([1 1]), tol);
assert(isequal([IN.npos IN.nneg IN.nzero],[2 1 0]) && IN.subspaceOK, 'PD case');
% Case indefinite: H=diag(1,-3), A=[1 1] -> reduced Hessian = (1-3)/2 = -1 < 0
IN2 = sosc_inertia(sparse(diag([1 -3])), sparse([1 1]), tol);
assert(~IN2.subspaceOK && IN2.nneg==2, 'indefinite case -> FAIL signature');
% Case rank-deficient A: two identical rows
IN3 = sosc_inertia(sparse(eye(3)), sparse([1 0 0; 1 0 0]), tol);
assert(IN3.nzero > 0, 'rank-deficient A -> nonzero nullity');
fprintf('test_sosc_inertia_qp PASSED\n');
```

- [ ] **Step 2: Run — FAIL.**

- [ ] **Step 3: Implement `verify/sosc/sosc_inertia.m`**

```matlab
function IN = sosc_inertia(H, A, tol)
% SOSC_INERTIA  Inertia of the KKT matrix [H A'; A 0] via sparse LDL^T; decide
% the subspace second-order condition (reduced Hessian PD on null(A)).
%
% INPUTS: H[nxn sparse] Lagrangian Hessian; A[m_a x n sparse] active Jacobian;
%         tol - sosc_defaults
% OUTPUTS: IN - struct .npos .nneg .nzero .expected[n m_a 0] .subspaceOK
%               .redMinEig (NaN placeholder; non-gating curvature margin)
% REFERENCES: process/DESIGN_sosc.md sec 4.5; Nocedal & Wright thm 16.3
%   (inertia(KKT)=(n,m_a,0) <=> reduced Hessian PD when A full row rank).
n = size(H,1);  ma = size(A,1);
K = [H, A.'; A, sparse(ma,ma)];
K = (K+K.')/2;                                  % symmetrize numerically
[~, D, ~] = ldl(K, 'vector');                   % D block-diagonal (1x1 & 2x2)
scale = max(1, normest(K));
zt = tol.inertiaZero * scale;
[npos, nneg, nzero] = deal(0,0,0);
i = 1; nD = size(D,1);
while i <= nD
    if i < nD && D(i+1,i) ~= 0                    % 2x2 block
        b = full(D(i:i+1,i:i+1)); ev = eig((b+b.')/2);
        for e = ev.'
            if e >  zt, npos=npos+1; elseif e < -zt, nneg=nneg+1; else, nzero=nzero+1; end
        end
        i = i + 2;
    else                                          % 1x1 block
        e = D(i,i);
        if e >  zt, npos=npos+1; elseif e < -zt, nneg=nneg+1; else, nzero=nzero+1; end
        i = i + 1;
    end
end
IN.npos=npos; IN.nneg=nneg; IN.nzero=nzero;
IN.expected = [n, ma, 0];
IN.subspaceOK = isequal([npos nneg nzero], [n ma 0]);
IN.redMinEig = NaN;   % optional non-gating margin; wired in a later enhancement
end
```

- [ ] **Step 4: Run — PASS** (all three cases).
- [ ] **Step 5: Commit** `git add verify/sosc/sosc_inertia.m tests/test_sosc_inertia_qp.m && git commit -m "feat(sosc): KKT-matrix inertia + subspace SOSC test (Task 7)"`

---

## Task 8: `verify_sosc_mee` orchestrator + verdict + tiered gate

**Files:**
- Create: `verify/sosc/verify_sosc_mee.m`
- Test: `tests/test_sosc_10N.m` (integration, solves) + `tests/test_sosc_verdict_logic.m` (synthetic, no solve)

**Interfaces:**
- Consumes: all of Tasks 1,3–7.
- Produces: `sosc = verify_sosc_mee(saved_or_path, opts)` → the verdict struct of DESIGN §5 (`.verdict .reason .status .drift .sign .kkt .active .inertia .redMinEig .thresholds .meta`). `opts` optional: `.tol` override.

- [ ] **Step 1: Write the failing tests.**

`tests/test_sosc_verdict_logic.m` (pure verdict mapping, no solve — call an exposed sub-helper `sosc_decide(K,AS,IN)`):
```matlab
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
mk = @(stat,pass) struct('signOK',stat,'pass',pass,'stat',0,'primalEq',0, ...
    'primalIneq',0,'dualFeas',0,'comp',0,'sign',1);
okK = mk(true,true);
% PASS: kkt ok, licq ok, no weak, inertia subspaceOK
AS1 = struct('licq',true,'nWeak',0,'m_active',1,'weakLabels',{{}},'nEq',1,'nStrong',0);
IN1 = struct('subspaceOK',true,'nzero',0,'npos',2,'nneg',1,'expected',[2 1 0]);
assert(strcmp(sosc_decide(okK,AS1,IN1).verdict,'PASS'));
% FAIL: inertia wrong, no weak, licq ok
IN2 = IN1; IN2.subspaceOK=false; IN2.nzero=0;
assert(strcmp(sosc_decide(okK,AS1,IN2).verdict,'FAIL'));
% INCONCLUSIVE: weak present
AS3 = AS1; AS3.nWeak=1; AS3.weakLabels={'thrHi, node 204'};
assert(strcmp(sosc_decide(okK,AS3,IN1).verdict,'INCONCLUSIVE'));
% ERROR: kkt not pass
badK = mk(false,false);
assert(strcmp(sosc_decide(badK,AS1,IN1).verdict,'ERROR'));
fprintf('test_sosc_verdict_logic PASSED\n');
```

`tests/test_sosc_10N.m` (integration):
```matlab
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));
sosc = verify_sosc_mee(fullfile(module_root(),'results','MEE_M2_10N.mat'));
fprintf('10 N verdict=%s stat=%.2e drift=%.2e inertia=[%d %d %d] expected=[%d %d %d]\n', ...
    sosc.verdict, sosc.kkt.stat, sosc.drift, sosc.inertia.npos, sosc.inertia.nneg, ...
    sosc.inertia.nzero, sosc.inertia.expected);
assert(strcmp(sosc.verdict,'PASS'), '10 N expected PASS, got %s (%s)', sosc.verdict, sosc.reason);
assert(strcmp(sosc.status,'certified-sosc'), 'status');
fprintf('test_sosc_10N PASSED\n');
```

- [ ] **Step 2: Run both — FAIL.**

- [ ] **Step 3: Implement `verify/sosc/verify_sosc_mee.m`** (with a nested/local `sosc_decide` exposed for the unit test — put `sosc_decide` in its own file `verify/sosc/sosc_decide.m` so the test can call it directly):

`verify/sosc/sosc_decide.m`:
```matlab
function v = sosc_decide(K, AS, IN)
% SOSC_DECIDE  Map KKT/active-set/inertia results to a verdict (DESIGN sec 5).
% OUTPUTS: v - struct .verdict .reason .status
if ~K.pass || ~K.signOK
    v.verdict='ERROR'; v.reason='KKT residual/sign check failed';
elseif ~AS.licq || AS.nWeak>0 || IN.nzero>0
    parts={}; if ~AS.licq, parts{end+1}='LICQ fails'; end
    if AS.nWeak>0, parts{end+1}=sprintf('%d weakly-active: %s',AS.nWeak,strjoin(AS.weakLabels,'; ')); end
    if IN.nzero>0, parts{end+1}=sprintf('%d zero KKT eigenvalue(s)',IN.nzero); end
    v.verdict='INCONCLUSIVE'; v.reason=strjoin(parts,'; ');
elseif IN.subspaceOK
    v.verdict='PASS'; v.reason='reduced Hessian PD on the critical cone (strict local min)';
else
    v.verdict='FAIL'; v.reason=sprintf('indefinite reduced Hessian: inertia [%d %d %d] != expected [%d %d %d]', ...
        IN.npos,IN.nneg,IN.nzero,IN.expected(1),IN.expected(2),IN.expected(3));
end
switch v.verdict
    case 'PASS',         v.status='certified-sosc';
    case 'FAIL',         v.status='feasible-only';
    otherwise,           v.status='certified-feasibility+sosc-inconclusive';
end
end
```

`verify/sosc/verify_sosc_mee.m`:
```matlab
function sosc = verify_sosc_mee(saved_or_path, opts)
% VERIFY_SOSC_MEE  NLP-level SOSC local-minimum certificate for a saved
% certified MEE min-fuel row. Orchestrates recover -> KKT re-check -> active
% set -> inertia -> verdict + tiered-gate status.
%
% INPUTS:
%   saved_or_path - a sosc_load_row struct OR a .mat path [char]
%   opts - optional struct: .tol (override sosc_defaults)
% OUTPUTS:
%   sosc - struct per process/DESIGN_sosc.md sec 5 (.verdict .reason .status
%          .drift .sign .kkt .active .inertia .redMinEig .thresholds .meta)
% REFERENCES: process/DESIGN_sosc.md secs 4-5.
if nargin<2, opts=struct(); end
tol = optdef(opts,'tol',sosc_defaults());
if ischar(saved_or_path)||isstring(saved_or_path), saved = sosc_load_row(char(saved_or_path));
else, saved = saved_or_path; end

R = sosc_recover_kkt(saved, tol);
sosc = struct('thresholds',tol,'drift',NaN,'sign',NaN, ...
    'kkt',[],'active',[],'inertia',[],'redMinEig',NaN, ...
    'meta',struct('thrustN',saved.thrustN,'tag',saved.tag,'when',datestr(now)));
if ~R.recoverOK
    sosc.verdict='ERROR'; sosc.reason=sprintf('warm re-solve failed: %s',R.ipoptStatus);
    sosc.status='certified-feasibility+sosc-inconclusive'; return;
end
sosc.drift = R.drift;
K  = sosc_kkt_residual(R, tol);   sosc.sign=K.sign;  sosc.kkt=K;
AS = sosc_active_set(R, K, tol);  sosc.active=AS;
IN = sosc_inertia(R.H, AS.A, tol);sosc.inertia=IN;
v  = sosc_decide(K, AS, IN);
sosc.verdict=v.verdict; sosc.reason=v.reason; sosc.status=v.status;
sosc.meta.n=R.n; sosc.meta.m=R.m; sosc.meta.m_active=AS.m_active;
if R.drift >= tol.drift
    warning('verify_sosc_mee:drift','warm re-solve drift %.2e >= %.1e (certifying the re-converged point)', R.drift, tol.drift);
end
end
```

- [ ] **Step 4: Run both tests.** Verdict-logic — PASS. 10 N integration — expect `verdict=PASS, status=certified-sosc`, small `stat`/`drift`. **This is the §6 tolerance-calibration point:** if `stat` or `drift` sits just above threshold on a clearly-good row, tighten/loosen the specific `tol` field in `sosc_defaults.m` and note why in a comment.

- [ ] **Step 5: Commit** `git add verify/sosc/sosc_decide.m verify/sosc/verify_sosc_mee.m tests/test_sosc_verdict_logic.m tests/test_sosc_10N.m && git commit -m "feat(sosc): orchestrator + verdict logic + tiered gate (Task 8)"`

---

## Task 9: Wire into the driver + reproducer (go-forward gate)

**Files:**
- Modify: `drivers/run_transfer_mee.m` (after `report`/`res` assembly, ~250–258)
- Modify: `reproduce/reproduce_row.m` (keep-best adoption)
- Test: `tests/test_sosc_gate_wiring.m`

**Interfaces:**
- Consumes: `verify_sosc_mee`.
- Produces: `res.sosc` on driver output; reproducer adopts a candidate only if `verdict ∈ {PASS, INCONCLUSIVE}` (never FAIL/ERROR-as-nonmin — ERROR is non-adopting here since we cannot certify).

**Context:** `run_transfer_mee` saves `res` only if `best.certified` (line 257). We attach `res.sosc` before that save and apply the tiered gate: a FAIL flips `report.certified=false` (demote) with a loud warning. Build a transient `saved`-shaped struct in-memory (no round-trip through disk) via `sosc_load_row`-equivalent fields.

- [ ] **Step 1: Write the failing test** `tests/test_sosc_gate_wiring.m`

```matlab
% The driver, run on the certified 10 N config (auto-reusing its cache path
% but forcing the SOSC attach), returns res.sosc with a verdict, and a FAIL
% verdict demotes report.certified. We test the gate mapping via a stub:
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
% Unit-level: apply_sosc_gate(report, sosc) demotes on FAIL, keeps otherwise.
rep = struct('certified',true);
gFail = apply_sosc_gate(rep, struct('verdict','FAIL','status','feasible-only','reason','x'));
assert(~gFail.certified, 'FAIL must demote');
gPass = apply_sosc_gate(rep, struct('verdict','PASS','status','certified-sosc','reason','y'));
assert(gPass.certified, 'PASS keeps certified');
gInc  = apply_sosc_gate(rep, struct('verdict','INCONCLUSIVE','status','certified-feasibility+sosc-inconclusive','reason','z'));
assert(gInc.certified, 'INCONCLUSIVE keeps certified (annotated)');
fprintf('test_sosc_gate_wiring PASSED\n');
```

- [ ] **Step 2: Run — FAIL** (`apply_sosc_gate` undefined).

- [ ] **Step 3: Implement.** Create `verify/sosc/apply_sosc_gate.m`:
```matlab
function report = apply_sosc_gate(report, sosc)
% APPLY_SOSC_GATE  Tiered gate: only a proven saddle (verdict FAIL) demotes a
% feasibility-certified row. PASS/INCONCLUSIVE/ERROR keep certified (annotated).
% INPUTS: report - run_transfer_mee report struct; sosc - verify_sosc_mee struct
% OUTPUTS: report - with .sosc attached and .certified possibly demoted.
report.sosc = sosc;
if strcmp(sosc.verdict,'FAIL')
    report.certified = false;
    warning('run_transfer_mee:soscFail', ...
        'SOSC FAIL (proven saddle) -> demoted to feasible-only: %s', sosc.reason);
end
end
```
Then in `run_transfer_mee.m`, after line 250 (`report` built) and before the `if best.certified` save at 257, when `best.certified` is true build the transient `saved` and certify:
```matlab
if best.certified
    savedT = struct('sigma',sigma,'X',best.X,'U',best.U,'dL',best.dL, ...
        'tfTarget',tf,'xf',xf,'thrustN',thrustN,'m0kg',m0kg,'ispS',ispS, ...
        'maxIter',maxIter,'tag',tag,'kind','MEE_M2');
    report = apply_sosc_gate(report, verify_sosc_mee(savedT));
    best.certified = report.certified;   % keep res.fuel + save gate consistent
end
```
In `reproduce/reproduce_row.m`, at the keep-best adoption of a fuel candidate, gate on the verdict — adopt the candidate into the best-mass pool only if `ismember(sosc.verdict,{'PASS','INCONCLUSIVE'})`; store `cand.sosc`. (Exact insertion at the `fuel_multistart` best-selection; consult the file — one added guard around the "is this candidate adoptable" test.)

- [ ] **Step 4: Run the unit test — PASS.** Then a driver smoke on 10 N to confirm `res.sosc` is attached and the row still certifies:
Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('$(pwd)'); setup_paths; addpath(fullfile(getenv('HOME'),'casadi-3.7.0')); r=run_transfer_mee(struct('thrustN',10,'ctf',1.5,'tfMinAnchor',22.2206,'seedThr',0.4,'betaMode','tangential','nodesPerRev',25)); disp(r.report.sosc.verdict)"`
Expected: prints `PASS` (or `INCONCLUSIVE`), row still certified.

- [ ] **Step 5: Commit** `git add verify/sosc/apply_sosc_gate.m drivers/run_transfer_mee.m reproduce/reproduce_row.m tests/test_sosc_gate_wiring.m && git commit -m "feat(sosc): wire SOSC gate into driver + reproducer (Task 9)"`

---

## Task 10: `recertify_table3` batch + sidecar report

**Files:**
- Create: `verify/sosc/recertify_table3.m`
- Test: `tests/test_sosc_recertify.m`

**Interfaces:**
- Produces: `T = recertify_table3(thrustList)` → struct array of `{thrustN, tag, verdict, drift, stat, inertia}`; writes `results/sosc/sosc_<tag>.mat` (`sosc` struct) per row + prints a summary table. **Never writes into `results/*.mat` campaign files.**

- [ ] **Step 1: Write the failing test** `tests/test_sosc_recertify.m`

```matlab
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));
campaign = fullfile(module_root(),'results','MEE_M2_10N.mat');
info0 = dir(campaign); before = info0.bytes;
T = recertify_table3(10);
assert(isscalar(T) && T(1).thrustN==10, 'one row for 10 N');
assert(ismember(T(1).verdict,{'PASS','INCONCLUSIVE','FAIL','ERROR'}), 'verdict set');
assert(isfile(fullfile(module_root(),'results','sosc','sosc_MEE_M2_10N.mat')), 'sidecar written');
info1 = dir(campaign);
assert(info1.bytes==before, 'campaign .mat must be untouched');
fprintf('test_sosc_recertify PASSED (10 N verdict=%s)\n', T(1).verdict);
```

- [ ] **Step 2: Run — FAIL.**

- [ ] **Step 3: Implement `verify/sosc/recertify_table3.m`**

```matlab
function T = recertify_table3(thrustList)
% RECERTIFY_TABLE3  Re-certify existing certified campaign rows with the SOSC
% certificate; write sidecar verdicts (results/sosc/) leaving campaign caches
% untouched, and print a summary table.
%
% INPUTS: thrustList - vector of thrust levels [N] (e.g. [10 5 2.5 1 0.5])
% OUTPUTS: T - struct array .thrustN .tag .verdict .drift .stat .inertia
% REFERENCES: process/DESIGN_sosc.md sec 8.
resDir  = fullfile(module_root(),'results');
sideDir = fullfile(resDir,'sosc'); if ~isfolder(sideDir), mkdir(sideDir); end
% tag map: the CERTIFIED HEADLINE row per rung. 10/5/2.5 N are the MEE_M2 fuel
% rows; 1 N and 0.5 N headline numbers (1371.44 kg, 1375.28 kg) are the
% PSR-refined solutions, so certify those PSR-final rows (sosc_load_row handles
% both the res and out shapes).
tagOf = containers.Map({10,5,2.5,1,0.5}, ...
    {'MEE_M2_10N','MEE_M2_5N','MEE_M2_2p5N', ...
     'MEE_M2_1N_PSR_psr_final','MEE_M2_0p5N_PSR_psr_final'});
T = struct('thrustN',{},'tag',{},'verdict',{},'drift',{},'stat',{},'inertia',{});
fprintf('\n  T[N]   tag                 verdict        drift      stat     inertia\n');
for Tn = thrustList(:).'
    tag = tagOf(Tn);
    mp  = fullfile(resDir,[tag '.mat']);
    if ~isfile(mp), fprintf('  %-5g  %-18s  MISSING (no cache)\n',Tn,tag); continue; end
    sosc = verify_sosc_mee(mp);
    save(fullfile(sideDir,['sosc_' tag '.mat']),'sosc');
    inb = sprintf('[%d %d %d]',sosc.inertia.npos,sosc.inertia.nneg,sosc.inertia.nzero);
    fprintf('  %-5g  %-18s  %-13s  %.2e  %.2e  %s\n', Tn, tag, sosc.verdict, sosc.drift, sosc.kkt.stat, inb);
    T(end+1) = struct('thrustN',Tn,'tag',tag,'verdict',sosc.verdict, ...
        'drift',sosc.drift,'stat',sosc.kkt.stat,'inertia',sosc.inertia); %#ok<AGROW>
end
end
```

- [ ] **Step 4: Run — PASS.** Then run the full ladder for the record:
`recertify_table3([10 5 2.5 1 0.5])` — capture the verdict table (a FAIL is a real finding; surface it).

- [ ] **Step 5: Commit** `git add verify/sosc/recertify_table3.m tests/test_sosc_recertify.m && git commit -m "feat(sosc): batch re-cert of campaign rows to sidecar (Task 10)"`

---

## Self-Review notes (author)

- **Spec coverage:** §3 files → Tasks 1–10 (all eight `verify/sosc/` files + the hook + two wiring points). §4 interfaces → Task steps carry the exact signatures. §5 verdict/gate → Tasks 8–9. §6 tolerances → Task 3 (+ calibration in Task 8.4). §7 error handling → ERROR path in Tasks 5/8, LICQ in 6/7, recon gate in 1. §8 plug-in → Tasks 9–10 (sidecar in 10). §9 testing → the QP unit (Task 7), active-set unit (6), recon (1), integration (8). §10 risks → Task 1 gate, CasADi-API note (Task 4.4), 2×2 inertia block handling (Task 7), MUMPS inherited (recover ERROR path).
- **Type consistency:** `saved`, `R`, `K`, `AS`, `IN`, `sosc` field names are used identically across tasks; `sosc_decide` is a separate file so its unit test can call it.
- **Deferred (YAGNI, per spec §1):** `redMinEig` is a NaN placeholder (non-gating); full critical-cone copositivity and OCP-level second-order are out of scope.
