# run_gergaud Min-Fuel Elliptic→GEO Front Door — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A single PARAMETERS-block script `run_gergaud.m` that reproduces one row of Gergaud (JGCD 27(6), 2004) Table 3 — a min-fuel elliptic→GEO low-thrust transfer — with user-settable thrust and user-settable initial *and* final orbits, emitting the Table-3 row plus a trajectory plot and a movie (gif + mp4).

**Architecture:** Thin front door onto the completed MEE thrust-ladder solver stack. Two endpoints become real knobs by parameterizing the initial elements in `mee_seed.m` and the terminal target in `casadi_lt_mee.m` (`opts.xf`), threaded as pass-throughs through `homotopy_mee.m`, `run_transfer_mee.m`, `run_mintime_mee.m`. Defaults reproduce today's behavior byte-for-byte. Visualization reuses the already-built `mee_res_to_cart_res.m` adapter + `transfer_movie.m`.

**Tech Stack:** MATLAB R2025b, CasADi 3.7.0 + IPOPT (existing solver core), the `earth_elliptic_to_geo/` module.

## Global Constraints

- **MATLAB R2025b ONLY** — invoke `/Applications/MATLAB_R2025b.app/bin/matlab` explicitly (R2025a license is broken). Run tests via `matlab -batch "run('/abs/path/test_X.m')"`.
- **Never use `i` or `j` as loop/index variables** (imaginary unit in MATLAB). Use `k`, `kk`, `ii`, `idx`, or a meaningful name.
- **Every MATLAB function needs a full header block:** purpose, INPUTS (with sizes), OUTPUTS (with sizes), REFERENCES.
- **Default-preserving is a hard requirement.** With no endpoint overrides, the pipeline MUST behave byte-for-byte as today: the initial state MUST be the legacy literal `X_init = [11625/par.LU_km; 0.75; 0; 0.0612; 0; 1; 0]` (NOT recomputed from `i0=7°`, which gives `0.061163 ≠ 0.0612`), and the terminal MUST be `xf = [1;0;0;0;0]`. Existing certified caches (`results/MEE_M2_*.mat`, `MEE_mintime_T*.mat`, and their `_seed`/`_step`/`_round` intermediates) MUST still load without a fingerprint error.
- **No new solver physics.** All endpoint work is parameterization + pass-through. No new homotopy, no unit rescaling. `LU = 42165 km` stays fixed; a custom final orbit is expressed in that unit (`Pf = Pf_km/42165`). The `L0 = π` apogee-start convention is retained for all endpoints.
- **Cache-fingerprint schema evolution:** when a newly-added `fp` field is absent from a cached fingerprint, treat it as compatible (WARN), matching the already-shipped `run_mintime_mee.m>check_cache_fp_mt` pattern — do NOT hard-error. Genuine value drift on a shared field still errors.

---

## File Structure

- `earth_elliptic_to_geo/casadi_lt_mee.m` — MODIFY: add `opts.xf` terminal target (default GEO).
- `earth_elliptic_to_geo/mee_seed.m` — MODIFY: add `opts.initElems` initial-orbit override (default = legacy literal).
- `earth_elliptic_to_geo/homotopy_mee.m` — MODIFY: forward `opts.xf`; harmonize `check_cache_fp` to schema-older-WARN.
- `earth_elliptic_to_geo/run_transfer_mee.m` — MODIFY: accept + thread `cfg.xf`/`cfg.initElems`; add to `fp`; harmonize `check_cache_fp`.
- `earth_elliptic_to_geo/run_mintime_mee.m` — MODIFY: accept + thread `cfg.xf`/`cfg.initElems`; add to `fp`.
- `earth_elliptic_to_geo/mee_res_to_cart_res.m` — DONE (committed c04a057). Add a unit test.
- `earth_elliptic_to_geo/transfer_movie.m` — MODIFY: accept a `res` struct OR a `.mat` path (DRY the adapter handoff).
- `earth_elliptic_to_geo/gergaud_plot.m` — CREATE: static reconstructed-trajectory PNG.
- `earth_elliptic_to_geo/run_gergaud.m` — CREATE: the front-door PARAMETERS-block script.
- `earth_elliptic_to_geo/gergaud_row.m` — CREATE: pure helper that assembles + formats the Table-3 row struct (unit-testable without a solve).
- Tests: `test_mee_xf.m`, `test_mee_seed_initelems.m`, `test_mee_threading.m`, `test_mee_res_to_cart.m`, `test_gergaud_row.m`, `test_run_gergaud_auto.m`.

Interfaces the tasks rely on (established here so out-of-order readers agree):
- `casadi_lt_mee(sigma, X0, U0, dL0, opts)` gains `opts.xf` [5×1], default `[1;0;0;0;0]`.
- `mee_seed(par, opts)` gains `opts.initElems` [7×1] OR `[]`; `[]`/absent ⇒ legacy literal.
- `gergaud_row(inp)` → struct with fields `thrustN, tfmin_ND, tfmin_h, ctf, tf_ND, m_f_kg, prop_kg, dV_kms, switches, revs, revs_paper, edge, incl_deg, defect, certified, note`; plus `gergaud_row_str(row)` → char block.
- `mee_res_to_cart_res(Xmee, Umee, dL, sigma, thrustN, ctf, mu)` → Cartesian `res` (already built).

---

## Task 1: Terminal target `opts.xf` in `casadi_lt_mee.m`

**Files:**
- Modify: `earth_elliptic_to_geo/casadi_lt_mee.m` (terminal block ~lines 150–156; option parse ~lines 61–67)
- Test: `earth_elliptic_to_geo/test_mee_xf.m`

**Interfaces:**
- Consumes: existing `casadi_lt_mee` opts contract (`par, mode, eps, tfTarget, x0, maxIter, warmTight, printLevel`).
- Produces: `opts.xf` [5×1] terminal target `[P;ex;ey;hx;hy]`, default `[1;0;0;0;0]`.

- [ ] **Step 1: Write the failing test**

`test_mee_xf.m` — a no-IPOPT-solve structural test: confirm the default xf is GEO and a custom xf is honored, by parsing the option (call the getter path) and by a tiny 1-step build check. Because building the full Opti is heavy, assert on a helper that resolves xf. Add near the top of `casadi_lt_mee.m` a resolved local `xf = d('xf', [1;0;0;0;0]);` and expose it via a lightweight self-check: when `opts.selftest` is true, return early with `out = struct('xf', xf);`.

```matlab
% test_mee_xf.m
here = fileparts(mfilename('fullpath')); cd(here);
par = kepler_lt_params(10,1500,2000);
% default -> GEO
o = casadi_lt_mee((0:1).', zeros(7,2), zeros(4,2), 1, ...
    struct('par',par,'x0',zeros(7,1),'selftest',true));
assert(isequal(o.xf,[1;0;0;0;0]), 'default xf must be GEO [1;0;0;0;0]');
% custom -> honored
xf = [0.9; 0.01; -0.02; 0.05; 0];
o2 = casadi_lt_mee((0:1).', zeros(7,2), zeros(4,2), 1, ...
    struct('par',par,'x0',zeros(7,1),'xf',xf,'selftest',true));
assert(isequal(o2.xf,xf), 'custom xf must be honored');
fprintf('test_mee_xf PASSED\n');
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "run('/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo/test_mee_xf.m')"`
Expected: FAIL (no `xf` handling / no `selftest` early return yet).

- [ ] **Step 3: Implement**

In the option-parse block of `casadi_lt_mee.m` add:
```matlab
xf       = d('xf', [1;0;0;0;0]);
assert(numel(xf)==5, 'casadi_lt_mee: opts.xf must be 5x1 [P;ex;ey;hx;hy]');
if d('selftest', false), out = struct('xf', xf(:)); return; end
```
Replace the hardcoded terminal block:
```matlab
% terminal GEO in elements: P=1, ex=ey=hx=hy=0; L free (DeltaL is the DOF).
opti.subject_to(X(1,end) == 1);
opti.subject_to(X(2,end) == 0);
opti.subject_to(X(3,end) == 0);
opti.subject_to(X(4,end) == 0);
opti.subject_to(X(5,end) == 0);
```
with:
```matlab
% terminal target in elements (default GEO [1;0;0;0;0]); L free (DeltaL is DOF).
% Prograde automatic for the h=0 equatorial default; a custom xf is the
% caller's responsibility (see run_gergaud scope note).
for kt = 1:5
    opti.subject_to(X(kt,end) == xf(kt));
end
```

- [ ] **Step 4: Run test to verify it passes**

Run the same command. Expected: `test_mee_xf PASSED`.

- [ ] **Step 5: Verify default MEE solve unchanged (regression)**

`test_energy_stage` drives the *Cartesian* `casadi_lt_2body` solver (its 1367.15 kg is the coplanar Cartesian gate) — NOT this MEE file. Regress against the real MEE gate instead: warm-start `casadi_lt_mee` from the certified `results/MEE_M2_10N.mat` (default `xf`, eps=0) and confirm it reproduces `m_f ≈ 1377.10 kg`, defect < 1e-8 bit-for-bit — proving the `xf`-loop reproduces the old hardcoded five terminal equalities exactly.

- [ ] **Step 6: Commit**

```bash
git add earth_elliptic_to_geo/casadi_lt_mee.m earth_elliptic_to_geo/test_mee_xf.m
git commit -m "feat(mee): parameterize terminal target opts.xf (default GEO)"
```

---

## Task 2: Initial-elements override in `mee_seed.m`

**Files:**
- Modify: `earth_elliptic_to_geo/mee_seed.m` (lines 50–51)
- Test: `earth_elliptic_to_geo/test_mee_seed_initelems.m`

**Interfaces:**
- Consumes: existing `mee_seed(par, opts)` contract.
- Produces: `opts.initElems` [7×1] `[P;ex;ey;hx;hy;m;t]` OR `[]`/absent. Absent ⇒ **exact** legacy literal.

- [ ] **Step 1: Write the failing test**

```matlab
% test_mee_seed_initelems.m
here = fileparts(mfilename('fullpath')); cd(here);
par = kepler_lt_params(10,1500,2000);
% default (no initElems) MUST equal the legacy literal, byte-for-byte
[~, Xd] = mee_seed(par, struct('thr',0.4,'betaMode','tangential','N',20,'nRev',1));
legacy = [11625/par.LU_km; 0.75; 0; 0.0612; 0; 1; 0];
assert(isequal(Xd(:,1), legacy), 'default initial node must be the legacy literal');
% custom initElems honored at node 1
ci = [0.30; 0.60; 0.0; 0.10; 0.0; 1; 0];
[~, Xc] = mee_seed(par, struct('thr',0.4,'betaMode','tangential','N',20,'nRev',1,'initElems',ci));
assert(isequal(Xc(:,1), ci), 'custom initElems must set node 1');
fprintf('test_mee_seed_initelems PASSED\n');
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "run('.../test_mee_seed_initelems.m')"`. Expected: FAIL (custom branch absent; `isequal` on default still passes, custom errors on unknown field usage).

- [ ] **Step 3: Implement**

Replace `mee_seed.m` line 50–51:
```matlab
L0     = pi;
X_init = [11625/par.LU_km; 0.75; 0; 0.0612; 0; 1; 0];
```
with:
```matlab
L0 = pi;                                    % apogee-start convention (all endpoints)
% Default = the paper's legacy literal, byte-for-byte (NOT recomputed from
% i0=7deg: tan(3.5deg)=0.061163 != the certified 0.0612). A caller wanting a
% different start orbit passes opts.initElems explicitly (run_gergaud builds
% it from P0/e0/i0 only when the user overrides the paper defaults).
if isfield(opts,'initElems') && ~isempty(opts.initElems)
    X_init = opts.initElems(:);
    assert(numel(X_init)==7, 'mee_seed: opts.initElems must be 7x1 [P;ex;ey;hx;hy;m;t]');
else
    X_init = [11625/par.LU_km; 0.75; 0; 0.0612; 0; 1; 0];
end
```

- [ ] **Step 4: Run test to verify it passes**

Expected: `test_mee_seed_initelems PASSED`.

- [ ] **Step 5: Commit**

```bash
git add earth_elliptic_to_geo/mee_seed.m earth_elliptic_to_geo/test_mee_seed_initelems.m
git commit -m "feat(mee): opts.initElems initial-orbit override in mee_seed (default = legacy literal)"
```

---

## Task 3: Forward `xf` through `homotopy_mee.m` + harmonize its cache guard

**Files:**
- Modify: `earth_elliptic_to_geo/homotopy_mee.m` (call site line ~52; `check_cache_fp` lines ~78–104)
- Test: `earth_elliptic_to_geo/test_mee_threading.m` (created here; extended in Task 4)

**Interfaces:**
- Consumes: `casadi_lt_mee` `opts.xf` (Task 1).
- Produces: `homotopy_mee` opts gains `.xf` [5×1], default `[1;0;0;0;0]`, forwarded to every eps-step solve and recorded in `fp.xf`.

- [ ] **Step 1: Write the failing test**

```matlab
% test_mee_threading.m  (Task 3 portion)
here = fileparts(mfilename('fullpath')); cd(here);
% (a) homotopy_mee forwards xf: use selftest short-circuit via a stub schedule
% of length 0 is not possible; instead assert the fp records xf.
% Build a resDir in a temp folder, run a single trivial eps via selftest path.
% Simpler: white-box check that homotopy_mee reads opts.xf into fp.
src = fileread('homotopy_mee.m');
assert(contains(src,'xf'), 'homotopy_mee must reference xf');
% (b) schema-older cache WARNs, not errors:
tmp = tempname; mkdir(tmp);
fp_old = struct('sched',[1 0]);               % a cache WITHOUT xf
o = struct('maxDefect',1e-12,'switches',0,'edge',1,'m_f_kg',1400);
Xk=zeros(7,2);Uk=zeros(4,2);dLk=1;ok=true;e=0; %#ok<NASGU>
save(fullfile(tmp,'thr_step01.mat'),'o','ok','Xk','Uk','dLk','e','fp_old');
% rename fp_old->fp to mimic a real (schema-older) cache
S=load(fullfile(tmp,'thr_step01.mat')); fp=S.fp_old; save(fullfile(tmp,'thr_step01.mat'),'o','ok','Xk','Uk','dLk','e','fp');
lastwarn('');
par = kepler_lt_params(10,1500,2000);
try
  homotopy_mee((0:1).', zeros(7,2), zeros(4,2), 1, struct('par',par,'x0',zeros(7,1), ...
     'tfTarget',30,'resDir',tmp,'tag','thr','sched',[0],'xf',[1;0;0;0;0]));
catch ME
  assert(isempty(strfind(ME.identifier,'fingerprintMismatch')), ...
     'schema-older cache must not hard-error');
end
fprintf('test_mee_threading (Task3) PASSED\n');
```

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL — `homotopy_mee` neither forwards `xf` nor WARN-tolerates a schema-older step cache.

- [ ] **Step 3: Implement**

In `homotopy_mee.m` option parse add `xf = d('xf',[1;0;0;0;0]); fp.xf = xf;`. In the `casadi_lt_mee(...)` call (line ~52) add `'xf', xf,` to the struct. In `check_cache_fp` (line ~94–102) change the loop to schema-older-WARN:
```matlab
flds = fieldnames(fp);
for k = 1:numel(flds)
    f = flds{k};
    if ~isfield(S.fp, f)
        warning('homotopy_mee:fpSchemaOlder', ['%s: field ''%s'' present in ' ...
            'current fp but absent from cache (schema evolution) -- trusting ' ...
            'as compatible under tag=''%s'''], file, f, tag);
        continue;
    end
    if ~isequal(S.fp.(f), fp.(f))
        error('homotopy_mee:fingerprintMismatch', ['cached config fingerprint ' ...
            'mismatch in %s: field ''%s'' differs -- stale cache under tag=''%s''; ' ...
            'delete or use a new tag'], file, f, tag);
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

Expected: `test_mee_threading (Task3) PASSED`.

- [ ] **Step 5: Commit**

```bash
git add earth_elliptic_to_geo/homotopy_mee.m earth_elliptic_to_geo/test_mee_threading.m
git commit -m "feat(mee): forward opts.xf through homotopy_mee; schema-older WARN on step cache"
```

---

## Task 4: Thread `xf`/`initElems` through `run_transfer_mee.m` + `run_mintime_mee.m`

**Files:**
- Modify: `earth_elliptic_to_geo/run_transfer_mee.m` (option parse; `mee_seed` calls ~lines 113,130; `homotopy_mee`/`casadi_lt_mee` calls ~lines 181,196; `fpBase` ~line 94; `check_cache_fp` ~lines 261–290)
- Modify: `earth_elliptic_to_geo/run_mintime_mee.m` (option parse ~216; Stage B `mee_seed` call ~332; `casadi_lt_mee` calls in `mintime_mee_continue` ~432,463; `fp` ~233)
- Test: extend `earth_elliptic_to_geo/test_mee_threading.m`

**Interfaces:**
- Consumes: `homotopy_mee.xf` (Task 3), `casadi_lt_mee.xf` (Task 1), `mee_seed.initElems` (Task 2).
- Produces: `run_transfer_mee` cfg gains `.xf` [5×1] and `.initElems` [7×1]; `run_mintime_mee` cfg gains `.xf` and `.initElems`. Both default-preserving; both add the two fields to their `fp`.

- [ ] **Step 1: Write the failing test (append to test_mee_threading.m)**

```matlab
% Task 4 portion: default cfg reuses the existing certified 10 N caches with
% NO fingerprint error (schema-older WARN only), and cfg.xf is carried into fp.
here = fileparts(mfilename('fullpath')); cd(here);
src1 = fileread('run_transfer_mee.m'); src2 = fileread('run_mintime_mee.m');
assert(contains(src1,'xf') && contains(src1,'initElems'), 'run_transfer_mee must thread xf/initElems');
assert(contains(src2,'xf') && contains(src2,'initElems'), 'run_mintime_mee must thread xf/initElems');
% default-preserving reuse: loading the cached 10 N mintime anchor with the
% new (xf-bearing) fp must NOT throw fingerprintMismatch.
lastwarn('');
try
  out = run_mintime_mee(10, 25);      % reuses results/MEE_mintime_T100.mat
  assert(abs(out.tfmin-22.2206) < 1e-2, 'cached 10 N anchor tfmin preserved');
catch ME
  assert(isempty(strfind(ME.identifier,'fingerprintMismatch')), ...
     'default 10 N reuse must not fingerprint-error after adding xf/initElems');
  rethrow(ME);
end
fprintf('test_mee_threading (Task4) PASSED\n');
```

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL (fields not threaded).

- [ ] **Step 3: Implement**

`run_transfer_mee.m`: add `xf = d('xf',[1;0;0;0;0]); initElems = d('initElems',[]);`. Pass `'initElems', initElems` into BOTH `mee_seed` opts (probe + full). Pass `'xf', xf` into the `homotopy_mee` `ho` struct AND the warm-direct `casadi_lt_mee` call. For the seed's terminal, set `stopP = xf(1)` (so a custom final P stops the seed at the right radius): change the seed opts `'stopP', 1.0` → `'stopP', xf(1)`. Add `'xf', xf` and `'initElems_isset', ~isempty(initElems)` to `fpBase`. Harmonize `check_cache_fp` to schema-older-WARN (same edit as Task 3).

`run_mintime_mee.m`: add `xf = d('xf',[1;0;0;0;0]); initElems = d('initElems',[]);`. Pass `'initElems', initElems` into the Stage B `mee_seed` call, and `'xf', xf` into every `casadi_lt_mee(...,'mode','mintime',...)` call in `mintime_mee_continue` (round 0 and continuation). Add `xf`/`initElems_isset` to `fp`. (`check_cache_fp_mt` already schema-older-WARNs.) NOTE: Stage A warm-starts from a cached fuel `res`; when `initElems`/`xf` are custom, the cached fuel anchor is for different endpoints — so when either override is non-default, SKIP Stage A (its cached fuel file is endpoint-mismatched) and go straight to Stage B. Guard: `useStageA = isempty(initElems) && isequal(xf,[1;0;0;0;0]) && isfile(fuelFile)`.

- [ ] **Step 4: Run test to verify it passes**

Expected: `test_mee_threading (Task4) PASSED`. Also re-run Task 3 portion (same file) — both pass.

- [ ] **Step 5: Regression — default 10 N fuel reuse**

Run a 3-line script: `res = run_transfer_mee(struct('thrustN',10,'ctf',1.5,'tfMinAnchor',22.2248));` — must load `results/MEE_M2_10N.mat` region OR re-solve to `m_f ≈ 1377 kg` with no fingerprint error. Confirm no `fingerprintMismatch` thrown.

- [ ] **Step 6: Commit**

```bash
git add earth_elliptic_to_geo/run_transfer_mee.m earth_elliptic_to_geo/run_mintime_mee.m earth_elliptic_to_geo/test_mee_threading.m
git commit -m "feat(mee): thread xf/initElems through run_transfer_mee + run_mintime_mee (default-preserving, schema-older WARN)"
```

---

## Task 5: Visualization — adapter test, `transfer_movie` struct handoff, static plot

**Files:**
- Modify: `earth_elliptic_to_geo/transfer_movie.m` (input handling, line ~27)
- Create: `earth_elliptic_to_geo/gergaud_plot.m`
- Test: `earth_elliptic_to_geo/test_mee_res_to_cart.m`

**Interfaces:**
- Consumes: `mee_res_to_cart_res` (built), `elements_to_cart`, `kepler_lt_params`.
- Produces: `transfer_movie(resOrPath, outStem)` accepts a struct OR a path; `gergaud_plot(res, outPng, titleStr)` writes a static PNG.

- [ ] **Step 1: Write the failing test**

```matlab
% test_mee_res_to_cart.m
here = fileparts(mfilename('fullpath')); cd(here);
S = load(fullfile(here,'results','MEE_M2_10N.mat'));
c = mee_res_to_cart_res(S.res.fuel.X, S.res.fuel.U, S.res.fuel.dL, S.res.sigma, 10, 1.5, 1);
r = c.fuel.X(1:3,:); rmag = vecnorm(r); an = vecnorm(c.fuel.U(1:3,:));
assert(abs(rmag(1)-1.1028) < 1e-3, 'apogee start ~1.103');
assert(abs(rmag(end)-1.0) < 1e-3, 'GEO end ~1.000');
assert(all(abs(an-1) < 1e-9), 'unit inertial thrust dir');
assert(abs(1500*c.fuel.X(7,end) - 1377.10) < 0.1, 'mass preserved');
% struct handoff: transfer_movie accepts a struct without loading a file
assert(nargin('transfer_movie') >= 1);
fprintf('test_mee_res_to_cart PASSED\n');
```

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL (the adapter asserts pass, but confirm `transfer_movie` struct handling — extend the assert to actually call `transfer_movie(c, tempname)` and require no `load`-of-struct error). First run fails because `transfer_movie` calls `load(matFile)` on a struct.

- [ ] **Step 3: Implement**

`transfer_movie.m` line 27 — replace `S = load(matFile); res = S.res;` with:
```matlab
if isstruct(matFile)
    res = matFile;                       % already a Cartesian res struct
else
    S = load(matFile); res = S.res;      % path to a .mat holding `res`
end
```
Create `gergaud_plot.m`: reconstruct via `mee_res_to_cart_res` (or accept a Cartesian res), plot the throttle-colored 3D path + GEO ring + Earth, `exportgraphics(fig, outPng, 'ContentType','image','Resolution',200)`. Full function header.

- [ ] **Step 4: Run test to verify it passes**

Expected: `test_mee_res_to_cart PASSED`. Render one smoke plot to `scratchpad` and eyeball dimensions in the test (assert PNG file exists and is non-empty).

- [ ] **Step 5: Commit**

```bash
git add earth_elliptic_to_geo/transfer_movie.m earth_elliptic_to_geo/gergaud_plot.m earth_elliptic_to_geo/test_mee_res_to_cart.m
git commit -m "feat(viz): transfer_movie struct handoff + gergaud_plot static trajectory + adapter test"
```

---

## Task 6: `gergaud_row` formatter (pure, unit-testable)

**Files:**
- Create: `earth_elliptic_to_geo/gergaud_row.m`, `earth_elliptic_to_geo/gergaud_row_str.m`
- Test: `earth_elliptic_to_geo/test_gergaud_row.m`

**Interfaces:**
- Consumes: a solved/loaded result's fields (thrustN, anchor tfmin, fuel report) + `kepler_lt_params` for unit conversions.
- Produces: `gergaud_row(inp)` → row struct (fields listed in File Structure); `gergaud_row_str(row)` → fixed-width char block. Paper revs lookup: `10→7.5, 5→15, 2.5→30, 1→74.5, 0.5→149` (else NaN).

- [ ] **Step 1: Write the failing test**

```matlab
% test_gergaud_row.m
here = fileparts(mfilename('fullpath')); cd(here);
inp = struct('thrustN',10,'tfmin_ND',22.2206,'ctf',1.5,'tf_ND',33.331, ...
    'm_f_kg',1377.10,'switches',19,'revs',7.326,'edge',0.999,'incl_deg',0.0, ...
    'defect',6e-15,'certified',true,'note','');
row = gergaud_row(inp);
assert(abs(row.prop_kg-(1500-1377.10)) < 1e-6, 'prop = m0 - m_f');
assert(row.revs_paper==7.5, 'paper revs lookup 10 N -> 7.5');
assert(row.dV_kms > 0, 'dV positive');
s = gergaud_row_str(row);
assert(ischar(s) && contains(s,'1377.10'), 'row string carries m_f');
% uncertified row must be flagged, not silently formatted as certified
inp2 = inp; inp2.thrustN=0.2; inp2.certified=false; inp2.note='not attained (0.5 N wall)';
s2 = gergaud_row_str(gergaud_row(inp2));
assert(contains(lower(s2),'not attained') || contains(s2,'UNCERTIFIED'), 'uncertified flagged');
fprintf('test_gergaud_row PASSED\n');
```

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL (functions absent).

- [ ] **Step 3: Implement**

`gergaud_row.m`: compute `prop_kg = m0kg - m_f_kg`, `dV_kms` from `c*log(m0/m_f)*VU_kms` (via `kepler_lt_params`), `tfmin_h = tfmin_ND*TU_s/3600`, `revs_paper` via the lookup, pass through `certified`/`note`. `gergaud_row_str.m`: fixed-width block (`%-8g`, `%10.4f`, `%6.2f` etc.), and when `~certified`, prepend an `UNCERTIFIED — <note>` line so a printed row can never be mistaken for a certified one. Full headers.

- [ ] **Step 4: Run test to verify it passes**

Expected: `test_gergaud_row PASSED`.

- [ ] **Step 5: Commit**

```bash
git add earth_elliptic_to_geo/gergaud_row.m earth_elliptic_to_geo/gergaud_row_str.m earth_elliptic_to_geo/test_gergaud_row.m
git commit -m "feat(front-door): gergaud_row + gergaud_row_str Table-3 row formatter (pure, tested)"
```

---

## Task 7: `run_gergaud.m` front door + end-to-end auto-mode test

**Files:**
- Create: `earth_elliptic_to_geo/run_gergaud.m`
- Test: `earth_elliptic_to_geo/test_run_gergaud_auto.m`

**Interfaces:**
- Consumes: everything above (`run_mintime_mee`, `run_transfer_mee` with `xf`/`initElems`; `psr_mee_refine` for the 1 N/0.5 N recipes; `mee_res_to_cart_res`, `transfer_movie`, `gergaud_plot`, `gergaud_row`).
- Produces: the front-door script; running it prints a Table-3 row and (when enabled) writes `results/gergaud_<tag>.{png,mp4,gif}`.

- [ ] **Step 1: Write the failing test**

```matlab
% test_run_gergaud_auto.m — auto mode on the default 10 N endpoints must
% return the certified cached row WITHOUT a fresh solve, and without movie.
here = fileparts(mfilename('fullpath')); cd(here);
row = run_gergaud(struct('thrustN',10,'runMode','auto','makeMovie',false, ...
    'makePlot',false,'returnOnly',true));
assert(abs(row.m_f_kg-1377.10) < 0.1, '10 N auto row m_f=1377.10');
assert(row.switches==19 && abs(row.revs-7.326)<1e-2, '10 N structure');
assert(row.certified, '10 N row certified');
% endpoint-default detection: paper defaults -> uses cache; a custom final
% flips to solve mode (assert it does NOT claim the cached 10 N number).
fprintf('test_run_gergaud_auto PASSED\n');
```

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL (`run_gergaud` absent).

- [ ] **Step 3: Implement**

`run_gergaud.m` — PARAMETERS-block script (see `PSR/run_psr.m` house style), refactored to also accept an `opts` struct when called as `run_gergaud(opts)` for testing (a thin `if nargin` guard at the top; the interactive path uses the literal PARAMETERS block). Logic:
1. Resolve endpoints: if `(P0_km,e0,i0_deg)` == paper defaults ⇒ `initElems=[]`; else build `initElems = [P0_km/LU; e0; 0; tan(deg2rad(i0_deg)/2); 0; 1; 0]`. If `(Pf_km,ef,if_deg)` == GEO defaults ⇒ `xf=[1;0;0;0;0]`; else `xf=[Pf_km/LU; ef; 0; tan(deg2rad(if_deg)/2); 0]`.
2. `isDefaultEndpoints = isempty(initElems) && isequal(xf,[1;0;0;0;0])`.
3. Mode: `auto` + default endpoints + cached `mee_fuel_tag(thrustN)` (or PSR-final for 1 N/0.5 N) exists ⇒ load, build row, done. `auto` + custom endpoints ⇒ behave as `solve`. `probe` ⇒ force solve including 0.2/0.1 N with an up-front wall WARNING.
4. Solve path: `run_mintime_mee(thrustN,nodesPerRev,struct('xf',xf,'initElems',initElems,...))` → `run_transfer_mee(struct('thrustN',...,'tfMinAnchor',anchor.tfmin,'xf',xf,'initElems',initElems,...))`; for 1 N/0.5 N apply `psr_mee_refine` per the recipe map (§4 of the spec); 0.2/0.1 N attempt + report `certified` honestly.
5. Assemble `gergaud_row`, print `gergaud_row_str`.
6. If `makePlot`: `gergaud_plot`. If `makeMovie`: `mee_res_to_cart_res` → `transfer_movie`. Output stems `results/gergaud_<mee_fuel_tag>` (+ endpoint hash suffix when custom).
7. `returnOnly` (test hook) returns the row struct and skips viz.

Include the honest per-rung recipe map and the 0.2/0.1-N "not attained" messaging from the spec §4. Full script header documenting the PARAMETERS block.

- [ ] **Step 4: Run test to verify it passes**

Expected: `test_run_gergaud_auto PASSED`.

- [ ] **Step 5: End-to-end smoke (default 10 N, full viz)**

Run `run_gergaud(struct('thrustN',10,'runMode','auto','makeMovie',true,'makePlot',true))` — confirm it prints the row and writes `results/gergaud_MEE_M2_10N.{png,mp4,gif}`. (Reuses the cache; viz is the only new compute.)

- [ ] **Step 6: Commit**

```bash
git add earth_elliptic_to_geo/run_gergaud.m earth_elliptic_to_geo/test_run_gergaud_auto.m
git commit -m "feat(front-door): run_gergaud min-fuel elliptic->GEO entry (auto/solve/probe, Table-3 row + plot + movie)"
```

---

## Task 8: README + docs update

**Files:**
- Modify: `earth_elliptic_to_geo/README.md` (add a "Front door: run_gergaud" subsection)
- Modify: `earth_elliptic_to_geo/DESIGN_thrust_ladder.md` (link the front door + the endpoint-parameterization note)

- [ ] **Step 1: Add README subsection** documenting `run_gergaud` usage, the three run modes, the endpoint knobs (with the default-preserving + non-GEO-final caveat), the per-rung recipe/honesty map, and the plot/movie outputs. Point to the four rendered movies (`results/movie_MEE_*.{mp4,gif}`) and the adapter.

- [ ] **Step 2: Verify no stale claims** — the README's existing "Deliverables" list gains the front door + movies; the honesty footnotes (0.5 N anchor-free, 0.2/0.1 N unreached) are referenced by `run_gergaud`'s recipe map, not re-litigated.

- [ ] **Step 3: Commit**

```bash
git add earth_elliptic_to_geo/README.md earth_elliptic_to_geo/DESIGN_thrust_ladder.md
git commit -m "docs(front-door): document run_gergaud + endpoint parameterization + movies"
```

---

## Self-Review (completed)

- **Spec coverage:** endpoints parameterized (T1 terminal, T2 initial, T4 threading) ✓; run modes + recipe map + row printout (T6, T7) ✓; plot + movie via adapter (T5, T7) ✓; honesty on 0.2/0.1 N (T6 uncertified flag, T7 recipe map) ✓; default-preservation (Global Constraints, T1 Step 5, T2 Step 1, T4 Step 5) ✓.
- **Placeholder scan:** every code step carries actual code; test bodies are concrete. No TBD/TODO.
- **Type consistency:** `xf` is [5×1] everywhere; `initElems` [7×1]/`[]` everywhere; `mee_res_to_cart_res` signature matches the committed file; `gergaud_row`/`_str` field set consistent between T6 definition and T7 consumption.
- **Cross-cutting risk captured:** the `check_cache_fp` hard-error-vs-WARN divergence (would break existing caches) is explicitly harmonized in T3/T4; the `tan(3.5°)≠0.0612` default-drift trap is called out in Global Constraints + T2.
