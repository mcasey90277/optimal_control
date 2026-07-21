# Thrust-Ladder Prep (P2) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make both CR3BP min-fuel pipelines (GTO_tulip sundman engine, GTO_ELFO freetf engine) thrust-ladder-ready: config fingerprints on every cached artifact, adaptive bounds + bound-saturation diagnostics, alias-free cross-rung chaining — validated by two live 20 mN pilot rungs.

**Architecture:** Three shared helpers land in `cr3bp_common/` (fingerprint build/check + thrust tag); each campaign then gets solver-level `boundSat` (output-only) + opt-in box overrides, fp wiring at its cache boundaries, and a chain helper (tulip: time-rescale + no-resample `sundman_seed_map` re-map; ELFO: pass-through — `cScale` already decouples the clock). Small review fixes (rF cascades, anchor gate, tf-sweep re-clean) fold in. Spec: `docs/superpowers/specs/2026-07-21-ladder-prep-design.md` (authoritative — read it).

**Tech Stack:** MATLAB R2025b (`/Applications/MATLAB_R2025b.app/bin/matlab -batch "<ONE line>"`), CasADi 3.7.0 (`~/casadi-3.7.0`), git.

## Global Constraints

- **Repo root** `/Users/msc/Desktop/optimal_control`; campaign dirs `orbit_transfer/GTO_tulip/direct/sundman_minfuel` (=`SMF` below) and `orbit_transfer/GTO_ELFO/direct/elfo` (=`ELFO` below); shared lib `orbit_transfer/cr3bp_common` (=`CC`).
- **Back-compat invariant (binding):** at nominal thrust (25 mN) with no new options, every path is byte-identical — same cache filenames, same NLP (bounds/values), same results. New solver behavior is opt-in via options; `boundSat` is output-only. Any test that shows a nominal-path change is a STOP.
- **MATLAB house style:** full comment headers (purpose/INPUTS with sizes/OUTPUTS with sizes/REFERENCES); never `i`/`j` as loop vars. One-line `-batch` commands only (zsh mangles multi-line). Filter license noise: `| grep -vE "License|academic|personal use|government"`.
- **Certified caches are read-only**: never overwrite `sundman_minfuel_certified.mat` or the ELFO seed bank; pilots write only new `_T20mN`-tagged artifacts. `.mat` files are gitignored — never `git add` one.
- **Staging discipline:** never `git add -A`; stage exactly the files each task names; `git status --short` before every commit (`papers/.DS_Store` + untracked clutter stay out).
- **Nominal thrust constant:** 0.025 N (`minfuel_config` default). `thrust_tag(0.025)` must return `''`.
- Commit trailer (exact): `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`
- On any unexpected failure (test fails, nominal path changes, git surprise): STOP, report BLOCKED with exact output.

---

### Task 1: Shared helpers — `cr3bp_fingerprint`, `check_cr3bp_fp`, `thrust_tag`

**Files:**
- Create: `CC/cr3bp_fingerprint.m`, `CC/check_cr3bp_fp.m`, `CC/thrust_tag.m`
- Modify: `CC/cr3bp_lt_params.m` — ADDITIVE input echo `p.thrustN = thrust_N; p.ispS = Isp_s;` (the struct today carries only `.Isp`/no thrust echo; the fingerprint reads `.thrustN/.ispS`, matching the earth `kepler_lt_params` convention; additive ⇒ all existing dependents unaffected). Update its header OUTPUTS comment accordingly. (Plan amendment 2026-07-21 after Task-1 BLOCKED escalation.)
- Test: `CC/test_cr3bp_fp.m`

**Interfaces:**
- Produces: `fp = cr3bp_fingerprint(p, extra)` — `p` a `cr3bp_lt_params` struct; `extra` optional struct merged in; returns struct with `thrustN,m0kg,ispS,Tmax,cEx,muStar,pSund` (+ every `extra` field). `check_cr3bp_fp(Scached, fpNow, file, tag)` — errors id `check_cr3bp_fp:mismatch` on differing shared field; warns ids `:noFingerprint` / `:schemaOlder` otherwise. `tag = thrust_tag(thrustN)` — `''` at 0.025, else e.g. `'_T20mN'` (0.020) / `'_T32p5mN'` (0.0325).

- [ ] **Step 1: Write the failing test** `CC/test_cr3bp_fp.m`:

```matlab
% TEST_CR3BP_FP  Unit tests: fingerprint build/check + thrust_tag.
here = fileparts(mfilename('fullpath'));  addpath(here);
p  = cr3bp_lt_params(0.025, 15, 2100);
fp = cr3bp_fingerprint(p, struct('tf', 7.23, 'insertion', 'campaign'));
assert(fp.thrustN==0.025 && fp.m0kg==15 && fp.ispS==2100, 'core fields');
assert(abs(fp.Tmax - p.Tmax) < 1e-15 && abs(fp.muStar - p.muStar) < 1e-15, 'derived fields');
assert(fp.tf==7.23 && strcmp(fp.insertion,'campaign'), 'extra fields merged');
% (a) match -> silent
S1 = struct('fp', fp);
check_cr3bp_fp(S1, fp, 'file.mat', 'tag');
% (b) legacy (no fp) -> warn, not error
wOld = warning('off','all'); lastwarn('');
check_cr3bp_fp(struct('X',1), fp, 'file.mat', 'tag');
[~, wid] = lastwarn; assert(strcmp(wid,'check_cr3bp_fp:noFingerprint'), 'legacy warn');
% (c) schema-older (cached fp missing a new field) -> warn
fpOld = rmfield(fp, 'insertion');  lastwarn('');
check_cr3bp_fp(struct('fp',fpOld), fp, 'file.mat', 'tag');
[~, wid] = lastwarn; assert(strcmp(wid,'check_cr3bp_fp:schemaOlder'), 'schema warn');
warning(wOld);
% (d) mismatch -> hard error naming the field
fpBad = fp;  fpBad.thrustN = 0.020;
ok = false;
try, check_cr3bp_fp(struct('fp',fpBad), fp, 'file.mat', 'tag');
catch err, ok = strcmp(err.identifier,'check_cr3bp_fp:mismatch') && contains(err.message,'thrustN'); end
assert(ok, 'mismatch must hard-error naming the field');
% (e) thrust_tag
assert(isempty(thrust_tag(0.025)), 'nominal tag must be empty');
assert(strcmp(thrust_tag(0.020), '_T20mN'), '20 mN tag');
assert(strcmp(thrust_tag(0.0325), '_T32p5mN'), 'fractional-mN tag');
fprintf('test_cr3bp_fp: ALL PASS\n');
```

- [ ] **Step 2: Run to verify it fails** —
`/Applications/MATLAB_R2025b.app/bin/matlab -batch "run('/Users/msc/Desktop/optimal_control/orbit_transfer/cr3bp_common/test_cr3bp_fp.m')"`
Expected: FAIL (`cr3bp_fingerprint` undefined).

- [ ] **Step 3: Implement the three helpers.**

`CC/cr3bp_fingerprint.m`:
```matlab
function fp = cr3bp_fingerprint(p, extra)
% CR3BP_FINGERPRINT  Build the config fingerprint that determines a solution.
%
% Captures the physics/config a cached artifact depends on, from the
% cr3bp_lt_params struct (so it cannot drift from the actual physics), plus
% caller-specific extras (tf, insertion, epsMin, ...). Consumed by
% check_cr3bp_fp at every cache read (2026-07-21 review triage C5/C6).
%
% INPUTS:
%   p     - cr3bp_lt_params struct (.thrustN .m0kg .ispS .Tmax .c .muStar) [struct]
%   extra - optional struct of run-specific fields to merge in [struct]
% OUTPUTS:
%   fp - fingerprint struct (.thrustN .m0kg .ispS .Tmax .cEx .muStar .pSund
%        + extras) [struct]
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-21-ladder-prep-design.md sec 2.
if nargin < 2 || isempty(extra), extra = struct(); end
fp = struct('thrustN', p.thrustN, 'm0kg', p.m0kg, 'ispS', p.ispS, ...
            'Tmax', p.Tmax, 'cEx', p.c, 'muStar', p.muStar);
if isfield(p, 'pSund'), fp.pSund = p.pSund; end
fn = fieldnames(extra);
for k = 1:numel(fn), fp.(fn{k}) = extra.(fn{k}); end
end
```

`CC/check_cr3bp_fp.m`:
```matlab
function check_cr3bp_fp(Scached, fpNow, file, tag)
% CHECK_CR3BP_FP  Fail-loud cache-fingerprint guard (earth-campaign pattern).
%
% Compares a loaded cache struct's .fp against the current fingerprint:
%   - no .fp at all       -> WARN (legacy cache, trusted under matching tag)
%   - field only in fpNow -> WARN (schema evolution, compatible)
%   - field on both sides with different values -> ERROR naming field + file
%
% INPUTS:
%   Scached - struct loaded from the cache file [struct]
%   fpNow   - current fingerprint (cr3bp_fingerprint) [struct]
%   file    - cache path, for messages [char]
%   tag     - run tag, for messages [char]
% OUTPUTS: (none) - warns or errors
% REFERENCES:
%   [1] earth_elliptic_to_geo/direct/core/homotopy_mee.m>check_cache_fp (the
%       precedent); [2] spec 2026-07-21-ladder-prep-design.md sec 2.
if ~isfield(Scached, 'fp')
    warning('check_cr3bp_fp:noFingerprint', ...
        '%s has no config fingerprint (legacy cache) -- trusting under tag ''%s''', file, tag);
    return;
end
fn = fieldnames(fpNow);
for k = 1:numel(fn)
    f = fn{k};
    if ~isfield(Scached.fp, f)
        warning('check_cr3bp_fp:schemaOlder', ...
            '%s: fingerprint field ''%s'' absent from cache (schema evolution) -- trusting', file, f);
        continue;
    end
    if ~isequal(Scached.fp.(f), fpNow.(f))
        error('check_cr3bp_fp:mismatch', ...
            ['fingerprint mismatch in %s: field ''%s'' differs from the current ' ...
             'config -- stale/foreign cache under tag ''%s''; delete it or use a new tag'], ...
            file, f, tag);
    end
end
end
```

`CC/thrust_tag.m`:
```matlab
function tag = thrust_tag(thrustN)
% THRUST_TAG  Artifact filename token for a thrust rung.
%
% '' at the nominal 25 mN (all existing cache names stay byte-identical);
% otherwise '_T<mN>mN' with '.' -> 'p' (0.020 -> '_T20mN', 0.0325 -> '_T32p5mN').
%
% INPUTS:  thrustN - thrust [N, scalar]
% OUTPUTS: tag     - filename token [char]
% REFERENCES: [1] spec 2026-07-21-ladder-prep-design.md sec 2.
if abs(thrustN - 0.025) < 1e-12
    tag = '';
else
    mN = thrustN * 1000;
    s  = strrep(sprintf('%g', mN), '.', 'p');
    tag = sprintf('_T%smN', s);
end
end
```

- [ ] **Step 4: Run test to verify PASS.** Same command; expect `test_cr3bp_fp: ALL PASS`.
- [ ] **Step 5: Commit**
```bash
cd /Users/msc/Desktop/optimal_control
git add orbit_transfer/cr3bp_common/cr3bp_fingerprint.m orbit_transfer/cr3bp_common/check_cr3bp_fp.m orbit_transfer/cr3bp_common/thrust_tag.m orbit_transfer/cr3bp_common/test_cr3bp_fp.m
git commit -m "feat(cr3bp_common): fingerprint build/check + thrust_tag helpers (ladder-prep T1)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: `minfuel_config` override + nominal byte-path guard

**Files:**
- Modify: `CC/minfuel_config.m` (signature `function cfg = minfuel_config()` → optional override)
- Test: `CC/test_minfuel_config_override.m`

**Interfaces:**
- Produces: `cfg = minfuel_config()` — byte-identical to today. `cfg = minfuel_config(struct('thrustN',0.020))` — same struct with the named fields overridden AFTER all defaults are built (so derived fields like schedules stay; `thrustN/m0kg/ispS` are pure scalars consumed later by `cr3bp_lt_params`, safe to override).

- [ ] **Step 1: Write the failing test** `CC/test_minfuel_config_override.m`:
```matlab
% TEST_MINFUEL_CONFIG_OVERRIDE  Override arg + nominal invariance.
here = fileparts(mfilename('fullpath'));  addpath(here);
c0 = minfuel_config();
c1 = minfuel_config(struct('thrustN', 0.020));
assert(c1.thrustN == 0.020, 'override applied');
c1b = rmfield(c1, 'thrustN');  c0b = rmfield(c0, 'thrustN');
assert(isequal(c1b, c0b), 'override must change ONLY the named field');
c2 = minfuel_config(struct());
assert(isequal(c2, c0), 'empty override = default');
assert(isequal(minfuel_config(), c0), 'no-arg call unchanged');
fprintf('test_minfuel_config_override: ALL PASS\n');
```
- [ ] **Step 2: Run — FAIL** (too many input arguments).
- [ ] **Step 3: Implement.** In `CC/minfuel_config.m`, change line 1 to `function cfg = minfuel_config(over)` and append before the final `end`:
```matlab
% --- optional override (ladder rungs): merge caller fields over the defaults ---
if nargin >= 1 && ~isempty(over)
    fo = fieldnames(over);
    for ko = 1:numel(fo), cfg.(fo{ko}) = over.(fo{ko}); end
end
```
(If the file has no trailing `end` because it is a script-style function, append at the bottom.) Also update the header comment: add an INPUTS block documenting `over`.
- [ ] **Step 4: Run — PASS.**
- [ ] **Step 5: Commit** (`git add orbit_transfer/cr3bp_common/minfuel_config.m orbit_transfer/cr3bp_common/test_minfuel_config_override.m`, message `feat(cr3bp_common): optional config override for thrust rungs (ladder-prep T2)` + trailer).

---

### Task 3: Tulip engine — solver opts + boundSat, fp wiring, chain helper

**Files:**
- Modify: `SMF/casadi_minfuel_sundman.m` (trailing `opts` + boundSat), `SMF/minfuel_at_tf.m`, `SMF/sundman_homotopy.m`, `SMF/gen_tulip_energy_2p.m`, `SMF/gen_tulip_mintime.m` (fp wiring)
- Create: `SMF/chain_rung_seed_tulip.m`
- Test: `SMF/test_ladder_prep_tulip.m`

**Interfaces:**
- `casadi_minfuel_sundman(..., warmTight, opts)` — NEW optional 15th arg `opts` struct: `.vBox` (velocity half-width, default 12), `.rBox` (position half-width, default 3). `out` gains `.boundSat = struct('minSlack', s, 'worst', label, 'hit', logical)` always.
- `[sigma, X0, U0, tauf0, fp] = chain_rung_seed_tulip(solSource, tfNew, pNew, extraFp)` — `solSource`: a solver `out` struct (`.X` 8×M with t=row 8, `.U` 4×M) or a saved-file path whose `out` field holds one, plus the source `fp` if present; errors id `chain_rung_seed_tulip:sameThrust` if source fp thrust equals `pNew.thrustN`.

- [ ] **Step 1: Write the failing test** `SMF/test_ladder_prep_tulip.m`:
```matlab
% TEST_LADDER_PREP_TULIP  boundSat fields + opts back-compat + chain helper.
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
addpath(fullfile(here,'..','..','..','cr3bp_common'));  setup_cr3bp_common();
cfg = minfuel_config();  p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
assert(isfile('sundman_minfuel_certified.mat'), 'need the certified cache');
C = load('sundman_minfuel_certified.mat');
% (a) chain helper from the certified solution to a 20 mN config
p20 = cr3bp_lt_params(0.020, cfg.m0kg, cfg.ispS);
tfNew = C.out.X(8,end);                     % same tf, new thrust (pilot pattern)
[sg, X0, U0, tauf0, fp] = chain_rung_seed_tulip(C, tfNew, p20, struct('note','test'));
assert(size(X0,1)==8 && size(U0,1)==4 && numel(sg)==size(X0,2), 'shapes');
assert(abs(X0(8,end)-tfNew) < 1e-9 && abs(X0(8,1)) < 1e-12, 'time row endpoints');
assert(tauf0 > 0 && issorted(sg) && abs(sg(1))<1e-12 && abs(sg(end)-1)<1e-12, 'sigma/tauf0');
assert(fp.thrustN==0.020 && isfield(fp,'chainedFrom'), 'fp thrust + provenance');
assert(max(abs(X0(1:6,1).'  - C.rv0(:).')) < 1e-9, 'rv0 pinned');
assert(max(abs(X0(1:6,end).'- C.rvf(:).')) < 1e-9, 'rvf pinned');
% (b) same-thrust chain refused
ok = false;
try, chain_rung_seed_tulip(C, tfNew, p, struct());
catch err, ok = strcmp(err.identifier,'chain_rung_seed_tulip:sameThrust'); end
assert(ok, 'same-thrust chain must error');
% (c) solver back-compat + boundSat: 3-iter probe, 14-arg call == 15-arg default
o14 = casadi_minfuel_sundman(C.sigma, tfNew, C.rv0, C.rvf, p.Tmax, p.c, p.muStar, ...
        C.out.X, C.out.U, C.tauf0, cfg.pSund, 3, 0, true);
o15 = casadi_minfuel_sundman(C.sigma, tfNew, C.rv0, C.rvf, p.Tmax, p.c, p.muStar, ...
        C.out.X, C.out.U, C.tauf0, cfg.pSund, 3, 0, true, struct());
assert(isfield(o14,'boundSat') && isfield(o14.boundSat,'minSlack') && ...
       isfield(o14.boundSat,'worst') && islogical(o14.boundSat.hit), 'boundSat fields');
assert(max(abs(o14.X(:)-o15.X(:))) < 1e-12, 'empty opts must be byte-compatible');
fprintf('test_ladder_prep_tulip: ALL PASS\n');
```
- [ ] **Step 2: Run — FAIL** (`chain_rung_seed_tulip` undefined).
- [ ] **Step 3a: Solver edit.** In `SMF/casadi_minfuel_sundman.m`: extend the signature `..., epsilon, warmTight, opts)`; right after the existing arg handling add:
```matlab
if nargin < 15 || isempty(opts), opts = struct(); end
vBox = 12;  if isfield(opts,'vBox') && ~isempty(opts.vBox), vBox = opts.vBox; end
rBox = 3;   if isfield(opts,'rBox') && ~isempty(opts.rBox), rBox = opts.rBox; end
```
Replace the bounds block (lines ~114–115) with:
```matlab
lbX = repmat([-rBox;-rBox;-rBox;-vBox;-vBox;-vBox;0.3;0], 1, nN);
ubX = repmat([ rBox; rBox; rBox; vBox; vBox; vBox;1.0; 2*tf], 1, nN);
```
After extraction (`Xs` available, before the `out = struct(...)`) add the boundSat diagnostic:
```matlab
% Bound-saturation diagnostic (2026-07-21 triage C4; output-only). Nonphysical
% boxes checked at INTERIOR nodes (BCs pin the endpoints by construction).
Xi = Xs(:,2:end-1);
slk = [ rBox - max(abs(Xi(1:3,:)),[],'all');            % position box
        vBox - max(abs(Xi(4:6,:)),[],'all');            % velocity box
        min(Xi(7,:),[],'all') - 0.3;                    % mass lower
        1.0 - max(Xi(7,:),[],'all') ];                  % mass upper
lbl = {'rBox','vBox','massLo','massHi'};
[minSlack, iw] = min(slk);
boundSat = struct('minSlack', minSlack, 'worst', lbl{iw}, 'hit', minSlack < 1e-4);
if boundSat.hit
    warning('casadi_minfuel_sundman:boundSaturation', ...
        'nonphysical box ''%s'' within %.2g of binding -- widen via opts before trusting', ...
        lbl{iw}, max(minSlack,0));
end
```
and add `'boundSat', boundSat, ...` into the `out = struct(...)` build.
- [ ] **Step 3b: Chain helper** `SMF/chain_rung_seed_tulip.m`:
```matlab
function [sigma, X0, U0, tauf0, fp] = chain_rung_seed_tulip(src, tfNew, pNew, extraFp)
% CHAIN_RUNG_SEED_TULIP  Alias-free cross-rung warm start for the sundman engine.
%
% Takes a previous rung's CONVERGED Sundman solution, rescales ONLY the time
% row to tfNew (same-mesh scalar rescale -- controls stay attached to their own
% spatial nodes; no sigma- or t-interpolation anywhere, so the phase-aliasing
% channel never opens), and re-maps through the house no-resample
% sundman_seed_map: fresh tauf0 for the new rung (never reuse the source's),
% endpoints pinned exactly. (2026-07-21 triage C5/C6 ladder-prep.)
%
% INPUTS:
%   src     - loaded cache struct with .out (solver struct, X 8xM with t=row 8,
%             U 4xM), .rv0 [1x6], .rvf [1x6], and optionally .fp  [struct]
%   tfNew   - target transfer time for the new rung [ND, scalar]
%   pNew    - cr3bp_lt_params struct for the NEW rung (thrust differs!) [struct]
%   extraFp - optional extra fingerprint fields (e.g. .note, .epsMin) [struct]
% OUTPUTS:
%   sigma [M'x1], X0 [8xM'], U0 [4xM'], tauf0 [scalar] - warm start for
%   casadi_minfuel_sundman at the new rung;  fp - fingerprint with
%   .chainedFrom provenance [struct]
% REFERENCES:
%   [1] sundman_seed_map.m (the no-resample map); [2] spec sec 4.
if nargin < 4 || isempty(extraFp), extraFp = struct(); end
srcThrust = NaN;
if isfield(src,'fp') && isfield(src.fp,'thrustN'), srcThrust = src.fp.thrustN; end
if isnan(srcThrust), srcThrust = 0.025; end          % legacy caches are nominal
assert(abs(srcThrust - pNew.thrustN) > 1e-12, 'chain_rung_seed_tulip:sameThrust', ...
    'source and target thrust are both %.4g N -- chaining to the same rung is a caller bug', srcThrust);
X = src.out.X;  U = src.out.U;
assert(size(X,1) >= 8, 'chain_rung_seed_tulip:badState', 'need 8-state X with t=row 8');
X(8,:) = X(8,:) * (tfNew / X(8,end));                % time row only; spatial mesh untouched
cfg = minfuel_config();
[sigma, X0, U0, tauf0] = sundman_seed_map(X(1:7,:), U, tfNew, X(8,:).', ...
                                          cfg.pSund, pNew.muStar, src.rv0, src.rvf);
extraFp.chainedFrom = sprintf('T=%.4gN', srcThrust);
extraFp.tf = tfNew;
fp = cr3bp_fingerprint(pNew, extraFp);
end
```
NOTE to implementer: check `sundman_seed_map`'s exact input contract (`Xseed [7xM]` `[r;v;m]`, `sgNorm` = the seed's own node parameter — the rescaled time row is the natural monotone choice). If it wants `[r;v;m]` only (7 rows), pass `X(1:7,:)` as written; run the test's stencil/shape asserts to confirm.
- [ ] **Step 3c: fp wiring (tulip).** Pattern at every save: build `fp = cr3bp_fingerprint(p, struct('tf',tf, ...))` and append `'fp'` to the `save(...)` list. At every load-for-warm-start: `check_cr3bp_fp(loadedStruct, fp, file, tag)`. Apply at: `minfuel_at_tf.m` (energy-seed load line ~52 + neighbor load ~57; result save line ~133), `sundman_homotopy.m` (saveFile save), `gen_tulip_energy_2p.m` (ckpt save/load + outFile save), `gen_tulip_mintime.m` (anchor save). All these files already have `p`/`cfg` in scope; `cr3bp_common` is on their path via `setup_paths`.
- [ ] **Step 4: Run the test — PASS** (`test_ladder_prep_tulip: ALL PASS`). Also rerun the existing no-solve guardrails: `matlab -batch "cd('<SMF>'); setup_paths; test_minfuel_lib"` — expect unchanged PASS.
- [ ] **Step 5: Commit** the six SMF files + test, message `feat(tulip): boundSat + opt-in boxes, cache fingerprints, cross-rung chain helper (ladder-prep T3)` + trailer.

---

### Task 4: ELFO engine — boundSat, fp wiring + seed filter, cBox scaling, chain helper

**Files:**
- Modify: `ELFO/casadi_energy_freetf.m`, `ELFO/casadi_mintime_freetf.m` (boundSat), `ELFO/elfo_find_energy_seed.m` (fp filter), `ELFO/gen_elfo_minfuel.m`, `ELFO/gen_elfo_energy_gravhom.m`, `ELFO/gen_elfo_energy_tfsweep.m`, `ELFO/elfo_run_one.m` (fp wiring + cBox scaling + thrust threading)
- Create: `ELFO/chain_rung_seed_elfo.m`
- Test: `ELFO/test_ladder_prep_elfo.m`

**Interfaces:**
- `out.boundSat` on both freetf solvers, same fields as Task 3 (+ a `cBox` slack row: `min(cScale-cBox(1), cBox(2)-cScale)`).
- `elfo_find_energy_seed(resDir, tfTarget, relTol, fpNow)` — NEW optional 4th arg: when given, seeds whose stored `fp.Tmax` mismatches `fpNow.Tmax` are skipped; seeds with no fp are eligible (legacy) with one aggregate warning.
- `[seedS, fp] = chain_rung_seed_elfo(seedS, pNew, extraFp)` — pass-through of `X,U,tauf0,sigma,rv0,rvf` (cScale decouples the clock; same nodes, no interp), same-thrust refusal, fp with `.chainedFrom`.
- Drivers accept `opts.thrustN` (default nominal) → `cfg = minfuel_config(struct('thrustN',opts.thrustN))`, `thrust_tag` appended to artifact names, `base.cBox` scaled: `Tfac = 0.025/thrustN; cBox = [0.15*min(1,Tfac), 6*max(1,Tfac)]` (tfsweep uses its own `[0.10 8]` nominal pair, same rule).

- [ ] **Step 1: Write the failing test** `ELFO/test_ladder_prep_elfo.m`:
```matlab
% TEST_LADDER_PREP_ELFO  chain helper, seed fp filter, cBox rule (no solves).
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
cfg = minfuel_config();  p20 = cr3bp_lt_params(0.020, cfg.m0kg, cfg.ispS);
% (a) chain helper pass-through + fp
S = load(fullfile(here,'results','energy_elfo_f1200.mat'));
[S2, fp] = chain_rung_seed_elfo(S, p20, struct('note','test'));
assert(isequal(S2.X, S.X) && isequal(S2.U, S.U), 'pass-through must not touch X/U');
assert(fp.thrustN==0.020 && isfield(fp,'chainedFrom'), 'fp');
ok=false; p25 = cr3bp_lt_params(0.025, cfg.m0kg, cfg.ispS);
try, chain_rung_seed_elfo(S, p25, struct());
catch err, ok = strcmp(err.identifier,'chain_rung_seed_elfo:sameThrust'); end
assert(ok, 'same-thrust refusal');
% (b) seed fp filter: legacy seeds eligible under a 25 mN fp, skipped under 20 mN
fp25 = cr3bp_fingerprint(p25);  fp20 = cr3bp_fingerprint(p20);
w = warning('off','all');
[sfA,~,~] = elfo_find_energy_seed(fullfile(here,'results'), S.X(8,end), 0.02, fp25);
[sfB,~,~] = elfo_find_energy_seed(fullfile(here,'results'), S.X(8,end), 0.02, fp20);
warning(w);
assert(~isempty(sfA), 'legacy seed eligible under nominal fp');
assert(isempty(sfB), 'legacy (nominal) seed must NOT satisfy a 20 mN request');
fprintf('test_ladder_prep_elfo: ALL PASS\n');
```
NOTE: (b) encodes a POLICY decision — under an off-nominal fp, legacy (fingerprint-less) seeds are treated as NOMINAL and therefore skipped; only the no-fp-given call keeps full legacy behavior. Implement `elfo_find_energy_seed` accordingly: legacy seed + `fpNow.thrustN==0.025` → eligible (warn once); legacy seed + off-nominal fpNow → skip.
- [ ] **Step 2: Run — FAIL.**
- [ ] **Step 3: Implement.** (a) boundSat in both freetf solvers: mirror Task 3's block, with rows `rBox=3` analog per their own `lbX/ubX` values (read the files; use their actual box constants), plus `cBoxSlack = min(Xs(9,end)-cBox(1), cBox(2)-Xs(9,end))` labeled `'cBox'`; add to `out`. (b) `chain_rung_seed_elfo.m` (complete):
```matlab
function [seedS, fp] = chain_rung_seed_elfo(seedS, pNew, extraFp)
% CHAIN_RUNG_SEED_ELFO  Cross-rung warm start for the freetf engine (trivial by
% design: the cScale slack state decouples the clock, so a thrust rung reuses
% the source X,U,tauf0 on the SAME nodes -- no interpolation, no aliasing).
% Adds the new rung's fingerprint + chainedFrom provenance; refuses a
% same-thrust chain. (2026-07-21 ladder-prep, spec sec 4.)
%
% INPUTS:
%   seedS   - loaded seed struct (.X [9xM] or [8xM], .U [4xM], .tauf0, .sigma,
%             .rv0, .rvf, optional .fp) [struct]
%   pNew    - cr3bp_lt_params for the NEW rung [struct]
%   extraFp - optional extra fingerprint fields [struct]
% OUTPUTS:
%   seedS - unchanged seed struct (pass-through); fp - new-rung fingerprint
% REFERENCES: [1] casadi_energy_freetf.m (cScale mechanics); [2] spec sec 4.
if nargin < 3 || isempty(extraFp), extraFp = struct(); end
srcThrust = 0.025;                                    % legacy caches are nominal
if isfield(seedS,'fp') && isfield(seedS.fp,'thrustN'), srcThrust = seedS.fp.thrustN; end
assert(abs(srcThrust - pNew.thrustN) > 1e-12, 'chain_rung_seed_elfo:sameThrust', ...
    'source and target thrust are both %.4g N', srcThrust);
extraFp.chainedFrom = sprintf('T=%.4gN', srcThrust);
if isfield(seedS,'X'), extraFp.tf = seedS.X(8,end); end
fp = cr3bp_fingerprint(pNew, extraFp);
end
```
(c) `elfo_find_energy_seed` 4th arg per the Step-1 NOTE (load `fp` alongside `X` in the scan; apply the policy). (d) fp wiring: every ELFO save adds `'fp'`; `gen_elfo_minfuel`'s checkpoint-resume check becomes `abs(C.tf0-tf0)<1e-9 && <fp check via check_cr3bp_fp>`; row saves in `elfo_run_one` include fp. (e) thrust threading + tags + cBox rule in the three drivers (`opts.thrustN` → config override → `p`; `thrust_tag(thrustN)` appended to `tag`/filenames; `cBox` scaled as specified — at nominal `Tfac=1` reproduces `[0.15 6]`/`[0.10 8]` exactly).
- [ ] **Step 4: Run the test — PASS.** Also re-run `test_cr3bp_fp` (unchanged PASS).
- [ ] **Step 5: Commit** the 8 ELFO files + test, message `feat(elfo): boundSat, cache fingerprints + seed fp filter, cBox scaling, chain helper (ladder-prep T4)` + trailer.

---

### Task 5: Folded review fixes (C4/C6): rF cascades, anchor gate, tf-sweep re-clean

**Files:**
- Modify: `SMF/gen_tulip_energy_2p.m`, `ELFO/gen_elfo_minfuel.m` (rF-after-rT cascades), `ELFO/gen_elfo_energy_tfsweep.m` (tight re-clean before banking), `ELFO/gen_elfo_mintime.m` (anchor certification gate)

- [ ] **Step 1: rF cascades.** In both step_solve helpers, the current flow is loose `rL` → (fail) tight-from-Xk `rF` → re-clean `rT`; when `rT` fails after a SUCCESSFUL loose probe, the step currently dies without trying tight-from-Xk. Change the `rT` failure branch (both files, same pattern — `gen_elfo_minfuel.m` step_solve tail shown; mirror in `gen_tulip_energy_2p.m`):
```matlab
if strcmp(rT.ipoptStatus,'Solve_Succeeded') && rT.maxDefect < 1e-6
    Xn = rT.X;  Un = rT.U;  ok = true;  info = rT;
else
    % re-clean failed after a good loose probe: try tight-from-Xk before
    % declaring the step dead (2026-07-21 triage C6 fallback cascade)
    oF2 = base;  oF2.maxIter = ctx.maxIter;  oF2.warmTight = true;
    rF2 = casadi_energy_freetf(ctx.sigma,ctx.rv0,ctx.rvf,ctx.Tmax,ctx.cEx,ctx.muStar,Xk,Uk,ctx.tauf0,oF2);
    if strcmp(rF2.ipoptStatus,'Solve_Succeeded') && rF2.maxDefect < 1e-6
        Xn = rF2.X;  Un = rF2.U;  ok = true;  info = rF2;
    else
        ok = false;  Xn = Xk;  Un = Uk;  info = rT;
    end
end
```
(For `gen_tulip_energy_2p.m` use its own solver call/argument pattern — read the file's step_solve and mirror the same branch structure with `casadi_energy_freetf`→its solver.)
- [ ] **Step 2: tf-sweep re-clean.** In `ELFO/gen_elfo_energy_tfsweep.m>solve_tf`, after the successful LOOSE branch (`rL` clean → currently returns immediately), insert a tight re-clean before returning: solve again with `maxIter=ctx.maxIter, warmTight=true` from `rL.X/rL.U`; return the re-cleaned result if it passes the gate, else fall through to the existing `rF` branch.
- [ ] **Step 3: anchor gate.** In `ELFO/gen_elfo_mintime.m`, wrap the `save(outFile, ...)` (line ~73) in:
```matlab
assert(strcmp(out.ipoptStatus,'Solve_Succeeded') && out.maxDefect < 1e-6, ...
    'gen_elfo_mintime:uncertified', ...
    'min-time anchor NOT certified (status=%s, defect=%.2g) -- refusing to save', ...
    out.ipoptStatus, out.maxDefect);
```
(also add `'fp'` to that save with `fp = cr3bp_fingerprint(p, struct('tf',out.tf,'insertion',insMeta.label))`).
- [ ] **Step 4: Lint all four** (`checkcode` — 0 parse errors) and re-run `test_ladder_prep_elfo` (PASS unchanged).
- [ ] **Step 5: Commit**, message `fix(gto): rF fallback cascades, tf-sweep re-clean, certified-only min-time anchor (ladder-prep T5)` + trailer.

---

### Task 6: Nominal regression + 20 mN pilot rungs + docs

**Files:**
- Create: `SMF/pilot_rung_20mN.m`, `ELFO/pilot_rung_20mN.m`
- Modify: `GTO_tulip/TODO.md`, `GTO_ELFO/TODO.md`, both `doc/reviews/2026-07-21_triage.md` (mark C4–C6 status), spec status line

**Interfaces:**
- Each pilot: chain from the named nominal source → energy re-clean at 20 mN → ε:1→0 sharpen (hardened gates) → save `_T20mN`-tagged artifacts + print PASS/FAIL line. Resume-safe (per-step caches under `results/`), background-runnable.

- [ ] **Step 1: Nominal byte-path regression (gate before pilots).** Run the fast no-solve suites: `test_cr3bp_fp`, `test_minfuel_config_override`, `test_ladder_prep_tulip`, `test_ladder_prep_elfo`, `test_minfuel_lib` — all PASS. Then one nominal smoke with defaults to confirm no numeric drift: `matlab -batch "cd('<ELFO>'); setup_paths; smoke_fixedtf"` — expect its historical PASS output.
- [ ] **Step 2: Write the two pilots.** `SMF/pilot_rung_20mN.m` (complete):
```matlab
function best = pilot_rung_20mN()
% PILOT_RUNG_20MN  Ladder-prep validation: one warm-chained 20 mN fuel rung.
%
% Chains from the certified nominal 1.15x solution (same t_f, new thrust) via
% chain_rung_seed_tulip, re-cleans the energy problem at 20 mN, then sharpens
% eps 1->0 through sundman_homotopy's hardened gates. PASS = certified=1 +
% clean boundSat + fp recorded. Artifacts under _T20mN tags; certified caches
% untouched. (2026-07-21 ladder-prep T6; spec sec 6.)
%
% OUTPUTS: best - sundman_homotopy best struct (.certified .epsReached ...)
% REFERENCES: [1] spec 2026-07-21-ladder-prep-design.md sec 6.
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
addpath(fullfile(here,'..','..','..','cr3bp_common'));  setup_cr3bp_common();
cfg = minfuel_config(struct('thrustN', 0.020));
p   = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
C   = load('sundman_minfuel_certified.mat');
tf  = C.out.X(8,end);
[sg, X0, U0, tauf0, fp] = chain_rung_seed_tulip(C, tf, p, struct('pilot','20mN')); %#ok<ASGLU>
tag = thrust_tag(cfg.thrustN);
saveFile = fullfile(here, 'results', sprintf('pilot_minfuel%s.mat', tag));
sched = [1 0.6 0.35 0.2 0.12 0.07 0.04 0.02 0.01 0.005 0.002 0.001 0];
[best, tbl] = sundman_homotopy(p, C.rv0, C.rvf, sg, X0, U0, tauf0, cfg.pSund, ...
                               sched, 3000, saveFile); %#ok<ASGLU>
sat = 'n/a'; if isfield(best,'boundSat'), sat = best.boundSat.worst; end
fprintf('\nPILOT 20mN TULIP: certified=%d epsReached=%.4g defect=%.2g sw=%d boundSatWorst=%s\n', ...
        best.certified, best.epsReached, best.maxDefect, best.switches, sat);
end
```
`ELFO/pilot_rung_20mN.m`: same shape on the freetf engine — `chain_rung_seed_elfo` from `results/energy_elfo_f1200.mat`, `cfg=minfuel_config(struct('thrustN',0.020))`, energy re-clean via `casadi_energy_freetf` (scaled `cBox=[0.15*min(1,1.25) 6*max(1,1.25)]=[0.15 7.5]`, `tfTarget=S.X(8,end)`), then `gen_elfo_minfuel(struct('seedFile',<recleaned 20mN seed saved with fp>,'target','ELFO','outFile',<pilot_T20mN path>,'maxIter',3000))` — write it following the same structure and printing the same PASS line (`PILOT 20mN ELFO: ...`).
NOTE: the energy re-clean at the new thrust must be gated (`Solve_Succeeded` + defect<1e-6) before the sharpening is attempted; if the re-clean itself fails, print `PILOT ... BLOCKED at energy re-clean` and return — that is a finding, not a crash.
- [ ] **Step 3: Launch both pilots in background** (separate processes, logs under the scratchpad), monitor to completion. PASS bar per spec: `certified=1`, `boundSat` clean, fp present in saved artifacts. **A pilot that honestly fails does NOT revert the package** — record the wall + attempt trajectory in the campaign TODO (house honesty rule) and report BLOCKED-as-finding.
- [ ] **Step 4: Docs.** Update both TODOs (ladder-prep trio items → done, with commit refs + pilot outcomes; add the deep-rung/free-span escalation note from the spec), both triage docs (C4–C6 status lines), and the spec's Status line (→ implemented). 
- [ ] **Step 5: Commit + push.**
```bash
cd /Users/msc/Desktop/optimal_control
git add orbit_transfer/GTO_tulip/direct/sundman_minfuel/pilot_rung_20mN.m orbit_transfer/GTO_ELFO/direct/elfo/pilot_rung_20mN.m orbit_transfer/GTO_tulip/TODO.md orbit_transfer/GTO_ELFO/TODO.md orbit_transfer/GTO_tulip/doc/reviews/2026-07-21_triage.md orbit_transfer/GTO_ELFO/doc/reviews/2026-07-21_triage.md docs/superpowers/specs/2026-07-21-ladder-prep-design.md
git commit -m "feat(gto): 20 mN pilot rungs + ladder-prep close-out docs (ladder-prep T6)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
git push origin main
```

---

## Self-Review notes

- **Spec coverage:** §2 fingerprints/tag/config-override → T1+T2; §3 boundSat + opt-in boxes + cBox rule → T3(a)/T4(a,e); §4 chain helpers (tulip re-map, ELFO pass-through, same-thrust refusal, chainedFrom) → T3(b)/T4(b); §5 folded fixes → T5; §6 validation (unit tests, nominal regression, two pilots + honest-failure rule) → T6; §7 file inventory matches the tasks; §8 risks → back-compat asserts in T3 test (o14==o15), read-only certified caches (pilots write `_T20mN` names only), legacy-cache policy encoded in T4 Step 1 NOTE.
- **Type consistency:** `fp` struct fields identical across producers/consumers; `check_cr3bp_fp(Scached, fpNow, file, tag)` signature used verbatim at every wiring site; both chain helpers share the same-thrust error-id naming pattern; `boundSat` fields (`minSlack/worst/hit`) identical on all three solvers.
- **Known judgment points for the implementer (explicitly delegated):** exact `sundman_seed_map` input contract (T3 Step 3b NOTE — the test's asserts are the acceptance); freetf solvers' actual box constants for their boundSat rows (T4 Step 3a says read the file); `gen_tulip_energy_2p`'s solver-call pattern for its rF cascade (T5 Step 1). Each is bounded by a named test or lint gate.
- **Placeholder scan:** all new files carry complete code; modification steps carry either exact replacement blocks or a named pattern + insertion code + file/line anchors.
