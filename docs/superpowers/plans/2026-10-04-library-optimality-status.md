# Optimality Status and Alternative Transfers — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Every entry of the 70 mN DRO -> tulip costate library carries an optimality status (4 sufficient / 3 necessary only / 2 conjugate point found / 1 neither) and its junction states, and the library holds the non-primary transfers the pipeline finds, so a mission designer can find, fly and judge any recorded transfer.

**Architecture:** One pure classifier (`optimality_status`) defines the status from a `certify_root` result; `certify_root` stamps it on every result (and runs the hypothesis gates after a conjugate FAIL, for the record). Builders keep every candidate above the floor; the packager writes per-cell `status`, per-entry `status_reason` + `junctions`, and a top-level `alternatives` table beside the unchanged primary grid. A new audit unit re-certifies every alternative and requires its status to reproduce. A backfill turns the current record into `library_70mN_24x48_v2` without a rebuild.

**Tech Stack:** MATLAB R2026a (`/Applications/MATLAB_R2026a.app/bin/matlab`), pumpkyn (`pumpkyn.cr3bp.tfMinProp`, `tfMinEoM`), Parallel Computing Toolbox for fences (`capped_pool`, `run_capped`). Tests are function files returning `ok` and printing PASS/FAIL per check (the folder's pattern).

**Spec:** `docs/superpowers/specs/2026-10-04-library-optimality-status-design.md` (read it first).

**Paths used below** (all under `/Users/msc/Desktop/optimal_control/orbit_transfer/`):
`IND` = `DRO_tulip/indirect/`, `CC` = `costate_common/`, `RES` = `DRO_tulip/indirect/results/`.

## Deviations from the spec (found while reading the code; flagged to Mike with the plan)

1. **`optimality_status` lives in `IND`, beside `certify_root`**, not in `costate_common`: it reads `certify_root`'s result struct, which is DRO_tulip's certifier (the costate_common admission rule is "generic by construction").
2. **`fly_transfer` already exists** (`CC/fly_transfer.m`, single-shot flight from z8). The spec's junction-based thrust program is a NEW function, `entry_thrust_program`.
3. **Harvested conjugate refusals must be RE-CERTIFIED, not just gated.** Legacy results store `C.conj = it.conj.pass`, so `conj == 0` means FAIL **or** UNDETERMINED; only a re-run can tell. The backfill re-certifies them (Task 11). Going forward `certify_root` stores the verdict string (`C.conjVerdict`).
4. **`rho <= 1e-6` refusals are below the floor**, not status 1: `certify_crossing` returns them before any normal-chart root exists (no z8, no flight).
5. **Status 2 also requires H6**, the conjugate instrument's own validity condition (`certify_root` step 6), besides the hypothesis gates of step 5.
6. **The audit is a new unit, `audit_status_layer`**, beside the unchanged `audit_phase_catalog`: primaries keep their full audit; the new unit checks primaries' status + junctions and re-certifies alternatives. Adoption requires both clean.

## Global Constraints

- MATLAB header on every new function: `%% Purpose:` / `%% Inputs:` (name, type, size) / `%% Outputs:` / `%% Revision History:` with `%  M. Casey ... (c) 10/04/2026` and `%  Copyright Coorbital Inc.`, then `%% ------------------------ Begin Code Sequence ---------------------------` (copy the layout of `IND/rib_as_catalog.m`).
- Never `i` or `j` as variables; never `%#ok` pragmas; in NEW code write vector lengths as `sqrt(sum(x.^2))`, not `norm`.
- Every solver/flight call outside a unit test runs in a `matlab -batch` process, fenced (`capped_pool` + `run_capped` via the existing callers), never in the shared interactive session.
- Tests that solve run as: `/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('<tests dir>'); exit(~<test>())"` (startup.m changes directory first; the `cd` inside the command is required). Pure tests may use the MATLAB MCP `evaluate_matlab_code`.
- Status codes are exactly: `4` sufficient, `3` necessary only, `2` conjugate point found, `1` neither, `0` empty cell, `-1` below the floor (never stored).
- Floor: a converged normal-chart root (finite z8) whose flight is admissible and arrives within `flyKm <= 100` and `flyVms <= 10`. Overridden (test-seam) results are never packaged.
- `sameRoot` rule (one definition, `IND/same_root.m`): the seven costates within 1e-6 of the larger norm (symmetric) and t_f within 1e-6 relative.
- Primaries of the backfilled library must be bit-identical to the record's (`catalog_content_key` unchanged).
- Each task ends with ONE commit made by the host after review (subagents never commit), using `/Library/Developer/CommandLineTools/usr/bin/git` from `/Users/msc/Desktop/optimal_control`, message ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. A legacy candidate whose `reason` text matches no known pattern (wording drifted between builds) must land in the lowest recordable tier with `inferred = true`, never in 4 — pinned in Task 1.
2. Catalogs WITHOUT the new fields (the 24 x 24 reference, the current record) passed to `transfers_at`, `merge_phase_catalogs` or `audit_status_layer` must work (primaries reported as "not annotated") or refuse by name — pinned in Tasks 8 and 10.
3. Junction arrays with K other than 24 must fly correctly — `entry_thrust_program` reads K from the array; pinned in Task 4 with K = 12.
4. A query at the phase wrap (s_D = 1 - 1e-9 vs 0) must find the s_D = 0 transfer — pinned in Task 10.
5. A re-certification that converges to a DIFFERENT root must be reported as moved, not silently relabelled — pinned in Task 9 (audit) and Task 11 (recertify).

---

### Task 1: The classifier `optimality_status` and the legend `status_key`

**Files:**
- Create: `IND/optimality_status.m`, `IND/status_key.m`
- Test: `IND/tests/test_optimality_status.m`

**Interfaces:**
- Produces: `[code, reason, inferred] = optimality_status(C, opts)` — `C` a `certify_root` result (new or legacy); `opts.gateKm` [100], `opts.gateVms` [10]; `code` double in {-1,1,2,3,4}; `reason` char; `inferred` logical (true when `C` has no `.stage` field).
- Produces: `K = status_key()` — struct with `.codes` [1x6 double] = [4 3 2 1 0 -1], `.names` {1x6 cellstr} = {'sufficient','necessary only','conjugate point found','neither','empty','below floor'}, `.meaning` {1x6 cellstr}, `.order` char (the certifier's check order), `.note` char.
- New-result fields read: `.ok .z .flyKm .flyVms .stage .conjFound .hypAfterConj .reason` (added to certify_root in Task 2).

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_optimality_status()
% TEST_OPTIMALITY_STATUS  The classifier maps every outcome of the gate stack
% to its status (spec 4.2), infers legacy records conservatively (4.3), and
% the legend names the codes.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
z = [(1:7).'; 4.1];
base = struct('ok', false, 'z', z, 'flyKm', 0.1, 'flyVms', 0.01, 'reason', '', ...
              'stage', 0, 'conjFound', false, 'hypAfterConj', '');
mk = @(varargin) setf(base, varargin{:});

% ---- new (stamped) results ------------------------------------------------
[c, r] = optimality_status(mk('ok', true, 'stage', 8, 'reason', 'certified'));
ok = chk(ok, c == 4 && strcmp(r, 'full stack passed'), 'certified -> 4 sufficient');
c = optimality_status(mk('z', nan(8, 1), 'reason', 'normal-chart polish did not converge'));
ok = chk(ok, c == -1, 'no converged root -> below the floor');
c = optimality_status(mk('stage', 1, 'flyKm', NaN, 'reason', 'flight inadmissible: hits the Moon'));
ok = chk(ok, c == -1, 'inadmissible flight (no flown miss) -> below the floor');
c = optimality_status(mk('stage', 1, 'flyKm', 150, 'reason', 'flown position miss 150.0 km > 100'));
ok = chk(ok, c == -1, 'misses the target -> below the floor');
[c, r] = optimality_status(mk('stage', 2, 'reason', 'Hamiltonian max|H| = 1.0e-03 > 1e-06'));
ok = chk(ok, c == 1 && contains(r, 'Hamiltonian'), 'pointwise PMP fails -> 1 neither, reason kept');
c = optimality_status(mk('stage', 3, 'reason', 'tfMin witness exceeded its 300 s cap'));
ok = chk(ok, c == 1, 'witness not confirmed -> 1 neither');
c = optimality_status(mk('stage', 4, 'reason', 'conjugate test verdict 0 (UNDETERMINED)'));
ok = chk(ok, c == 3, 'conjugate UNDETERMINED -> 3 necessary only');
c = optimality_status(mk('stage', 7, 'conjFound', true, 'hypAfterConj', 'held', 'reason', 'conjugate test verdict 0 (FAIL)'));
ok = chk(ok, c == 2, 'conjugate FAIL with gates + H6 held -> 2');
[c, r] = optimality_status(mk('stage', 4, 'conjFound', true, 'hypAfterConj', 'dim S = 2 (abnormal lift)', 'reason', 'conjugate test verdict 0 (FAIL)'));
ok = chk(ok, c == 3 && contains(r, 'hypotheses were not established') && contains(r, 'dim S'), ...
         'conjugate FAIL without its hypotheses -> 3, saying which failed');
c = optimality_status(mk('stage', 5, 'reason', 'H6 not established: margin 0.9'));
ok = chk(ok, c == 3, 'a sufficiency gate fails -> 3 necessary only');
c = optimality_status(mk('stage', 7, 'conjFound', true, 'hypAfterConj', 'held', ...
         'reason', 'dense conjugate scan not clear: 0 coarse sign change(s), 1 zero, 0 UNRESOLVED, 0 multiplicity (3 near-miss cleared)'));
ok = chk(ok, c == 2, 'dense scan finds a zero -> 2');
c = optimality_status(mk('stage', 7, 'reason', 'dense conjugate scan not clear: 0 coarse sign change(s), 0 zero, 1 UNRESOLVED, 0 multiplicity (3 near-miss cleared)'));
ok = chk(ok, c == 3, 'dense scan UNRESOLVED only -> 3');

% ---- legacy records (no .stage): inferred, conservative --------------------
L = rmfield(base, {'stage', 'conjFound', 'hypAfterConj'});
lg = @(varargin) setf(L, varargin{:});
[c, ~, inf_] = optimality_status(lg('ok', true, 'reason', 'certified'));
ok = chk(ok, c == 4 && inf_, 'legacy certified -> 4, inferred');
c = optimality_status(lg('reason', 'conjugate test verdict 0'));
ok = chk(ok, c == 3, 'legacy "verdict 0" is FAIL or UNDETERMINED -> 3 (re-certify to resolve)');
c = optimality_status(lg('reason', 'dense conjugate scan not clear: 0 coarse sign change(s), 1 zero, 0 UNRESOLVED, 0 multiplicity (3 near-miss cleared)'));
ok = chk(ok, c == 2, 'legacy dense-scan zero -> 2 (the dense scan runs after gates + H6)');
c = optimality_status(lg('reason', 'independent-field Hamiltonian residual 2.0e-05 > 1e-06'));
ok = chk(ok, c == 3, 'legacy gate failure that mentions "Hamiltonian" -> 3, not 1');
c = optimality_status(lg('reason', 'minimum-principle |full gap| 1.0e-09 > 1e-12'));
ok = chk(ok, c == 1, 'legacy pointwise failure -> 1');
[c, r, inf_] = optimality_status(lg('reason', 'some wording no build ever used'));
ok = chk(ok, c == 1 && inf_ && contains(r, 'inferred'), 'REVIEW FOCUS 1: unrecognised legacy reason -> lowest recordable tier, inferred');

% ---- the legend -------------------------------------------------------------
K = status_key();
ok = chk(ok, isequal(K.codes, [4 3 2 1 0 -1]) && numel(K.names) == 6 && strcmp(K.names{3}, 'conjugate point found'), ...
         'status_key names the six codes');

if ok, fprintf('test_optimality_status: ALL PASS\n'); else, fprintf('test_optimality_status: FAIL\n'); end
end

function s = setf(s, varargin)
% SETF  Set name/value pairs on a struct.  INPUTS: s; pairs.  OUTPUTS: s.
for k = 1:2:numel(varargin), s.(varargin{k}) = varargin{k+1}; end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
```

- [ ] **Step 2: Run it to verify it fails**

MCP `evaluate_matlab_code`: `cd('<IND>/tests'); test_optimality_status()`
Expected: error `Undefined function 'optimality_status'`.

- [ ] **Step 3: Implement `status_key.m`**

```matlab
function K = status_key()
%% Purpose:
%
%   The legend of the optimality STATUS carried by every library entry and
%   alternative: what each code means and the order of the checks that
%   define it. Stored in the catalog as .status_key so the file explains
%   itself to a recipient.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  K                        struct                  .codes [1x6] .names
%                                                   {1x6} .meaning {1x6}
%                                                   .order .note
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

K = struct();
K.codes = [4 3 2 1 0 -1];
K.names = {'sufficient', 'necessary only', 'conjugate point found', 'neither', 'empty', 'below floor'};
K.meaning = { ...
    'first- and second-order sufficient conditions hold (the full certify_root stack passed)', ...
    'first-order necessary conditions hold (pointwise Pontryagin checks + tfMin witness); sufficiency not established', ...
    'an extremal with a conjugate point on the arc, under the hypotheses (gates + H6) that make it a refutation: NOT locally optimal', ...
    'flies to the target, but the necessary conditions are not established', ...
    'no transfer recorded in this cell', ...
    'not a transfer of the library (no converged root, or it does not fly to the target); never stored'};
K.order = ['certify_root: 1 polish, 2 flight (100 km / 10 m/s), 3 pointwise Pontryagin, 4 tfMin witness, ' ...
           '5 coarse conjugate test, 6 hypothesis gates (min|lam_v|, Q_mt, dim S, X2, lift margin), 7 H6, 8 dense conjugate scan'];
K.note = ['A status records what was ESTABLISHED, not what might be true: a check that could not run gives the ' ...
          'highest tier actually established, with the reason. An orbit transfer does not have to be optimal to be useful.'];
end
```

- [ ] **Step 4: Implement `optimality_status.m`**

```matlab
function [code, reason, inferred] = optimality_status(C, opts)
%% Purpose:
%
%   THE status of one certify_root result: 4 sufficient, 3 necessary only,
%   2 conjugate point found, 1 neither, -1 below the floor (not a transfer
%   of the library). One function, called by certify_root (stamping) and by
%   the harvester (legacy records), so the two cannot disagree.
%
%  ASSUMPTIONS / NOTES:
%
% • A result stamped by the current certifier carries .stage (the last
%   stage passed: 1 polish, 2 flight, 3 pointwise, 4 witness, 5 coarse
%   conjugate, 6 gates, 7 H6, 8 dense scan) and .conjFound/.hypAfterConj.
% • A LEGACY result has none of these; its status is inferred from the
%   reason text, conservatively: anything unrecognised takes the lowest
%   recordable tier and inferred = true. A legacy "conjugate test verdict
%   0" is FAIL or UNDETERMINED and is 3 until re-certified.
% • The floor: a finite z8 and a flown miss within opts.gateKm / gateVms.
%
%% Inputs:
%
%  C                        struct                  certify_root result
%  opts                     struct (optional)       .gateKm [100] .gateVms [10]
%
%% Outputs:
%
%  code                     double                  4 | 3 | 2 | 1 | -1
%  reason                   char                    why, in the certifier's words
%  inferred                 logical                 true for a legacy record
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
gateKm = fieldd(opts, 'gateKm', 100);  gateVms = fieldd(opts, 'gateVms', 10);
r = '';  if isfield(C, 'reason') && ischar(C.reason), r = C.reason; end
inferred = ~isfield(C, 'stage');

% ---- the floor --------------------------------------------------------------
z = [];  if isfield(C, 'z'), z = C.z; end
if ~(isnumeric(z) && numel(z) == 8 && all(isfinite(z(:))))
    code = -1;  reason = sprintf('below the floor: no converged root (%s)', r);  return
end
km = fieldd(C, 'flyKm', NaN);  vms = fieldd(C, 'flyVms', NaN);
if ~(isscalar(km) && isscalar(vms) && km <= gateKm && vms <= gateVms)
    code = -1;  reason = sprintf('below the floor: does not fly to the target (%s)', r);  return
end
if isfield(C, 'ok') && isequal(C.ok, true)
    code = 4;  reason = 'full stack passed';  return
end

if ~inferred
    % ---- stamped by the current certifier -----------------------------------
    if C.stage < 4
        code = 1;  reason = ['neither: necessary conditions not established -- ' r];
    elseif isfield(C, 'conjFound') && isequal(C.conjFound, true)
        if strcmp(C.hypAfterConj, 'held')
            code = 2;  reason = ['conjugate point found (gates and H6 held): ' r];
        else
            code = 3;  reason = sprintf(['necessary only: conjugate point found, but its hypotheses were not ' ...
                                         'established (%s) -- %s'], C.hypAfterConj, r);
        end
    else
        code = 3;  reason = ['necessary only: ' r];
    end
    return
end

% ---- legacy: inferred from the reason text -------------------------------------
pointwise = '^(Hamiltonian max|transversality|adjoint equations|minimum-principle|applied throttle|pointwise checks|the tight-tolerance flight|tfMin witness|the witness flight|witness flight inadmissible|witness solution)';
gates = '^(X2:|gates |gate |hypothesis gates|min\|lam_v\||min Q_mt|dim S|accepted lift|independent-field|second gates|lift_margin|H6|H2|H3|DIAGNOSTIC)';
if ~isempty(regexp(r, pointwise, 'once'))
    code = 1;  reason = ['neither (inferred): ' r];
elseif startsWith(r, 'conjugate test verdict')
    code = 3;  reason = ['necessary only (inferred): legacy verdict 0 is FAIL or UNDETERMINED; re-certify to resolve -- ' r];
elseif startsWith(r, 'dense conjugate scan not clear')
    n = str2double(regexp(r, '(\d+) coarse sign change\(s\), (\d+) zero, (\d+) UNRESOLVED, (\d+) multiplicity', 'tokens', 'once'));
    if numel(n) == 4 && (n(1) + n(2) + n(4)) > 0
        code = 2;  reason = ['conjugate point found (inferred; the dense scan runs after gates and H6): ' r];
    else
        code = 3;  reason = ['necessary only (inferred): ' r];
    end
elseif ~isempty(regexp(r, gates, 'once')) || contains(r, 'lower-bound ESTIMATE')
    code = 3;  reason = ['necessary only (inferred): ' r];
else
    code = 1;  reason = ['neither (inferred: reason not recognised): ' r];
end
end

% ---------------------------------------------------------------------------
function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
```

- [ ] **Step 5: Run the test to verify it passes**

MCP: `cd('<IND>/tests'); clear functions; test_optimality_status()`
Expected: every line PASS, `test_optimality_status: ALL PASS`.

- [ ] **Step 6: Mutation check (house rule: every check must be able to fail)**

Copy `optimality_status.m` into three scratch folders with one mutation each, put a copy of the test beside each (as `<scratch>/mutN/tests/`, with `<scratch>/mutN/results` a symlink to `RES` if needed), `cd` into `mutN/tests`, run, expect FAIL:
1. `code = 2;  reason = ['conjugate point found (gates and H6 held): ' r];` → replace `strcmp(C.hypAfterConj, 'held')` with `true`.
2. In the legacy branch swap the order of the `pointwise` and `gates` tests.
3. Final `else` → `code = 4`.
Record which checks caught each.

- [ ] **Step 7: Commit (host)**

```bash
git add orbit_transfer/DRO_tulip/indirect/optimality_status.m orbit_transfer/DRO_tulip/indirect/status_key.m orbit_transfer/DRO_tulip/indirect/tests/test_optimality_status.m
git commit -m "optimality_status + status_key: the four statuses from a certify_root result (spec 4)"
```

---

### Task 2: `certify_root` stamps the status and runs the gates after a conjugate FAIL

**Files:**
- Modify: `IND/certify_root.m` (wrapper at the top, new fields at birth, stage marks, conjugate branch, H6 return, dense-scan flag, conj override seam)
- Test: `IND/tests/test_certify_status.m` (new); run existing `CC/tests/test_certify_schema.m`, `test_certify_enforcement.m`, `test_certify_caps.m`, `IND/tests/test_certify_crossing.m`

**Interfaces:**
- Consumes: `optimality_status` (Task 1).
- Produces: every `certify_root` result carries `.stage` (0-8), `.conjVerdict` (char: 'PASS' | 'FAIL' | 'UNDETERMINED' | ''), `.conjFound` (logical), `.conjReason` (char), `.hypAfterConj` (char: 'held' or the gate's reason), `.overridden` (logical), `.status` (double), `.status_reason` (char). Test seam: `opts.override.conj` struct overwrites fields of the coarse conjugate result (`.pass`, `.verdict`).

- [ ] **Step 1: Write the failing test** (`IND/tests/test_certify_status.m`)

```matlab
function ok = test_certify_status()
% TEST_CERTIFY_STATUS  certify_root stamps .status on every result; a coarse
% conjugate FAIL still runs the hypothesis gates and H6 (for the record) and
% is status 2 when they hold, 3 when they do not; UNDETERMINED is 3; a
% pointwise failure is 1; a dense-scan zero is 2; the untouched anchor is 4.
% Uses the test seam (.override), so every forced result is .overridden.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));                     % IND
cc = fullfile(fileparts(fileparts(here)), 'costate_common');
addpath(here, cc);
lib = dro_tulip_library();
E = lib(strcmp({lib.src}, 'anchor'));  E = E(1);
[B, ~] = arclength_arrival('setup');
rv0 = B.stateD(E.sD);  rvf = B.stateA(E.sA);
seed = seed_from_z8(E.z(:), rv0(1:6), 24, B.Tnd, B.cnd, B.mu);
base = struct('sA', E.sA, 'sD', E.sD, 'pool', capped_pool());
run = @(o) certify_root(seed, rv0(1:6), rvf(1:6), B, o);

C = run(base);
ok = chk(ok, C.ok && C.status == 4 && C.stage == 8 && strcmp(C.status_reason, 'full stack passed') && ~C.overridden, ...
         sprintf('anchor: status 4, stage 8 (%s)', C.reason));
ok = chk(ok, strcmp(C.conjVerdict, 'PASS'), 'the coarse verdict string is stored');

o = base;  o.override = struct('conj', struct('pass', false, 'verdict', 'FAIL'));
C = run(o);
ok = chk(ok, ~C.ok && C.conjFound && strcmp(C.hypAfterConj, 'held') && C.stage == 7 && C.status == 2 && C.overridden, ...
         sprintf('coarse FAIL, gates + H6 held -> 2 (stage %d, %s)', C.stage, C.status_reason));
ok = chk(ok, startsWith(C.reason, 'conjugate test verdict 0 (FAIL)') && ~isempty(C.g) && isfinite(C.h6Margin) && isempty(C.conjDense), ...
         'the reason is the conjugate one; gates and H6 were computed; the dense scan was not run');

o.override.g = struct('dimS', 2);
C = run(o);
ok = chk(ok, C.status == 3 && contains(C.hypAfterConj, 'dim S') && startsWith(C.reason, 'conjugate test verdict 0 (FAIL)'), ...
         sprintf('coarse FAIL, a gate fails -> 3 saying which (%s)', C.status_reason));

o = base;  o.override = struct('conj', struct('pass', false, 'verdict', 'UNDETERMINED'));
C = run(o);
ok = chk(ok, C.status == 3 && ~C.conjFound && C.stage == 4 && isempty(C.g), 'UNDETERMINED -> 3, gates not run');

o = base;  o.override = struct('PW', struct('Hmax', 1));
C = run(o);
ok = chk(ok, C.status == 1 && C.stage == 2, sprintf('pointwise failure -> 1 (%s)', C.reason));

o = base;  o.override = struct('CS', struct('nZero', 1, 'clear', false));
C = run(o);
ok = chk(ok, C.status == 2 && C.conjFound && C.stage == 7, sprintf('dense-scan zero -> 2 (%s)', C.reason));

if ok, fprintf('test_certify_status: ALL PASS\n'); else, fprintf('test_certify_status: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
```

- [ ] **Step 2: Run it to verify it fails** (batch: it solves)

```bash
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect/tests'); exit(~test_certify_status())"
```
Expected: FAIL at the first check (`Unrecognized field name "status"`).

- [ ] **Step 3: Add the public wrapper.** Rename line 1 `function C = certify_root(seed, rv0, rvf, B, opts)` to stay as is, and insert immediately after the line `%% ------------------------ Begin Code Sequence ---------------------------`:

```matlab
if nargin < 5, opts = struct(); end
C = certify_core(seed, rv0, rvf, B, opts);
% THE STATUS, stamped on every result (spec 4). A coarse conjugate FAIL ran
% the gates and H6 for the record: if they held the reason stays the
% conjugate one; if one failed, that gate's reason moves to .hypAfterConj and
% the conjugate reason comes back.
if C.conjFound
    if C.stage >= 7
        C.hypAfterConj = 'held';
    else
        C.hypAfterConj = C.reason;  C.reason = C.conjReason;
    end
end
[C.status, C.status_reason] = optimality_status(C);
end

function C = certify_core(seed, rv0, rvf, B, opts)
% CERTIFY_CORE  The gate stack itself (the body certify_root always had).
% INPUTS/OUTPUTS: as certify_root.
```
The existing body (starting with its own `if nargin < 5 ...` / option parsing) now belongs to `certify_core`.

- [ ] **Step 4: Fields at birth.** In the `C = struct('ok', false, ...)` initialisation, append before the closing `);`:

```matlab
           , 'stage', 0, 'conjVerdict', '', 'conjFound', false, 'conjReason', '', 'hypAfterConj', '', ...
           'overridden', ~isempty(fieldnames(ovr)), 'status', 0, 'status_reason', ''
```
(`ovr` is parsed above the initialisation; if it is parsed below, move the line `ovr = d('override', struct());` above `C = struct(...)`.)

- [ ] **Step 5: Stage marks.** Insert, each on its own line:
  - after `C.z = z(:);  C.Y = it.Y;  C.tfDays = z(8)*tStar/86400;` → `C.stage = 1;`
  - after the line `if ~(C.flyVms < gateVms), C.reason = ...` → `C.stage = 2;`
  - after the throttle check that closes the pointwise block (the `if` whose reason starts `'applied throttle is not 1'` and its `end`) → `C.stage = 3;`
  - after `if ~(isfinite(C.dz) && C.dz <= tolDz), ...` → `C.stage = 4;`
  - after the lift_margin block (the `if ~okL || ~certL ... end`), before `% ---- 6. H6, enforced` → `C.stage = 6;`
  - after the `if ~(h6m > h6MarginMin) ... end` block, before `% ---- 7. dense conjugate scan` → 
    ```matlab
    C.stage = 7;
    if C.conjFound, C.wallSec = toc(t0);  return, end   % gates + H6 were for the record only
    ```
  - after `C.fullStack = true;` → `C.stage = 8;`

- [ ] **Step 6: The conjugate branch.** Directly after the polish call succeeds (after `tick();` following `ms_tfmin`), add the seam:

```matlab
if isfield(it, 'conj') && isstruct(it.conj), it.conj = applyOverride(it.conj, ovr, 'conj'); end
```
Replace the whole `% ---- 4. conjugate test` block with:

```matlab
% ---- 4. conjugate test ----------------------------------------------------
if isfield(it, 'conj') && isfield(it.conj, 'pass')
    [okJ, jv] = scalar_verdict(it.conj.pass);
    if okJ, C.conj = jv; else, C.conj = -1; end
end
if isfield(it, 'conj') && isfield(it.conj, 'verdict'), C.conjVerdict = char(it.conj.verdict); end
if ~(C.conj == 1)
    C.reason = sprintf('conjugate test verdict %g (%s)', C.conj, C.conjVerdict);
    if ~strcmp(C.conjVerdict, 'FAIL'), C.wallSec = toc(t0);  return, end
    % A FAIL IS A FINDING, NOT A STOP: a conjugate point refutes optimality
    % only under the hypotheses the gates and H6 check, so they are run for
    % the record and the result returns after H6 (spec 2, Mike 2026-10-04).
    C.conjFound = true;  C.conjReason = C.reason;
else
    C.stage = 5;
end
```

- [ ] **Step 7: Dense-scan flag.** In the dense block, inside `if ~CS.clear`, before its `C.reason = sprintf(['dense conjugate scan not clear: ...`, add:

```matlab
        C.conjFound = (CS.nZero + CS.nInterior + CS.multiplicity) > 0;
        if C.conjFound, C.conjReason = ''; end
```
and after that `C.reason = sprintf(...)` line add `if C.conjFound, C.conjReason = C.reason; end`.

- [ ] **Step 8: Run the new test and the existing certifier suites** (batch, ~20-40 min)

```bash
M=/Applications/MATLAB_R2026a.app/bin/matlab; OT=/Users/msc/Desktop/optimal_control/orbit_transfer
$M -batch "cd('$OT/DRO_tulip/indirect/tests'); exit(~test_certify_status())"
for t in test_certify_schema test_certify_enforcement test_certify_caps; do $M -batch "cd('$OT/costate_common/tests'); exit(~$t())"; echo "$t exit $?"; done
$M -batch "cd('$OT/DRO_tulip/indirect/tests'); exit(~test_certify_crossing())"
```
Expected: all ALL PASS / exit 0. If `test_certify_schema` enumerates the result's fields, add the eight new names (`stage conjVerdict conjFound conjReason hypAfterConj overridden status status_reason`) to its expected list — that is the only allowed edit to existing tests.

- [ ] **Step 9: Mutation check.** Revert one edit at a time in a scratch copy (folder-precedence pattern of Task 1): (a) delete the `C.conjFound = true; ...` line; (b) remove the `if C.conjFound ... return` after H6; (c) drop the wrapper's `C.reason = C.conjReason` swap. Each must make `test_certify_status` FAIL.

- [ ] **Step 10: Commit (host)** — `certify_root: stamp the status; gates and H6 for the record after a conjugate FAIL; verdict string stored (spec 4, 5)`.

---

### Task 3: Shared units — `same_root`, `catalog_setup_request`, `make_alternative`, `dedup_alternatives`

**Files:**
- Create: `IND/same_root.m`, `IND/catalog_setup_request.m`, `IND/make_alternative.m`, `IND/dedup_alternatives.m`
- Modify: `IND/run_phase_torus.m` (local `sameRoot` delegates), `IND/fill_holes_direct.m` (local `physicsFromCatalog` delegates)
- Test: `IND/tests/test_alternative_rows.m`

**Interfaces:**
- Produces: `tf = same_root(zA, zB)` — the registry rule (Global Constraints).
- Produces: `so = catalog_setup_request(cat_, sD1)` — the `arclength_arrival('setup', so)` request from a catalog's own keys (closures only).
- Produces: `A = make_alternative(C, sD, sA, iD, iA, source)` — one alternatives row or `[]` when below the floor / overridden. Row fields, in this order: `sD sA iD iA z8 junctions tf_nd status status_reason inferred conj conjVerdict minLamV minQmt dimS h6Margin liftMargin flyKm flyVms source`. Status taken from `C.status`/`C.status_reason` when present, else `optimality_status(C)` (inferred).
- Produces: `A = dedup_alternatives(A, sheet)` — drops rows that are the same root (`same_root`) as the primary of their cell (`iD,iA > 0`) or as an earlier row at the same phases (1e-9, circular); returns rows sorted by (sD, sA, tf_nd).

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_alternative_rows()
% TEST_ALTERNATIVE_ROWS  make_alternative builds one row from a certify_root
% result (stamped or legacy), refuses what is below the floor or overridden;
% dedup_alternatives removes a primary's twin and repeated roots; same_root
% is symmetric; the delegates agree with the shared units.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
z = [(1:7).'; 4.1];  Y = rand(14, 24);
C = struct('ok', false, 'z', z, 'Y', Y, 'flyKm', 0.2, 'flyVms', 0.01, 'reason', 'H6 not established: x', ...
           'stage', 6, 'conjFound', false, 'hypAfterConj', '', 'overridden', false, 'status', 3, ...
           'status_reason', 'necessary only: H6 not established: x', 'conj', 1, 'conjVerdict', 'PASS', ...
           'g', struct('minLamV', 0.1, 'minQmt', 0.2, 'dimS', 1), 'h6Margin', 0.9, 'liftMargin', 30);
A = make_alternative(C, 0.25, 0.5, 7, 13, 'test');
ok = chk(ok, isstruct(A) && A.status == 3 && isequal(A.z8, z) && isequal(A.junctions, Y) && A.tf_nd == 4.1 && ~A.inferred ...
         && A.iD == 7 && A.dimS == 1 && strcmp(A.source, 'test'), 'a stamped result becomes one row with its numbers');
L = rmfield(C, {'stage', 'conjFound', 'hypAfterConj', 'status', 'status_reason', 'overridden', 'conjVerdict'});
L.reason = 'conjugate test verdict 0';
A2 = make_alternative(L, 0.25, 0.5, 7, 13, 'legacy');
ok = chk(ok, A2.status == 3 && A2.inferred, 'a legacy result is classified on the way in (inferred)');
Cb = C;  Cb.flyKm = 500;
ok = chk(ok, isempty(make_alternative(Cb, 0.25, 0.5, 7, 13, 'x')), 'below the floor -> no row');
Co = C;  Co.overridden = true;
ok = chk(ok, isempty(make_alternative(Co, 0.25, 0.5, 7, 13, 'x')), 'an overridden (test-seam) result -> no row');
Cc = C;  Cc.ok = true;  Cc.status = 4;
A4 = make_alternative(Cc, 0.25, 0.5, 7, 13, 'x');
ok = chk(ok, A4.status == 4, 'a certified slower root is a status-4 row');

% dedup: a primary's twin, a repeated root, a distinct root kept
sheet = struct('sD_frac', [0 0.25], 'sA_frac', [0.5 0.75], 'has_solution', [false false; true false], ...
               'entry_index', [0 0; 1 0], 'z8', z);
Ad = [make_alternative(C, 0.25, 0.5, 2, 1, 'twin of primary'), ...
      make_alternative(setz(C, z + [1; zeros(7, 1)]), 0.25, 0.5, 2, 1, 'distinct'), ...
      make_alternative(setz(C, z + [1; zeros(7, 1)]), 0.25, 0.5, 2, 1, 'repeat')];
Ad = dedup_alternatives(Ad, sheet);
ok = chk(ok, numel(Ad) == 1 && strcmp(Ad.source, 'distinct'), 'dedup drops the primary''s twin and the repeat');

% same_root and the delegates
ok = chk(ok, same_root(z, z + 1e-9) && same_root(z + 1e-9, z) && ~same_root(z, z + [0.01; zeros(7, 1)]), 'same_root: symmetric, 1e-6');
H = run_phase_torus('localfunctions');
ok = chk(ok, H.sameRoot(z, z + 1e-9) == same_root(z, z + 1e-9), 'run_phase_torus''s sameRoot delegates');
Lc = load(fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat'));  c = Lc.(char(fieldnames(Lc)));
Hf = fill_holes_direct('localfunctions');
ok = chk(ok, isequal(Hf.physicsFromCatalog(c, c.sheets(1), 0), catalog_setup_request(c, 0)), 'fill_holes_direct''s physicsFromCatalog delegates');

if ok, fprintf('test_alternative_rows: ALL PASS\n'); else, fprintf('test_alternative_rows: FAIL\n'); end
end

function C = setz(C, z)
% SETZ  Replace a result's z.  INPUTS: C; z.  OUTPUTS: C.
C.z = z;
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
```

- [ ] **Step 2: Run to verify it fails** (MCP) — expected `Undefined function 'make_alternative'`.

- [ ] **Step 3: `same_root.m`** — move the body of `run_phase_torus>sameRoot` (lines `function tf = sameRoot(zA, zB)` … `end`, including its comment) into `IND/same_root.m` with the house header, renamed `same_root`; replace the local function's body in `run_phase_torus.m` with:

```matlab
function tf = sameRoot(zA, zB)
% SAMEROOT  Delegate: the rule lives in same_root (one copy).
% INPUTS: zA, zB [8 x 1].  OUTPUTS: tf.
tf = same_root(zA, zB);
end
```

- [ ] **Step 4: `catalog_setup_request.m`** — move `fill_holes_direct>physicsFromCatalog`'s body into `IND/catalog_setup_request.m` as `so = catalog_setup_request(cat_, sD1)` (reading `sh = cat_.sheets(1)` itself); replace the local body with `so = catalog_setup_request(cat_, sD1);`.

- [ ] **Step 5: `make_alternative.m`**

```matlab
function A = make_alternative(C, sD, sA, iD, iA, source)
%% Purpose:
%
%   One row of a catalog's ALTERNATIVES table from a certify_root result: the
%   transfer (z8, junctions, t_f), its optimality status and the numbers
%   behind it, where it sits and where it came from. [] when the result is
%   not a transfer of the library (below the floor) or came through the
%   certifier's test seam.
%
%% Inputs:
%
%  C                        struct                  certify_root result (stamped or legacy)
%  sD, sA                   double                  exact phases
%  iD, iA                   double                  grid cell, 0 = off the grid
%  source                   char                    provenance
%
%% Outputs:
%
%  A                        struct or []            the row (fields: sD sA iD
%                                                   iA z8 junctions tf_nd
%                                                   status status_reason
%                                                   inferred conj conjVerdict
%                                                   minLamV minQmt dimS
%                                                   h6Margin liftMargin flyKm
%                                                   flyVms source)
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

A = [];
if isfield(C, 'overridden') && isequal(C.overridden, true), return, end
if isfield(C, 'status') && isfield(C, 'stage')
    code = C.status;  why = C.status_reason;  inferred = false;
else
    [code, why, inferred] = optimality_status(C);
end
if code < 1, return, end
Y = [];  if isfield(C, 'Y'), Y = C.Y; end
g = struct('minLamV', NaN, 'minQmt', NaN, 'dimS', NaN);
if isfield(C, 'g') && isstruct(C.g) && ~isempty(C.g)
    for f = fieldnames(g)', if isfield(C.g, f{1}), g.(f{1}) = C.g.(f{1}); end, end
end
A = struct('sD', mod(sD, 1), 'sA', mod(sA, 1), 'iD', iD, 'iA', iA, 'z8', C.z(:), 'junctions', Y, ...
           'tf_nd', C.z(8), 'status', code, 'status_reason', why, 'inferred', inferred, ...
           'conj', fieldd(C, 'conj', NaN), 'conjVerdict', fieldd(C, 'conjVerdict', ''), ...
           'minLamV', g.minLamV, 'minQmt', g.minQmt, 'dimS', g.dimS, ...
           'h6Margin', fieldd(C, 'h6Margin', NaN), 'liftMargin', fieldd(C, 'liftMargin', NaN), ...
           'flyKm', C.flyKm, 'flyVms', C.flyVms, 'source', source);
end

% ---------------------------------------------------------------------------
function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
```

- [ ] **Step 6: `dedup_alternatives.m`**

```matlab
function A = dedup_alternatives(A, sheet)
%% Purpose:
%
%   One root, one record: drop every alternative that is the same root
%   (same_root) as the PRIMARY of its cell, or as an earlier row at the same
%   phases; return the rest sorted by (sD, sA, t_f).
%
%% Inputs:
%
%  A                        struct array            alternatives rows (make_alternative)
%  sheet                    struct                  catalog sheet (.has_solution
%                                                   .entry_index .z8)
%
%% Outputs:
%
%  A                        struct array            deduplicated, sorted
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if isempty(A), return, end
keep = true(1, numel(A));
near = @(a, b) abs(mod(a - b + 0.5, 1) - 0.5) < 1e-9;
for k = 1:numel(A)
    a = A(k);
    if a.iD > 0 && a.iA > 0 && sheet.has_solution(a.iD, a.iA, 1)
        if same_root(a.z8, sheet.z8(:, sheet.entry_index(a.iD, a.iA, 1))), keep(k) = false;  continue, end
    end
    for m = 1:k-1
        if keep(m) && near(A(m).sD, a.sD) && near(A(m).sA, a.sA) && same_root(A(m).z8, a.z8)
            keep(k) = false;  break
        end
    end
end
A = A(keep);
[~, order] = sortrows([[A.sD].', [A.sA].', [A.tf_nd].']);
A = A(order);
end
```

- [ ] **Step 7: Run tests** — `test_alternative_rows` (MCP), then batch `test_fill_holes_physics`, `test_fill_holes_pool`, `test_run_phase_torus_p0` (the delegates). Expected ALL PASS.

- [ ] **Step 8: Mutation check** — (a) `make_alternative` without the overridden return; (b) `dedup_alternatives` without the primary test; (c) `same_root` tolerance 1e-3. Each must fail the test.

- [ ] **Step 9: Commit (host)** — `same_root, catalog_setup_request promoted (delegates left); make_alternative + dedup_alternatives (spec 3.2)`.

---

### Task 4: `entry_thrust_program` — fly an entry from its junctions

**Files:**
- Create: `IND/entry_thrust_program.m`
- Test: `IND/tests/test_entry_thrust_program.m`

**Interfaces:**
- Produces: `P = entry_thrust_program(z8, junctions, rv0, rvf, phys, opts)`; `phys` = `.Tnd .cnd .mu .lStar .tStar .m0kg`; `opts.nSample` [50 per segment]. Returns `P.t [n x 1]` (ND), `P.tDays`, `P.Y [n x 14]`, `P.u [n x 3]` thrust direction = `-lam_v/|lam_v|`, `P.massKg [n x 1]`, `P.defectKm [1 x K-1]` junction position defects, `P.defectVms [1 x K-1]`, `P.startErr` (junction 1 vs `[rv0; 1; z8(1:7)]`, relative), `P.flyKm`, `P.flyVms` (arrival miss), `P.K`.

- [ ] **Step 1: Write the failing test** — oracle: arrival state from the ORBIT at the entry's phase (`B.stateA`), and the thrust direction from the equations of motion (applied acceleration minus the ballistic one), not from the function's own formula.

```matlab
function ok = test_entry_thrust_program()
% TEST_ENTRY_THRUST_PROGRAM  Flying an entry segment by segment from its
% junctions reaches the arrival state the ORBIT gives at its phase (oracle:
% B.stateA, not the catalog), with small junction defects; the reported
% thrust direction is the one the equations of motion apply (oracle: EoM
% with thrust minus EoM without); K is read from the array (K = 12 works).
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
lib = dro_tulip_library();
E = lib(strcmp({lib.src}, 'anchor'));  E = E(1);
[B, ~] = arclength_arrival('setup');
rv0 = B.stateD(E.sD);  rv0 = rv0(1:6);  rvf = B.stateA(E.sA);  rvf = rvf(1:6);
phys = struct('Tnd', B.Tnd, 'cnd', B.cnd, 'mu', B.mu, 'lStar', 389703.264829278, 'tStar', 382981.289129055, 'm0kg', 150);
for K = [24 12]
    seed = seed_from_z8(E.z(:), rv0, K, B.Tnd, B.cnd, B.mu);
    P = entry_thrust_program(E.z(:), seed.Y(:, 1:K), rv0, rvf, phys);
    ok = chk(ok, P.K == K && P.flyKm < 1 && P.flyVms < 0.1, sprintf('K = %d: arrives at the orbit''s state (%.3f km, %.4f m/s)', K, P.flyKm, P.flyVms));
    ok = chk(ok, max(P.defectKm) < 1e-3 && P.startErr < 1e-12, sprintf('K = %d: junction defects %.1e km, start %.1e', K, max(P.defectKm), P.startErr));
end
k = round(numel(P.t)/2);
dyT = pumpkyn.cr3bp.tfMinEoM(P.t(k), P.Y(k, :).', B.Tnd, B.cnd, B.mu);
dy0 = pumpkyn.cr3bp.tfMinEoM(P.t(k), P.Y(k, :).', 0, B.cnd, B.mu);
a = dyT(4:6) - dy0(4:6);  a = a/sqrt(sum(a.^2));
ok = chk(ok, sqrt(sum((a(:) - P.u(k, :).').^2)) < 1e-9, 'thrust direction = the EoM''s applied acceleration direction');
ok = chk(ok, abs(P.massKg(end) - 150*(1 - B.Tnd*E.z(8)/B.cnd)) < 1e-6, 'final mass = the all-burn mass law');
if ok, fprintf('test_entry_thrust_program: ALL PASS\n'); else, fprintf('test_entry_thrust_program: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
```

- [ ] **Step 2: Run to verify it fails** (batch: it flies) — expected `Undefined function 'entry_thrust_program'`.

- [ ] **Step 3: Implement**

```matlab
function P = entry_thrust_program(z8, junctions, rv0, rvf, phys, opts)
%% Purpose:
%
%   A library entry's TRANSFER AND THRUST PROGRAM, flown segment by segment
%   from its multiple-shooting junction states -- accurate where a single
%   flight from z8 amplifies the costates' error over 17-26 days (FINDINGS
%   97). Reports the junction defects and the arrival miss so a recipient
%   sees how well the stored junctions represent the transfer.
%
%  ASSUMPTIONS / NOTES:
%
% • Junctions sit on the uniform grid linspace(0, t_f, K+1)
%   (flight_to_junctions); column k is the augmented state
%   [r; v; m; lam_r; lam_v; lam_m] at the start of segment k.
% • Minimum-time, all-burn: thrust direction u = -lam_v/|lam_v| (checked in
%   the test against the equations of motion), throttle 1.
%
%% Inputs:
%
%  z8                       [8 x 1]                 [lambda0(7); t_f] (ND)
%  junctions                [14 x K]                junction start states
%  rv0, rvf                 [6 x 1]                 departure / arrival states (ND)
%  phys                     struct                  .Tnd .cnd .mu .lStar .tStar .m0kg
%  opts                     struct (optional)       .nSample [50] per segment
%
%% Outputs:
%
%  P                        struct                  .t .tDays .Y .u .massKg
%                                                   .defectKm .defectVms
%                                                   .startErr .flyKm .flyVms .K
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 6, opts = struct(); end
nS = 50;  if isfield(opts, 'nSample') && ~isempty(opts.nSample), nS = opts.nSample; end
K = size(junctions, 2);
assert(size(junctions, 1) == 14 && K >= 1 && all(isfinite(junctions(:))), 'entry_thrust_program:junctions', ...
       'junctions must be a finite 14 x K array');
tf = z8(8);  h = tf/K;
y1 = [rv0(:); 1; z8(1:7)];
P = struct('K', K, 'startErr', sqrt(sum((junctions(:, 1) - y1).^2))/sqrt(sum(y1.^2)));
t = [];  Y = [];  dK = zeros(1, max(K - 1, 0));  dV = dK;
for k = 1:K
    [tk, Yk] = pumpkyn.cr3bp.tfMinProp(h, junctions(:, k), phys.Tnd, phys.cnd, phys.mu);
    tq = linspace(0, h, nS).';
    Yq = interp1(tk, Yk, tq, 'pchip');
    Yq(end, :) = Yk(end, :);
    if k < K
        dK(k) = sqrt(sum((Yk(end, 1:3) - junctions(1:3, k+1).').^2))*phys.lStar;
        dV(k) = sqrt(sum((Yk(end, 4:6) - junctions(4:6, k+1).').^2))*phys.lStar/phys.tStar*1000;
    end
    if k > 1, tq = tq(2:end);  Yq = Yq(2:end, :); end
    t = [t; (k - 1)*h + tq];  Y = [Y; Yq];
end
lv = Y(:, 11:13);
P.t = t;  P.tDays = t*phys.tStar/86400;  P.Y = Y;
P.u = -lv ./ sqrt(sum(lv.^2, 2));
P.massKg = Y(:, 7)*phys.m0kg;
P.defectKm = dK;  P.defectVms = dV;
P.flyKm = sqrt(sum((Y(end, 1:3) - rvf(1:3).').^2))*phys.lStar;
P.flyVms = sqrt(sum((Y(end, 4:6) - rvf(4:6).').^2))*phys.lStar/phys.tStar*1000;
end
```
(`t = [t; ...]` grows by K pieces only; K <= 48, so no preallocation is needed — say so if Code Analyzer flags it, do not add a pragma.)

- [ ] **Step 4: Run to verify it passes** (batch). Expected ALL PASS. If the thrust-direction check fails with the opposite sign, the convention is `+lam_v` — fix `P.u` and the header, do NOT weaken the test.

- [ ] **Step 5: Mutation check** — (a) `h = tf/(K+1)`; (b) `P.u = lv ./ ...`; (c) defects measured against `junctions(1:3, k)`. Each must fail.

- [ ] **Step 6: Commit (host)** — `entry_thrust_program: an entry's transfer and thrust program from its junctions (spec 8)`.

---

### Task 5: `rib_from_crossing` keeps its refusals; `walkPastConjugate` (opt-in)

**Files:**
- Modify: `IND/rib_from_crossing.m`
- Test: `IND/tests/test_rib_refused.m` (new); rerun `IND/tests/test_rib_from_crossing.m`

**Interfaces:**
- Consumes: `C.status` (Task 2).
- Produces: `R.refused` — struct array of certify_root results (every non-ok candidate with `status >= 1`, overridden ones included and flagged; the packager filters), each with added `.walked` (logical: the walk advanced through it). New option `opts.walkPastConjugate` [false]: when true a status-2 candidate advances the walk like a certified one (it goes to `R.refused` with `.walked = true`, never to `R.pts`). The checkpoint carries `refused`.

- [ ] **Step 1: Write the failing test** (uses the test seam: a forced coarse FAIL makes every trial status 2)

```matlab
function ok = test_rib_refused()
% TEST_RIB_REFUSED  A rib keeps every refused candidate above the floor,
% including the point it stopped at; with walkPastConjugate a status-2
% point advances the walk and is recorded, never as a certified point; off
% (the default) behaves as before.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
ot = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(fullfile(ot, 'costate_common'), fullfile(ot, 'DRO_tulip', 'indirect'));
[B, anc] = arclength_arrival('setup');
C0 = certify_crossing(anc.p, anc.sA, B, anc);
assert(C0.ok, 'fixture: the anchor must certify (%s)', C0.reason);
forced = struct('pool', capped_pool(), 'override', struct('conj', struct('pass', false, 'verdict', 'FAIL')));

R = rib_from_crossing(C0, B, anc, struct('nD', 12, 'direction', -1, 'nPts', 1, 'maxBisect', 1, 'copts', forced));
ok = chk(ok, isempty(R.pts) && startsWith(R.stop, 'stalled') && numel(R.refused) == 2, ...
         sprintf('default: stops, keeps both refused trials (%d)', numel(R.refused)));
ok = chk(ok, all([R.refused.status] == 2) && ~any([R.refused.walked]) && all([R.refused.overridden]), ...
         'refusals carry status 2, not walked, flagged overridden');

R = rib_from_crossing(C0, B, anc, struct('nD', 12, 'direction', -1, 'nPts', 1, 'maxBisect', 1, 'copts', forced, ...
                                         'walkPastConjugate', true));
ok = chk(ok, isempty(R.pts) && strcmp(R.stop, 'complete') && numel(R.refused) == 1 && R.refused(1).walked ...
         && abs(R.refused(1).sD - 11/12) < 1e-9, 'walkPastConjugate: advances through the status-2 point and records it');
if ok, fprintf('test_rib_refused: ALL PASS\n'); else, fprintf('test_rib_refused: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
```

- [ ] **Step 2: Run to verify it fails** (batch) — expected `Unrecognized field name "refused"`.

- [ ] **Step 3: Implement.** In `rib_from_crossing.m`:
  - header: document `.walkPastConjugate [false]` and the `R.refused` output.
  - after `progress = d('progress', []); ...` add `walkPast = d('walkPastConjugate', false);`
  - `R = struct('pts', struct([]), 'refused', struct([]), 'stop', '', ...` (add `'refused', struct([])`).
  - on resume, after `R.pts = Ck0.pts;` add `if isfield(Ck0, 'refused'), R.refused = Ck0.refused; end`.
  - replace `if Ct.ok` … through the end of its `else` branch with:

```matlab
        advance = Ct.ok || (walkPast && Ct.status == 2);
        if ~Ct.ok && Ct.status >= 1               % KEEP what the walk refused (spec 5)
            Ct.walked = advance;
            R.refused = append_point(R.refused, Ct);
        end
        if advance
            cur = trial;  zc = Ct.z;  Yc = Ct.Y;  Ck = Ct;  nGood = nGood + 1;
            lg('  rib sD %.4f: t_f %.4f d, %s', mod(cur,1), Ct.tfDays, Ct.reason);
            if nGood >= 2 && abs(step) < 1/nD
                step = 2*step;  nGood = 0;  nb = max(nb - 1, 0);
            end
        else
            if nb >= maxBisect
                R.stop = sprintf('stalled at sD = %.4f stepping to %.4f: %s', mod(cur,1), mod(trial,1), Ct.reason);
                lg('  rib STOP: %s', R.stop);
                return
            end
            step = step/2;  nb = nb + 1;  nGood = 0;  nHalved = nHalved + 1;
            lg('  rib sD %.4f -> %.4f refused (%s); halving to %.5f', mod(cur,1), mod(trial,1), Ct.reason, step);
        end
```
  - after the inner `while`: only certified lattice points join `R.pts`:

```matlab
    sD = mod(target, 1);  z = zc;  Y = Yc;
    Ck.sD = sD;
    if Ck.ok
        Ck.note = join_note(...);              % the existing note line, unchanged
        R.pts = append_point(R.pts, Ck);
    end
```
  - the checkpoint save struct gains `'refused', R.refused`.

- [ ] **Step 4: Run** `test_rib_refused` and `test_rib_from_crossing` (batch). Expected ALL PASS (the existing test still sees two certified points).

- [ ] **Step 5: Mutation check** — (a) keep only `Ct.status == 2` refusals; (b) let a walked status-2 point into `R.pts`; (c) default `walkPast = true`. Each must fail a test.

- [ ] **Step 6: Commit (host)** — `rib_from_crossing: keep refused candidates (stop point included); walkPastConjugate opt-in, off by default (spec 5)`.

---

### Task 6: `fill_holes_direct` keeps its candidates

**Files:**
- Modify: `IND/fill_holes_direct.m`
- Test: `IND/tests/test_fill_holes_pool.m` (extend with the helper and the wiring)

**Interfaces:**
- Produces: each cell record `rec(k)` gains `.others` — a struct array of certify_root results from that cell's seeds that are not the kept root: refused ones with `status >= 1`, and certified roots that lost (slower than the kept one, or not faster than the cell's existing entry). Local helper `others = keepCandidate(others, C)` appends `C` when `C.status >= 1` (via `append_point`, so differing field sets are harmonised).

- [ ] **Step 1: Extend the test** — append to `test_fill_holes_pool.m` before its final print:

```matlab
% ---- the filler keeps its candidates (optimality-status spec 5) -----------
C1 = struct('ok', false, 'status', 3, 'z', ones(8, 1));  C0 = struct('ok', false, 'status', -1, 'z', nan(8, 1));
o = H.keepCandidate(struct([]), C1);  o = H.keepCandidate(o, C0);
ok = chk(ok, numel(o) == 1 && o(1).status == 3, 'keepCandidate keeps a status >= 1 result, drops one below the floor');
ok = chk(ok, contains(loop{1}, 'others = keepCandidate(others, C)') && contains(src, '''others'''), ...
         'wiring: the seed loop keeps candidates and the cell record stores them');
```

- [ ] **Step 2: Run to verify it fails** (batch) — `Unrecognized field name "keepCandidate"`.

- [ ] **Step 3: Implement** in `fill_holes_direct.m`:
  - before the seed loop: `others = struct([]);`
  - in the seed loop, after `[C, dinfo] = direct_cell_solve(...)` and `tfD = ...`: `others = keepCandidate(others, C);` — then, where the race replaces `best` with a faster `C`, push the displaced `best` with `others = keepCandidate(others, best);` before reassigning; and in the `'not faster than the cell''s'` branch push `best` before `best = []`.
  - remove from `others` the kept `best` (if it was pushed) after the loop: `if ~isempty(best) && ~isempty(others), others = others(~arrayfun(@(o) same_root(o.z, best.z), others)); end`
  - `cellRec` gains a field: signature `cellRec(iD, iA, sDv, sAv, seed, tfD, tfC, ok, reason, wall, others)` and the struct gets `'others', {others}`; update both call sites (the no-neighbour call passes `struct([])`), and the empty `rec` initialiser gains `'others', {}`.
  - local helper:

```matlab
function others = keepCandidate(others, C)
% KEEPCANDIDATE  Keep a candidate the catalog may record as an alternative:
% any certify_root result with status >= 1 (the packager filters the test
% seam and the floor again).  INPUTS: others (struct array); C.  OUTPUTS: others.
if isstruct(C) && isfield(C, 'status') && C.status >= 1, others = append_point(others, C); end
end
```

- [ ] **Step 4: Run** `test_fill_holes_pool` and `test_fill_holes_physics` (batch). Expected ALL PASS.

- [ ] **Step 5: Commit (host)** — `fill_holes_direct: keep every candidate above the floor per cell (spec 5)`.

---

### Task 7: The packager writes status, junctions and alternatives

**Files:**
- Modify: `IND/sheet_to_catalog_file.m`, `CC/build_costate_catalog_family.m`, `IND/package_phase_catalog.m` (pass `extraAlternatives` through)
- Test: `IND/tests/test_package_status.m` (new); rerun `IND/tests/test_reproduce_library.m`, `CC/tests/test_build_costate_catalog_family*.m` if present (`ls CC/tests | grep -i catalog`)

**Interfaces:**
- Consumes: `make_alternative`, `dedup_alternatives` (Task 3); `R.refused` (Task 5).
- Produces (sheet file `Q`): `Q.STATUS int8 [nD x nA x 1]` (4 where OK), `Q.SREASON {nD x nA}`, `Q.JUNC {nD x nA}` (each `[14 x K]`), `Q.ALT` (alternatives struct array).
- Produces (catalog): `sheets(k).status` (= `Q.STATUS`), `sheets(k).status_reason {1 x nEntries}`, `sheets(k).junctions {1 x nEntries}` (z8 column order, like `entry_notes`), `cat_.alternatives` (all sheets' `Q.ALT`, each row gaining `.sheet = k`), `cat_.status_key = status_key()`.
- `sheet_to_catalog_file(S, ribs, outMat, opts)`: new `opts.extraAlternatives` (struct array of `make_alternative` rows, e.g. the filler's).

- [ ] **Step 1: Write the failing test** — a synthetic sheet using the REAL problem identity of the 96-phase sheet, two columns, fake z/Y.

```matlab
function ok = test_package_status()
% TEST_PACKAGE_STATUS  sheet_to_catalog_file + build_costate_catalog_family
% carry status, status_reason and junctions for every primary, and an
% alternatives table holding (a) the sheet's non-winning candidates above
% the floor, (b) a rib's refusals, (c) a certified rib point that lost its
% cell to a faster one -- deduplicated, overridden ones excluded.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
S96 = load(fullfile(here, 'results', 'sheet96_resolution_test', 'arrival_sheet_70mN_nA96.mat'));  P = S96.S.problem;
cert = @(z1, tf, sD, sA) struct('ok', true, 'status', 4, 'status_reason', 'full stack passed', 'stage', 8, ...
        'z', [z1; (2:7).'; tf], 'Y', z1*ones(14, 24), 'tfDays', tf*P.tStar/86400, 'sD', sD, 'sA', sA, ...
        'flyKm', 0.1, 'flyVms', 0.01, 'reason', 'certified', 'note', '', 'conj', 1, 'conjVerdict', 'PASS', ...
        'g', struct('minLamV', 1, 'minQmt', 1, 'dimS', 1), 'h6Margin', 5, 'liftMargin', 30, 'overridden', false);
ref = @(C, st) setf(C, 'ok', false, 'status', st, 'status_reason', sprintf('status %d', st), 'reason', 'x');
S = struct('problem', P, 'sA', [0.25 0.75], 'TF', [NaN NaN]);
c1 = [cert(1, 4.0, 0, 0.25), ref(cert(2, 4.5, 0, 0.25), 2), ref(setf(cert(3, 4.6, 0, 0.25), 'flyKm', 900), 3)];
c2 = [cert(4, 4.2, 0, 0.75), setf(ref(cert(5, 4.3, 0, 0.75), 3), 'overridden', true)];
S.cand = {c1, c2};  S.TF = [c1(1).tfDays, c2(1).tfDays];
rib = struct('sA', 0.25, 'pts', [cert(6, 4.1, 0.5, 0.25)], 'refused', [ref(cert(7, 4.4, 0.5, 0.25), 1)], 'stop', 'complete');
rib2 = struct('sA', 0.25, 'pts', [cert(8, 4.05, 0.5, 0.25)], 'refused', struct([]), 'stop', 'complete');
tmp = tempname;  mkdir(tmp);
save(fullfile(tmp, 'sheet.mat'), 'S');  R = rib;  save(fullfile(tmp, 'rib1.mat'), 'R');  R = rib2;  save(fullfile(tmp, 'rib2.mat'), 'R');
c = package_phase_catalog(fullfile(tmp, 'sheet.mat'), {fullfile(tmp, 'rib1.mat'), fullfile(tmp, 'rib2.mat')}, ...
        struct('nD', 2, 'outDir', tmp, 'name', 'cat_test', 'thrustN', P.thrustN, 'ispS', P.ispS, 'm0kg', P.m0kg));
sh = c.sheets(1);
ok = chk(ok, isequal(sh.status, int8([4 4; 4 0])) && numel(sh.junctions) == 3 && numel(sh.status_reason) == 3, ...
         'primaries: status 4 grid, junctions and reasons per entry');
k = sh.entry_index(2, 1);
ok = chk(ok, sh.z8(1, k) == 8 && isequal(sh.junctions{k}, 8*ones(14, 24)), 'the faster rib point is the primary and carries its junctions');
A = c.alternatives;  firsts = arrayfun(@(a) a.z8(1), A);
ok = chk(ok, isequal(sort(firsts), [2 6 7]), sprintf('alternatives: sheet refusal (2), rib refusal (7), the slower certified rib point (6) -- got %s', mat2str(sort(firsts))));
ok = chk(ok, ~any(firsts == 3) && ~any(firsts == 5), 'below the floor (3) and overridden (5) are not recorded');
ok = chk(ok, isfield(c, 'status_key') && isequal(c.status_key.codes, [4 3 2 1 0 -1]), 'the legend ships in the catalog');
rmdir(tmp, 's');
if ok, fprintf('test_package_status: ALL PASS\n'); else, fprintf('test_package_status: FAIL\n'); end
end

function s = setf(s, varargin)
% SETF  Set name/value pairs.  INPUTS: s; pairs.  OUTPUTS: s.
for k = 1:2:numel(varargin), s.(varargin{k}) = varargin{k+1}; end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
```

- [ ] **Step 2: Run to verify it fails** (MCP is enough: no solves) — expected missing `status` field.

- [ ] **Step 3: `sheet_to_catalog_file.m`.** After `Q.NOTE = repmat({''}, nD, nA);` add:

```matlab
% THE STATUS LAYER (optimality-status spec 3): every primary's status and
% reason and its junctions, and the ALTERNATIVES -- every other transfer the
% builders found above the floor, deduplicated against the primaries.
Q.STATUS = zeros(nD, nA, 1, 'int8');
Q.SREASON = repmat({''}, nD, nA);
Q.JUNC = cell(nD, nA);
alts = d('extraAlternatives', []);
```
In the spine loop: (a) immediately after the line `k = find([c.ok] & abs([c.tfDays] - S.TF(j)) < 1e-9, 1);` and BEFORE the line `if isempty(k), continue, end`, insert the loop below, which collects every candidate of the column except the winner `k` (so they are recorded whether or not the column has a winner); (b) after `Q.NOTE{iD0, j} = noteOf(...)` add `Q = putStatus(Q, iD0, j, c(k));`. The loop:

```matlab
    for m = 1:numel(c)
        if ~isempty(k) && m == k, continue, end
        alts = [alts, make_alternative(c(m), sD0, S.sA(j), iD0, j, sprintf('sheet candidate, column %d', j))];
    end
```
In the rib loop: before `if Q.OK(iD, iA, 1) && Q.TF(iD, iA, 1) <= Pt.z(8), continue, end` replace that line with:

```matlab
            if Q.OK(iD, iA, 1) && Q.TF(iD, iA, 1) <= Pt.z(8)
                alts = [alts, make_alternative(Pt, Pt.sD, Pt.sA, iD, iA, sprintf('rib point at sA %.4f, slower than the cell''s', R.sA))];
                continue
            end
            if Q.OK(iD, iA, 1)                     % the displaced primary is kept as an alternative
                alts = [alts, primaryAsAlternative(Q, iD, iA)];
            end
```
and after `Q.NOTE{iD, iA} = noteOf(Pt, F, ribCode);` add `Q = putStatus(Q, iD, iA, Pt);`. After the rib point loop (still inside `for k = 1:numel(ribs)`), add the rib's refusals:

```matlab
        if isfield(R, 'refused')
            for m = 1:numel(R.refused)
                Rf = R.refused(m);
                iD = idxOf(Q.sD, Rf.sD);  iA = idxOf(Q.sA, Rf.sA);
                if isempty(iD), iD = 0; end
                if isempty(iA), iA = 0; end
                alts = [alts, make_alternative(Rf, Rf.sD, Rf.sA, iD, iA, sprintf('rib refusal at sA %.4f (walked %d)', R.sA, Rf.walked))];
            end
        end
```
Before `% ---- meta`, build the temporary sheet view and deduplicate:

```matlab
view = struct('has_solution', Q.OK, 'entry_index', reshape(1:numel(Q.OK), size(Q.OK)), 'z8', reshape(Q.Z8, 8, []));
Q.ALT = dedup_alternatives(alts, view);
```
Add the local helpers:

```matlab
function Q = putStatus(Q, iD, iA, C)
% PUTSTATUS  A primary's status, reason and junctions.  INPUTS: Q; iD; iA;
% C (certify_root result; legacy results are classified).  OUTPUTS: Q.
if isfield(C, 'status') && isfield(C, 'stage'), st = C.status;  why = C.status_reason;
else, [st, why] = optimality_status(C); end
Q.STATUS(iD, iA, 1) = int8(st);  Q.SREASON{iD, iA} = why;
if isfield(C, 'Y'), Q.JUNC{iD, iA} = C.Y; end
end

function A = primaryAsAlternative(Q, iD, iA)
% PRIMARYASALTERNATIVE  The current primary of a cell as an alternatives row
% (it is being displaced by a faster root).  INPUTS: Q; iD; iA.  OUTPUTS: A.
C = struct('ok', true, 'stage', 8, 'status', double(Q.STATUS(iD, iA, 1)), 'status_reason', Q.SREASON{iD, iA}, ...
           'z', Q.Z8(:, iD, iA, 1), 'Y', Q.JUNC{iD, iA}, 'flyKm', 0, 'flyVms', 0, 'conj', double(Q.CONJ(iD, iA, 1)), ...
           'overridden', false, 'g', struct('minLamV', Q.MINLV(iD, iA, 1), 'minQmt', Q.MINQ(iD, iA, 1), 'dimS', Q.DIMS(iD, iA, 1)));
A = make_alternative(C, Q.sD(iD), Q.sA(iA), iD, iA, 'displaced primary (packaging)');
end
```
(`flyKm = 0` here: the displaced primary was certified, so its flight passed; the audit re-flies it.)

- [ ] **Step 4: `build_costate_catalog_family.m`.** In the per-entry loop add `reasons = repmat({''}, 1, nnz(Q.OK)); juncs = cell(1, nnz(Q.OK));` beside `notes`, and inside `if hasNotes && kr == 1 ...` add the parallel lines `if isfield(Q, 'SREASON') && kr == 1, reasons{n} = Q.SREASON{iD, iA}; end` and `if isfield(Q, 'JUNC') && kr == 1, juncs{n} = Q.JUNC{iD, iA}; end`. After the `family_index` line:

```matlab
    if isfield(Q, 'STATUS')
        sheets(nS,1).status        = Q.STATUS;
        sheets(nS,1).status_reason = reasons;
        sheets(nS,1).junctions     = juncs;
        if isfield(spec, 'statusKey') && ~isfield(cat_, 'status_key'), cat_.status_key = spec.statusKey; end
    end
    if isfield(Q, 'ALT') && ~isempty(Q.ALT)
        altK = Q.ALT;  [altK.sheet] = deal(nS);
        if ~isfield(cat_, 'alternatives') || isempty(cat_.alternatives), cat_.alternatives = altK;
        else, cat_.alternatives = [cat_.alternatives, altK]; end
    end
```
After the sheets are reordered (`cat_.sheets = sheets(ord);`) remap: `if isfield(cat_, 'alternatives'), inv(ord) = 1:numel(ord); for k = 1:numel(cat_.alternatives), cat_.alternatives(k).sheet = inv(cat_.alternatives(k).sheet); end, end` (name the variable `invOrd`, not `inv`).
`status_key` lives in `IND`; `build_costate_catalog_family` is in `CC` — add `addpath` of nothing: instead pass the legend in: in `package_phase_catalog.m` the `spec` struct for `build_costate_catalog_family` gains `'statusKey', status_key()`, and the family builder uses `spec.statusKey` when present (no costate_common -> campaign dependency).

- [ ] **Step 5: `package_phase_catalog.m`** — forward `opts.extraAlternatives` to `sheet_to_catalog_file` (it already passes `opts`) and add `'statusKey', status_key()` to the spec struct.

- [ ] **Step 6: Run** `test_package_status` (MCP), then batch `test_reproduce_library`, `test_score_interpolator`, `test_arrival_sheet_as_catalog`, `test_merge_phase_catalogs`. Expected ALL PASS.

- [ ] **Step 7: Mutation check** — (a) skip `make_alternative` for the sheet's losing candidates; (b) drop `primaryAsAlternative`; (c) omit `dedup_alternatives`. Each must fail `test_package_status`.

- [ ] **Step 8: Commit (host)** — `packager: status, status_reason, junctions per primary; alternatives table; status_key (spec 3, 5)`.

---

### Task 8: `merge_phase_catalogs` keeps a displaced primary and merges alternatives

**Files:**
- Modify: `IND/merge_phase_catalogs.m`
- Test: `IND/tests/test_merge_phase_catalogs.m` (extend)

**Interfaces:**
- Consumes: `make_alternative`, `dedup_alternatives`.
- Produces: `M.alternatives` = base alternatives + donor alternatives (cells remapped onto the base grid, off-grid -> 0) + every base primary the donor displaced (status from the base's `status` grid when present, `source = 'displaced primary (merge)'`), deduplicated against the merged primaries. Catalogs without the status layer merge exactly as today (Review Focus 2).

- [ ] **Step 1: Extend the test** — in `synth()` add `'status', int8(4*ones(nD, nA))`, `'status_reason', {repmat({'full stack passed'}, 1, n)}`, `'junctions', {repmat({ones(14, 24)}, 1, n)}` to the sheet, and a top-level `'alternatives', struct([])`; then append before the refusals:

```matlab
% ---- the status layer (optimality-status spec 5) ----------------------------
A = M.alternatives;
ok = chk(ok, numel(A) == 1 && A(1).iD == 1 && A(1).iA == 1 && A(1).status == 4 && strcmp(A(1).source, 'displaced primary (merge)') ...
         && isequal(A(1).z8, b.z8(:, b.entry_index(1, 1))), 'the displaced base primary is kept as an alternative');
ok = chk(ok, s.status(1, 3) == 4 && numel(s.junctions) == size(s.z8, 2), 'status grid and junctions travel with taken entries');
bOld = rmfield(base, 'alternatives');  bOld.sheets = rmfield(bOld.sheets, {'status', 'status_reason', 'junctions'});
dOld = rmfield(donor, 'alternatives');  dOld.sheets = rmfield(dOld.sheets, {'status', 'status_reason', 'junctions'});
[Mo, io] = merge_phase_catalogs(bOld, dOld);
ok = chk(ok, io.nFaster == 1 && io.nFilled == 1 && ~isfield(Mo, 'alternatives'), 'REVIEW FOCUS 2: catalogs without the layer merge as before');
```

- [ ] **Step 2: Run to verify it fails** (MCP).

- [ ] **Step 3: Implement.** In `merge_phase_catalogs.m`:
  - after the field discovery, `hasLayer = isfield(sb, 'status') && isfield(sd, 'status');` and `alts = struct([]);`
  - inside the take branch, BEFORE the entry fields are overwritten and only when `~filled && hasLayer`, record the displaced primary:

```matlab
        if ~filled && hasLayer
            kb0 = sb.entry_index(iD, iA, 1);
            Cp = struct('ok', true, 'stage', 8, 'status', double(sb.status(iD, iA, 1)), ...
                        'status_reason', sb.status_reason{kb0}, 'z', sb.z8(:, kb0), 'Y', sb.junctions{kb0}, ...
                        'flyKm', 0, 'flyVms', 0, 'overridden', false);
            alts = [alts, make_alternative(Cp, sb.sD_frac(iD), sb.sA_frac(iA), iD, iA, 'displaced primary (merge)')];
        end
```
  - after the loop, when `hasLayer`: remap donor alternatives (`iD = mapD(a.iD)` when `a.iD > 0`, else 0; same for `iA`), concatenate `[base.alternatives, remapped, alts]` (skipping empties), `M.alternatives = dedup_alternatives(all, M.sheets)`; `M.status_key = base.status_key` when present.

- [ ] **Step 4: Run** `test_merge_phase_catalogs` (MCP) — ALL PASS, including the 16 earlier checks.

- [ ] **Step 5: Mutation check** — (a) record the displaced primary AFTER overwriting (it would carry the donor's z8); (b) skip the remap. Each must fail.

- [ ] **Step 6: Commit (host)** — `merge_phase_catalogs: displaced primaries kept as alternatives; alternatives merged and remapped (spec 5)`.

---

### Task 9: `alternatives_content_key` and the audit `audit_status_layer`

**Files:**
- Create: `IND/alternatives_content_key.m`, `IND/audit_status_layer.m`
- Test: `IND/tests/test_audit_status_layer.m`

**Interfaces:**
- Produces: `key = alternatives_content_key(c)` — MD5 hex (32 chars) of every alternative's `sD sA status z8` in table order (same hashing as `catalog_content_key`: `java.security.MessageDigest`); `''` when the catalog has no alternatives.
- Produces: `A = audit_status_layer(catMat, opts)` — FAIL CLOSED. `opts`: `.idxAlt` [] (subset of alternatives; default all), `.skipPrimaries` [false], `.out` '' save path, `.pool` [capped_pool()], `.tolMove` [1e-6], `.junctionKm` [1] max junction defect, `.flyKm` [1] arrival miss for primaries, `.certifier` [] TEST SEAM `@(seed, rv0, rvf, B, opts) -> C` replacing `certify_root`. Returns `.primRows` (per primary: `ok`, `why`), `.altRows` (per alternative: `ok`, `why`, `statusNow`, `moved`), `.nOk .nBad .problems`, `.contentKey` (`catalog_content_key`), `.altContentKey`.
- Primary checks: `status == 4`; junctions present, `entry_thrust_program` defects `<= junctionKm`, arrival miss `<= flyKm`, `startErr < 1e-9`.
- Alternative checks: rebuild `rv0, rvf` from the row's phases via `B = arclength_arrival('setup', catalog_setup_request(c, sD))`; seed = `struct('tf', z8(8), 'tGrid', linspace(0, z8(8), K+1), 'Y', [junctions, junctions(:,end)])`; run the certifier; BAD if the root moved (`~same_root(C.z, row.z8)` — "moved to another root", Review Focus 5) or `C.status ~= row.status` ("status not reproduced: stored s, now t").

- [ ] **Step 1: Write the failing test** — fast, through the certifier seam (no solves) plus one real primary check.

```matlab
function ok = test_audit_status_layer()
% TEST_AUDIT_STATUS_LAYER  The status audit fails closed: a reproduced
% status passes; a different status, a root that moved, and corrupted
% junctions are each a BAD row naming why; the content keys are bound.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
L = load(fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat'));  c = L.(char(fieldnames(L)));
sh = c.sheets(1);  k = sh.entry_index(1, 1);
row = struct('sD', sh.sD_frac(1), 'sA', sh.sA_frac(1), 'iD', 0, 'iA', 0, 'z8', sh.z8(:, k), 'junctions', ones(14, 24), ...
             'tf_nd', sh.z8(8, k), 'status', 3, 'status_reason', 'x', 'inferred', false, 'conj', 1, 'conjVerdict', 'PASS', ...
             'minLamV', 1, 'minQmt', 1, 'dimS', 1, 'h6Margin', 1, 'liftMargin', 1, 'flyKm', 0, 'flyVms', 0, 'source', 't', 'sheet', 1);
c.alternatives = [row, setf(row, 'status', 2), row];
fake = @(st, dz) @(seed, rv0, rvf, B, o) struct('z', [seed.Y(8:14, 1); seed.tf] + dz, 'status', st, 'reason', 'fake');
tmp = [tempname '.mat'];  catalog = c;  save(tmp, 'catalog');
A = audit_status_layer(tmp, struct('skipPrimaries', true, 'idxAlt', 1, 'certifier', fake(3, 0), 'pool', []));
ok = chk(ok, A.nBad == 0 && A.altRows(1).ok, 'a reproduced status passes');
A = audit_status_layer(tmp, struct('skipPrimaries', true, 'idxAlt', 2, 'certifier', fake(3, 0), 'pool', []));
ok = chk(ok, A.nBad == 1 && contains(A.altRows(1).why, 'status not reproduced'), 'a different status is BAD, named');
A = audit_status_layer(tmp, struct('skipPrimaries', true, 'idxAlt', 3, 'certifier', fake(3, [0.01; zeros(7, 1)]), 'pool', []));
ok = chk(ok, A.nBad == 1 && A.altRows(1).moved && contains(A.altRows(1).why, 'moved'), 'REVIEW FOCUS 5: a root that moved is BAD, not relabelled');
ok = chk(ok, numel(A.altContentKey) == 32 && ~strcmp(A.altContentKey, alternatives_content_key(setf(c, 'alternatives', row))), ...
         'the alternatives key is bound to the table''s content');
delete(tmp);
if ok, fprintf('test_audit_status_layer: ALL PASS\n'); else, fprintf('test_audit_status_layer: FAIL\n'); end
end

function s = setf(s, varargin)
% SETF  Set name/value pairs.  INPUTS: s; pairs.  OUTPUTS: s.
for k = 1:2:numel(varargin), s.(varargin{k}) = varargin{k+1}; end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
```

- [ ] **Step 2: Run to verify it fails** (MCP).

- [ ] **Step 3: `alternatives_content_key.m`** — copy `catalog_content_key.m`'s hashing (read it first: `sed -n '/Begin Code/,$p' IND/catalog_content_key.m`), hashing `[A(k).sD; A(k).sA; A(k).status; A(k).z8(:)]` for every row in order, typecast to bytes exactly as `catalog_content_key` does.

- [ ] **Step 4: `audit_status_layer.m`** — implement the checks listed under Interfaces. Structure:

```matlab
L = load(catMat);  c = L.(char(fieldnames(L)));
assert(isscalar(c.sheets), 'audit_status_layer:sheets', 'one sheet expected');
sh = c.sheets(1);  P = struct();
problems = {};  primRows = struct('k', {}, 'ok', {}, 'why', {});  altRows = struct('k', {}, 'ok', {}, 'why', {}, 'statusNow', {}, 'moved', {});
% ---- primaries -------------------------------------------------------------
if ~skipPrimaries
    assert(isfield(sh, 'status') && isfield(sh, 'junctions'), 'audit_status_layer:noLayer', ...
           'the catalog has no status layer (status, junctions): nothing to audit');
    [B, phys] = setupAt(c, sh.sD_frac(1));
    for iA = 1:numel(sh.sA_frac)
        for iD = 1:numel(sh.sD_frac)
            if ~sh.has_solution(iD, iA, 1), continue, end
            k = sh.entry_index(iD, iA, 1);  why = '';
            if sh.status(iD, iA, 1) ~= 4, why = sprintf('primary status %d, not 4', sh.status(iD, iA, 1)); end
            Jn = sh.junctions{k};
            if isempty(why) && isempty(Jn), why = 'no junctions'; end
            if isempty(why)
                rv0 = B.stateD(sh.sD_frac(iD));  rvf = B.stateA(sh.sA_frac(iA));
                Pp = entry_thrust_program(sh.z8(:, k), Jn, rv0(1:6), rvf(1:6), phys);
                if max([Pp.defectKm 0]) > junctionKm, why = sprintf('junction defect %.2f km', max(Pp.defectKm));
                elseif Pp.flyKm > flyKm, why = sprintf('junction flight misses by %.2f km', Pp.flyKm);
                elseif Pp.startErr > 1e-9, why = sprintf('first junction differs from z8 (%.1e)', Pp.startErr); end
            end
            primRows(end+1) = struct('k', k, 'ok', isempty(why), 'why', why);
            if ~isempty(why), problems{end+1} = sprintf('primary (%d,%d): %s', iD, iA, why); end
        end
    end
end
% ---- alternatives -------------------------------------------------------------
A = struct([]);  if isfield(c, 'alternatives'), A = c.alternatives; end
idx = 1:numel(A);  if ~isempty(idxAlt), idx = idxAlt; end
for k = idx
    a = A(k);  K = size(a.junctions, 2);
    [B, ~] = setupAt(c, a.sD);
    rv0 = B.stateD(a.sD);  rvf = B.stateA(a.sA);
    seed = struct('tf', a.z8(8), 'tGrid', linspace(0, a.z8(8), K + 1), 'Y', [a.junctions, a.junctions(:, end)]);
    seed.Y(1:7, 1) = [rv0(1:6); 1];  seed.Y(8:14, 1) = a.z8(1:7);
    try
        C = certifier(seed, rv0(1:6), rvf(1:6), B, struct('sD', a.sD, 'sA', a.sA, 'pool', pool, 'allowUnfenced', isempty(pool)));
        moved = ~same_root(C.z, a.z8);
        if moved, why = 'the re-certification moved to another root';
        elseif C.status ~= a.status, why = sprintf('status not reproduced: stored %d, now %d (%s)', a.status, C.status, C.reason);
        else, why = ''; end
        sNow = C.status;
    catch ME
        moved = false;  sNow = NaN;  why = ['re-certification threw: ' ME.message];
    end
    altRows(end+1) = struct('k', k, 'ok', isempty(why), 'why', why, 'statusNow', sNow, 'moved', moved);
    if ~isempty(why), problems{end+1} = sprintf('alternative %d (%.4f, %.4f): %s', k, a.sD, a.sA, why); end
    if ~isempty(outFile), save(outFile, 'primRows', 'altRows', 'problems'); end    % resume-friendly
end
```
with `certifier = fieldd(opts, 'certifier', @certify_root)`, defaults as listed, and a local `setupAt(c, sD)` returning `[B, phys]` from `arclength_arrival('setup', catalog_setup_request(c, sD))` and `phys = struct('Tnd', B.Tnd, 'cnd', B.cnd, 'mu', B.mu, 'lStar', c.constants.lStar_km, 'tStar', c.constants.tStar_s, 'm0kg', c.thruster.m0_kg)`. The output struct: `nOk`, `nBad = numel(problems)`, `contentKey = catalog_content_key(c)`, `altContentKey = alternatives_content_key(c)`, saved to `opts.out` when given. Note `'allowUnfenced', isempty(pool)` exists only so the seam test can run without a pool; production calls always pass a pool.

- [ ] **Step 5: Run** `test_audit_status_layer` (MCP) — ALL PASS.

- [ ] **Step 6: Mutation check** — (a) drop the `moved` test; (b) compare `status` with `>=`; (c) key without `status`. Each must fail.

- [ ] **Step 7: Commit (host)** — `audit_status_layer + alternatives_content_key: status reproduction, fail closed (spec 7)`.

---

### Task 10: `transfers_at` — every transfer at a phase pair

**Files:**
- Create: `IND/transfers_at.m`
- Test: `IND/tests/test_transfers_at.m`

**Interfaces:**
- Produces: `T = transfers_at(c, sD, sA, opts)` — `opts.tol` [1e-6] circular phase tolerance. Struct array sorted by `tfDays`, fields: `kind` ('primary' | 'alternative'), `sD sA status statusName reason tfDays dvKms propellantKg z8 junctions source`. A catalog without the layer returns its primary with `status = NaN`, `statusName = 'not annotated'` (Review Focus 2). `dvKms = c_nd*log(1/mf)*lStar/tStar`, `mf = 1 - Tmax_nd*tf/c_nd`, `Tmax_nd = (rungs_N/m0_kg)*tStar^2/(lStar_km*1000)` (the catalog's `derive` formulas).

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_transfers_at()
% TEST_TRANSFERS_AT  Every transfer at a phase pair, ranked by flight time:
% the primary and the alternatives in its cell and off the grid; the wrap
% s_D = 1 - 1e-9 finds s_D = 0; a catalog without the status layer reports
% its primary as 'not annotated'; Delta-V follows the catalog's formula.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
L = load(fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat'));  c = L.(char(fieldnames(L)));
sh = c.sheets(1);  k = sh.entry_index(1, 1);  tS = c.constants.tStar_s/86400;
T = transfers_at(c, sh.sD_frac(1), sh.sA_frac(1));
ok = chk(ok, numel(T) == 1 && strcmp(T.kind, 'primary') && isnan(T.status) && strcmp(T.statusName, 'not annotated'), ...
         'REVIEW FOCUS 2: no layer -> the primary, not annotated');
Tmax = (c.rungs_N(1)/c.thruster.m0_kg)*c.constants.tStar_s^2/(c.constants.lStar_km*1000);
mf = 1 - Tmax*sh.z8(8, k)/c.thruster.c_nd;
ok = chk(ok, abs(T.dvKms - c.thruster.c_nd*log(1/mf)*c.constants.lStar_km/c.constants.tStar_s) < 1e-12, 'Delta-V from the catalog''s formula');
c.sheets(1).status = int8(4*sh.has_solution);  c.sheets(1).status_reason = repmat({'full stack passed'}, 1, size(sh.z8, 2));
c.sheets(1).junctions = repmat({[]}, 1, size(sh.z8, 2));
alt = @(z8, st, sDv) struct('sD', sDv, 'sA', sh.sA_frac(1), 'iD', 1, 'iA', 1, 'z8', z8, 'junctions', [], 'tf_nd', z8(8), ...
        'status', st, 'status_reason', 'r', 'inferred', false, 'conj', 0, 'conjVerdict', 'FAIL', 'minLamV', 1, 'minQmt', 1, ...
        'dimS', 1, 'h6Margin', 1, 'liftMargin', 1, 'flyKm', 0, 'flyVms', 0, 'source', 's', 'sheet', 1);
zS = sh.z8(:, k);  zF = zS;  zF(8) = zS(8) - 0.01;  zL = zS;  zL(8) = zS(8) + 0.02;
c.alternatives = [alt(zL, 2, sh.sD_frac(1)), alt(zF, 3, sh.sD_frac(1)), alt(zL, 1, 0.5)];
T = transfers_at(c, 1 - 1e-9, sh.sA_frac(1));
ok = chk(ok, numel(T) == 3 && issorted([T.tfDays]) && strcmp(T(2).kind, 'primary') && T(1).status == 3 && T(3).status == 2, ...
         'REVIEW FOCUS 4: the wrap finds s_D = 0; ranked by t_f: necessary-only, primary, conjugate');
ok = chk(ok, strcmp(T(3).statusName, 'conjugate point found'), 'status names from status_key');
if ok, fprintf('test_transfers_at: ALL PASS\n'); else, fprintf('test_transfers_at: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
```
(`sh.sD_frac(1)` is 0 in this catalog; the third alternative at s_D 0.5 must not appear.)

- [ ] **Step 2: Run to verify it fails** (MCP).

- [ ] **Step 3: Implement** `transfers_at.m` (house header) — loop the single sheet's cell matching `(sD, sA)` within `tol` (circular: `abs(mod(a - b + 0.5, 1) - 0.5) <= tol`), add the primary if present; add every alternative whose `sD`, `sA` match within `tol`; build each row with a local `row(kind, sD, sA, z8, status, reason, junctions, source)` computing `tfDays`, `dvKms`, `propellantKg = (1 - mf)*m0_kg`; `statusName` from `status_key()` (`names{codes == status}`), `'not annotated'` for NaN; sort by `tfDays`.

- [ ] **Step 4: Run** — ALL PASS. **Step 5: Mutation check** — (a) non-circular phase match; (b) sort removed. **Step 6: Commit (host)** — `transfers_at: every transfer at a phase pair, ranked (spec 8)`.

---

### Task 11: The backfill — `library_70mN_24x48_v2`

**Files:**
- Create: `IND/backfill_status_layer.m` (function, stages), `IND/recertify_candidates.m` (chunked, resumable), `IND/batch/backfill_v2_job.m`, `IND/batch/recertify_chunk_job.m`, `IND/batch/run_recertify.sh`
- Test: `IND/tests/test_backfill_status_layer.m`

**Interfaces:**
- Consumes: Tasks 1-10.
- `L = backfill_status_layer('harvest', recordMat, sources, outDir)` → writes `outDir/harvest.mat` with `cands` (struct array: certify_root results + `.sD .sA .source`), `primJ` ({1 x nEntries} junctions matched by `same_root`, [] when unmatched), `needRecert` (indices of cands whose status is inferred 3 from a legacy "conjugate test verdict" or whose re-solve is required: the rib stop point), `unmatched` (primary entry indices without a source `Y`).
- `recertify_candidates(harvestMat, chunk, nChunk, outMat)` — for its share of `needRecert` (and of `unmatched` primaries, re-polished from z8 via `seed_from_z8` + `certify_root`), certify under the new certifier, saving after every item (resume: skip items already in `outMat`). A primary whose re-polish leaves `same_root` is recorded `moved` and NOT written (Review Focus 5).
- `c2 = backfill_status_layer('assemble', recordMat, outDir)` → reads `harvest.mat` + every `recert_*.mat`, builds the v2 catalog: primaries unchanged, `status` 4 grid, `status_reason` 'full stack passed', `junctions` (matched or re-polished), `alternatives = dedup_alternatives(make_alternative rows, sheet)`, `status_key`; saves `outDir/costate_catalog_dro_tulip_70mN.mat`; asserts `catalog_content_key(c2) == catalog_content_key(record)`.
- Sources (from `RES`): `reproduce_70mN_24x48/round_01/arrival_sheet_70mN_nA48.mat`, `reproduce_70mN_24x48/round_01/fine_rib_col*.mat`, `reproduce_70mN_24x48/round_01/fine_rib_direct_holes.mat`, `reproduce_70mN_24x24/round_01/arrival_sheet_70mN_nA24.mat`, `reproduce_70mN_24x24/round_01/fine_rib_col*.mat`, `reproduce_70mN_24x24/round_01/fine_rib_direct_holes.mat`, `sheet96_resolution_test/arrival_sheet_70mN_nA96.mat`, `departure_rib96_test/rib96_col28.mat`, `departure_rib96_test/rib96_col28_ckpt.mat`, and the stop-point re-solve (the rib's last point at s_D 17/96 stepped to 0.1716 with `certify_root`).

- [ ] **Step 1: Write the failing test** — a 2 x 2 synthetic record plus synthetic source files in a temp folder; `harvest` must match junctions by root, classify a legacy refusal (inferred) and list a legacy "conjugate test verdict 0" in `needRecert`; `assemble` with a fake recert file must keep primaries bit-identical (content key), give every primary its junctions, and produce the expected alternatives.

```matlab
function ok = test_backfill_status_layer()
% TEST_BACKFILL_STATUS_LAYER  harvest matches primaries' junctions by root,
% classifies legacy candidates and lists the ambiguous ones for
% re-certification; assemble keeps the primaries bit-identical and builds
% the alternatives from harvest + re-certification.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
L = load(fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat'));  rec = L.(char(fieldnames(L)));
sh = rec.sheets(1);  k11 = sh.entry_index(1, 1);  k12 = sh.entry_index(1, 2);
tmp = tempname;  mkdir(tmp);
catalog = rec;  recMat = fullfile(tmp, 'record.mat');  save(recMat, 'catalog');
legacy = @(z, reason, ok_) struct('ok', ok_, 'z', z, 'Y', z(1)*ones(14, 24), 'flyKm', 0.1, 'flyVms', 0.01, 'reason', reason, ...
                                   'tfDays', z(8), 'sD', sh.sD_frac(1), 'sA', sh.sA_frac(1), 'conj', 0);
zAlt = sh.z8(:, k11) .* [1.01; ones(7, 1)];
S = struct('problem', struct('sD', 0), 'sA', sh.sA_frac(1:2), 'TF', [1 1], ...
           'cand', {{[legacy(sh.z8(:, k11), 'certified', true), legacy(zAlt, 'conjugate test verdict 0', false)], ...
                     legacy(sh.z8(:, k12), 'certified', true)}});
save(fullfile(tmp, 'sheet.mat'), 'S');
H = backfill_status_layer('harvest', recMat, {fullfile(tmp, 'sheet.mat')}, tmp);
ok = chk(ok, isequal(H.primJ{k11}, sh.z8(1, k11)*ones(14, 24)) && isequal(H.primJ{k12}, sh.z8(1, k12)*ones(14, 24)), ...
         'primaries get their junctions by root');
ok = chk(ok, numel(H.needRecert) == 1 && same_root(H.cands(H.needRecert).z, zAlt), 'the legacy verdict-0 candidate is listed for re-certification');
ok = chk(ok, numel(H.unmatched) == nnz(sh.has_solution) - 2, 'every other primary is listed as unmatched (needs a re-polish)');
% a fake re-certification: the verdict-0 candidate is a real FAIL with hypotheses held
Cr = H.cands(H.needRecert);  Cr.stage = 7;  Cr.status = 2;  Cr.status_reason = 'conjugate point found (gates and H6 held): x';
Cr.conjFound = true;  Cr.hypAfterConj = 'held';  Cr.overridden = false;
items = struct('kind', 'cand', 'index', H.needRecert, 'C', Cr, 'moved', false);
save(fullfile(tmp, 'recert_1.mat'), 'items');
c2 = backfill_status_layer('assemble', recMat, tmp, struct('allowUnmatched', true));
ok = chk(ok, strcmp(catalog_content_key(c2), catalog_content_key(rec)), 'primaries bit-identical (content key)');
ok = chk(ok, numel(c2.alternatives) == 1 && c2.alternatives.status == 2 && ~c2.alternatives.inferred, ...
         'the alternative carries the re-certified status, not the inferred one');
ok = chk(ok, all(c2.sheets(1).status(sh.has_solution) == 4), 'every primary is status 4');
rmdir(tmp, 's');
if ok, fprintf('test_backfill_status_layer: ALL PASS\n'); else, fprintf('test_backfill_status_layer: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
```
(`opts.allowUnmatched` lets `assemble` run in the test with primaries lacking junctions; in production it is false and an unmatched, un-repolished primary makes `assemble` refuse.)

- [ ] **Step 2: Run to verify it fails** (MCP).

- [ ] **Step 3: Implement `backfill_status_layer.m`** — one function, `mode` first argument (`'harvest'` | `'assemble'`), house header, local functions per mode:
  - **harvest:** load the record; for each source file load its struct (`S` for sheets: candidates from `S.cand{j}` with phases `(S.problem.sD, S.sA(j))`; `R` for rib files: `R.pts` and, if present, `R.refused`, phases from each point; the filler file's `R` the same way and `rec(k).others` when present; a checkpoint file's `pts`/`refused`); stamp each with `.source`; for every primary entry find the first candidate with `same_root(cand.z, z8)` and a non-empty `Y` → `primJ{k}`; `needRecert` = candidates whose `optimality_status` reason starts `'necessary only (inferred): legacy verdict 0'`; save `harvest.mat` with `cands primJ needRecert unmatched`.
  - **assemble:** load harvest + all `recert_*.mat` (`items` arrays: `kind` 'cand' | 'prim', `index`, `C`, `moved`); replace re-certified candidates; fill `primJ` from `'prim'` items that did not move; refuse if any primary still lacks junctions unless `opts.allowUnmatched`; build rows with `make_alternative` (iD/iA by matching the record's grid, 0 off-grid), `dedup_alternatives` against the record's sheet; write `status`, `status_reason`, `junctions`, `alternatives`, `status_key`, and `c2.status_layer = struct('built', date, 'sources', {sources}, 'nRecert', ..., 'nMoved', ...)`; assert the content key; save `outDir/costate_catalog_dro_tulip_70mN.mat` (variable name as the record's).

- [ ] **Step 4: Implement `recertify_candidates.m`** — items = `needRecert` (kind 'cand': seed from the candidate's own `Y` via the Task 9 seed construction, phases from the candidate) followed by `unmatched` (kind 'prim': seed from `seed_from_z8(z8, rv0, 24, ...)`); take items `chunk:nChunk:end`; skip those already in `outMat`; `C = certify_root(seed, rv0, rvf, B, struct('sD', sD, 'sA', sA, 'pool', capped_pool(1)))`; `moved = ~same_root(C.z, z8)`; append to `items` and save after every item (`publish_atomic` pattern: save to `[outMat '.part']`, then `movefile`). Log one line per item to `[outMat '.log']`.

- [ ] **Step 5: Jobs.** `batch/backfill_v2_job.m` runs `harvest` (sources above, `outDir = RES/library_70mN_24x48_v2`), then prints `numel(needRecert)` and `numel(unmatched)` and the command to launch the chunks. `batch/recertify_chunk_job.m` reads `CHUNK`, `NCHUNK` from the environment and calls `recertify_candidates`. `batch/run_recertify.sh N` launches N `nohup matlab -batch` processes (one per chunk, outputs `recert_<k>.mat`), each with its own log. A third invocation of `backfill_v2_job.m` with `BACKFILL_STAGE=assemble` assembles, runs `audit_status_layer` in N chunks (same chunking via `.idxAlt`) plus the primaries pass, runs `compare_phase_catalogs(v2, record)` and writes `~/BACKFILL_V2_VERDICT.txt` with: entries, alternatives by status, re-certified, moved, content key equal (yes/no), compare (agree / missing / slower / faster), audit ok/bad.

- [ ] **Step 6: Run** `test_backfill_status_layer` (MCP) — ALL PASS. **Mutation check** — (a) match junctions by t_f only; (b) let `assemble` keep the inferred status over a re-certified one. Each must fail.

- [ ] **Step 7: Commit (host)** — `backfill_status_layer + recertify_candidates + jobs: the status layer for the library of record without a rebuild (spec 6)`.

---

### Task 12: Run the backfill, document, hand the adoption decision to Mike

**Files:**
- Modify: `DRO_tulip/FINDINGS.md` (section 98), `DRO_tulip/README.md` (status-layer paragraph), `DRO_tulip/process/INTERPOLATION_STATUS.md` (pointer), `DRO_tulip/doc/reproduce_library_sdd.tex` (a "Status layer" section: data model of spec 3, classifier table of spec 4.2, audit of spec 7), `DRO_tulip/TODO.md`.

- [ ] **Step 1: Harvest** — `nohup matlab -batch "run('<IND>/batch/backfill_v2_job.m')" > ~/backfill_v2.out 2>&1 &`; record `needRecert` and `unmatched` counts.
- [ ] **Step 2: Re-certify** — `bash <IND>/batch/run_recertify.sh 8`; monitor each chunk's log age (stall = no new line for 40 min) and process liveness; resume by rerunning the same chunk (it skips done items).
- [ ] **Step 3: Assemble + audit + compare** — `BACKFILL_STAGE=assemble nohup matlab -batch "run('<IND>/batch/backfill_v2_job.m')" ...`. Acceptance (all required): content key equal; compare: 0 missing, 0 slower, 0 faster, 0 costate changes; `audit_status_layer` 0 bad; `audit_phase_catalog` on v2 still 1,152 ok / 0 bad.
- [ ] **Step 4: Docs** — FINDINGS 98 with the verdict numbers (alternatives by status, how many conjugate refusals were FAIL vs UNDETERMINED, moved count, where the status-2 transfers sit, the s_D 0.1716 stop point's status); README / SDD / status page / TODO pointers. Compile the SDD with `/Library/TeX/texbin/pdflatex` and clean aux files.
- [ ] **Step 5: Commit (host)** and report to Mike: the verdict line, the counts, and the adoption question (v2 replaces the merged library as the record, the FINDINGS 95 pattern).

---

## Self-review notes (kept for the executor)

- Spec coverage: 3.1/3.2/3.3 -> Tasks 7, 8, 11; 4 -> Tasks 1, 2; 5 -> Tasks 2, 5, 6, 7, 8; 6 -> Tasks 11, 12; 7 -> Task 9 (+ the unchanged `audit_phase_catalog`); 8 -> Tasks 4, 10; 9 -> every task's tests; 10 -> Task 12.
- Names used across tasks: `optimality_status`, `status_key`, `same_root`, `catalog_setup_request`, `make_alternative`, `dedup_alternatives`, `entry_thrust_program`, `alternatives_content_key`, `audit_status_layer`, `transfers_at`, `backfill_status_layer`, `recertify_candidates`; result fields `stage conjVerdict conjFound conjReason hypAfterConj overridden status status_reason`; rib field `refused` (+ `walked`); filler field `others`; sheet-file fields `STATUS SREASON JUNC ALT`; catalog fields `status status_reason junctions alternatives status_key status_layer`.
