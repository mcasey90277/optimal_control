# Costate-Library Tools and Deliverable 8 — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship the 70 mN DRO -> tulip library of record (`library_70mN_24x48_v2`, with the optimality-status layer) to Darin as deliverable 8, with a tool set that picks, lists, flies, checks and surveys its transfers, built in the repo and packaged by the builder.

**Architecture:** New fly-layer functions (`catalog_endpoints`, `fly_catalog_entry`, `tfmin_accepts`) wrap the existing engines (`get_family_orbit` + `phase_state`, `entry_thrust_program`, `fly_transfer`, `pumpkyn.cr3bp.tfMin`). The survey helpers (`costate_lib_describe`, `costate_catalog_extremes`) become status-aware with defaults unchanged; `plot_status_torus` is new. `export_deliverable_catalog` writes the cleaned catalog with the audited content keys. `build_dro_deliverable` assembles the folder from ONE file list, writes the README from computed facts, and refuses to zip without a CLEAN audit bound by content key. A recipient self-test proves the folder runs with only pumpkyn.

**Tech Stack:** MATLAB R2026a (`/Applications/MATLAB_R2026a.app/bin/matlab`), pumpkyn + pumpkynPie (`startup` at `/Users/msc/Desktop/proj7/external/pumpkynPie`), Parallel Computing Toolbox for the tfMin fence. Tests are function files returning `ok`, printing PASS/FAIL per check.

**Spec:** `docs/superpowers/specs/2026-10-13-costate-library-tools-deliverable8-design.md` (read it first).

**Paths** (under `/Users/msc/Desktop/optimal_control/orbit_transfer/`): `IND` = `DRO_tulip/indirect/`, `CC` = `costate_common/`, `V2` = `DRO_tulip/indirect/results/library_70mN_24x48_v2/costate_catalog_dro_tulip_70mN.mat`, `MERGED` = `DRO_tulip/indirect/results/library_70mN_24x48_merged/`, `DEL` = `DRO_tulip/deliverables/costate_catalog_70mN_status_deliverable8/`.

**Measured facts the plan relies on** (2026-10-13): V2 primaries content key `6a775f714b4cdf302a505aa56a8198ad`, alternatives key `3a1c997b6def0ebabc1bfe044cb2ccd3`; `MERGED/audit_70mN.mat` (variable `A`: `.nOk` 1152, `.nBad` 0, `.contentKey` = the primaries key, `.catMat` = the merged path); `V2`'s folder holds `audit_alt_1..8.mat` (variables `contentKey`, `altContentKey`, `audit` (`.nOk .nBad .primRows .altRows .problems`), `primRows`, `altRows`); the V2 catalog's top-level fields: `name description created provenance(char) schema constants thruster rungs_N notes_key families sheets n_entries env derive conj_test hyp_gates usage second_order merge alternatives status_key status_layer`; alternatives rows: `sD sA iD iA z8 junctions tf_nd status status_reason inferred conj conjVerdict minLamV minQmt dimS h6Margin liftMargin flyKm flyVms source sheet borderline`; `transfers_at` rows: `kind sD sA status statusName reason tfDays dvKms propellantKg z8 junctions source`; `costate_catalog_extremes` extreme struct: `tauDRO dep_family dep_params arr_family arr_params arrLab Np thrustN dep_frac dep_nd dep_days arr_frac arr_nd arr_days tf_nd tf_days deltaV_kms propellant_kg z8`; `print_transfer_summary(S)` reads `.tfDays .tfNd .dvKms .propellantKg .admissibility(.reason .moonKm .earthKm) .nSamples`; `plot_transfer_3d(T, B, opts)` reads `T.z T.sD T.sA`, `B.stateD B.stateA B.Tnd B.cnd B.problem(.lStar .tStar .muStar)`, `opts.flight(.t .Y)`; `pumpkyn.cr3bp.tfMin(rv0, rvf, z8, Tnd, cnd, mu)` returns z8.

## Global Constraints

- MATLAB header on every new function: `%% Purpose:` / `%% Inputs:` (name, type, size) / `%% Outputs:` / `%% Revision History:` with `%  M. Casey ... (c) 10/13/2026` and `%  Copyright Coorbital Inc.`, then `%% ------------------------ Begin Code Sequence ---------------------------` (copy `IND/rib_as_catalog.m`'s layout). Changed functions get a dated Revision History line.
- Never `i`/`j` as variables; never `%#ok` pragmas; `sqrt(sum(x.^2))` not `norm` in new code.
- Phase conventions: FRACTIONS of a period (wrapping, circular tolerance) in `transfers_at`, `catalog_endpoints`, `fly_catalog_entry`, `tfmin_accepts`, `plot_status_torus`; DAYS in `costate_catalog_pick` (unchanged). ND unless the field name says otherwise (`tfDays`, `dvKms`, `massKg`, `flyKm`, `flyVms`).
- Existing calls keep their behaviour: `costate_catalog_pick(cat, tauDep, Np, depDays, arrDays, thrustN[, warn])`, `costate_lib_describe(cat)`, `costate_catalog_extremes(cat, filt, doPlot, metric)` with no new filter fields -> identical output.
- Every solver/flight call outside a unit test runs in `matlab -batch`, fenced where a pool exists; never in the shared MATLAB session. Test command form: `/Applications/MATLAB_R2026a.app/bin/matlab -batch "addpath('<CC>'); capped_pool(); cd('<tests dir>'); exit(~<test>())"` (judge by MATLAB's exit code; `cd` inside the command is required). Pure tests may use the MATLAB MCP `evaluate_matlab_code`.
- The exported catalog's `catalog_content_key` and `alternatives_content_key` must equal V2's (asserted in code and in a test).
- Deliverable example scripts load the catalog by BARE filename; no absolute paths anywhere in `DEL`.
- The shipped folder may call only files inside `DEL` and pumpkyn/pumpkynPie (`+cr3bp`, `+util`, `pumpkynPie.*`, `movie2gifLite`); proven by the dependency scan in Task 10.
- Commits: one per reviewed task, by the controller, message ending with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`.

## Review Focus

1. A catalog WITHOUT the status layer (deliverables 3-7 files; the 24 x 24 reference) passed to the changed helpers must behave exactly as today — pinned in Tasks 5 and 6 (byte-equal printout / identical extremes).
2. A phase query at the wrap (s = 1 - 1e-9) in `catalog_endpoints` and `fly_catalog_entry` must give the s = 0 state — pinned in Task 2.
3. A `transfers_at` row with EMPTY junctions (a legacy row) passed to `fly_catalog_entry` with the default method must fall back to the z8 method and say so, not error — pinned in Task 3.
4. `tfmin_accepts` on a status-1 row must return `ok = false, expected = true` and never throw; a fence timeout must be a named refusal — pinned in Task 4.
5. The builder must refuse to zip when the audit files' keys do not match the export (stale audit) — pinned in Task 10.

---

### Task 1: `export_deliverable_catalog`

**Files:**
- Create: `IND/export_deliverable_catalog.m`
- Test: `IND/tests/test_export_deliverable_catalog.m`

**Interfaces:**
- Produces: `[c, info] = export_deliverable_catalog(recordMat, outMat, opts)`; `opts.varName` ['costate_catalog_dro_tulip_70mN']. Returns the exported struct and `info` (`.primKey .altKey .removed` (cellstr of top-level fields removed) `.nSourceRewritten`).
- Catalog changes: top-level `merge` and `status_layer` REMOVED; `provenance` (char) extended with one paragraph (build -> refill -> merge -> backfill -> re-certification -> audit CLEAN, with dates); new `build` struct: `.dates` (char), `.nCandidates`, `.nRecert`, `.nMoved`, `.nBorderline`, `.borderlineNote` (from `status_layer.borderlineNote`), `.auditVerdict` 'CLEAN', `.note`. Each alternative's `source` rewritten by `cleanSource` (local): any substring matching a path or file name (`[\w./-]*\.mat`, `reproduce_\S+`, `sheet96_resolution_test\S*`, `departure_rib96_test\S*`, `round_\d+\S*`) is dropped; the kind words are kept; collapse repeated spaces; e.g. `'reproduce_70mN_24x48/round_01/arrival_sheet_70mN_nA48 cand, column 7 (re-certified)'` -> `'sheet candidate, column 7 (re-certified)'`; `'rib refusal at sA 0.2837 (walked 0)'` unchanged; `'displaced primary (merge)'` unchanged. Everything else bit-identical.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_export_deliverable_catalog()
% TEST_EXPORT_DELIVERABLE_CATALOG  The deliverable export keeps every datum
% (both content keys equal, schema valid), removes the internal records
% (merge, status_layer), rewrites alternative sources without file paths,
% and loads with a bare load.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));            % IND
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
recordMat = fullfile(here, 'results', 'library_70mN_24x48_v2', 'costate_catalog_dro_tulip_70mN.mat');
tmp = [tempname '.mat'];
[c, info] = export_deliverable_catalog(recordMat, tmp);
R = load(recordMat);  r = R.(char(fieldnames(R)));
ok = chk(ok, strcmp(catalog_content_key(c), catalog_content_key(r)) && strcmp(alternatives_content_key(c), alternatives_content_key(r)), 'both content keys equal the record''s');
ok = chk(ok, strcmp(info.primKey, catalog_content_key(r)) && strcmp(info.altKey, alternatives_content_key(r)), 'info carries the keys');
ok = chk(ok, isempty(catalog_schema('validate', c)), 'the export validates');
ok = chk(ok, ~isfield(c, 'merge') && ~isfield(c, 'status_layer') && isfield(c, 'build') && ischar(c.provenance), 'merge/status_layer removed; build + provenance text present');
ok = chk(ok, isequal(c.sheets, r.sheets) && isequal(c.status_key, r.status_key) && isequal(c.constants, r.constants), 'sheets, legend and constants bit-identical');
src = {c.alternatives.source};
ok = chk(ok, ~any(contains(src, '.mat')) && ~any(contains(src, 'reproduce_')) && ~any(contains(src, '/')), 'no file paths in any alternative source');
ok = chk(ok, any(strcmp(src, 'displaced primary (merge)')) || any(startsWith(src, 'rib refusal')) || any(startsWith(src, 'sheet candidate')), 'the kind words survive');
ok = chk(ok, c.build.nBorderline == nnz([c.alternatives.borderline]) && strcmp(c.build.auditVerdict, 'CLEAN'), 'build counts come from the data');
L = load(tmp);
ok = chk(ok, isfield(L, 'costate_catalog_dro_tulip_70mN') && isequal(L.costate_catalog_dro_tulip_70mN.sheets(1).z8, r.sheets(1).z8), 'saved under the catalog''s variable name; loads with load alone');
delete(tmp);
if ok, fprintf('test_export_deliverable_catalog: ALL PASS\n'); else, fprintf('test_export_deliverable_catalog: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
```

- [ ] **Step 2: Run it to see it fail** (MCP, pure): `Undefined function 'export_deliverable_catalog'`.

- [ ] **Step 3: Implement** `IND/export_deliverable_catalog.m` (house header). Body:

```matlab
if nargin < 3, opts = struct(); end
varName = fieldd(opts, 'varName', 'costate_catalog_dro_tulip_70mN');
R = load(recordMat);  c = R.(char(fieldnames(R)));
primKey = catalog_content_key(c);  altKey = alternatives_content_key(c);
removed = {};
sl = struct();  if isfield(c, 'status_layer'), sl = c.status_layer; end
mg = struct();  if isfield(c, 'merge'), mg = c.merge; end
for f = {'merge', 'status_layer'}
    if isfield(c, f{1}), c = rmfield(c, f{1});  removed{end+1} = f{1}; end
end
c.build = struct( ...
    'dates', '24 x 48 build 2026-09-20/21; refill 2026-10-02/03; merge with the 24 x 24 record 2026-10-03; status-layer backfill, re-certification and audit 2026-10-04/05', ...
    'nCandidates', fieldd(sl, 'nCandidates', NaN), 'nRecert', fieldd(sl, 'nRecert', NaN), 'nMoved', fieldd(sl, 'nMoved', NaN), ...
    'nBorderline', nnz(arrayfun(@(a) isfield(a, 'borderline') && isequal(a.borderline, true), c.alternatives)), ...
    'borderlineNote', fieldd(sl, 'borderlineNote', ''), 'mergeNote', fieldd(mg, 'note', ''), ...
    'auditVerdict', 'CLEAN', ...
    'note', 'Every alternative was re-certified under the shipped certifier and audited fail-closed (1,808 ok / 0 bad); see README.');
c.provenance = sprintf(['%s Status layer (FINDINGS 98): candidates harvested from the builds, every inferred one ' ...
    're-certified, primaries matched to their junctions by root, assembled beside the unchanged primaries, ' ...
    'audited fail-closed (CLEAN), borderline statuses relabelled to the meet of two runs. Exported for deliverable 8 on %s.'], ...
    c.provenance, char(datetime('now', 'Format', 'yyyy-MM-dd')));
n = 0;
for k = 1:numel(c.alternatives)
    s0 = c.alternatives(k).source;  s1 = cleanSource(s0);
    if ~strcmp(s0, s1), n = n + 1; end
    c.alternatives(k).source = s1;
end
assert(strcmp(catalog_content_key(c), primKey) && strcmp(alternatives_content_key(c), altKey), 'export_deliverable_catalog:keys', 'the export changed a content key');
p = catalog_schema('validate', c);
assert(isempty(p), 'export_deliverable_catalog:schema', 'the export does not validate: %s', strjoin(p, '; '));
S.(varName) = c;  save(outMat, '-struct', 'S');
info = struct('primKey', primKey, 'altKey', altKey, 'removed', {removed}, 'nSourceRewritten', n);
```
with local `cleanSource(s)`: `s = regexprep(s, '[\w./-]*\.mat', '')`; `s = regexprep(s, '(reproduce_|sheet96_resolution_test|departure_rib96_test|round_\d+)\S*', '')`; `s = regexprep(s, '^\s*cand\b', 'sheet candidate')`; `s = strtrim(regexprep(s, '\s+', ' '))`; `s = regexprep(s, '^\s*,\s*', '')` — then assert the result is non-empty (fall back to `'build candidate'`). And `fieldd`.

- [ ] **Step 4: Run to see it pass** (MCP). Expected ALL PASS. Also print `info.nSourceRewritten`.
- [ ] **Step 5: Mutation check** — (a) drop the `rmfield` loop; (b) `cleanSource` returns its input; (c) skip the key assert and perturb one `z8` entry. Each must fail the test.
- [ ] **Step 6: Commit (controller)** — `export_deliverable_catalog: the cleaned deliverable catalog, same content keys (spec 3)`.

---

### Task 2: `catalog_endpoints`

**Files:**
- Create: `IND/catalog_endpoints.m`
- Test: `IND/tests/test_catalog_endpoints.m`

**Interfaces:**
- Produces: `[B, orbits] = catalog_endpoints(cat, sheetIdx)`; `sheetIdx` [1]. `B.stateD(sD)`, `B.stateA(sA)` closures (phase FRACTION, wrapping; return the 6 x 1 state; the underlying `phase_state` returns the full row of the orbit table — take `(1:6)`), `B.dstateA(sA)`, `B.Tnd`, `B.cnd`, `B.mu`, `B.problem` = struct(`thrustN`, `ispS`, `m0kg`, `tauDRO`, `NpTulip`, `pmTulip`, `muStar`, `lStar`, `tStar`, `Tnd`, `cnd`, `periodTulip`, `sD` NaN) — the fields `plot_transfer_3d` and the study tools read; `orbits` = struct(`tD`, `rvD`, `tT`, `rvT`) the one-period samples.
- Rule: `[tD, rvD] = get_family_orbit(sh.dep_family, p)` with `p = sh.dep_params; p.muStar = cat.constants.muStar`; same for arrival; `[B.stateD, ~] = phase_state(tD, rvD, struct('onMissing', 'notaknot'))` and `[B.stateA, ~, B.dstateA] = phase_state(tT, rvT, ...)` — exactly `ladder_endpoints` + `arclength_arrival`'s rule. `Tnd = (cat.rungs_N(1)/cat.thruster.m0_kg)*tStar^2/(lStar_km*1000)` (the catalog's `derive.thrust_nd` formula), `cnd = cat.thruster.c_nd`.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_catalog_endpoints()
% TEST_CATALOG_ENDPOINTS  The shipped endpoint rule is the certified one: on
% the library of record its departure and arrival closures equal
% arclength_arrival's (via catalog_setup_request) bitwise at every grid
% phase and at the wrap, and Tnd / cnd agree to 1e-15.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
L = load(fullfile(here, 'results', 'library_70mN_24x48_v2', 'costate_catalog_dro_tulip_70mN.mat'));  c = L.(char(fieldnames(L)));
[B, orb] = catalog_endpoints(c);
[Bref, ~] = arclength_arrival('setup', catalog_setup_request(c, 0));
sh = c.sheets(1);
eqD = true;  eqA = true;
for s = sh.sD_frac, r = B.stateD(s);  q = Bref.stateD(s);  eqD = eqD && isequal(r(1:6), q(1:6)); end
for s = sh.sA_frac, r = B.stateA(s);  q = Bref.stateA(s);  eqA = eqA && isequal(r(1:6), q(1:6)); end
ok = chk(ok, eqD && eqA, 'departure and arrival states bitwise-equal to the certified rule at every grid phase');
r = B.stateD(1 - 1e-12);  q = B.stateD(0);
ok = chk(ok, sqrt(sum((r(1:6) - q(1:6)).^2)) < 1e-9, 'REVIEW FOCUS 2: the wrap s = 1 - 1e-12 gives the s = 0 state');
ok = chk(ok, abs(B.Tnd - Bref.Tnd) <= 1e-15*Bref.Tnd && abs(B.cnd - Bref.cnd) <= 1e-15*Bref.cnd && B.mu == Bref.mu, 'Tnd, cnd, mu agree');
ok = chk(ok, isequal(size(B.stateA(0.3)), [6 1]) && all(isfield(B.problem, {'lStar', 'tStar', 'muStar', 'Tnd', 'cnd', 'thrustN', 'NpTulip', 'tauDRO'})), 'closures return 6 x 1; B.problem has the consumers'' fields');
ok = chk(ok, size(orb.rvD, 2) >= 6 && orb.tD(1) == 0 && size(orb.rvT, 2) >= 6, 'orbit samples returned for plotting');
if ok, fprintf('test_catalog_endpoints: ALL PASS\n'); else, fprintf('test_catalog_endpoints: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
```

- [ ] **Step 2: Run to see it fail** (batch: `arclength_arrival` propagates orbits; takes ~1 min).
- [ ] **Step 3: Implement** per Interfaces. If `phase_state`'s closure returns a row longer than 6 (the orbit table may carry extra columns), wrap: `B.stateD = @(s) col6(xD(s))`. If the bitwise test fails only in the last bits, find the cause (e.g. `onMissing`, the muStar field, column order) — do NOT loosen the test; the certified rule must ship.
- [ ] **Step 4: Run to see it pass** (batch).
- [ ] **Step 5: Mutation check** — (a) omit `p.muStar`; (b) use `'spline'` interp1 instead of `phase_state`. Each must fail.
- [ ] **Step 6: Commit** — `catalog_endpoints: the certified endpoint rule from the catalog's own keys (spec 4)`.

---

### Task 3: `fly_catalog_entry`

**Files:**
- Create: `IND/fly_catalog_entry.m`
- Test: `IND/tests/test_fly_catalog_entry.m`

**Interfaces:**
- Produces: `F = fly_catalog_entry(cat, entry, opts)`.
  `entry`: a `transfers_at` row (has `.z8`, `.sD`, `.sA`, optional `.junctions .status .statusName .reason .kind`), OR `struct('iD', iD, 'iA', iA)` (the primary of that cell), OR `struct('z8', z8, 'sD', sD, 'sA', sA)`.
  `opts.method` ['auto'] -> `'junctions'` when `entry.junctions` is a non-empty 14 x K array, else `'z8'`; explicit `'junctions'` with no junctions -> error `fly_catalog_entry:noJunctions`; `opts.nSample` [50] per segment; `opts.B` [] (a `catalog_endpoints` result to reuse).
  Returns `F` with: `method`, `sD`, `sA`, `z` (= z8, for `plot_transfer_3d`), `z8`, `t`, `tDays`, `Y` (n x 14), `u` (n x 3), `massKg`, `tfNd`, `tfDays`, `dvKms`, `propellantKg`, `finalMassKg`, `flyKm`, `flyVms`, `defectKm`, `defectVms` (empty for z8), `nSamples`, `admissibility` (`validate_flight` of the flight: `.ok .reason .moonKm .earthKm`), `status`, `statusName`, `reason` (NaN/'not annotated'/'' when the entry has none), `rv0`, `rvf`, `B`.
  dvKms etc. from the rocket equation with `cat.thruster.c_nd` and the mass law (`mf = Y(end,7)`), as `fly_transfer` does.
- Consumes: `catalog_endpoints` (Task 2), `entry_thrust_program`, `fly_transfer`, `validate_flight`, `status_key`.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_fly_catalog_entry()
% TEST_FLY_CATALOG_ENTRY  A library entry flies from its junctions and from
% its z8 to the catalog's own arrival state (< 1 km / 0.1 m/s), the thrust
% direction is the one the equations of motion apply, the summary and the
% 3D figure accept the output, a cell index selects the primary, a row with
% no junctions falls back to z8 and says so, and 'junctions' without any
% refuses by name.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
L = load(fullfile(here, 'results', 'library_70mN_24x48_v2', 'costate_catalog_dro_tulip_70mN.mat'));  c = L.(char(fieldnames(L)));
sh = c.sheets(1);  iD = 3;  iA = 7;
T = transfers_at(c, sh.sD_frac(iD), sh.sA_frac(iA));  prim = T(strcmp({T.kind}, 'primary'));
Fj = fly_catalog_entry(c, prim);
Fz = fly_catalog_entry(c, prim, struct('method', 'z8'));
ok = chk(ok, strcmp(Fj.method, 'junctions') && Fj.flyKm < 1 && Fj.flyVms < 0.1 && max([Fj.defectKm 0]) < 0.05, sprintf('junction flight arrives (%.3f km, %.4f m/s), defects %.1e km', Fj.flyKm, Fj.flyVms, max([Fj.defectKm 0])));
ok = chk(ok, strcmp(Fz.method, 'z8') && Fz.flyKm < 1 && Fz.flyVms < 0.1 && isempty(Fz.defectKm), sprintf('z8 flight arrives (%.3f km, %.4f m/s)', Fz.flyKm, Fz.flyVms));
ok = chk(ok, abs(Fj.tfDays - prim.tfDays) < 1e-9 && abs(Fj.dvKms - prim.dvKms) < 1e-9 && Fj.status == 4 && strcmp(Fj.statusName, 'sufficient'), 'numbers and status carried from the row');
B = Fj.B;
for k = round(linspace(2, numel(Fj.t) - 1, 5))
    dyT = pumpkyn.cr3bp.tfMinEoM(Fj.t(k), Fj.Y(k, :).', B.Tnd, B.cnd, B.mu);
    dy0 = pumpkyn.cr3bp.tfMinEoM(Fj.t(k), Fj.Y(k, :).', 0, B.cnd, B.mu);
    a = dyT(4:6) - dy0(4:6);  a = a/sqrt(sum(a.^2));
    ok = chk(ok, sqrt(sum((a(:) - Fj.u(k, :).').^2)) < 1e-9, sprintf('thrust direction = the EoM''s applied acceleration (sample %d)', k));
end
lines = print_transfer_summary(Fj);
ok = chk(ok, iscellstr(lines) || ischar(lines) || isstring(lines), 'print_transfer_summary accepts F');
P = plot_transfer_3d(Fj, B, struct('visible', false, 'flight', struct('t', Fj.t, 'Y', Fj.Y)));
ok = chk(ok, ~isempty(P), 'plot_transfer_3d accepts (F, F.B)');  close all
Fc = fly_catalog_entry(c, struct('iD', iD, 'iA', iA));
ok = chk(ok, isequal(Fc.z8, prim.z8) && strcmp(Fc.method, 'junctions'), 'a cell index selects the primary (with its junctions)');
rowNoJ = prim;  rowNoJ.junctions = [];
Fn = fly_catalog_entry(c, rowNoJ);
ok = chk(ok, strcmp(Fn.method, 'z8') && contains(Fn.reasonMethod, 'no junctions'), 'REVIEW FOCUS 3: a row without junctions falls back to z8 and says so');
try, fly_catalog_entry(c, rowNoJ, struct('method', 'junctions'));  bad = false; catch ME, bad = strcmp(ME.identifier, 'fly_catalog_entry:noJunctions'); end
ok = chk(ok, bad, 'explicit ''junctions'' with none refuses by name');
Fw = fly_catalog_entry(c, struct('z8', prim.z8, 'sD', 1 - 1e-12, 'sA', prim.sA), struct('method', 'z8'));
ok = chk(ok, Fw.flyKm < 1, 'REVIEW FOCUS 2: a z8-with-phases entry at the wrap flies to the same arrival');
if ok, fprintf('test_fly_catalog_entry: ALL PASS\n'); else, fprintf('test_fly_catalog_entry: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
```
(`sD_frac(3)` = 0.0833 is not 0, so the wrap check uses `1 - 1e-12` as the departure phase of the (iD=1) cell instead: use `prim1 = transfers_at(c, 0, sh.sA_frac(iA))` for `Fw`'s z8 — adjust the test so the wrap case flies the sD = 0 primary.)

- [ ] **Step 2: Run to see it fail** (batch).
- [ ] **Step 3: Implement.** Resolve the entry form; `B = opts.B` or `catalog_endpoints(cat)`; `rv0 = B.stateD(sD)`, `rvf = B.stateA(sA)`; `phys = struct('Tnd', B.Tnd, 'cnd', B.cnd, 'mu', B.mu, 'muStar', B.mu, 'lStar', cat.constants.lStar_km, 'tStar', cat.constants.tStar_s, 'm0kg', cat.thruster.m0_kg)`; method `'junctions'` -> `P = entry_thrust_program(z8, J, rv0, rvf, phys, struct('nSample', nSample))`, copy `t Y u massKg defectKm defectVms flyKm flyVms`; `'z8'` -> `Ff = fly_transfer(z8, rv0, rvf, phys)`, `u = -Y(:,11:13)./|.|`, `massKg = Y(:,7)*m0`, defects `[]`; both: `admissibility = validate_flight(t, Y, z8(8), B.Tnd, B.cnd, B.mu, lStar)`, `tfNd`, `tfDays`, `mf = Y(end,7)`, `dvKms = cnd*log(1/mf)*lStar/tStar`, `propellantKg = (1-mf)*m0`, `finalMassKg`, `nSamples = numel(t)`, `z = z8`, `reasonMethod` ('auto: junctions' | 'auto: no junctions stored, flown from z8' | 'requested'), status fields from the row via `status_key` names (NaN -> 'not annotated').
- [ ] **Step 4: Run to see it pass** (batch).
- [ ] **Step 5: Mutation check** — (a) `u = +lv`; (b) `rvf` taken at `sD` instead of `sA`; (c) auto-method ignores empty junctions (calls entry_thrust_program with []). Each must fail.
- [ ] **Step 6: Commit** — `fly_catalog_entry: fly a library entry from its junctions or z8; summary/plot-ready output (spec 4)`.

---

### Task 4: `tfmin_accepts`

**Files:**
- Create: `IND/tfmin_accepts.m`
- Test: `IND/tests/test_tfmin_accepts.m`

**Interfaces:**
- Produces: `[ok, dz, z8out, info] = tfmin_accepts(cat, entry, opts)`; `entry` as in Task 3; `opts.tolDz` [1e-6], `opts.pool` [`current_pool()`], `opts.capSec` [300], `opts.B` []. Calls `pumpkyn.cr3bp.tfMin(rv0, rvf, z8, B.Tnd, B.cnd, B.mu)` under `run_capped` when a pool exists (else direct, with `info.fenced = false`); `dz = max(abs(z8out - z8))` (the certifier's gate quantity — check `certify_root`'s `C.dz` definition and use the same); `ok = isfinite(dz) && dz <= tolDz`; `info.expected` = true when the entry's status is 1 ("the certification witness failed on this entry; tfMin is not expected to accept it"); `info.message` one sentence; a cap timeout -> `ok = false`, `info.timedOut = true`, `info.message` names the cap; a thrown tfMin -> `ok = false`, `info.threw = ME.message`. NEVER throws for a status-1 or a timed-out entry.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_tfmin_accepts()
% TEST_TFMIN_ACCEPTS  pumpkyn tfMin accepts a certified primary unchanged
% (|dz| < 1e-6); a status-1 alternative is reported as an expected refusal,
% not an error; the dz definition is the certifier's.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
L = load(fullfile(here, 'results', 'library_70mN_24x48_v2', 'costate_catalog_dro_tulip_70mN.mat'));  c = L.(char(fieldnames(L)));
sh = c.sheets(1);
T = transfers_at(c, sh.sD_frac(2), sh.sA_frac(5));  prim = T(strcmp({T.kind}, 'primary'));
[acc, dz, z8o, info] = tfmin_accepts(c, prim, struct('pool', capped_pool()));
ok = chk(ok, acc && dz < 1e-6 && numel(z8o) == 8 && ~info.expected, sprintf('a certified primary is accepted unchanged (|dz| = %.1e)', dz));
A = c.alternatives;  k1 = find([A.status] == 1, 1);
row1 = transfers_at(c, A(k1).sD, A(k1).sA);  row1 = row1([row1.status] == 1);  row1 = row1(1);
[acc1, dz1, ~, info1] = tfmin_accepts(c, row1, struct('pool', capped_pool()));
ok = chk(ok, ~acc1 && info1.expected && ischar(info1.message) && ~isempty(info1.message), sprintf('REVIEW FOCUS 4: a status-1 row -> ok = false, expected = true, no error (|dz| = %.2g)', dz1));
[acc2, ~, ~, info2] = tfmin_accepts(c, prim, struct('pool', capped_pool(), 'capSec', 1e-3));
ok = chk(ok, ~acc2 && isfield(info2, 'timedOut') && info2.timedOut, 'REVIEW FOCUS 4: a fence timeout is a named refusal');
if ok, fprintf('test_tfmin_accepts: ALL PASS\n'); else, fprintf('test_tfmin_accepts: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
```
(If a 1e-3 s cap cannot be made to fire reliably with `run_capped`, use the seam `opts.tfminFcn` (a handle replacing tfMin) with a function that `pause(2)`s and a `capSec` of 0.5.)

- [ ] **Step 2: Run to see it fail** (batch, pool).
- [ ] **Step 3: Implement** per Interfaces; read `certify_root.m`'s witness block (`C.dz`) and copy the dz definition; `opts.tfminFcn` [`@pumpkyn.cr3bp.tfMin`] as the test seam.
- [ ] **Step 4: Run to see it pass** (batch; ~2-5 min).
- [ ] **Step 5: Mutation check** — (a) `ok = true` regardless of dz; (b) `expected` never set. Each must fail.
- [ ] **Step 6: Commit** — `tfmin_accepts: the 'tfMin accepts it unchanged' check, status-1 expected refusals, fenced (spec 4)`.

---

### Task 5: `costate_lib_describe` status-aware

**Files:**
- Modify: `IND/costate_lib_describe.m`
- Test: `IND/tests/test_describe_status.m`

**Interfaces:**
- Produces: same signature `D = costate_lib_describe(libIn)`. When `isfield(lib, 'status_key')` and the sheet has `status`: appends, after the existing printout, a block:
  ```
  optimality status (status_key):
    4 sufficient              primaries 1152   alternatives 268
    3 necessary only          primaries    0   alternatives 110
    2 conjugate point found   primaries    0   alternatives 268
    1 neither                 primaries    0   alternatives  10
  alternatives: 656 (5 borderline: k21, k219, k419, k491, k544); junctions on every entry (K = 24)
  ```
  and `D` gains `.statusCounts` ([4 x 2] primaries/alternatives by code order 4 3 2 1), `.nAlternatives`, `.nBorderline`, `.K`. Without the layer: the printout and `D` are byte-identical to today's (the test captures `evalc` on the 24 x 24 reference before and after... the "before" is captured by running the ORIGINAL file from git: `git show HEAD:orbit_transfer/DRO_tulip/indirect/costate_lib_describe.m > <scratch>/old/costate_lib_describe.m` and comparing `evalc` outputs).

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_describe_status()
% TEST_DESCRIBE_STATUS  costate_lib_describe prints the status block on a
% catalog with the layer (counts from the data) and is byte-identical to
% its previous self on a catalog without it.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
v2 = load(fullfile(here, 'results', 'library_70mN_24x48_v2', 'costate_catalog_dro_tulip_70mN.mat'));  c = v2.(char(fieldnames(v2)));
r24 = load(fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat'));  r = r24.(char(fieldnames(r24)));
out = evalc('D = costate_lib_describe(c);');
A = c.alternatives;
ok = chk(ok, contains(out, 'optimality status') && isequal(D.statusCounts(:, 2).', [nnz([A.status]==4) nnz([A.status]==3) nnz([A.status]==2) nnz([A.status]==1)]) && D.statusCounts(1, 1) == nnz(c.sheets(1).has_solution), 'status block printed; counts equal the data');
ok = chk(ok, D.nBorderline == nnz([A.borderline]) && D.K == size(c.sheets(1).junctions{1}, 2) && contains(out, 'borderline'), 'borderline count and K reported');
% byte-identity without the layer: compare with the committed version
old = fullfile(tempdir, 'describe_old');  if ~isfolder(old), mkdir(old); end
[st, txt] = system(sprintf('cd "%s" && /Library/Developer/CommandLineTools/usr/bin/git show HEAD:orbit_transfer/DRO_tulip/indirect/costate_lib_describe.m > "%s/costate_lib_describe.m"', fileparts(fileparts(fileparts(here))), old));
assert(st == 0, 'could not extract the committed file: %s', txt);
outNew = evalc('Dn = costate_lib_describe(r);');
cur = pwd;  cd(old);  outOld = evalc('Do = costate_lib_describe(r);');  cd(cur);
ok = chk(ok, strcmp(outNew, outOld) && isequal(rmfield_if(Dn, {'statusCounts','nAlternatives','nBorderline','K'}), Do), 'REVIEW FOCUS 1: without the layer the printout and D are unchanged');
if ok, fprintf('test_describe_status: ALL PASS\n'); else, fprintf('test_describe_status: FAIL\n'); end
end

function s = rmfield_if(s, f)
% RMFIELD_IF  Remove fields that exist.  INPUTS: s; f cellstr.  OUTPUTS: s.
for k = 1:numel(f), if isfield(s, f{k}), s = rmfield(s, f{k}); end, end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
```
(The old-version comparison runs the committed file from a scratch folder that shadows the new one via `cd`; `D` without the layer must not gain the new fields at all — make the new fields present only when the layer is.)

- [ ] **Step 2: Run to see it fail** (MCP). **Step 3: Implement** the block (local function `printStatusBlock(c, sh)`), only inside the `fmt == 'catalog'` branch and only when `isfield(c, 'status_key') && isfield(sh, 'status')`. Borderline row indices: `find([A.borderline])`. **Step 4: Run to pass.** **Step 5: Mutation** — (a) counts from `status` of primaries used for alternatives too; (b) block printed without the layer. **Step 6: Commit** — `costate_lib_describe: the status block (spec 5)`.

---

### Task 6: `costate_catalog_extremes` status filters

**Files:**
- Modify: `IND/costate_catalog_extremes.m`
- Test: `IND/tests/test_extremes_status.m`

**Interfaces:**
- `filt.includeAlternatives` [false], `filt.status` [] (vector of codes; applies to alternatives AND primaries when given). With alternatives on, candidates include every `cat.alternatives` row of the filtered sheet(s) (match by `.sheet` index); `pack` gains fields `kind` ('primary'|'alternative'), `status` (double, NaN when no layer), `statusName`, `reason` (''), `junctions` ([] for primaries unless the sheet has them — fill from `sh.junctions{entry_index}` when present), `altIndex` (0 for primaries). `report` prints `status: <name> -- <reason>` for alternatives (and for primaries when the layer is present). With `includeAlternatives` false and no `status`: the returned structs equal today's PLUS the new fields (test: compare on the old fields).

- [ ] **Step 1: Write the failing test** — on V2: defaults -> `eMin`/`eMax` old fields identical to the committed version's (same `git show` trick as Task 5, compare `rmfield` of the new fields); `includeAlternatives = true, metric 'tf'` -> the minimum-t_f transfer over primaries+alternatives has `tfDays <= ` the primaries-only minimum and `kind` set; `status = 2, includeAlternatives = true` -> both extremes have `status == 2` and non-empty `reason`; on the 24 x 24 reference with `includeAlternatives = true` -> no error, `kind = 'primary'`, `status` NaN.
- [ ] **Step 2-6:** RED (MCP), implement (the alternatives loop beside the primaries loop; `filtstr` prints the new filters), GREEN, mutations ((a) alternatives included by default; (b) status filter ignored), commit — `costate_catalog_extremes: alternatives and status filters, defaults unchanged (spec 5)`.

---

### Task 7: `plot_status_torus`

**Files:**
- Create: `IND/plot_status_torus.m`
- Test: `IND/tests/test_plot_status_torus.m`

**Interfaces:**
- Produces: `H = plot_status_torus(cat, opts)`; `opts.out` '' (PNG path), `opts.visible` [true], `opts.panels` [{'tf','nalt','conj'}]. Panels on the (sD, sA) grid: `'tf'` primaries' t_f in days (imagesc, grid lines at the phases); `'nalt'` alternatives per cell (0 where none; off-grid alternatives drawn as markers at their phases, counted separately in the title); `'conj'` cells holding a status-2 transfer (binary), plus markers for off-grid status-2 rows. Returns `H` with `.fig`, `.nAltPerCell` (nD x nA), `.conjCell` (nD x nA logical), `.nOffGrid`, `.nOffGridConj`. Works on a catalog without alternatives (panels 'nalt'/'conj' empty with a note). 1280 x 720 figure; house polish (fixed axes, titles with the counts).

- [ ] **Step 1: Write the failing test** — synthetic catalog: 2 x 3 grid, primaries everywhere, alternatives: two in cell (1,2) (status 2 and 3), one in (2,3) (status 2), one off-grid status 2 at (0.17, 0.3); assert `H.nAltPerCell == [0 2 0; 0 0 1]`, `H.conjCell == [false true false; false false true]`, `H.nOffGrid == 1`, `H.nOffGridConj == 1`; `opts.out` writes a PNG that exists and is > 10 kB; then a smoke run on V2 (`visible = false`) returns `sum(H.nAltPerCell(:)) + H.nOffGrid == 656`.
- [ ] **Step 2-6:** RED, implement, GREEN, mutations ((a) off-grid rows counted into the nearest cell; (b) conj panel uses status >= 2), commit — `plot_status_torus: where the transfers and the conjugate points sit (spec 5)`.

---

### Task 8: `costate_catalog_extremes_movies` through the fly layer + the sample pair

**Files:**
- Modify: `IND/costate_catalog_extremes_movies.m`
- Test: `IND/tests/test_extremes_movies_short.m`

**Interfaces:**
- Signature kept: `outStem = costate_catalog_extremes_movies(cat_, filt, metric, outStem, opts)` with NEW optional `opts`: `.pair` [] — a 2-element cell of `transfers_at` rows to film instead of the extremes (the deliverable's sample: a cell's primary and its fastest status-2 alternative), `.labels` [] (two chars), `.nFrames` [] (cap, for tests), `.fly` ['auto'] passed as `fly_catalog_entry`'s method. Trajectories come from `fly_catalog_entry` (so an alternative's junction flight is filmed) instead of the inline `tfMinProp`; the title shows each panel's status name when known. Frame size and the 1280 x 720 rule unchanged.
- Also a helper `[prim, alt2] = sample_status_pair(cat)` (local to the builder or `IND/sample_status_pair.m`): the cell with the most status-2 alternatives (ties: smallest t_f gap); returns its primary row and the fastest status-2 row.

- [ ] **Step 1: Write the failing test** — `opts.nFrames = 6`, `opts.pair = {prim, alt2}` from `sample_status_pair(c)`, output to a temp stem: both `.mp4` and `.gif` exist; the mp4 frame size read back with `VideoReader` is 1280 x 720; the returned stem matches. (The default extremes path on V2 is exercised by the builder in Task 10, not here, to keep the test under 2 min.)
- [ ] **Step 2-6:** RED (batch; SatelliteAnimator needs pumpkynPie; guard the Aero-toolbox `ShowEarth` trap as the existing code does), implement, GREEN, mutation ((a) the fly layer replaced by a flight from z8 for a junction row -> the status-2 alternative's miss grows; assert on `F.method` recorded in the returned struct), commit — `extremes movies through the fly layer; the status sample pair (spec 6)`.

---

### Task 9: the two example scripts

**Files:**
- Create: `DRO_tulip/indirect/examples/costate_catalog_example.m`, `DRO_tulip/indirect/examples/status_layer_example.m`
- Test: `IND/tests/test_examples_run.m` (runs both from a scratch folder holding ONLY the exported catalog + the deliverable file list, with pumpkyn on the path — a rehearsal of the recipient self-test)

**Rules for both:** study-script style (numbered `%%` sections, every number printed with its units, a short `fprintf` verdict per section); load by bare filename `catFile = 'costate_catalog_dro_tulip_70mN.mat'; L = load(catFile); cat_ = L.costate_catalog_dro_tulip_70mN;`; no absolute paths; `close all` at the end when run non-interactively (`if ~usejava('desktop')`).

**`costate_catalog_example.m`** sections: 1 load + `costate_lib_describe`; 2 choose five coordinates (tauDep 1, Np 7, dep 0.5 d, arr 3.0 d, 0.07 N) -> `costate_catalog_pick` (watch the honesty warnings); 3 derive units (days, Delta-V, propellant via `cat_.derive`); 4 `catalog_endpoints` + the endpoint states; 5 `fly_catalog_entry` (z8 with phases from `info.delivered`) + `print_transfer_summary`; 6 `tfmin_accepts` -> |dz|; 7 `plot_transfer_3d`.

**`status_layer_example.m`** sections: 1 load; 2 `sample_status_pair` -> the cell; 3 `transfers_at` at that cell, the ranked table (kind, status name, t_f d, Delta-V, reason first 60 chars); 4 fly the primary and the status-2 alternative from junctions (`fly_catalog_entry`), print both summaries and junction defects; 5 plot both thrust programs (|u| components and mass vs days, two panels); 6 the status-2 row's `reason` in full + what status 2 means (from `status_key`); 7 `tfmin_accepts` on the alternative (accepted: a conjugate-point transfer is still a converged extremal); 8 `plot_status_torus`.

- [ ] **Step 1: Write the failing test** — builds the scratch folder (export via Task 1 + copy the file list of Task 10's table — define the list ONCE in `IND/deliverable_file_list.m` returning a cellstr, created in THIS task and consumed by Task 10), `cd`s there in a `matlab -batch` child process with ONLY pumpkyn + pumpkynPie startup and that folder on the path (`restoredefaultpath` is forbidden by the MCP rule — use a child `-batch` process launched with `system`, whose script does `startup` for pumpkynPie, `addpath(scratch)`, then `run('costate_catalog_example.m')` and `run('status_layer_example.m')`), and asserts exit 0 and that the output contains `|dz|` and `ALL SECTIONS DONE` (each script prints that at the end).
- [ ] **Step 2-6:** RED, write the scripts + `deliverable_file_list.m`, GREEN (batch, ~5-10 min), commit — `deliverable examples: first run + the status layer; deliverable_file_list (spec 6)`.

---

### Task 10: `build_dro_deliverable` for deliverable 8 + the recipient self-test

**Files:**
- Modify: `DRO_tulip/build_dro_deliverable.m`
- Create: `DRO_tulip/deliverable_selftest.m` (repo-only; not shipped)
- Test: `DRO_tulip/tests/test_build_dro_deliverable.m` (create `DRO_tulip/tests/` if absent)

**Interfaces:**
- `out = build_dro_deliverable(cfg)`; `cfg.catMat` [V2], `cfg.name` ['costate_catalog_70mN_status_deliverable8'], `cfg.outDir` [`DRO_tulip/deliverables/<name>`], `cfg.movies` [true], `cfg.zip` [true], `cfg.selftest` [true], `cfg.auditPrimariesMat` [`MERGED/audit_70mN.mat`], `cfg.auditAltDir` [V2's folder]. Stages, each logged:
  1. **Export** via `export_deliverable_catalog` into `outDir`.
  2. **Audit gate by CONTENT KEY** (replaces the path equality): `A = load(auditPrimariesMat).A`; require `A.nBad == 0` and `A.contentKey == info.primKey`; every `audit_alt_*.mat` in `auditAltDir`: require `contentKey == info.primKey`, `altContentKey == info.altKey`, sum `audit.nBad == 0`, and the union of `[altRows.k]` == 1..nAlternatives exactly once (coverage). Refuse otherwise with a named error (`build_dro_deliverable:audit`).
  3. **Copy** every file in `deliverable_file_list()` (Task 9) from `IND`, `CC` or `IND/examples` (the list gives the source folder per file).
  4. **Movie:** `sample_status_pair` + `costate_catalog_extremes_movies(c, struct(), 'tf', fullfile(outDir, 'status_pair_movie'), struct('pair', {{prim, alt2}}, 'labels', ...))` when `cfg.movies`.
  5. **Torus PNG:** `plot_status_torus(c, struct('out', fullfile(outDir, 'status_torus.png'), 'visible', false))`.
  6. **README** from computed facts (extend `writeReadme` with the sections of spec 7; every number from `c`, the audit structs and the file list; the status table from the data; the five borderline rows listed with their reasons' first clause).
  7. **Self-test** (`cfg.selftest`): `deliverable_selftest(outDir)` — launches a child `matlab -batch` whose script does pumpkynPie `startup`, `addpath(outDir)`, runs both examples, calls every helper's bare self-demo that has one (`costate_lib_describe`, `costate_catalog_extremes` with no args if supported — check each file's `nargin == 0` branch and list the ones that exist), and runs the dependency scan: for every `.m` in `outDir`, `matlab.codetools.requiredFilesAndProducts` must return only files under `outDir`, `matlabroot`, or the pumpkyn/pumpkynPie roots (`fileparts(which('pumpkyn.cr3bp.tfMin'))`'s parent, pumpkynPie's folder, movie2gifLite's folder). Any other file -> the self-test FAILS naming it. Returns `.ok .log`.
  8. **Zip** only if the self-test passed: `zip -qr <name>.zip <name>` from the parent; verify with `unzip -l` (no `.DS_Store`).
- `out`: `.dir .zip .facts .files .selftest .keys`.

- [ ] **Step 1: Write the failing test** (`DRO_tulip/tests/test_build_dro_deliverable.m`): (a) with `cfg.movies = false, cfg.zip = false, cfg.selftest = false, cfg.outDir = <temp>`: builds; every file of `deliverable_file_list()` exists in `outDir`; README contains the strings `'656 alternatives'`, `'conjugate point found'`, `'borderline'`, `'fractions of a period'`, `'DAYS'`; the exported catalog in `outDir` has both keys equal to V2's. (b) REVIEW FOCUS 5: with `cfg.auditAltDir` pointing at a temp copy of one `audit_alt_*.mat` whose `altContentKey` is altered -> error id `build_dro_deliverable:audit`. (c) with a temp audit dir missing one chunk -> the same error (coverage).
- [ ] **Step 2: Run to see it fail** (batch; the export + README is seconds).
- [ ] **Step 3: Implement** the stages; keep the old path for deliverable 3 working behind `cfg.catMat` pointing at the old catalog (the old audit gate branch stays for catalogs without `alternatives`; the key gate applies when the catalog has them).
- [ ] **Step 4: Run to see it pass** (batch). **Step 5: Mutation** — (a) key gate compares `nBad` only; (b) coverage check removed. Each must fail. **Step 6: Commit** — `build_dro_deliverable: deliverable 8 (export, key-bound audit gate, file list, movie, torus, README, self-test, zip) (spec 7)`.

---

### Task 11: Build, self-test as the recipient, zip, docs (controller-run, with one docs subagent)

- [ ] **Step 1:** `build_dro_deliverable()` with defaults in batch (movies on, self-test on, zip on); keep the log. Expected: self-test `ok`, zip written, README numbers = the known ones (1,152 / 656 / 268-110-268-10 / 5 borderline / audit 1,808 ok / 0 bad).
- [ ] **Step 2:** Independent recipient self-test by the controller: fresh `matlab -batch`, pumpkynPie `startup`, `cd` into a COPY of the deliverable folder unpacked from the zip in a scratch location, run both examples, read the printed |dz| (< 1e-6) and flown misses (< 1 km).
- [ ] **Step 3:** The library-deliverable skill checklist: `unzip -l` reviewed, README sections present, no absolute paths (`grep -rn "/Users/" <folder>` empty).
- [ ] **Step 4 (docs subagent):** FINDINGS 99 (the tools, the deliverable, the self-test numbers); `DRO_tulip/README.md` (deliverable 8 row + the tool table); `orbit_transfer/STATUS_AND_ROADMAP.md` status line; memory note.
- [ ] **Step 5:** Commit the deliverable folder's SOURCE-controlled parts as the earlier deliverables were committed (check `git ls-files DRO_tulip/deliverables | head` to see whether zips/mats are tracked; follow the same rule), push, report to Mike: zip path + size, contents table, self-test printout, caveats to say verbally (interpolation not included; status-1 rows; phases in days vs fractions).

---

## Self-review notes

- Spec coverage: §3 -> Tasks 1, 10; §4 -> Tasks 2-4; §5 -> Tasks 5-7; §6 -> Tasks 8-9; §7 -> Tasks 10-11; §8 tests in every task; §9 honoured (no interpolation file in `deliverable_file_list`).
- Names used across tasks: `export_deliverable_catalog`, `catalog_endpoints`, `fly_catalog_entry` (+ `reasonMethod`), `tfmin_accepts`, `plot_status_torus`, `sample_status_pair`, `deliverable_file_list`, `deliverable_selftest`; row fields of `transfers_at` unchanged; extremes gains `kind status statusName reason junctions altIndex`.
- Review Focus 1-5 each pinned (Tasks 5/6, 2, 3, 4, 10).
