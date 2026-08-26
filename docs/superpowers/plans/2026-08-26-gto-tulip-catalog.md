# GTO→Tulip Min-Time Costate Catalog Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build, gate, and package the GTO→tulip min-time costate catalog (deliverable-7 candidate): 16 sheets (4 GTO orientations × 4 petal counts), 12×6×9 grid per sheet, every entry three-gate + conjugate certified.

**Architecture:** One new library piece (`'gto'` pseudo-family in `costate_common/get_family_orbit`, algebraic fixed-orientation locus) makes the entire existing catalog machinery (`thrust_ladder_library`, packager, pickers, `conj_catalog_pass`) work on this pair. New campaign module `GTO_tulip/catalog/` clones the halo-catalog driver pattern. Pilot sheet gates the fleet.

**Tech Stack:** MATLAB R2025b (`/Applications/MATLAB_R2025b.app/bin/matlab`), pumpkyn/pumpkynPie (via `~/Desktop/proj7/external/pumpkynPie` startup), `costate_common` pipeline, matlab-campaign batching discipline (nohup + script files + OS-kill backstop + Monitor).

**Spec:** `docs/superpowers/specs/2026-08-26-gto-tulip-catalog-design.md` — read it first; every decision below argues from it.

## Global Constraints

- Propulsion: **m0 = 150 kg, Isp = 1710 s, rungs [15 12 10 7 5 3 2 1.5 1] N** (spec §2).
- Sheets: orientDeg ∈ {0, 90, 180, 270} × Np ∈ {5, 7, 9, 12}; per-sheet grid **12 anomaly × 6 arrival**; pm = −1 only (spec §2).
- orientDeg convention: **angle from the Earth→Moon line to the ellipse's perigee direction, in the frame's rotation sense**; the flagship geometry is orientDeg = −25° (spec §3.1). `sD_frac` = TIME fraction from perigee (mean-anomaly fraction), because every engine interpolates `frac*tau(end)` on a time-parametrized locus.
- MATLAB house style: pumpkyn headers (`%% Purpose/Inputs/Outputs/Revision History`), never `i`/`j`, no `%#ok` pragmas, `check_matlab_code` clean on new files.
- Long runs: matlab-campaign discipline — nohup + script file, file logging, resume from sidecar/progress files, OS kill at budget+grace, Monitor with stall detection. Never `matlab -batch` inline `&`.
- **Do NOT `git commit`** — Mike commits. Each task's final step = report the exact file list + a suggested commit message and STOP.
- Prefer the `matlab` MCP tools (shared session) for short verifications; `matlab -batch` (nohup) for anything over ~10 min.
- Order correction vs spec §4: **packaging (Task 6) precedes the conjugate sweep (Task 7)** — `conj_catalog_pass` consumes a packaged catalog .mat.

---

### Task 1: `'gto'` pseudo-family in `get_family_orbit` (TDD)

**Files:**
- Modify: `orbit_transfer/costate_common/get_family_orbit.m` (add one `case 'gto'` + header FAMILIES line)
- Create: `orbit_transfer/costate_common/tests/test_gto_family.m`

**Interfaces:**
- Consumes: `pumpkyn.cr3bp.orb2eci`, `pumpkyn.cr3bp.fromPCI` (read their headers in `~/Desktop/proj7/external/pumpkyn/src/+pumpkyn/+cr3bp/` before coding — the plan mirrors `mintime_params.m`'s exact call pattern; the 6th orbital element's anomaly convention and the trailing flag MUST be taken from those headers, not guessed).
- Produces: `[tau, rv, info] = get_family_orbit('gto', p)` with `p.orientDeg` (required for sheets; default −25 = flagship), `p.sma_km` (default (6728+42164)/2), `p.ecc` (default (42164−6728)/(6728+42164)); `tau` = time since perigee over one Kepler period [M×1, ND]; `rv` = rotating-frame states [M×6]; `info.periodND` = Kepler period ND. Every later task keys on exactly these names.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_gto_family()
% TEST_GTO_FAMILY  Unit test for the 'gto' pseudo-family (Stage B spec 3.1).
% Three checks: flagship equality, rotation equivariance, geometry.
% INPUTS: none  OUTPUTS: ok [logical]
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));                          % costate_common
addpath(fullfile(fileparts(fileparts(fileparts(here))), ...  % repo root
        'orbit_transfer', 'GTO_tulip', 'indirect', 'min_time'));
ok = true;

% 1. FLAGSHIP EQUALITY: orientDeg=-25 at perigee reproduces mintime_params rv0
[rv0Ref, ~, P] = mintime_params();
[tau, rv, info] = get_family_orbit('gto', struct('orientDeg', -25));
ok = check('perigee state == flagship rv0 (1e-10)', ...
           max(abs(rv(1,:) - rv0Ref(:).')) < 1e-10) && ok;
ok = check('tau starts at 0', tau(1) == 0) && ok;
ok = check('periodND positive and ~ Kepler', ...
           info.periodND > 0 && abs(tau(end) - info.periodND) < 1e-9) && ok;

% 2. ROTATION EQUIVARIANCE: rotating the orientation rotates the whole locus
%    (positions about the barycenter-frame z axis THROUGH EARTH at (-mu,0,0):
%    the ellipse is Earth-centered, so rotate the Earth-relative vector)
th = 37 * pi/180;
Rz = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
[~, rvA] = get_family_orbit('gto', struct('orientDeg', 10));
[~, rvB] = get_family_orbit('gto', struct('orientDeg', 10 + 37));
rE = [-P.muStar 0 0];
posErr = 0; velErr = 0;
for k = [1, 7, 13]                                  % spot rows
    pA = (Rz*(rvA(k,1:3) - rE).').' + rE;
    vA = (Rz*rvA(k,4:6).').';
    posErr = max(posErr, max(abs(pA - rvB(k,1:3))));
    velErr = max(velErr, max(abs(vA - rvB(k,4:6))));
end
ok = check(sprintf('rotation equivariance pos %.1e vel %.1e < 1e-10', ...
           posErr, velErr), posErr < 1e-10 && velErr < 1e-10) && ok;

% 3. GEOMETRY: radii from Earth at frac 0 / 0.5 are perigee / apogee
rp_km = (6378+350); ra_km = (6378+35786);
rPer = norm(rv(1,1:3) - rE) * P.lStar;
kAp = round((numel(tau)+1)/2);                      % M = pi at half period
rAp  = norm(rv(kAp,1:3) - rE) * P.lStar;
ok = check(sprintf('perigee radius %.1f ~ %d km', rPer, rp_km), ...
           abs(rPer - rp_km) < 1.0) && ok;
ok = check(sprintf('apogee radius %.1f ~ %d km', rAp, ra_km), ...
           abs(rAp - ra_km) < 1.0) && ok;
fprintf('test_gto_family: %s\n', string(ok));
end

function ok = check(name, cond)
% CHECK  One gate line.  INPUTS: name;cond  OUTPUTS: ok
ok = logical(cond);
if ok, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('  [%s] %s\n', tag, name);
end
```

Note for the implementer: the odd-count sampling makes row `(M+1)/2` land exactly on apogee — the reference implementation uses M=2001 uniform in ECCENTRIC anomaly (review fix; NOT 25, NOT uniform time). Additionally the test MUST include a NON-APSIDAL point (Gemini, critical — apsidal points mask a mean-vs-true-anomaly confusion): at mean anomaly π/2, compute the expected rotating-frame state INDEPENDENTLY inside the test (own 5-line Newton Kepler solve E−e·sinE=M → true anomaly → ellipse state in inertial → rotate by orientDeg → subtract ω×r for the rotating velocity) and require agreement < 1e-9. Also assert `tau` is strictly increasing and document (in the case comment) that tau is a locus parameter, never an epoch offset.

- [ ] **Step 2: Run test to verify it fails**

Via MCP `evaluate_matlab_code`, cd `costate_common/tests`: `ok = test_gto_family()`.
Expected: error `unknown family 'gto'` (or the case's equivalent from `get_family_orbit`).

- [ ] **Step 3: Read the two pumpkyn headers, then implement the case**

First read `orb2eci.m` and `fromPCI.m` headers (anomaly convention of element 6, meaning of the trailing flags) and `GTO_tulip/indirect/min_time/mintime_params.m:23-26` (the exact flagship call). Then add to `get_family_orbit`'s switch (adapt if the anomaly flag differs — the FLAGSHIP-EQUALITY test is the arbiter):

```matlab
case 'gto'
    % GTO pseudo-family (Stage B spec 2026-08-26): the locus of departure
    % states at FIXED orientation, parametrized by TIME fraction (mean
    % anomaly) over one Kepler period. Algebraic construction -- no
    % propagation; orientation is a sheet key, not a flow consequence.
    % orientDeg: Earth->Moon line to perigee direction, rotation sense.
    smaKm  = fieldd(p, 'sma_km', (6378+350 + 6378+35786)/2);
    eccGto = fieldd(p, 'ecc', (35786-350)/(2*((6378+350+6378+35786)/2)));
    orient = fieldd(p, 'orientDeg', -25);
    consts = local_constants();                 % muStar, lStar, tStar (one home)
    muE = 6.67384e-20*(1 - consts.muStar)*(5.9736E24 + 7.35E22);
    % REVIEW FIX (GPT+Gemini, critical): sample uniform in ECCENTRIC anomaly
    % (dense) so perigee is resolved on an e~0.72 ellipse; tau from Kepler.
    M  = 2001;                                  % odd: sample lands on apogee
    E_ = linspace(0, 2*pi, M);
    mAnom = E_ - eccGto*sin(E_);                % Kepler: M(E), monotone
    Tkep  = 2*pi*sqrt((smaKm)^3/muE);           % s
    rv  = zeros(M, 6);
    for km = 1:M
        % element-6 anomaly per orb2eci's own convention (checked in step 3
        % preamble); the flagship test pins correctness:
        [rr, vv] = pumpkyn.cr3bp.orb2eci(muE, ...
                     [smaKm, eccGto, 0, orient*pi/180, 0, mAnom(km)], 2);
        rv(km,:) = pumpkyn.cr3bp.fromPCI(0, [rr, vv], consts.muStar, ...
                     consts.tStar, consts.lStar, 1);
    end
    tau  = (mAnom.'/(2*pi)) * (Tkep/consts.tStar);
    info = struct('periodND', Tkep/consts.tStar, 'family', 'gto', ...
                  'sma_km', smaKm, 'ecc', eccGto, 'orientDeg', orient);
```

`local_constants` / `fieldd`: use whatever helpers `get_family_orbit` already has for the other families (it already resolves muStar/lStar/tStar for dro/tulip/halo — reuse that mechanism verbatim; if constants are inline per-case, inline the same three values `0.012150585609624 / 389703.264829278 / 382981.289129055`). If `orb2eci`'s element 6 is TRUE anomaly (flag-dependent), insert a Kepler solve `E - e*sin(E) = M`, `tan(nu/2) = sqrt((1+e)/(1-e))*tan(E/2)` before the call — 5 lines, Newton on E, tol 1e-14.

- [ ] **Step 4: Run test to verify it passes**

`ok = test_gto_family()` → all 6 checks PASS.

- [ ] **Step 5: Lint + regression**

`check_matlab_code` on `get_family_orbit.m` (no new warnings beyond pre-existing style notes); run `golden_cells()` → 20/20 (no engine touched, but this is the standing gate for any `costate_common` edit).

- [ ] **Step 6: Report for commit**

Report files (`get_family_orbit.m`, `tests/test_gto_family.m`) + suggested message: `costate_common: 'gto' pseudo-family -- fixed-orientation GTO departure locus (Stage B task 1)`. Do not commit.

---

### Task 2: Campaign module scaffold + endpoint/solver smoke (1×1 cell)

**Files:**
- Create: `orbit_transfer/GTO_tulip/catalog/setup_paths.m`
- Create: `orbit_transfer/GTO_tulip/catalog/smoke_gto_cell.m`
- Create: `orbit_transfer/GTO_tulip/catalog/results/` (dir; `.mat` gitignored already)

**Interfaces:**
- Consumes: `get_family_orbit('gto', ...)` from Task 1; `DRO_tulip/indirect/thrust_ladder_library.m` — READ ITS HEADER FIRST for the exact opts/meta contract (grid fields seen in `HALO_tulip/run_halo_catalog.m:103`: `'nD','nA','sD0','sA0',...`; mirror `run_halo_catalog`'s invocation precisely).
- Produces: a proven 1×1×[15 N] solve through the production engine — the cold-start path exercised exactly as the campaign will use it. Output `results/smoke_gto_cell.mat` with the engine's own OK/Z8/TF arrays.

- [ ] **Step 1: Write `setup_paths.m`**

```matlab
function setup_paths()
% SETUP_PATHS  GTO_tulip/catalog module paths: this dir + costate_common +
% DRO_tulip/indirect (ladder engine) + pumpkyn (via pumpkynPie startup if
% absent). INPUTS: none  OUTPUTS: none
here = fileparts(mfilename('fullpath'));
ot = fileparts(fileparts(here));
addpath(here, fullfile(ot, 'costate_common'), ...
        fullfile(ot, 'DRO_tulip', 'indirect'));
if isempty(which('pumpkyn.cr3bp.tfMin'))
    pp = '/Users/msc/Desktop/proj7/external/pumpkynPie';
    od = cd(pp); startup(); cd(od);
end
end
```

- [ ] **Step 2: Write `smoke_gto_cell.m`**

Read `HALO_tulip/run_halo_catalog.m` end-to-end first (it is <150 lines); then write the smoke as a 1×1-grid, single-rung [15] invocation of `thrust_ladder_library` with `depFamily='gto'`, `depParams=struct('orientDeg',0)`, `arrFamily='tulip'`, `arrParams=struct('Np',7,'pm',-1)`, `sD0` at apogee (0.5), `sA0` 0, m0 150, isp 1710, logFile `results/smoke.log`, outMat `results/smoke_gto_cell.mat` — copying the halo driver's exact option names. (The engine's own header is authoritative for field names; the halo driver is the working example.)

- [ ] **Step 3: Run the smoke via MCP** (minutes at 15 N)

Expected: 1/1 pair solved; engine log shows the three gates (ms residual ≤1e-10, flown <100 km, tfMin accept |Δz|<1e-6) passing; `nnz(OK)==1` in the saved .mat (census from the DATA file).

- [ ] **Step 4: If the cold 15 N start fails** — do NOT tune blindly: capture the log, check `preflight_screen` output, and try `sD0` at perigee (0.0) — one alternative, then STOP and report findings (this is the spec's early-risk surface). A smoke failure here is adjudication material, not a debugging spiral.

- [ ] **Step 5: Report for commit**

Files + message: `GTO_tulip/catalog: module scaffold + 1x1 smoke through the ladder engine (Stage B task 2)`.

---

### Task 3: `run_gto_catalog.m` driver + pilot sheet (orient 0°, Np 7)

**Files:**
- Create: `orbit_transfer/GTO_tulip/catalog/run_gto_catalog.m` (clone of `HALO_tulip/run_halo_catalog.m` with the sheet block below) — THE SWATH DRIVER: `run_gto_catalog()` = all 16 sheets; `run_gto_catalog(sheetSel)` = subset; resume/regeneration of unsolved cells is automatic (the engine revisits no-OK cells)
- Create: `orbit_transfer/GTO_tulip/catalog/gto_entry.m` — THE SINGLE-ENTRY DRIVER (Mike's product request 2026-08-26): pumpkyn style with `nargin==0` self-demo; signature `[entry, gates] = gto_entry(orientDeg, Np, depFrac, arrFrac, thrustN, opts)`; runs the SAME ladder engine on a 1×1 phasing grid with rungs = the standard ladder truncated at `thrustN` (identical solve path + three gates as campaign entries — never a bespoke solver); opts: `.conjTest` [false] (runs ms_tfmin conjugate test on the accepted entry), `.writeSheet` [false] (merge the entry into `results/catalog/gto_o<...>_Np<...>.mat`, creating it if absent). Returns `entry` (z8, tf_nd, dV_kms, mf, sheet/grid coordinates) + `gates` (msNormR, flownKm, acceptDz, conjPass when run). Self-demo: orient 0, Np 7, apogee departure, 15 N.
- Create: `orbit_transfer/GTO_tulip/catalog/run_gto_batched.sh` (clone of `HALO_tulip/run_halo_batched.sh`, names swapped)

**Interfaces:**
- Consumes: Task 2's proven invocation pattern.
- Produces: `run_gto_catalog(sheetSel)` — solves the selected sheets into `results/catalog/gto_o<ORIENT>_Np<NP>.mat` (one file per sheet, engine-native format: `OK/TF/Z8/ATT/rungs/sD/sA/meta`), resumable; `sheetSel` = index vector into the 16-sheet list (default 1:16). Sheet naming used by ALL later tasks: `gto_o000_Np5 … gto_o270_Np12`.

- [ ] **Step 1: Write the driver.** Copy `run_halo_catalog.m`; replace its sheet-list construction with:

```matlab
orients = [0 90 180 270];                    % deg -- sheet axis 1 (spec 2)
Nps     = [5 7 9 12];                        % sheet axis 2
sheets = {};
for od = orients
    for np = Nps
        sheets{end+1} = struct( ...
            'depFamily', 'gto', ...
            'depParams', struct('orientDeg', od), ...
            'arrFamily', 'tulip', ...
            'arrParams', struct('Np', np, 'pm', -1), ...
            'tag', sprintf('gto_o%03d_Np%d', od, np));
    end
end
```

and set the grid options: `nD = 12`, `nA = 6`, rungs `[15 12 10 7 5 3 2 1.5 1]`, m0 150 kg, isp 1710 s — otherwise field-for-field identical to the halo driver (progress files, resume, per-sheet outMat under `results/catalog/`). THREE REVIEW FIXES: (a) **the sheet meta's `tauDRO` field = orientDeg** (the sheet KEY — all three reviewers: the picker selects by nearest (tauDRO, Np) and the Kepler period is orientation-degenerate); the period rides in `depParams`. (b) **`sD` grid non-uniform**: if the engine accepts an explicit fraction vector (check its header; `Q.sD` is stored per sheet so it likely does), pass 12 values uniform in TRUE anomaly mapped through Kepler to time fractions: `nuG = (0:11)/12*2*pi; EG = 2*atan(sqrt((1-e)/(1+e))*tan(nuG/2)); sD = mod(EG - e*sin(EG), 2*pi)/(2*pi)` — else uniform with the caveat recorded in FINDINGS.md. (c) assert the wrapped fractions are unique (no sD=1 duplicate of 0).

- [ ] **Step 1b: Write `gto_entry.m`** per the Files description; its self-demo IS its test — run it (15 N apogee entry converges, three gates pass) and additionally call it once with `opts.conjTest=true`. Verify the returned entry's z8 round-trips through `pumpkyn.cr3bp.tfMin` unchanged (<1e-6).

- [ ] **Step 2: Lint** all new files with `check_matlab_code`.

- [ ] **Step 3: Launch the pilot** (sheet = orient 0 × Np 7 → index it from the list) via the batched shell script under nohup, with a Monitor (pass-boundary lines + ERROR lines + stall detection at 900 s log-age, per matlab-campaign). Pilot size: 12×6×9 = 648 attempts; measured catalog rates suggest hours-scale.

- [ ] **Step 4: Pilot census from the DATA file** (not the log):

```matlab
Q = load('results/catalog/gto_o000_Np7.mat');
fprintf('pairs solved %d/72, entries %d/648, ladder-complete pairs %d\n', ...
    nnz(any(Q.OK,3)), nnz(Q.OK), nnz(all(Q.OK,3)));
```

Census must ALSO report (review fixes): per-rung coverage (especially the 1 N rung), ladder-complete-pair count, minimum EARTH altitude over all flown entries (report-only Earth-clearance burn-in; Earth radius 6378 km at (−μ,0,0)), and a **deliberate kill-and-resume check**: kill one batch mid-sheet, resume, verify no lost/duplicated entries and monotone ATT counters from the DATA file.

**ADJUDICATION CHECKPOINT (spec §4.3): report the census + per-cell cost to Mike and STOP if pair solvability < 70%, OR 1 N-rung coverage < 50%, OR cost extrapolates the fleet beyond ~a week, OR any entry dips below 200 km Earth altitude.** Otherwise proceed.

- [ ] **Step 5: Report for commit** (driver + shell + pilot progress txt; .mat stays untracked): `GTO_tulip/catalog: driver + pilot sheet o000/Np7 (Stage B task 3) -- <X>% pairs, <Y> entries`.

---

### Task 4: Full campaign (remaining 15 sheets)

**Files:**
- Modify: none (operational task; uses Task 3's driver)
- Create: `results/catalog/gto_*_progress.txt` per sheet (engine-written)

**Interfaces:**
- Consumes: `run_gto_catalog(2:16)` (or the batched script's sheet loop).
- Produces: 16 complete sheet .mats — the packager's input.

- [ ] **Step 1: Launch** remaining sheets via `run_gto_batched.sh` (nohup, per-pass OS kill, resume). Expect days-scale; the batch driver loops passes until a no-progress pass.
- [ ] **Step 2: Monitor** (one persistent Monitor: pass boundaries, `CONJUGATE\|ERROR\|STALL`-class lines, stall detection). Record findings as they land (hard corners, wall patterns) in `GTO_tulip/catalog/FINDINGS.md` — created on first finding, not before.
- [ ] **Step 3: Densify pass** — rerun `run_gto_catalog(1:16)` once after the first full sweep (the engine revisits no-OK cells with fresh attempts; stop when a pass adds no entries).
- [ ] **Step 4: Full census from data files** (all 16 sheets): per-sheet pairs/entries table; compare against the ~85–92% expectation; identify the hard corners (expect shortest-departure × longest-tulip per every prior catalog).
- [ ] **Step 5: Report for commit** (progress txts + FINDINGS.md): `GTO_tulip/catalog: full 16-sheet campaign -- <N> entries, <pct>% pairs (Stage B task 4)`.

---

### Task 5: Spot-audit replay (3 entries)

**Files:**
- Create: `orbit_transfer/GTO_tulip/catalog/audit_gto_entries.m`

**Interfaces:**
- Consumes: sheet .mats; `get_family_orbit`; `seed_from_z8`; `pumpkyn.cr3bp.tfMin`.
- Produces: acceptance-criteria evidence (spec §5): ≥3 entries across different sheets replayed end-to-end — fly the stored z8, measure arrival, hand to tfMin, require unchanged.

- [ ] **Step 1: Write the audit** (pattern = the deliverable example scripts): pick entries PROGRAMMATICALLY (GPT review — hard-coded coordinates may be unsolved): from the packaged catalog select 3 solved entries spanning ≥2 different ORIENTATIONS and ≥2 rung classes (one low-rung), rebuild endpoints via `get_family_orbit` FROM THE SHEET'S OWN dep_params (this also proves sheet-local reconstruction — the packager must not homogenize orientations), fly via `pumpkyn.cr3bp.tfMinProp`, assert flown miss < 100 km, `norm(tfMin(...)-z8) < 1e-6`, and report each entry's minimum Earth altitude; print the table.
- [ ] **Step 2: Run it**; all three rows pass.
- [ ] **Step 3: Report for commit**: `GTO_tulip/catalog: 3-entry replay audit (Stage B task 5)`.

---

### Task 6: Package + validate + describe

**Files:**
- Create: `orbit_transfer/GTO_tulip/catalog/build_gto_catalog.m` (thin: one call)
- Output: `orbit_transfer/GTO_tulip/catalog/results/costate_catalog_gto_tulip.mat`

**Interfaces:**
- Consumes: `build_costate_catalog_family(catDir, outMat, spec)` (read its header; the halo campaign's builder call is the working example — find it via `grep -rn build_costate_catalog_family HALO_tulip/`).
- Produces: the schema-v2 catalog, single variable `costate_catalog_gto_tulip`; `dep_family='gto'`, `dep_params` = {sma_km, ecc, orientDeg} PER SHEET (sheet-local — audit proves it); **`tauDRO`/`tau_dep` sheet key = orientDeg** (review fix; period in dep_params/derive; README documents the override). This is the file Tasks 7–8 operate on.

- [ ] **Step 1: Write + run the builder** with `spec.glob='gto_o*_Np*.mat'`, `spec.name='costate_catalog_gto_tulip'`, description/provenance strings naming the orientation-axis convention, `spec.depReconstruction` = the one-line derive string for rebuilding the GTO locus from `dep_params`.
- [ ] **Step 2: Validate**: `catalog_schema('validate', cat)` → `{}` (empty problem list). If the validator objects to the non-tulip-style dep keys, the fix goes in the BUILDER call (v2 fields), not the validator — the validator only changes if a rule is genuinely wrong, reported first.
- [ ] **Step 3: Describe**: `costate_lib_describe('costate_catalog_gto_tulip.mat')` — sane grids/ranges/coverage; capture the output for the README.
- [ ] **Step 4: Report for commit**: `GTO_tulip/catalog: packaged schema-v2 catalog (<N> entries) (Stage B task 6)`.

---

### Task 7: Conjugate sweep + writeback

**Files:**
- Modify: none expected (`conj_catalog_pass` reads v2 recipes; `'gto'` flows through `get_family_orbit`)
- Output: verdicts in `costate_catalog_gto_tulip.mat` (+ `.bak_conj`, `_conjprog.mat` sidecar)

**Interfaces:**
- Consumes: `conj_catalog_pass(catMat, opts)` — the 2026-08-23 sweep engine (sidecar resume, honesty gate |z−z8|<1e-6, writeback on complete census).
- Produces: `conj_pass/conj_ncross/conj_atfinal` per sheet + `conj_test` provenance — the standing-rule requirement before any deliverable.

- [ ] **Step 1: Smoke 10 entries** (`maxEntries=10`, scratch sidecar): verdicts recorded, DZ ~1e-12..1e-10. If `'gto'` fails to flow (e.g. the fallback path triggers), fix the RECIPE (builder), not the sweep — report first if unclear.
- [ ] **Step 2: Full sweep** — batched nohup passes exactly like 2026-08-23 (900 s budget per pass, driver loop, Monitor with CONJUGATE POINT/UNVERIFIED/STALL lines); ~N×0.2 s ⇒ under an hour for ~10k entries.
- [ ] **Step 3: Writeback on complete census** (`writeback=true`); re-validate schema (conj fields checked); census printout pass/fail/notrun. Any FAIL entries: record coordinates in FINDINGS.md + README (the DPO precedent).
- [ ] **Step 4: Report for commit**: `GTO_tulip/catalog: conjugate verdicts stored -- <pass>/<fail>/<notrun> (Stage B task 7)`.

---

### Task 8: README + records wave

**Files:**
- Create: `orbit_transfer/GTO_tulip/catalog/README.md`
- Modify: `orbit_transfer/STATUS_AND_ROADMAP.md` (Stage B row → DONE with numbers), `orbit_transfer/TODO.md` (Stage B section), `orbit_transfer/GTO_tulip/README.md` (catalog product line added to the folder map), `orbit_transfer/OPTIMALITY_CERTIFICATION.md` §6 (campaign + sweep entries), `orbit_transfer/doc/transfer_problem_space.md` (GTO row), memory (new file + MEMORY.md line)

**Interfaces:**
- Consumes: every census number from Tasks 4–7.
- Produces: the record trail; the catalog README is the deliverable-facing document.

- [ ] **Step 1: `catalog/README.md`** — what it covers (sheets/grids/rungs/ranges), the **orientation-axis convention** (spec §3.1 wording, for Darin), quick start: the TWO-DRIVER product (`gto_entry` for one entry, `run_gto_catalog` for swaths/regeneration) + the picker call for lookup, coverage + hard corners, conjugate census, and the **identifiability rule** for the future mN extension (spec §2, verbatim: many-rev entries must ship full ms junction states, not bare z8). Note the promotion path: a generic `costate_common/catalog_entry` when a second campaign wants the single-entry driver (migration rule).
- [ ] **Step 2: The records wave** — each file listed above, dated entries in the established formats (register discipline: append §6, never rewrite history).
- [ ] **Step 3: Memory** — `gto-catalog-complete.md` (type: project): totals, coverage, hard corners, per-cell cost, anything that surprised; index line in MEMORY.md.
- [ ] **Step 4: Report for commit** (full session file list, grouped into the natural commits): `GTO_tulip/catalog: README + records (Stage B task 8)`.

---

## Self-review notes (done at write time)

- **Spec coverage:** §1→Tasks 3–8; §2 constants→Global Constraints; §3.1→Task 1; §3.2→Tasks 2–4; §3.3→Tasks 6–7; §4 sequencing→Tasks 1–8 (order fix: package before sweep, noted in Global Constraints); §5 acceptance→Tasks 1 (test+goldens), 3 (pilot gate), 5 (replay audit), 6 (validate/describe), 7 (conjugate 100%), 8 (records); §6 risks→Task 2 step 4, Task 3 checkpoint, Task 4 FINDINGS; §7 out-of-scope respected (no deliverable zip, no mN, no dual capture).
- **Known unknowns made explicit, not placeholder'd:** `orb2eci` anomaly-flag semantics (Task 1 step 3 preamble: read the header; the flagship-equality test is the arbiter) and `thrust_ladder_library`'s exact opts names (Task 2: the halo driver is the working example to mirror). Both are read-then-mirror steps against named working code, with tests pinning correctness.
- **Type consistency:** sheet tags `gto_o%03d_Np%d` used identically in Tasks 3/4/6; catalog variable `costate_catalog_gto_tulip` in 6/7/8; `p.orientDeg/sma_km/ecc` per Task 1's contract everywhere.
