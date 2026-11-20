function ok = test_build_dro_elfo_catalog()
%% Purpose:
%
%   Test of the chain build_dro_elfo_catalog (Task 7) on a REAL 2 x 2,
%   two-rung fixture sheet (cut once by the test, into a tempname folder,
%   from the matched-experiment sheet: departure rows [2 5], arrival
%   columns [2 5], rungs [15 12] N), in a sheet folder that also holds
%   decoys a glob would catch (<tag>_conj.mat, <tag>_direct.mat-style
%   names, a tulip control sheet). Real conjugate and gate sweeps run on
%   the 8 entries: a few minutes, run it in the background.
%     1  every stage runs (census, clearance, both repair modes with
%        maxRepair 0, package, conj, gates, audit, pictures, preflight; each
%        stage's log line is present); the catalog validates with
%        catalog_schema; exactly one packaged sheet (no decoy), carrying the
%        sheet's sheet_key; the audit replays 3 of 3 through dro_elfo_pick;
%        pictures written, no figure left open, DefaultFigureVisible
%        restored; lock released; package_stage removed; verdict file
%     2  RESUME: the same call again keeps the existing catalog (content
%        unchanged), flies no clearance entry and re-tests no verdict
%     3  an INJECTED REFUTED ENTRY (conj sidecar verdict 0 on a passing
%        entry) produces the 'conj:' blocker; allowConjBlockers waives it
%        from the verdict only (it stays in the returned list)
%     4  an INJECTED STALE CLEARANCE SIDECAR (z8 record of one OK entry
%        changed) produces 'stale:' and is NOT waivable (verdict false with
%        both waivers on); the clearance stage re-flies the entry and the
%        blocker clears
%     5  a changed sheet with the package stage off: 'stale:' (package)
%     5b DELIBERATE FAILURES of the other blocking gates: an OK entry whose
%        sidecar record is below the floor -> 'unverified:' (blocking with
%        every waiver on); a gate failure injected in the gates sidecar ->
%        'gates:' (not waivable); a sheet carrying the wrong sheet_key ->
%        'census:'; the packaged z8 of an audited entry off by 1e-9 ->
%        'audit:'; the folder is clean again after each is undone
%     5c entries below the Moon floor (floor raised to 8000 km on a copy):
%        the clearance stage and, with it off, the preflight report
%        'clearance-demoted:'; verdict false without allowClearanceDemotions,
%        true with it; the verdict file lists the demotion both times
%     6  refusals: a foreign fleet lock (:fleetRunning, lock survives), an
%        unknown runIn field (:runIn), an unknown stage (:stages); an empty
%        sheet folder returns a census blocker without throwing
%     7  a repair that changes a sheet (real engine: a 6 x 1 x 2 cut with
%        one cell emptied and retired) -> 2b re-sweeps the re-solved
%        entries and the package holds the recovered entry
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 true iff every check
%                                                   passed (throws on a
%                                                   failure)
%
%% Revision History:
%  M. Casey                                                   (c) 11/19/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here    = fileparts(mfilename('fullpath'));
oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(fileparts(here));
setup_dro_elfo_paths();
addpath(here, fullfile(fileparts(here), 'catalog'));
wd = tempname;  mkdir(wd);
cleanWd = onCleanup(@() rmdir(wd, 's'));
catDir = fullfile(wd, 'catalog');
sheet = fixture_dro_elfo_sheet(fullfile(catDir, 'dro_elfo_tau1_raan000.mat'), ...
                               struct('iD', [2 5], 'iA', [2 5], 'nR', 2));
copyfile(sheet, fullfile(catDir, 'dro_elfo_tau1_raan000_conj.mat'));     % decoys
copyfile(sheet, fullfile(catDir, 'dro_elfo_tau1_raan000_deep.mat'));
copyfile(sheet, fullfile(catDir, 'dro_tulip_tau1.mat'));
Q0 = load(sheet);
base = struct('catDir', catDir, 'outDir', wd, 'expectSheets', 1, 'maxRepair', 0);
catMat = fullfile(wd, 'costate_catalog_dro_elfo.mat');
visBefore  = get(0, 'DefaultFigureVisible');
figsBefore = findall(0, 'Type', 'figure');

%% 1. every stage
t1 = tic;
[b1, r1] = build_dro_elfo_catalog(base);
fprintf('  run 1 (%.0f s) blockers: [%s]\n', toc(t1), strjoin(b1, ' | '));
logTxt = fileread(r1.logFile);
for st = {'1. census', '2. clearance', '3/4. repair unsolved', '3/4. repair basin', '5. package', ...
          '6. conj', '7. gates', '8. audit', '9. pictures'}
    assert(contains(logTxt, st{1}), 'stage line "%s" missing from the log', st{1});
end
assert(isfile(catMat), 'catalog written');
C = load(catMat);  cat_ = C.costate_catalog_dro_elfo;
assert(isempty(catalog_schema('validate', cat_)), 'catalog validates');
assert(isscalar(cat_.sheets), 'one packaged sheet, no decoy (got %d)', numel(cat_.sheets));
assert(isequaln(cat_.sheets(1).sheet_key, Q0.sheet_key), 'packaged sheet_key = the sheet''s');
assert(isfield(cat_.sheets, 'conj_pass') && isfield(cat_.sheets, 'gate_min_lamv'), 'sweeps written back');
A = load(fullfile(wd, 'audit_dro_elfo.mat'));
assert(numel(A.results) == 3 && all([A.results.passPick] & [A.results.passFly] & ...
       [A.results.passAccept] & [A.results.passClear]), 'audit replays 3 of 3');
assert(~isempty(dir(fullfile(wd, 'pictures', '*.png'))), 'pictures');
assert(isequal(get(0, 'DefaultFigureVisible'), visBefore), 'DefaultFigureVisible restored');
assert(all(ismember(findall(0, 'Type', 'figure'), figsBefore)), 'no figure left open');
assert(~isfolder(fullfile(catDir, '.dro_elfo_driver.lock')), 'lock released');
assert(~isfolder(fullfile(wd, 'package_stage')), 'package_stage removed');
assert(isfile(r1.verdictFile), 'verdict file');
cp = cat_.sheets(1).conj_pass;
fprintf('  conj_pass per entry: %s; audit picks %s\n', mat2str(double(cp(cat_.sheets(1).has_solution)).'), ...
        strjoin(arrayfun(@(r) sprintf('(%d,%d,%gN)', r.iD, r.iA, r.rungN), A.results, 'UniformOutput', false), ' '));
for kk = 1:3
    fprintf('    audit %d: flown %.3f km, |dz| %.1e, Moon %.0f km, Earth %.0f km\n', kk, ...
            A.results(kk).flownKm, A.results(kk).acceptDz, A.results(kk).moonKm, A.results(kk).earthKm);
end

%% 2. resume
t2 = tic;
[b2, r2] = build_dro_elfo_catalog(base);
log2 = fileread(r2.logFile);
assert(contains(log2, 'content unchanged'), 'the existing catalog kept');
assert(~isempty(regexp(log2, '2\. clearance\s+dro_elfo_tau1_raan000\s[^\n]*\s0 flown now', 'once')), 'no clearance entry flown');
assert(~contains(log2, 'entries to test'), 'no verdict re-tested');
assert(isequal(sort(b2), sort(b1)), 'same blockers on resume');
fprintf('  run 2 (resume) %.0f s\n', toc(t2));

%% 3. injected refuted entry
conjSide = [catMat(1:end-4) '_conjprog.mat'];
P = load(conjSide);  P0 = P;
kPass = find(P.CONJ{1} == 1 & cat_.sheets(1).has_solution, 1);
assert(~isempty(kPass), 'a passing entry to inject into');
P.CONJ{1}(kPass) = 0;  save(conjSide, '-struct', 'P');
b3 = build_dro_elfo_catalog(stages(base, {'conj', 'preflight'}));
assert(any(startsWith(b3, 'conj:')), 'injected refuted entry -> conj: blocker');
r3w = base;  r3w.stages = {'preflight'};  r3w.allowConjBlockers = true;
[b3w, rep3w] = build_dro_elfo_catalog(r3w);
assert(any(startsWith(b3w, 'conj:')) && ~any(startsWith(rep3w.verdictBlockers, 'conj:')) && rep3w.nWaived >= 1, ...
       'allowConjBlockers waives conj: from the verdict only');
save(conjSide, '-struct', 'P0');
build_dro_elfo_catalog(stages(base, {'conj'}));
fprintf('  injected refuted: [%s]\n', strjoin(b3, ' | '));

%% 4. injected stale clearance sidecar: not waivable
clrSide = [sheet(1:end-4) '_clearance.mat'];
E = load(clrSide);  E0 = E;
kOK = find(Q0.OK, 1);
E.Z8SWEPT(1, kOK) = E.Z8SWEPT(1, kOK) + 1e-9;  save(clrSide, '-struct', 'E');
r4 = base;  r4.stages = {'preflight'};  r4.allowConjBlockers = true;  r4.allowClearanceDemotions = true;
[b4, rep4] = build_dro_elfo_catalog(r4);
assert(any(startsWith(b4, 'stale:')) && any(startsWith(rep4.verdictBlockers, 'stale:')) && ~rep4.verified, ...
       'stale clearance sidecar -> stale:, not waivable');
b4b = build_dro_elfo_catalog(stages(base, {'clearance', 'preflight'}));
assert(~any(startsWith(b4b, 'stale:')), 'the clearance stage re-flies the entry and clears stale:');
Ec = load(clrSide);
assert(isequaln(Ec.Z8SWEPT(:, kOK), E0.Z8SWEPT(:, kOK)), 'record re-made for the current z8');
fprintf('  injected stale sidecar: [%s]\n', strjoin(b4, ' | '));

%% 5. changed sheet, package off
Sx = load(sheet);  Sx0 = Sx;
Sx.Z8(1, kOK) = Sx.Z8(1, kOK) + 1e-9;  save(sheet, '-struct', 'Sx');
b5 = build_dro_elfo_catalog(stages(base, {'preflight'}));
assert(any(contains(b5, 'differ from the package')), 'changed sheet: stale: (package)');
save(sheet, '-struct', 'Sx0');

%% 5b. deliberate failures of every other blocking gate (on the verified folder)
allW = base;  allW.allowConjBlockers = true;  allW.allowClearanceDemotions = true;
% (b) an OK entry whose current sidecar record is below the floor: unverified:, not waivable
E = load(clrSide);  E0 = E;
E.MOONKM(kOK) = 100;  save(clrSide, '-struct', 'E');
[bB, repB] = build_dro_elfo_catalog(stages(allW, {'preflight'}));
assert(any(startsWith(bB, 'unverified:') & contains(bB, 'below a clearance floor')) && ~repB.verified && ...
       any(startsWith(repB.verdictBlockers, 'unverified:')), 'OK entry below the floor: unverified:, still blocking with every waiver');
save(clrSide, '-struct', 'E0');
fprintf('  (b) sidecar below floor: [%s]\n', strjoin(bB, ' | '));
% (c) an injected gate failure in the gates sidecar: gates:, not waivable
gatesSide = [catMat(1:end-4) '_gatesprog.mat'];
G = load(gatesSide);  G0 = G;
G.MINLV{1}(kOK) = 0;  save(gatesSide, '-struct', 'G');
[bC, repC] = build_dro_elfo_catalog(stages(allW, {'gates', 'preflight'}));
assert(any(startsWith(bC, 'gates:')) && ~repC.verified, 'injected gate failure: gates:, not waivable');
save(gatesSide, '-struct', 'G0');
build_dro_elfo_catalog(stages(base, {'gates'}));
fprintf('  (c) gate failure: [%s]\n', strjoin(bC, ' | '));
% (d) a sheet carrying the wrong sheet_key: census:
Sk = load(sheet);  Sk0 = Sk;
Sk.sheet_key.raanDeg = 90;  save(sheet, '-struct', 'Sk');
bD = build_dro_elfo_catalog(stages(base, {'census'}));
assert(any(startsWith(bD, 'census:') & contains(bD, 'sheet_key')), 'wrong sheet_key: census:');
save(sheet, '-struct', 'Sk0');
fprintf('  (d) wrong key: [%s]\n', strjoin(bD, ' | '));
% (e) an audit failure: the packaged z8 of the first audited entry off by 1e-9
catKeep = [tempname '.mat'];  copyfile(catMat, catKeep);
Ca = load(catMat);  ca = Ca.costate_catalog_dro_elfo;
col = ca.sheets(1).entry_index(A.results(1).iD, A.results(1).iA, A.results(1).kr);
ca.sheets(1).z8(8, col) = ca.sheets(1).z8(8, col) + 1e-9;
Ca.costate_catalog_dro_elfo = ca;  save(catMat, '-struct', 'Ca');
bE5 = build_dro_elfo_catalog(stages(base, {'audit'}));
assert(any(startsWith(bE5, 'audit:')), 'audit failure: audit:');
copyfile(catKeep, catMat);  delete(catKeep);
build_dro_elfo_catalog(stages(base, {'audit'}));
bOK = build_dro_elfo_catalog(stages(base, {'preflight'}));
assert(isempty(bOK), 'the folder is clean again after the injections (got [%s])', strjoin(bOK, ' | '));
fprintf('  (e) audit failure: [%s]\n', strjoin(bE5, ' | '));

%% 5c. (a) entries below the Moon floor: demoted, waivable, listed either way
catA = fullfile(wd, 'floor', 'catalog');  mkdir(catA);
copyfile(sheet, fullfile(catA, 'dro_elfo_tau1_raan000.mat'));
rA = struct('catDir', catA, 'outDir', fullfile(wd, 'floor'), 'expectSheets', 1, 'maxRepair', 0, ...
            'floorMoonKm', 8000);                 % the fixture's entries pass 7,933 to ~9,200 km
[bA1, repA1] = build_dro_elfo_catalog(rA);
stageMsg = startsWith(bA1, 'clearance-demoted:') & contains(bA1, 'demoted');
assert(any(stageMsg & contains(bA1, '(Moon')), 'the clearance stage reports the demotion');
assert(~repA1.verified, 'demotions block the verdict without allowClearanceDemotions');
VA1 = load(repA1.verdictFile);
assert(any(startsWith(VA1.blockers, 'clearance-demoted:')) && ~VA1.verified, 'verdict file lists the demotion (not waived)');
QA = load(fullfile(catA, 'dro_elfo_tau1_raan000.mat'), 'OK');
assert(nnz(QA.OK) < nnz(Q0.OK), 'entries demoted on the sheet');
bAp = build_dro_elfo_catalog(stages(rA, {'preflight'}));
assert(any(startsWith(bAp, 'clearance-demoted:') & contains(bAp, 'demoted by the clearance sweep')), ...
       'the preflight re-derives the demotion with the clearance stage off');
rA2 = rA;  rA2.allowClearanceDemotions = true;
[bA2, repA2] = build_dro_elfo_catalog(rA2);
VA2 = load(repA2.verdictFile);
assert(repA2.verified && repA2.nWaived >= 1 && any(startsWith(bA2, 'clearance-demoted:')) && ...
       any(startsWith(VA2.blockers, 'clearance-demoted:')) && VA2.verified, ...
       'allowClearanceDemotions: verdict true, the demotion still listed in the verdict file');
fprintf('  (a) floor 8000 km: %d of %d entries demoted; verdict %d -> %d with the waiver\n', ...
        nnz(Q0.OK) - nnz(QA.OK), nnz(Q0.OK), repA1.verified, repA2.verified);

%% 6. refusals, empty folder
lockDir = fullfile(catDir, '.dro_elfo_driver.lock');  mkdir(lockDir);
expect_err(@() build_dro_elfo_catalog(stages(base, {'clearance'})), 'build_dro_elfo_catalog:fleetRunning');
assert(isfolder(lockDir), 'foreign lock survives');  rmdir(lockDir);
expect_err(@() build_dro_elfo_catalog(struct('catdir', catDir)), 'build_dro_elfo_catalog:runIn');
expect_err(@() build_dro_elfo_catalog(stages(base, {'deliverable'})), 'build_dro_elfo_catalog:stages');
mkdir(fullfile(wd, 'empty', 'catalog'));
rE = base;  rE.catDir = fullfile(wd, 'empty', 'catalog');  rE.outDir = fullfile(wd, 'empty');
bE = build_dro_elfo_catalog(rE);
assert(any(startsWith(bE, 'census:')), 'empty folder: census blocker, no throw');

%% 7. a repair changes a sheet -> 2b re-sweeps it (real engine, 6 x 1 x 2 cut)
cat7 = fullfile(wd, 'repair', 'catalog');
s7 = fixture_dro_elfo_sheet(fullfile(cat7, 'dro_elfo_tau1_raan000.mat'), struct('iD', 1:6, 'iA', 2, 'nR', 2));
Q7 = load(s7);
Q7.TF(3, 1, :) = NaN;  Q7.FLYKM(3, 1, :) = NaN;  Q7.ACCDZ(3, 1, :) = NaN;  Q7.RES(3, 1, :) = NaN;
Q7.WALL(3, 1, :) = NaN;  Q7.OK(3, 1, :) = false;  Q7.Z8(:, 3, 1, :) = NaN;  Q7.ATT(3, 1) = 2;
save(s7, '-struct', 'Q7');
r7 = struct('catDir', cat7, 'outDir', fullfile(wd, 'repair'), 'expectSheets', 1, ...
            'stages', {{'clearance', 'repairUnsolved', 'package', 'preflight'}});
[~, rep7] = build_dro_elfo_catalog(r7);
log7 = fileread(rep7.logFile);
tok = regexp(log7, '2b\. clearance re-sweep[^\n]*?(\d+) flown now', 'tokens', 'once');
Q7a = load(s7);
assert(Q7a.OK(3, 1, 1) && ~isempty(tok) && str2double(tok{1}) > 0, ...
       'the repaired cell is back and 2b re-swept its new entries');
C7 = load(fullfile(wd, 'repair', 'costate_catalog_dro_elfo.mat'));
assert(C7.costate_catalog_dro_elfo.sheets(1).has_solution(3, 1, 1), 'the package holds the repaired entry');
fprintf('  repair + 2b: cell 3 recovered (15 N tf %.6f), %s entries re-swept\n', Q7a.TF(3, 1, 1), tok{1});

ok = true;
fprintf('test_build_dro_elfo_catalog: PASS\n');
end

% ------------------------------------------------------------------------
function r = stages(r, list)
% STAGES  runIn with only the listed stages on.  INPUTS: r [1 x 1] struct;
% list {1 x k} char.  OUTPUTS: r [1 x 1] struct.
r.stages = list;
end

function expect_err(fh, id)
% EXPECT_ERR  fh must throw identifier id.  INPUTS: fh [1 x 1] function
% handle; id [1 x n] char.  OUTPUTS: none (throws on a mismatch).
try
    fh();
catch err
    assert(strcmp(err.identifier, id), 'expected %s, got %s (%s)', id, err.identifier, err.message);
    fprintf('  refused as expected: %s\n', err.message);
    return
end
error('test_build_dro_elfo_catalog:noThrow', 'expected %s, nothing thrown', id);
end
