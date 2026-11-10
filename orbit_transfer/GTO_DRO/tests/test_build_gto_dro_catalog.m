function ok = test_build_gto_dro_catalog()
% TEST_BUILD_GTO_DRO_CATALOG  The chain on a one-sheet, one-cell catalog
% folder built from the demo fixture: census -> earth -> package -> conj ->
% gates -> audit -> pictures -> deliverable, every stage fenced, blockers
% returned; the repair stages are switched off here (Task 4 tests them).
%
%   1. full chain: package (sheet key = orientDeg, dep_params completed,
%      only the official sheet packaged although the catalog folder also
%      holds sidecars and a decoy), sweeps written back, audit
%      gates, pictures (no figure left open, DefaultFigureVisible restored),
%      and the refuted 15 N entry carried as a conj blocker that refuses the
%      zip
%   2. second run with allowConjBlockers (package off, sweep sidecars
%      reused): the zip is produced, the returned blockers still carry the
%      conj line
%   3. fleet lock: a foreign lock makes the chain refuse (and survives); a
%      normal run leaves no lock
%   4. deliverable preflight: a changed sheet with the package stage off ->
%      stale:; a fresh package without the sweeps -> unverified:; the old
%      sweep backups are kept under dated names; package_stage is removed
%   5. pickers resolve the DRO-arrival sheet by (orientation, DRO period)
%   6. an empty catalog folder returns a census blocker and does not throw
%   7. a FRESH folder, every stage on (repairs included), maxRepair = 1: the
%      basin repair re-solves the refuted 15 N entry, 2b re-sweeps it, the
%      package holds the repaired entry, and the chain ships with no waiver
%   8. the preflight re-derives, with the stages switched off: a census
%      shortfall, a gate failure corrupted into the catalog (not waivable)
%      and Earth demotions (impossible floor)
%
% INPUTS: none.  OUTPUTS: ok - logical.
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
addpath(here);  addpath(fullfile(fileparts(here), 'catalog'));
src = fixture_demo_sheet();
wd  = fullfile(tempdir, 'test_build_gto_dro');
if isfolder(wd), rmdir(wd, 's'); end
mkdir(fullfile(wd, 'catalog'));
copyfile(src, fullfile(wd, 'catalog', 'gto_dro_o000_tau1.mat'));
% a decoy that matches the packager's glob but is not an official sheet
copyfile(src, fullfile(wd, 'catalog', 'gto_dro_o000_tau1_deep.mat'));

visBefore  = get(0, 'DefaultFigureVisible');
figsBefore = findall(0, 'Type', 'figure');

run_ = struct('census', true, 'earth', true, 'repairUnsolved', false, 'repairBasin', false, ...
              'package', true, 'conj', true, 'gates', true, 'audit', true, 'pictures', true, ...
              'deliverable', true, 'catDir', fullfile(wd, 'catalog'), 'outDir', wd, ...
              'expectSheets', 1);

%% 1. full chain
blockers = build_gto_dro_catalog(run_);
catMat = fullfile(wd, 'costate_catalog_gto_dro.mat');
assert(isfile(catMat), 'package written');
C = load(catMat);  cat_ = C.costate_catalog_gto_dro;
assert(numel(cat_.sheets) == 1 && cat_.sheets(1).tauDRO == 0, ...
       'sheet key = orientDeg; the _deep decoy in the catalog folder was not packaged');
assert(cat_.sheets(1).dep_params.orientDeg == 0 && isfield(cat_.sheets(1).dep_params, 'sma_km'), 'dep_params completed');
assert(cat_.sheets(1).arr_params.tau == 1, 'arrival DRO period kept');
assert(isfield(cat_.sheets(1), 'conj_pass') && isfield(cat_.sheets(1), 'gate_min_lamv'), 'sweeps written back');
assert(isfile(fullfile(wd, 'audit_gto_dro.mat')), 'audit written');
A = load(fullfile(wd, 'audit_gto_dro.mat'));
assert(all([A.results.passFly]) && all([A.results.passAccept]) && all([A.results.passEarth]), 'audit gates');
assert(~isempty(dir(fullfile(wd, 'pictures', '*.png'))), 'pictures');
assert(isequal(get(0, 'DefaultFigureVisible'), visBefore), 'DefaultFigureVisible restored');
figsNow = findall(0, 'Type', 'figure');
assert(all(ismember(figsNow, figsBefore)), 'the picture stage left no figure open');
z = dir(fullfile(wd, 'deliverable', '*.zip'));
% the 15 N entry of the fixture is the refuted basin, so the conj sweep
% reports one FAIL and the chain must carry it as a blocker and refuse the zip
assert(any(startsWith(blockers, 'conj:')), 'refuted entry must be a blocker');
cpass = cat_.sheets(1).conj_pass;
assert(cpass(1, 1, 1) == 0 && cpass(1, 1, 2) == 1, 'conj verdict: FAIL at the 15 N entry, PASS at the 12 N entry');
assert(~isfolder(fullfile(wd, '.gtodrocat_driver.lock')), 'lock released after a normal run');
assert(~isfolder(fullfile(wd, 'package_stage')), 'package_stage removed after a successful package');
assert(isempty(z), 'deliverable refused while a blocker stands');
fprintf('  blockers: %s\n', strjoin(blockers, ' | '));

%% 2. allowConjBlockers: package off (sweep sidecars reused), conj on
run2 = run_;
run2.package = false;  run2.audit = false;  run2.pictures = false;
run2.allowConjBlockers = true;
blockers2 = build_gto_dro_catalog(run2);
assert(any(startsWith(blockers2, 'conj:')), 'the conj blocker stays in the returned list');
z = dir(fullfile(wd, 'deliverable', '*.zip'));
assert(numel(z) == 1, 'deliverable produced when the conj blocker is waived');
ex = fullfile(wd, 'unzipped');
unzip(fullfile(z.folder, z.name), ex);
inZip = dir(fullfile(ex, '**', '*'));
inZip = {inZip(~[inZip.isdir]).name};
assert(all(ismember({'costate_catalog_gto_dro.mat', 'costate_catalog_pick.m', ...
                     'costate_lib_describe.m', 'README.md'}, inZip)), 'zip contents');
fid = fopen(fullfile(fileparts(here), 'README.md'), 'r');  bytes = fread(fid, inf, 'uint8');  fclose(fid);
assert(all(bytes < 128), 'README is ASCII');

%% 3. fleet lock: a foreign lock refuses the chain and survives it
lockPath = fullfile(wd, '.gtodrocat_driver.lock');
mkdir(lockPath);
runL = run_;  runL.package = false;  runL.conj = false;  runL.gates = false;
runL.audit = false;  runL.pictures = false;  runL.deliverable = false;
gotId = '';
try
    build_gto_dro_catalog(runL);
catch err
    gotId = err.identifier;
end
assert(strcmp(gotId, 'build_gto_dro_catalog:fleetRunning'), 'a foreign lock refuses the chain');
assert(isfolder(lockPath), 'the foreign lock survives');
rmdir(lockPath);

%% 4. deliverable preflight
sheetPath = fullfile(wd, 'catalog', 'gto_dro_o000_tau1.mat');
sheetKeep = fullfile(wd, 'sheet_keep.mat');
copyfile(sheetPath, sheetKeep);
delete(fullfile(wd, 'deliverable', '*.zip'));
Sx = load(sheetPath);  Sx.Z8(1, 1, 1, 1) = Sx.Z8(1, 1, 1, 1) + 1e-9;  save(sheetPath, '-struct', 'Sx');
runS = run_;  runS.package = false;  runS.earth = false;  runS.conj = false;  runS.gates = false;
runS.audit = false;  runS.pictures = false;  runS.allowConjBlockers = true;
blockersS = build_gto_dro_catalog(runS);
assert(any(startsWith(blockersS, 'stale:')), 'a changed sheet with package off: stale blocker');
assert(isempty(dir(fullfile(wd, 'deliverable', '*.zip'))), 'stale package refused');
copyfile(sheetKeep, sheetPath);
runU = run_;  runU.earth = false;  runU.conj = false;  runU.gates = false;
runU.audit = false;  runU.pictures = false;  runU.allowConjBlockers = true;
blockersU = build_gto_dro_catalog(runU);
assert(any(startsWith(blockersU, 'unverified:')), 'a fresh package without the sweeps: unverified blocker');
assert(isempty(dir(fullfile(wd, 'deliverable', '*.zip'))), 'unverified package refused');
assert(~isempty(dir([catMat '.bak_conj_*'])), 'old conj backup kept under a dated name');
assert(~isfile([catMat '.bak_conj']), 'no undated conj backup left');
fprintf('  stale blockers: %s\n  unverified blockers: %s\n', strjoin(blockersS, ' | '), strjoin(blockersU, ' | '));

%% 5. pickers
Pd_days = get_dep_period_days(cat_);
[tf_nd, z8, info] = costate_catalog_pick(cat_, 0, 1, 0.10*Pd_days, 0, 12, false);
assert(isfinite(tf_nd) && numel(z8) == 8, 'picker returns the 12 N entry');
assert(~info.warned && abs(info.delivered.depDays - 0.10*Pd_days) < 1e-9 && info.delivered.arrKey == 1, ...
       'picker resolves the DRO sheet and phase without a warning');
D = costate_lib_describe(cat_);
assert(isstruct(D) && D.nEntries == 2, 'describe runs');

%% 6. an empty catalog folder: a census blocker, no throw
mkdir(fullfile(wd, 'empty', 'catalog'));
runE = run_;
runE.repairUnsolved = true;  runE.repairBasin = true;
runE.catDir = fullfile(wd, 'empty', 'catalog');  runE.outDir = fullfile(wd, 'empty');
blockersE = build_gto_dro_catalog(runE);
assert(any(startsWith(blockersE, 'census:')), 'empty folder: census blocker');
assert(isempty(dir(fullfile(wd, 'empty', 'deliverable', '*.zip'))), 'empty folder: nothing shipped');
assert(isequal(get(0, 'DefaultFigureVisible'), visBefore), 'DefaultFigureVisible restored (empty)');

%% 7. a fresh folder, every stage on, maxRepair = 1
wd3 = fullfile(tempdir, 'test_build_gto_dro_repair');
if isfolder(wd3), rmdir(wd3, 's'); end
mkdir(fullfile(wd3, 'catalog'));
copyfile(src, fullfile(wd3, 'catalog', 'gto_dro_o000_tau1.mat'));
run3 = run_;
run3.repairUnsolved = true;  run3.repairBasin = true;  run3.maxRepair = 1;
run3.catDir = fullfile(wd3, 'catalog');  run3.outDir = wd3;
blockers3 = build_gto_dro_catalog(run3);
fprintf('  third run blockers: [%s]\n', strjoin(blockers3, ' | '));
logFiles = dir(fullfile(wd3, 'chain_*.log'));
logTxt = fileread(fullfile(logFiles(1).folder, logFiles(1).name));
tok = regexp(logTxt, '2b\. earth re-sweep[^\n]*?(\d+) flown now', 'tokens', 'once');
assert(~isempty(tok) && str2double(tok{1}) > 0, '2b re-swept the repaired sheet (entries flown now > 0)');
C3 = load(fullfile(wd3, 'costate_catalog_gto_dro.mat'));  c3 = C3.costate_catalog_gto_dro;
assert(abs(c3.sheets(1).tf_nd(1, 1, 1) - 0.26682) < 1e-4, 'the package holds the repaired 15 N entry');
assert(all(c3.sheets(1).conj_pass(:) == 1), 'conj sweep: no FAIL after the repair');
assert(isempty(blockers3), 'the chain is clean');
assert(numel(dir(fullfile(wd3, 'deliverable', '*.zip'))) == 1, 'deliverable shipped with no waiver');
assert(~isfolder(fullfile(wd3, '.gtodrocat_driver.lock')), 'lock released (repair run)');

%% 8. the preflight re-derives what switched-off stages would report
% 8a. census shortfall, deliverable stage only, on the verified package of run 3
delete(fullfile(wd3, 'deliverable', '*.zip'));
runD = struct('census', false, 'earth', false, 'repairUnsolved', false, 'repairBasin', false, ...
              'package', false, 'conj', false, 'gates', false, 'audit', false, 'pictures', false, ...
              'deliverable', true, 'catDir', fullfile(wd3, 'catalog'), 'outDir', wd3, 'expectSheets', 2);
blockersD = build_gto_dro_catalog(runD);
assert(any(startsWith(blockersD, 'census:')), 'census shortfall re-derived with the census stage off');
assert(isempty(dir(fullfile(wd3, 'deliverable', '*.zip'))), 'census shortfall refused');
% 8b. the clean package still ships through the preflight alone
runD.expectSheets = 1;
blockersD = build_gto_dro_catalog(runD);
assert(isempty(blockersD) && numel(dir(fullfile(wd3, 'deliverable', '*.zip'))) == 1, ...
       'a verified package ships with only the deliverable stage on');
% 8c. corrupt one gate value so it fails H2; deliverable stage only
delete(fullfile(wd3, 'deliverable', '*.zip'));
catMat3 = fullfile(wd3, 'costate_catalog_gto_dro.mat');
G = load(catMat3);  gcat = G.costate_catalog_gto_dro;
gcat.sheets(1).gate_min_lamv(1, 1, 1) = 0;
Gout = struct('costate_catalog_gto_dro', gcat);  save(catMat3, '-struct', 'Gout');
blockersG = build_gto_dro_catalog(runD);
fprintf('  corrupted-gate blockers: %s\n', strjoin(blockersG, ' | '));
assert(any(startsWith(blockersG, 'gates:')), 'gate failure on disk re-derived with the gates stage off');
assert(isempty(dir(fullfile(wd3, 'deliverable', '*.zip'))), 'gate-failing package refused');
runD.allowConjBlockers = true;  runD.allowEarthDemotions = true;
blockersG = build_gto_dro_catalog(runD);
assert(any(startsWith(blockersG, 'gates:')) && isempty(dir(fullfile(wd3, 'deliverable', '*.zip'))), ...
       'the gates: blocker cannot be waived');
% 8d. Earth demotions: an impossible floor demotes everything; a later
% deliverable-only call must still report earth-demoted: (earth off)
wd4 = fullfile(tempdir, 'test_build_gto_dro_demote');
if isfolder(wd4), rmdir(wd4, 's'); end
mkdir(fullfile(wd4, 'catalog'));
copyfile(src, fullfile(wd4, 'catalog', 'gto_dro_o000_tau1.mat'));
runM = struct('census', false, 'earth', true, 'repairUnsolved', false, 'repairBasin', false, ...
              'package', false, 'conj', false, 'gates', false, 'audit', false, 'pictures', false, ...
              'deliverable', false, 'catDir', fullfile(wd4, 'catalog'), 'outDir', wd4, ...
              'expectSheets', 1, 'floorEarthKm', 1e6);
blockersM = build_gto_dro_catalog(runM);
assert(any(startsWith(blockersM, 'earth-demoted:')), 'impossible floor: entries demoted');
runM.earth = false;  runM.deliverable = true;
blockersM = build_gto_dro_catalog(runM);
fprintf('  demotion re-derived blockers: %s\n', strjoin(blockersM, ' | '));
assert(any(startsWith(blockersM, 'earth-demoted:')), 'earth-demoted: re-derived with the earth stage off');
assert(isempty(dir(fullfile(wd4, 'deliverable', '*.zip'))), 'nothing shipped');

rmdir(wd, 's');
rmdir(wd3, 's');
rmdir(wd4, 's');
ok = true;
fprintf('test_build_gto_dro_catalog: PASS\n');
end

function Pd = get_dep_period_days(cat_)
% GET_DEP_PERIOD_DAYS  The GTO Kepler period in days, from the family
% provider (first output of get_family_orbit, ND) and the catalog's tStar:
% the units costate_catalog_pick's departure phase uses.  INPUTS: cat_.
% OUTPUTS: Pd [scalar days].
tau = get_family_orbit('gto', cat_.sheets(1).dep_params);
Pd = tau(end) * cat_.constants.tStar_s / 86400;
end
