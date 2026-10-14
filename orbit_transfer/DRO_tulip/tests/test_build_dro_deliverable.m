function ok = test_build_dro_deliverable()
% TEST_BUILD_DRO_DELIVERABLE  Deliverable-8 builder: (a) a build into a
% temp folder (movies, zip and self-test off) holds every listed file, a
% README with the required strings and no machine paths, and an exported
% catalog whose two content keys equal the library of record's; the
% KEY-BOUND audit gate refuses (b) a chunk whose altContentKey differs,
% (c) an audit dir missing one chunk (coverage), (d) a primaries audit with
% nBad = 1, (e) a primaries audit with another contentKey, (f) a chunk with
% another contentKey, (g) a duplicated alternative index with the count
% preserved -- each with error id build_dro_deliverable:audit -- (h) a
% refused build leaves an earlier outDir in place, and (i) the USER profile
% (deliverable 9) builds exactly the user file list (no flight-layer file,
% no tfMin checker, no schema), a README that names every shipped file, the
% user example and the statuses, mentions no dropped tool, and the same two
% content keys; an unknown profile is refused.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
dro = fileparts(fileparts(mfilename('fullpath')));            % DRO_tulip
ind = fullfile(dro, 'indirect');
cc  = fullfile(fileparts(dro), 'costate_common');
addpath(dro, ind, cc);
v2Dir  = fullfile(ind, 'results', 'library_70mN_24x48_v2');
record = fullfile(v2Dir, 'costate_catalog_dro_tulip_70mN.mat');
if ~isfile(record)
    fprintf('  SKIP  the library of record is not on disk\n');  return
end
R = load(record);  fn = fieldnames(R);  rec = R.(fn{1});
nAlt = numel(rec.alternatives);
base = struct('movies', false, 'zip', false, 'selftest', false);

%% (a) a build
work = tempname;  mkdir(work);
cfg = base;  cfg.outDir = fullfile(work, 'costate_catalog_70mN_status_deliverable8');
try
    out = build_dro_deliverable(cfg);
    built = true;
catch ME
    built = false;
    fprintf('  build threw: %s (%s)\n', ME.message, ME.identifier);
end
ok = chk(ok, built, 'build (movies/zip/selftest off) completes');
if built
    L = deliverable_file_list();
    have = arrayfun(@(e) isfile(fullfile(cfg.outDir, e.file)), L);
    ok = chk(ok, all(have), sprintf('every listed file is in outDir (%d/%d)', nnz(have), numel(L)));
    ok = chk(ok, isfile(fullfile(cfg.outDir, 'status_torus.png')), 'status_torus.png written');
    rd = fullfile(cfg.outDir, 'README.md');
    ok = chk(ok, isfile(rd), 'README.md written');
    txt = '';
    if isfile(rd), txt = fileread(rd); end
    must = {sprintf('%d alternatives', nAlt), 'conjugate point found', 'borderline', ...
            'fractions of a period', 'DAYS'};
    for k = 1:numel(must)
        ok = chk(ok, contains(txt, must{k}), sprintf('README contains ''%s''', must{k}));
    end
    ok = chk(ok, ~contains(txt, '/Users/') && ~contains(txt, 'addpath('), ...
             'README has no machine path and no addpath(');
    bl = find([rec.alternatives.borderline]);
    rowsListed = arrayfun(@(k) contains(txt, sprintf('| %d |', k)), bl);
    ok = chk(ok, all(rowsListed), sprintf('README lists the %d borderline rows by index', numel(bl)));
    for k = 1:numel(L)
        ok = chk(ok, contains(txt, ['`' L(k).file '`']), sprintf('contents table names %s', L(k).file));
    end
    E = load(fullfile(cfg.outDir, 'costate_catalog_dro_tulip_70mN.mat'));
    fe = fieldnames(E);  ex = E.(fe{1});
    ok = chk(ok, strcmp(catalog_content_key(ex), catalog_content_key(rec)), ...
             'exported catalog_content_key equals the record''s');
    ok = chk(ok, strcmp(alternatives_content_key(ex), alternatives_content_key(rec)), ...
             'exported alternatives_content_key equals the record''s');
    ok = chk(ok, isfield(out, 'keys') && strcmp(out.keys.primKey, catalog_content_key(rec)) ...
             && strcmp(out.keys.altKey, alternatives_content_key(rec)), 'out.keys carries both keys');
    ok = chk(ok, all(isfield(out, {'dir', 'zip', 'facts', 'files', 'selftest', 'keys'})), ...
             'out has .dir .zip .facts .files .selftest .keys');
end

%% audit-dir copies for the refusals
D = dir(fullfile(v2Dir, 'audit_alt_*.mat'));
ok = chk(ok, numel(D) >= 2, sprintf('V2 holds %d audit_alt chunks', numel(D)));

% (b) one chunk's altContentKey altered, coverage intact
bDir = fullfile(work, 'audit_b');  mkdir(bDir);
for k = 1:numel(D), copyfile(fullfile(D(k).folder, D(k).name), bDir); end
X = load(fullfile(bDir, D(1).name));
X.altContentKey = repmat('0', 1, 32);
save(fullfile(bDir, D(1).name), '-struct', 'X');
cfgB = base;  cfgB.outDir = fullfile(work, 'deliv_b');  cfgB.auditAltDir = bDir;
ok = chk(ok, refusesAudit(cfgB), 'a chunk with an altered altContentKey: refused (build_dro_deliverable:audit)');

% (c) one chunk missing
cDir = fullfile(work, 'audit_c');  mkdir(cDir);
for k = 2:numel(D), copyfile(fullfile(D(k).folder, D(k).name), cDir); end
cfgC = base;  cfgC.outDir = fullfile(work, 'deliv_c');  cfgC.auditAltDir = cDir;
ok = chk(ok, refusesAudit(cfgC), 'an audit dir missing one chunk: refused (coverage)');

% (d) primaries audit with nBad = 1
P = load(fullfile(fileparts(fileparts(record)), 'library_70mN_24x48_merged', 'audit_70mN.mat'));
P.A.nBad = 1;
pMat = fullfile(work, 'audit_70mN_bad.mat');  save(pMat, '-struct', 'P');
cfgD = base;  cfgD.outDir = fullfile(work, 'deliv_d');  cfgD.auditPrimariesMat = pMat;
ok = chk(ok, refusesAudit(cfgD), 'a primaries audit with nBad = 1: refused');

% (e) primaries audit with an altered contentKey
P = load(fullfile(fileparts(fileparts(record)), 'library_70mN_24x48_merged', 'audit_70mN.mat'));
P.A.contentKey = repmat('f', 1, 32);
eMat = fullfile(work, 'audit_70mN_key.mat');  save(eMat, '-struct', 'P');
cfgE = base;  cfgE.outDir = fullfile(work, 'deliv_e');  cfgE.auditPrimariesMat = eMat;
ok = chk(ok, refusesAudit(cfgE), 'a primaries audit with an altered contentKey: refused');

% (f) one chunk's (primaries) contentKey altered, coverage intact
fDir = fullfile(work, 'audit_f');  mkdir(fDir);
for k = 1:numel(D), copyfile(fullfile(D(k).folder, D(k).name), fDir); end
X = load(fullfile(fDir, D(2).name));
X.contentKey = repmat('1', 1, 32);
save(fullfile(fDir, D(2).name), '-struct', 'X');
cfgF = base;  cfgF.outDir = fullfile(work, 'deliv_f');  cfgF.auditAltDir = fDir;
ok = chk(ok, refusesAudit(cfgF), 'a chunk with an altered contentKey: refused');

% (g) a duplicated k inside one chunk, row count preserved (one k missing, one twice)
gDir = fullfile(work, 'audit_g');  mkdir(gDir);
for k = 1:numel(D), copyfile(fullfile(D(k).folder, D(k).name), gDir); end
X = load(fullfile(gDir, D(3).name));
X.altRows(1).k = X.altRows(2).k;
save(fullfile(gDir, D(3).name), '-struct', 'X');
cfgG = base;  cfgG.outDir = fullfile(work, 'deliv_g');  cfgG.auditAltDir = gDir;
ok = chk(ok, refusesAudit(cfgG), 'a duplicated alternative index (count preserved): refused');

% (h) a refused build leaves an earlier deliverable folder untouched
keep = fullfile(work, 'deliv_keep');  mkdir(keep);
fid = fopen(fullfile(keep, 'README.md'), 'w');  fprintf(fid, 'earlier build\n');  fclose(fid);
cfgH = base;  cfgH.outDir = keep;  cfgH.auditAltDir = cDir;
refusesAudit(cfgH);
ok = chk(ok, isfile(fullfile(keep, 'README.md')), 'a refused audit leaves the earlier outDir in place');

%% (i) the user profile (deliverable 9)
cfgI = base;  cfgI.profile = 'user';  cfgI.outDir = fullfile(work, 'costate_library_70mN_user_deliverable9');
try
    outI = build_dro_deliverable(cfgI);
    builtI = true;
catch ME
    builtI = false;
    fprintf('  user build threw: %s (%s)\n', ME.message, ME.identifier);
end
ok = chk(ok, builtI, 'user-profile build (movies/zip/selftest off) completes');
if builtI
    LU = deliverable_file_list('user');
    shipped = dir(fullfile(cfgI.outDir, '*.m'));  shipped = {shipped.name};
    ok = chk(ok, isequal(sort(shipped(:)), sort({LU.file}.')), ...
             sprintf('outDir holds exactly the %d user-list .m files (%d found)', numel(LU), numel(shipped)));
    dropped = {'fly_catalog_entry.m', 'catalog_endpoints.m', 'fly_transfer.m', 'entry_thrust_program.m', ...
               'validate_flight.m', 'scalar_verdict.m', 'phase_state.m', 'periodic_pp.m', 'tfmin_accepts.m', ...
               'run_capped.m', 'current_pool.m', 'catalog_schema.m', 'sample_status_pair.m', ...
               'print_transfer_summary.m', 'plot_transfer_3d.m', 'status_layer_example.m', 'costate_catalog_example.m'};
    ok = chk(ok, ~any(ismember(dropped, shipped)), 'no flight-layer, tfMin-check, schema or developer-example file shipped');
    ok = chk(ok, isfile(fullfile(cfgI.outDir, 'status_torus.png')), 'user build writes status_torus.png');
    rdI = fullfile(cfgI.outDir, 'README.md');
    txt = '';
    if isfile(rdI), txt = fileread(rdI); end
    mustI = {'User Kit', sprintf('%d alternatives', nAlt), 'conjugate point found', 'borderline', ...
             'fractions of a period', 'DAYS', 'costate_library_example', 'pumpkyn.cr3bp.tfMin('};
    for k = 1:numel(mustI)
        ok = chk(ok, contains(txt, mustI{k}), sprintf('user README contains ''%s''', mustI{k}));
    end
    for k = 1:numel(LU)
        ok = chk(ok, contains(txt, ['`' LU(k).file '`']), sprintf('user README contents table names %s', LU(k).file));
    end
    notI = {'fly_catalog_entry', 'tfmin_accepts', 'catalog_endpoints', 'status_layer_example', 'Parallel Computing'};
    for k = 1:numel(notI)
        ok = chk(ok, ~contains(txt, notI{k}), sprintf('user README does not mention the dropped %s', notI{k}));
    end
    ok = chk(ok, ~contains(txt, '/Users/') && ~contains(txt, 'addpath('), 'user README has no machine path and no addpath(');
    EI = load(fullfile(cfgI.outDir, 'costate_catalog_dro_tulip_70mN.mat'));
    feI = fieldnames(EI);  exI = EI.(feI{1});
    ok = chk(ok, strcmp(catalog_content_key(exI), catalog_content_key(rec)) && ...
             strcmp(alternatives_content_key(exI), alternatives_content_key(rec)), ...
             'user build''s exported catalog carries the record''s two content keys');
    ok = chk(ok, strcmp(outI.profile, 'user'), 'out.profile = user');
end
cfgJ = base;  cfgJ.profile = 'customer';  cfgJ.outDir = fullfile(work, 'deliv_j');
refusedJ = false;
try, build_dro_deliverable(cfgJ); catch ME, refusedJ = strcmp(ME.identifier, 'build_dro_deliverable:profile'); end
ok = chk(ok, refusedJ, 'an unknown profile is refused (build_dro_deliverable:profile)');

if ok, fprintf('test_build_dro_deliverable: ALL PASS\n'); else, fprintf('test_build_dro_deliverable: FAIL\n'); end
end

function r = refusesAudit(cfg)
% REFUSESAUDIT  True when the build throws build_dro_deliverable:audit.
% INPUTS: cfg [struct].  OUTPUTS: r [logical].
r = false;
try
    build_dro_deliverable(cfg);
    fprintf('  (build did NOT refuse)\n');
catch ME
    r = strcmp(ME.identifier, 'build_dro_deliverable:audit');
    fprintf('  (threw %s: %s)\n', ME.identifier, ME.message);
end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
