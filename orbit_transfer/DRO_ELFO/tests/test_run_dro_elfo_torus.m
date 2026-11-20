function ok = test_run_dro_elfo_torus(topts)
%% Purpose:
%
%   Test of the DRO -> ELFO / DRO -> tulip fleet front door
%   run_dro_elfo_torus. Sections:
%     1  planOnly contract: 12 default ELFO sheets (tau 1, raanDeg 0:30:330),
%        tags, families, arrParams (raanDeg, IM oe), recipe, default grid
%        and rungs in pass A and pass B; a tulip sheet; option overrides;
%        selection strings; seeds tf0 0.15, tf0B 0.3; planOnly touches no disk
%     2  refusals by name: sel out of range, non-integer raanDeg, a tulip
%        sheet with a raanDeg, two sheets with one tag, a malformed sel
%        string, '3:1' / '1:0:3' / ';' (empty lists), an unknown target, a misspelt option, a one-rung or
%        unordered ladder
%     3  on-disk guards (fake sheet files, no solve, fence off): rung
%        mismatch, grid mismatch, sheet-key mismatch (key variable and
%        meta-derived key), a sheet of unverifiable identity, an
%        unreadable sheet (corruptSheet, file kept); sheet_key re-added to
%        a finished keyless sheet (engine resumes, solves nothing) and a
%        keyed sheet left untouched by a no-op call; the tagged SHEETS
%        SELECTED / COMPLETE lines, the all-sheets line only for the full
%        default list
%     4  the fence: requireFence refuses (run_dro_elfo_torus:noFence) when
%        the pool opener yields no pool (opts.poolFcn test seam, the
%        local no_pool); a zero-cell call with the default cap
%        opens the fence once, announces it and leaves the pool state as
%        it found it
%     5  REAL SOLVE (topts.solve, default true; cold about 4 to 6 min, run
%        it in the background): 1 x 1 sheet, rungs [15 12] N (one rung is
%        refused: the engine cannot resume a one-rung sheet), fenced as a
%        fleet runs, DRO tau 1 phase 0
%        -> ELFO raanDeg 0 phase 0 in a tempname folder; OK true, the
%        sheet_key variable matches, t_f matches 0.0900542599967 ND to
%        1e-6; a second call resumes without solving (file timestamp
%        unchanged, the engine logs the cell as resumed); the completion
%        marker is in combined_progress.txt
%   Nothing is cached; pass struct('solve', false) to run sections 1-4
%   alone (seconds).
%
%% Inputs:
%
%  topts                    struct (optional)       .solve [true] run
%                                                   section 5
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

if nargin < 1 || isempty(topts), topts = struct(); end
doSolve = true;
if isfield(topts, 'solve'), doSolve = topts.solve; end

here    = fileparts(mfilename('fullpath'));
oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(fileparts(here));
setup_dro_elfo_paths();
addpath(fullfile(fileparts(here), 'catalog'));

%% 1  planOnly contract:
P = run_dro_elfo_torus([], struct('planOnly', true));
assert(numel(P) == 12, '12 default sheets, got %d', numel(P));
assert(strcmp(P(1).fence, 'ms fence not opened (planOnly)'), 'planOnly opens no fence');
assert(numel(unique({P.tag})) == 12, 'distinct tags');
assert(strcmp(P(1).tag, 'dro_elfo_tau1_raan000') && strcmp(P(12).tag, 'dro_elfo_tau1_raan330'), ...
       'tag spelling: %s .. %s', P(1).tag, P(12).tag);
[~, ~, ~, oeIM] = im_elfo_optimum();
for k = 1:12
    for pass = {'optsA', 'optsB'}
        o = P(k).(pass{1});
        assert(strcmp(o.depFamily, 'dro') && o.depParams.tau == 1 && o.tauDRO == 1, 'departure DRO tau 1');
        assert(strcmp(o.arrFamily, 'elfo'), 'arrFamily elfo');
        assert(o.arrParams.raanDeg == 30*(k-1), 'raanDeg of sheet %d', k);
        assert(isequal(o.arrParams.oe, oeIM), 'oe from im_elfo_optimum');
        assert(isequal(o.rungs, [15 12 10 7 5 3 2 1.5 1 0.75 0.5]), 'default rungs');
        assert(o.nD == 12 && o.nA == 12 && o.sD0 == 0 && o.sA0 == 0, 'default grid');
        assert(o.keepDirect == true && o.msHardCapSec == 300 && o.thrLock == true, 'recipe flags');
        assert(o.maxCpuSec == 600 && o.maxIter == 6000 && o.N == 400, 'solver caps');
        assert(o.floorKm == 500 && isempty(o.floorEarthKm), 'floors (Moon 500 km, no Earth floor)');
        assert(isinf(o.directGateKm) && o.gateKm == 100, 'gates');
        assert(o.ispS == 1710 && o.m0kg == 150, 'catalog regime');
    end
    assert(abs(P(k).optsA.tf0 - 0.15) < 1e-12, 'pass A tf0 0.15');
    assert(abs(P(k).optsB.tf0 - 0.3) < 1e-12, 'pass B tf0B 0.3 (not the engine''s cold 4.0 ND)');
    assert(isempty(P(k).optsB.cells), 'pass B cells are scoped at run time');
    sk = P(k).sheet_key;
    assert(strcmp(sk.target, 'elfo') && sk.depSize == 1 && sk.raanDeg == 30*(k-1), 'sheet_key');
end
assert(strcmp(P(1).optsA.arrFamily, 'elfo') && P(1).optsA.arrParams.raanDeg == 0, 'sheet 1 contract');

shT = struct('target', {'tulip', 'elfo'}, 'droTau', {2, 0.5}, 'raanDeg', {NaN, 390});
PT  = run_dro_elfo_torus([1 2], struct('planOnly', true, 'sheets', shT, 'nD', 6, 'nA', 6, ...
          'rungs', [15 10 5 0.25], 'tf0', 0.27, 'tf0B', 0.12, 'msHardCapSec', 0, 'keepDirect', false));
assert(strcmp(PT(1).tag, 'dro_tulip_tau2') && strcmp(PT(1).optsA.arrFamily, 'tulip'), 'tulip sheet');
assert(PT(1).optsA.arrParams.Np == 7 && PT(1).optsA.arrParams.pm == -1, 'tulip Np 7, pm -1');
assert(isnan(PT(1).sheet_key.raanDeg) && strcmp(PT(1).sheet_key.target, 'tulip') && PT(1).sheet_key.depSize == 2, 'tulip key');
assert(strcmp(PT(2).tag, 'dro_elfo_tau0.5_raan030') && PT(2).optsA.arrParams.raanDeg == 30, 'raanDeg 390 wraps to 30');
assert(PT(2).optsA.nD == 6 && PT(2).optsB.nA == 6 && isequal(PT(2).optsB.rungs, [15 10 5 0.25]), 'grid + rung overrides');
assert(abs(PT(2).optsB.tf0 - 0.12) < 1e-12, 'tf0B override lands in pass B');
assert(abs(PT(2).optsA.tf0 - 0.27) < 1e-12 && PT(2).optsB.msHardCapSec == 0 && ~PT(2).optsB.keepDirect, 'recipe overrides');

PS = run_dro_elfo_torus('[1 3:5]', struct('planOnly', true));
assert(isequal([PS.raanDeg], [0 60 90 120]), 'selection string [1 3:5]');
PS = run_dro_elfo_torus("2,12", struct('planOnly', true));
assert(isequal([PS.raanDeg], [30 330]), 'selection string 2,12');

wd = tempname;
run_dro_elfo_torus(1, struct('planOnly', true, 'outDir', wd));
assert(~isfolder(wd), 'planOnly created the output folder');

%% 2  refusals by name:
expect_error(@() run_dro_elfo_torus(13, struct('planOnly', true)), 'run_dro_elfo_torus:sel');
expect_error(@() run_dro_elfo_torus([0 1], struct('planOnly', true)), 'run_dro_elfo_torus:sel');
expect_error(@() run_dro_elfo_torus('[1 x]', struct('planOnly', true)), 'run_dro_elfo_torus:sel');
for bad = {'3:1', '1:0:3', ';'}                       % parse to an empty list: never "all"
    expect_error(@() run_dro_elfo_torus(bad{1}, struct('planOnly', true)), 'run_dro_elfo_torus:sel');
end
assert(numel(run_dro_elfo_torus('  ', struct('planOnly', true))) == 12, 'a blank selection means all');
expect_error(@() run_dro_elfo_torus(1, struct('planOnly', true, ...
    'sheets', struct('target', 'elfo', 'droTau', 1, 'raanDeg', 12.5))), 'run_dro_elfo_torus:raanDeg');
expect_error(@() run_dro_elfo_torus(1, struct('planOnly', true, ...
    'sheets', struct('target', 'tulip', 'droTau', 1, 'raanDeg', 0))), 'run_dro_elfo_torus:raanDeg');
expect_error(@() run_dro_elfo_torus([], struct('planOnly', true, ...
    'sheets', struct('target', 'elfo', 'droTau', 1, 'raanDeg', {0, 360}))), 'run_dro_elfo_torus:dupSheet');
expect_error(@() run_dro_elfo_torus(1, struct('planOnly', true, ...
    'sheets', struct('target', 'halo', 'droTau', 1, 'raanDeg', NaN))), 'run_dro_elfo_torus:target');
% one rung: the engine's resume compares size(OK) with [nD nA nR], and MATLAB
% drops a trailing singleton, so a one-rung sheet is never resumed (measured
% 2026-11-19: re-solved on every call); refused, as are unordered rungs
expect_error(@() run_dro_elfo_torus(1, struct('planOnly', true, 'maxcells', 3)), 'run_dro_elfo_torus:opts');
expect_error(@() run_dro_elfo_torus(1, struct('planOnly', true, 'rungs', 15)), 'run_dro_elfo_torus:rungs');
expect_error(@() run_dro_elfo_torus(1, struct('planOnly', true, 'rungs', [12 15])), 'run_dro_elfo_torus:rungs');
expect_error(@() run_dro_elfo_torus(1, struct('planOnly', true, 'rungs', [15 12 12])), 'run_dro_elfo_torus:rungs');

%% 3  on-disk guards (fence off; every refusal fires before the engine):
mkdir(wd);
cleanWd = onCleanup(@() rm_dir(wd));
sheet1 = struct('target', 'elfo', 'droTau', 1, 'raanDeg', 0);
base   = struct('sheets', sheet1, 'nD', 1, 'nA', 1, 'rungs', [15 12], 'outDir', wd, ...
                'msHardCapSec', 0, 'maxCells', 0);
fake   = fullfile(wd, 'dro_elfo_tau1_raan000.mat');
rungs = [15 10];  OK = false(1, 1, 2);
save(fake, 'rungs', 'OK');
expect_error(@() run_dro_elfo_torus(1, base), 'run_dro_elfo_torus:rungMismatch');
rungs = [15 12];  OK = false(2, 1, 2);
save(fake, 'rungs', 'OK');
expect_error(@() run_dro_elfo_torus(1, base), 'run_dro_elfo_torus:gridMismatch');
OK = false(1, 1, 2);
meta = struct('arrFamily', 'elfo', 'arrParams', struct('oe', oeIM, 'raanDeg', 30), ...
              'depParams', struct('tau', 1), 'tauDRO', 1, 'sD0', 0, 'sA0', 0);
save(fake, 'rungs', 'OK', 'meta');
expect_error(@() run_dro_elfo_torus(1, base), 'run_dro_elfo_torus:sheetKey');   % meta-derived key
meta.arrParams.raanDeg = 0;
sheet_key = struct('target', 'elfo', 'depSize', 2, 'raanDeg', 0);
save(fake, 'rungs', 'OK', 'meta', 'sheet_key');
expect_error(@() run_dro_elfo_torus(1, base), 'run_dro_elfo_torus:sheetKey');   % stored key
meta.arrParams.oe = oeIM + [1 0 0 0];
save(fake, 'rungs', 'OK', 'meta');
expect_error(@() run_dro_elfo_torus(1, base), 'run_dro_elfo_torus:sheetKey');   % other ELFO shape
save(fake, 'rungs', 'OK');
expect_error(@() run_dro_elfo_torus(1, base), 'run_dro_elfo_torus:sheetKey');   % no identity at all
% unreadable sheet: named refusal, file left in place
fid = fopen(fake, 'w');  fprintf(fid, 'not a mat file');  fclose(fid);
expect_error(@() run_dro_elfo_torus(1, base), 'run_dro_elfo_torus:corruptSheet');
assert(isfile(fake), 'the corrupt sheet was deleted');
delete(fake);

% sheet_key re-add: an engine-valid finished sheet with matching meta and no
% sheet_key comes back with the key (the engine resumes it, solves nothing)
fabricate_sheet(fake, 0, oeIM);
assert(~ismember('sheet_key', who('-file', fake)), 'fixture already has a key');
bT1 = base;  bT1.runId = 'T1';
run_dro_elfo_torus(1, bT1);
S = load(fake, 'sheet_key');
assert(isfield(S, 'sheet_key') && isequaln(S.sheet_key, struct('target', 'elfo', 'depSize', 1, 'raanDeg', 0)), ...
       'sheet_key not re-added');
comb = fileread(fullfile(wd, 'combined_progress.txt'));
assert(contains(comb, 'SHEETS SELECTED [run T1]: dro_elfo_tau1_raan000') && ...
       contains(comb, 'SELECTED SHEETS COMPLETE [run T1]: dro_elfo_tau1_raan000'), ...
       'tagged SELECTED / COMPLETE lines missing:\n%s', comb);
d1 = dir(fake);  pause(1.1);
run_dro_elfo_torus(1, base);                     % key present: file untouched
d2 = dir(fake);
assert(d1.datenum == d2.datenum, 'a no-op call rewrote a keyed sheet');
fabricate_sheet(fake, 30, oeIM);                 % meta raanDeg 30, no key, at the raan000 tag
expect_error(@() run_dro_elfo_torus(1, base), 'run_dro_elfo_torus:sheetKey');
delete(fake);

% the all-sheets line only for the full DEFAULT 12-sheet list, with every tag
wdAll = tempname;  mkdir(wdAll);
cleanAll = onCleanup(@() rm_dir(wdAll));
for r = 0:30:330
    fabricate_sheet(fullfile(wdAll, sprintf('dro_elfo_tau1_raan%03d.mat', r)), r, oeIM);
end
bAll = rmfield(base, 'sheets');  bAll.outDir = wdAll;  bAll.runId = 'T2';
run_dro_elfo_torus([], bAll);
bAll.runId = 'T3';
run_dro_elfo_torus('[1 2]', bAll);
bAll.runId = 'T4';
bAll.sheets = struct('target', 'elfo', 'droTau', 1, 'raanDeg', {0, 30});
run_dro_elfo_torus([], bAll);                    % explicit list, all selected: still SELECTED
comb = fileread(fullfile(wdAll, 'combined_progress.txt'));
allTags = strjoin(arrayfun(@(r) sprintf('dro_elfo_tau1_raan%03d', r), 0:30:330, 'UniformOutput', false), ' ');
assert(contains(comb, ['CATALOG ALL SHEETS COMPLETE [run T2]: ' allTags]), 'all-sheets line:\n%s', comb);
assert(contains(comb, 'SELECTED SHEETS COMPLETE [run T3]: dro_elfo_tau1_raan000 dro_elfo_tau1_raan030'), 'T3 line');
assert(contains(comb, 'SELECTED SHEETS COMPLETE [run T4]: dro_elfo_tau1_raan000 dro_elfo_tau1_raan030') && ...
       ~contains(comb, 'CATALOG ALL SHEETS COMPLETE [run T4]'), 'T4 line');

%% 4  the fence:
noPool = @(n) no_pool(n);
expect_error(@() run_dro_elfo_torus(1, setfield(rmfield(base, 'msHardCapSec'), 'poolFcn', noPool)), ...
             'run_dro_elfo_torus:noFence');
comb = fileread(fullfile(wd, 'combined_progress.txt'));
assert(contains(comb, 'ms fence UNAVAILABLE -- running UNFENCED: fixture: no pool'), 'refusal not logged:\n%s', comb);
hadPool = ~isempty(current_pool());
PR = run_dro_elfo_torus(1, rmfield(base, 'msHardCapSec'));
comb = fileread(fullfile(wd, 'combined_progress.txt'));
assert(contains(comb, 'ms fence ACTIVE (pool of 1)'), 'fence line missing:\n%s', comb);
assert(strcmp(PR.fence, 'ms fence ACTIVE (pool of 1)'), 'P.fence = %s', PR.fence);
assert(~isempty(current_pool()) == hadPool, 'the front door left the pool state changed');
assert(~isfile(fake), 'zero cells must write no sheet');
fprintf('test_run_dro_elfo_torus: sections 1-4 PASS\n');

%% 5  real solve + resume:
if doSolve
    sv  = rmfield(base, {'maxCells', 'msHardCapSec'});      % fenced, as a fleet runs
    t0  = tic;
    run_dro_elfo_torus(1, sv);
    fprintf('  solve wall %.0f s\n', toc(t0));
    S = load(fake);
    assert(isfield(S, 'sheet_key') && strcmp(S.sheet_key.target, 'elfo') && ...
           S.sheet_key.depSize == 1 && S.sheet_key.raanDeg == 0, 'sheet_key variable');
    assert(S.OK(1, 1, 1), 'the 15 N cell is not OK');
    fprintf('  12 N rung: OK = %d, t_f = %.10f ND\n', S.OK(1, 1, 2), S.TF(1, 1, 2));
    fprintf('  t_f = %.13f ND (reference 0.0900542599967)\n', S.TF(1, 1, 1));
    assert(abs(S.TF(1, 1, 1) - 0.0900542599967) < 1e-6, 't_f off the reference');
    assert(strcmp(S.meta.arrFamily, 'elfo') && S.meta.arrParams.raanDeg == 0, 'meta arrival');
    d1 = dir(fake);
    logFile = fullfile(wd, 'dro_elfo_tau1_raan000_progress.txt');
    nLog = numel(splitlines(fileread(logFile)));
    pause(1.1);
    run_dro_elfo_torus(1, sv);
    d2 = dir(fake);
    assert(d1.datenum == d2.datenum, 'the resumed call rewrote the sheet');
    lines = splitlines(fileread(logFile));
    newL  = strjoin(lines(nLog:end), newline);
    assert(contains(newL, '0 cells x 2 rungs') && contains(newL, '(1 cells resumed/skipped)'), ...
           'resume line missing:\n%s', newL);
    comb = fileread(fullfile(wd, 'combined_progress.txt'));
    assert(contains(comb, 'SELECTED SHEETS COMPLETE: dro_elfo_tau1_raan000'), ...
           'completion marker missing:\n%s', comb);
    fprintf('test_run_dro_elfo_torus: section 5 (solve + resume) PASS\n');
else
    fprintf('test_run_dro_elfo_torus: section 5 SKIPPED (solve = false)\n');
end
ok = true;
fprintf('test_run_dro_elfo_torus: PASS\n');
end

% ------------------------------------------------------------------------
function expect_error(fcn, id)
% EXPECT_ERROR  fcn() must throw identifier id.
% INPUTS: fcn function handle; id char.  OUTPUTS: none.
try
    fcn();
catch err
    assert(strcmp(err.identifier, id), 'expected %s, got %s (%s)', id, err.identifier, err.message);
    return
end
error('test:noThrow', 'expected %s, nothing thrown', id);
end

% ------------------------------------------------------------------------
function fabricate_sheet(f, raanDeg, oe)
% FABRICATE_SHEET  An engine-valid FINISHED 1 x 1 x 2 sheet (cell OK at the
% top rung, ATT 1) whose meta names DRO tau 1 -> ELFO raanDeg, no sheet_key.
% INPUTS: f char; raanDeg double; oe [1 x 4].  OUTPUTS: none (file written).
TF = [0.09 NaN];  TF = reshape(TF, 1, 1, 2);
FLYKM = TF;  ACCDZ = TF;  RES = TF;  WALL = TF;
OK = reshape([true false], 1, 1, 2);
Z8 = nan(8, 1, 1, 2);  Z8(:, 1, 1, 1) = 0;
ATT = 1;  rungs = [15 12];  sD = 0;  sA = 0;
meta = struct('arrFamily', 'elfo', 'arrParams', struct('oe', oe, 'raanDeg', raanDeg), ...
              'depFamily', 'dro', 'depParams', struct('tau', 1), 'tauDRO', 1, 'sD0', 0, 'sA0', 0);
save(f, 'TF', 'FLYKM', 'ACCDZ', 'RES', 'WALL', 'OK', 'Z8', 'ATT', 'rungs', 'sD', 'sA', 'meta');
end

% ------------------------------------------------------------------------
function pool = no_pool(nWorkers)
% NO_POOL  Pool opener that yields no pool with a warning, as capped_pool
% reports an unavailable licence.  INPUTS: nWorkers.  OUTPUTS: pool [].
warning('test:no_pool', 'fixture: no pool for %d worker(s)', nWorkers);
pool = [];
end

% ------------------------------------------------------------------------
function rm_dir(wd)
% RM_DIR  Remove the test folder if present.  INPUTS: wd char.  OUTPUTS: none.
if isfolder(wd), rmdir(wd, 's'); end
end
