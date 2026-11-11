function ok = test_ladder_keepdirect()
%% Purpose:
%
%   Tests opts.keepDirect of thrust_ladder_library: the sidecar
%   <outMat minus .mat>_direct.mat that keeps every rung's DIRECT
%   (Hermite-Simpson) solution. Engine checks run on the stage 1 GTO -> DRO
%   demo cell (orient 0, DRO tau 1, departure phase 0.10, arrival phase 0,
%   Earth floor 200 km; the cell of test_ladder_seedtop):
%     (a) DEFAULT-OFF IDENTITY: rung 15 N solved once without the field and
%         once with keepDirect = false; every sheet variable is isequaln
%         except WALL (wall-clock seconds, never repeatable), and no
%         *_direct.mat is written
%     (b) keepDirect on, rungs [15 12]: DDONE true exactly where a direct
%         solve succeeded (here both rungs, so DDONE == isfinite(TF)); DTF
%         finite there; DX(1:6,1) equals the departure state to 1e-9; DU
%         columns unit norm to 1e-6; DEARTHKM <= the node-only Earth
%         altitude recomputed from DX (1e-6 km; DEARTHKM also covers the
%         midpoints) and DMOONKM likewise; identity fields; the sheet's
%         15 N rung equals run (a)'s except WALL (option on changes no
%         sheet value)
%     (c) RESUME leaves recorded slices untouched (cell holds OK rungs, so
%         it is skipped); a sidecar whose N does not match is replaced by a
%         fresh one, with one log line
%     (d) REDO clears a cell's slices BEFORE solving (a redo whose direct
%         solve fails at maxIter = 3 leaves DDONE false and DX NaN), then a
%         redo seeded from the kept direct solution rewrites them
%     (e) DELIBERATE FAILURE: the sidecar checker of (b) trips on a copy
%         with one DDONE entry corrupted
%   Two rungs from (b) on because the engine's resume check compares
%   size(OK) with [nD nA nR], which a 1 x 1 x 1 sheet never matches
%   (MATLAB drops the trailing singleton): a ONE-rung sheet is never
%   resumed (pre-existing, reported, not changed here). About 8 min: each
%   15 N cold solve with refinement is about 100 s.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 true when every check
%                                                   passes (errors otherwise)
%
%% Revision History:
%  M. Casey                                                   (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
if isempty(which('oc.ms_bvp'))
    addpath(fullfile(fileparts(fileparts(fileparts(here))), 'oclib'));
end
if isempty(which('pumpkyn.cr3bp.tfMin'))
    od = cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(od);
end
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));

outDir = [tempname '_keepdirect'];
mkdir(outDir);
cleanDir = onCleanup(@() rmdir_quiet(outDir));

N    = 400;
base = struct('rungs', 15, 'ispS', 1710, 'm0kg', 150, 'N', N, ...
    'floorKm', 500, 'floorEarthKm', 200, 'maxIter', 6000, 'maxCpuSec', 600, ...
    'gateKm', 100, 'accTol', 1e-6, 'nD', 1, 'nA', 1, 'sD0', 0.10, 'sA0', 0, ...
    'depFamily', 'gto', 'depParams', struct('orientDeg', 0), ...
    'arrFamily', 'dro', 'arrParams', struct('tau', 1), 'tauDRO', 0, ...
    'thrLock', true, 'tf0', 0.30, 'resume', false);
sheetVars = {'TF', 'FLYKM', 'ACCDZ', 'RES', 'OK', 'Z8', 'ATT', 'rungs', 'sD', 'sA', 'meta'};

% (a) default off: no field vs keepDirect = false
o = base;  o.logFile = fullfile(outDir, 'nofield.log');
matA = fullfile(outDir, 'nofield.mat');
thrust_ladder_library(matA, o);
o = base;  o.keepDirect = false;  o.logFile = fullfile(outDir, 'off.log');
matB = fullfile(outDir, 'off.mat');
thrust_ladder_library(matB, o);
QA = load(matA);  QB = load(matB);
assert(isequal(sort(fieldnames(QA)), sort(fieldnames(QB))), 'sheet variables differ with keepDirect = false');
for kv = 1:numel(sheetVars)
    assert(isequaln(QA.(sheetVars{kv}), QB.(sheetVars{kv})), ...
        'sheet variable %s differs between no field and keepDirect = false', sheetVars{kv});
end
assert(QA.OK(1, 1, 1), 'the 15 N demo rung must be OK');
d = dir(fullfile(outDir, '*_direct.mat'));
assert(isempty(d), 'option off wrote a sidecar: %s', strjoin({d.name}, ', '));
fprintf('  (a) option off: sheet isequaln to no-field run (WALL excluded), tf(15 N) = %.15f, no sidecar\n', QA.TF(1, 1, 1));

% (b) keepDirect on, rungs [15 12]
o = base;  o.rungs = [15 12];  o.keepDirect = true;  o.logFile = fullfile(outDir, 'on.log');
matC = fullfile(outDir, 'on.mat');
sideC = fullfile(outDir, 'on_direct.mat');
P = thrust_ladder_library(matC, o);
assert(isfile(sideC), 'keepDirect = true wrote no sidecar %s', sideC);
SC = load(sideC);
[tD, rvD] = get_family_orbit('gto', struct('orientDeg', 0));
depAt = phase_state(tD, rvD);
rv0 = depAt(0.10);
check_sidecar(SC, P, rv0, N, 'on.mat');
fprintf('  (b) DDONE %s, DTF %s ND (sheet TF %s), DEARTHKM %s km, DMOONKM %s km\n', ...
    mat2str(squeeze(SC.DDONE).'), mat2str(squeeze(SC.DTF).', 8), mat2str(squeeze(P.TF).', 8), ...
    mat2str(squeeze(SC.DEARTHKM).', 6), mat2str(squeeze(SC.DMOONKM).', 6));
QC = load(matC);
assert(isequal(sort(fieldnames(QA)), sort(fieldnames(QC))), 'sheet variables differ with keepDirect on');
for kv = {'TF', 'FLYKM', 'ACCDZ', 'RES', 'OK'}
    assert(isequaln(QC.(kv{1})(:, :, 1), QA.(kv{1})), '15 N %s differs with keepDirect on', kv{1});
end
assert(isequaln(QC.Z8(:, :, :, 1), QA.Z8) && isequaln(QC.ATT, QA.ATT), '15 N Z8 / ATT differ with keepDirect on');
assert(isequaln(QC.meta, QA.meta), 'meta differs with keepDirect on');
w = whos('-file', sideC);
fprintf('  (b) sidecar %.1f kB on disk for 1 cell x 2 rungs (%d variables)\n', dir(sideC).bytes/1024, numel(w));

% (c) resume: the cell holds OK rungs -> skipped, slices untouched
o.resume = true;  o.logFile = fullfile(outDir, 'resume.log');
thrust_ladder_library(matC, o);
SR = load(sideC);
assert(isequaln(SR, SC), 'resume changed the sidecar of a skipped cell');
% a sidecar whose N does not match the call -> a fresh one, logged
SX = SC;  SX.N = 200;  save(sideC, '-struct', 'SX');
o.logFile = fullfile(outDir, 'resume2.log');
thrust_ladder_library(matC, o);
logR = fileread(o.logFile);
assert(contains(logR, 'direct sidecar') && contains(logR, 'fresh'), ...
    'no "fresh direct sidecar" line on a mismatched sidecar:\n%s', logR);
SF = load(sideC);
assert(SF.N == N && ~any(SF.DDONE(:)) && all(isnan(SF.DX(:))), 'mismatched sidecar was not replaced by a fresh one');
fprintf('  (c) resume: skipped cell untouched; mismatched sidecar replaced (log: %s)\n', ...
    strtrim(regexp(logR, '[^\n]*direct sidecar[^\n]*', 'match', 'once')));
save(sideC, '-struct', 'SC');                    % put the recorded sidecar back

% (d) redo clears BEFORE solving: a redo whose direct solve fails
o = base;  o.rungs = [15 12];  o.keepDirect = true;  o.resume = true;
o.cells = [1 1];  o.redo = true;  o.maxIter = 3;  o.logFile = fullfile(outDir, 'redo_fail.log');
thrust_ladder_library(matC, o);
SD = load(sideC);
assert(~any(SD.DDONE(:)) && all(isnan(SD.DX(:))) && all(isnan(SD.DU(:))) && all(isnan(SD.DTF(:))) && ...
       all(isnan(SD.DEARTHKM(:))) && all(isnan(SD.DMOONKM(:))), ...
    'a redo whose direct solve failed left recorded slices (DDONE %s)', mat2str(squeeze(SD.DDONE).'));
logD = fileread(o.logFile);
assert(contains(logD, 'direct fail') || contains(logD, 'direct ERROR'), 'the maxIter = 3 redo did not fail:\n%s', logD);
% ... then a redo seeded from the kept 15 N direct solution rewrites them
seedTop = struct('X', SC.DX(:, :, 1, 1, 1), 'U', [SC.DU(:, :, 1, 1, 1); ones(1, N+1)], ...
                 'tf', SC.DTF(1, 1, 1));
o.maxIter = 6000;  o.seedTop = seedTop;  o.logFile = fullfile(outDir, 'redo.log');
P2 = thrust_ladder_library(matC, o);
SE = load(sideC);
check_sidecar(SE, P2, rv0, N, 'on.mat');
dTF = max(abs(SE.DTF(:) - SC.DTF(:)));
fprintf('  (d) redo from the kept direct solution: DDONE %s, max|dDTF| vs (b) = %.3g ND\n', ...
    mat2str(squeeze(SE.DDONE).'), dTF);
assert(dTF < 1e-6, 'the redo seeded from the kept 15 N solution moved DTF by %.3g ND', dTF);

% (e) deliberate failure: the checker must trip on a corrupted DDONE
SBad = SE;  SBad.DDONE(1, 1, 2) = false;
tripped = false;
try
    check_sidecar(SBad, P2, rv0, N, 'on.mat');
catch err
    tripped = true;
    fprintf('  (e) corrupted DDONE trips the checker: %s\n', err.message);
end
assert(tripped, 'the sidecar checker did not trip on a corrupted DDONE');

ok = true;
fprintf('test_ladder_keepdirect: PASS\n');
end

% ------------------------------------------------------------------------
function check_sidecar(S, P, rv0, N, sheetName)
% CHECK_SIDECAR  Assert the sidecar contract for a 1 x 1 sheet.
% INPUTS:  S sidecar struct; P engine output; rv0 [6 x 1] departure state;
%          N collocation intervals; sheetName char.  OUTPUTS: none (throws).
nR = numel(P.rungs);
assert(isequal(size(S.DX, 1:5), [7 N+1 1 1 nR]) && isequal(size(S.DU, 1:5), [3 N+1 1 1 nR]), 'DX/DU size');
assert(isequal(size(S.DTF, 1:3), [1 1 nR]) && isequal(size(S.DDONE, 1:3), [1 1 nR]) && islogical(S.DDONE), 'DTF/DDONE size');
assert(isequal(S.rungs, P.rungs) && isequal(S.sD, P.sD) && isequal(S.sA, P.sA) && S.N == N, 'identity fields');
assert(strcmp(S.sheetFile, sheetName), 'sheetFile %s, expected %s', S.sheetFile, sheetName);
expect = isfinite(P.TF);                          % every rung that ran passed the direct gate on this cell
assert(isequal(S.DDONE, expect), 'DDONE %s but direct solves succeeded at %s', ...
    mat2str(S.DDONE(:).'), mat2str(expect(:).'));
mu = P.meta.muStar;  lStar = P.meta.lStar;
for kr = 1:nR
    if ~S.DDONE(1, 1, kr)
        assert(all(isnan(S.DX(:, :, 1, 1, kr)), 'all') && isnan(S.DTF(1, 1, kr)), 'unrecorded rung %d not NaN', kr);
        continue
    end
    X = S.DX(:, :, 1, 1, kr);  U = S.DU(:, :, 1, 1, kr);
    assert(isfinite(S.DTF(1, 1, kr)) && S.DTF(1, 1, kr) > 0, 'DTF rung %d', kr);
    dx0 = max(abs(X(1:6, 1) - rv0(1:6)));
    assert(dx0 < 1e-9, 'DX(1:6,1) rung %d off the departure state by %.3g', kr, dx0);
    dU = max(abs(sqrt(sum(U.^2, 1)) - 1));
    assert(dU < 1e-6, 'DU columns rung %d not unit norm (%.3g)', kr, dU);
    eNode = min(sqrt(sum((X(1:3, :) - [-mu; 0; 0]).^2, 1)))*lStar - 6378.137;
    mNode = min(sqrt(sum((X(1:3, :) - [1-mu; 0; 0]).^2, 1)))*lStar - 1737.4;
    assert(S.DEARTHKM(1, 1, kr) <= eNode + 1e-6 && S.DEARTHKM(1, 1, kr) > eNode - 2000, ...
        'DEARTHKM %.6f vs node-only %.6f km, rung %d', S.DEARTHKM(1, 1, kr), eNode, kr);
    assert(S.DMOONKM(1, 1, kr) <= mNode + 1e-6 && S.DMOONKM(1, 1, kr) > mNode - 2000, ...
        'DMOONKM %.6f vs node-only %.6f km, rung %d', S.DMOONKM(1, 1, kr), mNode, kr);
end
end

function rmdir_quiet(d)
% RMDIR_QUIET  Remove a temp folder if present.  INPUTS: d char.  OUTPUTS: none.
if isfolder(d)
    try rmdir(d, 's'); catch, end
end
end
