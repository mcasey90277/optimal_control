function ok = test_dro_elfo_pick()
%% Purpose:
%
%   Test of dro_elfo_pick (Task 7) on a synthetic keyed catalog (no solve,
%   seconds): sheets (depSize 1, raanDeg 0 / 10 / 90) and (depSize 2,
%   raanDeg 0), 2 x 3 phases, 2 rungs [15 12] N, every z8 distinct.
%     1  refusals by name (deliberate failures of a nearest-sheet pick): a
%        raanDeg between two sheets (45) errors dro_elfo_pick:noSheet and
%        lists the available keys; a depSize off by 1e-6; a thrust that is
%        not a rung (:noRung); a rung with no entry on the sheet
%        (:noEntry); two sheets with one key (:dupSheet); a catalog
%        without sheet_key (:catalog)
%     2  angle wrapping: 370 and -350 deg pick the 10 deg sheet; a depSize
%        within 1e-12 picks
%     3  nearest cell on the torus (wrapped phases): sD 0.95, sA 0.70 ->
%        cell (1, 3) (sD 0, sA 2/3), with its own phases and distance; a
%        nearest cell with no entry at the rung is skipped for the next
%        nearest that has one, with e.nearestEmpty and the warning
%        dro_elfo_pick:nearestEmpty naming both cells; a sheet of another
%        target with the same depSize / raanDeg is not matched
%     4  verdicts returned (conj_pass 0 on a refuted entry, gate values);
%        -1 / NaN when the sheet carries none; the catalog loaded from a
%        .mat path gives the same entry
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
addpath(fullfile(fileparts(here), 'catalog'));

C = synth_catalog([1 0; 1 10; 1 90; 2 0]);

%% 1. refusals by name
expect_err(@() dro_elfo_pick(C, 1, 45, 0, 0, 15), 'dro_elfo_pick:noSheet', 'raanDeg 10');
expect_err(@() dro_elfo_pick(C, 1 + 1e-6, 0, 0, 0, 15), 'dro_elfo_pick:noSheet', 'depSize 1.000001');
expect_err(@() dro_elfo_pick(C, 1, 0, 0, 0, 11), 'dro_elfo_pick:noRung', '15');
Ce = C;  Ce.sheets(1).has_solution(:, :, 2) = false;
expect_err(@() dro_elfo_pick(Ce, 1, 0, 0, 0, 12), 'dro_elfo_pick:noEntry', '12 N');
Cd = C;  Cd.sheets(2).sheet_key.raanDeg = 0;
expect_err(@() dro_elfo_pick(Cd, 1, 0, 0, 0, 15), 'dro_elfo_pick:dupSheet', '');
Ct = C;  Ct.sheets(1).sheet_key.target = 'tulip';        % same depSize/raanDeg, other target
expect_err(@() dro_elfo_pick(Ct, 1, 0, 0, 0, 15), 'dro_elfo_pick:noSheet', 'target tulip, depSize 1, raanDeg 0');
Cn = rmfield(C, 'sheets');  Cn.sheets = rmfield(C.sheets, 'sheet_key');
expect_err(@() dro_elfo_pick(Cn, 1, 0, 0, 0, 15), 'dro_elfo_pick:catalog', '');

%% 2. wrapping
e370 = dro_elfo_pick(C, 1, 370, 0, 0, 15);
eNeg = dro_elfo_pick(C, 1, -350, 0, 0, 15);
assert(e370.sheet_key.raanDeg == 10 && e370.sheet == 2, '370 deg picks the 10 deg sheet');
assert(isequal(eNeg.z8, e370.z8), '-350 deg picks the 10 deg sheet');
eTol = dro_elfo_pick(C, 1 + 1e-12, 90, 0, 0, 15);
assert(eTol.sheet_key.raanDeg == 90, 'depSize within 1e-12 picks');

%% 3. nearest cell
e = dro_elfo_pick(C, 1, 0, 0.95, 0.70, 12);
assert(e.iD == 1 && e.iA == 3 && e.kr == 2 && abs(e.sA - 2/3) < 1e-15 && e.sD == 0, ...
       'nearest cell (1,3) on the torus (got (%d,%d))', e.iD, e.iA);
assert(abs(e.dPhase - sqrt(0.05^2 + (0.70 - 2/3)^2)) < 1e-12, 'distance on the torus');
assert(~e.nearestEmpty, 'nearest cell used: nearestEmpty false');
assert(isequal(e.z8, C.sheets(1).z8(:, C.sheets(1).entry_index(1, 3, 2))) && e.tf_nd == e.z8(8), 'the cell''s own z8');
Cs = C;  Cs.sheets(1).has_solution(1, 3, 2) = false;
Cs.sheets(1).entry_index(1, 3, 2) = 0;
ws = warning('off', 'dro_elfo_pick:nearestEmpty');  lastwarn('');
eS = dro_elfo_pick(Cs, 1, 0, 0.95, 0.70, 12);
[wMsg, wId] = lastwarn;  warning(ws);
assert(~(eS.iD == 1 && eS.iA == 3) && Cs.sheets(1).has_solution(eS.iD, eS.iA, 2), 'an empty nearest cell is skipped');
assert(eS.nearestEmpty && strcmp(wId, 'dro_elfo_pick:nearestEmpty') && contains(wMsg, '(1,3)') && ...
       contains(wMsg, sprintf('(%d,%d)', eS.iD, eS.iA)), 'nearestEmpty flag and named warning (got "%s")', wMsg);
fprintf('  warned as expected: %s\n', wMsg);

%% 4. verdicts, path input
eR = dro_elfo_pick(C, 2, 0, 0.5, 1/3, 15);
assert(eR.conj_pass == 0 && eR.gate_dimS == 1 && abs(eR.gate_min_lamv - 0.25) < 1e-15, 'verdicts returned');
Cv = C;  Cv.sheets = rmfield(C.sheets, {'conj_pass', 'gate_min_lamv', 'gate_min_qmt', 'gate_dimS'});
eV = dro_elfo_pick(Cv, 1, 0, 0, 0, 15);
assert(eV.conj_pass == -1 && isnan(eV.gate_min_lamv), 'no verdicts: -1 / NaN');
catMat = [tempname '.mat'];
S = struct('costate_catalog_dro_elfo', C);  save(catMat, '-struct', 'S');
eP = dro_elfo_pick(catMat, 1, 370, 0.95, 0.70, 12);
delete(catMat);
assert(isequal(eP.z8, dro_elfo_pick(C, 1, 10, 0.95, 0.70, 12).z8), 'path input gives the same entry');

ok = true;
fprintf('test_dro_elfo_pick: PASS\n');
end

% ------------------------------------------------------------------------
function C = synth_catalog(keys)
% SYNTH_CATALOG  A keyed catalog with one 2 x 3 x 2 sheet per key row; z8
% columns distinct (tf = sheet + cell/100); every entry passes but the
% sheet-4 entry (2,2,1), refuted.  INPUTS: keys [k x 2] (depSize, raanDeg).
% OUTPUTS: C [1 x 1] catalog struct (.rungs_N [1 x 2], .sheets [1 x k]).
C = struct('rungs_N', [15 12]);
nD = 2;  nA = 3;  nR = 2;
for ks = 1:size(keys, 1)
    has = true(nD, nA, nR);
    idx = reshape(1:nD*nA*nR, nD, nA, nR);
    z8  = [repmat((1:7).', 1, nD*nA*nR); ks + (1:nD*nA*nR)/100];
    cp  = ones(nD, nA, nR, 'int8');
    if ks == 4, cp(2, 2, 1) = 0; end
    C.sheets(ks) = struct('sheet_key', struct('target', 'elfo', 'depSize', keys(ks, 1), 'raanDeg', keys(ks, 2)), ...
        'sD_frac', [0 0.5], 'sA_frac', [0 1/3 2/3], 'has_solution', has, 'entry_index', idx, ...
        'z8', z8, 'conj_pass', cp, 'gate_min_lamv', 0.25*ones(nD, nA, nR), ...
        'gate_min_qmt', ones(nD, nA, nR), 'gate_dimS', ones(nD, nA, nR, 'int8'), ...
        'dep_family', 'dro', 'dep_params', struct('tau', keys(ks, 1)), 'arr_family', 'elfo', ...
        'arr_params', struct('oe', [12000 0.69 56.5 90], 'raanDeg', keys(ks, 2)));
end
end

function expect_err(fh, id, msgPart)
% EXPECT_ERR  fh must throw identifier id, its message containing msgPart.
% INPUTS: fh [1 x 1] function handle; id [1 x n] char; msgPart [1 x m] char.
% OUTPUTS: none (throws on a mismatch).
try
    fh();
catch err
    assert(strcmp(err.identifier, id), 'expected %s, got %s (%s)', id, err.identifier, err.message);
    assert(contains(err.message, msgPart), 'message of %s lacks "%s": %s', id, msgPart, err.message);
    fprintf('  refused as expected: %s\n', err.message);
    return
end
error('test_dro_elfo_pick:noThrow', 'expected %s, nothing thrown', id);
end
