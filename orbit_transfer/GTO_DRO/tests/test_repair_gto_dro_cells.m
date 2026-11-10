function ok = test_repair_gto_dro_cells()
% TEST_REPAIR_GTO_DRO_CELLS  Basin repair on the demo cell, restore on a
% failed repair, and the unsolved-mode neighbour search.
%
%   1. conjOnly: the fixture's 15 N entry is the refuted basin (stage 1)
%   2. basin repair: reverse warm start from the 12 N entry lands on the
%      faster basin; conjAfter is PRINTED (an open measurement, not asserted)
%   3. a second basin call repairs nothing (cell already in the sidecar)
%   4. restore on failure, on a FRESH copy of the fixture (refuted entry,
%      no sidecar) with a garbage seed: OK, TF, Z8, ATT identical to before
%   5. unsolved mode on a 2 x 1 sheet; the second cell sits at phase 0.60
%      (the engine owns the grid), which stage 1 did not measure: if it does
%      not solve at 15 N that is a measurement, reported by the printout
%
%   Fix round 1 adds: 6 engineOpts physics refusal (sheet untouched),
%   7 engine-throw path (wrong-shape seed -> restored (error), no backup
%   left), 8 interrupted-repair recovery from a hand-made backup + cleared
%   cell, 9 fleet-lock refusal, 10 conjOnly on mode unsolved refuses
%
% INPUTS:  none
% OUTPUTS: ok - logical
%
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
addpath(here);  addpath(fullfile(fileparts(here), 'catalog'));
src = fixture_demo_sheet();
wd  = fullfile(tempdir, 'test_repair_gto_dro');
if isfolder(wd), rmdir(wd, 's'); end
mkdir(wd);
sheetMat = fullfile(wd, 'demo_sheet.mat');
freshMat = fullfile(wd, 'fresh_sheet.mat');
copyfile(src, sheetMat);
copyfile(src, freshMat);
plan = run_gto_dro_catalog(2, [], [], [15 12], struct('planOnly', true));
eo = plan.optsA;  eo.logFile = fullfile(wd, 'repair.log');

% 1. conjOnly
R0 = repair_gto_dro_cells(sheetMat, 'basin', struct('engineOpts', eo, 'conjOnly', true));
assert(R0.conjPass(1, 1) == 0, 'the 15 N cold-start entry must be refuted (got %g)', R0.conjPass(1, 1));
assert(R0.nConjRefuted == 1 && R0.nRepaired == 0, 'conjOnly repairs nothing');
Q = load(sheetMat, 'TF', 'OK');
tfOld = Q.TF(1, 1, 1);
fprintf('  conjOnly: 15 N entry t_f %.5f, conjPass %g\n', tfOld, R0.conjPass(1, 1));

% 2. basin repair
R1 = repair_gto_dro_cells(sheetMat, 'basin', struct('engineOpts', eo));
assert(R1.nRepaired == 1 && R1.nRestored == 0, 'one cell repaired');
Q = load(sheetMat, 'TF', 'OK', 'ATT');
assert(Q.OK(1, 1, 1) && Q.TF(1, 1, 1) < 0.29 && Q.TF(1, 1, 1) < tfOld, 'repaired 15 N on the faster basin');
assert(Q.OK(1, 1, 2), 'the 12 N rung re-walked and OK');
S = load(R1.sideMat, 'rows');
assert(numel(S.rows) == 1 && strcmp(S.rows(1).outcome, 'repaired') && S.rows(1).conjBefore == 0, 'sidecar row');
assert(S.rows(1).directGateKm == 1000, 'sidecar records directGateKm');
fprintf('  basin repair: tf 15 N %.5f -> %.5f (12 N %.5f), conjBefore %g, conjAfter %g, %.0f s\n', ...
        tfOld, Q.TF(1, 1, 1), Q.TF(1, 1, 2), S.rows(1).conjBefore, S.rows(1).conjAfter, S.rows(1).wallS);

% 3. a second basin call does nothing (cell already repaired)
R2 = repair_gto_dro_cells(sheetMat, 'basin', struct('engineOpts', eo));
assert(R2.nCandidates == 0 || R2.nRepaired == 0, 'no second repair without force');
assert(R2.nRepaired == 0, 'no second repair without force');
assert(R2.conjPass(1, 1) == 1, 'cached conjAfter verdict reused');

% 4. restore on failure, fresh copy: garbage seed (t_f 1e-4 ND cannot reach the DRO)
Qb = load(freshMat);
seedBad = struct('X', repmat([0.1*ones(6, 1); 1], 1, 401), ...
                 'U', repmat([1; 0; 0; 1], 1, 401), 'tf', 1e-4);
eo4 = eo;  eo4.maxCpuSec = 60;
R3 = repair_gto_dro_cells(freshMat, 'basin', struct('engineOpts', eo4, 'seedOverride', seedBad));
assert(R3.nCandidates == 1, 'the fresh refuted entry must be a candidate (got %d)', R3.nCandidates);
assert(R3.nRestored == 1, 'a failed repair must restore');
Qa = load(freshMat);
assert(isequal(Qa.OK, Qb.OK), 'restored OK');
assert(max(abs(Qa.TF(:) - Qb.TF(:))) == 0, 'restored TF');
assert(isequaln(Qa.Z8, Qb.Z8), 'restored Z8');
assert(isequal(Qa.ATT, Qb.ATT), 'restored ATT (%g vs %g)', Qa.ATT, Qb.ATT);
S3 = load(R3.sideMat, 'rows');
assert(strncmp(S3.rows(end).outcome, 'restored', 8), 'sidecar outcome must start with restored (got %s)', S3.rows(end).outcome);
assert(isempty(dir([freshMat(1:end-4) '_repair_backup_*.mat'])), 'backup file removed after the row is saved');
assert(~isfolder(fullfile(tempdir, '.gtodrocat_driver.lock')), 'the tool released the lock it took');
snap = load([freshMat '.bak_repair'], '-mat');
assert(isequal(snap.OK, Qb.OK) && isequaln(snap.TF, Qb.TF) && isequaln(snap.Z8, Qb.Z8) && isequal(snap.ATT, Qb.ATT), ...
       'snapshot loads and equals the pre-call sheet');
fprintf('  failed repair: outcome "%s", cell restored\n', S3.rows(end).outcome);

% 6. physics guard: engineOpts.ispS differs from the sheet -> refuse, sheet untouched
Qpre = load(freshMat);
eoBad = eo;  eoBad.ispS = 900;
try
    repair_gto_dro_cells(freshMat, 'basin', struct('engineOpts', eoBad, 'conjOnly', true));
    error('test:noThrow', 'a physics mismatch did not refuse');
catch err
    assert(strcmp(err.identifier, 'repair_gto_dro_cells:engineOpts'), 'wrong error: %s', err.identifier);
    fprintf('  physics guard: %s\n', err.message);
end
assert(isequaln(load(freshMat), Qpre), 'refused call must leave the sheet untouched');
eoBad = eo;  eoBad.thrLock = false;
try
    repair_gto_dro_cells(freshMat, 'basin', struct('engineOpts', eoBad, 'conjOnly', true));
    error('test:noThrow', 'thrLock mismatch did not refuse');
catch err
    assert(strcmp(err.identifier, 'repair_gto_dro_cells:engineOpts'), 'wrong error: %s', err.identifier);
    fprintf('  physics guard: %s\n', err.message);
end
eoBad = eo;  eoBad.accTol = 1e-4;
try
    repair_gto_dro_cells(freshMat, 'basin', struct('engineOpts', eoBad, 'conjOnly', true));
    error('test:noThrow', 'accTol mismatch did not refuse');
catch err
    assert(strcmp(err.identifier, 'repair_gto_dro_cells:engineOpts'), 'wrong error: %s', err.identifier);
    fprintf('  physics guard: %s\n', err.message);
end
assert(isequaln(load(freshMat), Qpre), 'refused calls must leave the sheet untouched');

% 7. the engine throws (wrong-shape seed): restored (error ...), nFailed 1
errMat = fullfile(wd, 'err_sheet.mat');  copyfile(src, errMat);
Qe0 = load(errMat);
seedShape = struct('X', ones(7, 10), 'U', ones(4, 10), 'tf', 0.3);
R5 = repair_gto_dro_cells(errMat, 'basin', struct('engineOpts', eo, 'seedOverride', seedShape));
assert(R5.nCandidates == 1 && R5.nFailed == 1 && R5.nRestored == 1, 'engine throw counted');
Qe1 = load(errMat);
assert(isequal(Qe1.OK, Qe0.OK) && isequaln(Qe1.TF, Qe0.TF) && isequaln(Qe1.Z8, Qe0.Z8) && ...
       isequal(Qe1.ATT, Qe0.ATT), 'cell restored after an engine throw');
S5 = load(R5.sideMat, 'rows');
assert(strncmp(S5.rows(end).outcome, 'restored (error', 15), 'outcome %s', S5.rows(end).outcome);
assert(isempty(dir([errMat(1:end-4) '_repair_backup_*.mat'])), 'no backup left after an engine throw');
fprintf('  engine throw: outcome "%s"\n', S5.rows(end).outcome);

% 8. interrupted-repair recovery: cleared cell + pending backup, conjOnly call
catDir = fullfile(wd, 'catalog');  mkdir(catDir);
recMat = fullfile(catDir, 'rec_sheet.mat');  copyfile(src, recMat);
Qr = load(recMat);
bkR = struct('TF', squeeze(Qr.TF(1, 1, :)), 'FLYKM', squeeze(Qr.FLYKM(1, 1, :)), ...
             'ACCDZ', squeeze(Qr.ACCDZ(1, 1, :)), 'RES', squeeze(Qr.RES(1, 1, :)), ...
             'WALL', squeeze(Qr.WALL(1, 1, :)), 'OK', squeeze(Qr.OK(1, 1, :)), ...
             'Z8', squeeze(Qr.Z8(:, 1, 1, :)), 'ATT', Qr.ATT(1, 1));
rowR = struct('mode', 'basin', 'iD', 1, 'iA', 1, 'seedFrom', struct('iD', 1, 'iA', 1, 'kr', 2), ...
              'tfBefore', bkR.TF(1), 'tfAfter', NaN, 'okAfter', false, 'outcome', '', ...
              'conjBefore', 0, 'conjAfter', NaN, 'date', 'x', 'wallS', NaN, 'directGateKm', 1000);
bfileR = [recMat(1:end-4) '_repair_backup_i01_a01.mat'];
Bsave = struct('bk', bkR, 'mode', 'basin', 'iD', 1, 'iA', 1, 'seedFrom', rowR.seedFrom, ...
               'date', 'x', 'row', rowR);
save(bfileR, '-struct', 'Bsave');
Qc = Qr;                                              % the engine's cleared cell
Qc.TF(1, 1, :) = NaN;  Qc.FLYKM(1, 1, :) = NaN;  Qc.ACCDZ(1, 1, :) = NaN;
Qc.RES(1, 1, :) = NaN;  Qc.WALL(1, 1, :) = NaN;  Qc.OK(1, 1, :) = false;
Qc.Z8(:, 1, 1, :) = NaN;  Qc.ATT(1, 1) = 1;
save(recMat, '-struct', 'Qc');

% 9. fleet lock: refuses while a pending backup needs recovery (even conjOnly) ...
lockDir = fullfile(wd, '.gtodrocat_driver.lock');  mkdir(lockDir);
try
    repair_gto_dro_cells(recMat, 'basin', struct('engineOpts', eo, 'conjOnly', true));
    error('test:noThrow', 'the fleet lock did not refuse');
catch err
    assert(strcmp(err.identifier, 'repair_gto_dro_cells:fleetRunning'), 'wrong error: %s', err.identifier);
end
assert(isfile(bfileR), 'a refused call keeps the backup');
% ... and a repair call refuses under the lock
try
    repair_gto_dro_cells(recMat, 'basin', struct('engineOpts', eo));
    error('test:noThrow', 'the fleet lock did not refuse a repair');
catch err
    assert(strcmp(err.identifier, 'repair_gto_dro_cells:fleetRunning'), 'wrong error: %s', err.identifier);
end
assert(isfolder(lockDir), 'a foreign lock must not be removed');
rmdir(lockDir);

% 8 (continued): now recover
R6 = repair_gto_dro_cells(recMat, 'basin', struct('engineOpts', eo, 'conjOnly', true));
assert(R6.nRecovered == 1, 'one interrupted repair recovered');
assert(R6.nCandidates == 1, 'an interrupted repair is a candidate again (got %d)', R6.nCandidates);
Qd = load(recMat);
assert(isequal(Qd.OK, Qr.OK) && isequaln(Qd.TF, Qr.TF) && isequaln(Qd.Z8, Qr.Z8) && isequal(Qd.ATT, Qr.ATT), ...
       'recovered cell equals the original');
S6 = load(R6.sideMat, 'rows');
assert(strcmp(S6.rows(1).outcome, 'restored (interrupted)') && S6.rows(1).tfBefore == bkR.TF(1), 'interrupted row');
assert(~isfile(bfileR), 'backup deleted after recovery');
fprintf('  recovery: nRecovered %d, tf restored %.5f\n', R6.nRecovered, Qd.TF(1, 1, 1));

% conjOnly under the lock is allowed once nothing is pending
mkdir(lockDir);
R7 = repair_gto_dro_cells(recMat, 'basin', struct('engineOpts', eo, 'conjOnly', true));
assert(R7.nRecovered == 0 && R7.conjPass(1, 1) == 0, 'conjOnly runs under the lock');
rmdir(lockDir);

% 10. conjOnly on mode unsolved refuses
try
    repair_gto_dro_cells(recMat, 'unsolved', struct('engineOpts', eo, 'conjOnly', true));
    error('test:noThrow', 'conjOnly unsolved did not refuse');
catch err
    assert(strcmp(err.identifier, 'repair_gto_dro_cells:mode'), 'wrong error: %s', err.identifier);
end

% 5. unsolved mode on a 2 x 1 sheet: cell 2 empty and retired, seeded from cell 1
S2 = load(sheetMat);
S2.sD = mod(0.10 + (0:1)/2, 1);
f = {'TF', 'FLYKM', 'ACCDZ', 'RES', 'WALL'};
for kf = 1:numel(f), S2.(f{kf}) = cat(1, S2.(f{kf}), nan(1, 1, 2)); end
S2.OK  = cat(1, S2.OK, false(1, 1, 2));
S2.Z8  = cat(2, S2.Z8, nan(8, 1, 1, 2));
S2.ATT = [S2.ATT; 2];
sheet2 = fullfile(wd, 'two_cells.mat');
save(sheet2, '-struct', 'S2');
eo2 = eo;  eo2.logFile = fullfile(wd, 'repair2.log');
R4 = repair_gto_dro_cells(sheet2, 'unsolved', struct('engineOpts', eo2));
Q2 = load(sheet2, 'OK', 'TF', 'ATT');
assert(R4.nCandidates == 1, 'one unsolved candidate (got %d)', R4.nCandidates);
S4 = load(R4.sideMat, 'rows');
assert(numel(S4.rows) == 1 && S4.rows(end).seedFrom.iD == 1, 'seeded from cell 1');
fprintf('  phase continuation: cell 2 (sD 0.60) outcome "%s", OK = %d, tf 15 N = %.5f, ATT %d, %.0f s\n', ...
        S4.rows(end).outcome, Q2.OK(2, 1, 1), Q2.TF(2, 1, 1), Q2.ATT(2, 1), S4.rows(end).wallS);
assert(Q2.OK(2, 1, 1), 'phase continuation solved the neighbour cell at 15 N');

rmdir(wd, 's');
ok = true;
fprintf('test_repair_gto_dro_cells: PASS\n');
end
