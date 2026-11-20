function ok = test_repair_dro_elfo_cells(topts)
%% Purpose:
%
%   Test of repair_dro_elfo_cells (Task 7). Sections 1-8 use the test seams
%   (engineFcn: a fake engine that rewrites the cell the way
%   thrust_ladder_library does, sheet_key DROPPED as the real engine drops
%   it; conjFcn: a verdict by z8) on 1 x 1 x 3 cuts of the matched-
%   experiment sheet, seconds:
%     1  ASTRA C5: a refuted top rung, re-solved faster, the re-walk LOST
%        rung 2 -> accepted; rung 1 is the re-solve's, rungs 2 and 3 keep
%        their previous entries (keptOld [2 3], lostNew [2]); sheet_key back
%        on the sheet; the snapshot <sheet>.bak_repair beside it equals the
%        pre-call sheet; no cell backup left; conjAfter from the verdict
%     2  an empty lower rung takes the re-solve's entry when it is gated,
%        and stays empty when the re-solve's entry fails the campaign gate
%     3  deliberate failures of the acceptance rule: a re-solve that is not
%        faster, and one whose top rung fails the campaign gate, are
%        RESTORED (sheet identical to before, key kept)
%     3b CLEARANCE before acceptance (clearFcn seam): a faster, gated top
%        rung whose located minimum is 100 km (Moon) or 150 km (Earth) is
%        'restored (clearance)', the stored entry stays; a lower-rung
%        re-solve entry below the floor does not fill an empty rung
%     4  the engine throws -> 'restored (error ...)', nFailed 1, key kept
%     5  refusals by name: the fleet lock (sheet folder
%        .dro_elfo_driver.lock; a foreign lock survives), a physics
%        mismatch on the ELFO raanDeg (sheet untouched), a sheet with no
%        recoverable sheet_key, a sheet_key that disagrees with meta
%     6  an interrupted repair (cleared cell, keyless sheet, pending cell
%        backup carrying the key) is recovered: cell and key restored
%     7  REAL ENGINE (topts.solve, default true; a few minutes): 'unsolved'
%        on a 6 x 1 x 2 cut (departure phases 0:1/6:5/6, arrival 1/6) with
%        cell 3 emptied and retired: it is seeded from cell 2 and recovered
%        at 15 N; sheet_key on the sheet afterwards
%
%% Inputs:
%
%  topts                    struct (optional)       .solve [true] run
%                                                   section 7
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
addpath(here);
wd = tempname;  mkdir(wd);
cleanWd = onCleanup(@() rmdir(wd, 's'));

plan = run_dro_elfo_torus(1, struct('planOnly', true));
eo = plan.optsA;  eo.msHardCapSec = 0;  eo.keepDirect = false;
cutSpec = struct('iD', 2, 'iA', 2, 'nR', 3);
newSheet = @(name) fixture_dro_elfo_sheet(fullfile(wd, name, 'dro_elfo_tau1_raan000.mat'), cutSpec);

%% 1. astra C5
s1 = newSheet('c5');  Q0 = load(s1);
zTop0 = Q0.Z8(:, 1, 1, 1);  tf0 = Q0.TF(1, 1, 1);
spec = struct('tf', [tf0 - 1e-3, NaN, Q0.TF(1, 1, 3) - 5e-4], 'ok', [true false true], 'fly', [1 1 1]);
o1 = struct('engineOpts', eo, 'engineFcn', @(s, e) fake_engine(s, e, spec), ...
            'conjFcn', @(z8, varargin) double(~isequal(z8, zTop0)), ...
            'clearFcn', @(varargin) [1e4 1e5]);
R1 = repair_dro_elfo_cells(s1, 'basin', o1);
Q1 = load(s1);
assert(R1.nConjRefuted == 1 && R1.nCandidates == 1 && R1.nRepaired == 1 && R1.nRestored == 0, 'C5: accepted');
assert(abs(Q1.TF(1, 1, 1) - (tf0 - 1e-3)) < 1e-15 && Q1.OK(1, 1, 1), 'C5: the new top rung');
assert(Q1.OK(1, 1, 2) && isequal(Q1.Z8(:, 1, 1, 2), Q0.Z8(:, 1, 1, 2)) && Q1.TF(1, 1, 2) == Q0.TF(1, 1, 2), ...
       'C5: the lost rung 2 keeps its previous entry');
assert(isequal(Q1.Z8(:, 1, 1, 3), Q0.Z8(:, 1, 1, 3)), 'C5: rung 3 keeps its previous entry');
Sd = load(R1.sideMat);
assert(isequal(Sd.rows(end).keptOld, [2 3]) && isequal(Sd.rows(end).lostNew, 2) && ...
       strcmp(Sd.rows(end).outcome, 'repaired') && Sd.rows(end).conjAfter == 1 && R1.nKeptOld == 2, 'C5: sidecar row');
assert(isequaln(Q1.sheet_key, Q0.sheet_key), 'sheet_key restored after the engine dropped it');
snap = load([s1 '.bak_repair'], '-mat');
assert(isequaln(snap.Z8, Q0.Z8) && isequal(snap.OK, Q0.OK) && isfield(snap, 'sheet_key'), 'snapshot = pre-call sheet');
assert(isempty(dir([s1(1:end-4) '_repair_backup_*.mat'])), 'no cell backup left');
assert(~isfolder(fullfile(fileparts(s1), '.dro_elfo_driver.lock')), 'lock released');
R1b = repair_dro_elfo_cells(s1, 'basin', o1);
assert(R1b.nCandidates == 0 && R1b.conjPass(1, 1) == 1, 'second call: cached verdict, no candidate');
fprintf('  C5: tf %.6f -> %.6f, rungs kept from before [%s], lost by the re-walk [%s]\n', tf0, Q1.TF(1, 1, 1), ...
        num2str(Sd.rows(end).keptOld), num2str(Sd.rows(end).lostNew));

%% 2. empty lower rung: filled when gated, emptied when not
s2 = newSheet('fill');  Q = load(s2);
Q.OK(1, 1, 3) = false;  Q.TF(1, 1, 3) = NaN;  Q.Z8(:, 1, 1, 3) = NaN;  save(s2, '-struct', 'Q');
spec2 = struct('tf', [tf0 - 1e-3, Q0.TF(1, 1, 2), Q0.TF(1, 1, 3) + 1e-3], 'ok', [true true true], 'fly', [1 1 1]);
repair_dro_elfo_cells(s2, 'basin', with(o1, 'engineFcn', @(s, e) fake_engine(s, e, spec2)));
Q2 = load(s2);
assert(Q2.OK(1, 1, 3) && abs(Q2.TF(1, 1, 3) - spec2.tf(3)) < 1e-15, 'empty rung 3 takes the gated re-solve entry');
assert(isequal(Q2.Z8(:, 1, 1, 2), Q0.Z8(:, 1, 1, 2)), 'OK rung 2 keeps its previous entry over the re-solve''s');
s2b = newSheet('fillgate');  Q = load(s2b);
Q.OK(1, 1, 3) = false;  Q.TF(1, 1, 3) = NaN;  Q.Z8(:, 1, 1, 3) = NaN;  save(s2b, '-struct', 'Q');
spec2b = spec2;  spec2b.fly = [1 1 500];
repair_dro_elfo_cells(s2b, 'basin', with(o1, 'engineFcn', @(s, e) fake_engine(s, e, spec2b)));
Q2b = load(s2b);
assert(~Q2b.OK(1, 1, 3) && all(isnan(Q2b.Z8(:, 1, 1, 3))), 'an ungated re-solve entry does not fill rung 3');

%% 3. not faster / top rung ungated -> restored
for cs = {struct('tf', [tf0 + 1e-3, NaN, NaN], 'ok', [true false false], 'fly', [1 1 1], 'why', 'restored (not faster)'), ...
          struct('tf', [tf0 - 1e-3, NaN, NaN], 'ok', [true false false], 'fly', [500 1 1], 'why', 'restored (campaign gate)')}
    s3 = newSheet(['rest' num2str(cs{1}.fly(1))]);  Qb = load(s3);
    R3 = repair_dro_elfo_cells(s3, 'basin', with(o1, 'engineFcn', @(s, e) fake_engine(s, e, cs{1})));
    Qa = load(s3);  S3 = load(R3.sideMat);
    assert(R3.nRestored == 1 && R3.nRepaired == 0 && strcmp(S3.rows(end).outcome, cs{1}.why), ...
           'expected "%s", got "%s"', cs{1}.why, S3.rows(end).outcome);
    assert(isequal(Qa.OK, Qb.OK) && isequaln(Qa.TF, Qb.TF) && isequaln(Qa.Z8, Qb.Z8) && isequal(Qa.ATT, Qb.ATT) && ...
           isequaln(Qa.sheet_key, Qb.sheet_key), 'restored sheet identical, key kept');
end

%% 3b. clearance gate: a faster, gated re-solve that dips below the floor
% never replaces the stored entry; a lower-rung re-solve entry below the
% floor is not taken
s3c = newSheet('clear');  Qb = load(s3c);
specC = struct('tf', [tf0 - 1e-3, NaN, NaN], 'ok', [true false false], 'fly', [1 1 1]);
oC = with(with(o1, 'engineFcn', @(s, e) fake_engine(s, e, specC)), 'clearFcn', @(varargin) [100 1e5]);
RC = repair_dro_elfo_cells(s3c, 'basin', oC);
Qa = load(s3c);  SC = load(RC.sideMat);
assert(RC.nRestored == 1 && RC.nRepaired == 0 && strcmp(SC.rows(end).outcome, 'restored (clearance)') && ...
       isequal(SC.rows(end).clearKm, [100 1e5]), 'below-floor re-solve restored (got "%s")', SC.rows(end).outcome);
assert(isequaln(Qa.Z8, Qb.Z8) && isequal(Qa.OK, Qb.OK) && isequaln(Qa.sheet_key, Qb.sheet_key), ...
       'the stored entry that cleared the floor stays');
oE = with(o1, 'clearFcn', @(varargin) [1e4 150]);                 % Earth floor 200 km
s3e = newSheet('clearE');
RE = repair_dro_elfo_cells(s3e, 'basin', with(oE, 'engineFcn', @(s, e) fake_engine(s, e, specC)));
assert(RE.nRestored == 1, 'below the Earth floor: restored');
s3l = newSheet('clearLow');  Q = load(s3l);
Q.OK(1, 1, 3) = false;  Q.TF(1, 1, 3) = NaN;  Q.Z8(:, 1, 1, 3) = NaN;  save(s3l, '-struct', 'Q');
oL = with(with(o1, 'engineFcn', @(s, e) fake_engine(s, e, spec2)), ...
          'clearFcn', @(z8, varargin) [100 + 1e4*(abs(z8(1) - 0.3) > 1e-12), 1e5]);   % rung 3 dips
RL = repair_dro_elfo_cells(s3l, 'basin', oL);
QL = load(s3l);
assert(RL.nRepaired == 1 && ~QL.OK(1, 1, 3), 'a lower-rung re-solve entry below the floor is not taken');
fprintf('  clearance gate: top rung at 100 km -> restored; Earth 150 km -> restored; rung 3 at 100 km not taken\n');

%% 4. engine throws
s4 = newSheet('throw');  Qb = load(s4);
R4 = repair_dro_elfo_cells(s4, 'basin', with(o1, 'engineFcn', @(varargin) error('t:x', 'boom')));
Qa = load(s4);  S4 = load(R4.sideMat);
assert(R4.nFailed == 1 && strncmp(S4.rows(end).outcome, 'restored (error', 15), 'engine throw counted');
assert(isequaln(Qa.Z8, Qb.Z8) && isequaln(Qa.sheet_key, Qb.sheet_key), 'restored, key kept');
assert(isempty(dir([s4(1:end-4) '_repair_backup_*.mat'])), 'no backup left after a throw');

%% 5. refusals
s5 = newSheet('refuse');  Qb = load(s5);
lockDir = fullfile(fileparts(s5), '.dro_elfo_driver.lock');  mkdir(lockDir);
expect_err(@() repair_dro_elfo_cells(s5, 'basin', o1), 'repair_dro_elfo_cells:fleetRunning');
assert(isfolder(lockDir), 'a foreign lock survives');  rmdir(lockDir);
eoBad = eo;  eoBad.arrParams.raanDeg = 30;
expect_err(@() repair_dro_elfo_cells(s5, 'basin', with(o1, 'engineOpts', eoBad)), ...
           'repair_dro_elfo_cells:engineOpts');
assert(isequaln(load(s5), Qb), 'refused call leaves the sheet untouched');
Qk = rmfield(Qb, 'sheet_key');  save(s5, '-struct', 'Qk');
expect_err(@() repair_dro_elfo_cells(s5, 'basin', o1), 'repair_dro_elfo_cells:sheetKey');
Qk = Qb;  Qk.sheet_key.raanDeg = 90;  save(s5, '-struct', 'Qk');
expect_err(@() repair_dro_elfo_cells(s5, 'basin', o1), 'repair_dro_elfo_cells:sheetKey');

%% 6. interrupted repair: keyless cleared sheet + pending backup with the key
s6 = newSheet('interrupt');  Qr = load(s6);
bk = struct('TF', squeeze(Qr.TF(1, 1, :)), 'FLYKM', squeeze(Qr.FLYKM(1, 1, :)), ...
            'ACCDZ', squeeze(Qr.ACCDZ(1, 1, :)), 'RES', squeeze(Qr.RES(1, 1, :)), ...
            'WALL', squeeze(Qr.WALL(1, 1, :)), 'OK', squeeze(Qr.OK(1, 1, :)), ...
            'Z8', reshape(Qr.Z8(:, 1, 1, :), 8, []), 'ATT', Qr.ATT(1, 1));
row = struct('mode', 'basin', 'iD', 1, 'iA', 1, 'seedFrom', struct('iD', 1, 'iA', 1, 'kr', 2), ...
             'tfBefore', bk.TF(1), 'tfAfter', NaN, 'okAfter', false, 'outcome', '', 'conjBefore', 0, ...
             'conjAfter', NaN, 'date', 'x', 'wallS', NaN, 'directGateKm', inf, ...
             'keptOld', zeros(1, 0), 'lostNew', zeros(1, 0));
Bs = struct('bk', bk, 'mode', 'basin', 'iD', 1, 'iA', 1, 'seedFrom', row.seedFrom, 'date', 'x', ...
            'row', row, 'sheet_key', Qr.sheet_key);
save([s6(1:end-4) '_repair_backup_i01_a01.mat'], '-struct', 'Bs');
Qc = rmfield(Qr, 'sheet_key');
Qc.TF(1, 1, :) = NaN;  Qc.OK(1, 1, :) = false;  Qc.Z8(:, 1, 1, :) = NaN;  Qc.ATT(1, 1) = 1;
save(s6, '-struct', 'Qc');
R6 = repair_dro_elfo_cells(s6, 'basin', with(o1, 'conjOnly', true));
Qd = load(s6);
assert(R6.nRecovered == 1 && isequal(Qd.OK, Qr.OK) && isequaln(Qd.Z8, Qr.Z8) && ...
       isequaln(Qd.sheet_key, Qr.sheet_key), 'interrupted repair: cell and key restored');
fprintf('  seams: C5, fill/gate, restore x2, throw, refusals x4, recovery -- all as expected\n');

%% 7. real engine: unsolved cell recovered from its neighbour
if doSolve
    s7 = fixture_dro_elfo_sheet(fullfile(wd, 'real', 'dro_elfo_tau1_raan000.mat'), ...
                                struct('iD', 1:6, 'iA', 2, 'nR', 2));
    Q7 = load(s7);
    tfOrig = Q7.TF(3, 1, 1);
    Q7.TF(3, 1, :) = NaN;  Q7.FLYKM(3, 1, :) = NaN;  Q7.ACCDZ(3, 1, :) = NaN;  Q7.RES(3, 1, :) = NaN;
    Q7.WALL(3, 1, :) = NaN;  Q7.OK(3, 1, :) = false;  Q7.Z8(:, 3, 1, :) = NaN;  Q7.ATT(3, 1) = 2;
    save(s7, '-struct', 'Q7');
    t7 = tic;
    R7 = repair_dro_elfo_cells(s7, 'unsolved', struct('engineOpts', eo, 'logFile', fullfile(wd, 'real', 'repair.log')));
    Q7a = load(s7);  S7 = load(R7.sideMat);
    fprintf('  real unsolved repair: cell 3 seeded from cell %d, outcome "%s", 15 N tf %.6f (fleet %.6f), 12 N OK %d, %.0f s\n', ...
            S7.rows(end).seedFrom.iD, S7.rows(end).outcome, Q7a.TF(3, 1, 1), tfOrig, Q7a.OK(3, 1, 2), toc(t7));
    assert(R7.nCandidates == 1 && S7.rows(end).seedFrom.iD == 2, 'seeded from cell 2');
    assert(R7.nRepaired == 1 && Q7a.OK(3, 1, 1), 'the emptied cell is recovered at 15 N');
    assert(isequaln(Q7a.sheet_key, Q7.sheet_key), 'sheet_key on the sheet after a real engine call');
    assert(isequal(Q7a.OK([1 2 4 5 6], 1, :), Q7.OK([1 2 4 5 6], 1, :)) && ...
           isequaln(Q7a.Z8(:, [1 2 4 5 6], 1, :), Q7.Z8(:, [1 2 4 5 6], 1, :)), 'other cells untouched');
end

ok = true;
fprintf('test_repair_dro_elfo_cells: PASS\n');
end

% ------------------------------------------------------------------------
function fake_engine(sheetMat, e, spec)
% FAKE_ENGINE  Stand-in for thrust_ladder_library(sheetMat, e) on a redo of
% one cell: writes spec's rungs (OK, TF, flown miss; z8 distinct from the
% stored ones) and saves the engine's variables only (sheet_key dropped,
% as the real engine's full-file save drops it).
% INPUTS: sheetMat [1 x n] char; e [1 x 1] engine options (.cells [1 x 2]
% .redo); spec [1 x 1] struct .tf .ok .fly [1 x nR].  OUTPUTS: none (file
% written).
assert(e.redo && size(e.cells, 1) == 1, 'fake_engine: one redo cell');
Q = load(sheetMat);
iD = e.cells(1);  iA = e.cells(2);
for kr = 1:numel(spec.ok)
    if spec.ok(kr)
        Q.OK(iD, iA, kr) = true;  Q.TF(iD, iA, kr) = spec.tf(kr);  Q.FLYKM(iD, iA, kr) = spec.fly(kr);
        Q.ACCDZ(iD, iA, kr) = 1e-9;  Q.Z8(:, iD, iA, kr) = [0.1*kr*ones(7, 1); spec.tf(kr)];
    else
        Q.OK(iD, iA, kr) = false;  Q.TF(iD, iA, kr) = NaN;  Q.FLYKM(iD, iA, kr) = NaN;
        Q.ACCDZ(iD, iA, kr) = NaN;  Q.Z8(:, iD, iA, kr) = NaN;
    end
end
Q.ATT(iD, iA) = 1;
if isfield(Q, 'sheet_key'), Q = rmfield(Q, 'sheet_key'); end
save(sheetMat, '-struct', 'Q');
end

function o = with(o, f, v)
% WITH  Copy of struct o with field f set to v.  INPUTS: o [1 x 1] struct;
% f [1 x n] char; v [any].  OUTPUTS: o [1 x 1] struct.
o.(f) = v;
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
error('test_repair_dro_elfo_cells:noThrow', 'expected %s, nothing thrown', id);
end
