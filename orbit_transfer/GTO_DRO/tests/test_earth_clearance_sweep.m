function ok = test_earth_clearance_sweep()
% TEST_EARTH_CLEARANCE_SWEEP  Sidecar census, writeback and the demotion path.
% INPUTS: none.  OUTPUTS: ok - logical.
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
addpath(here);
src = fixture_demo_sheet();
wd  = fullfile(tempdir, 'test_earth_sweep');
if isfolder(wd), rmdir(wd, 's'); end
mkdir(wd);
sheetMat = fullfile(wd, 'demo_sheet.mat');
copyfile(src, sheetMat);

% census at the campaign floor: both entries clear it (the engine solved
% with the floor and the flown arc is within 1e-4 km of the direct one)
S = earth_clearance_sweep(sheetMat);
assert(S.done && S.nDone == 2 && S.nTodo == 0, 'census must cover both OK entries');
assert(isfile(S.sideMat), 'sidecar written');
E = load(S.sideMat, 'EARTHKM', 'floorKm');
assert(E.floorKm == 200, 'floor defaults to meta.floorEarthKm');
assert(all(E.EARTHKM(1, 1, :) > 200), 'both flown arcs clear 200 km');
assert(S.nBelow == 0, 'nothing below the floor');
fprintf('  min Earth altitude over the two entries: %.0f km\n', S.minAltKm);

% resume: a second call does no flights (done mask) and reports the same
S2 = earth_clearance_sweep(sheetMat);
assert(S2.nDone == 2 && abs(S2.minAltKm - S.minAltKm) < 1e-9, 'resume reproduces the census');

% writeback at the campaign floor: nothing demoted, the one-time backup is created, meta stamped
S3 = earth_clearance_sweep(sheetMat, struct('writeback', true));
Q = load(sheetMat, 'OK', 'meta');
assert(all(Q.OK(:)) && isfield(Q.meta, 'earthSweep') && Q.meta.earthSweep.nDemoted == 0, 'clean writeback');

% Z8 identity: perturb one costate of entry 2, only that entry is re-flown
Qz = load(sheetMat);  Z8 = Qz.Z8;  Z8(3, 1, 1, 2) = Z8(3, 1, 1, 2)*(1 + 1e-9);
save(sheetMat, '-append', 'Z8');
S5 = earth_clearance_sweep(sheetMat);
assert(S5.nFlown == 1 && S5.nDone == 2 && S5.done, 'exactly one changed entry re-flown');
S5b = earth_clearance_sweep(sheetMat);
assert(S5b.nFlown == 0, 'nothing re-flown when Z8 is unchanged');

% budget stop: maxEntries = 1 on a fresh sidecar, then a second call finishes
delete(S.sideMat);
S6 = earth_clearance_sweep(sheetMat, struct('maxEntries', 1));
assert(~S6.done && S6.nTodo == 1 && S6.nFlown == 1 && isfile(S6.sideMat), 'budget stop is partial and leaves a sidecar');
try
    earth_clearance_sweep(sheetMat, struct('maxEntries', 0, 'writeback', true));
    error('test:noThrow', 'writeback on an incomplete census must refuse');
catch err
    assert(strcmp(err.identifier, 'earth_clearance_sweep:incomplete'), 'wrong error: %s', err.identifier);
end
S7 = earth_clearance_sweep(sheetMat);
assert(S7.done && S7.nFlown == 1 && S7.nDone == 2, 'second call finishes the census');

% non-finite flight on an OK entry counts as below and is demoted
nanMat = fullfile(wd, 'nan_sheet.mat');
copyfile(sheetMat, nanMat);
Qn = load(nanMat);  Z8 = Qn.Z8;  Z8(1:7, 1, 1, 2) = NaN;  save(nanMat, '-append', 'Z8');
Sn = earth_clearance_sweep(nanMat);
assert(Sn.nNonFinite == 1 && Sn.nBelow == 1 && Sn.nDone == 2, 'NaN entry counted below');
Sn = earth_clearance_sweep(nanMat, struct('writeback', true));
Qn = load(nanMat, 'OK');
assert(Qn.OK(1, 1, 1) && ~Qn.OK(1, 1, 2) && isequal(Sn.demoted, [1 1 2]), 'NaN entry demoted, other kept');

% an existing backup is never overwritten
bakFile = [sheetMat(1:end-4) '.bak_earth'];
marker = 'kept';  save(bakFile, 'marker', '-mat');
earth_clearance_sweep(sheetMat, struct('floorKm', 1e6, 'writeback', true));
Bm = load(bakFile, '-mat');
assert(isfield(Bm, 'marker') && strcmp(Bm.marker, 'kept'), 'existing backup must not be overwritten');
delete(bakFile);
% restore the original sheet for the impossible-floor case below
copyfile(src, sheetMat);

% impossible floor: both demoted, backup kept, demoted list names them
delete(S.sideMat);
S4 = earth_clearance_sweep(sheetMat, struct('floorKm', 1e6, 'writeback', true));
assert(S4.nBelow == 2, 'both entries below an impossible floor');
Q = load(sheetMat, 'OK', 'meta');
assert(~any(Q.OK(:)), 'demoted entries are OK = false');
assert(isfile([sheetMat(1:end-4) '.bak_earth']), 'pre-sweep backup kept');
assert(isequal(sortrows(Q.meta.earthSweep.demoted), [1 1 1; 1 1 2]), 'demoted list');
bakInfo = dir([sheetMat(1:end-4) '.bak_earth']);
% a second writeback keeps the demotion record (the chain calls it every run)
S4b = earth_clearance_sweep(sheetMat, struct('floorKm', 1e6, 'writeback', true));
Q2  = load(sheetMat, 'OK', 'meta');
assert(S4b.nBelow == 2 && Q2.meta.earthSweep.nDemoted == 2, 'second writeback keeps nDemoted = 2');
assert(isequal(sortrows(Q2.meta.earthSweep.demoted), [1 1 1; 1 1 2]), 'demoted list unchanged');
bakInfo2 = dir([sheetMat(1:end-4) '.bak_earth']);
assert(bakInfo2.bytes == bakInfo.bytes && bakInfo2.datenum == bakInfo.datenum, 'backup untouched');
B = load([sheetMat(1:end-4) '.bak_earth'], '-mat', 'OK');
assert(all(B.OK(:)), 'the backup holds the original OK');

rmdir(wd, 's');
ok = true;
fprintf('test_earth_clearance_sweep: PASS\n');
end
