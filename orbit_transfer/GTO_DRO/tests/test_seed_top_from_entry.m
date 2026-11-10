function ok = test_seed_top_from_entry()
% TEST_SEED_TOP_FROM_ENTRY  Flown 12 N entry -> engine seedTop struct.
% INPUTS: none.  OUTPUTS: ok - logical.
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
addpath(here);
sheetMat = fixture_demo_sheet();
Q = load(sheetMat, 'Z8', 'rungs', 'meta', 'OK');
N = Q.meta.N;

s = seed_top_from_entry(sheetMat, 1, 1, 2, 15);
assert(isequal(size(s.X), [7 N+1]) && isequal(size(s.U), [4 N+1]), 'seedTop shape');
assert(all(s.U(4, :) == 1), 'thrust-locked seed must carry thr = 1');
assert(all(abs(vecnorm(s.U(1:3, :), 2, 1) - 1) < 1e-12), 'direction columns must be unit');
assert(abs(s.X(7, 1) - 1) < 1e-12, 'mass fraction starts at 1');
assert(s.X(7, end) < 1 && s.X(7, end) > 0.3, 'mass fraction ends in (0.3, 1)');
tf12 = Q.Z8(8, 1, 1, 2);
assert(abs(s.tf - tf12*(12/15)^0.6) < 1e-12, 'tf scaled by the engine rule (Tprev/TN)^0.6');
assert(s.from.iD == 1 && s.from.iA == 1 && s.from.kr == 2 && s.from.thrustN == 12, 'provenance');
% the resampled departure state is the sheet's own departure state
[tD, rvD] = get_family_orbit(Q.meta.depFamily, Q.meta.depParams);
dep = phase_state(tD, rvD);  rv0 = dep(0.10);
assert(max(abs(s.X(1:6, 1) - rv0(:))) < 1e-9, 'seed starts at the cell departure state');
% a non-OK cell refuses
badMat = fullfile(tempdir, 'seed_top_bad.mat');
S = load(sheetMat);  S.OK(1, 1, 2) = false;  save(badMat, '-struct', 'S');
try
    seed_top_from_entry(badMat, 1, 1, 2, 15);
    error('test:noThrow', 'a non-OK entry must refuse');
catch err
    assert(strcmp(err.identifier, 'seed_top_from_entry:notOK'), 'wrong error: %s', err.identifier);
end
delete(badMat);
% a NaN costate with OK forced true: refused, never a non-finite seed
nanMat = fullfile(tempdir, 'seed_top_nan.mat');
S = load(sheetMat);  S.Z8(1:7, 1, 1, 2) = NaN;  save(nanMat, '-struct', 'S');
try
    seed_top_from_entry(nanMat, 1, 1, 2, 15);
    error('test:noThrow', 'a NaN-costate entry must refuse');
catch err
    fprintf('  NaN-costate entry refused with %s\n', err.identifier);
    assert(~strcmp(err.identifier, 'test:noThrow'), 'a NaN-costate entry must refuse');
end
delete(nanMat);
ok = true;
fprintf('test_seed_top_from_entry: PASS\n');
end
