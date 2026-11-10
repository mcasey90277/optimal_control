function ok = test_run_gto_dro_catalog()
% TEST_RUN_GTO_DRO_CATALOG  planOnly contract and the on-disk rung guard.
% INPUTS: none.  OUTPUTS: ok - logical.
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
addpath(fullfile(fileparts(here), 'catalog'));

plan = run_gto_dro_catalog([], [], [], [], struct('planOnly', true));
assert(numel(plan) == 16, '16 sheets');
tags = {plan.tag};
assert(numel(unique(tags)) == 16, 'distinct tags');
assert(strcmp(plan(2).tag, 'gto_dro_o000_tau1') && plan(2).orientDeg == 0 && plan(2).tau == 1, 'sheet 2 is the pilot');
assert(strcmp(plan(5).tag, 'gto_dro_o090_tau0.5'), 'tau 0.5 tag spelling');
for k = 1:16
    oA = plan(k).optsA;
    assert(strcmp(oA.depFamily, 'gto') && strcmp(oA.arrFamily, 'dro'), 'families');
    assert(oA.depParams.orientDeg == plan(k).orientDeg, 'orientation');
    assert(oA.arrParams.tau == plan(k).tau, 'DRO period in arrParams');
    assert(oA.tauDRO == plan(k).orientDeg, 'sheet key rule: tauDRO = orientDeg');
    assert(oA.floorEarthKm == 200 && oA.floorKm == 500, 'floors');
    assert(oA.thrLock && abs(oA.tf0 - 0.30) < 1e-12, 'warm recipe');
    assert(oA.ispS == 1710 && oA.m0kg == 150 && oA.N == 400 && oA.maxIter == 6000 && oA.maxCpuSec == 600, 'regime + solver');
    assert(oA.nD == 12 && oA.nA == 6 && oA.sD0 == 0 && oA.sA0 == 0, 'grid');
    assert(isequal(oA.rungs, [15 12 10 7 5]), 'default rungs');
    oB = plan(k).optsB;
    assert(strcmp(oB.depFamily, 'gto') && strcmp(oB.arrFamily, 'dro'), 'B families');
    assert(oB.arrParams.tau == plan(k).tau && oB.depParams.orientDeg == plan(k).orientDeg, 'B params');
    assert(oB.floorEarthKm == 200 && oB.floorKm == 500, 'B floors');
    assert(oB.tauDRO == plan(k).orientDeg, 'B sheet key rule');
    assert(~isfield(plan(k).optsB, 'tf0') || isempty(plan(k).optsB.tf0), 'pass B is the cold tf0');
    assert(isempty(plan(k).optsB.cells), 'pass B cells are scoped at run time, empty in the plan');
end
plan9 = run_gto_dro_catalog(2, [], [], [15 12 10 7 5 3 2 1.5 1], struct('planOnly', true, 'tf0', 0.27));
assert(numel(plan9) == 1 && isequal(plan9.optsA.rungs, [15 12 10 7 5 3 2 1.5 1]) && abs(plan9.optsA.tf0 - 0.27) < 1e-12, 'overrides');

% planOnly touches no disk: an absent catDir stays absent
wd = fullfile(tempdir, 'test_run_gto_dro_catalog');
if isfolder(wd), rmdir(wd, 's'); end
run_gto_dro_catalog(2, [], [], [], struct('planOnly', true, 'catDir', wd));
assert(~isfolder(wd), 'planOnly created the catalog folder');

% on-disk rung guard: a fake sheet with rungs [15 12] refuses a 5-rung call
% (the refusal fires before the engine is reached, so a rungs-only file is enough)
mkdir(wd);
rungs = [15 12];
save(fullfile(wd, 'gto_dro_o000_tau1.mat'), 'rungs');
try
    run_gto_dro_catalog(2, 0, 1, [], struct('catDir', wd));
    error('test:noThrow', 'rung mismatch did not refuse');
catch err
    assert(strcmp(err.identifier, 'run_gto_dro_catalog:rungMismatch'), 'wrong error: %s', err.identifier);
end
delete(fullfile(wd, 'gto_dro_o000_tau1.mat'));
% a matching call with maxCells = 0 runs the engine on zero cells and returns
% (no sheet on disk: the engine's resume path needs an engine-valid sheet,
% which a rungs-only file is not)
run_gto_dro_catalog(2, 0, 1, [15 12], struct('catDir', wd));
assert(isfile(fullfile(wd, 'gto_dro_o000_tau1_progress.txt')), 'engine log written');
rmdir(wd, 's');
ok = true;
fprintf('test_run_gto_dro_catalog: PASS\n');
end
