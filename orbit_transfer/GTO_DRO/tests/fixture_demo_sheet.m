function sheetMat = fixture_demo_sheet(rebuild)
% FIXTURE_DEMO_SHEET  A cached 1 x 1 GTO -> DRO sheet for the stage 2 tests.
%
%   The stage 1 demo cell (orient 0, DRO tau 1, departure phase 0.10,
%   arrival phase 0, Earth floor 200 km) at rungs [15 12], solved through
%   thrust_ladder_library exactly as the fleet would solve it, saved under
%   results/test_fixture/ (git-ignored) and reused by every test that
%   needs a real entry. The 15 N cold start lands on the 0.3040 ND basin
%   (refuted by the conjugate test, stage 1) and the 12 N rung on 0.3070 ND.
%
% INPUTS:
%   rebuild  - true forces a fresh solve (default false)          [scalar]
%
% OUTPUTS:
%   sheetMat - path of the fixture sheet .mat                     [char]
%
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.

if nargin < 1, rebuild = false; end
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
fixDir   = fullfile(fileparts(here), 'results', 'test_fixture');
if ~isfolder(fixDir), mkdir(fixDir); end
sheetMat = fullfile(fixDir, 'demo_sheet.mat');
if isfile(sheetMat) && ~rebuild
    Q = load(sheetMat, 'OK', 'rungs');
    if isequal(Q.rungs, [15 12]) && all(Q.OK(:)), return, end
end
if isfile(sheetMat), delete(sheetMat); end
opts = struct('rungs', [15 12], 'ispS', 1710, 'm0kg', 150, 'N', 400, ...
    'floorKm', 500, 'floorEarthKm', 200, 'maxIter', 6000, 'maxCpuSec', 600, ...
    'gateKm', 100, 'accTol', 1e-6, 'nD', 1, 'nA', 1, 'sD0', 0.10, 'sA0', 0, ...
    'depFamily', 'gto', 'depParams', struct('orientDeg', 0), ...
    'arrFamily', 'dro', 'arrParams', struct('tau', 1), 'tauDRO', 0, ...
    'thrLock', true, 'tf0', 0.30, 'resume', false, ...
    'logFile', fullfile(fixDir, 'demo_sheet.log'));
P = thrust_ladder_library(sheetMat, opts);
assert(all(P.OK(:)), 'fixture_demo_sheet: the demo cell did not solve at both rungs');
end
