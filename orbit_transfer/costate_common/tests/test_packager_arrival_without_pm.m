function ok = test_packager_arrival_without_pm()
%% Purpose:
%
%   Pins build_costate_catalog_family's handling of the arrival branch sign
%   pm (2026-11-10): a non-tulip arrival whose arrParams carries NO pm (a DRO
%   arrival: the GTO -> DRO catalog) packages with sheets.pm = NaN instead of
%   throwing, a non-tulip arrival whose arrParams HAS pm (a halo arrival)
%   keeps its value, and the catalog passes catalog_schema validation. The
%   sheets are synthetic (hand-built minimal engine-schema sheets in a
%   tempdir), so the test needs no solver, no pumpkyn and no campaign.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every check passed
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
pathBefore = path;
addpath(fileparts(here));
wd = tempname;
mkdir(wd);
cleaner = onCleanup(@() cleanup(wd, pathBefore));
ok = true;

% two one-cell sheets, one arrival family ('halo': the packager only reads
% the family string), arrival 1 without pm, arrival 2 with pm = -1
write_sheet(fullfile(wd, 'sheet_nopm.mat'), 1, struct('tau', 1));
write_sheet(fullfile(wd, 'sheet_pm.mat'),   2, struct('tau', 2, 'pm', -1));

spec = struct('glob', 'sheet_*.mat', 'name', 'cat_pm_test', 'description', 'synthetic', ...
              'provenance', 'synthetic', 'depReconstruction', 'synthetic');
outMat = fullfile(wd, 'cat_pm_test.mat');
cat_ = build_costate_catalog_family(wd, outMat, spec);

tauArr = [cat_.sheets.tau_arr];
sNo  = cat_.sheets(tauArr == 1);
sYes = cat_.sheets(tauArr == 2);
ok = chk(ok, numel(cat_.sheets) == 2, 'both sheets packaged');
ok = chk(ok, isnan(sNo.pm), 'arrival without pm -> pm = NaN');
ok = chk(ok, sYes.pm == -1, 'arrival with pm keeps its value (-1)');
ok = chk(ok, isnan(sNo.Np) && isnan(sYes.Np), 'non-tulip arrivals carry Np = NaN');
ok = chk(ok, isempty(catalog_schema('validate', cat_)), 'catalog validates clean');
L = load(outMat);
ok = chk(ok, isfield(L, 'cat_pm_test'), 'catalog saved under its spec name');

if ok, fprintf('TEST_PACKAGER_ARRIVAL_WITHOUT_PM: ALL PASS\n');
else,  fprintf('TEST_PACKAGER_ARRIVAL_WITHOUT_PM: FAILURE\n');
end
end

function cleanup(wd, pathBefore)
% CLEANUP  Remove the temp folder and restore the path.  INPUTS: wd char;
% pathBefore char.  OUTPUTS: none.
if isfolder(wd), rmdir(wd, 's'); end
path(pathBefore);
end

function write_sheet(path, tauArr, arrParams)
% WRITE_SHEET  A minimal 1 x 1 x 2 engine-schema sheet (rungs [15 12], both
% rungs solved) with a gto-like departure and a non-tulip arrival.
% INPUTS: path char; tauArr scalar (arrival key); arrParams struct.
% OUTPUTS: none (writes path).
meta = struct('muStar', 0.012150585609624, 'lStar', 389703.264829278, ...
              'tStar', 382981.289129055, 'tauDRO', 0, 'NpTulip', 7, 'tauTulip', 5.2, ...
              'pmTulip', -1, 'depFamily', 'gto', ...
              'depParams', struct('sma_km', 24446, 'ecc', 0.72, 'orientDeg', 0), ...
              'arrFamily', 'halo', 'arrParams', arrParams, 'periodDRO', 0.0993, ...
              'periodTulip', tauArr, 'tauDep', 0.0993, 'tauArr', tauArr, 'sheetSchema', 2, ...
              'ispS', 1710, 'm0kg', 150, 'N', 400);
S = struct();
S.rungs = [15 12];
S.sD = 0.10;  S.sA = 0;
S.TF = reshape([0.30 0.31], 1, 1, 2);
S.OK = true(1, 1, 2);
S.Z8 = zeros(8, 1, 1, 2);
S.Z8(:, 1, 1, 1) = (1:8).' / 10;
S.Z8(:, 1, 1, 2) = (2:9).' / 10;
S.meta = meta;
save(path, '-struct', 'S');
end

function ok = chk(ok, cond, msg)
% CHK  PASS/FAIL line.  INPUTS: ok; cond logical; msg char.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
