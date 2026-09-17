function ok = test_thrust_ladder()
%% Purpose:
%
%   Tests thrust_ladder_library, the engine five campaigns solve through,
%   against a SHIPPED sheet: one stored cell's top rungs are re-solved from
%   that sheet's own meta and compared with what the campaign recorded.
%
%   Why this is a fair test rather than a re-run: each rung is warm-started
%   from the rung above, so a ladder truncated to the first nRung rungs
%   reproduces the stored run's first nRung rungs exactly. The comparison is
%   therefore against numbers produced months earlier, by the same engine in
%   its previous home (`DRO_tulip/indirect`), on a campaign this folder
%   knows nothing about.
%
%   Checks: every requested rung verified (OK matches the sheet), t_f within
%   1e-10 of the stored value, costates within 1e-8, and the returned
%   arrays/meta shaped as the packagers expect.
%
%   The fixture is a campaign results file, which is gitignored: with no
%   sheet on disk the test SKIPS (and says so) rather than failing.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Checks passed (or
%                                                   skipped: no fixture)
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
cc = fileparts(fileparts(mfilename('fullpath')));
ot = fileparts(cc);
addpath(cc, fullfile(fileparts(ot), 'oclib'));
if isempty(which('casadi.SX')), addpath(fullfile(getenv('HOME'), 'casadi-3.7.0')); end

sheetMat = fullfile(ot, 'HALO_tulip', 'direct', 'results', 'catalog', 'halo_tau1p75_Np5.mat');
if ~isfile(sheetMat)
    fprintf('  SKIP  no shipped sheet on disk (%s)\n', sheetMat);
    return
end
S = load(sheetMat);
m = S.meta;
nRung = 2;
[iD, iA] = find(S.OK(:,:,1), 1);

opts = struct('rungs', S.rungs(1:nRung), 'ispS', m.ispS, 'm0kg', m.m0kg, ...
    'N', m.N, 'floorKm', m.floorKm, 'gateKm', m.gateKm, ...
    'depFamily', m.depFamily, 'depParams', m.depParams, ...
    'arrFamily', m.arrFamily, 'arrParams', m.arrParams, ...
    'NpTulip', m.NpTulip, 'tauTulip', m.tauTulip, 'pmTulip', m.pmTulip, ...
    'tauDRO', m.tauDRO, 'nD', numel(S.sD), 'nA', numel(S.sA), ...
    'sD0', m.sD0, 'sA0', m.sA0, 'cells', [iD iA], 'maxCells', 1, ...
    'resume', false, 'maxIter', 3000, 'accTol', 1e-6, ...
    'logFile', fullfile(tempdir, 'test_thrust_ladder.log'));
outMat = fullfile(tempdir, sprintf('test_thrust_ladder_%d.mat', feature('getpid')));
cleaner = onCleanup(@() deleteIfThere(outMat));

P = thrust_ladder_library(outMat, opts);

okNew = squeeze(P.OK(iD,iA,1:nRung));   okRef = squeeze(S.OK(iD,iA,1:nRung));
tfNew = squeeze(P.TF(iD,iA,1:nRung));   tfRef = squeeze(S.TF(iD,iA,1:nRung));
z8New = squeeze(P.Z8(:,iD,iA,1:nRung)); z8Ref = squeeze(S.Z8(:,iD,iA,1:nRung));

ok = chk(ok, all(okNew) && isequal(okNew, okRef), ...
         sprintf('cell (%d,%d): both rungs verified, as the sheet records', iD, iA));
ok = chk(ok, max(abs(tfNew - tfRef)) < 1e-10, ...
         sprintf('t_f reproduces the shipped sheet: max |dt_f| = %.2e (%s)', ...
                 max(abs(tfNew - tfRef)), mat2str(tfNew.', 10)));
ok = chk(ok, max(abs(z8New(:) - z8Ref(:))) < 1e-8, ...
         sprintf('costates reproduce it: max |dz8| = %.2e', max(abs(z8New(:) - z8Ref(:)))));
ok = chk(ok, isequal(size(P.TF), [numel(S.sD) numel(S.sA) nRung]) ...
             && isequal(size(P.Z8), [8 numel(S.sD) numel(S.sA) nRung]), ...
         'arrays are [nD x nA x nRung] and [8 x nD x nA x nRung]');
ok = chk(ok, isfile(outMat), 'the sheet is written to disk (it saves after every rung)');
need = {'muStar','lStar','tStar','depFamily','depParams','arrFamily','arrParams','ispS','m0kg','N'};
ok = chk(ok, all(isfield(P.meta, need)), 'meta carries what the packager and the endpoint rebuild need');

if ok, fprintf('TEST_THRUST_LADDER: ALL PASS\n');
else,  fprintf('TEST_THRUST_LADDER: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end

function deleteIfThere(f)
%% Purpose:
%
%   Remove a temporary file if it exists.
%
if isfile(f), delete(f); end
end
