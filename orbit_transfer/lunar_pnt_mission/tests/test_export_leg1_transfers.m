function [ok, skipped] = test_export_leg1_transfers()
%% Purpose:
%
%   Verifies export_leg1_transfers on the SAVED leg-1 study results
%   (nothing re-polished): (1) ROUND TRIP: the sidecar (written to a
%   temporary name INSIDE results/, deleted afterwards) holds ONE variable,
%   leg1_transfers, with one entry per accepted transfer, in satellite
%   order, whose z8, t_f, endpoint states, epochs, orientation, phases and
%   seed cell equal the study's P1; every entry RE-FLIES (fly_transfer with
%   the stored regime) from its stored rv0 to its stored rvf with the
%   stored miss (1e-9 km, 1e-9 m/s) and that miss is below 1 km; the
%   verdict block carries the located clearances and the conjugate verdict
%   of the study's gates, and the not-checked list names the catalog
%   acceptance, the Earth floor in the optimisation and library admission;
%   the package identity equals the study's; (2) a transfer marked not
%   accepted (tampered copy of the results in a temporary folder) is not
%   exported; results with none accepted are refused (noTransfers); an
%   ACCEPTED transfer whose stored G.pass is false is refused (notPassed);
%   the note is built from the stored verdicts (its gate counts and
%   minimum clearances equal the entries');
%   (3) REFUSALS to write outside results/: an absolute path in a temporary
%   folder, a path that leaves results/ through '..', and a sibling folder
%   whose name starts with "results" -- each throws
%   export_leg1_transfers:outFile and creates nothing. Returns skipped =
%   true (and ok = true) when the results file does not exist.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 true if every assertion
%                                                   held (else it throws)
%  skipped                  logical                 true if the results file
%                                                   is missing (nothing run)
%
%% Revision History:
%  M. Casey                                                   (c) 11/11/2026
%  M. Casey  fix round 1: notPassed refusal, note from verdicts  (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
run(fullfile(root, 'setup_lunar_pnt_paths.m'));
resDir = fullfile(root, 'results');
resultsMat = fullfile(resDir, 'lunar_pnt_leg1_study.mat');
skipped = false;
if ~isfile(resultsMat)
    fprintf('  SKIPPED: no results file (%s); run lunar_pnt_leg1_study first\n', resultsMat);
    ok = true;  skipped = true;
    return
end
scratch = '/private/tmp/claude-501/-Users-msc-Desktop-optimal-control/70650963-d331-4c2e-9b00-6156c2bd5033/scratchpad/leg1';
if ~isfolder(scratch), mkdir(scratch); end
tmp = tempname(scratch);
mkdir(tmp);
cleanTmp = onCleanup(@() rmdir(tmp, 's'));
outTest = fullfile(resDir, sprintf('.leg1_export_test_%s.mat', char(java.util.UUID.randomUUID())));
cleanOut = onCleanup(@() delete_if(outTest));
R = load(resultsMat, 'P1', 'run_meta');
acc = find([R.P1.accepted]);

% ---- 1. round trip ----
E = export_leg1_transfers(resultsMat, struct('outFile', outTest, 'verbose', true));
w = whos('-file', outTest);
assert(isscalar(w) && strcmp(w.name, 'leg1_transfers'), 'sidecar must hold ONE variable leg1_transfers');
X = load(outTest);
X = X.leg1_transfers;
assert(isequaln(X.entries, E.entries), 'returned and written entries differ');
assert(numel(X.entries) == numel(acc), '%d entries for %d accepted transfers', numel(X.entries), numel(acc));
rg = X.regime;
pr = nd_propulsion(rg.thrustN, rg.ispS, rg.m0kg, rg.lStar, rg.tStar);
phys = struct('Tnd', pr.Tnd, 'cnd', pr.cnd, 'muStar', rg.muStar, 'lStar', rg.lStar, ...
              'tStar', rg.tStar, 'm0kg', rg.m0kg);
assert(strcmp(X.package.sha256, R.run_meta.identity.sha256), 'package identity differs from the study''s');
for ke = 1:numel(X.entries)
    e = X.entries(ke);
    q = R.P1(acc(ke));
    assert(e.sat == acc(ke), 'entry %d is satellite %d, expected %d', ke, e.sat, acc(ke));
    same = isequal(e.z8, q.z8(:)) && e.tf_nd == q.tf && isequal(e.rv0, q.rv0(:)) && ...
           isequal(e.rvf, q.rvf(:)) && e.tDep == q.tDep && e.tArr == q.tArr && ...
           e.orientDeg == q.orientDeg && e.sD == q.sD && e.arrPhase == q.arrPhase && ...
           isequal(e.seededFrom, [q.iD q.iA q.iS]);
    assert(same, 'entry %d does not carry the study''s transfer', ke);
    F = fly_transfer(e.z8, e.rv0, e.rvf, phys);
    dKm = abs(F.flyKm - e.flyKm);  dV = abs(F.flyVms - e.flyVms);
    assert(dKm < 1e-9 && dV < 1e-9, 'entry %d re-flies to %.6g km (stored %.6g)', ke, F.flyKm, e.flyKm);
    assert(e.flyKm < 1, 'entry %d stored miss %.3g km', ke, e.flyKm);
    v = e.verdicts;
    assert(v.earthMinKm == q.G.earth.hMinKm && v.moonMinKm == q.G.moon.hMinKm, 'entry %d clearances', ke);
    assert(strcmp(v.secondStatus, q.G.second.status), 'entry %d second-order status', ke);
    nc = strjoin(e.notChecked, ' ');
    assert(contains(nc, 'tfMin') && contains(nc, 'Earth floor') && contains(nc, 'library'), ...
           'entry %d: not-checked list incomplete', ke);
    fprintf('  entry %d (sat %d): re-flown %.4g km %.3g m/s (stored %.4g km), status %s, Earth %.1f km, Moon %.1f km\n', ...
            ke, e.sat, F.flyKm, F.flyVms, e.flyKm, v.secondStatus, v.earthMinKm, v.moonMinKm);
end

nE = numel(X.entries);
V = [X.entries.verdicts];
want = {sprintf('all gates %d of %d', nnz([V.gatesPass]), nE), ...
        sprintf('hand-off at T0 %d of %d', nnz([V.handoffPass]), nE), ...
        sprintf('min located %.1f km), Moon', min([V.earthMinKm])), ...
        sprintf('study verdict "%s"', R.run_meta.verdict)};
for kw = 1:numel(want)
    assert(contains(X.note, want{kw}), 'note does not carry "%s"', want{kw});
end
fprintf('  note built from the verdicts: %s\n', strjoin(want, '; '));

% ---- 2. not-accepted transfers are not exported ----
T = load(resultsMat);
T.P1(acc(1)).accepted = false;
tamp = fullfile(tmp, 'tampered.mat');
save(tamp, '-struct', 'T');
E2 = export_leg1_transfers(tamp, struct('outFile', outTest, 'verbose', false));
assert(numel(E2.entries) == numel(acc) - 1 && ~any([E2.entries.sat] == acc(1)), ...
       'a not-accepted transfer was exported');
fprintf('  satellite %d marked not accepted: %d entries, it is absent\n', acc(1), numel(E2.entries));
T2 = load(resultsMat);
T2.P1(acc(2)).G.pass = false;
save(tamp, '-struct', 'T2');
id = '';
try
    export_leg1_transfers(tamp, struct('outFile', outTest, 'verbose', false));
catch err
    id = err.identifier;
end
assert(strcmp(id, 'export_leg1_transfers:notPassed'), 'accepted with G.pass false: expected notPassed, got "%s"', id);
fprintf('  satellite %d accepted with G.pass false: refused (%s)\n', acc(2), id);
[T.P1.accepted] = deal(false);
save(tamp, '-struct', 'T');
id = '';
try
    export_leg1_transfers(tamp, struct('outFile', outTest, 'verbose', false));
catch err
    id = err.identifier;
end
assert(strcmp(id, 'export_leg1_transfers:noTransfers'), 'none accepted: expected noTransfers, got "%s"', id);
fprintf('  none accepted: refused (%s)\n', id);

% ---- 3. refusals outside results/ ----
sib = [resDir '_sibling_' char(java.util.UUID.randomUUID())];
bad = {fullfile(tmp, 'sub', 'x.mat'), fullfile(resDir, '..', 'escape_leg1.mat'), fullfile(sib, 'x.mat')};
for kb = 1:numel(bad)
    id = '';
    try
        export_leg1_transfers(resultsMat, struct('outFile', bad{kb}, 'verbose', false));
    catch err
        id = err.identifier;
    end
    assert(strcmp(id, 'export_leg1_transfers:outFile'), '%s: expected outFile refusal, got "%s"', bad{kb}, id);
    assert(~isfile(bad{kb}) && ~isfolder(fileparts(bad{kb})) || kb == 2, '%s: something was created', bad{kb});
    assert(~isfile(fullfile(root, 'escape_leg1.mat')), 'escape file created');
    fprintf('  refused: %s\n', bad{kb});
end
fprintf('  PASS\n');
ok = true;
end

% ------------------------------------------------------------------------
function delete_if(f)
% DELETE_IF  Delete a file if it exists.  INPUTS: f char.  OUTPUTS: none.
if isfile(f), delete(f); end
end
