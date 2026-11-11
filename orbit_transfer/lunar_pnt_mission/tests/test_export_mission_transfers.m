function [ok, skipped] = test_export_mission_transfers()
%% Purpose:
%
%   Verifies export_mission_transfers on the SAVED study results (nothing
%   re-polished): (0) catalog_identity reproduces the saved fingerprint's
%   identity fields (the one recipe the study and the export share);
%   (1) refusals: a non-accepted transfer (notAccepted) and a catalog
%   fingerprint mismatch (fingerprint), on tampered copies of the results;
%   an output path inside the kit folder (outFile), given absolute into a
%   new subfolder, relative from the kit's parent, and through '..', with
%   nothing created; (2) the sidecar holds ONE variable, mission_transfers,
%   whose .entries carry EXACTLY the field names of the library's
%   cat_.alternatives, in order, with the same class and shape; n = 6;
%   (3) HEADLINE PIN: every entry was certified by the library's own
%   certifier (statusRun 'certify_root (full)', certifyDz < 1e-6) and is
%   status 4; each carries the certifier's z8, its [14 x 24]
%   multiple-shooting junctions (column 1 = [rv0; 1; z8(1:7)], 1e-14) and
%   its flown miss (the library's convention), with the polish's root kept
%   in .mission.polishZ8; sD / sA reproduce the polish's rv0 / rvf through
%   the campaign's handles (1e-9); (4) the library junction convention on a
%   primary: the single-flight cut used for FALLBACK entries reproduces the
%   stored column 1 to 1e-14 and every column to 1e-5 relative (measured
%   3.3e-7 on (16,31), 7.1e-7 on (17,40), 3.0e-6 on (3,5)); (5) the
%   fullStatus = false fallback: status 1 by the library's rule, reason
%   naming the missing checks, the polish's z8 and the single-flight cut,
%   and the session path restored exactly (exercised with DRO_tulip/
%   indirect removed first, so the export has to add it); (6) DELIBERATE
%   FAILURE of the root-reproduction gate: one P(3).z8 costate moved by
%   1e-5 (sc 3 only, one certifier run) -> the certifier converges back to
%   the true root, 1e-5 from the stored z8, so the entry is a FALLBACK
%   (status 1) whose label and reason say so; (7) the KIT ROUND TRIP:
%   catalog_transfer(cat_, entry) lands on the tulip by both methods below
%   0.05 km (the library's gate is 100 km; measured max 2.5e-3 km, so this
%   catches an endpoint or spline regression with about a 10x margin);
%   (8) the shipped library file's MD5 is unchanged. Returns skipped = true
%   (and ok = true) with a message when the results file does not exist.
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
%  M. Casey                                                   (c) 11/10/2026
%  M. Casey  fix round 1: status-4 pin, root-reproduction failure,
%            certifier junctions, outFile refusals, path restore,
%            kit gate 0.05 km                                   (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
run(fullfile(root, 'setup_lunar_pnt_paths.m'));
resultsMat = fullfile(root, 'results', 'lunar_pnt_mission_study.mat');
skipped = false;
if ~isfile(resultsMat)
    fprintf('  SKIPPED: no results file (%s); run lunar_pnt_mission_study first\n', resultsMat);
    ok = true;  skipped = true;
    return
end
tmp = tempname;
mkdir(tmp);
cleanTmp = onCleanup(@() rmdir(tmp, 's'));
kit = fullfile(fileparts(root), 'DRO_tulip', 'deliverables', 'costate_library_70mN_user_deliverable9');
libFile = fullfile(kit, 'costate_catalog_dro_tulip_70mN.mat');
md5Before = file_md5(libFile);
M = mission_inputs();
cat_ = M.cat;
R = load(resultsMat, 'P', 'S', 'run_meta', 'fingerprint');
kmGate = 0.05;                                 % kit round trip [km], both methods

% ---- 0. one fingerprint recipe ----
fpI = catalog_identity(M);
fnI = fieldnames(fpI);
assert(all(cellfun(@(f) isequal(R.fingerprint.(f), fpI.(f)), fnI)), ...
       'catalog_identity does not reproduce the saved fingerprint');
fprintf('  catalog_identity reproduces the saved fingerprint (%d identity fields)\n', numel(fnI));

% ---- 1. refusals (cheap: before any status run) ----
bad = R;  bad.P(2).ok = false;  bad.P(2).reason = 'gate failed: conj (tampered)';
badMat = fullfile(tmp, 'tampered_ok.mat');
save(badMat, '-struct', 'bad');
thrown = refusal_id(badMat, struct('outFile', fullfile(tmp, 'never.mat')));
assert(strcmp(thrown, 'export_mission_transfers:notAccepted'), 'non-accepted transfer not refused (%s)', thrown);
assert(~isfile(fullfile(tmp, 'never.mat')), 'a refused export wrote its file');
fprintf('  refusal: non-accepted transfer -> %s (nothing written)\n', thrown);
bad = R;  bad.fingerprint.z8md5 = repmat('0', 1, 32);
badMat = fullfile(tmp, 'tampered_fp.mat');
save(badMat, '-struct', 'bad');
thrown = refusal_id(badMat, struct('outFile', fullfile(tmp, 'never.mat')));
assert(strcmp(thrown, 'export_mission_transfers:fingerprint'), 'fingerprint mismatch not refused (%s)', thrown);
fprintf('  refusal: catalog fingerprint mismatch -> %s\n', thrown);
[kitParent, kitName] = fileparts(kit);
od = pwd;
restoreCd = onCleanup(@() cd(od));
cd(kitParent);
targets = {fullfile(kit, 'newsub', 'x.mat'), fullfile(kitName, 'newsub', 'x.mat'), ...
           fullfile('..', 'deliverables', kitName, 'x.mat'), fullfile(kit, '..', kitName, 'sub2', 'y.mat')};
for kt = 1:numel(targets)
    thrown = refusal_id(resultsMat, struct('outFile', targets{kt}));
    assert(strcmp(thrown, 'export_mission_transfers:outFile'), 'output %s not refused (%s)', targets{kt}, thrown);
end
cd(od);
assert(~isfolder(fullfile(kit, 'newsub')) && ~isfile(fullfile(kit, 'x.mat')), 'a refused output created something');
fprintf('  refusal: %d output paths inside the kit (absolute new subfolder, relative, ''..'') -> outFile, nothing created\n', ...
        numel(targets));

% ---- 2. the full export: one variable, the library's field names and types ----
outFile = fullfile(tmp, 'mission_transfers_70mN.mat');
p0 = path;
t0 = tic;
E = export_mission_transfers(resultsMat, struct('outFile', outFile, 'verbose', false));
wallFull = toc(t0);
assert(isequal(path, p0), 'the full export changed the session path');
w = whos('-file', outFile);
assert(isscalar(w) && strcmp(w.name, 'mission_transfers'), 'sidecar variables: %s', strjoin({w.name}, ', '));
L = load(outFile);
X = L.mission_transfers;
assert(isequal(X.entries, E.entries), 'returned struct differs from the saved one');
A1 = cat_.alternatives(1);
fnLib = fieldnames(A1);
assert(isequal(fieldnames(X.entries), fnLib), 'entry field names differ from cat_.alternatives');
n = numel(X.entries);
assert(n == 6 && isequal(size(X.entries), [1 6]), 'expected [1 x 6] entries, got [%s]', num2str(size(X.entries)));
check_types(X.entries, A1);
for f = {'mission', 'catalog', 'provenance', 'usage'}
    assert(isfield(X, f{1}), 'sidecar has no .%s', f{1});
end
assert(isequal(X.catalog.fingerprint, R.fingerprint), '.catalog.fingerprint is not the results file''s');
assert(isequal(X.provenance.run_meta, R.run_meta), '.provenance.run_meta is not the results file''s');
assert(contains(X.provenance.note, 'NOT part of the library of record') && isequal(X.provenance.fenced, false) && ...
       contains(X.provenance.note, 'UNFENCED'), 'provenance does not disclose the sidecar status / the unfenced certifier');
fprintf('  sidecar: 1 variable, %d entries, %d library fields matched by name, class and shape (%.0f s, full status); path unchanged\n', ...
        n, numel(fnLib), wallFull);

% ---- 3. headline pin, the certifier's z8 / junctions / flight, endpoints ----
P = R.P;
for ke = 1:n
    e = X.entries(ke);  q = X.mission(ke);  k = q.sc;
    assert(strcmp(q.statusRun, 'certify_root (full)') && q.certifyDz < 1e-6 && e.status == 4 && ...
           strcmp(e.status_reason, 'full stack passed') && e.conj == 1 && strcmp(e.conjVerdict, 'PASS') && ...
           e.h6Margin > 1 && e.liftMargin > 10 && ~e.inferred && ~e.borderline, ...
           'entry %d: not a certified status 4 (%s, dz %.2e, status %g)', ke, q.statusRun, q.certifyDz, e.status);
    assert(e.iD == 0 && e.iA == 0 && e.sheet == 1 && q.slot == R.S.assign(k), 'entry %d: iD/iA/sheet/slot', ke);
    assert(max(abs(M.depState(e.sD) - P(k).rv0)) < 1e-9, 'entry %d: depState(sD) ~= rv0', ke);
    assert(max(abs(M.arrState(e.sA) - P(k).rvf)) < 1e-9, 'entry %d: arrState(sA) ~= rvf', ke);
    assert(isequal(q.polishZ8, P(k).z8) && abs(max(abs(e.z8 - P(k).z8)) - q.certifyDz) < 1e-15 && ...
           e.tf_nd == e.z8(8), 'entry %d: z8 is not the certifier''s root / polishZ8 not the polish''s', ke);
    assert(isequal(size(e.junctions), [14 24]) && max(abs(e.junctions(:, 1) - [P(k).rv0; 1; e.z8(1:7)])) < 1e-14, ...
           'entry %d: certifier junctions shape or column 1', ke);
    assert(e.flyKm == q.certify.flyKm && e.flyVms == q.certify.flyVms && e.flyKm < 100 && e.flyVms < 10, ...
           'entry %d: flown miss is not the certifier''s', ke);
    assert(q.reused == strcmp(R.run_meta.polishSource(k).origin, 'reused'), 'entry %d: reused flag', ke);
end
fprintf('  status: [%s], all certify_root (full), certifyDz max %.1e; certifier z8 + [14 x 24] junctions + miss stored\n', ...
        sprintf('%d ', [X.entries.status]), max([X.mission.certifyDz]));

% ---- 4. the library junction convention (the fallback cut), on a library primary ----
sh = M.sh;
cellLib = [16 31];
kL = sh.entry_index(cellLib(1), cellLib(2));
J  = sh.junctions{kL};
sd = seed_from_z8(sh.z8(:, kL), M.depState(M.sD(cellLib(1))), 24, M.Tnd, M.cnd, M.muStar);
Jc = sd.Y(:, 1:24);
c1 = max(abs(J(:, 1) - Jc(:, 1)));
rel = max(max(abs(J - Jc))./max(abs(J)));
assert(isequal(size(J), [14 24]) && c1 < 1e-14 && rel < 1e-5, ...
       'library convention: column 1 diff %.2e, relative %.2e', c1, rel);
fprintf('  fallback cut vs library primary (%d,%d): column 1 diff %.1e, max relative column diff %.2e (limit 1e-5)\n', ...
        cellLib, c1, rel);

% ---- 5. fallback, with the certifier's folder removed so the export must add it ----
indirect = fullfile(fileparts(root), 'DRO_tulip', 'indirect');
pFull = path;
restorePath = onCleanup(@() path(pFull));
if contains([pathsep path pathsep], [pathsep indirect pathsep]), rmpath(indirect); end
pBefore = path;
assert(isempty(which('optimality_status')), 'optimality_status still resolves with DRO_tulip/indirect removed');
Ef = export_mission_transfers(resultsMat, struct('outFile', fullfile(tmp, 'fallback.mat'), ...
                                                 'fullStatus', false, 'verbose', false));
assert(isequal(path, pBefore), 'the fallback export did not restore the path exactly');
path(pFull);
for ke = 1:n
    e = Ef.entries(ke);  k = Ef.mission(ke).sc;
    assert(e.status == 1 && contains(e.status_reason, 'neither') && contains(e.status_reason, 'dense conjugate scan') && ...
           contains(e.status_reason, 'lift_margin') && strcmp(Ef.mission(ke).statusRun, 'fallback: fullStatus = false'), ...
           'fallback entry %d: status %g reason %s', ke, e.status, e.status_reason);
    sdk = seed_from_z8(P(k).z8, P(k).rv0, 24, M.Tnd, M.cnd, M.muStar);
    assert(isequal(e.z8, P(k).z8) && isequal(e.junctions, sdk.Y(:, 1:24)) && e.flyKm == P(k).flight.flyKm, ...
           'fallback entry %d: not the polish''s z8 / single-flight cut / miss', ke);
end
check_types(Ef.entries, A1);
fprintf('  fallback (fullStatus false): status 1 on all %d, reason names the checks not run; polish z8 + cut; path restored exactly\n', n);

% ---- 6. deliberate failure of the root-reproduction gate ----
bad = R;  bad.P(3).z8(4) = bad.P(3).z8(4) + 1e-5;
badMat = fullfile(tmp, 'tampered_z8.mat');
save(badMat, '-struct', 'bad');
t0 = tic;
Et = export_mission_transfers(badMat, struct('outFile', fullfile(tmp, 'tampered_out.mat'), 'only', 3, 'verbose', false));
et = Et.entries(1);  qt = Et.mission(1);
assert(isscalar(Et.entries) && qt.sc == 3, 'only = 3 did not export exactly sc 3');
assert(startsWith(qt.statusRun, 'fallback: certify_root converged to a root') && qt.certifyDz > 1e-6 && ...
       et.status == 1 && contains(et.status_reason, 'converged to a root') && isequal(et.z8, bad.P(3).z8), ...
       'tampered z8: status %g, run "%s", dz %.2e', et.status, qt.statusRun, qt.certifyDz);
fprintf('  root-reproduction failure: z8(4) + 1e-5 -> certifier root %.2e away -> "%s", status %d (%.0f s)\n', ...
        qt.certifyDz, qt.statusRun, et.status, toc(t0));

% ---- 7. kit round trip ----
kmZ = nan(1, n);  kmJ = nan(1, n);
for ke = 1:n
    Fz = catalog_transfer(cat_, X.entries(ke), struct('method', 'z8'));
    Fj = catalog_transfer(cat_, X.entries(ke), struct('method', 'junctions'));
    kmZ(ke) = Fz.flyKm;  kmJ(ke) = Fj.flyKm;
    assert(strcmp(Fz.method, 'z8') && strcmp(Fj.method, 'junctions'), 'entry %d: kit method', ke);
end
fprintf('  kit round trip flyKm, method z8:        [%s] km (limit %g; library gate 100)\n', sprintf('%.3g ', kmZ), kmGate);
fprintf('  kit round trip flyKm, method junctions: [%s] km (limit %g)\n', sprintf('%.3g ', kmJ), kmGate);
assert(all(kmZ < kmGate), 'kit z8 flight misses %.3g km', max(kmZ));
assert(all(kmJ < kmGate), 'kit junction flight misses %.3g km', max(kmJ));

% ---- 8. the library file is never written ----
md5After = file_md5(libFile);
assert(strcmp(md5Before, md5After), 'the shipped library file changed');
fprintf('  library file MD5 unchanged (%s)\n', md5After);
ok = true;
clear cleanTmp restoreCd restorePath
end

% ------------------------------------------------------------------------
function check_types(entries, A1)
% CHECK_TYPES  Every field of every entry has the library's class and
%   shape (char fields: a row or empty, as the library stores '' too).
%   INPUTS: entries [1 x n] struct; A1 cat_.alternatives(1).  OUTPUTS: none.
fnLib = fieldnames(A1);
for kf = 1:numel(fnLib)
    f = fnLib{kf};
    for ke = 1:numel(entries)
        v = entries(ke).(f);
        assert(strcmp(class(v), class(A1.(f))), 'entry %d field %s class %s, library %s', ke, f, class(v), class(A1.(f)));
        if ischar(A1.(f))
            assert(isrow(v) || isempty(v), 'entry %d field %s is not a char row', ke, f);
        else
            assert(isequal(size(v), size(A1.(f))), 'entry %d field %s size [%s], library [%s]', ke, f, ...
                   num2str(size(v)), num2str(size(A1.(f))));
        end
    end
end
assert(all(arrayfun(@(e) ~isempty(e.source) && ~isempty(e.status_reason), entries)), 'empty source / status_reason');
end

% ------------------------------------------------------------------------
function id = refusal_id(resultsMat, opts)
% REFUSAL_ID  The error identifier export_mission_transfers throws ('' if none).
%   INPUTS: resultsMat char; opts struct.  OUTPUTS: id char.
id = '';
try
    export_mission_transfers(resultsMat, opts);
catch err
    id = err.identifier;
end
end

% ------------------------------------------------------------------------
function h = file_md5(file)
% FILE_MD5  Hex MD5 of a file's bytes.
%   INPUTS: file char.  OUTPUTS: h char [1 x 32].
fid = fopen(file, 'r');
assert(fid > 0, 'cannot open %s', file);
b = fread(fid, inf, '*uint8');
fclose(fid);
md = java.security.MessageDigest.getInstance('MD5');
md.update(b);
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end
