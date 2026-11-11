function ok = test_leg1_package_for()
%% Purpose:
%
%   Verifies leg1_package_for (and file_sha256) on SYNTHETIC package files in
%   a temporary folder: (1) the primary path matches -> it is returned
%   (source 'primary'); (2) the primary path holds ANOTHER package and a
%   snapshot folder holds a copy of the wanted one under an arbitrary name
%   (no "pilot", no sha8 in it) -> the snapshot is returned (source
%   'snapshot'); (3) nothing matches -> empty path and a reason naming the
%   wanted SHA-256 and every file checked; (4) HASH BEFORE LOAD: a matching
%   file that is NOT a MAT-file is returned without error, and a corrupt
%   .mat in the snapshot folder is hashed, never loaded; (5) a missing
%   primary file falls back to the snapshot; (6) snapshotsOnly skips the
%   primary; (6b) leg1_snapshot_package copies the package to
%   <dir>/costate_catalog_gto_dro_<sha8>.mat, reports 'existing' for any
%   file in <dir> with the wanted hash (any name), and never overwrites a
%   sha8-named file whose hash differs ('conflict'); (7) the saved leg-1 results, when present: their package is
%   found (default path or snapshot in results/inputs).
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 true if every assertion
%                                                   held (else it throws)
%
%% Revision History:
%  M. Casey  final fix wave: package lookup by SHA-256       (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
run(fullfile(root, 'setup_lunar_pnt_paths.m'));
tmp = tempname;
mkdir(tmp);
cleanTmp = onCleanup(@() rmdir(tmp, 's'));
snapDir = fullfile(tmp, 'inputs');
mkdir(snapDir);

pkgA = fullfile(tmp, 'pkgA.mat');  pkgVar = magic(4);      save(pkgA, 'pkgVar');
pkgB = fullfile(tmp, 'pkgB.mat');  pkgVar = magic(5);      save(pkgB, 'pkgVar');
shaA = file_sha256(pkgA);
shaB = file_sha256(pkgB);
assert(numel(shaA) == 64 && ~strcmp(shaA, shaB), 'file_sha256: bad hashes');
idA = struct('sha256', shaA, 'file', pkgA);

% ---- 1. primary matches ----
[f, why, info] = leg1_package_for(idA, struct('primary', pkgA, 'inputsDirs', {{snapDir}}));
assert(strcmp(f, pkgA) && strcmp(info.source, 'primary') && isempty(why), '1: got "%s" (%s)', f, why);
fprintf('  primary matches: %s (%s)\n', f, info.source);

% ---- 2. primary replaced, snapshot under an arbitrary name ----
snapA = fullfile(snapDir, 'anything_at_all.mat');
copyfile(pkgA, snapA);
[f, ~, info] = leg1_package_for(idA, struct('primary', pkgB, 'inputsDirs', {{snapDir}}));
assert(strcmp(f, snapA) && strcmp(info.source, 'snapshot'), '2: got "%s"', f);
fprintf('  primary replaced: snapshot %s found by its SHA-256\n', f);

% ---- 3. nothing matches ----
delete(snapA);
copyfile(pkgB, fullfile(snapDir, 'other.mat'));
[f, why, info] = leg1_package_for(idA, struct('primary', pkgB, 'inputsDirs', {{snapDir}}));
assert(isempty(f) && contains(why, shaA(1:12)) && size(info.checked, 1) == 2, ...
       '3: got "%s", why "%s", %d checked', f, why, size(info.checked, 1));
fprintf('  nothing matches: "%s"\n', why);

% ---- 4. hash before load ----
txt = fullfile(snapDir, 'not_a_mat_file.mat');
fid = fopen(txt, 'w');  fprintf(fid, 'this is not a MAT-file');  fclose(fid);
bad = fullfile(snapDir, 'corrupt.mat');
fid = fopen(bad, 'w');  fwrite(fid, uint8(1:200));  fclose(fid);
idT = struct('sha256', file_sha256(txt), 'file', pkgB);
[f, ~, info] = leg1_package_for(idT, struct('primary', pkgB, 'inputsDirs', {{snapDir}}));
assert(strcmp(f, txt) && strcmp(info.source, 'snapshot'), '4: got "%s"', f);
[f, why] = leg1_package_for(idA, struct('primary', pkgB, 'inputsDirs', {{snapDir}}));
assert(isempty(f) && ~isempty(why), '4: a corrupt file was not simply hashed');
fprintf('  hash before load: a non-MAT file with the wanted hash is returned; a corrupt .mat is only hashed\n');

% ---- 5. missing primary ----
copyfile(pkgA, fullfile(snapDir, 'zz_copy.mat'));
[f, ~, info] = leg1_package_for(idA, struct('primary', fullfile(tmp, 'gone.mat'), 'inputsDirs', {{snapDir}}));
assert(strcmp(f, fullfile(snapDir, 'zz_copy.mat')) && strcmp(info.source, 'snapshot'), '5: got "%s"', f);
fprintf('  missing primary: snapshot used\n');

% ---- 6. snapshotsOnly ----
copyfile(pkgA, fullfile(tmp, 'pkgA2.mat'));
[f, ~, info] = leg1_package_for(idA, struct('primary', fullfile(tmp, 'pkgA2.mat'), ...
                                            'inputsDirs', {{snapDir}}, 'snapshotsOnly', true));
assert(strcmp(f, fullfile(snapDir, 'zz_copy.mat')) && strcmp(info.source, 'snapshot'), '6: got "%s"', f);
fprintf('  snapshotsOnly: the primary is skipped\n');

% ---- 6b. leg1_snapshot_package: copied, existing (any name), never overwritten ----
sd = fullfile(tmp, 'snap2');
idA2 = struct('sha256', shaA, 'file', pkgA);
[sf, st] = leg1_snapshot_package(idA2, sd);
want = fullfile(sd, sprintf('costate_catalog_gto_dro_%s.mat', shaA(1:8)));
assert(strcmp(st, 'copied') && strcmp(sf, want) && strcmp(file_sha256(sf), shaA), 'snapshot: %s %s', st, sf);
[sf2, st2] = leg1_snapshot_package(idA2, sd);
assert(strcmp(st2, 'existing') && strcmp(sf2, want) && isscalar(dir(fullfile(sd, '*.mat'))), ...
       'snapshot second call: %s %s', st2, sf2);
movefile(want, fullfile(sd, 'renamed_by_hand.mat'));
[sf3, st3] = leg1_snapshot_package(idA2, sd);
assert(strcmp(st3, 'existing') && strcmp(sf3, fullfile(sd, 'renamed_by_hand.mat')) && ...
       isscalar(dir(fullfile(sd, '*.mat'))), 'snapshot by hash, any name: %s %s', st3, sf3);
delete(fullfile(sd, 'renamed_by_hand.mat'));
copyfile(pkgB, want);                                  % a sha8-named file that is NOT the package
[sf4, st4] = leg1_snapshot_package(idA2, sd);
assert(strcmp(st4, 'conflict') && isempty(sf4) && strcmp(file_sha256(want), shaB), ...
       'snapshot conflict: %s "%s" (file overwritten: %d)', st4, sf4, ~strcmp(file_sha256(want), shaB));
fprintf('  leg1_snapshot_package: copied, then existing, found by hash under any name, a sha8-named stranger never overwritten\n');

% ---- 7. the saved results ----
r1 = fullfile(root, 'results', 'lunar_pnt_leg1_study.mat');
if isfile(r1)
    R = load(r1, 'run_meta');
    [f, why, info] = leg1_package_for(R.run_meta.identity);
    assert(~isempty(f), 'saved results: %s', why);
    fprintf('  saved results (package %s...): %s (%s)\n', R.run_meta.identity.sha256(1:8), f, info.source);
else
    fprintf('  saved results: not present, skipped\n');
end
fprintf('  PASS\n');
ok = true;
end
