function [snapFile, status, msg] = leg1_snapshot_package(identity, snapDir)
%% Purpose:
%
%   KEEP A COPY OF THE GTO -> DRO PACKAGE A LEG-1 RUN USED, so the saved
%   results stay checkable after the package at the default path is
%   replaced. Any *.mat in snapDir whose own SHA-256 equals identity.sha256
%   IS the snapshot (the file name is not used): status 'existing'.
%   Otherwise identity.file is copied to
%   snapDir/costate_catalog_gto_dro_<first 8 hex of the SHA-256>.mat and the
%   copy's hash verified: status 'copied'. A file of that name that hashes
%   differently is NEVER overwritten: status 'conflict', snapFile ''.
%
%   Self-demo (nargin == 0): the package of the saved leg-1 results against
%   results/inputs (finds the existing snapshot; copies nothing new unless
%   the package is missing there).
%
%% Inputs:
%
%  identity                 struct                  .sha256 char [1 x 64]
%                                                   .file char (the package
%                                                   the run loaded)
%  snapDir                  char                    snapshot folder (created
%                                                   if needed)
%
%% Outputs:
%
%  snapFile                 char                    snapshot path ('' on
%                                                   conflict)
%  status                   char                    'existing' | 'copied' |
%                                                   'conflict'
%  msg                      char                    one line for the log
%
%% Revision History:
%  M. Casey  final fix wave: automatic package snapshot      (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    run(fullfile(here, 'setup_lunar_pnt_paths.m'));
    R = load(fullfile(here, 'results', 'lunar_pnt_leg1_study.mat'), 'run_meta');
    [snapFile, status, msg] = leg1_snapshot_package(R.run_meta.identity, fullfile(here, 'results', 'inputs'));
    fprintf('leg1_snapshot_package self-demo: %s\n', msg);
    return
end
want = lower(char(identity.sha256));
[snapFile, ~, info] = leg1_package_for(identity, struct('snapshotsOnly', true, 'inputsDirs', {{snapDir}}));
if ~isempty(snapFile)
    status = 'existing';
    msg = sprintf('package snapshot: existing %s (SHA-256 %s... matches; %d file(s) hashed)', snapFile, ...
                  want(1:8), size(info.checked, 1));
    return
end
target = fullfile(snapDir, sprintf('costate_catalog_gto_dro_%s.mat', want(1:8)));
if isfile(target)
    snapFile = '';
    status = 'conflict';
    msg = sprintf(['package snapshot: CONFLICT, %s exists but its SHA-256 is %s..., not %s...; ' ...
                   'not overwritten, no snapshot recorded'], target, file_sha256(target), want(1:12));
    return
end
src = char(identity.file);
if ~isfile(src) || ~strcmp(file_sha256(src), want)
    error('leg1_snapshot_package:source', 'the package file %s is missing or no longer has SHA-256 %s...', ...
          src, want(1:12));
end
if ~isfolder(snapDir), mkdir(snapDir); end
copyfile(src, target);
if ~strcmp(file_sha256(target), want)
    delete(target);
    error('leg1_snapshot_package:copy', 'the copy %s does not hash to %s...; deleted', target, want(1:12));
end
snapFile = target;
status = 'copied';
msg = sprintf('package snapshot: copied %s -> %s (SHA-256 %s... verified)', src, target, want(1:8));
end
