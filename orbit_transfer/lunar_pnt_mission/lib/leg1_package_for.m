function [pkgFile, why, info] = leg1_package_for(identity, opts)
%% Purpose:
%
%   THE GTO -> DRO PACKAGE FILE THAT A SAVED LEG-1 RESULT CAME FROM, found
%   by its SHA-256 alone (the package at the default path is replaced by
%   the fleet's; the leg-1 results stay tied to the package they were
%   seeded from). Candidates, in order, each HASHED BEFORE anything loads
%   it (a candidate is never loaded or validated here):
%     1. the primary path: opts.primary, else identity.file, else the
%        default GTO_DRO/results/costate_catalog_gto_dro.mat (skipped with
%        opts.snapshotsOnly);
%     2. opts.extra files (e.g. a snapshot path recorded in run_meta);
%     3. every *.mat in opts.inputsDirs (default: this campaign's
%        results/inputs), sorted by name; the file NAME is not used (a
%        snapshot is any file whose own SHA-256 equals the wanted one).
%   The first whose SHA-256 equals identity.sha256 is returned. None ->
%   pkgFile = '' and why says which hash was wanted and what was checked.
%
%   Self-demo (nargin == 0): the package of the saved leg-1 results.
%
%% Inputs:
%
%  identity                 struct                  .sha256 char [1 x 64]
%                                                   (run_meta.identity of
%                                                   the results); .file char
%                                                   (optional)
%
%  opts                     struct (optional)       .primary char ['']
%                                                   .extra {1 x m} cell [{}]
%                                                   .inputsDirs {1 x n} cell
%                                                   [{results/inputs}]
%                                                   .snapshotsOnly [false]
%
%% Outputs:
%
%  pkgFile                  char                    matching file ('' = none)
%  why                      char                    '' when found, else the
%                                                   reason
%  info                     struct                  .wantSha char; .source
%                                                   ('primary' | 'extra' |
%                                                   'snapshot' | ''); .checked
%                                                   {k x 2} cell (file, sha;
%                                                   sha '' if missing)
%
%% Revision History:
%  M. Casey  final fix wave: package lookup by SHA-256       (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(fileparts(mfilename('fullpath')));        % lunar_pnt_mission/
if nargin == 0
    run(fullfile(here, 'setup_lunar_pnt_paths.m'));
    R = load(fullfile(here, 'results', 'lunar_pnt_leg1_study.mat'), 'run_meta');
    [pkgFile, why, info] = leg1_package_for(R.run_meta.identity);
    fprintf('leg1_package_for self-demo: wanted %s...; found "%s" (%s)%s\n', ...
            info.wantSha(1:12), pkgFile, info.source, why);
    return
end
if nargin < 2 || isempty(opts), opts = struct(); end
known = {'primary', 'extra', 'inputsDirs', 'snapshotsOnly'};
for f = fieldnames(opts).'
    if ~any(strcmp(f{1}, known))
        error('leg1_package_for:opts', 'unknown option "%s"', f{1});
    end
end
if ~(isstruct(identity) && isfield(identity, 'sha256') && numel(identity.sha256) == 64)
    error('leg1_package_for:identity', 'identity must carry a 64-character sha256');
end
want = lower(char(identity.sha256));
defPkg = fullfile(fileparts(here), 'GTO_DRO', 'results', 'costate_catalog_gto_dro.mat');
primary = fdef(opts, 'primary', '');
if isempty(primary)
    if isfield(identity, 'file') && ~isempty(identity.file)
        primary = char(identity.file);
    else
        primary = defPkg;
    end
end
extra   = cellstr(fdef(opts, 'extra', {}));
inDirs  = cellstr(fdef(opts, 'inputsDirs', {fullfile(here, 'results', 'inputs')}));
snapOnly = fdef(opts, 'snapshotsOnly', false);

% candidates, in order
if snapOnly, pFiles = {}; else, pFiles = {primary}; end
eFiles = extra(~cellfun(@isempty, extra));
sFiles = cellfun(@snapshots_in, inDirs(:).', 'UniformOutput', false);
sFiles = [sFiles{:}];
cand = [pFiles, eFiles, sFiles];
srcs = [repmat({'primary'}, 1, numel(pFiles)), repmat({'extra'}, 1, numel(eFiles)), ...
        repmat({'snapshot'}, 1, numel(sFiles))];

pkgFile = '';  why = '';  src = '';
checked = cell(numel(cand), 2);  nC = 0;
for kk = 1:numel(cand)
    fc = cand{kk};
    if any(strcmp(checked(1:nC, 1), fc)), continue; end     % hashed already
    if isfile(fc), sh = file_sha256(fc); else, sh = ''; end
    nC = nC + 1;
    checked(nC, :) = {fc, sh};
    if strcmp(sh, want)
        pkgFile = fc;  src = srcs{kk};
        break
    end
end
checked = checked(1:nC, :);
if isempty(pkgFile)
    lst = strjoin(cellfun(@(a, b) sprintf('%s (%s)', a, short_sha(b)), checked(:, 1), checked(:, 2), ...
                          'UniformOutput', false), '; ');
    why = sprintf('no package with SHA-256 %s... found; checked: %s', want(1:12), lst);
end
info = struct('wantSha', want, 'source', src, 'checked', {checked});
end

% ------------------------------------------------------------------------
function c = snapshots_in(d)
% SNAPSHOTS_IN  Every *.mat file in a folder, sorted by name ({} if the
%   folder is missing).  INPUTS: d char.  OUTPUTS: c {1 x n} cell.
c = {};
if isempty(d) || ~isfolder(d), return; end
dd = dir(fullfile(d, '*.mat'));
dd = dd(~[dd.isdir]);
if isempty(dd), return; end
c  = sort(fullfile({dd.folder}, {dd.name}));
c  = reshape(c, 1, []);
end

% ------------------------------------------------------------------------
function s = short_sha(h)
% SHORT_SHA  First 8 hex of a hash, or 'missing'.  INPUTS: h char.
%   OUTPUTS: s char.
if isempty(h), s = 'missing'; else, s = h(1:8); end
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
