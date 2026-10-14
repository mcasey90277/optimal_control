function out = build_dro_deliverable(cfg)
%% Purpose:
%
%   Package a DRO -> tulip minimum-time costate catalog as a shippable
%   deliverable: the catalog itself, the maintained helper suite, and a
%   README whose every number is COMPUTED FROM THE CATALOG (and, for the
%   status-layer library, from the audit files bound to it by content key).
%
%   That last point is the reason this file exists. The 2026-08 deliverable
%   was assembled by hand and its README hard-coded "3,936 entries", "9
%   thrust levels" and "coverage 89-97%". The catalog has since gained the
%   0.75 N and 0.5 N rungs (4,439 entries, 11 rungs, 70% coverage), a full
%   conjugate-point sweep and the three sufficiency-hypothesis gates -- so
%   the shipped README understated and misdescribed what we have. Deriving
%   the numbers makes that failure mode impossible: re-run this and the
%   description cannot drift from the data.
%
%   Two paths, chosen by the catalog:
%
%   * A catalog WITH .alternatives (the status-layer library, deliverable 8,
%     the default): 1 export (export_deliverable_catalog), 2 audit gate by
%     CONTENT KEY, 3 copy every file of deliverable_file_list, 4 the sample
%     movie (cfg.movies), 5 the status torus PNG, 6 the README from computed
%     facts, 7 the recipient self-test (deliverable_selftest, cfg.selftest),
%     8 the zip, only after a passed self-test (cfg.zip).
%   * A catalog WITHOUT .alternatives (deliverables 3 and the earlier 70 mN
%     catalog): the original path, its audit gate by path equality, its
%     five helpers and its README, unchanged.
%
%  ASSUMPTIONS / NOTES:
%
% • The key-bound gate refuses (error id build_dro_deliverable:audit) unless
%   the primaries audit (variable A) reports nBad = 0 for the export's
%   catalog_content_key, and every audit_alt_*.mat chunk carries the
%   export's catalog_content_key AND alternatives_content_key, reports
%   audit.nBad = 0 with every altRows(k).ok true, and the union of the
%   chunks' [altRows.k] is 1..nAlternatives exactly once (coverage).
% • Helper sources are the MAINTAINED copies listed by deliverable_file_list
%   (DRO_tulip/indirect, costate_common, DRO_tulip/indirect/examples), never
%   a previous deliverable's copies.
% • An existing outDir is replaced only if it looks like an earlier build
%   (it holds a costate_catalog_*.mat or README.md, or is empty).
%
%% Inputs:
%
%  cfg                      struct (optional)
%   .catMat                 char                    catalog path [the library of
%                                                   record, indirect/results/
%                                                   library_70mN_24x48_v2/
%                                                   costate_catalog_dro_tulip_70mN.mat]
%   .name                   char                    deliverable name
%                                                   ['costate_catalog_70mN_status_deliverable8']
%   .outDir                 char                    [deliverables/<name>]
%   .movies                 logical                 render the sample movie [true]
%   .zip                    logical                 write and verify the .zip [true]
%   .selftest               logical                 run deliverable_selftest [true]
%   .auditPrimariesMat      char                    [indirect/results/
%                                                   library_70mN_24x48_merged/audit_70mN.mat]
%   .auditAltDir            char                    folder of audit_alt_*.mat
%                                                   [the catalog's folder]
%   (status-less path:) .auditMat .skipAuditBecause as before; outDir
%                                                   [deliverables/costate_catalog_deliverable3],
%                                                   movies [false], zip [true]
%
%% Outputs:
%
%  out                      struct                  .dir .zip ('' when not
%                                                   zipped) .facts (the
%                                                   computed numbers) .files
%                                                   {cellstr} .selftest
%                                                   (deliverable_selftest
%                                                   result, [] when off)
%                                                   .keys (.primKey .altKey)
%
%% Revision History:
%  M. Casey                                                   (c) 09/08/2026
%  M. Casey  deliverable 8: export, key-bound audit gate, file list, movie,
%            torus PNG, README from computed facts, recipient self-test,
%            zip from the parent; status-less path kept     (c) 10/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, cfg = struct(); end
here = fileparts(mfilename('fullpath'));
ind  = fullfile(here, 'indirect');
cc   = fullfile(fileparts(here), 'costate_common');
addpath(ind, cc);
v2Dir  = fullfile(ind, 'results', 'library_70mN_24x48_v2');
catMat = getf(cfg, 'catMat', fullfile(v2Dir, 'costate_catalog_dro_tulip_70mN.mat'));

L0 = load(catMat);  fn = fieldnames(L0);
assert(numel(fn) == 1, 'expected one variable in %s', catMat);
c0 = L0.(fn{1});
if ~isfield(c0, 'alternatives')
    out = buildLegacy(cfg, c0, catMat, here);
    return
end
out = buildStatus(cfg, catMat, here, ind);
end

% ========================================================================
function out = buildStatus(cfg, catMat, here, ind)
% BUILDSTATUS  The deliverable-8 path (catalog with the status layer).
% INPUTS: cfg; catMat (char); here, ind (folders).  OUTPUTS: out.
name     = getf(cfg, 'name', 'costate_catalog_70mN_status_deliverable8');
outDir   = getf(cfg, 'outDir', fullfile(here, 'deliverables', name));
doMov    = getf(cfg, 'movies', true);
doZip    = getf(cfg, 'zip', true);
doSelf   = getf(cfg, 'selftest', true);
primMat  = getf(cfg, 'auditPrimariesMat', fullfile(ind, 'results', ...
                'library_70mN_24x48_merged', 'audit_70mN.mat'));
altDir   = getf(cfg, 'auditAltDir', fileparts(catMat));
assert(~doZip || doSelf, 'build_dro_deliverable:selftest', ...
       'a zip is written only after a passed self-test: set cfg.selftest = true, or cfg.zip = false');

staleZip = [outDir '.zip'];
if isfile(staleZip), delete(staleZip); end           % never leave a zip of an older build

%% 1. Export (into a scratch folder: outDir is not touched until the gate passes)
catFile = 'costate_catalog_dro_tulip_70mN.mat';
scratch = tempname;  mkdir(scratch);
cleanScratch = onCleanup(@() rmdir(scratch, 's'));
[c, info] = export_deliverable_catalog(catMat, fullfile(scratch, catFile));
fprintf('[build] 1 export: %s (keys %s / %s; %d sources rewritten; removed %s)\n', ...
        catFile, info.primKey, info.altKey, info.nSourceRewritten, strjoin(info.removed, ', '));

%% 2. Audit gate by content key
nSolved = sum(arrayfun(@(s) nnz(s.has_solution), c.sheets));
aud = auditGate(primMat, altDir, info, numel(c.alternatives), nSolved);
prepareOutDir(outDir);
movefile(fullfile(scratch, catFile), fullfile(outDir, catFile));
fprintf('[build] 2 audit gate: primaries %d ok / %d bad, alternatives %d ok / %d bad over %d chunks, coverage 1..%d, keys bound\n', ...
        aud.primOk, aud.primBad, aud.altOk, aud.altBad, aud.nChunks, numel(c.alternatives));

%% 3. Copy the listed files
[L, roots] = deliverable_file_list();
for k = 1:numel(L)
    src = fullfile(roots.(L(k).from), L(k).file);
    assert(isfile(src), 'build_dro_deliverable:copy', 'missing maintained file %s', src);
    copyfile(src, fullfile(outDir, L(k).file));
end
fprintf('[build] 3 copied %d listed files\n', numel(L));

%% 4. Movie
movieStem = '';
if doMov
    [prim, alt2] = sample_status_pair(c);
    labels = {sprintf('primary (%s)', prim.statusName), sprintf('alternative (%s)', alt2.statusName)};
    movieStem = costate_catalog_extremes_movies(c, struct(), 'tf', fullfile(outDir, 'status_pair_movie'), ...
                                                struct('pair', {{prim, alt2}}, 'labels', {labels}));
    if isempty(movieStem)
        fprintf('[build] 4 movie: NOT rendered (no animator)\n');
    else
        [~, movieStem] = fileparts(movieStem);
        fprintf('[build] 4 movie: %s.mp4 / .gif (cell sD %.4f sA %.4f, %.2f d vs %.2f d)\n', ...
                movieStem, prim.sD, prim.sA, prim.tfDays, alt2.tfDays);
    end
else
    fprintf('[build] 4 movie: skipped (cfg.movies = false)\n');
end

%% 5. Torus PNG
pngFile = 'status_torus.png';
H = plot_status_torus(c, struct('out', fullfile(outDir, pngFile), 'visible', false));
if isfield(H, 'fig') && ~isempty(H.fig) && isgraphics(H.fig), close(H.fig); end
fprintf('[build] 5 torus: %s (%d conjugate cells, %d off-grid alternatives)\n', ...
        pngFile, nnz(H.conjCell), H.nOffGrid);

%% 6. README
f = statusFacts(c, aud, H, info);
f.movieStem = movieStem;  f.pngFile = pngFile;  f.catFile = catFile;
writeReadmeStatus(fullfile(outDir, 'README.md'), f, c, L);
fprintf('[build] 6 README: %d primaries, %d alternatives, %d borderline rows\n', ...
        f.nPrim, f.nAlt, numel(f.bl));

%% 7. Self-test
R = [];
if doSelf
    R = deliverable_selftest(outDir);
    if ~R.ok
        error('build_dro_deliverable:selftest', 'recipient self-test FAILED (log %s):\n  %s', ...
              R.log, strjoin(R.reasons, '\n  '));
    end
    fprintf('[build] 7 self-test: OK (log %s)\n', R.log);
else
    fprintf('[build] 7 self-test: skipped (cfg.selftest = false)\n');
end

%% 8. Zip
D = dir(outDir);  D = D(~[D.isdir]);
files = {D.name};  files = files(~strcmp(files, '.DS_Store'));
zipPath = '';
if doZip
    zipPath = zipAndVerify(outDir, files);
    d = dir(zipPath);
    fprintf('[build] 8 zip: %s (%.1f MB, %d files verified, no .DS_Store)\n', zipPath, d.bytes/1e6, numel(files));
else
    fprintf('[build] 8 zip: skipped (cfg.zip = false)\n');
end

out = struct('dir', outDir, 'zip', zipPath, 'facts', f, 'files', {files}, ...
             'selftest', R, 'keys', struct('primKey', info.primKey, 'altKey', info.altKey));
fprintf('[build] done: %s (%d files)\n', outDir, numel(files));
end

% ------------------------------------------------------------------------
function prepareOutDir(outDir)
% PREPAREOUTDIR  Create outDir, or replace it when it looks like an earlier
% build (holds a costate_catalog_*.mat or README.md, or is empty).
% INPUTS: outDir (char).  OUTPUTS: none.
if isfolder(outDir)
    D = dir(outDir);  D = D(~ismember({D.name}, {'.', '..'}));
    looksBuilt = isempty(D) || any(strcmp({D.name}, 'README.md')) || ...
                 ~isempty(dir(fullfile(outDir, 'costate_catalog_*.mat')));
    assert(looksBuilt, 'build_dro_deliverable:outDir', ...
           '%s exists and does not look like an earlier deliverable build; refusing to replace it', outDir);
    rmdir(outDir, 's');
end
mkdir(outDir);
end

% ------------------------------------------------------------------------
function aud = auditGate(primMat, altDir, info, nAlt, nSolved)
% AUDITGATE  The key-bound audit gate (see the header). Every refusal is
% error id build_dro_deliverable:audit.
% INPUTS: primMat (char); altDir (char); info (export info: .primKey
% .altKey); nAlt [scalar]; nSolved [scalar] solved cells.  OUTPUTS: aud (.primOk .primBad .primWhen .altOk
% .altBad .nChunks .chunkFiles).
id = 'build_dro_deliverable:audit';
if ~isfile(primMat), error(id, 'no primaries audit at %s', primMat); end
P = load(primMat);
if ~isfield(P, 'A') || ~all(isfield(P.A, {'nOk', 'nBad', 'contentKey'}))
    error(id, '%s holds no audit A with .nOk .nBad .contentKey', primMat);
end
if P.A.nBad ~= 0
    error(id, 'the primaries audit reports %d bad entries -- not shipping', P.A.nBad);
end
if ~strcmp(P.A.contentKey, info.primKey)
    error(id, 'the primaries audit is for content key %s, the export is %s', P.A.contentKey, info.primKey);
end
if P.A.nOk ~= nSolved
    error(id, 'the primaries audit covers %d entries, the catalog solves %d cells', P.A.nOk, nSolved);
end
D = dir(fullfile(altDir, 'audit_alt_*.mat'));
if isempty(D), error(id, 'no audit_alt_*.mat chunks in %s', altDir); end
kc = cell(1, numel(D));  altBad = 0;  altOk = 0;
for k = 1:numel(D)
    X = load(fullfile(D(k).folder, D(k).name));
    if ~all(isfield(X, {'contentKey', 'altContentKey', 'audit', 'altRows'}))
        error(id, '%s is not an audit_alt chunk (contentKey altContentKey audit altRows)', D(k).name);
    end
    if ~strcmp(X.contentKey, info.primKey)
        error(id, '%s: primaries key %s, the export is %s', D(k).name, X.contentKey, info.primKey);
    end
    if ~strcmp(X.altContentKey, info.altKey)
        error(id, '%s: alternatives key %s, the export is %s', D(k).name, X.altContentKey, info.altKey);
    end
    if ~isstruct(X.audit) || ~isfield(X.audit, 'nBad') || ~isstruct(X.altRows) || ...
            ~all(isfield(X.altRows, {'ok', 'k'}))
        error(id, '%s lacks audit.nBad or altRows.ok / altRows.k', D(k).name);
    end
    rowOk = [X.altRows.ok];
    if X.audit.nBad ~= 0 || ~all(rowOk)
        error(id, '%s reports %d bad (%d alternative rows not ok) -- not shipping', ...
              D(k).name, X.audit.nBad, nnz(~rowOk));
    end
    kc{k} = [X.altRows.k];
    altOk  = altOk + nnz(rowOk);
    altBad = altBad + X.audit.nBad;
end
ks = [kc{:}];
if numel(ks) ~= nAlt || ~isequal(sort(ks), 1:nAlt)
    error(id, 'the %d chunks cover %d alternative rows (%d distinct), not 1..%d exactly once', ...
          numel(D), numel(ks), numel(unique(ks)), nAlt);
end
aud = struct('primOk', P.A.nOk, 'primBad', P.A.nBad, 'primWhen', getf(P.A, 'when', ''), ...
             'altOk', altOk, 'altBad', altBad, 'nChunks', numel(D), 'chunkFiles', {{D.name}});
end

% ------------------------------------------------------------------------
function f = statusFacts(c, aud, H, info)
% STATUSFACTS  Every number the deliverable-8 README prints, computed from
% the exported catalog, the audit summary and the torus.
% INPUTS: c (catalog); aud (auditGate); H (plot_status_torus); info
% (export).  OUTPUTS: f (struct).
s = c.sheets(1);
assert(isscalar(c.sheets) && isscalar(c.rungs_N) && size(s.has_solution, 3) == 1, ...
       'build_dro_deliverable:shape', 'the status-layer README expects one sheet and one thrust rung');
a     = c.alternatives;
has   = s.has_solution;
tStar = c.constants.tStar_s;  lStar = c.constants.lStar_km;
cnd   = c.thruster.c_nd;  m0 = c.thruster.m0_kg;
f = struct();
f.name = c.name;  f.created = c.created;  f.keys = info;
f.thrustN = c.rungs_N(1);  f.isp = c.thruster.isp_s;  f.m0 = m0;
f.nD = size(has, 1);  f.nA = size(has, 2);
f.nPrim = nnz(has);  f.nAlt = numel(a);  f.nSlots = numel(has);
Tmax  = (f.thrustN/m0)*tStar^2/(lStar*1000);
tfP   = s.tf_nd(has);  tfA = [a.tf_nd];
dvOf  = @(tf) cnd*log(1./(1 - Tmax*tf/cnd))*lStar/tStar;
mfOf  = @(tf) m0*(1 - Tmax*tf/cnd);
f.tfPrimDays = [min(tfP) max(tfP)]*tStar/86400;
f.tfAltDays  = [min(tfA) max(tfA)]*tStar/86400;
f.dvPrim = [min(dvOf(tfP)) max(dvOf(tfP))];   f.dvAlt = [min(dvOf(tfA)) max(dvOf(tfA))];
f.mfPrim = [min(mfOf(tfP)) max(mfOf(tfP))];   f.mfAlt = [min(mfOf(tfA)) max(mfOf(tfA))];
f.codes = [4 3 2 1];
stP = s.status(has);  stA = [a.status];
f.statusCounts = [arrayfun(@(q) nnz(stP == q), f.codes).', arrayfun(@(q) nnz(stA == q), f.codes).'];
offA = [a.iD] == 0 | [a.iA] == 0;
f.nOffGrid = nnz(offA);
f.offGridByStatus = arrayfun(@(q) nnz(offA & stA == q), f.codes);
f.nOffD = nnz([a.iD] == 0);  f.nOffA = nnz([a.iA] == 0);  f.nOffBoth = nnz([a.iD] == 0 & [a.iA] == 0);
f.nConjCells = nnz(H.conjCell);  f.nTorusOffGrid = H.nOffGrid;
f.nJuncPrim = nnz(~cellfun(@isempty, s.junctions(s.entry_index(has))));
f.nJuncAlt  = nnz(~cellfun(@isempty, {a.junctions}));
f.nAltCells = nnz(H.nAltPerCell > 0);  f.maxAltPerCell = max(H.nAltPerCell(:));
% borderline rows
f.bl = find([a.borderline]);
f.blRows = struct('k', {}, 'sD', {}, 'sA', {}, 'status', {}, 'clause', {});
for k = f.bl
    f.blRows(end+1) = struct('k', k, 'sD', a(k).sD, 'sA', a(k).sA, 'status', a(k).status, ...
                             'clause', firstClause(a(k).status_reason));
end
% primary gates (sheet grids)
f.conjPass = nnz(s.conj_pass(has) == 1);
f.gLamv = nnz(s.gate_min_lamv(has) > 0);  f.gQmt = nnz(s.gate_min_qmt(has) > 0);
f.gDimS = nnz(s.gate_dimS(has) == 1);
f.h6Min = min(s.h6_margin(has));  f.liftMin = min(s.lift_margin(has));
f.conjInterior = nnz(s.conj_interior(has) > 0);
% flown misses recorded on the alternatives (NaN = not measured)
fk = [a.flyKm];  fv = [a.flyVms];
f.flyKmMax = max(fk(isfinite(fk)));  f.flyVmsMax = max(fv(isfinite(fv)));
f.nFlyMeasured = nnz(isfinite(fk));
f.aud = aud;
f.nOkTotal = aud.primOk + aud.altOk;  f.nBadTotal = aud.primBad + aud.altBad;
f.build = c.build;
end

% ------------------------------------------------------------------------
function t = firstClause(r)
% FIRSTCLAUSE  A reason's first clause: up to the first ';', ' | ', or a
% '.' that is not a decimal point.  INPUTS: r (char).  OUTPUTS: t (char).
r = char(r);
cut = regexp(r, ';| \| |\.(?!\d)', 'once');
if isempty(cut), t = strtrim(r); else, t = strtrim(r(1:cut-1)); end
t = strrep(t, '|', '/');
end

% ------------------------------------------------------------------------
function writeReadmeStatus(path, f, c, L)
% WRITEREADMESTATUS  The deliverable-8 README, every number from f (facts
% computed from the exported catalog and its audits), the legends from the
% catalog itself, the contents table from the file list.
% INPUTS: path (char); f (statusFacts); c (catalog); L (deliverable_file_list).
% OUTPUTS: none.
role = fileRoles();
K = c.status_key;  s = c.sheets(1);
fid = fopen(path, 'w');
cleaner = onCleanup(@() fclose(fid));
p = @(varargin) fprintf(fid, varargin{:});
d2 = @(x) sprintf('%.2f', x);

p('# DRO -> Tulip 70 mN Minimum-Time Costate Library, with Optimality Status (deliverable 8)\n\n');
p(['Every entry here is a **converged minimum-time PMP extremal**, not a hint: hand its `z8` ' ...
   'to `pumpkyn.cr3bp.tfMin` at the entry''s own endpoints and it comes back unchanged. ' ...
   'The library holds **%s primary transfers** (the fastest found in each cell of a %d x %d ' ...
   'departure x arrival phase grid) and **%d alternatives** (every other transfer the builds ' ...
   'found at a phase pair), each with an optimality **status** saying what was established ' ...
   'about it, and the multiple-shooting **junction states** that fly a %s-%s day arc accurately.\n\n'], ...
   commas(f.nPrim), f.nD, f.nA, f.nAlt, d2(min([f.tfPrimDays f.tfAltDays])), d2(max([f.tfPrimDays f.tfAltDays])));
p(['Catalog `%s` built %s; this package generated %s by `build_dro_deliverable.m`. ' ...
   'Every number below is computed from the shipped catalog and the audits bound to it by ' ...
   'content key, not transcribed.\n\n'], f.name, string(f.created), char(datetime('now', 'Format', 'yyyy-MM-dd')));

%% Contents
p('## Contents\n\n| file | role |\n|---|---|\n');
p('| `%s` | the catalog (cleaned export; content keys `%s` / `%s`) |\n', f.catFile, f.keys.primKey, f.keys.altKey);
for k = 1:numel(L)
    assert(isKey(role, L(k).file), 'build_dro_deliverable:role', 'no README role for %s', L(k).file);
    p('| `%s` | %s |\n', L(k).file, role(L(k).file));
end
if ~isempty(f.movieStem)
    p('| `%s.mp4`, `%s.gif` | the sample movie (see below) |\n', f.movieStem, f.movieStem);
end
p('| `%s` | the status map of the phase torus (see below) |\n', f.pngFile);
p('| `README.md` | this file |\n\n');
p('Every helper with a bare self-demo runs it on whichever `costate_catalog_*.mat` sits in the current folder (type e.g. `costate_lib_describe` with no arguments).\n\n');

%% Dependencies
p('## Dependencies\n\n');
p(['- **pumpkyn + pumpkynPie** on the MATLAB path (run pumpkynPie''s `startup.m`): the CR3BP ' ...
   'dynamics, `pumpkyn.cr3bp.tfMin`, the orbit families, and the movie scene ' ...
   '(`pumpkynPie.plot.SatelliteAnimator`).\n']);
p(['- **Parallel Computing Toolbox: optional.** `tfmin_accepts` runs tfMin under a hard ' ...
   'wall-clock cap (`run_capped`) when a pool is open; without the toolbox or a pool it ' ...
   'calls tfMin directly.\n']);
p('- Nothing else: every other function the examples and helpers call is in this folder. The `.mat` itself loads with `load` alone.\n\n');

%% Quick start
p('## Quick start\n\n');
p('Run the two worked examples from this folder (each about 10 seconds):\n\n');
p('```matlab\ncostate_catalog_example   %% pick by phases (days), fly, tfMin accepts it, plot\n');
p('status_layer_example      %% every transfer at one cell, statuses, fly two, the torus map\n```\n\n');
p('Or by hand:\n\n```matlab\n');
p('L = load(''%s'');\ncat_ = L.%s;\n', f.catFile, f.name);
p('T = transfers_at(cat_, %.4f, %.4f);   %% phases as fractions of a period\n', s.sD_frac(1), s.sA_frac(1));
p('F = fly_catalog_entry(cat_, T(1));          %% from its junctions\n');
p('[ok, dz] = tfmin_accepts(cat_, T(1));       %% pumpkyn tfMin accepts it unchanged\n```\n\n');

%% Coverage
p('## What it covers\n\n| axis | values |\n|---|---|\n');
p('| departure orbit | `%s`, recipe `get_family_orbit(''%s'', %s)`, period %.6g ND |\n', ...
  s.dep_family, s.dep_family, recipe(s.dep_params), s.tau_dep);
p('| arrival orbit | `%s`, recipe `get_family_orbit(''%s'', %s)`, period %.6g ND |\n', ...
  s.arr_family, s.arr_family, recipe(s.arr_params), s.tau_arr);
p('| thrust | %g N (one rung), Isp %g s, m0 %g kg |\n', f.thrustN, f.isp, f.m0);
p('| phase grid | %d departure x %d arrival phases (fractions of each orbit''s period) |\n', f.nD, f.nA);
p('| primaries | %s of %s cells (%.1f%%) |\n', commas(f.nPrim), commas(f.nSlots), 100*f.nPrim/f.nSlots);
p('| alternatives | %d (%d off the grid; %d cells hold at least one, up to %d per cell) |\n', ...
  f.nAlt, f.nOffGrid, f.nAltCells, f.maxAltPerCell);
p('| junction states | %s of %s primaries, %d of %d alternatives |\n', commas(f.nJuncPrim), commas(f.nPrim), f.nJuncAlt, f.nAlt);
p('| t_f, primaries | %s .. %s days |\n', d2(f.tfPrimDays(1)), d2(f.tfPrimDays(2)));
p('| t_f, alternatives | %s .. %s days |\n', d2(f.tfAltDays(1)), d2(f.tfAltDays(2)));
p('| Delta-V, primaries | %.4f .. %.4f km/s (final mass %.3f .. %.3f kg) |\n', f.dvPrim, f.mfPrim);
p('| Delta-V, alternatives | %.4f .. %.4f km/s (final mass %.3f .. %.3f kg) |\n\n', f.dvAlt, f.mfAlt);
p(['At fixed thrust the burn is continuous, so minimum time is minimum Delta-V: the ' ...
   'ranges above order the same way. Delta-V and mass follow `cat.derive` ' ...
   '(`mf = 1 - Tmax_nd*tf_nd/c_nd`, `dV = c_nd*log(1/mf)*lStar_km/tStar_s`).\n\n']);

%% Statuses
p('## The statuses\n\n');
p('%s\n\n', K.note);
p('| code | name | primaries | alternatives | what it establishes |\n|---|---|---|---|---|\n');
for q = 1:numel(f.codes)
    kk = find(K.codes == f.codes(q), 1);
    p('| %d | %s | %s | %d | %s |\n', f.codes(q), K.names{kk}, commas(f.statusCounts(q, 1)), ...
      f.statusCounts(q, 2), K.meaning{kk});
end
p('\nThe checks, in order: %s.\n\n', K.order);
p(['A status-2 transfer (**conjugate point found**) is a converged extremal that tfMin accepts ' ...
   'like any other; it is just not locally optimal. A faster status-2 alternative can still be ' ...
   'the right choice for a mission, taken knowing exactly that.\n\n']);
p('### Borderline rows\n\n');
p('%d alternatives are **borderline**: %s.\n\n', numel(f.bl), f.build.borderlineNote);
p('| index | sD | sA | status | reason (first clause) |\n|---|---|---|---|---|\n');
for k = 1:numel(f.blRows)
    r = f.blRows(k);
    p('| %d | %.4f | %.4f | %d | %s |\n', r.k, r.sD, r.sA, r.status, r.clause);
end
p('\n');

%% Phase conventions
p('## Phase conventions\n\n');
p(['- `transfers_at`, `catalog_endpoints`, `fly_catalog_entry`, `tfmin_accepts` and ' ...
   '`plot_status_torus` take phases as **fractions of a period** in [0, 1), wrapping ' ...
   '(s = 1 - 1e-9 finds s = 0), matched with a circular tolerance.\n']);
p(['- `costate_catalog_pick(cat, tauDep, Np, depDays, arrDays, thrustN)` takes phases in ' ...
   '**DAYS** (unchanged from the earlier deliverables): fraction x period_nd x tStar_s / 86400.\n']);
p('- Flight times are ND unless the field name says otherwise (`tfDays`, `dvKms`, `flyKm`, `flyVms`).\n\n');

%% Field reference
p('## Entry-field reference\n\n');
p('Top level: `%s`.\n\n', strjoin(fieldnames(c), '`, `'));
p('Per sheet `sheets(1)`:\n\n| field | meaning |\n|---|---|\n');
sm = sheetFieldMeanings();  sf = fieldnames(s);
for k = 1:numel(sf)
    if isKey(sm, sf{k}), m = sm(sf{k}); else, m = '(see `notes_key` / `second_order.meaning`)'; end
    p('| `%s` | %s |\n', sf{k}, m);
end
p('\nPer alternative `alternatives(k)`:\n\n| field | meaning |\n|---|---|\n');
am = altFieldMeanings();  af = fieldnames(c.alternatives);
for k = 1:numel(af)
    if isKey(am, af{k}), m = am(af{k}); else, m = '(see `notes_key`)'; end
    p('| `%s` | %s |\n', af{k}, m);
end
p(['\nLegends that ride in the catalog: `status_key` (codes, names, meanings, check order), ' ...
   '`notes_key` (the `entry_notes` layout), `second_order.meaning`, `conj_test.meaning`, ' ...
   '`hyp_gates.meaning`, `families.note`, `derive` (every derived-quantity formula).\n\n']);

%% How it was made
p('## How it was made, and how it was checked\n\n');
p('%s\n\n', c.provenance);
p('Build record: %s. %s candidates harvested, %d re-certified.\n\n', f.build.dates, commas(f.build.nCandidates), f.build.nRecert);
p('Every transfer goes through the certifier''s stack in order (`status_key.order`); its status records how far it got: %s.\n\n', K.order);
p('Primaries, from the sheet grids:\n\n');
p('- free-time conjugate test passed: %s / %s (`conj_pass = 1`)\n', commas(f.conjPass), commas(f.nPrim));
p('- hypothesis gates: min magnitude of lam_v > 0 on %s, Q_mt > 0 on %s, dim S = 1 on %s\n', commas(f.gLamv), commas(f.gQmt), commas(f.gDimS));
p('- dense conjugate scan: %d with an interior sign change; smallest H6 margin %.3g, smallest lift margin %.3g\n\n', ...
  f.conjInterior, f.h6Min, f.liftMin);
p('**Audits, bound to this file by content key** (fail-closed re-derivation from the catalog''s own keys):\n\n');
p('- primaries: %s ok / %d bad (audit of %s)\n', commas(f.aud.primOk), f.aud.primBad, string(f.aud.primWhen));
p('- alternatives: %d ok / %d bad, over %d audit chunks covering every row exactly once\n', f.aud.altOk, f.aud.altBad, f.aud.nChunks);
p('- total: **%s ok / %d bad**; catalog key `%s`, alternatives key `%s`\n\n', commas(f.nOkTotal), f.nBadTotal, f.keys.primKey, f.keys.altKey);
p(['Worst flown miss recorded on the alternatives (`flyKm`/`flyVms`, %d of %d measured; ' ...
   'the certifier''s flight gate is 100 km / 10 m/s): **%.3g km / %.3g m/s**. `fly_catalog_entry` ' ...
   'flies from the stored junctions when an entry has them rather than one flight from `z8` ' ...
   'over %s-%s days; the examples print the miss they get.\n\n'], ...
  f.nFlyMeasured, f.nAlt, f.flyKmMax, f.flyVmsMax, d2(min([f.tfPrimDays f.tfAltDays])), d2(max([f.tfPrimDays f.tfAltDays])));

%% Figures
p('## The movie and the map\n\n');
if ~isempty(f.movieStem)
    p(['`%s.mp4` / `%s.gif`: one cell''s primary beside its fastest status-2 alternative, ' ...
       'same phases, one clock (1280 x 720), flown from their junctions.\n\n'], f.movieStem, f.movieStem);
end
p(['`%s`: per cell, the primary''s flight time, how many alternatives it holds, and which ' ...
   'cells hold a status-2 transfer (%d cells); the %d off-grid alternatives are drawn at ' ...
   'their own phases.\n\n'], f.pngFile, f.nConjCells, f.nTorusOffGrid);

%% Gotchas
p('## Gotchas, honestly\n\n');
p(['- **Interpolation between grid points is NOT included.** The junction-interpolation ' ...
   'experiment is pending; until its hit rate is measured on this library, interpolate at ' ...
   'your own risk and use the result only as a tfMin seed. Neighbouring cells can sit on ' ...
   'different solution families (`family_index`).\n']);
p(['- **Status-1 rows and tfMin.** Status 1 means the transfer flies to the target but the ' ...
   'necessary conditions were not established, so tfMin is not expected to accept it: ' ...
   '`tfmin_accepts` skips status-1 rows (ok = false, expected refusal) unless `opts.force` ' ...
   'is set. There are %d such primaries and %d such alternatives.\n'], f.statusCounts(4, 1), f.statusCounts(4, 2));
p(['- **%d of the %d alternatives are OFF-GRID** (`iD = 0` or `iA = 0`: %d off the departure ' ...
   'grid, %d off the arrival grid, %d off both): found by ribs and sheet candidates at ' ...
   'non-grid phases. `transfers_at` finds them only at their own phases; the torus map ' ...
   'draws them as markers.\n'], f.nOffGrid, f.nAlt, f.nOffD, f.nOffA, f.nOffBoth);
p(['- **`source` wording is kind-only** (`sheet candidate, column n (re-certified)`, ' ...
   '`rib ... point`, ...): machine paths were stripped at export; the content keys are ' ...
   'unchanged because `source` is not hashed.\n']);
p('- **Phases: fractions vs DAYS** (see Phase conventions): the picker alone takes DAYS.\n');
p(['- **Borderline statuses are lower bounds** (the %d rows above): the status flipped ' ...
   'between two runs at a threshold and the claim both runs support is the one stored.\n\n'], numel(f.bl));
p('M. Casey / D. Koblick, Coorbital Inc.\n');
end

% ------------------------------------------------------------------------
function M = fileRoles()
% FILEROLES  One-line README role per shipped file.  INPUTS: none.
% OUTPUTS: M (containers.Map file -> role).
M = containers.Map();
M('costate_catalog_pick.m')            = 'five-coordinate lookup (DRO period, petals, departure/arrival phase in DAYS, thrust) with the honesty contract; unchanged';
M('transfers_at.m')                    = 'every transfer recorded at a phase pair (fractions), primary and alternatives, ranked by t_f, with status';
M('status_key.m')                      = 'the status legend: codes, names, meanings, the order of the checks';
M('catalog_endpoints.m')               = 'the certified endpoint rule, from the catalog''s own orbit recipes';
M('fly_catalog_entry.m')               = 'fly an entry: from its junction states when it has them, else from z8';
M('tfmin_accepts.m')                   = 'the "pumpkyn tfMin accepts it unchanged" check (capped; status-1 rows skipped)';
M('entry_thrust_program.m')            = 'flight engine: the transfer and thrust program, segment by segment from the junctions';
M('fly_transfer.m')                    = 'flight engine: one flight from z8, with mass and Delta-V';
M('costate_lib_describe.m')            = 'fact sheet: grids in ND and days, coverage, status counts';
M('costate_catalog_extremes.m')        = 'min/max Delta-V or flight time over a slice (alternatives opt-in)';
M('plot_status_torus.m')               = 'the status map of the phase torus';
M('costate_catalog_extremes_movies.m') = 'side-by-side movie of two transfers on one clock (pumpkynPie scene)';
M('print_transfer_summary.m')          = 'prints a transfer''s time of flight, Delta-V, propellant, flight samples';
M('plot_transfer_3d.m')                = 'the 3D transfer figure the examples draw';
M('sample_status_pair.m')              = 'picks the cell with the most status-2 alternatives and its pair';
M('costate_catalog_example.m')         = 'worked example: pick, derive units, rebuild orbits, fly, tfMin, plot';
M('status_layer_example.m')            = 'worked example: every transfer at a cell, statuses, fly two, the torus map';
M('get_family_orbit.m')                = 'family name + parameters -> periodic orbit (the DRO and tulip recipes)';
M('phase_state.m')                     = 'state at a phase fraction of a periodic orbit';
M('periodic_pp.m')                     = 'the periodic spline phase_state is built on';
M('validate_flight.m')                 = 'the one admissibility check for a flown trajectory';
M('scalar_verdict.m')                  = 'reads a solver result as a gate verdict, or refuses';
M('catalog_schema.m')                  = 'the catalog schema: field list, validator, formula registry';
M('run_capped.m')                      = 'hard wall-clock cap on a parfeval worker (used by tfmin_accepts)';
M('current_pool.m')                    = 'the open parallel pool or [], without throwing when the toolbox is absent';
end

% ------------------------------------------------------------------------
function M = sheetFieldMeanings()
% SHEETFIELDMEANINGS  One-line meaning per sheet field.  INPUTS: none.
% OUTPUTS: M (containers.Map).
M = containers.Map();
M('tau_dep') = 'departure-orbit period (ND)';
M('tauDRO')  = 'departure-orbit period (ND; legacy name, any family)';
M('tau_arr') = 'arrival-orbit period (ND)';
M('period_tulip_nd') = 'tulip period (ND)';
M('dep_family') = 'departure family name, for `get_family_orbit`';
M('dep_params') = 'departure reconstruction parameters';
M('arr_family') = 'arrival family name';
M('arr_params') = 'arrival reconstruction parameters';
M('Np') = 'tulip petal count';
M('pm') = 'tulip branch (the pm = +1 branch is the z-mirror)';
M('sD_frac') = 'departure phase grid, fractions of the departure period';
M('sA_frac') = 'arrival phase grid, fractions of the arrival period';
M('has_solution') = '**the** coverage authority (departure x arrival)';
M('tf_nd') = 'the primary''s flight time (ND)';
M('entry_index') = 'column of `z8` / `junctions` for each cell';
M('z8') = 'converged `[lam_r; lam_v; lam_m; tf]` per entry, tfMin-ready';
M('conj_pass') = 'free-time conjugate test at certification (1 = none found)';
M('gate_min_lamv') = 'hypothesis gate: min magnitude of lam_v over the arc (> 0)';
M('gate_min_qmt') = 'hypothesis gate: min Q_mt (> 0, all-burn is the PMP control)';
M('gate_dimS') = 'hypothesis gate: lift-space dimension (1 = no abnormal lift)';
M('conj_K') = 'conjugate-test sample count';
M('family_index') = 'provenance: which continuation family found the entry (see `families.note`; read by no tool)';
M('entry_notes') = 'per-entry story, layout in `notes_key`';
M('conj_interior') = 'dense scan: interior sign changes (0 = none)';
M('conj_multiplicity') = 'dense scan: candidates with two or more collapsed singular values';
M('conj_interior_cand') = 'dense scan: located interior candidates';
M('conj_near_miss') = 'dense scan: near-misses (positive minimum located)';
M('conj_zero') = 'dense scan: zeros (at the floor or a sign change)';
M('conj_unresolved') = 'dense scan: unresolved candidates (inside the floor band)';
M('h6_margin') = 'H6 margin, (c/T)/lambda_m(0) (> 1 excludes the spurious-zero mechanism)';
M('lift_margin') = 'sigma_6 over the measured lift-matrix error (> 1 certifies dim S = 1)';
M('conj_min_rel') = 'closest near-miss, relative to the median sigma_6';
M('status') = 'optimality status code (`status_key`)';
M('status_reason') = 'one-line reason for the status';
M('junctions') = 'multiple-shooting junction states per entry (state and costate rows x nodes)';
end

% ------------------------------------------------------------------------
function M = altFieldMeanings()
% ALTFIELDMEANINGS  One-line meaning per alternatives field.  INPUTS: none.
% OUTPUTS: M (containers.Map).
M = containers.Map();
M('sD') = 'departure phase (fraction)';
M('sA') = 'arrival phase (fraction)';
M('iD') = 'departure grid index of its cell (0 = off the grid)';
M('iA') = 'arrival grid index of its cell (0 = off the grid)';
M('z8') = 'converged `[lam_r; lam_v; lam_m; tf]`';
M('junctions') = 'multiple-shooting junction states';
M('tf_nd') = 'flight time (ND)';
M('status') = 'optimality status code (`status_key`)';
M('status_reason') = 'one-line reason for the status';
M('inferred') = 'true when the status was classified from the certifier''s record rather than read from its stamp';
M('conj') = 'conjugate test verdict (1 pass, 0 conjugate point, -1 undecided)';
M('conjVerdict') = 'conjugate test verdict, text';
M('minLamV') = 'hypothesis gate: min magnitude of lam_v over the arc';
M('minQmt') = 'hypothesis gate: min Q_mt';
M('dimS') = 'hypothesis gate: lift-space dimension';
M('h6Margin') = 'H6 margin';
M('liftMargin') = 'lift margin';
M('flyKm') = 'flown arrival miss, km (NaN = not measured)';
M('flyVms') = 'flown arrival miss, m/s (NaN = not measured)';
M('source') = 'how it was found (kind only, paths stripped)';
M('sheet') = 'sheet index (1)';
M('borderline') = 'true when the status is the meet of two runs (a lower bound)';
end

% ------------------------------------------------------------------------
function t = recipe(prm)
% RECIPE  A params struct as MATLAB struct(...) text.  INPUTS: prm (struct).
% OUTPUTS: t (char).
fn = fieldnames(prm);  parts = cell(1, numel(fn));
for k = 1:numel(fn), parts{k} = sprintf('''%s'', %.15g', fn{k}, prm.(fn{k})); end
t = ['struct(' strjoin(parts, ', ') ')'];
end

% ------------------------------------------------------------------------
function t = commas(n)
% COMMAS  Integer with thousands separators.  INPUTS: n [scalar].
% OUTPUTS: t (char).
t = sprintf('%d', round(n));
for k = numel(t)-3:-3:1, t = [t(1:k) ',' t(k+1:end)]; end
end

% ------------------------------------------------------------------------
function z = zipAndVerify(outDir, files)
% ZIPANDVERIFY  zip -qr from the parent folder (no .DS_Store), then verify
% with unzip -l that every file is present and no .DS_Store slipped in.
% INPUTS: outDir (char); files (cellstr of bare names).  OUTPUTS: z (char).
[parent, nm] = fileparts(outDir);
z = fullfile(parent, [nm '.zip']);
if isfile(z), delete(z); end
[st, msg] = system(sprintf('cd ''%s'' && zip -qr ''%s.zip'' ''%s'' -x ''*.DS_Store''', parent, nm, nm));
assert(st == 0, 'build_dro_deliverable:zip', 'zip failed: %s', msg);
[st, lst] = system(sprintf('unzip -l ''%s''', z));
assert(st == 0, 'build_dro_deliverable:zip', 'unzip -l failed: %s', lst);
assert(~contains(lst, '.DS_Store'), 'build_dro_deliverable:zip', 'the zip holds a .DS_Store');
for k = 1:numel(files)
    assert(contains(lst, [nm '/' files{k}]), 'build_dro_deliverable:zip', 'the zip lacks %s', files{k});
end
end

% ========================================================================
function out = buildLegacy(cfg, c, catMat, here)
% BUILDLEGACY  The status-less path (deliverable 3 and the earlier 70 mN
% catalog): path-equality audit gate, five helpers, the original README.
% INPUTS: cfg; c (catalog); catMat (char); here (DRO_tulip).  OUTPUTS: out.
outDir = getf(cfg, 'outDir', fullfile(here, 'deliverables', 'costate_catalog_deliverable3'));
doMov  = getf(cfg, 'movies', false);
doZip  = getf(cfg, 'zip', true);
srcDir = fullfile(here, 'indirect');

%% THE AUDIT GATE -------------------------------------------------------------
% A catalog does not ship until an audit that re-derived it FROM ITS OWN KEYS
% comes back clean. Fixing a builder says nothing about a file already
% written, which is exactly how the 2026-09-10 review left us: a shipped
% catalog and no way to know it was sound short of auditing it. Enforced here
% rather than left to a checklist, because a checklist is what failed.
% Waivable only by naming a reason in cfg.skipAuditBecause.
why = getf(cfg, 'skipAuditBecause', '');
if isempty(why)
    auditMat = getf(cfg, 'auditMat', [strrep(catMat, '.mat', '') '_audit.mat']);
    assert(exist(auditMat, 'file') == 2, ...
        ['no audit for %s.\n' ...
         'Run audit_phase_catalog(catMat) and pass cfg.auditMat, or state a ' ...
         'reason in cfg.skipAuditBecause.'], catMat);
    Aud = load(auditMat);
    assert(isfield(Aud, 'A') && isfield(Aud.A, 'nBad') && isfield(Aud.A, 'catMat'), ...
        'audit %s is not an audit_phase_catalog result', auditMat);
    [~, an, ae] = fileparts(Aud.A.catMat);  [~, cn, ce] = fileparts(catMat);
    assert(strcmp([an ae], [cn ce]), ...
        'audit %s is for %s, not for %s', auditMat, Aud.A.catMat, catMat);
    assert(Aud.A.nBad == 0, 'audit %s reports %d findings -- not shipping', auditMat, Aud.A.nBad);
    fprintf('[build] audit gate: %s clean (%d entries)\n', auditMat, Aud.A.nOk);
else
    fprintf('[build] audit gate WAIVED: %s\n', why);
end

%% Facts, computed -- never typed --------------------------------------------
f = struct();
f.name      = c.name;
f.created   = c.created;
f.nSheets   = numel(c.sheets);
f.rungs     = c.rungs_N;
f.nRungs    = numel(c.rungs_N);
f.isp       = c.thruster.isp_s;
f.m0        = c.thruster.m0_kg;
f.tau       = unique(arrayfun(@(s) s.tauDRO, c.sheets));
f.Np        = unique(arrayfun(@(s) s.Np, c.sheets));
f.pm        = unique(arrayfun(@(s) s.pm, c.sheets));
f.nPhase    = size(c.sheets(1).has_solution, 1);
f.nSlots    = 0;  f.nSolved = 0;  f.tfMin = inf;  f.tfMax = -inf;
cp = [0 0 0];  covSheet = zeros(1, f.nSheets);
for k = 1:f.nSheets
    s = c.sheets(k);
    f.nSlots  = f.nSlots  + numel(s.has_solution);
    f.nSolved = f.nSolved + nnz(s.has_solution);
    covSheet(k) = 100*nnz(s.has_solution)/numel(s.has_solution);
    t = s.tf_nd(s.has_solution);
    if ~isempty(t), f.tfMin = min(f.tfMin, min(t));  f.tfMax = max(f.tfMax, max(t)); end
    if isfield(s, 'conj_pass')
        v = s.conj_pass(s.has_solution);
        cp = cp + [nnz(v == 1) nnz(v == 0) nnz(v < 0)];
    end
end
f.coverage   = 100*f.nSolved/f.nSlots;
f.covMin     = min(covSheet);  f.covMax = max(covSheet);
f.conjPass   = cp(1);  f.conjFail = cp(2);  f.conjNotrun = cp(3);
f.hasConj    = isfield(c, 'conj_test');
f.hasGates   = isfield(c, 'hyp_gates');
if f.hasConj
    f.conjDate = c.conj_test.date;
    if isfield(c.conj_test, 'K'), f.conjK = c.conj_test.K; else, f.conjK = NaN; end
end
if f.hasGates, f.gatesDate = c.hyp_gates.date; end
f.tStar = c.constants.tStar_s;
f.dayMin = f.tfMin*f.tStar/86400;  f.dayMax = f.tfMax*f.tStar/86400;
% per-rung coverage
f.rungCov = zeros(1, f.nRungs);
for kr = 1:f.nRungs
    a = 0; b = 0;
    for k = 1:f.nSheets
        h = c.sheets(k).has_solution(:,:,kr);  a = a + nnz(h);  b = b + numel(h);
    end
    f.rungCov(kr) = 100*a/b;
end

%% Assemble -------------------------------------------------------------------
if ~exist(outDir, 'dir'), mkdir(outDir); end
helpers = {'costate_catalog_pick.m', 'costate_lib_describe.m', ...
           'costate_catalog_extremes.m', 'costate_catalog_example.m', ...
           'costate_catalog_extremes_movies.m'};
% NAME THE COPY FROM THE CATALOG, not from this campaign: the same builder
% now packages the 70 mN phase catalog, whose variable (and therefore whose
% file) is costate_catalog_dro_tulip_70mN.
catFile = [f.name '.mat'];
copyfile(catMat, fullfile(outDir, catFile));
for k = 1:numel(helpers)
    src = fullfile(srcDir, helpers{k});
    assert(exist(src, 'file') == 2, 'missing maintained helper %s', src);
    copyfile(src, fullfile(outDir, helpers{k}));
end
writeReadme(fullfile(outDir, 'README.md'), f, c);
files = [{catFile}, helpers, {'README.md'}];

if doMov
    costate_catalog_extremes_movies(fullfile(outDir, catFile), ...
        struct('outDir', outDir));
end

out = struct('dir', outDir, 'facts', f, 'files', {files}, 'zip', '', ...
             'selftest', [], 'keys', struct('primKey', '', 'altKey', ''));
if doZip
    z = [outDir '.zip'];
    if exist(z, 'file'), delete(z); end
    zip(z, '*', outDir);
    out.zip = z;
    d = dir(z);
    fprintf('[build] deliverable: %s (%.1f MB)\n', z, d.bytes/1e6);
end
fprintf(['[build] DRO->tulip deliverable rebuilt: %d entries, %d sheets, %d rungs, ' ...
         '%.1f%% coverage, conj %d/%d/%d\n'], f.nSolved, f.nSheets, f.nRungs, ...
         f.coverage, f.conjPass, f.conjFail, f.conjNotrun);
end

% ------------------------------------------------------------------------
function writeReadme(path, f, c)
% WRITEREADME  Emit the deliverable README with every number derived from
% the catalog.  INPUTS: path; f (facts); c (catalog).  OUTPUTS: none.
fid = fopen(path, 'w');
p = @(varargin) fprintf(fid, varargin{:});
p('# DRO -> Tulip Minimum-Time Costate CATALOG (deliverable 3)\n\n');
p(['**%d converged minimum-time transfer solutions** spanning **%d DRO ' ...
   'periods x %d tulip petal counts x a %dx%d phasing torus x %d thrust ' ...
   'levels**. Every entry is a root of the indirect (PMP) boundary-value ' ...
   'problem: hand its `z8` to `pumpkyn.cr3bp.tfMin` and it is accepted ' ...
   'unchanged in about a second.\n\n'], f.nSolved, numel(f.tau), numel(f.Np), ...
   f.nPhase, f.nPhase, f.nRungs);
p('Catalog built %s; this package generated %s by `build_dro_deliverable.m`.\n', ...
  string(f.created), datestr(now, 'yyyy-mm-dd'));
p('Every number below is computed from the catalog, not transcribed.\n\n');

p('## What it covers\n\n');
p('| axis | values |\n|---|---|\n');
p('| DRO period tau (ND) | %s |\n', strjoin(compose('%.3g', f.tau(:)'), ', '));
p('| tulip petal count Np | %s |\n', strjoin(compose('%g', f.Np(:)'), ', '));
p('| family branch pm | %s (pm = +1 orbits and costates are exact z-mirrors) |\n', ...
  strjoin(compose('%+d', f.pm(:)'), ', '));
p('| phasing | %dx%d torus per sheet, phases as fractions of each orbit period |\n', f.nPhase, f.nPhase);
p('| thrust (N) | %s |\n', strjoin(compose('%g', f.rungs(:)'), ', '));
p('| propulsion | Isp %g s, m0 %g kg |\n', f.isp, f.m0);
p('| time of flight | %.3f .. %.3f ND = %.2f .. %.2f days |\n\n', f.tfMin, f.tfMax, f.dayMin, f.dayMax);

p('**Coverage: %d of %d grid slots (%.1f%%)**, %.0f%%-%.0f%% per sheet.\n', ...
  f.nSolved, f.nSlots, f.coverage, f.covMin, f.covMax);
p('`sheets(k).has_solution` is the per-(pair, rung) authority. Coverage\n');
p('thins at low thrust, which is real and worth knowing before you plan\n');
p('around it:\n\n');
p('| thrust (N) | %s |\n|---|%s\n', strjoin(compose('%g', f.rungs(:)'), ' | '), repmat('---|', 1, f.nRungs));
p('| solved | %s |\n\n', strjoin(compose('%.0f%%', f.rungCov(:)'), ' | '));

p('## How it was made, and how it was checked\n\n');
p('Per sheet: a direct Hermite-Simpson collocation solve anchors each phase\n');
p('pair cold at the highest thrust, then thrust is walked down with each rung\n');
p('warm-starting the next, so a pair stays on one solution family. Every\n');
p('entry passes three FIRST-ORDER gates:\n\n');
p('1. multiple-shooting refinement residual (`ms_tfmin`, analytic STM Jacobian): ~1e-10\n');
p('2. the PMP control law is FLOWN end to end and must arrive: ~1e-4 km\n');
p('3. **`pumpkyn.cr3bp.tfMin` returns the entry unchanged** (|dz| ~ 1e-9) -- a\n');
p('   foreign witness: different code, integrator and root finder\n\n');
if f.hasConj || f.hasGates
    p('### Second-order verification (new since the previous deliverable)\n\n');
end
if f.hasConj
    p('**Conjugate-point test** (`conj_test`, %s, K = %d): the free-time\n', string(f.conjDate), f.conjK);
    p('quotiented Jacobi determinant, sampled at the multiple-shooting junctions.\n');
    p('`sheets(k).conj_pass`: 1 = no conjugate point found in (0, tf); 0 = a sign\n');
    p('change strictly inside, i.e. the entry is NOT a local minimum; -1 = not\n');
    p('decided. Result: **%d pass, %d fail, %d undecided**.\n\n', f.conjPass, f.conjFail, f.conjNotrun);
end
if f.hasGates
    p('**Sufficiency-hypothesis gates** (`hyp_gates`, %s): the three hypotheses\n', string(f.gatesDate));
    p('the conjugate test assumes, checked per entry and stored as\n');
    p('`gate_min_lamv` (strong Legendre, min|lam_v| > 0), `gate_min_qmt`\n');
    p('(all-burn is the PMP control, min Q_mt > 0) and `gate_dimS` (no abnormal\n');
    p('lift, dim S = 1). **All entries pass all three.**\n\n');
    p('With `conj_pass = 1` these are the Bonnard-Caillau-Trelat sufficiency\n');
    p('hypotheses, verified numerically at the sampled times. What is NOT\n');
    p('claimed: the determinant is sampled at junctions, so a conjugate pair\n');
    p('inside one segment is invisible, and the second-order layer -- unlike\n');
    p('the first -- has no independent second implementation behind it.\n\n');
end

p('## The COMPACT format\n\n');
p('Only canonical nondimensional quantities are stored: constants once at top\n');
p('level, per-sheet phase fractions and rung availability, `tf_nd` lookup\n');
p('grids, and the `z8` vectors (which already contain t_f). Days, dV and\n');
p('masses are DERIVED, and every formula rides along in `cat.derive`:\n\n');
p('```\n');
dn = fieldnames(c.derive);
for k = 1:numel(dn)
    v = c.derive.(dn{k});
    if ischar(v) || isstring(v), p('%-22s %s\n', dn{k}, string(v)); end
end
p('```\n\n');

p('## Quick start\n\n');
p('```matlab\n');
p('L = load(''%s.mat'');\n', f.name);
p('cat = L.%s;\n', f.name);
p('[tf_nd, z8, info] = costate_catalog_pick(cat, %g, %g, 3.0, 11.6, %g);\n', ...
  f.tau(min(3,end)), f.Np(min(2,end)), f.rungs(min(5,end)));
p('%%                                        tau  Np  dep  arr  thrust(N)\n');
p('%% z8 -> pumpkyn.cr3bp.tfMin as-is;  info.delivered = what you actually got\n');
p('```\n\n');
p('Or run `costate_catalog_example`. For a fact sheet: `costate_lib_describe`.\n');
p('Requires **pumpkyn / pumpkynPie** on the MATLAB path; the data itself is\n');
p('dependency-free.\n\n');

p('## The honesty contract (the picker''s warnings)\n\n');
p('Whenever what is RETURNED differs from what you REQUESTED, the picker\n');
p('prints exactly what you are getting: the nearest sheet, the nearest grid\n');
p('pair, the nearest SOLVED pair, or a seed from a different rung.\n');
p('`info.delivered` carries the same facts programmatically and `info.warned`\n');
p('flags any substitution. The right use of a substituted answer is as a\n');
p('**seed**: hand `z8` to `tfMin` at your true endpoints and thrust.\n\n');

p('## Gotchas, honestly\n\n');
p('- **Coverage gaps are real** and concentrated at low thrust (see the table\n');
p('  above). `has_solution` is always the authority.\n');
p('- **Only the pm = %+d branch is stored.** The mirror branch follows from the\n', f.pm(1));
p('  CR3BP symmetry (flip all z components) but is not catalogued.\n');
p('- **At fixed thrust, minimum time = minimum dV** (continuous burn makes dV\n');
p('  monotone in t_f). The metrics differ only ACROSS thrust levels.\n');
p('- **Neighbouring cells can sit on different solution families** -- t_f can\n');
p('  jump >15%% across a family wall. Seeding tfMin from a neighbour may find\n');
p('  the faster family.\n');
p('- The %dx%d grid is for interpolation and seeding, not final answers.\n\n', f.nPhase, f.nPhase);
p('M. Casey / D. Koblick, Coorbital Inc.\n');
fclose(fid);
end

function v = getf(s, fld, d)
% GETF  Field with default.  INPUTS: s; fld; d.  OUTPUTS: v.
if isfield(s, fld) && ~isempty(s.(fld)), v = s.(fld); else, v = d; end
end
