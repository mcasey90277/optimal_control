function out = run_phase_torus(spec)
%% Purpose:
%
%   Build a minimum-time costate library over a USER-CHOSEN set of
%   departure and arrival phases -- the phase torus -- by the round loop
%   that built the 70 mN DRO -> tulip 24 x 24 library by hand (FINDINGS
%   59-72, process/PHASE_TORUS_RUNBOOK.md, doc/phase_torus_methods.tex):
%
%     round k:  arcs for every anchor without them (batch jobs, both
%               arrival directions, partial saves)
%            -> the arrival sheet at the listed phases, from every arc and
%               every seed; columns whose spine changed since round k-1
%            -> ribs for those columns (the supervised queue), every
%               earlier round's ribs offered to the packager; package,
%               audit, sweep (the finalizer)
%            -> holes and the improve pass (fill_holes_direct), then a
%               re-package if anything was added
%            -> DISCOVERY: at every empty column and every column whose
%               spine is far slower than a neighbour's, a branch-blind
%               direct solve seeded from a DIFFERENT family's root nearby;
%               a certified root faster than the spine is seeded into the
%               next sheet and becomes a new anchor
%     until a round adds no anchor and changes no spine (or .maxRounds).
%
%   Every stage leaves its product on disk before the next starts, every
%   wait has a deadline and reads a verdict FILE, and the driver's own
%   state (which round, which anchors) is in <outDir>/torus_state.mat, so a
%   killed driver resumes where it was. .plan = true prints the next
%   round's plan and returns without launching anything.
%
%% Inputs:
%
%  spec                     struct
%   .sD, .sA                [1 x nD], [1 x nA]     phases in [0,1), any spacing,
%                                                   strictly increasing (give
%                                                   a lattice that wraps past 1
%                                                   as sort(mod(., 1)));
%                                                   sD(1) is the spine's
%                                                   departure phase
%   .orbits                 struct                  .tauDRO .NpTulip .pmTulip
%   .engine                 struct                  .thrustN .ispS .m0kg
%   .anchors                {name, file, sA, label; ...}  certified roots at
%                                                   sD(1) (best.z, best.it.Y)
%   .outDir                 char                    campaign root
%   .tag                    char                    names the arcs
%                                                   (arrival_arc_<tag>_*)
%                                                   ['torus']
%   .arc                    struct                  .nStep [4000] .deadlineSec
%                                                   [6*3600] .span [1.15]
%   .nWorkers               int                     rib workers [4]
%   .improveDays, .slowDays double                  [2], [2]
%   .maxRounds              int                     [6]
%   .discover               logical                 [true]
%   .librarySeeds           logical                 seed from the 70 mN
%                                                   library too [false]
%   .matlab                 char                    batch MATLAB binary
%                                                   [R2026a]
%   .plan                   logical                 print, launch nothing [false]
%
%% Outputs:
%
%  out                      struct                  .state ('done' | 'planned'
%                                                   | 'failed') .rounds
%                                                   (struct array) .final
%                                                   (folder) .anchors
%
%% Revision History:
%  M. Casey                                                   (c) 09/15/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
d = @(f, v) fieldd(spec, f, v);
sD = spec.sD(:).';  sA = spec.sA(:).';
assert(all(diff(sD) > 0) && all(diff(sA) > 0) && all(sD >= 0 & sD < 1) && all(sA >= 0 & sA < 1), ...
       'run_phase_torus: .sD and .sA must be strictly increasing lists in [0,1)');
orbits = spec.orbits;  engine = spec.engine;
outDir = absPath(spec.outDir);  if ~isfolder(outDir), mkdir(outDir); end
tag = d('tag', 'torus');
arc = fieldd(spec, 'arc', struct());  arc = withDefaults(arc, struct('nStep', 4000, 'deadlineSec', 6*3600, 'span', 1.15));
nWorkers = d('nWorkers', 4);  improveDays = d('improveDays', 2);  slowDays = d('slowDays', 2);
maxRounds = d('maxRounds', 6);  discover = d('discover', true);  plan = d('plan', false);
librarySeeds = d('librarySeeds', false);
matlabBin = d('matlab', '/Applications/MATLAB_R2026a.app/bin/matlab');
logF = fullfile(outDir, 'torus.log');
lg = @(varargin) logmsg(logF, sprintf(varargin{:}));
resDir = fullfile(here, 'results');                 % where the sheet builder globs arcs
seedFile = fullfile(outDir, 'direct_certified.mat'); % discovery's roots, seeded into every sheet
anchorDir = fullfile(outDir, 'anchors');  if ~isfolder(anchorDir), mkdir(anchorDir); end
jobDir = fullfile(outDir, 'jobs');  if ~isfolder(jobDir), mkdir(jobDir); end

% ---- the driver's state: resumable -----------------------------------------
stateF = fullfile(outDir, 'torus_state.mat');
if isfile(stateF)
    st = load(stateF);  st = st.st;
    lg('resuming: round %d done, %d anchor(s)', st.roundsDone, size(st.anchors, 1));
else
    st = struct('roundsDone', 0, 'anchors', {spec.anchors}, 'rounds', struct([]));
    assert(size(st.anchors, 2) == 4 && size(st.anchors, 1) >= 1, 'run_phase_torus: .anchors is {name, file, sA, label; ...}');
    for k = 1:size(st.anchors, 1)
        assert(isfile(st.anchors{k, 2}), 'anchor %s: file %s missing', st.anchors{k, 1}, st.anchors{k, 2});
    end
    save(stateF, 'st');
end
lg('run_phase_torus: %d x %d phases, tag %s, %d anchor(s), out %s%s', numel(sD), numel(sA), tag, ...
   size(st.anchors, 1), outDir, tern(plan, ' [PLAN ONLY]', ''));

common = struct('sD', sD, 'sA', sA, 'thrustN', engine.thrustN, 'ispS', engine.ispS, 'm0kg', engine.m0kg, ...
                'tauDRO', orbits.tauDRO, 'NpTulip', orbits.NpTulip, 'pmTulip', orbits.pmTulip, ...
                'arcPattern', sprintf('arrival_arc_%s_*.mat', tag), 'librarySeeds', librarySeeds);
if isfile(seedFile), common.seedFiles = {seedFile}; end

out = struct('state', 'running', 'rounds', st.rounds, 'final', '', 'anchors', {st.anchors});
for k = st.roundsDone + 1 : maxRounds
    V = fullfile(outDir, sprintf('round_%02d', k));
    lg('===== ROUND %d (%s) =====', k, V);

    % ---- 1. arcs: both directions for every anchor that lacks them --------
    need = {};
    for a = 1:size(st.anchors, 1)
        for dirn = [-1 +1]
            f = arcFile(resDir, tag, st.anchors{a, 1}, dirn);
            if ~isfile(f), need(end+1, :) = {st.anchors{a, 1}, st.anchors{a, 2}, st.anchors{a, 3}, dirn, f}; end
        end
    end
    if plan
        lg('PLAN round %d: %d arc(s) to walk: %s', k, size(need, 1), strjoin(cellfun(@(n, dn) sprintf('%s_%s', n, tern(dn < 0, 'dn', 'up')), need(:, 1), need(:, 4), 'UniformOutput', false), ' '));
    elseif ~isempty(need)
        jobs = {};
        for q = 1:size(need, 1)
            jobs{end+1} = spawnArc(need{q, 1}, need{q, 2}, need{q, 3}, need{q, 4}, need{q, 5}, ...
                                   engine, orbits, sD(1), sA, arc, jobDir, matlabBin, lg);
        end
        okArcs = waitForFiles(need(:, 5), arc.deadlineSec + 1800, 120, lg);
        assert(okArcs, 'run_phase_torus: %d arc(s) did not finish (see %s)', nnz(~cellfun(@isfile, need(:, 5))), jobDir);
    end

    % ---- 2. the sheet, and the columns whose spine changed ----------------
    if plan
        lg('PLAN round %d: sheet at %d arrival phases from arcs %s + seeds%s; then ribs for changed columns on %d workers; then package/audit/sweep; then holes + improve (%.1f d); then discovery (%s)', ...
           k, numel(sA), common.arcPattern, tern(isfile(seedFile), ' + discovery seeds', ''), nWorkers, improveDays, tern(discover, 'on', 'off'));
        out.state = 'planned';  return
    end
    o1 = run_costate_library(setfields(common, struct('outDir', V, 'launch', false, ...
             'run', struct('sheet', true, 'ribs', false, 'package', false, 'audit', false, 'sweep', false))));
    S = load(o1.sheet);  S = S.S;
    prevSheet = '';  if k > 1, prevSheet = fullfile(outDir, sprintf('round_%02d', k - 1), fileName(o1.sheet)); end
    changed = true(1, numel(sA));
    if isfile(prevSheet)
        Sp = load(prevSheet);  Sp = Sp.S;
        changed = ~(isfinite(S.TF) & isfinite(Sp.TF) & abs(S.TF - Sp.TF) <= 1e-6) & isfinite(S.TF);
        for j = find(~changed & isfinite(S.TF))
            src = fullfile(fileparts(prevSheet), sprintf('fine_rib_col%02d.mat', j));
            dst = fullfile(V, sprintf('fine_rib_col%02d.mat', j));
            if isfile(src) && ~isfile(dst), copyfile(src, dst); end
        end
    end
    lg('round %d sheet: %d of %d columns certified; %d to walk: %s', k, nnz(isfinite(S.TF)), numel(sA), ...
       nnz(changed), mat2str(find(changed)));

    % ---- 3. ribs, package, audit, sweep (the supervised campaign) --------
    extra = earlierRibs(outDir, k);
    o2 = run_costate_library(setfields(common, struct('outDir', V, 'nWorkers', nWorkers, 'launch', true, ...
             'extraRibFiles', {extra}, ...
             'run', struct('sheet', false, 'ribs', true, 'package', true, 'audit', true, 'sweep', true))));
    lg('round %d launched: state %s%s', k, o2.state, tern(isempty(o2.blockers), '', [' -- ' strjoin(o2.blockers, ' | ')]));
    verdict = waitForVerdict(fullfile(V, 'SUPERVISOR_VERDICT.txt'), 48*3600, 60, lg);
    assert(contains(verdict, 'FINISHED'), 'run_phase_torus: round %d did not finish: %s', k, verdict);
    catMat = fullfile(V, 'costate_catalog_dro_tulip_70mN.mat');
    assert(isfile(catMat), 'run_phase_torus: round %d left no catalog', k);

    % ---- 4. holes and the improve pass, then a re-package if needed -------
    fh = fill_holes_direct(catMat, struct('logFile', fullfile(V, 'fill_holes_direct.log'), ...
             'improveDays', improveDays, 'anchorMat', st.anchors{1, 2}));
    lg('round %d holes/improve: %d holes, %d tried, %d certified', k, fh.nHoles, fh.nTried, fh.nCert);
    if fh.nCert > 0
        o3 = run_costate_library(setfields(common, struct('outDir', V, 'launch', false, ...
                 'extraRibFiles', {[extra, {fh.file}]}, ...
                 'run', struct('sheet', false, 'ribs', true, 'package', true, 'audit', true, 'sweep', true))));
        lg('round %d re-packaged with the direct cells: state %s', k, o3.state);
    end

    % ---- 5. discovery ----------------------------------------------------
    newAnchors = 0;
    if discover
        [roots, newAnchors, st.anchors] = discoverRoots(S, st.anchors, sD(1), engine, orbits, tag, ...
                                                       slowDays, seedFile, anchorDir, lg);
        if ~isempty(roots), common.seedFiles = {seedFile}; end
    end

    % ---- bookkeeping and the stopping rule -------------------------------
    rec = struct('round', k, 'dir', V, 'nCert', nnz(isfinite(S.TF)), 'changed', find(changed), ...
                 'newAnchors', newAnchors, 'catalog', catMat);
    if isempty(st.rounds), st.rounds = rec; else, st.rounds(end+1) = rec; end
    st.roundsDone = k;  save(stateF, 'st');
    out.rounds = st.rounds;  out.anchors = st.anchors;
    if newAnchors == 0 && (k > 1 && ~any(changed))
        lg('round %d changed nothing and found no new family: done', k);  break
    end
    if newAnchors == 0 && k > 1 && all(isfinite(S.TF)) && fh.nCert == 0
        lg('round %d: every column certified, nothing added, no new family: done', k);  break
    end
end

% ---- the library of record -------------------------------------------
if ~isempty(st.rounds)
    last = st.rounds(end);
    fin = fullfile(outDir, 'final');  if ~isfolder(fin), mkdir(fin); end
    for f = {'costate_catalog_dro_tulip_70mN.mat', 'catalog_receipt.mat', 'second_order_progress_v3.mat', ...
             'fine_rib_direct_holes.mat', 'phase_torus_70mN.png', 'phase_torus_findings.png'}
        if isfile(fullfile(last.dir, f{1})), copyfile(fullfile(last.dir, f{1}), fin); end
    end
    sh = dir(fullfile(last.dir, 'arrival_sheet_*.mat'));
    for q = 1:numel(sh), copyfile(fullfile(sh(q).folder, sh(q).name), fin); end
    fid = fopen(fullfile(fin, 'README.txt'), 'w');
    fprintf(fid, 'Phase-torus library, tag %s: %d x %d phases, %d round(s), built by run_phase_torus on %s.\nRounds: %s\nAnchors: %s\n', ...
            tag, numel(sD), numel(sA), numel(st.rounds), char(datetime('now')), ...
            strjoin(arrayfun(@(r) sprintf('%d (%d certified spines, %d new anchors)', r.round, r.nCert, r.newAnchors), st.rounds, 'UniformOutput', false), '; '), ...
            strjoin(st.anchors(:, 1).', ' '));
    fclose(fid);
    out.final = fin;  out.state = 'done';
    lg('library of record: %s', fin);
end
end

% ==========================================================================
function [roots, nNew, anchors] = discoverRoots(S, anchors, sD0, engine, orbits, tag, slowDays, seedFile, anchorDir, lg)
% DISCOVERROOTS  Branch-blind direct solves at the spine's empty and slow
% columns, seeded from other families' roots nearby; certified roots faster
% than the spine are appended to the seed file and made anchors.
% INPUTS: S (sheet); anchors; sD0; engine; orbits; tag; slowDays;
% seedFile; anchorDir; lg.  OUTPUTS: roots (certify_root outputs); nNew;
% anchors (grown).
roots = struct([]);  nNew = 0;
nA = numel(S.sA);
F = family_map(S, struct('arcDir', fullfile(fileparts(mfilename('fullpath')), 'results')));
fam = [F.columns.family];
tf = S.TF;
targets = find(~isfinite(tf));
for j = find(isfinite(tf))
    nb = [mod(j - 2, nA) + 1, mod(j, nA) + 1];  nb = nb(isfinite(tf(nb)));
    if ~isempty(nb) && tf(j) - min(tf(nb)) > slowDays, targets(end+1) = j; end
end
if isempty(targets), lg('discovery: no empty or slow column'); return, end
lg('discovery: %d column(s) to probe: %s', numel(targets), mat2str(targets));
so = struct('thrustN', engine.thrustN, 'ispS', engine.ispS, 'm0kg', engine.m0kg, 'tauDRO', orbits.tauDRO, ...
            'NpTulip', orbits.NpTulip, 'pmTulip', orbits.pmTulip, 'sD', sD0, 'anchorMat', anchors{1, 2});
[B, ~] = arclength_arrival('setup', so);
pool = capped_pool();
rv0 = B.stateD(sD0);  rv0 = rv0(1:6);
if isfile(seedFile), L = load(seedFile);  direct = L.direct; else, direct = struct('sD', {}, 'sA', {}, 'tfDays', {}, 'z', {}, 'src', {}); end
for j = targets
    % seeds: the fastest certified spine roots of OTHER families within three columns
    cand = [];
    for dj = [1 -1 2 -2 3 -3]
        jj = mod(j - 1 + dj, nA) + 1;
        if isfinite(tf(jj)) && (fam(jj) ~= fam(j) || fam(j) <= 0), cand(end+1) = jj; end
    end
    if isempty(cand), lg('  column %d (sA %.4f): no seed of another family nearby', j, S.sA(j)); continue, end
    [~, ord] = sort(tf(cand));  cand = cand(ord(1:min(3, end)));
    best = [];
    for jj = cand
        c = S.cand{jj};  kk = find([c.ok] & abs([c.tfDays] - tf(jj)) < 1e-9, 1);
        [C, info] = direct_cell_solve(c(kk).z, rv0, sD0, S.sA(j), B, struct('pool', pool));
        lg('  column %d (sA %.4f) from column %d (%.3f d): direct %.3f d, %s', j, S.sA(j), jj, tf(jj), info.tfDirectDays, C.reason);
        if C.ok && (isempty(best) || C.tfDays < best.tfDays), best = C; end
    end
    if isempty(best) || (isfinite(tf(j)) && best.tfDays >= tf(j) - 1e-3), continue, end
    name = sprintf('d%02d', j);
    direct(end+1) = struct('sD', sD0, 'sA', S.sA(j), 'tfDays', best.tfDays, 'z', best.z(:), ...
        'src', sprintf('run_phase_torus discovery at column %d, %s', j, char(datetime('now'))));
    save(seedFile, 'direct');
    anc = fullfile(anchorDir, sprintf('mintime_%s_anchor_%s.mat', tag, name));
    saveAnchor(anc, struct('z', best.z(:), 'it', struct('Y', best.Y), 'sA', S.sA(j), 'sD', sD0, ...
                           'tfDays', best.tfDays, 'origin', direct(end).src), B.Tnd, B.cnd);
    anchors(end+1, :) = {name, anc, S.sA(j), name};
    nNew = nNew + 1;
    if isempty(roots), roots = best; else, roots(end+1) = best; end
    lg('  NEW ROOT at column %d: %.3f d (spine %s) -> anchor %s', j, best.tfDays, tern(isfinite(tf(j)), sprintf('%.3f d', tf(j)), 'none'), name);
end
end

function saveAnchor(anc, best, Tnd, cnd)
% SAVEANCHOR  An anchor file in the layout arclength_arrival reads (best.z,
% best.it.Y, Tnd, cnd).  INPUTS: anc; best; Tnd; cnd.  OUTPUTS: none.
save(anc, 'best', 'Tnd', 'cnd');
end

function jobFile = spawnArc(name, anchorMat, sA0, dirn, outMat, engine, orbits, sD0, sA, arc, jobDir, matlabBin, lg)
% SPAWNARC  Write and launch one arc job (a batch MATLAB process).
% INPUTS: as named.  OUTPUTS: jobFile.
dn = tern(dirn < 0, 'dn', 'up');
jobFile = fullfile(jobDir, sprintf('arc_%s_%s.m', name, dn));
lgF = strrep(outMat, '.mat', '.log');  partF = strrep(outMat, '.mat', '.partial.mat');
levels = [sA - 1, sA, sA + 1];
txt = sprintf([ ...
 '%%%% arc job written by run_phase_torus\n' ...
 'here = pwd; cd(''/Users/msc/Desktop/proj7/external/pumpkynPie''); startup(); cd(here);\n' ...
 'addpath(%s, %s);  cd(%s);  capped_pool(2);\n' ...
 'so = struct(''thrustN'', %.17g, ''ispS'', %.17g, ''m0kg'', %.17g, ''tauDRO'', %.17g, ''NpTulip'', %d, ''pmTulip'', %d, ''sD'', %.17g, ''sA0'', %.17g, ''anchorMat'', %s);\n' ...
 '[~, anc] = arclength_arrival(''setup'', so);\n' ...
 'ao = so;  ao.direction = %d;  ao.sAStop = %.17g;  ao.levels = %s;  ao.nStep = %d;  ao.deadlineSec = %d;\n' ...
 'ao.logFile = %s;  ao.partialFile = %s;  ao.saveEvery = 50;\n' ...
 'A = arclength_arrival(anc, ao);  save(%s, ''A'', ''-v7.3'');\n' ...
 'if isfile(ao.partialFile), delete(ao.partialFile); end\n' ...
 'fprintf(''ARC %s DONE: %%d roots, sA %%.4f -> %%.4f, %%d folds, %%d crossings, stop = %%s\\n'', numel(A.q), A.q(1), A.q(end), numel(A.folds), numel(A.crossings), A.stop);\n'], ...
 mlq(fileparts(mfilename('fullpath'))), mlq(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'costate_common')), ...
 mlq(fileparts(mfilename('fullpath'))), engine.thrustN, engine.ispS, engine.m0kg, orbits.tauDRO, orbits.NpTulip, orbits.pmTulip, ...
 sD0, sA0, mlq(anchorMat), dirn, sA0 + dirn*arc.span, mat2str(levels, 17), arc.nStep, round(arc.deadlineSec), ...
 mlq(lgF), mlq(partF), mlq(outMat), sprintf('%s_%s', name, dn));
fid = fopen(jobFile, 'w');  fprintf(fid, '%s', txt);  fclose(fid);
outLog = strrep(jobFile, '.m', '.out');
cmd = sprintf('nohup %s -batch "run(''%s'')" > %s 2>&1 &', matlabBin, jobFile, outLog);
[st, msg] = system(cmd);
assert(st == 0, 'run_phase_torus: could not launch %s: %s', jobFile, msg);
lg('  arc %s_%s launched (%s)', name, dn, jobFile);
end

function ok = waitForFiles(files, deadlineSec, pollSec, lg)
% WAITFORFILES  Wait until every file exists or the deadline passes.
% INPUTS: files (cell); deadlineSec; pollSec; lg.  OUTPUTS: ok.
t0 = tic;
while true
    have = cellfun(@isfile, files);
    if all(have), ok = true;  return, end
    if toc(t0) > deadlineSec, ok = false;  lg('  wait: deadline passed with %d of %d files', nnz(have), numel(files));  return, end
    pause(pollSec);
end
end

function v = waitForVerdict(f, deadlineSec, pollSec, lg)
% WAITFORVERDICT  Wait for the supervisor's verdict file.  INPUTS: f;
% deadlineSec; pollSec; lg.  OUTPUTS: v (its last line, or 'DEADLINE').
t0 = tic;
while ~isfile(f)
    if toc(t0) > deadlineSec, v = 'DEADLINE';  lg('  wait: no verdict in %.0f h', deadlineSec/3600);  return, end
    pause(pollSec);
end
txt = strtrim(fileread(f));  lines = strsplit(txt, newline);  v = strtrim(lines{end});
lg('  verdict: %s', v);
end

function x = earlierRibs(outDir, k)
% EARLIERRIBS  Every rib file of rounds 1..k-1 (and their direct-cell ribs).
% INPUTS: outDir; k.  OUTPUTS: x (cellstr).
x = {};
for q = 1:k-1
    V = fullfile(outDir, sprintf('round_%02d', q));
    f = [dir(fullfile(V, 'fine_rib_col*.mat')); dir(fullfile(V, 'fine_rib_direct_holes.mat'))];
    x = [x, fullfile(V, {f.name})];
end
end

function f = arcFile(resDir, tag, name, dirn)
% ARCFILE  The arc's file of record.  INPUTS: resDir; tag; name; dirn.
% OUTPUTS: f.
f = fullfile(resDir, sprintf('arrival_arc_%s_%s_%s_long.mat', tag, name, tern(dirn < 0, 'dn', 'up')));
end

function s = setfields(s, t)
% SETFIELDS  Copy t's fields onto s.  INPUTS: s; t.  OUTPUTS: s.
for f = fieldnames(t)', s.(f{1}) = t.(f{1}); end
end

function s = withDefaults(s, dflt)
% WITHDEFAULTS  Fill missing fields from dflt.  INPUTS: s; dflt.  OUTPUTS: s.
for f = fieldnames(dflt)', if ~isfield(s, f{1}), s.(f{1}) = dflt.(f{1}); end, end
end

function n = fileName(p)
% FILENAME  Name + extension of a path.  INPUTS: p.  OUTPUTS: n.
[~, a, b] = fileparts(p);  n = [a b];
end

function p = absPath(p)
% ABSPATH  Absolute path.  INPUTS: p.  OUTPUTS: p.
if ~startsWith(p, '/'), p = fullfile(pwd, p); end
end

function q = mlq(s)
% MLQ  A MATLAB single-quoted literal.  INPUTS: s.  OUTPUTS: q.
q = ['''' strrep(s, '''', '''''') ''''];
end

function s = tern(c, a, b)
% TERN  Ternary.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end

function logmsg(f, s)
% LOGMSG  Append to the log and echo.  INPUTS: f; s.  OUTPUTS: none.
line = sprintf('%s %s', char(datetime('now', 'Format', 'HH:mm:ss')), s);
fid = fopen(f, 'a');  fprintf(fid, '%s\n', line);  fclose(fid);  fprintf('%s\n', line);
end
