function out = run_phase_torus(spec)
%% Purpose:
%
%   Build a minimum-time costate library over a USER-CHOSEN set of
%   departure and arrival phases -- the phase torus -- by the round loop
%   that built the 70 mN DRO -> tulip 24 x 24 library by hand (FINDINGS
%   59-72, process/PHASE_TORUS_RUNBOOK.md, doc/phase_torus_methods.tex):
%
%     round k:  arcs for every anchor without them (batch jobs, both
%               arrival directions, partial saves, a verdict file each)
%            -> the arrival sheet at the listed phases, from every arc and
%               every registered root; columns whose spine ROOT changed
%            -> ribs for those columns (the supervised queue), every
%               earlier round's ribs offered to the packager; package,
%               audit, sweep (the finalizer)
%            -> holes and the improve pass (fill_holes_direct); a faster
%               root found ON THE SPINE is registered; re-package
%            -> DISCOVERY: at every empty column and every column whose
%               spine is far slower than a neighbour's (or every column,
%               .probeAll), direct solves seeded from other families'
%               certified roots nearby; every distinct certified root is
%               registered, and one faster than the spine that no known
%               family passes through becomes a new anchor
%     until a round registers no new root and adds no anchor (a fixed
%     point), or .maxRounds is reached (then the status is 'budget', not
%     'done', and the pending anchors are named).
%
%   Every stage leaves its product on disk before the next starts, every
%   wait has a deadline and reads a verdict FILE, the arcs of a campaign
%   live in the campaign's own folder, and the driver's state -- the
%   campaign manifest, the anchors, the rounds, the status -- is in
%   <outDir>/torus_state.mat, saved after every promotion, so a killed
%   driver resumes where it was and a finished one says so. .plan = true
%   prints what the next round would do and launches nothing.
%
%   Reviewed by GPT-6 Astra (xhigh) on 2026-09-15; the adjudication is in
%   FINDINGS 73 and reviews/run_phase_torus_astra_2026-09-15.md.
%
%% Inputs:
%
%  spec                     struct
%   .sD, .sA                [1 x nD], [1 x nA]     phases in [0,1), strictly
%                                                   increasing, no two closer
%                                                   than 1e-5 (give a lattice
%                                                   that wraps past 1 as
%                                                   sort(mod(., 1))); sD(1)
%                                                   is the spine's phase
%   .orbits                 struct                  .tauDRO [1] .NpTulip [7]
%                                                   .pmTulip [-1]
%   .engine                 struct                  .thrustN [0.070] .ispS
%                                                   [900] .m0kg [150]
%   .anchors                {name, file, sA, label; ...}  certified roots at
%                                                   sD(1) (best.z, best.it.Y)
%   .outDir                 char                    campaign root
%   .tag                    char                    names the arcs ['torus']
%   .arc                    struct                  .nStep [4000] .deadlineSec
%                                                   [6*3600] .span [1.15]
%   .rib                    struct                  .wallSec [900] per point
%   .nWorkers               int                     rib workers [4]
%   .improveDays, .slowDays double                  [2], [2]
%   .acceptDays             double                  a probe must beat the
%                                                   spine by this to be
%                                                   promoted [0.05]
%   .seedRadius             double                  circular phase distance
%                                                   within which other
%                                                   families' roots seed a
%                                                   probe [0.15]
%   .probeAll               logical                 probe every column, not
%                                                   only empty/slow [false]
%   .maxRounds              int                     [6]
%   .discover               logical                 [true]
%   .librarySeeds           logical                 seed from the 70 mN
%                                                   library too [false]
%   .matlab                 char                    batch MATLAB binary [R2026a]
%   .startup                char                    folder whose startup()
%                                                   builds the path [pumpkynPie]
%   .plan                   logical                 print, launch nothing [false]
%
%% Outputs:
%
%  out                      struct                  .state ('done' | 'budget'
%                                                   | 'planned' | 'nothing to
%                                                   do') .reason .rounds
%                                                   .final .anchors
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
assert(min([diff(sD), 1 - (sD(end) - sD(1))]) > 1e-5 && min([diff(sA), 1 - (sA(end) - sA(1))]) > 1e-5, ...
       'run_phase_torus: phases closer than 1e-5 of a period cannot be resolved');
orbits = withDefaults(fieldd(spec, 'orbits', struct()), struct('tauDRO', 1.0, 'NpTulip', 7, 'pmTulip', -1));
engine = withDefaults(fieldd(spec, 'engine', struct()), struct('thrustN', 0.070, 'ispS', 900, 'm0kg', 150));
outDir = absPath(spec.outDir);  if ~isfolder(outDir), mkdir(outDir); end
tag = d('tag', 'torus');
arc = withDefaults(fieldd(spec, 'arc', struct()), struct('nStep', 4000, 'deadlineSec', 6*3600, 'span', 1.15));
rib = withDefaults(fieldd(spec, 'rib', struct()), struct('wallSec', 900));
nWorkers = d('nWorkers', 4);  improveDays = d('improveDays', 2);  slowDays = d('slowDays', 2);
acceptDays = d('acceptDays', 0.05);  seedRadius = d('seedRadius', 0.15);  probeAll = d('probeAll', false);
maxRounds = d('maxRounds', 6);  discover = d('discover', true);  plan = d('plan', false);
librarySeeds = d('librarySeeds', false);
matlabBin = d('matlab', '/Applications/MATLAB_R2026a.app/bin/matlab');
startupDir = d('startup', '/Users/msc/Desktop/proj7/external/pumpkynPie');
logF = fullfile(outDir, 'torus.log');
lg = @(varargin) logmsg(logF, sprintf(varargin{:}));
arcDir = fullfile(outDir, 'arcs');  if ~isfolder(arcDir), mkdir(arcDir); end   % THIS campaign's arcs
seedFile = fullfile(outDir, 'direct_certified.mat');   % the registry: every distinct certified root found off the arcs
anchorDir = fullfile(outDir, 'anchors');  if ~isfolder(anchorDir), mkdir(anchorDir); end
jobDir = fullfile(outDir, 'jobs');  if ~isfolder(jobDir), mkdir(jobDir); end
manifest = struct('sD', sD, 'sA', sA, 'orbits', orbits, 'engine', engine, 'tag', tag);

% ---- the driver's state: identity-checked, resumable, terminal-aware ------
stateF = fullfile(outDir, 'torus_state.mat');
if isfile(stateF)
    st = load(stateF);  st = st.st;
    for f = fieldnames(manifest)'
        assert(isequaln(st.manifest.(f{1}), manifest.(f{1})), 'run_phase_torus:manifest', ...
               ['%s differs from this campaign''s manifest (%s): a torus is one grid, one engine, one orbit ' ...
                'pair, one tag -- use another outDir'], f{1}, outDir);
    end
    lg('resuming: status %s, %d round(s) done, %d anchor(s)', st.status, st.roundsDone, size(st.anchors, 1));
else
    anchors = spec.anchors;
    assert(size(anchors, 2) == 4 && size(anchors, 1) >= 1, 'run_phase_torus: .anchors is {name, file, sA, label; ...}');
    for k = 1:size(anchors, 1)
        anchors{k, 2} = absPath(anchors{k, 2});
        assert(isfile(anchors{k, 2}), 'anchor %s: file %s missing', anchors{k, 1}, anchors{k, 2});
    end
    assert(numel(unique(anchors(:, 1))) == size(anchors, 1), 'run_phase_torus: anchor names must be unique');
    st = struct('manifest', manifest, 'roundsDone', 0, 'anchors', {anchors}, 'rounds', struct([]), ...
                'status', 'running', 'reason', '', 'jobs', struct('pid', {}, 'file', {}));
    saveState(stateF, st);
end
out = struct('state', st.status, 'reason', st.reason, 'rounds', st.rounds, 'final', '', 'anchors', {st.anchors});
if any(strcmp(st.status, {'done', 'budget'}))
    lg('this campaign already ended: %s (%s). Nothing to do; use another outDir or delete torus_state.mat to rebuild.', st.status, st.reason);
    out.state = 'nothing to do';  out.final = fullfile(outDir, 'final');
    return
end
lg('run_phase_torus: %d x %d phases, tag %s, %d anchor(s), out %s%s', numel(sD), numel(sA), tag, ...
   size(st.anchors, 1), outDir, tern(plan, ' [PLAN ONLY]', ''));
if isfile(seedFile), L0 = load(seedFile);  lg('registry: %d certified root(s) registered off the arcs', numel(L0.direct)); end

common = struct('sD', sD, 'sA', sA, 'thrustN', engine.thrustN, 'ispS', engine.ispS, 'm0kg', engine.m0kg, ...
                'tauDRO', orbits.tauDRO, 'NpTulip', orbits.NpTulip, 'pmTulip', orbits.pmTulip, ...
                'arcPattern', sprintf('arrival_arc_%s_*.mat', tag), 'arcDir', arcDir, ...
                'librarySeeds', librarySeeds, 'ribWallSec', rib.wallSec);
if isfile(seedFile), common.seedFiles = {seedFile}; end

for k = st.roundsDone + 1 : maxRounds
    V = fullfile(outDir, sprintf('round_%02d', k));
    lg('===== ROUND %d (%s) =====', k, V);

    % ---- 1. arcs: both directions for every anchor that lacks them --------
    need = {};
    for a = 1:size(st.anchors, 1)
        for dirn = [-1 +1]
            f = arcFile(arcDir, tag, st.anchors{a, 1}, dirn);
            if ~isfile(f), need(end+1, :) = {st.anchors{a, 1}, st.anchors{a, 2}, st.anchors{a, 3}, dirn, f}; end
        end
    end
    if plan
        lg('PLAN round %d: %d arc(s) to walk: %s', k, size(need, 1), strjoin(cellfun(@(n, dn) sprintf('%s_%s', n, tern(dn < 0, 'dn', 'up')), need(:, 1), need(:, 4), 'UniformOutput', false), ' '));
        lg('PLAN round %d: sheet at %d arrival phases from %s%s; ribs for changed columns on %d workers (%.0f s per point); package/audit/sweep; holes + improve (%.1f d); discovery %s', ...
           k, numel(sA), common.arcPattern, tern(isfile(seedFile), ' + registered roots', ''), nWorkers, rib.wallSec, improveDays, tern(discover, tern(probeAll, '(every column)', '(empty and slow columns)'), 'off'));
        out.state = 'planned';  return
    end
    if ~isempty(need)
        st = runArcs(need, st, stateF, engine, orbits, sD(1), sA, arc, jobDir, matlabBin, startupDir, here, lg);
    end

    % ---- 2. the sheet, and the columns whose spine ROOT changed ------------
    o1 = run_costate_library(setfields(common, struct('outDir', V, 'launch', false, ...
             'run', struct('sheet', true, 'ribs', false, 'package', false, 'audit', false, 'sweep', false))));
    S = load(o1.sheet);  S = S.S;
    changed = isfinite(S.TF);
    prevSheet = '';  if k > 1, prevSheet = fullfile(outDir, sprintf('round_%02d', k - 1), fileName(o1.sheet)); end
    if isfile(prevSheet)
        Sp = load(prevSheet);  Sp = Sp.S;
        same = isfinite(S.TF) & isfinite(Sp.TF) & abs(S.TF - Sp.TF) <= 1e-6 & all(abs(S.Z8 - Sp.Z8) <= 1e-9*max(1, abs(Sp.Z8)), 1);
        changed = isfinite(S.TF) & ~same;
        for j = find(same)
            src = fullfile(fileparts(prevSheet), sprintf('fine_rib_col%02d.mat', j));
            dst = fullfile(V, sprintf('fine_rib_col%02d.mat', j));
            if isfile(src) && ~isfile(dst), copyfile(src, dst); end
        end
        lost = isfinite(Sp.TF) & ~isfinite(S.TF);
        if any(lost), lg('round %d: column(s) %s LOST their certified spine (regression) -- earlier ribs still package', k, mat2str(find(lost))); end
    end
    % a rib copied into THIS round by an earlier attempt for a spine that has
    % since changed must not satisfy the walk: set it aside
    for j = find(changed)
        f = fullfile(V, sprintf('fine_rib_col%02d.mat', j));
        if isfile(f) && ~ribMatchesSpine(f, S, j)
            movefile(f, sprintf('%s.stale_%s', f, char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
            lg('round %d: stale rib for column %d set aside (its spine changed)', k, j);
        end
    end
    lg('round %d sheet: %d of %d columns certified; %d to walk: %s', k, nnz(isfinite(S.TF)), numel(sA), nnz(changed), mat2str(find(changed)));

    % ---- 3. ribs, package, audit, sweep (the supervised campaign) --------
    extra = earlierRibs(outDir, k);
    verdict = runRound(V, common, nWorkers, extra, lg);
    assert(contains(verdict, 'FINISHED') || contains(verdict, 'packaged'), 'run_phase_torus: round %d did not finish: %s', k, verdict);
    catMat = fullfile(V, 'costate_catalog_dro_tulip_70mN.mat');
    assert(isfile(catMat), 'run_phase_torus: round %d left no catalog', k);

    % ---- 4. holes and the improve pass; spine roots registered; re-package -
    fh = fill_holes_direct(catMat, struct('logFile', fullfile(V, 'fill_holes_direct.log'), ...
             'improveDays', improveDays, 'anchorMat', st.anchors{1, 2}));
    lg('round %d holes/improve: %d holes, %d tried, %d certified', k, fh.nHoles, fh.nTried, fh.nCert);
    nRegistered = 0;
    if fh.nCert > 0
        [nReg, st] = registerSpineRoots(fh.file, S, sD(1), seedFile, st, stateF, arcDir, tag, anchorDir, acceptDays, lg);
        nRegistered = nRegistered + nReg;
        o3 = run_costate_library(setfields(common, struct('outDir', V, 'launch', false, ...
                 'extraRibFiles', {[extra, {fh.file}]}, ...
                 'run', struct('sheet', false, 'ribs', true, 'package', true, 'audit', true, 'sweep', true))));
        assert(strcmp(o3.state, 'packaged'), 'run_phase_torus: round %d re-package failed: %s', k, strjoin(o3.blockers, ' | '));
        lg('round %d re-packaged with the direct cells', k);
    end

    % ---- 5. discovery ----------------------------------------------------
    newAnchors = 0;
    if discover
        [nReg, newAnchors, st] = discoverRoots(S, st, stateF, sD(1), engine, orbits, tag, slowDays, acceptDays, ...
                                              seedRadius, probeAll, seedFile, anchorDir, arcDir, lg);
        nRegistered = nRegistered + nReg;
    end
    if isfile(seedFile), common.seedFiles = {seedFile}; end

    % ---- bookkeeping and the fixed-point rule ----------------------------
    rec = struct('round', k, 'dir', V, 'nCert', nnz(isfinite(S.TF)), 'changed', find(changed), ...
                 'newAnchors', newAnchors, 'newRoots', nRegistered, 'catalog', catMat);
    if isempty(st.rounds), st.rounds = rec; else, st.rounds(end+1) = rec; end
    st.roundsDone = k;
    if newAnchors == 0 && nRegistered == 0
        st.status = 'done';  st.reason = sprintf('round %d registered no new root and added no anchor: fixed point', k);
    elseif k == maxRounds
        st.status = 'budget';
        st.reason = sprintf('round %d was the last allowed (maxRounds = %d) and still added %d anchor(s) / %d root(s): raise .maxRounds and rerun', ...
                            k, maxRounds, newAnchors, nRegistered);
    end
    saveState(stateF, st);
    out.rounds = st.rounds;  out.anchors = st.anchors;
    lg('round %d: %d certified spines, %d walked, %d new root(s) registered, %d new anchor(s); status %s', ...
       k, rec.nCert, numel(rec.changed), nRegistered, newAnchors, st.status);
    if ~strcmp(st.status, 'running'), break, end
end

% ---- the library of record: a fresh folder, the last round's products --
if ~isempty(st.rounds)
    last = st.rounds(end);
    fin = fullfile(outDir, 'final');
    if isfolder(fin), rmdir(fin, 's'); end
    mkdir(fin);
    required = {'costate_catalog_dro_tulip_70mN.mat', 'catalog_receipt.mat'};
    optional = {'second_order_progress_v3.mat', 'fine_rib_direct_holes.mat', 'phase_torus_70mN.png', 'phase_torus_findings.png'};
    for f = required
        assert(isfile(fullfile(last.dir, f{1})), 'run_phase_torus: the last round lacks %s', f{1});
        copyfile(fullfile(last.dir, f{1}), fin);
    end
    missing = {};
    for f = optional
        if isfile(fullfile(last.dir, f{1})), copyfile(fullfile(last.dir, f{1}), fin); else, missing{end+1} = f{1}; end
    end
    sh = dir(fullfile(last.dir, 'arrival_sheet_*.mat'));
    for q = 1:numel(sh), copyfile(fullfile(sh(q).folder, sh(q).name), fin); end
    fid = fopen(fullfile(fin, 'README.txt'), 'w');
    fprintf(fid, ['Phase-torus library, tag %s: %d x %d phases, %d round(s), status %s (%s), built by run_phase_torus on %s.\n' ...
                  'Rounds: %s\nAnchors: %s\nNot produced: %s\n'], ...
            tag, numel(sD), numel(sA), numel(st.rounds), st.status, st.reason, char(datetime('now')), ...
            strjoin(arrayfun(@(r) sprintf('%d (%d certified spines, %d new roots, %d new anchors)', r.round, r.nCert, r.newRoots, r.newAnchors), st.rounds, 'UniformOutput', false), '; '), ...
            strjoin(st.anchors(:, 1).', ' '), tern(isempty(missing), 'nothing', strjoin(missing, ', ')));
    fclose(fid);
    out.final = fin;
    lg('library of record: %s (status %s)', fin, st.status);
end
out.state = st.status;  out.reason = st.reason;
end

% ==========================================================================
function st = runArcs(need, st, stateF, engine, orbits, sD0, sA, arc, jobDir, matlabBin, startupDir, here, lg)
% RUNARCS  Spawn one batch job per needed arc, record the PIDs, wait for
% each job's VERDICT file (.done / .fail), kill what is left at the deadline.
% INPUTS: as named.  OUTPUTS: st (jobs recorded).
files = need(:, 5);
for q = 1:size(need, 1)
    [pid, jobFile] = spawnArc(need{q, 1}, need{q, 2}, need{q, 3}, need{q, 4}, need{q, 5}, ...
                              engine, orbits, sD0, sA, arc, jobDir, matlabBin, startupDir, here, lg);
    st.jobs(end+1) = struct('pid', pid, 'file', jobFile);
end
saveState(stateF, st);
t0 = tic;  deadline = arc.deadlineSec + 1800;
while true
    done = cellfun(@(f) isfile([f '.done']), files);
    failed = cellfun(@(f) isfile([f '.fail']), files);
    if any(failed)
        f = files{find(failed, 1)};
        error('run_phase_torus:arc', 'arc %s FAILED: %s', fileName(f), strtrim(fileread([f '.fail'])));
    end
    if all(done), break, end
    if toc(t0) > deadline
        for j = 1:numel(st.jobs), system(sprintf('kill %d 2>/dev/null', st.jobs(j).pid)); end
        error('run_phase_torus:arc', '%d arc(s) did not finish within %.0f h; their jobs were killed', nnz(~done), deadline/3600);
    end
    pause(120);
end
st.jobs = struct('pid', {}, 'file', {});
saveState(stateF, st);
lg('  %d arc(s) done', numel(files));
end

function [pid, jobFile] = spawnArc(name, anchorMat, sA0, dirn, outMat, engine, orbits, sD0, sA, arc, jobDir, matlabBin, startupDir, here, lg)
% SPAWNARC  Write and launch one arc job. The job saves through a temp
% file and a rename, then writes <outMat>.done, or <outMat>.fail with the
% error.  INPUTS: as named.  OUTPUTS: pid; jobFile.
dn = tern(dirn < 0, 'dn', 'up');
jobFile = fullfile(jobDir, sprintf('arc_%s_%s.m', name, dn));
lgF = strrep(outMat, '.mat', '.log');  partF = strrep(outMat, '.mat', '.partial.mat');
kmin = floor(min(sA) + min(sA0 - arc.span, 0)) - 1;  kmax = ceil(max(sA) + max(sA0 + arc.span, 1)) + 1;
levels = reshape(sA(:) + (kmin:kmax), 1, []);            % every whole-period copy the walk can reach
txt = sprintf([ ...
 '%%%% arc job written by run_phase_torus\n' ...
 'here = pwd; cd(%s); startup(); cd(here);\n' ...
 'addpath(%s, %s);  cd(%s);  capped_pool(2);\n' ...
 'outMat = %s;\n' ...
 'try\n' ...
 '  so = struct(''thrustN'', %.17g, ''ispS'', %.17g, ''m0kg'', %.17g, ''tauDRO'', %.17g, ''NpTulip'', %d, ''pmTulip'', %d, ''sD'', %.17g, ''sA0'', %.17g, ''anchorMat'', %s);\n' ...
 '  [~, anc] = arclength_arrival(''setup'', so);\n' ...
 '  ao = so;  ao.direction = %d;  ao.sAStop = %.17g;  ao.levels = %s;  ao.nStep = %d;  ao.deadlineSec = %d;\n' ...
 '  ao.logFile = %s;  ao.partialFile = %s;  ao.saveEvery = 50;\n' ...
 '  A = arclength_arrival(anc, ao);\n' ...
 '  tmp = [outMat ''.part''];  save(tmp, ''A'', ''-v7.3'');  movefile(tmp, outMat, ''f'');\n' ...
 '  if isfile(ao.partialFile), delete(ao.partialFile); end\n' ...
 '  fid = fopen([outMat ''.done''], ''w'');  fprintf(fid, ''%%d roots, sA %%.4f -> %%.4f, %%d folds, %%d crossings, stop = %%s\\n'', numel(A.q), A.q(1), A.q(end), numel(A.folds), numel(A.crossings), A.stop);  fclose(fid);\n' ...
 '  fprintf(''ARC %s DONE\\n'');\n' ...
 'catch ME\n' ...
 '  fid = fopen([outMat ''.fail''], ''w'');  fprintf(fid, ''%%s\\n'', ME.message);  fclose(fid);\n' ...
 '  fprintf(''ARC %s FAILED: %%s\\n'', ME.message);\n' ...
 'end\n'], ...
 mlq(startupDir), mlq(here), mlq(fullfile(fileparts(fileparts(here)), 'costate_common')), mlq(here), mlq(outMat), ...
 engine.thrustN, engine.ispS, engine.m0kg, orbits.tauDRO, orbits.NpTulip, orbits.pmTulip, ...
 sD0, sA0, mlq(anchorMat), dirn, sA0 + dirn*arc.span, mat2str(levels, 17), arc.nStep, round(arc.deadlineSec), ...
 mlq(lgF), mlq(partF), sprintf('%s_%s', name, dn), sprintf('%s_%s', name, dn));
fid = fopen(jobFile, 'w');  fprintf(fid, '%s', txt);  fclose(fid);
outLog = strrep(jobFile, '.m', '.out');
cmd = sprintf('nohup %s -batch %s > %s 2>&1 & echo $!', shq(matlabBin), shq(sprintf('run(''%s'')', jobFile)), shq(outLog));
[st, msg] = system(cmd);
pid = str2double(strtrim(msg));
assert(st == 0 && isfinite(pid), 'run_phase_torus: could not launch %s: %s', jobFile, msg);
lg('  arc %s_%s launched (pid %d, %s)', name, dn, pid, jobFile);
end

function verdict = runRound(V, common, nWorkers, extra, lg)
% RUNROUND  Launch the supervised rib campaign for a round (or adopt one
% already running for this folder) and wait for a verdict newer than any on
% disk; a call that packages by itself needs no wait.  INPUTS: V; common;
% nWorkers; extra; lg.  OUTPUTS: verdict (char).
verdictF = fullfile(V, 'SUPERVISOR_VERDICT.txt');
before = 0;  if isfile(verdictF), before = dir(verdictF).datenum; end
[~, live] = system(sprintf('pgrep -f %s', shq(sprintf('campaign_supervisor.sh.*%s', V))));
if ~isempty(regexp(strtrim(live), '^\d+', 'once'))
    lg('  a supervisor for %s is still running (pid %s): adopting it instead of launching another', V, strtrim(live));
else
    o2 = run_costate_library(setfields(common, struct('outDir', V, 'nWorkers', nWorkers, 'launch', true, ...
             'extraRibFiles', {extra}, ...
             'run', struct('sheet', false, 'ribs', true, 'package', true, 'audit', true, 'sweep', true))));
    lg('  launched: state %s%s', o2.state, tern(isempty(o2.blockers), '', [' -- ' strjoin(o2.blockers, ' | ')]));
    if strcmp(o2.state, 'packaged'), verdict = 'packaged in this call';  return, end
    assert(strcmp(o2.state, 'launched'), 'run_phase_torus: round could not launch: %s', strjoin(o2.blockers, ' | '));
end
verdict = waitForVerdict(verdictF, before, 48*3600, 60, lg);
end

function [nReg, st] = registerSpineRoots(holesFile, S, sD0, seedFile, st, stateF, arcDir, tag, anchorDir, acceptDays, lg)
% REGISTERSPINEROOTS  Direct-solved cells ON THE SPINE are roots the next
% sheet must see; one faster than its column's spine by acceptDays that no
% known family passes through becomes an anchor.  INPUTS: as named.
% OUTPUTS: nReg (roots registered); st.
nReg = 0;
H = load(holesFile);
F = [];
for m = 1:numel(H.R)
    for q = 1:numel(H.R(m).pts)
        p = H.R(m).pts(q);
        if abs(mod(p.sD - sD0 + 0.5, 1) - 0.5) > 1e-8, continue, end
        j = find(abs(mod(S.sA - p.sA + 0.5, 1) - 0.5) < 1e-8, 1);
        if isempty(j), continue, end
        added = registerRoot(seedFile, p, sD0, sprintf('direct cell solve on the spine, round %d', st.roundsDone + 1));
        nReg = nReg + added;
        if added && (~isfinite(S.TF(j)) || p.tfDays < S.TF(j) - acceptDays)
            if isempty(F), F = family_map(S, struct('arcDir', arcDir)); end
            [famP, ~] = F.attach(S.sA(j), p.tfDays);
            if famP >= 1
                lg('  spine cell at column %d: %.3f d lies on family %s -- registered, not anchored', j, p.tfDays, F.families(famP).label);
            else
                st = promote(st, stateF, p, j, S, arcDir, tag, anchorDir, lg, 'spine cell');
            end
        end
    end
end
end

function [nReg, nNew, st] = discoverRoots(S, st, stateF, sD0, engine, orbits, tag, slowDays, acceptDays, seedRadius, probeAll, seedFile, anchorDir, arcDir, lg)
% DISCOVERROOTS  Branch-blind direct solves at the spine's target columns,
% seeded from other families' certified roots within a phase radius (every
% certified candidate of those columns, not only the winners); every
% distinct certified root is registered; one faster than the spine by
% acceptDays that no known family passes through is promoted.
% INPUTS: as named.  OUTPUTS: nReg; nNew; st.
nReg = 0;  nNew = 0;
nA = numel(S.sA);  tf = S.TF;
F = family_map(S, struct('arcDir', arcDir));
fam = [F.columns.family];
if probeAll, targets = 1:nA; else
    targets = find(~isfinite(tf));
    for j = find(isfinite(tf))
        nb = [mod(j - 2, nA) + 1, mod(j, nA) + 1];  nb = nb(isfinite(tf(nb)));
        if ~isempty(nb) && tf(j) - min(tf(nb)) > slowDays, targets(end+1) = j; end
    end
end
targets = unique(targets, 'stable');
if isempty(targets), lg('discovery: no target column'); return, end
lg('discovery: %d column(s) to probe: %s', numel(targets), mat2str(targets));
so = struct('thrustN', engine.thrustN, 'ispS', engine.ispS, 'm0kg', engine.m0kg, 'tauDRO', orbits.tauDRO, ...
            'NpTulip', orbits.NpTulip, 'pmTulip', orbits.pmTulip, 'sD', sD0, 'anchorMat', st.anchors{1, 2});
[B, ~] = arclength_arrival('setup', so);
pool = capped_pool();
rv0 = B.stateD(sD0);  rv0 = rv0(1:6);
for j = targets
    % the seed pool: every certified candidate at columns within seedRadius
    % (circular), of a family other than the target's (any family when the
    % target is empty or unattached), fastest first, distinct roots, three
    dist = abs(mod(S.sA - S.sA(j) + 0.5, 1) - 0.5);
    near = find(dist <= seedRadius & (1:nA) ~= j);
    seeds = struct('z', {}, 'tf', {}, 'col', {});
    for jj = near
        c = S.cand{jj};
        if isempty(c), continue, end
        for kk = find([c.ok])
            [famC, ~] = F.attach(S.sA(jj), c(kk).tfDays);
            if fam(j) >= 1 && famC == fam(j), continue, end
            if any(arrayfun(@(p) abs(p.tf - c(kk).tfDays) < 1e-3 && p.col == jj, seeds)), continue, end
            seeds(end+1) = struct('z', c(kk).z(:), 'tf', c(kk).tfDays, 'col', jj);
        end
    end
    if isempty(seeds), lg('  column %d (sA %.4f): no seed of another family within %.2f', j, S.sA(j), seedRadius); continue, end
    [~, ord] = sort([seeds.tf]);  seeds = seeds(ord(1:min(3, end)));
    best = [];
    for p = seeds
        try
            [C, info] = direct_cell_solve(p.z, rv0, sD0, S.sA(j), B, struct('pool', pool));
        catch ME
            C = struct('ok', false, 'reason', ['threw: ' ME.message]);  info = struct('tfDirectDays', NaN);
        end
        lg('  column %d (sA %.4f) from column %d (%.3f d): direct %.3f d, %s', j, S.sA(j), p.col, p.tf, info.tfDirectDays, C.reason);
        if ~C.ok, continue, end
        C.note = join_note(sprintf('discovery probe at column %d from column %d (%.3f d)', j, p.col, p.tf), C);
        added = registerRoot(seedFile, C, sD0, sprintf('discovery probe at column %d, round %d', j, st.roundsDone + 1));
        nReg = nReg + added;
        if isempty(best) || C.tfDays < best.tfDays, best = C; end
    end
    if isempty(best) || (isfinite(tf(j)) && best.tfDays >= tf(j) - acceptDays), continue, end
    [famB, ~] = F.attach(S.sA(j), best.tfDays);
    if famB >= 1
        lg('  column %d: the %.3f d root lies on family %s already -- registered, not anchored', j, best.tfDays, F.families(famB).label);
        continue
    end
    st = promote(st, stateF, best, j, S, arcDir, tag, anchorDir, lg, 'discovery');
    nNew = nNew + 1;
end
end

function st = promote(st, stateF, C, j, S, arcDir, tag, anchorDir, lg, how)
% PROMOTE  A certified root becomes an anchor with a unique, immutable name;
% the state is saved at once.  INPUTS: as named.  OUTPUTS: st.
base = sprintf('d%02d', j);
nSame = nnz(startsWith(st.anchors(:, 1), base));
name = sprintf('%s_%d', base, nSame + 1);
anc = fullfile(anchorDir, sprintf('mintime_%s_anchor_%s.mat', tag, name));
assert(~isfile(anc), 'run_phase_torus: anchor file %s exists already', anc);
Tnd = S.B.Tnd;  cnd = S.B.cnd;
best = struct('z', C.z(:), 'it', struct('Y', C.Y), 'sA', S.sA(j), 'sD', C.sD, 'tfDays', C.tfDays, ...
              'origin', sprintf('%s at column %d, %s', how, j, char(datetime('now'))));
save(anc, 'best', 'Tnd', 'cnd');
st.anchors(end+1, :) = {name, anc, S.sA(j), name};
saveState(stateF, st);
lg('  NEW ANCHOR %s: %.3f d at sA %.4f (%s; spine %s); arcs %s', name, C.tfDays, S.sA(j), how, ...
   tern(isfinite(S.TF(j)), sprintf('%.3f d', S.TF(j)), 'none'), fileName(arcFile(arcDir, tag, name, -1)));
end

function added = registerRoot(seedFile, C, sD0, src)
% REGISTERROOT  Append a certified root to the campaign's registry unless
% an equal one (same sA, t_f within 1e-3 d) is there; atomic write.
% INPUTS: seedFile; C (certify_root output); sD0; src.  OUTPUTS: added
% (logical).
if isfile(seedFile), L = load(seedFile);  direct = L.direct;
else, direct = struct('sD', {}, 'sA', {}, 'tfDays', {}, 'z', {}, 'src', {});
end
dup = ~isempty(direct) && any(abs(mod([direct.sA] - C.sA + 0.5, 1) - 0.5) < 1e-8 & abs([direct.tfDays] - C.tfDays) < 1e-3);
added = ~dup;
if added
    direct(end+1) = struct('sD', sD0, 'sA', C.sA, 'tfDays', C.tfDays, 'z', C.z(:), 'src', src);
    tmp = [seedFile '.part'];  save(tmp, 'direct');  movefile(tmp, seedFile, 'f');
end
end

function ok = ribMatchesSpine(ribFile, S, j)
% RIBMATCHESSPINE  Does a rib file's first point continue this sheet's spine
% root at column j (t_f within 0.5 d of the spine)?  INPUTS: ribFile; S; j.
% OUTPUTS: ok.
ok = false;
try
    L = load(ribFile);  r = L.R(1);
    if isempty(r.pts) || ~isfinite(S.TF(j)), return, end
    ok = abs(r.pts(1).tfDays - S.TF(j)) < 0.5;
catch
end
end

function v = waitForVerdict(f, before, deadlineSec, pollSec, lg)
% WAITFORVERDICT  Wait for a supervisor verdict file NEWER than `before`
% (a datenum; 0 = any).  INPUTS: f; before; deadlineSec; pollSec; lg.
% OUTPUTS: v (its last line, or 'DEADLINE').
t0 = tic;
while ~isfile(f) || dir(f).datenum <= before
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

function f = arcFile(arcDir, tag, name, dirn)
% ARCFILE  The arc's file of record.  INPUTS: arcDir; tag; name; dirn.
% OUTPUTS: f.
f = fullfile(arcDir, sprintf('arrival_arc_%s_%s_%s_long.mat', tag, name, tern(dirn < 0, 'dn', 'up')));
end

function saveState(stateF, st)
% SAVESTATE  Atomic write of the driver's state.  INPUTS: stateF; st.
tmp = [stateF '.part'];  save(tmp, 'st');  movefile(tmp, stateF, 'f');
end

function s = join_note(prefix, C)
% JOIN_NOTE  This producer's clause in front of the certifier's note.
% INPUTS: prefix; C.  OUTPUTS: s.
parts = {prefix};
if isfield(C, 'note') && ~isempty(C.note), parts{end+1} = C.note; end
s = strjoin(parts, ' | ');
end

function s = setfields(s, t)
% SETFIELDS  Copy t's fields onto s.  INPUTS: s; t.  OUTPUTS: s.
for f = fieldnames(t)', s.(f{1}) = t.(f{1}); end
end

function s = withDefaults(s, dflt)
% WITHDEFAULTS  Fill missing fields from dflt.  INPUTS: s; dflt.  OUTPUTS: s.
for f = fieldnames(dflt)', if ~isfield(s, f{1}) || isempty(s.(f{1})), s.(f{1}) = dflt.(f{1}); end, end
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

function q = shq(s)
% SHQ  A POSIX-shell single-quoted literal.  INPUTS: s.  OUTPUTS: q.
q = ['''' strrep(s, '''', '''\''''') ''''];
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
