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
%  spec                     struct                  every field below is read
%                                                   in section 0 (the input
%                                                   table); unknown fields
%                                                   are refused; call with no
%                                                   arguments for an example
%   .sD, .sA                [1 x nD], [1 x nA]     phases in [0,1), strictly
%                                                   increasing, no two closer
%                                                   than 1e-5 (give a lattice
%                                                   that wraps past 1 as
%                                                   sort(mod(., 1))); sD(1)
%                                                   is the spine's phase
%   .departure              struct                  the departure orbit:
%                                                   .family ['dro'] .tau [1.0]
%                                                   (period, ND)
%   .arrival                struct                  the arrival orbit:
%                                                   .family ['tulip'] .Np [7]
%                                                   .pm [-1]; its period is
%                                                   locked by Np
%   .engine                 struct                  the engine: .thrustN
%                                                   [0.070] N, .ispS [900] s,
%                                                   .m0kg [150] kg
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

%% ========================================================================
%  0. USER INPUTS -- everything a campaign IS, in one place. Every field of
%     `spec` is read, defaulted and checked here and nowhere else; a field
%     that is not in the table is refused by name (a typo must not become
%     a silently ignored setting). Call with no arguments for an example.
%% ========================================================================
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
if nargin == 0, out = exampleSpec(here);  return, end
% TEST SEAM: handles to this file's local functions, by name
% (tests/test_run_phase_torus_p0)
if ischar(spec) && strcmp(spec, 'localfunctions'), out = localHandles(localfunctions);  return, end
in = userInputs(spec);

% the names the body uses (read-only from here on)
sD = in.sD;  sA = in.sA;  orbits = in.orbits;  engine = in.engine;  tag = in.tag;
outDir = in.outDir;  arc = in.arc;  rib = in.rib;  nWorkers = in.nWorkers;
improveDays = in.improveDays;  slowDays = in.slowDays;  acceptDays = in.acceptDays;
seedRadius = in.seedRadius;  probeAll = in.probeAll;  maxRounds = in.maxRounds;
discover = in.discover;  plan = in.plan;  librarySeeds = in.librarySeeds;
matlabBin = in.matlab;  startupDir = in.startup;

%% ========================================================================
%  1. THE CAMPAIGN FOLDER -- its arcs, registry, anchors, jobs, log, state
%% ========================================================================
if ~isfolder(outDir), mkdir(outDir); end
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
    st = struct('manifest', manifest, 'roundsDone', 0, 'anchors', {in.anchors}, 'rounds', struct([]), ...
                'status', 'running', 'reason', '', 'jobs', struct('pid', {}, 'file', {}));
    saveState(stateF, st);
end
out = struct('state', st.status, 'reason', st.reason, 'rounds', st.rounds, 'final', '', 'anchors', {st.anchors});
[ended, whyEnded] = campaignEnded(st, maxRounds);
if ended
    lg('%s', whyEnded);
    out.state = 'nothing to do';  out.final = fullfile(outDir, 'final');
    return
end
if strcmp(st.status, 'budget')               % a larger .maxRounds reopens a budget-limited campaign
    lg('continuing a budget-limited campaign: %d round(s) done, .maxRounds is now %d', st.roundsDone, maxRounds);
    st.status = 'running';  st.reason = '';
    saveState(stateF, st);
end
lg('run_phase_torus: %d x %d phases, tag %s, %d anchor(s), out %s%s', numel(sD), numel(sA), tag, ...
   size(st.anchors, 1), outDir, tern(plan, ' [PLAN ONLY]', ''));
if isfile(seedFile), L0 = load(seedFile);  lg('registry: %d certified root(s) registered off the arcs', numel(L0.direct)); end

common = commonOpts(in, st.anchors, arcDir, sprintf('arrival_arc_%s_*.mat', tag));
if in.roundDeadlineHours > 0, common.roundDeadlineHours = in.roundDeadlineHours; end
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
        for jc = find(same)                         % jc: an arrival column
            src = fullfile(fileparts(prevSheet), sprintf('fine_rib_col%02d.mat', jc));
            dst = fullfile(V, sprintf('fine_rib_col%02d.mat', jc));
            if isfile(src) && ~isfile(dst), copyfile(src, dst); end
        end
        lost = isfinite(Sp.TF) & ~isfinite(S.TF);
        if any(lost), lg('round %d: column(s) %s LOST their certified spine (regression) -- earlier ribs still package', k, mat2str(find(lost))); end
    end
    % a rib copied into THIS round by an earlier attempt for a spine that has
    % since changed must not satisfy the walk: set it aside
    for jc = find(changed)
        f = fullfile(V, sprintf('fine_rib_col%02d.mat', jc));
        if isfile(f) && ~ribMatchesSpine(f, S, jc)
            movefile(f, sprintf('%s.stale_%s', f, char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
            lg('round %d: stale rib for column %d set aside (its spine changed)', k, jc);
        end
    end
    lg('round %d sheet: %d of %d columns certified; %d to walk: %s', k, nnz(isfinite(S.TF)), numel(sA), nnz(changed), mat2str(find(changed)));

    % ---- 3. ribs, package, audit, sweep (the supervised campaign) --------
    % earlier rounds' ribs, AND this round's own holes file when an earlier
    % attempt of this round already wrote one (a resumed round must package
    % the cells that attempt certified)
    extra = roundExtras(outDir, k);
    verdict = runRound(V, common, nWorkers, extra, lg);
    assert(contains(verdict, 'FINISHED') || contains(verdict, 'packaged'), 'run_phase_torus: round %d did not finish: %s', k, verdict);
    catMat = fullfile(V, 'costate_catalog_dro_tulip_70mN.mat');
    assert(isfile(catMat), 'run_phase_torus: round %d left no catalog', k);

    % ---- 4. holes and the improve pass; spine roots registered; re-package -
    fh = fill_holes_direct(catMat, struct('logFile', fullfile(V, 'fill_holes_direct.log'), ...
             'improveDays', improveDays));
    lg('round %d holes/improve: %d holes, %d tried, %d certified', k, fh.nHoles, fh.nTried, fh.nCert);
    % THE FILE DECIDES, not this call's counter. On a resumed round the
    % filler skips the cells an earlier attempt certified and reports
    % nCert = 0, yet those roots are on disk: they are registered here from
    % the file, and were offered to the packager above (roundExtras).
    nRegistered = 0;  nFillerAnchors = 0;
    if isfile(fh.file)
        [nReg, nFillerAnchors, st] = registerSpineRoots(fh.file, S, sD(1), seedFile, st, stateF, arcDir, tag, anchorDir, acceptDays, lg);
        nRegistered = nRegistered + nReg;
    end
    if fh.nCert > 0                              % new cells this call: the catalog must take them in
        o3 = run_costate_library(setfields(common, struct('outDir', V, 'launch', false, ...
                 'extraRibFiles', {unique([extra, {fh.file}], 'stable')}, ...
                 'run', struct('sheet', false, 'ribs', true, 'package', true, 'audit', true, 'sweep', true))));
        assertPackaged(o3, fileName(V), 're-package with the direct cells');
        lg('round %d re-packaged with the direct cells', k);
    end

    % ---- 5. discovery ----------------------------------------------------
    newAnchors = nFillerAnchors;                 % anchors come from the filler AND from discovery
    if discover
        [nReg, nDiscovered, st] = discoverRoots(S, st, stateF, sD(1), engine, orbits, tag, slowDays, acceptDays, ...
                                                seedRadius, probeAll, seedFile, anchorDir, arcDir, lg);
        nRegistered = nRegistered + nReg;  newAnchors = newAnchors + nDiscovered;
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
function in = userInputs(spec)
% USERINPUTS  The user-input table: name, default ([] = required), meaning.
% Unknown fields are refused by name, required ones must be present, the
% rest take their defaults; then the values are checked and a summary is
% printed.  INPUTS: spec.  OUTPUTS: in (every field of the table).
T = { ...
 % name           default                                         meaning
 'sD',            [],   'departure phases in [0,1), strictly increasing; sD(1) is the spine every arc and rib starts from (one phase = no ribs)'; ...
 'sA',            [],   'arrival phases in [0,1), strictly increasing (a lattice that wraps past 1: sort(mod(., 1)); one phase = one column)'; ...
 'anchors',       [],   '{name, file, sA, label; ...}: certified roots at sD(1) the arcs start from (best.z, best.it.Y)'; ...
 'outDir',        [],   'the campaign folder (rounds, arcs, registry, anchors, state, final)'; ...
 'departure',     struct('family', 'dro', 'tau', 1.0), ...
                        'the DEPARTURE orbit: pumpkyn family name and period (ND); a DRO of period tau. Only ''dro'' is wired today'; ...
 'arrival',       struct('family', 'tulip', 'Np', 7, 'pm', -1), ...
                        'the ARRIVAL orbit: family, petal count Np, branch pm (+1/-1); a tulip''s period is locked by Np, 2*pi*(Np-2)/(Np-1). Only ''tulip'' is wired today'; ...
 'engine',        struct('thrustN', 0.070, 'ispS', 900, 'm0kg', 150), ...
                        'the ENGINE: thrust thrustN (N, constant, always on -- minimum time), specific impulse ispS (s), initial mass m0kg (kg); 0.070 / 900 / 150 is the shipped 70 mN regime'; ...
 'tag',           'torus', 'names the campaign''s arcs (arrival_arc_<tag>_...)'; ...
 'arc',           struct('nStep', 4000, 'deadlineSec', 6*3600, 'span', 1.15), ...
                        'arc budget: steps, wall deadline, arrival-phase span from the anchor'; ...
 'rib',           struct('wallSec', 900), 'rib budget: wall cap per departure point'; ...
 'roundDeadlineHours', 0, 'how long to wait for a round''s ribs + finalizer; 0 = derive it from the grid and the workers (roundDeadlineSec)'; ...
 'nWorkers',      4,    'rib workers (one MATLAB each)'; ...
 'maxRounds',     6,    'rounds before the campaign stops with status ''budget'''; ...
 'improveDays',   2,    'improve pass: re-solve a filled cell slower than a column neighbour by more than this (days)'; ...
 'slowDays',      2,    'discovery: probe a column whose spine is slower than a neighbour''s by more than this (days)'; ...
 'acceptDays',    0.05, 'discovery: a probe must beat the spine by this (days) to become an anchor'; ...
 'seedRadius',    0.15, 'discovery: other families'' roots within this circular phase distance seed a probe'; ...
 'probeAll',      false, 'discovery: probe every column, not only the empty and slow ones'; ...
 'discover',      true, 'run the discovery step at all'; ...
 'librarySeeds',  false, 'also seed the sheet from the 70 mN library''s own roots (dro_tulip_library)'; ...
 'matlab',        '/Applications/MATLAB_R2026a.app/bin/matlab', 'the MATLAB that runs the batch jobs (R2026a: has the Parallel Computing Toolbox)'; ...
 'startup',       '/Users/msc/Desktop/proj7/external/pumpkynPie', 'folder whose startup() builds the path in every batch job'; ...
 'plan',          false, 'print what the next round would do and launch nothing'};
names = T(:, 1);
extra = setdiff(fieldnames(spec), names);
assert(isempty(extra), 'run_phase_torus:input', 'unknown input(s): %s. The inputs are: %s', ...
       strjoin(extra, ', '), strjoin(names, ', '));
in = struct();
for k = 1:size(T, 1)
    nm = T{k, 1};  dflt = T{k, 2};
    if isfield(spec, nm) && ~isempty(spec.(nm))
        v = spec.(nm);
        if isstruct(dflt), v = withDefaults(v, dflt); end       % a partial struct keeps the other defaults
        in.(nm) = v;
    else
        assert(~isempty(dflt), 'run_phase_torus:input', 'required input .%s is missing: %s', nm, T{k, 3});
        in.(nm) = dflt;
    end
end
% ---- checks --------------------------------------------------------------
in.sD = in.sD(:).';  in.sA = in.sA(:).';
for ax = {'sD', 'sA'}
    v = in.(ax{1});
    assert(isnumeric(v) && numel(v) >= 1 && all(v >= 0 & v < 1) && all(diff(v) > 0), 'run_phase_torus:input', ...
           '.%s must be one or more phases in [0,1), strictly increasing', ax{1});
    assert(numel(v) < 2 || min([diff(v), 1 - (v(end) - v(1))]) > 1e-5, 'run_phase_torus:input', ...
           '.%s has phases closer than 1e-5 of a period, which the sheet and rib matchers cannot resolve', ax{1});
end
assert(iscell(in.anchors) && size(in.anchors, 2) == 4 && size(in.anchors, 1) >= 1, 'run_phase_torus:input', ...
       '.anchors must be a cell table {name, file, sA, label; ...} with at least one row');
for k = 1:size(in.anchors, 1)
    in.anchors{k, 2} = absPath(in.anchors{k, 2});
    assert(isfile(in.anchors{k, 2}), 'run_phase_torus:input', 'anchor %s: file %s missing', in.anchors{k, 1}, in.anchors{k, 2});
end
assert(numel(unique(in.anchors(:, 1))) == size(in.anchors, 1), 'run_phase_torus:input', 'anchor names must be unique');
in.outDir = absPath(in.outDir);
assert(strcmpi(in.departure.family, 'dro') && strcmpi(in.arrival.family, 'tulip'), 'run_phase_torus:input', ...
       'the pipeline is wired for a DRO departure and a tulip arrival (got %s -> %s)', in.departure.family, in.arrival.family);
assert(in.departure.tau > 0 && in.arrival.Np >= 3 && in.arrival.Np == round(in.arrival.Np) && any(in.arrival.pm == [-1 1]), ...
       'run_phase_torus:input', '.departure.tau > 0; .arrival.Np an integer >= 3; .arrival.pm +1 or -1');
assert(in.engine.thrustN > 0 && in.engine.ispS > 0 && in.engine.m0kg > 0, 'run_phase_torus:input', ...
       '.engine.thrustN (N), .ispS (s) and .m0kg (kg) must be positive');
in.orbits = struct('tauDRO', in.departure.tau, 'NpTulip', in.arrival.Np, 'pmTulip', in.arrival.pm);   % what the chain reads
in.arrival.tau = 2*pi*(in.arrival.Np - 2)/(in.arrival.Np - 1);                                        % derived, never chosen
assert(in.nWorkers >= 1 && in.maxRounds >= 1 && in.arc.nStep >= 1 && in.arc.deadlineSec > 0 && in.rib.wallSec > 0, ...
       'run_phase_torus:input', 'budgets must be positive');
assert(in.acceptDays >= 0 && in.seedRadius > 0 && in.seedRadius <= 0.5, 'run_phase_torus:input', ...
       '.acceptDays >= 0 and 0 < .seedRadius <= 0.5');
% ---- the summary ---------------------------------------------------------
fprintf('PHASE TORUS\n');
fprintf('  grid      : %d departure x %d arrival phases (sD %s; sA %s)\n', numel(in.sD), numel(in.sA), phaseList(in.sD), phaseList(in.sA));
tStar = 382981.289129055;                                       % s per ND time unit (Earth-Moon)
fprintf('  departure : DRO, period %.4f ND (%.2f d)\n', in.departure.tau, in.departure.tau*tStar/86400);
fprintf('  arrival   : %d-petal tulip, branch %+d, period %.4f ND (%.2f d, locked by Np)\n', in.arrival.Np, in.arrival.pm, in.arrival.tau, in.arrival.tau*tStar/86400);
g0 = 9.80665;
fprintf('  engine    : thrust %.0f mN, Isp %g s (c = %.3f km/s), initial mass %g kg (%.4f mm/s^2 at m0)\n', ...
        in.engine.thrustN*1000, in.engine.ispS, in.engine.ispS*g0/1000, in.engine.m0kg, in.engine.thrustN/in.engine.m0kg*1000);
fprintf('  anchors   : %s\n', strjoin(cellfun(@(n, a) sprintf('%s (sA %.4f)', n, a), in.anchors(:, 1), in.anchors(:, 3), 'UniformOutput', false), ', '));
fprintf('  budgets   : arcs %d steps / %.1f h, ribs %.0f s per point, %d workers, %d rounds\n', in.arc.nStep, in.arc.deadlineSec/3600, in.rib.wallSec, in.nWorkers, in.maxRounds);
fprintf('  search    : improve > %.1f d, probe columns slower by > %.1f d%s, accept a gain > %.2f d, seeds within %.2f\n', ...
        in.improveDays, in.slowDays, tern(in.probeAll, ' (and every other column)', ''), in.acceptDays, in.seedRadius);
fprintf('  output    : %s (tag %s)%s\n', in.outDir, in.tag, tern(in.plan, '  [PLAN ONLY]', ''));
end

function s = phaseList(v)
% PHASELIST  A phase vector as text, abbreviated when long.  INPUTS: v.
% OUTPUTS: s.
if numel(v) <= 6, s = mat2str(v, 4);
else, s = sprintf('[%s ... %s] (%d)', strjoin(arrayfun(@(x) sprintf('%.4f', x), v(1:3), 'UniformOutput', false), ' '), ...
                  strjoin(arrayfun(@(x) sprintf('%.4f', x), v(end-1:end), 'UniformOutput', false), ' '), numel(v));
end
end

function spec = exampleSpec(here)
% EXAMPLESPEC  The 70 mN DRO -> tulip 24 x 24 torus as a spec, printed and
% returned (nothing runs).  INPUTS: here.  OUTPUTS: spec.
spec = struct( ...
    'sD',      (0:23)/24, ...                               % departure phases
    'sA',      sort(mod(0.0754 + (0:23)/24, 1)), ...        % arrival phases (the 24-lattice from the first anchor's phase)
    'departure', struct('family', 'dro', 'tau', 1.0), ...          % the DRO: period 1.0 ND
    'arrival',   struct('family', 'tulip', 'Np', 7, 'pm', -1), ...  % the tulip: 7 petals, branch -1
    'engine',  struct('thrustN', 0.070, 'ispS', 900, 'm0kg', 150), ...
    'anchors', {{'anchor', fullfile(here, 'results', 'mintime_70mN_anchor.mat'), 0.0754, 'fast'}}, ...
    'outDir',  fullfile(here, 'results', 'torus_70mN_24x24'), ...
    'tag',     '70mN', ...
    'nWorkers', 4, 'maxRounds', 6, 'plan', true);
fprintf('run_phase_torus: an example spec (plan mode). Edit and call run_phase_torus(spec).\n');
disp(spec);
end

% ==========================================================================
function st = runArcs(need, st, stateF, engine, orbits, sD0, sA, arc, jobDir, matlabBin, startupDir, here, lg)
% RUNARCS  Make every needed arc exist. A job an earlier driver left RUNNING
% is adopted (never a second writer on one file); the rest are spawned.
% Then wait on each arc's VERDICT file (.done / .fail): a .fail stops the
% campaign at once, and at the deadline whatever is left is killed.
% INPUTS: as named (need rows: name, anchor file, anchor phase, direction,
% arc file).  OUTPUTS: st (job record emptied once every arc is done).
files = need(:, 5);
[adopted, st] = adoptLiveJobs(need, st, lg);
for q = find(~adopted(:).')
    pid = spawnArc(need{q, 1}, need{q, 2}, need{q, 3}, need{q, 4}, need{q, 5}, ...
                   engine, orbits, sD0, sA, arc, jobDir, matlabBin, startupDir, here, lg);
    st.jobs(end+1) = struct('pid', pid, 'file', need{q, 5});
end
saveState(stateF, st);
t0 = tic;  deadline = arc.deadlineSec + 1800;
while true
    done = cellfun(@(f) isfile([f '.done']), files);
    failed = cellfun(@(f) isfile([f '.fail']), files);
    if any(failed)
        f = files{find(failed, 1)};
        why = strtrim(fileread([f '.fail']));
        st = killJobs(st, stateF);               % its siblings must not run on unowned
        error('run_phase_torus:arc', 'arc %s FAILED: %s', fileName(f), why);
    end
    if all(done), break, end
    if toc(t0) > deadline
        st = killJobs(st, stateF);
        error('run_phase_torus:arc', '%d arc(s) did not finish within %.0f h; their jobs were killed', nnz(~done), deadline/3600);
    end
    pause(120);
end
st.jobs = struct('pid', {}, 'file', {});
saveState(stateF, st);
lg('  %d arc(s) done', numel(files));
end

function [adopted, st] = adoptLiveJobs(need, st, lg)
% ADOPTLIVEJOBS  Which needed arcs already have a LIVE job from an earlier
% driver? Those are adopted -- waited on, not spawned again: two MATLABs
% walking one arc would race on one output file. A recorded job that is dead,
% or whose arc is no longer needed, leaves the record.
% INPUTS: need (rows as in runArcs); st (.jobs: pid, file = the arc file);
% lg.  OUTPUTS: adopted [nNeed x 1 logical]; st (dead jobs removed).
adopted = false(size(need, 1), 1);
keep = false(1, numel(st.jobs));
for q = 1:numel(st.jobs)
    job = st.jobs(q);
    hit = find(strcmp(need(:, 5), job.file), 1);
    alive = pidAlive(job.pid) || pidAlive(arcMatlabPid(job.file));
    if alive && ~isempty(hit)
        adopted(hit) = true;  keep(q) = true;
        lg('  arc %s: a job from an earlier driver is still running (pid %d) -- adopted, not respawned', fileName(job.file), job.pid);
    end
end
st.jobs = st.jobs(keep);
end

function st = killJobs(st, stateF)
% KILLJOBS  Stop every recorded arc job: the shell wrapper and the MATLAB it
% started (its pid is in <arc>.pid). The record is emptied and saved.
% INPUTS: st; stateF.  OUTPUTS: st.
for q = 1:numel(st.jobs)
    for pid = [st.jobs(q).pid, arcMatlabPid(st.jobs(q).file)]
        if isfinite(pid), system(sprintf('kill %d 2>/dev/null', pid)); end
    end
end
st.jobs = struct('pid', {}, 'file', {});
saveState(stateF, st);
end

function tf = pidAlive(pid)
% PIDALIVE  Is this process alive (signal 0 reaches it)?  INPUTS: pid.
% OUTPUTS: tf logical (false for NaN).
tf = isscalar(pid) && isfinite(pid) && pid > 0 && system(sprintf('kill -0 %d 2>/dev/null', pid)) == 0;
end

function pid = arcMatlabPid(arcFile)
% ARCMATLABPID  The pid of the MATLAB an arc's wrapper started, from
% <arcFile>.pid; NaN when there is none.  INPUTS: arcFile.  OUTPUTS: pid.
pid = NaN;
if isfile([arcFile '.pid']), pid = str2double(strtrim(fileread([arcFile '.pid']))); end
end

function [pid, jobFile] = spawnArc(name, anchorMat, sA0, dirn, outMat, engine, orbits, sD0, sA, arc, jobDir, matlabBin, startupDir, here, lg)
% SPAWNARC  Write and launch ONE arc job, as three readable files in jobDir:
%   arc_<name>_<dn|up>.m    the MATLAB job: walk the arc, save it through a
%                           temp file and a rename, then write <arc>.done --
%                           or <arc>.fail with the error
%   arc_<name>_<dn|up>.sh   a shell wrapper that starts MATLAB, records its
%                           pid in <arc>.pid, and GUARANTEES a verdict: if
%                           MATLAB exits having written neither file (a
%                           licence, path or pool failure before the job's own
%                           try/catch) the wrapper writes the .fail
%   arc_<name>_<dn|up>.out  MATLAB's stdout
% Verdict files of an EARLIER attempt are removed first: a stale .fail used
% to stop the retry on its first poll.
% INPUTS: as named.  OUTPUTS: pid (the wrapper's); jobFile.
dn = tern(dirn < 0, 'dn', 'up');
stem = fullfile(jobDir, sprintf('arc_%s_%s', name, dn));
jobFile = [stem '.m'];  shFile = [stem '.sh'];  outLog = [stem '.out'];
[arcFolder, arcBase] = fileparts(outMat);
lgF = fullfile(arcFolder, [arcBase '.log']);  partF = fullfile(arcFolder, [arcBase '.partial.mat']);
for ext = {'.done', '.fail', '.pid'}
    if isfile([outMat ext{1}]), delete([outMat ext{1}]); end
end

% ---- the MATLAB job ------------------------------------------------------
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
writeText(jobFile, txt);

% ---- the shell wrapper ---------------------------------------------------
sh = sprintf([ ...
 '#!/bin/sh\n' ...
 '# Written by run_phase_torus: run ONE arc job and guarantee a verdict file.\n' ...
 'ARC=%s\n' ...
 '%s -batch %s > %s 2>&1 &\n' ...
 'echo $! > "$ARC.pid"\n' ...
 'wait $!\n' ...
 'rc=$?\n' ...
 'if [ ! -f "$ARC.done" ] && [ ! -f "$ARC.fail" ]; then\n' ...
 '  echo "MATLAB exited (code $rc) without writing a verdict; see its log in the jobs folder" > "$ARC.fail"\n' ...
 'fi\n'], ...
 shq(outMat), shq(matlabBin), shq(sprintf('run(%s)', mlq(jobFile))), shq(outLog));
writeText(shFile, sh);

[rc, msg] = system(sprintf('nohup /bin/sh %s > /dev/null 2>&1 & echo $!', shq(shFile)));
pid = str2double(strtrim(msg));
assert(rc == 0 && isfinite(pid), 'run_phase_torus: could not launch %s: %s', shFile, msg);
lg('  arc %s_%s launched (pid %d, %s)', name, dn, pid, jobFile);
end

function writeText(file, txt)
% WRITETEXT  Write a text file, or say which one could not be written.
% INPUTS: file; txt.  OUTPUTS: none.
fid = fopen(file, 'w');
assert(fid >= 0, 'run_phase_torus: cannot write %s', file);
fprintf(fid, '%s', txt);  fclose(fid);
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
    if strcmp(o2.state, 'packaged')
        assertPackaged(o2, fileName(V), 'in-call package');
        verdict = 'packaged in this call';  return
    end
    assert(strcmp(o2.state, 'launched'), 'run_phase_torus: round could not launch: %s', strjoin(o2.blockers, ' | '));
end
verdict = waitForVerdict(verdictF, before, roundDeadlineSec(numel(common.sD), numel(common.sA), nWorkers, fieldd(common, 'roundDeadlineHours', [])), 60, lg);
end

function [nReg, nNew, st] = registerSpineRoots(holesFile, S, sD0, seedFile, st, stateF, arcDir, tag, anchorDir, acceptDays, lg)
% REGISTERSPINEROOTS  Direct-solved cells ON THE SPINE are roots the next
% sheet must see, so each is registered; one that passes promotionVerdict
% becomes an anchor. Works from the holes FILE, so cells certified by an
% earlier attempt of this round are registered too.
% INPUTS: as named.  OUTPUTS: nReg (roots registered); nNew (anchors
% added); st.
nReg = 0;  nNew = 0;
H = load(holesFile);
F = [];                                          % the family map: built only if a root needs it
for m = 1:numel(H.R)
    for q = 1:numel(H.R(m).pts)
        p = H.R(m).pts(q);
        onSpine = abs(mod(p.sD - sD0 + 0.5, 1) - 0.5) <= 1e-8;
        jc = find(abs(mod(S.sA - p.sA + 0.5, 1) - 0.5) < 1e-8, 1);      % its arrival column
        if ~onSpine || isempty(jc), continue, end
        added = registerRoot(seedFile, p, sD0, sprintf('direct cell solve on the spine, round %d', st.roundsDone + 1));
        nReg = nReg + added;
        if ~added, continue, end
        if isempty(F), F = family_map(S, struct('arcDir', arcDir)); end
        [famP, ~] = F.attach(S.sA(jc), p.tfDays, p.z);           % by flight time AND costates
        [promoteIt, why] = promotionVerdict(added, p.tfDays, S.TF(jc), acceptDays, famP);
        if promoteIt
            st = promote(st, stateF, p, jc, S, arcDir, tag, anchorDir, lg, 'spine cell');
            nNew = nNew + 1;
        else
            lg('  spine cell at column %d (%.3f d): registered, not anchored -- %s', jc, p.tfDays, why);
        end
    end
end
end

function [promoteIt, why] = promotionVerdict(added, tfDays, tfSpine, acceptDays, famCode)
% PROMOTIONVERDICT  THE rule for turning a certified root into an anchor,
% used by both paths (filler and discovery). All three must hold:
%   1. the root was NEWLY registered -- a root already in the registry was
%      judged when it arrived and is never anchored twice;
%   2. it beats its column's spine by acceptDays (any root beats an empty
%      column);
%   3. no known family passes through it -- a root on a walked family is a
%      seed, and walking its arcs again would buy nothing.
% INPUTS: added (logical); tfDays; tfSpine (NaN = empty column); acceptDays;
% famCode (family_map attachment: >= 1 is a known family).
% OUTPUTS: promoteIt (logical); why (char, the deciding clause).
promoteIt = false;
if ~added
    why = 'it was already registered';
elseif isfinite(tfSpine) && ~(tfDays < tfSpine - acceptDays)
    why = sprintf('it does not beat the %.3f d spine by %.2f d', tfSpine, acceptDays);
elseif famCode >= 1
    why = sprintf('it lies on known family %d, whose arcs are already walked', famCode);
else
    promoteIt = true;  why = 'new, faster than the spine, on no known family';
end
end

function [nReg, nNew, st] = discoverRoots(S, st, stateF, sD0, engine, orbits, tag, slowDays, acceptDays, seedRadius, probeAll, seedFile, anchorDir, arcDir, lg)
% DISCOVERROOTS  Branch-blind direct solves at the spine's target columns,
% seeded from other families' certified roots within a phase radius (every
% certified candidate of those columns, not only the winners). Every
% distinct certified root is registered; the fastest at each column is
% promoted if promotionVerdict says so.
% INPUTS: as named.  OUTPUTS: nReg; nNew; st.
nReg = 0;  nNew = 0;
nA = numel(S.sA);  tf = S.TF;
F = family_map(S, struct('arcDir', arcDir));
fam = [F.columns.family];

% ---- the target columns: empty ones, and ones far slower than a neighbour
if probeAll, targets = 1:nA; else
    targets = find(~isfinite(tf));
    for jc = find(isfinite(tf))
        nb = [mod(jc - 2, nA) + 1, mod(jc, nA) + 1];  nb = nb(isfinite(tf(nb)));
        if ~isempty(nb) && tf(jc) - min(tf(nb)) > slowDays, targets(end+1) = jc; end
    end
end
targets = unique(targets, 'stable');
if isempty(targets), lg('discovery: no target column'); return, end
lg('discovery: %d column(s) to probe: %s', numel(targets), mat2str(targets));

% the closures and the engine only: a probe has its own seed and needs no
% anchor polished (physicsOpts)
[B, ~] = arclength_arrival('setup', physicsOpts(engine, orbits, sD0));
pool = capped_pool();
rv0 = B.stateD(sD0);  rv0 = rv0(1:6);
for jc = targets
    % the seed pool: every certified candidate at columns within seedRadius
    % (circular), of a family other than the target's (any family when the
    % target is empty or unattached), fastest first, distinct roots, three
    dist = abs(mod(S.sA - S.sA(jc) + 0.5, 1) - 0.5);
    near = find(dist <= seedRadius & (1:nA) ~= jc);
    seeds = struct('z', {}, 'tf', {}, 'col', {});
    for jn = near                                % jn: a neighbouring column
        c = S.cand{jn};
        if isempty(c), continue, end
        for kk = find([c.ok])
            [famC, ~] = F.attach(S.sA(jn), c(kk).tfDays, c(kk).z);
            if fam(jc) >= 1 && famC == fam(jc), continue, end
            if any(arrayfun(@(p) p.col == jn && sameRoot(p.z, c(kk).z), seeds)), continue, end
            seeds(end+1) = struct('z', c(kk).z(:), 'tf', c(kk).tfDays, 'col', jn);
        end
    end
    if isempty(seeds), lg('  column %d (sA %.4f): no seed of another family within %.2f', jc, S.sA(jc), seedRadius); continue, end
    [~, ord] = sort([seeds.tf]);  seeds = seeds(ord(1:min(3, end)));

    % ---- probe from each seed; register every certified root ------------
    best = [];  bestAdded = false;
    for p = seeds
        try
            [C, info] = direct_cell_solve(p.z, rv0, sD0, S.sA(jc), B, struct('pool', pool));
        catch ME
            C = struct('ok', false, 'reason', ['threw: ' ME.message]);  info = struct('tfDirectDays', NaN);
        end
        lg('  column %d (sA %.4f) from column %d (%.3f d): direct %.3f d, %s', jc, S.sA(jc), p.col, p.tf, info.tfDirectDays, C.reason);
        if ~C.ok, continue, end
        C.note = join_note(sprintf('discovery probe at column %d from column %d (%.3f d)', jc, p.col, p.tf), C);
        added = registerRoot(seedFile, C, sD0, sprintf('discovery probe at column %d, round %d', jc, st.roundsDone + 1));
        nReg = nReg + added;
        if isempty(best) || C.tfDays < best.tfDays, best = C;  bestAdded = added; end
    end
    if isempty(best), continue, end

    % ---- the fastest one: anchor it, or say why not ----------------------
    [famB, ~] = F.attach(S.sA(jc), best.tfDays, best.z);         % by flight time AND costates
    [promoteIt, why] = promotionVerdict(bestAdded, best.tfDays, tf(jc), acceptDays, famB);
    if promoteIt
        st = promote(st, stateF, best, jc, S, arcDir, tag, anchorDir, lg, 'discovery');
        nNew = nNew + 1;
    else
        lg('  column %d: the %.3f d root is registered, not anchored -- %s', jc, best.tfDays, why);
    end
end
end

function common = commonOpts(in, anchors, arcDir, arcPattern)
% COMMONOPTS  What every run_costate_library call of this campaign shares:
% the two phase lists, the engine, the orbits, where this campaign's arcs
% live -- and the OPERATING POINT of the sheet: the first anchor's file and
% the arrival phase it was certified at. (The departure phase of the spine
% is sD(1), which travels in the list.) Without the anchor the sheet was set
% up from the shipped 70 mN anchor at 0.0754 whatever the campaign was.
% INPUTS: in (userInputs); anchors {name, file, sA, label; ...}; arcDir;
% arcPattern.  OUTPUTS: common struct.
common = struct('sD', in.sD, 'sA', in.sA, ...
                'thrustN', in.engine.thrustN, 'ispS', in.engine.ispS, 'm0kg', in.engine.m0kg, ...
                'tauDRO', in.orbits.tauDRO, 'NpTulip', in.orbits.NpTulip, 'pmTulip', in.orbits.pmTulip, ...
                'arcPattern', arcPattern, 'arcDir', arcDir, ...
                'anchorMat', anchors{1, 2}, 'anchorSA', anchors{1, 3}, ...
                'librarySeeds', in.librarySeeds, 'ribWallSec', in.rib.wallSec);
end

function so = physicsOpts(engine, orbits, sD0)
% PHYSICSOPTS  A setup request for the CLOSURES ONLY (endpoint states and
% propulsion constants) at the spine's departure phase -- no anchor is
% loaded or polished.  INPUTS: engine; orbits; sD0.  OUTPUTS: so struct for
% arclength_arrival('setup', so).
so = struct('thrustN', engine.thrustN, 'ispS', engine.ispS, 'm0kg', engine.m0kg, ...
            'tauDRO', orbits.tauDRO, 'NpTulip', orbits.NpTulip, 'pmTulip', orbits.pmTulip, ...
            'sD', sD0, 'physicsOnly', true);
end

function sec = roundDeadlineSec(nD, nA, nWorkers, userHours)
% ROUNDDEADLINESEC  How long to wait for one round's ribs + finalizer before
% giving up. It was a flat 48 h, sized for 24 x 24; a 48 x 48 round is about
% 55 h of honest work and the driver would have abandoned it. The wait is
% now THREE TIMES an estimate from the rates measured on the 2026-09-18/19
% rebuild (24 x 24, 4 workers): 235 worker-seconds per rib point (552 points
% in 9 h) and 28 s per entry for package + audit + sweep (519 entries in 4 h)
% -- and never less than 48 h. An explicit .roundDeadlineHours wins.
% INPUTS: nD, nA (grid); nWorkers; userHours ([] = derive).  OUTPUTS: sec.
if ~isempty(userHours), sec = userHours*3600;  return, end
ribSec   = nA*max(nD - 1, 0)*235/max(nWorkers, 1);
finalSec = nD*nA*28;
sec = max(48*3600, 3*(ribSec + finalSec));
end

function [ended, why] = campaignEnded(st, maxRounds)
% CAMPAIGNENDED  Is there nothing left for this call to do? A campaign at
% its fixed point ('done') has ended. One that stopped on its round budget
% ('budget') has ended only while .maxRounds has not been raised past the
% rounds already run -- its own stop message tells the user to raise it.
% INPUTS: st (.status .roundsDone .reason); maxRounds.  OUTPUTS: ended
% (logical); why (char, for the log).
ended = false;  why = '';
if strcmp(st.status, 'done')
    ended = true;
    why = sprintf('this campaign reached its fixed point (%s). Nothing to do; use another outDir to rebuild.', st.reason);
elseif strcmp(st.status, 'budget') && maxRounds <= st.roundsDone
    ended = true;
    why = sprintf('this campaign stopped on its budget after %d round(s) (%s). Raise .maxRounds above %d to continue it.', ...
                  st.roundsDone, st.reason, st.roundsDone);
end
end

function assertPackaged(o, roundName, what)
% ASSERTPACKAGED  A round's catalog counts only if the call PACKAGED it and
% every stage that ran behind it passed. 'packaged' alone says a catalog was
% written; an audit or sweep blocker used to ride along unnoticed.
% INPUTS: o (run_costate_library output: .state .blockers .stages);
% roundName; what (char, for the message).  OUTPUTS: none (throws).
assert(strcmp(o.state, 'packaged'), 'run_phase_torus:package', '%s: %s did not package (state %s): %s', ...
       roundName, what, o.state, strjoin(o.blockers, ' | '));
names = fieldnames(o.stages);
failed = names(structfun(@(v) isequal(v, false), o.stages));
assert(isempty(failed), 'run_phase_torus:stage', '%s: %s wrote a catalog but %s did not pass: %s', ...
       roundName, what, strjoin(failed, ', '), strjoin(o.blockers, ' | '));
end

function st = promote(st, stateF, C, jc, S, arcDir, tag, anchorDir, lg, how)
% PROMOTE  A certified root becomes an anchor with a unique, immutable name;
% the state is saved at once.  INPUTS: as named.  OUTPUTS: st.
base = sprintf('d%02d', jc);
nSame = nnz(startsWith(st.anchors(:, 1), base));
name = sprintf('%s_%d', base, nSame + 1);
anc = fullfile(anchorDir, sprintf('mintime_%s_anchor_%s.mat', tag, name));
assert(~isfile(anc), 'run_phase_torus: anchor file %s exists already', anc);
Tnd = S.B.Tnd;  cnd = S.B.cnd;
best = struct('z', C.z(:), 'it', struct('Y', C.Y), 'sA', S.sA(jc), 'sD', C.sD, 'tfDays', C.tfDays, ...
              'origin', sprintf('%s at column %d, %s', how, jc, char(datetime('now'))));
save(anc, 'best', 'Tnd', 'cnd');
st.anchors(end+1, :) = {name, anc, S.sA(jc), name};
saveState(stateF, st);
lg('  NEW ANCHOR %s: %.3f d at sA %.4f (%s; spine %s); arcs %s', name, C.tfDays, S.sA(jc), how, ...
   tern(isfinite(S.TF(jc)), sprintf('%.3f d', S.TF(jc)), 'none'), fileName(arcFile(arcDir, tag, name, -1)));
end

function added = registerRoot(seedFile, C, sD0, src)
% REGISTERROOT  Append a certified root to the campaign's registry unless
% the SAME root (same sA, same costates: sameRoot) is there; atomic write.
% INPUTS: seedFile; C (certify_root output); sD0; src.  OUTPUTS: added
% (logical).
if isfile(seedFile), L = load(seedFile);  direct = L.direct;
else, direct = struct('sD', {}, 'sA', {}, 'tfDays', {}, 'z', {}, 'src', {});
end
% A ROOT IS KNOWN BY ITS COSTATES. "The same arrival phase and a flight time
% within 1e-3 d" (86 s) merged distinct roots near a family crossing, which is
% where they matter; the same phase and the same costates is the same root.
assert(numel(C.z) >= 8 && all(isfinite(C.z(1:8))) && any(C.z(1:7) ~= 0), 'run_phase_torus:registry', ...
       'a root with non-finite or all-zero costates cannot be registered (sA %.4f, %s)', C.sA, src);
circ = @(x) abs(mod(x + 0.5, 1) - 0.5);
samePhase = circ([direct.sA] - C.sA) < 1e-8 & circ([direct.sD] - sD0) < 1e-8;      % BOTH phases: a record is a cell's root
dup = any(samePhase & arrayfun(@(r) sameRoot(r.z, C.z), direct));
added = ~dup;
if added
    direct(end+1) = struct('sD', sD0, 'sA', C.sA, 'tfDays', C.tfDays, 'z', C.z(:), 'src', src);
    tmp = [seedFile '.part'];  save(tmp, 'direct');  movefile(tmp, seedFile, 'f');
end
end

function tf = sameRoot(zA, zB)
% SAMEROOT  Are two solution vectors the same root? Both must be finite with
% non-zero costates; the seven initial costates must agree to 1e-6 of the
% larger norm (SYMMETRIC), and the flight times to 1e-6 relative. A re-polish
% moves a root by ~1e-9; this is a registry rule, not a proof of identity --
% near a fold distinct roots can be closer than any fixed tolerance, and an
% ill-conditioned polish can move one root further (review 2026-09-19).
% Normal chart, so there is no scale to quotient.
% INPUTS: zA, zB [8 x 1] (costates 1:7, t_f 8).  OUTPUTS: tf.
tf = false;
if numel(zA) < 8 || numel(zB) < 8 || ~all(isfinite(zA(1:8))) || ~all(isfinite(zB(1:8))), return, end
a = zA(1:7);  b = zB(1:7);
na = sqrt(sum(a(:).^2));  nb = sqrt(sum(b(:).^2));
if na == 0 || nb == 0, return, end
tf = sqrt(sum((a(:) - b(:)).^2)) <= 1e-6*max(na, nb) && abs(zA(8) - zB(8)) <= 1e-6*max(abs(zA(8)), abs(zB(8)));
end

function ok = ribMatchesSpine(ribFile, S, jc)
% RIBMATCHESSPINE  Does a rib file's first point continue this sheet's spine
% root at column jc (t_f within 0.5 d of the spine)?  INPUTS: ribFile; S; jc.
% OUTPUTS: ok.
ok = false;
try
    L = load(ribFile);  r = L.R(1);
    if isempty(r.pts) || ~isfinite(S.TF(jc)), return, end
    ok = abs(r.pts(1).tfDays - S.TF(jc)) < 0.5;
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

function x = roundExtras(outDir, k)
% ROUNDEXTRAS  The rib files round k offers the packager beside its own:
% every rib file of rounds 1..k-1 (and their direct-cell files), AND round
% k's own direct-cell file when an earlier attempt of the round left one.
% INPUTS: outDir; k.  OUTPUTS: x (cellstr).
x = {};
for q = 1:k-1
    V = fullfile(outDir, sprintf('round_%02d', q));
    f = [dir(fullfile(V, 'fine_rib_col*.mat')); dir(fullfile(V, 'fine_rib_direct_holes.mat'))];
    x = [x, fullfile(V, {f.name})];
end
holes = fullfile(outDir, sprintf('round_%02d', k), 'fine_rib_direct_holes.mat');
if isfile(holes), x{end+1} = holes; end
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

function H = localHandles(fh)
% LOCALHANDLES  This file's local functions as a struct of handles keyed by
% name -- the TEST SEAM: a test calls the real helper, not a copy of it.
% INPUTS: fh (cell of handles, from localfunctions).  OUTPUTS: H struct.
H = struct();
for k = 1:numel(fh), H.(func2str(fh{k})) = fh{k}; end
end
