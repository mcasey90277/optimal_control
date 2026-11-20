function plan = run_dro_elfo_torus(sel, opts)
%% Purpose:
%
%   FRONT DOOR for the DRO -> ELFO (and matched DRO -> tulip) phase-torus
%   fleet (DRO to ELFO plan, Task 5; spec docs/superpowers/specs/2026-11-19-
%   dro-elfo-design.md section E): the shared three-gate ladder engine
%   costate_common/thrust_ladder_library (direct solve -> ms_tfmin
%   refinement -> pumpkyn tfMin acceptance) walked down the thrust ladder
%   over the (DRO departure phase x arrival phase) torus of each selected
%   sheet. A copy of GTO_DRO/catalog/run_gto_dro_catalog.m with the
%   parameter block changed; kept from it unchanged: the resume from the
%   sheet on disk, the engine's attempt counters written before solving,
%   the maxCells / batchSec clean exit between cells, the two-point
%   top-rung multistart (pass A tf0 seed, pass B a second seed on the
%   cells pass A attempted and could not clear), the ms hard-cap fence opened
%   once per call with the requireFence refusal, the on-disk rung guard,
%   FAIL LOUD on an engine error, the open-cell sweep from disk and the
%   completion marker in combined_progress.txt that the shell driver
%   (run_dro_elfo_batched.sh) stops on.
%
%   PASS B SEED (differs from the GTO model, whose pass B falls back to the
%   engine's cold default tf0 = 4.0 ND): with thrust pinned at 15 N the
%   150 kg vehicle reaches the 5 percent mass floor at 0.416 ND, so a
%   4.0 ND seed is outside the feasible flight-time range and the retry is
%   nearly hopeless. Pass B seeds tf0B = 0.15 ND instead, the other side
%   of pass A's 0.3 ND around the 15 N reference flight (0.090 ND on the
%   stage 0 cell): a seed below and a seed above the expected t_f.
%
%   COMPLETION MARKER: written only when every selected sheet has no open
%   cell, and it names the sheet tags it covers --
%     SELECTED SHEETS COMPLETE[ [run <runId>]]: <tag> <tag> ...
%   or, only for the full default 12-sheet list, CATALOG ALL SHEETS
%   COMPLETE[ [run <runId>]]: <tags>. Every non-planOnly call first logs
%     SHEETS SELECTED[ [run <runId>]]: <tags>
%   so the shell driver (which passes its own runId) accepts only a marker
%   of ITS run whose tag list equals ITS selection; an ELFO fleet and a
%   tulip run sharing a folder cannot satisfy each other's marker.
%
%   SHEET IDENTITY (spec section E; nothing is overloaded): the engine's
%   tauDRO is the DRO period droTau, the arrival is arrFamily / arrParams,
%   the ELFO plane orientation lives in arrParams.raanDeg. The engine
%   rebuilds meta and rewrites the whole sheet file (save without -append)
%   at every save, so an extra meta field would not survive; the key is
%   therefore a separate top-level variable
%     sheet_key = struct('target', 'elfo'|'tulip', 'depSize', droTau,
%                        'raanDeg', raanDeg or NaN)
%   appended to the sheet (save -append) after each engine call when it is
%   not already there. Before the engine touches an existing sheet this
%   function checks, in order: its rungs (rungMismatch), its grid size and
%   phase origins (gridMismatch: the engine's resume test is a bare size
%   compare and would silently start over and overwrite), and its identity
%   (sheetKey: the stored sheet_key AND the key derived from meta, plus the
%   ELFO elements, must match the requested sheet; a sheet whose identity
%   cannot be read is refused too). Never overwritten on a mismatch.
%
%% Inputs:
%
%  sel                      [1 x k] | char | string Indices into
%                                                   opts.sheets; a string
%                                                   such as "[1 3:5]" or
%                                                   "2,12" (digits, ':',
%                                                   ',', spaces, brackets)
%                                                   is accepted for the
%                                                   shell driver; [] or ''
%                                                   = every sheet. Any
%                                                   index outside
%                                                   1..numel(sheets) is
%                                                   refused (:sel)
%  opts                     struct (optional)       Unknown fields are
%                                                   refused (:opts)
%   .sheets                 [s x 1] struct          .target 'elfo'|'tulip',
%                                                   .droTau (DRO period,
%                                                   ND), .raanDeg (integer
%                                                   deg, wrapped to
%                                                   [0,360); NaN for tulip)
%                                                   [12 ELFO sheets, droTau
%                                                   1, raanDeg 0:30:330]
%   .oe                     [1 x 4]                 ELFO elements [km - deg
%                                                   deg] [im_elfo_optimum]
%   .nD, .nA                double                  Torus resolution [12 12]
%   .sD0, .sA0              double                  Phase origins [0 0]
%   .rungs                  [1 x r]                 Thrust ladder, N, high
%                                                   to low, r >= 2,
%                                                   strictly decreasing
%                                                   (:rungs; one rung is
%                                                   refused, see the code)
%                                                   [15 12 10 7 5 3 2 1.5 1
%                                                   0.75 0.5]
%   .outDir                 char                    Sheets, logs, completion
%                                                   marker
%                                                   [DRO_ELFO/results/
%                                                   catalog]
%   .maxCells               double                  Fresh cells per sheet
%                                                   per pass [inf]
%   .batchSec               double                  Clean-exit wall budget
%                                                   (s) shared by the
%                                                   selected sheets [inf]
%   .planOnly               logical                 Return the plan, no
%                                                   solve, no disk [false]
%   .tf0                    double                  Pass A top-rung seed
%                                                   (ND) [0.3]
%   .tf0B                   double                  Pass B top-rung seed
%                                                   (ND), see PASS B SEED
%                                                   [0.15]
%   .runId                  char                    Run label carried in
%                                                   the SELECTED and
%                                                   COMPLETE log lines
%                                                   (the shell driver's)
%                                                   ['']
%   .maxCpuSec              double                  IPOPT CPU cap per rung
%                                                   (s) [600]
%   .maxIter                double                  IPOPT iteration cap
%                                                   [6000]
%   .msHardCapSec           double                  Hard cap (s) per
%                                                   ms_tfmin attempt; 0 =
%                                                   off [300]
%   .requireFence           logical                 Refuse (:noFence) when
%                                                   the cap is on and no
%                                                   pool can be had [true]
%   .floorKm                double                  Lunar altitude floor
%                                                   (km) [500]
%   .floorEarthKm           double                  Earth floor (km); [] =
%                                                   none [[]]
%   .directGateKm           double                  First flown gate on
%                                                   the direct control (km)
%                                                   [inf]
%   .gateKm                 double                  Refined flown gate (km)
%                                                   [100]
%   .keepDirect             logical                 Engine direct-solution
%                                                   sidecar <tag>_direct.mat
%                                                   [true]
%   .N                      double                  Collocation intervals
%                                                   [400]
%  (.poolFcn                function handle         TEST SEAM, not for
%                                                   campaigns: opens the
%                                                   fence pool, same
%                                                   contract as capped_pool
%                                                   [@capped_pool])
%
%% Outputs:
%
%  plan                     [k x 1] struct          .tag .target .droTau
%                                                   .raanDeg .sheet_key
%                                                   .fence .outMat .logFile
%                                                   .optsA (engine options,
%                                                   pass A) .optsB (pass B,
%                                                   tf0 = tf0B, .cells empty
%                                                   in the plan); batchSec
%                                                   and optsB.cells are set
%                                                   at run time. Files under
%                                                   outDir: <tag>.mat (with
%                                                   sheet_key),
%                                                   <tag>_progress.txt,
%                                                   <tag>_direct.mat,
%                                                   combined_progress.txt
%
%% Revision History:
%  M. Casey                                                   (c) 11/19/2026
%  M. Casey  tf0B pass B seed, tagged marker, corruptSheet    (c) 11/19/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, sel = []; end
if nargin < 2 || isempty(opts), opts = struct(); end
% A misspelt option would otherwise be ignored and a fleet run for hours
% on the default it meant to change: unknown fields are refused.
known = {'sheets', 'oe', 'nD', 'nA', 'sD0', 'sA0', 'rungs', 'outDir', 'maxCells', ...
         'batchSec', 'planOnly', 'tf0', 'maxCpuSec', 'maxIter', 'msHardCapSec', ...
         'requireFence', 'floorKm', 'floorEarthKm', 'directGateKm', 'gateKm', ...
         'keepDirect', 'N', 'poolFcn', 'tf0B', 'runId'};
bad = setdiff(fieldnames(opts), known);
if ~isempty(bad)
    error('run_dro_elfo_torus:opts', 'unknown option(s): %s', strjoin(bad, ', '));
end
planOnly = fdef(opts, 'planOnly', false);

%% ===================== ADJUSTABLE PARAMETERS ============================

%% Sheets (target x DRO period x ELFO plane orientation):
       sheets = fdef(opts, 'sheets', ...
                     struct('target', 'elfo', 'droTau', 1, 'raanDeg', num2cell(0:30:330)));

%% Phasing grid per sheet (departure x arrival; uniform in time fraction):
           nD = fdef(opts, 'nD', 12);
           nA = fdef(opts, 'nA', 12);
          sD0 = fdef(opts, 'sD0', 0);
          sA0 = fdef(opts, 'sA0', 0);     % the shipped tulip catalog's arrival grid starts at 0

%% Thrust ladder (N, high to low); must match any sheet already on disk:
        rungs = fdef(opts, 'rungs', [15 12 10 7 5 3 2 1.5 1 0.75 0.5]);
         ispS = 1710;                % specific impulse (s)
         m0kg = 150;                 % initial mass (kg)

%% Solver / gates (the stage 0 recipe, 18 of 18 accepted):
            N = fdef(opts, 'N', 400);              % collocation intervals
      floorKm = fdef(opts, 'floorKm', 500);        % lunar altitude floor (km)
 floorEarthKm = fdef(opts, 'floorEarthKm', []);    % Earth floor (km); [] = none
      maxIter = fdef(opts, 'maxIter', 6000);       % IPOPT iteration cap
    maxCpuSec = fdef(opts, 'maxCpuSec', 600);      % IPOPT CPU cap per rung (s)
       gateKm = fdef(opts, 'gateKm', 100);         % refined flown-arrival gate (km)
 directGateKm = fdef(opts, 'directGateKm', inf);   % first flown gate on the direct control; inf = none
       accTol = 1e-6;                              % tfMin acceptance
 requireFence = fdef(opts, 'requireFence', true);  % cap on but no pool: refuse (true) or run UNFENCED
 msHardCapSec = fdef(opts, 'msHardCapSec', 300);   % hard cap (s) per ms_tfmin attempt; 0 = off
      thrLock = true;                              % pin thr == 1 (classical min-time form)
   keepDirect = fdef(opts, 'keepDirect', true);    % keep each rung's direct solution (sidecar)
     tf0Fleet = fdef(opts, 'tf0', 0.3);            % pass A top-rung seed (ND)
      tf0Mop  = fdef(opts, 'tf0B', 0.15);          % pass B top-rung seed (ND), header PASS B SEED
        runId = char(fdef(opts, 'runId', ''));     % run label in the SELECTED / COMPLETE lines
     maxCells = fdef(opts, 'maxCells', inf);       % fresh cells per sheet per pass
     batchSec = fdef(opts, 'batchSec', inf);       % clean-exit wall budget (s)
      poolFcn = fdef(opts, 'poolFcn', @capped_pool); % TEST SEAM (header), not for campaigns

%% Arrival shapes:
   tulipArr   = struct('Np', 7, 'pm', -1);

%% ======================= END ADJUSTABLE PARAMETERS ======================

% Rungs: at least two, strictly decreasing, positive. ONE rung is refused:
% the engine's resume test is isequal(size(Q.OK), [nD nA nR]) and size()
% drops a trailing singleton, so a one-rung sheet is never recognised and
% is re-solved (and its attempt counter reset) on every call -- measured
% 2026-11-19, a 1 x 1 x 1 sheet re-solved on the second call; an
% unattended driver would loop on it.
if ~(isnumeric(rungs) && isvector(rungs) && numel(rungs) >= 2 && all(isfinite(rungs)) && ...
        all(rungs > 0) && all(diff(rungs(:)) < 0))
    error('run_dro_elfo_torus:rungs', ...
          ['rungs = %s: need at least two positive rungs, strictly decreasing (one rung ' ...
           'cannot be resumed by the engine)'], mat2str(rungs));
end
rungs = rungs(:)';

% Assert the departure-phase grid the engine will build internally has no
% duplicate phase (no sD=1 aliasing sD=0):
sDCheck = mod(sD0 + (0:nD-1)/nD, 1);
assert(numel(unique(sDCheck)) == nD, 'run_dro_elfo_torus:sDgrid', ...
    'departure-phase grid has a duplicate phase');

here = fileparts(mfilename('fullpath'));                 % DRO_ELFO/catalog
addpath(fileparts(here));
setup_dro_elfo_paths();
if isfield(opts, 'oe')
    oe = opts.oe;
else
    [~, ~, ~, oe] = im_elfo_optimum();                   % the IM constellation's ELFO
end
outDir = fdef(opts, 'outDir', fullfile(fileparts(here), 'results', 'catalog'));

% --- the sheet list: validated whole, before anything runs ---------------
sheets  = sheets(:);
nSheets = numel(sheets);
assert(nSheets >= 1, 'run_dro_elfo_torus:sheets', 'opts.sheets is empty');
resolved = cell(nSheets, 1);
for ks = 1:nSheets
    resolved{ks} = sheet_resolve(sheets(ks));
end
sheets  = vertcat(resolved{:});
tagsAll = {sheets.tag};
if numel(unique(tagsAll)) < nSheets
    error('run_dro_elfo_torus:dupSheet', ...
          'two sheets resolve to one tag (%s): the second would overwrite the first', ...
          strjoin(tagsAll, ', '));
end
sel = parse_sel(sel, nSheets);

if ~planOnly && ~isfolder(outDir), mkdir(outDir); end
logComb = fullfile(outDir, 'combined_progress.txt');
runSfx  = '';
if ~isempty(runId), runSfx = sprintf(' [run %s]', runId); end
selTags = strjoin(tagsAll(sel), ' ');
if ~planOnly
    log_line(logComb, sprintf('SHEETS SELECTED%s: %s', runSfx, selTags));
end

% THE FENCE (msHardCapSec > 0): one pool for the whole call, opened AFTER
% the paths are set so the workers inherit pumpkyn, costate_common, oclib
% and CasADi; every engine call finds it through current_pool and keeps
% it. A fleet must not run unfenced by accident: with requireFence
% (default) an unavailable pool is a refusal, not a fallback.
fenceLine = 'ms fence OFF (msHardCapSec = 0)';
if planOnly
    fenceLine = 'ms fence not opened (planOnly)';
elseif msHardCapSec > 0
    % WARNING FLOOD: a stuck flight emits one ode45 tolerance warning per
    % restart. Off for this call when the cap is on, restored on exit;
    % fenced workers silence it themselves.
    odeWarn = warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
    restoreOdeWarn = onCleanup(@() warning(odeWarn));
    hadPool = ~isempty(current_pool());
    fenceWhy = '';
    lastwarn('');
    try
        fencePool = poolFcn(1);
    catch ME
        fencePool = [];  fenceWhy = ME.message;
    end
    if isempty(fencePool)
        if isempty(fenceWhy), [fenceWhy, ~] = lastwarn; end
        fenceLine = sprintf('ms fence UNAVAILABLE -- running UNFENCED: %s', fenceWhy);
    else
        fenceLine = sprintf('ms fence ACTIVE (pool of %d)', fencePool.NumWorkers);
        if ~hadPool, closeFence = onCleanup(@() delete_pool(fencePool)); end
    end
    log_line(logComb, sprintf('=== %s  %s', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), fenceLine));
    fprintf('%s\n', fenceLine);
    if isempty(fencePool) && requireFence
        error('run_dro_elfo_torus:noFence', ...
              ['msHardCapSec = %g requested but no parallel pool could be opened (%s); refusing to run ' ...
               'unfenced (pass opts.requireFence = false to allow it)'], msHardCapSec, fenceWhy);
    end
end

plan  = struct('tag',{},'target',{},'droTau',{},'raanDeg',{},'sheet_key',{},'fence',{}, ...
               'outMat',{},'logFile',{},'optsA',{},'optsB',{});
planC = cell(numel(sel), 1);
kp    = 0;
tAll = tic;
allDone = true;
for ks = sel
    if ~planOnly && toc(tAll) > batchSec, allDone = false; break, end
    s = sheets(ks);
    outMat  = fullfile(outDir, [s.tag '.mat']);
    logFile = fullfile(outDir, [s.tag '_progress.txt']);
    if strcmp(s.target, 'elfo')
        arrParams = struct('oe', oe, 'raanDeg', s.raanDeg);
    else
        arrParams = tulipArr;
    end
    sheet_key = struct('target', s.target, 'depSize', s.droTau, 'raanDeg', s.raanDeg);
    optsCommon = struct('rungs',rungs, 'ispS',ispS, 'm0kg',m0kg, ...
        'N',N, 'floorKm',floorKm, 'floorEarthKm',floorEarthKm, ...
        'maxIter',maxIter, 'maxCpuSec',maxCpuSec, ...
        'gateKm',gateKm, 'directGateKm',directGateKm, 'msHardCapSec',msHardCapSec, ...
        'nD',nD, 'nA',nA, 'sD0',sD0, 'sA0',sA0, ...
        'depFamily','dro', 'depParams',struct('tau',s.droTau), ...
        'arrFamily',s.target, 'arrParams',arrParams, ...
        'tauDRO',s.droTau, ...             % the DRO period (spec section E: nothing overloaded)
        'accTol',accTol, 'maxCells',maxCells, ...
        'thrLock',thrLock, 'keepDirect',keepDirect, 'logFile',logFile);
    optsA = optsCommon;  optsA.tf0 = tf0Fleet;
    optsB = optsCommon;  optsB.tf0 = tf0Mop;  optsB.cells = zeros(0,2);
    kp = kp + 1;
    planC{kp} = struct('tag', s.tag, 'target', s.target, 'droTau', s.droTau, ...
        'raanDeg', s.raanDeg, 'sheet_key', sheet_key, 'fence', fenceLine, ...
        'outMat', outMat, 'logFile', logFile, 'optsA', optsA, 'optsB', optsB);
    if planOnly, continue, end
    fprintf('\n=== DRO -> %s sheet %d/%d: %s ===\n', upper(s.target), ks, nSheets, s.tag);
    % GUARDS before the engine ever touches an existing sheet (header):
    if isfile(outMat)
        guard_sheet(outMat, s.tag, rungs, nD, nA, sD0, sA0, sheet_key, oe);
    end
    % TWO-POINT MULTISTART, FIRST-SUCCESS-WINS (GTO_tulip/catalog/
    % run_gto_catalog.m): thrLock on both passes; only the top-rung seed
    % varies. Pass A -- warm tf0 seed:
    optsA.batchSec = max(0, batchSec - toc(tAll));
    run_engine(outMat, optsA, s.tag, 'A', logComb, sheet_key);
    % Pass B -- the second seed tf0B, SCOPED to the cells pass A attempted (ATT >= 1)
    % and could not clear (no OK rung), read fresh from disk.
    if toc(tAll) <= batchSec
        mopCells = zeros(0,2);
        if isfile(outMat)
            S = load(outMat, 'ATT', 'OK');
            [mD, mA] = find(S.ATT >= 1 & ~any(S.OK,3));
            mopCells = [mD, mA];
        end
        if ~isempty(mopCells)
            optsB.cells    = mopCells;
            optsB.batchSec = max(0, batchSec - toc(tAll));
            run_engine(outMat, optsB, s.tag, 'B', logComb, sheet_key);
        end
    end
    % OPEN-CELL SWEEP from disk (never in-memory state); no sheet on disk
    % (zero cells attempted) is trivially open.
    if isfile(outMat)
        Qsweep = load(outMat, 'OK', 'ATT');
        open_ = any(~any(Qsweep.OK, 3) & Qsweep.ATT < 2, 'all');
    else
        open_ = true;
    end
    if open_, allDone = false; end
end

plan = vertcat(plan, planC{1:kp});
if planOnly, return, end

if allDone && ~isfield(opts, 'sheets') && numel(sel) == nSheets
    log_line(logComb, sprintf('CATALOG ALL SHEETS COMPLETE%s: %s', runSfx, selTags));
elseif allDone
    log_line(logComb, sprintf('SELECTED SHEETS COMPLETE%s: %s', runSfx, selTags));
end
end

% ------------------------------------------------------------------------
function s = sheet_resolve(s)
% SHEET_RESOLVE  Validate one sheet entry, wrap raanDeg, attach its tag.
% INPUTS: s struct (.target .droTau .raanDeg).  OUTPUTS: s with .tag.
if ~all(isfield(s, {'target', 'droTau', 'raanDeg'}))
    error('run_dro_elfo_torus:sheets', 'each sheet needs .target, .droTau and .raanDeg');
end
if ~(isscalar(s.droTau) && isnumeric(s.droTau) && isfinite(s.droTau) && s.droTau > 0)
    error('run_dro_elfo_torus:droTau', 'droTau must be a positive finite scalar (ND)');
end
switch char(s.target)
    case 'elfo'
        r = s.raanDeg;
        if ~(isscalar(r) && isnumeric(r) && isfinite(r) && r == round(r))
            error('run_dro_elfo_torus:raanDeg', ...
                  'ELFO raanDeg must be an integer number of degrees (tags carry %%03d), got %s', mat2str(r));
        end
        s.raanDeg = mod(r, 360);
        s.tag     = sprintf('dro_elfo_tau%g_raan%03d', s.droTau, s.raanDeg);
    case 'tulip'
        if ~(isscalar(s.raanDeg) && isnan(s.raanDeg))
            error('run_dro_elfo_torus:raanDeg', 'a tulip sheet has no plane orientation: raanDeg must be NaN');
        end
        s.tag = sprintf('dro_tulip_tau%g', s.droTau);
    otherwise
        error('run_dro_elfo_torus:target', 'target must be ''elfo'' or ''tulip'', got ''%s''', char(s.target));
end
s = struct('target', char(s.target), 'droTau', s.droTau, 'raanDeg', s.raanDeg, 'tag', s.tag);
end

% ------------------------------------------------------------------------
function sel = parse_sel(sel, nSheets)
% PARSE_SEL  Sheet selection as a row of indices, every one in range.
% INPUTS: sel numeric | char | string; nSheets double.  OUTPUTS: sel [1 x k].
% Only a literally empty or blank selection means every sheet; a
% selection that parses to an empty list ('3:1', '1:0:3', ';') is refused.
if isstring(sel), sel = char(sel); end
if isempty(sel) || (ischar(sel) && isempty(strtrim(sel)))
    sel = 1:nSheets;
    return
end
if ischar(sel)
    txt = strtrim(sel);
    if isempty(regexp(txt, '^[\[\]\d\s:,;]*$', 'once'))
        error('run_dro_elfo_torus:sel', 'selection string "%s" is not an index list', txt);
    end
    toks = regexp(strrep(strrep(txt, '[', ' '), ']', ' '), '[^\s,;]+', 'match');
    parts = cell(1, numel(toks));
    for kt = 1:numel(toks)
        p = str2double(strsplit(toks{kt}, ':'));
        if any(isnan(p)) || numel(p) > 3
            error('run_dro_elfo_torus:sel', 'selection token "%s" is not an index or range', toks{kt});
        end
        switch numel(p)
            case 1, parts{kt} = p;
            case 2, parts{kt} = p(1):p(2);
            case 3, parts{kt} = p(1):p(2):p(3);
        end
    end
    sel = [parts{:}];
end
if isempty(sel)
    error('run_dro_elfo_torus:sel', 'the selection parses to no sheet at all');
end
sel = sel(:)';
if ~isnumeric(sel) || any(sel ~= round(sel)) || any(sel < 1) || any(sel > nSheets)
    error('run_dro_elfo_torus:sel', 'sel = %s: every index must lie in 1..%d', mat2str(sel), nSheets);
end
end

% ------------------------------------------------------------------------
function guard_sheet(outMat, tag, rungs, nD, nA, sD0, sA0, key, oe)
% GUARD_SHEET  Refuse an on-disk sheet whose rungs, grid or identity differ
% from the requested sheet (the engine would otherwise re-solve from
% scratch or mix two sheets in one file).
% INPUTS: outMat char; tag char; rungs [1 x r]; nD, nA, sD0, sA0 double;
%         key struct (target, depSize, raanDeg); oe [1 x 4].
% OUTPUTS: none (throws run_dro_elfo_torus:corruptSheet | :rungMismatch |
%          :gridMismatch | :sheetKey).
try
    vars = who('-file', outMat);
    want = vars(ismember(vars, {'rungs', 'OK', 'meta', 'sheet_key'}));
    Q = struct();
    if ~isempty(want), Q = load(outMat, want{:}); end
catch ME
    error('run_dro_elfo_torus:corruptSheet', ...
          'sheet file %s cannot be read (%s) -- refusing to run; the file is left in place', ...
          outMat, ME.message);
end
if ~isfield(Q, 'rungs') || ~isequal(Q.rungs(:)', rungs(:)')
    have = '(none)';
    if isfield(Q, 'rungs'), have = num2str(Q.rungs(:)'); end
    error('run_dro_elfo_torus:rungMismatch', ...
        ['sheet %s on disk has rungs=[%s] but this call requests rungs=[%s] -- refusing to run ' ...
         '(the engine''s resume check would silently re-solve from scratch and overwrite the ' ...
         'sheet''s accumulated entries). Remedy: pass opts.rungs matching the sheet, or move ' ...
         'the sheet aside deliberately.'], tag, have, num2str(rungs));
end
okGrid = isfield(Q, 'OK') && size(Q.OK, 1) == nD && size(Q.OK, 2) == nA;
if okGrid && isfield(Q, 'meta') && all(isfield(Q.meta, {'sD0', 'sA0'}))
    okGrid = abs(Q.meta.sD0 - sD0) < 1e-12 && abs(Q.meta.sA0 - sA0) < 1e-12;
end
if ~okGrid
    error('run_dro_elfo_torus:gridMismatch', ...
          ['sheet %s on disk has another torus grid than this call (nD = %d, nA = %d, sD0 = %g, ' ...
           'sA0 = %g) -- refusing to run (the engine would start the sheet over)'], tag, nD, nA, sD0, sA0);
end
keys = {};
if isfield(Q, 'sheet_key'), keys{end+1} = Q.sheet_key; end
if isfield(Q, 'meta') && all(isfield(Q.meta, {'arrFamily', 'arrParams', 'depParams'})) && ...
        isfield(Q.meta.depParams, 'tau')
    r = NaN;
    if strcmp(Q.meta.arrFamily, 'elfo')
        r = NaN(1, 2);                                   % unreadable until shown otherwise
        if isfield(Q.meta.arrParams, 'raanDeg'), r = Q.meta.arrParams.raanDeg; end
        if ~isfield(Q.meta.arrParams, 'oe') || ~isequal(Q.meta.arrParams.oe, oe)
            error('run_dro_elfo_torus:sheetKey', ...
                  'sheet %s on disk was built for other ELFO elements than [%s] -- refusing to run', ...
                  tag, num2str(oe));
        end
    end
    keys{end+1} = struct('target', Q.meta.arrFamily, 'depSize', Q.meta.depParams.tau, 'raanDeg', r);
end
if isempty(keys)
    error('run_dro_elfo_torus:sheetKey', ...
          'sheet %s on disk carries neither sheet_key nor a readable meta -- identity unverifiable, refusing', tag);
end
for kk = 1:numel(keys)
    if ~isequaln(keys{kk}, key)
        error('run_dro_elfo_torus:sheetKey', ...
              ['sheet %s on disk belongs to (target %s, depSize %g, raanDeg %s), not the requested ' ...
               '(target %s, depSize %g, raanDeg %g) -- refusing to run, the sheet is not overwritten'], ...
              tag, char(keys{kk}.target), keys{kk}.depSize, mat2str(keys{kk}.raanDeg), ...
              key.target, key.depSize, key.raanDeg);
    end
end
end

% ------------------------------------------------------------------------
function run_engine(outMat, o, tag, passName, logComb, sheet_key)
% RUN_ENGINE  One engine pass, FAIL LOUD on an error (logged to stderr and
% the combined log, then rethrown so the shell driver aborts); afterwards
% append sheet_key to the sheet when it is not already there (the engine
% rewrites the whole file at every save).
% INPUTS: outMat char; o struct (engine options); tag char; passName char;
%         logComb char; sheet_key struct.  OUTPUTS: none.
try
    thrust_ladder_library(outMat, o);
catch ME
    failMsg = sprintf('FATAL: sheet %s pass %s (thrust_ladder_library) errored -- %s', ...
        tag, passName, ME.message);
    fprintf(2, '%s\n', failMsg);
    log_line(logComb, failMsg);
    if isfile(outMat), append_key(outMat, sheet_key); end
    rethrow(ME);
end
if isfile(outMat), append_key(outMat, sheet_key); end
end

% ------------------------------------------------------------------------
function append_key(outMat, sheet_key)
% APPEND_KEY  save -append sheet_key unless the file already holds it, so a
% call that solved nothing leaves the sheet untouched.
% INPUTS: outMat char; sheet_key struct.  OUTPUTS: none (file written).
if ~ismember('sheet_key', who('-file', outMat))
    save(outMat, 'sheet_key', '-append');
end
end

% ------------------------------------------------------------------------
function log_line(f, s)
% LOG_LINE  Append one line to file f (silently skipped if it cannot open).
% INPUTS: f char; s char.  OUTPUTS: none.
fid = fopen(f, 'a');
if fid > 0, fprintf(fid, '%s\n', s); fclose(fid); end
end

% ------------------------------------------------------------------------
function delete_pool(pool)
% DELETE_POOL  Delete the fence pool this call opened, if still valid.
% INPUTS: pool parallel.Pool.  OUTPUTS: none.
if ~isempty(pool) && isvalid(pool)
    try delete(pool); catch, end
end
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
