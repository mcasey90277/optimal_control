function J = run_regime_map(cfg)
%% Purpose:
%
%   THE REGIME MAP (goal set by Mike 2026-09-07): run the eps, huber and
%   huberc continuation walks on the SAME (cell, gamma) records under ONE
%   identical budget, so that "this family arrives and that one walls" is a
%   measurement rather than an artefact of three differently-tuned campaigns.
%
%   The 2026-09-02 grid (eps vs huber) and the 2026-09-07 high-gamma race
%   (three families) used different schedules, tolerances and delta rules;
%   FINDINGS 23-28 therefore compare arms that were never given equal
%   budgets. This driver fixes every knob except the one thing that must
%   differ by construction:
%
%     * IDENTICAL: p schedule, tolR, wallSec, maxBisect, maxGaps, HdriftTol.
%     * FAMILY-CORRECT SEED: eps walks from the energy costates (its p = 1
%       minimiser is s* = Q/2 = the energy law); huber and huberc minimise
%       at s* = Q when p = 1, so lambda/2 is THEIR warm seed (Astra review
%       2026-09-05, FINDINGS 22). Handicapping them with seedScale = 1 would
%       measure the seed, not the family.
%     * FAMILY-CORRECT DELTA RULE: huberc walks p at FIXED delta = 0.03 and
%       then sharpens delta (FINDINGS 26: the fixed-delta p-walk is the
%       failure-free recipe; delta = p inherits eps's steep floor ramp).
%
%   One output file per (record, family) job, so the campaign is resumable
%   and can be partitioned across concurrent MATLAB processes.
%
%% Inputs:
%
%  cfg                      struct (optional)
%   .part                   [k n]                   run jobs k:n:end (for n
%                                                   concurrent processes)
%   .pick                   [1 x m]                 explicit job indices from
%                                                   a dryRun listing (wins
%                                                   over .part)
%   .dryRun                 logical                 list the jobs, run none
%   .outDir                 char                    default direct/results/regime
%   .logFile                char                    ''
%   .records                'all'|'pilot'|'hg'      which seed set [all]
%   .families               cellstr                 default all three
%   .redo                   logical                 recompute existing [false]
%
%% Outputs:
%
%  J                        struct array            the job list with .file,
%                                                   .done, .wallMin
%
%% Revision History:
%  M. Casey                                                   (c) 09/07/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, cfg = struct(); end
part     = getf(cfg, 'part', [1 1]);
recSel   = getf(cfg, 'records', 'all');
fams     = getf(cfg, 'families', {'eps', 'huber', 'huberc'});
redo     = getf(cfg, 'redo', false);
dryRun   = getf(cfg, 'dryRun', false);
logFile  = getf(cfg, 'logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'direct', 'results');
outDir = getf(cfg, 'outDir', fullfile(resDir, 'regime'));
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% THE ONE BUDGET ---------------------------------------------------------
B.sched      = [1 0.7 0.5 0.35 0.25 0.18 0.12 0.08 0.05 0.03 ...
                0.02 0.012 0.008 0.005 0.003 0.002 0.001];
B.wallSec    = 300;
B.tolR       = 3e-10;
B.maxBisect  = 3;
B.maxGaps    = 2;
B.HdriftTol  = 1e-6;
B.deltaFixed = 0.03;                                   % huberc p-walk
B.deltaSched = [0.02 0.012 0.008 0.005 0.003];         % huberc stage 2

%% Job list ---------------------------------------------------------------
J = struct([]);
srcs = {};
if any(strcmp(recSel, {'all', 'pilot'}))
    srcs{end+1} = struct('tag', 'grid', 'mat', fullfile(resDir, 'minenergy_pilot.mat'));
end
if any(strcmp(recSel, {'all', 'hg'}))
    srcs{end+1} = struct('tag', 'hg', 'mat', fullfile(resDir, 'minenergy_highgamma.mat'));
end
for ks = 1:numel(srcs)
    L = load(srcs{ks}.mat);
    for kr = 1:numel(L.R)
        r = L.R(kr);
        if ~usableRecord(r), continue, end
        for kf = 1:numel(fams)
            f = fams{kf};
            nm = sprintf('regime_%s_c%d%d_g%03.0f_%s.mat', srcs{ks}.tag, ...
                         r.iD, r.iA, round(r.gam*100), f);
            j = struct('src', srcs{ks}.tag, 'pilotMat', srcs{ks}.mat, ...
                       'cell', [r.iD r.iA], 'gamma', r.gam, 'family', f, ...
                       'tfDays', r.tf*4.3481, 'file', fullfile(outDir, nm), ...
                       'done', false, 'wallMin', NaN);
            if isempty(J), J = j; else, J(end+1) = j; end %#ok<AGROW>
        end
    end
end
pick = getf(cfg, 'pick', []);
if ~isempty(pick), J = J(pick); else, J = J(part(1):part(2):numel(J)); end
lg('[regime] %d jobs in this partition (%d of %d)', numel(J), part(1), part(2));
if dryRun
    for k = 1:numel(J)
        lg('  job %2d: %-34s tf %5.1f d -> %s', k, nameOf(J(k)), J(k).tfDays, ...
           tern(exist(J(k).file, 'file') > 0, 'EXISTS', 'todo'));
    end
    return
end

%% Run --------------------------------------------------------------------
for k = 1:numel(J)
    if exist(J(k).file, 'file') && ~redo
        J(k).done = true;
        lg('[%2d/%2d] SKIP (exists) %s', k, numel(J), nameOf(J(k)));
        continue
    end
    c = struct('cell', J(k).cell, 'gamma', J(k).gamma, 'families', {{J(k).family}}, ...
               'sched', B.sched, 'wallSec', B.wallSec, 'tolR', B.tolR, ...
               'maxBisect', B.maxBisect, 'maxGaps', B.maxGaps, ...
               'HdriftTol', B.HdriftTol, 'pilotMat', J(k).pilotMat, ...
               'outMat', J(k).file, 'logFile', logFile);
    if strcmp(J(k).family, 'eps')
        c.seedScale = 1;
    else
        c.seedScale = 0.5;                        % huber/huberc warm seed
    end
    if strcmp(J(k).family, 'huberc')
        c.delta      = B.deltaFixed;
        c.deltaSched = B.deltaSched;
    end
    t0 = tic;
    lg('[%2d/%2d] RUN %s (tf %.1f d)', k, numel(J), nameOf(J(k)), J(k).tfDays);
    try
        o = run_minfuel_race(c);
        A = o.arms.(J(k).family);
        J(k).done = true;  J(k).wallMin = toc(t0)/60;
        lg('[%2d/%2d] DONE %s: p_floor=%.4g mf=%.6f fails=%d bisects=%d retired=%d (%.1f min)', ...
           k, numel(J), nameOf(J(k)), lastOr(A.p, NaN), lastOr(A.mf, NaN), ...
           A.nFail, A.nBisect, A.retired, J(k).wallMin);
    catch ME
        lg('[%2d/%2d] ERROR %s: %s', k, numel(J), nameOf(J(k)), ME.message);
    end
end
lg('[regime] partition %d/%d complete: %d/%d jobs done', part(1), part(2), nnz([J.done]), numel(J));
end

% ------------------------------------------------------------------------
function ok = usableRecord(r)
% USABLERECORD  A seed record is usable when its multiple-shooting solution
% converged (the direct-fail rows carry no ms field content).
% INPUTS: r (pilot-format record).  OUTPUTS: ok [logical].
ok = isstruct(r.ms) && isfield(r.ms, 'z') && ~isempty(r.ms.z) && ...
     isfield(r.ms, 'converged') && r.ms.converged;
end

function s = nameOf(j)
% NAMEOF  Compact job label.  INPUTS: j.  OUTPUTS: s [char].
s = sprintf('%s (%d,%d)@%.4g %s', j.src, j.cell(1), j.cell(2), j.gamma, j.family);
end

function v = lastOr(x, d)
% LASTOR  Last element or a default.  INPUTS: x; d.  OUTPUTS: v.
if isempty(x), v = d; else, v = x(end); end
end

function s = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end

function v = getf(s, f, d)
% GETF  Field with default.  INPUTS: s; f; d.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function logmsg(f, s)
% LOGMSG  Append to a log file or stdout.  INPUTS: f; s.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
end
