function H = run_huber_race_grid(cfg)
% RUN_HUBER_RACE_GRID  The Huber (PLQ) continuation family over EVERY min-fuel
% grid record, with the saltation-correct propagator -- the robustness test
% that decides whether the 2026-09-05 single-cell result (17/17 rungs, 0
% fails) is real or one lucky cell.
%
% For each (cell, gamma) record in minfuel_grid.mat: runs the huber arm of
% run_minfuel_race (same 17-rung schedule, gates and hard caps as the eps
% arms already on disk), then re-solves the deepest huber solution with the
% corrected fixed-tf conjugate test (rows 1:7). One control run repeats
% (2,5)@1.2 with the seed costates halved (huber kappa=1 minimises at
% s* = Q, so lambda/2 is its warm seed). Aggregates huber vs the stored eps
% arms: depth, fails, bisections, wall, m_f, |dmf|, coast, conj verdict.
% Saves direct/results/huber_grid.mat and prints the table.
%
% INPUTS:
%   cfg - (optional) struct: .logFile [''], .wallSec [300] (per solve),
%         .control [true] (run the lambda/2 control cell)
%
% OUTPUTS:
%   H - struct array, one per record (+ the control): .cellIdx .gamma
%       .seedScale .pHuber .pEps .mfHuber .mfEps .dMf .failsH .failsE
%       .bisH .bisE .wallH .wallE .coastH .coastE .acceptH .conjPass
%       .nCross .firstFullRank
%
% REFERENCES:
%   [1] run_minfuel_race.m (the arm; seedScale added 2026-09-05).
%   [2] DRO_tulip/FINDINGS.md section 22 (the saltation fix and the one-cell
%       result this generalises).

if nargin < 1, cfg = struct(); end
logFile = '';  wallSec = 300;  doControl = true;
if isfield(cfg, 'logFile'),  logFile  = cfg.logFile;  end
if isfield(cfg, 'wallSec'),  wallSec  = cfg.wallSec;  end
if isfield(cfg, 'control'),  doControl = cfg.control; end
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'indirect'), fullfile(fileparts(here), 'costate_common'));
resDir = fullfile(here, 'direct', 'results');
Lg = load(fullfile(resDir, 'minfuel_grid.mat'));  G = Lg.G;
Lp = load(fullfile(resDir, 'minenergy_pilot.mat'));
outMat = fullfile(resDir, 'huber_grid.mat');

jobs = arrayfun(@(g) struct('cell', g.cellIdx, 'gamma', g.gamma, 'seedScale', 1), G);
if doControl
    % lambda/2 controls: the clean cell and the cell whose p = 1 rung failed cold
    jobs(end+1) = struct('cell', [2 5], 'gamma', 1.2, 'seedScale', 0.5);
    jobs(end+1) = struct('cell', [1 2], 'gamma', 1.2, 'seedScale', 0.5);
end
lg('HUBER GRID: %d runs (%d records%s)', numel(jobs), numel(G), tern(doControl, ' + 2 lambda/2 controls', ''));

H = struct([]);
tAll = tic;
for k = 1:numel(jobs)
    jb = jobs(k);
    tag = sprintf('c%d%d_g%03d', jb.cell, round(100*jb.gamma));
    if jb.seedScale ~= 1, tag = [tag '_half']; end
    om = fullfile(resDir, ['minfuel_huber_' tag '.mat']);
    lg('HUBER GRID: (%d,%d) gamma %.2f seedScale %.2g -> %s', jb.cell, jb.gamma, jb.seedScale, om);
    c2 = struct('cell', jb.cell, 'gamma', jb.gamma, 'families', {{'huber'}}, ...
                'outMat', om, 'logFile', logFile, 'wallSec', wallSec, 'seedScale', jb.seedScale);
    try
        if isfile(om)                      % resume: a finished arm is on disk
            Lo = load(om);  out = Lo.out;
            if isfield(out.arms, 'huber') && isfield(out.arms.huber, 'wallTotal')
                lg('HUBER GRID: (%d,%d) gamma %.2f already on disk -- reusing', jb.cell, jb.gamma);
            else
                out = run_minfuel_race(c2);
            end
        else
            out = run_minfuel_race(c2);
        end
        A = out.arms.huber;
    catch ME
        lg('HUBER GRID: (%d,%d) gamma %.2f ERROR: %s', jb.cell, jb.gamma, ME.message);
        continue
    end
    % the eps arm on disk for this record:
    E = [];
    ef = fullfile(resDir, sprintf('minfuel_eps_c%d%d_g%03d.mat', jb.cell, round(100*jb.gamma)));
    if isfile(ef), Le = load(ef);  E = Le.out.arms.eps; end
    % corrected fixed-tf conjugate test on the deepest huber solution:
    conjPass = -1;  nCross = NaN;  ffr = NaN;
    if ~isempty(A.p)
        rec = Lp.R(find(arrayfun(@(r) isequal([r.iD r.iA], jb.cell) && ...
                                 abs(r.gam - jb.gamma) < 1e-9, Lp.R), 1));
        sm = struct('family', 'huber', 'p', A.p(end));
        Yj = A.Y{end};  K = size(Yj, 2);
        yT = cr3bp_minfuel_prop(rec.tf/K, Yj(:,end), false, rec.Tmax, rec.c, rec.muStar, sm);
        seed = struct('tf', rec.tf, 'tGrid', linspace(0, rec.tf, K+1), 'Y', [Yj, yT]);
        try
            [~, it] = ms_minfuel(rec.rv0(:), rec.rvf(:), rec.tf, seed, rec.Tmax, rec.c, rec.muStar, sm, ...
                                 struct('conjTest', true, 'wallSec', wallSec));
            if it.converged
                conjPass = double(it.conj.pass);  nCross = it.conj.nCrossings;  ffr = it.conj.firstFullRank;
            end
            lg('  conj (huber p=%.4g): conv=%d normR=%.1e pass=%d nCross=%d firstFullRank=%d', ...
               A.p(end), it.converged, it.normR, conjPass, nCross, ffr);
        catch ME
            lg('  conj ERROR: %s', ME.message);
        end
    end
    h = struct('cellIdx', jb.cell, 'gamma', jb.gamma, 'seedScale', jb.seedScale, ...
        'pHuber', lastOr(A.p, NaN), 'pEps', lastOr(getf(E,'p'), NaN), ...
        'mfHuber', lastOr(A.mf, NaN), 'mfEps', lastOr(getf(E,'mf'), NaN), ...
        'failsH', A.nFail, 'failsE', getf(E,'nFail'), 'bisH', A.nBisect, 'bisE', getf(E,'nBisect'), ...
        'wallH', A.wallTotal/60, 'wallE', getf(E,'wallTotal')/60, ...
        'coastH', lastOr(A.coastFrac, NaN), 'coastE', lastOr(getf(E,'coastFrac'), NaN), ...
        'acceptH', lastOr(getf(A,'acceptOk'), false), 'conjPass', conjPass, 'nCross', nCross, 'firstFullRank', ffr, ...
        'nRungsH', numel(A.p), 'nRungsE', numel(getf(E,'p')));
    h.dMf = h.mfHuber - h.mfEps;
    if isempty(H), H = h; else, H(end+1) = h; end %#ok<AGROW>
    save(outMat, 'H');
end

lg('HUBER GRID SUMMARY (%.1f min):', toc(tAll)/60);
lg('cell   gam  seed   rungs H/E   p_min H / E        mf H / E             dmf       fails H/E  bis H/E  wall H/E[min]  coast H/E  conj');
for k = 1:numel(H)
    h = H(k);
    lg('(%d,%d) %.2f  %.2g   %2d/%2d   %.4g / %.4g   %.6f / %.6f  %+.1e   %d/%d      %d/%d    %4.1f/%4.1f     %.2f/%.2f   %d', ...
       h.cellIdx, h.gamma, h.seedScale, h.nRungsH, h.nRungsE, h.pHuber, h.pEps, ...
       h.mfHuber, h.mfEps, h.dMf, h.failsH, h.failsE, h.bisH, h.bisE, ...
       h.wallH, h.wallE, h.coastH, h.coastE, h.conjPass);
end
lg('HUBER GRID DONE -> %s', outMat);
end

% ------------------------------------------------------------------------
function v = lastOr(x, d)
% LASTOR  x(end) or d if empty.  INPUTS: x; d.  OUTPUTS: v.
if isempty(x), v = d; else, v = x(end); end
end

function v = getf(s, f)
% GETF  s.(f) or [] if s empty / field absent.  INPUTS: s; f.  OUTPUTS: v.
if isempty(s) || ~isfield(s, f), v = []; else, v = s.(f); end
end

function s = tern(c, a, b)
% TERN  a if c else b.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end

function logmsg(f, s)
% LOGMSG  Append to log file or stdout.  INPUTS: f; s.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
end
