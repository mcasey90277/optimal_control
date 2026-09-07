function V = run_conj_fixedtf_sweep(cfg)
% RUN_CONJ_FIXEDTF_SWEEP  Fixed-tf conjugate-point verdicts for every
% min-energy backbone record, every min-fuel grid record, and the
% gamma-continuation rewalk record -- each RE-SOLVED by multiple shooting
% with keepSTMs so the STMs belong to the exact solution being judged.
%
% This is the named form of the 2026-09-02 sweep (FINDINGS 20), re-run
% after the 2026-09-05 review: the instrument now monitors the FULL state
% block (rows 1:7, P0.3), skips initial-coast structural zeros (P0.4),
% samples through t_f and counts the last bracket (P1.1). Every verdict
% record stores the spec echo, the det samples and the initial costate it
% was computed on, so build_minfuel_catalog can bind verdicts to branches.
%
% INPUTS:
%   cfg - (optional) struct: .logFile [''], .wallSec [600] per solve
%
% OUTPUTS:
%   V - struct array (also saved to direct/results/conj_fixedtf_verdicts.mat):
%       .kind 'energy'|'fuel', .source 'pilot'|'grid'|'rewalk', .cellIdx,
%       .gamma, .p, .converged, .normR, .z [7x1], .conjPass, .nCross,
%       .atFinal, .firstFullRank, .dets [1 x K], .sigRatio, .stateRows,
%       .costateCols, .s1 (throttle at t_1), .mf
%
% REFERENCES:
%   [1] costate_common/ms_conjugate_test.m (the instrument; header = math).
%   [2] DRO_tulip/reviews/minfuel_code_review_2026-09-05.md (P0.3, P0.4, P1.1).

if nargin < 1, cfg = struct(); end
logFile = '';  wallSec = 600;
if isfield(cfg, 'logFile'), logFile = cfg.logFile; end
if isfield(cfg, 'wallSec'), wallSec = cfg.wallSec; end
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'indirect'));
addpath(fullfile(fileparts(here), 'costate_common'));
resDir = fullfile(here, 'direct', 'results');
Lp = load(fullfile(resDir, 'minenergy_pilot.mat'));
Lg = load(fullfile(resDir, 'minfuel_grid.mat'));
outMat = fullfile(resDir, 'conj_fixedtf_verdicts.mat');
mo = struct('conjTest', true, 'wallSec', wallSec);

V = struct([]);
tAll = tic;

%% Min-energy backbone records (p = 1; seed = the pilot's lambda0 flown):
sel = find(arrayfun(@(r) (ischar(r.verdict) && strcmp(r.verdict, 'PASS')) || ...
                         (~ischar(r.verdict) && ~isempty(r.gates) && all(r.gates)), Lp.R));
for ks = 1:numel(sel)
    rec = Lp.R(sel(ks));
    K = rec.ms.K;
    [~, ~, Tf, Yf] = cr3bp_minenergy_prop(rec.tf, [rec.rv0(:); 1; rec.ms.z(:)], false, ...
                                          rec.Tmax, rec.c, rec.muStar);
    tG = linspace(0, rec.tf, K+1);
    seed = struct('tf', rec.tf, 'tGrid', tG, 'Y', interp1(Tf, Yf, tG, 'pchip')');
    t0 = tic;
    [z, it] = ms_minenergy(rec.rv0(:), rec.rvf(:), rec.tf, seed, rec.Tmax, rec.c, rec.muStar, mo);
    V = addRec(V, 'energy', 'pilot', [rec.iD rec.iA], rec.gam, 1, z, it);
    lg('energy (%d,%d)@%.2f: conv=%d normR=%.1e -> %s nCross=%d firstFullRank=%d atFinal=%d (%.0fs)', ...
       rec.iD, rec.iA, rec.gam, it.converged, it.normR, it.conj.verdict, ...
       it.conj.nCrossings, it.conj.firstFullRank, it.conj.atFinal, toc(t0));
    save(outMat, 'V');
end

%% Min-fuel grid records (p = p_floor; seed = the banked junctions):
for k = 1:numel(Lg.G)
    g = Lg.G(k);
    rec = Lp.R(find(arrayfun(@(r) isequal([r.iD r.iA], g.cellIdx) && ...
                             abs(r.gam - g.gamma) < 1e-9, Lp.R), 1));
    V = fuelVerdict(V, 'grid', g.cellIdx, g.gamma, g.pDeepest, g.Y, rec, mo, lg);
    save(outMat, 'V');
end

%% The gamma-continuation rewalk record (FINDINGS 21), if present:
rwFile = fullfile(resDir, 'minfuel_rewalk_c25_g140.mat');
if isfile(rwFile)
    Lr = load(rwFile);  R = Lr.R;
    rec = Lp.R(find(arrayfun(@(r) isequal([r.iD r.iA], R.cellIdx) && ...
                             abs(r.gam - R.gamma) < 1e-9, Lp.R), 1));
    V = fuelVerdict(V, 'rewalk', R.cellIdx, R.gamma, R.pDeepest, R.Y, rec, mo, lg);
    save(outMat, 'V');
end

lg('SWEEP DONE: %d verdicts, %d PASS, %d FAIL, %d unconverged (%.1f min) -> %s', ...
   numel(V), nnz([V.conjPass] == 1 & [V.converged]), nnz([V.conjPass] == 0 & [V.converged]), ...
   nnz(~[V.converged]), toc(tAll)/60, outMat);
end

% ------------------------------------------------------------------------
function V = fuelVerdict(V, src, cell, gam, p, Yj, rec, mo, lg)
% FUELVERDICT  Re-solve one min-fuel record at its p from its junctions and
% record the verdict.  INPUTS: V; src; cell; gam; p; Yj [14 x K]; rec; mo;
% lg.  OUTPUTS: V.
sm = struct('family', 'eps', 'p', p);
K = size(Yj, 2);
yT = cr3bp_minfuel_prop(rec.tf/K, Yj(:,end), false, rec.Tmax, rec.c, rec.muStar, sm);
seed = struct('tf', rec.tf, 'tGrid', linspace(0, rec.tf, K+1), 'Y', [Yj, yT]);
t0 = tic;
[z, it] = ms_minfuel(rec.rv0(:), rec.rvf(:), rec.tf, seed, rec.Tmax, rec.c, rec.muStar, sm, mo);
V = addRec(V, 'fuel', src, cell, gam, p, z, it);
lg('fuel %-6s (%d,%d)@%.2f p=%.4g: conv=%d normR=%.1e -> %s nCross=%d firstFullRank=%d atFinal=%d s1=%.2f (%.0fs)', ...
   src, cell, gam, p, it.converged, it.normR, it.conj.verdict, ...
   it.conj.nCrossings, it.conj.firstFullRank, it.conj.atFinal, it.s(1), toc(t0));
end

function V = addRec(V, kind, src, cell, gam, p, z, it)
% ADDREC  Append one verdict record.  INPUTS: V; kind; src; cell; gam; p;
% z; it (ms info with .conj).  OUTPUTS: V.
c = it.conj;
% A verdict is only a verdict on a CONVERGED extremal (Astra review #2):
% unconverged -> -1 (unknown), regardless of what the instrument said.
cp = -1;  if it.converged && c.tested, cp = double(c.pass); end
mf = NaN;  if isfield(it, 'Yend') && ~isempty(it.Yend), mf = it.Yend(7); end   % FINAL mass, not the last junction start
r = struct('kind', kind, 'source', src, 'cellIdx', cell, 'gamma', gam, 'p', p, ...
    'converged', it.converged, 'normR', it.normR, 'z', z(:), ...
    'conjPass', cp, 'verdict', c.verdict, 'tested', c.tested, ...
    'nCross', c.nCrossings, 'nInterior', c.nInterior, 'atFinal', c.atFinal, ...
    'firstFullRank', c.firstFullRank, 'sampledThrough', c.sampledThrough, ...
    'dets', c.detScaled, 'sigRatio', c.sigRatio, ...
    'stateRows', c.stateRows, 'costateCols', c.costateCols, ...
    's1', it.s(1), 'mf', mf);
if isempty(V), V = r; else, V(end+1) = r; end
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
