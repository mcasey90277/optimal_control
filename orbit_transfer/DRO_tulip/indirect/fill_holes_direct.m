function out = fill_holes_direct(catMat, opts)
%% Purpose:
%
%   Fill the HOLES of a packaged phase catalog cell by cell: every grid
%   cell with no certified entry gets a direct HS+Sundman solve warm-started
%   from its nearest certified neighbour (the neighbour's PMP flight,
%   re-boundaried to the hole's departure and arrival states, clearance
%   enforced), the solution is harvested into a multiple-shooting seed and
%   put through the full gate stack (certify_root). Certified points are
%   written as ribs -- one rib per arrival column, in the rib_from_crossing
%   layout -- so the packager takes them like any other rib file.
%
%   This is the cell-by-cell route for the cells the continuation ribs
%   could not reach (a rib that stalls leaves every cell beyond it empty;
%   a column whose spine root's rib stalls at once leaves 22). A direct
%   solve does not care which family the neighbour is on, so it also finds
%   roots the walker's family cannot reach.
%
%% Inputs:
%
%  catMat                   char                    packaged catalog .mat
%                                                   (costate_catalog_*; one
%                                                   sheet)
%  opts                     struct (optional)
%   .out                    char                    rib file to write
%                                                   [<catalog dir>/fine_rib_direct_holes.mat]
%   .logFile                char                    ['' = stdout]
%   .cells                  [k x 2] int             explicit (iD, iA) cells;
%                                                   [] = every hole
%   .improveDays            double                  [0] also re-solve every
%                                                   FILLED cell whose t_f
%                                                   exceeds its faster
%                                                   column neighbour's by
%                                                   more than this (days),
%                                                   seeded from that
%                                                   neighbour; the packager
%                                                   keeps the faster root
%   .maxCells               int                     stop after this many
%                                                   cells [inf]
%   .N                      int                     collocation intervals [800]
%   .clearKm                double                  lunar clearance (radius)
%                                                   [1900]
%   .maxCpuSec              double                  per direct solve [300]
%                                                   (the good solves take
%                                                   1-3 min; the ones that
%                                                   ran 10 were junk)
%   .maxSeeds               int                     neighbours tried per
%                                                   cell [3]; the fastest
%                                                   certified root is kept
%   .wallSec                double                  per certification [600]
%   .K                      int                     shooting segments [24]
%   .anchorMat              char                    IGNORED since 2026-09-18 (the
%                                                   setup is closures-only); was a certified root file
%                                                   for arclength_arrival's
%                                                   setup (only its layout
%                                                   is used) [the 70 mN anchor]
%
%% Outputs:
%
%  out                      struct                  .nHoles .nTried .nCert
%                                                   .file .cells (struct
%                                                   array: iD iA sD sA seed
%                                                   tfDirect tfCert ok
%                                                   reason wall)
%
%% Revision History:
%  M. Casey                                                   (c) 09/15/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

% TEST SEAM: handles to this file's local functions (tests/test_fill_holes_physics)
if ischar(catMat) && strcmp(catMat, 'localfunctions'), out = localHandles(localfunctions);  return, end
if nargin < 2, opts = struct(); end
d = @(f, v) fieldd(opts, f, v);
here = fileparts(mfilename('fullpath'));
D = fullfile(fileparts(fileparts(here)), 'DRO_tulip', 'direct');
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'), fullfile(fileparts(fileparts(here)), 'campaign_common'), D, fullfile(D, 'lib'), fullfile(D, 'certify'));
logFile = d('logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));
N = d('N', 800);  clearKm = d('clearKm', 1900);  maxCpu = d('maxCpuSec', 300);
wallSec = d('wallSec', 600);  K = d('K', 24);  maxCells = d('maxCells', inf);
maxSeeds = d('maxSeeds', 3);  improveDays = d('improveDays', 0);
[catDir, ~] = fileparts(catMat);
outFile = d('out', fullfile(catDir, 'fine_rib_direct_holes.mat'));

% ---- the catalog: grid, holes, certified neighbours -----------------------
L = load(catMat);  fn = fieldnames(L);  cat_ = L.(fn{1});
assert(isscalar(cat_.sheets), 'fill_holes_direct: one sheet expected');
sh = cat_.sheets(1);
sD = sh.sD_frac(:).';  sA = sh.sA_frac(:).';  nD = numel(sD);  nA = numel(sA);
has = sh.has_solution(:, :, 1);  tfnd = sh.tf_nd(:, :, 1);  idx = sh.entry_index(:, :, 1);
z8 = sh.z8;
tStar = cat_.constants.tStar_s;  rMoonKm = 1737.4;
cells = d('cells', []);
if isempty(cells), [iDh, iAh] = find(~has);  cells = [iDh, iAh]; end
nHoles = nnz(~has);
lg('fill_holes_direct: %d x %d grid, %d holes, %d cells to try', nD, nA, nHoles, size(cells, 1));

% ---- the physics, FROM THE CATALOG (engine, orbits, spine) ---------------
[B, ~] = arclength_arrival('setup', physicsFromCatalog(cat_, sh, sD(1)));
pool = capped_pool();
floorKm = clearKm - rMoonKm;
problem = B.problem;
solveOpts = struct('N', N, 'K', K, 'clearKm', clearKm, 'maxCpuSec', maxCpu, 'wallSec', wallSec, 'pool', pool);

% ---- resume: keep what an earlier run certified --------------------------
R = struct('j', {}, 'sA', {}, 'pts', {}, 'stop', {}, 'nSolve', {});
rec = struct('iD', {}, 'iA', {}, 'sD', {}, 'sA', {}, 'seed', {}, 'tfDirect', {}, 'tfCert', {}, 'ok', {}, 'reason', {}, 'wall', {});
if isfile(outFile)
    Lo = load(outFile);
    if isfield(Lo, 'R'), R = Lo.R; end
    if isfield(Lo, 'rec'), rec = Lo.rec; end
    lg('  resuming: %d rib(s), %d cell record(s) on disk', numel(R), numel(rec));
end
% a cell certified earlier is done and becomes a SEED for its neighbours;
% a cell that failed is retried -- a neighbour certified since may seed it
done = false(nD, nA);
for k = 1:numel(rec)
    if rec(k).ok, done(rec(k).iD, rec(k).iA) = true; end
end
for m = 1:numel(R)
    for q = 1:numel(R(m).pts)
        [has, tfnd, idx, z8] = admit(has, tfnd, idx, z8, sD, sA, R(m).pts(q));
    end
end
rec = rec([rec.ok]);                          % failed records are re-made below

% the IMPROVE selection (AFTER the rib file's points are admitted, so a
% cell improved by an earlier run is the seed of the next one up): a filled cell far slower than the column neighbour
% below or above it (a rib of a slower family filled the column while a
% faster one stalled). Column-major order, so the chain runs upward.
if improveDays > 0
    tfDays = tfnd * tStar/86400;
    slow = false(nD, nA);
    for iA_ = 1:nA
        for iD_ = 1:nD
            if ~has(iD_, iA_), continue, end
            nbv = [mod(iD_ - 2, nD) + 1, mod(iD_, nD) + 1];
            nbv = nbv(has(nbv, iA_));
            if ~isempty(nbv) && tfDays(iD_, iA_) - min(tfDays(nbv, iA_)) > improveDays, slow(iD_, iA_) = true; end
        end
    end
    [iDs, iAs] = find(slow);
    cells = [cells; iDs, iAs];
    lg('  improve: %d filled cell(s) slower than a column neighbour by > %.1f d', nnz(slow), improveDays);
end

nTried = 0;  nCert = 0;
queued = false(nD, nA);  queued(sub2ind([nD nA], cells(:, 1), cells(:, 2))) = true;
kc = 0;
while kc < size(cells, 1)                       % the list grows as improvements chain
    kc = kc + 1;
    if nTried >= maxCells, break, end
    iD = cells(kc, 1);  iA = cells(kc, 2);
    if done(iD, iA), continue, end
    % the nearest certified neighbours: the same column's (they share the
    % arrival geometry -- a rib IS this continuation) before the same row's
    % (they share the departure state), each group fastest first; every
    % seed tried (up to maxSeeds) and the fastest certified root kept
    nbCol = [mod(iD - 2, nD) + 1, iA;  mod(iD, nD) + 1, iA];
    nbRow = [iD, mod(iA - 2, nA) + 1;  iD, mod(iA, nA) + 1];
    nbCol = nbCol(has(sub2ind([nD nA], nbCol(:, 1), nbCol(:, 2))), :);
    nbRow = nbRow(has(sub2ind([nD nA], nbRow(:, 1), nbRow(:, 2))), :);
    if has(iD, iA)                                  % an improve cell: faster seeds only
        nbCol = nbCol(tfnd(sub2ind([nD nA], nbCol(:, 1), nbCol(:, 2))) < tfnd(iD, iA), :);
        nbRow = nbRow(tfnd(sub2ind([nD nA], nbRow(:, 1), nbRow(:, 2))) < tfnd(iD, iA), :);
    end
    [~, oc] = sort(tfnd(sub2ind([nD nA], nbCol(:, 1), nbCol(:, 2))));
    [~, orw] = sort(tfnd(sub2ind([nD nA], nbRow(:, 1), nbRow(:, 2))));
    nb = [nbCol(oc, :); nbRow(orw, :)];
    nb = nb(1:min(end, maxSeeds), :);
    if isempty(nb)
        rec(end+1) = cellRec(iD, iA, sD(iD), sA(iA), [], NaN, NaN, false, 'no certified neighbour', 0);
        lg('  cell (%2d,%2d) sD %.4f sA %.4f: no certified neighbour', iD, iA, sD(iD), sA(iA));
        continue
    end
    t0 = tic;  okCell = false;  reason = '';  tfD = NaN;  tfC = NaN;  seedUsed = [];
    best = [];  reasons = {};  bestSeed = [];  bestSeedTf = NaN;  nRaced = 0;
    if has(iD, iA), reasons{end+1} = sprintf('improve: cell holds %.3f d', tfnd(iD, iA)*tStar/86400); end
    for kn = 1:size(nb, 1)
        nD_ = nb(kn, 1);  nA_ = nb(kn, 2);  seedUsed = [sD(nD_), sA(nA_)];
        zN = z8(:, idx(nD_, nA_));
        rvN = B.stateD(sD(nD_));  rvN = rvN(1:6);
        try
            % the neighbour's PMP flight as the warm start (direct_cell_solve)
            [C, dinfo] = direct_cell_solve(zN, rvN, sD(iD), sA(iA), B, solveOpts);
            tfD = dinfo.tfDirectDays;
            if ~dinfo.success
                reasons{end+1} = sprintf('direct solve failed (seed %.4f,%.4f)', seedUsed);  continue
            end
            if ~C.ok && contains(C.reason, 'clearance floor')
                reasons{end+1} = sprintf('%s (seed %.4f,%.4f)', C.reason, seedUsed);  continue
            end
            nRaced = nRaced + 1;
            if C.ok
                if isempty(best) || C.tfDays < best.tfDays
                    best = C;  tfC = C.tfDays;  bestSeed = seedUsed;  bestSeedTf = tfnd(nD_, nA_)*tStar/86400;
                end
                reasons{end+1} = sprintf('certified %.3f d (seed %.4f,%.4f)', C.tfDays, seedUsed);
            else
                reasons{end+1} = sprintf('%s (direct %.3f d, seed %.4f,%.4f)', C.reason, tfD, seedUsed);
            end
        catch ME
            reasons{end+1} = sprintf('ERROR %s (seed %.4f,%.4f)', ME.message, seedUsed);
        end
    end
    if ~isempty(best) && has(iD, iA) && best.z(8) >= tfnd(iD, iA)
        reasons{end+1} = sprintf('not faster than the cell''s %.3f d', tfnd(iD, iA)*tStar/86400);  best = [];
    end
    if ~isempty(best)
        clause = sprintf('direct cell solve seeded from (%.4f, %.4f) %.3f d, %d seed(s) raced', bestSeed, bestSeedTf, nRaced);
        if has(iD, iA), clause = sprintf('%s, replaced %.3f d', clause, tfnd(iD, iA)*tStar/86400); end
        best.note = join_note(clause, best);
        okCell = true;  R = putPoint(R, iA, sA(iA), best);  tfC = best.tfDays;
        % the new root seeds the cells still to come (its column's next
        % cell above all: that is how a rib chains)
        [has, tfnd, idx, z8] = admit(has, tfnd, idx, z8, sD, sA, best);
        % and a column neighbour now far slower than it joins the queue --
        % this is how the improve pass walks a whole column from one fast root
        if improveDays > 0
            for nbv = [mod(iD - 2, nD) + 1, mod(iD, nD) + 1]
                if has(nbv, iA) && ~queued(nbv, iA) && (tfnd(nbv, iA) - best.z(8))*tStar/86400 > improveDays
                    cells(end+1, :) = [nbv, iA];  queued(nbv, iA) = true;
                    lg('    queued (%2d,%2d): %.3f d against the new %.3f d', nbv, iA, tfnd(nbv, iA)*tStar/86400, best.tfDays);
                end
            end
        end
    end
    reason = strjoin(reasons, ' | ');
    nTried = nTried + 1;  nCert = nCert + okCell;
    rec(end+1) = cellRec(iD, iA, sD(iD), sA(iA), seedUsed, tfD, tfC, okCell, reason, toc(t0));
    lg('  cell (%2d,%2d) sD %.4f sA %.4f: %s  direct %.3f d  cert %.3f d  %s  (%.0f s)', iD, iA, sD(iD), sA(iA), ...
       pick(okCell, 'CERTIFIED', 'no'), tfD, tfC, reason, toc(t0));
    tmp = [outFile '.part'];
    save(tmp, 'R', 'rec', 'problem');  publish_atomic(tmp, outFile);
end
lg('fill_holes_direct: %d holes, %d tried, %d certified -> %s', nHoles, nTried, nCert, outFile);
out = struct('nHoles', nHoles, 'nTried', nTried, 'nCert', nCert, 'file', outFile, 'cells', rec);
end

function s = join_note(prefix, C)
% JOIN_NOTE  The entry's provenance note: this producer's clause in front
% of whatever the certifier (or an earlier producer) already wrote.
% INPUTS: prefix (char); C (struct, .note optional).  OUTPUTS: s.
parts = {prefix};
if isfield(C, 'note') && ~isempty(C.note), parts{end+1} = C.note; end
s = strjoin(parts, ' | ');
end

function [has, tfnd, idx, z8] = admit(has, tfnd, idx, z8, sD, sA, C)
% ADMIT  Enter a certified point into the seed grid.  INPUTS: has; tfnd;
% idx; z8; sD; sA; C (certify_root output with .sD .sA .z).  OUTPUTS: the
% four grids, updated.
iD = find(abs(mod(sD - C.sD + 0.5, 1) - 0.5) < 1e-8, 1);
iA = find(abs(mod(sA - C.sA + 0.5, 1) - 0.5) < 1e-8, 1);
if isempty(iD) || isempty(iA), return, end
if has(iD, iA) && tfnd(iD, iA) <= C.z(8), return, end
z8(:, end+1) = C.z(:);
has(iD, iA) = true;  tfnd(iD, iA) = C.z(8);  idx(iD, iA) = size(z8, 2);
end

function R = putPoint(R, j, sAj, C)
% PUTPOINT  Add a certified point to column j's direct rib.  INPUTS: R; j;
% sAj; C (certify_root output).  OUTPUTS: R.
m = find([R.j] == j, 1);
if isempty(m)
    R(end+1) = struct('j', j, 'sA', sAj, 'pts', C, 'stop', 'direct cells', 'nSolve', 1);
else
    R(m).pts = append_point(R(m).pts, C);  R(m).nSolve = R(m).nSolve + 1;
end
end

function r = cellRec(iD, iA, sDv, sAv, seed, tfD, tfC, ok, reason, wall)
% CELLREC  One cell's record.  INPUTS: as named.  OUTPUTS: r struct.
r = struct('iD', iD, 'iA', iA, 'sD', sDv, 'sA', sAv, 'seed', seed, 'tfDirect', tfD, ...
           'tfCert', tfC, 'ok', ok, 'reason', reason, 'wall', wall);
end

function v = pick(c, a, b)
% PICK  Ternary.  INPUTS: c; a; b.  OUTPUTS: v.
if c, v = a; else, v = b; end
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end

function logmsg(f, s)
% LOGMSG  Append to a log file or stdout.  INPUTS: f; s.  OUTPUTS: none.
line = sprintf('%s %s', char(datetime('now', 'Format', 'HH:mm:ss')), s);
if isempty(f), fprintf('%s\n', line);
else, fid = fopen(f, 'a');  fprintf(fid, '%s\n', line);  fclose(fid);  fprintf('%s\n', line);
end
end

function pts = append_point(pts, C)
% APPEND_POINT  Append a certified point to a rib, harmonising fields: a
% rib resumed from a checkpoint written before the .note field existed
% must still accept new points.  INPUTS: pts (struct array or []); C.
% OUTPUTS: pts.
if isempty(pts), pts = C;  return, end
for f = setdiff(fieldnames(C), fieldnames(pts))', [pts.(f{1})] = deal([]); end
for f = setdiff(fieldnames(pts), fieldnames(C))', C.(f{1}) = []; end
C = orderfields(C, pts);
pts(end+1) = C;
end

% ------------------------------------------------------------------------
function so = physicsFromCatalog(cat_, sh, sD1)
% PHYSICSFROMCATALOG  The setup request for this catalog's problem: engine,
% orbits and the spine's departure phase, CLOSURES ONLY. The filler needs
% the endpoint closures and the propulsion constants, never an anchor; asking
% for one made it re-polish the shipped 70 mN anchor, which fails for any
% other engine, orbit pair or departure phase (FINDINGS 78).
% INPUTS: cat_ (catalog struct); sh (its sheet); sD1 (spine departure phase).
% OUTPUTS: so struct for arclength_arrival('setup', so).
so = struct('thrustN', cat_.rungs_N(1), 'ispS', cat_.thruster.isp_s, 'm0kg', cat_.thruster.m0_kg, ...
            'tauDRO', sh.tauDRO, 'NpTulip', sh.Np, 'pmTulip', sh.pm, 'sD', sD1, 'physicsOnly', true);
end

% ------------------------------------------------------------------------
function H = localHandles(fh)
% LOCALHANDLES  This file's local functions as a struct of handles keyed by
% name -- the TEST SEAM: a test calls the real helper, not a copy of it.
% INPUTS: fh (cell of handles, from localfunctions).  OUTPUTS: H struct.
H = struct();
for k = 1:numel(fh), H.(func2str(fh{k})) = fh{k}; end
end
