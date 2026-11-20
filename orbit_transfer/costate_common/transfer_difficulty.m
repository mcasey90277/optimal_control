function D = transfer_difficulty(sheetMat, opts)
%% Purpose:
%
%   THE SCORECARD of a thrust-ladder sheet: how hard was this orbit pair
%   to solve? It is pure post-processing of a sheet written by
%   costate_common/thrust_ladder_library (and, when present, its
%   <tag>_direct.mat sidecar); no CasADi, no pumpkyn. Every statistic is
%   defined below exactly, so two pairs can be compared with numbers.
%
%   DEFINITIONS (per rung k of the sheet, over its nCells = nD*nA cells):
%
%     accepted        OK true. nAccepted, shareAccepted = nAccepted/nCells.
%     nDirect         accepted cells whose direct solution is stored
%                     (sidecar DDONE true); NaN with no sidecar.
%     tfMedianND, tfMinND, tfMaxND   over accepted TF (ND time).
%     tfSpread        (tfMax - tfMin)/tfMedian over accepted TF.
%     revMedian       median revolutions about the primary, from the
%                     DIRECT solution of each accepted cell with DDONE
%                     true: positions DX(1:3,:) (rotating frame, ND,
%                     barycentre origin) at t = linspace(0, DTF, 401) are
%                     rotated to the inertial frame by angle t about z;
%                     the primary's inertial position (rotating position
%                     [-muStar 0 0] Earth or [1-muStar 0 0] Moon, rotated
%                     the same way) is subtracted; the relative motion is
%                     projected on the plane normal to its mean
%                     angular-momentum direction (sum over samples of
%                     r_k x r_(k+1)); the polar angle in that plane is
%                     unwrapped and revs = |total angle|/(2 pi). NaN with
%                     no sidecar or no such cell.
%     dVmedianKms     median over accepted TF of the all-burn delta-V
%                     c ln(1/(1 - Tnd tf/cnd)) in km/s, the bookkeeping of
%                     pair_transfer's mass_of (nd_propulsion with ispS and
%                     m0kg from the sheet meta); NaN entries (mass
%                     exhausted) are skipped.
%     wallMedianS     median over accepted WALL (seconds).
%     roughness       for each accepted cell with at least two accepted
%                     torus neighbours among its four (wrap in both
%                     phases; distinct cells only, the cell itself
%                     excluded, so nD or nA <= 2 never double-counts): |TF - mean(neighbour TF)|/TF; the rung's
%                     roughness is the median over such cells, NaN if
%                     none.
%     shareRefuted    (only with opts.conjPass) count of conjPass == 0
%                     over accepted cells whose conjPass is 0 or 1;
%                     nUndecided counts accepted cells with any other
%                     conjPass value.
%
%   LADDER (once per sheet): lowestRungAllCells = lowest thrust whose
%   shareAccepted >= opts.shareAll (NaN if none); lowestRungAnyCell =
%   lowest thrust with any accepted cell (NaN if none); nFullLadders =
%   cells accepted at every rung; nEmptyCells = cells accepted at none;
%   attMedian = median of ATT over the cells with ATT > 0 (the engine
%   writes ATT before solving, so 0 = never attempted; NaN if none) and
%   nAttempted = number of such cells (ATT is per cell, not per rung, so
%   both are reported once here).
%
%   Self-demo (nargin == 0): a small synthetic sheet in tempname, its
%   scorecard printed, the folder removed.
%
%% Inputs:
%
%   sheetMat  char   path of the sheet .mat (fields rungs sD sA TF FLYKM
%                    ACCDZ RES WALL OK Z8 ATT meta). The sidecar is
%                    <sheetMat without .mat>_direct.mat.
%   opts      struct (optional)
%     .primary   'moon' | 'earth'  body the revolutions are about ['moon']
%     .print     logical           print the table [true]
%     .shareAll  scalar            share for lowestRungAllCells [0.9]
%     .conjPass  [nD x nA x nR]    conjugate verdicts 1/0/-1 [absent]
%
%% Outputs:
%
%   D  struct
%     .perRung  table, one row per rung in sheet order: thrustN nCells
%               nAccepted shareAccepted nDirect tfMedianND tfMinND
%               tfMaxND tfSpread revMedian dVmedianKms wallMedianS
%               roughness (+ shareRefuted nUndecided with conjPass)
%     .ladder   struct lowestRungAllCells lowestRungAnyCell nFullLadders
%               nEmptyCells nAttempted attMedian
%     .meta     struct families parameters grid (from the sheet meta)
%     .lines    cell of char, the printed table and summary lines
%
%% Revision History:
%
%   M. Casey (c) 11/19/2026
%
%   Copyright Coorbital Inc.
%

if nargin == 0
    tmp = tempname;  mkdir(tmp);
    cleaner = onCleanup(@() rmdir(tmp, 's'));
    D = transfer_difficulty(demo_sheet(tmp), struct('print', true));
    clear cleaner
    return
end
if nargin < 2 || isempty(opts),  opts = struct();  end
primary  = getopt(opts, 'primary',  'moon');
doPrint  = getopt(opts, 'print',    true);
shareAll = getopt(opts, 'shareAll', 0.9);
hasConj  = isfield(opts, 'conjPass') && ~isempty(opts.conjPass);
if ~any(strcmp(primary, {'moon', 'earth'}))
    error('transfer_difficulty:primary', 'opts.primary must be ''moon'' or ''earth'' (got "%s")', char(primary));
end

S      = load(sheetMat);
rungs  = S.rungs(:)';
nR     = numel(rungs);
nD     = numel(S.sD);
nA     = numel(S.sA);
nCells = nD*nA;
OK     = reshape(S.OK,   nD, nA, nR);                   % a one-rung sheet drops its trailing singleton
TF     = reshape(S.TF,   nD, nA, nR);
WALL   = reshape(S.WALL, nD, nA, nR);
meta   = S.meta;
muStar = fld(meta, 'muStar', 0.012150585609624);
lStar  = fld(meta, 'lStar',  389703.264829278);
tStar  = fld(meta, 'tStar',  382981.289129055);
ispS   = meta.ispS;
m0kg   = meta.m0kg;
if hasConj,  CP = reshape(opts.conjPass, nD, nA, nR);  end

% sidecar (optional)
[sdPath, sdName] = fileparts(sheetMat);
dirMat = fullfile(sdPath, [sdName '_direct.mat']);
hasDir = exist(dirMat, 'file') == 2;
if hasDir
    Q     = load(dirMat, 'DX', 'DTF', 'DDONE');
    DX    = reshape(Q.DX,    7, [], nD, nA, nR);
    DTF   = reshape(Q.DTF,   nD, nA, nR);
    DDONE = reshape(Q.DDONE, nD, nA, nR);
    if strcmp(primary, 'earth'),  xPrim = -muStar;  else,  xPrim = 1 - muStar;  end
end

nan1 = nan(nR, 1);
thrustN = rungs(:);  nCellsV = nCells*ones(nR, 1);
nAcc = nan1;  share = nan1;  nDir = nan1;
tfMed = nan1;  tfMin = nan1;  tfMax = nan1;  spread = nan1;
revMed = nan1;  dvMed = nan1;  wallMed = nan1;  rough = nan1;
shRef = nan1;  nUnd = nan1;
for k = 1:nR
    A  = OK(:,:,k);
    T  = TF(:,:,k);
    tf = T(A);
    nAcc(k)  = sum(A(:));
    share(k) = nAcc(k)/nCells;
    if ~isempty(tf)
        tfMed(k) = median(tf);  tfMin(k) = min(tf);  tfMax(k) = max(tf);
        spread(k) = (tfMax(k) - tfMin(k))/tfMed(k);
        ndp = nd_propulsion(rungs(k), ispS, m0kg, lStar, tStar);
        mf  = 1 - ndp.Tnd*tf/ndp.cnd;
        dv  = nan(size(tf));
        good = mf > 0;
        dv(good) = ndp.cnd*log(1./mf(good))*lStar/tStar;
        dvMed(k) = median(dv, 'omitnan');
        W = WALL(:,:,k);
        wallMed(k) = median(W(A), 'omitnan');
    end
    rough(k) = torus_roughness(A, T);
    if hasDir
        use = A & DDONE(:,:,k);
        nDir(k) = sum(use(:));
        idx = find(use)';
        rv  = nan(1, numel(idx));
        for q = 1:numel(idx)
            [a, b] = ind2sub([nD nA], idx(q));
            rv(q) = revolutions(DX(1:3,:,a,b,k), DTF(a,b,k), xPrim);
        end
        revMed(k) = median(rv, 'omitnan');
    end
    if hasConj
        C   = CP(:,:,k);
        dec = A & (C == 0 | C == 1);
        nUnd(k) = sum(A(:)) - sum(dec(:));
        if any(dec(:)),  shRef(k) = sum(C(dec) == 0)/sum(dec(:));  end
    end
end

D.perRung = table(thrustN, nCellsV, nAcc, share, nDir, tfMed, tfMin, tfMax, spread, ...
                  revMed, dvMed, wallMed, rough, ...
                  'VariableNames', {'thrustN','nCells','nAccepted','shareAccepted','nDirect', ...
                  'tfMedianND','tfMinND','tfMaxND','tfSpread','revMedian','dVmedianKms', ...
                  'wallMedianS','roughness'});
if hasConj
    D.perRung.shareRefuted = shRef;
    D.perRung.nUndecided   = nUnd;
end

qual = share >= shareAll;
anyR = nAcc > 0;
D.ladder = struct( ...
    'lowestRungAllCells', minorempty(rungs(qual)), ...
    'lowestRungAnyCell',  minorempty(rungs(anyR)), ...
    'nFullLadders',       sum(sum(all(OK, 3))), ...
    'nEmptyCells',        sum(sum(~any(OK, 3))), ...
    'nAttempted',         sum(S.ATT(:) > 0), ...
    'attMedian',          att_median(S.ATT));

D.meta = struct('families', struct('dep', fld(meta, 'depFamily', ''), 'arr', fld(meta, 'arrFamily', ''), ...
                                   'depParams', fld(meta, 'depParams', struct()), ...
                                   'arrParams', fld(meta, 'arrParams', struct())), ...
                'parameters', struct('ispS', ispS, 'm0kg', m0kg, 'N', fld(meta, 'N', NaN), ...
                                     'floorKm', fld(meta, 'floorKm', NaN), 'primary', primary, ...
                                     'shareAll', shareAll, 'hasSidecar', hasDir), ...
                'grid', struct('nD', nD, 'nA', nA, 'sD', S.sD, 'sA', S.sA, 'rungs', rungs));

D.lines = build_lines(D, hasConj, sheetMat);
if doPrint
    fprintf('%s\n', D.lines{:});
end
end

% ------------------------------------------------------------------------
function r = torus_roughness(A, T)
% TORUS_ROUGHNESS  Median over accepted cells with >= 2 accepted torus
%   neighbours of |TF - mean(neighbour TF)|/TF.  The neighbours of a cell
%   are the DISTINCT cells among its four torus neighbours (wrap in both
%   phases), the cell itself excluded; so on a grid with nD or nA <= 2 a
%   neighbour reached twice counts once and a wrap onto itself counts
%   not at all.  INPUTS: A logical [nD x nA] accepted; T [nD x nA] TF.
%   OUTPUT: scalar, NaN if no such cell.
[nD, nA] = size(A);
vals = [];
for c = find(A(:))'
    [a, b] = ind2sub([nD nA], c);
    ra = mod(a + [-1 1] - 1, nD) + 1;
    cb = mod(b + [-1 1] - 1, nA) + 1;
    nb = unique([sub2ind([nD nA], ra, [b b]), sub2ind([nD nA], [a a], cb)]);
    nb = nb(nb ~= c & A(nb));
    if numel(nb) >= 2
        vals(end+1) = abs(T(c) - mean(T(nb)))/T(c);
    end
end
if isempty(vals),  r = NaN;  else,  r = median(vals);  end
end

% ------------------------------------------------------------------------
function revs = revolutions(Xrot, tf, xPrim)
% REVOLUTIONS  Revolutions about the primary in the inertial frame.
%   INPUTS: Xrot [3 x nt] rotating-frame positions (ND); tf ND duration;
%   xPrim rotating-frame x of the primary.  OUTPUT: revs scalar (NaN when
%   the relative motion has no defined plane).
Xrot = reshape(Xrot, 3, []);
nt   = size(Xrot, 2);
t    = linspace(0, tf, nt);
c    = cos(t);  s = sin(t);
rel  = zeros(3, nt);
for kk = 1:nt
    Rz = [c(kk) -s(kk) 0; s(kk) c(kk) 0; 0 0 1];        % rotating -> inertial
    rel(:,kk) = Rz*Xrot(:,kk) - Rz*[xPrim; 0; 0];
end
H = sum(cross(rel(:,1:end-1), rel(:,2:end), 1), 2);
nH = sqrt(sum(H.^2));
if ~isfinite(nH) || nH == 0,  revs = NaN;  return,  end
n  = H/nH;
e1 = rel(:,1) - n*(n'*rel(:,1));
ne = sqrt(sum(e1.^2));
if ne == 0,  revs = NaN;  return,  end
e1 = e1/ne;
e2 = cross(n, e1);
ang = unwrap(atan2(e2'*rel, e1'*rel));
revs = abs(ang(end) - ang(1))/(2*pi);
end

% ------------------------------------------------------------------------
function L = build_lines(D, hasConj, sheetMat)
% BUILD_LINES  The aligned table (one row per rung) and summary lines.
P = D.perRung;
[~, nm] = fileparts(sheetMat);
L = {sprintf('Transfer difficulty: %s  (%s -> %s, %d x %d cells, primary %s)', nm, ...
             D.meta.families.dep, D.meta.families.arr, D.meta.grid.nD, D.meta.grid.nA, ...
             D.meta.parameters.primary)};
hdr = sprintf('%8s %6s %7s %6s %7s %9s %9s %9s %7s %8s %9s %9s %9s', 'thrustN', 'nCell', 'nAcc', ...
              'share', 'nDirect', 'tfMedND', 'tfMinND', 'tfMaxND', 'spread', 'revMed', 'dVmedKms', ...
              'wallMedS', 'roughness');
if hasConj,  hdr = [hdr sprintf(' %9s %6s', 'shRefuted', 'nUndec')];  end
L{end+1} = hdr;
for k = 1:height(P)
    row = sprintf('%8.4g %6d %7d %6.3f %7s %9.4f %9.4f %9.4f %7.3f %8.3f %9.4f %9.1f %9.4f', ...
                  P.thrustN(k), P.nCells(k), P.nAccepted(k), P.shareAccepted(k), nanstr(P.nDirect(k), '%d'), ...
                  P.tfMedianND(k), P.tfMinND(k), P.tfMaxND(k), P.tfSpread(k), P.revMedian(k), ...
                  P.dVmedianKms(k), P.wallMedianS(k), P.roughness(k));
    if hasConj
        row = [row sprintf(' %9.4f %6d', P.shareRefuted(k), P.nUndecided(k))];
    end
    L{end+1} = row;
end
Ld = D.ladder;
L{end+1} = sprintf('lowest rung with share >= %.2f of cells: %s N;  lowest rung with any cell: %s N', ...
                   D.meta.parameters.shareAll, nanstr(Ld.lowestRungAllCells, '%.4g'), nanstr(Ld.lowestRungAnyCell, '%.4g'));
L{end+1} = sprintf('full ladders (accepted at every rung): %d;  empty cells (accepted at none): %d;  median attempts over %d attempted cells: %.3g', ...
                   Ld.nFullLadders, Ld.nEmptyCells, Ld.nAttempted, Ld.attMedian);
end

% ------------------------------------------------------------------------
function s = nanstr(x, f)
% NANSTR  Format a scalar, 'NaN' text for NaN.
if isnan(x),  s = 'NaN';  else,  s = sprintf(f, x);  end
end

% ------------------------------------------------------------------------
function m = minorempty(x)
% MINOREMPTY  Minimum of a vector, NaN when empty.
if isempty(x),  m = NaN;  else,  m = min(x);  end
end

% ------------------------------------------------------------------------
function v = getopt(o, name, default)
% GETOPT  Option with default.
if isfield(o, name) && ~isempty(o.(name)),  v = o.(name);  else,  v = default;  end
end

% ------------------------------------------------------------------------
function v = fld(s, name, default)
% FLD  Struct field with default.
if isfield(s, name),  v = s.(name);  else,  v = default;  end
end

% ------------------------------------------------------------------------
function f = demo_sheet(tmp)
% DEMO_SHEET  Synthetic 6 x 5 x 4 sheet for the self-demo.
nD = 6;  nA = 5;  rungs = [15 10 5 2];
OK = false(nD, nA, 4);  TF = nan(nD, nA, 4);  WALL = nan(nD, nA, 4);
for k = 1:4
    OK(:,:,k) = mod((1:nD)'*3 + (1:nA)*7, 10) < 10 - 2*(k-1);   % fixed pattern
    TF(:,:,k) = (1 + 0.1*sin((1:nD)'*(1:nA)))*(0.5 + 0.4*(k-1));
    WALL(:,:,k) = 30*k;
end
TF(~OK) = NaN;
FLYKM = nan(nD, nA, 4);  Z8 = nan(8, nD, nA, 4);
ATT = 1 + mod((1:nD)' + (1:nA), 3);  sD = (0:nD-1)/nD;  sA = (0:nA-1)/nA;
meta = struct('muStar', 0.012150585609624, 'lStar', 389703.264829278, 'tStar', 382981.289129055, ...
              'depFamily', 'demo', 'arrFamily', 'demo', 'ispS', 1710, 'm0kg', 150, 'N', 400);
ACCDZ = FLYKM;  RES = FLYKM;
f = fullfile(tmp, 'demo_sheet.mat');
save(f, 'rungs', 'sD', 'sA', 'TF', 'FLYKM', 'ACCDZ', 'RES', 'WALL', 'OK', 'Z8', 'ATT', 'meta');
end

% ------------------------------------------------------------------------
function m = att_median(ATT)
% ATT_MEDIAN  Median attempts over the cells that were attempted (ATT > 0;
%   the engine writes ATT before solving, so 0 = never attempted).  NaN if none.
v = ATT(ATT > 0);
if isempty(v),  m = NaN;  else,  m = median(v);  end
end
