function P = ms_refine_catalog(sweepMat, cellsMat, outMat, opts)
% MS_REFINE_CATALOG  Wave 2: multiple-shooting indirect refinement of the map.
%
% For every green cell of the direct phase-sweep, build a multiple-shooting
% seed from the SAVED state+costate trajectory (X and lamDef at all N nodes
% -- stored for exactly this) and solve the PMP min-time problem with
% ms_tfmin. Single shooting (pumpkyn tfMin) was refuted for these seeds:
% collocation duals are ~1e-3 accurate and Lyapunov amplification makes the
% PMP flight miss by 36,000-560,000 km (see seed_quality_12x12.mat).
%
% Every cell is verified two ways after convergence: the multiple-shooting
% residual itself, and a single end-to-end PMP flight from the refined
% lambda0 (meaningful again once lambda0 is ~1e-12 accurate). Results are
% saved after EVERY cell with the sweep's own orbit metadata.
%
% INPUTS:
%   sweepMat - sweep result .mat holding S (map, meta with orbit definition)
%   cellsMat - per-cell products .mat holding CELLS (X, lamDef, tNodes, tf)
%   outMat   - output .mat path (written after every cell)
%   opts     - (optional) struct: .Kladder segment escalation [12 24 48],
%              .wallSec per cell [240],
%              .cells [k x 2] subset of (iD,iA) to run [all green],
%              .verbose [false]
%
% OUTPUTS:
%   P - struct of [nD x nA] arrays: .TFI refined tf, .RES ms residual,
%       .MISSKM flight miss (km), .ITERS, .WALL, .OKI verified; refined
%       costates .LAMI [8 x nD x nA]; .meta copied from the sweep.
%
% REFERENCES:
%   [1] ms_tfmin.m -- the multiple-shooting solver (pumpkyn conventions).
%   [2] sweep_phasing_direct.m -- source of the map and per-cell products.

if nargin < 4, opts = struct(); end
d = @(f,v) fielddef(opts, f, v);
Kladder = d('Kladder', [12 24 48]);   % escalate segmentation until verified
wallSec = d('wallSec', 240);
verbose = d('verbose', false);

M = load(sweepMat);   S = M.S;
CC = load(cellsMat);  CELLS = CC.CELLS;
ob = S.meta.orbit;
muStar = ob.muStar;
g0   = 9.80665*ob.tStar^2/(1000*ob.lStar);
Tmax = (S.meta.thrustN/S.meta.m0kg)*ob.tStar^2/(ob.lStar*1000);
c    = (S.meta.ispS/ob.tStar)*g0;

% endpoint states from the recorded orbit definition (never assumed)
[~, rvD0] = pumpkynPie.cr3bp.getDRO(ob.tauDRO);
rvD0 = pumpkyn.cr3bp.cont_np(rvD0, ob.tauDRO, muStar, 1e-12);
[tD, rvD] = pumpkyn.cr3bp.prop(ob.tauDRO, rvD0, muStar);
[~, rvT0] = pumpkyn.cr3bp.getTulip(ob.tauTulip, ob.NpTulip, ob.pmTulip);
rvT0 = pumpkyn.cr3bp.cont_np(rvT0, ob.tauTulip, muStar, 1e-12);
[tT, rvT] = pumpkyn.cr3bp.prop(ob.tauTulip, rvT0, muStar);

[nD, nA] = size(S.PASS);
todo = d('cells', []);
if isempty(todo)
    [iiD, iiA] = find(S.PASS);
    todo = [iiD, iiA];
end

TFI = nan(nD,nA); RES = TFI; MISSKM = TFI; ITERS = TFI; WALL = TFI; KUSED = TFI;
ACCDZ = TFI;
OKI = false(nD,nA);  LAMI = nan(8,nD,nA);
sD = S.sD; sA = S.sA; meta = S.meta; %#ok<NASGU>
fprintf('MS REFINE: %d cells, K-ladder [%s], %ds/cell budget\n', ...
    size(todo,1), num2str(Kladder), wallSec);
for kc = 1:size(todo,1)
    iD = todo(kc,1);  iA = todo(kc,2);
    cell_ = CELLS{iD,iA};
    if isempty(cell_) || ~S.PASS(iD,iA), continue, end
    rv0 = interp1(tD, rvD, mod(S.sD(iD),1)*tD(end), 'spline');
    rvf = interp1(tT, rvT, mod(S.sA(iA),1)*tT(end), 'spline');

    % seed trajectory per ladder rung: states at nodes, costates at defect
    % rows (node-k convention, validated at the anchor)
    tN = cell_.tNodes(:)';

    t0 = tic;
    try
        for K = Kladder
            tGrid = linspace(0, cell_.tf, K+1);
            Xg = interp1(tN, cell_.X', tGrid, 'pchip')';
            Lg = interp1(tN(1:end-1), cell_.lamDef(1:7,:)', tGrid, 'pchip', 'extrap')';
            seed = struct('tf', cell_.tf, 'tGrid', tGrid, 'Y', [Xg; Lg]);
            [z, info] = ms_tfmin(rv0(1:6), rvf(1:6), seed, Tmax, c, muStar, ...
                struct('wallSec', wallSec, 'verbose', verbose));
            [~, rvFly] = pumpkyn.cr3bp.tfMinProp(z(8), ...
                [rv0(1:6)'; 1; z(1:7)], Tmax, c, muStar);
            missKm = norm(rvFly(end,1:3) - rvf(1:3))*ob.lStar;
            ok = info.converged && missKm < 100;
            if ok
                % step-3 gate IN the driver: the production single-shooting
                % solver must accept the entry essentially unchanged
                evalc('zAcc = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), z, Tmax, c, muStar);');
                accDz = norm(zAcc - z);
                ok = accDz < 1e-6;
            else
                accDz = NaN;
            end
            if ok || toc(t0) > wallSec, break, end
        end
        TFI(iD,iA) = z(8);  RES(iD,iA) = info.normR;
        MISSKM(iD,iA) = missKm;  ITERS(iD,iA) = info.iters;
        WALL(iD,iA) = toc(t0);  OKI(iD,iA) = ok;  KUSED(iD,iA) = K;
        ACCDZ(iD,iA) = accDz;
        if ok, LAMI(:,iD,iA) = z; end
        fprintf('  (%2d,%2d) tfD=%.5f tfI=%.5f  R=%.1e  fly=%.2fkm  %s (%.0fs) [%d/%d]\n', ...
            iD, iA, S.TF(iD,iA), z(8), info.normR, missKm, okstr(ok), ...
            toc(t0), kc, size(todo,1));
    catch ME
        WALL(iD,iA) = toc(t0);
        fprintf('  (%2d,%2d) ERROR (%.0fs): %s\n', iD, iA, toc(t0), ME.message);
    end
    save(outMat, 'TFI','RES','MISSKM','ITERS','WALL','OKI','LAMI','KUSED','ACCDZ', ...
         'sD','sA','meta');                                    % EVERY cell
end
P = struct('TFI',TFI,'RES',RES,'MISSKM',MISSKM,'ITERS',ITERS,'WALL',WALL, ...
           'OKI',OKI,'LAMI',LAMI,'KUSED',KUSED,'ACCDZ',ACCDZ, ...
           'sD',S.sD,'sA',S.sA,'meta',S.meta);
nOK = nnz(OKI);
fprintf('MS REFINE DONE: %d/%d verified.', nOK, size(todo,1));
if nOK > 0
    fprintf('  median |tfD-tfI| = %.2e  fastest refined tf = %.6f', ...
        median(abs(S.TF(OKI) - TFI(OKI))), min(TFI(OKI)));
end
fprintf('\n');
end

% ------------------------------------------------------------------------
function v = fielddef(s, f, v0)
% FIELDDEF  s.(f) if present else v0.  INPUTS: s;f;v0. OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end

function s = okstr(t)
% OKSTR  OK/fail marker. INPUTS: t logical. OUTPUTS: s [char].
if t, s = 'OK'; else, s = 'fail'; end
end
