function P = densify_ladder(ladderMat, opts)
%% Purpose:
%
%   Fills the PER-RUNG GAPS in the thrust-ladder sheet: (phase pair, thrust)
%   combinations where the pair solves at OTHER thrusts but this rung is
%   missing -- the red cells of the coverage torus. Never-solved pairs are a
%   separate, harder campaign; this one only densifies where success nearby
%   is already proven.
%
%   Two fill routes per gap, tried in order:
%
%     A  SAME-CELL RUNG CONTINUATION  the nearest converged rung of the same
%        phase pair (in log-thrust) supplies its converged solution; the
%        trajectory is time-stretched onto the gap rung's expected flight
%        time (several stretch exponents tried, fastest convergent kept --
%        t_f ~ T^-p with p steepening from ~0.6 to ~1.7+ down the band) and
%        refined by ms_tfmin.
%
%     B  TORUS-NEIGHBOR SEED  the nearest phase pair (torus metric) that IS
%        converged at this rung supplies its trajectory; with this cell's
%        own endpoints it seeds ms_tfmin directly. Neighboring phasings on
%        one family differ by a small deformation, which multiple shooting
%        absorbs comfortably.
%
%   Every fill passes the same three gates as the rest of the library:
%   ms residual, flown arrival within gateKm, and pumpkyn.cr3bp.tfMin
%   accepting the entry unchanged. The sheet is updated IN PLACE (backup
%   written first) and saved after every fill.
%
%% Inputs:
%
%  ladderMat                char                    thrust_ladder_12x12.mat
%                                                   (modified in place)
%
%  opts                     struct                  (all optional)
%   .rungSel                vector                  thrust rungs (N) to
%                                                   densify [all]
%   .gateKm                 double                  flown gate, km    [100]
%   .accTol                 double                  tfMin acceptance  [1e-6]
%   .msWallS                double                  ms budget/attempt [120]
%   .K                     double                  ms segments       [24]
%   .tfExpList              vector                  stretch exponents
%                                                   [0.56 1.0 1.7 2.4]
%   .maxCells, .batchSec    double                  batch control     [inf]
%   .maxAtt                 double                  attempts per gap  [2]
%   .logFile                char                    append-mode log [stdout]
%
%% Outputs:
%
%  P                        struct                  the updated sheet (same
%                                                   schema, plus ATTD gap-
%                                                   attempt counter)
%
%% Revision History:
%  M. Casey                                                   (c) 08/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
d = @(f,v) fdef(opts, f, v);
rungSel  = d('rungSel', []);
gateKm   = d('gateKm', 100);
accTol   = d('accTol', 1e-6);
msWallS  = d('msWallS', 120);
K        = d('K', 24);
tfExpList= d('tfExpList', [0.56 1.0 1.7 2.4]);
maxCells = d('maxCells', inf);
batchSec = d('batchSec', inf);
maxAtt   = d('maxAtt', 2);
logFile  = d('logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

Q = load(ladderMat);
bak = strrep(ladderMat, '.mat', '_pre_densify.mat');
if ~isfile(bak), copyfile(ladderMat, bak); end
if ~isfield(Q, 'ATTD'), Q.ATTD = zeros(size(Q.OK)); end

ob = Q.meta;  muStar = ob.muStar;  lStar = ob.lStar;  tStar = ob.tStar;
g0  = 9.80665*tStar^2/(1000*lStar);
cnd = (ob.ispS/tStar)*g0;
ndT = @(TN) (TN/ob.m0kg)*tStar^2/(lStar*1000);

[~, rvD0] = pumpkynPie.cr3bp.getDRO(ob.tauDRO);
rvD0 = pumpkyn.cr3bp.cont_np(rvD0, ob.tauDRO, muStar, 1e-12);
[tD, rvD] = pumpkyn.cr3bp.prop(ob.tauDRO, rvD0, muStar);
[~, rvT0] = pumpkyn.cr3bp.getTulip(ob.tauTulip, ob.NpTulip, ob.pmTulip);
rvT0 = pumpkyn.cr3bp.cont_np(rvT0, ob.tauTulip, muStar, 1e-12);
[tT, rvT] = pumpkyn.cr3bp.prop(ob.tauTulip, rvT0, muStar);
dep = @(f) interp1(tD, rvD, mod(f,1)*tD(end), 'spline');
arr = @(f) interp1(tT, rvT, mod(f,1)*tT(end), 'spline');

[nD, nA, nR] = size(Q.OK);
if isempty(rungSel), rungSel = Q.rungs; end

%% Gap list: solves-somewhere pairs missing this rung, easiest rungs first
%  (descending thrust: the high-thrust interior gaps fill in seconds).
gaps = [];
for kr = 1:nR
    if ~any(abs(Q.rungs(kr) - rungSel) < 1e-9), continue, end
    for iD = 1:nD
        for iA = 1:nA
            if any(Q.OK(iD,iA,:)) && ~Q.OK(iD,iA,kr) ...
                    && Q.ATTD(iD,iA,kr) < maxAtt
                gaps(end+1,:) = [iD iA kr]; %#ok<AGROW>
            end
        end
    end
end
lg('DENSIFY: %d gaps among solving pairs (rungs [%s] N)', ...
   size(gaps,1), num2str(rungSel));

tAll = tic;  nFill = 0;
for kg = 1:min(size(gaps,1), maxCells)
    if toc(tAll) > batchSec
        lg('  [batch budget reached -- clean exit after %d gaps]', kg-1);
        break
    end
    iD = gaps(kg,1);  iA = gaps(kg,2);  kr = gaps(kg,3);
    TN = Q.rungs(kr);  Tnd = ndT(TN);
    Q.ATTD(iD,iA,kr) = Q.ATTD(iD,iA,kr) + 1;     % record BEFORE solving
    save(ladderMat, '-struct', 'Q');
    rv0 = dep(Q.sD(iD));  rvf = arr(Q.sA(iA));
    t0 = tic;  ok = false;  route = '-';

    %% Route A: same cell, nearest converged rung in log-thrust
    kSrc = find(squeeze(Q.OK(iD,iA,:)));
    if ~isempty(kSrc)
        [~, ks] = min(abs(log(Q.rungs(kSrc)) - log(TN)));
        kSrc = kSrc(ks);
        [ok, z, res, fly, dz] = fill_from(Q.Z8(:,iD,iA,kSrc), ...
            Q.rungs(kSrc), TN, rv0, rvf, rv0, tfExpList, K, ndT, cnd, ...
            muStar, lStar, msWallS, gateKm, accTol);
        if ok, route = sprintf('A(%.2gN)', Q.rungs(kSrc)); end
    end

    %% Route B: nearest torus neighbor converged at THIS rung
    if ~ok
        [nbD, nbA] = find(Q.OK(:,:,kr));
        if ~isempty(nbD)
            dtor = min(abs(nbD-iD),nD-abs(nbD-iD)).^2 ...
                 + min(abs(nbA-iA),nA-abs(nbA-iA)).^2;
            [~, q] = sort(dtor);
            for kn = q(1:min(2,end))'
                zN = Q.Z8(:,nbD(kn),nbA(kn),kr);
                rvN0 = dep(Q.sD(nbD(kn)));
                [ok, z, res, fly, dz] = fill_from(zN, TN, TN, rv0, rvf, ...
                    rvN0, 1, K, ndT, cnd, muStar, lStar, msWallS, ...
                    gateKm, accTol);
                if ok, route = sprintf('B(%d,%d)', nbD(kn), nbA(kn)); break, end
            end
        end
    end

    if ok
        nFill = nFill + 1;
        Q.TF(iD,iA,kr) = z(8);  Q.RES(iD,iA,kr) = res;
        Q.FLYKM(iD,iA,kr) = fly;  Q.ACCDZ(iD,iA,kr) = dz;
        Q.WALL(iD,iA,kr) = toc(t0);  Q.OK(iD,iA,kr) = true;
        Q.Z8(:,iD,iA,kr) = z;
        save(ladderMat, '-struct', 'Q');            % EVERY fill
        lg('  (%2d,%2d) T=%5.2f  tf=%.5f (%.3f d)  fly=%.2fkm  acc=%.1e  %s  OK (%.0fs) [%d/%d]', ...
           iD, iA, TN, z(8), z(8)*tStar/86400, fly, dz, route, toc(t0), ...
           kg, size(gaps,1));
    else
        lg('  (%2d,%2d) T=%5.2f  unfilled (%.0fs) [%d/%d]', ...
           iD, iA, TN, toc(t0), kg, size(gaps,1));
    end
end
P = Q;
lg('DENSIFY DONE: %d of %d gaps filled, %.1f min', ...
   nFill, min(size(gaps,1),maxCells), toc(tAll)/60);
end

% ------------------------------------------------------------------------
function [ok, z, res, fly, dz] = fill_from(zSrc, Tsrc, TN, rv0, rvf, ...
    rvSrc0, tfExpList, K, ndT, cnd, muStar, lStar, msWallS, gateKm, accTol)
% FILL_FROM  Seed ms_tfmin from a source solution (same cell at another
% rung, or a neighbor at the same rung) and run the acceptance gates.
% INPUTS:  zSrc [8x1] converged source; Tsrc its thrust; TN target thrust;
%          rv0/rvf THIS cell's endpoints; rvSrc0 the SOURCE's departure
%          state (differs on route B); stretch exponents; ms settings.
% OUTPUTS: ok; z [8x1]; ms residual; flown miss (km); tfMin |dz|.
z = zSrc;  res = NaN;  fly = inf;  dz = NaN;  ok = false;
try
    [tj, yj] = pumpkyn.cr3bp.tfMinProp(zSrc(8), ...
        [rvSrc0(1:6)'; 1; zSrc(1:7)], ndT(Tsrc), cnd, muStar);
    [tu, iu] = unique(tj);
    sG = linspace(0, 1, K+1);
    for p = tfExpList
        tfG = zSrc(8) * (Tsrc/TN)^p;
        Yg = interp1(tu/tu(end), yj(iu,1:14), sG, 'pchip')';
        Yg(7,:) = 1 + (Yg(7,:)-1)*(tfG/zSrc(8))*(TN/Tsrc);   % mdot ~ T (review)
        seed = struct('tf', tfG, 'tGrid', sG*tfG, 'Y', Yg);
        [zt, it] = ms_tfmin(rv0(1:6), rvf(1:6), seed, ndT(TN), cnd, ...
                            muStar, struct('wallSec', msWallS));
        [~, rvFly] = pumpkyn.cr3bp.tfMinProp(zt(8), ...
            [rv0(1:6)'; 1; zt(1:7)], ndT(TN), cnd, muStar);
        mk = norm(rvFly(end,1:3) - rvf(1:3))*lStar;
        if ~(it.converged && mk < gateKm), continue, end
        evalc('zC = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), zt, ndT(TN), cnd, muStar);');
        ad = norm(zC - zt);
        if ad >= accTol, continue, end
        if ~ok || zt(8) < z(8)               % keep the FASTEST valid fill
            ok = true;  z = zt;  res = it.normR;  fly = mk;  dz = ad;
        end
    end
catch
end
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end

function logmsg(f, s)
% LOGMSG  Append a line to the log file (or stdout if none).
% INPUTS: f path ('' = stdout); s text.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f,'a'); fprintf(fid,'%s\n',s); fclose(fid);
end
end
