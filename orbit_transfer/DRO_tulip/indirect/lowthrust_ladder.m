function P = lowthrust_ladder(outMat, opts)
%% Purpose:
%
%   Builds the LOW-THRUST sheet of the DRO -> tulip costate library by
%   walking thrust UPWARD from the converged 0.07 N solutions (the v1
%   library) toward the bottom of the high-thrust sheet at 0.5 N. This is
%   the "meet in the middle" strategy: the downward ladder from 15 N stalls
%   near one DRO revolution, but the 0.07 N regime is already solved -- so
%   climb out of it instead of descending into it.
%
%   Per phase pair, the steps are:
%
%     0. RE-ANCHOR   the v1 entry (0.07 N, Isp 900 s) is re-converged at the
%                    target Isp by indirect continuation -- measured to be a
%                    trivial walk (t_f moves ~0.4%).
%
%     1..R  CLIMB    for each thrust rung, the previous converged solution
%                    is propagated and becomes a DIRECT-solve warm start
%                    (states from the trajectory, control from the PMP law).
%                    Direct transcription deforms benignly through
%                    revolution-count changes where costate continuation
%                    diverges -- this is the lesson of the v2 ladder, now
%                    pointed upward. Each rung then runs the full pipeline:
%                    direct solve -> flown verification -> ms_tfmin
%                    refinement -> tfMin acceptance. All three gates, every
%                    rung.
%
%     JOIN           at the top rung, if the high-thrust sheet holds an
%                    entry for this phase pair, the two flight times are
%                    compared: agreement means one continuous solution
%                    family from 0.07 N to 15 N; disagreement is recorded as
%                    a family split (a finding, not a failure).
%
%   The mesh (N, K, Sundman on/off) comes from mesh_policy, scaled by the
%   expected transfer length in revolutions. In 'auto' mode a rung that
%   fails is retried once with the opposite Sundman setting, and the setting
%   that worked is recorded.
%
%% Inputs:
%
%  outMat                   char                    Output .mat, written
%                                                   after EVERY rung
%
%  opts                     struct                  (all optional)
%   .rungs                  [1 x R]                 Upward thrust rungs, N,
%                                                   ASCENDING
%                                                   [0.1 0.15 0.2 0.3 0.4 0.5]
%   .ispS                   double                  Target Isp, s     [1710]
%   .m0kg                   double                  Initial mass, kg  [150]
%   .meshMode               char                    'auto' | 'uniform' |
%                                                   'sundman'         [auto]
%   .gateKm                 double                  Flown-arrival gate, km
%                                                                     [100]
%   .accTol                 double                  tfMin acceptance  [1e-6]
%   .msWallS                double                  ms_tfmin budget/attempt,
%                                                   s                 [180]
%   .maxIter                double                  NLP iteration cap [3000]
%   .cpuSec                 double                  IPOPT wall cap, s [400]
%   .v1Mat                  char                    v1 refined catalog
%                                                   (ms_refined_12x12.mat)
%   .sheetAMat              char                    high-thrust ladder
%                                                   (thrust_ladder_12x12.mat)
%   .cells                  [k x 2]                 (iD,iA) subset [all v1]
%   .maxCells, .batchSec    double                  Batch control     [inf]
%   .maxAtt                 double                  Attempts/cell     [2]
%   .logFile                char                    Append-mode log [stdout]
%
%% Outputs:
%
%  P                        struct                  Arrays [nD x nA x R+1]
%                                                   (rung 1 = the 0.07 N
%                                                   re-anchor): TF, FLYKM,
%                                                   ACCDZ, RES, WALL, OK,
%                                                   SUNDUSED; Z8
%                                                   [8 x nD x nA x R+1];
%                                                   ATT, JOINDTF, rungs,
%                                                   sD, sA, meta
%
%% Revision History:
%  M. Casey                                                   (c) 08/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
d = @(f,v) fdef(opts, f, v);
rungs    = d('rungs', [0.1 0.15 0.2 0.3 0.4 0.5]);
ispS     = d('ispS', 1710);
m0kg     = d('m0kg', 150);
meshMode = d('meshMode', 'auto');
gateKm   = d('gateKm', 100);
accTol   = d('accTol', 1e-6);
msWallS  = d('msWallS', 180);
maxIter  = d('maxIter', 3000);
cpuSec   = d('cpuSec', 400);
maxCells = d('maxCells', inf);
batchSec = d('batchSec', inf);
maxAtt   = d('maxAtt', 2);
logFile  = d('logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, '..', 'direct', 'results');
v1Mat     = d('v1Mat',     fullfile(resDir, 'ms_refined_12x12.mat'));
sheetAMat = d('sheetAMat', fullfile(resDir, 'thrust_ladder_12x12.mat'));

assert(all(diff(rungs) > 0), 'rungs must ASCEND (this ladder climbs)');
assert(isfile(v1Mat), 'v1 catalog not found: %s', v1Mat);

%% The v1 catalog supplies the anchors; sheet A supplies the join test:
V1 = load(v1Mat);
ob = V1.meta.orbit;
muStar = ob.muStar;  lStar = ob.lStar;  tStar = ob.tStar;
g0   = 9.80665*tStar^2/(1000*lStar);
c900 = (V1.meta.ispS/tStar)*g0;                  % v1's exhaust velocity
cnd  = (ispS/tStar)*g0;                          % target exhaust velocity
ndT  = @(TN) (TN/m0kg)*tStar^2/(lStar*1000);
T0   = V1.meta.thrustN;                          % 0.07 N, the anchor thrust
haveA = isfile(sheetAMat);
if haveA, SA = load(sheetAMat); end

[~, rvD0] = pumpkynPie.cr3bp.getDRO(ob.tauDRO);
rvD0 = pumpkyn.cr3bp.cont_np(rvD0, ob.tauDRO, muStar, 1e-12);
[tD, rvD] = pumpkyn.cr3bp.prop(ob.tauDRO, rvD0, muStar);
[~, rvT0] = pumpkyn.cr3bp.getTulip(ob.tauTulip, ob.NpTulip, ob.pmTulip);
rvT0 = pumpkyn.cr3bp.cont_np(rvT0, ob.tauTulip, muStar, 1e-12);
[tT, rvT] = pumpkyn.cr3bp.prop(ob.tauTulip, rvT0, muStar);

[nD, nA] = size(V1.OKI);
rungsAll = [T0, rungs(:)'];                      % rung 1 = the re-anchor
nR = numel(rungsAll);

%% Result arrays (or resume from a previous run of this same campaign):
if isfile(outMat)
    Q = load(outMat);
    assert(isequal(Q.rungs, rungsAll), 'existing %s has different rungs', outMat);
else
    z = nan(nD,nA,nR);
    Q = struct('TF',z, 'FLYKM',z, 'ACCDZ',z, 'RES',z, 'WALL',z, ...
               'SUNDUSED',z, 'OK',false(nD,nA,nR), 'Z8',nan(8,nD,nA,nR), ...
               'ATT',zeros(nD,nA), 'JOINDTF',nan(nD,nA), ...
               'rungs',rungsAll, 'sD',V1.sD, 'sA',V1.sA);
    Q.meta = struct('orbit',ob, 'ispS',ispS, 'm0kg',m0kg, ...
                    'meshMode',meshMode, 'gateKm',gateKm, ...
                    'source_v1',v1Mat, 'source_sheetA',sheetAMat);
end

%% Work list: v1 cells not yet finished and not out of attempts:
only = d('cells', []);                % optional subset, e.g. the pilot cells
todo = [];
for iD = 1:nD
    for iA = 1:nA
        if ~isempty(only) && ~any(only(:,1)==iD & only(:,2)==iA), continue, end
        unfinished = V1.OKI(iD,iA) && ~Q.OK(iD,iA,nR) && Q.ATT(iD,iA) < maxAtt;
        if unfinished, todo(end+1,:) = [iD iA]; end %#ok<AGROW>
    end
end
lg('LOW-THRUST LADDER (%s mesh): %d cells, climbing [%s] N from %.3g N', ...
   meshMode, size(todo,1), num2str(rungs), T0);

tAll = tic;
for kc = 1:min(size(todo,1), maxCells)
    if toc(tAll) > batchSec
        lg('  [batch budget reached -- clean exit after %d cells]', kc-1);
        break
    end
    iD = todo(kc,1);  iA = todo(kc,2);
    Q.ATT(iD,iA) = Q.ATT(iD,iA) + 1;             % record BEFORE solving
    save(outMat, '-struct', 'Q');
    rv0 = interp1(tD, rvD, mod(Q.sD(iD),1)*tD(end), 'spline');
    rvf = interp1(tT, rvT, mod(Q.sA(iA),1)*tT(end), 'spline');

    %% Resume point: the highest rung this cell has already converged
    kStart = find(Q.OK(iD,iA,:), 1, 'last');
    if isempty(kStart)
        % ---- rung 1: Isp re-anchor by indirect continuation --------------
        t0 = tic;
        zPrev = V1.LAMI(:,iD,iA);
        [okA, zA, resA, flyA, dzA] = reanchor(zPrev, rv0, rvf, ndT(T0), ...
            c900, cnd, muStar, lStar, msWallS, gateKm, accTol);
        Q.TF(iD,iA,1) = zA(8);  Q.RES(iD,iA,1) = resA;
        Q.FLYKM(iD,iA,1) = flyA;  Q.ACCDZ(iD,iA,1) = dzA;
        Q.WALL(iD,iA,1) = toc(t0);  Q.OK(iD,iA,1) = okA;
        if okA, Q.Z8(:,iD,iA,1) = zA; end
        save(outMat, '-struct', 'Q');
        lg('  (%2d,%2d) reanchor %.3gN Isp%d  tf=%.4f  fly=%.2fkm  %s (%.0fs)', ...
           iD, iA, T0, ispS, zA(8), flyA, okstr(okA), toc(t0));
        if ~okA, continue, end
        kStart = 1;
    end
    zPrev = Q.Z8(:,iD,iA,kStart);
    Tprev = rungsAll(kStart);

    %% Climb the rungs above the resume point:
    for kr = kStart+1 : nR
        TN = rungsAll(kr);
        t0 = tic;
        tfGuess = zPrev(8) * (Tprev/TN)^0.7;     % t_f shrinks going up
        [N, K, sund1] = mesh_policy(tfGuess, meshMode);
        okR = false;
        for sund = unique([sund1, ~sund1 & strcmpi(meshMode,'auto')], 'stable')
            % (in forced modes this loop runs once; in auto it can retry
            %  once with the opposite mesh)
            if ~strcmpi(meshMode,'auto') && sund ~= sund1, break, end
            [okR, z, res, fly, dz] = climb_rung(zPrev, Tprev, TN, tfGuess, ...
                rv0, rvf, N, K, sund, ndT, cnd, muStar, lStar, tStar, ...
                maxIter, cpuSec, msWallS, gateKm, accTol);
            if okR, break, end
        end
        Q.TF(iD,iA,kr) = z(8);  Q.RES(iD,iA,kr) = res;
        Q.FLYKM(iD,iA,kr) = fly;  Q.ACCDZ(iD,iA,kr) = dz;
        Q.WALL(iD,iA,kr) = toc(t0);  Q.OK(iD,iA,kr) = okR;
        Q.SUNDUSED(iD,iA,kr) = double(sund);
        if okR, Q.Z8(:,iD,iA,kr) = z; end
        save(outMat, '-struct', 'Q');
        lg('  (%2d,%2d) T=%.2fN  tf=%.4f (%5.2f d, %.2f rev)  fly=%.2fkm  acc=%.1e  sund=%d  %s (%.0fs)', ...
           iD, iA, TN, z(8), z(8)*tStar/86400, z(8), fly, dz, sund, ...
           okstr(okR), toc(t0));
        if ~okR, break, end                      % ladder stops for this cell
        zPrev = z;  Tprev = TN;
    end

    %% Join test against the high-thrust sheet at the top rung:
    if Q.OK(iD,iA,nR) && haveA
        krA = find(abs(SA.rungs - rungsAll(nR)) < 1e-9, 1);
        if ~isempty(krA) && SA.OK(iD,iA,krA)
            Q.JOINDTF(iD,iA) = abs(Q.TF(iD,iA,nR) - SA.TF(iD,iA,krA));
            lg('  (%2d,%2d) JOIN at %.2fN: |tf_up - tf_down| = %.2e ND (%s)', ...
               iD, iA, rungsAll(nR), Q.JOINDTF(iD,iA), ...
               famstr(Q.JOINDTF(iD,iA)));
            save(outMat, '-struct', 'Q');
        end
    end
end
P = Q;
lg('LOW-THRUST LADDER DONE: %d verified entries, %d cells at the top rung, %.1f min', ...
   nnz(Q.OK), nnz(Q.OK(:,:,nR)), toc(tAll)/60);
end

% ------------------------------------------------------------------------
function [ok, z, res, fly, dz] = reanchor(z0, rv0, rvf, Tnd, cOld, cNew, ...
    muStar, lStar, wallSec, gateKm, accTol)
% REANCHOR  Re-converge a v1 entry at the target Isp (indirect continuation).
% INPUTS:  z0 [8x1] converged at cOld; endpoints; ND thrust; old/new exhaust
%          velocities; gates.  OUTPUTS: ok; z [8x1]; ms residual; flown miss
%          (km); tfMin |dz|.
z = z0;  res = NaN;  fly = inf;  dz = NaN;  ok = false;
try
    [tj, yj] = pumpkyn.cr3bp.tfMinProp(z0(8), [rv0(1:6)'; 1; z0(1:7)], ...
                                       Tnd, cOld, muStar);
    [tu, iu] = unique(tj);
    K = 24;
    sGrid = linspace(0, 1, K+1);
    Yg = interp1(tu/tu(end), yj(iu,1:14), sGrid, 'pchip')';
    seed = struct('tf', z0(8), 'tGrid', sGrid*z0(8), 'Y', Yg);
    [z, info] = ms_tfmin(rv0(1:6), rvf(1:6), seed, Tnd, cNew, muStar, ...
                         struct('wallSec', wallSec));
    res = info.normR;
    [~, rvFly] = pumpkyn.cr3bp.tfMinProp(z(8), [rv0(1:6)'; 1; z(1:7)], ...
                                         Tnd, cNew, muStar);
    fly = norm(rvFly(end,1:3) - rvf(1:3))*lStar;
    if info.converged && fly < gateKm
        evalc('zC = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), z, Tnd, cNew, muStar);');
        dz = norm(zC - z);
        ok = dz < accTol;
    end
catch
end
end

% ------------------------------------------------------------------------
function [ok, z, res, fly, dz] = climb_rung(zPrev, Tprev, TN, tfGuess, ...
    rv0, rvf, N, K, useSundman, ndT, cnd, muStar, lStar, tStar, ...
    maxIter, cpuSec, msWallS, gateKm, accTol)
% CLIMB_RUNG  One upward rung: direct warm start -> refine -> accept.
% INPUTS:  previous converged z and its thrust; the new rung; mesh settings;
%          physics; gates.  OUTPUTS: ok; z [8x1]; ms residual; flown miss
%          (km); tfMin |dz|.
z = zPrev;  res = NaN;  fly = inf;  dz = NaN;  ok = false;
try
    % direct warm start manufactured from the previous rung's solution:
    % states from the propagated trajectory, control from the PMP law
    [tj, yj] = pumpkyn.cr3bp.tfMinProp(zPrev(8), ...
        [rv0(1:6)'; 1; zPrev(1:7)], ndT(Tprev), cnd, muStar);
    [tu, iu] = unique(tj);
    sN = linspace(0, 1, N+1);
    Yn = interp1(tu/tu(end), yj(iu,1:14), sN, 'pchip')';
    X0 = Yn(1:7,:);
    lv = Yn(11:13,:);
    U0 = [-lv ./ max(vecnorm(lv,2,1), eps); ones(1, N+1)];
    o = casadi_mintime_dro(rv0(1:6), rvf(1:6), ndT(TN), cnd, muStar, N, ...
            X0, U0, tfGuess, struct('maxIter',maxIter, ...
            'scheme','hermite-simpson', 'sundman',useSundman, ...
            'returnModel',true, 'minAltKm',500, 'maxCpuSec',cpuSec));
    if ~(o.success && o.maxDefect < 1e-9), return, end
    % pre-flight screen before any unbounded integration
    rMoonKm = 1737.4;
    dMoon = vecnorm(o.X(1:3,:) - [1-muStar;0;0], 2, 1)*lStar - rMoonKm;
    if min(dMoon) < 250 || o.tf > 3*tfGuess || o.tf < 0.3*tfGuess, return, end
    C = certify_dro_mintime(o, struct('muStar',muStar,'lStar',lStar,'tStar',tStar), ...
            ndT(TN), cnd, struct('tfRef',[],'verbose',false,'posTolKm',inf));
    if C.globKm > gateKm, fly = C.globKm; return, end
    % refine: sign-vote the duals, seed multiple shooting at the midpoints
    lamD = o.lamDef(1:7,:);
    nv = min(10, size(lamD,2));  vote = 0;
    for kv = 1:nv
        aS = o.Um(1:3,kv)/max(norm(o.Um(1:3,kv)),eps);
        pv = -lamD(4:6,kv)/max(norm(lamD(4:6,kv)),eps);
        vote = vote + sign(pv.'*aS);
    end
    sg = 1;  if vote < 0, sg = -1; end
    if isempty(o.tNodes), tN_ = linspace(0, o.tf, size(o.X,2));
    else,                 tN_ = o.tNodes(:)';  end
    tMid = tN_(1:end-1) + diff(tN_)/2;
    for Kk = [K, min(2*K, 80)]
        tG = linspace(0, o.tf, Kk+1);
        Xg = interp1(tN_, o.X', tG, 'pchip')';
        Lg = interp1(tMid, (sg*lamD)', tG, 'pchip', 'extrap')';
        seed = struct('tf', o.tf, 'tGrid', tG, 'Y', [Xg; Lg]);
        [z, info] = ms_tfmin(rv0(1:6), rvf(1:6), seed, ndT(TN), cnd, ...
                             muStar, struct('wallSec', msWallS));
        res = info.normR;
        [~, rvFly] = pumpkyn.cr3bp.tfMinProp(z(8), [rv0(1:6)'; 1; z(1:7)], ...
                                             ndT(TN), cnd, muStar);
        fly = norm(rvFly(end,1:3) - rvf(1:3))*lStar;
        if info.converged && fly < gateKm
            evalc('zC = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), z, ndT(TN), cnd, muStar);');
            dz = norm(zC - z);
            if dz < accTol, ok = true; return, end
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

function s = okstr(t)
% OKSTR  OK/fail marker.  INPUTS: t logical.  OUTPUTS: s [char].
if t, s = 'OK'; else, s = 'fail'; end
end

function s = famstr(dtf)
% FAMSTR  Same-family verdict from a join |dtf|.  INPUTS: dtf.  OUTPUTS: s.
if dtf < 1e-5, s = 'SAME family'; else, s = 'DIFFERENT family'; end
end
