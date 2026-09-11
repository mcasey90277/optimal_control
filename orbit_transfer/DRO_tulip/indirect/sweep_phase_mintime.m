function S = sweep_phase_mintime(opts)
%% DEPRECATED (2026-09-09) -- superseded, kept for provenance only.
%
%   The fixed-step nearest-neighbour traversal in this file was judged the
%   WRONG METHOD for the phase sheet (FINDINGS 36): it steps OVER folds
%   instead of walking through them, so it stalls where the solution curve
%   turns, and six defects had to be fixed in it before it produced anything
%   trustworthy. It is superseded by the continuation route --
%   arclength_arrival on costate_common/arclength_ms for the arrival axis,
%   rib_from_crossing for the departure axis, both through certify_root.
%
%   It is kept because `results/sweep_phase_mintime.mat` is still a seed
%   source for dro_tulip_library, and this file is the record of how those
%   points were produced. DO NOT extend it; extend the continuation route.
%
%% Purpose:
%
%   The (departure phase x arrival phase) SHEET of the minimum-time
%   DRO -> tulip transfer at ONE thrust: breadth-first continuation from a
%   certified anchor, every point solved by MULTIPLE shooting and gated by
%   the flown arrival, the independent tfMin re-solve and the conjugate test.
%
%   Why a new sweep (sweep_phasing exists): that one solves each point with
%   SINGLE shooting and jumps a whole grid step at once. On the 6x6 grid an
%   arrival step is 1/6 of the tulip period = 3.8 days of phase, and its
%   stored run closed 1 of 36 points -- the BFS died on the first neighbour.
%   Everything learned since (FINDINGS 32-35) says the fix is multiple
%   shooting plus CONTINUATION IN PHASE: walk to a neighbour in small
%   increments, each seeded by the last, bisecting on failure.
%
%   The measurement this exists for: with orbits, engine and departure phase
%   fixed, arrival phase alone moves t_f 17.80 -> 26.44 d (+48%) at 70 mN
%   (FINDINGS 35). Arrival phase, not thrust, sets the cost of this
%   transfer, and for a constellation released into distinct phase slots
%   this sheet IS the deployment envelope.
%
%% Inputs:
%
%  opts                     struct (optional)
%   .nD, .nA                int                     grid sizes [12, 12]
%   .thrustN                double                  [0.070] N
%   .ispS, .m0kg            double                  [900] s, [150] kg
%   .tauDRO, .NpTulip       double                  [1.0], [7]
%   .K                      int                     ms segments [24]
%   .nSub                   int                     phase sub-steps per grid
%                                                   edge [4]
%   .maxBisect              int                     sub-step halvings [3]
%   .stepKm                 double                  target endpoint motion
%                                                   per sub-step, km [2000]
%                                                   (measured convergence
%                                                   limit ~ a few 1000 km)
%   .ptSec                  double                  clock backstop per grid
%                                                   point, seconds [2400]
%   .maxSolve               int                     hard sub-solve budget per
%                                                   grid point [800]
%   .gateKm                 double                  flown arrival gate [100]
%   .tolDz                  double                  tfMin acceptance [1e-6]
%   .conj                   logical                 run the conjugate test
%                                                   at each point [true]
%   .outMat, .logFile       char
%   .wallSec                double                  per solve [120]
%
%% Outputs:
%
%  S                        struct                  .sD .sA grids, .TF
%                                                   [nD x nA] ND, .DV, .MF,
%                                                   .Z8 [8 x nD x nA],
%                                                   .Yj {nD x nA} junctions,
%                                                   .FLYKM .ACCDZ .CONJ
%                                                   (1 pass / 0 fail / -1
%                                                   not run) .NSUB .meta
%
%% Revision History:
%  M. Casey                                                   (c) 09/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
nD = d('nD',12); nA = d('nA',12); thrustN = d('thrustN',0.070);
ispS = d('ispS',900); m0kg = d('m0kg',150);
tauDRO = d('tauDRO',1.0); NpTulip = d('NpTulip',7);
K = d('K',24); nSub = d('nSub',4); maxBis = d('maxBisect',3);
gateKm = d('gateKm',100); tolDz = d('tolDz',1e-6); doConj = d('conj',true);
wallSec = d('wallSec',120);
stepKm  = d('stepKm', 2000);      % target endpoint motion per sub-step, km
ptSec   = d('ptSec', 2400);       % clock backstop per grid point, seconds
maxSolve= d('maxSolve', 800);     % HARD sub-solve budget per grid point
here = fileparts(mfilename('fullpath'));
outMat = d('outMat', fullfile(here,'results','sweep_phase_mintime.mat'));
logFile = d('logFile','');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));
addpath(fullfile(fileparts(here),'..','costate_common'));

ob = struct('muStar',0.012150585609624,'lStar',389703.264829278, ...
            'tStar',382981.289129055,'tauDRO',tauDRO,'NpTulip',NpTulip, ...
            'tauTulip',5*2*pi/6,'pmTulip',-1,'ispS',ispS,'m0kg',m0kg);
lStar = ob.lStar; tStar = ob.tStar; mu = ob.muStar;
% THE shared propulsion conversion (costate_common/nd_propulsion)
ndp = nd_propulsion(thrustN, ispS, m0kg, lStar, tStar);
cnd = ndp.cnd;   Tnd = ndp.Tnd;
[tD, rvD, tT, rvT] = ladder_endpoints(ob);
% THE shared endpoint rule (costate_common/phase_state, FINDINGS 44)
stateD = phase_state(tD, rvD);
stateA = phase_state(tT, rvT);

%% Anchor: the certified demo phasing pair -------------------------------
Aanc = load(fullfile(here,'results','mintime_70mN_anchor.mat'));
sD0 = 0;  sA0 = 0.0754;
sD = mod(sD0 + (0:nD-1)/nD, 1);
sA = mod(sA0 + (0:nA-1)/nA, 1);
S = struct('sD',sD,'sA',sA,'TF',nan(nD,nA),'DV',nan(nD,nA),'MF',nan(nD,nA), ...
    'Z8',nan(8,nD,nA),'Yj',{cell(nD,nA)},'FLYKM',nan(nD,nA),'ACCDZ',nan(nD,nA), ...
    'CONJ',-ones(nD,nA),'NSUB',nan(nD,nA),'NCALL',nan(nD,nA), ...
    'meta',struct('thrustN',thrustN,'ispS',ispS,'m0kg',m0kg,'K',K,'ob',ob, ...
                  'Tnd',Tnd,'cnd',cnd,'anchor',[sD0 sA0]));
z0 = Aanc.z(:);
Y0 = Aanc.it.Y;
% The anchor's diagnostics are MEASURED here, never fabricated: the file may
% have been certified under different physics than this call requests.
% (Astra review 2026-09-09, defect 2.)
assert(abs(Aanc.Tnd - Tnd)/Tnd < 1e-10 && abs(Aanc.cnd - cnd)/cnd < 1e-10, ...
    ['anchor was certified at a different operating point (Tnd %.6g vs %.6g, ' ...
     'cnd %.6g vs %.6g) -- re-solve it before sweeping'], Aanc.Tnd, Tnd, Aanc.cnd, cnd);
[okA, z0, Y0, flyA, dzA, cjA] = solvePoint(z0, Y0, stateD(sD(1)), stateA(sA(1)), ...
    Tnd, cnd, mu, lStar, K, gateKm, tolDz, doConj, wallSec);
assert(okA, 'anchor did not pass its own gates at this operating point');
S = record(S,1,1,z0,Y0,flyA,dzA,cjA,stateD(sD(1)),stateA(sA(1)),Tnd,cnd,mu,lStar,tStar);
lg('sweep %dx%d at %.1f mN, Isp %g: anchor (%.4f, %.4f) tf = %.4f d', ...
   nD, nA, thrustN*1000, ispS, sD(1), sA(1), z0(8)*tStar/86400);

%% SPINE AND RIBS ---------------------------------------------------------
%  The two directions cost wildly different amounts. One grid step in
%  ARRIVAL phase moves the target ~40,900 km (the tulip is big) and needs
%  ~20 sub-solves; the same fraction of the DRO moves the departure point a
%  fraction of that and needs ~4. Breadth-first walks expensive arrival
%  edges from EVERY departure position -- ~12x the necessary work. So:
%  walk the arrival SPINE once at the anchor's departure phase, then walk
%  the cheap departure RIBS out of each spine point.
nSolved = 1;  tAll = tic;
walk = @(z, Y, a0, b0, a1, b1) walkPhase(z, Y, a0, b0, a1, b1, ...
    nD, nA, nSub, maxBis, stateD, stateA, Tnd, cnd, mu, lStar, ...
    K, gateKm, tolDz, doConj, wallSec, stepKm, ptSec, maxSolve);
put = @(S,i,j,z,Y,fly,dz,cj,used) recordFull(S,i,j,z,Y,fly,dz,cj,used, ...
    stateD(sD(i)),stateA(sA(j)),Tnd,cnd,mu,lStar,tStar);

lg('  -- arrival spine at sD = %.4f --', sD(1));
for j = 2:nA
    [okP, zN, YN, fly, dz, cj, used] = walk(squeeze(S.Z8(:,1,j-1)), S.Yj{1,j-1}, ...
        sD(1), sA(j-1), sD(1), sA(j));
    if ~okP
        lg('  spine (1,%2d) sA=%.4f: FAILED -- spine stops here', j, sA(j));
        break
    end
    S = put(S,1,j,zN,YN,fly,dz,cj,used);  nSolved = nSolved + 1;
    lg('  (%2d,%2d) sD=%.4f sA=%.4f: tf=%7.3f d dV=%.4f fly=%6.2f km dz=%.1e conj=%d [%d/%d, %d sub, %d calls] SPINE', ...
       1, j, sD(1), sA(j), S.TF(1,j)*tStar/86400, S.DV(1,j), fly, dz, cj, nSolved, nD*nA, used(1), used(2));
    save(outMat,'S');
end

for j = 1:nA
    if ~isfinite(S.TF(1,j)), continue, end
    for dir = [-1 1]
        zc = squeeze(S.Z8(:,1,j));  Yc = S.Yj{1,j};  ic = 1;
        for s = 1:floor(nD/2)
            i1 = mod(ic-1+dir, nD)+1;
            if isfinite(S.TF(i1,j)), break, end
            [okP, zN, YN, fly, dz, cj, used] = walk(zc, Yc, sD(ic), sA(j), sD(i1), sA(j));
            if ~okP
                lg('  (%2d,%2d) sD=%.4f sA=%.4f: FAILED from (%d,%d)', i1, j, sD(i1), sA(j), ic, j);
                break
            end
            S = put(S,i1,j,zN,YN,fly,dz,cj,used);  nSolved = nSolved + 1;
            lg('  (%2d,%2d) sD=%.4f sA=%.4f: tf=%7.3f d dV=%.4f fly=%6.2f km dz=%.1e conj=%d [%d/%d, %d sub, %d calls]', ...
               i1, j, sD(i1), sA(j), S.TF(i1,j)*tStar/86400, S.DV(i1,j), fly, dz, cj, nSolved, nD*nA, used(1), used(2));
            zc = zN;  Yc = YN;  ic = i1;
            save(outMat,'S');
        end
    end
end
S.meta.wallMin = toc(tAll)/60;
S.meta.nSolved = nSolved;
save(outMat,'S');
lg('sweep done: %d/%d solved in %.1f min -> %s', nSolved, nD*nA, S.meta.wallMin, outMat);
end

% ------------------------------------------------------------------------
function [ok, z, Y, fly, dz, cj, used] = walkPhase(z0, Y0, sD0, sA0, sD1, sA1, ...
    nD, nA, nSub, maxBis, stateD, stateA, Tnd, cnd, mu, lStar, K, gateKm, tolDz, doConj, wallSec, stepKm, ptSec, maxSolve)
% WALKPHASE  Continue one converged solution from phase (sD0,sA0) to
% (sD1,sA1) in nSub increments, each seeded by the last, halving the
% increment on failure. A whole grid edge in ONE jump is what killed the
% single-shooting sweep.
% INPUTS: see caller.  OUTPUTS: ok; z; Y; fly; dz; cj; used (sub-steps).
ok = false;  z = z0;  Y = Y0;  fly = NaN;  dz = NaN;  cj = -1;
dDs = wrapDiff(sD1, sD0);  dAs = wrapDiff(sA1, sA0);
% SUB-STEP COUNT FROM ENDPOINT MOTION, not a fixed number. Measured
% 2026-09-08 at the anchor: an arrival step converges at 2199 km of target
% motion and fails at 16512 km, so the limit is a few thousand km. One grid
% step in ARRIVAL phase moves the target 40871 km (the tulip is large);
% the same fraction of the DRO moves the departure point far less. That
% geometric asymmetry -- not t_f sensitivity, not seeding -- is why every
% arrival step failed at a fixed nSub = 4 while every departure step passed.
d0_ = stateD(sD0);  d1_ = stateD(sD1);
a0_ = stateA(sA0);  a1_ = stateA(sA1);
moveKm = (norm(d1_(1:3) - d0_(1:3)) + norm(a1_(1:3) - a0_(1:3)))*lStar;
n = max(nSub, ceil(moveKm/stepKm));  bis = 0;
% WALL BUDGET PER GRID POINT. Adaptive sizing (n ~ moveKm/stepKm, often 70+)
% MULTIPLIED by bisection doubling (maxBisect halvings) is a worst case of
% n*2^maxBis sub-solves: with n = 68 and maxBis = 5 that is 2176 solves at
% up to wallSec each = ~36 h on ONE point. Measured the hard way on
% 2026-09-08: the sweep burned 17 h of CPU on a single failed departure
% step and logged nothing. A point that cannot be reached inside its budget
% is a FAILED point, not a reason to keep the machine.
% Budget by SUB-SOLVE COUNT, with the clock only as a backstop. Measured
% costs of SUCCESSFUL points on 2026-09-08: 68-272 sub-solves for departure
% steps, 368 for the one arrival step that closed. A 7-minute clock cut all
% of those off; a count of ~800 covers them with margin and still bounds the
% worst case, which is what the 17-hour hang needed.
tPt = tic;  nSolve = 0;
while bis <= maxBis && toc(tPt) < ptSec && nSolve < maxSolve
    zc = z0;  Yc = Y0;  good = true;
    rvPrev = [stateD(sD0), stateA(sA0)];
    fDone = 0;                       % fraction of the edge already accepted
    for k = 1:n
        f = k/n;
        rv0 = stateD(sD0 + f*dDs);  rvf = stateA(sA0 + f*dAs);
        % TANGENT PREDICTOR. The single-shooting residual is
        % r = [y_f(1:6) - rvf; lam_m(tf); H(tf)], so dr/d(rvf) = [-I; 0; 0]
        % and dr/d(rv0) = B from the shooting kernel. Hence
        %   dz ~ J \ ([drvf; 0; 0] - B*drv0).
        % Without it the neighbour's t_f is the guess, and an ARRIVAL step of
        % one grid cell is 1.9 d of phase that moves t_f by a day or more --
        % which is why every arrival step failed and every departure step
        % (4.4 d orbit, t_f nearly flat) succeeded.
        zc = predictZ(zc, rvPrev, [rv0, rvf], Tnd, cnd, mu);
        [okS, zc, Yc, fly, dz, cj] = solvePoint(zc, Yc, rv0, rvf, Tnd, cnd, mu, ...
            lStar, K, gateKm, tolDz, doConj && k == n, wallSec);
        nSolve = nSolve + 1;
        if ~okS, good = false; break, end
        rvPrev = [rv0, rvf];  fDone = f;
        if k == n, break, end            % edge complete: accept BEFORE testing budget
        if toc(tPt) > ptSec || nSolve >= maxSolve, good = false; break, end
        rvPrev = [rv0, rvf];
    end
    if good, ok = true;  z = zc;  Y = Yc;  used = [n nSolve];  return, end
    n = 2*n;  bis = bis + 1;
end
used = [NaN nSolve];        % keep the call count on failure: it is the budget diagnostic
end

function z = predictZ(z, rvOld, rvNew, Tnd, cnd, mu)
% PREDICTZ  First-order predictor for the costates and t_f when the
% endpoints move: dz = J \ ([drvf; 0; 0] - B*drv0), with J and B from the
% shooting kernel. Returns z unchanged if the kernel or the solve fails --
% a bad predictor must never be worse than no predictor.
% INPUTS: z [8x1]; rvOld, rvNew [6x2] ([rv0, rvf]); Tnd; cnd; mu.
% OUTPUTS: z [8x1].
try
    [~, J, B] = sweep_phasing_shoot(z, rvOld(:,1)', rvOld(:,2)', Tnd, cnd, mu);
    drv0 = rvNew(:,1) - rvOld(:,1);
    drvf = rvNew(:,2) - rvOld(:,2);
    rhs = [drvf; 0; 0] - B*drv0;
    dz  = J \ rhs;
    if all(isfinite(dz)) && norm(dz) < 0.5*max(norm(z), 1), z = z + dz; end
catch
end
end

function [ok, z, Y, fly, dz, cj] = solvePoint(zSeed, YSeed, rv0, rvf, Tnd, cnd, mu, ...
    lStar, K, gateKm, tolDz, doConj, wallSec)
% SOLVEPOINT  One multiple-shooting solve at new endpoints, seeded by the
% neighbour's junctions, then the flown gate, the tfMin witness and (opt)
% the conjugate test.  INPUTS: see caller.  OUTPUTS: ok; z; Y; fly; dz; cj.
ok = false;  z = zSeed;  Y = YSeed;  fly = NaN;  dz = NaN;  cj = -1;
% KEEP THE NEIGHBOUR'S JUNCTIONS. Re-flying the whole trajectory from the
% initial costates (seed_from_z8) reintroduces exactly the long-horizon
% sensitivity multiple shooting exists to remove, and it is not needed: an
% MS seed does not have to satisfy the continuity defects -- removing them
% is the solver's job. This was almost certainly why an arrival grid step
% cost hundreds of sub-solves. (Astra review 2026-09-09, defect 3.)
K_ = size(YSeed, 2);
seed = struct('tf', zSeed(8), 'tGrid', linspace(0, zSeed(8), K_+1), ...
              'Y', [YSeed, YSeed(:,end)]);
seed.Y(8:14,1) = zSeed(1:7);          % predictor's costates on junction 1
seed.Y(1:7,1)  = [rv0(1:6); 1];       % the NEW departure state
o = struct('tolR', 3e-11, 'wallSec', wallSec);
if doConj, o.conjTest = true; end
try
    [zt, it] = ms_tfmin(rv0(1:6), rvf(1:6), seed, Tnd, cnd, mu, o);
catch
    return
end
if ~it.converged, return, end
[~, Yf] = pumpkyn.cr3bp.tfMinProp(zt(8), [rv0(1:6); 1; zt(1:7)], Tnd, cnd, mu);
fly = norm(Yf(end,1:3) - rvf(1:3)')*lStar;
if ~(fly < gateKm), return, end
try
    za = pumpkyn.cr3bp.tfMin(rv0(1:6)', rvf(1:6)', zt(:), Tnd, cnd, mu);
    dz = norm(za(:) - zt(:));
catch
    dz = NaN;
end
if doConj
    if isfield(it, 'conj') && isfield(it.conj, 'pass'), cj = double(it.conj.pass);
    else, cj = -1; end
end
% ENFORCE the gates this function's header promises. Until 2026-09-09 dz and
% cj were COMPUTED and then ignored -- tolDz was never referenced, a tfMin
% exception left dz = NaN and the point was still accepted, and a failed or
% missing conjugate verdict was accepted. The stored points happened to have
% dz = 0 and cj = 1, so the table was true; the harness's promise was not.
% (Astra review 2026-09-09, defect 1.)
witnessOK = isfinite(dz) && dz <= tolDz;
conjOK    = ~doConj || cj == 1;
if ~(witnessOK && conjOK), ok = false; return, end
z = zt;  Y = it.Y;  ok = true;
end

function S = recordFull(S, i, j, z, Y, fly, dz, cj, used, rv0, rvf, Tnd, cnd, mu, lStar, tStar)
% RECORDFULL  record() plus the sub-step count.  INPUTS: as record, + used.
% OUTPUTS: S.
S = record(S, i, j, z, Y, fly, dz, cj, rv0, rvf, Tnd, cnd, mu, lStar, tStar);
S.NSUB(i,j) = used(1);  S.NCALL(i,j) = used(2);
end

function S = record(S, i, j, z, Y, fly, dz, cj, rv0, rvf, Tnd, cnd, mu, lStar, tStar)
% RECORD  Store one solved grid point.  INPUTS: see caller.  OUTPUTS: S.
[~, Yf] = pumpkyn.cr3bp.tfMinProp(z(8), [rv0(1:6); 1; z(1:7)], Tnd, cnd, mu);
mf = Yf(end,7);
S.TF(i,j) = z(8);  S.Z8(:,i,j) = z(:);  S.Yj{i,j} = Y;
S.MF(i,j) = mf;    S.DV(i,j) = cnd*log(1/mf)*lStar/tStar;
S.FLYKM(i,j) = fly;  S.ACCDZ(i,j) = dz;  S.CONJ(i,j) = cj;
end

function d_ = wrapDiff(a, b)
% WRAPDIFF  Shortest signed phase difference a - b on the unit circle.
% INPUTS: a; b.  OUTPUTS: d_.
d_ = mod(a - b + 0.5, 1) - 0.5;
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end

function logmsg(f, s)
% LOGMSG  Append to a log file or stdout.  INPUTS: f; s.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f,'a'); fprintf(fid,'%s\n',s); fclose(fid);
end
end
