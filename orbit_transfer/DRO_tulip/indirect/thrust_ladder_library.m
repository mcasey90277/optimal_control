function P = thrust_ladder_library(outMat, opts)
% THRUST_LADDER_LIBRARY  Build a costate library with a THRUST axis.
%
% For every phase pair on the (departure x arrival) torus, anchor a direct
% minimum-time solve at the HIGH-thrust end -- where the transfer is nearly
% impulsive, sub-revolution, and converges cold in seconds -- then walk the
% thrust DOWN through a rung set, each rung warm-started from the rung above.
% Continuation down the ladder keeps every cell on ONE solution family;
% cold-solving each (cell, thrust) independently does not (measured: cold
% 15 N solves on two meshes land 50% apart in t_f, both feasible).
%
% Each rung runs the full pipeline: direct solve -> continuous verification
% (fly the control) -> ms_tfmin multiple-shooting refinement -> acceptance by
% pumpkyn.cr3bp.tfMin. Results save after EVERY rung.
%
% INPUTS:
%   outMat - output .mat path (rewritten after every rung)
%   opts   - (optional) struct:
%              .rungs     thrust rung set, N, high to low
%                         [15 12 10 7 5 3 2 1.5 1] -- 1.5 added
%                         because t_f(T) steepens below 2 N, so linear
%                         interpolation needs a tighter rung there
%              .ispS      specific impulse, s              [1710]
%              .m0kg      initial mass, kg                 [150]
%              .N         collocation intervals            [400]
%              .floorKm   lunar altitude floor, km         [500]
%              .gateKm    flown-arrival gate, km           [100]
%              .cells     [k x 2] (iD,iA) subset           [all 12x12]
%              .nD,.nA    torus resolution                 [12 12]
%              .sD0,.sA0  phase-grid origins; defaults MATCH the
%                         low-thrust library grid            [0, 0.075378]
%              .maxIter   NLP iteration cap                [3000]
%              .logFile   append-mode log path             [stdout]
%
% OUTPUTS:
%   P - struct with [nD x nA x nRung] arrays .TF (ND), .FLYKM, .ACCDZ,
%       .RES, .WALL, .OK, plus .Z8 [8 x nD x nA x nRung] refined costates,
%       .rungs, .sD, .sA, .meta (orbit + thruster definition).
%
% REFERENCES:
%   [1] casadi_mintime_dro.m -- the direct transcription.
%   [2] ms_tfmin.m -- multiple-shooting refinement.
%   [3] process/COSTATE_LIBRARY_PIPELINE.md -- the three-step process.

if nargin < 2, opts = struct(); end
d = @(f,v) fdef(opts, f, v);
rungs   = d('rungs', [15 12 10 7 5 3 2 1.5 1]);
ispS    = d('ispS', 1710);
m0kg    = d('m0kg', 150);
N       = d('N', 400);
floorKm = d('floorKm', 500);
gateKm  = d('gateKm', 100);
nD      = d('nD', 12);
nA      = d('nA', 12);
% Phase-grid ORIGINS. These must match the low-thrust library exactly or the
% two libraries sample different phase pairs and cannot be merged: that grid
% is offset in arrival phase by the demo's max-velocity-angle anchor.
sD0     = d('sD0', 0);
sA0     = d('sA0', 0.075378);
maxIter = d('maxIter', 3000);
% maxCells: process at most this many fresh cells then return. Lets a shell
% driver run the ladder in small batches under an OS-level timeout, so a
% hang inside an uninterruptible solver call costs minutes, not the run.
maxCells = d('maxCells', inf);
logFile = d('logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

% --- constants and orbits (same definitions as the low-thrust library) ----
muStar = 0.012150585609624;
lStar  = 389703.264829278;
tStar  = 382981.289129055;
tauDRO = 1.0;  NpT = 7;  tauT = 5*2*pi/6;  pmT = -1;
g0  = 9.80665*tStar^2/(1000*lStar);
cnd = (ispS/tStar)*g0;
ndT = @(TN) (TN/m0kg)*tStar^2/(lStar*1000);

[~, rvD0] = pumpkynPie.cr3bp.getDRO(tauDRO);
rvD0 = pumpkyn.cr3bp.cont_np(rvD0, tauDRO, muStar, 1e-12);
[tD, rvD] = pumpkyn.cr3bp.prop(tauDRO, rvD0, muStar);
[~, rvT0] = pumpkyn.cr3bp.getTulip(tauT, NpT, pmT);
rvT0 = pumpkyn.cr3bp.cont_np(rvT0, tauT, muStar, 1e-12);
[tT, rvT] = pumpkyn.cr3bp.prop(tauT, rvT0, muStar);

sD = mod(sD0 + (0:nD-1)/nD, 1);  sA = mod(sA0 + (0:nA-1)/nA, 1);
todo = d('cells', []);
if isempty(todo)
    [gA, gD] = meshgrid(1:nA, 1:nD);
    todo = [gD(:), gA(:)];
end

nR = numel(rungs);
TF = nan(nD,nA,nR); FLYKM = TF; ACCDZ = TF; RES = TF; WALL = TF;
OK = false(nD,nA,nR);  Z8 = nan(8,nD,nA,nR);
% ATT counts attempts per cell, written BEFORE the cell is solved so that a
% cell which HANGS (and is killed externally) still leaves a record. Without
% it a hung cell has no verified rung, resume would not skip it, and a
% batched driver would retry it forever.
ATT = zeros(nD,nA);
maxAtt = d('maxAtt', 2);
% RESUME: a prior run's products are data -- reload them and skip cells that
% already produced verified rungs (opts.resume, default true).
if d('resume', true) && isfile(outMat)
    Q = load(outMat);
    if isequal(size(Q.OK), [nD nA nR])
        TF=Q.TF; FLYKM=Q.FLYKM; ACCDZ=Q.ACCDZ; RES=Q.RES; WALL=Q.WALL;
        OK=Q.OK; Z8=Q.Z8;
        if isfield(Q,'ATT'), ATT = Q.ATT; end
        keep = false(size(todo,1),1);
        for kq = 1:size(todo,1)
            iq = todo(kq,1);  jq = todo(kq,2);
            keep(kq) = ~any(OK(iq,jq,:)) && ATT(iq,jq) < maxAtt;
        end
        nSkip = nnz(~keep);  todo = todo(keep,:);
    else
        nSkip = 0;
    end
else
    nSkip = 0;
end
meta = struct('muStar',muStar,'lStar',lStar,'tStar',tStar, ...
    'tauDRO',tauDRO,'NpTulip',NpT,'tauTulip',tauT,'pmTulip',pmT, ...
    'periodDRO',tD(end),'periodTulip',tT(end), ...
    'ispS',ispS,'m0kg',m0kg,'N',N,'floorKm',floorKm,'gateKm',gateKm, ...
    'sD0',sD0,'sA0',sA0);

lg('THRUST LADDER: %d cells x %d rungs [%s] N, Isp=%d s, m0=%d kg, N=%d (%d cells resumed/skipped)', ...
   size(todo,1), nR, num2str(rungs), ispS, m0kg, N, nSkip);
tAll = tic;
for kc = 1:min(size(todo,1), maxCells)
    iD = todo(kc,1);  iA = todo(kc,2);
    rv0 = interp1(tD, rvD, mod(sD(iD),1)*tD(end), 'spline');
    rvf = interp1(tT, rvT, mod(sA(iA),1)*tT(end), 'spline');
    seedX = [];  seedU = [];  seedTf = [];  Tprev = [];
    ATT(iD,iA) = ATT(iD,iA) + 1;      % record BEFORE solving (hang-proof)
    save(outMat, 'TF','FLYKM','ACCDZ','RES','WALL','OK','Z8','ATT', ...
         'rungs','sD','sA','meta');
    for kr = 1:nR
        TN = rungs(kr);  Tnd = ndT(TN);
        t0 = tic;
        if ~isempty(Tprev)
            % time-of-flight scales roughly as T^-0.6 in this band (pilot)
            seedTf = seedTf * (Tprev/TN)^0.6;
        end
        try
            o = casadi_mintime_dro(rv0(1:6), rvf(1:6), Tnd, cnd, muStar, N, ...
                    seedX, seedU, seedTf, ...
                    struct('maxIter',maxIter,'scheme','hermite-simpson', ...
                           'sundman',false,'returnModel',true, ...
                           'minAltKm',floorKm,'maxCpuSec',300));
        catch ME
            lg('  (%2d,%2d) T=%5.1f  direct ERROR: %s', iD, iA, TN, ME.message);
            break
        end
        if ~(o.success && o.maxDefect < 1e-9)
            lg('  (%2d,%2d) T=%5.1f  direct fail (defect %.1e) -- ladder stops', ...
               iD, iA, TN, o.maxDefect);
            break
        end
        % CHEAP SANITY PRE-CHECK. A solve can meet the discrete defect test
        % to 1e-9 and still be physically wild; integrating such a trajectory
        % crawls near the lunar singularity WITHOUT BOUND (measured: one cell
        % pinned the run 16 min at 100% CPU). Screen on the discrete nodes --
        % altitude floor and a plausible time of flight -- before ever
        % handing the solution to an integrator.
        rMoonKm = 1737.4;
        dMoon = vecnorm(o.X(1:3,:) - [1-muStar;0;0], 2, 1)*lStar - rMoonKm;
        tfBad = ~isempty(seedTf) && (o.tf > 3*seedTf || o.tf < 0.3*seedTf);
        if min(dMoon) < 0.5*floorKm || tfBad
            lg('  (%2d,%2d) T=%5.1f  rejected pre-flight (minAlt %.0f km, tf %.4f) -- ladder stops', ...
               iD, iA, TN, min(dMoon), o.tf);
            break
        end
        C = certify_dro_mintime(o, struct('muStar',muStar,'lStar',lStar,'tStar',tStar), ...
                Tnd, cnd, struct('tfRef',[],'verbose',false,'posTolKm',inf));
        flyKm = C.globKm;
        TF(iD,iA,kr) = o.tf;  FLYKM(iD,iA,kr) = flyKm;  WALL(iD,iA,kr) = toc(t0);
        if flyKm > gateKm
            lg('  (%2d,%2d) T=%5.1f  tf=%.5f flown %.1f km > gate -- ladder stops', ...
               iD, iA, TN, o.tf, flyKm);
            break
        end
        % --- indirect refinement + acceptance -----------------------------
        accDz = NaN;  okCell = false;
        if ~isempty(o.lamDef)
            lamD = o.lamDef(1:7,:);
            nv = min(10, size(lamD,2));  vote = 0;
            for kv = 1:nv
                aS = o.Um(1:3,kv)/max(norm(o.Um(1:3,kv)),eps);
                pv = -lamD(4:6,kv)/max(norm(lamD(4:6,kv)),eps);
                vote = vote + sign(pv.'*aS);
            end
            sg = 1;  if vote < 0, sg = -1; end
            % With Sundman off the mesh is uniform in time and the solver
            % leaves tNodes empty -- reconstruct it.
            if isempty(o.tNodes)
                tN = linspace(0, o.tf, size(o.X,2));
            else
                tN = o.tNodes(:)';
            end
            % Hermite-Simpson defect multipliers belong at interval
            % MIDPOINTS (review finding, Gemini 2026-08-04), not left nodes.
            tMid = tN(1:end-1) + diff(tN)/2;
            for K = [12 24]
                tG = linspace(0, o.tf, K+1);
                Xg = interp1(tN, o.X', tG, 'pchip')';
                Lg = interp1(tMid, (sg*lamD)', tG, 'pchip', 'extrap')';
                seed = struct('tf', o.tf, 'tGrid', tG, 'Y', [Xg; Lg]);
                [z, info] = ms_tfmin(rv0(1:6), rvf(1:6), seed, Tnd, cnd, muStar, ...
                        struct('wallSec', 120));
                [~, rvFly] = pumpkyn.cr3bp.tfMinProp(z(8), ...
                        [rv0(1:6)'; 1; z(1:7)], Tnd, cnd, muStar);
                msKm = norm(rvFly(end,1:3) - rvf(1:3))*lStar;
                if info.converged && msKm < 100
                    evalc('zA = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), z, Tnd, cnd, muStar);');
                    accDz = norm(zA - z);
                    if accDz < 1e-6
                        okCell = true;  Z8(:,iD,iA,kr) = z;
                        RES(iD,iA,kr) = info.normR;  TF(iD,iA,kr) = z(8);
                        FLYKM(iD,iA,kr) = msKm;
                    end
                    break
                end
            end
        end
        ACCDZ(iD,iA,kr) = accDz;  OK(iD,iA,kr) = okCell;
        WALL(iD,iA,kr) = toc(t0);
        lg('  (%2d,%2d) T=%5.1f  tf=%.5f (%.3f d)  fly=%.2fkm  acc=%.1e  %s (%.0fs)', ...
           iD, iA, TN, TF(iD,iA,kr), TF(iD,iA,kr)*tStar/86400, ...
           FLYKM(iD,iA,kr), accDz, okstr(okCell), toc(t0));
        % next rung warm-starts from THIS rung's direct solution. With the
        % uniform mesh the node grid already matches, so no resampling.
        seedX = o.X;
        seedU = o.U;
        seedU(1:3,:) = seedU(1:3,:) ./ max(vecnorm(seedU(1:3,:),2,1), eps);
        seedU(4,:)   = min(max(seedU(4,:),0),1);
        seedTf = o.tf;  Tprev = TN;
        rungs_ = rungs; %#ok<NASGU>
        save(outMat, 'TF','FLYKM','ACCDZ','RES','WALL','OK','Z8','ATT', ...
             'rungs','sD','sA','meta');                        % EVERY rung
    end
    if mod(kc,10) == 0
        lg('--- %d/%d cells, %d (cell,rung) entries verified, %.1f min elapsed', ...
           kc, size(todo,1), nnz(OK), toc(tAll)/60);
    end
end
P = struct('TF',TF,'FLYKM',FLYKM,'ACCDZ',ACCDZ,'RES',RES,'WALL',WALL, ...
           'OK',OK,'Z8',Z8,'rungs',rungs,'sD',sD,'sA',sA,'meta',meta);
lg('THRUST LADDER DONE: %d verified entries over %d cells, %.1f min', ...
   nnz(OK), size(todo,1), toc(tAll)/60);
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end

function logmsg(f, s)
% LOGMSG  Append a line to the log file (or stdout if none).
% INPUTS: f path ('' = stdout); s text.  OUTPUTS: none.
if isempty(f)
    fprintf('%s\n', s);
else
    fid = fopen(f, 'a');  fprintf(fid, '%s\n', s);  fclose(fid);
end
end

function s = okstr(t)
% OKSTR  OK/fail marker.  INPUTS: t logical.  OUTPUTS: s [char].
if t, s = 'OK'; else, s = 'fail'; end
end
