function R = probe_deep_rungs(opts)
%% Purpose:
%
%   Locates the deep-thrust CLOSURE WALL of the multiple-shooting catalog
%   pipeline: starting from a converged 0.5 N entry of the DRO fine sheet,
%   walks single-cell continuation down a geometric rung schedule toward
%   25 mN and records, per rung, whether ms_tfmin closes, whether the
%   flown arrival gate holds, and whether pumpkyn tfMin still accepts the
%   result -- the roadmap's "probe 0.1 N and 25 mN on one cell" item
%   (STATUS_AND_ROADMAP §5A / §6 step 4). The wall location IS the
%   deliverable; the walk stops after two consecutive unclosed rungs.
%
%   Per the identifiability rule (OPTIMALITY_CERTIFICATION §6, two-root
%   adjudication): every converged sub-1 N probe entry is recorded with
%   its FULL ms junction states (ms_bvp info.Y on the K+1 breakpoints),
%   never bare z8 -- bare z8 pins t_f only to ~1e-4 ND at ~40 revs.
%
%   Same engine, same gates as the catalogs (ms_tfmin + flown gate + tfMin
%   acceptance), orchestrated for a single cell with per-rung records; all
%   solver calls run under costate_common/run_capped hard timeouts.
%
%  ASSUMPTIONS / NOTES:
%
% • The probe works on a SCRATCH COPY of the fine sheet's meta/entries --
%   the sheet file itself is never modified.
% • "Closed" = ms converged AND flown arrival < gateKm. The tfMin
%   acceptance verdict is recorded separately: single shooting is expected
%   to die before multiple shooting as revolutions accumulate (the
%   direct-vs-indirect sensitivity axis), and that separation is itself a
%   finding.
%
%% Inputs:
%
%  opts                     struct                  Optional fields:
%                                                   .sheetMat  fine-sheet
%                                                     path [deliverable-2
%                                                     12x12 sheet]
%                                                   .rungs [0.375 0.28 0.21
%                                                     0.16 0.12 0.09 0.067
%                                                     0.05 0.0375 0.028
%                                                     0.025] N
%                                                   .K [24] .wallSec [300]
%                                                   .gateKm [100]
%                                                   .attemptSec [1200]
%                                                   .outMat [results/
%                                                     probe_deep_rungs.mat]
%                                                   .logFile ['']
%
%% Outputs:
%
%  R                        struct                  Per-rung record arrays:
%                                                   .rungs .closed .accOk
%                                                   .tf_nd .normR .flyKm
%                                                   .accDz .wall .revs
%                                                   .Y {K+1-junction states
%                                                    per closed rung}
%                                                   .cell (iD,iA) .z8seed
%
%% Revision History:
%  M. Casey                                                   (c) 09/01/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
here     = fileparts(mfilename('fullpath'));
sheetMat = d('sheetMat', fullfile(here, '..', 'direct', 'results', ...
                                  'thrust_ladder_12x12.mat'));
rungs    = d('rungs', [0.375 0.28 0.21 0.16 0.12 0.09 0.067 0.05 ...
                       0.0375 0.028 0.025]);
K        = d('K', 24);
wallSec  = d('wallSec', 300);
gateKm   = d('gateKm', 100);
attemptSec = d('attemptSec', 1200);
outMat   = d('outMat', fullfile(here, 'results', 'probe_deep_rungs.mat'));
logFile  = d('logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

Q  = load(sheetMat);
ob = Q.meta;  muStar = ob.muStar;  lStar = ob.lStar;  tStar = ob.tStar;
g0  = 9.80665*tStar^2/(1000*lStar);
cnd = (ob.ispS/tStar)*g0;
ndT = @(TN) (TN/ob.m0kg)*tStar^2/(lStar*1000);
[tD, rvD, tT, rvT] = ladder_endpoints(ob);

%% Start cell: the FASTEST converged entry at the sheet's deepest rung
%  (shortest transfer = least accumulated sensitivity, the best-behaved
%  launch point for a deep walk):
   kr0 = numel(Q.rungs);
   assert(abs(Q.rungs(kr0) - 0.5) < 1e-9, 'fine sheet floor is not 0.5 N');
    TFf = Q.TF(:,:,kr0);  TFf(~Q.OK(:,:,kr0)) = inf;
[~, kb] = min(TFf(:));
[iD, iA] = ind2sub(size(TFf), kb);
   zPrev = squeeze(Q.Z8(:,iD,iA,kr0));
   Tprev = Q.rungs(kr0);
     rv0 = interp1(tD, rvD, mod(Q.sD(iD),1)*tD(end), 'spline');
     rvf = interp1(tT, rvT, mod(Q.sA(iA),1)*tT(end), 'spline');
lg('PROBE: cell (%d,%d), 0.5 N seed tf=%.5f ND (%.2f d); rungs [%s] N', ...
   iD, iA, zPrev(8), zPrev(8)*tStar/86400, num2str(rungs));

pool = gcp;
lg('PROBE: hard per-call timeout ACTIVE (%d workers)', pool.NumWorkers);

nR = numel(rungs);
R = struct('rungs',rungs, 'cell',[iD iA], 'z8seed',zPrev, ...
    'closed',false(1,nR), 'accOk',false(1,nR), 'tf_nd',nan(1,nR), ...
    'normR',nan(1,nR), 'flyKm',nan(1,nR), 'accDz',nan(1,nR), ...
    'wall',nan(1,nR), 'revs',nan(1,nR), 'z8',nan(8,nR), 'Y',{cell(1,nR)});
tfExpList = [0.56 1.0 1.7 2.4];
nMiss = 0;
for kr = 1:nR
    TN = rungs(kr);
    t0 = tic;
    best = [];
    [tj, yj] = pumpkyn.cr3bp.tfMinProp(zPrev(8), ...
        [rv0(1:6)'; 1; zPrev(1:7)], ndT(Tprev), cnd, muStar);
    [tu, iu] = unique(tj);
    sGrid = linspace(0, 1, K+1);
    for tfExp = tfExpList
        if toc(t0) > attemptSec
            lg('  T=%.4f: [attemptSec %.0fs reached]', TN, attemptSec); break
        end
        tfGuess = zPrev(8) * (Tprev/TN)^tfExp;
        if tfGuess > 0.6*cnd/ndT(TN)
            lg('  T=%.4f guess tfExp=%.2f skipped (depletion margin)', TN, tfExp);
            continue
        end
        Yg = interp1(tu/tu(end), yj(iu,1:14), sGrid, 'pchip')';
        Yg(7,:) = 1 - ndT(TN)*(sGrid*tfGuess)/cnd;
        seed = struct('tf', tfGuess, 'tGrid', sGrid*tfGuess, 'Y', Yg);
        tG = tic;
        [okRun, zt, it] = run_capped(pool, @ms_tfmin, 2, wallSec + 90, ...
            rv0(1:6), rvf(1:6), seed, ndT(TN), cnd, muStar, ...
            struct('wallSec', wallSec));
        if ~okRun
            lg('  T=%.4f guess tfExp=%.2f: HARD TIMEOUT (>%ds)', TN, tfExp, wallSec+90);
            continue
        end
        if ~it.converged
            lg('  T=%.4f guess tfExp=%.2f: ms conv=0 normR=%.1e (%.0fs)', ...
               TN, tfExp, it.normR, toc(tG));
            continue
        end
        [~, rvFly] = pumpkyn.cr3bp.tfMinProp(zt(8), ...
            [rv0(1:6)'; 1; zt(1:7)], ndT(TN), cnd, muStar);
        mk = sqrt(sum((rvFly(end,1:3) - rvf(1:3)).^2))*lStar;
        lg('  T=%.4f guess tfExp=%.2f: ms conv=1 normR=%.1e fly=%.1fkm (%.0fs)', ...
           TN, tfExp, it.normR, mk, toc(tG));
        if mk >= gateKm, continue, end
        if isempty(best) || zt(8) < best.z(8)
            best = struct('z', zt, 'it', it, 'mk', mk);
        end
    end
    if ~isempty(best)
        % revolutions: TOTAL SWEPT Moon-centered angle of the flight (the
        % end-minus-start form under-reads petal geometries -- measured
        % 0.2 "revs" on a 1.3-rev transfer, fixed 2026-09-01):
        [~, yb] = pumpkyn.cr3bp.tfMinProp(best.z(8), ...
            [rv0(1:6)'; 1; best.z(1:7)], ndT(TN), cnd, muStar);
        ang = unwrap(atan2(yb(:,2), yb(:,1) - (1 - muStar)));
        R.revs(kr) = sum(abs(diff(ang)))/(2*pi);
        % tfMin acceptance under the hard cap (recorded, not required):
        [okA, zA] = run_capped(pool, ...
            @(a,b,z,T,c,mu) pumpkyn.cr3bp.tfMin(a,b,z,T,c,mu), 1, 180, ...
            rv0(1:6), rvf(1:6), best.z, ndT(TN), cnd, muStar);
        if okA, R.accDz(kr) = sqrt(sum((zA - best.z).^2)); end
        R.accOk(kr)  = okA && R.accDz(kr) < 1e-6;
        R.closed(kr) = true;
        R.tf_nd(kr)  = best.z(8);  R.normR(kr) = best.it.normR;
        R.flyKm(kr)  = best.mk;    R.z8(:,kr)  = best.z;
        R.Y{kr}      = best.it.Y;              % FULL ms junction states
        zPrev = best.z;  Tprev = TN;  nMiss = 0;
        lg('PROBE T=%.4f N CLOSED: tf=%.5f ND (%.2f d) revs=%.1f fly=%.2fkm accOk=%d accDz=%.1e (%.0fs)', ...
           TN, best.z(8), best.z(8)*tStar/86400, R.revs(kr), best.mk, ...
           R.accOk(kr), R.accDz(kr), toc(t0));
    else
        nMiss = nMiss + 1;
        R.wall(kr) = toc(t0);
        lg('PROBE T=%.4f N FAILED (%.0fs) [consecutive misses: %d]', ...
           TN, toc(t0), nMiss);
        if nMiss >= 2
            lg('PROBE WALL LOCATED: two consecutive failures at %.4f and %.4f N', ...
               rungs(kr-1), rungs(kr));
            break
        end
    end
    R.wall(kr) = toc(t0);
    save(outMat, 'R');                          % after EVERY rung
end
save(outMat, 'R');
if any(R.closed), deepest = min(rungs(R.closed)); else, deepest = NaN; end
lg('PROBE DONE: deepest closed rung %.4f N', deepest);
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end

function logmsg(f, s)
% LOGMSG  Append a line to the log file (or stdout if none).
% INPUTS: f path ('' = stdout); s text.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f,'a'); fprintf(fid,'%s\n',s); fclose(fid);
end
end
