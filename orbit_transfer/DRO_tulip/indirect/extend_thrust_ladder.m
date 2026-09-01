function P = extend_thrust_ladder(ladderMat, newRungs, opts)
% EXTEND_THRUST_LADDER  Continue an existing thrust ladder to lower thrust.
%
% Each phase pair resumes from its LOWEST already-converged rung and walks
% down the new rungs by pure indirect continuation: propagate that rung's
% costates, time-stretch the trajectory to the next rung's expected flight
% time, seed ms_tfmin with the stretched state+costate junctions, and gate the
% result on flown arrival plus acceptance by pumpkyn.cr3bp.tfMin. No direct
% solve is needed -- the existing entry already supplies a converged
% trajectory, which is a far better seed than any cold start.
%
% Valid while the transfer stays sub-revolution (above roughly 0.3-0.5 N for
% the tau=1 DRO): flight time scales about T^-0.56 there, and no
% revolution-count transitions intervene. Below that the transfer begins to
% wrap and this continuation is expected to fail -- that regime belongs to the
% low-thrust machinery (Sundman mesh, cold flood seeding).
%
% The ladder file is EXTENDED IN PLACE: its arrays grow along the rung
% dimension and existing results are preserved. A backup copy is written
% first.
%
% INPUTS:
%   ladderMat - engine-native ladder/sheet .mat to extend, modified in
%               place: the DRO fine sheet (thrust_ladder_12x12.mat) or any
%               catalog sheet from thrust_ladder_library (any family --
%               endpoints rebuild from the sheet's own meta recipe via
%               ladder_endpoints)
%   newRungs  - new thrust rungs, N, in DESCENDING order, all below the
%               existing minimum, e.g. [0.9 0.8 0.7 0.6 0.5]
%   opts      - (optional) struct:
%                 .K        multiple-shooting segments        [24]
%                 .wallSec  ms_tfmin budget per attempt, s    [120]
%                 .attemptSec  cumulative budget per (cell, rung), s:
%                           checked before each flight-time guess and
%                           before the tfMin acceptance call, so a
%                           grinding cell cannot stall the fleet (the
%                           2026-08-31 smoke test measured a >30 min
%                           single-cell grind at 0.45 N without it) [600]
%                 .hardCap  run each ms_tfmin guess and tfMin acceptance
%                           on a parfeval worker with a hard wall-clock
%                           cap (worker killed on timeout) -- the only
%                           bound that reaches a single crawling
%                           integration (near-primary iterate; measured
%                           twice, >4 h each). [true if Parallel toolbox
%                           licensed]
%                 .fromRungMax  only continue cells whose LOWEST converged
%                           rung is <= this thrust (N). Default: the
%                           original ladder's floor, so only full-depth
%                           cells are extended; cells stalled higher are
%                           densification work, and the continuation jump
%                           from a high rung measurably fails. inf =
%                           continue from anywhere (the old behavior)
%                           [min existing rung]
%                 .gateKm   flown-arrival gate, km            [100]
%                 .tfExp    t_f ~ T^-tfExp scaling for the
%                           flight-time guess                 [0.56]
%                 .maxCells cells to process this call        [inf]
%                 .batchSec stop cleanly between cells after
%                           this many seconds                 [inf]
%                 .logFile  append-mode log path              [stdout]
%
% OUTPUTS:
%   P - the updated ladder struct (also saved): TF, FLYKM, ACCDZ, RES, WALL,
%       OK, Z8 with the rung dimension extended, plus rungs, sD, sA, meta.
%
% REFERENCES:
%   [1] thrust_ladder_library.m -- the ladder this extends.
%   [2] ms_tfmin.m -- the multiple-shooting refinement.

if nargin < 3, opts = struct(); end
d = @(f,v) fdef(opts, f, v);
K        = d('K', 24);
wallSec  = d('wallSec', 120);
gateKm   = d('gateKm', 100);
% t_f ~ T^-tfExp. The exponent STEEPENS as thrust falls (measured: ~0.56
% between 2 and 1 N, ~1.7 between 1 and 0.9 N), so a single value makes a
% poor guess. Try several and keep the fastest solution that converges --
% this is a minimum-time problem, so smaller t_f is the better answer.
tfExpList = d('tfExpList', [0.56 1.0 1.7 2.4]);
attemptSec = d('attemptSec', 600);
% HARD per-call timeout (2026-08-31, pilot grind #2): ms_bvp's iterate
% guards bound tf and |p|, but a trust-region trial can still park a
% junction near a primary, where ONE segment integration crawls for hours
% inside every in-process fence (measured: cell 4 of the tau=2/Np=7 pilot,
% >4 h at a legal iterate). The only airtight bound is external: run each
% ms_tfmin guess and each tfMin acceptance on a parfeval worker and CANCEL
% the worker at the deadline. Requires Parallel Computing Toolbox; falls
% back to direct (unbounded) calls without it.
hardCap  = d('hardCap', license('test', 'Distrib_Computing_Toolbox'));
maxAtt   = d('maxAtt', 2);
% fromRungMax needs Q loaded -- resolved after the load below.
maxCells = d('maxCells', inf);
batchSec = d('batchSec', inf);
logFile  = d('logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

Q = load(ladderMat);
newRungs = newRungs(:)';
assert(all(diff(newRungs) < 0), 'new rungs must be in descending order');

nN = numel(newRungs);
[nD, nA] = size(Q.OK, [1 2]);
% Already extended with exactly these rungs? Then this is a resumed run:
% keep going, do not grow the arrays again and do not re-check ordering.
alreadyExtended = numel(Q.rungs) > nN && ...
    isequal(Q.rungs(end-nN+1:end), newRungs);
if alreadyExtended
    nR0 = numel(Q.rungs) - nN;
else
    nR0 = numel(Q.rungs);
    assert(all(newRungs < min(Q.rungs)), ...
        'new rungs must all lie below the existing minimum (%.3f N)', min(Q.rungs));
end
% Only continue cells converged down to this rung (default: the original
% ladder's floor -- see the todo-list comment below):
fromRungMax = d('fromRungMax', min(Q.rungs(1:nR0)));
if ~alreadyExtended                         % first extension: grow arrays
    bak = strrep(ladderMat, '.mat', sprintf('_pre_ext.mat'));
    if ~isfile(bak), copyfile(ladderMat, bak); end
    Q.rungs = [Q.rungs, newRungs];
    grow  = @(A) cat(3, A, nan(nD,nA,nN));
    Q.TF = grow(Q.TF);  Q.FLYKM = grow(Q.FLYKM);  Q.ACCDZ = grow(Q.ACCDZ);
    Q.RES = grow(Q.RES);  Q.WALL = grow(Q.WALL);
    Q.OK = cat(3, Q.OK, false(nD,nA,nN));
    Q.Z8 = cat(4, Q.Z8, nan(8,nD,nA,nN));
end
nR = numel(Q.rungs);
if ~isfield(Q,'EXT'), Q.EXT = zeros(nD,nA); end   % extension attempts
save(ladderMat, '-struct', 'Q');   % persist growth BEFORE any early return

ob = Q.meta;  muStar = ob.muStar;  lStar = ob.lStar;  tStar = ob.tStar;
g0  = 9.80665*tStar^2/(1000*lStar);
cnd = (ob.ispS/tStar)*g0;
ndT = @(TN) (TN/ob.m0kg)*tStar^2/(lStar*1000);

% Family-agnostic endpoints (2026-08-31, 0.5 N extension campaign): rebuilt
% from the sheet's OWN meta recipe via ladder_endpoints -- sheets carrying
% depFamily/arrFamily fields (halo/dpo/halo_halo/gto) route through
% get_family_orbit; legacy sheets (the DRO fine sheet, the DRO catalog)
% reproduce the original hardcoded construction bitwise
% (tests/test_ladder_endpoints.m).
[tD, rvD, tT, rvT] = ladder_endpoints(ob);

if hardCap
    pool = gcp;                    % warm once; workers inherit client path
    lg('EXTEND: hard per-call timeout ACTIVE (%d workers)', pool.NumWorkers);
else
    pool = [];
    lg('EXTEND: no Parallel toolbox -- solver calls UNBOUNDED (grind risk)');
end

% cells that have a converged rung to continue from, and no extension yet.
% By default only cells converged all the way down to fromRungMax (the
% existing ladder's floor) are continued: a cell whose ladder stalled at
% e.g. 1.5 N would face a 1.5 -> 0.45 continuation jump that measurably
% fails every guess (2026-08-31 smoke) -- such cells are densification
% work, not extension work. Pass opts.fromRungMax = inf for the old
% continue-from-anywhere behavior.
todo = [];
for iD = 1:nD
    for iA = 1:nA
        kLow = find(Q.OK(iD,iA,1:nR0), 1, 'last');
        hasOld = ~isempty(kLow) && Q.rungs(kLow) <= fromRungMax + 1e-12;
        hasNew = any(Q.OK(iD,iA,nR0+1:nR));
        if hasOld && ~hasNew && Q.EXT(iD,iA) < maxAtt
            todo(end+1,:) = [iD iA]; %#ok<AGROW>
        end
    end
end
lg('EXTEND: %d cells to continue onto rungs [%s] N (K=%d)', ...
   size(todo,1), num2str(newRungs), K);

tAll = tic;
for kc = 1:min(size(todo,1), maxCells)
    if toc(tAll) > batchSec
        lg('  [batch budget reached -- clean exit after %d cells]', kc-1);
        break
    end
    iD = todo(kc,1);  iA = todo(kc,2);
    Q.EXT(iD,iA) = Q.EXT(iD,iA) + 1;             % record BEFORE solving
    save(ladderMat, '-struct', 'Q');
    rv0 = interp1(tD, rvD, mod(Q.sD(iD),1)*tD(end), 'spline');
    rvf = interp1(tT, rvT, mod(Q.sA(iA),1)*tT(end), 'spline');

    % start from the LOWEST converged rung of this cell
    kPrev = find(Q.OK(iD,iA,1:nR0), 1, 'last');
    zPrev = Q.Z8(:,iD,iA,kPrev);
    Tprev = Q.rungs(kPrev);
    for kn = 1:nN
        kr = nR0 + kn;  TN = Q.rungs(kr);
        t0 = tic;
        okCell = false;  accDz = NaN;  z = zPrev;  info.normR = NaN;
        missKm = inf;
        try
            % propagate the previous rung's converged solution once; each
            % flight-time guess re-maps it onto a differently stretched arc
            [tj, yj] = pumpkyn.cr3bp.tfMinProp(zPrev(8), ...
                [rv0(1:6)'; 1; zPrev(1:7)], ndT(Tprev), cnd, muStar);
            lg('    prev-rung propagation (%.1f N, tf=%.4f): %.0fs', ...
               Tprev, zPrev(8), toc(t0));
            [tu, iu] = unique(tj);
            sGrid = linspace(0, 1, K+1);
            for tfExp = tfExpList
                % cumulative per-(cell,rung) budget -- a grinding cell must
                % not stall the fleet (measured 2026-08-31: >30 min on one
                % 0.45 N cell without this):
                if toc(t0) > attemptSec
                    lg('    [attemptSec budget %.0fs reached -- rung abandoned]', ...
                       attemptSec);
                    break
                end
                tfGuess = zPrev(8) * (Tprev/TN)^tfExp;
                % mass-depletion margin: the all-burn mass m = 1 - T t/c
                % hits zero at t = c/T; a guess near it lets fsolve wander
                % into the T/m blow-up where one integration grinds
                % unboundedly (the 2026-08-31 guess-4 grind). Skip it:
                if tfGuess > 0.6*cnd/ndT(TN)
                    lg('    guess tfExp=%.2f skipped: tf0=%.3f > 0.6 x depletion tf %.3f', ...
                       tfExp, tfGuess, cnd/ndT(TN));
                    continue
                end
                Yg = interp1(tu/tu(end), yj(iu,1:14), sGrid, 'pchip')';
                % mass: DERIVED from the all-burn identity m(t) = 1 - T t/c
                % rather than rescaled from the old trajectory. A derivation
                % from the invariant cannot carry a scaling mistake; the
                % previous rescale-based construction did (review finding).
                Yg(7,:) = 1 - ndT(TN)*(sGrid*tfGuess)/cnd;
                seed = struct('tf', tfGuess, 'tGrid', sGrid*tfGuess, 'Y', Yg);
                tG = tic;
                if hardCap
                    [okRun, zt, it] = run_capped(pool, @ms_tfmin, 2, ...
                        wallSec + 60, rv0(1:6), rvf(1:6), seed, ndT(TN), ...
                        cnd, muStar, struct('wallSec', wallSec));
                    if ~okRun
                        lg('    guess tfExp=%.2f (tf0=%.4f): HARD TIMEOUT -- worker killed (>%.0fs)', ...
                           tfExp, tfGuess, wallSec + 60);
                        continue
                    end
                else
                    [zt, it] = ms_tfmin(rv0(1:6), rvf(1:6), seed, ndT(TN), ...
                                        cnd, muStar, struct('wallSec', wallSec));
                end
                % NEVER fly an unconverged ms result: its tf can be garbage
                % (tens of ND), and the verification propagation below has
                % no wall cap -- flying it is the measured >70 min
                % single-cell grind of 2026-08-31, not the ms solve itself:
                if ~it.converged
                    lg('    guess tfExp=%.2f (tf0=%.4f): ms conv=0 normR=%.1e [not flown] (%.0fs)', ...
                       tfExp, tfGuess, it.normR, toc(tG));
                    continue
                end
                [~, rvFly] = pumpkyn.cr3bp.tfMinProp(zt(8), ...
                    [rv0(1:6)'; 1; zt(1:7)], ndT(TN), cnd, muStar);
                mk = norm(rvFly(end,1:3) - rvf(1:3))*lStar;
                lg('    guess tfExp=%.2f (tf0=%.4f): ms conv=1 normR=%.1e fly=%.1fkm (%.0fs)', ...
                   tfExp, tfGuess, it.normR, mk, toc(tG));
                if ~(mk < gateKm), continue, end
                if toc(t0) > attemptSec
                    lg('    [attemptSec budget reached before acceptance -- rung abandoned]');
                    break
                end
                tA = tic;
                if hardCap
                    [okRun, zA] = run_capped(pool, ...
                        @(a,b,z,T,c,mu) pumpkyn.cr3bp.tfMin(a,b,z,T,c,mu), ...
                        1, 120, rv0(1:6), rvf(1:6), zt, ndT(TN), cnd, muStar);
                    if ~okRun
                        lg('    acceptance tfMin: HARD TIMEOUT -- worker killed (>120s)');
                        continue
                    end
                else
                    evalc('zA = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), zt, ndT(TN), cnd, muStar);');
                end
                ad = norm(zA - zt);
                lg('    acceptance tfMin: |dz|=%.1e (%.0fs)', ad, toc(tA));
                if ad >= 1e-6, continue, end
                if ~okCell || zt(8) < z(8)          % keep the FASTEST valid
                    okCell = true;  z = zt;  info = it;  missKm = mk;  accDz = ad;
                end
            end
        catch ME
            lg('  (%2d,%2d) T=%5.2f  ERROR: %s', iD, iA, TN, ME.message);
            break
        end
        Q.TF(iD,iA,kr) = z(8);  Q.FLYKM(iD,iA,kr) = missKm;
        Q.RES(iD,iA,kr) = info.normR;  Q.ACCDZ(iD,iA,kr) = accDz;
        Q.WALL(iD,iA,kr) = toc(t0);  Q.OK(iD,iA,kr) = okCell;
        if okCell, Q.Z8(:,iD,iA,kr) = z; end
        save(ladderMat, '-struct', 'Q');                     % EVERY rung
        lg('  (%2d,%2d) T=%5.2f  tf=%.5f (%.3f d)  fly=%.2fkm  acc=%.1e  %s (%.0fs)', ...
           iD, iA, TN, z(8), z(8)*tStar/86400, missKm, accDz, ...
           okstr(okCell), toc(t0));
        if ~okCell, break, end                 % ladder stops at this cell
        zPrev = z;  Tprev = TN;
    end
end
P = Q;
nNew = nnz(Q.OK(:,:,nR0+1:nR));
lg('EXTEND DONE: %d new entries on rungs [%s] N, %.1f min', ...
   nNew, num2str(newRungs), toc(tAll)/60);
end

% run_capped now lives in costate_common (promoted 2026-09-01 on its
% second consumer, probe_deep_rungs -- consolidation queue item 4).

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
