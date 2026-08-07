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
%   ladderMat - thrust_ladder_12x12.mat to extend (modified in place)
%   newRungs  - new thrust rungs, N, in DESCENDING order, all below the
%               existing minimum, e.g. [0.9 0.8 0.7 0.6 0.5]
%   opts      - (optional) struct:
%                 .K        multiple-shooting segments        [24]
%                 .wallSec  ms_tfmin budget per attempt, s    [120]
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
maxAtt   = d('maxAtt', 2);
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

[~, rvD0] = pumpkynPie.cr3bp.getDRO(ob.tauDRO);
rvD0 = pumpkyn.cr3bp.cont_np(rvD0, ob.tauDRO, muStar, 1e-12);
[tD, rvD] = pumpkyn.cr3bp.prop(ob.tauDRO, rvD0, muStar);
[~, rvT0] = pumpkyn.cr3bp.getTulip(ob.tauTulip, ob.NpTulip, ob.pmTulip);
rvT0 = pumpkyn.cr3bp.cont_np(rvT0, ob.tauTulip, muStar, 1e-12);
[tT, rvT] = pumpkyn.cr3bp.prop(ob.tauTulip, rvT0, muStar);

% cells that have a converged rung to continue from, and no extension yet
todo = [];
for iD = 1:nD
    for iA = 1:nA
        hasOld = any(Q.OK(iD,iA,1:nR0));
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
            [tu, iu] = unique(tj);
            sGrid = linspace(0, 1, K+1);
            for tfExp = tfExpList
                tfGuess = zPrev(8) * (Tprev/TN)^tfExp;
                Yg = interp1(tu/tu(end), yj(iu,1:14), sGrid, 'pchip')';
                % mass rescale: consumed propellant scales with BOTH the time
                % stretch AND the thrust ratio (mdot ~ T) -- review finding
                Yg(7,:) = 1 + (Yg(7,:)-1)*(tfGuess/zPrev(8))*(TN/Tprev);
                seed = struct('tf', tfGuess, 'tGrid', sGrid*tfGuess, 'Y', Yg);
                [zt, it] = ms_tfmin(rv0(1:6), rvf(1:6), seed, ndT(TN), cnd, ...
                                    muStar, struct('wallSec', wallSec));
                [~, rvFly] = pumpkyn.cr3bp.tfMinProp(zt(8), ...
                    [rv0(1:6)'; 1; zt(1:7)], ndT(TN), cnd, muStar);
                mk = norm(rvFly(end,1:3) - rvf(1:3))*lStar;
                if ~(it.converged && mk < gateKm), continue, end
                evalc('zA = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), zt, ndT(TN), cnd, muStar);');
                ad = norm(zA - zt);
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
