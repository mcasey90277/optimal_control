function R = probe_abstract_case(opts)
%% PROBE (not production) -- the search that FOUND the 70 mN case.
%
%   An Isp ladder followed by a thrust walk, used once to locate the
%   abstract's operating point. Kept as the record of that search. To obtain
%   a 70 mN transfer now, use run_dro_tulip(sD, sA), which is tested and
%   puts every answer through the gate stack.
%
%% Purpose:
%
%   Solve the CISLUNAR ABSTRACT's operating point: minimum-time DRO -> 7-petal
%   tulip transfer for a 150 kg spacecraft with a 70 mN Hall thruster at
%   Isp 900 s ("MinTime Tulip Transfer Abstract.docx").
%
%   That point is NOT in the shipped catalog, which floors at 0.5 N and runs
%   Isp 1710 s. The abstract's quoted 18 days / 0.75 km/s follow from the
%   all-burn identity dV = (T/m0)*tf rather than from a converged transfer,
%   so the poster currently rests on an extrapolation seven times below the
%   lowest catalogued thrust. This driver produces the real number, or
%   establishes that the point is unreachable by continuation and says so.
%
%   Route: the deep probe (probe_deep_rungs, 2026-09-01) already closed
%   0.09 N at 13.11 d on cell (1,11) of the same sheet and BANKED its K+1
%   junction states; it then stalled at 0.067 N on a basin/winding wall
%   (FINDINGS section on the deep-thrust probe). 70 mN lies inside that
%   bracket, so the walk is short but crosses the approach to a known wall.
%
%   Two stages, easy one first:
%     1. Isp 1710 -> 900 at fixed 0.09 N. Mild: it only changes the mass
%        depletion rate, and FASTER depletion raises late-arc acceleration,
%        which should help rather than hurt the descent that follows.
%     2. Thrust 0.09 -> 0.07 N at Isp 900, fine steps with bisection.
%
%   Each accepted step re-seeds the next from its own junction states.
%
%% Inputs:
%
%  opts                     struct (optional)
%   .ispLadder              [1 x n]                 [1710 1400 1150 900] s
%   .thrustLadder           [1 x m]                 [0.09 0.085 0.08 0.075
%                                                   0.07] N
%   .K                      scalar                  segments [24]
%   .wallSec                scalar                  per-solve budget [300]
%   .gateKm                 scalar                  flown-arrival gate [100]
%   .maxBisect              scalar                  inserts per failed gap [3]
%   .tfExps                 [1 x k]                 t_f guess exponents per
%                                                   rung, tf ~ tf*(T/T')^e
%                                                   [1.0 1.3 1.7 2.2 3.0];
%                                                   the sweep is the route
%                                                   past a WINDING wall
%   .cell                   [iD iA]                 walk this fine-sheet cell
%                                                   from its 0.5 N entry
%   .resumeMat              char                    start from the last step
%                                                   of a previous run's .mat
%                                                   (skips stage 1)
%   .outMat, .logFile       char
%
%% Outputs:
%
%  R                        struct                  .steps (per accepted
%                                                   step: stage, isp, T_N,
%                                                   tf_nd, tf_days, dV_kms,
%                                                   prop_kg, normR, flyKm,
%                                                   accDz, revs, z8, Y),
%                                                   .reached, .wall
%
%% Revision History:
%  M. Casey                                                   (c) 09/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
here    = fileparts(mfilename('fullpath'));
ispLad  = d('ispLadder', [1710 1400 1150 900]);
thrLad  = d('thrustLadder', [0.09 0.085 0.08 0.075 0.07]);
K       = d('K', 24);
wallSec = d('wallSec', 300);
gateKm  = d('gateKm', 100);
maxBis  = d('maxBisect', 3);
tfExps  = d('tfExps', [1.0 1.3 1.7 2.2 3.0]);
outMat  = d('outMat', fullfile(here, 'results', 'probe_abstract_case.mat'));
logFile = d('logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

addpath(fullfile(fileparts(here), '..', 'costate_common'));
Q  = load(fullfile(here, '..', 'direct', 'results', 'thrust_ladder_12x12.mat'));
ob = Q.meta;  muStar = ob.muStar;  lStar = ob.lStar;  tStar = ob.tStar;
% THE shared propulsion conversion (costate_common/nd_propulsion)
ndp = nd_propulsion([], [], ob.m0kg, lStar, tStar);
cOf = ndp.ndC;                                    % ND exhaust speed at any Isp
ndT = ndp.ndT;
[tD, rvD, tT, rvT] = ladder_endpoints(ob);
% THE shared endpoint rule (costate_common/phase_state, FINDINGS 44)
depAt = phase_state(tD, rvD);   arrAt = phase_state(tT, rvT);

% RESUME: start from the last accepted step of a previous run instead of the
% banked 0.09 N rung (so a deeper push can change K or the guess sweep
% without re-walking the rungs that already closed).
% CELL: start from a chosen cell's 0.5 N fine-sheet entry instead of the
% banked deep rung. The deep wall is a BASIN property of the cell, so a cell
% other than (1,11) -- which was picked as the FASTEST 0.5 N entry, i.e. the
% least-wound one -- may continue further.
cellSel = d('cell', []);
resumeMat = d('resumeMat', '');
if ~isempty(cellSel)
    kr0 = numel(Q.rungs);
    iDs = cellSel(1);  iAs = cellSel(2);
    assert(Q.OK(iDs, iAs, kr0), 'cell (%d,%d) has no 0.5 N entry', iDs, iAs);
    z0 = squeeze(Q.Z8(:, iDs, iAs, kr0));
    [tj, yj] = pumpkyn.cr3bp.tfMinProp(z0(8), ...
        [depAt(Q.sD(iDs)); 1; z0(1:7)], ...
        ndT(Q.rungs(kr0)), cOf(ob.ispS), muStar);
    Y0 = flight_to_junctions(tj, yj, K);        % THE shared cut
    P = struct('R', struct('closed', true, 'rungs', Q.rungs(kr0), ...
                           'cell', [iDs iAs], 'tf_nd', z0(8), 'Y', {{Y0}}));
elseif ~isempty(resumeMat)
    Rp = load(resumeMat);  sp = Rp.R.steps(end);
    P = struct('R', struct('closed', true, 'rungs', sp.T_N, 'cell', Rp.R.cell, ...
                           'tf_nd', sp.tf_nd, 'Y', {{sp.Y}}));
    ob.ispS = sp.isp;                        % resume already carries the Isp
    ispLad  = sp.isp;                        % stage 1 is done
else
    P = load(fullfile(here, 'results', 'probe_deep_rungs.mat'));
end
if false
P = load(fullfile(here, 'results', 'probe_deep_rungs.mat'));
end
pr = P.R;
ic = find(pr.closed, 1, 'last');
assert(~isempty(ic), 'no closed rung to seed from');
iD = pr.cell(1);  iA = pr.cell(2);
rv0 = depAt(Q.sD(iD)).';            % row, as before
rvf = arrAt(Q.sA(iA)).';
lg(['ABSTRACT CASE: DRO tau=%g -> tulip Np=%g, cell (%d,%d), m0=%g kg.\n' ...
    'Seed: banked %.3f N solution, tf=%.4f ND (%.2f d), Isp %g s.\n' ...
    'Target: %.3f N at Isp %g s.'], ob.tauDRO, ob.NpTulip, iD, iA, ob.m0kg, ...
    pr.rungs(ic), pr.tf_nd(ic), pr.tf_nd(ic)*tStar/86400, ob.ispS, ...
    thrLad(end), ispLad(end));

pool = gcp;
Ycur = pr.Y{ic};                                   % [14 x K+1] junctions
tfCur = pr.tf_nd(ic);
if size(Ycur,2) ~= K+1
    sg = linspace(0,1,K+1);
    Ycur = interp1(linspace(0,1,size(Ycur,2))', Ycur', sg, 'pchip')';
end
ispCur = ob.ispS;  Tcur = pr.rungs(ic);
R = struct('steps', struct([]), 'reached', false, 'wall', NaN, ...
           'cell', [iD iA], 'meta', ob);

%% Stage 1: Isp, then Stage 2: thrust ------------------------------------
for stage = 1:2
    if stage == 1, ladder = ispLad(2:end); else, ladder = thrLad(2:end); end
    queue = ladder(:)';  nBis = 0;  good = [];
    if stage == 1, good = ispCur; else, good = Tcur; end
    while ~isempty(queue)
        x = queue(1);  queue(1) = [];
        if stage == 1, ispTry = x;  Ttry = Tcur; else, ispTry = ispCur;  Ttry = x; end
        cnd = cOf(ispTry);  Tnd = ndT(Ttry);
        % t_f GUESS SWEEP. Minimum time scales roughly as 1/acceleration, but
        % a single guess can only find the branch it lands in: the deep wall
        % of 2026-09-01 is a WINDING wall, where the solution past it needs
        % more revolutions than a linearly stretched seed carries. Sweeping
        % the exponent walks t_f upward at fixed thrust, which is the
        % recorded follow-up route. Several guesses may converge to DIFFERENT
        % branches, so keep the SHORTEST -- this is a minimum-time problem.
        t0 = tic;  z = [];  ok = false;  flyKm = NaN;  revs = NaN;
        best = inf;
        for tfExp = tfExps
            if stage == 1
                % t_f shrinks as Isp drops (faster depletion -> more
                % late-arc acceleration); scan a small bracket around it
                tfG = tfCur * (1 - 0.06*(tfExp - 1));
            else
                tfG = tfCur * (Tcur/Ttry)^tfExp;
            end
            if tfG > 0.6*cnd/Tnd, continue, end          % mass-depletion margin
            Yg = Ycur;  sg = linspace(0,1,size(Yg,2));
            Yg(7,:) = 1 - Tnd*(sg*tfG)/cnd;              % all-burn mass profile
            seed = struct('tf', tfG, 'tGrid', sg*tfG, 'Y', Yg);
            [okRun, zt, itt] = run_capped(pool, @ms_tfmin, 2, wallSec + 90, ...
                rv0(1:6), rvf(1:6), seed, Tnd, cnd, muStar, struct('wallSec', wallSec));
            if ~(okRun && itt.converged), continue, end
            [fk, rv] = flownGate(zt, rv0, rvf, Tnd, cnd, muStar, lStar);
            if fk >= gateKm, continue, end
            if zt(8) < best
                best = zt(8);  z = zt;  it = itt;  ok = true;  flyKm = fk;  revs = rv;
            end
            % (no early break: t_f moves with Isp too, and at 0.5 N that
            % move is large enough that a single guess fails -- measured on
            % cells (3,1) and (1,9), 2026-09-08)
        end
        if ok
            [accOk, accDz] = tfMinAccept(z, rv0, rvf, Tnd, cnd, muStar);
            mf = 1 - Tnd*z(8)/cnd;
            s = struct('stage', stage, 'isp', ispTry, 'T_N', Ttry, ...
                'tf_nd', z(8), 'tf_days', z(8)*tStar/86400, ...
                'dV_kms', cnd*log(1/mf)*lStar/tStar, ...
                'prop_kg', ob.m0kg*(1-mf), 'normR', it.normR, 'flyKm', flyKm, ...
                'accOk', accOk, 'accDz', accDz, 'revs', revs, 'z8', z(:), ...
                'Y', it.Y, 'wallSec', toc(t0));
            if isempty(R.steps), R.steps = s; else, R.steps(end+1) = s; end
            Ycur = it.Y;  tfCur = z(8);
            if stage == 1, ispCur = ispTry;  good = ispTry;
            else,          Tcur = Ttry;      good = Ttry; end
            lg(['  [s%d] Isp=%g T=%.4f N: tf=%.4f ND (%.2f d) dV=%.4f km/s ' ...
                'prop=%.2f kg revs=%.2f normR=%.1e fly=%.3f km tfMin=%d (%.0fs)'], ...
                stage, ispTry, Ttry, z(8), s.tf_days, s.dV_kms, s.prop_kg, ...
                revs, it.normR, flyKm, accOk, toc(t0));
            save(outMat, 'R');
        else
            lg('  [s%d] Isp=%g T=%.4f N: FAIL (no guess in [%s] converged+arrived) (%.0fs)', ...
               stage, ispTry, Ttry, num2str(tfExps), toc(t0));
            if ~isempty(good) && nBis < maxBis && abs(good - x)/abs(good) > 0.01
                xMid = sqrt(good*x);
                queue = [xMid, x, queue]; %#ok<AGROW>
                nBis = nBis + 1;
                lg('    bisect -> %.5g', xMid);
            else
                R.wall = x;
                lg('  [s%d] WALL at %.5g (bisection budget %d used)', stage, x, nBis);
                break
            end
        end
    end
    if ~isnan(R.wall), break, end
end

R.reached = ~isempty(R.steps) && abs(R.steps(end).T_N - thrLad(end)) < 1e-12 && ...
            abs(R.steps(end).isp - ispLad(end)) < 1e-9;
save(outMat, 'R');
if R.reached
    s = R.steps(end);
    lg(['ABSTRACT CASE REACHED: %.0f mN, Isp %g s, %.2f DAYS, dV %.4f km/s, ' ...
        '%.2f kg propellant (%.1f%% of wet mass), %.2f revs'], ...
        s.T_N*1000, s.isp, s.tf_days, s.dV_kms, s.prop_kg, ...
        100*s.prop_kg/ob.m0kg, s.revs);
else
    lg('ABSTRACT CASE NOT REACHED: wall at %.5g', R.wall);
end
end

% ------------------------------------------------------------------------
function [flyKm, revs] = flownGate(z, rv0, rvf, Tnd, cnd, muStar, lStar)
% FLOWNGATE  Fly the PMP control end to end; arrival error in km and the
% swept revolutions about the Moon.  INPUTS: z; rv0; rvf; Tnd; cnd; muStar;
% lStar.  OUTPUTS: flyKm; revs.
[~, Y] = pumpkyn.cr3bp.tfMinProp(z(8), [rv0(1:6)'; 1; z(1:7)], Tnd, cnd, muStar);
flyKm = sqrt(sum((Y(end,1:3) - rvf(1:3)).^2))*lStar;
th = unwrap(atan2(Y(:,2), Y(:,1) - (1 - muStar)));
revs = (max(th) - min(th))/(2*pi);
end

function [ok, dz] = tfMinAccept(z, rv0, rvf, Tnd, cnd, muStar)
% TFMINACCEPT  Hand z8 to pumpkyn tfMin and measure how far it moves -- the
% foreign witness. INPUTS: z; rv0; rvf; Tnd; cnd; muStar. OUTPUTS: ok; dz.
ok = false;  dz = NaN;
try
    zt = pumpkyn.cr3bp.tfMin(rv0(1:6)', rvf(1:6)', z(:), Tnd, cnd, muStar);
    dz = norm(zt(:) - z(:));  ok = dz < 1e-6;
catch
end
end

function v = fieldd(s, f, d)
% FIELDD  Field with default.  INPUTS: s; f; d.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function logmsg(f, s)
% LOGMSG  Append to a log file or stdout.  INPUTS: f; s.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
end
