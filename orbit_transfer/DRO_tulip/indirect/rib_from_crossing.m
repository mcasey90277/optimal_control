function R = rib_from_crossing(C0, B, anc, opts)
%% Purpose:
%
%   Hang a DEPARTURE-phase rib off one certified arrival-phase crossing:
%   hold the arrival phase and walk the departure phase in 1/nD steps,
%   re-solving at each step from the previous point's junctions and putting
%   every candidate through the full gate stack (certify_root) before it is
%   kept. Departure phase moves only the first junction's fixed state, so
%   the steps are cheap and a fixed walker with bisection beats an
%   arclength arc here; the arrival direction is where the sensitivity and
%   the folds live, and that is arclength_arrival's job.
%
%   A step that fails is BISECTED (up to maxBisect halvings). A rib that
%   cannot advance stops and says so in R.stop -- it never returns a short
%   list as if it had finished.
%
%% Inputs:
%
%  C0                       struct                  certified crossing
%                                                   (certify_crossing /
%                                                   certify_root output:
%                                                   .z .Y .sA)
%  B, anc                   struct                  arclength_arrival('setup')
%  opts                     struct (optional)
%   .nD [12] grid size, .direction [+1], .nPts [nD-1] points to walk,
%   .maxBisect [8] halvings allowed per grid step (the 2026-09-09 sweep
%   needed ~1/1728 of a period on the hardest step), .wallSec [600],
%   .copts (certify_root options),
%   .logFile ''
%
%% Outputs:
%
%  R                        struct                  .pts (struct array of
%                                                   certify_root outputs,
%                                                   certified only) .stop
%                                                   .nSolve .sA .nD
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 4, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
nD = d('nD', 12);  dirn = d('direction', +1);  nPts = d('nPts', nD - 1);
maxBisect = d('maxBisect', 8);  wallSec = d('wallSec', 600);
copts = d('copts', struct());  copts.wallSec = wallSec;  copts.sA = C0.sA;
logFile = d('logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

R = struct('pts', struct([]), 'stop', '', 'nSolve', 0, 'sA', C0.sA, 'nD', nD);
rvf = B.stateA(C0.sA);
sD = anc.sD;  z = C0.z;  Y = C0.Y;  K = size(Y, 2);

for k = 1:nPts
    target = sD + dirn/nD;
    step = dirn/nD;  nb = 0;  cur = sD;  zc = z;  Yc = Y;  Ck = [];  nGood = 0;
    while abs(wrapDiff(cur, target)) > 1e-12
        % never overshoot the grid point, and RESTORE the step after two
        % clean sub-steps -- a walk that only ever halves crawls at the
        % depth of its worst patch for the rest of the arc
        if abs(step) > abs(wrapDiff(cur, target)), step = wrapDiff(cur, target); end
        trial = cur + step;
        seed = struct('tf', zc(8), 'tGrid', linspace(0, zc(8), K+1), 'Y', [Yc, Yc(:,end)]);
        seed.Y(8:14, 1) = zc(1:7);
        rv0 = B.stateD(trial);
        seed.Y(1:7, 1) = [rv0(1:6); 1];
        copts.sD = mod(trial, 1);
        Ct = certify_root(seed, rv0, rvf, B, copts);
        R.nSolve = R.nSolve + 1;
        if Ct.ok
            cur = trial;  zc = Ct.z;  Yc = Ct.Y;  Ck = Ct;  nGood = nGood + 1;
            lg('  rib sD %.4f: t_f %.4f d, %s', mod(cur,1), Ct.tfDays, Ct.reason);
            if nGood >= 2 && abs(step) < 1/nD
                step = 2*step;  nGood = 0;  nb = max(nb - 1, 0);
            end
        else
            if nb >= maxBisect
                R.stop = sprintf('stalled at sD = %.4f stepping to %.4f: %s', mod(cur,1), mod(trial,1), Ct.reason);
                lg('  rib STOP: %s', R.stop);
                return
            end
            step = step/2;  nb = nb + 1;  nGood = 0;
            lg('  rib sD %.4f -> %.4f refused (%s); halving to %.5f', mod(cur,1), mod(trial,1), Ct.reason, step);
        end
    end
    sD = mod(target, 1);  z = zc;  Y = Yc;
    Ck.sD = sD;
    if isempty(R.pts), R.pts = Ck; else, R.pts(end+1) = Ck; end
end
R.stop = 'complete';
end

function d_ = wrapDiff(a, b)
% WRAPDIFF  Signed difference b - a wrapped to [-0.5, 0.5).
% INPUTS: a; b.  OUTPUTS: d_.
d_ = mod(b - a + 0.5, 1) - 0.5;
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end

function logmsg(f, s)
% LOGMSG  Print, and append to a log file when one is named.  INPUTS: f; s.
fprintf('%s\n', s);
if ~isempty(f), fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid); end
end
