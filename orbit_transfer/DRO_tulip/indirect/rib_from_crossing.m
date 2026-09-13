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
%   .targets [] EXPLICIT unwrapped departure offsets from anc.sD to walk to,
%   in order -- use these when the caller needs exact phases; otherwise
%   .nD [12] grid size, .direction [+1], .nPts [nD-1] points to walk,
%   .maxBisect [8] halvings allowed per grid step (the 2026-09-09 sweep
%   needed ~1/1728 of a period on the hardest step), .wallSec [600],
%   .copts (certify_root options),
%   .checkpoint '' path of the unit's resumable checkpoint (walk_checkpoint):
%   written after every accepted point, resumed from when it matches this
%   walk's identity; .problem [] identity struct stored in it,
%   .progress [] function handle called with NO arguments after EVERY
%   solve (about every 2-3 minutes) -- the campaign worker passes its
%   heartbeat here. Without it a claim on a multi-hour column went stale
%   after 30 minutes and could be reclaimed mid-walk (Astra 2026-09-13).
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
% EXPLICIT TARGETS beat a derived step count. `round(1/|delta|)` does not
% reproduce an arbitrary phase -- 5/12 becomes round(2.4) = 2 and walks a
% half period -- and a caller that then labels the answer with what it asked
% for has certified one phase and reported another. (Astra 2026-09-10.)
targets = d('targets', []);
if ~isempty(targets)
    targets = targets(:).';
    nPts = numel(targets);
else
    targets = dirn*(1:nPts)/nD;
end
maxBisect = d('maxBisect', 8);  wallSec = d('wallSec', 600);
copts = d('copts', struct());  copts.wallSec = wallSec;  copts.sA = C0.sA;
logFile = d('logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));
progress = d('progress', []);
if isempty(progress), progress = @() []; end
copts.progress = progress;                   % certify_root ticks it after every capped stage

R = struct('pts', struct([]), 'stop', '', 'nSolve', 0, 'sA', C0.sA, 'nD', nD);
rvf = B.stateA(C0.sA);
sD = anc.sD;  z = C0.z;  Y = C0.Y;  K = size(Y, 2);

% RESUME from the unit's checkpoint if one matches this exact walk: same
% column identity, lattice, direction, targets and problem. A killed
% attempt used to cost the whole column (up to nine hours).
ckptFile = d('checkpoint', '');
ident = struct('sA', C0.sA, 'nD', nD, 'dirn', dirn, 'targets', targets, 'sD0', anc.sD, ...
               'problem', d('problem', struct()));
kStart = 1;
if ~isempty(ckptFile)
    [Ck0, why] = walk_checkpoint('load', ckptFile, ident);
    if ~isempty(Ck0)
        kStart = Ck0.k + 1;  sD = Ck0.sD;  z = Ck0.z;  Y = Ck0.Y;
        R.pts = Ck0.pts;  R.nSolve = Ck0.nSolve;
        lg('  rib RESUMED at point %d of %d from %s (saved %s)', kStart, nPts, ckptFile, Ck0.saved);
        if kStart > nPts, R.stop = 'complete';  return, end
    elseif isfile(ckptFile)
        lg('  rib checkpoint %s IGNORED: %s', ckptFile, why);
    end
end

for k = kStart:nPts
    target = anc.sD + targets(k);          % unwrapped, absolute
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
        try progress(); catch, end         % a beat that fails must not stop the walk
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
    if ~isempty(ckptFile)              % after every ACCEPTED point, atomically
        walk_checkpoint('save', ckptFile, struct('identity', ident, 'k', k, 'sD', sD, 'z', z, 'Y', Y, ...
                                                  'pts', R.pts, 'nSolve', R.nSolve));
    end
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
