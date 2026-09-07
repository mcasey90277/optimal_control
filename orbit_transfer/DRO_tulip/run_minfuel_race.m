function out = run_minfuel_race(cfg)
% RUN_MINFUEL_RACE  Step-5 Task 3: the energy->fuel continuation RACE on one
% (cell, gamma) -- discrete Bertrand-Epenoy 'eps' ladder vs the PLQ 'huber'
% kappa walk (Mike's pinned experiment, 2026-08-31), both from the SAME
% min-energy seed, both through fixed-tf ms_minfuel + the pilot gates.
%
% Both arms walk the same decreasing parameter schedule. A failed step
% bisects toward the last good parameter (up to maxBisect inserts); an arm
% retires after two ABANDONED gaps. Every accepted step records the final
% mass (from the arm's own flight), coast fraction, Hdrift, iterations and
% wall time; the final step also runs the single-shooting acceptance gate.
%
% PRE-REGISTERED SCORING (design 2026-09-01): failed/bisected steps, basin
% flips (m_f discontinuities along the walk), total wall, and -- decisive --
% final m_f agreement between arms (same limit problem). Structural family
% differences already on record in cr3bp_minfuel_pmp's header (huber's
% switch jump; huber never coasts exactly).
%
% INPUTS:
%   cfg - (optional) struct overriding ADJUSTABLE PARAMETERS below.
%
% OUTPUTS:
%   out - struct: .arms (per-family step records), .compare, .seedRec info;
%         saved to direct/results/minfuel_race.mat after every step.
%
% REFERENCES:
%   [1] costate_common/cr3bp_minfuel_pmp (families + recorded theory).
%   [2] run_minenergy_pilot.m (the seed records and gate conventions).
%   [3] GTO_tulip/direct/lib/sundman_homotopy.m (the direct-side eps walk).

%% ADJUSTABLE PARAMETERS ---------------------------------------------------
P.cell     = [2 5];              % (iD, iA) of the seed record
P.gamma    = 1.2;                % t_f multiplier of the seed record
P.sched    = [1 0.7 0.5 0.35 0.25 0.18 0.12 0.08 0.05 0.03 ...
              0.02 0.012 0.008 0.005 0.003 0.002 0.001];
P.families = {'eps', 'huber'};
P.wallSec  = 300;                % ms budget per solve (hard-capped +90 s)
P.maxBisect = 3;                 % midpoint inserts per failed gap
P.maxGaps  = 2;                  % abandoned gaps before the arm retires
P.HdriftTol = 1e-6;              % absolute first-integral gate (pilot rule)
P.outMat   = '';                 % '' = direct/results/minfuel_race.mat
P.logFile  = '';
P.rungTolR  = [];                % LOOSE-RUNG gate (MfMax review 2026-09-06,
                                 % idea 1.1): [] = a rung must CONVERGE to
                                 % ms tolR (the 09-02 protocol); a number
                                 % accepts a rung when normR < rungTolR
                                 % (MfMax's homCI uses 1e-3) -- continuation
                                 % points are predictors, not deliverables.
                                 % The FLOOR is always refined to full tolR.
P.rungHdriftTol = [];            % Hdrift gate for LOOSE rungs ([] = HdriftTol)
P.delta     = [];                % huberc ramp width during the p-walk: [] =
                                 % the field default delta = p; a scalar =
                                 % FIXED delta; a function handle @(p) ->
                                 % delta (e.g. @(p) 2*p to match eps's width).
                                 % Ignored by eps/huber. (FINDINGS 25 knob.)
P.deltaSched = [];               % optional SECOND STAGE: after the p-walk,
                                 % hold p at its deepest value and walk delta
                                 % down through this decreasing schedule
                                 % (same gates, bisection on delta).
P.seedScale = 1;                 % multiply the seed COSTATES (rows 8:14) by
                                 % this before the first rung. 1 = the energy
                                 % costates verbatim (the 09-02 protocol).
                                 % huber kappa=1 minimises at s* = Q, not Q/2,
                                 % so 0.5 is its warm seed (Astra, 2026-09-05).
if nargin >= 1 && ~isempty(cfg)
    fn = fieldnames(cfg);
    for k = 1:numel(fn), P.(fn{k}) = cfg.(fn{k}); end
end
lg = @(varargin) logmsg(P.logFile, sprintf(varargin{:}));

%% Paths and the seed record ----------------------------------------------
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'indirect'), fullfile(fileparts(here), 'costate_common'));
if isempty(which('pumpkyn.cr3bp.tfMinEoM'))
    run(fullfile(fileparts(here), 'cr3bp_common', 'setup_cr3bp_common.m'));
end
if isempty(which('casadi.SX')), addpath(fullfile(getenv('HOME'), 'casadi-3.7.0')); end
outMat = P.outMat;
if isempty(outMat)
    outMat = fullfile(here, 'direct', 'results', 'minfuel_race.mat');
end
L = load(fullfile(here, 'direct', 'results', 'minenergy_pilot.mat'));
ki = find(arrayfun(@(r) isequal([r.iD r.iA], P.cell) && ...
                        abs(r.gam - P.gamma) < 1e-9, L.R), 1);
assert(~isempty(ki), 'no minenergy_pilot record for cell (%d,%d) gamma %.2f', ...
       P.cell(1), P.cell(2), P.gamma);
rec = L.R(ki);
tf = rec.tf;  Tmax = rec.Tmax;  c = rec.c;  muStar = rec.muStar;
rv0 = rec.rv0(:);  rvf = rec.rvf(:);

% Seed junctions: the pilot record stores only the converged lambda0 and K
% (compact summary), so rebuild the junction states by flying the ENERGY
% field from [rv0; 1; lambda0] and sampling the K+1 breakpoints -- smooth
% field, few revs, and the eps = 1 first step re-polishes in ~1 iteration:
z0 = rec.ms.z(:);
K  = rec.ms.K;
[~, ~, Tf, Yf] = cr3bp_minenergy_prop(tf, [rv0; 1; z0], false, Tmax, c, muStar);
tG = linspace(0, tf, K+1);
Yj = interp1(Tf, Yf, tG, 'pchip')';
seed0 = struct('tf', tf, 'tGrid', tG, 'Y', Yj);
lg('RACE: cell (%d,%d) gamma %.2f, tf = %.4f ND (%.2f d), K = %d, sched %d steps', ...
   P.cell(1), P.cell(2), P.gamma, tf, tf*382981.289129055/86400, K, numel(P.sched));

pool = gcp;
lg('RACE: hard per-call timeout ACTIVE (%d workers)', pool.NumWorkers);

%% The two arms ------------------------------------------------------------
out = struct('cellIdx', P.cell, 'gamma', P.gamma, 'sched', P.sched, 'arms', struct());
for kf = 1:numel(P.families)
    fam = P.families{kf};
    lg('RACE ARM %s: start', fam);
    tArm = tic;
    seed = seed0;
    seed.Y(8:14,:) = P.seedScale * seed.Y(8:14,:);
    if P.seedScale ~= 1, lg('RACE ARM %s: seed costates scaled by %.3g', fam, P.seedScale); end
    A = struct('p', [], 'delta', [], 'mf', [], 'coastFrac', [], 'Hdrift', [], ...
               'iters', [], 'wall', [], 'nBisect', 0, 'nFail', 0, ...
               'retired', false, 'z', [], 'Y', {{}}, ...
               'tight', logical([]), 'condJ', [], ...
               'failP', [], 'failDelta', [], 'failCondJ', [], 'failNormR', []);
    % delta rule for the p-walk (huberc only; NaN = field default delta = p)
    if isempty(P.delta),               deltaOf = @(p) NaN(size(p));
    elseif isa(P.delta, 'function_handle'), deltaOf = P.delta;
    else,                              deltaOf = @(p) P.delta*ones(size(p));
    end
    % Rungs are (p, delta) pairs; stage 1 walks p (delta from the rule),
    % stage 2 (optional) holds p at its deepest value and walks delta.
    for stage = 1:2
    if stage == 1
        queue = [P.sched; deltaOf(P.sched)];  mv = 1;
    else
        if isempty(P.deltaSched) || isempty(A.p) || ~strcmp(fam, 'huberc'), break, end
        dLast = A.delta(end);  if isnan(dLast), dLast = A.p(end); end
        ds = P.deltaSched(P.deltaSched < dLast);
        if isempty(ds), break, end
        queue = [A.p(end)*ones(1, numel(ds)); ds];  mv = 2;
        lg('RACE ARM %s: STAGE 2 -- p held at %.4g, walking delta %s', fam, A.p(end), mat2str(ds, 3));
    end
    xGood = [];  gaps = 0;
    while ~isempty(queue)
        pTry = queue(1,1);  dTry = queue(2,1);  queue(:,1) = [];
        sm = struct('family', fam, 'p', pTry);
        if ~isnan(dTry), sm.delta = dTry; end
        xTry = pTry;  if mv == 2, xTry = dTry; end
        t0 = tic;
        [okRun, z, it] = run_capped(pool, @ms_minfuel, 2, P.wallSec + 90, ...
            rv0, rvf, tf, seed, Tmax, c, muStar, sm, ...
            struct('wallSec', P.wallSec));
        tight = okRun && it.converged && it.Hdrift < P.HdriftTol;
        loose = okRun && ~isempty(P.rungTolR) && ~tight && ...
                it.normR < P.rungTolR && it.Hdrift < dflt(P.rungHdriftTol, P.HdriftTol);
        okStep = tight || loose;
        cj = NaN;  if okRun && isfield(it, 'condJ'), cj = it.condJ; end
        if okStep
            A.tight(end+1) = tight;  A.condJ(end+1) = cj;
            % final mass from this arm's own flight of its own solution:
            [~, ~, ~, Yfly] = cr3bp_minfuel_prop(tf, it.Y(:,1), false, ...
                                                 Tmax, c, muStar, sm);
            mf = Yfly(end, 7);
            A.p(end+1) = pTry;  A.delta(end+1) = dTry;  A.mf(end+1) = mf;
            A.coastFrac(end+1) = it.coastFrac;  A.Hdrift(end+1) = it.Hdrift;
            A.iters(end+1) = it.iters;  A.wall(end+1) = toc(t0);
            A.z = z;  A.Y{end+1} = it.Y;
            yT2 = cr3bp_minfuel_prop(it.tGrid(end) - it.tGrid(end-1), ...
                                     it.Y(:,end), false, Tmax, c, muStar, sm);
            seed = struct('tf', tf, 'tGrid', it.tGrid(:)', 'Y', [it.Y, yT2]);
            xGood = xTry;
            lg('  %s p=%.4g%s %s: mf=%.6f coast=%.2f normR=%.1e Hdrift=%.1e condJ=%.2e iters=%d (%.0fs)', ...
               fam, pTry, dtag(dTry), tern(tight, 'OK', 'ok-LOOSE'), mf, it.coastFrac, it.normR, ...
               it.Hdrift, cj, it.iters, toc(t0));
        else
            A.nFail = A.nFail + 1;
            A.failP(end+1) = pTry;  A.failDelta(end+1) = dTry;  A.failCondJ(end+1) = cj;
            A.failNormR(end+1) = NaN;
            if okRun
                A.failNormR(end) = it.normR;
                lg('  %s p=%.4g%s FAIL: conv=%d normR=%.1e Hdrift=%.1e condJ=%.2e (%.0fs)', ...
                   fam, pTry, dtag(dTry), it.converged, it.normR, it.Hdrift, cj, toc(t0));
            else
                lg('  %s p=%.4g%s FAIL: HARD TIMEOUT/worker error (%.0fs)', ...
                   fam, pTry, dtag(dTry), toc(t0));
            end
            if ~isempty(xGood) && A.nBisect < P.maxBisect*numel(P.sched) && ...
               abs(xGood - xTry)/xGood > 0.02
                xMid = sqrt(xGood*xTry);          % geometric midpoint on the moving coordinate
                if mv == 1, ins = [xMid, pTry; deltaOf(xMid), dTry];
                else,       ins = [pTry, pTry; xMid, dTry];
                end
                queue = [ins, queue]; %#ok<AGROW>
                A.nBisect = A.nBisect + 1;
                lg('  %s bisect -> %s=%.4g', fam, tern(mv == 1, 'p', 'delta'), xMid);
            else
                gaps = gaps + 1;
                lg('  %s gap ABANDONED at p=%.4g%s (gap %d/%d)', ...
                   fam, pTry, dtag(dTry), gaps, P.maxGaps);
                if gaps >= P.maxGaps, A.retired = true; break, end
            end
        end
        out.arms.(fam) = A;  save(outMat, 'out');   % after EVERY step
    end
    if A.retired && stage == 1 && ~isempty(P.deltaSched)
        A.retired = false;                          % stage 2 still gets its chance
    end
    end                                             % stage loop
    % FLOOR REFINEMENT (MfMax Fullpath pattern): if the deepest accepted rung
    % was a LOOSE one, re-solve it at full tolerance from its own junctions;
    % a loose point that cannot be tightened is dropped back to the deepest
    % tight rung. Loose rungs are never reported as the arm's answer.
    if ~isempty(A.p) && ~A.tight(end)
        sm = smOf(fam, A.p(end), A.delta(end));
        yT2 = cr3bp_minfuel_prop(seed.tGrid(end) - seed.tGrid(end-1), ...
                                 A.Y{end}(:,end), false, Tmax, c, muStar, sm);
        seedR = struct('tf', tf, 'tGrid', seed.tGrid, 'Y', [A.Y{end}, yT2]);
        [okR, zR, itR] = run_capped(pool, @ms_minfuel, 2, P.wallSec + 90, ...
            rv0, rvf, tf, seedR, Tmax, c, muStar, sm, struct('wallSec', P.wallSec));
        if okR && itR.converged && itR.Hdrift < P.HdriftTol
            [~, ~, ~, Yfly] = cr3bp_minfuel_prop(tf, itR.Y(:,1), false, Tmax, c, muStar, sm);
            A.mf(end) = Yfly(end, 7);  A.z = zR;  A.Y{end} = itR.Y;
            A.tight(end) = true;  A.Hdrift(end) = itR.Hdrift;  A.condJ(end) = itR.condJ;
            lg('  %s FLOOR REFINED at p=%.4g: mf=%.6f normR=%.1e (tight)', fam, A.p(end), A.mf(end), itR.normR);
        else
            kT = find(A.tight, 1, 'last');
            lg('  %s floor p=%.4g could NOT be tightened (normR %.1e); falling back to deepest TIGHT rung p=%s', ...
               fam, A.p(end), tern(okR, itR.normR, NaN), tern(isempty(kT), 'NONE', sprintf('%.4g', A.p(max(kT,1)))));
            if isempty(kT)
                A.p = []; A.mf = []; A.coastFrac = []; A.Hdrift = []; A.iters = []; A.wall = [];
                A.tight = logical([]); A.condJ = []; A.Y = {}; A.z = [];
            else
                A.p = A.p(1:kT); A.mf = A.mf(1:kT); A.coastFrac = A.coastFrac(1:kT);
                A.Hdrift = A.Hdrift(1:kT); A.iters = A.iters(1:kT); A.wall = A.wall(1:kT);
                A.tight = A.tight(1:kT); A.condJ = A.condJ(1:kT); A.Y = A.Y(1:kT);
            end
        end
    end
    % acceptance gate at the arm's deepest converged parameter:
    if ~isempty(A.p)
        sm = smOf(fam, A.p(end), A.delta(end));
        [okA, ~, itA] = run_capped(pool, @ms_minfuel, 2, 2*P.wallSec, ...
            rv0, rvf, tf, seed, Tmax, c, muStar, sm, ...
            struct('wallSec', P.wallSec, 'accept', true));
        if okA && isfield(itA, 'accept')
            A.acceptDz = itA.accept.dz;  A.acceptOk = itA.accept.accepted;
        else
            A.acceptDz = NaN;  A.acceptOk = false;
        end
        lg('RACE ARM %s DONE: deepest p=%.4g%s, mf=%.6f, coast=%.2f, accept ok=%d dz=%.1e, %d fails/%d bisects, %.1f min', ...
           fam, A.p(end), dtag(A.delta(end)), A.mf(end), A.coastFrac(end), A.acceptOk, ...
           A.acceptDz, A.nFail, A.nBisect, toc(tArm)/60);
    else
        lg('RACE ARM %s DONE: NO converged steps', fam);
    end
    A.wallTotal = toc(tArm);
    out.arms.(fam) = A;  save(outMat, 'out');
end

%% Compare -----------------------------------------------------------------
if isfield(out.arms, 'eps') && isfield(out.arms, 'huber') && ...
   ~isempty(out.arms.eps.p) && ~isempty(out.arms.huber.p)
    out.compare = struct( ...
        'deepestP',  [out.arms.eps.p(end), out.arms.huber.p(end)], ...
        'mfFinal',   [out.arms.eps.mf(end), out.arms.huber.mf(end)], ...
        'mfAgree',   abs(out.arms.eps.mf(end) - out.arms.huber.mf(end)), ...
        'fails',     [out.arms.eps.nFail, out.arms.huber.nFail], ...
        'bisects',   [out.arms.eps.nBisect, out.arms.huber.nBisect], ...
        'wallMin',   [out.arms.eps.wallTotal, out.arms.huber.wallTotal]/60);
    lg('RACE COMPARE: deepest p [%.4g %.4g], mf [%.6f %.6f] (|d| %.2e), fails [%d %d], bisects [%d %d], wall [%.1f %.1f] min', ...
       out.compare.deepestP, out.compare.mfFinal, out.compare.mfAgree, ...
       out.compare.fails, out.compare.bisects, out.compare.wallMin);
end
save(outMat, 'out');
end

% ------------------------------------------------------------------------
function v = dflt(x, d)
% DFLT  x, or d when x is empty.  INPUTS: x; d.  OUTPUTS: v.
if isempty(x), v = d; else, v = x; end
end

function sm = smOf(fam, p, delta)
% SMOF  Smoothing spec for a (family, p, delta) rung; NaN delta = field
% default.  INPUTS: fam; p; delta.  OUTPUTS: sm struct.
sm = struct('family', fam, 'p', p);
if ~isnan(delta), sm.delta = delta; end
end

function s = dtag(delta)
% DTAG  ' delta=...' log tag, empty for the field default.  INPUTS: delta.
% OUTPUTS: s char.
if isnan(delta), s = ''; else, s = sprintf(' delta=%.4g', delta); end
end

function s = tern(c, a, b)
% TERN  a if c else b.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end

function logmsg(f, s)
% LOGMSG  Append to log file or stdout.  INPUTS: f path; s.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f,'a'); fprintf(fid,'%s\n',s); fclose(fid);
end
end
