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
P.logFile  = '';
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
outMat = fullfile(here, 'direct', 'results', 'minfuel_race.mat');
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
    A = struct('p', [], 'mf', [], 'coastFrac', [], 'Hdrift', [], ...
               'iters', [], 'wall', [], 'nBisect', 0, 'nFail', 0, ...
               'retired', false, 'z', [], 'Y', {{}});
    pGood = [];  gaps = 0;
    queue = P.sched;
    while ~isempty(queue)
        pTry = queue(1);  queue(1) = [];
        sm = struct('family', fam, 'p', pTry);
        t0 = tic;
        [okRun, z, it] = run_capped(pool, @ms_minfuel, 2, P.wallSec + 90, ...
            rv0, rvf, tf, seed, Tmax, c, muStar, sm, ...
            struct('wallSec', P.wallSec));
        okStep = okRun && it.converged && it.Hdrift < P.HdriftTol;
        if okStep
            % final mass from this arm's own flight of its own solution:
            [~, ~, ~, Yfly] = cr3bp_minfuel_prop(tf, it.Y(:,1), false, ...
                                                 Tmax, c, muStar, sm);
            mf = Yfly(end, 7);
            A.p(end+1) = pTry;  A.mf(end+1) = mf;
            A.coastFrac(end+1) = it.coastFrac;  A.Hdrift(end+1) = it.Hdrift;
            A.iters(end+1) = it.iters;  A.wall(end+1) = toc(t0);
            A.z = z;  A.Y{end+1} = it.Y;
            yT2 = cr3bp_minfuel_prop(it.tGrid(end) - it.tGrid(end-1), ...
                                     it.Y(:,end), false, Tmax, c, muStar, sm);
            seed = struct('tf', tf, 'tGrid', it.tGrid(:)', 'Y', [it.Y, yT2]);
            pGood = pTry;
            lg('  %s p=%.4g OK: mf=%.6f coast=%.2f Hdrift=%.1e iters=%d (%.0fs)', ...
               fam, pTry, mf, it.coastFrac, it.Hdrift, it.iters, toc(t0));
        else
            A.nFail = A.nFail + 1;
            if okRun
                lg('  %s p=%.4g FAIL: conv=%d normR=%.1e Hdrift=%.1e (%.0fs)', ...
                   fam, pTry, it.converged, it.normR, it.Hdrift, toc(t0));
            else
                lg('  %s p=%.4g FAIL: HARD TIMEOUT/worker error (%.0fs)', ...
                   fam, pTry, toc(t0));
            end
            if ~isempty(pGood) && A.nBisect < P.maxBisect*numel(P.sched) && ...
               abs(pGood - pTry)/pGood > 0.02
                pMid = sqrt(pGood*pTry);          % geometric midpoint
                queue = [pMid, pTry, queue]; %#ok<AGROW>
                A.nBisect = A.nBisect + 1;
                lg('  %s bisect -> %.4g', fam, pMid);
            else
                gaps = gaps + 1;
                lg('  %s gap ABANDONED at p=%.4g (gap %d/%d)', ...
                   fam, pTry, gaps, P.maxGaps);
                if gaps >= P.maxGaps, A.retired = true; break, end
            end
        end
        out.arms.(fam) = A;  save(outMat, 'out');   % after EVERY step
    end
    % acceptance gate at the arm's deepest converged parameter:
    if ~isempty(A.p)
        sm = struct('family', fam, 'p', A.p(end));
        [okA, ~, itA] = run_capped(pool, @ms_minfuel, 2, 2*P.wallSec, ...
            rv0, rvf, tf, seed, Tmax, c, muStar, sm, ...
            struct('wallSec', P.wallSec, 'accept', true));
        if okA && isfield(itA, 'accept')
            A.acceptDz = itA.accept.dz;  A.acceptOk = itA.accept.accepted;
        else
            A.acceptDz = NaN;  A.acceptOk = false;
        end
        lg('RACE ARM %s DONE: deepest p=%.4g, mf=%.6f, coast=%.2f, accept ok=%d dz=%.1e, %d fails/%d bisects, %.1f min', ...
           fam, A.p(end), A.mf(end), A.coastFrac(end), A.acceptOk, ...
           A.acceptDz, A.nFail, A.nBisect, toc(tArm)/60);
    else
        lg('RACE ARM %s DONE: NO converged steps', fam);
    end
    A.wallTotal = toc(tArm);
    out.arms.(fam) = A;  save(outMat, 'out');
end

%% Compare -----------------------------------------------------------------
if all(isfield(out.arms, P.families)) && ...
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
function logmsg(f, s)
% LOGMSG  Append to log file or stdout.  INPUTS: f path; s.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f,'a'); fprintf(fid,'%s\n',s); fclose(fid);
end
end
