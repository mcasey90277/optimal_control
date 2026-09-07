function R = run_gamma_walk(cfg)
% RUN_GAMMA_WALK  Boundary-data homotopy into the HIGH-GAMMA band: walk the
% fixed-tf MIN-ENERGY problem (p = 1, smooth) in tf = gamma * tfMin from the
% deepest converged pilot record of a cell up to the gammas where the direct
% collocation solve fails (FINDINGS 19: gamma >= 1.7 on (2,5), >= 1.4 on
% (6,8) and (1,2)). Each step is seeded by the previous solution's junction
% states at the SAME FRACTIONAL TIMES (the K+1 breakpoints of the new,
% longer grid), solved by ms_minenergy, gated by convergence + absolute H
% conservation, with geometric bisection on gamma. This is the MfMax homCI
% idea (continuation in the boundary data, no direct solve anywhere) in the
% form this program already trusts (FINDINGS 21's gamma-continuation).
%
% Records are written in the minenergy_pilot.mat record format so
% run_minfuel_race can consume them unchanged (cfg.pilotMat).
%
% INPUTS:
%   cfg - (optional) struct: .cells [3x2: (2,5);(6,8);(1,2)], .gammas
%         [1.4 1.55 1.7 1.85 2.0] (targets; each cell starts above its own
%         deepest PASS record), .maxBisect [4], .wallSec [600], .logFile
%         [''], .outMat [direct/results/minenergy_highgamma.mat]
%
% OUTPUTS:
%   R - struct array of records (pilot format: iD iA gam tfMin tf Tmax c
%       muStar rv0 rvf direct.mf ms.z ms.K ms.converged ms.normR ms.Hdrift
%       verdict 'PASS'|'FAIL', plus .gammaFrom, .conjPass, .Y)
%
% REFERENCES:
%   [1] doc/mfmax_ideas_review.md section 1.5 (the seed route that needs no
%       direct solve); FINDINGS 21 (gamma-continuation at fixed depth).
%   [2] run_minenergy_pilot.m (record format), ms_minenergy.m (the solve).

if nargin < 1, cfg = struct(); end
P.cells     = [2 5; 6 8; 1 2];
P.gammas    = [1.4 1.55 1.7 1.85 2.0];
P.maxBisect = 4;
P.wallSec   = 600;
P.logFile   = '';
P.outMat    = '';
fn = fieldnames(cfg);  for k = 1:numel(fn), P.(fn{k}) = cfg.(fn{k}); end
lg = @(varargin) logmsg(P.logFile, sprintf(varargin{:}));

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'indirect'), fullfile(fileparts(here), 'costate_common'));
resDir = fullfile(here, 'direct', 'results');
if isempty(P.outMat), P.outMat = fullfile(resDir, 'minenergy_highgamma.mat'); end
Lp = load(fullfile(resDir, 'minenergy_pilot.mat'));

R = struct([]);
for kc = 1:size(P.cells, 1)
    cell_ = P.cells(kc, :);
    % deepest PASS pilot record of this cell = the walk's start
    ok = arrayfun(@(r) isequal([r.iD r.iA], cell_) && ischar(r.verdict) && strcmp(r.verdict, 'PASS'), Lp.R);
    cand = Lp.R(ok);  [~, im] = max([cand.gam]);  rec = cand(im);
    K = rec.ms.K;
    [~, ~, Tf, Yf] = cr3bp_minenergy_prop(rec.tf, [rec.rv0(:); 1; rec.ms.z(:)], false, rec.Tmax, rec.c, rec.muStar);
    tG = linspace(0, rec.tf, K+1);
    Yj = interp1(Tf, Yf, tG, 'pchip')';                       % [14 x K+1] at fractional times
    gamGood = rec.gam;
    lg('GAMMA WALK (%d,%d): start gamma %.2f (tf %.4f, K %d) -> targets %s', ...
       cell_, gamGood, rec.tf, K, mat2str(P.gammas(P.gammas > gamGood), 3));
    queue = P.gammas(P.gammas > gamGood);  nB = 0;
    while ~isempty(queue)
        gamTry = queue(1);  queue(1) = [];
        tf = gamTry * rec.tfMin;
        seed = struct('tf', tf, 'tGrid', linspace(0, tf, K+1), 'Y', Yj);
        t0 = tic;
        try
            [z, it] = ms_minenergy(rec.rv0(:), rec.rvf(:), tf, seed, rec.Tmax, rec.c, rec.muStar, ...
                                   struct('conjTest', true, 'wallSec', P.wallSec));
            okStep = it.converged && it.Hdrift < 1e-6;
        catch ME
            okStep = false;  it = struct('normR', NaN, 'Hdrift', NaN);  z = [];
            lg('  gamma %.3f ERROR: %s', gamTry, ME.message);
        end
        if okStep
            [~, ~, ~, Yfly] = cr3bp_minenergy_prop(tf, it.Y(:,1), false, rec.Tmax, rec.c, rec.muStar);
            yT = cr3bp_minenergy_prop(it.tGrid(end) - it.tGrid(end-1), it.Y(:,end), false, rec.Tmax, rec.c, rec.muStar);
            Yj = [it.Y, yT];  gamGood = gamTry;
            r = struct('iD', cell_(1), 'iA', cell_(2), 'gam', gamTry, 'tfMin', rec.tfMin, 'tf', tf, ...
                'Tmax', rec.Tmax, 'c', rec.c, 'muStar', rec.muStar, 'rv0', rec.rv0, 'rvf', rec.rvf, ...
                'direct', struct('mf', Yfly(end, 7), 'source', 'gamma-walk (no direct solve)'), ...
                'ms', struct('z', z(:), 'K', K, 'converged', it.converged, 'normR', it.normR, ...
                             'Hdrift', it.Hdrift, 'iters', it.iters, 'condJ', it.condJ), ...
                'verdict', 'PASS', 'gammaFrom', gamGood, 'conjPass', double(it.conj.pass), ...
                'conjVerdict', it.conj.verdict, 'Y', it.Y);
            if isempty(R), R = r; else, R(end+1) = r; end %#ok<AGROW>
            save(P.outMat, 'R');
            lg('  gamma %.3f OK: mf=%.6f normR=%.1e Hdrift=%.1e iters=%d condJ=%.2e conj=%s (%.0fs)', ...
               gamTry, Yfly(end,7), it.normR, it.Hdrift, it.iters, it.condJ, it.conj.verdict, toc(t0));
        else
            lg('  gamma %.3f FAIL: normR=%.1e Hdrift=%.1e (%.0fs)', gamTry, it.normR, it.Hdrift, toc(t0));
            if nB < P.maxBisect && (gamTry - gamGood)/gamGood > 0.01
                gMid = sqrt(gamGood * gamTry);
                queue = [gMid, gamTry, queue];  nB = nB + 1;             %#ok<AGROW>
                lg('  bisect -> gamma %.4f', gMid);
            else
                lg('  gamma %.3f ABANDONED (%d bisections); walk for (%d,%d) stops at %.3f', gamTry, nB, cell_, gamGood);
                break
            end
        end
    end
end
lg('GAMMA WALK DONE: %d records -> %s', numel(R), P.outMat);
end

% ------------------------------------------------------------------------
function logmsg(f, s)
% LOGMSG  Append to log file or stdout.  INPUTS: f; s.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
end
