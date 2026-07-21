function [best, tbl] = sundman_homotopy(p, rv0, rvf, sigma, X0, U0, tauf0, pSund, epsSched, maxIter, saveFile)
% SUNDMAN_HOMOTOPY  Guarded energy->fuel homotopy sweep over a Sundman-
% regularized min-fuel collocation problem.
%
% Solves CASADI_MINFUEL_SUNDMAN at each epsilon in a decreasing schedule,
% warm-starting each solve from the previous one, so the throttle sharpens
% continuously from the smooth energy control (eps=1) to bang-bang fuel
% (eps=0). GUARDED: a step that does not converge tight (success=0 or
% maxDefect>1e-6) is DISCARDED -- it neither advances the warm start nor
% overwrites the best-so-far solution -- so a wedged near-zero-eps step cannot
% poison the chain or clobber a good result. This single routine replaces the
% earlier run_sundman_homotopy + run_sundman_tail pair.
%
% INPUTS:
%   p        - parameter struct (CR3BP_LT_PARAMS); uses .c,.Tmax,.muStar,.m0kg,
%              .lStar,.tStar
%   rv0,rvf  - departure / arrival states [1x6]
%   sigma    - Sundman nodes [ (N+1) x1 ] (from SUNDMAN_SEED_MAP)
%   X0,U0    - warm-start states/controls [8x(N+1)] / [4x(N+1)]
%   tauf0    - fixed regularized length [scalar]
%   pSund    - Sundman power [scalar]
%   epsSched - decreasing homotopy schedule, e.g. [1 .6 .35 ... .001 0]
%   maxIter  - IPOPT max iters per step [scalar, default 1500]
%   saveFile - (optional) .mat path; the best solution is saved after each
%              improving step
%
% OUTPUTS:
%   best - best (smallest-eps, tight) solver struct from CASADI_MINFUEL_SUNDMAN,
%          carrying a .certified logical: true iff at least one schedule step
%          converged tight (success & maxDefect<1e-6). When false, `best` is the
%          last UNCERTIFIED iterate and must not be treated as a solution.
%   tbl  - per-step table, columns [eps defect switches edge% prop_kg dV_kms]

if nargin < 10 || isempty(maxIter), maxIter = 1500; end
if nargin < 11, saveFile = ''; end

tf = X0(8,end);                          % transfer time, pinned by the seed map
tbl = zeros(numel(epsSched), 6);
best = []; bestEps = NaN;  out = [];
for ie = 1:numel(epsSched)
    epsH = epsSched(ie);
    fprintf('\n---- eps=%.4g (step %d/%d) ----\n', epsH, ie, numel(epsSched));
    out = casadi_minfuel_sundman(sigma, tf, rv0, rvf, p.Tmax, p.c, p.muStar, ...
                                 X0, U0, tauf0, pSund, maxIter, epsH);
    dV = p.c*log(1/out.mf)*p.lStar/p.tStar;
    tbl(ie,:) = [epsH, out.maxDefect, out.switches, 100*out.edge, p.m0kg*(1-out.mf), dV];
    % 2026-07-21 triage C1: require full convergence -- `success` alone also
    % covers Solved_To_Acceptable_Level, whose 1e-5-grade duals must not enter
    % the warm chain or PMP checks.
    ok = strcmp(out.ipoptStatus,'Solve_Succeeded') && out.maxDefect < 1e-6;
    verdict = 'DISCARD'; if ok, verdict = 'KEEP'; end
    fprintf('eps=%.4g: success=%d defect=%.2g switches=%d edge=%.1f%% prop=%.4f kg dV=%.4f  %s\n', ...
            epsH, out.success, out.maxDefect, out.switches, 100*out.edge, p.m0kg*(1-out.mf), dV, verdict);
    if ok
        X0 = out.X; U0 = out.U; best = out; bestEps = epsH;
        if ~isempty(saveFile)
            eps = epsH; %#ok<NASGU>  (saved as 'eps' for .mat back-compat)
            save(saveFile, 'out','sigma','rv0','rvf','tauf0','pSund','eps','tbl');
        end
    else
        fprintf('   (loose step: warm start not advanced; best kept at eps=%.4g)\n', bestEps);
    end
end
anyClean = ~isempty(best);
if ~anyClean
    warning('sundman_homotopy:noCleanStep', ...
        ['no homotopy step converged tight (Solve_Succeeded & maxDefect<1e-6); ' ...
         'returning the last UNCERTIFIED iterate -- do NOT treat it as a solution']);
    best = out;  bestEps = epsSched(end);
end
% 2026-07-21 triage C2: `certified` now requires the REQUESTED homotopy
% endpoint (last scheduled eps), not merely "some clean step" -- an eps>0
% best must not carry the min-fuel certification flag.
best.epsReached = bestEps;
best.certified  = anyClean && abs(bestEps - epsSched(end)) < 1e-12;
certified = best.certified;
if anyClean && ~certified
    warning('sundman_homotopy:endpointNotReached', ...
        ['homotopy stalled at eps=%.4g (requested endpoint eps=%.4g): result is ' ...
         'a clean INTERMEDIATE solution, NOT certified at the requested objective'], ...
        bestEps, epsSched(end));
end
dVb = p.c*log(1/best.mf)*p.lStar/p.tStar;
tag = 'BEST'; if ~certified, tag = 'UNCERTIFIED (no tight step)'; end
fprintf('\n=== HOMOTOPY %s: eps=%.4g defect=%.2g switches=%d edge=%.1f%% prop=%.4f kg dV=%.4f km/s ===\n', ...
        tag, bestEps, best.maxDefect, best.switches, 100*best.edge, p.m0kg*(1-best.mf), dVb);
end
