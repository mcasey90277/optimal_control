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
%   best - best (smallest-eps, tight) solver struct from CASADI_MINFUEL_SUNDMAN
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
    ok = out.success && out.maxDefect < 1e-6;
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
if isempty(best), best = out; bestEps = epsSched(end); end
dVb = p.c*log(1/best.mf)*p.lStar/p.tStar;
fprintf('\n=== HOMOTOPY BEST: eps=%.4g defect=%.2g switches=%d edge=%.1f%% prop=%.4f kg dV=%.4f km/s ===\n', ...
        bestEps, best.maxDefect, best.switches, 100*best.edge, p.m0kg*(1-best.mf), dVb);
end
