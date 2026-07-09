function best = run_sundman_tail(startFile, epsSched, maxIter)
% RUN_SUNDMAN_TAIL  Guarded fine-step tail of the energy->fuel homotopy.
%
% Continues the Sundman-regularized homotopy from an already-converged Sundman
% solution (in Sundman coordinates, correct mesh -- no remap needed) toward
% eps=0, using small eps steps near the bang-bang limit. GUARDED: a step that
% fails to converge tight (success=0 or defect>1e-6) is discarded -- it neither
% becomes the next warm start nor overwrites the best-so-far solution -- so a
% wedged near-zero-eps step cannot poison the chain or clobber a good result.
%
% INPUTS:
%   startFile - .mat holding .out (Sundman solver output), .sigma, .tf, .rv0,
%               .rvf, .tauf0, .pSund  [default sundman_minfuel_certified.mat]
%   epsSched  - decreasing eps schedule to append below the start point
%   maxIter   - IPOPT max iters per step [default 2000]
%
% OUTPUTS:
%   best - best (smallest-eps, tight) solver struct reached; saved to
%          sundman_minfuel_certified.mat

here=fileparts(mfilename('fullpath')); addpath(here); run(fullfile(here,'setup_paths.m'));
if nargin<1||isempty(startFile), startFile=fullfile(here,'sundman_minfuel_certified.mat'); end
if nargin<2||isempty(epsSched), epsSched=[0.03 0.022 0.016 0.012 0.009 0.006 0.004 0.0025 0.0015 0.001 0.0005 0]; end
if nargin<3||isempty(maxIter), maxIter=2000; end
muStar=0.012150585609624; lStar=389703.264829278; tStar=382981.289129055;
m0kg=15; g0=9.80665*tStar^2/(1000*lStar);
Tmax=(0.025/m0kg)*tStar^2/(lStar*1000); c=(2100/tStar)*g0;

S=load(startFile); sigma=S.sigma; tf=S.tf; rv0=S.rv0; rvf=S.rvf; tauf0=S.tauf0; pSund=S.pSund;
X0=S.out.X; U0=S.out.U; best=S.out; bestEps=S.eps;
fprintf('TAIL start: eps=%.4g defect=%.2g switches=%d edge=%.1f%% prop=%.4f kg\n', ...
        bestEps, best.maxDefect, best.switches, 100*best.edge, m0kg*(1-best.mf));

for ie=1:numel(epsSched)
  epsH=epsSched(ie);                            % epsH, not eps (eps is builtin)
  fprintf('\n---- TAIL eps=%.4g (%d/%d) ----\n', epsH, ie, numel(epsSched));
  out=casadi_minfuel_sundman(sigma, tf, rv0, rvf, Tmax, c, muStar, X0, U0, tauf0, pSund, maxIter, epsH);
  dV=c*log(1/out.mf)*lStar/tStar;
  ok = out.success && out.maxDefect<1e-6;
  verdict = 'DISCARD'; if ok, verdict = 'KEEP'; end
  fprintf('eps=%.4g: success=%d defect=%.2g switches=%d edge=%.1f%% prop=%.4f kg dV=%.4f  %s\n', ...
          epsH, out.success, out.maxDefect, out.switches, 100*out.edge, m0kg*(1-out.mf), dV, verdict);
  if ok
    X0=out.X; U0=out.U; best=out; bestEps=epsH;     % advance only on success
    eps=epsH; %#ok<NASGU>  (saved as 'eps' for .mat back-compat with loaders)
    save(fullfile(here,'sundman_minfuel_certified.mat'), 'out','sigma','tf','rv0','rvf','tauf0','pSund','eps');
  else
    fprintf('   (kept eps=%.4g as best; not advancing)\n', bestEps);
  end
end
dVb=c*log(1/best.mf)*lStar/tStar;
fprintf('\n=== TAIL BEST: eps=%.4g defect=%.2g switches=%d edge=%.1f%% prop=%.4f kg dV=%.4f km/s ===\n', ...
        bestEps, best.maxDefect, best.switches, 100*best.edge, m0kg*(1-best.mf), dVb);
end
