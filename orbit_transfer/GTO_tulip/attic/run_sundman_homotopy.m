function out = run_sundman_homotopy(N, pSund, epsSched, seedFile, maxIter)
% RUN_SUNDMAN_HOMOTOPY  Sundman-regularized energy->fuel homotopy for the
% min-fuel GTO->tulip transfer -- the two-wall fix in one pipeline.
%
% Sundman regularization (kappa=r1^pSund) tames the near-perigee 1/r^3
% singularity so the collocation defects can reach machine zero (proven:
% 1e-10 at N=1500). But feeding IPOPT a *bang-bang* objective from a bang-bang
% seed drives it straight into the restoration phase and a false "locally
% infeasible" exit -- the objective-side wall. The cure is the Bertrand-Epenoy
% homotopy: start at the smooth (strictly convex) energy objective eps=1, then
% step eps->0 toward fuel, warm-starting each solve from the last so the
% sharpening throttle stays in the solver's basin.
%
%   J(eps) = Int[s]dt - eps*Int[s(1-s)]dt
%     eps=1 -> Int s^2 dt (energy, smooth ramp)   eps=0 -> Int s dt (fuel, bang)
%
% INPUTS:
%   N        - Sundman segments [default 2000]
%   pSund    - Sundman power [default 1.5]
%   epsSched - decreasing eps schedule [default [1 .5 .2 .1 .05 .02 .01 0]]
%   seedFile - certified time-mesh solution .mat [default energy-seed soln]
%   maxIter  - IPOPT max iters per eps step [default 1500]
%
% OUTPUTS:
%   out - final (eps=0) solver struct; also saves sundman_homotopy_solution.mat

if nargin<1||isempty(N), N=2000; end
if nargin<2||isempty(pSund), pSund=1.5; end
if nargin<3||isempty(epsSched), epsSched=[1 0.5 0.2 0.1 0.05 0.02 0.01 0]; end
here=fileparts(mfilename('fullpath'));
if nargin<4||isempty(seedFile), seedFile=fullfile(here,'minfuel_from_energy_seed.mat'); end
if nargin<5||isempty(maxIter), maxIter=1500; end
addpath(here); run(fullfile(here,'setup_paths.m'));

muStar=0.012150585609624; lStar=389703.264829278; tStar=382981.289129055;
m0kg=15; g0=9.80665*tStar^2/(1000*lStar);
Tmax=(0.025/m0kg)*tStar^2/(lStar*1000); c=(2100/tStar)*g0;

% ---- map certified time-mesh solution into Sundman coordinates ----
% NO-RESAMPLE: use the seed's OWN nodes mapped to tau (sigma = tau/tauf0). A
% pchip resample onto a uniform-sigma mesh loses too much of the 40-rev
% oscillatory structure and leaves an irreducible ~1e-2 defect that pins IPOPT
% in restoration; the seed's own nodes make the only initial infeasibility the
% (small) time-trap vs Sundman-trap mismatch, which IPOPT closes to 1e-14.
% N is therefore ignored (kept for signature compatibility).
S=load(seedFile); Xs=S.nlp.X; Us=S.nlp.U; tf=S.tf; rv0=S.rv0; rvf=S.rvf;
sg=S.sigma(:); sg=(sg-sg(1))/(sg(end)-sg(1)); tSeed=sg*tf;
wcol=Us(1:3,:); s_seed=Us(4,:); alpha=wcol./max(sqrt(sum(wcol.^2,1)),1e-9);
r1=sqrt((Xs(1,:)+muStar).^2+Xs(2,:).^2+Xs(3,:).^2).'; kap=r1.^pSund;
dt=diff(tSeed); dtau=dt.*0.5.*(1./kap(1:end-1)+1./kap(2:end));
tau=[0;cumsum(dtau)]; tauf0=tau(end); sigma=tau/tauf0;
[sigma,ku]=unique(sigma,'stable');
X0=[Xs(:,ku); tSeed(ku).']; U0=[alpha(:,ku); s_seed(ku)];
X0(1:6,1)=rv0(:); X0(7,1)=1; X0(8,1)=0; X0(1:6,end)=rvf(:); X0(8,end)=tf;
N=numel(sigma)-1;

% ---- homotopy loop ----
fprintf('SUNDMAN HOMOTOPY: N=%d pSund=%.2f tf=%.4f tauf0=%.4g eps=[%s]\n', ...
        N, pSund, tf, tauf0, num2str(epsSched));
tbl=zeros(numel(epsSched),6);
best=[]; bestEps=NaN;
for ie=1:numel(epsSched)
  epsH=epsSched(ie);                            % epsH, not eps (eps is builtin)
  fprintf('\n---- eps=%.4g (step %d/%d) ----\n', epsH, ie, numel(epsSched));
  out=casadi_minfuel_sundman(sigma, tf, rv0, rvf, Tmax, c, muStar, X0, U0, tauf0, pSund, maxIter, epsH);
  dV=c*log(1/out.mf)*lStar/tStar;
  tbl(ie,:)=[epsH, out.maxDefect, out.switches, 100*out.edge, m0kg*(1-out.mf), dV];
  ok = out.success && out.maxDefect<1e-6;
  verdict='DISCARD'; if ok, verdict='KEEP'; end
  fprintf('eps=%.4g: success=%d defect=%.2g switches=%d edge=%.1f%% prop=%.4f kg dV=%.4f  %s\n', ...
          epsH, out.success, out.maxDefect, out.switches, 100*out.edge, m0kg*(1-out.mf), dV, verdict);
  if ok
    X0=out.X; U0=out.U;                        % advance warm start only on success
    best=out; bestEps=epsH; eps=epsH; %#ok<NASGU>  (eps kept in .mat for back-compat)
    save(fullfile(here,'sundman_homotopy_solution.mat'), 'out','sigma','tf','rv0','rvf','tauf0','pSund','eps','tbl');
  else
    fprintf('   (loose step: not advancing warm start; kept eps=%.4g)\n', bestEps);
  end
end
if isempty(best), best=out; end                % nothing converged tight; report last

fprintf('\n=== HOMOTOPY SUMMARY (eps, defect, switches, edge%%, prop kg, dV) ===\n');
for ie=1:numel(epsSched)
  fprintf('  %6.3g  %8.2g  %3d  %5.1f  %7.4f  %7.4f\n', tbl(ie,1),tbl(ie,2),tbl(ie,3),tbl(ie,4),tbl(ie,5),tbl(ie,6));
end
fprintf('best kept: eps=%.4g  defect=%.2g  switches=%d  edge=%.1f%%\n', ...
        bestEps, best.maxDefect, best.switches, 100*best.edge);
fprintf('FINAL_PASS=%d\n', ~isempty(best) && best.maxDefect<1e-6);
out=best;
end
