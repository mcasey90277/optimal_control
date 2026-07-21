function out = test_sundman_eps1_noresample(pSund, epsVal, maxIter)
% Quick probe: map the certified min-fuel solution into Sundman tau using its
% OWN nodes (no pchip resample) and solve at a single eps. Tests whether a
% dynamically-minimal-error seed lets the energy (eps=1) problem converge tight
% in Sundman coordinates -- isolating "resample loss" from "conditioning".
if nargin<1||isempty(pSund), pSund=1.5; end
if nargin<2||isempty(epsVal), epsVal=1; end
if nargin<3||isempty(maxIter), maxIter=1500; end
here=fileparts(mfilename('fullpath')); addpath(here); run(fullfile(here,'setup_paths.m'));
muStar=0.012150585609624; lStar=389703.264829278; tStar=382981.289129055;
m0kg=15; g0=9.80665*tStar^2/(1000*lStar);
Tmax=(0.025/m0kg)*tStar^2/(lStar*1000); c=(2100/tStar)*g0;
S=load(fullfile(here,'minfuel_from_energy_seed.mat'));
Xs=S.nlp.X; Us=S.nlp.U; tf=S.tf; rv0=S.rv0; rvf=S.rvf;
sg=S.sigma(:); sg=(sg-sg(1))/(sg(end)-sg(1)); tSeed=sg*tf;
wcol=Us(1:3,:); s_seed=Us(4,:); alpha=wcol./max(sqrt(sum(wcol.^2,1)),1e-9);
r1=sqrt((Xs(1,:)+muStar).^2+Xs(2,:).^2+Xs(3,:).^2).'; kap=r1.^pSund;
dt=diff(tSeed); dtau=dt.*0.5.*(1./kap(1:end-1)+1./kap(2:end));
tau=[0;cumsum(dtau)]; tauf0=tau(end);
sigma=tau/tauf0;                              % seed's own tau nodes, no resample
[sigma,ku]=unique(sigma,'stable');
X0=[Xs(:,ku); tSeed(ku).']; U0=[alpha(:,ku); s_seed(ku)];
X0(1:6,1)=rv0(:); X0(7,1)=1; X0(8,1)=0; X0(1:6,end)=rvf(:); X0(8,end)=tf;
fprintf('NORESAMPLE eps=%.3g: N=%d pSund=%.2f tauf0=%.4g\n', epsVal, numel(sigma)-1, pSund, tauf0);
out=casadi_minfuel_sundman(sigma, tf, rv0, rvf, Tmax, c, muStar, X0, U0, tauf0, pSund, maxIter, epsVal);
dV=c*log(1/out.mf)*lStar/tStar;
fprintf('\n=== NORESAMPLE RESULT eps=%.3g ===\n', epsVal);
fprintf('success=%d defect=%.2g unit=%.2g switches=%d edge=%.1f%% prop=%.4f kg dV=%.4f\n', ...
        out.success, out.maxDefect, out.maxUnit, out.switches, 100*out.edge, m0kg*(1-out.mf), dV);
end
