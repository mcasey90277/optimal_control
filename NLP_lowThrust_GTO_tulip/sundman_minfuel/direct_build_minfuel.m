function direct_build_minfuel(factor, outFile, maxIter)
% DIRECT_BUILD_MINFUEL  Solve the min-fuel problem at one t_f FROM SCRATCH, with
% no continuation from any other t_f.
%
% Builds a fresh burn+coast warm start AT this t_f (min-time burn to t_f^min,
% then coast for the remaining t_f - t_f^min), maps it into Sundman coordinates
% (no-resample) constrained to the REAL endpoints, solves min-energy, then the
% energy->fuel homotopy to min-fuel. Because the warm start is built at this
% t_f, it is dynamically consistent -- unlike rescaling a solution from another
% t_f. The quality of the warm start scales with how SHORT the coast is: near
% min-time (small factor) the coast is tiny and the burn+coast is nearly the
% exact solution, so this is expected to work best close to min-time and to
% degrade (hit the restoration wall) as the coast grows.
%
% INPUTS: factor - t_f/t_f^min ; outFile - output .mat ; maxIter [default 2000]
% OUTPUT FILE: out (full solver struct incl. X,U,lamDef), factor, tf, dV,...

here=fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','..','lowThrust_GTO_tulip'));   % lt_pmp_eom
addpath(fullfile(here,'..'));                              % lt_dynamics_dirthrottle
addpath(here);                                             % sundman_minfuel LAST (priority)
if nargin<3||isempty(maxIter), maxIter=2000; end
pSund=1.5;
p=cr3bp_lt_params(0.025,15,2100); Tmax=p.Tmax; c=p.c; muStar=p.muStar;
S=load(fullfile(here,'minfuel_from_energy_seed.mat')); rv0=S.rv0; rvf=S.rvf; tfMin=S.tf/1.15;
tf=factor*tfMin;
zMin=[190.4760481;-79.7060409;-0.4298691037;0.3011592775;0.5866700046;-0.007117348902;4.329378839];
opt=odeset('RelTol',1e-12,'AbsTol',1e-14);
[tB,yB]=ode113(@lt_pmp_eom,[0 tfMin],[rv0(:);1;zMin],opt,Tmax,c,muStar); [tB,k]=unique(tB,'stable'); yB=yB(k,:);
cf=@(t,x) lt_dynamics_dirthrottle(x,[1;0;0;0],Tmax,c,muStar);
[tC,yC]=ode113(cf,[0, tf-tfMin], yB(end,1:7).', opt); [tC,k]=unique(tC,'stable'); yC=yC(k,:);
xAll=[yB(:,1:7); yC(2:end,:)]; tAll=[tB; tfMin+tC(2:end)];
sAll=[ones(numel(tB),1); zeros(numel(tC)-1,1)];
alB=-yB(:,11:13)./sqrt(sum(yB(:,11:13).^2,2)); alC=repmat([1 0 0],numel(tC)-1,1); alAll=[alB;alC];
fprintf('DIRECT %.3fx: burn+coast %d nodes, coast=%.1f%% of tf\n', factor, numel(tAll), 100*(tf-tfMin)/tf);
[sigma,X0,U0,tauf0]=sundman_seed_map(xAll.', [alAll.';sAll.'], tf, tAll, pSund, muStar, rv0, rvf);
oE=casadi_minfuel_sundman(sigma,tf,rv0,rvf,Tmax,c,muStar,X0,U0,tauf0,pSund,maxIter,1,false);
fprintf('DIRECT %.3fx ENERGY: ok=%d defect=%.2g edge=%.1f%%\n', factor, oE.success, oE.maxDefect, 100*oE.edge);
if ~(oE.success && oE.maxDefect<1e-6)
    fprintf('DIRECT %.3fx: energy did not converge -- not sharpening.\n', factor); return;
end
oT=casadi_minfuel_sundman(sigma,tf,rv0,rvf,Tmax,c,muStar,oE.X,oE.U,tauf0,pSund,maxIter,1,true);
Xk=oT.X; Uk=oT.U; best=oT;
for e=[0.6 0.35 0.2 0.12 0.07 0.04 0.025 0.015 0.008 0.004 0.002 0.001]
    o=casadi_minfuel_sundman(sigma,tf,rv0,rvf,Tmax,c,muStar,Xk,Uk,tauf0,pSund,maxIter,e,true);
    if o.success && o.maxDefect<1e-6, Xk=o.X; Uk=o.U; best=o; end
end
dV=c*log(1/best.mf)*p.lStar/p.tStar;
out=best; out.factor=factor; out.tf=tf; out.tf_days=tf*p.tStar/86400; out.dV=dV; out.prop_kg=p.m0kg*(1-best.mf);
save(outFile,'out','sigma','tauf0','rv0','rvf','pSund','factor');
fprintf('DIRECT %.3fx MIN-FUEL: dV=%.4f sw=%d edge=%.1f%% def=%.2g primer=%.3f -> %s\n', ...
        factor, dV, best.switches, 100*best.edge, best.maxDefect, best.primerAlignDeg, outFile);
end
