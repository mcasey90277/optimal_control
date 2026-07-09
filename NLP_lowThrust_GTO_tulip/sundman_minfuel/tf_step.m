function tf_step(seedMat, targetFactor, outMat)
% TF_STEP  One t_f-continuation step, in an ISOLATED MATLAB process.
%
% Loads a bang-bang solution from seedMat, rescales its time state to
% t_f = targetFactor * t_f^min, re-sharpens (fine schedule), and saves the
% result to outMat. Designed to be called once per point from a shell loop so a
% CasADi/IPOPT MEX FATAL crash (uncatchable, kills the process -- and it reliably
% fires on t_f-DECREASING steps) is isolated to a single point instead of taking
% out the whole continuation chain.
%
% INPUTS:
%   seedMat      - .mat with X [8xnN], U [4xnN], factor (the seed's t_f/t_f^min)
%   targetFactor - t_f/t_f^min for this step [scalar]
%   outMat       - output .mat path; on success holds X,U,factor,dV,switches,
%                  edge,defect,primerAlignDeg,lamDef,tf_days,prop_kg,success

here=fileparts(mfilename('fullpath')); addpath(here);
pSund=1.5; sched=[0.05 0.02 0.008 0.003 0.001 0]; maxIter=900;
p=cr3bp_lt_params(0.025,15,2100);
C=load(fullfile(here,'sundman_minfuel_certified.mat'));
sigma=C.sigma; rv0=C.rv0; rvf=C.rvf; tauf0=C.tauf0; tfMin=C.out.X(8,end)/1.15;
S=load(seedMat); X0=S.X; U0=S.U; fSeed=S.factor;

tf=targetFactor*tfMin; tfPrev=fSeed*tfMin; X0(8,:)=X0(8,:)*(tf/tfPrev);
fprintf('TF_STEP: %.3fx -> %.3fx (t_f=%.4f, %.1f d)\n', fSeed, targetFactor, tf, tf*p.tStar/86400);
best=[]; Xk=X0; Uk=U0; o=[];
for ie=1:numel(sched)
    e=sched(ie); tight=ie>1;
    o=casadi_minfuel_sundman(sigma,tf,rv0,rvf,p.Tmax,p.c,p.muStar,Xk,Uk,tauf0,pSund,maxIter,e,tight);
    if o.success && o.maxDefect<1e-6, Xk=o.X; Uk=o.U; best=o; end
end
adv=~isempty(best); if ~adv, best=o; end
dV=p.c*log(1/best.mf)*p.lStar/p.tStar;
X=best.X; U=best.U; factor=targetFactor; %#ok<NASGU>
out=struct('factor',targetFactor,'tf_days',tf*p.tStar/86400,'dV',dV, ...
    'prop_kg',p.m0kg*(1-best.mf),'switches',best.switches,'edge',best.edge, ...
    'defect',best.maxDefect,'primerAlignDeg',best.primerAlignDeg,'success',adv, ...
    'tf',tf,'X',best.X,'U',best.U,'lamDef',best.lamDef); %#ok<NASGU>
save(outMat,'X','U','factor','out');
fprintf('TF_STEP done: dV=%.4f km/s switches=%d edge=%.1f%% defect=%.2g %s -> %s\n', ...
    dV,best.switches,100*best.edge,best.maxDefect,string(adv),outMat);
end
