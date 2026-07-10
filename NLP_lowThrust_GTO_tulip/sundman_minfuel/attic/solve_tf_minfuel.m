function solve_tf_minfuel(factor, outFile, maxIter)
% SOLVE_TF_MINFUEL  Phase 2 of the down-sweep: sharpen one backbone ENERGY
% solution into a certified min-fuel solution at its t_f.
%
% Consumes energy_<factor>.mat from BUILD_ENERGY_BACKBONE (Phase 1). Each t_f is
% INDEPENDENT here, so many instances run in parallel. The recipe (validated
% 2026-07-09, the first successful down-step in the campaign):
%   1. TIGHT re-clean: re-solve eps=1 (energy) AT this t_f with warmTight=true.
%      The backbone's energy was converged with the LOOSE warm start (needed for
%      the continuation move), which leaves the KKT multipliers inconsistent --
%      sharpening it directly blows up (inf_du ~1e10 on the first eps=0.6 step).
%      A tight re-solve at the SAME t_f (no move, so no wedge) cleans the duals.
%   2. FINE sharpen: energy->fuel homotopy eps 0.6->0.001 (warmTight=true), each
%      step re-solving at a near-bang-bang point.
% Saves ALL data products (state X, control U, discrete costates lamDef, primer,
% transversality) so the result is PMP-verifiable (verify_tf_front.m) as a
% genuine first-order extremal, no re-solve needed. Monotonicity check
% (dV must rise as t_f drops below 1.15x) is done at aggregation.
%
% INPUTS:
%   factor  - t_f/t_f^min (must have a matching energy_<factor>.mat)
%   outFile - output .mat path
%   maxIter - IPOPT max iters per step [default 1500]

here=fileparts(mfilename('fullpath')); addpath(here);
if nargin<3||isempty(maxIter), maxIter=1500; end
pSund=1.5; sched=[0.6 0.35 0.2 0.12 0.07 0.04 0.025 0.015 0.008 0.004 0.002 0.001];
p=cr3bp_lt_params(0.025,15,2100);
E=load(fullfile(here,sprintf('energy_%.2f.mat',factor)));
sigma=E.sigma; rv0=E.rv0; rvf=E.rvf; tauf0=E.tauf0;
Sm=load(fullfile(here,'minfuel_from_energy_seed.mat'),'tf'); tfMin=Sm.tf/1.15; tf=factor*tfMin;

fprintf('SOLVE_TF_MINFUEL: factor=%.3f  t_f=%.4f (%.1f d)\n', factor, tf, tf*p.tStar/86400);
% (1) tight re-clean of the backbone energy (same t_f, so no wedge)
oT=casadi_minfuel_sundman(sigma,tf,rv0,rvf,p.Tmax,p.c,p.muStar,E.X,E.U,tauf0,pSund,maxIter,1,true);
fprintf('  re-clean energy: ok=%d defect=%.2g\n', oT.success, oT.maxDefect);
% (2) fine sharpen
Xk=oT.X; Uk=oT.U; best=oT; nrow=numel(sched); tbl=zeros(nrow,4);
for ie=1:nrow
    e=sched(ie);
    o=casadi_minfuel_sundman(sigma,tf,rv0,rvf,p.Tmax,p.c,p.muStar,Xk,Uk,tauf0,pSund,maxIter,e,true);
    ok=o.success && o.maxDefect<1e-6; tbl(ie,:)=[e,o.maxDefect,o.switches,100*o.edge];
    fprintf('  eps=%.4g: ok=%d defect=%.2g sw=%d\n', e, ok, o.maxDefect, o.switches);
    if ok, Xk=o.X; Uk=o.U; best=o; end
end
dV=p.c*log(1/best.mf)*p.lStar/p.tStar;
out=best; out.factor=factor; out.tf=tf; out.tf_days=tf*p.tStar/86400;
out.dV=dV; out.prop_kg=p.m0kg*(1-best.mf);
save(outFile,'out','sigma','tauf0','rv0','rvf','pSund','factor','tbl');
fprintf('SOLVE_TF_MINFUEL done: factor=%.3f dV=%.4f km/s sw=%d edge=%.1f%% defect=%.2g primer=%.3f -> %s\n', ...
        factor, dV, best.switches, 100*best.edge, best.maxDefect, best.primerAlignDeg, outFile);
end
