function energy_step(seedFile, targetFactor, outFile)
% ENERGY_STEP  One energy->energy continuation step, in an ISOLATED process.
%
% Loads an energy solution (energy_<f>.mat), rescales its time state to
% t_f = targetFactor*t_f^min, and re-solves eps=1 (energy) with the LOOSE warm
% start (a genuine continuation move). Run once per step from a shell loop so a
% sporadic CasADi/IPOPT MEX FATAL crash (uncatchable) is isolated to one step.
% Saves energy_<targetFactor>.mat for solve_tf_minfuel (Phase 2 sharpen).
%
% INPUTS: seedFile - .mat with X,U,factor,tauf0,sigma,rv0,rvf
%         targetFactor - t_f/t_f^min for this step
%         outFile - output .mat

here=fileparts(mfilename('fullpath')); addpath(here); pSund=1.5;
p=cr3bp_lt_params(0.025,15,2100);
Sm=load(fullfile(here,'minfuel_from_energy_seed.mat'),'tf'); tfMin=Sm.tf/1.15;
E=load(seedFile); sigma=E.sigma; rv0=E.rv0; rvf=E.rvf; tauf0=E.tauf0;
tf=targetFactor*tfMin; tfp=E.factor*tfMin;
Xk=E.X; Xk(8,:)=Xk(8,:)*(tf/tfp);
% (a) loose continuation to the new t_f (a genuine move)
o=casadi_minfuel_sundman(sigma,tf,rv0,rvf,p.Tmax,p.c,p.muStar,Xk,E.U,tauf0,pSund,1500,1,false);
fprintf('ENERGY_STEP: %.3f -> %.3f  loose ok=%d defect=%.2g\n', ...
        E.factor, targetFactor, o.success, o.maxDefect);
% (b) TIGHT re-clean at the SAME t_f (no move -> no wedge): keeps the backbone
%     multiplier-consistent so the NEXT continuation starts clean (chaining
%     loose continuations without this amplifies inf_du to ~1e14 and diverges)
if o.success && o.maxDefect<1e-6
    oT=casadi_minfuel_sundman(sigma,tf,rv0,rvf,p.Tmax,p.c,p.muStar,o.X,o.U,tauf0,pSund,1500,1,true);
    fprintf('  re-clean ok=%d defect=%.2g edge=%.1f%%\n', oT.success, oT.maxDefect, 100*oT.edge);
    if oT.success && oT.maxDefect<1e-6
        X=oT.X; U=oT.U; factor=targetFactor; %#ok<NASGU>
        save(outFile,'X','U','factor','tauf0','sigma','rv0','rvf');
        fprintf('  saved (clean) %s\n', outFile);
    else
        fprintf('  re-clean failed -- not saved\n');
    end
else
    fprintf('  loose continuation failed -- not saved; next step continues from last good\n');
end
end
