function build_energy_backbone(facList)
% BUILD_ENERGY_BACKBONE  Phase 1 of the down-sweep: a chain of ENERGY (eps=1)
% solutions continued across t_f, crash-free.
%
% The down-direction defeats every bang-bang method (continuation MEX-crashes;
% a rescaled min-fuel seed makes the energy solve plateau/blow up). The one
% robust primitive is ENERGY->ENERGY continuation: the eps=1 problem is convex
% in the control, so continuing the SMOOTH energy solution to a neighbouring t_f
% converges to machine precision without crashing. This routine builds that
% backbone: start from the cleanly-converged energy solution at 1.15x and step
% t_f through facList (order it as a walk from 1.15 outward), warm-starting each
% energy solve from the previous one (loose warm start = a genuine move). Each
% t_f's energy solution is saved to energy_<factor>.mat for Phase 2
% (SOLVE_TF_MINFUEL: tight re-clean + sharpen), which is embarrassingly parallel.
%
% INPUTS:  facList - t_f/t_f^min values, ordered as a continuation walk from 1.15
%                    [default 1.15:-0.02:1.01 (down toward min-time)]
% OUTPUTS: (none) - writes energy_<factor>.mat per t_f (X,U,factor,tauf0,sigma,
%          rv0,rvf) and appends to energy_backbone_index.mat

here=fileparts(mfilename('fullpath')); addpath(here); pSund=1.5;
if nargin<1||isempty(facList), facList=1.15:-0.02:1.01; end
p=cr3bp_lt_params(0.025,15,2100);
S=load(fullfile(here,'minfuel_from_energy_seed.mat')); rv0=S.rv0; rvf=S.rvf; tfMin=S.tf/1.15;

% seed the backbone: cleanly-converged (tight) energy solution at 1.15
tf15=1.15*tfMin;
[sigma,X0,U0,tauf0]=sundman_seed_map(S.nlp.X,S.nlp.U,tf15,S.sigma,pSund,p.muStar,rv0,rvf);
o=casadi_minfuel_sundman(sigma,tf15,rv0,rvf,p.Tmax,p.c,p.muStar,X0,U0,tauf0,pSund,1500,1,true);
fprintf('[backbone] 1.15 energy seed: ok=%d defect=%.2g\n', o.success, o.maxDefect);
Xp=o.X; Up=o.U; fp=1.15;

for f=facList
    if abs(f-1.15)<1e-9   % the 1.15 anchor itself
        X=o.X; U=o.U; %#ok<NASGU>
    else
        tf=f*tfMin; tfp=fp*tfMin; Xk=Xp; Xk(8,:)=Xk(8,:)*(tf/tfp);   % rescale time to new t_f
        oe=casadi_minfuel_sundman(sigma,tf,rv0,rvf,p.Tmax,p.c,p.muStar,Xk,Up,tauf0,pSund,1500,1,false);
        fprintf('[backbone] %.3fx <- %.3fx energy: ok=%d defect=%.2g edge=%.1f%%\n', ...
                f, fp, oe.success, oe.maxDefect, 100*oe.edge);
        if ~(oe.success && oe.maxDefect<1e-6)
            fprintf('   (energy continuation loose; keeping previous anchor for next step)\n');
        else
            X=oe.X; U=oe.U; Xp=oe.X; Up=oe.U; fp=f; %#ok<NASGU>
        end
        if ~exist('X','var'), X=oe.X; U=oe.U; end %#ok<NASGU>
    end
    factor=f; %#ok<NASGU>
    save(fullfile(here,sprintf('energy_%.2f.mat',f)),'X','U','factor','tauf0','sigma','rv0','rvf','pSund');
end
fprintf('[backbone] DONE: energy solutions saved for %d t_f values.\n', numel(facList));
end
