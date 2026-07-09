function results = run_tf_front(upTo, downTo, step, maxIter)
% RUN_TF_FRONT  Accurate minimum-fuel Delta-V vs transfer-time front by
% SMALL-STEP continuation of the certified good-basin solution.
%
% The energy-continuation sweep (run_tf_sweep) scatters because it re-derives a
% fresh bang-bang basin at each t_f (dense local minima). This driver instead
% starts from the CERTIFIED 25-switch solution (the best-known basin) and
% continues IT in small t_f steps, with a light re-sharpen per step, so the
% solve stays in a good family and Delta-V tracks smoothly. Small steps are
% essential: a large t_f jump breaks the family (a 1.15->1.60 jump wedges).
%
% INPUTS:
%   upTo    - max t_f factor (t_f/t_f^min) to continue UP to   [default 1.85]
%   downTo  - min t_f factor to continue DOWN to               [default 1.05]
%   step    - multiplicative t_f step per continuation move     [default 0.05]
%   maxIter - IPOPT max iters per solve                         [default 900]
%
% OUTPUTS:
%   results - struct array (sorted by factor): .factor .tf_days .dV .prop_kg
%             .switches .edge .defect .primerAlignDeg .success .tf .X .U .lamDef
%             (saved to tf_front_results.mat; figure tf_front.png)

here = fileparts(mfilename('fullpath'));  addpath(here);
if nargin<1||isempty(upTo),   upTo=1.85;  end
if nargin<2||isempty(downTo), downTo=1.05; end
if nargin<3||isempty(step),   step=0.05;  end
if nargin<4||isempty(maxIter),maxIter=900; end
pSund = 1.5;  sched = [0.02 0.006 0.002 0];   % light re-sharpen per step

p = cr3bp_lt_params(0.025, 15, 2100);
C = load(fullfile(here,'sundman_minfuel_certified.mat'));
sigma=C.sigma; rv0=C.rv0; rvf=C.rvf; tauf0=C.tauf0;
Xc=C.out.X; Uc=C.out.U; tfA=Xc(8,end); tfMin=tfA/1.15;

results = struct('factor',{},'tf_days',{},'dV',{},'prop_kg',{},'switches',{}, ...
                 'edge',{},'defect',{},'primerAlignDeg',{},'success',{},'tf',{},'X',{},'U',{},'lamDef',{});
saveP = fullfile(here,'tf_front_results.mat');

% factor grids (anchor 1.15 solved once, in the up pass)
upF   = 1.15:step:upTo;
downF = (1.15-step):-step:downTo;

% ---- up pass (anchor first, continuing certified) ----
Xk=Xc; Uk=Uc; tfPrev=tfA;
for f = upF
    [rec,Xk,Uk,tfPrev] = step_tf(f,tfMin,tfPrev,Xk,Uk,sigma,tauf0,rv0,rvf,pSund,p,sched,maxIter);
    results(end+1)=rec; save(saveP,'results','tfMin'); %#ok<AGROW>
end
% ---- down pass (restart from certified anchor) ----
Xk=Xc; Uk=Uc; tfPrev=tfA;
for f = downF
    [rec,Xk,Uk,tfPrev] = step_tf(f,tfMin,tfPrev,Xk,Uk,sigma,tauf0,rv0,rvf,pSund,p,sched,maxIter);
    results(end+1)=rec; save(saveP,'results','tfMin'); %#ok<AGROW>
end

[~,ord]=sort([results.factor]); results=results(ord); save(saveP,'results','tfMin');
fprintf('\n=== dV-TIME FRONT (factor, days, dV, prop, switches, primer, ok) ===\n');
for k=1:numel(results)
    r=results(k);
    fprintf('  %.2f  %6.2f d  %7.4f km/s  %7.4f kg  %3d sw  %.2f deg  %d\n', ...
        r.factor,r.tf_days,r.dV,r.prop_kg,r.switches,r.primerAlignDeg,r.success);
end
plot_it(results,tfMin,p.tStar,here);
end

% -------------------------------------------------------------------------
function [rec,Xout,Uout,tfOut] = step_tf(f,tfMin,tfPrev,Xk,Uk,sigma,tauf0,rv0,rvf,pSund,p,sched,maxIter)
tf=f*tfMin;  Xk(8,:)=Xk(8,:)*(tf/tfPrev);
fprintf('\n==== t_f=%.3f (%.2fx, %.1f d) ====\n', tf, f, tf*p.tStar/86400);
X0=Xk; U0=Uk; best=[]; o=[];
for ie=1:numel(sched)
    e=sched(ie); tight=ie>1;
    o=casadi_minfuel_sundman(sigma,tf,rv0,rvf,p.Tmax,p.c,p.muStar,X0,U0,tauf0,pSund,maxIter,e,tight);
    if o.success && o.maxDefect<1e-6, X0=o.X; U0=o.U; best=o; end
end
adv=~isempty(best); if ~adv, best=o; end
dV=p.c*log(1/best.mf)*p.lStar/p.tStar;
rec=struct('factor',f,'tf_days',tf*p.tStar/86400,'dV',dV,'prop_kg',p.m0kg*(1-best.mf), ...
    'switches',best.switches,'edge',best.edge,'defect',best.maxDefect, ...
    'primerAlignDeg',best.primerAlignDeg,'success',adv,'tf',tf,'X',best.X,'U',best.U,'lamDef',best.lamDef);
fprintf('  dV=%.4f km/s  switches=%d  edge=%.1f%%  defect=%.2g  %s\n', ...
    dV,best.switches,100*best.edge,best.maxDefect,string(adv));
% advance the continuation only on a clean solve (else keep the prior good state)
if adv, Xout=best.X; Uout=best.U; tfOut=tf; else, Xout=Xk; Uout=Uk; tfOut=tfPrev; end
end

% -------------------------------------------------------------------------
function plot_it(results,tfMin,tStar,here)
ok=[results.success]==1; d=[results.tf_days]; v=[results.dV];
fig=figure('Color','w','Position',[100 100 780 460],'Visible','off');
try
    theme(fig,'light');
catch
end
hold on; grid on; box on;
plot(d(ok),v(ok),'-o','Color',[0.60 0.10 0.10],'MarkerFaceColor',[0.60 0.10 0.10],'LineWidth',1.8);
if any(~ok), plot(d(~ok),v(~ok),'x','Color',[0.5 0.5 0.5],'MarkerSize',9); end
plot(tfMin*tStar/86400,4.4665,'ks','MarkerFaceColor','k','MarkerSize',8);
text(tfMin*tStar/86400,4.4665,'  min-time (4.4665, 0 sw)','FontSize',9,'Color',[0.2 0.2 0.2]);
xlabel('transfer time t_f (days)'); ylabel('\DeltaV (km/s)');
title('Minimum-fuel \DeltaV--time front (certified good-basin continuation)');
exportgraphics(fig,fullfile(here,'tf_front.png'),'Resolution',150); close(fig);
fprintf('WROTE %s\n', fullfile(here,'tf_front.png'));
end
