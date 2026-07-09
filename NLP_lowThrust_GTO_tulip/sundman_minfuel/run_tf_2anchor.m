function run_tf_2anchor(maxIter)
% RUN_TF_2ANCHOR  Two-anchor continuation of the min-fuel dV-time front.
%
% The single certified basin only stays optimal (PMP-certified) over t_f ~
% 1.15-1.25x; past that it drifts. But a t_f-sweep found a DIFFERENT, lower
% basin at 1.75x (23-switch, 2.52 km/s). This driver continues BOTH good
% families in small t_f steps with a finer re-sharpen, so each has its best
% chance to stay a genuine PMP extremal:
%   anchor A = certified 25-switch solution  -> down to 1.00x, up to 1.35x
%   anchor B = the 1.75x 23-switch solution   -> down to 1.30x, up to 1.90x
% The lower envelope of the PMP-certified points from the two threads is the
% front; the key question is whether B's low family extends DOWN into the
% 1.30-1.70x band where A drifts. Results saved per anchor for verify_tf_front.
%
% INPUTS:  maxIter - IPOPT max iters per solve [default 900]
% OUTPUTS: (none) - writes tf_2anchor_A.mat, tf_2anchor_B.mat

here=fileparts(mfilename('fullpath')); addpath(here);
if nargin<1||isempty(maxIter), maxIter=900; end
pSund=1.5;  sched=[0.05 0.02 0.008 0.003 0.001 0];   % finer re-sharpen -> certify
p=cr3bp_lt_params(0.025,15,2100);
C=load(fullfile(here,'sundman_minfuel_certified.mat'));
sigma=C.sigma; rv0=C.rv0; rvf=C.rvf; tauf0=C.tauf0; tfMin=C.out.X(8,end)/1.15;

% anchor A = certified; anchor B = the 1.75x low-basin solution from the sweep
F=load(fullfile(here,'tf_front_results.mat')); Rf=F.results;
[~,jb]=min(abs([Rf.factor]-1.75)); XB=Rf(jb).X; UB=Rf(jb).U; fB=Rf(jb).factor;

% A: down toward min-time, then up through the drift onset
continue_anchor('A', C.out.X, C.out.U, 1.15, [1.10 1.05 1.00 1.20 1.25 1.30 1.35], ...
                sigma,rv0,rvf,tauf0,pSund,p,tfMin,sched,maxIter,fullfile(here,'tf_2anchor_A.mat'));
% B: down into the gap (the interesting direction), then up
continue_anchor('B', XB, UB, fB, [1.70 1.65 1.60 1.55 1.50 1.45 1.40 1.35 1.30 1.80 1.85 1.90], ...
                sigma,rv0,rvf,tauf0,pSund,p,tfMin,sched,maxIter,fullfile(here,'tf_2anchor_B.mat'));
fprintf('\nDONE. Verify with: verify_tf_front(''tf_2anchor_A.mat''), verify_tf_front(''tf_2anchor_B.mat'')\n');
end

% -------------------------------------------------------------------------
function continue_anchor(tag, Xa, Ua, fa, facList, sigma,rv0,rvf,tauf0,pSund,p,tfMin,sched,maxIter,saveFile)
% Continue one bang-bang anchor (Xa,Ua at factor fa) through facList. facList
% is walked in the order given; each entry continues from the nearest already-
% solved factor (so order it outward from the anchor in each direction).
results=struct('factor',{},'tf_days',{},'dV',{},'prop_kg',{},'switches',{},'edge',{}, ...
               'defect',{},'primerAlignDeg',{},'success',{},'tf',{},'X',{},'U',{},'lamDef',{});
solX=containers.Map('KeyType','double','ValueType','any'); solU=containers.Map('KeyType','double','ValueType','any');
solX(fa)=Xa; solU(fa)=Ua;
for f=facList
    known=cell2mat(solX.keys); [~,kk]=min(abs(known-f)); fp=known(kk);
    Xk=solX(fp); Uk=solU(fp); tfPrev=fp*tfMin; tf=f*tfMin;
    Xk(8,:)=Xk(8,:)*(tf/tfPrev);
    fprintf('\n==[%s]== t_f=%.3f (%.2fx, %.1f d)  from %.2fx ==\n', tag,tf,f,tf*p.tStar/86400,fp);
    X0=Xk; U0=Uk; best=[]; o=[];
    for ie=1:numel(sched)
        e=sched(ie); tight=ie>1;
        o=casadi_minfuel_sundman(sigma,tf,rv0,rvf,p.Tmax,p.c,p.muStar,X0,U0,tauf0,pSund,maxIter,e,tight);
        if o.success && o.maxDefect<1e-6, X0=o.X; U0=o.U; best=o; end
    end
    adv=~isempty(best); if ~adv, best=o; end
    dV=p.c*log(1/best.mf)*p.lStar/p.tStar;
    results(end+1)=struct('factor',f,'tf_days',tf*p.tStar/86400,'dV',dV,'prop_kg',p.m0kg*(1-best.mf), ...
        'switches',best.switches,'edge',best.edge,'defect',best.maxDefect, ...
        'primerAlignDeg',best.primerAlignDeg,'success',adv,'tf',tf,'X',best.X,'U',best.U,'lamDef',best.lamDef); %#ok<AGROW>
    fprintf('  dV=%.4f km/s  switches=%d  edge=%.1f%%  defect=%.2g  %s\n', ...
        dV,best.switches,100*best.edge,best.maxDefect,string(adv));
    if adv, solX(f)=best.X; solU(f)=best.U; end
    save(saveFile,'results');    % incremental (survives a MEX crash)
end
end
