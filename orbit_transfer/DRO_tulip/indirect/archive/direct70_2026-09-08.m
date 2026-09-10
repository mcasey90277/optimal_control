% 70 mN by DIRECT collocation, seeded with the FULL 75.5 mN trajectory
% (Astra 2026-09-08: initialize with the whole state/control history, a
% modest time stretch and an exactly consistent mass history -- NOT with
% lambda0 and a longer tf).
addpath('/Users/msc/Desktop/optimal_control/orbit_transfer/costate_common');
addpath('/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect');
addpath('/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/direct/lib');
if isempty(which('casadi.Opti')), addpath(fullfile(getenv('HOME'),'casadi-3.7.0')); end
S='/private/tmp/claude-501/-Users-msc-Desktop-optimal-control/0a56dbd5-11b8-4df7-b2f2-e22812d68425/scratchpad/';
A=load([S 'abstract_080.mat']); B=load([S 'abstract_K48.mat']);
st=[A.R.steps B.R.steps]; st=st([st.isp]==900); s0=st(end);      % the 75.5 mN root
Q=load('/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/direct/results/thrust_ladder_12x12.mat');
ob=Q.meta; lStar=ob.lStar; tStar=ob.tStar; muStar=ob.muStar;
g0=9.80665*tStar^2/(1000*lStar); cnd=(900/tStar)*g0;
ndT=@(TN)(TN/ob.m0kg)*tStar^2/(lStar*1000);
[tD,rvD,tT,rvT]=ladder_endpoints(ob);
rv0=interp1(tD,rvD,mod(Q.sD(1),1)*tD(end),'spline')';
rvf=interp1(tT,rvT,mod(Q.sA(11),1)*tT(end),'spline')';

Ttar=0.070; Tnd=ndT(Ttar); N=800;
fprintf('seed: %.1f mN, tf=%.4f ND (%.2f d); target %.1f mN\n', s0.T_N*1000, s0.tf_nd, s0.tf_days, Ttar*1000);

% Fly the 75.5 mN extremal densely, resample to N+1 nodes
[tau,Y]=pumpkyn.cr3bp.tfMinProp(s0.tf_nd,[rv0(1:6);1;s0.z8(1:7)],ndT(s0.T_N),cnd,muStar);
[tu,iu]=unique(tau); Y=Y(iu,:);
sN=linspace(0,1,N+1);
Ys=interp1(tu/tu(end),Y,sN,'pchip')';                 % 14 x (N+1)
tf0=s0.tf_nd*(s0.T_N/Ttar);                            % modest stretch
X0=zeros(7,N+1); X0(1:6,:)=Ys(1:6,:);
X0(7,:)=1-Tnd*(sN*tf0)/cnd;                            % mass consistent at 70 mN
lv=Ys(11:13,:); nlv=sqrt(sum(lv.^2,1));
U0=[-lv./max(nlv,1e-300); ones(1,N+1)];                % primer direction, full throttle
fprintf('tf0 = %.4f ND (%.2f d); N=%d; mf0=%.4f\n', tf0, tf0*tStar/86400, N, X0(7,end));

o=casadi_mintime_dro(rv0(1:6),rvf(1:6),Tnd,cnd,muStar,N,X0,U0,tf0, ...
    struct('scheme','hermite-simpson','maxIter',4000,'printLevel',5,'thrLock',true,'returnModel',true));
fprintf('\nDIRECT 70 mN: success=%d tf=%.6f ND (%.3f d)\n', o.success, o.tf, o.tf*tStar/86400);
save([S 'direct70.mat'],'o','X0','U0','tf0','rv0','rvf','Tnd','cnd','muStar','N');
