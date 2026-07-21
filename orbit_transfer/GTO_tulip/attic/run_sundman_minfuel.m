function out = run_sundman_minfuel(N, pSund, tfFactor, maxIter)
% Sundman-regularized min-fuel GTO->tulip test (nominal 25 mN, full spiral).
if nargin<1||isempty(N), N=1500; end
if nargin<2||isempty(pSund), pSund=1.5; end
if nargin<3||isempty(tfFactor), tfFactor=1.15; end
if nargin<4||isempty(maxIter), maxIter=3000; end
here=fileparts(mfilename('fullpath'));
addpath(here); run(fullfile(here,'setup_paths.m'));
addpath(fullfile(here,'..','indirect','lowThrust_GTO_tulip'));
muStar=0.012150585609624; lStar=389703.264829278; tStar=382981.289129055;
m0kg=15; g0=9.80665*tStar^2/(1000*lStar);
Tmax=(0.025/m0kg)*tStar^2/(lStar*1000); c=(2100/tStar)*g0;
muEarth=6.67384e-20*(1-muStar)*(5.9736E24+7.35E22);
sma=(6378+350+6378+35786)/2; ecc=(35786-350)/(2*sma);
[r0,v0]=pumpkyn.cr3bp.orb2eci(muEarth,[sma,ecc,0,-25*pi/180,0,0],2);
rv0=pumpkyn.cr3bp.fromPCI(0,[r0,v0],muStar,tStar,lStar,1);
tfMin=6.2906939607;
zMin=[190.4760481;-79.7060409;-0.4298691037;0.3011592775;0.5866700046;-0.007117348902;4.329378839];
optsInt=odeset('RelTol',1e-12,'AbsTol',1e-14);
tf=tfFactor*tfMin;
% burn+coast warm start (TIME)
[tauB,yB]=ode113(@lt_pmp_eom,[0 tfMin],[rv0(:);1;zMin],optsInt,Tmax,c,muStar); [tauB,k]=unique(tauB,'stable'); yB=yB(k,:);
coastFun=@(t,x) lt_dynamics_dirthrottle(x,[1;0;0;0],Tmax,c,muStar);
[tauC,yC]=ode113(coastFun,[0, tf-tfMin], yB(end,1:7).', optsInt); [tauC,k]=unique(tauC,'stable'); yC=yC(k,:);
rvfT=yC(end,1:6);
tAll=[tauB; tfMin+tauC(2:end)]; xAll=[yB(:,1:7); yC(2:end,:)];
sAll=[ones(numel(tauB),1); zeros(numel(tauC)-1,1)];
lamVB=yB(:,11:13); alB=-lamVB./sqrt(sum(lamVB.^2,2)); alC=repmat([1 0 0],numel(tauC)-1,1);
alAll=[alB;alC];
% --- map time -> Sundman tau: dtau = dt/kappa, kappa=r1^p ---
r1=sqrt((xAll(:,1)+muStar).^2 + xAll(:,2).^2 + xAll(:,3).^2);
kap=r1.^pSund;
dt=diff(tAll); dtau=dt.*0.5.*(1./kap(1:end-1)+1./kap(2:end));
tau=[0; cumsum(dtau)]; tauf0=tau(end); sig_i=tau/tauf0;
[sig_i,ku]=unique(sig_i,'stable'); xAll=xAll(ku,:); tAll=tAll(ku); sAll=sAll(ku); alAll=alAll(ku,:);
sigma=linspace(0,1,N+1).';
Xrvm=interp1(sig_i,xAll,sigma,'pchip');           % r,v,m
tW=interp1(sig_i,tAll,sigma,'pchip');
X0=[Xrvm.'; tW.'];                                 % 8 x (N+1)
sG=min(max(interp1(sig_i,sAll,sigma,'pchip').',0),1);
aG=interp1(sig_i,alAll,sigma,'pchip').'; aG=aG./sqrt(sum(aG.^2,1));
U0=[aG; sG];
X0(1:6,1)=rv0(:); X0(7,1)=1; X0(8,1)=0; X0(1:6,end)=rvfT(:); X0(8,end)=tf;
fprintf('SUNDMAN test: N=%d, pSund=%.2f, tf=%.4f, tauf0=%.4g\n', N, pSund, tf, tauf0);
out=casadi_minfuel_sundman(sigma, tf, rv0, rvfT, Tmax, c, muStar, X0, U0, tauf0, pSund, maxIter);
dV=c*log(1/out.mf)*lStar/tStar;
fprintf('\n=== SUNDMAN RESULT (pSund=%.2f) ===\n', pSund);
fprintf('success=%d status=%s\n', out.success, out.ipoptStatus);
fprintf('maxDefect=%.2g maxUnit=%.2g tauf=%.4g\n', out.maxDefect, out.maxUnit, out.tauf);
fprintf('switches=%d bang-bang=%.1f%% prop=%.4f kg dV=%.4f km/s\n', out.switches, 100*out.edge, m0kg*(1-out.mf), dV);
fprintf('SUNDMAN_PASS=%d\n', out.maxDefect<1e-6);
end
