% VERIFY_ELFO_SEED  Confirm energy_elfo_freetf.mat is a genuine GTO->ELFO
% min-energy solution: endpoints hit, terminal really at the ELFO, dynamics
% residual small at full gravity + two-primary clock.
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
addpath(fullfile(here,'..','PSR'));
cfg = minfuel_config();  p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
if ~exist('SEEDFILE','var') || isempty(SEEDFILE)
    SEEDFILE = fullfile(here,'results','energy_elfo_freetf.mat');
end
S = load(SEEDFILE);

muGainS = 1;  if isfield(S,'muGain'), muGainS = S.muGain; end
epsS = NaN;   if isfield(S,'epsilon'), epsS = S.epsilon; end
fprintf('seed: X is %dx%d, U %dx%d, moonZone=%.3f muGain=%.2f eps=%.2f\n', ...
        size(S.X,1),size(S.X,2),size(S.U,1),size(S.U,2),S.moonZone,muGainS,epsS);
r0err = norm(S.X(1:6,1)   - S.rv0(:));
rferr = norm(S.X(1:6,end) - S.rvf(:));
mMoon = [1-p.muStar,0,0];
dMoon0 = norm(S.X(1:3,1)  - mMoon.')*p.lStar;
dMoonf = norm(S.X(1:3,end)- mMoon.')*p.lStar;
dEarthf= norm(S.X(1:3,end)- [-p.muStar,0,0].')*p.lStar;
spdf   = norm(S.X(4:6,end));
fprintf('  ||X(:,1)-rv0||   = %.2e   (GTO departure)\n', r0err);
fprintf('  ||X(:,end)-rvf|| = %.2e   (ELFO rendezvous hit)\n', rferr);
fprintf('  terminal: dMoon=%.0f km  dEarth=%.0f km  speed=%.4f ND\n', dMoonf, dEarthf, spdf);
fprintf('  start   : dMoon=%.0f km\n', dMoon0);
fprintf('  tf=%.4f ND (%.2f d)  mf=%.4f (prop %.1f%%)  cScale=%.4f\n', ...
        S.X(8,end), S.X(8,end)*p.tStar/86400, S.X(7,end), 100*(1-S.X(7,end)), S.X(9,end));

% independent (solver-free) defect: recompute the two-primary Sundman dynamics
% and the trapezoidal defect straight from the saved X,U (full gravity, eps=1).
mu = p.muStar;  q = S.qSund;  pw = S.pSund;  D = S.moonZone;
Tmax = p.Tmax;  cEx = p.c;  tauf = S.tauf0;
X = S.X;  U = S.U;  nN = size(X,2);
r = X(1:3,:);  v = X(4:6,:);  m = X(7,:);  cs = X(9,:);
al = U(1:3,:);  s = U(4,:);
dd = [r(1,:)+mu; r(2,:); r(3,:)];     rr = [r(1,:)-1+mu; r(2,:); r(3,:)];
d2 = sum(dd.^2,1)+1e-12;  e2 = sum(rr.^2,1)+1e-12;
r1e = sqrt(d2);  r2e = sqrt(e2);  d3 = d2.^1.5;  r3e = e2.^1.5;
gr = [r(1,:);r(2,:);zeros(1,nN)] - (1-mu)*dd./d3 - mu*rr./r3e;   % muGain=1
hv = [2*v(2,:); -2*v(1,:); zeros(1,nN)];
accel = gr + hv + (s.*Tmax./m).*al;
mdot  = -(Tmax/cEx)*s;
kappa = ( r1e.^(-q) + (r2e/D).^(-q) ).^(-pw/q);
F = [cs.*kappa.*[v; accel; mdot; ones(1,nN)]; zeros(1,nN)];   % 9 x nN
dsig = diff(S.sigma(:)).';
Dd = X(:,2:end) - X(:,1:end-1) - tauf*(repmat(dsig,9,1)/2).*(F(:,1:end-1)+F(:,2:end));
maxDef = max(abs(Dd(:)));
maxUnit = max(abs(sum(al.^2,1)-1));
fprintf('  INDEP defect (pure MATLAB, full gravity) = %.2e   maxUnit=%.2e\n', maxDef, maxUnit);
ok = (r0err<1e-6)&&(rferr<1e-6)&&(maxDef<1e-6)&&(maxUnit<1e-6);
fprintf('VERIFY: %s\n', repmat('PASS', 1, ok));
