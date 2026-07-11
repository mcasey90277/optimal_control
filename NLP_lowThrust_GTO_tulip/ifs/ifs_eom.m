function [dY, S] = ifs_eom(~, Y, Tmax, c, muStar, pSund, uArc)
% IFS_EOM  Hard-throttle min-fuel PMP dynamics, Sundman tau-domain (16-dim).
%
% Fixed-structure form of the min-fuel PMP EOM: the throttle is a KNOWN
% constant uArc in {0,1} for the whole arc (NOT computed from S), so there is
% no smoothing and no epsilon. Independent variable tau with dt/dtau = kappa =
% r1^pSund; physical time t carried as a state. Identical gravity / costate
% structure to SMS_EOM with the entropy term removed. Complex-step safe.
%
% INPUTS:
%   ~       - tau (unused; autonomous) [scalar]
%   Y       - augmented state [16x1]: [r;v;m;t;lamR;lamV;lamM;lamT]
%   Tmax    - max thrust accel at m=1 (ND) [scalar]
%   c       - exhaust velocity (ND) [scalar]
%   muStar  - Earth-Moon mass ratio [scalar]
%   pSund   - Sundman exponent (1.5) [scalar]
%   uArc    - fixed arc throttle in {0,1} [scalar]
%
% OUTPUTS:
%   dY - dY/dtau [16x1]
%   S  - min-fuel switching function 1 - ||lamV||c/m - lamM [scalar]
%
% REFERENCES:
%   [1] Zhang et al., JGCD 38(8), 2015 (indirect min-fuel CR3BP).
%   [2] docs/superpowers/specs/2026-07-11-ifs-design.md

r = Y(1:3);  v = Y(4:6);  m = Y(7);
lamR = Y(9:11);  lamV = Y(12:14);  lamM = Y(15);  lamT = Y(16);

dd = [r(1)+muStar; r(2); r(3)];
rr = [r(1)-1+muStar; r(2); r(3)];
d1 = sqrt(sum(dd.^2));
d3 = d1^3;  r3 = sqrt(sum(rr.^2))^3;
gr = [r(1); r(2); 0] - (1-muStar)*dd./d3 - muStar*rr./r3;
hv = [2*v(2); -2*v(1); 0];
d5 = d1^5;  r5 = sqrt(sum(rr.^2))^5;
G  = diag([1,1,0]) ...
     - (1-muStar)*(eye(3)./d3 - 3*(dd*dd.')./d5) ...
     -      muStar*(eye(3)./r3 - 3*(rr*rr.')./r5);
Hc = [0 2 0; -2 0 0; 0 0 0];

lamvMag = sqrt(sum(lamV.^2));
S = 1 - lamvMag*c/m - lamM;
% primer direction only enters on a burn arc; guard on the real part for CS
if uArc ~= 0
    if real(lamvMag) < 1e-8
        error('ifs_eom:primerSingular', '||lamV||=%.2e too small on a burn arc', real(lamvMag));
    end
    alpha = -lamV./lamvMag;
else
    alpha = zeros(3,1);
end

kap    = d1^pSund;
dkapdr = pSund*d1^(pSund-2)*dd;
Ht     = lamR.'*v + lamV.'*(gr+hv) + (Tmax/c)*uArc*S;
Hval   = Ht + lamT;

dY = [ kap*v; ...
       kap*(gr + hv + uArc*Tmax/m.*alpha); ...
       kap*(-uArc*Tmax/c); ...
       kap; ...
       kap*(-G.'*lamV) - dkapdr*Hval; ...
       kap*(-lamR - Hc.'*lamV); ...
       kap*(-lamvMag*uArc*Tmax/m^2); ...
       0 ];
end
