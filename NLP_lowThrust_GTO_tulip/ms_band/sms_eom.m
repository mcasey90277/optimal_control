function [dY, Ht, S, u] = sms_eom(~, Y, Tmax, c, muStar, epsSmooth, pSund)
% SMS_EOM  Sundman-domain min-fuel PMP dynamics (16-dim augmented system).
%
% Independent variable sigma with dt/dsigma = kappa(r) = r1^pSund,
% r1 = ||r - rE||, rE = [-muStar; 0; 0]. Time t is carried as a state with
% costate lamT. The sigma-domain Hamiltonian is
%     H_sigma = kappa(r) * ( Ht + lamT ),
% where Ht is the TIME-domain smoothed min-fuel Hamiltonian INCLUDING the
% entropy running cost that generates the tanh throttle law:
%     Lear = u*log(u) + (1-u)*log(1-u)
%     Ht   = lamR'*v + lamV'*(g + h) + (Tmax/c)*( u*S + epsSmooth*Lear ).
% The entropy term is required for H_sigma conservation: the smoothed
% throttle u = (1 - tanh(S/(2 eps)))/2 = logistic(-S/eps) is the interior
% argmin of the entropy-smoothed Hamiltonian (envelope theorem then leaves
% the costate equations formally the hard-throttle ones evaluated at u).
% Lear is computed via the CS-safe softplus identity: with z = -S/eps and
% u = logistic(z), the binary entropy obeys H_b(u) = softplus(z) - z*u, so
% Lear = -H_b(u) = z*u - softplus(z); softplus branches on real(z) only
% (branch locally constant, each branch analytic and overflow-free), which
% preserves complex-step derivatives. No abs/max/norm on complex anywhere.
%
% INPUTS:
%   ~         - sigma (unused; autonomous) [scalar]
%   Y         - augmented state [16x1]:
%               [r(1:3); v(4:6); m(7); t(8); lamR(9:11); lamV(12:14);
%                lamM(15); lamT(16)]
%   Tmax      - max thrust acceleration at m = 1 (ND) [scalar]
%   c         - exhaust velocity (ND) [scalar]
%   muStar    - Earth-Moon mass ratio [scalar]
%   epsSmooth - throttle smoothing parameter [scalar]
%   pSund     - Sundman exponent (campaign: 1.5) [scalar]
%
% OUTPUTS:
%   dY - d Y / d sigma [16x1]
%   Ht - time-domain smoothed Hamiltonian value (entropy cost included)
%        [scalar]; H_sigma = kappa*(Ht + lamT)
%   S  - min-fuel switching function 1 - ||lamV||c/m - lamM [scalar]
%   u  - smoothed throttle in (0,1) [scalar]
%
% REFERENCES:
%   [1] Bertrand, Epenoy, OCAM 23(4), 2002 (throttle smoothing).
%   [2] Zhang et al., JGCD 38(8), 2015 (min-fuel CR3BP PMP).
%   [3] .superpowers/sdd/task-S1-brief.md (this exact 16-dim system).

r    = Y(1:3);
v    = Y(4:6);
m    = Y(7);
lamR = Y(9:11);
lamV = Y(12:14);
lamM = Y(15);
lamT = Y(16);

% CR3BP gravity, Coriolis, gravity gradient (same expressions as
% lt_pmp_eom_minfuel; dd = r - rE is the Earth-relative vector)
dd = [r(1) + muStar;     r(2); r(3)];
rr = [r(1) - 1 + muStar; r(2); r(3)];
d1 = sqrt(sum(dd.^2));
d3 = d1^3;
r3 = sqrt(sum(rr.^2))^3;
gr = [r(1); r(2); 0] - (1 - muStar)*dd./d3 - muStar*rr./r3;
hv = [2*v(2); -2*v(1); 0];
d5 = d1^5;
r5 = sqrt(sum(rr.^2))^5;
G  = diag([1, 1, 0]) ...
     - (1 - muStar)*(eye(3)./d3 - 3*(dd*dd.')./d5) ...
     -      muStar *(eye(3)./r3 - 3*(rr*rr.')./r5);
Hc = [0 2 0; -2 0 0; 0 0 0];

% Smoothed min-fuel control law (tanh form: CS-safe, no exp overflow)
lamvMag = sqrt(sum(lamV.^2));
alpha   = -lamV./lamvMag;
S       = 1 - lamvMag*c/m - lamM;
u       = (1 - tanh(S/(2*epsSmooth)))/2;

% Entropy running-cost term via the softplus identity (CS-safe)
z = -S/epsSmooth;
if real(z) > 0
    sp = z + log(1 + exp(-z));
else
    sp = log(1 + exp(z));
end
Lear = z*u - sp;

% Time-domain Hamiltonian value (entropy-smoothed) + Sundman factors
Ht     = lamR.'*v + lamV.'*(gr + hv) + (Tmax/c)*(u*S + epsSmooth*Lear);
kap    = d1^pSund;
dkapdr = pSund * d1^(pSund - 2) * dd;      % gradient of kappa [3x1]
Hval   = Ht + lamT;

dY = [kap*v; ...
      kap*(gr + hv + u*Tmax/m.*alpha); ...
      kap*(-u*Tmax/c); ...
      kap; ...
      kap*(-G.'*lamV) - dkapdr*Hval; ...
      kap*(-lamR - Hc.'*lamV); ...
      kap*(-lamvMag*u*Tmax/m^2); ...
      0];
end
