function test_ifs_eom()
% TEST_IFS_EOM  Hard EOM: matches sms_eom as eps->0, conserves H_sigma.
%
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths();
p = cr3bp_lt_params(0.025, 15, 2100);
pSund = 1.5;
% a representative augmented state near perigee, burn arc (nonzero primer)
Y = [ 0.9; 0.02; -0.01; 0.1; 1.2; -0.3; 0.98; 0.5; ...
      -0.4; 0.2; 0.1;  -0.6; 0.3; 0.05;  -1.3; -0.5 ];

% (a) hard EOM (u=1) matches smoothed sms_eom at tiny eps on a burn state (S<0)
[dY_hard, S] = ifs_eom(0, Y, p.Tmax, p.c, p.muStar, pSund, 1);
assert(S < 0, 'test state should be a burn (S<0), got S=%.3f', S);
dY_sms = sms_eom(0, Y, p.Tmax, p.c, p.muStar, 1e-6, pSund);
assert(max(abs(dY_hard - dY_sms)) < 1e-4, ...
       'hard(u=1) must match sms_eom at eps=1e-6 on a burn state, gap %.2e', ...
       max(abs(dY_hard - dY_sms)));

% (b) coast arc (u=0): thrust/mdot terms vanish
dY_coast = ifs_eom(0, Y, p.Tmax, p.c, p.muStar, pSund, 0);
assert(abs(dY_coast(7)) < 1e-14, 'coast mdot must be 0');

% (c) H_sigma = kappa*(Ht+lamT) conserved along an integrated burn arc
odeOpts = odeset('RelTol',1e-13,'AbsTol',1e-15);
[~, Yout] = ode113(@(s,y) ifs_eom(s,y,p.Tmax,p.c,p.muStar,pSund,1), [0 0.05], Y, odeOpts);
Hs = zeros(size(Yout,1),1);
rE = [-p.muStar;0;0];
for q = 1:size(Yout,1)
    yq = Yout(q,:).';
    r1 = sqrt(sum((yq(1:3)-rE).^2));
    Htq = ifs_Ht(yq, p.Tmax, p.c, p.muStar, 1);
    Hs(q) = r1^pSund * (Htq + yq(16));
end
assert(max(abs(Hs - Hs(1))) < 1e-8, 'H_sigma must be conserved, drift %.2e', max(abs(Hs-Hs(1))));

% (d) complex-step safety: a complex perturbation yields a finite imag derivative
Yc = Y;  Yc(11) = Yc(11) + 1i*1e-20;
dYc = ifs_eom(0, Yc, p.Tmax, p.c, p.muStar, pSund, 1);
assert(all(isfinite(imag(dYc))) && any(imag(dYc)~=0), 'EOM must be complex-step safe');
fprintf('ALL PASS (S=%.3f)\n', S);
end

function Ht = ifs_Ht(Y, Tmax, c, muStar, u)
% local: hard min-fuel time-domain Hamiltonian value for the conservation check
r=Y(1:3); v=Y(4:6); m=Y(7); lamR=Y(9:11); lamV=Y(12:14); lamM=Y(15);
dd=[r(1)+muStar;r(2);r(3)]; rr=[r(1)-1+muStar;r(2);r(3)];
d3=sqrt(sum(dd.^2))^3; r3=sqrt(sum(rr.^2))^3;
gr=[r(1);r(2);0]-(1-muStar)*dd./d3-muStar*rr./r3; hv=[2*v(2);-2*v(1);0];
S = 1 - sqrt(sum(lamV.^2))*c/m - lamM;
Ht = lamR.'*v + lamV.'*(gr+hv) + (Tmax/c)*u*S;
end
