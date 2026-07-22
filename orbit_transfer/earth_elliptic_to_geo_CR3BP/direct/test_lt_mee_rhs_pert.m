% TEST_LT_MEE_RHS_PERT  Opt-in third-body term: back-compat + exact oracle.
here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths;
par = kepler_lt_params(10, 1500, 2000);
Xr = [1.3; 0.2; -0.1; 0.05; 0.02; 0.9; 2.0];   % generic elliptic 3D state
Ur = [0.36; 0.48; 0.80; 0.7];                   % unit beta, mid throttle
par.L = 1.1;
[d0, L0] = lt_mee_rhs(Xr, Ur, par);
% (a) gate 1: pert ABSENT vs gain=0 vs gain=1e-30 -- identical / continuous
parG0 = par;  parG0.pert = lunar_params(par, 0, 0);
[dg0, Lg0] = lt_mee_rhs(Xr, Ur, parG0);
assert(isequal(d0, dg0) && isequal(L0, Lg0), 'gain=0 takes nominal branch: BITWISE identical');
parGt = par;  parGt.pert = lunar_params(par, 0, 1e-12);
[dgt, ~] = lt_mee_rhs(Xr, Ur, parGt);
assert(max(abs(d0-dgt)) < 1e-9, 'cross-branch continuity at tiny gain (FP noise << 1e-9)');
% (b) exact oracle: equatorial circular state at L=0, Moon at phi0=0, t=0.
%     r=[1;0;0], Rhat=[1;0;0], That=[0;1;0], Nhat=[0;0;1]; Moon on +x =>
%     a_M = muM*(1/(DM-1)^2 - 1/DM^2) purely RADIAL (+x). Gauss response:
%     dP/dt += 0; dex/dt += sqrt(P)*sin(L)*aR = 0; dey/dt += -sqrt(P)*cos(L)*aR
%     = -aR; dhx,dhy,Ldot,mdot unchanged.
Xe = [1; 0; 0; 0; 0; 1; 0];  Ue = [0;1;0; 0];   % zero throttle isolates pert
pe = par;  pe.L = 0;  pe.pert = lunar_params(par, 0, 1);
[dp, Lp] = lt_mee_rhs(Xe, Ue, pe);
pe0 = par; pe0.L = 0;
[dq, Lq] = lt_mee_rhs(Xe, Ue, pe0);
aR = pe.pert.muM * (1/(pe.pert.DM-1)^2 - 1/pe.pert.DM^2);
assert(abs(Lp - Lq) < 1e-15, 'radial pert does not touch Ldot');   % check FIRST
delta = (dp - dq) * Lp;                          % back to d/dt (Lp==Lq just proven)
assert(abs(delta(1)) < 1e-14, 'dP unchanged under radial accel');
assert(abs(delta(3) - (-aR)) < 1e-12, 'dey/dt == -aR (radial oracle)');
assert(all(abs(delta([2 4 5 6 7])) < 1e-14), 'ex,hx,hy,m,t untouched');
% (c) frame identities at the generic 3D state (analytic: |Nhat|=1, R.N=0)
%     -- exercised inside the RHS; here we just confirm no NaN and mdot clean
parP = par;  parP.pert = lunar_params(par, 0.7, 1);
[dfull, Lfull] = lt_mee_rhs(Xr, Ur, parP);
assert(all(isfinite(dfull)) && isfinite(Lfull), 'finite with pert on');
assert(abs(dfull(6)*Lfull - d0(6)*L0) < 1e-15, 'PHYSICAL mdot has NO perturbation coupling (d/dL*Ldot roundtrip)');
% (d) transverse oracle: equatorial circular at L=pi/2, Moon on +x, t=0.
%     r=[0;1;0], Rhat=[0;1;0], That=[-1;0;0], Nhat=[0;0;1].
%     a_M = muM*((rM-r)/|rM-r|^3 - rM/DM^3), rM=[DM;0;0] => in-plane only.
%     Gauss response at this state (P=1, Z=1, cL=0, sL=1, A1=0, A2=2):
%     dP/dt += 2*aT; dex/dt += aR; dey/dt += 2*aT; dhx,dhy += 0; Ldot += 0.
Xq = [1; 0; 0; 0; 0; 1; 0];  Uq = [0;1;0; 0];   % zero throttle isolates pert
pq = par;  pq.L = pi/2;  pq.pert = lunar_params(par, 0, 1);
[dqp, Lqp] = lt_mee_rhs(Xq, Uq, pq);
pq0 = par; pq0.L = pi/2;
[dq0, Lq0] = lt_mee_rhs(Xq, Uq, pq0);
DM = pq.pert.DM;  muM = pq.pert.muM;
dvec = [DM; -1; 0];  d3 = (dvec.'*dvec)^1.5;
aM = muM*(dvec/d3 - [DM;0;0]/DM^3);            % inertial accel, z==0
aRq = aM(2);  aTq = -aM(1);                     % Rhat=[0;1;0], That=[-1;0;0]
assert(abs(Lqp - Lq0) < 1e-15, 'in-plane pert (aN=0) does not touch Ldot');
dq = (dqp - dq0) * Lq0;                         % to d/dt (Ldot equal, just proven)
assert(abs(dq(1) - 2*aTq) < 1e-12, 'dP/dt == 2*aT (transverse oracle)');
assert(abs(dq(2) - aRq)   < 1e-12, 'dex/dt == aR at L=pi/2');
assert(abs(dq(3) - 2*aTq) < 1e-12, 'dey/dt == 2*aT at L=pi/2');
assert(all(abs(dq([4 5 6 7])) < 1e-14), 'hx,hy,m,t untouched (aN=0)');
% (e) 3D frame check: closed-form Nhat vs numeric orbit normal r(L) x r(L+dL)
%     at random 3D elliptic states (unit, orthogonal to Rhat, right orientation).
rng(7);
for kk = 1:5
    Pv = 0.5 + rand;  exv = 0.6*(rand-0.5);  eyv = 0.6*(rand-0.5);
    hxv = 0.8*(rand-0.5);  hyv = 0.8*(rand-0.5);  Lv = 2*pi*rand;
    r1 = mee_pos_local(Pv, exv, eyv, hxv, hyv, Lv);
    r2 = mee_pos_local(Pv, exv, eyv, hxv, hyv, Lv + 1e-5);
    nNum = cross(r1, r2);  nNum = nNum / norm(nNum);
    Xhv = 1 + hxv^2 + hyv^2;
    nForm = [2*hyv; -2*hxv; 1 - hxv^2 - hyv^2] / Xhv;
    assert(abs(norm(nForm)-1) < 1e-14, 'Nhat unit');
    assert(abs(dot(nForm, r1/norm(r1))) < 1e-13, 'Nhat orthogonal to Rhat');
    assert(dot(nForm, nNum) > 0.999999, 'Nhat matches numeric orbit normal (orientation)');
end
fprintf('test_lt_mee_rhs_pert: ALL PASS\n');

function r = mee_pos_local(P, ex, ey, hx, hy, L)
% Position from MEE (the same closed forms the RHS pert branch uses).
cL = cos(L); sL = sin(L);
Z = 1 + ex*cL + ey*sL;  rmag = P/Z;
Xh = 1 + hx^2 + hy^2;  a2 = hx^2 - hy^2;
r = [ (rmag/Xh)*(cL + a2*cL + 2*hx*hy*sL);
      (rmag/Xh)*(sL - a2*sL + 2*hx*hy*cL);
      (2*rmag/Xh)*(hx*sL - hy*cL) ];
end
