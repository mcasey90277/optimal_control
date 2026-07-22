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
fprintf('test_lt_mee_rhs_pert: ALL PASS\n');
