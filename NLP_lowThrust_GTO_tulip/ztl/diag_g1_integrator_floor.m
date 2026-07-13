% DIAG_G1_INTEGRATOR_FLOOR  Is G1's 9.3e-7 terminal miss a cross-integrator
% floor or a ztl_eom/map defect? Integrate the SAME legacy EOM with ode113
% and ode89 (both 1e-13/1e-15) over the same 13-rev arc and compare; then
% compare ztl_flow against the ode89-legacy result (like-for-like).

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
A = load(fullfile(here, 'results', 'p0i_fd_finish.mat'));  a = A.anchor;
[rv0, ~, P0] = ztl_endpoints();
Tmax = 3*P0.Tmax25;  tfL = a.tf;  lamL = a.lam0(:);
y0L = [rv0(:); 1; lamL];

optsI = odeset('RelTol', 1e-13, 'AbsTol', 1e-15);
[~, y113] = ode113(@lt_pmp_eom_energy, [0 tfL], y0L, optsI, Tmax, P0.c, P0.muStar);
[~, y89]  = ode89(@(t,y) lt_pmp_eom_energy(t, y, Tmax, P0.c, P0.muStar), ...
                  [0 tfL], y0L, optsI);
eInt = max(abs(y113(end,1:7) - y89(end,1:7)));
fprintf('legacy EOM, ode113 vs ode89 terminal diff (states 1:7): %.3e\n', eInt);

P = struct('muStar', P0.muStar, 'c', P0.c, 'Tmax', Tmax, 'eps', 1, ...
           'odeRelTol', 1e-13, 'odeAbsTol', 1e-15);
o = ztl_flow([rv0(:); 1; (2*Tmax/P0.c)*lamL], [0 tfL], P, false);
eZtl = max(abs(o.yf(1:7) - y89(end,1:7).'));
fprintf('ztl_flow(eps=1, mapped) vs ode89-legacy terminal diff: %.3e\n', eZtl);
if eZtl < 1e-9
    fprintf('VERDICT: G1 miss was the cross-integrator floor (%.1e); EOM+map exact.\n', eInt);
else
    fprintf('VERDICT: real discrepancy beyond integrator floor -- investigate.\n');
end
