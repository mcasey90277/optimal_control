% PROBE_PERIGEE_M  Worst-arc amplification vs node count at p=1 perigee placement.
here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
warning('off', 'MATLAB:ode113:IntegrationTolNotMet');
T = load('results/z0_accept2_trace.mat');
A = load('results/p0i_fd_finish.mat');  a = A.anchor;
[rv0, rvf, P0] = ztl_endpoints();  Tmax = 3*P0.Tmax25;
P = struct('muStar', P0.muStar, 'c', P0.c, 'Tmax', Tmax, 'eps', 1, ...
           'odeRelTol', 1e-13, 'odeAbsTol', 1e-15);
lam0 = T.lam(:);
for M = [104 156 208 312]
    tN = ztl_ms_nodes(lam0, rv0, a.tf, P, M, 'perigee1');
    Y = [rv0(:); 1; lam0];  amp = zeros(1, M);
    for k = 1:M
        o = ztl_flow(Y, [tN(k) tN(k+1)], P, true);  amp(k) = norm(o.PHI);  Y = o.yf;
    end
    fprintf('M=%3d p=1: worst||Phi||=%.2e  max/med=%.0f  unknowns=%d\n', ...
            M, max(amp), max(amp)/median(amp), 14*M-7);
end
